"""Build a target dossier from Open Targets, gnomAD, cBioPortal, and DepMap.

Open Targets supplies identity, disease associations, drugs, tractability
and safety; gnomAD supplies population constraint and loss-of-function
carrier counts; cBioPortal supplies the gene's mutations in tumours
alongside their MSI scores; the MSI dose-response is delegated to
depmap-local, which owns access to the local DepMap release.

The Open Targets table functions are pure: they shape a response fetched
elsewhere (see SKILL.md). The gnomAD and cBioPortal fetch functions do
their own HTTP, following literature-evidence. Every table carries the
underlying ids (Ensembl, ChEMBL, EFO, gnomAD, cBioPortal study/sample)
alongside display names, so a later skill can look up the same records
without re-fetching.
"""

import json
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd


TARGET_URL = "https://platform.opentargets.org/target/{id}"
DRUG_URL = "https://platform.opentargets.org/drug/{id}"
DISEASE_URL = "https://platform.opentargets.org/disease/{id}"

# The datatypeScores id that carries genetic-evidence strength.
GENETIC_DATATYPE = "genetic_association"

GNOMAD_API_URL = "https://gnomad.broadinstitute.org/api"
GNOMAD_GENE_URL = "https://gnomad.broadinstitute.org/gene/{id}?dataset={dataset}"
GNOMAD_VARIANT_URL = "https://gnomad.broadinstitute.org/variant/{id}?dataset={dataset}"

# LOFTEE's high-confidence class. Only these variants are counted as
# loss-of-function; "LC" calls are annotation artefacts often enough that
# including them inflates the carrier count.
GNOMAD_LOF_HIGH_CONFIDENCE = "HC"

CBIOPORTAL_API_URL = "https://www.cbioportal.org/api"
CBIOPORTAL_SAMPLE_URL = ("https://www.cbioportal.org/patient?studyId={study}"
                         "&sampleId={sample}")

# The clinical attributes cBioPortal studies use to record MSI burden, and
# the published call thresholds for each. These are the values to PASS to
# msi_high_samples, not a default it applies on its own: MANTIS >= 0.4
# (Kautto 2017) and MSIsensor >= 10 (Niu 2014) are the authors' cutoffs for
# those specific scores and do not transfer between them.
MSI_ATTRIBUTE_THRESHOLDS = {
    "MSI_SCORE_MANTIS": 0.4,
    "MSI_SENSOR_SCORE": 10.0,
}


def resolve_target(data: dict) -> dict:
    """The target dict, or a clear error if the ensemblId didn't resolve."""
    target = data.get("target")
    if target is None:
        raise ValueError(
            "Open Targets returned no target for this id - check the "
            "Ensembl id resolved correctly before re-running"
        )
    return target


def identity_table(data: dict) -> pd.DataFrame:
    """One row: ensembl id, approved symbol/name, biotype, platform URL."""
    target = resolve_target(data)
    ensembl_id = target.get("id")
    return pd.DataFrame([{
        "ensembl_id": ensembl_id,
        "approved_symbol": target.get("approvedSymbol"),
        "approved_name": target.get("approvedName"),
        "biotype": target.get("biotype"),
        "open_targets_url": TARGET_URL.format(id=ensembl_id) if ensembl_id else None,
    }])


def known_drugs_table(data: dict) -> pd.DataFrame:
    """One row per candidate drug x indication: molecule, phase, mechanism, ids.

    A candidate lists several indications, so a drug appears once per
    disease; rows whose disease is null (unmapped indication) are kept
    with empty disease columns rather than dropped. mechanism_of_action
    is drug-level in Open Targets - it covers every target the molecule
    hits, not only this one - so multiple mechanisms are joined with
    "; " into a single cell.
    """
    columns = ["drug_id", "drug_name", "drug_type", "max_clinical_stage",
               "action_type", "mechanism_of_action",
               "disease_id", "disease_name", "drug_url"]
    rows = ((resolve_target(data).get("drugAndClinicalCandidates") or {}).get("rows") or [])
    out = []
    for row in rows:
        drug = row.get("drug") or {}
        moa_rows = ((drug.get("mechanismsOfAction") or {}).get("rows") or [])
        base = {
            "drug_id": drug.get("id"),
            "drug_name": drug.get("name"),
            "drug_type": drug.get("drugType"),
            "max_clinical_stage": row.get("maxClinicalStage"),
            "action_type": "; ".join(sorted({m["actionType"] for m in moa_rows
                                             if m.get("actionType")})) or None,
            "mechanism_of_action": "; ".join(m["mechanismOfAction"] for m in moa_rows
                                             if m.get("mechanismOfAction")) or None,
            "drug_url": DRUG_URL.format(id=drug["id"]) if drug.get("id") else None,
        }
        for entry in (row.get("diseases") or [{}]):
            disease = entry.get("disease") or {}
            out.append({**base,
                        "disease_id": disease.get("id"),
                        "disease_name": disease.get("name")})
    return pd.DataFrame(out, columns=columns)


def tractability_table(data: dict) -> pd.DataFrame:
    """One row per modality: label, modality code, value. Keeps False rows."""
    columns = ["label", "modality", "value"]
    rows = resolve_target(data).get("tractability") or []
    return pd.DataFrame(rows, columns=columns)


def association_base_rows(data: dict) -> list[dict]:
    """Shared base rows for associations_table and genetic_evidence_table."""
    rows = ((resolve_target(data).get("associatedDiseases") or {}).get("rows") or [])
    out = []
    for row in rows:
        disease = row.get("disease") or {}
        scores = {s["id"]: s.get("score") for s in (row.get("datatypeScores") or [])}
        out.append({
            "disease_id": disease.get("id"),
            "disease_name": disease.get("name"),
            "overall_score": row.get("score"),
            "disease_url": DISEASE_URL.format(id=disease["id"]) if disease.get("id") else None,
            "_datatype_scores": scores,
        })
    return out


def associations_table(data: dict) -> pd.DataFrame:
    """Up to 20 rows, one per disease, plus one column per datatype score
    (wide). Missing datatype for a disease renders as NaN, not 0.
    """
    base = association_base_rows(data)
    datatype_ids = sorted({dt for row in base for dt in row["_datatype_scores"]})
    out = []
    for row in base:
        entry = {k: v for k, v in row.items() if k != "_datatype_scores"}
        for dt in datatype_ids:
            entry[dt] = row["_datatype_scores"].get(dt)
        out.append(entry)
    columns = ["disease_id", "disease_name", "overall_score", "disease_url", *datatype_ids]
    return pd.DataFrame(out, columns=columns)


def genetic_evidence_table(data: dict) -> pd.DataFrame:
    """Disease id/name plus the genetic_association datatype score only."""
    columns = ["disease_id", "disease_name", "genetic_association_score", "disease_url"]
    out = [{
        "disease_id": row["disease_id"],
        "disease_name": row["disease_name"],
        "genetic_association_score": row["_datatype_scores"].get(GENETIC_DATATYPE),
        "disease_url": row["disease_url"],
    } for row in association_base_rows(data)]
    return pd.DataFrame(out, columns=columns)


def safety_table(data: dict) -> pd.DataFrame:
    """One row per reported safety liability. Often empty."""
    columns = ["event", "event_id", "datasource", "literature", "url"]
    rows = resolve_target(data).get("safetyLiabilities") or []
    out = [{
        "event": row.get("event"),
        "event_id": row.get("eventId"),
        "datasource": row.get("datasource"),
        "literature": row.get("literature"),
        "url": row.get("url"),
    } for row in rows]
    return pd.DataFrame(out, columns=columns)



def http_json(url: str, payload: dict | None = None, timeout: int = 60):
    """GET url, or POST payload as JSON to it, and parse the JSON response."""
    if payload is None:
        request = urllib.request.Request(url)
    else:
        request = urllib.request.Request(
            url, data=json.dumps(payload).encode(),
            headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(request, timeout=timeout) as fh:
        return json.load(fh)


# ---------------------------------------------------------------------------
# gnomAD: population constraint and loss-of-function carriers.
#
# This is the number Open Targets' safety field cannot give for a recessive
# gene. LOEUF (oe_lof_upper) says how strongly the population depletes
# loss-of-function alleles; the observed LoF count and the homozygote count
# say how many people already carry them and how many carry two.
# ---------------------------------------------------------------------------

GNOMAD_GENE_QUERY = """
query Gene($symbol: String!, $referenceGenome: ReferenceGenomeId!, $dataset: DatasetId!) {
  gene(gene_symbol: $symbol, reference_genome: $referenceGenome) {
    gene_id
    symbol
    name
    gnomad_constraint {
      exp_lof obs_lof oe_lof oe_lof_lower oe_lof_upper
      exp_mis obs_mis oe_mis lof_z mis_z pLI flags
    }
    variants(dataset: $dataset) {
      variant_id
      consequence
      hgvsp
      lof
      lof_filter
      flags
      exome { ac an homozygote_count }
      genome { ac an homozygote_count }
    }
  }
}
"""


def fetch_gnomad_gene(symbol: str, dataset: str = "gnomad_r4",
                      reference_genome: str = "GRCh38",
                      timeout: int = 120) -> dict:
    """Raw gnomAD gene record: constraint plus every variant in the gene."""
    body = http_json(GNOMAD_API_URL, {
        "query": GNOMAD_GENE_QUERY,
        "variables": {"symbol": symbol, "dataset": dataset,
                      "referenceGenome": reference_genome},
    }, timeout=timeout)
    if body.get("errors"):
        raise ValueError("gnomAD returned errors for %s: %s"
                         % (symbol, body["errors"]))
    gene = (body.get("data") or {}).get("gene")
    if gene is None:
        raise ValueError("gnomAD has no gene %r on %s - check the symbol"
                         % (symbol, reference_genome))
    gene["dataset"] = dataset
    return gene


def gnomad_constraint_table(gene: dict) -> pd.DataFrame:
    """One row of constraint metrics. LOEUF is the oe_lof_upper column.

    Raises when the gene has no constraint record rather than returning an
    empty row: gnomAD omits constraint for genes with too little coverage,
    and a blank row reads as "unconstrained", which is the opposite claim.
    """
    constraint = gene.get("gnomad_constraint")
    if not constraint:
        raise ValueError("gnomAD has no constraint record for %s - it is not "
                         "an unconstrained gene, it is an unmeasured one"
                         % gene.get("symbol"))
    flags = constraint.get("flags") or []
    return pd.DataFrame([{
        "gene_id": gene.get("gene_id"),
        "symbol": gene.get("symbol"),
        "dataset": gene.get("dataset"),
        "loeuf": constraint.get("oe_lof_upper"),
        "oe_lof": constraint.get("oe_lof"),
        "oe_lof_lower": constraint.get("oe_lof_lower"),
        "obs_lof": constraint.get("obs_lof"),
        "exp_lof": constraint.get("exp_lof"),
        "pli": constraint.get("pLI"),
        "lof_z": constraint.get("lof_z"),
        "oe_mis": constraint.get("oe_mis"),
        "mis_z": constraint.get("mis_z"),
        "constraint_flags": "; ".join(flags) or None,
        "gnomad_url": GNOMAD_GENE_URL.format(id=gene.get("gene_id"),
                                             dataset=gene.get("dataset")),
    }])


def gnomad_lof_variants_table(gene: dict) -> pd.DataFrame:
    """One row per high-confidence loss-of-function variant, with carrier counts.

    Exome and genome callsets are separate cohorts in gnomAD, so their allele
    counts are summed and their allele numbers summed alongside - the two
    columns are what the frequency is derived from, and both are kept so a
    reader can see which callset carries the signal.
    """
    columns = ["variant_id", "consequence", "hgvsp", "lof", "lof_filter",
               "flags", "allele_count", "allele_number", "allele_frequency",
               "homozygote_count", "variant_url"]
    rows = []
    for variant in gene.get("variants") or []:
        if variant.get("lof") != GNOMAD_LOF_HIGH_CONFIDENCE:
            continue
        ac = an = hom = 0
        for callset in ("exome", "genome"):
            block = variant.get(callset) or {}
            ac += block.get("ac") or 0
            an += block.get("an") or 0
            hom += block.get("homozygote_count") or 0
        rows.append({
            "variant_id": variant.get("variant_id"),
            "consequence": variant.get("consequence"),
            "hgvsp": variant.get("hgvsp"),
            "lof": variant.get("lof"),
            "lof_filter": variant.get("lof_filter"),
            "flags": "; ".join(variant.get("flags") or []) or None,
            "allele_count": ac,
            "allele_number": an,
            "allele_frequency": (ac / an) if an else None,
            "homozygote_count": hom,
            "variant_url": GNOMAD_VARIANT_URL.format(
                id=variant.get("variant_id"), dataset=gene.get("dataset")),
        })
    out = pd.DataFrame(rows, columns=columns)
    if len(out):
        out = out.sort_values("allele_count", ascending=False).reset_index(drop=True)
    return out


def gnomad_carrier_summary_table(gene: dict) -> pd.DataFrame:
    """One row: how many people carry a loss-of-function allele, and how many two.

    homozygote_count is per variant, so total_homozygotes counts people
    homozygous for one LoF variant and misses compound heterozygotes, who
    carry two different LoF alleles. It is a floor on the two-hit count, not
    an estimate of it.
    """
    variants = gnomad_lof_variants_table(gene)
    constraint = gene.get("gnomad_constraint") or {}
    if variants.empty:
        return pd.DataFrame([{
            "symbol": gene.get("symbol"), "dataset": gene.get("dataset"),
            "n_lof_variants": 0, "total_lof_allele_count": 0,
            "total_homozygotes": 0, "max_allele_frequency": None,
            "loeuf": constraint.get("oe_lof_upper"),
        }])
    return pd.DataFrame([{
        "symbol": gene.get("symbol"),
        "dataset": gene.get("dataset"),
        "n_lof_variants": int(len(variants)),
        "total_lof_allele_count": int(variants.allele_count.sum()),
        "total_homozygotes": int(variants.homozygote_count.sum()),
        "max_allele_frequency": float(variants.allele_frequency.max()),
        "most_common_lof_variant": variants.iloc[0]["variant_id"],
        "loeuf": constraint.get("oe_lof_upper"),
    }])


# ---------------------------------------------------------------------------
# cBioPortal: is the resistance residue already there before treatment?
#
# Hypermutated tumours carry large numbers of passenger mutations, so a
# residue that confers resistance to a target's inhibitor may already be
# present subclonally at baseline. That is checkable: pull the gene's
# mutations from MSI-scored studies, keep the samples called MSI-high, and
# ask which mutations fall in the residue ranges of interest.
# ---------------------------------------------------------------------------


def cbioportal_gene_id(symbol: str, timeout: int = 60) -> int:
    """Entrez gene id for a HUGO symbol."""
    record = http_json(f"{CBIOPORTAL_API_URL}/genes/{urllib.parse.quote(symbol)}",
                       timeout=timeout)
    entrez = record.get("entrezGeneId")
    if entrez is None:
        raise ValueError("cBioPortal has no gene %r" % symbol)
    return int(entrez)


def cbioportal_msi_studies(thresholds: dict, timeout: int = 120) -> pd.DataFrame:
    """Every public study carrying one of the MSI attributes in `thresholds`.

    Discovered from the live attribute list rather than a hardcoded study
    list, so a study added after this was written is still found. A study
    recording two of the attributes appears once per attribute.
    """
    attributes = http_json(f"{CBIOPORTAL_API_URL}/clinical-attributes",
                           timeout=timeout)
    rows = [{"study_id": a["studyId"], "attribute_id": a["clinicalAttributeId"],
             "threshold": thresholds[a["clinicalAttributeId"]]}
            for a in attributes if a.get("clinicalAttributeId") in thresholds]
    return pd.DataFrame(rows, columns=["study_id", "attribute_id", "threshold"]) \
        .sort_values(["study_id", "attribute_id"]).reset_index(drop=True)


def fetch_cbioportal_clinical(study_id: str, attribute_id: str,
                              page_size: int = 2000, max_pages: int = 20,
                              timeout: int = 120) -> list[dict]:
    """Sample-level values of one clinical attribute across a study."""
    out = []
    for page in range(max_pages):
        params = urllib.parse.urlencode({
            "clinicalDataType": "SAMPLE", "attributeId": attribute_id,
            "pageSize": page_size, "pageNumber": page,
        })
        rows = http_json(f"{CBIOPORTAL_API_URL}/studies/"
                         f"{urllib.parse.quote(study_id)}/clinical-data?{params}",
                         timeout=timeout)
        out.extend(rows)
        if len(rows) < page_size:
            break
    return out


def fetch_cbioportal_mutations(study_id: str, entrez_gene_id: int,
                               timeout: int = 120) -> list[dict]:
    """Every mutation in one gene across a study's mutation profile.

    The profile id is the study's own `<study>_mutations`; a study without
    one has no mutation data and yields an empty list rather than an error,
    since the caller sweeps many studies at once.
    """
    profiles = http_json(f"{CBIOPORTAL_API_URL}/studies/"
                         f"{urllib.parse.quote(study_id)}/molecular-profiles",
                         timeout=timeout)
    mutation_profiles = [p["molecularProfileId"] for p in profiles
                         if p.get("molecularAlterationType") == "MUTATION_EXTENDED"]
    if not mutation_profiles:
        return []
    return http_json(f"{CBIOPORTAL_API_URL}/mutations/fetch?projection=DETAILED",
                     {"molecularProfileIds": mutation_profiles,
                      "entrezGeneIds": [entrez_gene_id]},
                     timeout=timeout)


def msi_high_samples(clinical_rows: list[dict], threshold: float) -> dict:
    """{sample_id: score} for samples whose MSI score meets `threshold`.

    Non-numeric values (a study recording "NA" or a text call in a score
    column) are dropped, not coerced.
    """
    out = {}
    for row in clinical_rows:
        try:
            score = float(row.get("value"))
        except (TypeError, ValueError):
            continue
        if score >= threshold:
            out[row["sampleId"]] = score
    return out


def residue_range_hit(protein_position, residue_ranges: dict) -> str | None:
    """Name of the first range in `residue_ranges` containing the position.

    residue_ranges maps a name to an inclusive (start, end) residue pair and
    is always caller-supplied: domain boundaries differ between isoforms and
    annotation sources, so they belong with the report that cites one.
    """
    if protein_position is None:
        return None
    try:
        position = int(protein_position)
    except (TypeError, ValueError):
        return None
    for name, bounds in residue_ranges.items():
        start, end = bounds
        if int(start) <= position <= int(end):
            return name
    return None


def tumour_mutations_table(mutations: list[dict], msi_scores: dict,
                           residue_ranges: dict,
                           attribute_id: str) -> pd.DataFrame:
    """One row per mutation of the gene in an MSI-high sample.

    Rows are restricted to samples in `msi_scores`; a mutation in a sample
    the study never scored is not evidence about MSI-high tumours. VAF comes
    from the read counts when both are present, and is left empty otherwise -
    it is the number that says whether an allele is clonal or a subclonal
    passenger, so an imputed one would be worse than none.
    """
    columns = ["study_id", "sample_id", "patient_id", "protein_change",
               "mutation_type", "protein_pos_start", "residue_range",
               "variant_allele_frequency", "tumor_alt_count", "tumor_ref_count",
               "msi_attribute", "msi_score", "sample_url"]
    rows = []
    for mutation in mutations:
        sample_id = mutation.get("sampleId")
        if sample_id not in msi_scores:
            continue
        alt = mutation.get("tumorAltCount")
        ref = mutation.get("tumorRefCount")
        depth = (alt + ref) if (alt is not None and ref is not None) else None
        rows.append({
            "study_id": mutation.get("studyId"),
            "sample_id": sample_id,
            "patient_id": mutation.get("patientId"),
            "protein_change": mutation.get("proteinChange"),
            "mutation_type": mutation.get("mutationType"),
            "protein_pos_start": mutation.get("proteinPosStart"),
            "residue_range": residue_range_hit(mutation.get("proteinPosStart"),
                                               residue_ranges),
            "variant_allele_frequency": (alt / depth) if depth else None,
            "tumor_alt_count": alt,
            "tumor_ref_count": ref,
            "msi_attribute": attribute_id,
            "msi_score": msi_scores[sample_id],
            "sample_url": CBIOPORTAL_SAMPLE_URL.format(
                study=mutation.get("studyId"), sample=sample_id),
        })
    return pd.DataFrame(rows, columns=columns)


def cbioportal_cohort_table(cohort_rows: list[dict]) -> pd.DataFrame:
    """One row per study: the denominators the mutation table is read against.

    Without n_msi_high_samples, a count of mutated samples cannot be turned
    into a rate, and the mutation table alone invites reading "12 samples
    carry it" as common or rare depending on nothing.
    """
    columns = ["study_id", "attribute_id", "threshold", "n_samples_scored",
               "n_msi_high_samples", "n_msi_high_mutated",
               "n_msi_high_in_range"]
    return pd.DataFrame(cohort_rows, columns=columns)


def fetch_cbioportal_msi_cohort(symbol: str, residue_ranges: dict,
                                thresholds: dict, study_ids: list = None,
                                timeout: int = 120) -> dict:
    """Sweep MSI-scored studies for the gene's mutations in MSI-high samples.

    Returns {"mutations": DataFrame, "cohort": DataFrame}. Pass study_ids to
    restrict the sweep; omitting it covers every public study that records
    one of the attributes in `thresholds`, which is slow but is the only way
    to state a rate over all MSI-high tumours rather than over the studies
    someone thought to name.

    A study that records two MSI attributes is swept once per attribute,
    because the two score different sample sets and neither is a substitute
    for the other. The same mutation therefore appears twice, under
    different msi_attribute and msi_score values. Deduplicate on study_id +
    sample_id + protein_change before counting mutations across attributes.
    """
    if not residue_ranges:
        raise ValueError("residue_ranges is required: give the residue windows "
                         "to test, with the annotation source they came from")
    entrez = cbioportal_gene_id(symbol, timeout=timeout)
    studies = cbioportal_msi_studies(thresholds, timeout=timeout)
    if study_ids is not None:
        studies = studies[studies.study_id.isin(study_ids)]
    mutation_frames = []
    cohort_rows = []
    for study in studies.itertuples():
        clinical = fetch_cbioportal_clinical(study.study_id, study.attribute_id,
                                             timeout=timeout)
        high = msi_high_samples(clinical, study.threshold)
        if not high:
            cohort_rows.append({
                "study_id": study.study_id, "attribute_id": study.attribute_id,
                "threshold": study.threshold, "n_samples_scored": len(clinical),
                "n_msi_high_samples": 0, "n_msi_high_mutated": 0,
                "n_msi_high_in_range": 0,
            })
            continue
        mutations = fetch_cbioportal_mutations(study.study_id, entrez,
                                               timeout=timeout)
        table = tumour_mutations_table(mutations, high, residue_ranges,
                                       study.attribute_id)
        mutation_frames.append(table)
        cohort_rows.append({
            "study_id": study.study_id,
            "attribute_id": study.attribute_id,
            "threshold": study.threshold,
            "n_samples_scored": len(clinical),
            "n_msi_high_samples": len(high),
            "n_msi_high_mutated": int(table.sample_id.nunique()),
            "n_msi_high_in_range": int(
                table[table.residue_range.notna()].sample_id.nunique()),
        })
    mutations_table = (pd.concat(mutation_frames, ignore_index=True)
                       if mutation_frames else tumour_mutations_table([], {}, {}, ""))
    return {"mutations": mutations_table,
            "cohort": cbioportal_cohort_table(cohort_rows)}


# ---------------------------------------------------------------------------
# DepMap: does dependency scale with instability burden, or switch on?
# ---------------------------------------------------------------------------


def msi_dose_response(gene: str, root=None) -> dict:
    """{"stats": DataFrame, "bins": DataFrame} from depmap-local.

    Delegates rather than reimplements: depmap-local owns the local release's
    matrix reads and its MSI score column.
    """
    from lib.common.skills import load_sibling_kernel
    depmap = load_sibling_kernel("depmap-local")
    return {"stats": depmap.depmap_msi_dose_response(gene, root=root),
            "bins": depmap.depmap_msi_dose_bins(gene, root=root)}


def write_run(data: dict, out_dir: Path, subject: str, summary: str,
              steps: list[dict], source: str | None = None,
              limits: list[str] | None = None, gnomad: dict | None = None,
              cbioportal: dict | None = None,
              msi_dose: dict | None = None) -> dict[str, Path]:
    """Write one CSV per table into out_dir, then write README.md.

    The six Open Targets CSVs are always written. gnomad (a fetch_gnomad_gene
    record), cbioportal and msi_dose (the dicts their fetch functions return)
    each add their own CSVs when supplied and are silently absent when not:
    depmap-ot-fusion writes an Open Targets-only sub-run through this same
    function and must not be made to fetch four sources to do it.
    """
    from lib.common.reporting import write_readme
    written = {}
    for name, table_fn in (
        ("identity", identity_table),
        ("known_drugs", known_drugs_table),
        ("tractability", tractability_table),
        ("associations", associations_table),
        ("genetic_evidence", genetic_evidence_table),
        ("safety", safety_table),
    ):
        path = out_dir / f"{name}.csv"
        table_fn(data).to_csv(path, index=False)
        written[name] = path

    tables = {}
    if gnomad is not None:
        tables["gnomad_constraint"] = gnomad_constraint_table(gnomad)
        tables["gnomad_lof_variants"] = gnomad_lof_variants_table(gnomad)
        tables["gnomad_carrier_summary"] = gnomad_carrier_summary_table(gnomad)
    if cbioportal is not None:
        tables["tumour_mutations"] = cbioportal["mutations"]
        tables["tumour_cohort"] = cbioportal["cohort"]
    if msi_dose is not None:
        tables["msi_dose_response"] = msi_dose["stats"]
        tables["msi_dose_bins"] = msi_dose["bins"]
    for name, table in tables.items():
        path = out_dir / f"{name}.csv"
        table.to_csv(path, index=False)
        written[name] = path

    written["readme"] = write_readme(out_dir / "README.md", subject, summary,
                                     steps, source=source, limits=limits)
    return written
