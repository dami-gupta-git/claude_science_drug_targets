"""Build a target dossier from Open Targets, gnomAD, cBioPortal, and DepMap.

Open Targets supplies identity, disease associations, drugs, tractability
and safety; gnomAD supplies population constraint, loss-of-function
carrier counts and common missense variants; cBioPortal supplies the gene's mutations in tumours
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

# Open Targets maxClinicalStage values, highest first, mapped to the tier the
# drug summary reports. APPROVAL and PHASE_4 (a post-marketing trial) share
# the APPROVED tier: both mean the drug is approved somewhere. WITHDRAWAL and
# UNKNOWN are real values that sit off this ladder; any other value raises.
CLINICAL_STAGE_TIERS = {
    "APPROVAL": "APPROVED", "PHASE_4": "APPROVED",
    "PHASE_3": "PHASE_3", "PHASE_2_3": "PHASE_2_3", "PHASE_2": "PHASE_2",
    "PHASE_1_2": "PHASE_1_2", "PHASE_1": "PHASE_1",
}
UNRANKED_CLINICAL_STAGES = ["WITHDRAWAL", "UNKNOWN"]


def markdown_table(rows, headers=None):
    """Markdown table from a list of dicts, resolved lazily.

    Bound at module scope so the published sidecar still imports: the repo
    implementation lives in lib.common.reporting, which only exists inside
    the repo, and importing it at module scope made the whole file
    unloadable elsewhere.
    """
    from lib.common.reporting import markdown_table as _impl
    return _impl(rows, headers)


def write_readme(path, subject, summary, steps, source="", limits=None):
    """Write a run README, resolved lazily. See markdown_table."""
    from lib.common.reporting import write_readme as _impl
    return _impl(path, subject, summary, steps, source=source, limits=limits)


GNOMAD_API_URL = "https://gnomad.broadinstitute.org/api"
GNOMAD_GENE_URL = "https://gnomad.broadinstitute.org/gene/{id}?dataset={dataset}"
GNOMAD_VARIANT_URL = "https://gnomad.broadinstitute.org/variant/{id}?dataset={dataset}"

# LOFTEE's high-confidence class. Only these variants are counted as
# loss-of-function; "LC" calls are annotation artefacts often enough that
# including them inflates the carrier count.
GNOMAD_LOF_HIGH_CONFIDENCE = "HC"

# gnomAD's filter for "no high-quality genotype carries this allele".
GNOMAD_AC0_FILTER = "AC0"

# gnomAD's consequence label for a missense variant.
GNOMAD_MISSENSE_CONSEQUENCE = "missense_variant"

# gnomAD chrom values on which every person carries two copies.
AUTOSOMES = {str(n) for n in range(1, 23)}

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
    "; " into a single cell. parent_drug_id/name are ChEMBL's parent
    molecule, empty for a parent itself.
    """
    columns = ["drug_id", "drug_name", "drug_type", "max_clinical_stage",
               "action_type", "mechanism_of_action",
               "disease_id", "disease_name", "drug_url",
               "parent_drug_id", "parent_drug_name"]
    rows = ((resolve_target(data).get("drugAndClinicalCandidates") or {}).get("rows") or [])
    out = []
    for row in rows:
        drug = row.get("drug") or {}
        parent = drug.get("parentMolecule") or {}
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
            "parent_drug_id": parent.get("id"),
            "parent_drug_name": parent.get("name"),
        }
        for entry in (row.get("diseases") or [{}]):
            disease = entry.get("disease") or {}
            out.append({**base,
                        "disease_id": disease.get("id"),
                        "disease_name": disease.get("name")})
    return pd.DataFrame(out, columns=columns)


def is_salt_form(drug_name, parent_name) -> bool:
    """True when ChEMBL's parent link is a salt or hydrate, not a distinct drug.

    ChEMBL names a salt as its parent plus a counter-ion (TOFACITINIB CITRATE
    -> TOFACITINIB). The same parent link also joins distinct deuterated
    drugs to their non-deuterated parent (DEURUXOLITINIB PHOSPHATE ->
    RUXOLITINIB, DEUCRAVACITINIB -> CRAVACITINIB); those names do not start
    with the parent's, so they stay separate drugs.
    """
    if not drug_name or not parent_name:
        return False
    return drug_name.upper().startswith(parent_name.upper() + " ")


def rank_clinical_stages(stages) -> list[str]:
    """Distinct stages, highest first; WITHDRAWAL/UNKNOWN last. Raises on an unseen value."""
    present = {s for s in stages if isinstance(s, str)}
    unseen = present - set(CLINICAL_STAGE_TIERS) - set(UNRANKED_CLINICAL_STAGES)
    if unseen:
        raise ValueError(f"Unrecognised Open Targets clinical stage(s) {sorted(unseen)}; "
                         "add them to CLINICAL_STAGE_TIERS or UNRANKED_CLINICAL_STAGES")
    return [s for s in [*CLINICAL_STAGE_TIERS, *UNRANKED_CLINICAL_STAGES] if s in present]


def join_distinct(values) -> str | None:
    """Join "; "-separated cells into one, keeping first-seen order without repeats."""
    seen = []
    for cell in values:
        if isinstance(cell, str):
            seen += [p for p in cell.split("; ") if p and p not in seen]
    return "; ".join(seen) or None


def drug_summary_table(data: dict) -> pd.DataFrame:
    """One row per drug, salt forms folded into their parent, highest stage first.

    known_drugs.csv has one row per drug x indication and lists a salt as its
    own drug. This collapses both: rows group on the parent when
    is_salt_form holds, otherwise on the drug itself. max_clinical_stage is
    the tier (CLINICAL_STAGE_TIERS) of the highest stage across the group's
    rows, empty when every row is WITHDRAWAL or UNKNOWN; all_stages lists
    every raw stage seen, so a withdrawal stays visible. forms lists every
    ChEMBL id folded in.
    """
    columns = ["drug_id", "drug_name", "drug_type", "max_clinical_stage",
               "all_stages", "action_type", "mechanism_of_action",
               "n_indications", "forms", "drug_url"]
    drugs = known_drugs_table(data)
    if drugs.empty:
        return pd.DataFrame(columns=columns)
    salt = [is_salt_form(n, p) for n, p in zip(drugs["drug_name"], drugs["parent_drug_name"])]
    drugs["group_id"] = drugs["parent_drug_id"].where(salt, drugs["drug_id"])
    drugs["group_name"] = drugs["parent_drug_name"].where(salt, drugs["drug_name"])

    out = []
    for (group_id, group_name), rows in drugs.groupby(["group_id", "group_name"], sort=False,
                                                     dropna=False):
        stages = rank_clinical_stages(rows["max_clinical_stage"])
        ranked = [s for s in stages if s in CLINICAL_STAGE_TIERS]
        out.append({
            "drug_id": group_id,
            "drug_name": group_name,
            "drug_type": join_distinct(rows["drug_type"]),
            "max_clinical_stage": CLINICAL_STAGE_TIERS[ranked[0]] if ranked else None,
            "all_stages": "; ".join(stages) or None,
            "action_type": join_distinct(rows["action_type"]),
            "mechanism_of_action": join_distinct(rows["mechanism_of_action"]),
            "n_indications": rows["disease_id"].dropna().nunique(),
            "forms": join_distinct(rows["drug_id"]),
            "drug_url": DRUG_URL.format(id=group_id),
        })
    rank = {tier: i for i, tier in enumerate(dict.fromkeys(CLINICAL_STAGE_TIERS.values()))}
    out.sort(key=lambda r: (rank.get(r["max_clinical_stage"], len(rank)), r["drug_name"]))
    return pd.DataFrame(out, columns=columns)


def tractability_table(data: dict) -> pd.DataFrame:
    """One row per modality: label, modality code, value. Keeps False rows."""
    columns = ["label", "modality", "value"]
    rows = resolve_target(data).get("tractability") or []
    return pd.DataFrame(rows, columns=columns)


def association_base_rows(data: dict) -> list[dict]:
    """Shared base rows for associations_table and top_association_genetic_scores_table."""
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


def association_count(data: dict) -> dict:
    """How many disease associations exist versus how many were returned.

    The query pages `associatedDiseases` at size 20, so a gene with a thousand
    associations yields twenty rows and no indication in the tables that the
    rest exist. Open Targets reports the true total in `count`; this surfaces
    it so a run can state the cap rather than presenting the top slice as the
    whole picture. Ordering is by descending overall score, so the returned
    rows are the strongest - the cap loses the tail, not a random sample.
    """
    block = ((data.get("target") or {}).get("associatedDiseases") or {})
    rows = block.get("rows") or []
    total = block.get("count")
    return {
        "n_returned": len(rows),
        "n_total": total,
        "truncated": bool(total is not None and total > len(rows)),
        "lowest_returned_score": min((r.get("score") for r in rows
                                      if r.get("score") is not None),
                                     default=None),
    }


def associations_table(data: dict) -> pd.DataFrame:
    """One row per returned disease, plus one column per datatype score (wide).

    Capped at the query's page size, currently 20, so this is the strongest
    associations and not all of them. Call `association_count` for the true
    total and quote it alongside; a table of 20 rows for a gene with a
    thousand associations is a top slice, not a survey. Missing datatype for a
    disease renders as NaN, not 0.
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


def top_association_genetic_scores_table(data: dict) -> pd.DataFrame:
    """Genetic_association score for the top associations that have one.

    Drawn from the same capped, overall-score-ranked slice as
    associations_table, so it is not a ranking of the gene's genetic
    evidence: a disease with strong genetic support but a lower overall
    score is absent. Diseases in the slice with no genetic score are dropped.
    """
    columns = ["disease_id", "disease_name", "genetic_association_score", "disease_url"]
    out = [{
        "disease_id": row["disease_id"],
        "disease_name": row["disease_name"],
        "genetic_association_score": row["_datatype_scores"][GENETIC_DATATYPE],
        "disease_url": row["disease_url"],
    } for row in association_base_rows(data)
        if row["_datatype_scores"].get(GENETIC_DATATYPE) is not None]
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
    chrom
    canonical_transcript_id
    gnomad_constraint {
      exp_lof obs_lof oe_lof oe_lof_lower oe_lof_upper
      exp_mis obs_mis oe_mis lof_z mis_z pLI flags
    }
    variants(dataset: $dataset) {
      variant_id
      consequence
      transcript_id
      hgvsp
      lof
      lof_filter
      flags
      exome { ac an homozygote_count filters }
      genome { ac an homozygote_count filters }
      joint { ac an homozygote_count filters }
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
    # Recorded on the record, not just requested, so the downstream tables can
    # state it. Variant ids are build-specific coordinates: 8-31081132-C-T means
    # nothing without knowing which assembly it is on, and a reader who assumes
    # the wrong one looks up a different position.
    gene["reference_genome"] = reference_genome
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
        "reference_genome": gene.get("reference_genome"),
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


GNOMAD_VARIANT_COLUMNS = ["variant_id", "reference_genome", "consequence", "hgvsp",
                          "lof", "lof_filter", "flags", "callsets_counted",
                          "allele_count", "allele_number", "allele_frequency",
                          "homozygote_count", "variant_url"]


def passes_gnomad_qc(block: dict) -> bool:
    """True when a callset block carries no gnomAD quality filter.

    AC0 is not a quality failure: it says no high-quality genotype carries the
    allele, so the block still contributes the people called at the site.
    """
    return not set(block.get("filters") or []) - {GNOMAD_AC0_FILTER}


def gnomad_variant_counts(variant: dict) -> tuple:
    """(ac, an, homozygotes, callsets) counted only from callsets that passed QC.

    When every callset passes, gnomAD's joint block is used: it keeps the
    people called in a callset where the variant was not seen. The joint block
    does not carry exome or genome quality failures (TYK2 19-10352476-A-C fails
    exome VQSR with a clean joint block), so when a callset fails, only the
    passing callsets are summed and the failing one is dropped from both the
    count and the denominator. No passing callset gives all-None counts.
    """
    blocks = {c: variant[c] for c in ("exome", "genome") if variant.get(c)}
    joint = variant.get("joint")
    if joint and passes_gnomad_qc(joint) and all(map(passes_gnomad_qc, blocks.values())):
        return joint.get("ac"), joint.get("an"), joint.get("homozygote_count"), "joint"
    passing = {c: b for c, b in blocks.items() if passes_gnomad_qc(b)}
    if not passing:
        return None, None, None, None
    total = lambda key: sum(b.get(key) or 0 for b in passing.values())
    return total("ac"), total("an"), total("homozygote_count"), "; ".join(passing)


def gnomad_variant_row(variant: dict, gene: dict) -> dict:
    """One variant with its QC-passing carrier counts; see gnomad_variant_counts."""
    ac, an, hom, callsets = gnomad_variant_counts(variant)
    return {
        "variant_id": variant.get("variant_id"),
        # variant_id is a build-specific coordinate; carried per row so a
        # reader who opens only this file can still resolve the position.
        "reference_genome": gene.get("reference_genome"),
        "consequence": variant.get("consequence"),
        "hgvsp": variant.get("hgvsp"),
        "lof": variant.get("lof"),
        "lof_filter": variant.get("lof_filter"),
        "flags": "; ".join(variant.get("flags") or []) or None,
        "callsets_counted": callsets,
        "allele_count": ac,
        "allele_number": an,
        "allele_frequency": (ac / an) if an and ac is not None else None,
        "homozygote_count": hom,
        "variant_url": GNOMAD_VARIANT_URL.format(
            id=variant.get("variant_id"), dataset=gene.get("dataset")),
    }


def gnomad_carried_variant_rows(gene: dict, keep) -> list[dict]:
    """Rows for variants where keep(variant) holds and a QC-passing carrier exists.

    gnomAD lists variants whose every call failed QC or that no high-quality
    genotype carries (filter AC0); counting them inflates a variant count
    with alleles nobody is observed to carry.
    """
    rows = [gnomad_variant_row(v, gene) for v in gene.get("variants") or [] if keep(v)]
    return [r for r in rows if r["allele_count"]]


def gnomad_lof_variants_table(gene: dict) -> pd.DataFrame:
    """One row per carried high-confidence loss-of-function variant, most carriers first."""
    rows = gnomad_carried_variant_rows(
        gene, lambda v: v.get("lof") == GNOMAD_LOF_HIGH_CONFIDENCE)
    out = pd.DataFrame(rows, columns=GNOMAD_VARIANT_COLUMNS)
    if len(out):
        out = out.sort_values("allele_count", ascending=False).reset_index(drop=True)
    return out


def minor_allele_carriers(row: dict, chrom: str) -> int | None:
    """People carrying at least one copy of the minor allele; None off the autosomes.

    On an autosome everyone called at the site has two copies, so an_people is
    an/2 and an alt-allele carrier count is ac - homozygotes (each homozygote
    holds two of the ac copies). Where the reference allele is the minor one,
    its carriers are everyone not homozygous for the alt. On X, Y and the
    mitochondrial genome copy number differs by sex or is not two, so the
    count is left empty rather than computed on the two-copy assumption.
    """
    if chrom not in AUTOSOMES:
        return None
    ac, an, hom = row["allele_count"], row["allele_number"], row["homozygote_count"]
    if ac is None or not an or hom is None:
        return None
    return an // 2 - hom if row["reference_is_minor"] else ac - hom


def gnomad_missense_variants_table(gene: dict, min_allele_frequency: float) -> pd.DataFrame:
    """Missense variants in the canonical transcript whose minor allele frequency
    is at least the cutoff.

    min_allele_frequency is caller-supplied with no default. A common missense
    variant is a candidate natural experiment (TYK2 P1104A), but gnomAD says
    only that it is common, not that it changes protein function.

    gnomAD reports each variant's consequence on one transcript; a row is kept
    only when that is the gene's canonical transcript, so hgvsp is in canonical
    numbering and a change to a minor isoform is not listed as a change to the
    protein. gnomAD's consequence_in_canonical_transcript flag is null on this
    endpoint, so transcript_id is compared instead.

    The cutoff applies to the minor allele, min(AF, 1 - AF). Where the
    reference genome carries the rarer allele (PELO L221M, AF 0.91), the
    gnomAD "variant" is the population's usual form; reference_is_minor marks
    these rows, and their homozygote_count counts people with the usual form.
    minor_allele_carrier_fraction is the share of people, not of gene copies;
    see minor_allele_carriers. A variant whose allele number is zero has no
    frequency and is excluded.
    """
    if not 0 < min_allele_frequency <= 0.5:
        raise ValueError("min_allele_frequency is a minor allele frequency and "
                         "must be in (0, 0.5], got %r" % min_allele_frequency)
    canonical = gene.get("canonical_transcript_id")
    if not canonical:
        raise ValueError("gnomAD gave no canonical transcript for %s"
                         % gene.get("symbol"))
    columns = GNOMAD_VARIANT_COLUMNS[:-1] + [
        "transcript_id", "minor_allele_frequency", "reference_is_minor",
        "minor_allele_carriers", "minor_allele_carrier_fraction", "variant_url"]
    rows = []
    for row in gnomad_carried_variant_rows(
            gene, lambda v: v.get("consequence") == GNOMAD_MISSENSE_CONSEQUENCE
            and v.get("transcript_id") == canonical):
        af = row["allele_frequency"]
        if af is None or min(af, 1 - af) < min_allele_frequency:
            continue
        row["transcript_id"] = canonical
        row["minor_allele_frequency"] = min(af, 1 - af)
        row["reference_is_minor"] = af > 0.5
        carriers = minor_allele_carriers(row, gene.get("chrom"))
        row["minor_allele_carriers"] = carriers
        row["minor_allele_carrier_fraction"] = (
            carriers / (row["allele_number"] // 2) if carriers is not None else None)
        rows.append(row)
    out = pd.DataFrame(rows, columns=columns)
    if len(out):
        out = out.sort_values("minor_allele_frequency", ascending=False).reset_index(drop=True)
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
        "reference_genome": gene.get("reference_genome"),
            "n_lof_variants": 0, "total_lof_allele_count": 0,
            "total_homozygotes": 0, "max_allele_frequency": None,
            "loeuf": constraint.get("oe_lof_upper"),
        }])
    return pd.DataFrame([{
        "symbol": gene.get("symbol"),
        "dataset": gene.get("dataset"),
        "reference_genome": gene.get("reference_genome"),
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
    """Entrez gene id for a HUGO symbol.

    The symbol is upper-cased before the lookup. cBioPortal indexes HUGO
    symbols in upper case and returns 404 for any other casing, so a gene
    whose approved symbol carries lower-case letters - C9orf72, C4orf46, the
    whole `orf` family - fails on the symbol the rest of the pipeline uses.
    Open Targets, gnomAD and DepMap all accept the mixed-case form, so the
    casing difference surfaces only here.
    """
    lookup = (symbol or "").strip().upper()
    if not lookup:
        raise ValueError("symbol is required")
    try:
        record = http_json(
            f"{CBIOPORTAL_API_URL}/genes/{urllib.parse.quote(lookup)}",
            timeout=timeout)
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            raise ValueError(
                "cBioPortal has no gene %r (looked up as %r). Check the "
                "approved HUGO symbol; an alias or a withdrawn symbol will "
                "404 here while still resolving in Open Targets"
                % (symbol, lookup)) from exc
        raise
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


def clonality_floor(mutations_table) -> dict:
    """The VAF distribution of what was called, which bounds an absence claim.

    A finding that residue X is unmutated across a bulk-sequenced cohort is
    not evidence that no cell carries it. Bulk calls have a detection floor:
    below some allele fraction a variant is not called at all, and
    pre-existing resistance is characteristically subclonal, so the honest
    claim is that no such allele is clonal or near-clonal - not that none
    exists. Sample size is the looser bound and the one people reach for
    first; the floor is usually tighter.

    Returns the observed floor and quartiles over mutations that carry a VAF,
    for quoting in the run's limits. `min_vaf` is the lowest fraction this
    cohort actually called, which is an upper bound on the true floor rather
    than a measurement of it - the pipeline's real limit is lower and
    unstated. Deduplicate the table before calling: a study scored on two MSI
    attributes lists each mutation twice.
    """
    values = pd.to_numeric(mutations_table.get("variant_allele_frequency"),
                           errors="coerce").dropna()
    if values.empty:
        return {"n_with_vaf": 0, "min_vaf": None, "p5_vaf": None,
                "median_vaf": None, "frac_below_0.10": None}
    return {
        "n_with_vaf": int(len(values)),
        "min_vaf": round(float(values.min()), 4),
        "p5_vaf": round(float(values.quantile(0.05)), 3),
        "median_vaf": round(float(values.median()), 3),
        "frac_below_0.10": round(float((values < 0.10).mean()), 3),
    }


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
                                test_positions: dict = None,
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

    Pass `test_positions` - a {label: residue_position} map of the specific
    residues the run is asking about - whenever the question is "is residue X
    mutated here". It is checked against `residue_ranges` and raises if any
    position falls outside every range, because that combination silently
    produces the wrong answer: the mutation table labels such a mutation as
    outside the named ranges, the in-range count omits it, and the run reads
    as evidence of absence when the position was never queried. Choosing
    `residue_ranges` from a protein's annotated domains is the easy way to
    hit this - annotated domains do not tile the sequence, so a residue in an
    unannotated stretch between two domains is invisible.
    """
    if not residue_ranges:
        raise ValueError("residue_ranges is required: give the residue windows "
                         "to test, with the annotation source they came from")
    if test_positions:
        uncovered = {}
        for label, position in test_positions.items():
            if residue_range_hit(position, residue_ranges) is None:
                uncovered[label] = position
        if uncovered:
            raise ValueError(
                "residue_ranges does not cover %d tested position(s): %s. A "
                "mutation there would be reported as outside the named ranges "
                "and omitted from n_msi_high_in_range, so the run would read as "
                "absence when the position was never queried. Add a window "
                "around each, or pass ranges that tile the whole sequence."
                % (len(uncovered), ", ".join("%s=%s" % kv for kv in
                                             sorted(uncovered.items()))))
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
    # Empty frames are dropped before concatenation rather than passed to it. A
    # study with no mutations in MSI-high samples contributes an all-NA frame,
    # and pandas warns that such frames will stop being excluded from dtype
    # inference in a future version - at which point numeric columns here would
    # silently become object, and variant_allele_frequency would stop comparing
    # as a number. Dropping them keeps the dtypes fixed either way.
    populated = [frame for frame in mutation_frames if not frame.empty]
    mutations_table = (pd.concat(populated, ignore_index=True)
                       if populated else tumour_mutations_table([], {}, {}, ""))
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
              missense_min_allele_frequency: float | None = None,
              cbioportal: dict | None = None,
              msi_dose: dict | None = None) -> dict[str, Path]:
    """Write one CSV per table into out_dir, then write README.md.

    The seven Open Targets CSVs are always written. gnomad (a fetch_gnomad_gene
    record), cbioportal and msi_dose (the dicts their fetch functions return)
    each add their own CSVs when supplied and are silently absent when not:
    depmap-ot-fusion writes an Open Targets-only sub-run through this same
    function and must not be made to fetch four sources to do it. gnomad
    requires missense_min_allele_frequency, so the missense table cannot be
    skipped by omission.
    """
    if gnomad is not None and missense_min_allele_frequency is None:
        raise ValueError("gnomad was supplied without missense_min_allele_frequency")
    written = {}
    for name, table_fn in (
        ("identity", identity_table),
        ("known_drugs", known_drugs_table),
        ("drug_summary", drug_summary_table),
        ("tractability", tractability_table),
        ("associations", associations_table),
        ("top_association_genetic_scores", top_association_genetic_scores_table),
        ("safety", safety_table),
    ):
        path = out_dir / f"{name}.csv"
        table_fn(data).to_csv(path, index=False)
        written[name] = path

    # Written unconditionally: it records that associations.csv is a capped top
    # slice and how much was left behind. A run that omitted this read as a
    # complete survey of a gene's diseases when it held 20 of a thousand.
    coverage = out_dir / "association_coverage.csv"
    pd.DataFrame([association_count(data)]).to_csv(coverage, index=False)
    written["association_coverage"] = coverage

    tables = {}
    if gnomad is not None:
        tables["gnomad_constraint"] = gnomad_constraint_table(gnomad)
        tables["gnomad_lof_variants"] = gnomad_lof_variants_table(gnomad)
        tables["gnomad_carrier_summary"] = gnomad_carrier_summary_table(gnomad)
        tables["gnomad_missense_variants"] = gnomad_missense_variants_table(
            gnomad, missense_min_allele_frequency)
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
