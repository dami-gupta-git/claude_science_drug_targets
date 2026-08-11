"""Fuse DepMap dependency data and Open Targets evidence into a per-gene dossier."""

import json
import math
import re
from datetime import datetime, timedelta
from pathlib import Path

from lib.common.paths import project_root, results_dir
from lib.common.reporting import write_readme

OT_GRAPHQL_URL = "https://api.platform.opentargets.org/api/v4/graphql"
OT_DEFAULT_MAX_AGE_DAYS = 90

DEPMAP_REQUIRED = {"selectivity.csv", "lineage_enrichment.csv"}
OT_REQUIRED = {
    "identity.csv", "associations.csv", "tractability.csv",
    "known_drugs.csv", "safety.csv", "genetic_evidence.csv",
}
CONTRAST_FILES = [
    "msi_contrast.csv", "msi_lineage_contrast.csv",
    "mutation_contrast_overall.csv", "mutation_contrast_by_lineage.csv",
    "fusion_contrast.csv",
]

# Maps each DepMap OncotreeLineage to keywords found in Open Targets disease names.
# Longer keywords match preferentially over shorter ones (see lineage_for_disease).
LINEAGE_KEYWORDS = {
    "Bowel": ["colorectal carcinoma", "colon carcinoma", "rectal carcinoma",
              "colorectal", "colon cancer", "bowel", "intestinal"],
    "Lung": ["lung carcinoma", "non-small cell lung", "small cell lung",
             "lung adenocarcinoma", "lung squamous", "lung cancer"],
    "Breast": ["breast carcinoma", "breast cancer", "breast neoplasm"],
    "Skin": ["cutaneous melanoma", "skin carcinoma", "melanoma", "skin cancer"],
    "Eye": ["uveal melanoma", "retinoblastoma", "ocular melanoma"],
    "Pancreas": ["pancreatic carcinoma", "pancreatic ductal", "pancreatic cancer"],
    "Liver": ["hepatocellular carcinoma", "liver carcinoma", "hepatic carcinoma"],
    "Ovary/Fallopian Tube": ["ovarian carcinoma", "fallopian tube", "ovarian cancer"],
    "Bladder/Urinary Tract": ["bladder carcinoma", "urothelial carcinoma", "bladder cancer"],
    "Kidney": ["renal cell carcinoma", "kidney carcinoma", "renal carcinoma"],
    "Prostate": ["prostate carcinoma", "prostate cancer"],
    "Head and Neck": ["head and neck carcinoma", "oral cavity", "pharyngeal carcinoma",
                      "laryngeal carcinoma", "head and neck squamous"],
    "Esophagus/Stomach": ["esophageal carcinoma", "gastric carcinoma",
                          "stomach carcinoma", "esophag"],
    "CNS/Brain": ["glioblastoma", "glioma", "brain carcinoma", "brain cancer"],
    "Myeloid": ["acute myeloid leukemia", "myeloid leukemia", "multiple myeloma",
                "myeloid malignancy"],
    "Lymphoid": ["lymphocytic leukemia", "lymphoid leukemia", "lymphoma",
                 "hodgkin", "burkitt"],
    "Uterus": ["uterine carcinoma", "endometrial carcinoma", "endometrium"],
    "Cervix": ["cervical carcinoma", "cervix uteri"],
    "Thyroid": ["thyroid gland carcinoma", "thyroid carcinoma", "thyroid cancer"],
    "Biliary Tract": ["cholangiocarcinoma", "bile duct carcinoma",
                      "gallbladder carcinoma", "biliary"],
    "Bone": ["osteosarcoma", "ewing sarcoma", "bone carcinoma"],
    "Soft Tissue": ["soft tissue sarcoma", "leiomyosarcoma", "liposarcoma", "sarcoma"],
    "Pleura": ["mesothelioma", "pleural mesothelioma"],
    "Adrenal Gland": ["adrenal carcinoma", "pheochromocytoma", "adrenocortical"],
    "Peripheral Nervous System": ["neuroblastoma"],
    "Testis": ["testicular germ cell", "testicular carcinoma"],
    "Vulva/Vagina": ["vulvar carcinoma", "vaginal carcinoma"],
    "Ampulla of Vater": ["ampullary carcinoma", "ampulla of vater"],
}

# A literal tuple rather than frozenset(...): the kernel.py sidecar loader
# only permits literal constant assignment at module scope, and a function
# call there (even one wrapping literals) is rejected, defining none of this
# module's helpers. Membership tests (`in`) work the same on a tuple.
CANCER_TERMS = (
    "cancer", "carcinoma", "tumor", "tumour", "sarcoma", "lymphoma",
    "leukemia", "leukaemia", "melanoma", "glioma", "myeloma", "blastoma",
    "adenocarcinoma", "mesothelioma", "neuroblastoma", "malignant",
    "neoplasm", "neoplasia",
)

OT_SEARCH_QUERY = """
query Search($q: String!) {
  search(queryString: $q, entityNames: ["target"], page: {index: 0, size: 5}) {
    hits { id name }
  }
}
"""


# ---------------------------------------------------------------------------
# Readers
# ---------------------------------------------------------------------------

def read_depmap_run(run_dir):
    """Read all DepMap result CSVs from run_dir.

    Returns dict with 'frames' (name -> DataFrame) and 'optional_files' (list
    of optional file names that were present).
    """
    import pandas as pd
    run_dir = Path(run_dir)
    frames = {
        "selectivity": pd.read_csv(run_dir / "selectivity.csv"),
        "lineage_enrichment": pd.read_csv(run_dir / "lineage_enrichment.csv"),
    }
    optional = []
    for name in CONTRAST_FILES + ["prism_compounds.csv"]:
        p = run_dir / name
        if p.exists():
            frames[name.replace(".csv", "")] = pd.read_csv(p)
            optional.append(name)
    return {"frames": frames, "optional_files": optional}


def read_open_targets_run(run_dir):
    """Read all Open Targets result CSVs from run_dir.

    Returns dict with 'frames' (name -> DataFrame).
    """
    import pandas as pd
    run_dir = Path(run_dir)
    frames = {}
    for name in ["identity", "associations", "tractability",
                 "known_drugs", "safety", "genetic_evidence"]:
        p = run_dir / f"{name}.csv"
        if p.exists():
            frames[name] = pd.read_csv(p)
        else:
            frames[name] = pd.DataFrame()
    return {"frames": frames}


# ---------------------------------------------------------------------------
# Completeness and staleness
# ---------------------------------------------------------------------------

def depmap_run_is_complete(run_dir):
    """True if selectivity.csv and lineage_enrichment.csv are present with >= 1 data row."""
    run_dir = Path(run_dir)
    for fname in DEPMAP_REQUIRED:
        p = run_dir / fname
        if not p.exists():
            return False
        # header + at least one data row
        lines = [l for l in p.read_text().splitlines() if l.strip()]
        if len(lines) < 2:
            return False
    return True


def open_targets_run_is_complete(run_dir, gene):
    """True if all six OT CSVs are present and identity.csv has exactly one row for gene."""
    import pandas as pd
    run_dir = Path(run_dir)
    for fname in OT_REQUIRED:
        if not (run_dir / fname).exists():
            return False
    identity = pd.read_csv(run_dir / "identity.csv")
    if len(identity) != 1:
        return False
    symbol = str(identity.iloc[0].get("approved_symbol", "")).upper()
    return symbol == gene.upper()


def parse_run_timestamp(run_dir):
    """Extract datetime from a run directory named YYYY-MM-DD_HHMMSS, or None."""
    m = re.match(r"(\d{4}-\d{2}-\d{2}_\d{6})$", Path(run_dir).name)
    if not m:
        return None
    try:
        return datetime.strptime(m.group(1), "%Y-%m-%d_%H%M%S")
    except ValueError:
        return None


# ---------------------------------------------------------------------------
# Run selection
# ---------------------------------------------------------------------------

def find_depmap_run(gene, max_age_days=None, results_base=None):
    """Newest complete DepMap run for gene passing staleness check, or None."""
    base = Path(results_base) if results_base else project_root() / "results"
    gene_dir = base / "depmap_local" / gene.lower()
    if not gene_dir.is_dir():
        return None
    candidates = sorted(
        (d for d in gene_dir.iterdir() if d.is_dir()),
        key=lambda d: d.name,
        reverse=True,
    )
    for d in candidates:
        if not depmap_run_is_complete(d):
            continue
        if max_age_days is not None:
            ts = parse_run_timestamp(d)
            if ts is None or (datetime.now() - ts) > timedelta(days=max_age_days):
                continue
        return d
    return None


def find_open_targets_run(gene, max_age_days=None, results_base=None):
    """Newest complete non-stale Open Targets run for gene, or None.

    `max_age_days` defaults to OT_DEFAULT_MAX_AGE_DAYS through an explicit
    None check rather than in the signature: the kernel.py sidecar loader
    rejects a non-literal default, and a rejected file defines none of this
    module's helpers.
    """
    max_age_days = OT_DEFAULT_MAX_AGE_DAYS if max_age_days is None else max_age_days
    base = Path(results_base) if results_base else project_root() / "results"
    gene_dir = base / "open_targets" / gene.lower()
    if not gene_dir.is_dir():
        return None
    candidates = sorted(
        (d for d in gene_dir.iterdir() if d.is_dir()),
        key=lambda d: d.name,
        reverse=True,
    )
    for d in candidates:
        if not open_targets_run_is_complete(d, gene):
            continue
        if max_age_days is not None:
            ts = parse_run_timestamp(d)
            if ts is None or (datetime.now() - ts) > timedelta(days=max_age_days):
                continue
        return d
    return None


# ---------------------------------------------------------------------------
# Fusion helpers
# ---------------------------------------------------------------------------

def lineage_for_disease(disease_name):
    """Best-matching DepMap lineage for an OT disease name (longest keyword wins), or None."""
    if not disease_name:
        return None
    disease_lower = str(disease_name).lower()
    best_lineage = None
    best_kw_len = 0
    for lineage, keywords in LINEAGE_KEYWORDS.items():
        for kw in keywords:
            if kw in disease_lower and len(kw) > best_kw_len:
                best_lineage = lineage
                best_kw_len = len(kw)
    return best_lineage


def stage_to_num(stage_str):
    """Convert 'PHASE_1', 'APPROVAL' etc. to an int (APPROVAL=4), or None.

    Open Targets' maxClinicalStage field uses "APPROVAL", not "APPROVED" --
    confirmed against live known_drugs.csv output for multiple genes. The
    original "APPROVED" in s check never matched, so every approved drug fell
    through to the digit regex, found no digit, and returned None: max_stage
    in verdict.csv silently dropped every approved drug and reported the
    highest *unapproved* phase instead, understating drug-development status
    for any target with an approved therapy.
    """
    if stage_str is None:
        return None
    s = str(stage_str).upper()
    if "APPROV" in s:  # matches APPROVAL and APPROVED, whichever the API uses
        return 4
    m = re.search(r"(\d+)", s)
    return int(m.group(1)) if m else None


def pick_genotype_marker(frames):
    """Return (marker_name, detail_str) for the most informative contrast, or (None, None)."""
    for key, pval_col in [
        ("msi_contrast", "p_msi_more_dependent"),
        ("mutation_contrast_overall", "p_mutant_more_dependent"),
    ]:
        df = frames.get(key)
        if df is None or len(df) == 0:
            continue
        row = df.iloc[0]
        pval = row.get(pval_col)
        if pval is not None and not (isinstance(pval, float) and math.isnan(pval)):
            return (key.replace("_contrast", "").replace("_overall", ""),
                    f"p={float(pval):.2e}")
    df = frames.get("fusion_contrast")
    if df is not None and len(df) > 0:
        status = str(df.iloc[0].get("status", ""))
        if status and status != "untestable":
            return ("fusion", status)
    return (None, None)


def build_verdict(depmap_frames, ot_frames):
    """One-row DataFrame summarising the dependency and clinical picture."""
    import pandas as pd
    sel = depmap_frames["selectivity"]
    lin = depmap_frames["lineage_enrichment"]

    gene = str(sel.iloc[0]["gene"])
    classification = str(sel.iloc[0]["classification"])
    frac_dependent = float(sel.iloc[0]["frac_dependent"])

    sig_lin = lin[lin["p_bh"] < 0.05] if len(lin) > 0 else pd.DataFrame()
    top_lineage_significant = len(sig_lin) > 0
    if len(sig_lin) > 0:
        top = sig_lin.sort_values("p_bh").iloc[0]
    elif len(lin) > 0:
        top = lin.sort_values("p_bh").iloc[0]
    else:
        top = None
    top_lineage = str(top["lineage"]) if top is not None else None
    top_lineage_p_bh = float(top["p_bh"]) if top is not None else None

    marker_name, marker_detail = pick_genotype_marker(depmap_frames)

    tract = ot_frames.get("tractability", pd.DataFrame())
    if len(tract) > 0:
        true_labels = tract[tract["value"] == True]["label"].tolist()
        tractability_modalities = "; ".join(true_labels) if true_labels else None
    else:
        tractability_modalities = None

    drugs = ot_frames.get("known_drugs", pd.DataFrame())
    max_stage = None
    if len(drugs) > 0 and "max_clinical_stage" in drugs.columns:
        nums = drugs["max_clinical_stage"].dropna().apply(stage_to_num)
        nums = nums.dropna()
        if len(nums) > 0:
            max_stage = int(nums.max())

    safety = ot_frames.get("safety", pd.DataFrame())

    return pd.DataFrame([{
        "gene": gene,
        "classification": classification,
        "frac_dependent": frac_dependent,
        "top_lineage": top_lineage,
        "top_lineage_p_bh": top_lineage_p_bh,
        "top_lineage_significant": top_lineage_significant,
        "genotype_marker": marker_name,
        "genotype_marker_detail": marker_detail,
        "tractability_modalities": tractability_modalities,
        "max_clinical_stage": max_stage,
        "has_safety_liabilities": bool(len(safety) > 0),
    }])


def find_conflicts(depmap_frames, ot_frames):
    """Rows describing disagreements between DepMap and Open Targets."""
    import pandas as pd
    sel = depmap_frames["selectivity"]
    lin = depmap_frames["lineage_enrichment"]
    prism = depmap_frames.get("prism_compounds")
    assoc = ot_frames.get("associations", pd.DataFrame())
    tract = ot_frames.get("tractability", pd.DataFrame())
    drugs = ot_frames.get("known_drugs", pd.DataFrame())

    gene = str(sel.iloc[0]["gene"])
    classification = str(sel.iloc[0]["classification"])
    sig_lin = lin[lin["p_bh"] < 0.05] if len(lin) > 0 else pd.DataFrame()
    assoc_scored = assoc.dropna(subset=["overall_score"]) if len(assoc) > 0 else assoc
    rows = []

    # 1. High-scoring OT disease with no significant DepMap lineage enrichment.
    if len(assoc_scored) > 0 and len(sig_lin) == 0:
        top = assoc_scored.nlargest(1, "overall_score").iloc[0]
        if float(top["overall_score"]) > 0.5:
            rows.append({
                "conflict_type": "disease_without_lineage_enrichment",
                "observation_depmap": "no lineage passes BH-adjusted p < 0.05",
                "observation_ot": (
                    f"{top['disease_name']} (score {float(top['overall_score']):.2f})"
                ),
                "description": (
                    f"Open Targets ranks {top['disease_name']} as the top disease "
                    f"association for {gene} (score {float(top['overall_score']):.2f}), "
                    "but no cancer lineage shows significant selective dependency in DepMap "
                    "after multiple-testing correction."
                ),
            })

    # 2. Common-essential with small-molecule tractability.
    if classification == "common essential" and len(tract) > 0:
        sm_true = tract[(tract["modality"] == "SM") & (tract["value"] == True)]
        if len(sm_true) > 0:
            labels = "; ".join(sm_true["label"].tolist())
            rows.append({
                "conflict_type": "common_essential_with_sm_tractability",
                "observation_depmap": "common essential (dependent in >=90% of lines)",
                "observation_ot": f"small-molecule tractability: {labels}",
                "description": (
                    f"{gene} is common-essential in DepMap: loss of function reduces "
                    "viability across nearly all cancer cell lines. Small-molecule "
                    "tractability indicates a druggable pocket exists, but a narrow "
                    "therapeutic window is expected."
                ),
            })

    # 3. PRISM compounds with activity but no drug at any clinical phase in OT.
    if prism is not None and len(prism) > 0 and "median_lfc" in prism.columns:
        active = prism[prism["median_lfc"].notna() & (prism["median_lfc"] < -0.5)]
        if len(active) > 0 and len(drugs) == 0:
            drug_names = "; ".join(active["drug"].head(3).tolist())
            rows.append({
                "conflict_type": "prism_activity_without_clinical_drug",
                "observation_depmap": f"PRISM compounds with activity: {drug_names}",
                "observation_ot": "no drug at any clinical phase",
                "description": (
                    f"PRISM compounds annotated to {gene} reduce viability in cancer "
                    "cell lines (median LFC < -0.5), but Open Targets records no drug "
                    "against this target at any clinical phase."
                ),
            })

    # 4. Selective dependency in one lineage, but top OT disease is not a cancer.
    if classification == "selective" and len(sig_lin) > 0 and len(assoc_scored) > 0:
        top_disease = assoc_scored.nlargest(1, "overall_score").iloc[0]
        disease_lower = str(top_disease["disease_name"]).lower()
        if not any(t in disease_lower for t in CANCER_TERMS):
            top_dep = sig_lin.sort_values("p_bh")["lineage"].head(2).tolist()
            rows.append({
                "conflict_type": "non_cancer_top_disease_with_selective_dependency",
                "observation_depmap": "selective, enriched in " + ", ".join(top_dep),
                "observation_ot": (
                    f"top disease: {top_disease['disease_name']} "
                    f"(score {float(top_disease['overall_score']):.2f})"
                ),
                "description": (
                    f"{gene} shows selective dependency in DepMap, enriched in "
                    f"{', '.join(top_dep)}, but the highest-scoring disease in Open "
                    f"Targets is {top_disease['disease_name']}, which is not a cancer. "
                    "The DepMap signal may reflect a synthetic lethal relationship "
                    "distinct from the gene's primary disease context."
                ),
            })

    cols = ["conflict_type", "observation_depmap", "observation_ot", "description"]
    return pd.DataFrame(rows, columns=cols) if rows else pd.DataFrame(columns=cols)


def build_disease_lineage(assoc, lin):
    """Align OT diseases to DepMap lineages using keyword matching on disease names."""
    import pandas as pd
    lin_index = (
        lin.set_index("lineage")
        if len(lin) > 0 and "lineage" in lin.columns
        else pd.DataFrame()
    )
    rows = []
    for _, row in assoc.iterrows():
        matched = lineage_for_disease(row.get("disease_name"))
        lin_effect = None
        lin_p_bh = None
        # Three outcomes share empty effect columns and are distinguished here:
        # the disease name matched no lineage keyword; it matched a lineage that
        # DepMap did not screen for this gene; or it matched and has data.
        if matched is None:
            status = "no_lineage_keyword"
        elif len(lin_index) > 0 and matched in lin_index.index:
            lr = lin_index.loc[matched]
            lin_effect = float(lr["mean_effect"]) if "mean_effect" in lr else None
            lin_p_bh = float(lr["p_bh"]) if "p_bh" in lr else None
            status = "matched"
        else:
            status = "lineage_not_in_depmap_run"
        rows.append({
            "disease_id": row.get("disease_id"),
            "disease_name": row.get("disease_name"),
            "ot_overall_score": row.get("overall_score"),
            "depmap_lineage": matched,
            "lineage_status": status,
            "lineage_mean_effect": lin_effect,
            "lineage_p_bh": lin_p_bh,
        })
    cols = ["disease_id", "disease_name", "ot_overall_score",
            "depmap_lineage", "lineage_status", "lineage_mean_effect", "lineage_p_bh"]
    return pd.DataFrame(rows, columns=cols) if rows else pd.DataFrame(columns=cols)


def build_drugs_joined(ot_drugs, prism):
    """Join OT known drugs to PRISM compounds on drug name (case-insensitive)."""
    import pandas as pd
    out_cols = ["drug_name", "drug_id", "max_clinical_stage", "action_type",
                "compound_id", "median_lfc", "min_lfc"]
    no_ot = len(ot_drugs) == 0
    no_prism = prism is None or len(prism) == 0

    if no_ot and no_prism:
        return pd.DataFrame(columns=out_cols)

    if not no_ot:
        ot_copy = ot_drugs.copy()
        ot_copy["_stage_num"] = ot_copy["max_clinical_stage"].apply(stage_to_num)
        ot_dedup = (
            ot_copy
            .sort_values("_stage_num", ascending=False, na_position="last")
            .drop_duplicates(subset=["drug_id"])
            [["drug_id", "drug_name", "max_clinical_stage", "action_type"]]
            .copy()
        )
        ot_dedup["_key"] = ot_dedup["drug_name"].str.lower().str.strip()
    else:
        ot_dedup = pd.DataFrame(columns=["drug_id", "drug_name", "max_clinical_stage",
                                          "action_type", "_key"])

    if not no_prism:
        prism_copy = prism[["drug", "median_lfc", "min_lfc"]].copy()
        if "compound_id" in prism.columns:
            prism_copy["compound_id"] = prism["compound_id"]
        else:
            prism_copy["compound_id"] = None
        prism_copy["_key"] = prism_copy["drug"].str.lower().str.strip()
    else:
        prism_copy = pd.DataFrame(columns=["drug", "compound_id", "median_lfc",
                                            "min_lfc", "_key"])

    merged = ot_dedup.merge(prism_copy, on="_key", how="outer")
    merged["drug_name"] = merged["drug_name"].fillna(merged.get("drug", pd.Series(dtype=str)))
    for col in out_cols:
        if col not in merged.columns:
            merged[col] = None
    return merged[out_cols].reset_index(drop=True)


# ---------------------------------------------------------------------------
# Fusion entry point
# ---------------------------------------------------------------------------

def fuse(depmap_dir, ot_dir):
    """Read both source run directories and produce the four fused tables.

    Returns dict with keys: verdict, conflicts, disease_lineage, drugs_joined.
    Each value is a DataFrame.
    """
    import pandas as pd
    depmap_data = read_depmap_run(depmap_dir)
    ot_data = read_open_targets_run(ot_dir)
    df = depmap_data["frames"]
    of = ot_data["frames"]
    return {
        "verdict": build_verdict(df, of),
        "conflicts": find_conflicts(df, of),
        "disease_lineage": build_disease_lineage(
            of.get("associations", pd.DataFrame()),
            df.get("lineage_enrichment", pd.DataFrame()),
        ),
        "drugs_joined": build_drugs_joined(
            of.get("known_drugs", pd.DataFrame()),
            df.get("prism_compounds"),
        ),
        "lineage_enrichment": df.get("lineage_enrichment", pd.DataFrame()),
    }


# ---------------------------------------------------------------------------
# Write
# ---------------------------------------------------------------------------

FUSION_CSV_TABLES = {"verdict", "conflicts", "disease_lineage", "drugs_joined"}


def write_fusion_run(tables, sources, out_dir):
    """Write verdict, conflicts, disease_lineage, drugs_joined as CSVs; write sources.json."""
    out_dir = Path(out_dir)
    written = {}
    for name, df in tables.items():
        if name not in FUSION_CSV_TABLES:
            continue
        path = out_dir / f"{name}.csv"
        df.to_csv(path, index=False)
        written[name] = path
    sources_path = out_dir / "sources.json"
    sources_path.write_text(json.dumps(sources, indent=2, default=str))
    written["sources"] = sources_path
    return written


def write_fusion_readme(out_dir, gene, tables, sources):
    """Write README.md with a plain-language intro paragraph and per-section findings."""
    import pandas as pd
    from lib.common.reporting import markdown_table

    verdict = tables["verdict"]
    conflicts = tables["conflicts"]
    disease_lineage = tables["disease_lineage"]
    drugs_joined = tables["drugs_joined"]

    v = verdict.iloc[0]
    classification = str(v["classification"])
    top_lineage = v["top_lineage"]
    top_lineage_significant = bool(v["top_lineage_significant"])
    max_stage = v["max_clinical_stage"]
    has_safety = bool(v["has_safety_liabilities"])
    tractability = v.get("tractability_modalities")

    # Intro paragraph
    if top_lineage and top_lineage_significant:
        lineage_clause = f", with the strongest selective dependency in the {top_lineage} lineage"
    elif top_lineage:
        lineage_clause = (
            f", with the {top_lineage} lineage showing the highest-ranked but "
            "not statistically significant enrichment"
        )
    else:
        lineage_clause = ""
    intro = f"{gene} is {classification} in DepMap" + lineage_clause + ". "
    if max_stage is not None and not (isinstance(max_stage, float) and pd.isna(max_stage)):
        s = int(max_stage)
        intro += (
            f"An approved drug against {gene} exists in Open Targets."
            if s >= 4
            else f"At least one drug against {gene} has reached Phase {s} clinical trials in Open Targets."
        )
    else:
        intro += f"No drug against {gene} appears in clinical trials in Open Targets."
    if len(conflicts) > 0:
        n_c = len(conflicts)
        word = "conflict" if n_c == 1 else "conflicts"
        intro += f" {n_c} {word} between the two sources are noted under Conflicts below."

    lin_enrich = tables.get("lineage_enrichment", pd.DataFrame())

    # SM tractability headline: highest tier that is True
    SM_TIERS = [
        ("Approved Drug", "an approved drug"),
        ("Advanced Clinical", "a drug at advanced clinical stage"),
        ("Phase 1 Clinical", "a drug at phase 1"),
        ("Structure with Ligand", "a structural binding site with a ligand"),
        ("High-Quality Ligand", "a high-quality ligand"),
        ("Medium-Quality Ligand", "a medium-quality ligand"),
        ("Small Molecule Binder", "a small molecule binder"),
    ]
    tract_str = str(tractability) if tractability else ""
    sm_headline = next((desc for label, desc in SM_TIERS if label in tract_str), None)

    steps = []

    # Dependency
    frac = float(v["frac_dependent"])
    dep_finding = f"{gene} is {classification} (dependent in {frac:.1%} of screened lines)."
    if top_lineage:
        p_bh = v.get("top_lineage_p_bh")
        p_bh_str = (
            f" (BH p = {float(p_bh):.2e})"
            if p_bh is not None and not (isinstance(p_bh, float) and pd.isna(p_bh))
            else ""
        )
        if top_lineage_significant:
            dep_finding += f" The {top_lineage} lineage shows the strongest enrichment{p_bh_str}."
        else:
            dep_finding += (
                f" The {top_lineage} lineage ranks highest{p_bh_str}, "
                "but does not reach BH-adjusted significance."
            )
    marker = v.get("genotype_marker")
    if marker:
        detail = v.get("genotype_marker_detail") or ""
        dep_finding += f" A {marker} genotype contrast was run ({detail})."

    lin_step = {"name": "Dependency", "finding": dep_finding}
    if len(lin_enrich) > 0:
        sig = lin_enrich[lin_enrich["p_bh"] < 0.05].sort_values("p_bh")
        display = sig if len(sig) > 0 else lin_enrich.sort_values("p_bh").head(5)
        lin_rows = []
        for _, r in display.iterrows():
            lin_rows.append({
                "lineage": r["lineage"],
                "n": int(r["n"]),
                "mean_effect": f"{float(r['mean_effect']):.3f}",
                "frac_dependent": f"{float(r['frac_dependent']):.2f}",
                "p_bh": f"{float(r['p_bh']):.2e}",
            })
        if lin_rows:
            lin_step["table"] = lin_rows
    steps.append(lin_step)

    # Disease associations
    if len(disease_lineage) > 0:
        top5 = disease_lineage.nlargest(5, "ot_overall_score")
        lead = str(top5.iloc[0]["disease_name"])
        assoc_finding = f"The highest-scoring association in Open Targets is {lead}."
        top_rows = []
        for _, r in top5.iterrows():
            score = r["ot_overall_score"]
            effect = r["lineage_mean_effect"]
            top_rows.append({
                "disease_name": r["disease_name"],
                "ot_score": ("" if pd.isna(score) else f"{float(score):.2f}"),
                "depmap_lineage": r["depmap_lineage"] if r["depmap_lineage"] else "",
                "lineage_status": r["lineage_status"],
                "lineage_mean_effect": ("" if pd.isna(effect) else f"{float(effect):.3f}"),
            })
        steps.append({
            "name": "Disease associations",
            "finding": assoc_finding,
            "table": top_rows,
        })
    else:
        steps.append({
            "name": "Disease associations",
            "skipped": "No disease associations returned by Open Targets.",
        })

    # Tractability and drugs
    if sm_headline:
        tract_finding = f"Open Targets records {sm_headline} for {gene}."
    else:
        tract_finding = f"No small-molecule tractability is recorded in Open Targets for {gene}."
    ot_drugs = drugs_joined[drugs_joined["drug_id"].notna()] if len(drugs_joined) > 0 else pd.DataFrame()
    if len(ot_drugs) > 0:
        drug_rows = []
        for _, r in drugs_joined.head(10).iterrows():
            stage = r.get("max_clinical_stage")
            lfc = r.get("median_lfc")
            drug_rows.append({
                "drug_name": r.get("drug_name", ""),
                "max_clinical_stage": "" if pd.isna(stage) else stage,
                "compound_id": r.get("compound_id", "") or "",
                "median_lfc": "" if pd.isna(lfc) else f"{float(lfc):.3f}",
            })
        steps.append({
            "name": "Tractability and drugs",
            "finding": tract_finding,
            "table": drug_rows,
        })
    else:
        steps.append({"name": "Tractability and drugs", "finding": tract_finding})

    # Safety
    steps.append({
        "name": "Safety",
        "finding": (
            "Safety liabilities are recorded in Open Targets. See safety.csv for details."
            if has_safety
            else "No safety liabilities are recorded in Open Targets."
        ),
    })

    # Conflicts
    if len(conflicts) > 0:
        steps.append({
            "name": "Conflicts",
            "finding": "The DepMap and Open Targets data disagree as follows.",
            "table": conflicts[["conflict_type", "description"]].to_dict("records"),
        })
    else:
        steps.append({
            "name": "Conflicts",
            "finding": "No conflicts detected between the two sources.",
        })

    # Provenance
    prov = []
    for src_name, src in sources.items():
        action = "reused" if src.get("reused") else "generated"
        ts = src.get("timestamp", "unknown")
        prov.append(f"{src_name} run {ts} ({action})")

    write_readme(
        path=Path(out_dir) / "README.md",
        subject=gene,
        summary=intro,
        steps=steps,
        source="DepMap + Open Targets",
        limits=[
            "Disease-to-lineage mapping uses keyword matching on disease names; "
            "diseases without a recognised keyword are unassigned.",
            "Drug join uses case-insensitive name matching. "
            "compound_id in drugs_joined.csv identifies the specific PRISM compound.",
        ] + (["Provenance: " + "; ".join(prov)] if prov else []),
    )


# ---------------------------------------------------------------------------
# Ensure upstream runs exist
# ---------------------------------------------------------------------------

def resolve_ensembl_id(gene):
    """Resolve a gene symbol to an Ensembl ID via the Open Targets search API."""
    import requests
    r = requests.post(
        OT_GRAPHQL_URL,
        json={"query": OT_SEARCH_QUERY, "variables": {"q": gene}},
        timeout=30,
    )
    r.raise_for_status()
    hits = r.json().get("data", {}).get("search", {}).get("hits", [])
    for h in hits:
        if str(h.get("name", "")).upper() == gene.upper():
            return h["id"]
    if hits:
        return hits[0]["id"]
    raise ValueError(f"no Open Targets target found for gene symbol {gene!r}")


def ensure_depmap_run(gene, reuse=True, max_age_days=None, root=None,
                      dry_run=False, results_base=None):
    """Return a valid DepMap run directory for gene, running the analysis if needed.

    reuse=True checks for an existing complete run first. dry_run=True reports
    what would be run without executing.
    """
    if reuse:
        existing = find_depmap_run(gene, max_age_days=max_age_days,
                                   results_base=results_base)
        if existing is not None:
            return existing

    if dry_run:
        print(f"[dry_run] Would run depmap-local for {gene}")
        return None

    import importlib.util, sys as _sys
    kernel_path = Path(__file__).resolve().parent.parent / "depmap-local" / "kernel.py"
    spec = importlib.util.spec_from_file_location("depmap_local_kernel", kernel_path)
    dk = importlib.util.module_from_spec(spec)
    _sys.modules[spec.name] = dk
    spec.loader.exec_module(dk)

    out = results_dir("depmap-local", gene.lower())
    tables = {}
    tables["selectivity"] = dk.depmap_selectivity(gene, root=root)
    tables["lineage_enrichment"] = dk.depmap_lineage_enrichment(gene, root=root)
    for fn, call in [
        ("fusion_contrast", lambda: dk.depmap_fusion_contrast(gene, root=root)),
        ("msi_contrast", lambda: dk.depmap_msi_contrast(gene, root=root)),
    ]:
        try:
            tables[fn] = call()
        except (FileNotFoundError, KeyError):
            pass
    prism = dk.depmap_prism_compounds(gene, root=root)
    if len(prism) > 0:
        tables["prism_compounds"] = prism
    dk.write_run(tables, out)
    return out


def ensure_open_targets_run(gene, reuse=True, max_age_days=None,
                             dry_run=False, results_base=None):
    """Return a valid Open Targets run directory for gene, fetching data if needed.

    `max_age_days` defaults to OT_DEFAULT_MAX_AGE_DAYS through an explicit
    None check rather than in the signature; see find_open_targets_run.
    """
    max_age_days = OT_DEFAULT_MAX_AGE_DAYS if max_age_days is None else max_age_days
    if reuse:
        existing = find_open_targets_run(gene, max_age_days=max_age_days,
                                         results_base=results_base)
        if existing is not None:
            return existing

    if dry_run:
        print(f"[dry_run] Would run open-targets for {gene}")
        return None

    import importlib.util, sys as _sys, requests
    kernel_path = Path(__file__).resolve().parent.parent / "open-targets" / "kernel.py"
    spec = importlib.util.spec_from_file_location("open_targets_kernel", kernel_path)
    ok = importlib.util.module_from_spec(spec)
    _sys.modules[spec.name] = ok
    spec.loader.exec_module(ok)

    ensembl_id = resolve_ensembl_id(gene)

    query_path = Path(__file__).resolve().parent.parent / "open-targets" / "query.graphql"
    query = query_path.read_text()

    r = requests.post(
        OT_GRAPHQL_URL,
        json={"query": query, "variables": {"ensemblId": ensembl_id}},
        timeout=60,
    )
    r.raise_for_status()
    data = r.json().get("data", {})

    out = results_dir("open-targets", gene.lower())
    ok.write_run(data, out)
    return out


# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------

def run_depmap_ot_fusion(gene, reuse=True, max_age_days=None,
                          root=None, dry_run=False):
    """Run the full DepMap + Open Targets fusion for gene. Returns the output directory.

    With reuse=True, existing upstream runs are reused when complete and fresh.
    With dry_run=True, prints planned actions and returns None.

    `max_age_days` defaults to OT_DEFAULT_MAX_AGE_DAYS through an explicit
    None check rather than in the signature; see find_open_targets_run.
    """
    max_age_days = OT_DEFAULT_MAX_AGE_DAYS if max_age_days is None else max_age_days
    gene_upper = gene.upper()

    # Whether each source will be reused is resolved before the ensure_* calls
    # so sources.json can record what happened rather than what was requested:
    # reuse=True with nothing cached generates a fresh run.
    depmap_reused = reuse and find_depmap_run(gene_upper) is not None
    ot_reused = reuse and find_open_targets_run(
        gene.lower(), max_age_days=max_age_days) is not None

    depmap_dir = ensure_depmap_run(
        gene_upper, reuse=reuse, max_age_days=None, root=root, dry_run=dry_run,
    )
    ot_dir = ensure_open_targets_run(
        gene.lower(), reuse=reuse, max_age_days=max_age_days, dry_run=dry_run,
    )

    if dry_run:
        return None

    tables = fuse(depmap_dir, ot_dir)
    out = results_dir("depmap-ot-fusion", gene.lower())

    depmap_data = read_depmap_run(depmap_dir)
    sources = {
        "depmap": {
            "run_dir": str(depmap_dir),
            "timestamp": depmap_dir.name,
            "reused": depmap_reused,
            "optional_files": depmap_data["optional_files"],
        },
        "open_targets": {
            "run_dir": str(ot_dir),
            "timestamp": ot_dir.name,
            "reused": ot_reused,
        },
    }

    write_fusion_run(tables, sources, out)
    write_fusion_readme(out, gene_upper, tables, sources)
    return out
