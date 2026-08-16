"""Fetch and shape gene expression across normal tissue and tumor types.

kernel.py loads with this skill: resolve_gtex_gene, fetch_gtex_expression,
tissue_expression_table, resolve_cbioportal_gene, fetch_cbioportal_expression,
tumor_expression_table, classify_level, write_run, and the constants
GTEX_GENE_URL, GTEX_EXPRESSION_URL, GTEX_DATASET_ID, CBIOPORTAL_BASE,
CANCER_TYPE_STUDIES, EXPRESSION_FOLD_THRESHOLD, TISSUE_EXPRESSION_COLUMNS,
TUMOR_EXPRESSION_COLUMNS.
"""
import json
import statistics
import urllib.error
import urllib.parse
import urllib.request

import pandas as pd

from lib.common.reporting import write_readme

GTEX_GENE_URL = "https://gtexportal.org/api/v2/reference/gene"
GTEX_EXPRESSION_URL = "https://gtexportal.org/api/v2/expression/medianGeneExpression"
GTEX_DATASET_ID = "gtex_v8"

CBIOPORTAL_BASE = "https://www.cbioportal.org/api"

# Fixed set of major TCGA PanCancer Atlas studies, the same role DepMap's
# lineage column plays in depmap-local: broad tissue-of-origin coverage
# without a live per-run study lookup against cBioPortal's full list.
CANCER_TYPE_STUDIES = {
    "brca_tcga_pan_can_atlas_2018": "Breast invasive carcinoma",
    "luad_tcga_pan_can_atlas_2018": "Lung adenocarcinoma",
    "coadread_tcga_pan_can_atlas_2018": "Colorectal adenocarcinoma",
    "paad_tcga_pan_can_atlas_2018": "Pancreatic adenocarcinoma",
    "skcm_tcga_pan_can_atlas_2018": "Skin cutaneous melanoma",
    "prad_tcga_pan_can_atlas_2018": "Prostate adenocarcinoma",
    "ov_tcga_pan_can_atlas_2018": "Ovarian serous cystadenocarcinoma",
    "gbm_tcga_pan_can_atlas_2018": "Glioblastoma multiforme",
    "kirc_tcga_pan_can_atlas_2018": "Kidney renal clear cell carcinoma",
    "hnsc_tcga_pan_can_atlas_2018": "Head and neck squamous cell carcinoma",
    "lihc_tcga_pan_can_atlas_2018": "Liver hepatocellular carcinoma",
    "blca_tcga_pan_can_atlas_2018": "Bladder urothelial carcinoma",
}

# A tissue/cancer type is "high" at >= 4x the gene's own median across the
# table, "low" at <= 1/4x. Fixed here rather than decided per run, the same
# role depmap-local's DEP_THRESHOLD plays for dependency calls.
EXPRESSION_FOLD_THRESHOLD = 4.0

TISSUE_EXPRESSION_COLUMNS = [
    "tissue", "median_tpm", "fold_vs_gene_median", "level",
]

TUMOR_EXPRESSION_COLUMNS = [
    "cancer_type", "study_id", "median_expression", "n_samples",
    "fold_vs_gene_median", "level",
]


def classify_level(value: float, reference_median: float,
                   threshold: float = EXPRESSION_FOLD_THRESHOLD) -> tuple:
    """Classify value against reference_median as ("level", fold).

    level is "high" at fold >= threshold, "low" at fold <= 1/threshold,
    "typical" otherwise, "unknown" when value or reference_median is
    missing or reference_median is zero (fold is undefined).
    """
    if value is None or reference_median in (None, 0):
        return "unknown", None
    fold = value / reference_median
    if fold >= threshold:
        return "high", fold
    if fold <= 1 / threshold:
        return "low", fold
    return "typical", fold


def resolve_gtex_gene(gene: str, timeout: int = 60) -> str:
    """The GTEx gencodeId for a gene symbol, or raise if not found."""
    params = {"geneId": gene, "format": "json"}
    url = "%s?%s" % (GTEX_GENE_URL, urllib.parse.urlencode(params))
    with urllib.request.urlopen(url, timeout=timeout) as fh:
        data = json.load(fh)
    for record in data.get("data") or []:
        if record.get("geneSymbol", "").upper() == gene.upper():
            return record["gencodeId"]
    raise ValueError("no GTEx gene record found for %r" % gene)


def fetch_gtex_expression(gencode_id: str, dataset_id: str = GTEX_DATASET_ID,
                          timeout: int = 60) -> list:
    """Raw GTEx per-tissue median expression records for gencode_id."""
    params = {"gencodeId": gencode_id, "datasetId": dataset_id, "format": "json"}
    url = "%s?%s" % (GTEX_EXPRESSION_URL, urllib.parse.urlencode(params))
    with urllib.request.urlopen(url, timeout=timeout) as fh:
        data = json.load(fh)
    return data.get("data") or []


def tissue_expression_table(records: list) -> pd.DataFrame:
    """One row per tissue, classified against the gene's own median here."""
    values = [r["median"] for r in records if r.get("median") is not None]
    reference = statistics.median(values) if values else None
    rows = []
    for r in records:
        level, fold = classify_level(r.get("median"), reference)
        rows.append({
            "tissue": r.get("tissueSiteDetailId"),
            "median_tpm": r.get("median"),
            "fold_vs_gene_median": fold,
            "level": level,
        })
    return pd.DataFrame(rows, columns=TISSUE_EXPRESSION_COLUMNS)


def resolve_cbioportal_gene(gene: str, timeout: int = 60) -> int:
    """The cBioPortal Entrez gene id for a gene symbol."""
    url = "%s/genes/%s" % (CBIOPORTAL_BASE, urllib.parse.quote(gene))
    with urllib.request.urlopen(url, timeout=timeout) as fh:
        data = json.load(fh)
    return data["entrezGeneId"]


def fetch_cbioportal_expression(entrez_gene_id: int,
                                studies: dict = CANCER_TYPE_STUDIES,
                                timeout: int = 60) -> dict:
    """Per-sample expression values for entrez_gene_id in each study.

    Returns {study_id: [value, ...]}. A study whose RNA-seq profile or
    sample list isn't found under the expected naming is omitted rather
    than raising, since PanCancer Atlas naming isn't guaranteed forever -
    the caller sees fewer studies in the output, not a crashed run.
    """
    results = {}
    for study_id in studies:
        profile_id = "%s_rna_seq_v2_mrna" % study_id
        sample_list_id = "%s_all" % study_id
        url = "%s/molecular-profiles/%s/molecular-data/fetch?projection=SUMMARY" % (
            CBIOPORTAL_BASE, profile_id)
        body = json.dumps({
            "entrezGeneIds": [entrez_gene_id],
            "sampleListId": sample_list_id,
        }).encode()
        req = urllib.request.Request(
            url, data=body, method="POST",
            headers={"Content-Type": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=timeout) as fh:
                data = json.load(fh)
        except urllib.error.HTTPError:
            continue
        values = [d["value"] for d in data if d.get("value") is not None]
        if values:
            results[study_id] = values
    return results


def tumor_expression_table(study_values: dict,
                           studies: dict = CANCER_TYPE_STUDIES) -> pd.DataFrame:
    """One row per cancer type with data, classified against the gene's
    own median across the studies present in study_values."""
    medians = {sid: statistics.median(vals) for sid, vals in study_values.items()}
    reference = statistics.median(medians.values()) if medians else None
    rows = []
    for study_id, median_value in medians.items():
        level, fold = classify_level(median_value, reference)
        rows.append({
            "cancer_type": studies.get(study_id, study_id),
            "study_id": study_id,
            "median_expression": median_value,
            "n_samples": len(study_values[study_id]),
            "fold_vs_gene_median": fold,
            "level": level,
        })
    return pd.DataFrame(rows, columns=TUMOR_EXPRESSION_COLUMNS)


def write_run(tissue_records: list, study_values: dict, out_dir, subject: str,
              summary: str, steps: list, source: str = "GTEx / cBioPortal",
              limits: list = None) -> dict:
    """Write tissue_expression.csv, tumor_expression.csv and README.md."""
    written = {}
    tissue_path = out_dir / "tissue_expression.csv"
    tissue_expression_table(tissue_records).to_csv(tissue_path, index=False)
    written["tissue_expression"] = tissue_path

    tumor_path = out_dir / "tumor_expression.csv"
    tumor_expression_table(study_values).to_csv(tumor_path, index=False)
    written["tumor_expression"] = tumor_path

    written["readme"] = write_readme(out_dir / "README.md", subject, summary,
                                     steps, source=source, limits=limits)
    return written
