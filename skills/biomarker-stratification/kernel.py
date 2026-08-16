"""Find biomarkers used to stratify patients for a disease's clinical trials.

Reuses literature-evidence's Europe PMC access and ai-origination-audit's
ClinicalTrials.gov bulk fetch rather than duplicating them - see
load_sibling_kernel. kernel.py loads with this skill: load_sibling_kernel,
search_ctgov_by_condition, extract_eligibility_row, trial_eligibility_table,
biomarkers_table, write_run, and the constants CTGOV_STUDIES_API,
TRIAL_ELIGIBILITY_COLUMNS, BIOMARKER_EVIDENCE_TYPES, BIOMARKER_COLUMNS.
"""
import json
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

from lib.common.reporting import write_readme

# Same ClinicalTrials.gov API v2 endpoint ai-origination-audit's CTGOV_API
# points at. Kept as its own constant here (not imported) because only the
# URL is shared - the query shape (condition, not intervention/sponsor) is
# specific to this skill.
CTGOV_STUDIES_API = "https://clinicaltrials.gov/api/v2/studies"

TRIAL_ELIGIBILITY_COLUMNS = [
    "nct_id", "title", "phase", "conditions", "eligibility_criteria",
]

BIOMARKER_EVIDENCE_TYPES = ("registry_criterion", "literature_claim")

BIOMARKER_COLUMNS = [
    "biomarker_id", "disease", "biomarker", "evidence_type", "source_id",
    "evidence_text", "url",
]


def load_sibling_kernel(skill_dir: str):
    """Dynamically import another skill's kernel.py as <skill_dir>_kernel.

    Skill directories are hyphenated, so they can't be imported as normal
    Python packages; this loads a sibling kernel.py by file path, the same
    pattern depmap-ot-fusion uses to reuse depmap-local and open-targets.
    """
    import importlib.util
    import sys
    here = Path(__file__).resolve().parent
    kernel_path = here.parent / skill_dir / "kernel.py"
    module_name = skill_dir.replace("-", "_") + "_kernel"
    spec = importlib.util.spec_from_file_location(module_name, kernel_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def search_ctgov_by_condition(condition: str, page_size: int = 100,
                              max_pages: int = None, timeout: int = 60) -> list:
    """NCT ids of studies matching condition, paginated via pageToken.

    ai-origination-audit's search_ctgov searches by intervention/sponsor
    only; this fills the condition-search gap it doesn't cover.
    """
    nct_ids = []
    page_token = None
    pages = 0
    while True:
        params = {"query.cond": condition, "pageSize": page_size,
                 "fields": "NCTId"}
        if page_token:
            params["pageToken"] = page_token
        url = "%s?%s" % (CTGOV_STUDIES_API, urllib.parse.urlencode(params))
        with urllib.request.urlopen(url, timeout=timeout) as fh:
            data = json.load(fh)
        for study in data.get("studies") or []:
            nct_id = (study.get("protocolSection", {})
                          .get("identificationModule", {})
                          .get("nctId"))
            if nct_id:
                nct_ids.append(nct_id)
        pages += 1
        page_token = data.get("nextPageToken")
        if not page_token or (max_pages and pages >= max_pages):
            break
    return nct_ids


def extract_eligibility_row(study: dict) -> dict:
    """One row of eligibility criteria text from a full ClinicalTrials.gov
    API v2 study record (as returned by ai-origination-audit's
    fetch_ctgov_study/fetch_ctgov_studies)."""
    protocol = study.get("protocolSection", {})
    ident = protocol.get("identificationModule", {})
    design = protocol.get("designModule", {})
    conditions = protocol.get("conditionsModule", {}).get("conditions") or []
    eligibility = protocol.get("eligibilityModule", {})
    return {
        "nct_id": ident.get("nctId"),
        "title": ident.get("briefTitle"),
        "phase": ", ".join(design.get("phases") or []),
        "conditions": "; ".join(conditions),
        "eligibility_criteria": eligibility.get("eligibilityCriteria"),
    }


def trial_eligibility_table(studies: dict) -> pd.DataFrame:
    """One row per trial, from {nct_id: study} as returned by
    ai-origination-audit's fetch_ctgov_studies. Entries under the special
    "_errors" key (failed fetches) are skipped, not fabricated as rows."""
    rows = [extract_eligibility_row(study)
            for nct_id, study in studies.items() if nct_id != "_errors"]
    return pd.DataFrame(rows, columns=TRIAL_ELIGIBILITY_COLUMNS)


def biomarkers_table(biomarker_rows: list) -> pd.DataFrame:
    """One row per biomarker, hand-extracted from a trial's eligibility
    criteria text or a literature claim.

    Every row must carry evidence_text (the exact criterion sentence or
    literature quote) and source_id (the NCT id or PMID it came from) - a
    biomarker without a locatable source sentence must not be recorded.
    evidence_type must be one of BIOMARKER_EVIDENCE_TYPES.
    """
    for row in biomarker_rows:
        missing = [c for c in ("disease", "biomarker", "evidence_type",
                               "source_id", "evidence_text")
                  if not row.get(c)]
        if missing:
            raise ValueError("biomarker row %r missing required field(s): %s"
                             % (row, missing))
        if row["evidence_type"] not in BIOMARKER_EVIDENCE_TYPES:
            raise ValueError("evidence_type %r not in %s"
                             % (row["evidence_type"], BIOMARKER_EVIDENCE_TYPES))
    return pd.DataFrame(biomarker_rows, columns=BIOMARKER_COLUMNS)


def write_run(publications_df: pd.DataFrame, studies: dict, biomarker_rows: list,
              out_dir: Path, subject: str, summary: str, steps: list,
              source: str = "Europe PMC / ClinicalTrials.gov",
              limits: list = None) -> dict:
    """Write literature_publications.csv, trial_eligibility.csv,
    biomarkers.csv, and README.md into out_dir."""
    written = {}
    pubs_path = out_dir / "literature_publications.csv"
    publications_df.to_csv(pubs_path, index=False)
    written["literature_publications"] = pubs_path

    trials_path = out_dir / "trial_eligibility.csv"
    trial_eligibility_table(studies).to_csv(trials_path, index=False)
    written["trial_eligibility"] = trials_path

    biomarkers_path = out_dir / "biomarkers.csv"
    biomarkers_table(biomarker_rows).to_csv(biomarkers_path, index=False)
    written["biomarkers"] = biomarkers_path

    written["readme"] = write_readme(out_dir / "README.md", subject, summary,
                                     steps, source=source, limits=limits)
    return written
