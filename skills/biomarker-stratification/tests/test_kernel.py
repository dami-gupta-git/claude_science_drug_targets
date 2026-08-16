import pandas as pd
import pytest

import kernel

STUDY = {
    "protocolSection": {
        "identificationModule": {"nctId": "NCT00000001", "briefTitle": "A trial"},
        "designModule": {"phases": ["PHASE2"]},
        "conditionsModule": {"conditions": ["Non-small cell lung cancer"]},
        "eligibilityModule": {
            "eligibilityCriteria": "Inclusion: EGFR exon 19 deletion or L858R "
                                   "mutation confirmed by tissue biopsy.",
        },
    },
}

STUDIES = {"NCT00000001": STUDY, "_errors": ["NCT99999999"]}


def test_load_sibling_kernel_imports_literature_evidence():
    lit = kernel.load_sibling_kernel("literature-evidence")
    assert hasattr(lit, "build_query")
    assert hasattr(lit, "search_europepmc")


def test_load_sibling_kernel_imports_ai_origination_audit():
    audit = kernel.load_sibling_kernel("ai-origination-audit")
    assert hasattr(audit, "fetch_ctgov_studies")


def test_extract_eligibility_row():
    row = kernel.extract_eligibility_row(STUDY)
    assert row["nct_id"] == "NCT00000001"
    assert row["phase"] == "PHASE2"
    assert row["conditions"] == "Non-small cell lung cancer"
    assert "EGFR exon 19 deletion" in row["eligibility_criteria"]


def test_trial_eligibility_table_skips_errors_key():
    df = kernel.trial_eligibility_table(STUDIES)
    assert len(df) == 1
    assert list(df.columns) == kernel.TRIAL_ELIGIBILITY_COLUMNS
    assert df.iloc[0]["nct_id"] == "NCT00000001"


def test_trial_eligibility_table_empty():
    df = kernel.trial_eligibility_table({})
    assert len(df) == 0
    assert list(df.columns) == kernel.TRIAL_ELIGIBILITY_COLUMNS


VALID_BIOMARKER = {
    "biomarker_id": "B1", "disease": "NSCLC", "biomarker": "EGFR exon 19 deletion",
    "evidence_type": "registry_criterion", "source_id": "NCT00000001",
    "evidence_text": "Inclusion: EGFR exon 19 deletion or L858R mutation "
                     "confirmed by tissue biopsy.",
    "url": "https://clinicaltrials.gov/study/NCT00000001",
}


def test_biomarkers_table_valid_row():
    df = kernel.biomarkers_table([VALID_BIOMARKER])
    assert len(df) == 1
    assert list(df.columns) == kernel.BIOMARKER_COLUMNS


def test_biomarkers_table_empty():
    df = kernel.biomarkers_table([])
    assert len(df) == 0
    assert list(df.columns) == kernel.BIOMARKER_COLUMNS


def test_biomarkers_table_missing_field_raises():
    bad = {**VALID_BIOMARKER, "evidence_text": ""}
    with pytest.raises(ValueError, match="missing required field"):
        kernel.biomarkers_table([bad])


def test_biomarkers_table_bad_evidence_type_raises():
    bad = {**VALID_BIOMARKER, "evidence_type": "vibes"}
    with pytest.raises(ValueError, match="not in"):
        kernel.biomarkers_table([bad])


SAMPLE_STEPS = [
    {"name": "Trial eligibility", "finding": "One trial reviewed."},
    {"name": "Biomarkers", "finding": "One registry-criterion biomarker recorded."},
]


def test_write_run(tmp_path):
    publications_df = pd.DataFrame(
        [{"pmid": "1", "title": "x"}], columns=["pmid", "title"])
    written = kernel.write_run(publications_df, STUDIES, [VALID_BIOMARKER],
                               tmp_path, "NSCLC", "Summary.", SAMPLE_STEPS)
    assert set(written) == {
        "literature_publications", "trial_eligibility", "biomarkers", "readme",
    }
    for path in written.values():
        assert path.is_file()
    readme = written["readme"].read_text()
    assert "# NSCLC — Europe PMC / ClinicalTrials.gov" in readme
