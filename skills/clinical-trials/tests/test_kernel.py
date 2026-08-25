import pytest

import kernel

STUDY_PHASE2 = {
    "hasResults": True,
    "protocolSection": {
        "identificationModule": {
            "nctId": "NCT01000001",
            "briefTitle": "Drug X in KRAS G12C NSCLC",
        },
        "statusModule": {
            "overallStatus": "COMPLETED",
            "startDateStruct": {"date": "2020-03-01", "type": "ACTUAL"},
            "primaryCompletionDateStruct": {"date": "2022-06-01",
                                            "type": "ACTUAL"},
        },
        "sponsorCollaboratorsModule": {"leadSponsor": {"name": "Acme Pharma"}},
        "conditionsModule": {"conditions": ["Non-Small Cell Lung Cancer"]},
        "designModule": {
            "studyType": "INTERVENTIONAL",
            "phases": ["PHASE1", "PHASE2"],
            "enrollmentInfo": {"count": 300, "type": "ACTUAL"},
        },
        "armsInterventionsModule": {
            "armGroups": [
                {"label": "Combo", "interventionNames":
                    ["Drug: Drug X", "Drug: Pembrolizumab"]},
            ],
            "interventions": [
                {"type": "DRUG", "name": "Drug X"},
                {"type": "DRUG", "name": "Pembrolizumab"},
            ],
        },
        "eligibilityModule": {
            "eligibilityCriteria":
                "Inclusion: KRAS G12C mutation. Prior platinum chemotherapy.",
        },
    },
}

STUDY_TERMINATED = {
    "hasResults": False,
    "protocolSection": {
        "identificationModule": {
            "nctId": "NCT01000002",
            "briefTitle": "Drug X monotherapy",
        },
        "statusModule": {
            "overallStatus": "TERMINATED",
            "whyStopped": "Development change",
            "startDateStruct": {"date": "2021-10-01", "type": "ACTUAL"},
        },
        "designModule": {
            "studyType": "INTERVENTIONAL",
            "phases": ["PHASE3"],
            "enrollmentInfo": {"count": 50, "type": "ESTIMATED"},
        },
    },
}

STUDY_UNPHASED = {
    "hasResults": False,
    "protocolSection": {
        "identificationModule": {
            "nctId": "NCT01000003",
            "briefTitle": "Observational registry",
        },
        "statusModule": {"overallStatus": "RECRUITING"},
        "designModule": {"studyType": "OBSERVATIONAL"},
    },
}

STUDY_EMPTY = {"protocolSection": {}}


def test_trial_row_maps_fields():
    row = kernel.trial_row(STUDY_PHASE2)
    assert row["nct_id"] == "NCT01000001"
    assert row["phase"] == "PHASE1/PHASE2"
    assert row["sponsor"] == "Acme Pharma"
    assert row["enrollment"] == 300
    assert row["enrollment_type"] == "ACTUAL"
    assert row["arms"] == "Combo: Drug: Drug X + Drug: Pembrolizumab"
    assert row["interventions"] == "Drug X; Pembrolizumab"
    assert row["has_results"] is True
    assert row["url"] == "https://clinicaltrials.gov/study/NCT01000001"


def test_trial_row_terminated_carries_why_stopped():
    row = kernel.trial_row(STUDY_TERMINATED)
    assert row["overall_status"] == "TERMINATED"
    assert row["why_stopped"] == "Development change"
    assert row["enrollment_type"] == "ESTIMATED"


def test_trial_row_unphased_has_null_phase():
    row = kernel.trial_row(STUDY_UNPHASED)
    assert row["phase"] is None
    assert row["study_type"] == "OBSERVATIONAL"


def test_trial_row_empty_study_gives_no_url():
    row = kernel.trial_row(STUDY_EMPTY)
    assert row["nct_id"] is None
    assert row["url"] is None
    assert row["why_stopped"] is None


def test_trials_table_columns_and_rows():
    df = kernel.trials_table([STUDY_PHASE2, STUDY_TERMINATED])
    assert len(df) == 2
    assert list(df.columns) == kernel.TRIAL_COLUMNS


def test_trials_table_empty():
    df = kernel.trials_table([])
    assert len(df) == 0
    assert list(df.columns) == kernel.TRIAL_COLUMNS


def test_phase_rank_orders_phases():
    assert kernel.phase_rank("PHASE4") > kernel.phase_rank("PHASE3")
    assert kernel.phase_rank("PHASE1/PHASE2") == kernel.phase_rank("PHASE2")
    assert kernel.phase_rank("EARLY_PHASE1") < kernel.phase_rank("PHASE1")
    assert kernel.phase_rank("NA") is None
    assert kernel.phase_rank(None) is None


def test_rank_trials_orders_by_phase_and_keeps_all_rows():
    df = kernel.trials_table([STUDY_UNPHASED, STUDY_PHASE2, STUDY_TERMINATED])
    ranked = kernel.rank_trials(df)
    assert len(ranked) == 3
    assert list(ranked["nct_id"]) == ["NCT01000002", "NCT01000001",
                                      "NCT01000003"]
    assert ranked["phase_rank"].isna().iloc[-1]  # unphased sorts last


def test_phase_summary_counts_and_keeps_unphased():
    df = kernel.trials_table([STUDY_PHASE2, STUDY_TERMINATED, STUDY_UNPHASED])
    summary = kernel.phase_summary(df)
    assert list(summary.columns) == ["phase", "overall_status", "n_trials"]
    assert len(summary) == 3
    assert summary["n_trials"].sum() == 3
    assert summary.iloc[0]["phase"] == "PHASE3"  # highest phase first


def test_terminated_table_subsets_stopped():
    df = kernel.trials_table([STUDY_PHASE2, STUDY_TERMINATED])
    stopped = kernel.terminated_table(df)
    assert len(stopped) == 1
    assert stopped.iloc[0]["nct_id"] == "NCT01000002"
    assert stopped.iloc[0]["why_stopped"] == "Development change"


VALID_SIGNAL = {
    "signal_id": "S1", "drug": "Drug X",
    "signal_type": "biomarker_requirement", "value": "KRAS G12C mutation",
    "evidence_text": "Inclusion: KRAS G12C mutation.",
    "nct_id": "NCT01000001",
    "url": "https://clinicaltrials.gov/study/NCT01000001",
}


def test_signals_table_valid_row():
    df = kernel.signals_table([VALID_SIGNAL])
    assert len(df) == 1
    assert list(df.columns) == kernel.SIGNAL_COLUMNS


def test_signals_table_empty():
    df = kernel.signals_table([])
    assert len(df) == 0
    assert list(df.columns) == kernel.SIGNAL_COLUMNS


def test_signals_table_missing_field_raises():
    bad = {**VALID_SIGNAL, "evidence_text": ""}
    with pytest.raises(ValueError, match="missing required field"):
        kernel.signals_table([bad])


def test_signals_table_bad_signal_type_raises():
    bad = {**VALID_SIGNAL, "signal_type": "vibes"}
    with pytest.raises(ValueError, match="not in"):
        kernel.signals_table([bad])


SAMPLE_STEPS = [
    {"name": "Trials", "finding": "Three trials found; highest phase 3."},
    {"name": "Terminations",
     "finding": "One phase 3 trial terminated (development change)."},
]


def test_write_run(tmp_path):
    written = kernel.write_run([STUDY_PHASE2, STUDY_TERMINATED],
                               [VALID_SIGNAL], tmp_path, "Drug X",
                               "Summary.", SAMPLE_STEPS)
    assert set(written) == {"trials", "phase_summary", "terminated",
                            "signals", "readme"}
    for path in written.values():
        assert path.is_file()
    readme = written["readme"].read_text()
    assert "# Drug X — ClinicalTrials.gov" in readme


def test_write_run_bad_signal_stops_before_signals_csv(tmp_path):
    # trials.csv etc. are written before signals are validated, matching
    # the single-writer/table-per-call pattern used elsewhere in this repo.
    bad = {**VALID_SIGNAL, "signal_type": "vibes"}
    with pytest.raises(ValueError):
        kernel.write_run([STUDY_PHASE2], [bad], tmp_path, "Drug X",
                         "Summary.", SAMPLE_STEPS)
    assert (tmp_path / "trials.csv").is_file()
    assert not (tmp_path / "population_signals.csv").is_file()
