import math

import pytest

import kernel

NULL_TARGET = {"target": None}

SAMPLE = {
    "target": {
        "id": "ENSG00000165392",
        "approvedSymbol": "WRN",
        "approvedName": "WRN RecQ like helicase",
        "biotype": "protein_coding",
        "tractability": [
            {"label": "Approved Drug", "modality": "SM", "value": False},
            {"label": "Advanced Clinical", "modality": "AB", "value": True},
        ],
        "drugAndClinicalCandidates": {
            "count": 1,
            "rows": [{
                "maxClinicalStage": 1,
                "drug": {
                    "id": "CHEMBL123",
                    "name": "Example Drug",
                    "drugType": "Small molecule",
                    "mechanismsOfAction": {
                        "rows": [
                            {"actionType": "INHIBITOR", "mechanismOfAction": "WRN helicase inhibitor"},
                        ],
                    },
                },
                "diseases": [
                    {"disease": {"id": "EFO_0000305", "name": "breast carcinoma"}},
                ],
            }],
        },
        "associatedDiseases": {
            "count": 2,
            "rows": [
                {
                    "score": 0.8,
                    "disease": {"id": "EFO_0000305", "name": "breast carcinoma"},
                    "datatypeScores": [
                        {"id": "genetic_association", "score": 0.5},
                        {"id": "known_drug", "score": 0.9},
                    ],
                },
                {
                    "score": 0.3,
                    "disease": {"id": "EFO_0000311", "name": "cancer"},
                    "datatypeScores": [
                        {"id": "literature", "score": 0.3},
                    ],
                },
            ],
        },
        "safetyLiabilities": [
            {"event": "cardiotoxicity", "eventId": "EFO_123", "datasource": "AOP-Wiki",
             "literature": "12345", "url": "https://example.org/12345"},
        ],
    }
}


def test_identity_table():
    df = kernel.identity_table(SAMPLE)
    assert df.iloc[0]["ensembl_id"] == "ENSG00000165392"
    assert df.iloc[0]["approved_symbol"] == "WRN"
    assert "ENSG00000165392" in df.iloc[0]["open_targets_url"]


def test_known_drugs_table():
    df = kernel.known_drugs_table(SAMPLE)
    assert len(df) == 1
    row = df.iloc[0]
    assert row["drug_id"] == "CHEMBL123"
    assert row["max_clinical_stage"] == 1
    assert row["action_type"] == "INHIBITOR"
    assert row["mechanism_of_action"] == "WRN helicase inhibitor"
    assert row["disease_id"] == "EFO_0000305"


def test_known_drugs_table_multiple_diseases_one_row_each():
    data = {"target": {"drugAndClinicalCandidates": {"rows": [{
        "maxClinicalStage": 2,
        "drug": {"id": "CHEMBL1", "name": "D1", "drugType": "SM",
                 "mechanismsOfAction": {"rows": []}},
        "diseases": [
            {"disease": {"id": "EFO_1", "name": "disease one"}},
            {"disease": {"id": "EFO_2", "name": "disease two"}},
        ],
    }]}}}
    df = kernel.known_drugs_table(data)
    assert len(df) == 2
    assert set(df["disease_id"]) == {"EFO_1", "EFO_2"}


def test_known_drugs_table_null_disease_kept():
    data = {"target": {"drugAndClinicalCandidates": {"rows": [{
        "maxClinicalStage": 1,
        "drug": {"id": "CHEMBL1", "name": "D1", "drugType": "SM",
                 "mechanismsOfAction": {"rows": []}},
        "diseases": [{"disease": None}],
    }]}}}
    df = kernel.known_drugs_table(data)
    assert len(df) == 1
    assert df.iloc[0]["disease_id"] is None


def test_known_drugs_table_empty():
    empty = {"target": {"drugAndClinicalCandidates": {"count": 0, "rows": []}}}
    df = kernel.known_drugs_table(empty)
    assert len(df) == 0
    assert list(df.columns) == [
        "drug_id", "drug_name", "drug_type", "max_clinical_stage",
        "action_type", "mechanism_of_action",
        "disease_id", "disease_name", "drug_url",
    ]


def test_tractability_table_keeps_false_rows():
    df = kernel.tractability_table(SAMPLE)
    assert len(df) == 2
    assert df.iloc[0]["value"] == False  # noqa: E712 - pandas gives np.bool_, not bool


def test_associations_table_pivots_datatypes_and_fills_absent_with_nan():
    df = kernel.associations_table(SAMPLE)
    assert len(df) == 2
    assert set(df.columns) == {
        "disease_id", "disease_name", "overall_score", "disease_url",
        "genetic_association", "known_drug", "literature",
    }
    row = df[df["disease_id"] == "EFO_0000305"].iloc[0]
    assert row["genetic_association"] == 0.5
    assert row["known_drug"] == 0.9
    assert math.isnan(row["literature"])


def test_genetic_evidence_table():
    df = kernel.genetic_evidence_table(SAMPLE)
    assert len(df) == 2
    row = df[df["disease_id"] == "EFO_0000305"].iloc[0]
    assert row["genetic_association_score"] == 0.5
    row2 = df[df["disease_id"] == "EFO_0000311"].iloc[0]
    assert math.isnan(row2["genetic_association_score"])


def test_safety_table():
    df = kernel.safety_table(SAMPLE)
    assert len(df) == 1
    assert df.iloc[0]["event"] == "cardiotoxicity"
    assert df.iloc[0]["event_id"] == "EFO_123"
    assert df.iloc[0]["url"] == "https://example.org/12345"


def test_safety_table_empty():
    empty = {"target": {}}
    df = kernel.safety_table(empty)
    assert len(df) == 0
    assert list(df.columns) == ["event", "event_id", "datasource", "literature", "url"]


@pytest.mark.parametrize("table_fn", [
    kernel.identity_table,
    kernel.known_drugs_table,
    kernel.tractability_table,
    kernel.associations_table,
    kernel.genetic_evidence_table,
    kernel.safety_table,
])
def test_null_target_raises_clear_error(table_fn):
    with pytest.raises(ValueError, match="no target"):
        table_fn(NULL_TARGET)


SAMPLE_STEPS = [
    {"name": "Disease associations", "finding": "Two diseases associated.", "table": [
        {"disease": "breast carcinoma", "score": "0.80"},
    ]},
    {"name": "Safety", "skipped": "No liabilities reported in this dataset."},
]


def test_write_readme_assembles_structure(tmp_path):
    path = kernel.write_readme(tmp_path / "README.md", "WRN", "Summary text.", SAMPLE_STEPS,
                               source="Open Targets")
    text = path.read_text()
    assert "# WRN — Open Targets" in text
    assert "Summary text." in text
    assert "## Disease associations" in text
    assert "breast carcinoma" in text
    assert "## Safety" in text
    assert "Not run." in text


def test_write_readme_requires_summary(tmp_path):
    with pytest.raises(ValueError, match="summary is required"):
        kernel.write_readme(tmp_path / "README.md", "WRN", "", SAMPLE_STEPS)


def test_write_readme_requires_steps(tmp_path):
    with pytest.raises(ValueError, match="at least one step"):
        kernel.write_readme(tmp_path / "README.md", "WRN", "Summary.", [])


def test_write_readme_step_missing_finding_raises(tmp_path):
    steps = [{"name": "Drugs"}]
    with pytest.raises(ValueError, match="finding.*skipped"):
        kernel.write_readme(tmp_path / "README.md", "WRN", "Summary.", steps)


def test_write_readme_limits_section(tmp_path):
    text = kernel.write_readme(
        tmp_path / "README.md", "WRN", "Summary.", SAMPLE_STEPS,
        limits=["Capped at 20 associations."]
    ).read_text()
    assert "## Limits" in text
    assert "Capped at 20" in text


def test_markdown_table():
    rows = [{"a": 1, "b": 2}, {"a": 3, "b": 4}]
    out = kernel.markdown_table(rows)
    assert "| a | b |" in out
    assert "| 1 | 2 |" in out


def test_write_run(tmp_path):
    written = kernel.write_run(SAMPLE, tmp_path)
    assert set(written) == {
        "identity", "known_drugs", "tractability",
        "associations", "genetic_evidence", "safety",
    }
    for path in written.values():
        assert path.is_file()
