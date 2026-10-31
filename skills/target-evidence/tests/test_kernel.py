import math

import pandas as pd
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
                "maxClinicalStage": "PHASE_1",
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
    assert row["max_clinical_stage"] == "PHASE_1"
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
        "parent_drug_id", "parent_drug_name",
    ]


def drug_row(drug_id, name, stage, parent=None, disease_id="EFO_1", moa="JAK inhibitor"):
    return {
        "maxClinicalStage": stage,
        "drug": {"id": drug_id, "name": name, "drugType": "Small molecule",
                 "parentMolecule": parent,
                 "mechanismsOfAction": {"rows": [
                     {"actionType": "INHIBITOR", "mechanismOfAction": moa}]}},
        "diseases": [{"disease": {"id": disease_id, "name": disease_id}}],
    }


def drugs_data(*rows):
    return {"target": {"drugAndClinicalCandidates": {"rows": list(rows)}}}


def test_drug_summary_folds_salt_into_parent():
    # TYK2 shape: the salt carries APPROVAL, the parent PHASE_4.
    df = kernel.drug_summary_table(drugs_data(
        drug_row("CHEMBL221959", "TOFACITINIB", "PHASE_4", disease_id="EFO_1"),
        drug_row("CHEMBL2103743", "TOFACITINIB CITRATE", "APPROVAL",
                 parent={"id": "CHEMBL221959", "name": "TOFACITINIB"}, disease_id="EFO_2"),
    ))
    assert len(df) == 1
    row = df.iloc[0]
    assert row["drug_id"] == "CHEMBL221959"
    assert row["max_clinical_stage"] == "APPROVED"
    assert row["n_indications"] == 2
    assert row["forms"] == "CHEMBL221959; CHEMBL2103743"


def test_drug_summary_keeps_deuterated_drug_separate():
    df = kernel.drug_summary_table(drugs_data(
        drug_row("CHEMBL1789941", "RUXOLITINIB", "PHASE_4"),
        drug_row("CHEMBL4594382", "DEURUXOLITINIB PHOSPHATE", "APPROVAL",
                 parent={"id": "CHEMBL1789941", "name": "RUXOLITINIB"}),
    ))
    assert set(df["drug_name"]) == {"RUXOLITINIB", "DEURUXOLITINIB PHOSPHATE"}


def test_drug_summary_salt_without_listed_parent_uses_parent_identity():
    df = kernel.drug_summary_table(drugs_data(
        drug_row("CHEMBL4298167", "FILGOTINIB MALEATE", "APPROVAL",
                 parent={"id": "CHEMBL3301607", "name": "FILGOTINIB"}),
    ))
    assert df.iloc[0]["drug_id"] == "CHEMBL3301607"
    assert df.iloc[0]["forms"] == "CHEMBL4298167"


def test_drug_summary_orders_approved_first():
    df = kernel.drug_summary_table(drugs_data(
        drug_row("C1", "ALPHA", "PHASE_2"),
        drug_row("C2", "BETA", "APPROVAL"),
        drug_row("C3", "GAMMA", "PHASE_4"),
        drug_row("C4", "DELTA", "UNKNOWN"),
    ))
    assert list(df["drug_name"]) == ["BETA", "GAMMA", "ALPHA", "DELTA"]


def test_drug_summary_phase_4_only_counts_as_approved():
    # TYK2 shape: deucravacitinib has only post-marketing (PHASE_4) rows.
    df = kernel.drug_summary_table(drugs_data(drug_row("C1", "ALPHA", "PHASE_4")))
    assert df.iloc[0]["max_clinical_stage"] == "APPROVED"
    assert df.iloc[0]["all_stages"] == "PHASE_4"


def test_drug_summary_unranked_stages_kept_visible():
    df = kernel.drug_summary_table(drugs_data(
        drug_row("C1", "ALPHA", "WITHDRAWAL", disease_id="EFO_1"),
        drug_row("C1", "ALPHA", "PHASE_3", disease_id="EFO_2"),
        drug_row("C2", "BETA", "UNKNOWN"),
    ))
    by_name = df.set_index("drug_name")
    assert by_name.loc["ALPHA", "max_clinical_stage"] == "PHASE_3"
    assert by_name.loc["ALPHA", "all_stages"] == "PHASE_3; WITHDRAWAL"
    assert pd.isna(by_name.loc["BETA", "max_clinical_stage"])


def test_drug_summary_unrecognised_stage_raises():
    with pytest.raises(ValueError, match="Unrecognised"):
        kernel.drug_summary_table(drugs_data(drug_row("C1", "ALPHA", "PHASE_9")))


def test_drug_summary_empty():
    df = kernel.drug_summary_table(drugs_data())
    assert len(df) == 0
    assert "max_clinical_stage" in df.columns


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


def test_top_association_genetic_scores_table_drops_diseases_without_score():
    df = kernel.top_association_genetic_scores_table(SAMPLE)
    assert list(df["disease_id"]) == ["EFO_0000305"]
    assert df.iloc[0]["genetic_association_score"] == 0.5


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
    kernel.drug_summary_table,
    kernel.tractability_table,
    kernel.associations_table,
    kernel.top_association_genetic_scores_table,
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


def test_markdown_table_formats_floats():
    rows = [{"p": 2.0763817748514768e-11, "q": 0.819959670562054,
             "x": 17.36, "flag": True}]
    out = kernel.markdown_table(rows)
    assert "| 2.08e-11 | 0.82 | 17.36 | True |" in out


def test_write_run(tmp_path):
    written = kernel.write_run(SAMPLE, tmp_path, "WRN", "Summary.", SAMPLE_STEPS)
    assert set(written) == {
        "identity", "known_drugs", "drug_summary", "tractability",
        "associations", "top_association_genetic_scores", "safety",
        "association_coverage", "readme",
    }
    for path in written.values():
        assert path.is_file()


def test_write_run_records_association_truncation(tmp_path):
    # associations.csv is a capped top slice; the coverage file is what says so.
    written = kernel.write_run(SAMPLE, tmp_path, "WRN", "Summary.",
                               SAMPLE_STEPS)
    coverage = pd.read_csv(written["association_coverage"])
    assert list(coverage.columns) == ["n_returned", "n_total", "truncated",
                                      "lowest_returned_score"]
    assert len(coverage) == 1


def test_association_count_flags_truncation():
    capped = {"target": {"associatedDiseases": {
        "count": 1002, "rows": [{"score": 0.5}] * 20}}}
    out = kernel.association_count(capped)
    assert out["truncated"] is True
    assert out["n_total"] == 1002 and out["n_returned"] == 20
    complete = {"target": {"associatedDiseases": {
        "count": 3, "rows": [{"score": 0.5}] * 3}}}
    assert kernel.association_count(complete)["truncated"] is False


# --- gnomAD -----------------------------------------------------------------

GNOMAD_GENE = {
    "gene_id": "ENSG00000165392",
    "symbol": "WRN",
    "dataset": "gnomad_r4",
    "gnomad_constraint": {
        "oe_lof": 0.82, "oe_lof_lower": 0.72, "oe_lof_upper": 0.94,
        "obs_lof": 149, "exp_lof": 180.7, "pLI": 0.0, "lof_z": -1.2,
        "oe_mis": 0.98, "mis_z": 0.1, "flags": [],
    },
    "variants": [
        {"variant_id": "8-1-A-T", "consequence": "frameshift_variant",
         "hgvsp": "p.Ser2fs", "lof": "HC", "lof_filter": None,
         "flags": ["lof_flag"],
         "exome": {"ac": 5, "an": 1000, "homozygote_count": 1},
         "genome": {"ac": 1, "an": 500, "homozygote_count": 0}},
        {"variant_id": "8-2-G-C", "consequence": "missense_variant",
         "hgvsp": "p.Ala3Gly", "lof": None, "lof_filter": None, "flags": [],
         "exome": {"ac": 90, "an": 1000, "homozygote_count": 4},
         "genome": None},
        {"variant_id": "8-3-C-T", "consequence": "stop_gained",
         "hgvsp": "p.Arg9Ter", "lof": "LC", "lof_filter": "END_TRUNC",
         "flags": [], "exome": {"ac": 2, "an": 1000, "homozygote_count": 0},
         "genome": None},
    ],
}


def test_gnomad_constraint_table():
    row = kernel.gnomad_constraint_table(GNOMAD_GENE).iloc[0]
    assert row["loeuf"] == 0.94
    assert row["obs_lof"] == 149
    assert "ENSG00000165392" in row["gnomad_url"]


def test_gnomad_constraint_table_raises_without_constraint():
    gene = {**GNOMAD_GENE, "gnomad_constraint": None}
    with pytest.raises(ValueError, match="unmeasured"):
        kernel.gnomad_constraint_table(gene)


def test_gnomad_lof_variants_keeps_high_confidence_only():
    table = kernel.gnomad_lof_variants_table(GNOMAD_GENE)
    assert list(table.variant_id) == ["8-1-A-T"]
    assert table.iloc[0]["allele_count"] == 6
    assert table.iloc[0]["allele_number"] == 1500
    assert table.iloc[0]["homozygote_count"] == 1


def test_gnomad_missense_variants_filters_by_frequency():
    table = kernel.gnomad_missense_variants_table(GNOMAD_GENE, 0.05)
    assert list(table.variant_id) == ["8-2-G-C"]
    assert table.iloc[0]["allele_frequency"] == 0.09
    assert table.iloc[0]["homozygote_count"] == 4
    assert not table.iloc[0]["reference_is_minor"]
    assert kernel.gnomad_missense_variants_table(GNOMAD_GENE, 0.1).empty


def test_gnomad_missense_variants_uses_minor_allele_when_reference_is_rare():
    variant = {"variant_id": "8-4-C-A", "consequence": "missense_variant",
               "hgvsp": "p.Leu5Met", "lof": None, "lof_filter": None, "flags": [],
               "exome": {"ac": 910, "an": 1000, "homozygote_count": 820},
               "genome": None}
    table = kernel.gnomad_missense_variants_table(
        {**GNOMAD_GENE, "variants": [variant]}, 0.05)
    row = table.iloc[0]
    assert row["allele_frequency"] == 0.91
    assert math.isclose(row["minor_allele_frequency"], 0.09)
    assert row["reference_is_minor"]
    assert kernel.gnomad_missense_variants_table(
        {**GNOMAD_GENE, "variants": [variant]}, 0.1).empty


def test_gnomad_missense_variants_rejects_bad_frequency():
    for bad in (0, -0.01, 0.6):
        with pytest.raises(ValueError, match="min_allele_frequency"):
            kernel.gnomad_missense_variants_table(GNOMAD_GENE, bad)


def test_write_run_requires_missense_frequency_with_gnomad(tmp_path):
    with pytest.raises(ValueError, match="missense_min_allele_frequency"):
        kernel.write_run(SAMPLE, tmp_path, "WRN", "Summary.", SAMPLE_STEPS,
                         gnomad=GNOMAD_GENE)


def test_gnomad_carrier_summary():
    row = kernel.gnomad_carrier_summary_table(GNOMAD_GENE).iloc[0]
    assert row["n_lof_variants"] == 1
    assert row["total_lof_allele_count"] == 6
    assert row["total_homozygotes"] == 1
    assert row["loeuf"] == 0.94


def test_gnomad_carrier_summary_with_no_lof_variants():
    row = kernel.gnomad_carrier_summary_table({**GNOMAD_GENE, "variants": []}).iloc[0]
    assert row["n_lof_variants"] == 0
    assert row["max_allele_frequency"] is None


# --- cBioPortal -------------------------------------------------------------

RANGES = {"helicase": (558, 724)}

MUTATIONS = [
    {"studyId": "s1", "sampleId": "S-1", "patientId": "P-1",
     "proteinChange": "R600C", "mutationType": "Missense_Mutation",
     "proteinPosStart": 600, "tumorAltCount": 10, "tumorRefCount": 90},
    {"studyId": "s1", "sampleId": "S-2", "patientId": "P-2",
     "proteinChange": "L44I", "mutationType": "Missense_Mutation",
     "proteinPosStart": 44, "tumorAltCount": 20, "tumorRefCount": 20},
    {"studyId": "s1", "sampleId": "S-9", "patientId": "P-9",
     "proteinChange": "R700H", "mutationType": "Missense_Mutation",
     "proteinPosStart": 700, "tumorAltCount": None, "tumorRefCount": None},
]


def test_residue_range_hit():
    assert kernel.residue_range_hit(600, RANGES) == "helicase"
    assert kernel.residue_range_hit(558, RANGES) == "helicase"
    assert kernel.residue_range_hit(724, RANGES) == "helicase"
    assert kernel.residue_range_hit(44, RANGES) is None
    assert kernel.residue_range_hit(None, RANGES) is None
    assert kernel.residue_range_hit("", RANGES) is None


def test_msi_high_samples_filters_and_ignores_non_numeric():
    rows = [{"sampleId": "S-1", "value": "25.0"},
            {"sampleId": "S-2", "value": "3.1"},
            {"sampleId": "S-3", "value": "NA"}]
    assert kernel.msi_high_samples(rows, 10.0) == {"S-1": 25.0}


def test_tumour_mutations_table_restricts_to_scored_msi_high_samples():
    table = kernel.tumour_mutations_table(
        MUTATIONS, {"S-1": 25.0, "S-2": 30.0}, RANGES, "MSI_SENSOR_SCORE")
    assert list(table.sample_id) == ["S-1", "S-2"]
    assert list(table.residue_range) == ["helicase", None]
    assert table.iloc[0]["variant_allele_frequency"] == 0.1


def test_tumour_mutations_table_leaves_vaf_empty_without_read_counts():
    table = kernel.tumour_mutations_table(
        MUTATIONS, {"S-9": 40.0}, RANGES, "MSI_SENSOR_SCORE")
    assert table.iloc[0]["variant_allele_frequency"] is None


def test_write_run_adds_supplied_sections(tmp_path):
    cohort = {"mutations": kernel.tumour_mutations_table(
                  MUTATIONS, {"S-1": 25.0}, RANGES, "MSI_SENSOR_SCORE"),
              "cohort": kernel.cbioportal_cohort_table([
                  {"study_id": "s1", "attribute_id": "MSI_SENSOR_SCORE",
                   "threshold": 10.0, "n_samples_scored": 100,
                   "n_msi_high_samples": 12, "n_msi_high_mutated": 1,
                   "n_msi_high_in_range": 1}])}
    written = kernel.write_run(SAMPLE, tmp_path, "WRN", "Summary.", SAMPLE_STEPS,
                               gnomad=GNOMAD_GENE, missense_min_allele_frequency=0.05,
                               cbioportal=cohort)
    assert "gnomad_constraint" in written
    assert "gnomad_missense_variants" in written
    assert "tumour_cohort" in written
    assert "msi_dose_response" not in written
    for path in written.values():
        assert path.is_file()
