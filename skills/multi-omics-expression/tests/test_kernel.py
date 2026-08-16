import kernel

TISSUE_RECORDS = [
    {"tissueSiteDetailId": "Liver", "median": 40.0},
    {"tissueSiteDetailId": "Skin", "median": 10.0},
    {"tissueSiteDetailId": "Brain - Cortex", "median": 1.0},
    {"tissueSiteDetailId": "Muscle - Skeletal", "median": 10.0},
]


def test_classify_level_high():
    level, fold = kernel.classify_level(40.0, 10.0)
    assert level == "high"
    assert fold == 4.0


def test_classify_level_low():
    level, fold = kernel.classify_level(2.5, 10.0)
    assert level == "low"
    assert fold == 0.25


def test_classify_level_typical():
    level, fold = kernel.classify_level(10.0, 10.0)
    assert level == "typical"
    assert fold == 1.0


def test_classify_level_unknown_when_reference_zero():
    level, fold = kernel.classify_level(5.0, 0)
    assert level == "unknown"
    assert fold is None


def test_classify_level_unknown_when_value_missing():
    level, fold = kernel.classify_level(None, 10.0)
    assert level == "unknown"
    assert fold is None


def test_tissue_expression_table_shape_and_flags():
    df = kernel.tissue_expression_table(TISSUE_RECORDS)
    assert list(df.columns) == kernel.TISSUE_EXPRESSION_COLUMNS
    assert len(df) == 4
    # median of [40, 10, 1, 10] is 10 -> Liver at 40 is 4x -> high
    liver = df[df["tissue"] == "Liver"].iloc[0]
    assert liver["level"] == "high"
    brain = df[df["tissue"] == "Brain - Cortex"].iloc[0]
    assert brain["level"] == "low"
    skin = df[df["tissue"] == "Skin"].iloc[0]
    assert skin["level"] == "typical"


def test_tissue_expression_table_empty():
    df = kernel.tissue_expression_table([])
    assert len(df) == 0
    assert list(df.columns) == kernel.TISSUE_EXPRESSION_COLUMNS


STUDY_VALUES = {
    "brca_tcga_pan_can_atlas_2018": [40.0, 44.0, 36.0],
    "luad_tcga_pan_can_atlas_2018": [10.0, 9.0, 11.0],
    "gbm_tcga_pan_can_atlas_2018": [1.0, 1.2, 0.8],
}


def test_tumor_expression_table_shape_and_flags():
    df = kernel.tumor_expression_table(STUDY_VALUES)
    assert list(df.columns) == kernel.TUMOR_EXPRESSION_COLUMNS
    assert len(df) == 3
    brca = df[df["study_id"] == "brca_tcga_pan_can_atlas_2018"].iloc[0]
    assert brca["cancer_type"] == "Breast invasive carcinoma"
    assert brca["n_samples"] == 3
    assert brca["level"] == "high"
    gbm = df[df["study_id"] == "gbm_tcga_pan_can_atlas_2018"].iloc[0]
    assert gbm["level"] == "low"


def test_tumor_expression_table_empty():
    df = kernel.tumor_expression_table({})
    assert len(df) == 0
    assert list(df.columns) == kernel.TUMOR_EXPRESSION_COLUMNS


def test_tumor_expression_table_unknown_study_falls_back_to_id():
    df = kernel.tumor_expression_table({"made_up_study": [5.0, 5.0]})
    row = df.iloc[0]
    assert row["cancer_type"] == "made_up_study"


SAMPLE_STEPS = [
    {"name": "Tissue expression", "finding": "Liver is high (4.0x median)."},
    {"name": "Tumor expression", "finding": "GBM is low (0.1x median)."},
]


def test_write_run(tmp_path):
    written = kernel.write_run(TISSUE_RECORDS, STUDY_VALUES, tmp_path, "WRN",
                               "Summary.", SAMPLE_STEPS)
    assert set(written) == {"tissue_expression", "tumor_expression", "readme"}
    for path in written.values():
        assert path.is_file()
    readme = written["readme"].read_text()
    assert "# WRN — GTEx / cBioPortal" in readme
