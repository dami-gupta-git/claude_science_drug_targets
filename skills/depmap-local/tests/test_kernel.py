import sys

import pytest

import kernel


def make_root(tmp_path):
    root = tmp_path / "depmap_root"
    root.mkdir()
    (root / "CRISPRGeneEffect.csv").write_text(
        "ModelID,GeneA (123),GeneB (456)\n"
        "ACH-1,-1.2,0.1\n"
        "ACH-2,0.05,-0.9\n"
        "ACH-3,,-0.3\n"
    )
    (root / "OmicsSomaticMutationsMatrixHotspot.csv").write_text(
        "ModelID,IsDefaultEntryForModel,GeneA (123)\n"
        "ACH-1,Yes,1\n"
        "ACH-1,No,1\n"
        "ACH-2,Yes,0\n"
    )
    (root / "CRISPRInferredCommonEssentials.csv").write_text(
        "Gene\nGeneA (123)\nGeneB (456)\n"
    )
    return root


def test_depmap_root_requires_env_or_arg(monkeypatch):
    monkeypatch.delenv("DEPMAP_ROOT", raising=False)
    with pytest.raises(FileNotFoundError, match="No DepMap root configured"):
        kernel.depmap_root()


def test_depmap_root_requires_existing_dir():
    with pytest.raises(FileNotFoundError, match="not found"):
        kernel.depmap_root("/no/such/path")


def test_depmap_root_accepts_explicit_arg(tmp_path):
    assert kernel.depmap_root(str(tmp_path)) == str(tmp_path)


@pytest.mark.parametrize("raw,clean", [
    ("KRAS (3845)", "KRAS"),
    ("KRAS", "KRAS"),
    ("GENE (1)", "GENE"),
])
def test_strip_entrez(raw, clean):
    assert kernel.strip_entrez(raw) == clean


def test_depmap_inventory_lists_files_and_cache_flag(tmp_path):
    root = make_root(tmp_path)
    inv = kernel.depmap_inventory(str(root))
    assert "CRISPRGeneEffect.csv" in inv
    assert inv["_cache_built"] is False


def test_depmap_matrix_backend_prefers_parquet_then_csv(tmp_path):
    root = make_root(tmp_path)
    cache = root / kernel.CACHE_DIRNAME
    cache.mkdir()
    pq = cache / "crispr_gene_effect.parquet"
    pq.write_text("stub")

    backend, path = kernel.depmap_matrix_backend("effect", root=str(root))
    assert backend == "parquet"
    assert path == str(pq)

    pq.unlink()
    backend, path = kernel.depmap_matrix_backend("effect", root=str(root))
    assert backend == "csv"
    assert path.endswith("CRISPRGeneEffect.csv")


def test_depmap_matrix_backend_raises_when_no_source(tmp_path):
    root = tmp_path / "empty"
    root.mkdir()
    with pytest.raises(FileNotFoundError):
        kernel.depmap_matrix_backend("effect", root=str(root))


def test_depmap_matrix_backend_unknown_kind(tmp_path):
    root = make_root(tmp_path)
    with pytest.raises(ValueError, match="unknown matrix"):
        kernel.depmap_matrix_backend("bogus", root=str(root))


def test_depmap_read_matrix_csv_keeps_order_and_reports_missing(tmp_path):
    root = make_root(tmp_path)
    df, missing = kernel.depmap_read_matrix("effect", ["GeneB", "GeneA", "GeneC"], root=str(root))
    assert list(df.columns) == ["GeneB", "GeneA"]
    assert missing == ["GeneC"]
    assert df.index.name == "ModelID"
    assert df.loc["ACH-1", "GeneA"] == -1.2


def test_depmap_read_matrix_raises_when_all_missing(tmp_path):
    root = make_root(tmp_path)
    with pytest.raises(KeyError):
        kernel.depmap_read_matrix("effect", ["NoSuchGene"], root=str(root))


def test_depmap_read_matrix_collapses_to_default_entry(tmp_path):
    root = make_root(tmp_path)
    df, missing = kernel.depmap_read_matrix("hotspot", ["GeneA"], root=str(root))
    assert missing == []
    assert set(df.index) == {"ACH-1", "ACH-2"}
    assert df.loc["ACH-1", "GeneA"] == 1


def test_depmap_common_essentials_strips_entrez(tmp_path):
    root = make_root(tmp_path)
    assert kernel.depmap_common_essentials(str(root)) == {"GeneA", "GeneB"}


def test_depmap_prism_releases_sorted(tmp_path):
    root = make_root(tmp_path)
    (root / "Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv").write_text("x")
    (root / "Repurposing_Public_23Q2_Extended_Primary_Data_Matrix.csv").write_text("x")
    assert kernel.depmap_prism_releases(str(root)) == ["23Q2", "24Q2"]


def test_depmap_prism_releases_empty(tmp_path):
    root = make_root(tmp_path)
    assert kernel.depmap_prism_releases(str(root)) == []


def test_floor_pvalue_clamps_zero():
    assert kernel.floor_pvalue(0.0) == sys.float_info.min


def test_floor_pvalue_passes_through_nonzero():
    assert kernel.floor_pvalue(0.5) == 0.5


def test_bh_adjust_equal_pvalues():
    out = kernel.bh_adjust([0.01, 0.02, 0.03])
    assert out == pytest.approx([0.03, 0.03, 0.03])


def test_bh_adjust_aligns_to_input_order():
    out = kernel.bh_adjust([0.001, 0.2, 0.05, 0.01])
    assert out == pytest.approx([0.004, 0.2, 0.0667, 0.02], abs=1e-4)


def test_write_run_writes_csvs(tmp_path):
    import pandas as pd
    tables = {"foo": pd.DataFrame([{"a": 1}])}
    written = kernel.write_run(tables, tmp_path)
    assert set(written) == {"foo"}
    assert written["foo"].is_file()
    assert written["foo"].read_text().strip().splitlines() == ["a", "1"]
