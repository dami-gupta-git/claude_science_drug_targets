import sys
from pathlib import Path

import pandas as pd
import pytest

import depmap_ot_fusion_kernel as kernel


# ---------------------------------------------------------------------------
# read_depmap_run / read_open_targets_run
# ---------------------------------------------------------------------------

def test_read_depmap_run_required_frames(wrn_depmap_dir):
    data = kernel.read_depmap_run(wrn_depmap_dir)
    assert "selectivity" in data["frames"]
    assert "lineage_enrichment" in data["frames"]


def test_read_depmap_run_optional_files(wrn_depmap_dir):
    data = kernel.read_depmap_run(wrn_depmap_dir)
    assert "msi_contrast.csv" in data["optional_files"]
    assert "fusion_contrast.csv" in data["optional_files"]
    assert "prism_compounds.csv" not in data["optional_files"]


def test_read_depmap_run_prism_present(kras_depmap_dir):
    data = kernel.read_depmap_run(kras_depmap_dir)
    assert "prism_compounds.csv" in data["optional_files"]
    assert "compound_id" in data["frames"]["prism_compounds"].columns


def test_read_open_targets_run(wrn_ot_dir):
    data = kernel.read_open_targets_run(wrn_ot_dir)
    assert "identity" in data["frames"]
    assert len(data["frames"]["identity"]) == 1
    assert "known_drugs" in data["frames"]


# ---------------------------------------------------------------------------
# Completeness checks
# ---------------------------------------------------------------------------

def test_depmap_run_is_complete(wrn_depmap_dir):
    assert kernel.depmap_run_is_complete(wrn_depmap_dir)


def test_depmap_run_is_incomplete_missing_file(tmp_path, wrn_depmap_files):
    d = tmp_path / "run"
    d.mkdir()
    (d / "selectivity.csv").write_text(wrn_depmap_files["selectivity.csv"])
    assert not kernel.depmap_run_is_complete(d)


def test_depmap_run_is_incomplete_empty_file(tmp_path):
    d = tmp_path / "run"
    d.mkdir()
    (d / "selectivity.csv").write_text("gene,classification\n")
    (d / "lineage_enrichment.csv").write_text("lineage,n\n")
    assert not kernel.depmap_run_is_complete(d)


def test_ot_run_is_complete(wrn_ot_dir):
    assert kernel.open_targets_run_is_complete(wrn_ot_dir, "WRN")


def test_ot_run_is_complete_case_insensitive(wrn_ot_dir):
    assert kernel.open_targets_run_is_complete(wrn_ot_dir, "wrn")


def test_ot_run_is_incomplete_missing_file(tmp_path, wrn_ot_files):
    d = tmp_path / "run"
    d.mkdir()
    for fname in ("identity.csv", "associations.csv", "tractability.csv"):
        (d / fname).write_text(wrn_ot_files[fname])
    assert not kernel.open_targets_run_is_complete(d, "WRN")


# ---------------------------------------------------------------------------
# Run selection
# ---------------------------------------------------------------------------

def test_find_depmap_run_picks_complete(results_tree):
    found = kernel.find_depmap_run("wrn", results_base=results_tree)
    assert found is not None
    assert kernel.depmap_run_is_complete(found)


def test_find_depmap_run_skips_incomplete(results_tree):
    found = kernel.find_depmap_run("wrn", results_base=results_tree)
    assert found.name == "2026-01-01_120000"


def test_find_depmap_run_missing_gene(results_tree):
    assert kernel.find_depmap_run("unknown_gene", results_base=results_tree) is None


def test_find_ot_run_picks_fresh(results_tree):
    found = kernel.find_open_targets_run("wrn", max_age_days=90, results_base=results_tree)
    assert found is not None
    assert found.name == "2026-07-01_120000"


def test_find_ot_run_rejects_stale(results_tree, wrn_ot_files):
    # Only the stale run exists; require freshness of 1 day
    stale_only = results_tree.parent / "stale_only"
    d = stale_only / "open_targets" / "wrn" / "2025-01-01_120000"
    d.mkdir(parents=True)
    for fname, content in wrn_ot_files.items():
        (d / fname).write_text(content)
    found = kernel.find_open_targets_run("wrn", max_age_days=1, results_base=stale_only)
    assert found is None


# ---------------------------------------------------------------------------
# fuse — verdict
# ---------------------------------------------------------------------------

def test_fuse_wrn_verdict_one_row(wrn_depmap_dir, wrn_ot_dir):
    tables = kernel.fuse(wrn_depmap_dir, wrn_ot_dir)
    assert len(tables["verdict"]) == 1


def test_fuse_wrn_verdict_fields(wrn_depmap_dir, wrn_ot_dir):
    v = kernel.fuse(wrn_depmap_dir, wrn_ot_dir)["verdict"].iloc[0]
    assert v["gene"] == "WRN"
    assert v["classification"] == "selective"
    assert v["top_lineage"] == "Ovary/Fallopian Tube"
    assert v["genotype_marker"] == "msi"
    assert v["max_clinical_stage"] is None or pd.isna(v["max_clinical_stage"])
    assert v["has_safety_liabilities"] == False


def test_fuse_kras_verdict_fields(kras_depmap_dir, kras_ot_dir):
    v = kernel.fuse(kras_depmap_dir, kras_ot_dir)["verdict"].iloc[0]
    assert v["gene"] == "KRAS"
    assert v["classification"] == "selective"
    assert v["genotype_marker"] == "mutation"
    assert v["max_clinical_stage"] == 4  # APPROVED


# ---------------------------------------------------------------------------
# fuse — conflicts
# ---------------------------------------------------------------------------

def test_fuse_wrn_produces_conflicts(wrn_depmap_dir, wrn_ot_dir):
    tables = kernel.fuse(wrn_depmap_dir, wrn_ot_dir)
    assert len(tables["conflicts"]) >= 1


def test_fuse_wrn_non_cancer_disease_conflict(wrn_depmap_dir, wrn_ot_dir):
    conflicts = kernel.fuse(wrn_depmap_dir, wrn_ot_dir)["conflicts"]
    types = conflicts["conflict_type"].tolist()
    assert "non_cancer_top_disease_with_selective_dependency" in types


def test_fuse_no_prism_no_drugs_no_prism_conflict(wrn_depmap_dir, wrn_ot_dir):
    # WRN has no PRISM data and no drugs; should NOT produce prism_activity conflict
    conflicts = kernel.fuse(wrn_depmap_dir, wrn_ot_dir)["conflicts"]
    assert "prism_activity_without_clinical_drug" not in conflicts["conflict_type"].tolist()


# ---------------------------------------------------------------------------
# fuse — disease_lineage
# ---------------------------------------------------------------------------

def test_fuse_disease_lineage_columns(wrn_depmap_dir, wrn_ot_dir):
    dl = kernel.fuse(wrn_depmap_dir, wrn_ot_dir)["disease_lineage"]
    for col in ["disease_id", "disease_name", "ot_overall_score",
                "depmap_lineage", "lineage_mean_effect", "lineage_p_bh"]:
        assert col in dl.columns


def test_fuse_gastric_carcinoma_maps_to_stomach_lineage(wrn_depmap_dir, wrn_ot_dir):
    dl = kernel.fuse(wrn_depmap_dir, wrn_ot_dir)["disease_lineage"]
    row = dl[dl["disease_name"] == "gastric carcinoma"]
    assert len(row) == 1
    assert row.iloc[0]["depmap_lineage"] == "Esophagus/Stomach"


def test_fuse_werner_syndrome_has_no_lineage(wrn_depmap_dir, wrn_ot_dir):
    dl = kernel.fuse(wrn_depmap_dir, wrn_ot_dir)["disease_lineage"]
    row = dl[dl["disease_name"] == "Werner syndrome"]
    assert len(row) == 1
    assert pd.isna(row.iloc[0]["depmap_lineage"]) or row.iloc[0]["depmap_lineage"] is None


# ---------------------------------------------------------------------------
# fuse — drugs_joined
# ---------------------------------------------------------------------------

def test_fuse_wrn_drugs_joined_empty_ot_drugs(wrn_depmap_dir, wrn_ot_dir):
    dj = kernel.fuse(wrn_depmap_dir, wrn_ot_dir)["drugs_joined"]
    # WRN has no known drugs and no PRISM data
    assert len(dj) == 0 or dj["drug_id"].isna().all()


def test_fuse_kras_drugs_joined_has_sotorasib(kras_depmap_dir, kras_ot_dir):
    dj = kernel.fuse(kras_depmap_dir, kras_ot_dir)["drugs_joined"]
    assert len(dj) > 0
    ot_drugs = dj[dj["drug_id"].notna()]
    assert "SOTORASIB" in ot_drugs["drug_name"].str.upper().tolist()


# ---------------------------------------------------------------------------
# write_fusion_run
# ---------------------------------------------------------------------------

def test_write_fusion_run_creates_files(wrn_depmap_dir, wrn_ot_dir, tmp_path):
    tables = kernel.fuse(wrn_depmap_dir, wrn_ot_dir)
    sources = {"depmap": {"reused": True, "timestamp": "2026-01-01_120000"},
               "open_targets": {"reused": False, "timestamp": "2026-07-01_120000"}}
    out = tmp_path / "out"
    out.mkdir()
    written = kernel.write_fusion_run(tables, sources, out)
    for name in ("verdict", "conflicts", "disease_lineage", "drugs_joined"):
        assert (out / f"{name}.csv").exists()
    assert (out / "sources.json").exists()
    import json
    src = json.loads((out / "sources.json").read_text())
    assert "depmap" in src
    assert "open_targets" in src


# ---------------------------------------------------------------------------
# write_fusion_readme
# ---------------------------------------------------------------------------

def test_write_fusion_readme_creates_file(wrn_depmap_dir, wrn_ot_dir, tmp_path):
    tables = kernel.fuse(wrn_depmap_dir, wrn_ot_dir)
    sources = {"depmap": {"reused": True, "timestamp": "2026-01-01_120000"},
               "open_targets": {"reused": False, "timestamp": "2026-07-01_120000"}}
    out = tmp_path / "out"
    out.mkdir()
    kernel.write_fusion_readme(out, "WRN", tables, sources)
    readme = out / "README.md"
    assert readme.exists()
    text = readme.read_text()
    assert "WRN" in text
    assert "selective" in text


# ---------------------------------------------------------------------------
# lineage_for_disease
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# lineage_status distinguishes the three empty-effect cases
# ---------------------------------------------------------------------------

def test_disease_lineage_status_matched(kras_depmap_dir, kras_ot_dir):
    """A disease whose lineage DepMap screened carries effect values."""
    dl = kernel.fuse(kras_depmap_dir, kras_ot_dir)["disease_lineage"]
    row = dl[dl["disease_name"] == "pancreatic carcinoma"].iloc[0]
    assert row["lineage_status"] == "matched"
    assert row["depmap_lineage"] == "Pancreas"
    assert row["lineage_mean_effect"] == -2.141


def test_disease_lineage_status_no_keyword(wrn_depmap_dir, wrn_ot_dir):
    """Werner syndrome is not a cancer and matches no lineage keyword."""
    dl = kernel.fuse(wrn_depmap_dir, wrn_ot_dir)["disease_lineage"]
    row = dl[dl["disease_name"] == "Werner syndrome"].iloc[0]
    assert row["lineage_status"] == "no_lineage_keyword"
    assert pd.isna(row["depmap_lineage"]) or row["depmap_lineage"] is None


def test_disease_lineage_status_lineage_absent_from_run(wrn_depmap_dir, wrn_ot_dir):
    """Gastric carcinoma maps to a lineage the WRN enrichment run does not cover."""
    dl = kernel.fuse(wrn_depmap_dir, wrn_ot_dir)["disease_lineage"]
    row = dl[dl["disease_name"] == "gastric carcinoma"].iloc[0]
    assert row["lineage_status"] == "lineage_not_in_depmap_run"
    assert row["depmap_lineage"] == "Esophagus/Stomach"
    assert pd.isna(row["lineage_mean_effect"])


def test_readme_omits_raw_score_precision(wrn_depmap_dir, wrn_ot_dir, tmp_path):
    """Association scores are rounded in prose; full precision stays in the CSV."""
    tables = kernel.fuse(wrn_depmap_dir, wrn_ot_dir)
    sources = {"depmap": {"reused": True, "timestamp": "2026-01-01_120000"},
               "open_targets": {"reused": True, "timestamp": "2026-01-01_120000"}}
    out = tmp_path / "out"
    out.mkdir()
    kernel.write_fusion_readme(out, "WRN", tables, sources)
    text = (out / "README.md").read_text()
    assert "0.828645145754454" not in text
    assert "0.83" in text


# ---------------------------------------------------------------------------
# sources.json records observed reuse
# ---------------------------------------------------------------------------

def _stub_fusion_run(monkeypatch, tmp_path, wrn_depmap_dir, wrn_ot_dir,
                     depmap_cached, ot_cached):
    """Run run_depmap_ot_fusion with upstream lookup and generation stubbed.

    depmap_cached / ot_cached control whether find_* reports an existing run,
    which is what distinguishes a reused source from a freshly generated one.
    """
    monkeypatch.setattr(kernel, "find_depmap_run",
                        lambda *a, **k: wrn_depmap_dir if depmap_cached else None)
    monkeypatch.setattr(kernel, "find_open_targets_run",
                        lambda *a, **k: wrn_ot_dir if ot_cached else None)
    monkeypatch.setattr(kernel, "ensure_depmap_run",
                        lambda *a, **k: wrn_depmap_dir)
    monkeypatch.setattr(kernel, "ensure_open_targets_run",
                        lambda *a, **k: wrn_ot_dir)
    out = tmp_path / "fusion_out"
    out.mkdir()
    monkeypatch.setattr(kernel, "results_dir", lambda *a, **k: out)
    kernel.run_depmap_ot_fusion("WRN")
    import json
    return json.loads((out / "sources.json").read_text())


def test_sources_records_reuse_when_runs_cached(monkeypatch, tmp_path,
                                                wrn_depmap_dir, wrn_ot_dir):
    src = _stub_fusion_run(monkeypatch, tmp_path, wrn_depmap_dir, wrn_ot_dir,
                           depmap_cached=True, ot_cached=True)
    assert src["depmap"]["reused"] is True
    assert src["open_targets"]["reused"] is True


def test_sources_records_generation_when_nothing_cached(monkeypatch, tmp_path,
                                                        wrn_depmap_dir, wrn_ot_dir):
    """reuse=True with no cached run generates fresh data; reused must be False."""
    src = _stub_fusion_run(monkeypatch, tmp_path, wrn_depmap_dir, wrn_ot_dir,
                           depmap_cached=False, ot_cached=False)
    assert src["depmap"]["reused"] is False
    assert src["open_targets"]["reused"] is False


def test_sources_records_each_source_independently(monkeypatch, tmp_path,
                                                   wrn_depmap_dir, wrn_ot_dir):
    """A cached DepMap run beside a stale Open Targets run is recorded as such."""
    src = _stub_fusion_run(monkeypatch, tmp_path, wrn_depmap_dir, wrn_ot_dir,
                           depmap_cached=True, ot_cached=False)
    assert src["depmap"]["reused"] is True
    assert src["open_targets"]["reused"] is False


@pytest.mark.parametrize("disease,expected_lineage", [
    ("colorectal carcinoma", "Bowel"),
    ("pancreatic carcinoma", "Pancreas"),
    ("uveal melanoma", "Eye"),
    ("cutaneous melanoma", "Skin"),
    ("gastric carcinoma", "Esophagus/Stomach"),
    ("glioblastoma", "CNS/Brain"),
    ("Werner syndrome", None),
    ("cancer", None),
])
def test_lineage_for_disease(disease, expected_lineage):
    assert kernel.lineage_for_disease(disease) == expected_lineage


def test_lineage_uveal_melanoma_beats_skin_melanoma():
    # 'uveal melanoma' must map to Eye, not Skin
    assert kernel.lineage_for_disease("uveal melanoma") == "Eye"


# ---------------------------------------------------------------------------
# stage_to_num
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("stage,expected", [
    ("PHASE_1", 1),
    ("PHASE_3", 3),
    ("APPROVED", 4),
    # Open Targets' live maxClinicalStage field uses "APPROVAL", not
    # "APPROVED" -- confirmed against real known_drugs.csv output. The
    # original code checked only for "APPROVED", which never matched, so
    # every approved drug fell through to the digit regex and returned None
    # instead of 4; verdict.csv's max_clinical_stage silently understated
    # drug-development status for any gene with an approved therapy.
    ("APPROVAL", 4),
    ("approval", 4),
    (None, None),
    ("UNKNOWN", None),
])
def test_stage_to_num(stage, expected):
    assert kernel.stage_to_num(stage) == expected
