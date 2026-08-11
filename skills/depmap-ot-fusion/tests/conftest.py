import importlib.util
import sys
from pathlib import Path

# Load this skill's kernel under a unique module name so it does not collide
# with 'kernel' modules registered by other skills' conftest.py files.
_skill_dir = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("depmap_ot_fusion_kernel",
                                                _skill_dir / "kernel.py")
_mod = importlib.util.module_from_spec(_spec)
sys.modules["depmap_ot_fusion_kernel"] = _mod
_spec.loader.exec_module(_mod)

import pytest

# ---------------------------------------------------------------------------
# Fixture CSV content (trimmed from real WRN / KRAS runs)
# ---------------------------------------------------------------------------

WRN_DEPMAP = {
    "selectivity.csv": (
        "gene,n_lines,mean_effect,median_effect,min_effect,pct05_effect,"
        "frac_dependent,n_dependent,sd_effect,frac_positive,p95_effect,classification\n"
        "WRN,1208,-0.18,-0.121,-4.18,-0.565,0.0679,82,0.342,0.0025,0.115,selective\n"
    ),
    "lineage_enrichment.csv": (
        "lineage,n,mean_effect,mean_effect_other,frac_dependent,p_enriched,p_bh\n"
        "Ovary/Fallopian Tube,59,-0.287,-0.175,0.119,0.001877,0.041\n"
        "Uterus,34,-0.573,-0.169,0.382,0.006772,0.074\n"
        "Bowel,63,-0.435,-0.166,0.190,0.361,0.739\n"
        "Lung,126,-0.117,-0.188,0.024,0.947,0.999\n"
    ),
    "msi_contrast.csv": (
        "gene,msi_threshold,n_msi_high,n_mss,mean_effect_msi_high,mean_effect_mss,"
        "frac_dependent_msi_high,frac_dependent_mss,p_msi_more_dependent\n"
        "WRN,20.0,75,1133,-0.854,-0.136,0.52,0.038,1.985e-12\n"
    ),
    "fusion_contrast.csv": (
        "gene,n_fusion_models,n_fusion_with_crispr,n_fusion_negative,"
        "n_unprofiled_excluded,canonical_only,mean_effect_fusion_pos,"
        "mean_effect_fusion_neg,frac_dependent_fusion_pos,p_fusion_more_dependent,status\n"
        "WRN,5,3,1137,68,True,-0.271,-0.184,0.0,0.0796,not-dependent-despite-fusion\n"
    ),
}

WRN_OT = {
    "identity.csv": (
        "ensembl_id,approved_symbol,approved_name,biotype,open_targets_url\n"
        "ENSG00000165392,WRN,WRN RecQ like helicase,protein_coding,"
        "https://platform.opentargets.org/target/ENSG00000165392\n"
    ),
    "associations.csv": (
        "disease_id,disease_name,overall_score,disease_url,"
        "affected_pathway,genetic_association,genetic_literature,literature,somatic_mutation\n"
        "MONDO_0010196,Werner syndrome,0.829,https://platform.opentargets.org/disease/MONDO_0010196,"
        ",0.886,0.846,0.527,\n"
        "MONDO_0004992,cancer,0.555,https://platform.opentargets.org/disease/MONDO_0004992,"
        "0.865,,,0.953,\n"
        "MONDO_0004950,gastric carcinoma,0.438,https://platform.opentargets.org/disease/MONDO_0004950,"
        "0.449,,,0.049,0.608\n"
    ),
    "tractability.csv": (
        "label,modality,value\n"
        "Approved Drug,SM,False\n"
        "Advanced Clinical,SM,False\n"
        "Structure with Ligand,SM,True\n"
        "High-Quality Ligand,SM,True\n"
    ),
    "known_drugs.csv": (
        "drug_id,drug_name,drug_type,max_clinical_stage,action_type,"
        "mechanism_of_action,disease_id,disease_name,drug_url\n"
    ),
    "safety.csv": "event,event_id,datasource,literature,url\n",
    "genetic_evidence.csv": (
        "disease_id,disease_name,genetic_association_score,disease_url\n"
        "MONDO_0010196,Werner syndrome,0.886,"
        "https://platform.opentargets.org/disease/MONDO_0010196\n"
    ),
}

KRAS_DEPMAP = {
    "selectivity.csv": (
        "gene,n_lines,mean_effect,median_effect,min_effect,pct05_effect,"
        "frac_dependent,n_dependent,sd_effect,frac_positive,p95_effect,classification\n"
        "KRAS,1208,-0.725,-0.517,-4.461,-2.177,0.524,633,0.648,0.0,-0.157,selective\n"
    ),
    "lineage_enrichment.csv": (
        "lineage,n,mean_effect,mean_effect_other,frac_dependent,p_enriched,p_bh\n"
        "Pancreas,48,-2.141,-0.666,0.958,4.966e-22,1.093e-20\n"
        "Bowel,63,-1.43,-0.686,0.857,4.314e-14,4.746e-13\n"
        "Lung,126,-0.882,-0.707,0.579,0.00190,0.01390\n"
    ),
    "mutation_contrast_overall.csv": (
        "dep_gene,marker_gene,marker_kind,n_mutant,n_wt,mean_effect_mutant,"
        "mean_effect_wt,delta,cohens_d,p_mutant_more_dependent\n"
        "KRAS,KRAS,hotspot,177,1031,-1.915,-0.521,-1.395,-3.323,2.225e-308\n"
    ),
    "fusion_contrast.csv": (
        "gene,n_fusion_models,n_fusion_with_crispr,n_fusion_negative,"
        "n_unprofiled_excluded,canonical_only,mean_effect_fusion_pos,"
        "mean_effect_fusion_neg,frac_dependent_fusion_pos,p_fusion_more_dependent,status\n"
        "KRAS,11,8,1132,68,True,-1.473,-0.728,1.0,0.0008,dependent-in-fusion-positive\n"
    ),
    "prism_compounds.csv": (
        "drug,compound_id,moa,dose_uM,screens,release,n_lines,"
        "median_lfc,min_lfc,frac_killed_lfc_lt_1\n"
        "K-RAS(G12C)-INHIBITOR-12,BRD-K00001,K-RAS INHIBITOR,2.5,REP.1M,24Q2,"
        "892,-0.65,-5.185,0.341\n"
        "LONAFARNIB,BRD-K00002,FARNESYLTRANSFERASE INHIBITOR,2.5,REP.PRIMARY,24Q2,"
        "562,-0.169,-2.777,0.146\n"
    ),
}

# Minimal OT data for KRAS: approved drug exists (sotorasib) and lineage-matching diseases.
KRAS_OT = {
    "identity.csv": (
        "ensembl_id,approved_symbol,approved_name,biotype,open_targets_url\n"
        "ENSG00000133703,KRAS,KRAS proto-oncogene GTPase,protein_coding,"
        "https://platform.opentargets.org/target/ENSG00000133703\n"
    ),
    "associations.csv": (
        "disease_id,disease_name,overall_score,disease_url\n"
        "EFO_0002618,pancreatic carcinoma,0.750,https://platform.opentargets.org/disease/EFO_0002618\n"
        "EFO_0005575,colorectal carcinoma,0.680,https://platform.opentargets.org/disease/EFO_0005575\n"
        "EFO_0001071,lung carcinoma,0.620,https://platform.opentargets.org/disease/EFO_0001071\n"
    ),
    "tractability.csv": (
        "label,modality,value\n"
        "Approved Drug,SM,True\n"
        "Structure with Ligand,SM,True\n"
    ),
    "known_drugs.csv": (
        "drug_id,drug_name,drug_type,max_clinical_stage,action_type,"
        "mechanism_of_action,disease_id,disease_name,drug_url\n"
        "CHEMBL4523582,SOTORASIB,Small molecule,APPROVED,INHIBITOR,"
        "KRAS G12C inhibitor,EFO_0001071,lung carcinoma,"
        "https://platform.opentargets.org/drug/CHEMBL4523582\n"
    ),
    "safety.csv": "event,event_id,datasource,literature,url\n",
    "genetic_evidence.csv": (
        "disease_id,disease_name,genetic_association_score,disease_url\n"
        "EFO_0002618,pancreatic carcinoma,0.720,"
        "https://platform.opentargets.org/disease/EFO_0002618\n"
    ),
}


@pytest.fixture
def wrn_depmap_files():
    return WRN_DEPMAP


@pytest.fixture
def wrn_ot_files():
    return WRN_OT


def _write_fixture_dir(tmp_path, subpath, files):
    d = tmp_path / subpath
    d.mkdir(parents=True, exist_ok=True)
    for fname, content in files.items():
        (d / fname).write_text(content)
    return d


@pytest.fixture
def wrn_depmap_dir(tmp_path):
    return _write_fixture_dir(tmp_path, "depmap", WRN_DEPMAP)


@pytest.fixture
def wrn_ot_dir(tmp_path):
    return _write_fixture_dir(tmp_path, "ot", WRN_OT)


@pytest.fixture
def kras_depmap_dir(tmp_path):
    return _write_fixture_dir(tmp_path, "depmap", KRAS_DEPMAP)


@pytest.fixture
def kras_ot_dir(tmp_path):
    return _write_fixture_dir(tmp_path, "ot", KRAS_OT)


# ---------------------------------------------------------------------------
# Run-selection fixture tree
# ---------------------------------------------------------------------------

@pytest.fixture
def results_tree(tmp_path):
    """
    results/
      depmap_local/wrn/
        2026-01-01_120000/  <- complete
        2026-01-01_110000/  <- incomplete (missing lineage_enrichment.csv)
      open_targets/wrn/
        2026-07-01_120000/  <- complete, fresh
        2025-01-01_120000/  <- complete, stale (>90 days old)
        2026-07-01_110000/  <- incomplete (missing safety.csv)
    """
    base = tmp_path / "results"

    # DepMap complete
    dm_complete = base / "depmap_local" / "wrn" / "2026-01-01_120000"
    dm_complete.mkdir(parents=True)
    for fname, content in WRN_DEPMAP.items():
        (dm_complete / fname).write_text(content)

    # DepMap incomplete
    dm_incomplete = base / "depmap_local" / "wrn" / "2026-01-01_110000"
    dm_incomplete.mkdir(parents=True)
    (dm_incomplete / "selectivity.csv").write_text(WRN_DEPMAP["selectivity.csv"])

    # OT complete, fresh
    ot_fresh = base / "open_targets" / "wrn" / "2026-07-01_120000"
    ot_fresh.mkdir(parents=True)
    for fname, content in WRN_OT.items():
        (ot_fresh / fname).write_text(content)

    # OT complete, stale
    ot_stale = base / "open_targets" / "wrn" / "2025-01-01_120000"
    ot_stale.mkdir(parents=True)
    for fname, content in WRN_OT.items():
        (ot_stale / fname).write_text(content)

    # OT incomplete
    ot_incomplete = base / "open_targets" / "wrn" / "2026-07-01_110000"
    ot_incomplete.mkdir(parents=True)
    for fname in ("identity.csv", "associations.csv", "tractability.csv"):
        (ot_incomplete / fname).write_text(WRN_OT[fname])

    return base
