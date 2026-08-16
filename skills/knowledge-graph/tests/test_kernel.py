import csv

import pytest

import kernel


def write_csv(path, rows, headers):
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=headers)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def test_node_id_is_deterministic_and_normalized():
    assert kernel.node_id("gene", " WRN ") == "GENE:wrn"
    assert kernel.node_id("gene", "wrn") == kernel.node_id("gene", "WRN")


def test_make_node_rejects_bad_type():
    with pytest.raises(ValueError, match="not in"):
        kernel.make_node("protein", "WRN")


def test_find_latest_run_picks_newest_timestamp(tmp_path):
    subject_dir = tmp_path / "open_targets" / "wrn"
    (subject_dir / "2024-01-01_000000").mkdir(parents=True)
    (subject_dir / "2024-06-01_000000").mkdir(parents=True)
    latest = kernel.find_latest_run("open_targets", "wrn", base=tmp_path)
    assert latest.name == "2024-06-01_000000"


def test_find_latest_run_missing_subject_returns_none(tmp_path):
    assert kernel.find_latest_run("open_targets", "nope", base=tmp_path) is None


def test_ingest_open_targets(tmp_path):
    run_dir = tmp_path / "2024-01-01_000000"
    run_dir.mkdir()
    write_csv(run_dir / "associations.csv",
              [{"disease_id": "EFO:1", "disease_name": "MSI cancer",
                "overall_score": "0.8", "disease_url": "http://x/EFO:1"}],
              ["disease_id", "disease_name", "overall_score", "disease_url"])
    write_csv(run_dir / "known_drugs.csv",
              [{"drug_id": "CHEMBL1", "drug_name": "Drugolol",
                "drug_url": "http://x/CHEMBL1"}],
              ["drug_id", "drug_name", "drug_url"])

    nodes, edges = kernel.ingest_open_targets("WRN", run_dir)
    node_ids = {n["node_id"] for n in nodes}
    assert "GENE:wrn" in node_ids
    assert "DISEASE:msi cancer" in node_ids
    assert "DRUG:drugolol" in node_ids
    predicates = {e["predicate"] for e in edges}
    assert predicates == {"associated_with", "targeted_by_drug"}
    assoc = next(e for e in edges if e["predicate"] == "associated_with")
    assert assoc["score"] == 0.8
    assert assoc["source_run"] == "2024-01-01_000000"


def test_ingest_depmap_local_filters_by_direction(tmp_path):
    run_dir = tmp_path / "2024-01-01_000000"
    run_dir.mkdir()
    write_csv(run_dir / "lineage_enrichment.csv", [
        {"lineage": "Bowel", "mean_effect": "-0.6", "frac_dependent": "0.4", "p_bh": "0.01"},
        {"lineage": "Lymphoid", "mean_effect": "0.1", "frac_dependent": "0.0", "p_bh": "0.03"},
    ], ["lineage", "mean_effect", "frac_dependent", "p_bh"])

    nodes, edges = kernel.ingest_depmap_local("RB1", run_dir)
    lineage_labels = {n["label"] for n in nodes if n["node_type"] == "lineage"}
    assert lineage_labels == {"Bowel"}
    assert len(edges) == 1
    assert edges[0]["predicate"] == "dependency_in_lineage"


def test_ingest_literature_evidence(tmp_path):
    run_dir = tmp_path / "2024-01-01_000000"
    run_dir.mkdir()
    write_csv(run_dir / "claims.csv", [
        {"claim_type": "disease_association", "entity": "MSI cancer", "pmid": "111"},
        {"claim_type": "interaction", "entity": "NCK-interacting kinase", "pmid": "222"},
        {"claim_type": "", "entity": "unused", "pmid": "333"},
    ], ["claim_type", "entity", "pmid"])

    nodes, edges = kernel.ingest_literature_evidence("TRAF2", run_dir)
    types_by_label = {n["label"]: n["node_type"] for n in nodes}
    assert types_by_label["MSI cancer"] == "disease"
    assert types_by_label["NCK-interacting kinase"] == "gene"
    assert len(edges) == 2


def test_nodes_table_dedupes_by_id():
    nodes = [kernel.make_node("gene", "WRN"), kernel.make_node("gene", "wrn")]
    df = kernel.nodes_table(nodes)
    assert len(df) == 1
    assert list(df.columns) == kernel.NODE_COLUMNS


def test_edges_table_dedupes_by_id():
    e1 = kernel.make_edge("GENE:wrn", "associated_with", "DISEASE:x",
                          "open-targets", "wrn", "run1")
    e2 = dict(e1)
    df = kernel.edges_table([e1, e2])
    assert len(df) == 1
    assert list(df.columns) == kernel.EDGE_COLUMNS


def test_merge_into_store_replaces_by_source(tmp_path):
    e_old = kernel.make_edge("GENE:wrn", "associated_with", "DISEASE:old",
                             "open-targets", "wrn", "run1")
    n_old = [kernel.make_node("gene", "WRN"), kernel.make_node("disease", "old")]
    kernel.merge_into_store(kernel.nodes_table(n_old), kernel.edges_table([e_old]),
                            "open-targets", "wrn", store_dir=tmp_path)

    e_new = kernel.make_edge("GENE:wrn", "associated_with", "DISEASE:new",
                             "open-targets", "wrn", "run2")
    n_new = [kernel.make_node("gene", "WRN"), kernel.make_node("disease", "new")]
    result = kernel.merge_into_store(kernel.nodes_table(n_new), kernel.edges_table([e_new]),
                                     "open-targets", "wrn", store_dir=tmp_path)

    _, edges_df = kernel.load_store(store_dir=tmp_path)
    assert result["edges_added"] == 1
    assert "DISEASE:old" not in edges_df["object_id"].values
    assert "DISEASE:new" in edges_df["object_id"].values


def test_merge_into_store_keeps_other_sources(tmp_path):
    e_ot = kernel.make_edge("GENE:wrn", "associated_with", "DISEASE:x",
                            "open-targets", "wrn", "run1")
    kernel.merge_into_store(kernel.nodes_table([kernel.make_node("gene", "WRN")]),
                            kernel.edges_table([e_ot]), "open-targets", "wrn",
                            store_dir=tmp_path)
    e_dm = kernel.make_edge("GENE:wrn", "dependency_in_lineage", "LINEAGE:bowel",
                            "depmap-local", "wrn", "run1")
    kernel.merge_into_store(kernel.nodes_table([kernel.make_node("gene", "WRN")]),
                            kernel.edges_table([e_dm]), "depmap-local", "wrn",
                            store_dir=tmp_path)
    _, edges_df = kernel.load_store(store_dir=tmp_path)
    assert len(edges_df) == 2


def test_query_neighbors():
    edges = kernel.edges_table([
        kernel.make_edge("GENE:wrn", "associated_with", "DISEASE:x", "open-targets", "wrn", "r1"),
        kernel.make_edge("GENE:other", "associated_with", "DISEASE:y", "open-targets", "other", "r1"),
    ])
    neighbors = kernel.query_neighbors("GENE:wrn", edges)
    assert len(neighbors) == 1


def test_query_path_direct_and_multihop():
    edges = kernel.edges_table([
        kernel.make_edge("GENE:a", "interacts_with", "GENE:b", "literature-evidence", "a", "r1"),
        kernel.make_edge("GENE:b", "associated_with", "DISEASE:c", "open-targets", "b", "r1"),
    ])
    assert kernel.query_path("GENE:a", "DISEASE:c", edges) == ["GENE:a", "GENE:b", "DISEASE:c"]
    assert kernel.query_path("GENE:a", "GENE:a", edges) == ["GENE:a"]


def test_query_path_unreachable_returns_none():
    edges = kernel.edges_table([
        kernel.make_edge("GENE:a", "interacts_with", "GENE:b", "literature-evidence", "a", "r1"),
    ])
    assert kernel.query_path("GENE:a", "DISEASE:nowhere", edges) is None


def test_write_sync_report(tmp_path, monkeypatch):
    monkeypatch.setenv("DRUG_TARGETS_ROOT", str(tmp_path))
    path = kernel.write_sync_report(
        "WRN", "Synced open-targets and depmap-local for WRN.",
        [{"name": "Sync", "finding": "2 edges added, 0 removed."}])
    assert path.is_file()
    assert "# WRN — Knowledge graph sync" in path.read_text()
