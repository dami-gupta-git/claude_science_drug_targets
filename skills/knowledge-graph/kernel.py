"""Build and query a persistent gene/disease/drug/lineage graph from other
skills' results/ output.

Unlike the other skills, this one does not call an external API - it reads
the CSVs already written by open-targets, depmap-local, and
literature-evidence runs. It is also the only skill with a store that
accumulates across calls rather than one dossier per run: nodes.csv and
edges.csv live in graph_store_dir(), not inside a timestamped
results_dir() run. See merge_into_store for the replace-by-source update
rule. kernel.py loads with this skill: read_csv_rows, find_latest_run,
ingest_open_targets, ingest_depmap_local, ingest_literature_evidence,
nodes_table, edges_table, graph_store_dir, load_store, merge_into_store,
query_neighbors, query_path, write_sync_report, and the constants
NODE_COLUMNS, EDGE_COLUMNS, NODE_TYPES, CLAIM_PREDICATES.
"""
import csv
from pathlib import Path

import pandas as pd


NODE_COLUMNS = ["node_id", "node_type", "label"]

EDGE_COLUMNS = [
    "edge_id", "subject_id", "predicate", "object_id",
    "source_skill", "source_subject", "source_run", "evidence_ref", "score",
    "support_p", "support_kind",
]

# `score` holds a different, incomparable quantity per predicate: an Open
# Targets association score for associated_with, a dependent fraction for
# dependency_in_lineage, and nothing at all for targeted_by_drug. It is
# therefore not orderable across predicates and must never be used to rank or
# compare edges of different kinds.
#
# `support_p` is the statistical support where one exists - a BH-adjusted
# p-value for dependency edges - promoted out of the free-text `evidence_ref`
# field, which also carries URLs and PMIDs and so cannot be parsed as a
# number. `support_kind` names what the p-value is, or the reason there is
# none. Query functions filter on support_p; nothing should reach into
# evidence_ref to recover it.
SUPPORT_KIND_BH = "bh_adjusted_p"
SUPPORT_KIND_NONE = "no_statistical_test"
SUPPORT_ALPHA = 0.05

# depmap-local surfaces lineage (tissue/cancer-type) dependency, not pathway
# membership - "pathway" from the original four-skill plan is realized as
# "lineage" here since no pathway source is wired up yet.
NODE_TYPES = ("gene", "disease", "drug", "lineage")

# literature-evidence's claim_type -> this graph's edge predicate.
CLAIM_PREDICATES = {
    "disease_association": "associated_with",
    "interaction": "interacts_with",
    "mechanism": "mechanism_of",
    "biomarker": "biomarker_for",
}


def normalize_label(label: str) -> str:
    return str(label).strip().lower()


def node_id(node_type: str, label: str) -> str:
    return "%s:%s" % (node_type.upper(), normalize_label(label))


def make_node(node_type: str, label: str) -> dict:
    if node_type not in NODE_TYPES:
        raise ValueError("node_type %r not in %s" % (node_type, NODE_TYPES))
    return {"node_id": node_id(node_type, label), "node_type": node_type, "label": label}


def make_edge(subject_id: str, predicate: str, object_id: str, source_skill: str,
             source_subject: str, source_run: str, evidence_ref: str = None,
             score: float = None, support_p=None,
             support_kind: str = None) -> dict:
    """One edge row. `support_p` is the statistical support, not `score`.

    Pass `support_p` with a `support_kind` naming what it is whenever the
    source run computed a test. Leave both unset when it did not, and the
    edge is recorded as untested rather than as supported - an edge with no
    p-value must not read as a strong one.
    """
    edge_id = "|".join([subject_id, predicate, object_id, source_run])
    if support_p is not None and not support_kind:
        raise ValueError("support_p given without support_kind: a bare number "
                         "cannot be interpreted downstream")
    try:
        support_p = None if support_p is None else float(support_p)
    except (TypeError, ValueError):
        support_p = None
    return {
        "edge_id": edge_id, "subject_id": subject_id, "predicate": predicate,
        "object_id": object_id, "source_skill": source_skill,
        "source_subject": source_subject, "source_run": source_run,
        "evidence_ref": evidence_ref, "score": score,
        "support_p": support_p,
        "support_kind": support_kind or SUPPORT_KIND_NONE,
    }


def read_csv_rows(path: Path) -> list:
    """Rows of a CSV file as dicts, or [] if the file doesn't exist."""
    path = Path(path)
    if not path.is_file():
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def find_latest_run(skill_domain: str, subject_snake: str, base: Path = None) -> Path:
    """Newest results/<skill_domain>/<subject_snake>/<timestamp>/ dir, or None."""
    from lib.common.paths import project_root
    base = Path(base) if base else project_root() / "results"
    subject_dir = base / skill_domain / subject_snake
    if not subject_dir.is_dir():
        return None
    candidates = sorted(
        (d for d in subject_dir.iterdir() if d.is_dir()),
        key=lambda d: d.name, reverse=True)
    return candidates[0] if candidates else None


def ingest_open_targets(gene: str, run_dir: Path) -> tuple:
    """Gene-disease and gene-drug edges from an open-targets run directory."""
    run_dir = Path(run_dir)
    nodes = [make_node("gene", gene)]
    edges = []
    source_run = run_dir.name
    gene_id = node_id("gene", gene)

    for row in read_csv_rows(run_dir / "associations.csv"):
        if not row.get("disease_name"):
            continue
        disease_node = make_node("disease", row["disease_name"])
        nodes.append(disease_node)
        edges.append(make_edge(
            gene_id, "associated_with", disease_node["node_id"],
            "open-targets", gene.lower(), source_run,
            evidence_ref=row.get("disease_url"),
            score=float(row["overall_score"]) if row.get("overall_score") else None))

    for row in read_csv_rows(run_dir / "known_drugs.csv"):
        if not row.get("drug_name"):
            continue
        drug_node = make_node("drug", row["drug_name"])
        nodes.append(drug_node)
        edges.append(make_edge(
            gene_id, "targeted_by_drug", drug_node["node_id"],
            "open-targets", gene.lower(), source_run,
            evidence_ref=row.get("drug_url")))

    return nodes, edges


def ingest_depmap_local(gene: str, run_dir: Path) -> tuple:
    """Gene-lineage dependency edges from a depmap-local run directory.

    Only lineages with mean_effect < 0 (knockout reduces viability) and
    frac_dependent > 0 are edges - the same direction filter
    docs/depmap_ot_fusion_lineage_direction.md fixed for depmap-ot-fusion,
    applied here for the same reason: a smaller positive mean_effect is not
    a dependency in any direction.
    """
    run_dir = Path(run_dir)
    nodes = [make_node("gene", gene)]
    edges = []
    source_run = run_dir.name
    gene_id = node_id("gene", gene)

    for row in read_csv_rows(run_dir / "lineage_enrichment.csv"):
        if not row.get("lineage"):
            continue
        mean_effect = float(row["mean_effect"]) if row.get("mean_effect") else None
        frac_dependent = float(row["frac_dependent"]) if row.get("frac_dependent") else None
        if mean_effect is None or frac_dependent is None:
            continue
        if not (mean_effect < 0 and frac_dependent > 0):
            continue
        lineage_node = make_node("lineage", row["lineage"])
        nodes.append(lineage_node)
        edges.append(make_edge(
            gene_id, "dependency_in_lineage", lineage_node["node_id"],
            "depmap-local", gene.lower(), source_run,
            evidence_ref=row.get("p_bh"), score=frac_dependent,
            support_p=row.get("p_bh"), support_kind=SUPPORT_KIND_BH))

    return nodes, edges


def ingest_literature_evidence(gene: str, run_dir: Path) -> tuple:
    """Gene-entity claim edges from a literature-evidence run directory."""
    run_dir = Path(run_dir)
    nodes = [make_node("gene", gene)]
    edges = []
    source_run = run_dir.name
    gene_id = node_id("gene", gene)

    for row in read_csv_rows(run_dir / "claims.csv"):
        claim_type = row.get("claim_type")
        predicate = CLAIM_PREDICATES.get(claim_type)
        if not predicate or not row.get("entity"):
            continue
        # A literature claim's entity is typed only when the claim_type
        # itself implies a type (disease_association -> disease);
        # interaction/mechanism/biomarker entities could be a gene, protein,
        # or pathway that claims.csv doesn't distinguish, so they're
        # recorded as "gene" nodes rather than guessed further.
        entity_type = "disease" if claim_type == "disease_association" else "gene"
        entity_node = make_node(entity_type, row["entity"])
        nodes.append(entity_node)
        edges.append(make_edge(
            gene_id, predicate, entity_node["node_id"],
            "literature-evidence", gene.lower(), source_run,
            evidence_ref=row.get("pmid")))

    return nodes, edges


def nodes_table(nodes: list) -> pd.DataFrame:
    """One row per distinct node_id (last label wins on duplicates)."""
    by_id = {}
    for n in nodes:
        by_id[n["node_id"]] = n
    return pd.DataFrame(list(by_id.values()), columns=NODE_COLUMNS)


def edges_table(edges: list) -> pd.DataFrame:
    """One row per distinct edge_id."""
    by_id = {}
    for e in edges:
        by_id[e["edge_id"]] = e
    return pd.DataFrame(list(by_id.values()), columns=EDGE_COLUMNS)


def graph_store_dir(root: str = None) -> Path:
    """results/knowledge_graph/, created if missing.

    Not a results_dir() run: this store accumulates across calls rather
    than getting a fresh timestamped directory each time - see
    merge_into_store for how updates are reconciled.
    """
    from lib.common.paths import project_root
    base = Path(root) if root else project_root() / "results"
    store = base / "knowledge_graph"
    store.mkdir(parents=True, exist_ok=True)
    return store


def load_store(store_dir: Path = None) -> tuple:
    """(nodes_df, edges_df) from the persistent store, empty if not yet built."""
    store_dir = Path(store_dir) if store_dir else graph_store_dir()
    nodes_path = store_dir / "nodes.csv"
    edges_path = store_dir / "edges.csv"
    nodes_df = pd.read_csv(nodes_path) if nodes_path.is_file() else pd.DataFrame(columns=NODE_COLUMNS)
    edges_df = pd.read_csv(edges_path) if edges_path.is_file() else pd.DataFrame(columns=EDGE_COLUMNS)
    return nodes_df, edges_df


def merge_into_store(nodes_df: pd.DataFrame, edges_df: pd.DataFrame, source_skill: str,
                     source_subject: str, store_dir: Path = None) -> dict:
    """Merge a fresh ingest into the persistent store.

    Update rule: every existing edge with this exact (source_skill,
    source_subject) is dropped and replaced by the new edges - a newer run
    for the same gene/skill supersedes the old one rather than
    accumulating alongside it, since re-ingesting reflects a re-run, not a
    second independent source. Nodes are unioned by node_id (new label
    wins on conflict); nodes are never deleted, since another still-current
    edge may reference one this ingest didn't touch.
    """
    store_dir = Path(store_dir) if store_dir else graph_store_dir()
    existing_nodes, existing_edges = load_store(store_dir)

    if len(existing_edges):
        keep = ~((existing_edges["source_skill"] == source_skill) &
                (existing_edges["source_subject"] == source_subject))
        existing_edges = existing_edges[keep]
    merged_edges = pd.concat([existing_edges, edges_df], ignore_index=True)
    merged_edges = merged_edges.drop_duplicates(subset="edge_id", keep="last")

    merged_nodes = pd.concat([existing_nodes, nodes_df], ignore_index=True)
    merged_nodes = merged_nodes.drop_duplicates(subset="node_id", keep="last")

    merged_nodes.to_csv(store_dir / "nodes.csv", index=False)
    merged_edges.to_csv(store_dir / "edges.csv", index=False)
    return {"nodes": len(merged_nodes), "edges": len(merged_edges),
           "edges_added": len(edges_df)}


def query_neighbors(node_id_: str, edges_df: pd.DataFrame,
                    alpha: float = None) -> pd.DataFrame:
    """Edges touching node_id_ on either side.

    Unfiltered by default, because listing a node's neighbourhood is a
    different question from asserting a connection - but read `support_p`
    before quoting any row: most dependency edges in a typical store do not
    survive multiple-testing correction. Pass `alpha` to drop those, and note
    that `score` is not comparable across predicates, so the returned rows
    cannot be ranked by it.
    """
    hits = edges_df[(edges_df["subject_id"] == node_id_) |
                    (edges_df["object_id"] == node_id_)]
    return hits if alpha is None else supported_edges(hits, alpha=alpha)


def supported_edges(edges_df: pd.DataFrame, alpha: float = 0.05,
                    keep_untested: bool = True) -> pd.DataFrame:
    """Edges whose statistical support clears `alpha`, plus untested ones.

    An edge with `support_p` above alpha failed the test its own source run
    performed, and a route through it is noise presented as a connection. An
    edge with no `support_p` was never tested - Open Targets association and
    drug edges carry no p-value - and is kept by default because dropping it
    would remove most of the graph; pass keep_untested=False to see only
    edges with positive statistical support.
    """
    if "support_p" not in edges_df.columns:
        raise ValueError(
            "edges_df has no support_p column: it predates support tracking. "
            "Re-ingest the source runs, or filter on the p-value in "
            "evidence_ref yourself and say so in the report")
    support = pd.to_numeric(edges_df["support_p"], errors="coerce")
    tested = support.notna()
    return edges_df[(tested & (support <= alpha)) |
                    (~tested if keep_untested else False)]


def query_path(from_id: str, to_id: str, edges_df: pd.DataFrame,
               max_hops: int = 4, alpha: float = 0.05) -> list:
    """Shortest undirected node_id path from from_id to to_id, or None.

    Breadth-first search over edges_df treated as an undirected graph -
    the predicate direction matters for reading an edge, not for whether
    two nodes are connected.

    Edges failing their own significance test at `alpha` are excluded before
    the search. Without that filter the shortest path is the wrong answer
    rather than a weak one: on a store where most dependency edges do not
    survive multiple-testing correction, BFS returns routes through edges at
    p near 1.0, and hop count makes such a route look preferable to a longer
    well-supported one. Pass alpha=1.0 for the unfiltered topology, knowing
    the result is connectivity and not evidence.
    """
    edges_df = supported_edges(edges_df, alpha=alpha)
    adjacency = {}
    for _, row in edges_df.iterrows():
        adjacency.setdefault(row["subject_id"], set()).add(row["object_id"])
        adjacency.setdefault(row["object_id"], set()).add(row["subject_id"])

    if from_id == to_id:
        return [from_id]
    frontier = [[from_id]]
    visited = {from_id}
    for _ in range(max_hops):
        next_frontier = []
        for path in frontier:
            for neighbor in adjacency.get(path[-1], ()):
                if neighbor in visited:
                    continue
                new_path = path + [neighbor]
                if neighbor == to_id:
                    return new_path
                visited.add(neighbor)
                next_frontier.append(new_path)
        frontier = next_frontier
        if not frontier:
            break
    return None


def write_sync_report(subject: str, summary: str, steps: list,
                      source: str = "Knowledge graph sync", limits: list = None) -> Path:
    """Write this ingest's own README into results/knowledge_graph/<subject>/<ts>/.

    The persistent store (nodes.csv/edges.csv in graph_store_dir()) is not
    duplicated here; this is only the audit trail for one sync call.
    """
    from lib.common.paths import results_dir
    from lib.common.reporting import write_readme
    out_dir = results_dir("knowledge-graph", subject.lower())
    return write_readme(out_dir / "README.md", subject, summary, steps,
                        source=source, limits=limits)
