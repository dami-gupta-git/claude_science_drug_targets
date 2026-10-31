---
name: knowledge-graph
description: Build and query a persistent gene/disease/drug/lineage graph from target-evidence, depmap-local, and literature-evidence run output, and answer path/neighbor questions like how a gene connects to a disease. Use when asked to connect evidence across those skills for a gene, or to trace a path between a gene and a disease/drug.
---

# Knowledge graph

## Workflow

- No external API is called. `find_latest_run(skill_domain, subject_snake)`
  locates the newest run under `results/<skill_domain>/<subject>/` for
  `target_evidence`, `depmap_local`, and `literature_evidence` — run those
  skills first for a gene before syncing it into the graph.
- `ingest_open_targets` / `ingest_depmap_local` / `ingest_literature_evidence`
  each read one skill's CSVs into `(nodes, edges)` lists.
  `nodes_table`/`edges_table` dedupe them into DataFrames.
- `merge_into_store(nodes_df, edges_df, source_skill, source_subject)`
  writes into the persistent store at `graph_store_dir()`
  (`results/knowledge_graph/nodes.csv` and `edges.csv`) — **not** a
  `results_dir()` run. See its docstring for the replace-by-source update
  rule: re-syncing a gene from the same skill replaces that skill's edges
  for that gene, it does not accumulate duplicates.
- `write_sync_report(subject, summary, steps)` writes this call's own
  audit-trail README into a normal timestamped
  `results/knowledge-graph/<subject>/<ts>/` run — separate from the store
  itself. Invoke `docs-style` before writing `summary`/`steps` prose.
- `query_neighbors(node_id, edges_df)` and `query_path(from_id, to_id,
  edges_df)` answer questions against `load_store()`'s current `edges_df`.
  Build a `node_id` with `node_id(node_type, label)` — types are `gene`,
  `disease`, `drug`, `lineage` (`NODE_TYPES`).

## Principles

- This skill reads other skills' `results/` CSVs; it never re-derives their
  numbers. If a needed run doesn't exist, say so and point at which skill
  to run, rather than approximating from stale or partial data.
- `ingest_depmap_local` applies the same dependency-direction filter
  (`mean_effect < 0` and `frac_dependent > 0`) that
  `docs/depmap_ot_fusion_lineage_direction.md` fixed for depmap-ot-fusion —
  do not remove it; a positive mean effect is not a dependency edge in any
  direction.
- An edge's `evidence_ref` and `source_run` are what make it checkable —
  never write an edge without both.
- `score` is NOT comparable across predicates and must never be used to rank
  or compare edges of different kinds. It holds an Open Targets association
  score for `associated_with`, a dependent fraction for
  `dependency_in_lineage`, and nothing at all for `targeted_by_drug`. Three
  incomparable quantities in one column.
- Statistical support lives in `support_p`, with `support_kind` naming what it
  is. Pass both to `make_edge` whenever the source run computed a test; omit
  both when it did not, and the edge records as untested rather than as
  supported. Never recover a p-value by parsing `evidence_ref` — that field
  also carries URLs and PMIDs, so it is not numeric.
- `query_path` filters out edges failing their own test at `alpha` before
  searching, and this changes answers rather than trimming them. Most
  dependency edges in a typical store do not survive multiple-testing
  correction, so unfiltered BFS returns routes through edges at p near 1.0,
  and hop count makes such a route look preferable to a longer well-supported
  one. `alpha=1.0` gives the raw topology; say so if you report it.
- A store written before `support_p` existed will make the query functions
  raise. That is intended: re-ingest the source runs rather than filtering by
  hand, so the support travels with the edge.
- A literature-evidence claim's entity is typed `disease` only when
  `claim_type` is `disease_association`; every other claim type's entity is
  recorded as a generic `gene` node, since `claims.csv` doesn't distinguish
  further. Don't guess a more specific type from the entity string.
- Single writer to the persistent store: only `merge_into_store` writes
  `nodes.csv`/`edges.csv`; nothing else in this skill or another touches
  them directly.

## Writing findings

Follow `docs-style`. For a path query, name every hop and the edge's
`source_skill`/`evidence_ref` — a path is only as trustworthy as its
weakest edge. For a neighbor query, group by predicate, not by a raw edge
count.

## Verify

- Every `edge_id` in `edges.csv` is unique and matches
  `subject_id|predicate|object_id|source_run`.
- Every `node_id` referenced by an edge (`subject_id` or `object_id`)
  exists in `nodes.csv`.
- After `merge_into_store` for a given `(source_skill, source_subject)`,
  no edge with that pair predates the new `source_run`.
