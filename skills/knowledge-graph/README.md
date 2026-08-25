# knowledge-graph

Builds and queries a persistent gene/disease/drug/lineage graph assembled
from `target-evidence`, `depmap-local`, and `literature-evidence` run output.
No external API is called — this skill reads the CSVs those skills already
wrote to `results/`. The graph store lives at `results/knowledge_graph/`
(`nodes.csv`, `edges.csv`) and accumulates across syncs; each sync also
writes its own audit-trail README to
`results/knowledge-graph/<subject>/<timestamp>/`.

| File | Purpose |
|---|---|
| `kernel.py` | `find_latest_run` locates a sibling skill's newest run for a subject. `ingest_open_targets`/`ingest_depmap_local`/`ingest_literature_evidence` read that run's CSVs into nodes/edges. `nodes_table`/`edges_table` dedupe them. `merge_into_store` writes them into the persistent store, replacing that source's prior edges for the subject. `query_neighbors`/`query_path` answer questions against the store. `write_sync_report` writes the per-sync README. |
| `SKILL.md` | Workflow, principles, and verification steps. |

Module constants: `NODE_COLUMNS`, `EDGE_COLUMNS`, `NODE_TYPES`,
`CLAIM_PREDICATES`.

## Running it

`python` cell — once per session, `pip install -e .` from the repo root.
Run `target-evidence`, `depmap-local`, and/or `literature-evidence` for the
gene first; this skill only reads their output.

```python
from kernel import (find_latest_run, ingest_open_targets, ingest_depmap_local,
                    ingest_literature_evidence, nodes_table, edges_table,
                    merge_into_store, load_store, query_path, write_sync_report)

nodes, edges = [], []
for skill_domain, ingest in [
    ("target_evidence", ingest_open_targets),
    ("depmap_local", ingest_depmap_local),
    ("literature_evidence", ingest_literature_evidence),
]:
    run_dir = find_latest_run(skill_domain, gene.lower())
    if run_dir:
        n, e = ingest(gene, run_dir)
        nodes += n
        edges += e
        merge_into_store(nodes_table(n), edges_table(e), skill_domain, gene.lower())

write_sync_report(gene, summary, steps)

# querying:
_, edges_df = load_store()
query_path("GENE:%s" % gene.lower(), "DISEASE:some disease", edges_df)
```
