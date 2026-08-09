# open-targets

Fetches a target's identity, top 20 disease associations, genetic
evidence, known drugs and clinical status, tractability, and safety
liabilities from the Open Targets Platform GraphQL API for a given
gene, and writes them as CSVs under `results/open-targets/<gene>/`.
The query text lives in `query.graphql`; the response is shaped into
tables by `kernel.py`.

| File | Purpose |
|---|---|
| `kernel.py` | `identity_table`, `associations_table`, `genetic_evidence_table`, `known_drugs_table`, `tractability_table`, `safety_table` build one table each from a raw Open Targets response; `write_run` calls all six and writes their CSVs. |
| `query.graphql` | The GraphQL query run in the `repl` cell. |
| `SKILL.md` | Workflow, principles, and verification steps. |

## Running it

`repl` cell — repl can't import this repo, so the handoff write is inline:

```python
data = host.mcp("clinical-genomics", "open_targets_graphql",
                query=open("query.graphql").read(),
                variables={"ensemblId": ensembl_id})
import json, os
os.makedirs("handoff", exist_ok=True)
json.dump(data, open(f"handoff/open_targets_{gene}.json", "w"))
```

`python` cell — once per session, `pip install -e .` from the repo root:

```python
from lib.common.handoff import read_handoff
from lib.common.paths import results_dir
from kernel import write_run

data = read_handoff(f"open_targets_{gene}")
out = results_dir("open-targets", gene)
write_run(data, out)
```
