---
name: open-targets
description: Fetch target identity, top disease associations, genetic evidence, known drugs/clinical status, tractability, safety, and raw evidence ids for a gene from the Open Targets Platform GraphQL API. Use when asked for a target's druggability profile, its top associated diseases, known drugs against it, or an audit-ready set of Open Targets ids for a gene.
---

# Open Targets

## Workflow

- Resolve the gene symbol to an Ensembl id first (a plain symbol is not
  a valid `ensemblId` argument).
- `repl`: run `query.graphql` against `open_targets_graphql`, write the
  raw response to `handoff/open_targets_<gene>.json`.
- `python`: `read_handoff` → `results_dir("open-targets", gene)` →
  `write_run(data, out)` writes six CSVs.
- Read the CSVs back. Write findings as `summary` and `steps` derived
  strictly from what the CSVs contain — do not state values not present
  in the saved data. Invoke `docs-style` for prose tone.
- Call `write_readme(out / "README.md", gene, summary, steps)` to write
  the README. See README.md for the exact call signature.

## Writing findings

- **Associations**: put scores in a table (disease, score), not in prose.
  Name the top diseases; do not summarise as a count. Do not quote scores
  in the summary paragraph — the scale is undefined to the reader.
- **Tractability**: name the specific labels where `value` is True (e.g.
  "Small-molecule structural evidence, PROTAC criteria"). Do not report
  True/False counts or row totals — they carry no meaning to a reader.
- **Drugs**: name drugs and their highest clinical phase; a count alone
  is not a finding.
- **Safety**: if liabilities are recorded, name the events; do not
  summarise as a count.

## Principles

- Single writer: only `python` writes to `results/`, `results_dir`
  called once per run.
- Every table carries ids (Ensembl, ChEMBL, EFO) next to display names.
- Scoping decisions (20-association cap, `NaN` vs 0, genetic-evidence
  scope, `tractability`'s `False` rows) are documented on each table
  function in `kernel.py` — not restated here.

## Verify

- Six CSVs land in the run directory, one per table function.
- `identity.csv` has exactly one row matching the requested gene;
  `associations.csv` has at most 20.
- Some `query.graphql` fields are unconfirmed against the live schema
  (marked inline) — spot-check the first real response before trusting
  those columns.
