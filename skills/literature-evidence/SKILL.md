---
name: literature-evidence
description: Search Europe PMC for a gene's literature and extract disease-association, interaction, mechanism, and biomarker claims with a source quote and PMID for each. Use when asked for literature evidence on a gene, to check what's published about a target beyond Open Targets and DepMap, or to build a citable claims table for a gene.
---

# Literature evidence

## Workflow

- `build_query(gene, extra_terms=None)` → `search_europepmc(query)` — network
  call, run in the `python` cell directly (Europe PMC is a plain REST API,
  no MCP connector needed). Narrow a broad symbol with `extra_terms`
  (e.g. `"drug target"`, `"biomarker"`) rather than raising `max_pages`.
- Europe PMC indexes bioRxiv/medRxiv preprints (source `PPR`) alongside
  peer-reviewed literature, so `build_query` includes them by default.
  `publications_table`'s `is_preprint` column marks which rows these are.
  Pass `include_preprints=False` to exclude them instead of filtering the
  table after the fact.
- Read the returned abstracts. For each claim worth recording, write a row
  with `claim_type` (one of `CLAIM_TYPES`), `entity` (the disease, gene,
  pathway, or biomarker the claim is about), `direction` where applicable,
  and `evidence_text` set to the exact sentence the claim came from — not
  a paraphrase. Do not record a claim you cannot point to a sentence for.
- `claims_table(claim_rows)` raises if a required field is missing or
  `claim_type` isn't recognized — fix the row, don't work around it.
- `write_run(hits, claim_rows, out, subject, summary, steps)` writes
  `publications.csv`, `claims.csv`, and `README.md` in one call. Invoke
  `docs-style` before writing `summary`/`steps` prose.

## Principles

- A claim is a claim about a *sentence in a specific publication*, not a
  property of the gene. Two publications making the same claim are two
  rows, not one — a reader needs to see how many independent sources agree.
- Never summarise an abstract into a claim the abstract doesn't state
  outright. If a claim is implied but not stated, do not record it.
- A claim sourced from a preprint (`is_preprint` true) has not been
  peer-reviewed. Say so wherever the claim is reported in findings prose;
  do not present it with the same weight as a peer-reviewed claim.
- Single writer: only the `python` cell writes to `results/`.
- `search_europepmc` returns publications; it does not extract claims.
  Claim extraction is a reading step, done against the actual abstract
  text in `publications_table`, not inferred from title alone.

## Writing findings

Follow `docs-style`. Literature-specific guidance:

- Group claims by `claim_type` in the findings, not by publication.
- Name the entity and cite the PMID next to it; do not summarise as a count
  of claims found.
- If the search returned publications but none yielded a recordable claim,
  say so — that is a finding, not a skipped step.

## Verify

- `claims.csv` has no row whose `pmid` is absent from `publications.csv`.
- Every claim's `evidence_text` is a substring of that publication's
  `abstract` in `publications.csv` — spot-check a sample before trusting
  the table.
