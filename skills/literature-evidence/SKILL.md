---
name: literature-evidence
description: Search Europe PMC for a gene's literature and extract disease-association, interaction, mechanism, and biomarker claims with a source quote and PMID for each. Use when asked for literature evidence on a gene, to check what's published about a target beyond Open Targets and DepMap, or to build a citable claims table for a gene.
---

# Literature evidence

## Workflow

- Invoke the `literature-search-terms` skill for the gene before
  searching, every run, and follow it to the end; it writes the terms
  file. Then
  `terms = load_sibling_kernel("literature-search-terms").latest_search_terms(gene)`
  (`lib.common.skills`). Never hand-write the terms dict or reuse another
  session's file without re-invoking the skill.
- `retrieval = search_balanced(gene, terms)` is the default entry point
  for a dossier — network call, run in the `python` cell directly (Europe
  PMC is a plain REST API, no MCP connector needed). It returns a dict,
  not a tuple: `hits = retrieval["hits"]`. Pass `retrieval` to
  `write_run`; the terms are recorded in `retrieval.json`. State them and
  their reasons in the Retrieval finding. Narrow a broad symbol with
  `extra_terms` (e.g. `"drug target"`, `"biomarker"`) rather than
  requesting more records.
- Three passes run: `recent` (newest first), `landmark` (most cited ever,
  title match) and `established` (most cited per year in the last
  `established_years` years, re-ranked from a larger candidate pool). See
  `search_balanced` for why each is needed.
- **Do not use `search_europepmc` alone for a dossier on an established
  gene.** The Europe PMC default order is newest first, not relevance, so
  a single pass returns only the last year or two. On WRN it returned
  nothing before 2025 at either 12 or 30 records, omitting the 2019 papers
  that established the target — and raising the cap does not help, because
  the order, not the count, is what excludes them. A run that reports the
  current front of a field as if it were the field is the failure this
  guards against.
- Citation sorting alone does not fix it either. On a title-or-abstract
  query the founding WRN paper ranks 15th by citations, behind papers that
  merely mention the gene; half the top 30 are off topic. `search_balanced`
  therefore narrows the gene term to the title for the landmark pass only,
  where that paper ranks 2nd, and keeps abstract matching for the recent
  pass so current work that does not name the gene in its title is still
  seen. Verified across WRN, KRAS and TNIK.
- Each row's `retrieval_pass` names every pass that returned it, joined
  with "; ". Report the year span of what was retrieved — if it starts at
  last year, the retrieval is the finding.
- `build_query(gene, extra_terms=None)` → `search_europepmc(query)` remains
  available for a single targeted pass. Pass `sort=` with one of the
  `SORT_*` constants; an unrecognized sort string returns no results
  rather than erroring at the API, so the helper re-runs the query
  unsorted and raises if that has hits.
- Europe PMC indexes bioRxiv/medRxiv preprints (source `PPR`) alongside
  peer-reviewed literature, so `build_query` includes them by default.
  `publications_table`'s `is_preprint` column marks which rows these are.
  Pass `include_preprints=False` to exclude them instead of filtering the
  table after the fact.
- Read the returned abstracts. For each claim worth recording, write a row
  with `claim_type` (one of `CLAIM_TYPES`), `entity` (the disease, gene,
  pathway, or biomarker the claim is about), `direction` where applicable
  (one of `DIRECTIONS`, read as a statement about the gene; check the
  definition, not just the word — a paper where blocking the gene helps is
  `therapeutic target`, not `increases risk`),
  `source_id` (the `<source>/<id>` of the publication it came from, e.g.
  `MED/30971823`), and `evidence_text` set to the exact sentence the claim
  came from — not a paraphrase. Do not record a claim you cannot point to a
  sentence for. Do not hand-copy `pmid`, `doi` or `url`; `claims_table`
  fills those in from the publication.
- `claims_table(claim_rows, hits)` raises if a required field is missing,
  if `claim_type` or `direction` isn't recognized, if `source_id` names a publication that
  was not retrieved, or if `evidence_text` is not a verbatim substring of
  that publication's abstract. Fix the row, don't work around it.
- `write_run(hits, claim_rows, out, subject, summary, steps,
  retrieval=retrieval)` writes `publications.csv`, `claims.csv`,
  `claims_preprint.csv` where preprint claims exist, `retrieval.json`, and
  `README.md` in one call. Invoke `docs-style` before writing
  `summary`/`steps` prose.

## Principles

- A claim is a claim about a *sentence in a specific publication*, not a
  property of the gene. Two publications making the same claim are two
  rows, not one — a reader needs to see how many independent sources agree.
  Check `duplicate_titles` in `retrieval.json` before counting independent
  sources: a preprint and its published version share a title and would
  otherwise count twice.
- Never summarise an abstract into a claim the abstract doesn't state
  outright. If a claim is implied but not stated, do not record it.
- A claim sourced from a preprint (`is_preprint` true) has not been
  peer-reviewed. Say so wherever the claim is reported in findings prose;
  do not present it with the same weight as a peer-reviewed claim. The
  split into `claims_preprint.csv` marks review status only — both files
  pass the same schema and verbatim checks.
- `pub_types` distinguishes a review from a research article. A mechanism
  claim read from a review is the field restating a finding, not the paper
  that established it; attribute it accordingly rather than treating the
  retrieved record as the primary source.
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
- Report the retrieved year span and the sort used for each pass. Both are
  in `retrieval.json`; a dossier that does not state what the retrieval
  could not have reached is overstating its coverage.

## Verify

`claims_table` now enforces what this section used to ask you to check by
hand: every claim's `source_id` appears in `publications.csv`, and every
`evidence_text` is a verbatim substring of that publication's abstract
after HTML tags and typographic punctuation are folded. A run that writes
a claims table has passed both. What still needs a human eye:

- `publications.csv` spans more than the last two years. A run whose
  earliest record is last year has not seen the field's founding work; if
  that is genuinely a young target, say so explicitly rather than leaving
  the span unremarked.
- Skim the landmark titles. If several are not about the gene, or the
  list is thin for a well-known gene, re-invoke `literature-search-terms`
  with that finding and search again before extracting claims.
- All three passes contributed. If no `retrieval_pass` names `landmark`
  or `established`, that query returned nothing — check the symbol
  appears in titles rather than assuming no such papers exist.
- The year coverage has no gap between the landmark span and the recent
  span that the established pass should have filled.
- `duplicate_titles` in `retrieval.json` is empty, or the duplicates are
  accounted for when claims are counted.
- The landmark pass ranks by citation count, which favours older work and
  admits neighbours that merely mention the gene. Screen its records for
  topicality and record what was excluded.
