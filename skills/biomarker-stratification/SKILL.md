---
name: biomarker-stratification
description: Find biomarkers already used to select patients for a disease's clinical trials, from ClinicalTrials.gov eligibility criteria and Europe PMC literature. Use when asked what biomarkers stratify patients in a given disease, or to build a citable biomarker table for a disease before designing a trial or evaluating a target's biomarker fit.
---

# Biomarker stratification

## Workflow

- This skill is scoped by disease, not gene.
- `load_sibling_kernel("literature-evidence")` → `lit.build_query(disease,
  extra_terms="biomarker")` → `lit.search_europepmc(query)` →
  `lit.publications_table(hits)`. Reuses `literature-evidence`'s Europe PMC
  access rather than duplicating it.
- **Finding the trials: prefer the `clinical-trials` connector.** In a
  `repl` cell call `search_trials(condition=..., page_size=..., count_total=True)`
  and keep the `nct_id` values. It returns a `total` alongside the page,
  so a run can state what fraction of the disease's trials it read — the
  direct path cannot, and reporting a biomarker table from an unstated
  fraction of a cohort is the failure this avoids. It also accepts phase
  and status filters and boolean condition syntax. Write the ids to
  `./handoff/` for the `python` cell.
- Do NOT use `search_by_eligibility` to find trials for this skill: it
  defaults to recruiting trials only, and a biomarker table built from
  recruiting trials alone silently omits every completed one. It also
  returns no criteria text despite the name.
- **Fetching the eligibility text: keep the direct concurrent path.**
  `load_sibling_kernel("ai-origination-audit").fetch_ctgov_studies(nct_ids)`
  → `trial_eligibility_table(studies)`. The connector carries eligibility
  text only on `get_trial_details`, one trial per call and serial:
  measured at 0.24 s per trial against 0.04 s for the 8-worker bulk
  fetch, a 6.6x difference that scales with the cohort — roughly 2 s
  versus 12 s at 50 trials, 18 s versus 118 s at 500. This skill is
  disease-scoped, so cohorts are large and that gap decides it. The two
  routes return byte-identical `eligibility_criteria`; verified on
  NCT05838768, 2010 characters from both.
- `search_ctgov_by_condition` remains as a fallback for when the connector
  is unavailable. It paginates exhaustively and returns ids only.
- Read the eligibility criteria text and the publication abstracts. For
  each biomarker used to select or stratify patients, write a row naming
  the biomarker, `evidence_type` (`registry_criterion` from eligibility
  text, or `literature_claim` from an abstract), the source id (NCT or
  PMID), and `evidence_text` set to the exact source sentence. Do not
  record a biomarker you cannot point to a sentence for.
- `biomarkers_table(biomarker_rows)` raises if a required field is missing
  or `evidence_type` isn't recognized.
- `write_run(publications_df, studies, biomarker_rows, out, subject,
  summary, steps)` writes all three CSVs and README.md in one call. Invoke
  `docs-style` before writing `summary`/`steps` prose.

## Principles

- A biomarker row is a claim about one trial's eligibility text or one
  publication's abstract, not a property of the disease. The same
  biomarker named by five trials is five rows, not one - the count of
  independent sources is itself a finding.
- No keyword scanner infers a biomarker from criteria/abstract text
  automatically; extraction is a reading step against the actual text
  captured in `trial_eligibility.csv`/the literature table, same discipline
  as `literature-evidence`'s claims.
- Single writer: only the `python` cell writes to `results/`.

## Writing findings

Follow `docs-style`. Group biomarkers by how many independent sources
(trials + publications) name them; do not report a biomarker's presence
without naming which sources support it. If no biomarkers are found for a
disease, say so - that is a finding, not a skipped step.

## Verify

- `biomarkers.csv` has no row whose `source_id` is absent from both
  `trial_eligibility.csv` (`nct_id`) and `literature_publications.csv`
  (`pmid`).
- Every biomarker's `evidence_text` is a substring of that source's
  `eligibility_criteria` or `abstract` - spot-check a sample before
  trusting the table.
