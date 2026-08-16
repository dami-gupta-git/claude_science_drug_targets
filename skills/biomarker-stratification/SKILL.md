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
- `search_ctgov_by_condition(disease)` → `nct_ids` →
  `load_sibling_kernel("ai-origination-audit").fetch_ctgov_studies(nct_ids)`
  → `trial_eligibility_table(studies)`. Reuses `ai-origination-audit`'s
  concurrent bulk fetch; only the condition search is new here, since that
  skill's `search_ctgov` covers intervention/sponsor, not condition.
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
