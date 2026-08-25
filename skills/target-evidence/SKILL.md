---
name: target-evidence
description: Build a target dossier from four sources - Open Targets (identity, disease associations, genetic evidence, drugs, tractability, safety), gnomAD (loss-of-function constraint and carrier counts), cBioPortal (the gene's mutations in MSI-high tumours, with variant allele frequencies), and DepMap via depmap-local (whether dependency scales with MSI burden). Use when asked for a target's druggability profile, whether heterozygous loss of the gene is tolerated in the population, whether resistance residues are already present in treatment-naive tumours, or whether an MSI dependency is graded or a threshold effect.
---

# Target evidence

## Workflow

- Every section is optional except Open Targets. Run only the sections the
  question needs, and record the ones skipped in `steps` so a reader can
  tell "not checked" from "checked and empty".
- Open Targets: resolve the gene symbol to an Ensembl id first (a plain
  symbol is not a valid `ensemblId`). `repl`: run `query.graphql` against
  `open_targets_graphql`, write the raw response to
  `handoff/open_targets_<gene>.json`. `python`: `read_handoff` →
  `results_dir("target-evidence", gene)`.
- gnomAD: `fetch_gnomad_gene(symbol)` → pass the record to `write_run`.
  Read `gnomad_constraint_table` for LOEUF and observed-vs-expected LoF,
  and `gnomad_carrier_summary_table` for how many carriers exist.
- cBioPortal: `fetch_cbioportal_msi_cohort(symbol, residue_ranges,
  thresholds=MSI_ATTRIBUTE_THRESHOLDS)`. `residue_ranges` is required and
  caller-supplied - look the domain boundaries up (UniProt or the paper
  being tested) and cite the source in the findings. Sweeping every
  MSI-scored study takes minutes; pass `study_ids` to restrict it, and say
  in the findings which studies were covered.
- DepMap MSI dose-response: `msi_dose_response(gene)`, which delegates to
  `depmap-local`. Needs `$DEPMAP_ROOT`.
- Invoke `docs-style` before writing findings.
- `write_run(data, out, subject, summary, steps, gnomad=..., cbioportal=...,
  msi_dose=...)` writes the six Open Targets CSVs plus a CSV per supplied
  section, then README.md.

## Writing findings

Follow the prose-vs-tables rule in `docs-style`. Source-specific guidance:

- **Associations**: put scores in a table (disease, score), not in prose.
  Name the top diseases; do not summarise as a count. The finding must
  define the score scale: `overall_score` is 0-1, aggregated across
  evidence datatypes by Open Targets; higher means stronger evidence.
  Do not quote scores in the summary paragraph - define the scale once
  in the finding text, then let the table carry the values.
- **Tractability**: name the specific labels where `value` is True (e.g.
  "Small-molecule structural evidence, PROTAC criteria"). Do not report
  True/False counts or row totals - they carry no meaning to a reader.
- **Drugs**: name drugs and their highest clinical phase; a count alone
  is not a finding.
- **Safety**: if liabilities are recorded, name the events; do not
  summarise as a count. An empty Open Targets safety table is not
  evidence of safety; say what gnomAD constraint shows instead.
- **Constraint**: LOEUF is the upper bound of the observed/expected LoF
  ratio; lower means the population tolerates loss of the gene less.
  State the scale, then give observed and expected LoF counts. For a
  recessive gene, say explicitly that carrier counts describe
  heterozygotes and that `total_homozygotes` is a floor (see the
  docstring), so no claim of a two-hit frequency rests on it.
- **Tumour mutations**: a mutation count is only readable against the
  MSI-high sample count in `tumour_cohort.csv` - give both. Report
  variant allele frequency when arguing a variant is clonal or
  subclonal, and name the sample. Zero mutations in a residue range,
  over a stated number of MSI-high samples, is a finding.
- **Dose-response**: report the correlation within the MSI-high arm as
  well as over all lines; the all-lines correlation is partly the
  MSI-high/MSS split restated. A flat within-arm correlation is a
  finding, not a failed test - say which of the two the data supports.

## Principles

- Single writer: only `python` writes to `results/`, `results_dir`
  called once per run.
- Every table carries ids (Ensembl, ChEMBL, EFO, gnomAD variant,
  cBioPortal study/sample) next to display names.
- Nothing scientific is defaulted: MSI call thresholds and residue
  ranges are passed in by the caller with a cited source, and a missing
  VAF or constraint record is left empty or raises rather than filled.
- Scoping decisions (20-association cap, `NaN` vs 0, LOFTEE
  high-confidence-only filter, exome/genome pooling, MSI-high-only
  mutation rows) are documented on each function in `kernel.py` - not
  restated here.

## Verify

- Six Open Targets CSVs land in the run directory, plus two or three
  per supplied section.
- `identity.csv` has exactly one row matching the requested gene;
  `associations.csv` has at most 20.
- Every `sample_id` in `tumour_mutations.csv` belongs to a study listed
  in `tumour_cohort.csv`, and `n_msi_high_mutated` never exceeds
  `n_msi_high_samples`.
- `gnomad_lof_variants.csv` contains only `lof == HC` rows.
- A study scored by two MSI attributes contributes one mutation row per
  attribute; deduplicate on study, sample and protein change before
  counting mutations across attributes.
- Some `query.graphql` fields are unconfirmed against the live schema
  (marked inline) - spot-check the first real response before trusting
  those columns.
