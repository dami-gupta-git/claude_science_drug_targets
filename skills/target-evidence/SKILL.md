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
  `write_run` also needs `missense_min_allele_frequency`; it has no
  default, so state the cutoff chosen and why in the findings.
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
  msi_dose=...)` writes the seven Open Targets CSVs plus a CSV per supplied
  section, then README.md.

## Writing findings

Follow the prose-vs-tables rule in `docs-style`.

- Each section states its finding in one or two sentences, plus the one
  definition the reader needs to read its table (a score scale, a
  cutoff). Method caveats go in Limits, once each, never in a section
  body: the association cap, the genetic-score slice, salt folding, the
  approved-tier merge, why the two LoF counts differ, ancestry pooling
  and the genome build.

Source-specific guidance:

- **Associations**: put scores in a table (disease, score), not in prose.
  Name the top diseases; do not summarise as a count. The finding must
  define the score scale: `overall_score` is 0-1, aggregated across
  evidence datatypes by Open Targets; higher means stronger evidence.
  Do not quote scores in the summary paragraph - define the scale once
  in the finding text, then let the table carry the values.
- **Genetic scores**: `top_association_genetic_scores.csv` covers only
  the top associations (see `top_association_genetic_scores_table`). Do
  not present it as where the gene's genetic support lies.
- **Tractability**: name the specific labels where `value` is True (e.g.
  "Small-molecule structural evidence, PROTAC criteria"). Do not report
  True/False counts or row totals - they carry no meaning to a reader.
- **Drugs**: name drugs and their highest clinical phase; a count alone
  is not a finding. Build the README drug table from `drug_summary.csv`
  in its row order (see `drug_summary_table`), never by hand from
  `known_drugs.csv`. If trimmed, cut from the bottom and state how many
  drugs were left out. Write "approved per Open Targets", not "approved":
  its labels lag regulators and miss some regional approvals.
- **Safety**: if liabilities are recorded, name the events; do not
  summarise as a count. An empty Open Targets safety table is not
  evidence of safety; say what gnomAD constraint shows instead.
- **Constraint**: LOEUF is the upper bound of the observed/expected LoF
  ratio; lower means the population tolerates loss of the gene less.
  State the scale, then give observed and expected LoF counts. For a
  recessive gene, say explicitly that carrier counts describe
  heterozygotes and that `total_homozygotes` is a floor (see the
  docstring), so no claim of a two-hit frequency rests on it. For an X
  or Y gene, also give `total_hemizygotes`: those XY people have no
  working copy from a single LoF allele.
- **Missense variants**: name each common variant by protein change
  (hgvsp) with its carrier share and homozygote count, and state the cutoff.
  gnomAD shows only that a variant is common, not that it changes
  function or disease risk; say so, and cite a source for any functional
  claim. An empty table means no missense variant cleared the cutoff.
  Where `reference_is_minor` is true, the listed allele is the usual form
  in the population; describe the reference allele as the rare one.
  When saying how many people carry a variant, quote
  `minor_allele_carrier_fraction`, never the allele frequency, which
  counts gene copies. It is empty on X, Y and mitochondrial genes; say so
  instead of estimating it (see `minor_allele_carriers`).
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

- Seven Open Targets CSVs land in the run directory, plus two to four
  per supplied section.
- Every APPROVED drug in `drug_summary.csv` appears in the README drug
  table.
- `identity.csv` has exactly one row matching the requested gene;
  `associations.csv` has at most 20.
- Every `sample_id` in `tumour_mutations.csv` belongs to a study listed
  in `tumour_cohort.csv`, and `n_msi_high_mutated` never exceeds
  `n_msi_high_samples`.
- `associations.csv` is capped at the query's page size, currently 20, ordered
  by descending score — the strongest associations, not all of them. WRN has
  1,002. Every run writes `association_coverage.csv` with the true total and a
  `truncated` flag; quote it whenever you characterise what a gene is
  associated with, or the table reads as a survey when it is a top slice. The
  cap loses the tail rather than a random sample, so "the top associations are
  X" is supportable and "the gene is associated with X diseases" is not.
- Both gnomAD variant tables hold only canonical-transcript rows with a
  nonzero `allele_count` from QC-passing callsets (see
  `gnomad_carried_variant_rows`); `transcript_id` matches the gene's
  canonical transcript on every row.
- `gnomad_lof_variants.csv` contains only `lof == HC` rows.
- `gnomad_missense_variants.csv` contains only missense rows whose
  `minor_allele_frequency` is at or above the stated cutoff.
- The two gnomAD LoF counts differ and both are correct: `obs_lof` in
  `gnomad_constraint.csv` is the constraint model's count, which also
  excludes variants LOFTEE flags, while `gnomad_lof_variants.csv` keeps
  flagged HC calls and marks them in `flags`. Say which you mean when
  quoting one; adjacent in a sentence they read as a contradiction. Both tables carry
  `reference_genome`, and variant ids are coordinates on that build — state it,
  because the same id means a different position on another assembly.
- When the question is "is residue X mutated in these tumours", pass
  `test_positions={label: position}` to `fetch_cbioportal_msi_cohort`. It raises
  if any tested position falls outside every `residue_range`. Without it the
  wrong answer is silent rather than absent: an uncovered position is labelled
  outside the named ranges, dropped from `n_msi_high_in_range`, and the run reads
  as evidence of absence. Annotated domains are the trap — they do not tile the
  sequence, so a residue in an unannotated stretch between two domains is
  invisible. Either tile the full sequence with explicit unannotated blocks, or
  query resistance windows and domains as two separate passes and write both
  tables.
- `residue_range_hit` returns the FIRST matching range, so overlapping ranges
  hide all but one label. Two passes, not one merged dict, when a residue belongs
  to both a resistance window and a domain.
- An absence result over a bulk-sequenced cohort is bounded by detection, not by
  sample size. Call `clonality_floor(mutations_table)` and quote it: the claim
  supported is that no allele is clonal or near-clonal, since a subclonal variant
  falls below the calling floor and pre-existing resistance is characteristically
  subclonal. Sample size is the looser bound and the tempting one to cite.
- The cohort table's `n_msi_high_in_range` counts SAMPLES; the mutation table
  counts MUTATIONS, and a sample can carry several. Label whichever appears in
  prose, and never compare one against the other.
- A study scored by two MSI attributes contributes one mutation row per
  attribute; deduplicate on study, sample and protein change before
  counting mutations across attributes.
- Some `query.graphql` fields are unconfirmed against the live schema
  (marked inline) - spot-check the first real response before trusting
  those columns.
