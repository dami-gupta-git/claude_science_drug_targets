---
name: depmap-ot-fusion
description: Fuse a DepMap dependency run and an Open Targets run into a per-gene dossier. Outputs verdict, conflicts, disease-lineage alignment, and a drug join table. Reuses cached upstream runs when they are complete and fresh.
---

# DepMap + Open Targets Fusion

## When to use

Use when you have already run (or want to run) both `depmap-local` and `open-targets` for a gene and want a side-by-side view: does the cancer dependency pattern match the disease landscape? Does clinical drug precedent exist? Where do the two sources disagree?

## Workflow

- `python`: `run_depmap_ot_fusion(gene)` — the single entry point. Ensures both upstream runs exist, fuses them, and writes output to `results/depmap_ot_fusion/<gene>/<timestamp>/`.
- Set `reuse=False` to regenerate both upstream runs regardless of what is cached.
- Set `dry_run=True` to print planned actions without writing anything.
- Invoke `docs-style` before writing findings prose.

Step-by-step when you need finer control:

1. `ensure_depmap_run(gene, reuse, root)` → depmap run directory
2. `ensure_open_targets_run(gene, reuse, max_age_days)` → OT run directory
3. `fuse(depmap_dir, ot_dir)` → dict of DataFrames (verdict, conflicts, disease_lineage, drugs_joined)
4. `write_fusion_run(tables, sources, out_dir)` → CSVs + `sources.json`
5. `write_fusion_readme(out_dir, gene, tables, sources)` → `README.md`

## Outputs

- `verdict.csv` — one row: dependency classification, top lineage (BH p < 0.05), genotype marker if tested, tractability modalities where True, max clinical stage, safety flag.
- `conflicts.csv` — one row per disagreement, with the two observations and a plain-language description. See conflict types below.
- `disease_lineage.csv` — OT diseases aligned to DepMap lineages by keyword matching. `lineage_status` distinguishes the three cases behind an empty effect column: `matched` (lineage found and screened), `no_lineage_keyword` (disease name matched no lineage), `lineage_not_in_depmap_run` (lineage identified but absent from this gene's enrichment run).
- `drugs_joined.csv` — OT known drugs joined to PRISM compounds on drug name (case-insensitive). `compound_id` identifies the specific PRISM entry.
- `sources.json` — both source run directories, their timestamps, and whether each was reused.
- `README.md` — plain-language summary with one intro paragraph, then per-section findings.

## Conflict types

- `disease_without_lineage_enrichment` — a high-scoring OT disease (score > 0.5) with no DepMap lineage passing BH correction.
- `common_essential_with_sm_tractability` — common-essential classification alongside small-molecule tractability (indicates a narrow therapeutic window).
- `prism_activity_without_clinical_drug` — PRISM compounds with median LFC < -0.5 where OT records no drug at any clinical phase.
- `non_cancer_top_disease_with_selective_dependency` — selective dependency enriched in a cancer lineage, but the top OT disease is not a cancer (WRN / Werner syndrome is the canonical example).

## Reuse and staleness

- DepMap: treated as fresh unless `max_age_days` is set. The DepMap release on disk is static; passing `max_age_days` enforces a time limit.
- Open Targets: 90-day default. Scores move between OT platform releases, so older runs may not reflect current data.
- Completeness is checked before reuse: a run directory missing required CSVs, or whose `identity.csv` does not match the requested gene, is skipped.

## Decisions and constraints

- Disease-to-lineage mapping uses keyword matching on disease names. The keyword list (`LINEAGE_KEYWORDS` in `kernel.py`) covers the DepMap OncotreeLineage vocabulary. Keyword specificity determines priority (longer keyword wins); "uveal melanoma" maps to Eye, not Skin.
- Drug join is on normalised drug name. One drug name can span multiple compound IDs with different behaviour; `compound_id` in `drugs_joined.csv` identifies the specific PRISM record.
- `ensure_depmap_run` runs selectivity, lineage enrichment, MSI contrast, and fusion contrast. Mutation contrast requires a marker gene and is not run automatically; run `depmap-local` directly to add it.
- No composite score is produced. `verdict.csv` copies or labels fields from source columns without weighting.

## Verify

- `verdict.csv` has exactly one row.
- WRN produces at least one row in `conflicts.csv` with type `non_cancer_top_disease_with_selective_dependency`.
- A gene with no PRISM data and no known drugs has an empty or all-null `drugs_joined.csv`, not an error.
- `sources.json` contains both `depmap` and `open_targets` keys with `run_dir` and `timestamp`.
