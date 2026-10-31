---
name: depmap-ot-fusion
description: Fuse a DepMap dependency run and an Open Targets run into a per-gene dossier. Outputs verdict, conflicts, disease-lineage alignment, and a drug join table. Reuses cached upstream runs when they are complete and fresh.
---

# DepMap + Open Targets Fusion

## When to use

Use when you have already run (or want to run) both `depmap-local` and `target-evidence` for a gene and want a side-by-side view: does the cancer dependency pattern match the disease landscape? Does clinical drug precedent exist? Where do the two sources disagree?

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

- `verdict.csv` — one row: dependency classification, direction of effect, top lineage (BH p < 0.05), genotype contrast (first significant one, else first tested; see `genotype_signal`), tractability modalities where True, max clinical stage, safety flag. `effect_direction`, `mean_effect` and `frac_positive` sit next to `frac_dependent` so the fraction is never read without the sign that gives it meaning; `dependency_framework_applies` is False for a growth-suppressing gene.
- `conflicts.csv` — one row per disagreement, with the two observations and a plain-language description. See conflict types below.
- `disease_lineage.csv` — OT diseases aligned to DepMap lineages by keyword matching. `lineage_status` distinguishes the three cases behind an empty effect column: `matched` (lineage found and screened), `no_lineage_keyword` (disease name matched no lineage), `lineage_not_in_depmap_run` (lineage identified but absent from this gene's enrichment run).
- `drugs_joined.csv` — OT known drugs joined to PRISM compounds on drug name (case-insensitive). `compound_id` identifies the specific PRISM entry.
- `sources.json` — both source run directories, their timestamps, and whether each was reused.
- `README.md` — plain-language summary with one intro paragraph, then per-section findings.

## Conflict types

- `growth_suppressing_gene_not_a_dependency` — knockout increases proliferation, so the dependency framework does not apply. Takes precedence over `disease_without_lineage_enrichment`, which fires on the same input but attributes the missing dependency to weak evidence; for a gene whose loss is selected for that is the wrong diagnosis. TP53 is the canonical example: dependent in 0.4% of lines, with 41.9% of lines above +0.3.
- `disease_without_lineage_enrichment` — a high-scoring OT disease (score > 0.5) with no DepMap lineage passing BH correction. Suppressed for a growth-suppressing gene.
- `essential_with_sm_tractability` — common or broadly essential alongside small-molecule tractability (narrow therapeutic window). States whether a genotype deepens the dependency.
- `prism_activity_without_clinical_drug` — PRISM compounds with median LFC < -0.5 where OT records no drug at any clinical phase.
- `non_cancer_top_disease_with_cancer_dependency` — a cancer-context dependency, but the top OT disease is not a cancer. For a selective gene a significant lineage or genotype counts; for an essential gene only a genotype (WRN / Werner syndrome and PELO / MSI are the canonical examples).

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
- WRN produces at least one row in `conflicts.csv` with type `non_cancer_top_disease_with_cancer_dependency`.
- PELO produces `essential_with_sm_tractability` naming the MSI contrast, and its README intro leads with MSI, not a lineage.
- An unknown classification label raises; `CLASSIFICATION_ROLES` keys match `CLASSIFICATIONS` in depmap-local.
- TP53 classifies as `growth-suppressing`, raises `growth_suppressing_gene_not_a_dependency`, names no top lineage despite four lineages passing BH correction, and its README does not contain the string "TP53 is non-essential".
- For any gene, check that the README's dependency sentence agrees in sign with `mean_effect` in `verdict.csv`. A positive mean effect described as a low dependency fraction is the failure this guards.
- A gene with no PRISM data and no known drugs has an empty or all-null `drugs_joined.csv`, not an error.
- `sources.json` contains both `depmap` and `open_targets` keys with `run_dir` and `timestamp`.
