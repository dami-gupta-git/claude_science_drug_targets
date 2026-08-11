---
name: depmap-local
description: Query a locally downloaded DepMap release - CRISPR gene effect/dependency (Chronos), cell-line metadata, hotspot mutations, and PRISM Repurposing compound sensitivity. Use when asked whether cancer cells depend on a gene, which lineages or genotypes are selectively dependent, whether a gene is common-essential or selective, for synthetic-lethal / biomarker contrasts (e.g. WRN in MSI, KRAS in KRAS-mutant lines), or to find PRISM compounds annotated against a target. Handles the wide 18k-column matrices without loading them whole.
---

# DepMap (local release)

## Requirements

- The `claude_science_drug_targets` repo must be installed in the active environment: `pip install -e /path/to/claude_science_drug_targets`. The kernel imports `lib.common.reporting`; without this install the helpers will not load.
- Set `$DEPMAP_ROOT` to the release directory, or pass `root=` to every helper.

## Workflow

- `depmap_inventory()` — confirm files present and whether the cache is built. Run this first; call `build_depmap_cache()` if the cache is absent (one-off, ~15 s, needs `pyarrow`).
- `depmap_selectivity(gene)` — common-essential / selective / non-essential.
- `depmap_lineage_enrichment(gene)` — which lineages are selectively dependent (BH-adjusted).
- `depmap_mutation_contrast(dep_gene, marker_gene)` — genotype-driven dependency.
- `depmap_fusion_contrast(gene, canonical_only=True)` — dependency in fusion-positive lines; distinguishes *untestable* (too few lines) from a real negative.
- `depmap_prism_compounds(target)` — compounds annotated to the target. Defaults to newest release on disk; `depmap_prism_releases()` lists them.
- `write_run(tables, out_dir, subject, summary, steps)` — write CSVs and README.md in one call. Defaults to `source="DepMap"`.
- Invoke `docs-style` before writing findings prose.

## Writing findings

Follow the prose-vs-tables rule in `docs-style`. DepMap-specific guidance:

- **Selectivity**: state the classification and what it implies ("selective — most lines unaffected, a subset dependent"). 
- **Lineage enrichment**: name which lineages show enrichment and whether any survive multiple-testing correction.
- **Contrast sections** (mutation, MSI, fusion): state direction and biological meaning. Explain acronyms on first use (e.g., "MSI — microsatellite instability").
- **PRISM**: name compounds and whether they show activity. If none, note the target may lack annotated compounds in this library.
- **Fusion**: state whether the result is testable; if not, give the reason in plain terms.

## Interpreting gene effect

- `< -1.0` strong dependency; `< -0.5` conventional dependent cut (`DEP_THRESHOLD`).
- `frac_dependent >= 0.90` → common essential: usually a poor therapeutic window.
- Positive effect means knockout helps growth — typical of tumour suppressors. Never read it as a weak dependency.

## Pitfalls

- `CRISPRInferredCommonEssentials.csv` is not a frequency threshold — it includes KRAS (dependent in ~half of lines). Use `frac_dependent >= 0.90`; call `depmap_common_essentials()` for cross-reference only.
- PRISM and CRISPR are different releases with ~three-fifths overlap; coverage is uneven by lineage. Always report the joint n.
- The secondary screen (`prism-repurposing-20q2-secondary-*`) is a separate, older release covering only primary-screen hits. AUC polarity is inverted vs primary LFC (higher AUC = more resistant); prefer AUC when the compound is present.
- Never aggregate PRISM by `Drug.Name`; one name can span multiple compound IDs with different behaviour. The compound ID is the unit.
