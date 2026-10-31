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
- `depmap_selectivity(gene)` — common-essential / selective / non-essential / growth-suppressing.
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
- `0.50 <= frac_dependent < 0.90` → broadly essential: most lines dependent, a minority spared. Also a narrow window; check whether a genotype deepens the effect.
- Positive effect means knockout helps growth — typical of tumour suppressors. Never read it as a weak dependency.

### The dependency question is one-sided

`frac_dependent` counts only the negative tail, so on a gene whose loss is
selected for it returns a true but misleading null. TP53 is dependent in 0.4% of
lines; reported on its own that reads as "not a target", when the finding is that
knockout confers a proliferation advantage (mean effect +0.42, with 41.9% of
lines above +0.3). The sign of the effect is the result, not a detail.

`depmap_selectivity` therefore returns `growth-suppressing` when the positive
tail clears `POS_THRESHOLD` (10% of lines above `POS_EFFECT` = +0.3) and the mean
effect is positive. Both values were calibrated against a 579-gene random
sample, where the median gene has 0.75% of lines above +0.3 and only about 2%
of genes reach 10%.

Two consequences when reporting such a gene:

- Lead with the positive tail and the direction. Quote `frac_dependent` only to
  say it is the expected result for this class of gene, never as evidence about
  druggability.
- Do not name a top dependent lineage. For a growth-suppressing gene every
  lineage mean effect is typically positive, and a Mann-Whitney test will still
  return significant lineages — those are lineages where loss helps growth
  *less*, not dependencies. TP53 has four lineages passing BH correction and no
  dependency among them.

The classification is decided on the data, not on whether the gene is annotated
as a tumour suppressor. APC, STK11 and VHL are tumour suppressors that are real
dependencies in this screen (VHL is below `DEP_THRESHOLD` in 77% of lines) and
correctly classify as dependencies (`selective` or `broadly essential`). Conversely a gene can show a positive tail
with few screened lines — check `n_lines` before drawing conclusions.

## Pitfalls

- `CRISPRInferredCommonEssentials.csv` is not a frequency threshold — it includes KRAS (dependent in ~half of lines). Use `frac_dependent >= 0.90`; call `depmap_common_essentials()` for cross-reference only.
- PRISM and CRISPR are different releases with ~three-fifths overlap; coverage is uneven by lineage. Always report the joint n.
- The secondary screen (`prism-repurposing-20q2-secondary-*`) is a separate, older release covering only primary-screen hits. AUC polarity is inverted vs primary LFC (higher AUC = more resistant); prefer AUC when the compound is present.
- Never aggregate PRISM by `Drug.Name`; one name can span multiple compound IDs with different behaviour. The compound ID is the unit.
- `depmap_msi_dose_response()` returns three arms and the pooled one is
  deliberately labelled `not interpretable`. Roughly 94% of lines sit below the
  MSI threshold, so a pooled rank correlation mixes the between-group step with
  any within-group gradient and reports neither: on WRN it comes out at -0.086,
  significant only because n=1208, explaining 0.7% of variance against 51.8%
  inside the MSI-high arm. Read the `mss_only` and `msi_high_only` rows, and use
  `r_squared` rather than `p_value` to compare them — the arms differ ~15x in n.
- Read the MSS arm as the control it is. A gradient inside the MSI-high arm means
  little unless the stable arm is flat, and the two together distinguish claims
  that are easy to conflate: flat-then-graded is a threshold for entry followed
  by a dose-response, which is neither "graded" nor "switch" on its own.
- `depmap_msi_dose_bins()` bins within each arm, not across the panel. Binning
  across it cannot answer the question: whole-panel quantiles put every MSI-high
  line in the top bin, so the output is flat-then-step whatever the true shape,
  including when a clean gradient exists above the threshold. If you compare
  against an older run, that is why the bin table changed.
- The conventional MSI cut-off is not necessarily where dependency begins. On WRN
  the lowest MSI-high quartile (scores 24-55) is indistinguishable from stable
  lines at p = 0.69, while above ~55 the contrast is p = 3e-17. Check the lowest
  quartile against the stable arm before treating threshold membership as the
  selecting variable.
