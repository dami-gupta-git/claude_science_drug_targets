# depmap-ot-fusion: lineage direction filter

## What changed

`build_verdict`, `find_conflicts`, and `write_fusion_readme` in
`skills/depmap-ot-fusion/kernel.py` selected a gene's "top lineage" from
`lineage_enrichment.csv` using only `p_bh` (Benjamini-Hochberg-adjusted
significance). They did not check whether the lineage's effect was actually
in the dependency direction.

All three now first filter to rows where `mean_effect < 0` (knockout reduces
viability) and `frac_dependent > 0` (at least one line in the lineage crosses
the dependency threshold), and only rank by `p_bh` within that filtered set.
If no lineage passes the filter, no top lineage is reported.

## Why this was wrong

`depmap_lineage_enrichment` (in `skills/depmap-local/kernel.py`) runs a
Mann-Whitney test comparing each lineage's knockout effect against all other
lineages. This test only asks "is this lineage's effect different from the
rest," not "is this lineage dependent on the gene." For a tumor-suppressor
gene, knocking the gene out helps cell growth almost everywhere, so every
lineage's `mean_effect` is positive. The test can still find one lineage
whose positive effect is *smaller* than the others' and call that
"significant" - but a smaller positive effect is not a dependency in any
direction; it is a lineage that benefits less than most from losing the
gene.

## How it surfaced

Running the fusion on RB1 (a canonical tumor suppressor) produced a verdict
naming "Lymphoid" as the lineage with "the strongest selective dependency"
(BH p = 0.030), worded identically to how the tool reports a real dependency
finding (e.g. EGFR in Head and Neck lines). The Lymphoid row's own
`frac_dependent` was 0.000 - no line in it was actually dependent - and its
`mean_effect` was +0.126, positive like every other lineage in RB1's table.
The gene's own selectivity call was correctly "non-essential" (0.17% of
lines dependent); the lineage claim contradicted it.

## Scope

Any tumor-suppressor-like gene (uniformly positive knockout effect across
lineages) is subject to the same false positive, since a smaller-magnitude
positive effect can still reach significance by chance in a 20+ lineage
comparison. Genes with genuine negative-effect dependencies (ALK, BRCA1,
KRAS, WRN, EGFR - all run earlier in this project) were not affected, because
their top-ranked lineage by `p_bh` already had `mean_effect < 0`.

## Files affected

- `skills/depmap-ot-fusion/kernel.py`: `build_verdict`, `find_conflicts`,
  `write_fusion_readme`.
- No change to `skills/depmap-local/kernel.py` - `depmap_lineage_enrichment`
  itself is unchanged and correctly documents the interpretation trap
  ("Positive effect means knockout helps growth... never read it as a weak
  dependency") in its own `SKILL.md`. The bug was in `depmap-ot-fusion`
  ignoring that guidance when it consumed the table.
