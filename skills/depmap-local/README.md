# depmap-local

Queries a locally downloaded DepMap release for gene-level dependency
analysis: CRISPR gene effect (Chronos), lineage enrichment, genotype
and fusion contrasts, and PRISM compound sensitivity. No network access;
all data files are downloaded manually from depmap.org.

| File | Purpose |
|---|---|
| `kernel.py` | All helpers: `depmap_selectivity`, `depmap_lineage_enrichment`, `depmap_mutation_contrast`, `depmap_fusion_contrast`, `depmap_msi_contrast`, `depmap_prism_compounds`, and `write_run`. |
| `SKILL.md` | Workflow, interpretation rules, and pitfalls. |

## Setup

Set `$DEPMAP_ROOT` to the directory holding the DepMap release files,
then install the repo so `lib.common.reporting` is importable:

```bash
pip install -e /path/to/claude_science_drug_targets
export DEPMAP_ROOT=/path/to/depmap/release
```

## Running it

All analysis runs in a `python` cell. A typical single-gene run:

```python
from lib.common.paths import results_dir
from lib.common.reporting import write_readme
from kernel import (depmap_selectivity, depmap_lineage_enrichment,
                    depmap_mutation_contrast, write_run)

gene = "KRAS"
sel   = depmap_selectivity(gene)
lin   = depmap_lineage_enrichment(gene)
mut   = depmap_mutation_contrast(gene, "KRAS", kind="hotspot")

out = results_dir("depmap-local", gene.lower())
write_run({"selectivity": sel, "lineage_enrichment": lin,
           "mutation_contrast": mut}, out)
write_readme(out / "README.md", gene, summary, steps, source="DepMap")
```

Call `build_depmap_cache()` once after a new release download to build
the parquet cache (~15 s, needs `pyarrow`).
