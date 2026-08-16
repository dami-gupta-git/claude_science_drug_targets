# multi-omics-expression

Fetches a gene's median expression across GTEx normal tissues and a fixed
set of 12 major TCGA cancer types (via cBioPortal), and flags tissues or
cancer types where the gene is markedly over- or under-expressed relative
to its own median, to `results/multi-omics-expression/<gene>/`.

| File | Purpose |
|---|---|
| `kernel.py` | `resolve_gtex_gene`/`fetch_gtex_expression` and `resolve_cbioportal_gene`/`fetch_cbioportal_expression` fetch raw records from each API. `tissue_expression_table`/`tumor_expression_table` shape them, using the shared `classify_level` to flag high/low expression. `write_run` calls both and writes their CSVs plus the README. |
| `SKILL.md` | Workflow, principles, and verification steps. |

Module constants: `GTEX_GENE_URL`, `GTEX_EXPRESSION_URL`, `GTEX_DATASET_ID`,
`CBIOPORTAL_BASE`, `CANCER_TYPE_STUDIES`, `EXPRESSION_FOLD_THRESHOLD`,
`TISSUE_EXPRESSION_COLUMNS`, `TUMOR_EXPRESSION_COLUMNS`.

## Running it

`python` cell — once per session, `pip install -e .` from the repo root:

```python
from lib.common.paths import results_dir
from kernel import (resolve_gtex_gene, fetch_gtex_expression,
                    resolve_cbioportal_gene, fetch_cbioportal_expression,
                    write_run)

tissue_records = fetch_gtex_expression(resolve_gtex_gene(gene))
study_values = fetch_cbioportal_expression(resolve_cbioportal_gene(gene))

out = results_dir("multi-omics-expression", gene.lower())
write_run(tissue_records, study_values, out, gene, summary, steps)
```
