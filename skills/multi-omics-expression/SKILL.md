---
name: multi-omics-expression
description: Fetch a gene's median expression across normal GTEx tissues and a fixed set of major TCGA cancer types via cBioPortal, and flag tissues/cancer types where it is markedly over- or under-expressed. Use when asked whether a target is biologically active in a given tissue or tumor type, as a companion check to DepMap dependency evidence.
---

# Multi-omics expression

## Workflow

The two halves of this skill reach their data by different routes, because
only one of them has a connector that covers it.

- **Normal tissue, via the `expression` connector.** In a `repl` cell call
  `gtex_expression_summary(gene=...)`, which returns all 54 tissues ranked
  with a `median` per tissue plus the resolved GENCODE id, then write the
  rows to `./handoff/` and build the table with
  `tissue_expression_table(records)` in the `python` cell. The records it
  needs are dicts with `tissueSiteDetailId` and `median`, which is the
  connector's shape already. Use `gtex_resolve_genes` when you need the
  GENCODE id alone.
- Prefer the connector for GTEx: it reaches releases the kernel cannot
  (the kernel pins `gtex_v8`; the connector reports `gtex_v10` available)
  and offers sample-level distributions via `gtex_gene_expression`, which
  takes `gencode_id`, not a symbol. `resolve_gtex_gene` and
  `fetch_gtex_expression` remain as a fallback and were verified
  equivalent on WRN — same 54 tissues, zero difference in median TPM,
  identical classification.
- **Tumour expression, direct HTTP — no connector covers this.** The
  `cancer-models` connector exposes cBioPortal studies, mutations,
  copy-number and clinical attributes, but has no expression tool, so
  `resolve_cbioportal_gene(gene)` → `fetch_cbioportal_expression(entrez_id)`
  → `tumor_expression_table(study_values)` stays a direct POST to the
  cBioPortal molecular-data route from the `python` cell. Do not replace
  it with a mutation or copy-number tool: those measure something else.
- Say which release each half came from in the run's limits. The two
  halves are on different scales and are each compared only against their
  own table's median, never against each other.
- `write_run(tissue_records, study_values, out, subject, summary, steps)`
  writes both CSVs and README.md in one call. Invoke `docs-style` before
  writing `summary`/`steps` prose.

## Principles

- `CANCER_TYPE_STUDIES` is a fixed list of 12 major TCGA PanCancer Atlas
  studies, not a live per-run lookup — broad coverage, not exhaustive. If
  a study you need isn't on the list, say so in `limits` rather than
  querying cBioPortal's full study catalog ad hoc.
- `classify_level` is the single place "high"/"low"/"typical" is decided,
  for both tissue and tumor tables — do not duplicate the fold-change
  logic in an ad hoc line elsewhere.
- `EXPRESSION_FOLD_THRESHOLD = 4.0`: a tissue/cancer type is "high" at
  ≥4× the gene's own median across the table it appears in, "low" at
  ≤¼×. The two tables use separate reference medians (tissue vs. tumor) -
  they are different assays on different scales and are not compared to
  each other directly.
- `fetch_cbioportal_expression` silently omits a study whose RNA-seq
  profile or sample list isn't found under the expected naming, rather
  than failing the whole run. Check the returned dict's keys against
  `CANCER_TYPE_STUDIES` before writing findings, and name any missing
  study in `limits`.
- Single writer: only the `python` cell writes to `results/`.

## Writing findings

Follow `docs-style`. Table the tissue/cancer type, level, and fold value;
do not restate the raw median TPM/expression figures in prose. Name every
tissue and cancer type flagged "high" or "low" — do not summarise as a
count. If a gene has no tissue or cancer type at "high" or "low", say so;
"typical everywhere" is itself a finding, not an empty result.

## Verify

- `tissue_expression.csv` has one row per GTEx tissue returned.
- `tumor_expression.csv` has at most 12 rows (`CANCER_TYPE_STUDIES`), fewer
  if any study was omitted by `fetch_cbioportal_expression` - cross-check
  against `limits` before treating a short table as complete coverage.
- Every non-"unknown" `level` has a `fold_vs_gene_median` consistent with
  `EXPRESSION_FOLD_THRESHOLD`.
