# literature-evidence

Searches Europe PMC for a gene's literature — including bioRxiv/medRxiv
preprints, flagged via `is_preprint` — and writes the matching publications
plus a table of disease-association, interaction, mechanism, and biomarker
claims extracted from their abstracts, each carrying the source sentence
and PMID, to `results/literature-evidence/<gene>/`.

| File | Purpose |
|---|---|
| `kernel.py` | `build_query` builds the Europe PMC query, including bioRxiv/medRxiv preprints (source `PPR`) by default — pass `include_preprints=False` to exclude them. `search_europepmc` fetches and pages through results. `publications_table` shapes the raw hits, flagging preprints via `is_preprint`. `claims_table` validates hand-extracted claim rows. `write_run` calls both and writes their CSVs plus the README. |
| `SKILL.md` | Workflow, principles, and verification steps. |

Module constants: `EUROPEPMC_SEARCH_URL`, `EUROPEPMC_ARTICLE_URL`,
`PUBLICATION_COLUMNS`, `CLAIM_TYPES`, `CLAIM_COLUMNS`, `PREPRINT_SOURCE`.

## Running it

`python` cell — once per session, `pip install -e .` from the repo root:

```python
from lib.common.paths import results_dir
from lib.common.skills import load_sibling_kernel
from kernel import search_balanced, write_run

# run the literature-search-terms skill for the gene first
terms = load_sibling_kernel("literature-search-terms").latest_search_terms(gene)
retrieval = search_balanced(gene, terms)
hits = retrieval["hits"]

# read hits' abstracts, then build claim_rows by hand:
claim_rows = [
    {"claim_id": "C1", "gene": gene, "claim_type": "disease_association",
     "entity": "...", "direction": "...", "evidence_text": "...",
     "pmid": "...", "url": "..."},
]

out = results_dir("literature-evidence", gene.lower())
write_run(hits, claim_rows, out, gene, summary, steps, retrieval=retrieval)
```
