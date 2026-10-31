---
name: literature-search-terms
description: Choose the extra names and context words literature-evidence searches with for a gene, from the gene symbol and its approved name only, and save them with a reason for each. Use before literature-evidence for any gene, and especially when papers use another name for it (TP53 as p53) or the symbol is also an ordinary word (PELO).
---

# Literature search terms

## Workflow

- `name = approved_name(gene)`. Use only the symbol and this name; do not
  search the literature or look anything else up.
- Decide two lists:
  - `aliases`: other names papers put in titles for this gene (TP53 →
    `["p53"]`, MAP3K13 → `["LZK"]`). Empty if the symbol is what the field
    uses. Never a name shared with another gene.
  - `context`: only if the symbol is also an ordinary word or another
    gene's name, up to `MAX_CONTEXT_TERMS` words a real paper's abstract
    would contain (PELO → `["ribosome", "mRNA", "Hbs1", "Dom34"]`).
    Otherwise empty; it narrows every search pass.
- Give each term a one-line reason in `reasons`.
- `terms = search_terms(gene, name, aliases, context, reasons)` validates;
  `write_search_terms(terms, results_dir("literature-search-terms", gene.lower()))`
  saves `search_terms.json`.

## Principles

- Empty lists are a valid answer and the usual one. Add a term only when
  the symbol alone would miss or mis-hit the field.
- A wrong alias pulls in another gene's papers; a wrong context word
  silently drops real ones. When unsure, leave it out.

## Verify

- `search_terms.json` exists for the gene and every term has a reason.
- After literature-evidence runs, skim its landmark titles. Off-topic
  titles mean a collision `context` should fix; a thin landmark list for a
  well-known gene means a missing alias. Rerun this skill, then
  literature-evidence.
