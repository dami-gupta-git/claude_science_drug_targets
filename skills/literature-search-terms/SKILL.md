---
name: literature-search-terms
description: Choose the extra names and context words literature-evidence searches with for a gene, from the gene symbol and its approved name only, and save them with a reason for each. Use before literature-evidence for any gene, and especially when papers use another name for it (TP53 as p53) or the symbol is also an ordinary word (PELO).
---

# Literature search terms

## Workflow

- `name = approved_name(gene)`. Use only the symbol and this name; do not
  search the literature or look anything else up.
- Decide three lists:
  - `aliases`: other names papers put in titles for this gene (TP53 →
    `["p53"]`, MAP3K13 → `["LZK"]`). Empty if the symbol is what the field
    uses. Never a name shared with another gene.
  - `context`: only if the symbol or an alias is also an ordinary word,
    another gene's name or some other label (MAP3K13's alias LZK also names
    bacterial strains), up to `MAX_CONTEXT_TERMS` words a real paper's
    title or abstract would contain (PELO → `["ribosome", "Hbs1", "Dom34"]`).
    Otherwise empty.
    Each word must be specific to this gene's biology. A word that any
    molecular-biology abstract might contain does not disambiguate:
    "mRNA" let a goat-pregnancy proteomics paper and a rat-artery paper
    through as PELO hits, because both are Portuguese-language papers
    using "pelo" whose abstracts happen to mention mRNA. Prefer a binding
    partner, a complex, the homologue's name, or a process named after
    the gene; reject "mRNA", "protein", "gene", "cell", "expression".
  - `ambiguous`: the symbol and/or aliases that collide (PELO →
    `["PELO", "pelota"]`, MAP3K13 → `["LZK"]`). Only papers matched on
    these names must carry a context word; the others match freely, so
    real papers that never use the context words are kept. Set together
    with `context` or not at all.
- Give each term a one-line reason in `reasons`.
- `terms = search_terms(gene, name, aliases, context, ambiguous, reasons)` validates;
  `write_search_terms(terms, results_dir("literature-search-terms", gene.lower()))`
  saves `search_terms.json`.

## Principles

- Empty lists are a valid answer and the usual one. Add a term only when
  the symbol alone would miss or mis-hit the field.
- A wrong alias pulls in another gene's papers. A context word that is
  too narrow silently drops real ones; one that is too generic lets
  unrelated papers back in. When unsure, leave it out.

## Verify

- `search_terms.json` exists for the gene and every term has a reason.
- After literature-evidence runs, skim the titles from every pass, not
  just the landmark one. Off-topic titles mean a collision `context`
  should fix, or a context word too generic to exclude them; a thin
  landmark list for a well-known gene means a missing alias. Rerun this
  skill, then literature-evidence.
