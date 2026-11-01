# literature-search-terms

Chooses the extra names and context words that `literature-evidence`
searches with for a gene, and writes them with a reason for each to
`results/literature_search_terms/<gene>/<timestamp>/search_terms.json`.
The model running the skill decides the terms from the gene symbol and
its approved name only.

| File | Purpose |
|---|---|
| `kernel.py` | `approved_name` fetches the gene's full name from Open Targets. `search_terms` validates the chosen terms. `write_search_terms` saves them; `latest_search_terms` reads the newest record for `literature-evidence`. |
| `SKILL.md` | Workflow, principles, and verification steps. |

## Running it

`python` cell — once per session, `pip install -e .` from the repo root:

```python
from lib.common.paths import results_dir
from kernel import approved_name, search_terms, write_search_terms

name = approved_name(gene)
# decide aliases and context from gene and name (see SKILL.md), e.g. for PELO:
terms = search_terms(gene, name, aliases=["pelota"],
                     context=["ribosome", "Hbs1", "Dom34"],
                     ambiguous=["PELO", "pelota"],
                     reasons={"pelota": "older name, still used in titles",
                              "ribosome": "PELO is also a Portuguese word", ...})
write_search_terms(terms, results_dir("literature-search-terms", gene.lower()))
```

## Sources

| Source | Endpoint | Auth |
|---|---|---|
| Open Targets | `api.platform.opentargets.org` GraphQL | none |
