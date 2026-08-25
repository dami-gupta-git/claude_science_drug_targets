# biomarker-stratification

Finds biomarkers already used to select or stratify patients for a
disease's clinical trials, from ClinicalTrials.gov eligibility criteria and
Europe PMC literature, and writes them with source citations to
`results/biomarker-stratification/<disease>/`. Scoped by disease, not gene.

| File | Purpose |
|---|---|
| `kernel.py` | `load_sibling_kernel` (a thin wrapper over `lib.common.skills`) imports `literature-evidence` and `ai-origination-audit`'s kernels so their Europe PMC and ClinicalTrials.gov access isn't duplicated. `search_ctgov_by_condition` fills the condition-search gap those skills don't cover. `trial_eligibility_table` shapes fetched trial records; `biomarkers_table` validates hand-extracted biomarker rows. `write_run` writes all three CSVs plus the README. |
| `SKILL.md` | Workflow, principles, and verification steps. |

Module constants: `CTGOV_STUDIES_API`, `TRIAL_ELIGIBILITY_COLUMNS`,
`BIOMARKER_EVIDENCE_TYPES`, `BIOMARKER_COLUMNS`.

## Running it

`python` cell — once per session, `pip install -e .` from the repo root:

```python
from lib.common.paths import results_dir
from kernel import load_sibling_kernel, search_ctgov_by_condition, write_run

lit = load_sibling_kernel("literature-evidence")
hits = lit.search_europepmc(lit.build_query(disease, extra_terms="biomarker"))
publications_df = lit.publications_table(hits)

audit = load_sibling_kernel("ai-origination-audit")
nct_ids = search_ctgov_by_condition(disease)
studies = audit.fetch_ctgov_studies(nct_ids)

# read eligibility criteria + abstracts, then build biomarker_rows by hand:
biomarker_rows = [
    {"biomarker_id": "B1", "disease": disease, "biomarker": "...",
     "evidence_type": "registry_criterion", "source_id": "NCT...",
     "evidence_text": "...", "url": "..."},
]

out = results_dir("biomarker-stratification", disease.lower().replace(" ", "-"))
write_run(publications_df, studies, biomarker_rows, out, disease, summary, steps)
```
