# clinical-trials

Fetches a drug's trials from the ClinicalTrials.gov v2 API and writes the
ranked trial list, a phase × status summary, the stopped trials with the
sponsor's stated reason, and a table of population signals (biomarker
requirements, prior-therapy requirements, combination partners) extracted
from eligibility text, each carrying the source passage and NCT ID, to
`results/clinical_trials/<drug>/`.

| File | Purpose |
|---|---|
| `kernel.py` | `search_trials` fetches and pages through studies for an intervention. `trials_table` shapes the raw records; `rank_trials` orders by phase, enrollment, recency. `phase_summary` counts trials per phase × status. `terminated_table` subsets the stopped trials. `signals_table` validates hand-extracted signal rows. `write_run` calls them all and writes the CSVs plus the README. |
| `SKILL.md` | Workflow, principles, and verification steps. |

Module constants: `CTGOV_STUDIES_URL`, `CTGOV_STUDY_URL`, `TRIAL_COLUMNS`,
`STOPPED_STATUSES`, `SIGNAL_TYPES`, `SIGNAL_COLUMNS`, `PHASE_ORDER`.

## Running it

`python` cell — once per session, `pip install -e .` from the repo root:

```python
from lib.common.paths import results_dir
from kernel import search_trials, trials_table, rank_trials, write_run

hits = search_trials(drug)
df = rank_trials(trials_table(hits))

# read eligibility_criteria and arms for top-ranked + stopped trials,
# then build signal_rows by hand:
signal_rows = [
    {"signal_id": "S1", "drug": drug, "signal_type": "biomarker_requirement",
     "value": "...", "evidence_text": "...", "nct_id": "...", "url": "..."},
]

out = results_dir("clinical-trials", drug.lower())
write_run(hits, signal_rows, out, drug, summary, steps)
```
