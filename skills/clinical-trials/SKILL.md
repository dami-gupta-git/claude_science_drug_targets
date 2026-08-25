---
name: clinical-trials
description: Fetch a drug's clinical trials from ClinicalTrials.gov, analyze phase distribution and terminations, and extract population signals (biomarker requirements, prior-therapy requirements, combination partners) from eligibility text with a source quote and NCT ID for each. Use when asked what trials exist for a drug, how far it got in the clinic, why trials stopped, or what population sponsors selected.
---

# Clinical trials

## Workflow

- `search_trials(drug, condition=None, max_records=100)` is the entry
  point — network call, run in the `python` cell directly
  (ClinicalTrials.gov v2 is a plain REST API, no key, no MCP connector).
  The API matches registered synonyms and brand names itself; do not
  pre-expand the drug name. Pass `condition` only when the question is
  about one disease.
- `rank_trials(trials_table(hits))` orders by importance: phase, then
  enrollment, then recency. Unphased trials (observational, expanded
  access) sort last with a null `phase_rank` — report them separately,
  never bucket them into a phase.
- `phase_summary(df)` gives trial counts per phase × status; report the
  highest phase reached and whether any late-phase trial is still active.
- `terminated_table(df)` is the stopped subset. Read every row — a
  termination is a finding regardless of rank. Classify the reason
  (safety, efficacy, business, enrollment) only when `why_stopped` text
  says so; an empty reason is reported as "reason not given", never
  inferred (see `terminated_table` docstring).
- Read `eligibility_criteria` and `arms` for the top-ranked trials plus
  all stopped trials — eligibility text runs long, so state the read
  cutoff in the run's limits. For each population signal, write a row
  with `signal_type` (one of `SIGNAL_TYPES`), `value` (the biomarker and
  threshold, the required prior therapy, or the combination partner
  drug), and `evidence_text` set to the exact passage it came from — not
  a paraphrase.
- `signals_table(signal_rows)` raises if a required field is missing or
  `signal_type` isn't recognized — fix the row, don't work around it.
- `write_run(hits, signal_rows, out, subject, summary, steps)` writes
  `trials.csv`, `phase_summary.csv`, `terminated.csv`,
  `population_signals.csv`, and `README.md` in one call. Invoke
  `docs-style` before writing `summary`/`steps` prose.

## Principles

- A signal is a fact about a *passage in a specific trial's registration*,
  not a property of the drug. Two trials requiring the same biomarker are
  two rows — convergence across independent trials is the evidence.
- Eligibility criteria, prior-therapy requirements, combination arms and
  enrollment show what sponsors think the population is — independent
  evidence on the biomarker question, separate from literature or omics.
- Never record a signal from a trial's title alone; titles advertise,
  eligibility text binds.
- An `enrollment_type` of ESTIMATED is a plan, not a result — say which
  it is when citing enrollment.
- Single writer: only the `python` cell writes to `results/`.

## Writing findings

Follow `docs-style`. Trial-specific guidance:

- Group signals by `signal_type`, not by trial; cite the NCT ID next to
  each value.
- Compare early- vs late-phase eligibility: a sponsor that narrowed to a
  biomarker-selected population after phase 1 is a finding, as is one
  that never selected.
- Flag enrollment inconsistent with the claimed population — a large
  all-comers phase 3 after a biomarker-selected phase 2 needs remarking.
- If trials were fetched but no signal was recordable, say so — that is
  a finding, not a skipped step.

## Verify

- `population_signals.csv` has no row whose `nct_id` is absent from
  `trials.csv`.
- Every signal's `evidence_text` is a substring of that trial's
  `eligibility_criteria` or `arms` in `trials.csv` — spot-check a sample
  before trusting the table.
- Every `terminated.csv` row is mentioned in the findings; a stopped
  trial left unremarked is a gap.
- If `max_records` was hit, the fetch is truncated — say so in limits
  and consider narrowing with `condition` rather than trusting a
  partial phase distribution.
