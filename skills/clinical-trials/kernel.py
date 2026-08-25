"""Fetch and shape a drug's clinical trials from the ClinicalTrials.gov v2 API.

kernel.py loads with this skill: search_trials, trial_row, trials_table,
rank_trials, phase_summary, terminated_table, signals_table, write_run,
phase_rank, and the constants CTGOV_STUDIES_URL, CTGOV_STUDY_URL,
TRIAL_COLUMNS, STOPPED_STATUSES, SIGNAL_TYPES, SIGNAL_COLUMNS, PHASE_ORDER.

The API matches drug synonyms and brand names on the intervention query
itself (e.g. "sotorasib" also matches trials registered under "AMG 510"),
so no separate synonym expansion is done here.
"""
import json
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

CTGOV_STUDIES_URL = "https://clinicaltrials.gov/api/v2/studies"
CTGOV_STUDY_URL = "https://clinicaltrials.gov/study/{nct_id}"

TRIAL_COLUMNS = [
    "nct_id", "title", "phase", "study_type", "overall_status", "why_stopped",
    "conditions", "interventions", "arms", "enrollment", "enrollment_type",
    "sponsor", "start_date", "primary_completion_date", "has_results",
    "eligibility_criteria", "url",
]

# Registry statuses meaning the trial stopped before finishing as planned.
STOPPED_STATUSES = ("TERMINATED", "WITHDRAWN", "SUSPENDED")

# Registry phase codes, highest first. A trial spanning two phases (the API
# returns both, e.g. PHASE1 and PHASE2) ranks by the highest one. NA marks
# trials the registry classes as not applicable (e.g. observational).
PHASE_ORDER = ("PHASE4", "PHASE3", "PHASE2", "PHASE1", "EARLY_PHASE1", "NA")

SIGNAL_TYPES = ("biomarker_requirement", "prior_therapy", "combination_partner")

SIGNAL_COLUMNS = [
    "signal_id", "drug", "signal_type", "value",
    "evidence_text", "nct_id", "url",
]


def search_trials(drug: str, condition: str | None = None,
                  max_records: int = 100, timeout: int = 60) -> list[dict]:
    """Raw ClinicalTrials.gov study records for an intervention, paginated.

    drug goes to the API's intervention query, which matches registered
    synonyms and brand names on its own. condition, if given, restricts to
    trials in that disease. Stops at max_records or when the API returns
    no further page token.
    """
    hits = []
    page_token = None
    while len(hits) < max_records:
        params = {
            "query.intr": drug,
            "pageSize": min(100, max_records - len(hits)),
            "format": "json",
        }
        if condition:
            params["query.cond"] = condition
        if page_token:
            params["pageToken"] = page_token
        url = "%s?%s" % (CTGOV_STUDIES_URL, urllib.parse.urlencode(params))
        with urllib.request.urlopen(url, timeout=timeout) as fh:
            data = json.load(fh)
        hits.extend(data.get("studies") or [])
        page_token = data.get("nextPageToken")
        if not page_token:
            break
    return hits


def trial_row(study: dict) -> dict:
    """One row from a raw v2 study record."""
    protocol = study.get("protocolSection") or {}
    ident = protocol.get("identificationModule") or {}
    status = protocol.get("statusModule") or {}
    design = protocol.get("designModule") or {}
    arms_mod = protocol.get("armsInterventionsModule") or {}
    elig = protocol.get("eligibilityModule") or {}
    sponsor = ((protocol.get("sponsorCollaboratorsModule") or {})
               .get("leadSponsor") or {})
    enrollment = design.get("enrollmentInfo") or {}

    nct_id = ident.get("nctId")
    phases = design.get("phases") or []
    conditions = (protocol.get("conditionsModule") or {}).get("conditions") or []
    interventions = arms_mod.get("interventions") or []
    arms = arms_mod.get("armGroups") or []
    arm_texts = [
        "%s: %s" % (arm.get("label"),
                    " + ".join(arm.get("interventionNames") or []))
        for arm in arms
    ]
    return {
        "nct_id": nct_id,
        "title": ident.get("briefTitle"),
        "phase": "/".join(phases) if phases else None,
        "study_type": design.get("studyType"),
        "overall_status": status.get("overallStatus"),
        "why_stopped": status.get("whyStopped"),
        "conditions": "; ".join(conditions) if conditions else None,
        "interventions": "; ".join(
            i.get("name") for i in interventions if i.get("name")) or None,
        "arms": "; ".join(arm_texts) if arm_texts else None,
        "enrollment": enrollment.get("count"),
        # ACTUAL vs ESTIMATED - an estimated count is a plan, not a result
        "enrollment_type": enrollment.get("type"),
        "sponsor": sponsor.get("name"),
        "start_date": (status.get("startDateStruct") or {}).get("date"),
        "primary_completion_date":
            (status.get("primaryCompletionDateStruct") or {}).get("date"),
        "has_results": study.get("hasResults"),
        "eligibility_criteria": elig.get("eligibilityCriteria"),
        "url": CTGOV_STUDY_URL.format(nct_id=nct_id) if nct_id else None,
    }


def trials_table(hits: list[dict]) -> pd.DataFrame:
    """One row per trial, in TRIAL_COLUMNS order."""
    rows = [trial_row(study) for study in hits]
    return pd.DataFrame(rows, columns=TRIAL_COLUMNS)


def phase_rank(phase: str | None) -> int | None:
    """Ordering rank for a trials_table phase value; higher is later-stage.

    A joined value like "PHASE1/PHASE2" ranks by its highest component.
    Returns None when the phase is absent or NA - the registry did not
    assign one, so the trial has no place in a phase ordering.
    """
    if not phase or not isinstance(phase, str):
        return None
    ranks = [len(PHASE_ORDER) - 1 - PHASE_ORDER.index(p)
             for p in phase.split("/") if p in PHASE_ORDER and p != "NA"]
    return max(ranks) if ranks else None


def rank_trials(df: pd.DataFrame) -> pd.DataFrame:
    """trials_table sorted by importance: phase, then enrollment, then recency.

    Adds a phase_rank column (None for unphased trials, which sort last
    within the table but are not dropped - being unphased is information,
    not noise). Sorting is presentation only; every input row is kept.
    """
    out = df.copy()
    out["phase_rank"] = out["phase"].map(phase_rank)
    return out.sort_values(
        by=["phase_rank", "enrollment", "start_date"],
        ascending=False, na_position="last",
    ).reset_index(drop=True)


def phase_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Trial counts per (phase, overall_status), unphased trials included."""
    counts = (df.groupby(["phase", "overall_status"], dropna=False)
                .size().reset_index(name="n_trials"))
    counts["phase_rank"] = counts["phase"].map(phase_rank)
    counts = counts.sort_values(
        by=["phase_rank", "overall_status"],
        ascending=[False, True], na_position="last",
    ).drop(columns="phase_rank")
    return counts.reset_index(drop=True)


def terminated_table(df: pd.DataFrame) -> pd.DataFrame:
    """The stopped subset (STOPPED_STATUSES) with the sponsor's stated reason.

    why_stopped is free text the sponsor optionally supplies; when absent
    it stays empty. Do not infer a reason - report it as not given.
    """
    cols = ["nct_id", "title", "phase", "overall_status", "why_stopped",
            "enrollment", "start_date", "url"]
    stopped = df[df["overall_status"].isin(STOPPED_STATUSES)]
    return stopped[cols].reset_index(drop=True)


def signals_table(signal_rows: list[dict]) -> pd.DataFrame:
    """One row per population signal read from a trial's eligibility text.

    Every signal must carry an nct_id and evidence_text - the exact
    passage it was read from in that trial's eligibility criteria or arm
    descriptions. A signal without a locatable source passage must not be
    recorded. signal_type must be one of SIGNAL_TYPES.
    """
    for row in signal_rows:
        missing = [c for c in ("drug", "signal_type", "value",
                               "evidence_text", "nct_id")
                   if not row.get(c)]
        if missing:
            raise ValueError(
                "signal %r missing required field(s): %s" % (row, missing))
        if row["signal_type"] not in SIGNAL_TYPES:
            raise ValueError(
                "signal_type %r not in %s" % (row["signal_type"], SIGNAL_TYPES))
    return pd.DataFrame(signal_rows, columns=SIGNAL_COLUMNS)


def markdown_table_local(rows: list[dict], headers: list[str] | None = None) -> str:
    """A markdown table from a list of dicts, used by the fallback writer."""
    if not rows:
        return ""
    cols = headers or list(rows[0].keys())
    out = ["| " + " | ".join(cols) + " |",
           "| " + " | ".join("---" for _ in cols) + " |"]
    for row in rows:
        out.append("| " + " | ".join(str(row.get(c, "")) for c in cols) + " |")
    return "\n".join(out)


def write_readme_local(path, subject: str, summary: str, steps: list[dict],
                       source: str = "", limits: list[str] | None = None):
    """Standalone README writer, used when lib.common.reporting is absent.

    Mirrors the repo writer's structure and its validation: a summary and at
    least one named step are required, so a run cannot be written without
    them whichever writer is in use.
    """
    if not steps:
        raise ValueError("write_readme requires at least one step")
    if not summary or not str(summary).strip():
        raise ValueError("summary is required")
    title = "# %s — %s" % (subject, source) if source else "# %s" % subject
    parts = [title, "", str(summary).strip(), ""]
    for step in steps:
        if not step.get("name"):
            raise ValueError("step %r missing 'name'" % step)
        if not step.get("finding"):
            raise ValueError("step %r missing 'finding'" % step)
        parts += ["## %s" % step["name"], "", str(step["finding"]).strip(), ""]
        if step.get("table"):
            parts += [markdown_table_local(step["table"], step.get("headers")), ""]
    if limits:
        parts += ["## Limits", ""] + ["- %s" % lim for lim in limits] + [""]
    path = Path(path)
    path.write_text("\n".join(parts))
    return path


def resolve_readme_writer():
    """The repo's write_readme when importable, else the local fallback."""
    try:
        from lib.common.reporting import write_readme
        return write_readme
    except ImportError:
        return write_readme_local


def write_run(hits: list[dict], signal_rows: list[dict], out_dir: Path,
              subject: str, summary: str, steps: list[dict],
              source: str = "ClinicalTrials.gov",
              limits: list[str] | None = None) -> dict[str, Path]:
    """Write trials.csv, phase_summary.csv, terminated.csv,
    population_signals.csv and README.md into out_dir.

    The README writer is imported here rather than at module scope. In this
    repo it comes from lib.common.reporting, which is only importable with
    the repo root on sys.path; the published skill loads its sidecar from
    the skill directory instead, where `lib` does not exist. Importing at
    module scope therefore broke the whole sidecar outside the repo, taking
    the fetch and table helpers down with it. The fallback keeps every
    other function usable and fails only on this call.
    """
    write_readme = resolve_readme_writer()
    written = {}
    df = trials_table(hits)

    trials_path = out_dir / "trials.csv"
    rank_trials(df).to_csv(trials_path, index=False)
    written["trials"] = trials_path

    summary_path = out_dir / "phase_summary.csv"
    phase_summary(df).to_csv(summary_path, index=False)
    written["phase_summary"] = summary_path

    terminated_path = out_dir / "terminated.csv"
    terminated_table(df).to_csv(terminated_path, index=False)
    written["terminated"] = terminated_path

    signals_path = out_dir / "population_signals.csv"
    signals_table(signal_rows).to_csv(signals_path, index=False)
    written["signals"] = signals_path

    written["readme"] = write_readme(out_dir / "README.md", subject, summary,
                                     steps, source=source, limits=limits)
    return written
