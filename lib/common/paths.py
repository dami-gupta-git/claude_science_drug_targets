from __future__ import annotations

import os
from datetime import datetime
from pathlib import Path

ROOT_ENV = "DRUG_TARGETS_ROOT"

# Used only when DRUG_TARGETS_ROOT is unset. The Claude Science sandbox
# cannot read this repo's .env, so the env var never arrives there by
# itself. This repo lives on one machine, so a hardcoded fallback is
# simpler than a setup step every session needs to repeat.
_FALLBACK_ROOT = "/Users/dgupta/code/portfolio/claude_science/claude_science_drug_targets"


def project_root() -> Path:
    """Repo root from DRUG_TARGETS_ROOT, falling back to _FALLBACK_ROOT when unset."""
    value = os.environ.get(ROOT_ENV) or _FALLBACK_ROOT
    root = Path(value).expanduser()
    if not root.is_dir():
        raise RuntimeError(f"{ROOT_ENV}={root} does not exist or is not a directory")
    return root.resolve()


def results_dir(skill_name: str, subject_name: str) -> Path:
    """results/<skill_name>/<subject_name>/<YYYY-MM-DD_HHMMSS>/, created if missing.

    skill_name matches the skills/ folder naming (kebab-case, e.g.
    "open-targets") and is converted to snake_case for the results path.
    """
    timestamp = datetime.now().strftime("%Y-%m-%d_%H%M%S")
    domain_name = skill_name.replace("-", "_")
    path = project_root() / "results" / domain_name / subject_name / timestamp
    path.mkdir(parents=True, exist_ok=True)
    return path
