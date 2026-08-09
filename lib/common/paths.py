"""Resolve and create results/<skill_name>/<subject>/<timestamp>/ run directories.

Project root comes from the DRUG_TARGETS_ROOT env var, falling back to
the value in .env.example (see _read_env_example) when unset.
"""

from __future__ import annotations

import os
import re
from datetime import datetime
from pathlib import Path

ROOT_ENV = "DRUG_TARGETS_ROOT"

# snake_case only: blocks a "/" or ".." component from escaping results/.
_PART_RE = re.compile(r"[a-z0-9][a-z0-9_]*\Z")


def _read_env_example() -> str | None:
    """DRUG_TARGETS_ROOT out of .env.example, next to this package's checkout.

    Used only when the env var itself is unset. The Claude Science sandbox
    cannot read this repo's .env, but it can read .env.example, so the
    fallback value lives there instead of a literal in this file. Located
    via __file__ rather than cwd: pip install -e . keeps this file at its
    real source location even when invoked from elsewhere.
    """
    env_file = Path(__file__).resolve().parents[2] / ".env.example"
    if not env_file.is_file():
        return None
    for line in env_file.read_text().splitlines():
        line = line.strip()
        if line.startswith(f"{ROOT_ENV}="):
            return line.split("=", 1)[1].strip()
    return None


def project_root() -> Path:
    """Repo root from DRUG_TARGETS_ROOT, falling back to .env.example when unset."""
    value = os.environ.get(ROOT_ENV) or _read_env_example()
    if not value:
        raise RuntimeError(f"{ROOT_ENV} is not set and .env.example has no fallback")
    root = Path(value).expanduser()
    if not root.is_dir():
        raise RuntimeError(f"{ROOT_ENV}={root} does not exist or is not a directory")
    return root.resolve()


def results_dir(skill_name: str, subject_name: str) -> Path:
    """results/<skill_name>/<subject_name>/<YYYY-MM-DD_HHMMSS>/, created if missing.

    skill_name matches the skills/ folder naming (kebab-case, e.g.
    "open-targets") and is converted to snake_case for the results path.
    subject_name must already be snake_case.
    """
    domain_name = skill_name.replace("-", "_")
    if not _PART_RE.fullmatch(domain_name):
        raise ValueError(f"invalid skill_name {skill_name!r}: use kebab-case or snake_case")
    if not _PART_RE.fullmatch(subject_name):
        raise ValueError(f"invalid subject_name {subject_name!r}: use snake_case")
    timestamp = datetime.now().strftime("%Y-%m-%d_%H%M%S")
    path = project_root() / "results" / domain_name / subject_name / timestamp
    path.mkdir(parents=True, exist_ok=True)
    return path
