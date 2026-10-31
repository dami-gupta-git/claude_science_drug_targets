"""Shared README and markdown-table utilities for skill kernels."""

from __future__ import annotations

from pathlib import Path


def format_cell(value) -> str:
    """Render a table cell; floats are shortened for display, CSVs keep full precision.

    |x| < 1e-3 → scientific with 3 significant figures (keeps tiny p-values
    distinct from zero); |x| < 1 → 3 significant figures; otherwise 3 decimals.
    """
    if isinstance(value, bool) or not isinstance(value, float):
        return str(value)
    if value != value or value == 0:
        return str(value)
    if abs(value) < 1e-3:
        return f"{value:.2e}"
    if abs(value) < 1:
        return f"{value:.3g}"
    return str(round(value, 3))


def markdown_table(rows: list[dict], headers: list[str] | None = None) -> str:
    """Render a list of dicts as a markdown table."""
    if not rows:
        return ""
    keys = headers or list(rows[0].keys())
    sep = " | ".join("---" for _ in keys)
    header = " | ".join(keys)
    lines = [f"| {header} |", f"| {sep} |"]
    for row in rows:
        cells = " | ".join(format_cell(row.get(k, "")) for k in keys)
        lines.append(f"| {cells} |")
    return "\n".join(lines)


def write_readme(path: Path, subject: str, summary: str, steps: list[dict],
                 source: str | None = None, limits: list[str] | None = None) -> Path:
    """Write a run README from caller-supplied findings.

    steps: [{"name": str, "finding": str, "table": list[dict] (opt),
             "skipped": str (opt)}]. Use "skipped" when a section could not
    be run so the reader can distinguish it from "not checked".
    source: appended to the title as "# {subject} — {source}" when given.
    """
    if not steps:
        raise ValueError("write_readme requires at least one step")
    if not summary or not str(summary).strip():
        raise ValueError("summary is required")

    title = f"# {subject} — {source}" if source else f"# {subject}"
    parts = [title, "", str(summary).strip(), ""]

    for step in steps:
        name = step.get("name")
        if not name:
            raise ValueError(f"step {step!r} missing 'name'")
        finding = step.get("finding")
        skipped = step.get("skipped")
        if not finding and not skipped:
            raise ValueError(
                f"step '{name}' needs 'finding' or 'skipped'; "
                "state what it showed or why it could not be run"
            )
        parts += [f"## {name}", ""]
        if skipped:
            parts += [f"Not run. {str(skipped).strip()}", ""]
        else:
            parts += [str(finding).strip(), ""]
            table = step.get("table")
            if table:
                parts += [markdown_table(table, step.get("headers")), ""]

    if limits:
        parts += ["## Limits", ""]
        parts += [f"- {s}" for s in limits]
        parts.append("")

    path = Path(path)
    path.write_text("\n".join(parts).rstrip() + "\n")
    return path
