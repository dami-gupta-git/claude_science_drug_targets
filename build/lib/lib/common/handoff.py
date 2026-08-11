"""Read a repl->python handoff file.

repl can't import this module (no filesystem grant), so it writes
handoff/<name>.json inline: os.makedirs("handoff", exist_ok=True);
json.dump(data, open(f"handoff/{name}.json", "w")). Path is relative —
the workspace path has a per-session ID.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

_HANDOFF_DIR = Path("handoff")


def read_handoff(name: str) -> Any:
    """Read back handoff/<name>.json written by a repl cell."""
    path = _HANDOFF_DIR / f"{name}.json"
    if not path.exists():
        raise FileNotFoundError(
            f"{path.resolve()} not found — did the repl cell that writes "
            f"'{name}' run in this same session?"
        )
    return json.loads(path.read_text())
