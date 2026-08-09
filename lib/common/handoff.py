"""Read data a repl cell wrote to hand off to a python cell.

repl (has connectors, no filesystem grant) can't import this module, so
it writes the handoff file inline:

    import json, os
    os.makedirs("handoff", exist_ok=True)
    json.dump(data, open(f"handoff/{name}.json", "w"))

python (has the filesystem grant, no connectors) reads it back with
read_handoff. The path is relative to the session's working directory,
never a hardcoded absolute path, because the workspace path contains a
per-session ID.
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
