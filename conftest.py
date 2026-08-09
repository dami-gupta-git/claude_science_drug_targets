"""Root conftest: makes `lib` importable and keeps collection off `.env`.

The `.env` exclusion matters because some sandboxed environments deny
`stat` on it, which fails collection before any test runs.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

collect_ignore_glob = [".env*"]
