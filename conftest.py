import sys
from pathlib import Path

# Not redundant with `pip install -e .`: this repo isn't actually
# installed in the environment tests run in (confirmed - `import lib`
# fails without this line), only in whatever kernel a skill runs in.
sys.path.insert(0, str(Path(__file__).resolve().parent))

# A local .env (secrets), if one is ever added, should not be collected.
collect_ignore_glob = [".env*"]
