import sys
from pathlib import Path

# skills/biomarker-stratification/ has hyphens, so it can't be a Python
# package - add it to sys.path directly so tests can `import kernel`.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
