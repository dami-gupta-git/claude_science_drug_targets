import importlib.util
import sys
from pathlib import Path

# Load under a unique module name so this kernel does not collide with the
# 'kernel' modules other skills' conftest.py files register.
_spec = importlib.util.spec_from_file_location(
    "literature_search_terms_kernel", Path(__file__).resolve().parents[1] / "kernel.py")
_mod = importlib.util.module_from_spec(_spec)
sys.modules["literature_search_terms_kernel"] = _mod
_spec.loader.exec_module(_mod)
