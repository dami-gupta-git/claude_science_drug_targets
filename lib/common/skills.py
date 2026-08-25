"""Import one skill's kernel.py from another skill."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

from lib.common.paths import project_root


def load_sibling_kernel(skill_dir: str):
    """Import skills/<skill_dir>/kernel.py as the module <skill_dir>_kernel.

    Skill directories are hyphenated, so they are not importable as normal
    Python packages; this loads one by file path. Repeated calls re-use the
    module already in sys.modules.
    """
    module_name = skill_dir.replace("-", "_") + "_kernel"
    if module_name in sys.modules:
        return sys.modules[module_name]
    kernel_path = project_root() / "skills" / skill_dir / "kernel.py"
    if not kernel_path.is_file():
        raise FileNotFoundError(f"no kernel.py for skill {skill_dir!r} at {kernel_path}")
    spec = importlib.util.spec_from_file_location(module_name, kernel_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module
