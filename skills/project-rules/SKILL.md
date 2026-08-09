---
name: project-rules
description: All standing rules for the drug-targets project — not just one topic. Load at the start of any session in this project and before writing any script, saving any result, or deciding where a file goes. New project-wide conventions are added here as sections rather than split into separate skills.
---

# Project rules

- The standing rulebook for the drug-targets project. Any convention the
  user settles that applies project-wide goes here as a new section
  rather than into a topic-specific skill.
- Content should be expressed as bullet points, not prose paragraphs.

## Path resolution

- `lib/common/paths.py` owns where things live. Call `project_root()`
  and `results_dir(skill_name, subject_name)` to resolve paths.
- Do not compose paths by hand, and do not restate their shape here —
  the module's docstrings are the specification.
- `lib/` cannot be imported until it's on `sys.path`, so every `python`
  cell that touches `results_dir` starts with:
  ```python
  import sys
  sys.path.insert(0, "/Users/dgupta/code/portfolio/claude_science/claude_science_drug_targets")
  from lib.common.paths import results_dir
  ```
- Writes go straight to `results_dir(...)`'s path from that same
  `python` cell — the sandbox has a filesystem grant on this repo, no
  `bash` or `save_artifacts` needed.
- Exception: MCP connectors (`host.mcp`) only run in `repl`, which has
  no filesystem grant. A connector fetch writes to the shared workspace
  first, then a `python` cell copies it into the run directory.

## Skill layout

- Skills live in the top-level `skills/` folder, one directory per
  skill, named in kebab-case (e.g. `skills/open-targets/`).

## Coding rules

- Do not duplicate code.
- Code docstrings should be brief and tight.

## SKILL.md files

- Should only have: high-level workflow and steps; principles, rules,
  constraints; decision criteria and judgement calls; how to verify the
  result.
- Should be lean. Move detailed reference material, long examples, or
  schemas into separate files.
- Content should be expressed as bullet points, not prose paragraphs.
- They should be brief and tight.
