# claude_science_drug_targets

Claude Science skills for drug target discovery: gathering evidence on a
gene from public databases and local data releases, and writing it up as
CSVs plus a findings README.

## Skills

| Skill | Purpose |
|---|---|
| [`skills/open-targets`](skills/open-targets/README.md) | Target identity, disease associations, genetic evidence, known drugs, tractability, and safety liabilities from the Open Targets Platform GraphQL API. |
| [`skills/depmap-local`](skills/depmap-local/README.md) | Gene dependency analysis against a locally downloaded DepMap release: CRISPR gene effect, lineage enrichment, genotype/fusion contrasts, PRISM compound sensitivity. |
| [`skills/depmap-ot-fusion`](skills/depmap-ot-fusion/SKILL.md) | Fuses a DepMap dependency run and an Open Targets run into a per-gene dossier: verdict, conflicts between sources, disease-lineage alignment, and a drug join table. |
| [`skills/ai-origination-audit`](skills/ai-origination-audit/SKILL.md) | Audits how much of a drug candidate's discovery was actually done by AI, tiering claims against clinical trial registry evidence with provenance. |
| [`skills/literature-evidence`](skills/literature-evidence/README.md) | Searches Europe PMC (including bioRxiv/medRxiv preprints) for a gene's literature and extracts disease-association, interaction, mechanism, and biomarker claims with a source quote and PMID for each. |
| [`skills/multi-omics-expression`](skills/multi-omics-expression/README.md) | Median gene expression across GTEx normal tissues and a fixed set of major TCGA cancer types (via cBioPortal), flagging markedly over- or under-expressed tissues/cancer types. |
| [`skills/biomarker-stratification`](skills/biomarker-stratification/SKILL.md) | Biomarkers used to select or stratify patients in a disease's clinical trials, from ClinicalTrials.gov eligibility criteria and Europe PMC literature. Scoped by disease, not gene. |
| [`skills/knowledge-graph`](skills/knowledge-graph/README.md) | Persistent gene/disease/drug/lineage graph assembled from `open-targets`, `depmap-local`, and `literature-evidence` run output; answers path/neighbor queries. |
| [`skills/docs-style`](skills/docs-style/SKILL.md) | Shared tone and writing conventions for findings prose and other documents; not a data skill, invoked by the others. |

## Setup

```bash
pip install -e .
```

Each data skill's own README has skill-specific setup (API access,
`$DEPMAP_ROOT`, etc.).

## Tests

```bash
pytest
```
