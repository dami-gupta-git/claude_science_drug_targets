# claude_science_drug_targets

Claude Science skills for drug target discovery: gathering evidence on a
gene from public databases and local data releases, and writing it up as
CSVs plus a findings README.

## Skills

| Skill | Purpose |
|---|---|
| [`skills/target-evidence`](skills/target-evidence/README.md) | Target dossier: identity, disease associations, genetic evidence, drugs, tractability and safety from Open Targets; population constraint and loss-of-function carrier counts from gnomAD; the gene's mutations in MSI-high tumours from cBioPortal; and the MSI dose-response from DepMap. |
| [`skills/depmap-local`](skills/depmap-local/README.md) | Gene dependency analysis against a locally downloaded DepMap release: CRISPR gene effect, lineage enrichment, genotype/fusion contrasts, PRISM compound sensitivity. |
| [`skills/depmap-ot-fusion`](skills/depmap-ot-fusion/SKILL.md) | Fuses a DepMap dependency run and an Open Targets run into a per-gene dossier: verdict, conflicts between sources, disease-lineage alignment, and a drug join table. |
| [`skills/ai-origination-audit`](skills/ai-origination-audit/SKILL.md) | Audits how much of a drug candidate's discovery was actually done by AI, tiering claims against clinical trial registry evidence with provenance. |
| [`skills/literature-evidence`](skills/literature-evidence/README.md) | Searches Europe PMC (including bioRxiv/medRxiv preprints) for a gene's literature and extracts disease-association, interaction, mechanism, and biomarker claims with a source quote and PMID for each. |
| [`skills/multi-omics-expression`](skills/multi-omics-expression/README.md) | Median gene expression across GTEx normal tissues and a fixed set of major TCGA cancer types (via cBioPortal), flagging markedly over- or under-expressed tissues/cancer types. |
| [`skills/biomarker-stratification`](skills/biomarker-stratification/SKILL.md) | Biomarkers used to select or stratify patients in a disease's clinical trials, from ClinicalTrials.gov eligibility criteria and Europe PMC literature. Scoped by disease, not gene. |
| [`skills/knowledge-graph`](skills/knowledge-graph/README.md) | Persistent gene/disease/drug/lineage graph assembled from `target-evidence`, `depmap-local`, and `literature-evidence` run output; answers path/neighbor queries. |
| [`skills/docs-style`](skills/docs-style/SKILL.md) | Shared tone and writing conventions for findings prose and other documents; not a data skill, invoked by the others. |

## Workflow

Ten live skills (plus `open-targets`, which is stale). Which apply depends on what you're starting from:

- Start from a gene — `target-evidence`, `depmap-local`, `literature-evidence`, `multi-omics-expression`, `depmap-ot-fusion`
- Start from a drug — `clinical-trials`, `ai-origination-audit`
- Start from a disease — `biomarker-stratification`
- Across runs — `knowledge-graph`, and `docs-style` for anything written

### Triage sequence for a gene

1. **`target-evidence`** — run this first, always. It's the broadest single call: identity, whether the protein is chemically druggable, whether losing it is survivable in the population, whether resistance mutations already exist in tumours, and whether the dependency is graded or a threshold. If the answer here is discouraging, stop — one run is enough to know.
2. **`depmap-local`** — the fuller dependency picture. `target-evidence` borrows only the dose-response from it, so this adds which cancer lineages depend on the gene, the mismatch-repair contrast, and whether any existing compound hits it. Run it when step 1 shows a dependency worth characterising.
3. **`literature-evidence`** — what the field already claims, with a quoted sentence and a PMID behind each claim. Run this before forming a view of your own; it's how you find out the result was published in 2019.
4. **`multi-omics-expression`** — a cheap sanity check: is the gene even expressed in the tissues where the dependency showed up. Fast, and occasionally contradicts steps 2 and 3.
5. **`depmap-ot-fusion`** — reconciles the dependency and Open Targets views and reports where they disagree. It reuses the runs from steps 1 and 2 rather than refetching, so it only makes sense after them.
6. **`knowledge-graph`** — connects the runs so you can ask which of your genes share a disease. Worth it once you have several genes; for a single gene the individual runs are simpler.

### Then the drug side

- **`clinical-trials`**, if a compound exists — who is running trials, what phase, and how they're selecting patients. On WRN this found four sponsors taking three different positions on patient selection, which no other source showed.
- **`biomarker-stratification`**, if you're asking how patients would be chosen for an indication rather than about one drug.
- **`ai-origination-audit`**, only when the provenance of a discovery matters. It audits a drug candidate, not a gene — hand it a gene and it has nothing to work with until you resolve that gene to actual compounds first.

### Caveats

- `biomarker-stratification` is the one skill in the suite that hasn't actually been run. Its code was read and fixed, but the recommendation above rests on its documentation rather than on observed output.
- For `knowledge-graph`, the existing store was read and repaired, but never freshly synced. The store on disk predates the `support_p` fix, so the query functions will now refuse it by design — it needs re-ingesting from the source runs before path queries will work.
- `target-evidence` needs `DEPMAP_ROOT` set for its dose-response section, and `literature-evidence` should be called through `search_balanced` rather than the single-pass search, or it returns only recent papers.

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
