# Literature and multi-omics skills: plan

Four new skills extend the existing gene- and drug-target evidence base into
literature, expression, and cross-source graph queries. Each follows the
shape already used by `open-targets`, `depmap-local`, and `depmap-ot-fusion`:
a `SKILL.md` for workflow and principles, a `README.md` describing the
functions in `kernel.py`, pure table-building functions plus a `write_run`
that writes CSVs and a README into `results/<skill>/<subject>/`, and a
`tests/` folder. All four use `lib.common.paths.results_dir` and
`lib.common.reporting.write_readme`, and invoke `docs-style` before writing
findings prose.

Build order: `literature-evidence` first, since `biomarker-stratification`
reuses its Europe PMC access rather than duplicating it, and
`knowledge-graph` ingests its output. `multi-omics-expression` has no
dependency on the others and can be built in any order relative to them.
`knowledge-graph` is last, because it consumes the run output the other
three produce.

## literature-evidence

Status: built.

Subject: a gene symbol. Source: the Europe PMC REST API, which needs no key
and is called directly from the `python` cell with the standard library's
HTTP client, the same pattern `ai-origination-audit` already uses for
ClinicalTrials.gov. No MCP connector is involved.

The skill searches Europe PMC for the gene and writes every matching
publication to a table. Claim extraction is a reading step: for each
publication worth recording, a row is written naming the claim type -
disease association, interaction, mechanism, or biomarker - the entity the
claim concerns, and the exact sentence from the abstract that states it. A
claim with no locatable source sentence is not recorded. This keeps every
claim traceable to a specific publication and a specific sentence, so two
publications making the same claim produce two rows, not one.

Output: a publications table and a claims table, one row per publication and
one row per claim respectively, in `results/literature-evidence/<gene>/`.

## multi-omics-expression

Status: built. Source ended up as GTEx's own portal API for normal-tissue
expression (not Expression Atlas as originally sketched) paired with
cBioPortal for tumor expression - GTEx returns per-tissue medians for a gene
in one call, which fits the same direct-call pattern more directly than
Expression Atlas's bulk-download-oriented API. A tissue/cancer type is
flagged "high" at >=4x the gene's own median across the table it appears
in, "low" at <=1/4x - fixed as `EXPRESSION_FOLD_THRESHOLD` in `kernel.py`,
the same role `depmap-local`'s `DEP_THRESHOLD` plays.

Subject: a gene symbol. Source: a public gene-expression API - Expression
Atlas for normal-tissue expression and cBioPortal for tumor expression by
cancer type - both called directly, following the same no-MCP pattern as
`literature-evidence`.

The skill writes two tables: expression level by normal tissue, and
expression level by cancer type. The finding is which tissues or cancer
types show expression markedly above or below the gene's typical level.
What counts as "markedly" is a threshold fixed in `SKILL.md`, not decided
per run, so two runs of the same gene reach the same call.

This skill has no dependency on the other three and sits next to
`depmap-local` as a second check on whether a target is biologically active
where the disease of interest occurs.

Output: a normal-tissue expression table and a tumor expression table, in
`results/multi-omics-expression/<gene>/`.

## biomarker-stratification

Status: built. `search_ctgov_by_condition` fills the condition-search gap
`ai-origination-audit`'s `search_ctgov` doesn't cover (it searches
intervention/sponsor only); trial full-record fetching still goes through
that skill's `fetch_ctgov_studies`, loaded by file path via
`load_sibling_kernel` since skill directories are hyphenated and can't be
imported as normal Python packages.

Subject: a disease name, not a gene - the only skill in this set scoped by
indication rather than target. Source: the same Europe PMC access used by
`literature-evidence` and the same ClinicalTrials.gov access used by
`ai-origination-audit`, imported and called from this skill's `kernel.py`
rather than reimplemented.

The skill searches literature and trial records for biomarkers already used
to select patients for trials in the given disease, and writes one row per
biomarker naming its evidence type - a registry-listed inclusion criterion,
or a literature claim - and its source citation.

This skill depends on `literature-evidence` and `ai-origination-audit`
existing as importable modules, not on their run output, so it can be built
any time after `literature-evidence` is in place.

Output: a biomarkers table, in `results/biomarker-stratification/<disease>/`.

## knowledge-graph

Status: built. Node types ended up as gene/disease/drug/lineage rather than
gene/disease/drug/pathway as originally sketched: `depmap-local` surfaces
lineage (tissue/cancer-type) dependency, not pathway membership, and no
pathway source is wired up yet. The update rule settled on is
replace-by-source: re-syncing a gene from a given upstream skill drops that
skill's prior edges for that gene and writes the new ones, rather than
accumulating duplicates across re-runs.

Subject: none fixed - this skill answers queries rather than producing a
single gene or disease dossier. It is the only one of the four that reads
other skills' `results/` output as input rather than calling an external
source itself, and the only one that needs a persistent store - accumulated
nodes and edges - rather than one dossier per run.

Nodes are genes, diseases, drugs, and pathways. Edges are drawn from
`open-targets` (target-disease, target-drug), `depmap-local`
(gene-dependency-lineage), and `literature-evidence` (claim-based edges,
each retaining its source PMID). The skill answers path and neighbor
queries for a gene or disease against the accumulated graph.

Because the store persists and accumulates across runs, `SKILL.md` for this
skill needs an explicit rule for whether a new run of an upstream skill
replaces or merges into existing edges, which the other three skills do not
need since each of their runs is self-contained. This is worked out when
the skill is built, after the other three have produced real output to
ingest.
