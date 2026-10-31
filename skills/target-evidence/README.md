# target-evidence

Builds a target dossier from four sources and writes it as CSVs under
`results/target_evidence/<gene>/`. Open Targets supplies identity, the
top 20 disease associations, genetic evidence, known drugs and clinical
status, tractability and safety liabilities. gnomAD supplies
loss-of-function constraint (LOEUF, observed vs expected LoF), the
carrier counts behind it, and missense variants above a caller-supplied
frequency cutoff. cBioPortal supplies the gene's mutations in
MSI-high tumour samples with their variant allele frequencies. DepMap,
through `depmap-local`, supplies the MSI dose-response. Only the Open
Targets section is required; each other section is written when its data
is passed to `write_run`.

| File | Purpose |
|---|---|
| `kernel.py` | Open Targets table functions (`identity_table`, `associations_table`, `top_association_genetic_scores_table`, `known_drugs_table`, `drug_summary_table`, `tractability_table`, `safety_table`); gnomAD fetch and tables (`fetch_gnomad_gene`, `gnomad_constraint_table`, `gnomad_lof_variants_table`, `gnomad_carrier_summary_table`, `gnomad_missense_variants_table`); cBioPortal fetch and tables (`fetch_cbioportal_msi_cohort` and the per-step functions it composes); `msi_dose_response` delegating to `depmap-local`; `write_run` writes every supplied section plus the README. |
| `query.graphql` | The GraphQL query run in the `repl` cell. |
| `SKILL.md` | Workflow, principles, and verification steps. |

## Running it

`repl` cell — repl can't import this repo, so the handoff write is inline:

```python
data = host.mcp("clinical-genomics", "open_targets_graphql",
                query=open("query.graphql").read(),
                variables={"ensemblId": ensembl_id})
import json, os
os.makedirs("handoff", exist_ok=True)
json.dump(data, open(f"handoff/open_targets_{gene}.json", "w"))
```

`python` cell — once per session, `pip install -e .` from the repo root:

```python
from lib.common.handoff import read_handoff
from lib.common.paths import results_dir
from kernel import (write_run, fetch_gnomad_gene, fetch_cbioportal_msi_cohort,
                    msi_dose_response, MSI_ATTRIBUTE_THRESHOLDS)

data = read_handoff(f"open_targets_{gene}")
gnomad = fetch_gnomad_gene(gene)

# Residue ranges are caller-supplied and must be cited in the findings.
ranges = {"RecQ helicase ATP-binding": (558, 724)}   # UniProt Q14191
cbio = fetch_cbioportal_msi_cohort(gene, ranges, MSI_ATTRIBUTE_THRESHOLDS,
                                   study_ids=["coadread_tcga_pan_can_atlas_2018"])

out = results_dir("target-evidence", gene.lower())
write_run(data, out, subject, summary, steps,
          gnomad=gnomad, missense_min_allele_frequency=min_af,  # caller-chosen, cite it
          cbioportal=cbio, msi_dose=msi_dose_response(gene))
```

## Sources

| Source | Endpoint | Auth |
|---|---|---|
| Open Targets | `api.platform.opentargets.org` GraphQL | none |
| gnomAD | `gnomad.broadinstitute.org/api` GraphQL | none |
| cBioPortal | `www.cbioportal.org/api` REST | none |
| DepMap | local release via `$DEPMAP_ROOT` | n/a |
