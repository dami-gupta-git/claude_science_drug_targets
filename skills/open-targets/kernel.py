"""Shape a raw Open Targets target-dossier response into result tables.

Pure functions: no network, no filesystem, except write_run's own CSV
writes. Every table carries the underlying Ensembl/ChEMBL/EFO ids
alongside display names, so a later skill can look up the same records
without re-fetching.
"""

from pathlib import Path

import pandas as pd

from lib.common.reporting import markdown_table, write_readme

TARGET_URL = "https://platform.opentargets.org/target/{id}"
DRUG_URL = "https://platform.opentargets.org/drug/{id}"
DISEASE_URL = "https://platform.opentargets.org/disease/{id}"

# The datatypeScores id that carries genetic-evidence strength.
GENETIC_DATATYPE = "genetic_association"


def resolve_target(data: dict) -> dict:
    """The target dict, or a clear error if the ensemblId didn't resolve."""
    target = data.get("target")
    if target is None:
        raise ValueError(
            "Open Targets returned no target for this id - check the "
            "Ensembl id resolved correctly before re-running"
        )
    return target


def identity_table(data: dict) -> pd.DataFrame:
    """One row: ensembl id, approved symbol/name, biotype, platform URL."""
    target = resolve_target(data)
    ensembl_id = target.get("id")
    return pd.DataFrame([{
        "ensembl_id": ensembl_id,
        "approved_symbol": target.get("approvedSymbol"),
        "approved_name": target.get("approvedName"),
        "biotype": target.get("biotype"),
        "open_targets_url": TARGET_URL.format(id=ensembl_id) if ensembl_id else None,
    }])


def known_drugs_table(data: dict) -> pd.DataFrame:
    """One row per candidate drug x indication: molecule, phase, mechanism, ids.

    A candidate lists several indications, so a drug appears once per
    disease; rows whose disease is null (unmapped indication) are kept
    with empty disease columns rather than dropped. mechanism_of_action
    is drug-level in Open Targets - it covers every target the molecule
    hits, not only this one - so multiple mechanisms are joined with
    "; " into a single cell.
    """
    columns = ["drug_id", "drug_name", "drug_type", "max_clinical_stage",
               "action_type", "mechanism_of_action",
               "disease_id", "disease_name", "drug_url"]
    rows = ((resolve_target(data).get("drugAndClinicalCandidates") or {}).get("rows") or [])
    out = []
    for row in rows:
        drug = row.get("drug") or {}
        moa_rows = ((drug.get("mechanismsOfAction") or {}).get("rows") or [])
        base = {
            "drug_id": drug.get("id"),
            "drug_name": drug.get("name"),
            "drug_type": drug.get("drugType"),
            "max_clinical_stage": row.get("maxClinicalStage"),
            "action_type": "; ".join(sorted({m["actionType"] for m in moa_rows
                                             if m.get("actionType")})) or None,
            "mechanism_of_action": "; ".join(m["mechanismOfAction"] for m in moa_rows
                                             if m.get("mechanismOfAction")) or None,
            "drug_url": DRUG_URL.format(id=drug["id"]) if drug.get("id") else None,
        }
        for entry in (row.get("diseases") or [{}]):
            disease = entry.get("disease") or {}
            out.append({**base,
                        "disease_id": disease.get("id"),
                        "disease_name": disease.get("name")})
    return pd.DataFrame(out, columns=columns)


def tractability_table(data: dict) -> pd.DataFrame:
    """One row per modality: label, modality code, value. Keeps False rows."""
    columns = ["label", "modality", "value"]
    rows = resolve_target(data).get("tractability") or []
    return pd.DataFrame(rows, columns=columns)


def association_base_rows(data: dict) -> list[dict]:
    """Shared base rows for associations_table and genetic_evidence_table."""
    rows = ((resolve_target(data).get("associatedDiseases") or {}).get("rows") or [])
    out = []
    for row in rows:
        disease = row.get("disease") or {}
        scores = {s["id"]: s.get("score") for s in (row.get("datatypeScores") or [])}
        out.append({
            "disease_id": disease.get("id"),
            "disease_name": disease.get("name"),
            "overall_score": row.get("score"),
            "disease_url": DISEASE_URL.format(id=disease["id"]) if disease.get("id") else None,
            "_datatype_scores": scores,
        })
    return out


def associations_table(data: dict) -> pd.DataFrame:
    """Up to 20 rows, one per disease, plus one column per datatype score
    (wide). Missing datatype for a disease renders as NaN, not 0.
    """
    base = association_base_rows(data)
    datatype_ids = sorted({dt for row in base for dt in row["_datatype_scores"]})
    out = []
    for row in base:
        entry = {k: v for k, v in row.items() if k != "_datatype_scores"}
        for dt in datatype_ids:
            entry[dt] = row["_datatype_scores"].get(dt)
        out.append(entry)
    columns = ["disease_id", "disease_name", "overall_score", "disease_url", *datatype_ids]
    return pd.DataFrame(out, columns=columns)


def genetic_evidence_table(data: dict) -> pd.DataFrame:
    """Disease id/name plus the genetic_association datatype score only."""
    columns = ["disease_id", "disease_name", "genetic_association_score", "disease_url"]
    out = [{
        "disease_id": row["disease_id"],
        "disease_name": row["disease_name"],
        "genetic_association_score": row["_datatype_scores"].get(GENETIC_DATATYPE),
        "disease_url": row["disease_url"],
    } for row in association_base_rows(data)]
    return pd.DataFrame(out, columns=columns)


def safety_table(data: dict) -> pd.DataFrame:
    """One row per reported safety liability. Often empty."""
    columns = ["event", "event_id", "datasource", "literature", "url"]
    rows = resolve_target(data).get("safetyLiabilities") or []
    out = [{
        "event": row.get("event"),
        "event_id": row.get("eventId"),
        "datasource": row.get("datasource"),
        "literature": row.get("literature"),
        "url": row.get("url"),
    } for row in rows]
    return pd.DataFrame(out, columns=columns)



def write_run(data: dict, out_dir: Path, subject: str, summary: str,
              steps: list[dict], source: str | None = None,
              limits: list[str] | None = None) -> dict[str, Path]:
    """Write one CSV per table into out_dir, then write README.md."""
    written = {}
    for name, table_fn in (
        ("identity", identity_table),
        ("known_drugs", known_drugs_table),
        ("tractability", tractability_table),
        ("associations", associations_table),
        ("genetic_evidence", genetic_evidence_table),
        ("safety", safety_table),
    ):
        path = out_dir / f"{name}.csv"
        table_fn(data).to_csv(path, index=False)
        written[name] = path
    written["readme"] = write_readme(out_dir / "README.md", subject, summary,
                                     steps, source=source, limits=limits)
    return written
