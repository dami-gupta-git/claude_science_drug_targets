"""Fetch and shape gene literature evidence from the Europe PMC REST API.

kernel.py loads with this skill: build_query, search_europepmc,
publication_row, publications_table, claims_table, write_run, and the
constants EUROPEPMC_SEARCH_URL, EUROPEPMC_ARTICLE_URL, PUBLICATION_COLUMNS,
CLAIM_TYPES, CLAIM_COLUMNS, PREPRINT_SOURCE.

Europe PMC indexes bioRxiv and medRxiv preprints under source "PPR"
alongside peer-reviewed literature, so no separate preprint source is
queried - build_query includes them by default.
"""
import json
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

from lib.common.reporting import write_readme

EUROPEPMC_SEARCH_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
EUROPEPMC_ARTICLE_URL = "https://europepmc.org/article/{source}/{id}"

PUBLICATION_COLUMNS = [
    "pmid", "pmcid", "doi", "title", "journal", "pub_year", "first_author",
    "citation_count", "is_open_access", "abstract", "url",
    "source", "is_preprint",
]

# Europe PMC's source code for bioRxiv/medRxiv (and other server) preprints.
PREPRINT_SOURCE = "PPR"

CLAIM_TYPES = ("disease_association", "interaction", "mechanism", "biomarker")

CLAIM_COLUMNS = [
    "claim_id", "gene", "claim_type", "entity", "direction",
    "evidence_text", "pmid", "url",
]


def build_query(gene: str, extra_terms: str | None = None,
                include_preprints: bool = True) -> str:
    """Europe PMC query for a gene symbol in title or abstract.

    extra_terms, if given, is AND-ed onto the gene term as free text
    (e.g. "drug target", "biomarker") to narrow a broad symbol. Europe PMC
    includes bioRxiv/medRxiv preprints (source PPR) by default; set
    include_preprints=False to restrict to peer-reviewed sources only.
    """
    if '"' in gene:
        raise ValueError('gene %r contains a double quote, which breaks '
                          'the Europe PMC query syntax' % gene)
    query = '(TITLE:"%s" OR ABSTRACT:"%s")' % (gene, gene)
    if extra_terms:
        query += " AND (%s)" % extra_terms
    if not include_preprints:
        query += " NOT SRC:%s" % PREPRINT_SOURCE
    return query


def search_europepmc(query: str, page_size: int = 25, max_pages: int = 4,
                     timeout: int = 60) -> list[dict]:
    """Raw Europe PMC result records for query, paginated via cursorMark.

    Stops after max_pages pages, or sooner if the API returns no further
    cursor or an empty page.
    """
    hits = []
    cursor = "*"
    for _ in range(max_pages):
        params = {
            "query": query, "format": "json", "resultType": "core",
            "pageSize": page_size, "cursorMark": cursor,
        }
        url = "%s?%s" % (EUROPEPMC_SEARCH_URL, urllib.parse.urlencode(params))
        with urllib.request.urlopen(url, timeout=timeout) as fh:
            data = json.load(fh)
        page_hits = (data.get("resultList") or {}).get("result") or []
        hits.extend(page_hits)
        next_cursor = data.get("nextCursorMark")
        if not next_cursor or next_cursor == cursor or not page_hits:
            break
        cursor = next_cursor
    return hits


def publication_row(hit: dict) -> dict:
    """One row from a raw Europe PMC result record."""
    source = hit.get("source")
    ident = hit.get("id")
    journal = ((hit.get("journalInfo") or {}).get("journal") or {}).get("title")
    return {
        "pmid": hit.get("pmid"),
        "pmcid": hit.get("pmcid"),
        "doi": hit.get("doi"),
        "title": hit.get("title"),
        "journal": journal,
        "pub_year": hit.get("pubYear"),
        "first_author": hit.get("authorString"),
        "citation_count": hit.get("citedByCount"),
        "is_open_access": hit.get("isOpenAccess"),
        "abstract": hit.get("abstractText"),
        "url": EUROPEPMC_ARTICLE_URL.format(source=source, id=ident)
               if source and ident else None,
        "source": source,
        "is_preprint": source == PREPRINT_SOURCE,
    }


def publications_table(hits: list[dict]) -> pd.DataFrame:
    """One row per publication, in PUBLICATION_COLUMNS order."""
    rows = [publication_row(hit) for hit in hits]
    return pd.DataFrame(rows, columns=PUBLICATION_COLUMNS)


def claims_table(claim_rows: list[dict]) -> pd.DataFrame:
    """One row per extracted literature claim.

    Every claim must carry a pmid and evidence_text (the exact sentence
    it was read from in that publication's abstract) - a claim without a
    locatable source sentence must not be recorded. claim_type must be one
    of CLAIM_TYPES.
    """
    for row in claim_rows:
        missing = [c for c in ("gene", "claim_type", "entity", "evidence_text", "pmid")
                  if not row.get(c)]
        if missing:
            raise ValueError("claim %r missing required field(s): %s" % (row, missing))
        if row["claim_type"] not in CLAIM_TYPES:
            raise ValueError(
                "claim_type %r not in %s" % (row["claim_type"], CLAIM_TYPES))
    return pd.DataFrame(claim_rows, columns=CLAIM_COLUMNS)


def write_run(hits: list[dict], claim_rows: list[dict], out_dir: Path,
              subject: str, summary: str, steps: list[dict],
              source: str = "Europe PMC",
              limits: list[str] | None = None) -> dict[str, Path]:
    """Write publications.csv, claims.csv and README.md into out_dir."""
    written = {}
    pubs_path = out_dir / "publications.csv"
    publications_table(hits).to_csv(pubs_path, index=False)
    written["publications"] = pubs_path

    claims_path = out_dir / "claims.csv"
    claims_table(claim_rows).to_csv(claims_path, index=False)
    written["claims"] = claims_path

    written["readme"] = write_readme(out_dir / "README.md", subject, summary,
                                     steps, source=source, limits=limits)
    return written
