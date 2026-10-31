"""Fetch and shape gene literature evidence from the Europe PMC REST API.

kernel.py loads with this skill: approved_name, build_query, search_europepmc,
search_balanced, publication_row, publications_table, publication_index,
normalise_text, verify_claims, claims_table, write_run, and the constants
EUROPEPMC_SEARCH_URL, EUROPEPMC_ARTICLE_URL, PUBLICATION_COLUMNS,
CLAIM_TYPES, CLAIM_COLUMNS, DIRECTIONS, PREPRINT_SOURCE, SORT_*.

Europe PMC indexes bioRxiv and medRxiv preprints under source "PPR"
alongside peer-reviewed literature, so no separate preprint source is
queried - build_query includes them by default.

Claims are keyed on source_id ("<source>/<id>", e.g. "MED/30971823"),
not on pmid. A preprint has no pmid, so keying on pmid made preprints
unrepresentable in the validated table and pushed them into a hand-made
file with no schema and no verbatim check - the least reviewed claims
got the least validation. source_id identifies the record; is_preprint,
derived from the source code, carries review status. write_run still
emits two files, but both are built and checked the same way.
"""
import datetime
import html
import json
import re
import time
import unicodedata
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd


EUROPEPMC_SEARCH_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
EUROPEPMC_ARTICLE_URL = "https://europepmc.org/article/{source}/{id}"

PUBLICATION_COLUMNS = [
    "source_id", "pmid", "pmcid", "doi", "title", "journal", "pub_year",
    "first_author", "authors", "pub_types", "citation_count",
    "is_open_access", "abstract", "url", "source", "is_preprint",
    "retrieval_pass",
]

# Europe PMC's source code for bioRxiv/medRxiv (and other server) preprints.
PREPRINT_SOURCE = "PPR"

CLAIM_TYPES = ("disease_association", "interaction", "mechanism", "biomarker")

# Europe PMC sort strings, confirmed against the live API. Omitting sort
# gives newest-first, NOT relevance - see search_europepmc. Only these
# values are accepted; anything else returns a response with no
# resultList rather than an error.
SORT_CITED = "CITED desc"
SORT_DATE_DESC = "P_PDATE_D desc"
SORT_DATE_ASC = "P_PDATE_D asc"
# Spelled out as literals rather than referencing the constants above:
# the sidecar gate allows only literal constants at module scope, and a
# tuple of names is not one. Keep in step with the SORT_* values.
SORT_OPTIONS = ("CITED desc", "P_PDATE_D desc", "P_PDATE_D asc",
                "CITED asc", "AUTH_FIRST asc")

# claim_id, gene, claim_type, entity, direction and evidence_text are
# written by hand. source_id is written by hand. Everything after it is
# filled in by claims_table from the matching publication - a hand-copied
# pmid or url is a transcription error waiting to happen.
CLAIM_COLUMNS = [
    "claim_id", "gene", "claim_type", "entity", "direction",
    "evidence_text", "source_id", "pmid", "doi", "is_preprint", "url",
]

CLAIM_REQUIRED = ("gene", "claim_type", "entity", "evidence_text", "source_id")

# The only accepted `direction` values, each read as a statement about the
# gene. A free-text field let a paper showing that blocking TYK2 rescues
# Alzheimer's-type toxicity be tagged "risk", the opposite reading.
# direction may be left empty where none of these applies.
DIRECTIONS = {
    "increases risk": "the gene, its activity or a variant in it raises risk or severity of the entity",
    "decreases risk": "the gene, its activity or a variant in it lowers risk of the entity",
    "loss causes": "losing the gene causes the entity (e.g. an inherited deficiency)",
    "therapeutic target": "inhibiting or modulating the gene treats, or is proposed to treat, the entity",
    "dependency": "cells of the entity need the gene to survive or grow",
    "suppresses": "the gene restrains the entity; losing or blocking it promotes the entity",
    "required for": "the gene is needed for the entity (a pathway, process or response)",
    "activates": "the gene activates or raises the entity",
    "activated by": "the entity activates the gene",
    "inhibited by": "the entity inhibits, dephosphorylates or blocks the gene",
    "binds": "the entity physically binds the gene, with no sign of effect stated",
    "elevated": "the gene's expression or activity is raised in the entity",
    "reduced": "the gene's expression or activity is lowered in, or by, the entity",
}

OPEN_TARGETS_GRAPHQL_URL = "https://api.platform.opentargets.org/api/v4/graphql"


def approved_name(symbol: str, timeout: int = 60) -> str:
    """The gene's approved full name from Open Targets (TYK2 -> "tyrosine kinase 2").

    Raises when no target's approved symbol matches exactly, rather than
    taking the top search hit, which for a mistyped symbol is another gene.
    """
    query = ('query($q:String!){search(queryString:$q,entityNames:["target"],'
             'page:{index:0,size:10}){hits{object{... on Target{approvedSymbol approvedName}}}}}')
    request = urllib.request.Request(
        OPEN_TARGETS_GRAPHQL_URL,
        data=json.dumps({"query": query, "variables": {"q": symbol}}).encode(),
        headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(request, timeout=timeout) as fh:
        body = json.load(fh)
    hits = (((body.get("data") or {}).get("search") or {}).get("hits")) or []
    for hit in hits:
        target = hit.get("object") or {}
        if (target.get("approvedSymbol") or "").upper() == symbol.upper():
            return target["approvedName"]
    raise ValueError("Open Targets has no target with approved symbol %r" % symbol)


def build_query(gene: str, extra_terms: str | None = None,
                include_preprints: bool = True,
                title_only: bool = False,
                full_name: str | None = None,
                years: tuple[int, int] | None = None,
                aliases: list[str] | None = None,
                context: list[str] | None = None) -> str:
    """Europe PMC query for a gene symbol in title or abstract.

    aliases are other names the literature uses, matched wherever the
    symbol is (TP53 papers mostly say p53; MAP3K13 papers say LZK). context
    words, when given, must appear in the abstract: they narrow a symbol
    that is also an ordinary word (PELO is Portuguese for "by the"). Both
    come from the "Choose search terms" step in SKILL.md.

    extra_terms, if given, is AND-ed onto the gene term as free text
    (e.g. "drug target", "biomarker") to narrow a broad symbol. Europe PMC
    includes bioRxiv/medRxiv preprints (source PPR) by default; set
    include_preprints=False to restrict to peer-reviewed sources only.

    title_only drops the ABSTRACT clause, so only papers whose title
    carries the symbol match. Combined with citation sorting this is what
    surfaces a field's founding papers; used alone it loses current work
    that discusses a gene without naming it in the title. search_balanced
    applies it to the landmark and established passes.

    full_name widens a title_only query to titles that spell the gene out
    (TYK2's most-cited paper is titled "Human tyrosine kinase 2
    deficiency"), but only when the abstract also names the symbol: the
    bare phrase "tyrosine kinase 2" also titles papers on proline-rich
    tyrosine kinase 2 (PYK2). The title-or-abstract query already reaches
    such papers through the abstract, so full_name does not change it.

    years, an inclusive (first, last) pair, restricts by publication year.
    """
    names = [gene, *(aliases or [])]
    for term in [*names, full_name, *(context or [])]:
        if term and '"' in term:
            raise ValueError('%r contains a double quote, which breaks '
                             'the Europe PMC query syntax' % term)
    if title_only:
        clauses = ['TITLE:"%s"' % n for n in names]
        if full_name:
            clauses.append('(TITLE:"%s" AND ABSTRACT:"%s")' % (full_name, gene))
    else:
        clauses = ['TITLE:"%s" OR ABSTRACT:"%s"' % (n, n) for n in names]
    query = "(%s)" % " OR ".join(clauses)
    if context:
        query += " AND (%s)" % " OR ".join('ABSTRACT:"%s"' % c for c in context)
    if years:
        query += " AND PUB_YEAR:[%d TO %d]" % years
    if extra_terms:
        query += " AND (%s)" % extra_terms
    if not include_preprints:
        query += " NOT SRC:%s" % PREPRINT_SOURCE
    return query


def fetch_json_retrying(url: str, timeout: int, attempts: int = 3,
                backoff: float = 2.0) -> dict:
    """GET url and parse JSON, retrying transient network failures.

    A dossier makes several calls; without this, one timed-out page loses
    the whole run. HTTP 4xx is not retried - a malformed query will fail
    the same way three times.
    """
    last = None
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as fh:
                return json.load(fh)
        except urllib.error.HTTPError as exc:
            if exc.code < 500:
                raise
            last = exc
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            last = exc
        if attempt < attempts - 1:
            time.sleep(backoff * (2 ** attempt))
    raise RuntimeError("Europe PMC request failed after %d attempts: %s"
                       % (attempts, last))


def search_europepmc(query: str, page_size: int = 25, max_pages: int = 4,
                     timeout: int = 60, sort: str | None = None) -> list[dict]:
    """Raw Europe PMC result records for query, paginated via cursorMark.

    Stops after max_pages pages, or sooner if the API returns no further
    cursor or an empty page.

    sort selects the API's ordering. Passing None uses the API default,
    which is NEWEST FIRST, not relevance - on a well-studied gene that
    returns only the last year or two and silently omits the founding
    papers. SORT_CITED reaches those. Prefer search_balanced() over
    calling this directly for a dossier; see its docstring.

    An invalid sort string is not an error at the API: the response comes
    back with no resultList, which would otherwise read here as an empty
    result set. When a sorted first page comes back empty, the same query
    is retried unsorted; if that has hits, the sort string is at fault and
    this raises rather than returning [].
    """
    hits = []
    cursor = "*"
    for page in range(max_pages):
        params = {
            "query": query, "format": "json", "resultType": "core",
            "pageSize": page_size, "cursorMark": cursor,
        }
        if sort:
            params["sort"] = sort
        url = "%s?%s" % (EUROPEPMC_SEARCH_URL, urllib.parse.urlencode(params))
        data = fetch_json_retrying(url, timeout)
        page_hits = (data.get("resultList") or {}).get("result") or []
        if page == 0 and sort and not page_hits:
            raise_if_sort_at_fault(query, sort, page_size, timeout)
        hits.extend(page_hits)
        next_cursor = data.get("nextCursorMark")
        if not next_cursor or next_cursor == cursor or not page_hits:
            break
        cursor = next_cursor
    return hits


def raise_if_sort_at_fault(query: str, sort: str, page_size: int,
                            timeout: int) -> None:
    """Raise if query has hits unsorted but none with this sort string.

    Distinguishes a rejected sort value, which the API reports as a
    response with no resultList, from a query that genuinely matches
    nothing. Costs one extra call, and only on the failing path.
    """
    params = {
        "query": query, "format": "json", "resultType": "idlist",
        "pageSize": page_size, "cursorMark": "*",
    }
    url = "%s?%s" % (EUROPEPMC_SEARCH_URL, urllib.parse.urlencode(params))
    unsorted = (fetch_json_retrying(url, timeout).get("resultList") or {}).get("result") or []
    if unsorted:
        raise ValueError(
            "Europe PMC returned nothing for sort=%r but %d hit(s) for the "
            "same query unsorted - the sort string is not accepted. Use one "
            "of: %s" % (sort, len(unsorted),
                        ", ".join(repr(s) for s in SORT_OPTIONS))
        )


# A retrieval pass and a whole retrieval are plain dicts, not NamedTuples:
# the sidecar gate allows no class definitions at module scope. Dict keys
# are the former field names, so retrieval.json is unchanged.
#
# pass keys:      label, query, sort, requested, returned, year_min, year_max
# retrieval keys: hits, passes, retrieved_at


def year_span(hits: list[dict]) -> tuple[int | None, int | None]:
    """Earliest and latest pubYear among hits, or (None, None) if none parse."""
    years = []
    for hit in hits:
        try:
            years.append(int(hit.get("pubYear")))
        except (TypeError, ValueError):
            continue
    return (min(years), max(years)) if years else (None, None)


def pass_info(label, query, sort, requested, returned, year_min, year_max,
              reranked_from=None):
    """One retrieval pass: what it asked for and what came back.

    reranked_from is the candidate pool size when the pass re-ranked the
    API's top records by citations per year; None when the API order stood.
    """
    return {"label": label, "query": query, "sort": sort,
            "reranked_from": reranked_from,
            "requested": requested, "returned": returned,
            "year_min": year_min, "year_max": year_max}


def citations_per_year(hit: dict, this_year: int) -> float | None:
    """citedByCount over years since publication, counting the publication
    year as one; None when either value is missing, so the record ranks last
    rather than as zero."""
    try:
        cited, year = int(hit.get("citedByCount")), int(hit.get("pubYear"))
    except (TypeError, ValueError):
        return None
    return cited / max(this_year - year + 1, 1)


def retrieval_record(hits, passes, retrieved_at, search_terms=None):
    """Records plus the provenance of how they were retrieved.

    search_balanced returns this so write_run can write retrieval.json
    without the caller having to remember what it asked for. Access it by
    key:

        result = search_balanced("WRN", approved_name("WRN"))
        hits = result["hits"]

    Note it is not a 2-tuple: an earlier version of this docstring showed
    `hits, retrieval = search_balanced(...)`, which raised ValueError
    because the record carries three fields, not two.
    """
    return {"hits": hits, "passes": passes, "retrieved_at": retrieved_at,
            "search_terms": search_terms}


def search_balanced(gene: str, full_name: str | None,
                    aliases: list[str] | None = None,
                    context: list[str] | None = None,
                    extra_terms: str | None = None,
                    recent: int = 20, landmark: int = 15,
                    established: int = 15, established_years: int = 5,
                    established_pool: int = 100,
                    include_preprints: bool = True,
                    timeout: int = 60) -> dict:
    """Records from three passes - newest first, most cited, and most cited
    per year in the last established_years years - merged.

    full_name is required (approved_name gives it) and may be None only for
    a gene with no full name; omitting it loses papers whose title spells
    the gene out. aliases and context apply to every pass and are recorded
    in retrieval.json; see build_query.

    Returns a dict with keys hits, passes and retrieved_at (see
    retrieval_record); hand the whole dict to write_run, which writes it to
    retrieval.json. Recording the query,
    sort and year span alongside the records is what lets a later reader
    see what the retrieval could not have reached; describing it in prose
    means it is only there when someone remembers to write it.

    A single pass cannot cover a well-studied gene, for three separate
    reasons that have to be fixed separately.

    The API default order is newest first, NOT relevance. On WRN that
    returned only 2025-2026 records, so the 2019 papers establishing the
    target were absent whether 12 or 30 records were requested; widening
    the cap cannot reach them.

    Sorting by citation count reaches old work but not, on its own, the
    right old work. On the title-or-abstract query the founding WRN paper
    ranks 15th by citations, behind papers on SIRT6, RECQL4 and yeast
    Srs2 that merely mention the gene. Half the top 30 are not about it.
    So the landmark pass narrows the gene term to the TITLE, where the
    same paper ranks 2nd. Narrowing is right for this pass and wrong for
    the recent pass, which needs abstract matches to see current work
    that does not lead with the gene name.

    The two passes still leave a hole. Newest-first stops at the last few
    months, and all-time citation counts favour papers that have had
    decades to accrue them, so on TYK2 nothing from 2020-2025 was
    retrieved: the deucravacitinib trials and approval were absent. The
    established pass searches the last established_years years on the
    landmark pass's title query. Ranked by total citations it still took
    only the window's oldest years (TYK2: 2022-2023, nothing from
    2024-2025), so it fetches the top established_pool records by citations
    and keeps the established with the most citations per year (see
    citations_per_year). Europe PMC has no per-year sort.

    Records are deduplicated on (source, id). Each carries a
    `retrieval_pass` key naming every pass that returned it, joined with
    "; " (e.g. "landmark; established"). A record in more than one pass is
    current and well cited, which is worth knowing when weighting claims.

    recent, landmark and established are record counts, not pages; 0 skips
    a pass.
    """
    out = {}
    infos = []
    this_year = datetime.date.today().year
    window = (this_year - established_years + 1, this_year)
    terms = {"aliases": aliases, "context": context}
    passes = (
        ("recent", build_query(gene, extra_terms, include_preprints, **terms),
         None, recent, None),
        ("landmark", build_query(gene, extra_terms, include_preprints,
                                 title_only=True, full_name=full_name, **terms),
         SORT_CITED, landmark, None),
        ("established", build_query(gene, extra_terms, include_preprints,
                                    title_only=True, full_name=full_name,
                                    years=window, **terms),
         SORT_CITED, established, established_pool),
    )
    for label, query, sort_key, n, pool in passes:
        if n <= 0:
            continue
        found = search_europepmc(query, page_size=max(n, pool or 0), max_pages=1,
                                 timeout=timeout, sort=sort_key)
        if pool:
            rate = lambda h: citations_per_year(h, this_year)
            found = sorted(found, key=lambda h: (rate(h) is not None, rate(h) or 0),
                           reverse=True)[:n]
        year_min, year_max = year_span(found)
        infos.append(pass_info(label, query, sort_key, n, len(found),
                              year_min, year_max, reranked_from=pool))
        for hit in found:
            key = (hit.get("source"), hit.get("id"))
            if key in out:
                out[key]["retrieval_pass"] += "; " + label
            else:
                hit["retrieval_pass"] = label
                out[key] = hit
    stamp = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
    return retrieval_record(list(out.values()), infos, stamp,
                            {"full_name": full_name, **terms})


# --- text normalisation and verbatim checking -------------------------

# Patterns are kept as literal strings rather than compiled objects: the
# skill sidecar gate only allows literal constants at module scope, and
# re caches compiled patterns internally so re.sub on a string costs the
# same after the first call.
TAG_PATTERN = r"<[^>]+>"
WS_PATTERN = r"\s+"
# Abstracts arrive with typographic punctuation; a sentence copied out of
# one and back in can pick up ASCII equivalents. Fold both to ASCII before
# comparing so a verbatim quote does not fail on a dash.
PUNCT_MAP = {
    0x2010: "-", 0x2011: "-", 0x2012: "-", 0x2013: "-", 0x2014: "-",
    0x2018: "'", 0x2019: "'", 0x201c: '"', 0x201d: '"', 0x00a0: " ",
    0x2032: "'", 0x2033: '"',
}


def normalise_text(text: str | None) -> str:
    """Strip HTML markup and fold punctuation and whitespace.

    Europe PMC wraps gene symbols in italic tags inside abstractText, so
    a genuinely verbatim sentence fails a naive substring test. This is
    the failure the check has to survive, not work around.
    """
    if not text:
        return ""
    text = re.sub(TAG_PATTERN, "", text)
    text = html.unescape(text)
    text = unicodedata.normalize("NFKC", text)
    text = text.translate(PUNCT_MAP)
    return re.sub(WS_PATTERN, " ", text).strip()


def source_id(hit: dict) -> str | None:
    """Stable "<source>/<id>" key for a raw Europe PMC record."""
    source, ident = hit.get("source"), hit.get("id")
    return "%s/%s" % (source, ident) if source and ident else None


def publication_index(hits: list[dict]) -> dict[str, dict]:
    """source_id -> publication row, for claim validation and enrichment."""
    index = {}
    for hit in hits:
        key = source_id(hit)
        if key:
            index[key] = publication_row(hit)
    return index


def verify_claims(claim_rows: list[dict], hits: list[dict]) -> list[str]:
    """Problems found checking each claim against its publication.

    Returns a list of human-readable problems; empty means every claim
    names a retrieved publication and its evidence_text appears verbatim
    in that publication's abstract. claims_table raises on a non-empty
    result, so a run cannot write a claims table it has not checked.
    """
    index = publication_index(hits)
    problems = []
    for row in claim_rows:
        cid = row.get("claim_id") or "<no claim_id>"
        key = row.get("source_id")
        if not key:
            problems.append("%s: no source_id" % cid)
            continue
        pub = index.get(key)
        if pub is None:
            problems.append("%s: source_id %r is not in the retrieved "
                            "publications" % (cid, key))
            continue
        abstract = normalise_text(pub.get("abstract"))
        if not abstract:
            problems.append("%s: publication %s has no abstract to check "
                            "the quote against" % (cid, key))
            continue
        if normalise_text(row.get("evidence_text")) not in abstract:
            problems.append("%s: evidence_text is not a verbatim substring "
                            "of the abstract of %s" % (cid, key))
    return problems


def publication_row(hit: dict) -> dict:
    """One row from a raw Europe PMC result record."""
    source = hit.get("source")
    ident = hit.get("id")
    journal = ((hit.get("journalInfo") or {}).get("journal") or {}).get("title")
    authors = hit.get("authorString") or ""
    pub_types = (hit.get("pubTypeList") or {}).get("pubType") or []
    if isinstance(pub_types, str):
        pub_types = [pub_types]
    return {
        "source_id": source_id(hit),
        "pmid": hit.get("pmid"),
        "pmcid": hit.get("pmcid"),
        "doi": hit.get("doi"),
        "title": hit.get("title"),
        "journal": journal,
        "pub_year": hit.get("pubYear"),
        # authorString is the whole author list; first_author is the first
        # name in it. The old column held the full list under the singular
        # name, which broke anything joining on it.
        "first_author": authors.split(",")[0].strip().rstrip(".") or None,
        "authors": authors or None,
        # "review" vs "research-article" vs preprint distinguishes a paper
        # that established a finding from one restating it.
        "pub_types": "; ".join(pub_types) or None,
        "citation_count": hit.get("citedByCount"),
        "is_open_access": hit.get("isOpenAccess"),
        "abstract": hit.get("abstractText"),
        "url": EUROPEPMC_ARTICLE_URL.format(source=source, id=ident)
               if source and ident else None,
        "source": source,
        "is_preprint": source == PREPRINT_SOURCE,
        # set by search_balanced; absent when search_europepmc was called
        # directly, in which case the pass is unknown rather than missing
        "retrieval_pass": hit.get("retrieval_pass"),
    }


def publications_table(hits: list[dict]) -> pd.DataFrame:
    """One row per publication, in PUBLICATION_COLUMNS order."""
    rows = [publication_row(hit) for hit in hits]
    return pd.DataFrame(rows, columns=PUBLICATION_COLUMNS)


def duplicate_titles(hits: list[dict]) -> list[list[str]]:
    """Groups of source_ids sharing a normalised title.

    Deduplication is on (source, id), so a preprint and its published
    version are two records. Claims are counted per publication on the
    principle that two publications agreeing is stronger evidence than
    one; a preprint and its own journal version would inflate that count.
    """
    seen: dict[str, list[str]] = {}
    for hit in hits:
        title = normalise_text(hit.get("title")).lower().rstrip(".")
        key = source_id(hit)
        if title and key:
            seen.setdefault(title, []).append(key)
    return [ids for ids in seen.values() if len(ids) > 1]


def claims_table(claim_rows: list[dict], hits: list[dict]) -> pd.DataFrame:
    """One row per extracted literature claim, checked against its source.

    Every claim must carry a source_id naming a publication in hits, and
    an evidence_text that appears verbatim in that publication's abstract
    - a claim without a locatable source sentence must not be recorded.
    claim_type must be one of CLAIM_TYPES. pmid, doi, url and is_preprint
    are filled in from the publication rather than copied by hand.

    Raises on any failure. A claims table that has not passed this check
    should not exist on disk.
    """
    for row in claim_rows:
        missing = [c for c in CLAIM_REQUIRED if not row.get(c)]
        if missing:
            raise ValueError("claim %r missing required field(s): %s"
                             % (row.get("claim_id") or row, missing))
        if row["claim_type"] not in CLAIM_TYPES:
            raise ValueError(
                "claim_type %r not in %s" % (row["claim_type"], CLAIM_TYPES))
        if row.get("direction") and row["direction"] not in DIRECTIONS:
            raise ValueError("direction %r not in DIRECTIONS: %s"
                             % (row["direction"], sorted(DIRECTIONS)))

    problems = verify_claims(claim_rows, hits)
    if problems:
        raise ValueError("%d claim(s) failed verification:\n  %s"
                         % (len(problems), "\n  ".join(problems)))

    index = publication_index(hits)
    enriched = []
    for row in claim_rows:
        pub = index[row["source_id"]]
        enriched.append({**row,
                         "pmid": pub.get("pmid"),
                         "doi": pub.get("doi"),
                         "url": pub.get("url"),
                         "is_preprint": pub.get("is_preprint")})
    return pd.DataFrame(enriched, columns=CLAIM_COLUMNS)


def retrieval_manifest(hits: list[dict], retrieval: dict | None,
                       gene: str) -> dict:
    """What was asked of Europe PMC and what came back."""
    year_min, year_max = year_span(hits)
    preprints = sum(1 for h in hits if h.get("source") == PREPRINT_SOURCE)
    by_pass: dict[str, int] = {}
    for hit in hits:
        label = hit.get("retrieval_pass") or "unknown"
        by_pass[label] = by_pass.get(label, 0) + 1
    return {
        "gene": gene,
        "retrieved_at": retrieval["retrieved_at"] if retrieval else None,
        "search_terms": retrieval.get("search_terms") if retrieval else None,
        "passes": list(retrieval["passes"]) if retrieval else [],
        "records": {
            "total": len(hits),
            "year_min": year_min,
            "year_max": year_max,
            "preprints": preprints,
            "by_retrieval_pass": by_pass,
        },
        "duplicate_titles": duplicate_titles(hits),
    }


def write_run(hits: list[dict], claim_rows: list[dict], out_dir: Path,
              subject: str, summary: str, steps: list[dict],
              source: str = "Europe PMC",
              limits: list[str] | None = None,
              retrieval: dict | None = None) -> dict[str, Path]:
    """Write publications.csv, claims.csv, retrieval.json and README.md.

    Claims from preprints go to claims_preprint.csv. Both files are built
    by claims_table and so carry the same schema check and the same
    verbatim check; the split marks review status, it is not a fallback
    for rows that could not be validated.
    """
    from lib.common.reporting import write_readme
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written = {}

    pubs_path = out_dir / "publications.csv"
    publications_table(hits).to_csv(pubs_path, index=False)
    written["publications"] = pubs_path

    # Built and checked in one call, then split - so a preprint claim
    # cannot skip validation by virtue of being a preprint.
    claims = claims_table(claim_rows, hits)
    peer_reviewed = claims[~claims["is_preprint"].astype(bool)]
    preprint = claims[claims["is_preprint"].astype(bool)]

    claims_path = out_dir / "claims.csv"
    peer_reviewed.to_csv(claims_path, index=False)
    written["claims"] = claims_path

    if not preprint.empty:
        preprint_path = out_dir / "claims_preprint.csv"
        preprint.to_csv(preprint_path, index=False)
        written["claims_preprint"] = preprint_path

    manifest_path = out_dir / "retrieval.json"
    manifest_path.write_text(
        json.dumps(retrieval_manifest(hits, retrieval, subject), indent=2))
    written["retrieval"] = manifest_path

    written["readme"] = write_readme(out_dir / "README.md", subject, summary,
                                     steps, source=source, limits=limits)
    return written
