"""Record the search terms literature-evidence uses for a gene.

The terms are chosen by the model running this skill from two facts only,
the gene symbol and its approved name (approved_name). This file fetches
the name, checks the chosen terms, and writes them to search_terms.json,
which literature-evidence's search_balanced takes as input.
"""
import datetime
import json
import urllib.request
from pathlib import Path


OPEN_TARGETS_GRAPHQL_URL = "https://api.platform.opentargets.org/api/v4/graphql"

SEARCH_TERMS_FILE = "search_terms.json"

# More context words narrow a search past what the field publishes; one is
# a single point of failure for a symbol that collides with a word.
MAX_CONTEXT_TERMS = 4


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


def search_terms(gene: str, full_name: str | None, aliases: list[str],
                 context: list[str], reasons: dict[str, str]) -> dict:
    """The validated terms record. Raises on any problem rather than dropping a term.

    aliases and context may be empty lists, but every term in them needs a
    one-line entry in reasons, so a reader of search_terms.json can see why
    the search was widened or narrowed.
    """
    for label, terms in (("aliases", aliases), ("context", context)):
        if not isinstance(terms, list):
            raise ValueError("%s must be a list, got %r" % (label, terms))
        for term in terms:
            if not isinstance(term, str) or not term.strip():
                raise ValueError("%s has an empty or non-text term: %r" % (label, term))
            if '"' in term:
                raise ValueError("%s term %r contains a double quote, which breaks "
                                 "the Europe PMC query syntax" % (label, term))
    if len(context) > MAX_CONTEXT_TERMS:
        raise ValueError("context has %d terms; at most %d" % (len(context), MAX_CONTEXT_TERMS))
    if gene.upper() in (a.upper() for a in aliases):
        raise ValueError("aliases repeat the symbol %r" % gene)
    missing = [t for t in aliases + context if not (reasons.get(t) or "").strip()]
    if missing:
        raise ValueError("no reason given for term(s): %s" % missing)
    return {"gene": gene, "full_name": full_name, "aliases": aliases,
            "context": context, "reasons": {t: reasons[t] for t in aliases + context}}


def write_search_terms(terms: dict, out_dir: Path) -> Path:
    """Write a search_terms record (see search_terms) to out_dir/search_terms.json."""
    path = Path(out_dir) / SEARCH_TERMS_FILE
    stamp = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
    path.write_text(json.dumps({**terms, "written_at": stamp}, indent=2))
    return path


def latest_search_terms(gene: str, results_base: Path | None = None) -> dict:
    """The newest search_terms.json for gene under results/literature_search_terms/.

    Raises FileNotFoundError when none exists: literature-evidence must not
    search on the bare symbol by default.
    """
    if results_base is None:
        from lib.common.paths import project_root
        results_base = project_root() / "results"
    runs = sorted((Path(results_base) / "literature_search_terms" / gene.lower()).glob(
        "*/" + SEARCH_TERMS_FILE))
    if not runs:
        raise FileNotFoundError("no search terms for %s; run literature-search-terms first" % gene)
    return json.loads(runs[-1].read_text())
