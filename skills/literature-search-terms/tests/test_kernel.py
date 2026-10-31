import json

import pytest

import literature_search_terms_kernel as kernel


PELO = dict(gene="PELO", full_name="pelota mRNA surveillance and ribosome rescue factor",
            aliases=["pelota"], context=["ribosome", "mRNA"],
            reasons={"pelota": "older name used in titles",
                     "ribosome": "PELO is a Portuguese word",
                     "mRNA": "PELO is a Portuguese word"})


def test_search_terms_valid():
    terms = kernel.search_terms(**PELO)
    assert terms["aliases"] == ["pelota"]
    assert set(terms["reasons"]) == {"pelota", "ribosome", "mRNA"}


def test_search_terms_empty_lists_allowed():
    terms = kernel.search_terms("TYK2", "tyrosine kinase 2", [], [], {})
    assert terms["aliases"] == [] and terms["context"] == []


def test_search_terms_needs_a_reason_per_term():
    with pytest.raises(ValueError, match="no reason"):
        kernel.search_terms(**{**PELO, "reasons": {"pelota": "older name"}})


def test_search_terms_rejects_quote_and_symbol_repeat_and_too_many_context():
    with pytest.raises(ValueError, match="double quote"):
        kernel.search_terms(**{**PELO, "aliases": ['pe"lota']})
    with pytest.raises(ValueError, match="repeat the symbol"):
        kernel.search_terms(**{**PELO, "aliases": ["pelo"]})
    with pytest.raises(ValueError, match="at most"):
        kernel.search_terms(**{**PELO, "context": list("abcde"),
                               "reasons": {**PELO["reasons"], **{c: "x" for c in "abcde"}}})


def test_write_and_read_latest(tmp_path):
    base = tmp_path / "results"
    for stamp in ("2026-01-01_000000", "2026-02-01_000000"):
        run = base / "literature_search_terms" / "pelo" / stamp
        run.mkdir(parents=True)
        kernel.write_search_terms({**kernel.search_terms(**PELO), "stamp": stamp}, run)
    latest = kernel.latest_search_terms("PELO", results_base=base)
    assert latest["stamp"] == "2026-02-01_000000"
    assert "written_at" in latest


def test_latest_raises_when_missing(tmp_path):
    with pytest.raises(FileNotFoundError, match="run literature-search-terms"):
        kernel.latest_search_terms("PELO", results_base=tmp_path)


def test_approved_name_requires_exact_symbol(monkeypatch):
    class Response:
        def __enter__(self):
            return self
        def __exit__(self, *args):
            return False
        def read(self):
            return json.dumps({"data": {"search": {"hits": [{"object": {
                "approvedSymbol": "TYK2", "approvedName": "tyrosine kinase 2"}}]}}}).encode()

    monkeypatch.setattr(kernel.urllib.request, "urlopen", lambda *a, **k: Response())
    assert kernel.approved_name("tyk2") == "tyrosine kinase 2"
    with pytest.raises(ValueError, match="no target"):
        kernel.approved_name("TYK3")
