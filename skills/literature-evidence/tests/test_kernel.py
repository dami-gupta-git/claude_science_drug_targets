import pytest

import kernel

HIT_MED = {
    "id": "38012345",
    "source": "MED",
    "pmid": "38012345",
    "pmcid": "PMC1111111",
    "doi": "10.1000/example.001",
    "title": "WRN dependency in microsatellite-instable cancer",
    "authorString": "Smith J, Doe A.",
    "journalInfo": {"journal": {"title": "Nat Genet"}},
    "pubYear": "2023",
    "citedByCount": 42,
    "isOpenAccess": "Y",
    "abstractText": "WRN is selectively essential in MSI cancer cell lines.",
}

HIT_NO_IDS = {
    "title": "An untitled record",
    "abstractText": "No source or id given.",
}

HIT_PPR = {
    "id": "2024.01.01.500000",
    "source": "PPR",
    "doi": "10.1101/2024.01.01.500000",
    "title": "WRN synthetic lethality: a bioRxiv preprint",
    "authorString": "Lee K.",
    "pubYear": "2024",
    "citedByCount": 0,
    "isOpenAccess": "Y",
    "abstractText": "Preliminary evidence for WRN synthetic lethality in MSI models.",
}


def test_build_query_plain():
    q = kernel.build_query("WRN")
    assert 'TITLE:"WRN"' in q
    assert 'ABSTRACT:"WRN"' in q
    assert "AND" not in q


def test_build_query_with_extra_terms():
    q = kernel.build_query("WRN", extra_terms="biomarker")
    assert "AND (biomarker)" in q


def test_build_query_includes_preprints_by_default():
    q = kernel.build_query("WRN")
    assert "SRC:" not in q


def test_build_query_can_exclude_preprints():
    q = kernel.build_query("WRN", include_preprints=False)
    assert "NOT SRC:PPR" in q


def test_build_query_full_name_title_match_needs_symbol_in_abstract():
    q = kernel.build_query("TYK2", title_only=True, full_name="tyrosine kinase 2")
    assert q == '(TITLE:"TYK2" OR (TITLE:"tyrosine kinase 2" AND ABSTRACT:"TYK2"))'


def test_build_query_full_name_leaves_title_or_abstract_query_unchanged():
    assert (kernel.build_query("TYK2", full_name="tyrosine kinase 2")
            == kernel.build_query("TYK2"))


def test_build_query_year_window():
    q = kernel.build_query("TYK2", title_only=True, years=(2022, 2026))
    assert q.endswith("AND PUB_YEAR:[2022 TO 2026]")


def test_build_query_rejects_quote_in_full_name():
    with pytest.raises(ValueError, match="double quote"):
        kernel.build_query("TYK2", title_only=True, full_name='tyrosine "kinase" 2')


def test_search_balanced_runs_three_passes_and_joins_labels(monkeypatch):
    calls = []

    def fake_search(query, page_size, max_pages, timeout, sort):
        calls.append((query, sort))
        return [dict(HIT_MED)] if sort == kernel.SORT_CITED else []

    monkeypatch.setattr(kernel, "search_europepmc", fake_search)
    result = kernel.search_balanced("TYK2", "tyrosine kinase 2", established_years=5)
    assert [p["label"] for p in result["passes"]] == ["recent", "landmark", "established"]
    assert "PUB_YEAR:[" in calls[2][0]
    assert 'TITLE:"tyrosine kinase 2"' in calls[1][0]
    assert result["hits"][0]["retrieval_pass"] == "landmark; established"


def test_approved_name_requires_exact_symbol(monkeypatch):
    class Response:
        def __init__(self, body):
            self.body = body
        def __enter__(self):
            return self
        def __exit__(self, *args):
            return False
        def read(self):
            return self.body

    body = (b'{"data":{"search":{"hits":[{"object":{"approvedSymbol":"TYK2",'
            b'"approvedName":"tyrosine kinase 2"}}]}}}')
    monkeypatch.setattr(kernel.urllib.request, "urlopen", lambda *a, **k: Response(body))
    assert kernel.approved_name("tyk2") == "tyrosine kinase 2"
    with pytest.raises(ValueError, match="no target"):
        kernel.approved_name("TYK3")


def test_publication_row_maps_fields():
    row = kernel.publication_row(HIT_MED)
    assert row["pmid"] == "38012345"
    assert row["journal"] == "Nat Genet"
    assert row["citation_count"] == 42
    assert row["url"] == "https://europepmc.org/article/MED/38012345"


def test_publication_row_missing_source_or_id_gives_no_url():
    row = kernel.publication_row(HIT_NO_IDS)
    assert row["url"] is None
    assert row["pmid"] is None
    assert row["is_preprint"] is False


def test_publication_row_flags_preprint():
    row = kernel.publication_row(HIT_PPR)
    assert row["source"] == "PPR"
    assert row["is_preprint"] is True
    assert row["journal"] is None


def test_publications_table_columns_and_rows():
    df = kernel.publications_table([HIT_MED, HIT_NO_IDS])
    assert len(df) == 2
    assert list(df.columns) == kernel.PUBLICATION_COLUMNS


def test_publications_table_empty():
    df = kernel.publications_table([])
    assert len(df) == 0
    assert list(df.columns) == kernel.PUBLICATION_COLUMNS


VALID_CLAIM = {
    "claim_id": "C1", "gene": "WRN", "claim_type": "disease_association",
    "entity": "microsatellite-instable cancer", "direction": "dependency",
    "evidence_text": "WRN is selectively essential in MSI cancer cell lines.",
    "source_id": "MED/38012345",
}


def test_claims_table_valid_row():
    df = kernel.claims_table([VALID_CLAIM], [HIT_MED])
    assert len(df) == 1
    assert list(df.columns) == kernel.CLAIM_COLUMNS
    assert df.iloc[0]["claim_type"] == "disease_association"
    assert df.iloc[0]["pmid"] == "38012345"
    assert df.iloc[0]["url"] == "https://europepmc.org/article/MED/38012345"


def test_claims_table_empty():
    df = kernel.claims_table([], [])
    assert len(df) == 0
    assert list(df.columns) == kernel.CLAIM_COLUMNS


def test_claims_table_missing_field_raises():
    bad = {**VALID_CLAIM, "evidence_text": ""}
    with pytest.raises(ValueError, match="missing required field"):
        kernel.claims_table([bad], [HIT_MED])


def test_claims_table_bad_claim_type_raises():
    bad = {**VALID_CLAIM, "claim_type": "vibes"}
    with pytest.raises(ValueError, match="not in"):
        kernel.claims_table([bad], [HIT_MED])


def test_claims_table_bad_direction_raises():
    bad = {**VALID_CLAIM, "direction": "risk"}
    with pytest.raises(ValueError, match="direction"):
        kernel.claims_table([bad], [HIT_MED])


def test_claims_table_allows_empty_direction():
    df = kernel.claims_table([{**VALID_CLAIM, "direction": None}], [HIT_MED])
    assert len(df) == 1


def test_claims_table_unknown_source_raises():
    stray = {**VALID_CLAIM, "source_id": "MED/99999999"}
    with pytest.raises(ValueError, match="not in the retrieved"):
        kernel.claims_table([stray], [HIT_MED])


def test_claims_table_non_verbatim_evidence_raises():
    paraphrased = {**VALID_CLAIM,
                   "evidence_text": "WRN is essential in MSI cancers."}
    with pytest.raises(ValueError, match="verbatim"):
        kernel.claims_table([paraphrased], [HIT_MED])


SAMPLE_STEPS = [
    {"name": "Publications", "finding": "One publication found.", "table": [
        {"title": HIT_MED["title"], "pmid": HIT_MED["pmid"]},
    ]},
    {"name": "Claims", "finding": "One disease-association claim recorded."},
]


def test_write_run(tmp_path):
    written = kernel.write_run([HIT_MED], [VALID_CLAIM], tmp_path, "WRN",
                               "Summary.", SAMPLE_STEPS)
    assert set(written) == {"publications", "claims", "retrieval", "readme"}
    for path in written.values():
        assert path.is_file()
    readme = written["readme"].read_text()
    assert "# WRN — Europe PMC" in readme


def test_write_run_bad_claim_writes_nothing_after_it_but_publications_already_written(tmp_path):
    # publications.csv is written before claims are validated, matching the
    # single-writer/table-per-call pattern used elsewhere in this repo.
    bad_claim = {**VALID_CLAIM, "claim_type": "vibes"}
    with pytest.raises(ValueError):
        kernel.write_run([HIT_MED], [bad_claim], tmp_path, "WRN",
                         "Summary.", SAMPLE_STEPS)
    assert (tmp_path / "publications.csv").is_file()
    assert not (tmp_path / "claims.csv").is_file()
