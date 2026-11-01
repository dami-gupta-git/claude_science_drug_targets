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


def terms(gene, full_name, aliases=(), context=(), ambiguous=()):
    """A search_terms record as literature-search-terms writes it."""
    return {"gene": gene, "full_name": full_name, "aliases": list(aliases),
            "context": list(context), "ambiguous": list(ambiguous), "reasons": {}}


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


def test_build_query_aliases_widen_every_name_clause():
    assert (kernel.build_query("TP53", aliases=["p53"])
            == '(TITLE:"TP53" OR ABSTRACT:"TP53" OR TITLE:"p53" OR ABSTRACT:"p53")')
    assert (kernel.build_query("TP53", title_only=True, aliases=["p53"])
            == '(TITLE:"TP53" OR TITLE:"p53")')


def test_build_query_context_applies_only_to_ambiguous_names():
    q = kernel.build_query("MAP3K13", title_only=True, aliases=["LZK"],
                           context=["kinase", "JNK"], ambiguous=["LZK"])
    assert q == ('(TITLE:"MAP3K13" OR ((TITLE:"LZK") AND '
                 '(TITLE:"kinase" OR ABSTRACT:"kinase" OR TITLE:"JNK" OR ABSTRACT:"JNK")))')


def test_build_query_context_and_ambiguous_go_together():
    with pytest.raises(ValueError, match="both"):
        kernel.build_query("PELO", context=["ribosome"])
    with pytest.raises(ValueError, match="both"):
        kernel.build_query("PELO", ambiguous=["PELO"])


def test_search_balanced_records_search_terms(monkeypatch):
    monkeypatch.setattr(kernel, "search_europepmc", lambda *a, **k: [])
    pelo = terms("PELO", "pelota factor", aliases=["pelota"], context=["ribosome"],
                 ambiguous=["PELO"])
    r = kernel.search_balanced("PELO", pelo)
    assert r["search_terms"] == pelo
    # no records, so no seeds and no snowball pass: every pass is a search
    assert all('ABSTRACT:"ribosome"' in p["query"] for p in r["passes"])
    assert all('TITLE:"pelota"' in p["query"] for p in r["passes"])
    assert "snowball" not in [p["label"] for p in r["passes"]]


def test_search_balanced_rejects_terms_for_another_gene(monkeypatch):
    monkeypatch.setattr(kernel, "search_europepmc", lambda *a, **k: [])
    with pytest.raises(ValueError, match="not 'TYK2'"):
        kernel.search_balanced("TYK2", terms("PELO", None))


def test_build_query_rejects_quote_in_full_name():
    with pytest.raises(ValueError, match="double quote"):
        kernel.build_query("TYK2", title_only=True, full_name='tyrosine "kinase" 2')


def test_search_balanced_runs_every_pass_and_joins_labels(monkeypatch):
    calls = []

    def fake_search(query, page_size, max_pages, timeout, sort):
        calls.append((query, sort))
        return [dict(HIT_MED)] if sort == kernel.SORT_CITED else []

    monkeypatch.setattr(kernel, "search_europepmc", fake_search)
    monkeypatch.setattr(kernel, "snowball",
                        lambda seeds, link, timeout=60: [dict(HIT_PPR)])
    result = kernel.search_balanced("TYK2", terms("TYK2", "tyrosine kinase 2"))
    assert [p["label"] for p in result["passes"]] == (
        ["relevant", "landmark", "recent_years"] + list(kernel.EVIDENCE_TOPICS)
        + ["referenced", "citing"])
    assert "PUB_YEAR:[" in calls[2][0]
    assert 'TITLE:"tyrosine kinase 2"' in calls[1][0]
    # every evidence topic is searched on every run, whatever the gene
    for topic_terms in kernel.EVIDENCE_TOPICS.values():
        assert any(topic_terms in q for q, _ in calls)
    by_id = {h["id"]: h for h in result["hits"]}
    # the fixture answers both citation-ranked passes, so its labels join
    assert by_id[HIT_MED["id"]]["retrieval_pass"] == "landmark; recent_years"
    assert by_id[HIT_PPR["id"]]["retrieval_pass"] == "referenced; citing"


def test_search_balanced_snowballs_from_the_best_cited_records(monkeypatch):
    low = {**HIT_MED, "id": "low", "citedByCount": 1}
    high = {**HIT_MED, "id": "high", "citedByCount": 900}
    seen = {}

    monkeypatch.setattr(kernel, "search_europepmc",
                        lambda *a, **k: [dict(low), dict(high)])

    def fake_snowball(seeds, link, timeout=60):
        seen["ids"] = [s["id"] for s in seeds]
        return []

    monkeypatch.setattr(kernel, "snowball", fake_snowball)
    kernel.search_balanced("TYK2", terms("TYK2", None), snowball_seeds=1)
    assert seen["ids"] == ["high"]


def test_snowball_refetches_linked_records_for_their_abstracts(monkeypatch):
    linked = {}

    def fake_linked(source, ident, link, limit=100, timeout=60):
        linked[link] = (source, ident)
        return [("MED", "111")] if link == kernel.LINK_REFERENCES else [("MED", "222")]

    def fake_by_key(keys, batch_size=25, timeout=60):
        return [{**HIT_MED, "id": i, "abstractText": "text"} for _, i in keys]

    monkeypatch.setattr(kernel, "fetch_linked_keys", fake_linked)
    monkeypatch.setattr(kernel, "fetch_records_by_key", fake_by_key)
    seed = [{"source": "MED", "id": "seed"}]
    refs = kernel.snowball(seed, kernel.LINK_REFERENCES)
    citing = kernel.snowball(seed, kernel.LINK_CITATIONS)
    assert linked == {kernel.LINK_REFERENCES: ("MED", "seed"),
                      kernel.LINK_CITATIONS: ("MED", "seed")}
    assert [h["id"] for h in refs] == ["111"]
    assert [h["id"] for h in citing] == ["222"]
    assert all(h.get("abstractText") for h in refs + citing)


def test_snowball_drops_the_seeds_themselves(monkeypatch):
    monkeypatch.setattr(kernel, "fetch_linked_keys",
                        lambda *a, **k: [("MED", "seed"), ("MED", "other")])
    monkeypatch.setattr(kernel, "fetch_records_by_key",
                        lambda keys, **k: [{"source": s, "id": i} for s, i in keys])
    out = kernel.snowball([{"source": "MED", "id": "seed"}], kernel.LINK_REFERENCES)
    assert [h["id"] for h in out] == ["other"]


def test_fetch_records_by_key_batches_per_source(monkeypatch):
    queries = []
    monkeypatch.setattr(kernel, "search_europepmc",
                        lambda query, **k: queries.append(query) or [])
    kernel.fetch_records_by_key([("MED", "1"), ("MED", "2"), ("PPR", "3")],
                                batch_size=2)
    assert len(queries) == 2
    assert 'EXT_ID:"1"' in queries[0] and "AND SRC:MED" in queries[0]
    assert 'EXT_ID:"3"' in queries[1] and "AND SRC:PPR" in queries[1]



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
