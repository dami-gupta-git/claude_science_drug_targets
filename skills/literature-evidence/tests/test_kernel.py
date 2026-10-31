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
    "entity": "microsatellite-instable cancer", "direction": "positive",
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
