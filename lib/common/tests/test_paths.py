from pathlib import Path

import pytest

from lib.common import paths


def test_project_root_from_env(tmp_path, monkeypatch):
    monkeypatch.setenv(paths.ROOT_ENV, str(tmp_path))
    assert paths.project_root() == tmp_path.resolve()


def test_project_root_falls_back_to_env_example(monkeypatch):
    monkeypatch.delenv(paths.ROOT_ENV, raising=False)
    expected = Path(__file__).resolve().parents[3]
    assert paths.project_root() == expected


def test_project_root_raises_for_missing_dir(tmp_path, monkeypatch):
    monkeypatch.setenv(paths.ROOT_ENV, str(tmp_path / "does_not_exist"))
    with pytest.raises(RuntimeError):
        paths.project_root()


def test_results_dir_creates_path(tmp_path, monkeypatch):
    monkeypatch.setenv(paths.ROOT_ENV, str(tmp_path))
    out = paths.results_dir("open-targets", "wrn")
    assert out.is_dir()
    assert out.parent == tmp_path / "results" / "open_targets" / "wrn"


def test_results_dir_rejects_bad_subject_name(tmp_path, monkeypatch):
    monkeypatch.setenv(paths.ROOT_ENV, str(tmp_path))
    with pytest.raises(ValueError):
        paths.results_dir("open-targets", "../escape")


def test_results_dir_rejects_bad_skill_name(tmp_path, monkeypatch):
    monkeypatch.setenv(paths.ROOT_ENV, str(tmp_path))
    with pytest.raises(ValueError):
        paths.results_dir("../escape", "wrn")
