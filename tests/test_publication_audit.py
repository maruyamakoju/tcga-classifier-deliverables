"""Robustness tests for publication-readiness metadata checks."""
import audit_publication_readiness as publication


def test_json_metadata_check_does_not_crash_when_version_is_missing(tmp_path, monkeypatch):
    monkeypatch.setattr(publication, "ROOT", tmp_path)
    (tmp_path / ".zenodo.json").write_text("{}\n", encoding="utf-8")
    (tmp_path / "codemeta.json").write_text("{}\n", encoding="utf-8")
    messages = []

    publication.check_json_metadata(messages)

    assert messages == []


def test_release_artifacts_non_object_is_reported_without_opening_zip(tmp_path, monkeypatch):
    monkeypatch.setattr(publication, "ROOT", tmp_path)
    (tmp_path / "RELEASE_ARTIFACTS.json").write_text("[]\n", encoding="utf-8")
    (tmp_path / publication.ZIP_NAME).write_bytes(b"not-a-zip")
    messages = []

    publication.check_release_artifacts(messages)

    assert {message["code"] for message in messages} == {"release_artifacts_not_object"}


def test_corrupt_release_zip_is_a_structured_error(tmp_path, monkeypatch):
    monkeypatch.setattr(publication, "ROOT", tmp_path)
    (tmp_path / "RELEASE_ARTIFACTS.json").write_text("{}\n", encoding="utf-8")
    (tmp_path / publication.ZIP_NAME).write_bytes(b"not-a-zip")
    messages = []

    publication.check_release_artifacts(messages)

    assert {message["code"] for message in messages} == {"release_zip_invalid"}
