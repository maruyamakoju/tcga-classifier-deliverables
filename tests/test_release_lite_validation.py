"""Unit tests for release-lite bundle validation."""
import json
import zipfile

import pytest

from validate_release_lite import (
    EXPECTED_BUNDLE_NAME,
    EXPECTED_MANIFEST_SCHEMA_VERSION,
    FORBIDDEN_NAMES,
    validate_manifest_metadata,
    validate_release_dir,
    validate_source_parity,
    validate_zip,
)
from validate_zip_bundle import validate_zip_members


def write_release_metadata(path, version="v-test", release_date="2026-07-09"):
    (path / "VERSION").write_text(version + "\n", encoding="utf-8")
    (path / "RELEASE_METADATA.json").write_text(
        json.dumps({"version": version, "release_date": release_date}) + "\n",
        encoding="utf-8",
    )


def base_manifest(file_count):
    return {
        "schema_version": EXPECTED_MANIFEST_SCHEMA_VERSION,
        "bundle_name": EXPECTED_BUNDLE_NAME,
        "version": "v-test",
        "release_date": "2026-07-09",
        "file_count_excluding_manifest_and_checksums": file_count,
        "forbidden_artifact_names": sorted(FORBIDDEN_NAMES),
        "files": [],
    }


def test_validate_manifest_metadata_accepts_current_contract(tmp_path):
    write_release_metadata(tmp_path)
    manifest_files = {"a.txt": {}, "b.txt": {}}

    errors = validate_manifest_metadata(tmp_path, base_manifest(2), manifest_files)

    assert errors == []


def test_validate_manifest_metadata_rejects_stale_top_level_fields(tmp_path):
    write_release_metadata(tmp_path)
    manifest = base_manifest(99)
    manifest.update(
        {
            "schema_version": "0.9",
            "bundle_name": "wrong-bundle",
            "version": "v-stale",
            "release_date": "2025-01-01",
            "forbidden_artifact_names": ["unexpected.pkl"],
        }
    )

    errors = "\n".join(validate_manifest_metadata(tmp_path, manifest, {"a.txt": {}}))

    assert "schema_version mismatch" in errors
    assert "bundle_name mismatch" in errors
    assert "file_count_excluding_manifest_and_checksums mismatch" in errors
    assert "forbidden_artifact_names mismatch" in errors
    assert "expected 'v-test' from VERSION" in errors
    assert "expected '2026-07-09' from RELEASE_METADATA.json" in errors


def test_validate_manifest_metadata_rejects_non_integer_file_count(tmp_path):
    write_release_metadata(tmp_path)
    manifest = base_manifest(1)
    manifest["file_count_excluding_manifest_and_checksums"] = True

    errors = validate_manifest_metadata(tmp_path, manifest, {"a.txt": {}})

    assert (
        "release_manifest.json file_count_excluding_manifest_and_checksums "
        "must be an integer"
    ) in errors


def test_validate_release_dir_rejects_non_object_manifest_without_crashing(tmp_path):
    (tmp_path / "release_manifest.json").write_text("[]\n", encoding="utf-8")

    errors, warnings, summary = validate_release_dir(tmp_path, max_file_bytes=5_000_000)

    assert "release_manifest.json top-level value must be an object" in errors
    assert warnings == []
    assert summary["has_manifest"] is True


def test_validate_release_dir_rejects_null_manifest_without_fallback(tmp_path):
    (tmp_path / "release_manifest.json").write_text("null\n", encoding="utf-8")

    errors, warnings, summary = validate_release_dir(tmp_path, max_file_bytes=5_000_000)

    assert "release_manifest.json top-level value must be an object" in errors
    assert warnings == []
    assert summary["has_manifest"] is True


def test_validate_release_dir_rejects_non_integer_manifest_bytes(tmp_path):
    payload = tmp_path / "payload.txt"
    payload.write_text("payload\n", encoding="utf-8")
    manifest = base_manifest(1)
    manifest["files"] = [
        {"path": "payload.txt", "bytes": "not-an-int", "sha256": "0" * 64}
    ]
    (tmp_path / "release_manifest.json").write_text(
        json.dumps(manifest) + "\n",
        encoding="utf-8",
    )

    errors, _, _ = validate_release_dir(tmp_path, max_file_bytes=5_000_000)

    assert "Manifest byte size for payload.txt must be an integer" in errors


def test_validate_source_parity_rejects_non_list_manifest_files(tmp_path):
    errors = validate_source_parity(tmp_path, tmp_path, {"files": {}})

    assert errors == ["release_manifest.json files must be a list"]


def test_validate_zip_rejects_duplicate_file_entries(tmp_path):
    release_dir = tmp_path / "release-lite"
    release_dir.mkdir()
    (release_dir / "payload.txt").write_text("payload\n", encoding="utf-8")
    zip_path = tmp_path / "bundle.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.writestr("payload.txt", "payload\n")
        with pytest.warns(UserWarning, match="Duplicate name"):
            zf.writestr("payload.txt", "payload\n")

    errors, summary = validate_zip(zip_path, release_dir)

    assert "Zip contains duplicate member path: payload.txt" in errors
    assert summary["zip_entries"] == 2


def test_validate_zip_bundle_members_rejects_duplicate_file_entries(tmp_path):
    zip_path = tmp_path / "bundle.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.writestr("payload.txt", "payload\n")
        with pytest.warns(UserWarning, match="Duplicate name"):
            zf.writestr("payload.txt", "payload\n")

    errors = validate_zip_members(zip_path)

    assert errors == ["Zip contains duplicate member path: payload.txt"]
