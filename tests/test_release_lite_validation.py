"""Unit tests for release-lite bundle validation."""
import json
import stat
import zipfile

import pytest

from validate_release_lite import (
    EXPECTED_BUNDLE_NAME,
    EXPECTED_MANIFEST_SCHEMA_VERSION,
    FORBIDDEN_NAMES,
    main as validate_release_lite_main,
    validate_manifest_metadata,
    validate_artifacts,
    validate_release_dir,
    validate_source_parity,
    validate_zip,
)
from validate_zip_bundle import (
    main as validate_zip_bundle_main,
    normalize_expected_sha,
    validate_zip_members,
)
from release_tools.common import (
    RELEASE_SCHEMA_VERSION,
    RELEASE_VALIDATION_COMMAND,
    RELEASE_ZIP_NAME,
    ZIP_ACCEPTANCE_COMMAND,
    canonical_zip_datetime,
    canonical_zip_info,
    sha256_file,
)


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
        "builder": "build_release_lite.py",
        "validation_command": RELEASE_VALIDATION_COMMAND,
        "files": [],
    }


def write_canonical_zip(zip_path, release_dir, release_date="2026-07-09"):
    timestamp = canonical_zip_datetime(release_date)
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_STORED) as zf:
        paths = sorted(
            (p for p in release_dir.rglob("*") if p.is_file()),
            key=lambda path: path.relative_to(release_dir).as_posix(),
        )
        for path in paths:
            rel = path.relative_to(release_dir).as_posix()
            zf.writestr(canonical_zip_info(rel, timestamp), path.read_bytes())


def artifact_contract(release_dir, zip_path, release_summary, zip_summary):
    digest = zip_summary["zip_sha256"]
    return {
        "schema_version": RELEASE_SCHEMA_VERSION,
        "version": "v-test",
        "release_date": "2026-07-09",
        "release_dir": release_dir.name,
        "release_file_count": release_summary["file_count"],
        "release_total_bytes": release_summary["total_bytes"],
        "zip_path": zip_path.name,
        "zip_entries": zip_summary["zip_entries"],
        "zip_bytes": zip_summary["zip_bytes"],
        "zip_sha256": digest,
        "validation_command": RELEASE_VALIDATION_COMMAND,
        "zip_acceptance_command": (
            f"{ZIP_ACCEPTANCE_COMMAND} --expected-sha256 {digest}"
        ),
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


def test_validate_release_dir_rejects_malformed_manifest_without_crashing(tmp_path):
    (tmp_path / "release_manifest.json").write_text("{not json\n", encoding="utf-8")

    errors, warnings, summary = validate_release_dir(tmp_path, max_file_bytes=5_000_000)

    assert any(error.startswith("Could not parse release_manifest.json:") for error in errors)
    assert warnings == []
    assert summary["has_manifest"] is True


def test_validate_release_lite_main_reports_malformed_manifest_without_crashing(tmp_path):
    release_dir = tmp_path / "release-lite"
    release_dir.mkdir()
    (release_dir / "release_manifest.json").write_text("{not json\n", encoding="utf-8")

    code = validate_release_lite_main([
        "--release-dir",
        str(release_dir),
        "--no-zip",
    ])

    assert code == 1


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


def test_validate_release_dir_rejects_empty_checksums(tmp_path):
    write_release_metadata(tmp_path)
    (tmp_path / "release_manifest.json").write_text(
        json.dumps(base_manifest(0)) + "\n", encoding="utf-8"
    )
    (tmp_path / "SHA256SUMS.txt").write_text("", encoding="ascii")

    errors, _, summary = validate_release_dir(tmp_path, max_file_bytes=5_000_000)

    assert any("at least one checksum" in error for error in errors)
    assert summary["checksum_count"] == 0


def test_validate_release_dir_reports_non_ascii_checksums_without_crashing(tmp_path):
    write_release_metadata(tmp_path)
    (tmp_path / "release_manifest.json").write_text(
        json.dumps(base_manifest(0)) + "\n", encoding="utf-8"
    )
    (tmp_path / "SHA256SUMS.txt").write_bytes(b"\xff\n")

    errors, _, _ = validate_release_dir(tmp_path, max_file_bytes=5_000_000)

    assert any("Could not parse SHA256SUMS.txt" in error for error in errors)


def test_validate_release_dir_rejects_incomplete_checksums(tmp_path):
    write_release_metadata(tmp_path)
    payload = tmp_path / "payload.txt"
    payload.write_text("payload\n", encoding="utf-8")
    manifest = base_manifest(1)
    manifest["files"] = [
        {"path": "payload.txt", "bytes": payload.stat().st_size, "sha256": sha256_file(payload)}
    ]
    manifest_path = tmp_path / "release_manifest.json"
    manifest_path.write_text(json.dumps(manifest) + "\n", encoding="utf-8")
    (tmp_path / "SHA256SUMS.txt").write_text(
        f"{sha256_file(manifest_path)}  release_manifest.json\n", encoding="ascii"
    )

    errors, _, _ = validate_release_dir(tmp_path, max_file_bytes=5_000_000)

    assert any("Files missing from SHA256SUMS" in error for error in errors)


def test_validate_source_parity_rejects_non_list_manifest_files(tmp_path):
    errors = validate_source_parity(tmp_path, tmp_path, {"files": {}})

    assert errors == ["release_manifest.json files must be a list"]


def test_validate_zip_accepts_canonical_archive(tmp_path):
    release_dir = tmp_path / "release-lite"
    release_dir.mkdir()
    write_release_metadata(release_dir)
    (release_dir / "payload.txt").write_text("payload\n", encoding="utf-8")
    zip_path = tmp_path / RELEASE_ZIP_NAME
    write_canonical_zip(zip_path, release_dir)

    errors, summary = validate_zip(zip_path, release_dir)

    assert errors == []
    assert summary["zip_entries"] == 3


def test_validate_zip_rejects_noncanonical_metadata_and_order(tmp_path):
    release_dir = tmp_path / "release-lite"
    release_dir.mkdir()
    write_release_metadata(release_dir)
    (release_dir / "payload.txt").write_text("payload\n", encoding="utf-8")
    zip_path = tmp_path / RELEASE_ZIP_NAME
    paths = sorted((p for p in release_dir.iterdir() if p.is_file()), reverse=True)
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.comment = b"noncanonical"
        for path in paths:
            info = zipfile.ZipInfo(path.name, date_time=(2001, 2, 3, 4, 5, 6))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.extra = b"\x01\x00\x00\x00"
            info.comment = b"member-comment"
            zf.writestr(info, path.read_bytes())

    errors, _ = validate_zip(zip_path, release_dir)
    joined = "\n".join(errors)

    assert "archive comment" in joined
    assert "member order" in joined
    assert "timestamp" in joined
    assert "compression is not canonical" in joined
    assert "extra field" in joined
    assert "member comment" in joined


def test_validate_zip_rejects_special_file_and_resource_limit(tmp_path):
    release_dir = tmp_path / "release-lite"
    release_dir.mkdir()
    write_release_metadata(release_dir)
    payload = release_dir / "payload.txt"
    payload.write_text("payload-too-large\n", encoding="utf-8")
    zip_path = tmp_path / RELEASE_ZIP_NAME
    timestamp = canonical_zip_datetime("2026-07-09")
    with zipfile.ZipFile(zip_path, "w") as zf:
        for path in sorted(release_dir.iterdir()):
            info = canonical_zip_info(path.name, timestamp)
            if path.name == "payload.txt":
                info.external_attr = (stat.S_IFLNK | 0o777) << 16
            zf.writestr(info, path.read_bytes())

    errors, _ = validate_zip(zip_path, release_dir, max_member_bytes=5)
    joined = "\n".join(errors)

    assert "special-file member" in joined
    assert "exceeds 5 uncompressed bytes" in joined


def test_validate_zip_rejects_entry_total_and_compression_ratio_limits(tmp_path):
    release_dir = tmp_path / "release-lite"
    release_dir.mkdir()
    write_release_metadata(release_dir)
    (release_dir / "payload.txt").write_text("x" * 10_000, encoding="utf-8")
    zip_path = tmp_path / RELEASE_ZIP_NAME
    timestamp = canonical_zip_datetime("2026-07-09")
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in sorted(release_dir.iterdir(), key=lambda item: item.name):
            info = canonical_zip_info(path.name, timestamp)
            info.compress_type = zipfile.ZIP_DEFLATED
            zf.writestr(info, path.read_bytes())

    errors, _ = validate_zip(
        zip_path,
        release_dir,
        max_entries=2,
        max_total_bytes=1_000,
        max_compression_ratio=2,
    )
    joined = "\n".join(errors)

    assert "entries; limit is 2" in joined
    assert "expands to" in joined
    assert "compression ratio exceeds 2" in joined


def test_validate_zip_rejects_raw_archive_above_limit_before_hashing(
    tmp_path, monkeypatch
):
    release_dir = tmp_path / "release-lite"
    release_dir.mkdir()
    write_release_metadata(release_dir)
    zip_path = tmp_path / RELEASE_ZIP_NAME
    write_canonical_zip(zip_path, release_dir)

    def forbidden_hash(_path):
        raise AssertionError("oversized archive must not be hashed")

    monkeypatch.setattr("validate_release_lite.sha256_file", forbidden_hash)
    errors, summary = validate_zip(
        zip_path, release_dir, max_archive_bytes=zip_path.stat().st_size - 1
    )

    assert "Zip archive is" in errors[0]
    assert summary["zip_entries"] == 0
    assert summary["zip_sha256"] is None


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


def test_validate_zip_reports_bad_zip_without_crashing(tmp_path):
    release_dir = tmp_path / "release-lite"
    release_dir.mkdir()
    zip_path = tmp_path / "bad.zip"
    zip_path.write_text("not a zip\n", encoding="utf-8")

    errors, summary = validate_zip(zip_path, release_dir)

    assert any(error.startswith("Could not read zip archive:") for error in errors)
    assert summary is None


def test_validate_artifacts_rejects_non_object_without_crashing(tmp_path):
    path = tmp_path / "RELEASE_ARTIFACTS.json"
    path.write_text("[]\n", encoding="utf-8")

    errors = validate_artifacts(path, {"unused": True}, {"unused": True})

    assert errors == ["RELEASE_ARTIFACTS.json top-level value must be an object"]


def test_validate_artifacts_enforces_complete_contract(tmp_path):
    release_dir = tmp_path / "release-lite"
    release_dir.mkdir()
    write_release_metadata(release_dir)
    zip_path = tmp_path / RELEASE_ZIP_NAME
    zip_path.write_bytes(b"zip-bytes")
    release_summary = {
        "release_dir": str(release_dir),
        "file_count": 2,
        "total_bytes": 42,
    }
    zip_summary = {
        "zip_path": str(zip_path),
        "zip_entries": 2,
        "zip_bytes": zip_path.stat().st_size,
        "zip_sha256": sha256_file(zip_path),
    }
    artifacts_path = tmp_path / "RELEASE_ARTIFACTS.json"
    artifact = artifact_contract(release_dir, zip_path, release_summary, zip_summary)
    artifacts_path.write_text(json.dumps(artifact) + "\n", encoding="utf-8")
    assert validate_artifacts(artifacts_path, release_summary, zip_summary) == []

    artifact.update(
        {
            "schema_version": "bogus",
            "version": "v-wrong",
            "release_date": "not-a-date",
            "release_dir": "wrong",
            "zip_path": "wrong.zip",
            "validation_command": "false",
            "zip_acceptance_command": "false",
        }
    )
    artifacts_path.write_text(json.dumps(artifact) + "\n", encoding="utf-8")

    joined = "\n".join(validate_artifacts(artifacts_path, release_summary, zip_summary))
    for field in [
        "schema_version",
        "version",
        "release_date",
        "release_dir",
        "zip_path",
        "validation_command",
        "zip_acceptance_command",
    ]:
        assert field in joined


def test_validate_zip_bundle_members_rejects_duplicate_file_entries(tmp_path):
    zip_path = tmp_path / "bundle.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.writestr("payload.txt", "payload\n")
        with pytest.warns(UserWarning, match="Duplicate name"):
            zf.writestr("payload.txt", "payload\n")

    errors = validate_zip_members(zip_path)

    assert "Zip contains duplicate member path: payload.txt" in errors


def test_validate_zip_bundle_members_reports_bad_zip_without_crashing(tmp_path):
    zip_path = tmp_path / "bad.zip"
    zip_path.write_text("not a zip\n", encoding="utf-8")

    errors = validate_zip_members(zip_path)

    assert any(error.startswith("Could not read zip archive:") for error in errors)


def test_validate_zip_bundle_members_accepts_canonical_archive(tmp_path):
    release_dir = tmp_path / "release-lite"
    release_dir.mkdir()
    write_release_metadata(release_dir)
    (release_dir / "payload.txt").write_text("payload\n", encoding="utf-8")
    zip_path = tmp_path / RELEASE_ZIP_NAME
    write_canonical_zip(zip_path, release_dir)

    assert validate_zip_members(zip_path) == []


def test_expected_sha_normalization_and_mismatch(tmp_path, capsys):
    assert normalize_expected_sha("sha256:" + "A" * 64) == "a" * 64
    with pytest.raises(ValueError, match="64-character"):
        normalize_expected_sha("bad")
    zip_path = tmp_path / "bad-digest.zip"
    zip_path.write_bytes(b"not the expected bytes")

    code = validate_zip_bundle_main(
        [str(zip_path), "--expected-sha256", "0" * 64, "--skip-acceptance"]
    )

    assert code == 1
    assert "Zip SHA256 mismatch" in capsys.readouterr().err


def test_validate_zip_bundle_main_reports_bad_zip_without_crashing(tmp_path):
    zip_path = tmp_path / "bad.zip"
    zip_path.write_text("not a zip\n", encoding="utf-8")

    code = validate_zip_bundle_main([str(zip_path), "--skip-acceptance"])

    assert code == 1
