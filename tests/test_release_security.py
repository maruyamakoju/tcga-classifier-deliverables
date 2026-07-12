"""Security regressions for release path and trusted-ZIP handling."""
import json
import zipfile

import pytest

import run_release_acceptance as acceptance
import validate_zip_bundle as zip_bundle
from release_tools.common import (
    RELEASE_SCHEMA_VERSION,
    RELEASE_ZIP_NAME,
    ZIP_ACCEPTANCE_COMMAND,
    canonical_zip_datetime,
    canonical_zip_info,
    normalize_release_path,
)


def write_structurally_canonical_zip(path):
    timestamp = canonical_zip_datetime("2026-07-09")
    metadata = json.dumps(
        {"version": "v-test", "release_date": "2026-07-09"}, sort_keys=True
    ).encode("utf-8")
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_STORED) as archive:
        archive.writestr(
            canonical_zip_info("RELEASE_METADATA.json", timestamp), metadata
        )
        archive.writestr(
            canonical_zip_info("validate_release_lite.py", timestamp),
            b"raise RuntimeError('must never execute')\n",
        )


@pytest.mark.parametrize(
    "unsafe_path",
    [
        "C:/payload.txt",
        "C:payload.txt",
        "payload.txt:stream",
        'quote".txt',
        "less-than<.txt",
        "greater-than>.txt",
        "pipe|.txt",
        "question?.txt",
        "star*.txt",
        "CON",
        "aux.txt",
        "nested/NUL.json",
        "nested/COM1.csv",
        "trailing-dot.",
        "trailing-space ",
        "control\x01.txt",
        "cafe\u0301.txt",
    ],
)
def test_normalize_release_path_rejects_cross_platform_ambiguous_names(unsafe_path):
    with pytest.raises(ValueError):
        normalize_release_path(unsafe_path)


def test_zip_members_reject_case_insensitive_path_collision(tmp_path):
    zip_path = tmp_path / "collision.zip"
    with zipfile.ZipFile(zip_path, "w") as archive:
        archive.writestr("Payload.txt", "first")
        archive.writestr("payload.txt", "second")

    errors = zip_bundle.validate_zip_members(zip_path)

    assert any("case-insensitive member collision" in error for error in errors)


def test_zip_members_reject_raw_archive_above_limit_before_parsing(tmp_path):
    zip_path = tmp_path / "oversized.zip"
    zip_path.write_bytes(b"not-a-zip")

    errors = zip_bundle.validate_zip_members(zip_path, max_archive_bytes=4)

    assert errors == [f"Zip archive is {zip_path.stat().st_size} bytes; limit is 4"]


def test_acceptance_without_trusted_digest_never_extracts_or_executes(
    tmp_path, monkeypatch, capsys
):
    zip_path = tmp_path / RELEASE_ZIP_NAME
    write_structurally_canonical_zip(zip_path)

    def forbidden(*_args, **_kwargs):
        raise AssertionError("untrusted archive content must not be extracted or executed")

    monkeypatch.setattr(zip_bundle, "safe_extract", forbidden)
    monkeypatch.setattr(zip_bundle, "run_step", forbidden)

    code = zip_bundle.main([str(zip_path)])

    assert code == 1
    assert "--expected-sha256 is required" in capsys.readouterr().err


def test_digest_free_skip_acceptance_is_structure_only(tmp_path, monkeypatch):
    zip_path = tmp_path / RELEASE_ZIP_NAME
    report_path = tmp_path / "structure-report.json"
    write_structurally_canonical_zip(zip_path)

    def forbidden(*_args, **_kwargs):
        raise AssertionError("structure-only validation must not extract or execute")

    monkeypatch.setattr(zip_bundle, "safe_extract", forbidden)
    monkeypatch.setattr(zip_bundle, "run_step", forbidden)

    assert (
        zip_bundle.main(
            [str(zip_path), "--skip-acceptance", "--output", str(report_path)]
        )
        == 0
    )
    report = json.loads(report_path.read_text(encoding="utf-8"))
    assert report["validation_scope"] == "archive-structure-only"
    assert report["zip_sha256"] == zip_bundle.sha256(zip_path)
    assert report["steps"] == []


def test_archive_size_limit_is_enforced_before_hashing(tmp_path, monkeypatch):
    zip_path = tmp_path / RELEASE_ZIP_NAME
    zip_path.write_bytes(b"too-large")

    def forbidden_hash(_path):
        raise AssertionError("oversized archive must not be hashed")

    monkeypatch.setattr(zip_bundle, "sha256", forbidden_hash)

    code = zip_bundle.main(
        [
            str(zip_path),
            "--expected-sha256",
            "0" * 64,
            "--skip-acceptance",
            "--max-archive-bytes",
            "4",
        ]
    )

    assert code == 1


def artifact_metadata(digest):
    return {
        "schema_version": RELEASE_SCHEMA_VERSION,
        "zip_path": RELEASE_ZIP_NAME,
        "zip_bytes": 1,
        "zip_sha256": digest,
        "zip_acceptance_command": (
            f"{ZIP_ACCEPTANCE_COMMAND} --expected-sha256 {digest}"
        ),
    }


def test_acceptance_loads_strict_artifact_digest(tmp_path):
    digest = "a" * 64
    path = tmp_path / acceptance.ARTIFACTS_NAME
    path.write_text(json.dumps(artifact_metadata(digest)), encoding="utf-8")

    assert acceptance.load_trusted_zip_sha256(path) == digest

    metadata = artifact_metadata(digest)
    metadata["zip_acceptance_command"] = ZIP_ACCEPTANCE_COMMAND
    path.write_text(json.dumps(metadata), encoding="utf-8")
    with pytest.raises(ValueError, match="must include its trusted digest"):
        acceptance.load_trusted_zip_sha256(path)


def test_top_level_acceptance_passes_artifact_digest_to_zip_validator(
    tmp_path, monkeypatch
):
    digest = "b" * 64
    (tmp_path / "release-lite").mkdir()
    (tmp_path / RELEASE_ZIP_NAME).write_bytes(b"x")
    (tmp_path / acceptance.ARTIFACTS_NAME).write_text(
        json.dumps(artifact_metadata(digest)), encoding="utf-8"
    )
    commands = []

    def passing_step(label, cmd, required=True, timeout_seconds=300):
        commands.append((label, cmd))
        return {
            "label": label,
            "command": cmd,
            "cwd": str(tmp_path),
            "required": required,
            "returncode": 0,
            "status": "PASS",
            "duration_seconds": 0.0,
            "stdout": "",
            "stderr": "",
        }

    monkeypatch.setattr(acceptance, "ROOT", tmp_path)
    monkeypatch.setattr(acceptance, "run_step", passing_step)

    assert acceptance.main([]) == 0
    zip_commands = [command for label, command in commands if label == "zip_bundle"]
    assert zip_commands == [
        [
            acceptance.sys.executable,
            "validate_zip_bundle.py",
            RELEASE_ZIP_NAME,
            "--expected-sha256",
            digest,
        ]
    ]


def test_malformed_artifact_digest_fails_before_any_acceptance_subprocess(
    tmp_path, monkeypatch
):
    (tmp_path / "release-lite").mkdir()
    (tmp_path / RELEASE_ZIP_NAME).write_bytes(b"x")
    (tmp_path / acceptance.ARTIFACTS_NAME).write_text("[]\n", encoding="utf-8")

    def forbidden(*_args, **_kwargs):
        raise AssertionError("preflight failure must not start subprocesses")

    monkeypatch.setattr(acceptance, "ROOT", tmp_path)
    monkeypatch.setattr(acceptance, "run_step", forbidden)

    assert acceptance.main([]) == 1
