"""Regression tests for deterministic staging and rollback-safe release builds."""
import os
import stat
import zipfile
from pathlib import Path

import pytest

import build_release_lite as builder
from release_tools.common import CANONICAL_ZIP_CREATE_SYSTEM, sha256_file


def test_write_zip_is_byte_reproducible_and_host_neutral(tmp_path):
    release_dir = tmp_path / "release-lite"
    release_dir.mkdir()
    (release_dir / "z.txt").write_text("z\n", encoding="utf-8")
    (release_dir / "a.txt").write_text("a\n", encoding="utf-8")
    first = tmp_path / "first.zip"
    second = tmp_path / "second.zip"

    builder.write_zip(release_dir, first, "2026-07-09")
    builder.write_zip(release_dir, second, "2026-07-09")

    assert sha256_file(first) == sha256_file(second)
    with zipfile.ZipFile(first) as zf:
        assert zf.namelist() == ["a.txt", "z.txt"]
        for info in zf.infolist():
            assert info.create_system == CANONICAL_ZIP_CREATE_SYSTEM
            assert info.compress_type == zipfile.ZIP_STORED
            assert stat.S_IFMT(info.external_attr >> 16) == stat.S_IFREG
            assert stat.S_IMODE(info.external_attr >> 16) == 0o644
            assert info.date_time == (2026, 7, 9, 0, 0, 0)
            assert info.extra == b""
            assert info.comment == b""


def _prepare_old_and_staged_outputs(tmp_path, monkeypatch):
    monkeypatch.setattr(builder, "ROOT", tmp_path)
    monkeypatch.setattr(builder, "RELEASE_DIR", tmp_path / "release-lite")
    monkeypatch.setattr(builder, "ZIP_PATH", tmp_path / "bundle.zip")
    monkeypatch.setattr(builder, "ARTIFACTS_PATH", tmp_path / "RELEASE_ARTIFACTS.json")
    builder.RELEASE_DIR.mkdir()
    (builder.RELEASE_DIR / "payload.txt").write_text("old-release\n", encoding="utf-8")
    builder.ZIP_PATH.write_bytes(b"old-zip")
    builder.ARTIFACTS_PATH.write_text("old-artifacts\n", encoding="utf-8")
    stage = tmp_path / "stage"
    stage.mkdir()
    staged_release = stage / "release-lite"
    staged_release.mkdir()
    (staged_release / "payload.txt").write_text("new-release\n", encoding="utf-8")
    staged_zip = stage / "bundle.zip"
    staged_zip.write_bytes(b"new-zip")
    staged_artifacts = stage / "RELEASE_ARTIFACTS.json"
    staged_artifacts.write_text("new-artifacts\n", encoding="utf-8")
    staged = {
        "release_dir": staged_release,
        "zip_path": staged_zip,
        "artifacts_path": staged_artifacts,
    }
    return stage, staged


def test_atomic_replace_rolls_back_every_output_on_install_failure(tmp_path, monkeypatch):
    stage, staged = _prepare_old_and_staged_outputs(tmp_path, monkeypatch)
    real_replace = os.replace
    calls = 0

    def fail_fourth_replace(source, destination):
        nonlocal calls
        calls += 1
        if calls == 4:
            raise OSError("simulated install failure")
        return real_replace(source, destination)

    monkeypatch.setattr(builder.os, "replace", fail_fourth_replace)

    with pytest.raises(RuntimeError, match="atomically"):
        builder.atomic_replace_outputs(staged, stage)

    assert (builder.RELEASE_DIR / "payload.txt").read_text(encoding="utf-8") == "old-release\n"
    assert builder.ZIP_PATH.read_bytes() == b"old-zip"
    assert builder.ARTIFACTS_PATH.read_text(encoding="utf-8") == "old-artifacts\n"


def test_stage_validation_failure_leaves_existing_outputs_unchanged(tmp_path, monkeypatch):
    _prepare_old_and_staged_outputs(tmp_path, monkeypatch)

    def fail_validation(*args, **kwargs):
        raise RuntimeError("staging validation failed")

    monkeypatch.setattr(builder, "build_stage", fail_validation)

    with pytest.raises(RuntimeError, match="staging validation failed"):
        builder.build()

    assert (builder.RELEASE_DIR / "payload.txt").read_text(encoding="utf-8") == "old-release\n"
    assert builder.ZIP_PATH.read_bytes() == b"old-zip"
    assert builder.ARTIFACTS_PATH.read_text(encoding="utf-8") == "old-artifacts\n"


def test_check_mode_never_installs_staged_outputs(tmp_path, monkeypatch):
    _prepare_old_and_staged_outputs(tmp_path, monkeypatch)

    def fake_build_stage(stage_root, **kwargs):
        stage_root = Path(stage_root)
        release = stage_root / "release-lite"
        release.mkdir()
        (release / "payload.txt").write_text("candidate\n", encoding="utf-8")
        zip_path = stage_root / "bundle.zip"
        zip_path.write_bytes(b"candidate")
        artifacts = stage_root / "RELEASE_ARTIFACTS.json"
        artifacts.write_text("candidate\n", encoding="utf-8")
        return {
            "release_dir": release,
            "zip_path": zip_path,
            "artifacts_path": artifacts,
            "manifest_count": 1,
            "checksum_count": 1,
            "zip_entries": 1,
        }

    monkeypatch.setattr(builder, "build_stage", fake_build_stage)
    monkeypatch.setattr(
        builder,
        "atomic_replace_outputs",
        lambda *args, **kwargs: pytest.fail("check mode attempted to install outputs"),
    )

    assert builder.build(check=True) is False
    assert (builder.RELEASE_DIR / "payload.txt").read_text(encoding="utf-8") == "old-release\n"
    assert builder.ZIP_PATH.read_bytes() == b"old-zip"
