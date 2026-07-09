"""Regression tests for subprocess timeout reporting helpers."""
import subprocess

import run_release_acceptance
import validate_zip_bundle


def timeout_run(*args, **kwargs):
    raise subprocess.TimeoutExpired(
        cmd=args[0],
        timeout=kwargs.get("timeout"),
        output=b"partial stdout\n",
        stderr=b"partial stderr\n",
    )


def test_acceptance_run_step_decodes_timeout_bytes(monkeypatch, capsys):
    monkeypatch.setattr(run_release_acceptance.subprocess, "run", timeout_run)

    result = run_release_acceptance.run_step("timeout", ["demo"], timeout_seconds=1)

    captured = capsys.readouterr()
    assert "FAIL timeout" in captured.err
    assert result["returncode"] == 124
    assert result["status"] == "FAIL"
    assert result["stdout"] == "partial stdout\n"
    assert result["stderr"] == "partial stderr\nTimed out after 1s"


def test_zip_bundle_run_step_decodes_timeout_bytes(monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(validate_zip_bundle.subprocess, "run", timeout_run)

    result = validate_zip_bundle.run_step(
        "timeout",
        ["demo"],
        cwd=tmp_path,
        timeout_seconds=2,
    )

    captured = capsys.readouterr()
    assert "FAIL timeout" in captured.err
    assert result["returncode"] == 124
    assert result["status"] == "FAIL"
    assert result["stdout"] == "partial stdout\n"
    assert result["stderr"] == "partial stderr\nTimed out after 2s"


def test_subprocess_output_text_accepts_none_bytes_and_strings():
    assert run_release_acceptance.subprocess_output_text(None) == ""
    assert run_release_acceptance.subprocess_output_text(b"abc") == "abc"
    assert run_release_acceptance.subprocess_output_text("abc") == "abc"
    assert validate_zip_bundle.subprocess_output_text(None) == ""
    assert validate_zip_bundle.subprocess_output_text(b"abc") == "abc"
    assert validate_zip_bundle.subprocess_output_text("abc") == "abc"


def test_append_timeout_message_avoids_extra_blank_line():
    assert (
        run_release_acceptance.append_timeout_message("stderr\n", 3)
        == "stderr\nTimed out after 3s"
    )
    assert (
        validate_zip_bundle.append_timeout_message("", 4)
        == "Timed out after 4s"
    )
