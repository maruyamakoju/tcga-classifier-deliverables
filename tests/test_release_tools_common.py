"""Unit tests for release_tools.common, shared by the audit/validate/build scripts."""
import json
import subprocess

from release_tools.common import (
    add_message,
    append_timeout_message,
    exit_code_for_status,
    load_manifest_paths,
    release_python_files,
    release_target_root,
    run_subprocess_step,
    sha256_file,
    status_from_levels,
    subprocess_output_text,
    write_json_report,
)


def test_add_message_with_and_without_path():
    messages = []
    add_message(messages, "ERROR", "some_code", "some message")
    add_message(messages, "WARNING", "other_code", "other message", path="a/b.txt")
    assert messages[0] == {"level": "ERROR", "code": "some_code", "message": "some message"}
    assert messages[1] == {
        "level": "WARNING", "code": "other_code", "message": "other message", "path": "a/b.txt",
    }


def test_status_from_levels_precedence():
    assert status_from_levels([]) == "PASS"
    assert status_from_levels([{"level": "WARNING"}]) == "WARN"
    assert status_from_levels([{"level": "WARNING"}, {"level": "ERROR"}]) == "FAIL"


def test_exit_code_for_status():
    assert exit_code_for_status("PASS", strict=False) == 0
    assert exit_code_for_status("WARN", strict=False) == 0
    assert exit_code_for_status("WARN", strict=True) == 1
    assert exit_code_for_status("FAIL", strict=False) == 1
    assert exit_code_for_status("FAIL", strict=True) == 1


def test_sha256_file_matches_hashlib(tmp_path):
    import hashlib

    path = tmp_path / "data.bin"
    path.write_bytes(b"some bytes to hash")
    assert sha256_file(path) == hashlib.sha256(b"some bytes to hash").hexdigest()


def test_write_json_report_creates_parent_dirs_and_trailing_newline(tmp_path):
    out_path = tmp_path / "nested" / "deep" / "report.json"
    write_json_report(out_path, {"status": "PASS"})
    assert out_path.exists()
    text = out_path.read_text(encoding="utf-8")
    assert text.endswith("\n")
    assert json.loads(text) == {"status": "PASS"}


def test_write_json_report_resolves_relative_path_against_root(tmp_path):
    write_json_report("report.json", {"a": 1}, root=tmp_path)
    assert (tmp_path / "report.json").exists()


def test_release_target_root_prefers_nested_release_lite(tmp_path):
    assert release_target_root(tmp_path) == tmp_path

    nested = tmp_path / "release-lite"
    nested.mkdir()
    (nested / "release_manifest.json").write_text("{}", encoding="utf-8")
    (nested / "SHA256SUMS.txt").write_text("", encoding="utf-8")
    assert release_target_root(tmp_path) == nested

    (tmp_path / "release_manifest.json").write_text("{}", encoding="utf-8")
    (tmp_path / "SHA256SUMS.txt").write_text("", encoding="utf-8")
    assert release_target_root(tmp_path) == tmp_path


def test_load_manifest_paths_missing_and_malformed(tmp_path):
    assert load_manifest_paths(tmp_path) is None

    (tmp_path / "release_manifest.json").write_text("{not json", encoding="utf-8")
    assert load_manifest_paths(tmp_path) is None

    (tmp_path / "release_manifest.json").write_text(
        json.dumps({"files": [{"path": "a.py"}, {"path": "b.md"}]}), encoding="utf-8"
    )
    assert load_manifest_paths(tmp_path) == ["a.py", "b.md"]


def test_release_python_files_uses_manifest_when_present(tmp_path):
    (tmp_path / "a.py").write_text("", encoding="utf-8")
    (tmp_path / "b.py").write_text("", encoding="utf-8")
    (tmp_path / "release_manifest.json").write_text(
        json.dumps({"files": [{"path": "a.py"}, {"path": "not_here.py"}, {"path": "c.md"}]}),
        encoding="utf-8",
    )
    files = release_python_files(tmp_path)
    assert [p.name for p in files] == ["a.py"]


def test_release_python_files_falls_back_to_glob(tmp_path):
    (tmp_path / "a.py").write_text("", encoding="utf-8")
    (tmp_path / "b.py").write_text("", encoding="utf-8")
    files = release_python_files(tmp_path)
    assert sorted(p.name for p in files) == ["a.py", "b.py"]


def test_subprocess_output_text_accepts_none_bytes_and_str():
    assert subprocess_output_text(None) == ""
    assert subprocess_output_text(b"abc") == "abc"
    assert subprocess_output_text("abc") == "abc"


def test_append_timeout_message_avoids_extra_blank_line():
    assert append_timeout_message("", 4) == "Timed out after 4s"
    assert append_timeout_message("stderr\n", 3) == "stderr\nTimed out after 3s"


def test_run_subprocess_step_reports_pass_and_fail(tmp_path):
    ok = run_subprocess_step("ok", ["python", "-c", "print('hi')"], tmp_path)
    assert ok["status"] == "PASS"
    assert ok["returncode"] == 0
    assert ok["required"] is True
    assert "hi" in ok["stdout"]

    bad = run_subprocess_step("bad", ["python", "-c", "import sys; sys.exit(3)"], tmp_path)
    assert bad["status"] == "FAIL"
    assert bad["returncode"] == 3


def test_run_subprocess_step_not_required_warns_instead_of_fails(tmp_path):
    result = run_subprocess_step(
        "optional", ["python", "-c", "import sys; sys.exit(1)"], tmp_path, required=False,
    )
    assert result["status"] == "WARN"
    assert result["required"] is False


def test_run_subprocess_step_handles_timeout(monkeypatch, tmp_path):
    def timeout_run(*args, **kwargs):
        raise subprocess.TimeoutExpired(
            cmd=args[0], timeout=kwargs.get("timeout"),
            output=b"partial stdout\n", stderr=b"partial stderr\n",
        )

    import release_tools.common as common
    monkeypatch.setattr(common.subprocess, "run", timeout_run)

    result = run_subprocess_step("timeout", ["demo"], tmp_path, timeout_seconds=1)
    assert result["returncode"] == 124
    assert result["status"] == "FAIL"
    assert result["stdout"] == "partial stdout\n"
    assert result["stderr"] == "partial stderr\nTimed out after 1s"
