#!/usr/bin/env python3
"""Run end-to-end acceptance checks for the lightweight release."""
import argparse
import json
import shutil
import subprocess  # noqa: F401 -- re-exported: tests monkeypatch subprocess.run via this module
import sys
from datetime import datetime, timezone
from pathlib import Path

from release_tools.common import (
    RELEASE_SCHEMA_VERSION,
    RELEASE_ZIP_NAME,
    ZIP_ACCEPTANCE_COMMAND,
    append_timeout_message,  # noqa: F401 -- re-exported for tests/test_subprocess_reporting.py
    run_subprocess_step,
    subprocess_output_text,  # noqa: F401 -- re-exported for tests/test_subprocess_reporting.py
    write_json_report as _write_json_report,
)


ROOT = Path(__file__).resolve().parent
ZIP_NAME = RELEASE_ZIP_NAME
ARTIFACTS_NAME = "RELEASE_ARTIFACTS.json"


def run_step(label, cmd, required=True, timeout_seconds=300):
    return run_subprocess_step(label, cmd, ROOT, timeout_seconds=timeout_seconds,
                               required=required, prefix="acceptance")


def cleanup_transient_files():
    root = ROOT.resolve()
    targets = []
    for pattern in ("_smoke_*", "_acceptance_*", "*.pyc"):
        targets.extend(path for path in ROOT.rglob(pattern) if path.is_file())
    targets.extend(path for path in ROOT.rglob("__pycache__") if path.is_dir())
    for path in sorted(targets, key=lambda item: len(item.parts), reverse=True):
        resolved = path.resolve()
        if resolved != root and root not in resolved.parents:
            raise RuntimeError(f"Refusing to remove outside acceptance root: {resolved}")
        if path.is_dir():
            shutil.rmtree(path)
        elif path.exists():
            path.unlink()


def choose_release_validation():
    in_release_dir = (ROOT / "release_manifest.json").exists() and (ROOT / "SHA256SUMS.txt").exists()
    parent_zip = ROOT.parent / ZIP_NAME
    if in_release_dir and parent_zip.exists():
        return ["validate_release_lite.py", "--release-dir", ".", "--zip", str(parent_zip)]
    if in_release_dir:
        return ["validate_release_lite.py", "--release-dir", ".", "--no-zip"]
    if (ROOT / "release-lite").is_dir() and (ROOT / ZIP_NAME).exists():
        return [
            "validate_release_lite.py", "--release-dir", "release-lite",
            "--zip", ZIP_NAME, "--source-root", ".",
            "--artifacts", "RELEASE_ARTIFACTS.json",
        ]
    return None


def load_trusted_zip_sha256(path=None):
    """Load the canonical trusted digest used before executing ZIP payload code."""
    artifacts_path = Path(path) if path is not None else ROOT / ARTIFACTS_NAME
    if artifacts_path.is_symlink():
        raise ValueError(f"{ARTIFACTS_NAME} must not be a symbolic link")
    try:
        artifacts = json.loads(artifacts_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"Could not read {ARTIFACTS_NAME}: {exc}") from exc
    if not isinstance(artifacts, dict):
        raise ValueError(f"{ARTIFACTS_NAME} top-level value must be an object")
    if artifacts.get("schema_version") != RELEASE_SCHEMA_VERSION:
        raise ValueError(f"{ARTIFACTS_NAME} schema_version is not canonical")
    if artifacts.get("zip_path") != ZIP_NAME:
        raise ValueError(f"{ARTIFACTS_NAME} zip_path must be {ZIP_NAME!r}")
    zip_bytes = artifacts.get("zip_bytes")
    if not isinstance(zip_bytes, int) or isinstance(zip_bytes, bool) or zip_bytes < 0:
        raise ValueError(f"{ARTIFACTS_NAME} zip_bytes must be a non-negative integer")
    digest = artifacts.get("zip_sha256")
    if (
        not isinstance(digest, str)
        or len(digest) != 64
        or any(character not in "0123456789abcdef" for character in digest)
    ):
        raise ValueError(
            f"{ARTIFACTS_NAME} zip_sha256 must be 64 lowercase hexadecimal characters"
        )
    expected_command = f"{ZIP_ACCEPTANCE_COMMAND} --expected-sha256 {digest}"
    if artifacts.get("zip_acceptance_command") != expected_command:
        raise ValueError(
            f"{ARTIFACTS_NAME} zip_acceptance_command must include its trusted digest"
        )
    return digest


def preflight_failure_step(message):
    return {
        "label": "trusted_zip_digest",
        "command": [],
        "cwd": str(ROOT),
        "required": True,
        "returncode": 2,
        "status": "FAIL",
        "duration_seconds": 0.0,
        "stdout": "",
        "stderr": message,
    }


def write_json_report(path, report):
    _write_json_report(path, report, root=ROOT, prefix="acceptance")


def write_markdown_report(path, report):
    path = Path(path)
    if not path.is_absolute():
        path = ROOT / path
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = [
        "# Release acceptance report",
        "",
        f"- Status: `{report['status']}`",
        f"- Generated UTC: `{report['generated_utc']}`",
        f"- Root: `{report['root']}`",
        "",
        "| Step | Status | Required | Seconds |",
        "|---|---:|---:|---:|",
    ]
    for step in report["steps"]:
        rows.append(
            f"| {step['label']} | `{step['status']}` | {step['required']} | "
            f"{step['duration_seconds']:.1f} |"
        )
    if report["skipped"]:
        rows.extend(["", "## Skipped"])
        rows.extend(f"- {item}" for item in report["skipped"])
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    print(f"[acceptance] wrote {path}")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Run environment, smoke, safety, and release-integrity checks."
    )
    parser.add_argument("-o", "--output",
                        help="optional JSON report path; prefer a path outside release-lite/")
    parser.add_argument("--markdown",
                        help="optional Markdown summary path; prefer a path outside release-lite/")
    parser.add_argument("--timeout-seconds", type=int, default=300,
                        help="per-step subprocess timeout (default: 300)")
    args = parser.parse_args(argv)

    steps = []
    skipped = []
    commands = [
        ("environment", [sys.executable, "check_environment.py", "--self-test"]),
        ("dependency_audit", [sys.executable, "audit_lightweight_dependencies.py"]),
        ("cli_audit", [sys.executable, "audit_cli_entrypoints.py"]),
        ("docs_audit", [sys.executable, "audit_release_docs.py"]),
        ("output_contracts", [sys.executable, "validate_output_contracts.py"]),
        ("smoke", [sys.executable, "run_smoke_tests.py"]),
        ("safety", [sys.executable, "run_safety_tests.py"]),
    ]
    release_validation = choose_release_validation()
    if release_validation:
        commands.append(("release_validation", [sys.executable] + release_validation))
    else:
        skipped.append("release validation: no release-lite directory/zip or in-place manifest found")
    preflight_error = None
    if (ROOT / "release-lite").is_dir() and (ROOT / ZIP_NAME).exists():
        try:
            trusted_digest = load_trusted_zip_sha256()
        except ValueError as exc:
            preflight_error = str(exc)
        else:
            commands.append(
                (
                    "zip_bundle",
                    [
                        sys.executable,
                        "validate_zip_bundle.py",
                        ZIP_NAME,
                        "--expected-sha256",
                        trusted_digest,
                    ],
                )
            )

    if preflight_error is not None:
        print(f"[acceptance] ERROR: {preflight_error}", file=sys.stderr)
        steps.append(preflight_failure_step(preflight_error))
    else:
        for label, cmd in commands:
            if label in {"release_validation", "zip_bundle"}:
                cleanup_transient_files()
            steps.append(run_step(label, cmd, timeout_seconds=args.timeout_seconds))
            if steps[-1]["returncode"] != 0:
                break

    failed = [step for step in steps if step["required"] and step["returncode"] != 0]
    report = {
        "schema_version": "1.0",
        "status": "FAIL" if failed else "PASS",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "root": str(ROOT),
        "python": sys.executable,
        "steps": steps,
        "skipped": skipped,
    }

    if args.output:
        write_json_report(args.output, report)
    if args.markdown:
        write_markdown_report(args.markdown, report)

    print(f"[acceptance] {report['status']}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
