#!/usr/bin/env python3
"""Run end-to-end acceptance checks for the lightweight release."""
import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent
ZIP_NAME = "tcga-tumor-normal-release-lite.zip"


def run_step(label, cmd, required=True, timeout_seconds=300):
    print(f"[acceptance] {label}: {' '.join(str(x) for x in cmd)}")
    started = time.perf_counter()
    env = dict(os.environ)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    try:
        result = subprocess.run(
            cmd, cwd=ROOT, text=True, capture_output=True, env=env,
            timeout=timeout_seconds,
        )
    except subprocess.TimeoutExpired as exc:
        duration = time.perf_counter() - started
        stdout = exc.stdout or ""
        stderr = exc.stderr or ""
        print(f"[acceptance] {label}: FAIL timeout after {duration:.1f}s", file=sys.stderr)
        return {
            "label": label,
            "command": [str(x) for x in cmd],
            "required": required,
            "returncode": 124,
            "status": "FAIL",
            "duration_seconds": round(duration, 3),
            "stdout": stdout,
            "stderr": stderr + f"\nTimed out after {timeout_seconds}s",
        }
    duration = time.perf_counter() - started
    status = "PASS" if result.returncode == 0 else "FAIL"
    if not required and result.returncode != 0:
        status = "WARN"
    if result.stdout:
        print(result.stdout.rstrip())
    if result.stderr:
        print(result.stderr.rstrip(), file=sys.stderr)
    print(f"[acceptance] {label}: {status} ({duration:.1f}s)")
    return {
        "label": label,
        "command": [str(x) for x in cmd],
        "required": required,
        "returncode": result.returncode,
        "status": status,
        "duration_seconds": round(duration, 3),
        "stdout": result.stdout,
        "stderr": result.stderr,
    }


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


def write_json_report(path, report):
    path = Path(path)
    if not path.is_absolute():
        path = ROOT / path
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[acceptance] wrote {path}")


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
    if (ROOT / "release-lite").is_dir() and (ROOT / ZIP_NAME).exists():
        commands.append(("zip_bundle", [sys.executable, "validate_zip_bundle.py", ZIP_NAME]))

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
