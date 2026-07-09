#!/usr/bin/env python3
"""Extract the release zip into a clean temp directory and run acceptance checks."""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
import zipfile
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent
ZIP_NAME = "tcga-tumor-normal-release-lite.zip"


def subprocess_output_text(value):
    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return str(value)


def append_timeout_message(stderr, timeout_seconds):
    message = f"Timed out after {timeout_seconds}s"
    if stderr:
        return stderr.rstrip("\n") + "\n" + message
    return message


def sha256(path):
    import hashlib

    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run_step(label, cmd, cwd, timeout_seconds=300):
    print(f"[zip-bundle] {label}: {' '.join(str(x) for x in cmd)}")
    started = time.perf_counter()
    env = dict(os.environ)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    try:
        result = subprocess.run(
            cmd, cwd=cwd, text=True, capture_output=True, env=env,
            timeout=timeout_seconds,
        )
    except subprocess.TimeoutExpired as exc:
        duration = time.perf_counter() - started
        stdout = subprocess_output_text(exc.stdout)
        stderr = subprocess_output_text(exc.stderr)
        print(f"[zip-bundle] {label}: FAIL timeout after {duration:.1f}s", file=sys.stderr)
        return {
            "label": label,
            "command": [str(x) for x in cmd],
            "cwd": str(cwd),
            "returncode": 124,
            "status": "FAIL",
            "duration_seconds": round(duration, 3),
            "stdout": stdout,
            "stderr": append_timeout_message(stderr, timeout_seconds),
        }
    duration = time.perf_counter() - started
    if result.stdout:
        print(result.stdout.rstrip())
    if result.stderr:
        print(result.stderr.rstrip(), file=sys.stderr)
    status = "PASS" if result.returncode == 0 else "FAIL"
    print(f"[zip-bundle] {label}: {status} ({duration:.1f}s)")
    return {
        "label": label,
        "command": [str(x) for x in cmd],
        "cwd": str(cwd),
        "returncode": result.returncode,
        "status": status,
        "duration_seconds": round(duration, 3),
        "stdout": result.stdout,
        "stderr": result.stderr,
    }


def validate_zip_members(zip_path):
    errors = []
    with zipfile.ZipFile(zip_path) as zf:
        bad = zf.testzip()
        if bad is not None:
            errors.append(f"Zip archive is corrupt at {bad}")
        for info in zf.infolist():
            name = info.filename.replace("\\", "/")
            path = Path(name)
            if name.startswith("/") or path.is_absolute():
                errors.append(f"Zip contains absolute path: {info.filename}")
            if any(part == ".." for part in path.parts):
                errors.append(f"Zip contains parent traversal path: {info.filename}")
            if info.is_dir():
                continue
            if not name or name.endswith("/"):
                errors.append(f"Zip contains malformed file path: {info.filename}")
    return errors


def safe_extract(zip_path, extract_dir):
    extract_root = extract_dir.resolve()
    extract_dir.mkdir(parents=True)
    with zipfile.ZipFile(zip_path) as zf:
        for info in zf.infolist():
            target = (extract_dir / info.filename).resolve()
            if target != extract_root and extract_root not in target.parents:
                raise RuntimeError(f"Refusing to extract outside target: {info.filename}")
        zf.extractall(extract_dir)


def cleanup_temp(temp_root):
    temp_root = temp_root.resolve()
    temp_parent = Path(tempfile.gettempdir()).resolve()
    if temp_root.parent != temp_parent or not temp_root.name.startswith("tcga_zip_bundle_"):
        raise RuntimeError(f"Refusing to remove unexpected temp directory: {temp_root}")
    shutil.rmtree(temp_root, ignore_errors=True)


def write_report(path, report):
    out_path = Path(path)
    if not out_path.is_absolute():
        out_path = ROOT / out_path
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")
    print(f"[zip-bundle] wrote {out_path}")


def resolve_zip_path(path_arg):
    zip_path = Path(path_arg)
    if zip_path.is_absolute():
        return zip_path.resolve()
    for base in (ROOT, ROOT.parent):
        candidate = (base / zip_path).resolve()
        if candidate.exists():
            return candidate
    return (ROOT / zip_path).resolve()


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Validate a release zip after extraction into a clean temp directory."
    )
    parser.add_argument("zip_path", nargs="?", default=ZIP_NAME)
    parser.add_argument("-o", "--output", help="optional JSON report path")
    parser.add_argument("--keep-extracted", action="store_true",
                        help="do not delete the temporary extraction directory")
    parser.add_argument("--skip-acceptance", action="store_true",
                        help="only validate zip structure and release manifest/checksums")
    parser.add_argument("--timeout-seconds", type=int, default=300,
                        help="per-step subprocess timeout (default: 300)")
    args = parser.parse_args(argv)

    zip_path = resolve_zip_path(args.zip_path)
    steps = []
    errors = []
    temp_root = None
    extract_dir = None
    zip_copy = None
    try:
        if not zip_path.exists():
            errors.append(f"Zip archive not found: {zip_path}")
        else:
            errors.extend(validate_zip_members(zip_path))

        if not errors:
            temp_root = Path(tempfile.mkdtemp(prefix="tcga_zip_bundle_"))
            extract_dir = temp_root / "release-lite"
            zip_copy = temp_root / ZIP_NAME
            shutil.copy2(zip_path, zip_copy)
            safe_extract(zip_copy, extract_dir)
            steps.append(run_step(
                "release_validation",
                [
                    sys.executable,
                    str(extract_dir / "validate_release_lite.py"),
                    "--release-dir",
                    str(extract_dir),
                    "--zip",
                    str(zip_copy),
                ],
                cwd=temp_root,
                timeout_seconds=args.timeout_seconds,
            ))
            if steps[-1]["returncode"] != 0:
                errors.append("Release validation failed after zip extraction.")
            if not args.skip_acceptance and not errors:
                steps.append(run_step(
                    "extracted_acceptance",
                    [sys.executable, "run_release_acceptance.py"],
                    cwd=extract_dir,
                    timeout_seconds=args.timeout_seconds,
                ))
                if steps[-1]["returncode"] != 0:
                    errors.append("Acceptance failed inside extracted zip bundle.")

        report = {
            "schema_version": "1.0",
            "status": "FAIL" if errors else "PASS",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "zip_path": str(zip_path),
            "zip_bytes": zip_path.stat().st_size if zip_path.exists() else None,
            "zip_sha256": sha256(zip_path) if zip_path.exists() else None,
            "temp_root": str(temp_root) if temp_root else None,
            "extract_dir": str(extract_dir) if extract_dir else None,
            "steps": steps,
            "errors": errors,
        }
        for error in errors:
            print(f"[zip-bundle] ERROR: {error}", file=sys.stderr)
        if args.output:
            write_report(args.output, report)
        print(f"[zip-bundle] {report['status']}")
        return 1 if errors else 0
    finally:
        if temp_root is not None and not args.keep_extracted:
            cleanup_temp(temp_root)


if __name__ == "__main__":
    raise SystemExit(main())
