#!/usr/bin/env python3
"""Validate and safely extract the canonical lightweight release ZIP."""
import argparse
import json
import shutil
import subprocess  # noqa: F401 -- tests monkeypatch subprocess.run via this module
import sys
import tempfile
import zipfile
from datetime import datetime, timezone
from pathlib import Path

from release_tools.common import (
    DEFAULT_MAX_ZIP_ARCHIVE_BYTES,
    DEFAULT_MAX_ZIP_COMPRESSION_RATIO,
    DEFAULT_MAX_ZIP_ENTRIES,
    DEFAULT_MAX_ZIP_MEMBER_BYTES,
    DEFAULT_MAX_ZIP_TOTAL_BYTES,
    RELEASE_ZIP_NAME,
    append_timeout_message,  # noqa: F401 -- backwards-compatible test API
    canonical_zip_datetime,
    canonical_zip_errors,
    normalize_release_path,
    run_subprocess_step,
    sha256_file as sha256,
    subprocess_output_text,  # noqa: F401 -- backwards-compatible test API
    write_json_report,
    zip_safety_errors,
)


ROOT = Path(__file__).resolve().parent
ZIP_NAME = RELEASE_ZIP_NAME


def run_step(label, cmd, cwd, required=True, timeout_seconds=300):
    return run_subprocess_step(
        label,
        cmd,
        cwd,
        timeout_seconds=timeout_seconds,
        required=required,
        prefix="zip-bundle",
    )


def _zip_release_datetime(zf, infos, errors):
    matches = [info for info in infos if info.filename == "RELEASE_METADATA.json"]
    if len(matches) != 1:
        errors.append("Zip must contain exactly one RELEASE_METADATA.json file")
        return None
    try:
        metadata = json.loads(zf.read(matches[0]).decode("utf-8"))
    except (UnicodeError, json.JSONDecodeError, RuntimeError, zipfile.BadZipFile) as exc:
        errors.append(f"Could not parse zipped RELEASE_METADATA.json: {exc}")
        return None
    if not isinstance(metadata, dict):
        errors.append("Zipped RELEASE_METADATA.json top-level value must be an object")
        return None
    try:
        return canonical_zip_datetime(metadata.get("release_date"))
    except ValueError as exc:
        errors.append(f"Zipped RELEASE_METADATA.json {exc}")
        return None


def validate_zip_members(
    zip_path,
    *,
    max_archive_bytes=DEFAULT_MAX_ZIP_ARCHIVE_BYTES,
    max_entries=DEFAULT_MAX_ZIP_ENTRIES,
    max_member_bytes=DEFAULT_MAX_ZIP_MEMBER_BYTES,
    max_total_bytes=DEFAULT_MAX_ZIP_TOTAL_BYTES,
    max_compression_ratio=DEFAULT_MAX_ZIP_COMPRESSION_RATIO,
):
    errors = []
    zip_arg = Path(zip_path)
    if zip_arg.is_symlink():
        return [f"Zip archive must not be a symbolic link: {zip_arg}"]
    if not zip_arg.is_file():
        return [f"Zip archive not found: {zip_arg}"]
    try:
        archive_bytes = zip_arg.stat().st_size
    except OSError as exc:
        return [f"Could not stat zip archive: {exc}"]
    if archive_bytes > max_archive_bytes:
        return [f"Zip archive is {archive_bytes} bytes; limit is {max_archive_bytes}"]
    try:
        zf = zipfile.ZipFile(zip_arg)
    except (OSError, zipfile.BadZipFile) as exc:
        return [f"Could not read zip archive: {exc}"]
    with zf:
        infos = zf.infolist()
        safety_errors = zip_safety_errors(
            infos,
            max_entries=max_entries,
            max_member_bytes=max_member_bytes,
            max_total_bytes=max_total_bytes,
            max_compression_ratio=max_compression_ratio,
        )
        errors.extend(safety_errors)
        if safety_errors:
            return errors
        expected_datetime = _zip_release_datetime(zf, infos, errors)
        if expected_datetime is not None:
            errors.extend(canonical_zip_errors(zf, infos, expected_datetime))
        try:
            bad = zf.testzip()
        except (OSError, RuntimeError, zipfile.BadZipFile) as exc:
            errors.append(f"Could not test zip archive: {exc}")
            bad = None
        if bad is not None:
            errors.append(f"Zip archive is corrupt at {bad}")
    return errors


def safe_extract(
    zip_path,
    extract_dir,
    *,
    max_member_bytes=DEFAULT_MAX_ZIP_MEMBER_BYTES,
    max_total_bytes=DEFAULT_MAX_ZIP_TOTAL_BYTES,
):
    """Extract regular files with a second path/size check and no extractall."""
    extract_dir = Path(extract_dir)
    extract_root = extract_dir.resolve()
    extract_dir.mkdir(parents=True)
    declared_total = 0
    written_total = 0
    try:
        with zipfile.ZipFile(zip_path) as zf:
            for info in zf.infolist():
                rel = normalize_release_path(info.filename)
                if info.is_dir():
                    raise RuntimeError(f"Refusing to extract directory entry: {info.filename}")
                if info.file_size > max_member_bytes:
                    raise RuntimeError(f"Refusing oversized member: {info.filename}")
                declared_total += info.file_size
                if declared_total > max_total_bytes:
                    raise RuntimeError("Refusing archive above total extraction limit")
                target = (extract_dir / rel).resolve()
                if target != extract_root and extract_root not in target.parents:
                    raise RuntimeError(f"Refusing to extract outside target: {info.filename}")
                target.parent.mkdir(parents=True, exist_ok=True)
                with zf.open(info) as source, open(target, "xb") as destination:
                    member_written = 0
                    while True:
                        chunk = source.read(1024 * 1024)
                        if not chunk:
                            break
                        member_written += len(chunk)
                        written_total += len(chunk)
                        if member_written > max_member_bytes:
                            raise RuntimeError(
                                f"Member exceeded extraction limit while reading: {info.filename}"
                            )
                        if written_total > max_total_bytes:
                            raise RuntimeError(
                                "Archive exceeded total extraction limit while reading"
                            )
                        destination.write(chunk)
                    if member_written != info.file_size:
                        raise RuntimeError(
                            "Member size changed while extracting: "
                            f"{info.filename} (declared {info.file_size}, "
                            f"read {member_written})"
                        )
    except (OSError, ValueError, zipfile.BadZipFile) as exc:
        raise RuntimeError(f"Could not extract zip archive: {exc}") from exc


def copy_archive_bounded(source, destination, max_archive_bytes):
    """Copy an archive once, rejecting growth past the pre-hash size ceiling."""
    copied = 0
    try:
        with open(source, "rb") as input_handle, open(destination, "xb") as output_handle:
            while True:
                chunk = input_handle.read(1024 * 1024)
                if not chunk:
                    break
                copied += len(chunk)
                if copied > max_archive_bytes:
                    raise RuntimeError(
                        f"Zip archive exceeded {max_archive_bytes} bytes while copying"
                    )
                output_handle.write(chunk)
    except OSError as exc:
        raise RuntimeError(f"Could not copy zip archive: {exc}") from exc
    return copied


def cleanup_temp(temp_root):
    temp_root = Path(temp_root).resolve()
    temp_parent = Path(tempfile.gettempdir()).resolve()
    if temp_root.parent != temp_parent or not temp_root.name.startswith("tcga_zip_bundle_"):
        raise RuntimeError(f"Refusing to remove unexpected temp directory: {temp_root}")
    shutil.rmtree(temp_root)


def write_report(path, report):
    write_json_report(path, report, root=ROOT, prefix="zip-bundle")


def resolve_zip_path(path_arg):
    zip_path = Path(path_arg)
    if zip_path.is_absolute():
        return zip_path.absolute()
    for base in (ROOT, ROOT.parent):
        candidate = (base / zip_path).absolute()
        if candidate.exists():
            return candidate
    return (ROOT / zip_path).absolute()


def normalize_expected_sha(value):
    if value is None:
        return None
    digest = value.lower().removeprefix("sha256:")
    if len(digest) != 64 or any(ch not in "0123456789abcdef" for ch in digest):
        raise ValueError("--expected-sha256 must be a 64-character hexadecimal digest")
    return digest


def report_validation_scope(expected_sha, actual_sha, skip_acceptance):
    if expected_sha is not None and actual_sha == expected_sha:
        return "trusted-extracted-release"
    if expected_sha is not None and actual_sha is not None:
        return "digest-rejected"
    if skip_acceptance and actual_sha is not None:
        return "archive-structure-only"
    return "not-validated"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("zip_path", nargs="?", default=ZIP_NAME)
    parser.add_argument("-o", "--output", help="optional JSON report path")
    parser.add_argument(
        "--expected-sha256",
        help="trusted archive digest required before any extracted code may run",
    )
    parser.add_argument(
        "--keep-extracted",
        action="store_true",
        help="retain the trusted temporary archive copy and extraction directory",
    )
    parser.add_argument(
        "--skip-acceptance",
        action="store_true",
        help=(
            "skip extracted acceptance; without a trusted digest, perform only "
            "non-extracting structural validation"
        ),
    )
    parser.add_argument("--timeout-seconds", type=int, default=300)
    parser.add_argument(
        "--max-archive-bytes", type=int, default=DEFAULT_MAX_ZIP_ARCHIVE_BYTES
    )
    parser.add_argument("--max-entries", type=int, default=DEFAULT_MAX_ZIP_ENTRIES)
    parser.add_argument(
        "--max-member-bytes", type=int, default=DEFAULT_MAX_ZIP_MEMBER_BYTES
    )
    parser.add_argument("--max-total-bytes", type=int, default=DEFAULT_MAX_ZIP_TOTAL_BYTES)
    parser.add_argument(
        "--max-compression-ratio",
        type=float,
        default=DEFAULT_MAX_ZIP_COMPRESSION_RATIO,
    )
    args = parser.parse_args(argv)

    zip_path = resolve_zip_path(args.zip_path)
    steps = []
    errors = []
    temp_root = None
    extract_dir = None
    zip_copy = None
    expected_sha = None
    actual_sha = None
    zip_bytes = None
    try:
        try:
            expected_sha = normalize_expected_sha(args.expected_sha256)
        except ValueError as exc:
            errors.append(str(exc))
        if not args.skip_acceptance and args.expected_sha256 is None:
            errors.append(
                "--expected-sha256 is required before extracted release code may run; "
                "use --skip-acceptance for structure-only validation"
            )

        if zip_path.is_symlink():
            errors.append(f"Zip archive must not be a symbolic link: {zip_path}")
        elif not zip_path.is_file():
            errors.append(f"Zip archive not found: {zip_path}")
        else:
            try:
                zip_bytes = zip_path.stat().st_size
            except OSError as exc:
                errors.append(f"Could not stat zip archive: {exc}")
            if zip_bytes is not None and zip_bytes > args.max_archive_bytes:
                errors.append(
                    f"Zip archive is {zip_bytes} bytes; limit is {args.max_archive_bytes}"
                )

        if not errors:
            validation_path = zip_path
            temp_root = Path(tempfile.mkdtemp(prefix="tcga_zip_bundle_"))
            zip_copy = temp_root / ZIP_NAME
            try:
                zip_bytes = copy_archive_bounded(
                    zip_path, zip_copy, args.max_archive_bytes
                )
            except RuntimeError as exc:
                errors.append(str(exc))
            if not errors:
                actual_sha = sha256(zip_copy)
                if expected_sha is not None and actual_sha != expected_sha:
                    errors.append(
                        "Zip SHA256 mismatch: "
                        f"expected {expected_sha}, found {actual_sha}"
                    )
                validation_path = zip_copy

        if not errors:
            errors.extend(
                validate_zip_members(
                    validation_path,
                    max_archive_bytes=args.max_archive_bytes,
                    max_entries=args.max_entries,
                    max_member_bytes=args.max_member_bytes,
                    max_total_bytes=args.max_total_bytes,
                    max_compression_ratio=args.max_compression_ratio,
                )
            )

        # A digest-free --skip-acceptance run is deliberately structure-only:
        # never extract or execute Python supplied by an untrusted archive.
        if not errors and expected_sha is not None:
            if temp_root is None or zip_copy is None:
                raise RuntimeError("Internal error: trusted archive copy was not created")
            extract_dir = temp_root / "release-lite"
            try:
                safe_extract(
                    zip_copy,
                    extract_dir,
                    max_member_bytes=args.max_member_bytes,
                    max_total_bytes=args.max_total_bytes,
                )
            except RuntimeError as exc:
                errors.append(str(exc))
            if not errors:
                steps.append(
                    run_step(
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
                    )
                )
                if steps[-1]["returncode"] != 0:
                    errors.append("Release validation failed after zip extraction.")
            if not args.skip_acceptance and not errors:
                steps.append(
                    run_step(
                        "extracted_acceptance",
                        [sys.executable, "run_release_acceptance.py"],
                        cwd=extract_dir,
                        timeout_seconds=args.timeout_seconds,
                    )
                )
                if steps[-1]["returncode"] != 0:
                    errors.append("Acceptance failed inside extracted zip bundle.")

        report = {
            "schema_version": "1.0",
            "status": "FAIL" if errors else "PASS",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "zip_path": str(zip_path),
            "zip_bytes": zip_bytes,
            "zip_sha256": actual_sha,
            "expected_sha256": expected_sha,
            "validation_scope": report_validation_scope(
                expected_sha, actual_sha, args.skip_acceptance
            ),
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
