#!/usr/bin/env python3
"""Validate the lightweight release directory, canonical ZIP, and sidecar."""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import zipfile
from pathlib import Path

from release_tools.common import (
    DEFAULT_MAX_ZIP_ARCHIVE_BYTES,
    DEFAULT_MAX_ZIP_COMPRESSION_RATIO,
    DEFAULT_MAX_ZIP_ENTRIES,
    DEFAULT_MAX_ZIP_MEMBER_BYTES,
    DEFAULT_MAX_ZIP_TOTAL_BYTES,
    FORBIDDEN_NAMES,
    RELEASE_BUNDLE_NAME,
    RELEASE_FILES,
    RELEASE_SCHEMA_VERSION,
    RELEASE_VALIDATION_COMMAND,
    RELEASE_ZIP_NAME,
    ZIP_ACCEPTANCE_COMMAND,
    canonical_zip_datetime,
    canonical_zip_errors,
    normalize_release_path,
    release_path_collision_key,
    sha256_file,
    zip_safety_errors,
)


REQUIRED_FALLBACK = {
    "audit_cli_entrypoints.py",
    "audit_lightweight_dependencies.py",
    "audit_release_docs.py",
    "calibrate_threshold.py",
    "check_environment.py",
    "cohort_adapt_score.py",
    "DATA_DICTIONARY.md",
    "deployable_lr_weights.npz",
    "EXECUTIVE_SUMMARY.md",
    ".zenodo.json",
    "codemeta.json",
    "example_input.csv",
    "example_labels.csv",
    "example_output.csv",
    "example_workflow_output/workflow_report.md",
    "explain_scores.py",
    "inspect_expression_input.py",
    "model_gene_metadata.csv",
    "model_qc_reference.json",
    "README.md",
    "RELEASE_BUNDLE.md",
    "RELEASE_METADATA.json",
    "REPRODUCIBILITY.md",
    "requirements-light.txt",
    "release_manifest.json",
    "run_release_acceptance.py",
    "run_safety_tests.py",
    "run_smoke_tests.py",
    "run_tumor_normal_workflow.py",
    "score_tumor_normal.py",
    "SHA256SUMS.txt",
    "templates/input_matrix_template.csv",
    "templates/labels_template.csv",
    "TROUBLESHOOTING.md",
    "USER_GUIDE.md",
    "validate_output_contracts.py",
    "VERSION",
    "validate_zip_bundle.py",
}
GENERATED_RELEASE_FILES = {"release_manifest.json", "SHA256SUMS.txt"}
assert REQUIRED_FALLBACK - GENERATED_RELEASE_FILES <= set(RELEASE_FILES), (
    "REQUIRED_FALLBACK has entries release_tools.common.RELEASE_FILES does not ship"
)

EXPECTED_MANIFEST_SCHEMA_VERSION = RELEASE_SCHEMA_VERSION
EXPECTED_BUNDLE_NAME = RELEASE_BUNDLE_NAME
TRANSIENT_NAMES = {"__pycache__"}
TRANSIENT_PREFIXES = ("_smoke_", "_acceptance_")
TRANSIENT_SUFFIXES = (".pyc",)
BINARY_SUFFIXES = {".npy", ".npz", ".pkl", ".png", ".zip"}
TEXT_NAMES = {"LICENSE", "VERSION"}


def sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()


def sha256_zip_member(zf, info):
    digest = hashlib.sha256()
    with zf.open(info) as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def release_files(release_dir):
    release_dir = Path(release_dir)
    return {
        path.relative_to(release_dir).as_posix(): path
        for path in release_dir.rglob("*")
        if path.is_file() or path.is_symlink()
    }


def is_text_release_path(rel):
    path = Path(rel)
    if path.name in TEXT_NAMES:
        return True
    return path.suffix.lower() not in BINARY_SUFFIXES


def parse_sha256sums(path):
    expected = {}
    collision_keys = {}
    text = Path(path).read_text(encoding="ascii")
    for lineno, line in enumerate(text.splitlines(), start=1):
        if not line.strip():
            continue
        if "  " not in line:
            raise ValueError(f"Malformed SHA256SUMS line {lineno}: {line!r}")
        digest, rel = line.split("  ", 1)
        rel = normalize_release_path(rel)
        if len(digest) != 64 or any(ch not in "0123456789abcdef" for ch in digest):
            raise ValueError(f"Malformed SHA256 digest on line {lineno}: {digest!r}")
        if rel in expected:
            raise ValueError(f"Duplicate SHA256SUMS entry on line {lineno}: {rel}")
        collision_key = release_path_collision_key(rel)
        previous = collision_keys.get(collision_key)
        if previous is not None:
            raise ValueError(
                "Case-insensitive SHA256SUMS path collision on line "
                f"{lineno}: {previous} and {rel}"
            )
        collision_keys[collision_key] = rel
        expected[rel] = digest
    if not expected:
        raise ValueError("SHA256SUMS.txt must contain at least one checksum entry")
    return expected


def load_manifest(release_dir):
    path = Path(release_dir) / "release_manifest.json"
    if not path.exists():
        return None
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def load_manifest_for_validation(release_dir):
    try:
        return load_manifest(release_dir), []
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        return None, [f"Could not parse release_manifest.json: {exc}"]


def read_release_metadata(release_dir):
    errors = []
    release_dir = Path(release_dir)
    version = None
    metadata = None
    version_path = release_dir / "VERSION"
    metadata_path = release_dir / "RELEASE_METADATA.json"
    if version_path.is_symlink():
        errors.append("VERSION must not be a symbolic link")
    if metadata_path.is_symlink():
        errors.append("RELEASE_METADATA.json must not be a symbolic link")
    if errors:
        return None, None, errors
    try:
        version = version_path.read_text(encoding="utf-8").strip()
        if not version:
            errors.append("VERSION must be non-empty")
    except (OSError, UnicodeError) as exc:
        errors.append(f"Could not read VERSION: {exc}")
    try:
        metadata = json.loads(
            metadata_path.read_text(encoding="utf-8")
        )
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        errors.append(f"Could not parse RELEASE_METADATA.json: {exc}")
        return version, None, errors
    if not isinstance(metadata, dict):
        errors.append("RELEASE_METADATA.json top-level value must be an object")
        return version, None, errors
    if version is not None and metadata.get("version") != version:
        errors.append(
            "RELEASE_METADATA.json version mismatch: "
            f"expected {version!r} from VERSION, found {metadata.get('version')!r}"
        )
    try:
        canonical_zip_datetime(metadata.get("release_date"))
    except ValueError as exc:
        errors.append(f"RELEASE_METADATA.json {exc}")
    return version, metadata, errors


def validate_manifest_metadata(release_dir, manifest, manifest_files):
    errors = []
    if not isinstance(manifest, dict):
        return ["release_manifest.json top-level value must be an object"]
    if manifest.get("schema_version") != EXPECTED_MANIFEST_SCHEMA_VERSION:
        errors.append(
            "release_manifest.json schema_version mismatch: "
            f"expected {EXPECTED_MANIFEST_SCHEMA_VERSION!r}, "
            f"found {manifest.get('schema_version')!r}"
        )
    if manifest.get("bundle_name") != EXPECTED_BUNDLE_NAME:
        errors.append(
            "release_manifest.json bundle_name mismatch: "
            f"expected {EXPECTED_BUNDLE_NAME!r}, found {manifest.get('bundle_name')!r}"
        )
    expected_count = len(manifest_files)
    manifest_count = manifest.get("file_count_excluding_manifest_and_checksums")
    if not isinstance(manifest_count, int) or isinstance(manifest_count, bool):
        errors.append(
            "release_manifest.json file_count_excluding_manifest_and_checksums "
            "must be an integer"
        )
    elif manifest_count != expected_count:
        errors.append(
            "release_manifest.json file_count_excluding_manifest_and_checksums "
            f"mismatch: expected {expected_count}, found {manifest_count}"
        )
    if manifest.get("forbidden_artifact_names") != sorted(FORBIDDEN_NAMES):
        errors.append(
            "release_manifest.json forbidden_artifact_names mismatch: "
            "expected validator deny-list"
        )
    if manifest.get("builder") != "build_release_lite.py":
        errors.append("release_manifest.json builder must be 'build_release_lite.py'")
    if manifest.get("validation_command") != RELEASE_VALIDATION_COMMAND:
        errors.append("release_manifest.json validation_command is not canonical")

    version, metadata, metadata_errors = read_release_metadata(release_dir)
    errors.extend(metadata_errors)
    if version is not None and manifest.get("version") != version:
        errors.append(
            "release_manifest.json version mismatch: "
            f"expected {version!r} from VERSION, found {manifest.get('version')!r}"
        )
    if metadata is not None and manifest.get("release_date") != metadata.get("release_date"):
        errors.append(
            "release_manifest.json release_date mismatch: "
            f"expected {metadata.get('release_date')!r} from RELEASE_METADATA.json, "
            f"found {manifest.get('release_date')!r}"
        )
    return errors


def validate_release_dir(release_dir, max_file_bytes):
    errors = []
    warnings = []
    release_arg = Path(release_dir)
    if release_arg.is_symlink():
        return [f"Release directory must not be a symbolic link: {release_arg}"], warnings, None
    release_dir = release_arg.resolve()
    if not release_dir.exists():
        return [f"Release directory not found: {release_dir}"], warnings, None
    if not release_dir.is_dir():
        return [f"Release path is not a directory: {release_dir}"], warnings, None

    files = release_files(release_dir)
    file_collision_keys = {}
    for rel in sorted(files):
        try:
            normalized = normalize_release_path(rel)
        except ValueError as exc:
            errors.append(f"Release payload path {rel!r}: {exc}")
            continue
        collision_key = release_path_collision_key(normalized)
        previous = file_collision_keys.get(collision_key)
        if previous is not None:
            errors.append(
                "Release payload contains a case-insensitive path collision: "
                f"{previous} and {rel}"
            )
        else:
            file_collision_keys[collision_key] = rel
    manifest_path = release_dir / "release_manifest.json"
    manifest_present = manifest_path.is_file() and not manifest_path.is_symlink()
    if manifest_path.is_symlink():
        errors.append("release_manifest.json must not be a symbolic link")
        manifest, manifest_errors = None, []
    else:
        manifest, manifest_errors = load_manifest_for_validation(release_dir)
        errors.extend(manifest_errors)
    manifest_files = {}
    manifest_collision_keys = {}
    if not manifest_present:
        errors.append("release_manifest.json is required")
    elif not manifest_errors:
        if not isinstance(manifest, dict):
            errors.append("release_manifest.json top-level value must be an object")
        else:
            manifest_items = manifest.get("files")
            if not isinstance(manifest_items, list):
                errors.append("release_manifest.json files must be a list")
                manifest_items = []
            for idx, item in enumerate(manifest_items, start=1):
                if not isinstance(item, dict):
                    errors.append(f"Manifest file entry {idx} is not an object")
                    continue
                try:
                    rel = normalize_release_path(item.get("path"))
                except ValueError as exc:
                    errors.append(f"Manifest file entry {idx}: {exc}")
                    continue
                if rel in manifest_files:
                    errors.append(f"Duplicate release_manifest.json file entry: {rel}")
                    continue
                collision_key = release_path_collision_key(rel)
                previous = manifest_collision_keys.get(collision_key)
                if previous is not None:
                    errors.append(
                        "Case-insensitive release_manifest.json path collision: "
                        f"{previous} and {rel}"
                    )
                    continue
                manifest_collision_keys[collision_key] = rel
                manifest_files[rel] = item
            errors.extend(validate_manifest_metadata(release_dir, manifest, manifest_files))

    missing_required = sorted(path for path in REQUIRED_FALLBACK if path not in files)
    if missing_required:
        errors.append("Missing required files: " + ", ".join(missing_required))

    for rel, path in sorted(files.items()):
        name = Path(rel).name
        if path.is_symlink():
            errors.append(f"Symbolic links are not allowed in release payloads: {rel}")
            continue
        if name in FORBIDDEN_NAMES:
            errors.append(f"Forbidden training/full-artifact file present: {rel}")
        if any(part in TRANSIENT_NAMES for part in Path(rel).parts):
            errors.append(f"Transient cache path present: {rel}")
        if name.startswith(TRANSIENT_PREFIXES) or name.endswith(TRANSIENT_SUFFIXES):
            errors.append(f"Transient test/cache file present: {rel}")
        if path.stat().st_size > max_file_bytes:
            errors.append(f"File exceeds max size {max_file_bytes} bytes: {rel}")
        if is_text_release_path(rel) and b"\r" in path.read_bytes():
            errors.append(f"Text release file contains CR newline bytes: {rel}")

    checksum_path = release_dir / "SHA256SUMS.txt"
    expected_hashes = None
    if checksum_path.is_symlink():
        errors.append("SHA256SUMS.txt must not be a symbolic link")
    elif not checksum_path.is_file():
        errors.append("Missing SHA256SUMS.txt")
    else:
        try:
            expected_hashes = parse_sha256sums(checksum_path)
        except (OSError, UnicodeError, ValueError) as exc:
            errors.append(f"Could not parse SHA256SUMS.txt: {exc}")
    if expected_hashes is not None:
        listed = set(expected_hashes)
        actual_for_hash = set(files) - {"SHA256SUMS.txt"}
        extra_listed = sorted(listed - actual_for_hash)
        unlisted = sorted(actual_for_hash - listed)
        if extra_listed:
            errors.append("SHA256SUMS lists missing files: " + ", ".join(extra_listed))
        if unlisted:
            errors.append("Files missing from SHA256SUMS: " + ", ".join(unlisted))
        for rel, expected in sorted(expected_hashes.items()):
            path = release_dir / rel
            if path.is_file() and sha256_file(path) != expected:
                errors.append(f"SHA256 mismatch for {rel}")

    if manifest_present:
        manifest_actual = set(files) - {"SHA256SUMS.txt", "release_manifest.json"}
        missing_from_manifest = sorted(manifest_actual - set(manifest_files))
        stale_manifest = sorted(set(manifest_files) - manifest_actual)
        if missing_from_manifest:
            errors.append(
                "Files missing from release_manifest.json: " + ", ".join(missing_from_manifest)
            )
        if stale_manifest:
            errors.append("Manifest lists missing files: " + ", ".join(stale_manifest))
        for rel, item in sorted(manifest_files.items()):
            path = release_dir / rel
            if not path.is_file():
                continue
            size = path.stat().st_size
            digest = sha256_file(path)
            manifest_bytes = item.get("bytes")
            if not isinstance(manifest_bytes, int) or isinstance(manifest_bytes, bool):
                errors.append(f"Manifest byte size for {rel} must be an integer")
            elif manifest_bytes != size:
                errors.append(f"Manifest byte size mismatch for {rel}")
            manifest_digest = item.get("sha256")
            if (
                not isinstance(manifest_digest, str)
                or len(manifest_digest) != 64
                or any(ch not in "0123456789abcdef" for ch in manifest_digest)
            ):
                errors.append(f"Manifest SHA256 for {rel} is malformed")
            elif manifest_digest != digest:
                errors.append(f"Manifest SHA256 mismatch for {rel}")

    summary = {
        "release_dir": str(release_dir),
        "file_count": len(files),
        "checksum_count": len(expected_hashes or {}),
        "has_manifest": manifest_present,
        "total_bytes": sum(path.stat().st_size for path in files.values() if path.is_file()),
    }
    return errors, warnings, summary


def validate_source_parity(release_dir, source_root, manifest):
    errors = []
    release_dir = Path(release_dir).resolve()
    source_root = Path(source_root).resolve()
    if not source_root.exists() or not source_root.is_dir():
        return [f"Source root not found or not a directory: {source_root}"]
    if not isinstance(manifest, dict):
        return ["Cannot validate source parity without a valid release_manifest.json"]
    manifest_items = manifest.get("files")
    if not isinstance(manifest_items, list):
        return ["release_manifest.json files must be a list"]
    for idx, item in enumerate(manifest_items, start=1):
        if not isinstance(item, dict):
            errors.append(f"Manifest file entry {idx} is not an object")
            continue
        try:
            rel = normalize_release_path(item.get("path"))
        except ValueError as exc:
            errors.append(f"Manifest file entry {idx}: {exc}")
            continue
        src = source_root / rel
        dst = release_dir / rel
        if not src.is_file():
            errors.append(f"Source file missing for release payload: {rel}")
        elif not dst.is_file():
            errors.append(f"Release file missing for source parity: {rel}")
        elif sha256_file(src) != sha256_file(dst):
            errors.append(f"Release file is stale relative to source: {rel}")
    return errors


def validate_zip(
    zip_path,
    release_dir,
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
        return [f"Zip archive must not be a symbolic link: {zip_arg}"], None
    zip_path = zip_arg.resolve()
    release_dir = Path(release_dir).resolve()
    if not zip_path.is_file():
        return [f"Zip archive not found: {zip_path}"], None
    zip_bytes = zip_path.stat().st_size
    if zip_bytes > max_archive_bytes:
        return [f"Zip archive is {zip_bytes} bytes; limit is {max_archive_bytes}"], {
            "zip_path": str(zip_path),
            "zip_entries": 0,
            "zip_bytes": zip_bytes,
            "zip_sha256": None,
        }
    release = release_files(release_dir)
    try:
        zf = zipfile.ZipFile(zip_path)
    except (OSError, zipfile.BadZipFile) as exc:
        return [f"Could not read zip archive: {exc}"], None
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
        _, metadata, metadata_errors = read_release_metadata(release_dir)
        errors.extend(metadata_errors)
        expected_datetime = None
        if metadata is not None:
            try:
                expected_datetime = canonical_zip_datetime(metadata.get("release_date"))
            except ValueError:
                pass
        if expected_datetime is not None:
            errors.extend(canonical_zip_errors(zf, infos, expected_datetime))

        zip_records = []
        zip_names = []
        for info in infos:
            try:
                rel = normalize_release_path(info.filename)
            except ValueError:
                continue
            zip_names.append(rel)
            zip_records.append((rel, info))
        release_names = sorted(release)
        if zip_names != release_names:
            missing = sorted(set(release_names) - set(zip_names))
            extra = sorted(set(zip_names) - set(release_names))
            if missing:
                errors.append("Zip is missing files: " + ", ".join(missing))
            if extra:
                errors.append("Zip has extra files: " + ", ".join(extra))
            if not missing and not extra:
                errors.append("Zip member order differs from the release directory")

        if not safety_errors:
            archive_readable = True
            try:
                bad = zf.testzip()
            except (OSError, RuntimeError, zipfile.BadZipFile) as exc:
                errors.append(f"Could not test zip archive: {exc}")
                bad = None
                archive_readable = False
            if bad is not None:
                errors.append(f"Zip archive is corrupt at {bad}")
                archive_readable = False
            if archive_readable:
                for rel, info in zip_records:
                    if rel not in release:
                        continue
                    try:
                        member_digest = sha256_zip_member(zf, info)
                    except (OSError, RuntimeError, zipfile.BadZipFile) as exc:
                        errors.append(f"Could not read zip member {info.filename}: {exc}")
                        continue
                    if member_digest != sha256_file(release[rel]):
                        errors.append(f"Zip content differs from release dir: {info.filename}")
    summary = {
        "zip_path": str(zip_path),
        "zip_entries": len(infos),
        "zip_bytes": zip_bytes,
        "zip_sha256": sha256_file(zip_path),
    }
    return errors, summary


def _relative_artifact_path(target, base, label, errors):
    try:
        return Path(target).resolve().relative_to(Path(base).resolve()).as_posix()
    except ValueError:
        errors.append(f"{label} must be located below the RELEASE_ARTIFACTS.json directory")
        return None


def validate_artifacts(path, release_summary, zip_summary):
    errors = []
    path_arg = Path(path)
    if path_arg.is_symlink():
        return [f"RELEASE_ARTIFACTS.json must not be a symbolic link: {path_arg}"]
    path = path_arg.resolve()
    if not path.is_file():
        return [f"Release artifact metadata not found: {path}"]
    try:
        artifact = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        return [f"Could not read release artifact metadata: {exc}"]
    if not isinstance(artifact, dict):
        return ["RELEASE_ARTIFACTS.json top-level value must be an object"]
    if not release_summary or not zip_summary:
        return ["Release and zip summaries are required to validate RELEASE_ARTIFACTS.json"]

    release_dir = Path(release_summary["release_dir"])
    version, metadata, metadata_errors = read_release_metadata(release_dir)
    errors.extend(metadata_errors)
    release_rel = _relative_artifact_path(release_dir, path.parent, "release_dir", errors)
    zip_rel = _relative_artifact_path(zip_summary["zip_path"], path.parent, "zip_path", errors)
    zip_digest = zip_summary["zip_sha256"]
    expected = {
        "schema_version": RELEASE_SCHEMA_VERSION,
        "version": version,
        "release_date": metadata.get("release_date") if metadata else None,
        "release_dir": release_rel,
        "release_file_count": release_summary["file_count"],
        "release_total_bytes": release_summary["total_bytes"],
        "zip_path": zip_rel,
        "zip_entries": zip_summary["zip_entries"],
        "zip_bytes": zip_summary["zip_bytes"],
        "zip_sha256": zip_digest,
        "validation_command": RELEASE_VALIDATION_COMMAND,
        "zip_acceptance_command": (
            f"{ZIP_ACCEPTANCE_COMMAND} --expected-sha256 {zip_digest}"
        ),
    }
    required_keys = set(expected)
    missing = sorted(required_keys - set(artifact))
    if missing:
        errors.append("RELEASE_ARTIFACTS.json missing keys: " + ", ".join(missing))
    integer_keys = {
        "release_file_count",
        "release_total_bytes",
        "zip_entries",
        "zip_bytes",
    }
    for key, expected_value in expected.items():
        actual = artifact.get(key)
        if key in integer_keys and (
            not isinstance(actual, int) or isinstance(actual, bool) or actual < 0
        ):
            errors.append(f"RELEASE_ARTIFACTS.json {key} must be a non-negative integer")
            continue
        if actual != expected_value:
            errors.append(
                f"RELEASE_ARTIFACTS.json {key} mismatch: "
                f"expected {expected_value!r}, found {actual!r}"
            )
    return errors


def run_smoke(release_dir, timeout_seconds):
    env = dict(os.environ)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env["PYTHONUTF8"] = "1"
    try:
        result = subprocess.run(
            [sys.executable, "run_smoke_tests.py"],
            cwd=release_dir,
            text=True,
            env=env,
            timeout=timeout_seconds,
        )
    except subprocess.TimeoutExpired:
        return 124
    return result.returncode


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release-dir", default="release-lite")
    parser.add_argument("--zip", dest="zip_path", default=RELEASE_ZIP_NAME)
    parser.add_argument("--no-zip", action="store_true")
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--max-file-bytes", type=int, default=5_000_000)
    parser.add_argument("--max-zip-entries", type=int, default=DEFAULT_MAX_ZIP_ENTRIES)
    parser.add_argument(
        "--max-zip-archive-bytes", type=int, default=DEFAULT_MAX_ZIP_ARCHIVE_BYTES
    )
    parser.add_argument(
        "--max-zip-member-bytes", type=int, default=DEFAULT_MAX_ZIP_MEMBER_BYTES
    )
    parser.add_argument(
        "--max-zip-total-bytes", type=int, default=DEFAULT_MAX_ZIP_TOTAL_BYTES
    )
    parser.add_argument(
        "--max-zip-compression-ratio",
        type=float,
        default=DEFAULT_MAX_ZIP_COMPRESSION_RATIO,
    )
    parser.add_argument("--source-root")
    parser.add_argument("--artifacts")
    parser.add_argument("--timeout-seconds", type=int, default=300)
    args = parser.parse_args(argv)

    release_dir = Path(args.release_dir)
    errors, warnings, release_summary = validate_release_dir(
        release_dir, args.max_file_bytes
    )
    manifest = None
    if release_dir.exists() and not errors:
        manifest, manifest_errors = load_manifest_for_validation(release_dir.resolve())
        errors.extend(manifest_errors)
    if args.source_root and not errors:
        errors.extend(validate_source_parity(release_dir, Path(args.source_root), manifest))
    zip_summary = None
    if not args.no_zip:
        zip_errors, zip_summary = validate_zip(
            Path(args.zip_path),
            release_dir,
            max_archive_bytes=args.max_zip_archive_bytes,
            max_entries=args.max_zip_entries,
            max_member_bytes=args.max_zip_member_bytes,
            max_total_bytes=args.max_zip_total_bytes,
            max_compression_ratio=args.max_zip_compression_ratio,
        )
        errors.extend(zip_errors)
    if args.artifacts and not errors:
        errors.extend(validate_artifacts(Path(args.artifacts), release_summary, zip_summary))
    if args.smoke and not errors:
        smoke_code = run_smoke(release_dir, args.timeout_seconds)
        if smoke_code != 0:
            errors.append(f"Smoke test failed with exit code {smoke_code}")

    for warning in warnings:
        print(f"[validate] WARNING: {warning}", file=sys.stderr)
    for error in errors:
        print(f"[validate] ERROR: {error}", file=sys.stderr)
    if release_summary:
        print(f"[validate] release files: {release_summary['file_count']}")
        print(f"[validate] checksum entries: {release_summary['checksum_count']}")
        print(f"[validate] release bytes: {release_summary['total_bytes']}")
    if zip_summary:
        print(f"[validate] zip entries: {zip_summary['zip_entries']}")
        print(f"[validate] zip bytes: {zip_summary['zip_bytes']}")
    if errors:
        print("[validate] FAIL", file=sys.stderr)
        return 1
    print("[validate] PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
