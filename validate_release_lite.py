#!/usr/bin/env python3
"""Validate the lightweight release bundle and optional zip archive."""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import zipfile
from pathlib import Path, PurePosixPath


REQUIRED_FALLBACK = {
    "audit_cli_entrypoints.py",
    "audit_lightweight_dependencies.py",
    "audit_release_docs.py",
    "calibrate_threshold.py",
    "check_environment.py",
    "DATA_DICTIONARY.md",
    "deployable_lr_weights.npz",
    "EXECUTIVE_SUMMARY.md",
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

FORBIDDEN_NAMES = {
    "deployable_pipeline.pkl",
    "feature_selection.pkl",
    "final_model_results.pkl",
    "gene_id_to_name.pkl",
    "groups_full.pkl",
    "model_lr.pkl",
    "model_rf.pkl",
    "model_xgb.pkl",
    "projects_full.pkl",
    "sample_metadata.pkl",
    "selected_files.csv",
    "train_classifier.py",
    "train_idx.npy",
    "test_idx.npy",
    "X_full_filtered.pkl",
    "y_full.pkl",
}

TRANSIENT_NAMES = {"__pycache__"}
TRANSIENT_PREFIXES = ("_smoke_", "_acceptance_")
TRANSIENT_SUFFIXES = (".pyc",)


def sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalize_release_path(rel):
    if not isinstance(rel, str) or not rel.strip():
        raise ValueError(f"Invalid empty release path: {rel!r}")
    if "\\" in rel:
        raise ValueError(f"Release paths must use forward slashes: {rel!r}")
    path = PurePosixPath(rel)
    if path.is_absolute() or rel.startswith("/"):
        raise ValueError(f"Release path must be relative: {rel!r}")
    if any(part in {"", ".", ".."} for part in path.parts):
        raise ValueError(f"Release path contains unsafe component: {rel!r}")
    return path.as_posix()


def release_files(release_dir):
    return {
        path.relative_to(release_dir).as_posix(): path
        for path in release_dir.rglob("*")
        if path.is_file()
    }


def parse_sha256sums(path):
    expected = {}
    for lineno, line in enumerate(path.read_text(encoding="ascii").splitlines(), start=1):
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
        expected[rel] = digest
    return expected


def load_manifest(release_dir):
    path = release_dir / "release_manifest.json"
    if not path.exists():
        return None
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_release_dir(release_dir, max_file_bytes):
    errors = []
    warnings = []
    release_dir = release_dir.resolve()
    if not release_dir.exists():
        return [f"Release directory not found: {release_dir}"], warnings, None
    if not release_dir.is_dir():
        return [f"Release path is not a directory: {release_dir}"], warnings, None

    files = release_files(release_dir)
    manifest = load_manifest(release_dir)
    required = set(REQUIRED_FALLBACK)
    manifest_files = {}
    if manifest:
        seen_manifest_paths = set()
        for idx, item in enumerate(manifest.get("files", []), start=1):
            if not isinstance(item, dict):
                errors.append(f"Manifest file entry {idx} is not an object")
                continue
            try:
                rel = normalize_release_path(item.get("path"))
            except ValueError as exc:
                errors.append(f"Manifest file entry {idx}: {exc}")
                continue
            if rel in seen_manifest_paths:
                errors.append(f"Duplicate release_manifest.json file entry: {rel}")
                continue
            seen_manifest_paths.add(rel)
            manifest_files[rel] = item
        required.update(manifest_files)
        required.update({"release_manifest.json", "SHA256SUMS.txt"})
    missing_required = sorted(path for path in required if path not in files)
    if missing_required:
        errors.append("Missing required files: " + ", ".join(missing_required))

    for rel, path in sorted(files.items()):
        name = Path(rel).name
        if name in FORBIDDEN_NAMES:
            errors.append(f"Forbidden training/full-artifact file present: {rel}")
        if any(part in TRANSIENT_NAMES for part in Path(rel).parts):
            errors.append(f"Transient cache path present: {rel}")
        if name.startswith(TRANSIENT_PREFIXES) or name.endswith(TRANSIENT_SUFFIXES):
            errors.append(f"Transient test/cache file present: {rel}")
        if path.stat().st_size > max_file_bytes:
            errors.append(f"File exceeds max size {max_file_bytes} bytes: {rel}")

    checksum_path = release_dir / "SHA256SUMS.txt"
    if not checksum_path.exists():
        errors.append("Missing SHA256SUMS.txt")
        expected_hashes = {}
    else:
        try:
            expected_hashes = parse_sha256sums(checksum_path)
        except ValueError as exc:
            errors.append(str(exc))
            expected_hashes = {}

    if expected_hashes:
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
            if path.exists():
                actual = sha256_file(path)
                if actual != expected:
                    errors.append(f"SHA256 mismatch for {rel}")

    if manifest:
        manifest_actual = set(files) - {"SHA256SUMS.txt", "release_manifest.json"}
        missing_from_manifest = sorted(manifest_actual - set(manifest_files))
        stale_manifest = sorted(set(manifest_files) - manifest_actual)
        if missing_from_manifest:
            errors.append("Files missing from release_manifest.json: "
                          + ", ".join(missing_from_manifest))
        if stale_manifest:
            errors.append("Manifest lists missing files: " + ", ".join(stale_manifest))
        for rel, item in sorted(manifest_files.items()):
            path = release_dir / rel
            if path.exists():
                size = path.stat().st_size
                digest = sha256_file(path)
                if int(item.get("bytes", -1)) != size:
                    errors.append(f"Manifest byte size mismatch for {rel}")
                if item.get("sha256") != digest:
                    errors.append(f"Manifest SHA256 mismatch for {rel}")
    else:
        warnings.append("release_manifest.json not found; using fallback required-file checks")

    summary = {
        "release_dir": str(release_dir),
        "file_count": len(files),
        "checksum_count": len(expected_hashes),
        "has_manifest": manifest is not None,
        "total_bytes": sum(path.stat().st_size for path in files.values()),
    }
    return errors, warnings, summary


def validate_source_parity(release_dir, source_root, manifest):
    errors = []
    source_root = source_root.resolve()
    if not source_root.exists() or not source_root.is_dir():
        return [f"Source root not found or not a directory: {source_root}"]
    if manifest is None:
        return ["Cannot validate source parity without release_manifest.json"]

    for idx, item in enumerate(manifest.get("files", []), start=1):
        try:
            rel = normalize_release_path(item.get("path"))
        except ValueError as exc:
            errors.append(f"Manifest file entry {idx}: {exc}")
            continue
        src = source_root / rel
        dst = release_dir / rel
        if not src.exists():
            errors.append(f"Source file missing for release payload: {rel}")
            continue
        if not dst.exists():
            errors.append(f"Release file missing for source parity: {rel}")
            continue
        if sha256_file(src) != sha256_file(dst):
            errors.append(f"Release file is stale relative to source: {rel}")
    return errors


def validate_zip(zip_path, release_dir):
    errors = []
    zip_path = zip_path.resolve()
    if not zip_path.exists():
        return [f"Zip archive not found: {zip_path}"], None
    release = release_files(release_dir.resolve())
    with zipfile.ZipFile(zip_path) as zf:
        bad = zf.testzip()
        if bad is not None:
            errors.append(f"Zip archive is corrupt at {bad}")
        infos = [info for info in zf.infolist() if not info.is_dir()]
        zip_names = {info.filename for info in infos}
        release_names = set(release)
        if zip_names != release_names:
            missing = sorted(release_names - zip_names)
            extra = sorted(zip_names - release_names)
            if missing:
                errors.append("Zip is missing files: " + ", ".join(missing))
            if extra:
                errors.append("Zip has extra files: " + ", ".join(extra))
        for info in infos:
            if info.filename in release:
                with zf.open(info) as handle:
                    digest = sha256_bytes(handle.read())
                if digest != sha256_file(release[info.filename]):
                    errors.append(f"Zip content differs from release dir: {info.filename}")
    summary = {
        "zip_path": str(zip_path),
        "zip_entries": len(zip_names),
        "zip_bytes": zip_path.stat().st_size,
    }
    return errors, summary


def validate_artifacts(path, release_summary, zip_summary):
    errors = []
    path = path.resolve()
    if not path.exists():
        return [f"Release artifact metadata not found: {path}"]
    try:
        artifact = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return [f"Could not read release artifact metadata: {exc}"]

    if release_summary:
        expected_release = {
            "release_file_count": release_summary["file_count"],
            "release_total_bytes": release_summary["total_bytes"],
        }
        for key, expected in expected_release.items():
            if artifact.get(key) != expected:
                errors.append(
                    f"RELEASE_ARTIFACTS.json {key} mismatch: "
                    f"expected {expected}, found {artifact.get(key)}"
                )
    if zip_summary:
        zip_path = Path(zip_summary["zip_path"])
        expected_zip = {
            "zip_entries": zip_summary["zip_entries"],
            "zip_bytes": zip_summary["zip_bytes"],
            "zip_sha256": sha256_file(zip_path),
        }
        for key, expected in expected_zip.items():
            if artifact.get(key) != expected:
                errors.append(
                    f"RELEASE_ARTIFACTS.json {key} mismatch: "
                    f"expected {expected}, found {artifact.get(key)}"
                )
    return errors


def run_smoke(release_dir, timeout_seconds):
    env = dict(os.environ)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
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
    parser = argparse.ArgumentParser(description="Validate release-lite and optional zip.")
    parser.add_argument("--release-dir", default="release-lite")
    parser.add_argument("--zip", dest="zip_path",
                        default="tcga-tumor-normal-release-lite.zip")
    parser.add_argument("--no-zip", action="store_true")
    parser.add_argument("--smoke", action="store_true",
                        help="run run_smoke_tests.py inside the release directory")
    parser.add_argument("--max-file-bytes", type=int, default=5_000_000)
    parser.add_argument("--source-root",
                        help="optional full deliverables root; fail if payload files differ")
    parser.add_argument("--artifacts",
                        help="optional RELEASE_ARTIFACTS.json sidecar to validate")
    parser.add_argument("--timeout-seconds", type=int, default=300,
                        help="subprocess timeout for --smoke (default: 300)")
    args = parser.parse_args(argv)

    release_dir = Path(args.release_dir)
    errors, warnings, release_summary = validate_release_dir(release_dir, args.max_file_bytes)
    manifest = load_manifest(release_dir.resolve()) if release_dir.exists() else None
    if args.source_root and not errors:
        errors.extend(validate_source_parity(release_dir.resolve(), Path(args.source_root), manifest))
    zip_summary = None
    if not args.no_zip:
        zip_errors, zip_summary = validate_zip(Path(args.zip_path), release_dir)
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
