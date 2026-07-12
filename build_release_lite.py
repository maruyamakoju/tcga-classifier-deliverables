#!/usr/bin/env python3
"""Build the lightweight release bundle without exposing partial outputs."""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

from release_tools.common import (
    CANONICAL_ZIP_COMPRESSION,
    FORBIDDEN_NAMES,
    RELEASE_BUNDLE_NAME,
    RELEASE_FILES,
    RELEASE_SCHEMA_VERSION,
    RELEASE_VALIDATION_COMMAND,
    RELEASE_ZIP_NAME,
    ZIP_ACCEPTANCE_COMMAND,
    canonical_zip_datetime,
    canonical_zip_info,
    sha256_file,
)


ROOT = Path(__file__).resolve().parent
RELEASE_DIR = ROOT / "release-lite"
ZIP_PATH = ROOT / RELEASE_ZIP_NAME
ARTIFACTS_PATH = ROOT / "RELEASE_ARTIFACTS.json"

BINARY_SUFFIXES = {".npy", ".npz", ".pkl", ".png", ".zip"}
TEXT_NAMES = {"LICENSE", "VERSION"}


def assert_inside_root(path, root=None):
    resolved = Path(path).resolve()
    resolved_root = Path(ROOT if root is None else root).resolve()
    if resolved != resolved_root and resolved_root not in resolved.parents:
        raise RuntimeError(f"Refusing to operate outside project root: {resolved}")
    return resolved


def write_text_lf(path, text, encoding="utf-8"):
    Path(path).write_bytes(text.encode(encoding))


def is_text_release_file(rel):
    path = Path(rel)
    if path.name in TEXT_NAMES:
        return True
    return path.suffix.lower() not in BINARY_SUFFIXES


def sorted_release_paths(release_dir):
    release_dir = Path(release_dir)
    return sorted(
        (path for path in release_dir.rglob("*") if path.is_file()),
        key=lambda path: path.relative_to(release_dir).as_posix(),
    )


def load_release_metadata(source_root):
    source_root = Path(source_root)
    metadata_path = source_root / "RELEASE_METADATA.json"
    version_path = source_root / "VERSION"
    try:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"Could not read RELEASE_METADATA.json: {exc}") from exc
    if not isinstance(metadata, dict):
        raise RuntimeError("RELEASE_METADATA.json top-level value must be an object")
    try:
        version = version_path.read_text(encoding="utf-8").strip()
    except OSError as exc:
        raise RuntimeError(f"Could not read VERSION: {exc}") from exc
    if not version:
        raise RuntimeError("VERSION must be non-empty")
    if metadata.get("version") != version:
        raise RuntimeError("VERSION and RELEASE_METADATA.json version disagree")
    canonical_zip_datetime(metadata.get("release_date"))
    return metadata


def write_release_file(src, dst, rel):
    if is_text_release_file(rel):
        data = src.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
        dst.write_bytes(data)
    else:
        shutil.copy2(src, dst)


def copy_release_files(source_root, release_dir):
    source_root = Path(source_root)
    release_dir = Path(release_dir)
    missing = []
    for rel in RELEASE_FILES:
        src = source_root / rel
        if src.is_symlink():
            raise RuntimeError(f"Release source file must not be a symbolic link: {rel}")
        if not src.is_file():
            missing.append(rel)
            continue
        dst = release_dir / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        write_release_file(src, dst, rel)
    if missing:
        raise FileNotFoundError("Missing release source files:\n" + "\n".join(missing))


def manifest_file_records(release_dir):
    records = []
    for path in sorted_release_paths(release_dir):
        if path.name in {"SHA256SUMS.txt", "release_manifest.json"}:
            continue
        records.append(
            {
                "path": path.relative_to(release_dir).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    return records


def write_manifest(release_dir, release_metadata):
    release_dir = Path(release_dir)
    records = manifest_file_records(release_dir)
    manifest = {
        "schema_version": RELEASE_SCHEMA_VERSION,
        "bundle_name": RELEASE_BUNDLE_NAME,
        "version": release_metadata["version"],
        "release_date": release_metadata["release_date"],
        "intended_input": (
            "GDC STAR-Counts-style log2(TPM+1), rows=samples, columns=Ensembl genes."
        ),
        "builder": "build_release_lite.py",
        "validation_command": RELEASE_VALIDATION_COMMAND,
        "file_count_excluding_manifest_and_checksums": len(records),
        "forbidden_artifact_names": sorted(FORBIDDEN_NAMES),
        "files": records,
    }
    write_text_lf(
        release_dir / "release_manifest.json",
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
    )
    return len(records)


def write_checksums(release_dir):
    release_dir = Path(release_dir)
    rows = []
    for path in sorted_release_paths(release_dir):
        if path.name == "SHA256SUMS.txt":
            continue
        rel = path.relative_to(release_dir).as_posix()
        rows.append(f"{sha256_file(path)}  {rel}")
    if not rows:
        raise RuntimeError("Refusing to write an empty SHA256SUMS.txt")
    write_text_lf(release_dir / "SHA256SUMS.txt", "\n".join(rows) + "\n", encoding="ascii")
    return len(rows)


def write_zip(release_dir, zip_path, release_date):
    release_dir = Path(release_dir)
    zip_path = Path(zip_path)
    fixed_datetime = canonical_zip_datetime(release_date)
    with zipfile.ZipFile(
        zip_path,
        "x",
        compression=CANONICAL_ZIP_COMPRESSION,
        allowZip64=False,
    ) as zf:
        zf.comment = b""
        for path in sorted_release_paths(release_dir):
            rel = path.relative_to(release_dir).as_posix()
            info = canonical_zip_info(rel, fixed_datetime)
            zf.writestr(info, path.read_bytes(), compress_type=CANONICAL_ZIP_COMPRESSION)
    with zipfile.ZipFile(zip_path) as zf:
        bad = zf.testzip()
        if bad is not None:
            raise RuntimeError(f"Zip validation failed at {bad}")
        return len(zf.infolist())


def write_artifact_metadata(
    release_dir,
    zip_path,
    artifacts_path,
    release_metadata,
    zip_entries,
):
    release_dir = Path(release_dir)
    zip_path = Path(zip_path)
    artifacts_path = Path(artifacts_path)
    base = artifacts_path.parent.resolve()
    files = sorted_release_paths(release_dir)
    zip_digest = sha256_file(zip_path)
    artifact = {
        "schema_version": RELEASE_SCHEMA_VERSION,
        "version": release_metadata["version"],
        "release_date": release_metadata["release_date"],
        "release_dir": release_dir.resolve().relative_to(base).as_posix(),
        "release_file_count": len(files),
        "release_total_bytes": sum(path.stat().st_size for path in files),
        "zip_path": zip_path.resolve().relative_to(base).as_posix(),
        "zip_entries": zip_entries,
        "zip_bytes": zip_path.stat().st_size,
        "zip_sha256": zip_digest,
        "validation_command": RELEASE_VALIDATION_COMMAND,
        "zip_acceptance_command": f"{ZIP_ACCEPTANCE_COMMAND} --expected-sha256 {zip_digest}",
    }
    write_text_lf(artifacts_path, json.dumps(artifact, indent=2, sort_keys=True) + "\n")


def cleanup_transient_files(release_dir):
    release_dir = Path(release_dir).resolve()
    patterns = ["_smoke_*", "_acceptance_*", "*.pyc"]
    for pattern in patterns:
        for path in release_dir.rglob(pattern):
            if path.is_file():
                path.unlink()
    for path in sorted(release_dir.rglob("__pycache__"), reverse=True):
        if path.is_dir():
            shutil.rmtree(path)


def _run_checked(cmd, cwd, timeout_seconds):
    env = dict(os.environ)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env["PYTHONUTF8"] = "1"
    try:
        result = subprocess.run(cmd, cwd=cwd, text=True, env=env, timeout=timeout_seconds)
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(
            f"Command timed out after {timeout_seconds}s: {' '.join(map(str, cmd))}"
        ) from exc
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed with exit code {result.returncode}: {' '.join(map(str, cmd))}"
        )


def run_smoke_test(release_dir, timeout_seconds):
    for script in ["run_smoke_tests.py", "run_safety_tests.py"]:
        _run_checked([sys.executable, script], release_dir, timeout_seconds)


def run_release_validation(
    release_dir,
    zip_path,
    artifacts_path,
    source_root,
    timeout_seconds,
):
    _run_checked(
        [
            sys.executable,
            str(Path(source_root) / "validate_release_lite.py"),
            "--release-dir",
            str(release_dir),
            "--zip",
            str(zip_path),
            "--source-root",
            str(source_root),
            "--artifacts",
            str(artifacts_path),
        ],
        source_root,
        timeout_seconds,
    )


def run_zip_validation(zip_path, source_root, timeout_seconds):
    _run_checked(
        [
            sys.executable,
            str(Path(source_root) / "validate_zip_bundle.py"),
            str(zip_path),
            "--expected-sha256",
            sha256_file(zip_path),
            "--skip-acceptance",
            "--timeout-seconds",
            str(timeout_seconds),
        ],
        source_root,
        timeout_seconds,
    )


def build_stage(stage_root, source_root=ROOT, run_smoke=False, timeout_seconds=300):
    stage_root = Path(stage_root)
    source_root = Path(source_root).resolve()
    release_dir = stage_root / "release-lite"
    zip_path = stage_root / RELEASE_ZIP_NAME
    artifacts_path = stage_root / "RELEASE_ARTIFACTS.json"
    release_dir.mkdir(parents=True)
    metadata = load_release_metadata(source_root)
    copy_release_files(source_root, release_dir)
    if run_smoke:
        run_smoke_test(release_dir, timeout_seconds)
    cleanup_transient_files(release_dir)
    manifest_count = write_manifest(release_dir, metadata)
    checksum_count = write_checksums(release_dir)
    zip_entries = write_zip(release_dir, zip_path, metadata["release_date"])
    write_artifact_metadata(
        release_dir,
        zip_path,
        artifacts_path,
        metadata,
        zip_entries,
    )
    run_release_validation(
        release_dir,
        zip_path,
        artifacts_path,
        source_root,
        timeout_seconds,
    )
    run_zip_validation(zip_path, source_root, timeout_seconds)
    return {
        "release_dir": release_dir,
        "zip_path": zip_path,
        "artifacts_path": artifacts_path,
        "manifest_count": manifest_count,
        "checksum_count": checksum_count,
        "zip_entries": zip_entries,
    }


def directory_hashes(root):
    root = Path(root)
    if not root.is_dir():
        return None
    return {
        path.relative_to(root).as_posix(): sha256_file(path)
        for path in sorted_release_paths(root)
    }


def compare_outputs(staged):
    differences = []
    expected_tree = directory_hashes(staged["release_dir"])
    current_tree = directory_hashes(RELEASE_DIR)
    if current_tree != expected_tree:
        current_names = set(current_tree or {})
        expected_names = set(expected_tree or {})
        for rel in sorted(expected_names - current_names):
            differences.append(f"release-lite missing: {rel}")
        for rel in sorted(current_names - expected_names):
            differences.append(f"release-lite extra: {rel}")
        for rel in sorted(current_names & expected_names):
            if current_tree[rel] != expected_tree[rel]:
                differences.append(f"release-lite differs: {rel}")
        if current_tree is None:
            differences.append("release-lite directory is missing")
    for label, current, candidate in [
        ("release ZIP", ZIP_PATH, staged["zip_path"]),
        ("RELEASE_ARTIFACTS.json", ARTIFACTS_PATH, staged["artifacts_path"]),
    ]:
        if not current.is_file():
            differences.append(f"{label} is missing")
        elif sha256_file(current) != sha256_file(candidate):
            differences.append(f"{label} differs")
    return differences


def remove_path(path):
    path = Path(path)
    if path.is_dir() and not path.is_symlink():
        shutil.rmtree(path)
    elif path.exists() or path.is_symlink():
        path.unlink()


def atomic_replace_outputs(staged, backup_root):
    """Install all staged outputs as a rollback-capable transaction."""
    backup_root = Path(backup_root)
    pairs = [
        (Path(staged["release_dir"]), assert_inside_root(RELEASE_DIR)),
        (Path(staged["zip_path"]), assert_inside_root(ZIP_PATH)),
        (Path(staged["artifacts_path"]), assert_inside_root(ARTIFACTS_PATH)),
    ]
    backups = []
    installed = []
    try:
        for index, (candidate, target) in enumerate(pairs):
            backup = backup_root / f"backup-{index}-{target.name}"
            if target.exists() or target.is_symlink():
                os.replace(target, backup)
                backups.append((target, backup))
            os.replace(candidate, target)
            installed.append(target)
    except Exception as exc:
        rollback_errors = []
        for target in reversed(installed):
            try:
                remove_path(target)
            except OSError as rollback_exc:
                rollback_errors.append(str(rollback_exc))
        for target, backup in reversed(backups):
            try:
                if backup.exists() or backup.is_symlink():
                    os.replace(backup, target)
            except OSError as rollback_exc:
                rollback_errors.append(str(rollback_exc))
        detail = f"; rollback errors: {'; '.join(rollback_errors)}" if rollback_errors else ""
        raise RuntimeError(f"Could not install release outputs atomically{detail}") from exc


def print_build_summary(staged, prefix):
    print(f"[release] files in manifest: {staged['manifest_count']}")
    print(f"[release] files with checksums: {staged['checksum_count']}")
    print(f"[release] zip entries: {staged['zip_entries']}")
    print(f"[release] zip bytes: {Path(staged['zip_path']).stat().st_size}")
    print(f"[release] zip sha256: {sha256_file(staged['zip_path'])}")
    print(f"[release] {prefix}")


def build(run_smoke=False, timeout_seconds=300, check=False):
    # Stage on the project filesystem so os.replace never crosses volumes.
    with tempfile.TemporaryDirectory(prefix=".release-stage-", dir=ROOT) as temp:
        stage_root = Path(temp)
        staged = build_stage(
            stage_root,
            source_root=ROOT,
            run_smoke=run_smoke,
            timeout_seconds=timeout_seconds,
        )
        if check:
            differences = compare_outputs(staged)
            if differences:
                for difference in differences:
                    print(f"[release] CHECK FAIL: {difference}", file=sys.stderr)
                return False
            print_build_summary(staged, "CHECK PASS: committed outputs are reproducible")
            return True
        # Print before moving the staged files; their paths cease to exist after
        # successful installation.
        print_build_summary(staged, "staging validation passed")
        atomic_replace_outputs(staged, stage_root)
    print(f"[release] wrote {RELEASE_DIR}")
    print(f"[release] wrote {ZIP_PATH}")
    print(f"[release] wrote {ARTIFACTS_PATH}")
    return True


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--smoke",
        action="store_true",
        help="run release-lite smoke and safety tests in staging",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="build and validate in staging, then fail on drift without modifying outputs",
    )
    parser.add_argument(
        "--timeout-seconds",
        type=int,
        default=300,
        help="per-step subprocess timeout (default: 300)",
    )
    args = parser.parse_args(argv)
    try:
        success = build(
            run_smoke=args.smoke,
            timeout_seconds=args.timeout_seconds,
            check=args.check,
        )
    except (OSError, RuntimeError, ValueError, zipfile.BadZipFile) as exc:
        print(f"[release] ERROR: {exc}", file=sys.stderr)
        return 1
    return 0 if success else 1


if __name__ == "__main__":
    raise SystemExit(main())
