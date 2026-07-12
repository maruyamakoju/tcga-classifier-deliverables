"""Shared integrity, provenance, and transaction helpers for model training.

This module is deliberately separate from :mod:`tcga_rnaseq`: the public
scoring package remains dependency-light, while the development-only training
programs can enforce stronger filesystem and environment contracts.
"""

from __future__ import annotations

import contextlib
import hashlib
import json
import os
import platform
import shutil
import stat
import subprocess
import sys
import tempfile
import uuid
from pathlib import Path

import numpy as np
import pandas as pd


CANONICAL_TRAINING_ENVIRONMENT = {
    "python": "3.11",
    "numpy": "1.26.4",
    "pandas": "2.3.3",
    "scipy": "1.15.3",
    "scikit_learn": "1.8.0",
}

FEATURE_MANIFEST_NAME = "X_full.export.json"
FEATURE_EXPORT_LOCK_PATH = Path(__file__).with_name("feature_export_lock.json")


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    """Return the SHA-256 of *path* without loading it all into memory."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_record(path: Path, *, relative_name: str | None = None) -> dict:
    """Build a JSON-safe size/hash record for one existing regular file."""
    path = Path(path).resolve()
    if not path.is_file():
        raise ValueError(f"expected a regular file: {path}")
    return {
        "path": relative_name if relative_name is not None else str(path),
        "bytes": int(path.stat().st_size),
        "sha256": sha256_file(path),
    }


def snapshot_inputs(paths: dict[str, Path]) -> dict[str, dict]:
    """Hash all inputs before computation so later mutation can be detected."""
    records = {}
    for name, raw_path in paths.items():
        path = Path(raw_path).resolve()
        record = file_record(path)
        record["path"] = str(path)
        records[name] = record
    return records


def verify_input_snapshot(records: dict[str, dict]) -> None:
    """Fail if any snapshotted input changed size or content."""
    for name, expected in records.items():
        path = Path(expected["path"])
        if not path.is_file():
            raise ValueError(f"input changed during the run ({name} disappeared): {path}")
        size = int(path.stat().st_size)
        digest = sha256_file(path)
        if size != int(expected["bytes"]) or digest != expected["sha256"]:
            raise ValueError(f"input changed during the run ({name}): {path}")


def _load_json_object(path: Path) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"could not read feature export manifest {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"feature export manifest must contain a JSON object: {path}")
    return value


def load_feature_export_lock(path: Path = FEATURE_EXPORT_LOCK_PATH) -> dict:
    """Load the committed, path-neutral canonical feature-export contract."""
    path = Path(path)
    lock = _load_json_object(path)
    if lock.get("schema_version") != "1.0":
        raise ValueError(f"feature export lock must use schema_version 1.0: {path}")
    if lock.get("artifact") != "canonical-tcga-feature-export":
        raise ValueError(f"unexpected feature export lock artifact identity: {path}")
    for field in ("source", "converter_environment", "outputs"):
        if not isinstance(lock.get(field), dict):
            raise ValueError(f"feature export lock {field} must be an object: {path}")
    return lock


def canonical_feature_export_issues(manifest: dict) -> list[str]:
    """Return path-neutral differences from the committed feature-export lock."""
    lock = load_feature_export_lock()
    issues = []
    lock_record = manifest.get("canonical_lock")
    expected_lock_record = file_record(
        FEATURE_EXPORT_LOCK_PATH,
        relative_name="training_tools/feature_export_lock.json",
    )
    if not isinstance(lock_record, dict) or any(
        lock_record.get(field) != expected_lock_record[field]
        for field in ("path", "bytes", "sha256")
    ):
        issues.append("canonical lock identity/hash differs")

    source = manifest.get("source")
    expected_source = lock["source"]
    if not isinstance(source, dict):
        issues.append("source record is missing")
    else:
        source_name = source.get("path")
        if (
            not isinstance(source_name, str)
            or Path(source_name).name != source_name
            or Path(source_name).is_absolute()
        ):
            issues.append("source path is not a path-neutral basename")
        for field in ("bytes", "sha256"):
            if source.get(field) != expected_source.get(field):
                issues.append(f"source {field} differs")

    environment = manifest.get("converter_environment")
    if environment != lock["converter_environment"]:
        issues.append("converter environment differs")

    outputs = manifest.get("outputs")
    if not isinstance(outputs, dict):
        issues.append("outputs record is missing")
        return issues
    expected_outputs = lock["outputs"]
    for role, record in outputs.items():
        expected = expected_outputs.get(role)
        if expected is None:
            issues.append(f"unexpected output role {role}")
            continue
        if not isinstance(record, dict):
            issues.append(f"output {role} is not an object")
            continue
        for field in ("path", "bytes", "sha256", "dtype", "shape"):
            if record.get(field) != expected.get(field):
                issues.append(f"output {role} {field} differs")
    return issues


def validate_feature_manifest(
    features_path: Path,
    genes_path: Path,
    samples_path: Path,
    *,
    expected_dtype,
    allow_unverified: bool,
    require_canonical_lock: bool = True,
) -> tuple[Path | None, dict | None]:
    """Verify that a feature matrix and its axes belong to one export generation.

    Canonical runs require the schema-v3 manifest emitted by
    ``export_features_npy.py``.  The explicit escape hatch skips only manifest
    verification; canonical dtype validation remains the caller's obligation.
    """
    if allow_unverified:
        return None, None

    features_path = Path(features_path).resolve()
    genes_path = Path(genes_path).resolve()
    samples_path = Path(samples_path).resolve()
    manifest_path = features_path.with_name(FEATURE_MANIFEST_NAME)
    if not manifest_path.is_file():
        raise ValueError(
            f"verified feature export manifest is required: {manifest_path}; "
            "use --allow-unverified-features only for deliberate development inputs"
        )
    manifest = _load_json_object(manifest_path)
    if manifest.get("schema_version") != "3.0":
        raise ValueError(
            f"feature export manifest must use schema_version 3.0: {manifest_path}"
        )
    if require_canonical_lock:
        issues = canonical_feature_export_issues(manifest)
        if manifest.get("canonical_lock_verified") is not True or issues:
            detail = "; ".join(issues) if issues else "manifest is marked noncanonical"
            raise ValueError(
                "feature export is not verified by the committed canonical lock: "
                + detail
            )
    outputs = manifest.get("outputs")
    if not isinstance(outputs, dict):
        raise ValueError("feature export manifest outputs must be an object")

    required = {
        "features": (features_path, np.dtype(expected_dtype)),
        "genes": (genes_path, None),
        "samples": (samples_path, None),
    }
    chosen_records = {}
    for role, (path, dtype) in required.items():
        matches = [
            record
            for record in outputs.values()
            if isinstance(record, dict) and record.get("path") == path.name
        ]
        if len(matches) != 1:
            raise ValueError(
                f"feature export manifest must bind exactly one {role} file named {path.name}"
            )
        record = matches[0]
        for field in ("bytes", "sha256", "dtype", "shape"):
            if field not in record:
                raise ValueError(f"feature export manifest {role} record lacks {field}")
        actual = file_record(path, relative_name=path.name)
        if (
            int(record["bytes"]) != actual["bytes"]
            or record["sha256"] != actual["sha256"]
        ):
            raise ValueError(f"feature export manifest hash mismatch for {role}: {path}")
        if dtype is not None and np.dtype(record["dtype"]) != dtype:
            raise ValueError(
                f"feature export manifest dtype mismatch for {role}: "
                f"expected {dtype}, found {record['dtype']}"
            )
        chosen_records[role] = record

    feature_shape = chosen_records["features"]["shape"]
    gene_shape = chosen_records["genes"]["shape"]
    sample_shape = chosen_records["samples"]["shape"]
    if (
        not isinstance(feature_shape, list)
        or len(feature_shape) != 2
        or gene_shape != [feature_shape[1]]
        or sample_shape != [feature_shape[0]]
    ):
        raise ValueError("feature export manifest contains inconsistent axis shapes")
    return manifest_path, manifest


def _run_git(root: Path, args: list[str]) -> bytes | None:
    try:
        result = subprocess.run(
            ["git", *args],
            cwd=root,
            check=True,
            capture_output=True,
            timeout=20,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return result.stdout


def code_provenance(root: Path, source_paths: list[Path]) -> dict:
    """Record code hashes plus an honest clean/dirty Git attribution."""
    root = Path(root).resolve()
    commit_raw = _run_git(root, ["rev-parse", "HEAD"])
    status_raw = _run_git(root, ["status", "--porcelain=v1", "--untracked-files=all"])
    diff_raw = _run_git(root, ["diff", "--binary", "HEAD"])
    sources = {}
    for raw_path in source_paths:
        path = Path(raw_path).resolve()
        try:
            label = path.relative_to(root).as_posix()
        except ValueError:
            label = str(path)
        sources[label] = sha256_file(path)
    return {
        "git_commit": None if commit_raw is None else commit_raw.decode().strip() or None,
        "git_dirty": None if status_raw is None else bool(status_raw.strip()),
        "git_diff_sha256": (
            None if diff_raw is None else hashlib.sha256(diff_raw).hexdigest()
        ),
        "source_sha256": sources,
    }


def training_environment() -> dict:
    """Return actual versions and the immutable canonical training contract."""
    import scipy
    import sklearn

    actual = {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "scikit_learn": sklearn.__version__,
    }
    comparisons = {
        "python": actual["python"].split(".")[:2]
        == CANONICAL_TRAINING_ENVIRONMENT["python"].split("."),
        **{
            key: actual[key] == expected
            for key, expected in CANONICAL_TRAINING_ENVIRONMENT.items()
            if key != "python"
        },
    }
    return {
        "actual": actual,
        "canonical": dict(CANONICAL_TRAINING_ENVIRONMENT),
        "canonical_match": all(comparisons.values()),
        "canonical_component_match": comparisons,
    }


def _remove_path(path: Path) -> None:
    if path.is_symlink() or path.is_file():
        path.unlink(missing_ok=True)
    elif path.exists():
        shutil.rmtree(path)


def _is_reparse_point(path: Path) -> bool:
    """Return whether *path* is a Windows reparse point (including junctions)."""
    try:
        attributes = getattr(os.lstat(path), "st_file_attributes", 0)
    except (FileNotFoundError, OSError):
        return False
    return bool(attributes & getattr(stat, "FILE_ATTRIBUTE_REPARSE_POINT", 0))


def validate_output_directory_path(output_dir: Path) -> Path:
    """Reject non-directory and link-like output locations before resolution."""
    output = Path(output_dir)
    if output.is_symlink() or _is_reparse_point(output):
        raise ValueError(f"output directory must not be a symlink or reparse point: {output}")
    if output.exists() and not output.is_dir():
        raise ValueError(
            f"output directory path exists but is not a directory: {output}"
        )
    return output


def _output_lock_path(target: Path) -> Path:
    """Map one resolved output identity to a persistent per-user OS lock file."""
    identity = os.path.normcase(str(Path(target).resolve(strict=False)))
    digest = hashlib.sha256(identity.encode("utf-8")).hexdigest()
    user_component = str(os.getuid()) if hasattr(os, "getuid") else "windows-user"
    lock_root = (
        Path(tempfile.gettempdir())
        / f"tcga-classifier-training-locks-{user_component}"
    )
    lock_root.mkdir(mode=0o700, parents=True, exist_ok=True)
    if lock_root.is_symlink() or _is_reparse_point(lock_root) or not lock_root.is_dir():
        raise ValueError(f"unsafe training lock directory: {lock_root}")
    if hasattr(os, "getuid"):
        root_stat = lock_root.stat()
        if root_stat.st_uid != os.getuid():
            raise ValueError(f"training lock directory has an unexpected owner: {lock_root}")
        lock_root.chmod(0o700)
    return lock_root / f"{digest}.lock"


@contextlib.contextmanager
def exclusive_output_lock(target: Path):
    """Hold a non-blocking, cross-platform single-writer lock for *target*."""
    lock_path = _output_lock_path(Path(target))
    if lock_path.is_symlink() or _is_reparse_point(lock_path):
        raise ValueError(f"unsafe training lock file: {lock_path}")
    flags = os.O_RDWR | os.O_CREAT | getattr(os, "O_BINARY", 0)
    flags |= getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(lock_path, flags, 0o600)
    with os.fdopen(descriptor, "r+b") as handle:
        handle.seek(0, os.SEEK_END)
        if handle.tell() == 0:
            handle.write(b"\0")
            handle.flush()
        handle.seek(0)
        try:
            if os.name == "nt":
                import msvcrt

                msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                import fcntl

                fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as exc:
            raise ValueError(
                f"another process is already writing this output generation: {target}"
            ) from exc
        try:
            yield
        finally:
            handle.seek(0)
            if os.name == "nt":
                import msvcrt

                msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                import fcntl

                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def ensure_output_isolated(output_dir: Path, protected_paths: list[Path]) -> Path:
    """Reject output directories that contain an input or source file."""
    output = Path(output_dir).resolve()
    for raw_path in protected_paths:
        protected = Path(raw_path).resolve()
        try:
            protected.relative_to(output)
        except ValueError:
            continue
        raise ValueError(
            f"output directory must not contain an input or source path: "
            f"{output} contains {protected}"
        )
    return output


@contextlib.contextmanager
def staged_output_directory(
    output_dir: Path,
    *,
    force: bool,
    protected_paths: list[Path],
):
    """Stage one complete directory and promote it with single-writer rollback."""
    raw_output = Path(output_dir)
    with exclusive_output_lock(raw_output):
        validate_output_directory_path(raw_output)
        output = ensure_output_isolated(raw_output, protected_paths)
        output.parent.mkdir(parents=True, exist_ok=True)
        if output.exists() and not force:
            raise ValueError(
                f"output directory already exists: {output}; pass --force to replace it"
            )
        stage = Path(
            tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=output.parent)
        )
        committed = False
        backup = output.parent / f".{output.name}.backup.{uuid.uuid4().hex}"
        try:
            yield stage
            if output.exists():
                os.replace(output, backup)
            try:
                os.replace(stage, output)
                committed = True
            except BaseException:
                if backup.exists():
                    if output.exists():
                        _remove_path(output)
                    os.replace(backup, output)
                raise
            if backup.exists():
                _remove_path(backup)
        finally:
            if not committed and stage.exists():
                _remove_path(stage)
            if backup.exists() and not committed:
                if output.exists():
                    _remove_path(output)
                os.replace(backup, output)
            elif backup.exists():
                _remove_path(backup)


def commit_file_generation(
    stage_dir: Path,
    output_dir: Path,
    *,
    staged_names: list[str],
    managed_names: list[str],
    manifest_name: str,
    force: bool,
) -> None:
    """Commit a managed file set with rollback and manifest-last visibility.

    This variant is for feature exports whose destination also contains source
    files, so replacing the whole directory would be destructive.  Consumers
    reject the generation while its manifest is absent, and Python-level
    failures restore every previous managed file.
    """
    stage = Path(stage_dir).resolve()
    raw_output = Path(output_dir)
    if manifest_name not in staged_names or staged_names[-1] != manifest_name:
        raise ValueError("the generation manifest must be staged and committed last")
    if len(set(staged_names)) != len(staged_names):
        raise ValueError("staged generation names must be unique")
    for name in staged_names:
        if Path(name).name != name or not (stage / name).is_file():
            raise ValueError(f"invalid or missing staged generation file: {name}")

    with exclusive_output_lock(raw_output):
        validate_output_directory_path(raw_output)
        output = raw_output.resolve()
        output.mkdir(parents=True, exist_ok=True)
        existing = [name for name in managed_names if (output / name).exists()]
        if existing and not force:
            raise ValueError(
                "feature generation already exists ("
                + ", ".join(existing)
                + "); pass --force to replace it"
            )
        backup = Path(
            tempfile.mkdtemp(prefix=".feature-generation.backup.", dir=output.parent)
        )
        backed_up = []
        moved_new = []
        completed = False
        rolled_back = False
        try:
            backup_order = [manifest_name] + [
                name for name in managed_names if name != manifest_name
            ]
            for name in backup_order:
                source = output / name
                if source.exists():
                    os.replace(source, backup / name)
                    backed_up.append(name)
            for name in staged_names:
                os.replace(stage / name, output / name)
                moved_new.append(name)
            completed = True
        except BaseException:
            for name in reversed(moved_new):
                (output / name).unlink(missing_ok=True)
            for name in backed_up:
                old = backup / name
                if old.exists():
                    os.replace(old, output / name)
            rolled_back = True
            raise
        finally:
            if backup.exists() and (completed or rolled_back):
                _remove_path(backup)
            if stage.exists():
                _remove_path(stage)


def output_records(paths: dict[str, Path]) -> dict[str, dict]:
    """Create manifest records for a completed staged generation."""
    return {
        name: file_record(path, relative_name=Path(path).name)
        for name, path in paths.items()
    }
