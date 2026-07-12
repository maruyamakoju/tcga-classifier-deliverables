"""Cache integrity and run-provenance helpers for external validation.

This module deliberately has no scikit-learn dependency.  External-validation
downloads are mutable research inputs, so a cache is accepted only when its
semantic fingerprint, byte hash, and tabular axes all match the current run.
"""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import math
import os
import platform
import shutil
import subprocess
import sys
import tempfile
from contextlib import contextmanager
from collections.abc import Callable, Mapping, Sequence
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


CACHE_SCHEMA_VERSION = "3"
RUN_MANIFEST_SCHEMA_VERSION = "1"
PROVENANCE_SOURCE_PATH = Path(__file__).resolve()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def _jsonable(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Mapping):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, set):
        return sorted(_jsonable(item) for item in value)
    if isinstance(value, np.ndarray):
        return _jsonable(value.tolist())
    if isinstance(value, np.generic):
        return _jsonable(value.item())
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("JSON provenance must not contain NaN or infinite values")
        return value
    if isinstance(value, (str, int, bool)) or value is None:
        return value
    return str(value)


def canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(
        _jsonable(value), sort_keys=True, separators=(",", ":"), ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")


def _reject_nonfinite_json_constant(token: str):
    raise ValueError(f"non-finite JSON constant {token}")


def semantic_fingerprint(components: Mapping[str, Any]) -> str:
    """Hash every declared semantic input to a derived cache object."""
    return hashlib.sha256(canonical_json_bytes(components)).hexdigest()


def sequence_hash(values: Sequence[Any]) -> str:
    return hashlib.sha256(canonical_json_bytes([str(value) for value in values])).hexdigest()


def sha256_file(path: os.PathLike[str] | str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def cache_meta_path(cache_path: os.PathLike[str] | str) -> Path:
    path = Path(cache_path)
    return path.with_suffix(path.suffix + ".meta.json")


def safe_cache_key(identifier: Any, namespace: str) -> str:
    """Return a path-safe, namespace-bound key for an external identifier."""
    value = str(identifier)
    if not value or value != value.strip():
        raise ValueError(f"{namespace} identifier must be non-empty and unpadded")
    return hashlib.sha256(
        canonical_json_bytes({"namespace": namespace, "identifier": value})
    ).hexdigest()


def contained_cache_path(
    cache_root: os.PathLike[str] | str,
    identifier: Any,
    *,
    namespace: str,
    suffix: str,
) -> Path:
    """Build a hashed cache path and prove that it remains below ``cache_root``."""
    root = Path(cache_root).resolve(strict=False)
    key = f"{safe_cache_key(identifier, namespace)}{suffix}"
    # safe_cache_key returns a bare sha256 hex digest, so ``key`` is a single path
    # component with no separators or traversal, and ``root / key`` is contained by
    # construction. Verify that directly instead of re-resolving the child and
    # calling relative_to(): on Windows, Path.resolve() adds a \\?\ extended-length
    # prefix to only one side once the full path crosses MAX_PATH, which made
    # relative_to() spuriously reject legitimate long cache directories.
    if key != os.path.basename(key) or key in {"", ".", ".."} or "/" in key or "\\" in key:
        raise ValueError(f"cache path escapes cache root: {root / key}")
    return root / key


def validate_source_revision(source_revision: Any, *, live: bool) -> str:
    revision = str(source_revision).strip()
    if not revision:
        raise ValueError("--source-revision must be non-empty")
    if live and revision.casefold() == "unversioned":
        raise ValueError(
            "live external validation requires a provider snapshot/revision; "
            "replace --source-revision unversioned with a concrete value"
        )
    return revision


def _paths_alias(left: Path, right: Path) -> bool:
    if left.resolve(strict=False) == right.resolve(strict=False):
        return True
    try:
        return left.exists() and right.exists() and os.path.samefile(left, right)
    except OSError:
        return False


def validate_managed_paths(
    *,
    protected_inputs: Mapping[str, os.PathLike[str] | str],
    managed_files: Mapping[str, os.PathLike[str] | str],
    managed_directories: Mapping[str, os.PathLike[str] | str] | None = None,
) -> None:
    """Reject aliases between read-only inputs and every mutable run path."""
    protected = {name: Path(path) for name, path in protected_inputs.items()}
    managed = {name: Path(path) for name, path in managed_files.items()}
    directories = {
        name: Path(path) for name, path in (managed_directories or {}).items()
    }

    managed_items = list(managed.items())
    for index, (left_name, left_path) in enumerate(managed_items):
        for right_name, right_path in managed_items[index + 1:]:
            if _paths_alias(left_path, right_path):
                raise ValueError(
                    f"managed paths collide: {left_name} and {right_name}: "
                    f"{left_path.resolve(strict=False)}"
                )

    for input_name, input_path in protected.items():
        for managed_name, managed_path in managed.items():
            if _paths_alias(input_path, managed_path):
                raise ValueError(
                    f"protected input {input_name} collides with {managed_name}: "
                    f"{input_path.resolve(strict=False)}"
                )
        resolved_input = input_path.resolve(strict=False)
        for directory_name, directory in directories.items():
            resolved_directory = directory.resolve(strict=False)
            try:
                resolved_input.relative_to(resolved_directory)
            except ValueError:
                continue
            raise ValueError(
                f"protected input {input_name} is inside mutable directory "
                f"{directory_name}: {resolved_directory}"
            )


@contextmanager
def staged_output_directory(output_dir: os.PathLike[str] | str):
    """Create a sibling staging directory that is always removed on exit."""
    destination = Path(output_dir).resolve(strict=False)
    destination.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(
        dir=destination.parent,
        prefix=f".{destination.name}.external-validation-run-",
    ))
    try:
        yield stage
    finally:
        shutil.rmtree(stage, ignore_errors=True)


def publish_staged_files(
    staged_files: Mapping[str, os.PathLike[str] | str],
    final_files: Mapping[str, os.PathLike[str] | str],
    *,
    manifest_key: str = "run_manifest",
) -> None:
    """Publish a complete output set with rollback and the manifest moved last."""
    if set(staged_files) != set(final_files):
        raise ValueError("staged and final output keys must match")
    if manifest_key not in staged_files:
        raise ValueError(f"staged output set must contain {manifest_key!r}")
    staged = {name: Path(path) for name, path in staged_files.items()}
    final = {name: Path(path) for name, path in final_files.items()}
    missing = [name for name, path in staged.items() if not path.is_file()]
    if missing:
        raise ValueError(f"staged outputs are missing: {missing}")

    validate_managed_paths(protected_inputs={}, managed_files=final)
    for path in final.values():
        path.parent.mkdir(parents=True, exist_ok=True)
    backup_root = Path(tempfile.mkdtemp(
        dir=next(iter(final.values())).parent,
        prefix=".external-validation-backup-",
    ))
    backups: dict[str, Path] = {}
    published: list[str] = []
    order = [name for name in staged if name != manifest_key] + [manifest_key]
    try:
        for index, name in enumerate(order):
            destination = final[name]
            if destination.exists():
                backup = backup_root / f"{index:04d}.bak"
                os.replace(destination, backup)
                backups[name] = backup
            os.replace(staged[name], destination)
            published.append(name)
    except Exception:
        for name in reversed(published):
            try:
                final[name].unlink()
            except FileNotFoundError:
                pass
        for name, backup in backups.items():
            if backup.exists():
                os.replace(backup, final[name])
        raise
    finally:
        shutil.rmtree(backup_root, ignore_errors=True)


def atomic_write_path(path: os.PathLike[str] | str, write_fn: Callable[[Path], Any]) -> None:
    """Write through a unique sibling temp file and atomically replace ``path``."""
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=destination.parent,
        prefix=f".{destination.name}.",
        suffix=".tmp",
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        write_fn(temporary)
        os.replace(temporary, destination)
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def atomic_write_json(path: os.PathLike[str] | str, value: Mapping[str, Any]) -> None:
    def _write(temporary: Path) -> None:
        temporary.write_text(
            json.dumps(
                _jsonable(value), indent=2, sort_keys=True, ensure_ascii=False,
                allow_nan=False,
            ) + "\n",
            encoding="utf-8",
        )

    atomic_write_path(path, _write)


def atomic_write_text(path: os.PathLike[str] | str, text: str) -> None:
    atomic_write_path(path, lambda temporary: temporary.write_text(text, encoding="utf-8"))


def atomic_write_csv(path: os.PathLike[str] | str, frame: pd.DataFrame) -> None:
    text = frame.to_csv(index=False, lineterminator="\n")
    atomic_write_path(path, lambda temporary: temporary.write_bytes(text.encode("utf-8")))


def _axis_metadata(frame: pd.DataFrame) -> dict[str, Any]:
    return {
        "shape": [int(frame.shape[0]), int(frame.shape[1])],
        "index_count": int(len(frame.index)),
        "index_sha256": sequence_hash(frame.index.tolist()),
        "column_count": int(len(frame.columns)),
        "columns_sha256": sequence_hash(frame.columns.tolist()),
    }


def _series_metadata(series: pd.Series) -> dict[str, Any]:
    return {
        "length": int(len(series)),
        "index_sha256": sequence_hash(series.index.tolist()),
        "name": None if series.name is None else str(series.name),
    }


def _base_cache_metadata(
    cache_path: Path,
    fingerprint: str,
    fingerprint_inputs: Mapping[str, Any],
    cache_kind: str,
    object_type: str,
    extra_metadata: Mapping[str, Any] | None,
) -> dict[str, Any]:
    metadata: dict[str, Any] = {
        "cache_schema_version": CACHE_SCHEMA_VERSION,
        "cache_kind": cache_kind,
        "object_type": object_type,
        "fingerprint": fingerprint,
        "fingerprint_inputs": _jsonable(fingerprint_inputs),
        "content_sha256": sha256_file(cache_path),
        "written_at_utc": utc_now(),
    }
    if extra_metadata:
        metadata["extra"] = _jsonable(extra_metadata)
    return metadata


def _require_parquet_path(cache_path: os.PathLike[str] | str) -> Path:
    path = Path(cache_path)
    if path.suffix.lower() != ".parquet":
        raise ValueError(f"executable pickle caches are unsupported; expected .parquet: {path}")
    return path


def write_dataframe_cache(
    cache_path: os.PathLike[str] | str,
    frame: pd.DataFrame,
    *,
    fingerprint: str,
    fingerprint_inputs: Mapping[str, Any],
    cache_kind: str,
    extra_metadata: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    path = _require_parquet_path(cache_path)
    atomic_write_path(path, lambda temporary: frame.to_parquet(temporary, index=True))
    metadata = _base_cache_metadata(
        path, fingerprint, fingerprint_inputs, cache_kind, "dataframe_parquet", extra_metadata
    )
    metadata.update(_axis_metadata(frame))
    atomic_write_json(cache_meta_path(path), metadata)
    return metadata


def write_series_cache(
    cache_path: os.PathLike[str] | str,
    series: pd.Series,
    *,
    fingerprint: str,
    fingerprint_inputs: Mapping[str, Any],
    cache_kind: str,
    extra_metadata: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    path = _require_parquet_path(cache_path)
    stored = series.rename("__value__").to_frame()
    atomic_write_path(path, lambda temporary: stored.to_parquet(temporary, index=True))
    metadata = _base_cache_metadata(
        path, fingerprint, fingerprint_inputs, cache_kind, "series_parquet", extra_metadata
    )
    metadata.update(_series_metadata(series))
    atomic_write_json(cache_meta_path(path), metadata)
    return metadata


def read_cache_metadata(cache_path: os.PathLike[str] | str) -> dict[str, Any] | None:
    path = cache_meta_path(cache_path)
    if not path.exists():
        return None
    try:
        value = json.loads(
            path.read_text(encoding="utf-8"),
            parse_constant=_reject_nonfinite_json_constant,
        )
        value = _jsonable(value)
    except (OSError, json.JSONDecodeError, ValueError):
        return None
    return value if isinstance(value, dict) else None


def _metadata_header_is_valid(
    cache_path: Path,
    metadata: Mapping[str, Any] | None,
    fingerprint: str,
    object_type: str,
) -> bool:
    if not metadata:
        return False
    if metadata.get("cache_schema_version") != CACHE_SCHEMA_VERSION:
        return False
    if metadata.get("fingerprint") != fingerprint:
        return False
    if metadata.get("object_type") != object_type:
        return False
    expected_hash = metadata.get("content_sha256")
    if not isinstance(expected_hash, str):
        return False
    try:
        return sha256_file(cache_path) == expected_hash
    except OSError:
        return False


def load_dataframe_cache(
    cache_path: os.PathLike[str] | str,
    *,
    fingerprint: str,
    expected_index: Sequence[Any] | None = None,
    expected_columns: Sequence[Any] | None = None,
    columns_may_be_subset: bool = False,
    reject_all_missing_columns: bool = False,
) -> pd.DataFrame | None:
    path = Path(cache_path)
    if path.suffix.lower() != ".parquet":
        return None
    if not path.exists():
        return None
    metadata = read_cache_metadata(path)
    if not _metadata_header_is_valid(path, metadata, fingerprint, "dataframe_parquet"):
        return None
    try:
        frame = pd.read_parquet(path)
    except Exception:  # Backend-specific corrupt-Parquet exceptions vary by engine.
        return None
    if not isinstance(frame, pd.DataFrame):
        return None
    if not frame.index.is_unique or not frame.columns.is_unique:
        return None
    axes = _axis_metadata(frame)
    for key, value in axes.items():
        if metadata.get(key) != value:
            return None
    if expected_index is not None:
        if [str(value) for value in frame.index] != [str(value) for value in expected_index]:
            return None
    if expected_columns is not None:
        actual = [str(value) for value in frame.columns]
        expected = [str(value) for value in expected_columns]
        if columns_may_be_subset:
            expected_set = set(expected)
            if any(column not in expected_set for column in actual):
                return None
            if actual != [column for column in expected if column in set(actual)]:
                return None
        elif actual != expected:
            return None
    if reject_all_missing_columns and bool(frame.isna().all(axis=0).any()):
        return None
    return frame


def load_series_cache(
    cache_path: os.PathLike[str] | str,
    *,
    fingerprint: str,
    expected_name: str | None = None,
    allowed_index: Sequence[Any] | None = None,
) -> pd.Series | None:
    path = Path(cache_path)
    if path.suffix.lower() != ".parquet":
        return None
    if not path.exists():
        return None
    metadata = read_cache_metadata(path)
    if not _metadata_header_is_valid(path, metadata, fingerprint, "series_parquet"):
        return None
    try:
        stored = pd.read_parquet(path)
    except Exception:  # Backend-specific corrupt-Parquet exceptions vary by engine.
        return None
    if not isinstance(stored, pd.DataFrame) or stored.columns.tolist() != ["__value__"]:
        return None
    series = stored["__value__"].rename(metadata.get("name"))
    if not isinstance(series, pd.Series) or not series.index.is_unique:
        return None
    for key, value in _series_metadata(series).items():
        if metadata.get(key) != value:
            return None
    if expected_name is not None and str(series.name) != str(expected_name):
        return None
    if allowed_index is not None:
        allowed = {str(value) for value in allowed_index}
        if any(str(value) not in allowed for value in series.index):
            return None
    return series


def write_csv_cache(
    cache_path: os.PathLike[str] | str,
    frame: pd.DataFrame,
    *,
    fingerprint: str,
    fingerprint_inputs: Mapping[str, Any],
    cache_kind: str,
    extra_metadata: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    path = Path(cache_path)
    atomic_write_csv(path, frame)
    metadata = _base_cache_metadata(
        path, fingerprint, fingerprint_inputs, cache_kind, "dataframe_csv", extra_metadata
    )
    metadata.update({
        "shape": [int(frame.shape[0]), int(frame.shape[1])],
        "column_count": int(len(frame.columns)),
        "columns_sha256": sequence_hash(frame.columns.tolist()),
    })
    atomic_write_json(cache_meta_path(path), metadata)
    return metadata


def load_csv_cache(
    cache_path: os.PathLike[str] | str,
    *,
    fingerprint: str,
    required_columns: Sequence[str] | None = None,
) -> pd.DataFrame | None:
    path = Path(cache_path)
    if not path.exists():
        return None
    metadata = read_cache_metadata(path)
    if not _metadata_header_is_valid(path, metadata, fingerprint, "dataframe_csv"):
        return None
    try:
        frame = pd.read_csv(path)
    except (OSError, pd.errors.ParserError, UnicodeError):
        return None
    if not frame.columns.is_unique:
        return None
    if metadata.get("shape") != [int(frame.shape[0]), int(frame.shape[1])]:
        return None
    if metadata.get("columns_sha256") != sequence_hash(frame.columns.tolist()):
        return None
    if required_columns and not set(required_columns).issubset(frame.columns):
        return None
    return frame


def validate_identifier_column(frame: pd.DataFrame, column: str, context: str) -> pd.Series:
    if column not in frame.columns:
        raise ValueError(f"{context} must contain {column!r}")
    values = frame[column]
    if values.isna().any():
        raise ValueError(f"{context} contains missing {column} values")
    raw = values.astype(str)
    normalized = raw.str.strip()
    if normalized.eq("").any():
        raise ValueError(f"{context} contains blank {column} values")
    if not raw.eq(normalized).all():
        raise ValueError(f"{context} contains padded {column} values")
    if normalized.duplicated().any():
        examples = normalized[normalized.duplicated(keep=False)].unique().tolist()[:5]
        raise ValueError(f"{context} contains duplicate {column} values: {examples}")
    return normalized


def validate_unique_strings(values: Sequence[Any], context: str) -> list[str]:
    normalized = [str(value).strip() for value in values]
    if not normalized or any(not value for value in normalized):
        raise ValueError(f"{context} must contain non-empty values")
    if len(set(normalized)) != len(normalized):
        duplicates = sorted({value for value in normalized if normalized.count(value) > 1})[:5]
        raise ValueError(f"{context} contains duplicate values: {duplicates}")
    return normalized


def group_audit(frame: pd.DataFrame, id_column: str, group_column: str) -> dict[str, Any]:
    identifiers = validate_identifier_column(frame, id_column, "sample manifest")
    if group_column not in frame.columns:
        return {
            "id_column": id_column,
            "n_rows": int(len(frame)),
            "n_unique_ids": int(identifiers.nunique()),
            "group_column": group_column,
            "n_unique_groups": None,
            "n_repeated_groups": None,
            "max_rows_per_group": None,
        }
    groups = frame[group_column].astype(str).str.strip()
    usable = groups.ne("") & frame[group_column].notna()
    counts = groups[usable].value_counts()
    return {
        "id_column": id_column,
        "n_rows": int(len(frame)),
        "n_unique_ids": int(identifiers.nunique()),
        "group_column": group_column,
        "n_unique_groups": int(counts.size),
        "n_repeated_groups": int((counts > 1).sum()),
        "max_rows_per_group": int(counts.max()) if not counts.empty else 0,
    }


def scored_dataframe(
    sample_ids: Sequence[Any], probabilities: Sequence[float], threshold: float
) -> pd.DataFrame:
    """Build an external-validation score table without display rounding."""
    samples = []
    for value in sample_ids:
        if value is None or bool(pd.isna(value)):
            raise ValueError("sample IDs must not be missing")
        text = str(value)
        if not text or text != text.strip():
            raise ValueError("sample IDs must be non-empty and have no surrounding whitespace")
        samples.append(text)
    if len(set(samples)) != len(samples):
        raise ValueError("sample IDs must be unique")
    scores = np.asarray(probabilities, dtype=float)
    if scores.ndim != 1 or len(scores) != len(samples):
        raise ValueError("sample_ids and probabilities must be same-length one-dimensional data")
    if not np.isfinite(scores).all():
        raise ValueError("probabilities must contain only finite values")
    if np.any((scores < 0.0) | (scores > 1.0)):
        raise ValueError("probabilities must be in [0, 1]")
    try:
        threshold = float(threshold)
    except (TypeError, ValueError) as exc:
        raise ValueError("threshold must be a finite number in [0, 1]") from exc
    if not math.isfinite(threshold) or not 0.0 <= threshold <= 1.0:
        raise ValueError("threshold must be a finite number in [0, 1]")
    return pd.DataFrame({
        "sample": samples,
        "tumor_probability": scores,
        "call": np.where(scores >= threshold, "tumor", "normal"),
    })


def git_state(root: os.PathLike[str] | str) -> dict[str, Any]:
    root = str(Path(root).resolve())
    try:
        sha = subprocess.run(
            ["git", "-C", root, "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
            timeout=5,
        ).stdout.strip()
        dirty = bool(subprocess.run(
            ["git", "-C", root, "status", "--porcelain"],
            check=True,
            capture_output=True,
            text=True,
            timeout=5,
        ).stdout.strip())
        return {"commit": sha, "dirty": dirty}
    except (OSError, subprocess.CalledProcessError, subprocess.TimeoutExpired):
        return {"commit": None, "dirty": None}


def runtime_state() -> dict[str, Any]:
    versions = {}
    for package in ("numpy", "pandas", "requests"):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = None
    return {
        "python": sys.version.split()[0],
        "implementation": platform.python_implementation(),
        "platform": platform.platform(),
        "packages": versions,
    }


def file_inventory(
    paths: Mapping[str, os.PathLike[str] | str],
    *,
    display_paths: Mapping[str, os.PathLike[str] | str] | None = None,
) -> dict[str, Any]:
    inventory = {}
    for name, raw_path in paths.items():
        path = Path(raw_path)
        display_path = Path((display_paths or {}).get(name, raw_path))
        if path.is_file():
            inventory[name] = {
                "path": str(display_path.resolve(strict=False)),
                "size_bytes": int(path.stat().st_size),
                "sha256": sha256_file(path),
            }
        else:
            inventory[name] = {
                "path": str(display_path.resolve(strict=False)),
                "missing": True,
            }
    return inventory


def write_run_manifest(
    path: os.PathLike[str] | str,
    *,
    root: os.PathLike[str] | str,
    run_kind: str,
    started_at_utc: str,
    argv: Sequence[str],
    parameters: Mapping[str, Any],
    model_path: os.PathLike[str] | str,
    sources: Mapping[str, Any],
    inputs: Mapping[str, os.PathLike[str] | str],
    outputs: Mapping[str, os.PathLike[str] | str],
    alignment: Mapping[str, Any],
    cohort_audit: Mapping[str, Any],
    cache_details: Mapping[str, Any] | None = None,
    source_code: Mapping[str, os.PathLike[str] | str] | None = None,
    output_display_paths: Mapping[str, os.PathLike[str] | str] | None = None,
) -> dict[str, Any]:
    model_path = Path(model_path)
    manifest = {
        "schema_version": RUN_MANIFEST_SCHEMA_VERSION,
        "run_kind": run_kind,
        "started_at_utc": started_at_utc,
        "finished_at_utc": utc_now(),
        "argv": list(argv),
        "parameters": _jsonable(parameters),
        "git": git_state(root),
        "runtime": runtime_state(),
        "source_code": file_inventory({
            "provenance": PROVENANCE_SOURCE_PATH,
            **(source_code or {}),
        }),
        "model": {
            "path": str(model_path.resolve()),
            "size_bytes": int(model_path.stat().st_size),
            "sha256": sha256_file(model_path),
        },
        "sources": _jsonable(sources),
        "inputs": file_inventory(inputs),
        "outputs": file_inventory(outputs, display_paths=output_display_paths),
        "alignment": _jsonable(alignment),
        "cohort_audit": _jsonable(cohort_audit),
        "cache_details": _jsonable(cache_details or {}),
    }
    atomic_write_json(path, manifest)
    return manifest
