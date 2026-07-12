#!/usr/bin/env python3
"""Export the trusted training DataFrame to version-neutral NumPy arrays.

``X_full_filtered.pkl`` was written by pandas 3 and is a trusted, local
development artifact; the required CLI acknowledgment does not make an
untrusted pickle safe.  Canonical conversion is locked to the source/output
hashes and exact environment in ``training_tools/feature_export_lock.json``.
The default writes both representations: ``X_full.npy`` as float32 for exact
reproduction of the cancer-type model's historical training input, and
``X_full_float64.npy`` for the binary/LOCO analyses that used the original
float64 DataFrame.
"""

from __future__ import annotations

import argparse
import json
import os
import pickle
import platform
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

from training_tools import (  # noqa: E402
    FEATURE_EXPORT_LOCK_PATH,
    FEATURE_MANIFEST_NAME,
    canonical_feature_export_issues,
    commit_file_generation,
    file_record,
    load_feature_export_lock,
    snapshot_inputs,
    validate_output_directory_path,
    verify_input_snapshot,
)

MANAGED_EXPORT_NAMES = [
    "X_full.npy",
    "X_full_float64.npy",
    "X_genes.npy",
    "X_samples.npy",
    FEATURE_MANIFEST_NAME,
]


def atomic_save(path: Path, array) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as handle:
            np.save(handle, array, allow_pickle=False)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def _converter_environment() -> dict[str, str]:
    return {
        "implementation": platform.python_implementation(),
        "python": ".".join(platform.python_version_tuple()[:2]),
        "numpy": np.__version__,
        "pandas": pd.__version__,
    }


def export(
    source: Path,
    output_dir: Path,
    dtype: str,
    *,
    force: bool = False,
    trusted_source_pickle: bool = False,
    allow_noncanonical_export: bool = False,
) -> dict:
    if not trusted_source_pickle:
        raise ValueError(
            "trusted_source_pickle=True is required because loading a pickle can execute code"
        )
    source = source.resolve()
    if not source.is_file():
        raise ValueError(f"trusted source pickle not found: {source}")
    if dtype not in {"both", "float32", "float64"}:
        raise ValueError("dtype must be one of: both, float32, float64")
    output_dir = Path(output_dir)
    validate_output_directory_path(output_dir)
    resolved_output = output_dir.resolve(strict=False)
    managed_destinations = {
        (resolved_output / name).resolve(strict=False) for name in MANAGED_EXPORT_NAMES
    }
    if source in managed_destinations:
        raise ValueError(
            "trusted source pickle must not collide with a managed export destination: "
            f"{source}"
        )
    if pd.__version__.split(".", 1)[0].isdigit() and int(pd.__version__.split(".", 1)[0]) < 3:
        raise ValueError("this pickle requires pandas>=3; use a compatible conversion environment")
    lock = load_feature_export_lock()
    source_snapshot = snapshot_inputs({"source": source})
    source_record = dict(source_snapshot["source"])
    source_record["path"] = source.name
    converter_environment = _converter_environment()
    preflight_issues = []
    for field in ("bytes", "sha256"):
        if source_record[field] != lock["source"][field]:
            preflight_issues.append(f"source {field} differs")
    if converter_environment != lock["converter_environment"]:
        preflight_issues.append("converter environment differs")
    if preflight_issues and not allow_noncanonical_export:
        raise ValueError(
            "canonical feature export lock verification failed before pickle load: "
            + "; ".join(preflight_issues)
        )
    with source.open("rb") as handle:
        frame = pickle.load(handle)  # noqa: S301 - explicitly trusted development artifact.
    if not isinstance(frame, pd.DataFrame):
        raise ValueError("trusted source pickle must contain a pandas DataFrame")
    if frame.empty or frame.shape[1] == 0:
        raise ValueError("training feature matrix must contain samples and genes")
    if frame.index.hasnans or frame.columns.hasnans:
        raise ValueError("training feature matrix has missing sample or gene identifiers")
    if getattr(frame.index, "nlevels", 1) != 1 or getattr(frame.columns, "nlevels", 1) != 1:
        raise ValueError("training feature matrix axes must be one-dimensional indexes")
    samples = np.asarray([str(value) for value in frame.index], dtype=str)
    genes = np.asarray([str(value) for value in frame.columns], dtype=str)
    for name, values in [("sample", samples), ("gene", genes)]:
        if any(not value.strip() or value != value.strip() for value in values):
            raise ValueError(
                f"training feature matrix has blank or padded {name} identifiers"
            )
        if len(set(values)) != len(values):
            raise ValueError(f"training feature matrix has duplicate {name} identifiers")
    source_values = frame.to_numpy(dtype=np.float64, copy=True)
    if not np.isfinite(source_values).all():
        raise ValueError("training feature matrix contains NaN or infinite values")

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(
        tempfile.mkdtemp(
            prefix=f".{output_dir.name}.feature-staging.", dir=output_dir.parent
        )
    )
    paths = {
        "genes": stage / "X_genes.npy",
        "samples": stage / "X_samples.npy",
    }
    try:
        if dtype in {"both", "float32"}:
            paths["features_float32"] = stage / "X_full.npy"
            atomic_save(paths["features_float32"], source_values.astype(np.float32))
        if dtype in {"both", "float64"}:
            paths["features_float64"] = stage / "X_full_float64.npy"
            atomic_save(paths["features_float64"], source_values)
        atomic_save(paths["genes"], genes)
        atomic_save(paths["samples"], samples)

        array_metadata = {
            "genes": (genes.dtype, genes.shape),
            "samples": (samples.dtype, samples.shape),
        }
        if "features_float32" in paths:
            array_metadata["features_float32"] = (
                np.dtype(np.float32),
                source_values.shape,
            )
        if "features_float64" in paths:
            array_metadata["features_float64"] = (
                np.dtype(np.float64),
                source_values.shape,
            )
        output_manifest = {}
        for name, path in paths.items():
            record = file_record(path, relative_name=path.name)
            record["dtype"] = str(np.dtype(array_metadata[name][0]))
            record["shape"] = [int(value) for value in array_metadata[name][1]]
            output_manifest[name] = record
        manifest = {
            "schema_version": "3.0",
            "source": source_record,
            "shape": [int(source_values.shape[0]), int(source_values.shape[1])],
            "requested_export": dtype,
            "scientific_contract": {
                "features_float32": "historical cancer-type classifier input (lossy)",
                "features_float64": "exact binary and LOCO analysis input",
            },
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "converter_environment": converter_environment,
            "canonical_lock": file_record(
                FEATURE_EXPORT_LOCK_PATH,
                relative_name="training_tools/feature_export_lock.json",
            ),
            "outputs": output_manifest,
        }
        canonical_issues = canonical_feature_export_issues(manifest)
        manifest["canonical_lock_verified"] = not canonical_issues
        if canonical_issues:
            manifest["canonical_lock_issues"] = canonical_issues
            if not allow_noncanonical_export:
                raise ValueError(
                    "canonical feature export lock verification failed: "
                    + "; ".join(canonical_issues)
                )
        manifest_path = stage / FEATURE_MANIFEST_NAME
        atomic_write_text(
            manifest_path,
            json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n",
        )
        verify_input_snapshot(source_snapshot)
        staged_names = [path.name for path in paths.values()] + [FEATURE_MANIFEST_NAME]
        commit_file_generation(
            stage,
            output_dir,
            staged_names=staged_names,
            managed_names=MANAGED_EXPORT_NAMES,
            manifest_name=FEATURE_MANIFEST_NAME,
            force=force,
        )
        return manifest
    finally:
        if stage.exists():
            shutil.rmtree(stage, ignore_errors=True)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=ROOT / "X_full_filtered.pkl")
    parser.add_argument("--output-dir", type=Path, default=HERE)
    parser.add_argument(
        "--dtype",
        choices=["both", "float64", "float32"],
        default="both",
        help="write both historical representations by default",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="replace an existing managed feature generation transactionally",
    )
    parser.add_argument(
        "--trusted-source-pickle",
        action="store_true",
        help="acknowledge that the pickle was independently verified as trusted",
    )
    parser.add_argument(
        "--allow-noncanonical-export",
        action="store_true",
        help=(
            "allow a development export that does not match the committed source, "
            "converter, and output lock"
        ),
    )
    args = parser.parse_args(argv)
    if not args.trusted_source_pickle:
        parser.error(
            "--trusted-source-pickle is required because loading a pickle can execute code"
        )
    try:
        manifest = export(
            args.source,
            args.output_dir,
            args.dtype,
            force=args.force,
            trusted_source_pickle=True,
            allow_noncanonical_export=args.allow_noncanonical_export,
        )
    except (
        OSError,
        ValueError,
        TypeError,
        AttributeError,
        NotImplementedError,
        pickle.UnpicklingError,
    ) as exc:
        parser.error(str(exc))
    print(
        f"exported {manifest['requested_export']} feature arrays "
        f"{tuple(manifest['shape'])} to {args.output_dir.resolve()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
