#!/usr/bin/env python3
"""Reproduce leave-one-cancer-type-out (LOCO) validation from versioned arrays.

The original LOCO tables were delivered without their generating program.
This runner closes that provenance gap.  Every sample is scored by a pipeline
whose feature selection, scaling, and logistic regression fit excluded the
sample's entire cancer type.

The full feature matrix is intentionally not committed.  Generate the default
``cancer-type-classifier/X_full_float64.npy`` plus its ``X_genes.npy`` and
``X_samples.npy`` siblings with ``export_features_npy.py`` first.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_FEATURES = ROOT / "cancer-type-classifier" / "X_full_float64.npy"
DEFAULT_METADATA = ROOT / "selected_files.csv"
DEFAULT_OUTPUT_DIR = ROOT / "cross-cancer-holdout"

sys.path.insert(0, str(ROOT))
from tcga_rnaseq import metrics as M  # noqa: E402
from training_tools import (  # noqa: E402
    code_provenance,
    output_records,
    snapshot_inputs,
    staged_output_directory,
    training_environment,
    validate_feature_manifest,
    verify_input_snapshot,
)


SOURCE_PATHS = [
    Path(__file__).resolve(),
    ROOT / "training_tools" / "__init__.py",
    ROOT / "tcga_rnaseq" / "metrics.py",
]


def _atomic_write_text(path: Path, text: str) -> None:
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


def _atomic_write_csv(frame: pd.DataFrame, path: Path) -> None:
    _atomic_write_text(path, frame.to_csv(index=False, float_format="%.17g", lineterminator="\n"))


def load_inputs(
    features_path: Path,
    metadata_path: Path,
    *,
    allow_unverified_features: bool = False,
):
    features_path = features_path.resolve()
    genes_path = features_path.with_name("X_genes.npy")
    samples_path = features_path.with_name("X_samples.npy")
    paths = {
        "features": features_path,
        "genes": genes_path,
        "samples": samples_path,
        "metadata": metadata_path.resolve(),
    }
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise ValueError("missing LOCO input files: " + ", ".join(missing))

    manifest_path, _ = validate_feature_manifest(
        paths["features"],
        paths["genes"],
        paths["samples"],
        expected_dtype=np.float64,
        allow_unverified=allow_unverified_features,
    )
    if manifest_path is not None:
        paths["feature_manifest"] = manifest_path
    input_snapshot = snapshot_inputs(paths)

    X = np.load(features_path, mmap_mode="r", allow_pickle=False)
    genes = np.load(genes_path, allow_pickle=False).astype(str)
    samples = np.load(samples_path, allow_pickle=False).astype(str)
    if X.dtype != np.float64:
        raise ValueError(
            f"canonical LOCO requires the exact float64 matrix, found {X.dtype}"
        )
    if genes.ndim != 1 or samples.ndim != 1:
        raise ValueError("X_genes.npy and X_samples.npy must be one-dimensional arrays")
    if X.ndim != 2 or X.shape != (len(samples), len(genes)):
        raise ValueError(
            "feature matrix shape must equal (number of samples, number of genes): "
            f"X={X.shape}, samples={len(samples)}, genes={len(genes)}"
        )
    if X.shape[0] == 0 or X.shape[1] == 0:
        raise ValueError("feature matrix must contain samples and genes")
    if not np.isfinite(X).all():
        raise ValueError("feature matrix contains NaN or infinite values")
    for name, values in [("sample", samples), ("gene", genes)]:
        if any(not value.strip() or value != value.strip() for value in values):
            raise ValueError(f"{name} identifiers must be non-empty and unpadded")
        if len(set(values)) != len(values):
            raise ValueError(f"{name} identifiers must be unique")

    metadata = pd.read_csv(metadata_path, keep_default_na=False, dtype=str)
    required_columns = {"file_id", "project", "label", "case_id"}
    absent = sorted(required_columns - set(metadata.columns))
    if absent:
        raise ValueError("metadata is missing required columns: " + ", ".join(absent))
    if metadata["file_id"].duplicated().any():
        raise ValueError("metadata contains duplicate file_id values")
    for field in required_columns:
        values = metadata[field].astype(str)
        if values.str.strip().eq("").any() or values.ne(values.str.strip()).any():
            raise ValueError(f"metadata column {field} contains blank or padded values")
    aligned = metadata.set_index("file_id").reindex(samples)
    metadata_fields = list(required_columns - {"file_id"})
    if (
        aligned[metadata_fields].isna().any(axis=None)
        or aligned[metadata_fields].eq("").any(axis=None)
    ):
        raise ValueError("metadata is missing project, label, or case_id for feature samples")
    unknown_labels = sorted(set(aligned["label"]) - {"tumor", "normal"})
    if unknown_labels:
        raise ValueError(f"metadata contains unsupported labels: {unknown_labels[:5]}")
    y = (aligned["label"].to_numpy() == "tumor").astype(int)
    projects = aligned["project"].astype(str).to_numpy()
    cases = aligned["case_id"].astype(str).to_numpy()
    if np.unique(projects).size < 2:
        raise ValueError("LOCO requires at least two cancer types")
    if aligned.groupby("case_id")["project"].nunique().gt(1).any():
        raise ValueError("each case_id must belong to exactly one cancer type")
    verify_input_snapshot(input_snapshot)
    return X, genes, samples, y, projects, cases, input_snapshot


def run_loco(X, samples, y, projects, cases, n_features=2000, c_value=0.1):
    from sklearn.exceptions import ConvergenceWarning
    from sklearn.feature_selection import SelectKBest, f_classif
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    X = np.asarray(X)
    samples = np.asarray(samples)
    y = np.asarray(y)
    projects = np.asarray(projects)
    cases = np.asarray(cases)
    if X.ndim != 2 or X.dtype != np.float64 or not np.isfinite(X).all():
        raise ValueError("LOCO X must be a finite two-dimensional float64 matrix")
    for name, values in [
        ("samples", samples),
        ("y", y),
        ("projects", projects),
        ("cases", cases),
    ]:
        if values.ndim != 1 or len(values) != len(X):
            raise ValueError(f"LOCO {name} must be a 1D array aligned to X rows")
    if len(set(samples.astype(str))) != len(samples):
        raise ValueError("LOCO sample identifiers must be unique")
    for name, values in [("sample", samples), ("project", projects), ("case", cases)]:
        strings = values.astype(str)
        if any(not value.strip() or value != value.strip() for value in strings):
            raise ValueError(f"LOCO {name} identifiers must be non-empty and unpadded")
    try:
        y_numeric = y.astype(float)
    except (TypeError, ValueError) as exc:
        raise ValueError("LOCO labels must contain only binary 0/1 values") from exc
    if not np.isfinite(y_numeric).all() or not np.isin(y_numeric, [0.0, 1.0]).all():
        raise ValueError("LOCO labels must contain only binary 0/1 values")
    y = y_numeric.astype(int)
    if isinstance(n_features, (bool, np.bool_)) or not isinstance(
        n_features, (int, np.integer)
    ):
        raise ValueError("n_features must be an integer")
    if not 1 <= n_features <= X.shape[1]:
        raise ValueError(f"n_features must be between 1 and {X.shape[1]}")
    if not np.isfinite(c_value) or c_value <= 0:
        raise ValueError("c_value must be finite and positive")

    prediction_parts = []
    metric_rows = []
    for held_out in sorted(np.unique(projects)):
        test_mask = projects == held_out
        train_mask = ~test_mask
        if np.unique(y[train_mask]).size != 2 or np.unique(y[test_mask]).size != 2:
            raise ValueError(f"{held_out} does not contain both classes in train and test")
        if set(cases[train_mask]) & set(cases[test_mask]):
            raise ValueError(f"patient overlap detected for held-out type {held_out}")

        pipeline = Pipeline([
            ("select", SelectKBest(f_classif, k=n_features)),
            ("scale", StandardScaler()),
            ("classifier", LogisticRegression(
                max_iter=5000,
                C=float(c_value),
                class_weight="balanced",
                random_state=42,
                solver="lbfgs",
                tol=1e-4,
                fit_intercept=True,
            )),
        ])
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            pipeline.fit(X[train_mask], y[train_mask])
        convergence = [
            item for item in caught if issubclass(item.category, ConvergenceWarning)
        ]
        if convergence:
            raise RuntimeError(
                f"LOCO model for {held_out} did not converge: {convergence[0].message}"
            )
        for item in caught:
            warnings.warn_explicit(
                str(item.message), item.category, item.filename, item.lineno
            )
        probability = pipeline.predict_proba(X[test_mask])[:, 1]
        if not np.all(np.isfinite(probability)):
            raise ValueError(f"non-finite LOCO probabilities for {held_out}")
        binary = y[test_mask]
        metrics = M.classification_metrics(binary, probability, threshold=0.5)
        metric_rows.append({
            "held_out_type": held_out,
            "model": "logistic_regression",
            "n_test": int(test_mask.sum()),
            "n_normal": int((binary == 0).sum()),
            "n_tumor": int((binary == 1).sum()),
            "auc": metrics["auc"],
            "avg_precision": M.average_precision(binary, probability),
            "accuracy": metrics["accuracy"],
            "f1": metrics["f1"],
            "precision": metrics["precision"],
            "recall": metrics["recall"],
        })
        prediction_parts.append(pd.DataFrame({
            "sample": samples[test_mask],
            "case_id": cases[test_mask],
            "held_out_type": held_out,
            "label": binary,
            "tumor_probability": probability,
            "call": np.where(probability >= 0.5, "tumor", "normal"),
        }))

    predictions = pd.concat(prediction_parts, ignore_index=True)
    if len(predictions) != len(samples) or predictions["sample"].duplicated().any():
        raise RuntimeError("LOCO predictions do not cover every sample exactly once")
    per_type = pd.DataFrame(metric_rows).sort_values("held_out_type").reset_index(drop=True)
    pooled = M.classification_metrics(
        predictions["label"].to_numpy(),
        predictions["tumor_probability"].to_numpy(),
        threshold=0.5,
    )
    pooled_summary = pd.DataFrame([{
        "model": "logistic_regression",
        "pooled_auc": pooled["auc"],
        "pooled_avg_precision": M.average_precision(
            predictions["label"].to_numpy(),
            predictions["tumor_probability"].to_numpy(),
        ),
        "pooled_accuracy": pooled["accuracy"],
        "pooled_f1": pooled["f1"],
        "macro_mean_auc": float(per_type["auc"].mean()),
        "macro_min_auc": float(per_type["auc"].min()),
    }])
    return per_type, pooled_summary, predictions.sort_values("sample").reset_index(drop=True)


def verify_existing(per_type, pooled, expected_dir: Path, tolerance: float) -> None:
    if not np.isfinite(tolerance) or tolerance < 0:
        raise ValueError("LOCO verification tolerance must be finite and non-negative")
    expected_dir = Path(expected_dir).resolve()
    expected_per_type = pd.read_csv(
        expected_dir / "loco_per_cancer_metrics.csv", keep_default_na=False
    )
    expected_pooled = pd.read_csv(
        expected_dir / "loco_pooled_summary.csv", keep_default_na=False
    )
    per_type_columns = [
        "held_out_type",
        "model",
        "n_test",
        "n_normal",
        "n_tumor",
        "auc",
        "avg_precision",
        "accuracy",
        "f1",
        "precision",
        "recall",
    ]
    pooled_columns = [
        "model",
        "pooled_auc",
        "pooled_avg_precision",
        "pooled_accuracy",
        "pooled_f1",
        "macro_mean_auc",
        "macro_min_auc",
    ]
    for name, frame, columns in [
        ("actual per-type", per_type, per_type_columns),
        ("expected per-type", expected_per_type, per_type_columns),
        ("actual pooled", pooled, pooled_columns),
        ("expected pooled", expected_pooled, pooled_columns),
    ]:
        if list(frame.columns) != columns:
            raise ValueError(
                f"LOCO {name} schema differs; expected columns {columns}, "
                f"found {list(frame.columns)}"
            )
        if frame.empty:
            raise ValueError(f"LOCO {name} table must not be empty")
    if len(pooled) != 1 or len(expected_pooled) != 1:
        raise ValueError("LOCO pooled tables must contain exactly one row")

    actual = per_type.sort_values("held_out_type").reset_index(drop=True)
    expected = expected_per_type.sort_values("held_out_type").reset_index(drop=True)
    if actual["held_out_type"].duplicated().any() or expected["held_out_type"].duplicated().any():
        raise ValueError("LOCO held_out_type rows must be unique")
    if not actual[["held_out_type", "model"]].equals(
        expected[["held_out_type", "model"]]
    ):
        raise ValueError("LOCO cancer-type identities or model names differ")
    if not actual["model"].eq("logistic_regression").all():
        raise ValueError("LOCO per-type model must be logistic_regression")
    if not pooled["model"].eq("logistic_regression").all() or not expected_pooled[
        "model"
    ].eq("logistic_regression").all():
        raise ValueError("LOCO pooled model must be logistic_regression")

    per_numeric = per_type_columns[2:]
    pooled_numeric = pooled_columns[1:]
    try:
        actual_values = actual[per_numeric].to_numpy(dtype=float)
        expected_values = expected[per_numeric].to_numpy(dtype=float)
        actual_pooled = pooled[pooled_numeric].to_numpy(dtype=float)
        expected_pooled_values = expected_pooled[pooled_numeric].to_numpy(dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("LOCO verification metric columns must be numeric") from exc
    for name, values in [
        ("actual per-type", actual_values),
        ("expected per-type", expected_values),
        ("actual pooled", actual_pooled),
        ("expected pooled", expected_pooled_values),
    ]:
        if not np.isfinite(values).all():
            raise ValueError(f"LOCO {name} metrics contain NaN or infinite values")
    if actual_values.shape != expected_values.shape:
        raise ValueError("LOCO per-type row counts differ")
    delta = np.abs(actual_values - expected_values)
    if not np.allclose(
        actual_values,
        expected_values,
        rtol=0,
        atol=tolerance,
        equal_nan=False,
    ):
        raise ValueError(
            f"LOCO per-type metrics differ from committed artifacts; max delta={delta.max():.3g}"
        )
    pooled_delta = np.abs(actual_pooled - expected_pooled_values)
    if not np.allclose(
        actual_pooled,
        expected_pooled_values,
        rtol=0,
        atol=tolerance,
        equal_nan=False,
    ):
        raise ValueError(
            "LOCO pooled metrics differ from committed artifacts; "
            f"max delta={pooled_delta.max():.3g}"
        )


def write_outputs(
    output_dir,
    per_type,
    pooled,
    predictions,
    input_snapshot,
    config,
    code,
    environment,
    verification_reference,
):
    output_dir.mkdir(parents=True, exist_ok=True)
    output_paths = {
        "per_type": output_dir / "loco_per_cancer_metrics.csv",
        "pooled": output_dir / "loco_pooled_summary.csv",
        "predictions": output_dir / "loco_oof_predictions.csv",
    }
    _atomic_write_csv(per_type, output_paths["per_type"])
    _atomic_write_csv(pooled, output_paths["pooled"])
    _atomic_write_csv(predictions, output_paths["predictions"])
    manifest = {
        "schema_version": "2.0",
        "analysis": "leave-one-cancer-type-out",
        "git_commit": code["git_commit"],
        "code": code,
        "environment": environment,
        "config": config,
        "inputs": input_snapshot,
        "feature_bundle_verified": "feature_manifest" in input_snapshot,
        "verification_reference": verification_reference,
        "outputs": output_records(output_paths),
    }
    manifest_path = output_dir / "loco_run_manifest.json"
    _atomic_write_text(
        manifest_path,
        json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n",
    )
    verify_input_snapshot(input_snapshot)
    if verification_reference is not None:
        verify_input_snapshot(verification_reference)
    return output_paths | {"manifest": manifest_path}


def loco_config(n_features, c_value):
    from sklearn.linear_model import LogisticRegression

    classifier = LogisticRegression(
        max_iter=5000,
        C=float(c_value),
        class_weight="balanced",
        random_state=42,
        solver="lbfgs",
        tol=1e-4,
        fit_intercept=True,
    )
    return {
        "pipeline_order": [
            "SelectKBest(f_classif)",
            "StandardScaler",
            "LogisticRegression",
        ],
        "feature_selection": {"score_func": "f_classif", "k": int(n_features)},
        "scaler": {"with_mean": True, "with_std": True},
        "classifier": classifier.get_params(deep=False),
        "held_out_group": "project",
        "patient_identifier": "case_id",
        "threshold": 0.5,
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", type=Path, default=DEFAULT_FEATURES)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument(
        "--force", action="store_true", help="replace an existing output generation"
    )
    parser.add_argument(
        "--allow-unverified-features",
        action="store_true",
        help="allow canonical-dtype arrays without an export manifest (development only)",
    )
    parser.add_argument("--n-features", type=int, default=2000)
    parser.add_argument("--c", dest="c_value", type=float, default=0.1)
    parser.add_argument(
        "--verify-existing",
        type=Path,
        default=None,
        metavar="DIR",
        help="compare recomputed metrics with committed LOCO CSVs in DIR",
    )
    parser.add_argument("--tolerance", type=float, default=1e-10)
    args = parser.parse_args(argv)
    if not np.isfinite(args.tolerance) or args.tolerance < 0:
        parser.error("--tolerance must be finite and non-negative")

    try:
        X, genes, samples, y, projects, cases, input_snapshot = load_inputs(
            args.features,
            args.metadata,
            allow_unverified_features=args.allow_unverified_features,
        )
        verification_reference = None
        if args.verify_existing:
            expected_dir = args.verify_existing.resolve()
            verification_reference = snapshot_inputs({
                "expected_per_type": expected_dir / "loco_per_cancer_metrics.csv",
                "expected_pooled": expected_dir / "loco_pooled_summary.csv",
            })
        per_type, pooled, predictions = run_loco(
            X,
            samples,
            y,
            projects,
            cases,
            n_features=args.n_features,
            c_value=args.c_value,
        )
        if args.verify_existing:
            verify_existing(per_type, pooled, args.verify_existing, args.tolerance)
        code = code_provenance(ROOT, SOURCE_PATHS)
        environment = training_environment()
        if not environment["canonical_match"]:
            warnings.warn(
                "training environment differs from the canonical stack; "
                "exact LOCO verification may fail",
                RuntimeWarning,
                stacklevel=2,
            )
        config = loco_config(args.n_features, args.c_value)
        if args.output_dir:
            protected = [Path(record["path"]) for record in input_snapshot.values()]
            protected.extend(SOURCE_PATHS)
            if verification_reference is not None:
                protected.extend(
                    Path(record["path"]) for record in verification_reference.values()
                )
            with staged_output_directory(
                args.output_dir,
                force=args.force,
                protected_paths=protected,
            ) as stage:
                staged_paths = write_outputs(
                    stage,
                    per_type,
                    pooled,
                    predictions,
                    input_snapshot,
                    config,
                    code,
                    environment,
                    verification_reference,
                )
            final_output = args.output_dir.resolve()
            paths = {
                name: final_output / path.name for name, path in staged_paths.items()
            }
            for name, path in paths.items():
                print(f"[loco] {name}: {path}")
        else:
            verify_input_snapshot(input_snapshot)
            if verification_reference is not None:
                verify_input_snapshot(verification_reference)
    except (OSError, ValueError, RuntimeError, TypeError, KeyError) as exc:
        parser.error(str(exc))

    print(pooled.to_string(index=False))
    print(f"[loco] verified {len(per_type)} held-out cancer types / {len(predictions)} samples")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
