#!/usr/bin/env python3
"""Reproduce the deployed tumor-vs-normal logistic-regression model.

This is the canonical training path for the public NumPy model.  It consumes
the exact float64 feature export, validates the patient-disjoint split, fits
feature selection and scaling on training samples only, evaluates the held-out
test set and grouped CV, and optionally writes a provenance-rich NPZ release
candidate.  Historical RF/XGBoost comparisons remain reported artifacts but
are intentionally outside the deployed model's canonical reproduction path.
"""

from __future__ import annotations

import argparse
import json
import os
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.exceptions import ConvergenceWarning
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.preprocessing import StandardScaler

from tcga_rnaseq import metrics as M
from tcga_rnaseq import write_dataframe_csv, write_json
from training_tools import (
    code_provenance,
    output_records,
    snapshot_inputs,
    staged_output_directory,
    training_environment,
    validate_feature_manifest,
    verify_input_snapshot,
)


ROOT = Path(__file__).resolve().parent
DEFAULT_FEATURES = ROOT / "cancer-type-classifier" / "X_full_float64.npy"
SOURCE_PATHS = [
    Path(__file__).resolve(),
    ROOT / "training_tools" / "__init__.py",
    ROOT / "tcga_rnaseq" / "metrics.py",
]


def load_training_inputs(
    features_path: Path,
    metadata_path: Path,
    train_index_path: Path,
    test_index_path: Path,
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
        "train_index": train_index_path.resolve(),
        "test_index": test_index_path.resolve(),
    }
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise ValueError("missing training inputs: " + ", ".join(missing))

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
    train_index = np.load(train_index_path, allow_pickle=False)
    test_index = np.load(test_index_path, allow_pickle=False)
    if X.dtype != np.float64:
        raise ValueError(
            f"canonical binary training requires the exact float64 matrix, found {X.dtype}"
        )
    if genes.ndim != 1 or samples.ndim != 1:
        raise ValueError("X_genes.npy and X_samples.npy must be one-dimensional arrays")
    if X.ndim != 2 or X.shape != (len(samples), len(genes)):
        raise ValueError(
            f"feature shape mismatch: X={X.shape}, samples={len(samples)}, genes={len(genes)}"
        )
    if X.shape[0] == 0 or X.shape[1] == 0 or not np.isfinite(X).all():
        raise ValueError("feature matrix must be non-empty and finite")
    for name, values in [("sample", samples), ("gene", genes)]:
        if any(not value.strip() or value != value.strip() for value in values):
            raise ValueError(f"{name} identifiers must be non-empty and unpadded")
        if len(set(values)) != len(values):
            raise ValueError(f"{name} identifiers must be unique")

    for name, index in [("train", train_index), ("test", test_index)]:
        if index.ndim != 1 or index.dtype.kind not in "iu":
            raise ValueError(f"{name} index must be a one-dimensional integer array")
        if len(np.unique(index)) != len(index):
            raise ValueError(f"{name} index contains duplicates")
        if len(index) == 0 or index.min() < 0 or index.max() >= len(samples):
            raise ValueError(f"{name} index is empty or out of range")
    if set(train_index) & set(test_index):
        raise ValueError("train and test indices overlap")
    if set(train_index) | set(test_index) != set(range(len(samples))):
        raise ValueError("train and test indices must partition every sample exactly once")

    metadata = pd.read_csv(metadata_path, keep_default_na=False, dtype=str)
    required = {"file_id", "project", "label", "case_id", "submitter_id"}
    absent = sorted(required - set(metadata.columns))
    if absent:
        raise ValueError("metadata is missing columns: " + ", ".join(absent))
    if metadata["file_id"].duplicated().any():
        raise ValueError("metadata contains duplicate file_id values")
    for field in required:
        values = metadata[field].astype(str)
        if values.str.strip().eq("").any() or values.ne(values.str.strip()).any():
            raise ValueError(f"metadata column {field} contains blank or padded values")
    aligned = metadata.set_index("file_id").reindex(samples)
    fields = ["project", "label", "case_id", "submitter_id"]
    if aligned[fields].isna().any(axis=None) or aligned[fields].eq("").any(axis=None):
        raise ValueError("metadata does not completely cover the feature matrix")
    unknown_labels = sorted(set(aligned["label"]) - {"tumor", "normal"})
    if unknown_labels:
        raise ValueError(f"unsupported labels: {unknown_labels[:5]}")
    y = (aligned["label"].to_numpy() == "tumor").astype(int)
    projects = aligned["project"].astype(str).to_numpy()
    cases = aligned["case_id"].astype(str).to_numpy()
    submitters = aligned["submitter_id"].astype(str).to_numpy()
    case_submitter = aligned[["case_id", "submitter_id"]].drop_duplicates()
    if (
        case_submitter.groupby("case_id")["submitter_id"].nunique().gt(1).any()
        or case_submitter.groupby("submitter_id")["case_id"].nunique().gt(1).any()
    ):
        raise ValueError("case_id and submitter_id must have a one-to-one patient mapping")
    if aligned.groupby("case_id")["project"].nunique().gt(1).any():
        raise ValueError("each case_id must belong to exactly one project")
    for group_name, groups in [("case_id", cases), ("submitter_id", submitters)]:
        overlap = set(groups[train_index]) & set(groups[test_index])
        if overlap:
            raise ValueError(
                f"patient leakage across train/test by {group_name}: {sorted(overlap)[:5]}"
            )
    for split_name, index in [("train", train_index), ("test", test_index)]:
        if np.unique(y[index]).size != 2:
            raise ValueError(f"{split_name} split must contain both binary classes")
    verify_input_snapshot(input_snapshot)
    return (
        X,
        genes,
        samples,
        y,
        projects,
        cases,
        train_index,
        test_index,
        input_snapshot,
    )


def fit_checked(estimator, X, y, context):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        estimator.fit(X, y)
    convergence = [warning for warning in caught if issubclass(warning.category, ConvergenceWarning)]
    if convergence:
        raise RuntimeError(f"{context} did not converge: {convergence[0].message}")
    for warning in caught:
        warnings.warn_explicit(
            str(warning.message), warning.category, warning.filename, warning.lineno
        )
    return estimator


def fit_release_model(X, y, train_index, n_features=2000, c_value=0.1, seed=42):
    selector = SelectKBest(f_classif, k=n_features)
    selected = selector.fit_transform(X[train_index], y[train_index])
    scaler = StandardScaler().fit(selected)
    classifier = LogisticRegression(
        max_iter=5000,
        C=c_value,
        class_weight="balanced",
        random_state=seed,
        solver="lbfgs",
        tol=1e-4,
        fit_intercept=True,
    )
    fit_checked(
        classifier,
        scaler.transform(selected),
        y[train_index],
        "release logistic regression",
    )
    return selector, scaler, classifier


def heldout_evaluation(
    selector, scaler, classifier, X, samples, y, projects, test_index
):
    transformed = scaler.transform(selector.transform(X[test_index]))
    probability = classifier.predict_proba(transformed)[:, 1]
    if not np.isfinite(probability).all():
        raise RuntimeError("held-out evaluation produced non-finite probabilities")
    metrics = M.classification_metrics(y[test_index], probability, threshold=0.5)
    metrics["average_precision"] = M.average_precision(y[test_index], probability)
    rows = []
    for project in sorted(np.unique(projects[test_index])):
        mask = projects[test_index] == project
        labels = y[test_index][mask]
        scores = probability[mask]
        if np.unique(labels).size != 2:
            raise ValueError(f"held-out project {project} must contain both classes")
        project_metrics = M.classification_metrics(labels, scores, threshold=0.5)
        rows.append({
            "project": project,
            "n": int(mask.sum()),
            "n_tumor": int((labels == 1).sum()),
            "n_normal": int((labels == 0).sum()),
            "auc": project_metrics["auc"],
            "accuracy": project_metrics["accuracy"],
        })
    predictions = pd.DataFrame({
        "sample": samples[test_index],
        "label": y[test_index],
        "project": projects[test_index],
        "tumor_probability": probability,
        "call": np.where(probability >= 0.5, "tumor", "normal"),
    })
    return metrics, pd.DataFrame(rows), predictions


def grouped_cv(
    X,
    y,
    projects,
    cases,
    train_index,
    n_features,
    c_value,
    *,
    n_folds=5,
    cv_seed=1,
    model_seed=42,
):
    splitter = StratifiedGroupKFold(
        n_splits=n_folds, shuffle=True, random_state=cv_seed
    )
    stratification = np.array(
        [
            f"{project}_{label}"
            for project, label in zip(
                projects[train_index], y[train_index], strict=True
            )
        ]
    )
    rows = []
    for fold, (train_rel, validation_rel) in enumerate(
        splitter.split(
            np.empty(len(train_index)), stratification, groups=cases[train_index]
        ),
        start=1,
    ):
        fold_train_index = train_index[train_rel]
        fold_validation_index = train_index[validation_rel]
        if set(cases[fold_train_index]) & set(cases[fold_validation_index]):
            raise RuntimeError(f"patient leakage in grouped CV fold {fold}")
        selector = SelectKBest(f_classif, k=n_features)
        train_selected = selector.fit_transform(
            X[fold_train_index], y[fold_train_index]
        )
        validation_selected = selector.transform(X[fold_validation_index])
        scaler = StandardScaler().fit(train_selected)
        classifier = LogisticRegression(
            max_iter=5000,
            C=c_value,
            class_weight="balanced",
            random_state=model_seed,
            solver="lbfgs",
            tol=1e-4,
            fit_intercept=True,
        )
        fit_checked(
            classifier,
            scaler.transform(train_selected),
            y[fold_train_index],
            f"grouped CV fold {fold}",
        )
        probability = classifier.predict_proba(scaler.transform(validation_selected))[:, 1]
        if np.unique(y[fold_validation_index]).size != 2:
            raise ValueError(f"grouped CV fold {fold} does not contain both classes")
        if not np.isfinite(probability).all():
            raise RuntimeError(f"grouped CV fold {fold} produced non-finite probabilities")
        rows.append({
            "fold": fold,
            "n_train": int(len(train_rel)),
            "n_validation": int(len(validation_rel)),
            "auc": M.roc_auc(y[fold_validation_index], probability),
        })
    return pd.DataFrame(rows)


def build_export_arrays(selector, scaler, classifier, genes, provenance):
    selected_indices = np.flatnonzero(selector.get_support())
    class_order = np.asarray(classifier.classes_)
    if not np.array_equal(class_order, np.array([0, 1])):
        raise RuntimeError(f"unexpected binary classifier class order: {class_order.tolist()}")
    return {
        "selected_genes": genes[selected_indices].astype(str),
        "selected_gene_index": selected_indices.astype(np.int32),
        "scaler_mean": scaler.mean_,
        "scaler_scale": scaler.scale_,
        "coef": classifier.coef_[0],
        "intercept": classifier.intercept_[0],
        "class_order": class_order.astype(np.int8),
        "model_metadata_json": json.dumps(provenance, sort_keys=True),
        "notes": (
            "Pure NumPy logistic-regression export. Input must be GDC STAR-Counts "
            "log2(TPM+1), rows=samples, columns=Ensembl genes."
        ),
    }


def write_npz(path: Path, arrays) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        with temporary.open("wb") as handle:
            np.savez(handle, **arrays)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _validated_binary_model_arrays(values, *, source: str) -> dict[str, np.ndarray]:
    required = {
        "selected_genes",
        "scaler_mean",
        "scaler_scale",
        "coef",
        "intercept",
        "class_order",
    }
    available = set(values.files) if hasattr(values, "files") else set(values)
    missing = sorted(required - available)
    if missing:
        raise ValueError(f"{source} model is missing keys: {', '.join(missing)}")
    arrays = {key: np.asarray(values[key]) for key in required}
    genes = arrays["selected_genes"]
    if genes.ndim != 1 or genes.size == 0:
        raise ValueError(f"{source} selected_genes must be a non-empty 1D array")
    genes = genes.astype(str)
    if (
        len(set(genes)) != len(genes)
        or any(not value.strip() or value != value.strip() for value in genes)
    ):
        raise ValueError(f"{source} selected_genes must be unique, non-blank, and unpadded")
    n_features = len(genes)
    for key in ["scaler_mean", "scaler_scale", "coef"]:
        array = arrays[key]
        if array.shape != (n_features,) or array.dtype.kind not in "fiu":
            raise ValueError(
                f"{source} {key} must have numeric shape ({n_features},), found {array.shape}"
            )
        if not np.isfinite(array).all():
            raise ValueError(f"{source} {key} contains NaN or infinite values")
    intercept = arrays["intercept"]
    if intercept.shape != () or intercept.dtype.kind not in "fiu" or not np.isfinite(intercept):
        raise ValueError(f"{source} intercept must be one finite numeric scalar")
    if np.any(arrays["scaler_scale"] <= 0):
        raise ValueError(f"{source} scaler_scale must be strictly positive")
    class_order = arrays["class_order"]
    if class_order.shape != (2,) or not np.array_equal(
        class_order.astype(float), np.array([0.0, 1.0])
    ):
        raise ValueError(f"{source} class_order must be exactly [0, 1]")
    arrays["selected_genes"] = genes
    if "selected_gene_index" in available:
        indices = np.asarray(values["selected_gene_index"])
        if (
            indices.shape != (n_features,)
            or indices.dtype.kind not in "iu"
            or np.any(indices < 0)
            or len(np.unique(indices)) != n_features
            or np.any(np.diff(indices) <= 0)
        ):
            raise ValueError(
                f"{source} selected_gene_index must be a strictly increasing integer vector"
            )
        arrays["selected_gene_index"] = indices
    return arrays


def verify_shipped(arrays, shipped_path: Path, tolerance):
    if not np.isfinite(tolerance) or tolerance < 0:
        raise ValueError("model verification tolerance must be finite and non-negative")
    candidate = _validated_binary_model_arrays(arrays, source="retrained")
    try:
        with np.load(shipped_path, allow_pickle=False) as loaded:
            shipped = _validated_binary_model_arrays(loaded, source="shipped")
    except (OSError, ValueError, KeyError) as exc:
        raise ValueError(f"could not validate shipped model {shipped_path}: {exc}") from exc
    if not np.array_equal(candidate["selected_genes"], shipped["selected_genes"]):
        raise ValueError("selected genes differ from the shipped model")
    if not np.array_equal(candidate["class_order"], shipped["class_order"]):
        raise ValueError("class order differs from the shipped model")
    if "selected_gene_index" in shipped:
        if "selected_gene_index" not in candidate or not np.array_equal(
            candidate["selected_gene_index"], shipped["selected_gene_index"]
        ):
            raise ValueError("selected gene indices differ from the shipped model")
    deltas = {}
    for key in ["scaler_mean", "scaler_scale", "coef", "intercept"]:
        if candidate[key].shape != shipped[key].shape:
            raise ValueError(
                f"{key} shape differs from shipped model: "
                f"{candidate[key].shape} != {shipped[key].shape}"
            )
        delta = np.abs(candidate[key].astype(float) - shipped[key].astype(float))
        if not np.isfinite(delta).all():
            raise ValueError(f"non-finite {key} comparison against shipped model")
        deltas[key] = float(delta.max()) if delta.size else float(delta)
    worst = max(deltas.values())
    if not np.isfinite(worst) or worst > tolerance:
        raise ValueError(
            f"retrained model differs from shipped weights; max delta={worst:.3g}, "
            f"tolerance={tolerance:.3g}"
        )
    return deltas


def training_config(classifier, *, n_features, c_value, skip_cv):
    return {
        "pipeline_order": [
            "SelectKBest(f_classif)",
            "StandardScaler",
            "LogisticRegression",
        ],
        "feature_selection": {"score_func": "f_classif", "k": int(n_features)},
        "scaler": {"with_mean": True, "with_std": True},
        "classifier": classifier.get_params(deep=False),
        "heldout_threshold": 0.5,
        "grouped_cv": {
            "enabled": not skip_cv,
            "splitter": "StratifiedGroupKFold",
            "n_splits": 5,
            "shuffle": True,
            "random_state": 1,
            "group": "case_id",
            "stratification": "project + binary label",
            "feature_selection_fit_inside_fold": True,
        },
        "requested_C": float(c_value),
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", type=Path, default=DEFAULT_FEATURES)
    parser.add_argument("--metadata", type=Path, default=ROOT / "selected_files.csv")
    parser.add_argument("--train-index", type=Path, default=ROOT / "train_idx.npy")
    parser.add_argument("--test-index", type=Path, default=ROOT / "test_idx.npy")
    parser.add_argument("--n-features", type=int, default=2000)
    parser.add_argument("--c", dest="c_value", type=float, default=0.1)
    parser.add_argument("--skip-cv", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument(
        "--force", action="store_true", help="replace an existing output generation"
    )
    parser.add_argument(
        "--allow-unverified-features",
        action="store_true",
        help="allow canonical-dtype arrays without an export manifest (development only)",
    )
    parser.add_argument("--verify-shipped", type=Path, default=None, metavar="NPZ")
    parser.add_argument("--weight-tolerance", type=float, default=1e-7)
    args = parser.parse_args(argv)
    if args.n_features < 1 or not np.isfinite(args.c_value) or args.c_value <= 0:
        parser.error("--n-features must be positive and --c must be finite and positive")
    if not np.isfinite(args.weight_tolerance) or args.weight_tolerance < 0:
        parser.error("--weight-tolerance must be finite and non-negative")

    try:
        X, genes, samples, y, projects, cases, train_index, test_index, input_snapshot = (
            load_training_inputs(
                args.features,
                args.metadata,
                args.train_index,
                args.test_index,
                allow_unverified_features=args.allow_unverified_features,
            )
        )
        if args.n_features > X.shape[1]:
            raise ValueError(f"--n-features exceeds input gene count {X.shape[1]}")
        verification_snapshot = None
        if args.verify_shipped:
            verification_snapshot = snapshot_inputs(
                {"shipped_model": args.verify_shipped.resolve()}
            )
        selector, scaler, classifier = fit_release_model(
            X, y, train_index, args.n_features, args.c_value
        )
        heldout, per_project, predictions = heldout_evaluation(
            selector, scaler, classifier, X, samples, y, projects, test_index
        )
        cv = None if args.skip_cv else grouped_cv(
            X, y, projects, cases, train_index, args.n_features, args.c_value
        )
        code = code_provenance(ROOT, SOURCE_PATHS)
        environment = training_environment()
        if not environment["canonical_match"]:
            warnings.warn(
                "training environment differs from the canonical stack; "
                "strict shipped-model verification may fail",
                RuntimeWarning,
                stacklevel=2,
            )
        provenance = {
            "schema_version": "3.0",
            "task": "TCGA/GDC tumor-vs-adjacent-normal classification",
            "git_commit": code["git_commit"],
            "code": code,
            "config": training_config(
                classifier,
                n_features=args.n_features,
                c_value=args.c_value,
                skip_cv=args.skip_cv,
            ),
            "inputs": {
                name: record["sha256"] for name, record in input_snapshot.items()
            },
            "input_files": input_snapshot,
            "feature_bundle_verified": "feature_manifest" in input_snapshot,
            "environment": environment,
        }
        arrays = build_export_arrays(selector, scaler, classifier, genes, provenance)
        deltas = None
        if args.verify_shipped:
            deltas = verify_shipped(arrays, args.verify_shipped.resolve(), args.weight_tolerance)
        if verification_snapshot is not None:
            verify_input_snapshot(verification_snapshot)

        summary = {
            "schema_version": "3.0",
            "n_train": int(len(train_index)),
            "n_test": int(len(test_index)),
            "heldout": heldout,
            "grouped_cv": None if cv is None else {
                "folds": cv.to_dict(orient="records"),
                "mean_auc": float(cv["auc"].mean()),
                "std_auc": float(cv["auc"].std(ddof=0)),
            },
            "weight_deltas_vs_shipped": deltas,
            "verification_reference": verification_snapshot,
            "provenance": provenance,
        }
        if args.output_dir:
            protected = [Path(record["path"]) for record in input_snapshot.values()]
            protected.extend(SOURCE_PATHS)
            if verification_snapshot is not None:
                protected.extend(
                    Path(record["path"]) for record in verification_snapshot.values()
                )
            with staged_output_directory(
                args.output_dir,
                force=args.force,
                protected_paths=protected,
            ) as output:
                output_paths = {
                    "weights": output / "deployable_lr_weights.npz",
                    "per_cancer_type": (
                        output / "binary_lr_per_cancer_type_performance.csv"
                    ),
                    "heldout_predictions": output / "binary_lr_heldout_predictions.csv",
                    "summary": output / "binary_lr_training_summary.json",
                }
                if cv is not None:
                    output_paths["grouped_cv"] = output / "binary_lr_grouped_cv_metrics.csv"
                write_npz(output_paths["weights"], arrays)
                write_dataframe_csv(per_project, output_paths["per_cancer_type"])
                write_dataframe_csv(predictions, output_paths["heldout_predictions"])
                if cv is not None:
                    write_dataframe_csv(cv, output_paths["grouped_cv"])
                write_json(summary, output_paths["summary"], sort_keys=True)
                verify_input_snapshot(input_snapshot)
                if verification_snapshot is not None:
                    verify_input_snapshot(verification_snapshot)
                run_manifest = {
                    "schema_version": "1.0",
                    "analysis": "canonical-binary-logistic-regression",
                    "code": code,
                    "config": provenance["config"],
                    "environment": environment,
                    "inputs": input_snapshot,
                    "verification_reference": verification_snapshot,
                    "outputs": output_records(output_paths),
                }
                write_json(
                    run_manifest,
                    output / "binary_lr_run_manifest.json",
                    sort_keys=True,
                )
        else:
            verify_input_snapshot(input_snapshot)
    except (OSError, ValueError, RuntimeError, TypeError, KeyError) as exc:
        parser.error(str(exc))

    print(
        f"heldout AUC={heldout['auc']:.6f} accuracy={heldout['accuracy']:.6f} "
        f"n={len(test_index)}"
    )
    if cv is not None:
        print(f"grouped CV AUC={cv['auc'].mean():.6f} +/- {cv['auc'].std(ddof=0):.6f}")
    if deltas is not None:
        print(f"shipped weight reproduction max_delta={max(deltas.values()):.3g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
