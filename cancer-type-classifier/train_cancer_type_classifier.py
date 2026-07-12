#!/usr/bin/env python3
"""Train, patient-held-out-evaluate, and export the cancer-type classifier.

Every out-of-fold prediction is produced by a pipeline whose scaler and feature
selector were fitted without that patient's samples.  The final deployable NPZ
is fitted on all tumor samples and records hashes of its exact training inputs,
configuration, environment, and source commit.
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
from sklearn.exceptions import ConvergenceWarning
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
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


def _atomic_path(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    os.close(descriptor)
    return Path(temporary)


def atomic_write_text(path: Path, text: str) -> None:
    temporary = _atomic_path(path)
    try:
        with temporary.open("w", encoding="utf-8", newline="\n") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def atomic_write_csv(frame: pd.DataFrame, path: Path, index=False) -> None:
    atomic_write_text(
        path,
        frame.to_csv(index=index, float_format="%.17g", lineterminator="\n"),
    )


def atomic_savez(path: Path, **arrays) -> None:
    temporary = _atomic_path(path)
    try:
        with temporary.open("wb") as handle:
            np.savez(handle, **arrays)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def load_training_data(
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
        raise ValueError("missing cancer-type training inputs: " + ", ".join(missing))

    manifest_path, _ = validate_feature_manifest(
        paths["features"],
        paths["genes"],
        paths["samples"],
        expected_dtype=np.float32,
        allow_unverified=allow_unverified_features,
    )
    if manifest_path is not None:
        paths["feature_manifest"] = manifest_path
    input_snapshot = snapshot_inputs(paths)

    X = np.load(features_path, mmap_mode="r", allow_pickle=False)
    genes = np.load(genes_path, allow_pickle=False).astype(str)
    samples = np.load(samples_path, allow_pickle=False).astype(str)
    if X.dtype != np.float32:
        raise ValueError(
            f"canonical cancer-type training requires the exact float32 matrix, found {X.dtype}"
        )
    if genes.ndim != 1 or samples.ndim != 1:
        raise ValueError("X_genes.npy and X_samples.npy must be one-dimensional arrays")
    if X.ndim != 2 or X.shape != (len(samples), len(genes)):
        raise ValueError(
            f"feature shape mismatch: X={X.shape}, samples={len(samples)}, genes={len(genes)}"
        )
    if X.shape[0] == 0 or X.shape[1] == 0 or not np.isfinite(X).all():
        raise ValueError("feature matrix must be non-empty and contain only finite values")
    for name, values in [("sample", samples), ("gene", genes)]:
        if any(not value.strip() or value != value.strip() for value in values):
            raise ValueError(f"{name} identifiers must be non-empty and unpadded")
        if len(set(values)) != len(values):
            raise ValueError(f"{name} identifiers must be unique")

    metadata = pd.read_csv(metadata_path, keep_default_na=False, dtype=str)
    required_columns = {"file_id", "project", "label", "case_id"}
    missing_columns = sorted(required_columns - set(metadata.columns))
    if missing_columns:
        raise ValueError("metadata is missing columns: " + ", ".join(missing_columns))
    if metadata["file_id"].duplicated().any():
        raise ValueError("metadata contains duplicate file_id values")
    for field in required_columns:
        values = metadata[field].astype(str)
        if values.str.strip().eq("").any() or values.ne(values.str.strip()).any():
            raise ValueError(f"metadata column {field} contains blank or padded values")
    aligned = metadata.set_index("file_id").reindex(samples)
    fields = ["project", "label", "case_id"]
    if aligned[fields].isna().any(axis=None) or aligned[fields].eq("").any(axis=None):
        raise ValueError("metadata does not completely cover feature samples")

    unknown_labels = sorted(set(aligned["label"]) - {"tumor", "normal"})
    if unknown_labels:
        raise ValueError(f"metadata contains unsupported labels: {unknown_labels[:5]}")
    tumor_mask = aligned["label"].to_numpy() == "tumor"
    if not tumor_mask.any():
        raise ValueError("metadata contains no tumor samples")
    tumor_projects = aligned.loc[tumor_mask, "project"].astype(str)
    if not tumor_projects.str.startswith("TCGA-").all() or tumor_projects.str.len().le(5).any():
        raise ValueError("tumor project identifiers must use the TCGA-<TYPE> form")
    project = tumor_projects.str.removeprefix("TCGA-").to_numpy(dtype=str)
    groups = aligned.loc[tumor_mask, "case_id"].to_numpy(dtype=str)
    tumor_samples = samples[tumor_mask]
    if any(not value for value in project) or any(not value for value in groups):
        raise ValueError("tumor project and case identifiers must be non-empty")
    group_targets = pd.DataFrame({"case_id": groups, "project": project})
    if group_targets.groupby("case_id")["project"].nunique().gt(1).any():
        raise ValueError("each tumor case_id must map to exactly one cancer type")
    verify_input_snapshot(input_snapshot)
    return (
        X[tumor_mask],
        project,
        groups,
        tumor_samples,
        genes,
        input_snapshot,
    )


def make_pipeline(n_features, c_value, max_iter, seed):
    return Pipeline([
        ("scale", StandardScaler()),
        ("select", SelectKBest(f_classif, k=n_features)),
        ("classifier", LogisticRegression(
            C=c_value,
            max_iter=max_iter,
            random_state=seed,
            solver="lbfgs",
            tol=1e-4,
            fit_intercept=True,
        )),
    ])


def fit_without_convergence_warning(estimator, X, y, context):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        estimator.fit(X, y)
    convergence = [item for item in caught if issubclass(item.category, ConvergenceWarning)]
    if convergence:
        raise RuntimeError(f"{context} did not converge: {convergence[0].message}")
    for item in caught:
        warnings.warn_explicit(
            str(item.message), item.category, item.filename, item.lineno
        )
    return estimator


def cross_validated_predictions(X, y, groups, n_features, c_value, n_folds, seed, max_iter):
    classes = np.unique(y)
    if len(classes) < 2:
        raise ValueError("cancer-type training requires at least two classes")
    for label in classes:
        n_groups = np.unique(groups[y == label]).size
        if n_groups < n_folds:
            raise ValueError(
                f"class {label} has only {n_groups} patients; {n_folds} folds requested"
            )

    splitter = StratifiedGroupKFold(n_splits=n_folds, shuffle=True, random_state=seed)
    predictions = np.empty(len(y), dtype=object)
    confidence = np.full(len(y), np.nan, dtype=float)
    assigned = np.zeros(len(y), dtype=bool)
    fold_index = np.full(len(y), -1, dtype=int)
    fold_rows = []
    for fold, (train, test) in enumerate(splitter.split(X, y, groups), start=1):
        overlap = set(groups[train]) & set(groups[test])
        if overlap:
            raise RuntimeError(f"patient leakage in fold {fold}: {sorted(overlap)[:5]}")
        if assigned[test].any():
            raise RuntimeError(f"samples assigned to more than one fold at fold {fold}")
        pipeline = make_pipeline(n_features, c_value, max_iter, seed)
        fit_without_convergence_warning(pipeline, X[train], y[train], f"fold {fold}")
        probabilities = pipeline.predict_proba(X[test])
        if not np.isfinite(probabilities).all():
            raise RuntimeError(f"fold {fold} produced non-finite probabilities")
        fold_predictions = pipeline.classes_[probabilities.argmax(axis=1)]
        predictions[test] = fold_predictions
        confidence[test] = probabilities.max(axis=1)
        fold_index[test] = fold
        assigned[test] = True
        fold_rows.append({
            "fold": fold,
            "n_train": int(len(train)),
            "n_test": int(len(test)),
            "n_train_patients": int(np.unique(groups[train]).size),
            "n_test_patients": int(np.unique(groups[test]).size),
            "accuracy": M.accuracy(y[test], fold_predictions),
            "macro_f1": M.macro_f1(y[test], fold_predictions, labels=classes),
        })
    if not assigned.all() or not np.isfinite(confidence).all():
        raise RuntimeError("out-of-fold predictions do not cover every tumor exactly once")
    return predictions.astype(str), confidence, fold_index, pd.DataFrame(fold_rows)


def fit_final_model(X, y, n_features, c_value, seed, max_iter):
    scaler = StandardScaler().fit(X)
    scaled = scaler.transform(X)
    selector = SelectKBest(f_classif, k=n_features).fit(scaled, y)
    selected_mask = selector.get_support()
    classifier = LogisticRegression(
        C=c_value,
        max_iter=max_iter,
        random_state=seed,
        solver="lbfgs",
        tol=1e-4,
        fit_intercept=True,
    )
    fit_without_convergence_warning(
        classifier, scaled[:, selected_mask], y, "final cancer-type model"
    )
    return scaler, selector, classifier


def load_gene_symbols(path: Path) -> dict:
    if not path.is_file():
        raise ValueError(f"gene symbol table not found: {path}")
    table = pd.read_csv(path, keep_default_na=False, dtype=str)
    if not {"gene_id", "symbol"}.issubset(table.columns):
        raise ValueError("gene symbol table must contain gene_id and symbol columns")
    if table["gene_id"].duplicated().any():
        raise ValueError("gene symbol table contains duplicate gene_id values")
    for field in ["gene_id", "symbol"]:
        values = table[field].astype(str)
        if values.str.strip().eq("").any() or values.ne(values.str.strip()).any():
            raise ValueError(f"gene symbol table column {field} has blank or padded values")
    return table.set_index("gene_id")["symbol"].astype(str).to_dict()


def _validated_cancer_model_arrays(values, *, source: str) -> dict[str, np.ndarray]:
    required = {
        "selected_genes",
        "selected_gene_index",
        "scaler_mean",
        "scaler_scale",
        "coef",
        "intercept",
        "classes",
    }
    available = set(values.files) if hasattr(values, "files") else set(values)
    missing = sorted(required - available)
    if missing:
        raise ValueError(f"{source} model is missing keys: {', '.join(missing)}")
    arrays = {key: np.asarray(values[key]) for key in required}
    genes = arrays["selected_genes"]
    classes = arrays["classes"]
    if genes.ndim != 1 or genes.size == 0 or classes.ndim != 1 or classes.size < 2:
        raise ValueError(f"{source} selected_genes/classes have invalid shapes")
    genes = genes.astype(str)
    classes = classes.astype(str)
    for name, entries in [("selected_genes", genes), ("classes", classes)]:
        if (
            len(set(entries)) != len(entries)
            or any(not value.strip() or value != value.strip() for value in entries)
        ):
            raise ValueError(f"{source} {name} must be unique, non-blank, and unpadded")
    n_features = len(genes)
    n_classes = len(classes)
    indices = arrays["selected_gene_index"]
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
    expected_shapes = {
        "scaler_mean": (n_features,),
        "scaler_scale": (n_features,),
        "coef": (n_classes, n_features),
        "intercept": (n_classes,),
    }
    for key, shape in expected_shapes.items():
        array = arrays[key]
        if array.shape != shape or array.dtype.kind not in "fiu":
            raise ValueError(f"{source} {key} must have numeric shape {shape}")
        if not np.isfinite(array).all():
            raise ValueError(f"{source} {key} contains NaN or infinite values")
    if np.any(arrays["scaler_scale"] <= 0):
        raise ValueError(f"{source} scaler_scale must be strictly positive")
    arrays["selected_genes"] = genes
    arrays["classes"] = classes
    return arrays


def verify_shipped(arrays, shipped_path: Path, tolerance: float) -> dict[str, float]:
    if not np.isfinite(tolerance) or tolerance < 0:
        raise ValueError("model verification tolerance must be finite and non-negative")
    candidate = _validated_cancer_model_arrays(arrays, source="retrained")
    try:
        with np.load(shipped_path, allow_pickle=False) as loaded:
            shipped = _validated_cancer_model_arrays(loaded, source="shipped")
    except (OSError, ValueError, KeyError) as exc:
        raise ValueError(f"could not validate shipped model {shipped_path}: {exc}") from exc
    for key in ["selected_genes", "selected_gene_index", "classes"]:
        if not np.array_equal(candidate[key], shipped[key]):
            raise ValueError(f"{key} differs from the shipped cancer-type model")
    deltas = {}
    for key in ["scaler_mean", "scaler_scale", "coef", "intercept"]:
        if candidate[key].shape != shipped[key].shape:
            raise ValueError(f"{key} shape differs from the shipped cancer-type model")
        delta = np.abs(candidate[key].astype(float) - shipped[key].astype(float))
        if not np.isfinite(delta).all():
            raise ValueError(f"non-finite {key} comparison against shipped model")
        deltas[key] = float(delta.max())
    worst = max(deltas.values())
    if not np.isfinite(worst) or worst > tolerance:
        raise ValueError(
            "retrained cancer-type model differs from shipped weights; "
            f"max delta={worst:.3g}, tolerance={tolerance:.3g}"
        )
    return deltas


def train_and_write(
    X,
    y,
    groups,
    samples,
    genes,
    input_snapshot,
    output_dir,
    n_features=1000,
    c_value=2.0,
    n_folds=5,
    seed=0,
    max_iter=5000,
    gene_symbols_path=None,
    code=None,
    environment=None,
    verification_reference=None,
    verify_shipped_path=None,
    weight_tolerance=1e-10,
):
    if not 1 <= n_features <= X.shape[1]:
        raise ValueError(f"n_features must be between 1 and {X.shape[1]}")
    if not np.isfinite(c_value) or c_value <= 0:
        raise ValueError("C must be finite and positive")
    if n_folds < 2 or max_iter < 1:
        raise ValueError("n_folds must be >=2 and max_iter must be >=1")

    output_dir.mkdir(parents=True, exist_ok=True)
    classes = np.unique(y)
    oof, confidence, fold_index, fold_metrics = cross_validated_predictions(
        X, y, groups, n_features, c_value, n_folds, seed, max_iter
    )
    per_class = pd.DataFrame(M.per_class_prf(y, oof, labels=classes)).rename(
        columns={"label": "cancer_type"}
    )
    matrix, matrix_labels = M.confusion_matrix(y, oof, labels=classes)
    scaler, selector, classifier = fit_final_model(
        X, y, n_features, c_value, seed, max_iter
    )
    selected_indices = np.flatnonzero(selector.get_support())
    selected_genes = genes[selected_indices]
    config = {
        "pipeline_order": [
            "StandardScaler",
            "SelectKBest(f_classif)",
            "LogisticRegression",
        ],
        "feature_selection": {"score_func": "f_classif", "k": int(n_features)},
        "scaler": {"with_mean": True, "with_std": True},
        "classifier": classifier.get_params(deep=False),
        "cross_validation": {
            "splitter": "StratifiedGroupKFold",
            "n_splits": int(n_folds),
            "shuffle": True,
            "random_state": int(seed),
            "group": "case_id",
            "pipeline_fit_inside_fold": True,
        },
    }
    code = code if code is not None else code_provenance(ROOT, SOURCE_PATHS)
    environment = environment if environment is not None else training_environment()

    output_paths = {
        "oof_predictions": output_dir / "cancer_type_oof_predictions.csv",
        "per_class_metrics": output_dir / "cancer_type_per_class_metrics.csv",
        "confusion_matrix": output_dir / "cancer_type_confusion_matrix.csv",
        "fold_metrics": output_dir / "cancer_type_fold_metrics.csv",
        "summary": output_dir / "cancer_type_summary.json",
        "top_genes": output_dir / "cancer_type_top_genes.csv",
        "weights": output_dir / "cancer_type_lr_weights.npz",
    }
    atomic_write_csv(pd.DataFrame({
        "file_id": samples,
        "case_id": groups,
        "fold": fold_index,
        "true": y,
        "pred": oof,
        "correct": (oof == y).astype(int),
        "confidence": confidence,
    }), output_paths["oof_predictions"])
    atomic_write_csv(
        per_class.sort_values("f1", ascending=False), output_paths["per_class_metrics"]
    )
    confusion_frame = pd.DataFrame(matrix, index=matrix_labels, columns=matrix_labels)
    confusion_frame.index.name = "true"
    atomic_write_text(
        output_paths["confusion_matrix"],
        confusion_frame.to_csv(float_format="%.17g", lineterminator="\n"),
    )
    atomic_write_csv(fold_metrics, output_paths["fold_metrics"])

    symbols = load_gene_symbols(gene_symbols_path) if gene_symbols_path else {}
    top_rows = []
    for class_index, class_name in enumerate(classifier.classes_):
        order = np.argsort(-classifier.coef_[class_index], kind="mergesort")[:15]
        for rank, selected_index in enumerate(order, start=1):
            gene = str(selected_genes[selected_index])
            top_rows.append({
                "cancer_type": class_name,
                "rank": rank,
                "gene": gene,
                "symbol": symbols.get(gene, ""),
                "coef": float(classifier.coef_[class_index, selected_index]),
            })
    atomic_write_csv(pd.DataFrame(top_rows), output_paths["top_genes"])

    model_metadata = {
        "schema_version": "3.0",
        "task": "TCGA cancer type classification",
        "training_config": config,
        "git_commit": code["git_commit"],
        "code": code,
        "input_sha256": {
            name: record["sha256"] for name, record in input_snapshot.items()
        },
        "feature_bundle_verified": "feature_manifest" in input_snapshot,
        "environment": environment,
    }
    model_arrays = {
        "selected_genes": selected_genes.astype(str),
        "selected_gene_index": selected_indices.astype(np.int32),
        "scaler_mean": scaler.mean_[selected_indices],
        "scaler_scale": scaler.scale_[selected_indices],
        "coef": classifier.coef_,
        "intercept": classifier.intercept_,
        "classes": np.asarray(classifier.classes_, dtype=str),
        "model_metadata_json": json.dumps(
            model_metadata, sort_keys=True, allow_nan=False
        ),
        "notes": (
            "Pure-NumPy multinomial LR cancer-type classifier. Input: GDC "
            "STAR-Counts log2(TPM+1), columns=Ensembl gene IDs."
        ),
    }
    weight_deltas = None
    if verify_shipped_path is not None:
        weight_deltas = verify_shipped(
            model_arrays, verify_shipped_path, weight_tolerance
        )
    atomic_savez(output_paths["weights"], **model_arrays)

    summary = {
        "schema_version": "3.0",
        "n_tumors": int(len(y)),
        "n_classes": int(len(classes)),
        "n_patients": int(np.unique(groups).size),
        "n_input_genes": int(X.shape[1]),
        "selected_genes": int(n_features),
        "model": (
            f"StandardScaler+SelectKBest(f_classif,k={n_features})+"
            f"LogisticRegression(C={c_value})"
        ),
        "evaluation": f"{n_folds}-fold StratifiedGroupKFold by case_id (patient-held-out)",
        "accuracy": M.accuracy(y, oof),
        "balanced_accuracy": M.balanced_accuracy(y, oof),
        "macro_f1": M.macro_f1(y, oof, labels=classes),
        "weighted_f1": float(
            np.average(per_class.set_index("cancer_type").loc[classes, "f1"], weights=[
                int((y == label).sum()) for label in classes
            ])
        ),
        "config": config,
        "git_commit": code["git_commit"],
        "code": code,
        "environment": environment,
        "inputs": input_snapshot,
        "feature_bundle_verified": "feature_manifest" in input_snapshot,
        "weight_deltas_vs_shipped": weight_deltas,
        "verification_reference": verification_reference,
    }
    atomic_write_text(
        output_paths["summary"],
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
    )
    manifest_path = output_dir / "cancer_type_run_manifest.json"
    manifest = {
        "schema_version": "1.0",
        "git_commit": code["git_commit"],
        "code": code,
        "config": config,
        "environment": environment,
        "inputs": input_snapshot,
        "verification_reference": verification_reference,
        "outputs": output_records(output_paths),
    }
    atomic_write_text(
        manifest_path,
        json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n",
    )
    verify_input_snapshot(input_snapshot)
    if verification_reference is not None:
        verify_input_snapshot(verification_reference)
    return summary, output_paths | {"manifest": manifest_path}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", type=Path, default=HERE / "X_full.npy")
    parser.add_argument("--metadata", type=Path, default=ROOT / "selected_files.csv")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--force", action="store_true", help="replace an existing output generation"
    )
    parser.add_argument(
        "--allow-unverified-features",
        action="store_true",
        help="allow canonical-dtype arrays without an export manifest (development only)",
    )
    parser.add_argument("--n-features", type=int, default=1000)
    parser.add_argument("--c", dest="c_value", type=float, default=2.0)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--max-iter", type=int, default=5000)
    symbols = parser.add_mutually_exclusive_group()
    symbols.add_argument(
        "--gene-symbols", type=Path, default=HERE / "gene_id_to_name.csv"
    )
    symbols.add_argument(
        "--no-gene-symbols",
        action="store_true",
        help="deliberately omit optional gene-symbol annotations",
    )
    parser.add_argument("--verify-shipped", type=Path, default=None, metavar="NPZ")
    parser.add_argument("--weight-tolerance", type=float, default=1e-10)
    args = parser.parse_args(argv)
    if not np.isfinite(args.weight_tolerance) or args.weight_tolerance < 0:
        parser.error("--weight-tolerance must be finite and non-negative")
    try:
        X, y, groups, samples, genes, input_snapshot = load_training_data(
            args.features,
            args.metadata,
            allow_unverified_features=args.allow_unverified_features,
        )
        gene_symbols_path = None if args.no_gene_symbols else args.gene_symbols.resolve()
        if gene_symbols_path is not None:
            input_snapshot.update(
                snapshot_inputs({"gene_symbols": gene_symbols_path})
            )
        verification_reference = None
        if args.verify_shipped:
            verification_reference = snapshot_inputs(
                {"shipped_model": args.verify_shipped.resolve()}
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
            summary, staged_paths = train_and_write(
                X,
                y,
                groups,
                samples,
                genes,
                input_snapshot,
                stage,
                n_features=args.n_features,
                c_value=args.c_value,
                n_folds=args.folds,
                seed=args.seed,
                max_iter=args.max_iter,
                gene_symbols_path=gene_symbols_path,
                code=code,
                environment=environment,
                verification_reference=verification_reference,
                verify_shipped_path=(
                    None if args.verify_shipped is None else args.verify_shipped.resolve()
                ),
                weight_tolerance=args.weight_tolerance,
            )
        final_output = args.output_dir.resolve()
        paths = {
            name: final_output / path.name for name, path in staged_paths.items()
        }
    except (OSError, ValueError, RuntimeError, TypeError, KeyError) as exc:
        parser.error(str(exc))
    print(
        f"tumors={summary['n_tumors']} classes={summary['n_classes']} "
        f"patients={summary['n_patients']} genes={summary['n_input_genes']}"
    )
    print(
        f"accuracy={summary['accuracy']:.4f} "
        f"balanced_accuracy={summary['balanced_accuracy']:.4f} "
        f"macro_f1={summary['macro_f1']:.4f}"
    )
    for name, path in paths.items():
        print(f"[cancer-type] {name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
