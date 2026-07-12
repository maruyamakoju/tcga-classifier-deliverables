"""Contract tests for the cancer-type training/provenance pipeline."""

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "cancer-type-classifier" / "train_cancer_type_classifier.py"
SPEC = importlib.util.spec_from_file_location("train_cancer_type_classifier", SCRIPT)
TRAIN = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(TRAIN)


def synthetic_multiclass():
    rng = np.random.default_rng(7)
    labels = np.repeat(["A", "B", "C"], 9)
    groups = np.array([f"{label}-{i}" for label in ["A", "B", "C"] for i in range(9)])
    X = rng.normal(size=(27, 8))
    X[labels == "A", 0] += 3
    X[labels == "B", 1] += 3
    X[labels == "C", 2] += 3
    return X, labels, groups


def write_cancer_inputs(tmp_path, *, dtype=np.float32, labels=None, projects=None, cases=None):
    labels = labels or ["tumor", "tumor", "normal"]
    projects = projects or ["TCGA-A", "TCGA-B", "TCGA-A"]
    cases = cases or ["c1", "c2", "c3"]
    np.save(tmp_path / "X_full.npy", np.ones((3, 3), dtype=dtype))
    np.save(tmp_path / "X_genes.npy", np.array(["g1", "g2", "g3"]))
    np.save(tmp_path / "X_samples.npy", np.array(["s1", "s2", "s3"]))
    metadata = tmp_path / "selected_files.csv"
    pd.DataFrame({
        "file_id": ["s1", "s2", "s3"],
        "project": projects,
        "label": labels,
        "case_id": cases,
    }).to_csv(metadata, index=False)
    return tmp_path / "X_full.npy", metadata


def test_cross_validated_predictions_cover_every_sample_deterministically():
    X, labels, groups = synthetic_multiclass()
    first = TRAIN.cross_validated_predictions(
        X, labels, groups, n_features=5, c_value=1.0, n_folds=3, seed=11, max_iter=2000
    )
    second = TRAIN.cross_validated_predictions(
        X, labels, groups, n_features=5, c_value=1.0, n_folds=3, seed=11, max_iter=2000
    )
    assert np.array_equal(first[0], second[0])
    assert np.array_equal(first[1], second[1])
    assert np.array_equal(first[2], second[2])
    assert set(first[2]) == {1, 2, 3}
    assert np.isfinite(first[1]).all()
    pd.testing.assert_frame_equal(first[3], second[3])


def test_cross_validation_rejects_too_few_patient_groups():
    X, labels, groups = synthetic_multiclass()
    groups[labels == "C"] = "one-patient"
    with pytest.raises(ValueError, match="only 1 patients"):
        TRAIN.cross_validated_predictions(
            X, labels, groups, n_features=5, c_value=1.0, n_folds=3, seed=0, max_iter=100
        )


def test_load_training_data_rejects_incomplete_metadata(tmp_path):
    np.save(tmp_path / "X_full.npy", np.ones((2, 3), dtype=np.float32))
    np.save(tmp_path / "X_genes.npy", np.array(["g1", "g2", "g3"]))
    np.save(tmp_path / "X_samples.npy", np.array(["s1", "s2"]))
    metadata = tmp_path / "selected_files.csv"
    metadata.write_text(
        "file_id,project,label,case_id\ns1,TCGA-A,tumor,c1\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="completely cover"):
        TRAIN.load_training_data(
            tmp_path / "X_full.npy",
            metadata,
            allow_unverified_features=True,
        )


def test_load_training_data_rejects_unknown_labels(tmp_path):
    paths = write_cancer_inputs(tmp_path, labels=["tumor", "typo", "normal"])
    with pytest.raises(ValueError, match="unsupported labels"):
        TRAIN.load_training_data(*paths, allow_unverified_features=True)


def test_load_training_data_requires_canonical_float32(tmp_path):
    paths = write_cancer_inputs(tmp_path, dtype=np.float64)
    with pytest.raises(ValueError, match="exact float32"):
        TRAIN.load_training_data(*paths, allow_unverified_features=True)


def test_load_training_data_rejects_one_patient_with_multiple_targets(tmp_path):
    paths = write_cancer_inputs(
        tmp_path,
        projects=["TCGA-A", "TCGA-B", "TCGA-A"],
        cases=["same", "same", "other"],
    )
    with pytest.raises(ValueError, match="exactly one cancer type"):
        TRAIN.load_training_data(*paths, allow_unverified_features=True)


def test_verify_shipped_cancer_model_rejects_nan_and_shape_broadcast(tmp_path):
    candidate = {
        "selected_genes": np.array(["g1", "g2"]),
        "selected_gene_index": np.array([0, 1], dtype=np.int32),
        "scaler_mean": np.zeros(2),
        "scaler_scale": np.ones(2),
        "coef": np.ones((2, 2)),
        "intercept": np.zeros(2),
        "classes": np.array(["A", "B"]),
    }
    malformed = dict(candidate)
    malformed["coef"] = np.ones((2, 1))
    malformed["scaler_mean"] = np.array([np.nan, 0.0])
    path = tmp_path / "malformed.npz"
    np.savez(path, **malformed)
    with pytest.raises(ValueError):
        TRAIN.verify_shipped(candidate, path, tolerance=0)
