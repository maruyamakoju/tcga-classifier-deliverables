"""Contract tests for the canonical binary-model training pipeline."""

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "train_classifier.py"
SPEC = importlib.util.spec_from_file_location("train_classifier", SCRIPT)
TRAIN = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(TRAIN)


def write_inputs(tmp_path, *, dtype=np.float64, leaked_case=False):
    rng = np.random.default_rng(91)
    samples = np.array([f"sample-{index:02d}" for index in range(20)])
    genes = np.array([f"gene-{index:02d}" for index in range(10)])
    labels = np.tile(["normal", "tumor"], 10)
    X = rng.normal(size=(len(samples), len(genes))).astype(dtype)
    X[labels == "tumor", :3] += 4
    cases = np.array([f"case-{index:02d}" for index in range(len(samples))])
    if leaked_case:
        cases[12] = cases[0]

    features = tmp_path / "X_full_float64.npy"
    np.save(features, X)
    np.save(tmp_path / "X_genes.npy", genes)
    np.save(tmp_path / "X_samples.npy", samples)
    np.save(tmp_path / "train.npy", np.arange(12))
    np.save(tmp_path / "test.npy", np.arange(12, 20))
    metadata = pd.DataFrame(
        {
            "file_id": samples,
            "project": np.repeat("TCGA-SYNTHETIC", len(samples)),
            "label": labels,
            "case_id": cases,
            "submitter_id": [f"submitter-{index:02d}" for index in range(len(samples))],
        }
    ).sample(frac=1, random_state=3)
    metadata_path = tmp_path / "metadata.csv"
    metadata.to_csv(metadata_path, index=False)
    return features, metadata_path, tmp_path / "train.npy", tmp_path / "test.npy"


def test_load_fit_and_evaluate_preserve_canonical_sample_ids(tmp_path):
    paths = write_inputs(tmp_path)
    X, genes, samples, y, projects, _, train_index, test_index, _ = (
        TRAIN.load_training_inputs(*paths, allow_unverified_features=True)
    )
    selector, scaler, classifier = TRAIN.fit_release_model(
        X, y, train_index, n_features=5
    )
    metrics, per_project, predictions = TRAIN.heldout_evaluation(
        selector,
        scaler,
        classifier,
        X,
        samples,
        y,
        projects,
        test_index,
    )

    assert len(genes) == 10
    assert predictions["sample"].tolist() == samples[test_index].tolist()
    assert predictions["sample"].is_unique
    assert predictions["tumor_probability"].between(0, 1).all()
    assert np.isfinite(list(metrics.values())).all()
    assert per_project["n"].tolist() == [len(test_index)]


def test_load_training_inputs_requires_exact_float64_matrix(tmp_path):
    paths = write_inputs(tmp_path, dtype=np.float32)
    with pytest.raises(ValueError, match="requires the exact float64 matrix"):
        TRAIN.load_training_inputs(*paths, allow_unverified_features=True)


def test_load_training_inputs_rejects_patient_leakage(tmp_path):
    paths = write_inputs(tmp_path, leaked_case=True)
    with pytest.raises(ValueError, match="one-to-one patient mapping|patient leakage"):
        TRAIN.load_training_inputs(*paths, allow_unverified_features=True)


def test_verify_shipped_rejects_selected_gene_drift(tmp_path):
    shipped = tmp_path / "shipped.npz"
    np.savez(
        shipped,
        selected_genes=np.array(["g1", "g2"]),
        scaler_mean=np.zeros(2),
        scaler_scale=np.ones(2),
        coef=np.ones(2),
        intercept=np.array(0.0),
        class_order=np.array([0, 1]),
    )
    arrays = {
        "selected_genes": np.array(["g1", "different"]),
        "scaler_mean": np.zeros(2),
        "scaler_scale": np.ones(2),
        "coef": np.ones(2),
        "intercept": np.array(0.0),
        "class_order": np.array([0, 1]),
    }
    with pytest.raises(ValueError, match="selected genes differ"):
        TRAIN.verify_shipped(arrays, shipped, tolerance=0)


@pytest.mark.parametrize("malformation", ["broadcast", "nan", "wrong_classes"])
def test_verify_shipped_rejects_malformed_model_schema(tmp_path, malformation):
    arrays = {
        "selected_genes": np.array(["g1", "g2"]),
        "selected_gene_index": np.array([0, 1], dtype=np.int32),
        "scaler_mean": np.zeros(2),
        "scaler_scale": np.ones(2),
        "coef": np.ones(2),
        "intercept": np.array(0.0),
        "class_order": np.array([0, 1]),
    }
    shipped = dict(arrays)
    if malformation == "broadcast":
        shipped["coef"] = np.ones(1)
    elif malformation == "nan":
        shipped["scaler_mean"] = np.array([np.nan, 0.0])
    else:
        shipped["class_order"] = np.array([1, 0])
    path = tmp_path / "malformed.npz"
    np.savez(path, **shipped)
    with pytest.raises(ValueError):
        TRAIN.verify_shipped(arrays, path, tolerance=0)


def test_load_training_inputs_requires_manifest_by_default(tmp_path):
    paths = write_inputs(tmp_path)
    with pytest.raises(ValueError, match="export manifest is required"):
        TRAIN.load_training_inputs(*paths)


def test_load_training_inputs_rejects_non_bijective_patient_ids(tmp_path):
    paths = write_inputs(tmp_path)
    metadata = pd.read_csv(paths[1], dtype=str)
    metadata.loc[1, "submitter_id"] = metadata.loc[0, "submitter_id"]
    metadata.to_csv(paths[1], index=False)
    with pytest.raises(ValueError, match="one-to-one patient mapping"):
        TRAIN.load_training_inputs(*paths, allow_unverified_features=True)


def test_load_training_inputs_rejects_multidimensional_axis(tmp_path):
    paths = write_inputs(tmp_path)
    np.save(tmp_path / "X_genes.npy", np.array([["g"]] * 10))
    with pytest.raises(ValueError, match="one-dimensional arrays"):
        TRAIN.load_training_inputs(*paths, allow_unverified_features=True)


def test_binary_output_is_fresh_hash_bound_generation(tmp_path):
    paths = write_inputs(tmp_path)
    output = tmp_path / "binary-output"
    argv = [
        "--features",
        str(paths[0]),
        "--metadata",
        str(paths[1]),
        "--train-index",
        str(paths[2]),
        "--test-index",
        str(paths[3]),
        "--n-features",
        "5",
        "--skip-cv",
        "--allow-unverified-features",
        "--output-dir",
        str(output),
    ]
    assert TRAIN.main(argv) == 0
    assert (output / "binary_lr_run_manifest.json").is_file()
    per_type = pd.read_csv(output / "binary_lr_per_cancer_type_performance.csv")
    assert list(per_type.columns) == [
        "project",
        "n",
        "n_tumor",
        "n_normal",
        "auc",
        "accuracy",
    ]
    (output / "stale.txt").write_text("stale", encoding="utf-8")
    assert TRAIN.main([*argv, "--force"]) == 0
    assert not (output / "stale.txt").exists()
    assert not (output / "binary_lr_grouped_cv_metrics.csv").exists()
