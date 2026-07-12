"""Tests for the reproducible LOCO analysis runner."""

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "cross-cancer-holdout" / "run_loco.py"
SPEC = importlib.util.spec_from_file_location("run_loco", SCRIPT)
LOCO = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(LOCO)


def synthetic_loco_inputs():
    rng = np.random.default_rng(42)
    projects = np.repeat(["A", "B", "C"], 8)
    y = np.tile([0, 0, 0, 0, 1, 1, 1, 1], 3)
    X = rng.normal(size=(24, 6))
    X[:, 0] += y * 4
    X[:, 1] += y * 2
    samples = np.array([f"s{i:02d}" for i in range(24)])
    cases = np.array([f"case{i:02d}" for i in range(24)])
    return X, samples, y, projects, cases


def test_run_loco_covers_each_sample_once_and_is_deterministic():
    inputs = synthetic_loco_inputs()
    first = LOCO.run_loco(*inputs, n_features=4, c_value=0.1)
    second = LOCO.run_loco(*inputs, n_features=4, c_value=0.1)

    per_type, pooled, predictions = first
    assert per_type["held_out_type"].tolist() == ["A", "B", "C"]
    assert len(predictions) == 24
    assert predictions["sample"].is_unique
    assert set(predictions["call"]) <= {"tumor", "normal"}
    pd.testing.assert_frame_equal(per_type, second[0])
    pd.testing.assert_frame_equal(pooled, second[1])
    pd.testing.assert_frame_equal(predictions, second[2])


def test_run_loco_rejects_patient_overlap_across_projects():
    X, samples, y, projects, cases = synthetic_loco_inputs()
    cases[8] = cases[0]
    with pytest.raises(ValueError, match="patient overlap"):
        LOCO.run_loco(X, samples, y, projects, cases, n_features=4)


def test_verify_existing_detects_metric_drift(tmp_path):
    per_type, pooled, _ = LOCO.run_loco(*synthetic_loco_inputs(), n_features=4)
    per_type.to_csv(tmp_path / "loco_per_cancer_metrics.csv", index=False)
    pooled.to_csv(tmp_path / "loco_pooled_summary.csv", index=False)
    LOCO.verify_existing(per_type, pooled, tmp_path, tolerance=1e-12)

    changed = per_type.copy()
    changed.loc[0, "auc"] -= 0.1
    with pytest.raises(ValueError, match="differ"):
        LOCO.verify_existing(changed, pooled, tmp_path, tolerance=1e-12)


def test_verify_existing_rejects_missing_schema_and_nan_drift(tmp_path):
    per_type, pooled, _ = LOCO.run_loco(*synthetic_loco_inputs(), n_features=4)
    pd.DataFrame({"held_out_type": ["A"], "bogus": ["x"]}).to_csv(
        tmp_path / "loco_per_cancer_metrics.csv", index=False
    )
    pd.DataFrame({"bogus": ["x"]}).to_csv(
        tmp_path / "loco_pooled_summary.csv", index=False
    )
    with pytest.raises(ValueError, match="schema differs"):
        LOCO.verify_existing(per_type, pooled, tmp_path, tolerance=0)

    per_type.to_csv(tmp_path / "loco_per_cancer_metrics.csv", index=False)
    pooled.to_csv(tmp_path / "loco_pooled_summary.csv", index=False)
    changed = per_type.copy()
    changed.loc[0, "auc"] = np.nan
    with pytest.raises(ValueError, match="NaN or infinite"):
        LOCO.verify_existing(changed, pooled, tmp_path, tolerance=0)


def test_load_inputs_requires_canonical_float64_even_with_escape(tmp_path):
    X, samples, y, projects, cases = synthetic_loco_inputs()
    np.save(tmp_path / "X_full.npy", X.astype(np.float32))
    np.save(tmp_path / "X_genes.npy", np.array([f"g{i}" for i in range(X.shape[1])]))
    np.save(tmp_path / "X_samples.npy", samples)
    pd.DataFrame({
        "file_id": samples,
        "project": projects,
        "label": np.where(y == 1, "tumor", "normal"),
        "case_id": cases,
    }).to_csv(tmp_path / "metadata.csv", index=False)
    with pytest.raises(ValueError, match="exact float64"):
        LOCO.load_inputs(
            tmp_path / "X_full.npy",
            tmp_path / "metadata.csv",
            allow_unverified_features=True,
        )


def test_run_loco_rejects_non_integer_feature_count():
    with pytest.raises(ValueError, match="must be an integer"):
        LOCO.run_loco(*synthetic_loco_inputs(), n_features=4.5)
