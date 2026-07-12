"""Fail-closed and numerical invariants for public scoring APIs."""
import io

import numpy as np
import pandas as pd
import pytest

from explain_scores import explain_dataframe
from score_tumor_normal import load_lr_weights
from tcga_rnaseq import score_binary_dataframe
from tcga_rnaseq.score import predict_proba


def binary_model(coef=(1.0, -1.0), classes=(0, 1)):
    return {
        "genes": np.array(["g1", "g2"]),
        "mean": np.array([0.0, 0.0]),
        "scale": np.array([1.0, 1.0]),
        "coef": np.asarray(coef, dtype=float),
        "intercept": 0.0,
        "classes": np.asarray(classes),
        "kind": "binary",
    }


def test_core_rejects_zero_gene_coverage_by_default():
    frame = pd.DataFrame({"unrelated": [1.0]}, index=["s1"])

    with pytest.raises(ValueError, match="low model-gene coverage"):
        predict_proba(binary_model(), frame)

    probabilities, report = predict_proba(
        binary_model(),
        frame,
        allow_low_gene_coverage=True,
        return_alignment_report=True,
    )
    assert probabilities.tolist() == pytest.approx([0.5])
    assert report["n_matched_genes"] == 0


@pytest.mark.parametrize("threshold", [np.nan, np.inf, -0.01, 1.01, True])
def test_core_rejects_invalid_threshold(threshold):
    frame = pd.DataFrame({"g1": [1.0], "g2": [0.0]}, index=["s1"])
    with pytest.raises(ValueError, match="threshold"):
        score_binary_dataframe(binary_model(), frame, threshold=threshold)


def test_nonfinite_model_and_numeric_overflow_fail_closed():
    frame = pd.DataFrame({"g1": [1.0], "g2": [1.0]}, index=["s1"])
    with pytest.raises(ValueError, match="coefficients and intercept must be finite"):
        score_binary_dataframe(binary_model(coef=(np.nan, 0.0)), frame)

    huge_model = binary_model(coef=(1e308, 1e308))
    huge_frame = pd.DataFrame({"g1": [1e308], "g2": [1e308]}, index=["s1"])
    with pytest.raises(ValueError, match="overflow|non-finite"):
        score_binary_dataframe(huge_model, huge_frame)


def test_tumor_score_requires_normal_then_tumor_class_order():
    frame = pd.DataFrame({"g1": [10.0], "g2": [0.0]}, index=["s1"])

    with pytest.raises(ValueError, match="binary class order"):
        score_binary_dataframe(
            binary_model(classes=("tumor", "normal")), frame
        )


def test_public_weight_adapter_preserves_binary_class_meaning(tmp_path):
    path = tmp_path / "reversed.npz"
    np.savez(
        path,
        selected_genes=np.array(["g1", "g2"]),
        scaler_mean=np.zeros(2),
        scaler_scale=np.ones(2),
        coef=np.array([1.0, -1.0]),
        intercept=np.array(0.0),
        classes=np.array(["tumor", "normal"]),
    )

    with pytest.raises(ValueError, match="binary class order"):
        load_lr_weights(path)


def test_serialized_probability_and_call_remain_consistent():
    probability = 0.4999996
    logit = float(np.log(probability / (1.0 - probability)))
    frame = pd.DataFrame({"g1": [logit], "g2": [0.0]}, index=["edge"])
    scored, _, _ = score_binary_dataframe(binary_model(coef=(1.0, 0.0)), frame)

    buffer = io.StringIO()
    scored.to_csv(buffer, index=False)
    observed = pd.read_csv(io.StringIO(buffer.getvalue()))

    assert observed.loc[0, "tumor_probability"] < 0.5
    assert observed.loc[0, "call"] == "normal"
    assert observed.loc[0, "call"] == (
        "tumor" if observed.loc[0, "tumor_probability"] >= 0.5 else "normal"
    )


@pytest.mark.parametrize("index", [[""], [" s1 "], ["s1", "s1"]])
def test_core_rejects_invalid_sample_ids(index):
    frame = pd.DataFrame(
        {"g1": [1.0] * len(index), "g2": [0.0] * len(index)}, index=index
    )
    with pytest.raises(ValueError, match="sample identifiers|duplicate sample"):
        score_binary_dataframe(binary_model(), frame)


def test_explanations_use_true_contribution_sign_and_no_overlap():
    weights = {
        "selected_genes": np.array(["g1", "g2"]),
        "scaler_mean": np.zeros(2),
        "scaler_scale": np.ones(2),
        "coef": np.array([1.0, 2.0]),
        "intercept": 0.0,
    }
    frame = pd.DataFrame(
        {"g1": [1.0, -1.0], "g2": [1.0, -1.0]}, index=["positive", "negative"]
    )

    explanations, _, _ = explain_dataframe(frame, weights, top_n=2)

    tumor = explanations[explanations["direction"] == "tumor"]
    normal = explanations[explanations["direction"] == "normal"]
    assert (tumor["contribution_logit"] > 0).all()
    assert (normal["contribution_logit"] < 0).all()
    assert not explanations.duplicated(["sample", "gene_id"]).any()
    assert set(tumor["sample"]) == {"positive"}
    assert set(normal["sample"]) == {"negative"}
