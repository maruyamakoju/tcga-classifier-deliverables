"""Regression tests for cohort_adapt_score.py's metrics/output consistency."""
import json

import numpy as np
import pandas as pd
import pytest

import cohort_adapt_score


def write_tiny_binary_model(path, coef=1.0, intercept=0.0):
    np.savez(
        path,
        selected_genes=np.array(["ENSG_TEST"]),
        scaler_mean=np.array([0.0]),
        scaler_scale=np.array([1.0]),
        coef=np.array([coef]),
        intercept=np.array(intercept),
        classes=np.array(["normal", "tumor"]),
    )


def test_metrics_use_raw_probability_not_rounded_call_boundary(tmp_path, capsys):
    """A sample whose raw probability rounds across the threshold must still
    be scored consistently between the 'call' column and the metrics dict.

    Regression test: score_binary_dataframe derives `call` from the raw
    (unrounded) probability, but rounds `tumor_probability` to 6 decimals for
    display. cohort_adapt_score.py used to recompute metrics by
    re-thresholding that *rounded* column, which could disagree with `call`
    for a probability that rounds across the threshold (e.g. 0.4999995 -> a
    true "normal" call, but round(.., 6) == 0.5 >= 0.5 threshold => "tumor").
    """
    weights_path = tmp_path / "weights.npz"
    write_tiny_binary_model(weights_path)

    # x such that sigmoid(x) == 0.4999995 (< 0.5, rounds to exactly 0.5 at 6dp).
    x = float(np.log((0.5 - 5e-7) / (0.5 + 5e-7)))
    assert round(1 / (1 + np.exp(-x)), 6) == 0.5
    assert (1 / (1 + np.exp(-x))) < 0.5

    # A second, unambiguous true-tumor sample so both classes are present in
    # the labeled subset (classification_metrics only runs with >1 class).
    input_csv = tmp_path / "input.csv"
    pd.DataFrame(
        {"ENSG_TEST": [x, 10.0]}, index=["sample_a", "sample_b"]
    ).to_csv(input_csv, index_label="sample")

    labels_csv = tmp_path / "labels.csv"
    pd.DataFrame(
        {"sample": ["sample_a", "sample_b"], "label": ["normal", "tumor"]}
    ).to_csv(labels_csv, index=False)

    out_csv = tmp_path / "out.csv"
    code = cohort_adapt_score.main([
        str(input_csv),
        "--weights", str(weights_path),
        "--adapt", "none",
        "--threshold", "0.5",
        "--labels", str(labels_csv),
        "--out", str(out_csv),
    ])
    assert code == 0

    out_df = pd.read_csv(out_csv).set_index("sample")
    assert out_df.loc["sample_a", "call"] == "normal"
    assert out_df.loc["sample_a", "tumor_probability"] == pytest.approx(0.5)
    assert out_df.loc["sample_b", "call"] == "tumor"

    report = json.loads(capsys.readouterr().out)
    assert report["normal_calls"] == 1
    assert report["tumor_calls"] == 1
    # sample_a's rounded probability lands exactly on the threshold, but its
    # raw probability is < 0.5 and its call is "normal" -- metrics must not
    # score it as a false positive just because the *displayed* probability
    # rounds up to the threshold.
    assert report["metrics"]["specificity"] == 1.0
    assert report["metrics"]["accuracy"] == 1.0
