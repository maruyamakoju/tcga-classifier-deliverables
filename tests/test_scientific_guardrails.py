"""Conservative opt-in contracts for adaptation and calibration."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import calibrate_threshold
import cohort_adapt_score


ROOT = Path(__file__).resolve().parents[1]


def write_tiny_model(path):
    np.savez(
        path,
        selected_genes=np.array(["g1"]),
        scaler_mean=np.array([0.0]),
        scaler_scale=np.array([1.0]),
        coef=np.array([1.0]),
        intercept=np.array(0.0),
        classes=np.array(["normal", "tumor"]),
    )


def write_expression(path, n_samples):
    pd.DataFrame(
        {"g1": np.linspace(-2.0, 2.0, n_samples)},
        index=[f"s{i:02d}" for i in range(n_samples)],
    ).to_csv(path, index_label="sample")


def test_adaptation_is_opt_in_and_default_is_inductive(tmp_path, capsys):
    model = tmp_path / "model.npz"
    expression = tmp_path / "expression.csv"
    output = tmp_path / "scores.csv"
    write_tiny_model(model)
    write_expression(expression, 2)

    code = cohort_adapt_score.main(
        [str(expression), "--weights", str(model), "--out", str(output)]
    )

    assert code == 0
    captured = capsys.readouterr()
    report = json.loads(captured.out)
    assert report["adapt"] == "none"
    assert not any("transductive" in warning for warning in report["warnings"])


def test_small_adapted_cohort_hard_fails_before_output(tmp_path, capsys):
    model = tmp_path / "model.npz"
    expression = tmp_path / "expression.csv"
    output = tmp_path / "scores.csv"
    write_tiny_model(model)
    write_expression(expression, 2)

    with pytest.raises(SystemExit) as exc_info:
        cohort_adapt_score.main(
            [
                str(expression),
                "--weights",
                str(model),
                "--adapt",
                "cohort_zscore",
                "--min-samples",
                "20",
                "--out",
                str(output),
            ]
        )

    assert exc_info.value.code == 2
    assert "requires at least --min-samples=20" in capsys.readouterr().err
    assert not output.exists()


def test_adapted_report_discloses_composition_and_transductive_behavior(
    tmp_path, capsys
):
    model = tmp_path / "model.npz"
    expression = tmp_path / "expression.csv"
    output = tmp_path / "scores.csv"
    write_tiny_model(model)
    write_expression(expression, 20)

    code = cohort_adapt_score.main(
        [
            str(expression),
            "--weights",
            str(model),
            "--adapt",
            "cohort_zscore",
            "--out",
            str(output),
        ]
    )

    assert code == 0
    captured = capsys.readouterr()
    report = json.loads(captured.out)
    assert any("composition" in warning for warning in report["warnings"])
    assert any("transductive" in warning for warning in report["warnings"])
    assert "transductive" in captured.err
    assert report["score_interpretation"] == "model score, not calibrated clinical risk"


def test_malformed_adaptation_labels_fail_before_output(tmp_path, capsys):
    model = tmp_path / "model.npz"
    expression = tmp_path / "expression.csv"
    labels = tmp_path / "labels.csv"
    output = tmp_path / "scores.csv"
    write_tiny_model(model)
    write_expression(expression, 2)
    pd.DataFrame(
        {"sample": ["s00", "s01"], "label": ["normal", "not-a-label"]}
    ).to_csv(labels, index=False)

    with pytest.raises(SystemExit) as exc_info:
        cohort_adapt_score.main(
            [
                str(expression),
                "--weights",
                str(model),
                "--labels",
                str(labels),
                "--out",
                str(output),
            ]
        )

    assert exc_info.value.code == 2
    assert "Unrecognized label" in capsys.readouterr().err
    assert not output.exists()


def test_calibration_summary_discloses_resubstitution_and_small_classes(
    tmp_path, capsys
):
    scores = tmp_path / "scores.csv"
    labels = tmp_path / "labels.csv"
    metrics = tmp_path / "thresholds.csv"
    summary_path = tmp_path / "calibration.json"
    pd.DataFrame(
        {"sample": ["n1", "n2", "t1", "t2"],
         "tumor_probability": [0.1, 0.2, 0.8, 0.9]}
    ).to_csv(scores, index=False)
    pd.DataFrame(
        {"sample": ["n1", "n2", "t1", "t2"],
         "label": ["normal", "normal", "tumor", "tumor"]}
    ).to_csv(labels, index=False)

    code = calibrate_threshold.main(
        [
            str(scores),
            str(labels),
            "-o",
            str(metrics),
            "--json-output",
            str(summary_path),
        ]
    )

    assert code == 0
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    assert summary["evaluation_type"] == "apparent_resubstitution"
    assert summary["warnings"]
    assert "fewer than 10" in summary["warnings"][0]
    assert "apparent/resubstitution" in capsys.readouterr().err


def test_full_research_reports_retain_scientific_scope_caveats():
    adaptation = (
        ROOT / "cross-platform-adaptation" / "CROSS_PLATFORM_ADAPTATION.md"
    ).read_text(encoding="utf-8")
    loco = (ROOT / "cross-cancer-holdout" / "LOCO_REPORT.md").read_text(
        encoding="utf-8"
    )
    cancer_type = (
        ROOT / "cancer-type-classifier" / "CANCER_TYPE_CLASSIFIER.md"
    ).read_text(encoding="utf-8")
    adaptation = " ".join(adaptation.split())
    loco = " ".join(loco.split())
    cancer_type = " ".join(cancer_type.split())

    assert "transductive, composition-dependent" in adaptation
    assert "does not" in adaptation and "restore probability calibration" in adaptation
    assert "does not rule out shared GDC/procurement confounding" in loco
    assert "cross-fitted aggregates across 17 different" in loco
    assert "per-type AUCs and their macro summary as the primary" in loco
    assert "not an external or prospective generalization estimate" in cancer_type
    assert "not causal interpretation" in cancer_type

    adaptation_html = (
        ROOT / "cross-platform-adaptation" / "cross_platform_adaptation.html"
    ).read_text(encoding="ascii")
    loco_html = (ROOT / "cross-cancer-holdout" / "loco_report.html").read_text(
        encoding="utf-8"
    )
    cancer_html = (
        ROOT / "cancer-type-classifier" / "cancer_type_classifier.html"
    ).read_text(encoding="ascii")
    assert "transductive, composition-dependent" in adaptation_html
    assert "does not exclude GDC-wide" in loco_html
    assert "not global metrics from one frozen model" in loco_html
    assert "not causal proof" in cancer_html
