"""Unit tests for the tcga_rnaseq core primitives."""
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

from tcga_rnaseq import score as S
from tcga_rnaseq import metrics as M
from tcga_rnaseq import load_lr_model
from tcga_rnaseq.align import align_to_genes, align_to_genes_with_report
from calibrate_threshold import (
    choose_youden_threshold,
    load_scores_and_labels,
    metrics_at_threshold,
    rank_auc,
    validate_threshold,
)
from score_tumor_normal import load_lr_weights, score_dataframe_lr_weights
from run_tumor_normal_workflow import build_report
from explain_scores import EXPLANATION_COLUMNS, explain_dataframe
from cohort_adapt_score import load_label_vector


def test_sigmoid_stable_no_overflow(recwarn):
    assert S.sigmoid(0) == pytest.approx(0.5)
    # large magnitudes must not overflow and must saturate correctly
    out = S.sigmoid(np.array([-1000.0, 1000.0]))
    assert np.allclose(out, [0.0, 1.0], atol=1e-12)
    assert not any(issubclass(w.category, RuntimeWarning) for w in recwarn.list)


def test_softmax():
    p = S.softmax(np.array([[1.0, 1.0, 1.0]]))
    assert np.allclose(p, 1 / 3)
    assert S.softmax(np.array([[0.0, 100.0]]))[0, 1] == pytest.approx(1.0)
    # invariance to constant shift
    a = np.array([[0.2, 0.5, 0.3]])
    assert np.allclose(S.softmax(a), S.softmax(a + 7.0))


def test_align_imputes_missing_at_mean():
    genes = np.array(["g1", "g2", "g3"])
    mean = np.array([10.0, 20.0, 30.0])
    X = pd.DataFrame({"g1": [1.0, 2.0], "g3": [3.0, 4.0]})  # g2 missing
    v, n_matched, missing = align_to_genes(X, genes, impute_mean=mean)
    assert n_matched == 2 and missing == ["g2"]
    assert v[:, 1].tolist() == [20.0, 20.0]  # imputed at mean -> standardized 0
    assert v[:, 0].tolist() == [1.0, 2.0]


def test_align_matches_across_version_suffix():
    """Regression: a versioned model must align to an unversioned input CSV
    (and vice-versa) instead of silently NaN-imputing every column."""
    mean = np.array([99.0, 99.0])
    # model expects versioned IDs; input provides unversioned
    v, n_matched, missing = align_to_genes(
        pd.DataFrame({"ENSG0001": [1.0], "ENSG0002": [2.0]}),
        np.array(["ENSG0001.7", "ENSG0002.3"]), impute_mean=mean)
    assert n_matched == 2 and missing == []
    assert v.tolist() == [[1.0, 2.0]]
    # reverse: unversioned model, versioned input
    v2, n2, _ = align_to_genes(
        pd.DataFrame({"ENSG0001.7": [5.0]}), np.array(["ENSG0001"]), impute_mean=np.array([0.0]))
    assert n2 == 1 and v2.tolist() == [[5.0]]


def test_align_coerces_nonnumeric():
    v, _, _ = align_to_genes(pd.DataFrame({"g1": ["1.5", "oops"]}),
                             np.array(["g1"]), impute_mean=np.array([7.0]))
    assert v[:, 0].tolist() == [1.5, 7.0]  # non-numeric -> imputed at mean


def test_align_imputes_nonfinite_and_validates_mean_length():
    X = pd.DataFrame({"g1": [1.0, np.inf, -np.inf, np.nan]})
    v, n_matched, missing = align_to_genes(X, np.array(["g1"]), impute_mean=np.array([7.0]))
    assert n_matched == 1 and missing == []
    assert v[:, 0].tolist() == [1.0, 7.0, 7.0, 7.0]
    with pytest.raises(ValueError, match="impute_mean length"):
        align_to_genes(X, np.array(["g1", "g2"]), impute_mean=np.array([7.0]))


def test_align_report_counts_invalid_matched_values():
    X = pd.DataFrame(
        {"g1": ["bad", "also_bad"], "g2": [1.0, np.nan], "not_model": [5.0, 6.0]},
        index=["s1", "s2"],
    )
    values, report = align_to_genes_with_report(
        X,
        np.array(["g1", "g2", "g3"]),
        impute_mean=np.array([10.0, 20.0, 30.0]),
    )
    assert values.tolist() == [[10.0, 1.0, 30.0], [10.0, 20.0, 30.0]]
    assert report["n_matched_genes"] == 2
    assert report["missing_genes"] == ["g3"]
    assert report["matched_cells"] == 4
    assert report["invalid_matched_cells"] == 3
    assert report["n_genes_with_invalid_values"] == 2
    assert report["n_genes_with_all_invalid_values"] == 1
    assert report["first_genes_with_all_invalid_values"] == ["g1"]
    assert report["n_samples_with_invalid_values"] == 2
    assert report["max_invalid_matched_cell_fraction_per_sample"] == pytest.approx(1.0)


def test_align_rejects_duplicate_and_version_colliding_columns():
    with pytest.raises(ValueError, match="Duplicate gene columns"):
        align_to_genes(
            pd.DataFrame([[1.0, 2.0]], columns=["g1", "g1"]),
            np.array(["g1"]),
            impute_mean=np.array([0.0]),
        )
    with pytest.raises(ValueError, match="Ambiguous gene columns"):
        align_to_genes(
            pd.DataFrame({"ENSG1.1": [1.0], "ENSG1.2": [2.0]}),
            np.array(["ENSG1"]),
            impute_mean=np.array([0.0]),
        )


def test_standardize_modes_differ():
    model = {"mean": np.array([0.0, 0.0]), "scale": np.array([1.0, 1.0])}
    v = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
    z_none = S.standardize(v, model, "none")
    z_cohort = S.standardize(v, model, "cohort_zscore")
    assert np.allclose(z_none, v)                      # identity when mean=0,scale=1
    assert np.allclose(z_cohort.mean(axis=0), 0, atol=1e-9)  # cohort-centered
    assert np.allclose(z_cohort.std(axis=0), 1, atol=1e-9)   # cohort unit-variance
    with pytest.raises(ValueError):
        S.standardize(v, model, "bogus")


def test_roc_auc_matches_sklearn():
    sk = pytest.importorskip("sklearn.metrics")
    rng = np.random.default_rng(0)
    for _ in range(5):
        y = rng.integers(0, 2, 200)
        if len(set(y)) < 2:
            continue
        s = rng.random(200) + 0.3 * y
        assert M.roc_auc(y, s) == pytest.approx(sk.roc_auc_score(y, s), abs=1e-9)


def test_roc_auc_handles_ties():
    y = np.array([0, 0, 1, 1])
    s = np.array([0.5, 0.5, 0.5, 0.5])  # all tied -> AUC 0.5
    assert M.roc_auc(y, s) == pytest.approx(0.5)


def test_calibration_wrappers_preserve_release_columns():
    y = np.array([0, 0, 1, 1])
    s = np.array([0.1, 0.2, 0.8, 0.9])
    row = metrics_at_threshold(y, s, 0.5, "default")
    assert list(row) == [
        "threshold_name", "threshold", "accuracy", "f1", "precision",
        "recall", "specificity", "tn", "fp", "fn", "tp",
    ]
    assert row["accuracy"] == pytest.approx(1.0)
    assert rank_auc(y, s) == pytest.approx(M.roc_auc(y, s))
    best = choose_youden_threshold(y, s)
    assert best["threshold_name"] == "youden_j"
    assert best["youden_j"] == pytest.approx(1.0)


def test_youden_tie_break_matches_core_lower_threshold():
    y = np.array([0, 1, 0, 1])
    s = np.array([0.2, 0.4, 0.6, 0.8])
    best = choose_youden_threshold(y, s)
    core = M.youden_threshold(y, s)
    assert best["threshold"] == pytest.approx(0.4)
    assert best["threshold"] == pytest.approx(core["threshold"])


def test_calibration_rejects_bad_inputs(tmp_path):
    scores = tmp_path / "scores.csv"
    labels = tmp_path / "labels.csv"
    scores.write_text(
        "sample,tumor_probability\ns1,0.1\ns1,0.2\n",
        encoding="utf-8",
    )
    labels.write_text("sample,label\ns1,normal\n", encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate sample"):
        load_scores_and_labels(scores, labels, "sample", "label")

    scores.write_text("sample,tumor_probability\ns1,1.2\ns2,0.2\n", encoding="utf-8")
    labels.write_text("sample,label\ns1,tumor\ns2,normal\n", encoding="utf-8")
    with pytest.raises(ValueError, match="between 0 and 1"):
        load_scores_and_labels(scores, labels, "sample", "label")

    with pytest.raises(ValueError, match="between 0 and 1"):
        validate_threshold(2.0)


def test_calibration_rejects_accidental_label_subset(tmp_path):
    scores = tmp_path / "scores.csv"
    labels = tmp_path / "labels.csv"
    scores.write_text(
        "sample,tumor_probability\ns1,0.1\ns2,0.9\n",
        encoding="utf-8",
    )
    labels.write_text("sample,label\ns1,normal\n", encoding="utf-8")
    with pytest.raises(ValueError, match="below --min-match-fraction"):
        load_scores_and_labels(scores, labels, "sample", "label")


def test_cohort_adapt_labels_preserve_numeric_strings_and_missing(tmp_path):
    labels = tmp_path / "labels.csv"
    labels.write_text(
        "sample,label\ns1,1\ns2,0\nextra,tumor\n",
        encoding="utf-8",
    )
    y, matched, stats = load_label_vector(labels, pd.Index(["s1", "s2", "s3"]))
    assert y.tolist()[:2] == [1.0, 0.0]
    assert matched.tolist() == [True, True, False]
    assert stats == {
        "n_labels": 3,
        "n_labeled": 2,
        "n_unmatched_samples": 1,
        "n_extra_labels": 1,
    }


def test_accuracy_balanced_and_prf():
    y = np.array(["a", "a", "b", "b", "b"])
    p = np.array(["a", "b", "b", "b", "a"])
    assert M.accuracy(y, p) == pytest.approx(3 / 5)
    # class a: recall 1/2; class b: recall 2/3 -> balanced = (0.5+0.6667)/2
    assert M.balanced_accuracy(y, p) == pytest.approx((0.5 + 2 / 3) / 2)
    prf = {r["label"]: r for r in M.per_class_prf(y, p)}
    assert prf["a"]["support"] == 2 and prf["b"]["support"] == 3


def test_model_loading_shapes(binary_model, cancer_type_model):
    assert binary_model["kind"] == "binary"
    assert binary_model["coef"].ndim == 1
    assert binary_model["genes"].shape == (2000,)
    assert cancer_type_model["kind"] == "multiclass"
    assert cancer_type_model["coef"].shape == (17, 1000)
    assert len(cancer_type_model["classes"]) == 17


def test_score_binary_dataframe_contract():
    model = {
        "genes": np.array(["g1", "g2"]),
        "mean": np.array([0.0, 0.0]),
        "scale": np.array([1.0, 1.0]),
        "coef": np.array([1.0, -1.0]),
        "intercept": 0.0,
        "classes": np.array([0, 1]),
        "kind": "binary",
    }
    X = pd.DataFrame({"g1": [2.0, 0.0], "g2": [0.0, 2.0]}, index=["s1", "s2"])
    scored, n_matched, missing = S.score_binary_dataframe(model, X, threshold=0.5)
    assert list(scored.columns) == ["sample", "tumor_probability", "call"]
    assert scored["sample"].tolist() == ["s1", "s2"]
    assert scored["call"].tolist() == ["tumor", "normal"]
    assert n_matched == 2 and missing == []
    scored2, _, _, report = S.score_binary_dataframe(
        model, pd.DataFrame({"g1": ["bad"], "g2": [2.0]}, index=["s3"]),
        return_alignment_report=True,
    )
    assert list(scored2.columns) == ["sample", "tumor_probability", "call"]
    assert report["invalid_matched_cells"] == 1


def test_load_lr_model_rejects_inconsistent_shapes(tmp_path):
    bad = tmp_path / "bad_model.npz"
    np.savez(
        bad,
        selected_genes=np.array(["g1", "g2"]),
        scaler_mean=np.array([0.0]),
        scaler_scale=np.array([1.0, 1.0]),
        coef=np.array([1.0, 2.0]),
        intercept=np.array(0.0),
    )
    with pytest.raises(ValueError, match="one value per selected gene"):
        load_lr_model(bad)


def test_binary_cli_rejects_multiclass_weights(root, cancer_type_model):
    with pytest.raises(ValueError, match="requires a binary LR weights file"):
        load_lr_weights(f"{root}/cancer-type-classifier/cancer_type_lr_weights.npz")
    with pytest.raises(ValueError, match="requires a binary model"):
        score_dataframe_lr_weights(pd.DataFrame({"g1": [1.0]}), cancer_type_model)


def test_cancer_type_cli_rejects_invalid_matched_values(tmp_path, root, cancer_type_model):
    gene = str(cancer_type_model["genes"][0])
    input_path = tmp_path / "invalid_cancer_type.csv"
    output_path = tmp_path / "predictions.csv"
    input_path.write_text(f"sample,{gene}\ns1,not_numeric\n", encoding="utf-8")
    result = subprocess.run(
        [
            sys.executable,
            "cancer-type-classifier/predict_cancer_type.py",
            str(input_path),
            "--out",
            str(output_path),
        ],
        cwd=root,
        text=True,
        capture_output=True,
    )
    assert result.returncode != 0
    assert "invalid matched values" in result.stderr
    assert "Refusing to write predictions" in result.stderr
    assert not output_path.exists()

    result = subprocess.run(
        [
            sys.executable,
            "cancer-type-classifier/predict_cancer_type.py",
            str(input_path),
            "--out",
            str(output_path),
            "--allow-invalid-values",
        ],
        cwd=root,
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0
    assert output_path.exists()


def test_workflow_qc_fail_report_does_not_claim_zero_scored_samples():
    qc = {
        "status": "FAIL",
        "messages": [],
        "shape": {"samples": 3, "input_genes": 2, "model_genes": 2},
        "gene_match": {"matched_model_genes": 0, "match_rate": 0.0},
        "value_summary": {"p50": None, "p99": None, "max": None},
        "distribution_summary": {
            "abs_z_gt_6_fraction": None,
            "cohort_gene_mean_abs_z": {"p99": None},
        },
        "score_summary": {
            "tumor_probability": {"p50": None, "p90": None, "max": None},
        },
    }
    report = build_report("bad.csv", qc, None, {"qc_json": "qc.json"})
    assert "Scoring was not run" in report
    assert "Input samples inspected: 3" in report
    assert "Samples: 0" not in report


def test_explain_empty_input_preserves_output_columns():
    weights = {
        "selected_genes": np.array(["g1"]),
        "scaler_mean": np.array([0.0]),
        "scaler_scale": np.array([1.0]),
        "coef": np.array([1.0]),
        "intercept": 0.0,
    }
    df = pd.DataFrame({"g1": []})
    explanations, n_matched, missing = explain_dataframe(df, weights, top_n=1)
    assert list(explanations.columns) == EXPLANATION_COLUMNS
    assert explanations.empty
    assert n_matched == 1 and missing == []


def test_explain_dataframe_can_return_alignment_report():
    weights = {
        "selected_genes": np.array(["g1"]),
        "scaler_mean": np.array([0.0]),
        "scaler_scale": np.array([1.0]),
        "coef": np.array([1.0]),
        "intercept": 0.0,
    }
    df = pd.DataFrame({"g1": ["bad"]}, index=["s1"])
    explanations, n_matched, missing, report = explain_dataframe(
        df,
        weights,
        top_n=1,
        return_alignment_report=True,
    )
    assert list(explanations.columns) == EXPLANATION_COLUMNS
    assert n_matched == 1 and missing == []
    assert report["invalid_matched_cells"] == 1
