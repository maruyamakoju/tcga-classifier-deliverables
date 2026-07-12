"""Unit tests for the tcga_rnaseq core primitives."""
import json
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

from tcga_rnaseq import score as S
from tcga_rnaseq import metrics as M
from tcga_rnaseq import load_lr_model, read_matrix, write_json
from tcga_rnaseq.io import load_pipeline, _XgbStub
from tcga_rnaseq.align import (
    align_to_genes,
    align_to_genes_with_report,
    build_gene_column_lookups,
    format_alignment_issues,
    format_gene_match_issues,
    print_invalid_alignment_summary,
    validate_alignment_report,
    validate_gene_match_report,
)
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


def test_gene_match_report_flags_low_coverage():
    assert validate_gene_match_report(
        {"n_model_genes": 4, "n_matched_genes": 1},
        min_match_rate=0.5,
    )
    assert validate_gene_match_report(
        {"n_model_genes": 4, "n_matched_genes": 2},
        min_match_rate=0.5,
    ) == []


def test_align_coerces_nonnumeric():
    with pytest.raises(ValueError, match="invalid matched values"):
        align_to_genes(
            pd.DataFrame({"g1": ["1.5", "oops"]}),
            np.array(["g1"]),
            impute_mean=np.array([7.0]),
        )
    v, _, _ = align_to_genes(
        pd.DataFrame({"g1": ["1.5", "oops"]}),
        np.array(["g1"]),
        impute_mean=np.array([7.0]),
        allow_invalid_values=True,
    )
    assert v[:, 0].tolist() == [1.5, 7.0]  # non-numeric -> imputed at mean


def test_align_imputes_nonfinite_and_validates_mean_length():
    X = pd.DataFrame({"g1": [1.0, np.inf, -np.inf, np.nan]})
    with pytest.raises(ValueError, match="invalid matched values"):
        align_to_genes(X, np.array(["g1"]), impute_mean=np.array([7.0]))
    v, n_matched, missing = align_to_genes(
        X,
        np.array(["g1"]),
        impute_mean=np.array([7.0]),
        allow_invalid_values=True,
    )
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
    # s2 is invalid across every matched gene (g1="also_bad", g2=NaN); s1 only
    # has one of its two matched genes invalid, so it isn't "all invalid".
    assert report["n_samples_with_all_invalid_values"] == 1
    assert report["first_samples_with_all_invalid_values"] == ["s2"]


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


def test_build_gene_column_lookups_returns_exact_and_stripped():
    exact, stripped = build_gene_column_lookups(["g1", "g2.3", "g3"])
    assert exact == {"g1": "g1", "g2.3": "g2.3", "g3": "g3"}
    assert stripped == {"g1": "g1", "g2": "g2.3", "g3": "g3"}


def test_build_gene_column_lookups_rejects_exact_duplicates():
    with pytest.raises(ValueError, match="Duplicate gene columns"):
        build_gene_column_lookups(["g1", "g2", "g1"])


def test_build_gene_column_lookups_rejects_version_suffix_collisions():
    with pytest.raises(ValueError, match="Ambiguous gene columns"):
        build_gene_column_lookups(["g1.1", "g1.2"])


def test_validate_alignment_report_flags_all_invalid_genes():
    report = {
        "invalid_matched_cells": 5,
        "n_genes_with_all_invalid_values": 2,
        "n_samples_with_all_invalid_values": 0,
        "invalid_matched_fraction": 0.0,
        "max_invalid_matched_cell_fraction_per_sample": 0.0,
        "first_genes_with_all_invalid_values": ["g1", "g2"],
        "first_samples_with_all_invalid_values": [],
    }
    issues = validate_alignment_report(report, max_invalid_cell_fraction=0.0)
    assert issues == ["2 matched model genes have no finite values. Examples: g1, g2."]


def test_validate_alignment_report_flags_all_invalid_samples():
    report = {
        "invalid_matched_cells": 5,
        "n_genes_with_all_invalid_values": 0,
        "n_samples_with_all_invalid_values": 1,
        "invalid_matched_fraction": 0.0,
        "max_invalid_matched_cell_fraction_per_sample": 0.0,
        "first_genes_with_all_invalid_values": [],
        "first_samples_with_all_invalid_values": ["s2"],
    }
    issues = validate_alignment_report(report, max_invalid_cell_fraction=0.0)
    assert issues == ["1 samples have no finite matched model-gene values. Examples: s2."]


def test_validate_alignment_report_flags_overall_fraction_exceeded():
    report = {
        "invalid_matched_cells": 3,
        "n_genes_with_all_invalid_values": 0,
        "n_samples_with_all_invalid_values": 0,
        "invalid_matched_fraction": 0.02,
        "max_invalid_matched_cell_fraction_per_sample": 0.0,
        "first_genes_with_all_invalid_values": [],
        "first_samples_with_all_invalid_values": [],
    }
    issues = validate_alignment_report(report, max_invalid_cell_fraction=0.01)
    assert issues == [
        "Invalid matched-value fraction 2.000% exceeds --max-invalid-cell-fraction 1.000%."
    ]


def test_validate_alignment_report_flags_worst_sample_fraction_exceeded():
    report = {
        "invalid_matched_cells": 3,
        "n_genes_with_all_invalid_values": 0,
        "n_samples_with_all_invalid_values": 0,
        "invalid_matched_fraction": 0.0,
        "max_invalid_matched_cell_fraction_per_sample": 0.5,
        "first_genes_with_all_invalid_values": [],
        "first_samples_with_all_invalid_values": [],
    }
    issues = validate_alignment_report(report, max_invalid_cell_fraction=0.1)
    assert issues == [
        "Worst-sample invalid matched-value fraction 50.000% exceeds "
        "--max-invalid-cell-fraction 10.000%."
    ]


def test_format_gene_match_issues_clean_and_low_coverage():
    assert format_gene_match_issues({"n_model_genes": 4, "n_matched_genes": 4}) == ""
    message = format_gene_match_issues(
        {"n_model_genes": 4, "n_matched_genes": 1}, min_match_rate=0.5
    )
    assert message.startswith("low model-gene coverage: ")
    assert "25.0%" in message


def test_format_alignment_issues_clean_and_invalid():
    assert format_alignment_issues({"invalid_matched_cells": 0}) == ""
    bad_report = {
        "invalid_matched_cells": 5,
        "n_genes_with_all_invalid_values": 1,
        "n_samples_with_all_invalid_values": 0,
        "invalid_matched_fraction": 0.0,
        "max_invalid_matched_cell_fraction_per_sample": 0.0,
        "first_genes_with_all_invalid_values": ["g1"],
        "first_samples_with_all_invalid_values": [],
    }
    message = format_alignment_issues(bad_report)
    assert message.startswith("invalid matched values: ")
    assert "g1" in message


def test_print_invalid_alignment_summary_includes_gene_and_sample_examples(capsys):
    report = {
        "invalid_matched_cells": 4,
        "matched_cells": 10,
        "invalid_matched_fraction": 0.4,
        "n_genes_with_invalid_values": 2,
        "n_samples_with_invalid_values": 2,
        "first_genes_with_invalid_values": [
            {"gene": "g1", "invalid_cells": 2, "total_cells": 5},
        ],
        "first_samples_with_invalid_values": [
            {"sample": "s1", "invalid_cells": 2, "matched_genes": 2},
        ],
    }
    print_invalid_alignment_summary(report, sys.stdout)
    out = capsys.readouterr().out
    assert "g1:2/5" in out
    assert "s1:2/2" in out


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

    # cohort_center: location-only adaptation -> (x - cohort_mean) / train_scale
    scaled_model = {"mean": np.array([0.0, 0.0]), "scale": np.array([2.0, 4.0])}
    z_center = S.standardize(v, scaled_model, "cohort_center")
    cohort_mean = v.mean(axis=0)
    assert np.allclose(z_center, (v - cohort_mean) / scaled_model["scale"])


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


def test_external_labels_override_same_named_scores_column_without_merge_collision(tmp_path):
    scores = tmp_path / "scores.csv"
    labels = tmp_path / "labels.csv"
    scores.write_text(
        "sample,tumor_probability,label\ns1,0.1,stale\ns2,0.9,stale\n",
        encoding="utf-8",
    )
    labels.write_text("sample,label\ns1,normal\ns2,tumor\n", encoding="utf-8")

    observed = load_scores_and_labels(scores, labels, "sample", "label")

    assert observed["label_binary"].tolist() == [0, 1]


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


def test_cohort_adapt_labels_reject_missing_sample_ids(tmp_path):
    labels = tmp_path / "labels.csv"
    labels.write_text(
        "sample,label\n,normal\ns2,tumor\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="sample identifiers must be non-empty"):
        load_label_vector(labels, pd.Index(["s1", "s2"]))


def test_accuracy_balanced_and_prf():
    y = np.array(["a", "a", "b", "b", "b"])
    p = np.array(["a", "b", "b", "b", "a"])
    assert M.accuracy(y, p) == pytest.approx(3 / 5)
    # class a: recall 1/2; class b: recall 2/3 -> balanced = (0.5+0.6667)/2
    assert M.balanced_accuracy(y, p) == pytest.approx((0.5 + 2 / 3) / 2)
    prf = {r["label"]: r for r in M.per_class_prf(y, p)}
    assert prf["a"]["support"] == 2 and prf["b"]["support"] == 3
    # confusion: true a -> {a:1, b:1}; true b -> {a:1, b:2}
    # class a: tp=1, fp=1 (one "a" prediction was actually b), fn=1 -> prec=rec=f1=0.5
    assert prf["a"]["precision"] == pytest.approx(0.5)
    assert prf["a"]["recall"] == pytest.approx(0.5)
    assert prf["a"]["f1"] == pytest.approx(0.5)
    # class b: tp=2, fp=1, fn=1 -> prec=rec=f1=2/3
    assert prf["b"]["precision"] == pytest.approx(2 / 3)
    assert prf["b"]["recall"] == pytest.approx(2 / 3)
    assert prf["b"]["f1"] == pytest.approx(2 / 3)


def test_confusion_matrix_explicit_labels_ordering():
    y = np.array(["a", "b", "a", "c"])
    p = np.array(["a", "a", "b", "c"])
    m, labels = M.confusion_matrix(y, p, labels=["c", "b", "a"])
    assert labels == ["c", "b", "a"]
    # rows/cols follow the explicit ["c","b","a"] order, not sorted order
    assert m.tolist() == [
        [1, 0, 0],  # true=c: predicted c once
        [0, 0, 1],  # true=b: predicted a once
        [0, 1, 1],  # true=a: predicted b once, predicted a once
    ]


def test_macro_f1_matches_hand_computed_value():
    y = np.array(["a", "a", "b", "b", "b"])
    p = np.array(["a", "b", "b", "b", "a"])
    # per-class f1: a=0.5, b=2/3 (see test_accuracy_balanced_and_prf)
    assert M.macro_f1(y, p) == pytest.approx((0.5 + 2 / 3) / 2)


def test_threshold_sweep_returns_metrics_per_threshold():
    y = np.array([0, 0, 1, 1])
    s = np.array([0.1, 0.4, 0.6, 0.9])
    thresholds = [0.2, 0.5, 0.8]
    results = M.threshold_sweep(y, s, thresholds)
    assert len(results) == len(thresholds)
    assert [r["threshold"] for r in results] == thresholds
    expected_keys = set(M.classification_metrics(y, s, thresholds[0]).keys())
    for r in results:
        assert set(r.keys()) == expected_keys


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
        allow_invalid_values=True,
        return_alignment_report=True,
    )
    assert list(scored2.columns) == ["sample", "tumor_probability", "call"]
    assert report["invalid_matched_cells"] == 1


def test_core_scoring_rejects_invalid_matched_values_by_default():
    model = {
        "genes": np.array(["g1", "g2"]),
        "mean": np.array([0.0, 0.0]),
        "scale": np.array([1.0, 1.0]),
        "coef": np.array([1.0, -1.0]),
        "intercept": 0.0,
        "classes": np.array([0, 1]),
        "kind": "binary",
    }
    X = pd.DataFrame({"g1": ["bad"], "g2": [2.0]}, index=["s1"])
    with pytest.raises(ValueError, match="invalid matched values"):
        S.predict_proba(model, X)
    with pytest.raises(ValueError, match="invalid matched values"):
        S.score_binary_dataframe(model, X)

    p, report = S.predict_proba(model, X, allow_invalid_values=True, return_alignment_report=True)
    assert p.shape == (1,)
    assert report["invalid_matched_cells"] == 1
    scored, _, _, report = S.score_binary_dataframe(
        model,
        X,
        allow_invalid_values=True,
        return_alignment_report=True,
    )
    assert list(scored.columns) == ["sample", "tumor_probability", "call"]
    assert report["invalid_matched_cells"] == 1


def test_predict_proba_from_aligned_validates_shape_and_finiteness():
    model = {
        "genes": np.array(["g1", "g2"]),
        "mean": np.array([0.0, 0.0]),
        "scale": np.array([1.0, 1.0]),
        "coef": np.array([1.0, -1.0]),
        "intercept": 0.0,
        "classes": np.array([0, 1]),
        "kind": "binary",
    }
    with pytest.raises(ValueError, match="expected 2 model genes"):
        S.predict_proba_from_aligned(model, np.array([[1.0]]))
    with pytest.raises(ValueError, match="non-finite"):
        S.predict_proba_from_aligned(model, np.array([[1.0, np.nan]]))


def test_predict_returns_binary_class_labels():
    model = {
        "genes": np.array(["g1", "g2"]),
        "mean": np.array([0.0, 0.0]),
        "scale": np.array([1.0, 1.0]),
        "coef": np.array([1.0, -1.0]),
        "intercept": 0.0,
        "classes": np.array([0, 1]),
        "kind": "binary",
    }
    # s1: g1=2,g2=0 -> logit=2  -> proba~0.88 -> class 1
    # s2: g1=0,g2=2 -> logit=-2 -> proba~0.12 -> class 0
    X = pd.DataFrame({"g1": [2.0, 0.0], "g2": [0.0, 2.0]}, index=["s1", "s2"])
    calls = S.predict(model, X, threshold=0.5)
    assert calls.tolist() == [1, 0]


def test_predict_uses_default_binary_classes_after_model_canonicalization():
    model = {
        "genes": np.array(["g1"]),
        "mean": np.array([0.0]),
        "scale": np.array([1.0]),
        "coef": np.array([1.0]),
        "intercept": 0.0,
    }
    X = pd.DataFrame({"g1": [1.0, -1.0]}, index=["s1", "s2"])
    assert S.predict(model, X).tolist() == [1, 0]


def test_predict_returns_multiclass_class_labels():
    model = {
        "genes": np.array(["g1", "g2"]),
        "mean": np.array([0.0, 0.0]),
        "scale": np.array([1.0, 1.0]),
        "coef": np.array([[1.0, 0.0], [0.0, 1.0], [-1.0, -1.0]]),
        "intercept": np.array([0.0, 0.0, 0.0]),
        "classes": np.array(["A", "B", "C"]),
        "kind": "multiclass",
    }
    # s1 favors class A, s2 favors class B, s3 favors class C
    X = pd.DataFrame(
        {"g1": [5.0, 0.0, -5.0], "g2": [0.0, 5.0, -5.0]}, index=["s1", "s2", "s3"]
    )
    calls = S.predict(model, X)
    assert calls.tolist() == ["A", "B", "C"]


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


def _binary_npz_kwargs(n_genes=3):
    return dict(
        selected_genes=np.array([f"g{i}" for i in range(n_genes)]),
        scaler_mean=np.zeros(n_genes),
        scaler_scale=np.ones(n_genes),
        coef=np.ones(n_genes),
        intercept=np.array(0.0),
    )


def _multiclass_npz_kwargs(n_genes=3, n_classes=2):
    return dict(
        selected_genes=np.array([f"g{i}" for i in range(n_genes)]),
        scaler_mean=np.zeros(n_genes),
        scaler_scale=np.ones(n_genes),
        coef=np.zeros((n_classes, n_genes)),
        intercept=np.zeros(n_classes),
        classes=np.array([f"c{i}" for i in range(n_classes)]),
    )


def test_load_lr_model_rejects_nonfinite_scaler_values(tmp_path):
    nan_mean = tmp_path / "nan_mean.npz"
    kwargs = _binary_npz_kwargs()
    kwargs["scaler_mean"] = np.array([0.0, np.nan, 1.0])
    np.savez(nan_mean, **kwargs)
    with pytest.raises(ValueError, match="must be finite"):
        load_lr_model(nan_mean)

    inf_scale = tmp_path / "inf_scale.npz"
    kwargs = _binary_npz_kwargs()
    kwargs["scaler_scale"] = np.array([1.0, np.inf, 1.0])
    np.savez(inf_scale, **kwargs)
    with pytest.raises(ValueError, match="must be finite"):
        load_lr_model(inf_scale)


def test_load_lr_model_rejects_nonpositive_scale(tmp_path):
    bad = tmp_path / "nonpositive_scale.npz"
    kwargs = _binary_npz_kwargs()
    kwargs["scaler_scale"] = np.array([1.0, 0.0, 1.0])
    np.savez(bad, **kwargs)
    with pytest.raises(ValueError, match="must be positive"):
        load_lr_model(bad)


def test_load_lr_model_rejects_binary_coef_shape_mismatch(tmp_path):
    bad = tmp_path / "binary_coef_mismatch.npz"
    kwargs = _binary_npz_kwargs(n_genes=3)
    kwargs["coef"] = np.array([1.0, 2.0])  # length 2, but 3 genes
    np.savez(bad, **kwargs)
    with pytest.raises(ValueError, match="does not match genes="):
        load_lr_model(bad)


def test_load_lr_model_rejects_multiclass_coef_gene_mismatch(tmp_path):
    bad = tmp_path / "multiclass_coef_mismatch.npz"
    kwargs = _multiclass_npz_kwargs(n_genes=3, n_classes=2)
    kwargs["coef"] = np.zeros((2, 4))  # coef.shape[1]=4 != n_genes=3
    np.savez(bad, **kwargs)
    with pytest.raises(ValueError, match="does not match genes="):
        load_lr_model(bad)


def test_load_lr_model_rejects_multiclass_intercept_length_mismatch(tmp_path):
    bad = tmp_path / "multiclass_intercept_mismatch.npz"
    kwargs = _multiclass_npz_kwargs(n_genes=3, n_classes=2)
    kwargs["intercept"] = np.zeros(3)  # coef.shape[0]=2, intercept length 3
    np.savez(bad, **kwargs)
    with pytest.raises(
        ValueError, match="Multiclass intercept length must match number of coefficient rows"
    ):
        load_lr_model(bad)


def test_load_lr_model_rejects_multiclass_class_count_mismatch(tmp_path):
    bad = tmp_path / "multiclass_class_count_mismatch.npz"
    kwargs = _multiclass_npz_kwargs(n_genes=3, n_classes=2)
    kwargs["classes"] = np.array(["c0", "c1", "c2"])  # coef.shape[0]=2, 3 classes
    np.savez(bad, **kwargs)
    with pytest.raises(
        ValueError, match="Class label count must match number of coefficient rows"
    ):
        load_lr_model(bad)


def test_load_pipeline_stubs_xgboost_objects_without_importing(tmp_path):
    """A legacy pickle that references an xgboost.* class must unpickle as
    _XgbStub via _SafeUnpickler.find_class, without ever importing xgboost
    (importing it segfaults in the project's conda env). Hand-crafted at the
    pickle-opcode level (not pickle.dump) because the real Pickler refuses to
    reference a class from a module it can't actually import."""
    payload = b"".join([
        b"(",                            # MARK (outer dict)
        b"S'model'\n",                   # key "model"
        b"cxgboost.sklearn\n",           # GLOBAL module
        b"XGBClassifier\n",              # GLOBAL qualname
        b"(t",                           # EMPTY_TUPLE (mark + TUPLE)
        b"R",                            # REDUCE -> calls stubbed class() -> instance
        b"(S'note'\nS'stub-state'\nd",   # build state dict {"note": "stub-state"}
        b"b",                            # BUILD -> instance.__setstate__(state)
        b"d",                            # DICT -> {"model": instance}
        b".",                            # STOP
    ])
    path = tmp_path / "legacy_pipeline.pkl"
    path.write_bytes(payload)

    loaded = load_pipeline(path, trusted=True)

    assert "xgboost" not in sys.modules
    assert isinstance(loaded, dict)
    assert isinstance(loaded["model"], _XgbStub)
    assert loaded["model"]._state == {"note": "stub-state"}


def test_read_matrix_rejects_pickle_by_default(tmp_path):
    path = tmp_path / "input.pkl"
    expected = pd.DataFrame({"g1": [1.0]}, index=["s1"])
    expected.to_pickle(path)

    with pytest.raises(ValueError, match="Pickle expression inputs are disabled"):
        read_matrix(path)

    observed = read_matrix(path, allow_pickle=True)
    pd.testing.assert_frame_equal(observed, expected)


def test_read_matrix_reports_missing_file_as_valueerror(tmp_path):
    missing = tmp_path / "missing.csv"
    with pytest.raises(ValueError, match="expression matrix file not found"):
        read_matrix(missing)


def test_read_matrix_rejects_directory_path(tmp_path):
    with pytest.raises(ValueError, match="path is a directory"):
        read_matrix(tmp_path)


def test_read_matrix_tsv_and_txt_round_trip(tmp_path):
    df = pd.DataFrame({"g1": [1.0, 2.0], "g2": [3.0, 4.0]}, index=["s1", "s2"])

    tsv_path = tmp_path / "matrix.tsv"
    df.to_csv(tsv_path, sep="\t")
    pd.testing.assert_frame_equal(read_matrix(tsv_path), df)

    txt_path = tmp_path / "matrix.txt"
    df.to_csv(txt_path, sep="\t")
    pd.testing.assert_frame_equal(read_matrix(txt_path), df)


def test_read_matrix_parquet_round_trip(tmp_path):
    pytest.importorskip("pyarrow")
    df = pd.DataFrame({"g1": [1.0, 2.0], "g2": [3.0, 4.0]}, index=["s1", "s2"])
    path = tmp_path / "matrix.parquet"
    df.to_parquet(path)
    pd.testing.assert_frame_equal(read_matrix(path), df)


def test_write_json_creates_parent_dirs_and_trailing_newline(tmp_path):
    path = tmp_path / "nested" / "deep" / "report.json"
    payload = {"a": 1, "b": [1, 2, 3]}
    result = write_json(payload, path)

    assert result == path
    assert path.parent.is_dir()
    content = path.read_text(encoding="utf-8")
    assert content.endswith("\n")
    assert json.loads(content) == payload


def test_scoring_cli_reports_missing_input_without_traceback(tmp_path, root):
    output = tmp_path / "missing.scored.csv"
    result = subprocess.run(
        [
            sys.executable,
            "score_tumor_normal.py",
            str(tmp_path / "missing.csv"),
            "-o",
            str(output),
        ],
        cwd=root,
        text=True,
        capture_output=True,
    )

    assert result.returncode != 0
    assert "expression matrix file not found" in result.stderr
    assert "Traceback" not in result.stderr
    assert not output.exists()


def test_binary_cli_rejects_multiclass_weights(root, cancer_type_model):
    with pytest.raises(ValueError, match="requires a binary LR weights file"):
        load_lr_weights(f"{root}/cancer-type-classifier/cancer_type_lr_weights.npz")
    with pytest.raises(ValueError, match="requires a binary model"):
        score_dataframe_lr_weights(pd.DataFrame({"g1": [1.0]}), cancer_type_model)


def test_cancer_type_cli_rejects_missing_weights_file(tmp_path, root):
    result = subprocess.run(
        [
            sys.executable,
            "cancer-type-classifier/predict_cancer_type.py",
            f"{root}/example_input.csv",
            "--weights",
            str(tmp_path / "does_not_exist.npz"),
        ],
        cwd=root,
        text=True,
        capture_output=True,
    )
    assert result.returncode != 0
    assert "weights file not found" in result.stderr
    assert "Traceback" not in result.stderr


def test_cancer_type_cli_rejects_malformed_weights_file(tmp_path, root):
    weights_path = tmp_path / "malformed.npz"
    np.savez(
        weights_path,
        selected_genes=np.array(["a", "b"]),
        scaler_mean=np.array([0.0]),
        scaler_scale=np.array([1.0, 1.0]),
        coef=np.array([[1.0, 1.0], [1.0, 1.0]]),
        intercept=np.array([0.0, 0.0]),
        classes=np.array(["x", "y"]),
    )
    result = subprocess.run(
        [
            sys.executable,
            "cancer-type-classifier/predict_cancer_type.py",
            f"{root}/example_input.csv",
            "--weights",
            str(weights_path),
        ],
        cwd=root,
        text=True,
        capture_output=True,
    )
    assert result.returncode != 0
    assert "must have one value per selected gene" in result.stderr
    assert "Traceback" not in result.stderr


def test_cancer_type_cli_rejects_invalid_matched_values(tmp_path, root, cancer_type_model):
    genes = [str(gene) for gene in cancer_type_model["genes"][:600]]
    input_path = tmp_path / "invalid_cancer_type.csv"
    output_path = tmp_path / "predictions.csv"
    values = {gene: [1.0] for gene in genes}
    values[genes[0]] = ["not_numeric"]
    pd.DataFrame(values, index=["s1"]).to_csv(input_path, index_label="sample")
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


def test_cancer_type_cli_rejects_low_gene_coverage_unless_allowed(tmp_path, root):
    input_path = tmp_path / "no_model_genes.csv"
    output_path = tmp_path / "predictions.csv"
    pd.DataFrame({"NOT_A_MODEL_GENE": [1.0]}, index=["s1"]).to_csv(input_path, index_label="sample")

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
    assert "low model-gene coverage" in result.stderr
    assert "Refusing to write predictions" in result.stderr
    assert not output_path.exists()

    result = subprocess.run(
        [
            sys.executable,
            "cancer-type-classifier/predict_cancer_type.py",
            str(input_path),
            "--out",
            str(output_path),
            "--allow-low-gene-coverage",
        ],
        cwd=root,
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0
    assert output_path.exists()


def test_cancer_type_cli_validates_topk_bounds(root, cancer_type_model):
    n_classes = len(cancer_type_model["classes"])
    result = subprocess.run(
        [sys.executable, "cancer-type-classifier/predict_cancer_type.py",
         f"{root}/example_input.csv", "--topk", "0"],
        cwd=root, text=True, capture_output=True,
    )
    assert result.returncode != 0
    assert "--topk must be >= 1" in result.stderr

    result = subprocess.run(
        [sys.executable, "cancer-type-classifier/predict_cancer_type.py",
         f"{root}/example_input.csv", "--topk", str(n_classes + 1)],
        cwd=root, text=True, capture_output=True,
    )
    assert result.returncode != 0
    assert f"--topk must be <= number of classes ({n_classes})" in result.stderr


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


def test_explain_empty_input_is_rejected():
    weights = {
        "selected_genes": np.array(["g1"]),
        "scaler_mean": np.array([0.0]),
        "scaler_scale": np.array([1.0]),
        "coef": np.array([1.0]),
        "intercept": 0.0,
    }
    df = pd.DataFrame({"g1": []})
    with pytest.raises(ValueError, match="at least one sample"):
        explain_dataframe(df, weights, top_n=1)


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
