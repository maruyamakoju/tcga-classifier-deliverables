"""Failure contracts and differential tests for dependency-light metrics."""

import numpy as np
import pytest

from tcga_rnaseq import metrics as M


@pytest.mark.parametrize(
    ("function", "args", "message"),
    [
        (M.accuracy, ([0, 1], [0]), "same length"),
        (M.confusion_matrix, ([0, 1], [0]), "same length"),
        (M.roc_auc, ([0, 1], [0.2]), "same length"),
        (M.roc_auc, ([], []), "at least one"),
        (M.roc_auc, ([0, 2], [0.2, 0.8]), "binary 0/1"),
        (M.roc_auc, ([0, 1], [0.2, np.nan]), "finite numeric"),
        (M.confusion_at, ([0, 1], [0.2, 0.8], np.nan), "finite number"),
    ],
)
def test_metrics_reject_malformed_arrays(function, args, message):
    with pytest.raises(ValueError, match=message):
        function(*args)


def test_multiclass_metrics_reject_missing_labels_and_duplicate_label_order():
    with pytest.raises(ValueError, match="missing"):
        M.accuracy(["a", None], ["a", "b"])
    with pytest.raises(ValueError, match="duplicates"):
        M.confusion_matrix(["a"], ["a"], labels=["a", "a"])


def test_youden_requires_both_classes():
    with pytest.raises(ValueError, match="both binary classes"):
        M.youden_threshold([0, 0], [0.1, 0.2])


def test_youden_evaluates_probability_boundaries_for_tied_scores():
    best = M.youden_threshold([0, 0, 0, 1], [0.5, 0.5, 0.5, 0.5])

    assert best["threshold"] == pytest.approx(1.0)
    assert best["accuracy"] == pytest.approx(0.75)
    assert best["specificity"] == pytest.approx(1.0)
    assert best["sensitivity"] == pytest.approx(0.0)


def test_youden_rejects_non_probability_scores():
    with pytest.raises(ValueError, match=r"probabilities in \[0, 1\]"):
        M.youden_threshold([0, 1], [-1.0, 2.0])


def test_single_class_auc_and_average_precision_are_explicitly_undefined():
    assert np.isnan(M.roc_auc([0, 0], [0.1, 0.2]))
    assert np.isnan(M.average_precision([0, 0], [0.1, 0.2]))


def test_single_class_balanced_accuracy_has_one_consistent_definition():
    y = np.ones(3, dtype=int)
    scores = np.array([0.6, 0.7, 0.4])
    predictions = (scores >= 0.5).astype(int)

    expected = M.balanced_accuracy(y, predictions)
    assert expected == pytest.approx(2 / 3)
    assert M.classification_metrics(y, scores)["balanced_accuracy"] == pytest.approx(
        expected
    )


def test_average_precision_handles_ties():
    sklearn_metrics = pytest.importorskip("sklearn.metrics")
    y = np.array([0, 1, 0, 1, 1, 0])
    scores = np.array([0.9, 0.9, 0.4, 0.4, 0.1, 0.1])
    assert M.average_precision(y, scores) == pytest.approx(
        sklearn_metrics.average_precision_score(y, scores), abs=1e-12
    )


def test_seeded_binary_metrics_match_sklearn():
    sklearn_metrics = pytest.importorskip("sklearn.metrics")
    rng = np.random.default_rng(20260712)
    for n_samples in [2, 3, 7, 25, 100]:
        for _ in range(25):
            y = rng.integers(0, 2, n_samples)
            if np.unique(y).size != 2:
                continue
            scores = np.round(rng.random(n_samples), int(rng.integers(0, 7)))
            threshold = float(rng.random())
            pred = (scores >= threshold).astype(int)
            actual = M.classification_metrics(y, scores, threshold)

            assert actual["auc"] == pytest.approx(
                sklearn_metrics.roc_auc_score(y, scores), abs=1e-12
            )
            assert M.average_precision(y, scores) == pytest.approx(
                sklearn_metrics.average_precision_score(y, scores), abs=1e-12
            )
            assert actual["accuracy"] == pytest.approx(
                sklearn_metrics.accuracy_score(y, pred), abs=1e-12
            )
            assert actual["precision"] == pytest.approx(
                sklearn_metrics.precision_score(y, pred, zero_division=0), abs=1e-12
            )
            assert actual["recall"] == pytest.approx(
                sklearn_metrics.recall_score(y, pred, zero_division=0), abs=1e-12
            )
            assert actual["f1"] == pytest.approx(
                sklearn_metrics.f1_score(y, pred, zero_division=0), abs=1e-12
            )
