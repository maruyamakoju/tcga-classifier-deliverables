"""Dependency-light classification metrics with strict input contracts.

The release deliberately keeps scoring independent of scikit-learn.  These
helpers therefore implement the small set of metrics used by calibration and
external validation, while failing closed on malformed arrays instead of
silently truncating, sorting NaNs, or accepting non-binary labels.
"""

from __future__ import annotations

from typing import Any

import numpy as np


def _as_1d(values: Any, name: str) -> np.ndarray:
    array = np.asarray(values)
    if array.ndim != 1:
        raise ValueError(f"{name} must be a one-dimensional array")
    if array.size == 0:
        raise ValueError(f"{name} must contain at least one value")
    return array


def _validate_same_length(left: Any, right: Any, left_name: str = "y_true", right_name: str = "y_pred") -> tuple[np.ndarray, np.ndarray]:
    left = _as_1d(left, left_name)
    right = _as_1d(right, right_name)
    if left.size != right.size:
        raise ValueError(
            f"{left_name} and {right_name} must have the same length "
            f"({left.size} != {right.size})"
        )
    return left, right


def _has_missing_label(values: np.ndarray) -> bool:
    for value in values.tolist():
        if value is None:
            return True
        try:
            if bool(np.asarray(value != value).item()):  # NaN/NaT without pandas.
                return True
        except (TypeError, ValueError):
            # pandas.NA and non-scalar label objects have no unambiguous truth
            # value; neither is a valid class label for these metrics.
            return True
    return False


def _classification_arrays(y_true: Any, y_pred: Any) -> tuple[np.ndarray, np.ndarray]:
    y, pred = _validate_same_length(y_true, y_pred)
    if _has_missing_label(y) or _has_missing_label(pred):
        raise ValueError("class labels must not contain missing values")
    return y, pred


def _binary_arrays(y_true: Any, scores: Any, require_both_classes: bool = False) -> tuple[np.ndarray, np.ndarray]:
    y_raw, score_raw = _validate_same_length(y_true, scores, "y_true", "scores")
    try:
        y_numeric = y_raw.astype(float)
    except (TypeError, ValueError) as exc:
        raise ValueError("y_true must contain only binary 0/1 labels") from exc
    if not np.all(np.isfinite(y_numeric)) or not np.all(np.isin(y_numeric, [0.0, 1.0])):
        raise ValueError("y_true must contain only binary 0/1 labels")
    try:
        score_values = score_raw.astype(float)
    except (TypeError, ValueError) as exc:
        raise ValueError("scores must contain only finite numeric values") from exc
    if not np.all(np.isfinite(score_values)):
        raise ValueError("scores must contain only finite numeric values")
    y = y_numeric.astype(int)
    if require_both_classes and np.unique(y).size != 2:
        raise ValueError("y_true must contain both binary classes")
    return y, score_values


def _finite_threshold(threshold: Any) -> float:
    try:
        value = float(threshold)
    except (TypeError, ValueError) as exc:
        raise ValueError("threshold must be a finite number") from exc
    if not np.isfinite(value):
        raise ValueError("threshold must be a finite number")
    return value


def roc_auc(y_true: Any, scores: Any) -> float:
    """Binary ROC AUC via the rank statistic (Mann-Whitney U).

    A single-class target has no defined ROC AUC and returns ``nan``.  Invalid
    labels, non-finite scores, shape mismatches, and empty inputs raise a clear
    ``ValueError``.
    """
    y, s = _binary_arrays(y_true, scores)
    n_pos = int((y == 1).sum())
    n_neg = int((y == 0).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    order = np.argsort(s, kind="mergesort")
    ranks = np.empty(len(s), dtype=float)
    s_sorted = s[order]
    ranks_sorted = np.arange(1, len(s) + 1, dtype=float)
    i = 0
    while i < len(s_sorted):
        j = i
        while j + 1 < len(s_sorted) and s_sorted[j + 1] == s_sorted[i]:
            j += 1
        if j > i:
            ranks_sorted[i:j + 1] = (i + 1 + j + 1) / 2.0
        i = j + 1
    ranks[order] = ranks_sorted
    return float((ranks[y == 1].sum() - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg))


def average_precision(y_true: Any, scores: Any) -> float:
    """Binary average precision, grouping tied scores at one threshold.

    This matches ``sklearn.metrics.average_precision_score`` for a target that
    contains both classes.  A target with no positive samples has no useful AP
    for this project and returns ``nan``.
    """
    y, s = _binary_arrays(y_true, scores)
    n_pos = int((y == 1).sum())
    if n_pos == 0:
        return float("nan")

    order = np.argsort(-s, kind="mergesort")
    sorted_y = y[order]
    sorted_scores = s[order]
    tp = np.cumsum(sorted_y == 1)
    fp = np.cumsum(sorted_y == 0)
    group_ends = np.r_[np.flatnonzero(np.diff(sorted_scores)), len(sorted_scores) - 1]
    precision = tp[group_ends] / (tp[group_ends] + fp[group_ends])
    recall = tp[group_ends] / n_pos
    recall_increments = np.diff(np.r_[0.0, recall])
    return float(np.sum(recall_increments * precision))


def accuracy(y_true: Any, y_pred: Any) -> float:
    y, pred = _classification_arrays(y_true, y_pred)
    return float(np.mean(y == pred))


def balanced_accuracy(y_true: Any, y_pred: Any) -> float:
    y, pred = _classification_arrays(y_true, y_pred)
    recalls = [np.mean(pred[y == label] == label) for label in np.unique(y)]
    return float(np.mean(recalls))


def confusion_matrix(y_true: Any, y_pred: Any, labels: Any = None) -> tuple[np.ndarray, list[Any]]:
    y, pred = _classification_arrays(y_true, y_pred)
    if labels is None:
        labels = np.unique(np.concatenate([y, pred])).tolist()
    else:
        labels = list(labels)
        if not labels:
            raise ValueError("labels must contain at least one class")
        if len(set(labels)) != len(labels):
            raise ValueError("labels must not contain duplicates")
        if _has_missing_label(np.asarray(labels, dtype=object)):
            raise ValueError("labels must not contain missing values")
    index = {label: i for i, label in enumerate(labels)}
    matrix = np.zeros((len(labels), len(labels)), dtype=int)
    for truth, prediction in zip(y, pred, strict=True):
        if truth in index and prediction in index:
            matrix[index[truth], index[prediction]] += 1
    return matrix, labels


def per_class_prf(y_true: Any, y_pred: Any, labels: Any = None) -> list[dict[str, Any]]:
    """Per-class precision/recall/F1/support. Returns a list of dicts."""
    matrix, labels = confusion_matrix(y_true, y_pred, labels)
    output = []
    for i, label in enumerate(labels):
        tp = matrix[i, i]
        fp = matrix[:, i].sum() - tp
        fn = matrix[i, :].sum() - tp
        precision = tp / (tp + fp) if (tp + fp) else 0.0
        recall = tp / (tp + fn) if (tp + fn) else 0.0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
        output.append({
            "label": label,
            "support": int(matrix[i, :].sum()),
            "precision": precision,
            "recall": recall,
            "f1": f1,
        })
    return output


def macro_f1(y_true: Any, y_pred: Any, labels: Any = None) -> float:
    return float(np.mean([row["f1"] for row in per_class_prf(y_true, y_pred, labels)]))


def confusion_at(y_true: Any, scores: Any, threshold: Any) -> tuple[int, int, int, int]:
    """Binary confusion counts ``(tn, fp, fn, tp)`` at a threshold."""
    y, score_values = _binary_arrays(y_true, scores)
    threshold = _finite_threshold(threshold)
    pred = (score_values >= threshold).astype(int)
    tp = int(((pred == 1) & (y == 1)).sum())
    tn = int(((pred == 0) & (y == 0)).sum())
    fp = int(((pred == 1) & (y == 0)).sum())
    fn = int(((pred == 0) & (y == 1)).sum())
    return tn, fp, fn, tp


def classification_metrics(y_true: Any, scores: Any, threshold: Any = 0.5, name: str | None = None) -> dict[str, Any]:
    """Binary threshold metrics plus threshold-independent ROC AUC."""
    y, score_values = _binary_arrays(y_true, scores)
    threshold = _finite_threshold(threshold)
    pred = (score_values >= threshold).astype(int)
    tn, fp, fn, tp = confusion_at(y, score_values, threshold)
    n = tn + fp + fn + tp
    acc = (tp + tn) / n
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    specificity = tn / (tn + fp) if (tn + fp) else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
    output = {
        "threshold": threshold,
        "accuracy": acc,
        "precision": precision,
        "recall": recall,
        "sensitivity": recall,
        "specificity": specificity,
        "f1": f1,
        # Match balanced_accuracy(): average recall across classes actually
        # present in y_true. ROC AUC remains undefined for a single class, but
        # threshold metrics still have an unambiguous value.
        "balanced_accuracy": balanced_accuracy(y, pred),
        "auc": roc_auc(y, score_values),
        "tn": tn,
        "fp": fp,
        "fn": fn,
        "tp": tp,
    }
    if name is not None:
        output["name"] = name
    return output


def youden_threshold(y_true: Any, scores: Any) -> dict[str, Any]:
    """Choose a probability cutoff by J, accuracy, then lower-threshold tie-break.

    Candidate cutoffs include the closed probability-domain boundaries.  This
    matters when every score is tied: under the public ``score >= threshold``
    rule, threshold 1 can be a better all-negative classifier than the sole
    observed score even though both have the same Youden J.
    """
    y, score_values = _binary_arrays(y_true, scores, require_both_classes=True)
    if np.any((score_values < 0.0) | (score_values > 1.0)):
        raise ValueError("scores must be probabilities in [0, 1]")
    best_key = None
    best_metrics = None
    candidates = np.unique(np.concatenate([score_values, np.array([0.0, 1.0])]))
    for threshold in candidates:
        metrics = classification_metrics(y, score_values, threshold)
        key = (
            metrics["sensitivity"] + metrics["specificity"] - 1,
            metrics["accuracy"],
            -threshold,
        )
        if best_key is None or key > best_key:
            best_key = key
            best_metrics = metrics
    if best_metrics is None:  # _binary_arrays already rejects an empty input.
        raise RuntimeError("could not choose a threshold")
    return best_metrics


def threshold_sweep(y_true: Any, scores: Any, thresholds: Any) -> list[dict[str, Any]]:
    """Return classification metrics across the requested thresholds."""
    return [classification_metrics(y_true, scores, threshold) for threshold in thresholds]
