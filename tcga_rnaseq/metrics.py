"""Dependency-light classification metrics (numpy only, no scikit-learn)."""
import numpy as np


def roc_auc(y_true, scores):
    """Binary ROC AUC via the rank statistic (Mann-Whitney U)."""
    y = np.asarray(y_true).astype(int)
    s = np.asarray(scores, dtype=float)
    n_pos = int((y == 1).sum())
    n_neg = int((y == 0).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    order = np.argsort(s, kind="mergesort")
    ranks = np.empty(len(s), dtype=float)
    # average ranks for ties
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


def accuracy(y_true, y_pred):
    return float(np.mean(np.asarray(y_true) == np.asarray(y_pred)))


def balanced_accuracy(y_true, y_pred):
    y = np.asarray(y_true)
    p = np.asarray(y_pred)
    recalls = [np.mean(p[y == c] == c) for c in np.unique(y)]
    return float(np.mean(recalls))


def confusion_matrix(y_true, y_pred, labels=None):
    y = np.asarray(y_true)
    p = np.asarray(y_pred)
    if labels is None:
        labels = np.unique(np.concatenate([y, p]))
    idx = {c: i for i, c in enumerate(labels)}
    m = np.zeros((len(labels), len(labels)), dtype=int)
    for t, pr in zip(y, p):
        if t in idx and pr in idx:
            m[idx[t], idx[pr]] += 1
    return m, list(labels)


def per_class_prf(y_true, y_pred, labels=None):
    """Per-class precision/recall/F1/support. Returns list of dicts."""
    m, labels = confusion_matrix(y_true, y_pred, labels)
    out = []
    for i, c in enumerate(labels):
        tp = m[i, i]
        fp = m[:, i].sum() - tp
        fn = m[i, :].sum() - tp
        prec = tp / (tp + fp) if (tp + fp) else 0.0
        rec = tp / (tp + fn) if (tp + fn) else 0.0
        f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
        out.append({"label": c, "support": int(m[i, :].sum()),
                    "precision": prec, "recall": rec, "f1": f1})
    return out


def macro_f1(y_true, y_pred, labels=None):
    return float(np.mean([r["f1"] for r in per_class_prf(y_true, y_pred, labels)]))


def confusion_at(y_true, scores, threshold):
    """Binary confusion counts (tn, fp, fn, tp) at a probability threshold."""
    y = np.asarray(y_true).astype(int)
    pred = (np.asarray(scores, dtype=float) >= threshold).astype(int)
    tp = int(((pred == 1) & (y == 1)).sum())
    tn = int(((pred == 0) & (y == 0)).sum())
    fp = int(((pred == 1) & (y == 0)).sum())
    fn = int(((pred == 0) & (y == 1)).sum())
    return tn, fp, fn, tp


def classification_metrics(y_true, scores, threshold=0.5, name=None):
    """Binary metrics at a threshold: accuracy, precision, recall/sensitivity,
    specificity, f1, plus AUC (threshold-independent) and confusion counts."""
    tn, fp, fn, tp = confusion_at(y_true, scores, threshold)
    n = tn + fp + fn + tp
    acc = (tp + tn) / n if n else float("nan")
    prec = tp / (tp + fp) if (tp + fp) else 0.0
    rec = tp / (tp + fn) if (tp + fn) else 0.0
    spec = tn / (tn + fp) if (tn + fp) else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
    out = {"threshold": float(threshold), "accuracy": acc, "precision": prec,
           "recall": rec, "sensitivity": rec, "specificity": spec, "f1": f1,
           "balanced_accuracy": (rec + spec) / 2, "auc": roc_auc(y_true, scores),
           "tn": tn, "fp": fp, "fn": fn, "tp": tp}
    if name is not None:
        out["name"] = name
    return out


def youden_threshold(y_true, scores):
    """Single canonical Youden-J threshold: scan unique scores, maximize
    (sensitivity + specificity - 1), tie-break on accuracy then lower threshold."""
    y = np.asarray(y_true).astype(int)
    s = np.asarray(scores, dtype=float)
    best = None
    for thr in np.unique(s):
        m = classification_metrics(y, s, thr)
        key = (m["sensitivity"] + m["specificity"] - 1, m["accuracy"], -thr)
        if best is None or key > best[0]:
            best = (key, m)
    return best[1] if best else {"threshold": 0.5}


def threshold_sweep(y_true, scores, thresholds):
    """Return a list of classification_metrics dicts across thresholds."""
    return [classification_metrics(y_true, scores, t) for t in thresholds]
