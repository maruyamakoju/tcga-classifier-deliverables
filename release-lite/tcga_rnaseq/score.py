"""Standardization (with cross-platform domain adaptation) and logistic scoring."""
import numpy as np

from .align import align_to_genes, align_to_genes_with_report

ADAPT_MODES = ("none", "cohort_zscore", "cohort_center")


def sigmoid(x):
    """Numerically stable logistic sigmoid (no overflow at large |x|)."""
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    pos = x >= 0
    out[pos] = 1.0 / (1.0 + np.exp(-x[pos]))
    e = np.exp(x[~pos])
    out[~pos] = e / (1.0 + e)
    return out


def softmax(logits, axis=1):
    logits = np.asarray(logits, dtype=float)
    z = logits - logits.max(axis=axis, keepdims=True)
    e = np.exp(z)
    return e / e.sum(axis=axis, keepdims=True)


def standardize(values, model, adapt="none"):
    """Standardize aligned feature values before the linear model.

    none           z = (x - train_mean) / train_scale        (deployed)
    cohort_zscore  z = (x - cohort_mean) / cohort_std         (domain adaptation)
    cohort_center  z = (x - cohort_mean) / train_scale        (location-only DA)

    Cohort modes realign a foreign-pipeline batch onto the training marginal;
    see cross-platform-adaptation/ for when to use them.
    """
    v = np.asarray(values, dtype=float)
    if adapt == "none":
        return (v - model["mean"]) / model["scale"]
    if adapt not in ADAPT_MODES:
        raise ValueError(f"unknown adapt mode {adapt!r}; choose from {ADAPT_MODES}")
    mu = np.nanmean(v, axis=0)
    if adapt == "cohort_center":
        return (v - mu) / model["scale"]
    sd = np.nanstd(v, axis=0)
    sd[sd == 0] = 1.0
    return (v - mu) / sd  # cohort_zscore


def predict_proba_from_aligned(model, values, adapt="none"):
    """Score already-aligned feature values.

    ``values`` must be samples x model genes in ``model["genes"]`` order.
    This helper lets CLIs align once, then reuse the same matrix for reporting
    matched/missing genes and probabilities.
    """
    z = standardize(values, model, adapt=adapt)
    coef = model["coef"]
    if model["kind"] == "binary":
        return sigmoid(z @ coef + model["intercept"])
    return softmax(z @ coef.T + model["intercept"], axis=1)


def predict_proba(model, X, adapt="none"):
    """Score an expression DataFrame.

    Returns P(positive class) as a 1-D array for a binary model, or an
    (n_samples, n_classes) probability matrix for a multi-class model.
    Missing model genes are imputed at the training mean.
    """
    values, _n_matched, _missing = align_to_genes(X, model["genes"], impute_mean=model["mean"])
    return predict_proba_from_aligned(model, values, adapt=adapt)


def score_binary_dataframe(
    model,
    X,
    threshold=0.5,
    adapt="none",
    positive_label="tumor",
    negative_label="normal",
    round_digits=6,
    return_alignment_report=False,
):
    """Return the stable public scoring CSV shape for a binary model.

    Returns ``(dataframe, n_matched, missing)``. The DataFrame columns are
    ``sample,tumor_probability,call`` to preserve the release contract.
    When ``return_alignment_report`` is true, returns
    ``(dataframe, n_matched, missing, report)``.
    """
    if model["kind"] != "binary":
        raise ValueError("score_binary_dataframe requires a binary model")
    values, report = align_to_genes_with_report(
        X, model["genes"], impute_mean=model["mean"]
    )
    n_matched = report["n_matched_genes"]
    missing = report["missing_genes"]
    proba = predict_proba_from_aligned(model, values, adapt=adapt)
    calls = np.where(proba >= threshold, positive_label, negative_label)
    import pandas as pd

    result = pd.DataFrame({
        "sample": X.index,
        "tumor_probability": np.round(proba, round_digits),
        "call": calls,
    })
    if return_alignment_report:
        return result, n_matched, missing, report
    return result, n_matched, missing


def predict(model, X, adapt="none", threshold=0.5):
    """Return hard class calls.

    Binary: array of the two class labels using `threshold` on P(positive).
    Multi-class: array of argmax class labels.
    """
    proba = predict_proba(model, X, adapt=adapt)
    classes = model["classes"]
    if model["kind"] == "binary":
        return np.where(proba >= threshold, classes[1], classes[0])
    return classes[proba.argmax(axis=1)]
