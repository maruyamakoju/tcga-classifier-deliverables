"""Standardization (with cross-platform domain adaptation) and logistic scoring."""
from __future__ import annotations

from typing import Any, cast

import numpy as np
import pandas as pd

from .align import (
    align_to_genes_with_report,
    format_alignment_issues,
    format_gene_match_issues,
)
from .io import validate_lr_model
from .validation import validate_expression_matrix, validate_threshold

ADAPT_MODES = ("none", "cohort_zscore", "cohort_center")


def sigmoid(x: Any) -> np.ndarray:
    """Numerically stable logistic sigmoid (no overflow at large |x|)."""
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    pos = x >= 0
    out[pos] = 1.0 / (1.0 + np.exp(-x[pos]))
    e = np.exp(x[~pos])
    out[~pos] = e / (1.0 + e)
    return out


def softmax(logits: Any, axis: int = 1) -> np.ndarray:
    logits = np.asarray(logits, dtype=float)
    if logits.ndim < 1:
        raise ValueError("softmax logits must have at least one dimension")
    if not np.all(np.isfinite(logits)):
        raise ValueError("softmax logits must be finite")
    z = logits - logits.max(axis=axis, keepdims=True)
    e = np.exp(z)
    return e / e.sum(axis=axis, keepdims=True)


def standardize(values: Any, model: dict[str, Any], adapt: str = "none") -> np.ndarray:
    """Standardize aligned feature values before the linear model.

    none           z = (x - train_mean) / train_scale        (deployed)
    cohort_zscore  z = (x - cohort_mean) / cohort_std         (domain adaptation)
    cohort_center  z = (x - cohort_mean) / train_scale        (location-only DA)

    Cohort modes are experimental, transductive transforms: every score depends
    on the other samples in the batch.  They do not restore calibration or
    guarantee that a decision threshold transfers across platforms.
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


def _validate_aligned_values(model: dict[str, Any], values: np.ndarray) -> None:
    if values.ndim != 2:
        raise ValueError("aligned values must be a 2-D samples x genes matrix")
    n_genes = len(model["genes"])
    if values.shape[1] != n_genes:
        raise ValueError(
            f"aligned values have {values.shape[1]} columns; expected {n_genes} model genes"
        )
    if not np.isfinite(values).all():
        raise ValueError("aligned values contain non-finite values")
    if values.shape[0] == 0:
        raise ValueError("aligned values must contain at least one sample")


def _binary_role(value: Any) -> str | None:
    text = str(value).strip().lower()
    if text in {"0", "0.0", "false", "normal", "healthy", "negative", "neg"}:
        return "normal"
    if text in {"1", "1.0", "true", "tumor", "tumour", "cancer", "positive", "pos"}:
        return "tumor"
    return None


def validate_tumor_binary_model(model: dict[str, Any]) -> dict[str, Any]:
    """Validate that binary class 1 really denotes tumor probability."""
    model = validate_lr_model(model)
    if model["kind"] != "binary":
        raise ValueError("tumor-vs-normal scoring requires a binary model")
    roles = [_binary_role(value) for value in model["classes"]]
    if roles != ["normal", "tumor"]:
        raise ValueError(
            "binary class order must identify normal as class 0 and tumor as class 1; "
            f"received {model['classes'].tolist()}"
        )
    return model


def predict_proba_from_aligned(model: dict[str, Any], values: Any, adapt: str = "none", validate_values: bool = True) -> np.ndarray:
    """Score already-aligned feature values.

    ``values`` must be samples x model genes in ``model["genes"]`` order.
    This helper lets CLIs align once, then reuse the same matrix for reporting
    matched/missing genes and probabilities.
    """
    model = validate_lr_model(model)
    values = np.asarray(values, dtype=float)
    if validate_values:
        _validate_aligned_values(model, values)
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            z = standardize(values, model, adapt=adapt)
            if not np.all(np.isfinite(z)):
                raise ValueError("standardized expression values are non-finite")
            coef = model["coef"]
            if model["kind"] == "binary":
                logits = z @ coef + model["intercept"]
            else:
                logits = z @ coef.T + model["intercept"]
    except FloatingPointError as exc:
        raise ValueError(f"numeric overflow while scoring expression values: {exc}") from exc
    if not np.all(np.isfinite(logits)):
        raise ValueError("model logits are non-finite")
    proba = sigmoid(logits) if model["kind"] == "binary" else softmax(logits, axis=1)
    if not np.all(np.isfinite(proba)) or np.any((proba < 0) | (proba > 1)):
        raise ValueError("model probabilities must be finite and between 0 and 1")
    if model["kind"] == "multiclass" and not np.allclose(
        proba.sum(axis=1), 1.0, rtol=0.0, atol=1e-12
    ):
        raise ValueError("multiclass probability rows must sum to 1")
    return proba


def _raise_for_invalid_alignment(report: dict[str, Any], max_invalid_cell_fraction: float) -> None:
    message = format_alignment_issues(
        report,
        max_invalid_cell_fraction=max_invalid_cell_fraction,
    )
    if message:
        raise ValueError(message)


def predict_proba(
    model: dict[str, Any],
    X: pd.DataFrame,
    adapt: str = "none",
    max_invalid_cell_fraction: float = 0.0,
    allow_invalid_values: bool = False,
    return_alignment_report: bool = False,
    min_model_gene_match_rate: float = 0.5,
    allow_low_gene_coverage: bool = False,
) -> np.ndarray | tuple[np.ndarray, dict[str, Any]]:
    """Score an expression DataFrame.

    Returns P(positive class) as a 1-D array for a binary model, or an
    (n_samples, n_classes) probability matrix for a multi-class model.
    Missing model genes are imputed at the training mean. Invalid values in
    matched model-gene cells raise ``ValueError`` by default; pass
    ``allow_invalid_values=True`` only after reviewing mean imputation.
    Inputs matching fewer than ``min_model_gene_match_rate`` model genes also
    fail by default; bypassing that guard requires ``allow_low_gene_coverage``.
    """
    model = validate_lr_model(model)
    X = validate_expression_matrix(X)
    min_model_gene_match_rate = validate_threshold(
        min_model_gene_match_rate, "min_model_gene_match_rate"
    )
    max_invalid_cell_fraction = validate_threshold(
        max_invalid_cell_fraction, "max_invalid_cell_fraction"
    )
    values, report = align_to_genes_with_report(
        X, model["genes"], impute_mean=model["mean"]
    )
    if not allow_low_gene_coverage:
        message = format_gene_match_issues(
            report, min_match_rate=min_model_gene_match_rate
        )
        if message:
            raise ValueError(message)
    if not allow_invalid_values:
        _raise_for_invalid_alignment(report, max_invalid_cell_fraction)
    proba = predict_proba_from_aligned(model, values, adapt=adapt)
    if return_alignment_report:
        return proba, report
    return proba


def score_binary_dataframe(
    model: dict[str, Any],
    X: pd.DataFrame,
    threshold: float = 0.5,
    adapt: str = "none",
    positive_label: str = "tumor",
    negative_label: str = "normal",
    round_digits: int | None = None,
    max_invalid_cell_fraction: float = 0.0,
    allow_invalid_values: bool = False,
    return_alignment_report: bool = False,
    min_model_gene_match_rate: float = 0.5,
    allow_low_gene_coverage: bool = False,
) -> tuple[Any, ...]:
    """Return the stable public scoring CSV shape for a binary model.

    Returns ``(dataframe, n_matched, missing)``. The DataFrame columns are
    ``sample,tumor_probability,call`` to preserve the release contract. Public
    probabilities retain full float precision so the serialized probability
    and hard call cannot disagree at the threshold.
    When ``return_alignment_report`` is true, returns
    ``(dataframe, n_matched, missing, report)``.
    """
    model = validate_tumor_binary_model(model)
    X = validate_expression_matrix(X)
    threshold = validate_threshold(threshold)
    min_model_gene_match_rate = validate_threshold(
        min_model_gene_match_rate, "min_model_gene_match_rate"
    )
    max_invalid_cell_fraction = validate_threshold(
        max_invalid_cell_fraction, "max_invalid_cell_fraction"
    )
    values, report = align_to_genes_with_report(
        X, model["genes"], impute_mean=model["mean"]
    )
    n_matched = report["n_matched_genes"]
    missing = report["missing_genes"]
    if not allow_low_gene_coverage:
        message = format_gene_match_issues(
            report, min_match_rate=min_model_gene_match_rate
        )
        if message:
            raise ValueError(message)
    if not allow_invalid_values:
        _raise_for_invalid_alignment(report, max_invalid_cell_fraction)
    proba = predict_proba_from_aligned(model, values, adapt=adapt)
    calls = np.where(proba >= threshold, positive_label, negative_label)

    displayed_proba = proba if round_digits is None else np.round(proba, round_digits)
    if round_digits is not None:
        calls = np.where(displayed_proba >= threshold, positive_label, negative_label)
    result = pd.DataFrame({
        "sample": X.index,
        "tumor_probability": displayed_proba,
        "call": calls,
    })
    if return_alignment_report:
        return result, n_matched, missing, report
    return result, n_matched, missing


def predict(
    model: dict[str, Any],
    X: pd.DataFrame,
    adapt: str = "none",
    threshold: float = 0.5,
    max_invalid_cell_fraction: float = 0.0,
    allow_invalid_values: bool = False,
    min_model_gene_match_rate: float = 0.5,
    allow_low_gene_coverage: bool = False,
) -> np.ndarray:
    """Return hard class calls.

    Binary: array of the two class labels using `threshold` on P(positive).
    Multi-class: array of argmax class labels.
    """
    model = validate_lr_model(model)
    threshold = validate_threshold(threshold)
    proba = cast(np.ndarray, predict_proba(
        model,
        X,
        adapt=adapt,
        max_invalid_cell_fraction=max_invalid_cell_fraction,
        allow_invalid_values=allow_invalid_values,
        min_model_gene_match_rate=min_model_gene_match_rate,
        allow_low_gene_coverage=allow_low_gene_coverage,
    ))
    classes = model["classes"]
    if model["kind"] == "binary":
        return np.where(proba >= threshold, classes[1], classes[0])
    return classes[proba.argmax(axis=1)]
