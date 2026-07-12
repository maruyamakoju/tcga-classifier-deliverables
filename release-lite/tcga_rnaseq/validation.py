"""Shared public input validation for :mod:`tcga_rnaseq`.

These were previously duplicated (or cross-imported from one leaf CLI script,
calibrate_threshold.py, into five others) across the scoring CLI suite.
"""
from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd


def validate_threshold(value: Any, name: str = "threshold") -> float:
    """Validate a 0-1 probability/fraction argument shared by the scoring CLIs.

    Used for --threshold, --max-invalid-cell-fraction, --min-model-gene-match-rate,
    and similar 0-1 range arguments across the CLI suite.
    """
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be a numeric value between 0 and 1")
    try:
        value = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a numeric value between 0 and 1") from exc
    if not np.isfinite(value) or not 0 <= value <= 1:
        raise ValueError(f"{name} must be between 0 and 1")
    return value


def normalize_label(value: Any) -> int:
    """Map a tumor/normal label value (string or 0/1) to a binary int."""
    if isinstance(value, (int, float, np.integer, np.floating, bool, np.bool_)):
        numeric = float(value)
        if np.isfinite(numeric) and numeric in {0.0, 1.0}:
            return int(numeric)
    text = str(value).strip().lower()
    if text in {"1", "tumor", "tumour", "primary tumor", "cancer", "positive", "pos", "true"}:
        return 1
    if text in {"0", "normal", "solid tissue normal", "healthy", "negative", "neg", "false"}:
        return 0
    raise ValueError(f"Unrecognized label: {value!r}")


def sample_key(series: Any, source_name: str = "sample identifiers") -> pd.Series:
    """Return canonical string sample identifiers.

    Public outputs promise non-empty identifiers without leading/trailing
    whitespace.  Reject padding instead of silently trimming it: silently
    changing an identifier can join the wrong label or overwrite another
    sample's result.
    """
    raw = pd.Series(series, dtype=object)
    keys = raw.astype(str).str.strip()
    if raw.isna().any() or (keys == "").any():
        raise ValueError(f"{source_name} must be non-empty")
    padded = raw.astype(str) != keys
    if padded.any():
        examples = ", ".join(repr(value) for value in raw[padded].head(5))
        raise ValueError(
            f"{source_name} must not contain leading or trailing whitespace"
            + (f": {examples}" if examples else "")
        )
    return keys


def _reject_duplicate_keys(keys: pd.Series, source_name: str, kind: str) -> None:
    """Raise ValueError if ``keys`` (a canonical string Series) has duplicates."""
    duplicated = sorted(keys[keys.duplicated()].unique())
    if duplicated:
        preview = ", ".join(duplicated[:5])
        raise ValueError(f"{source_name} contains duplicate {kind}: {preview}")


def require_unique_samples(df: pd.DataFrame, sample_col: str, source_name: str) -> pd.Series:
    """sample_key(df[sample_col]), raising ValueError on any duplicate."""
    if sample_col not in df.columns:
        raise ValueError(f"{source_name} must contain {sample_col!r}")
    keys = sample_key(df[sample_col], f"{source_name} sample identifiers")
    _reject_duplicate_keys(keys, source_name, "sample IDs")
    return keys


def validate_expression_matrix(df: Any, source_name: str = "expression matrix") -> pd.DataFrame:
    """Validate and canonicalize a samples-by-genes expression DataFrame.

    The returned frame is a shallow copy with string sample/gene axes.  Empty
    matrices, hierarchical axes, invalid sample IDs, duplicate samples, and
    blank/padded/duplicate gene columns are rejected before any scoring occurs.
    Version-colliding Ensembl columns are checked by the alignment layer, where
    the model's matching rules are available.
    """
    if not isinstance(df, pd.DataFrame):
        raise ValueError(f"{source_name} must be a pandas DataFrame")
    if isinstance(df.index, pd.MultiIndex) or isinstance(df.columns, pd.MultiIndex):
        raise ValueError(f"{source_name} must use flat sample and gene identifiers")
    if df.shape[0] == 0:
        raise ValueError(f"{source_name} must contain at least one sample")
    if df.shape[1] == 0:
        raise ValueError(f"{source_name} must contain at least one gene column")

    sample_keys = sample_key(df.index, f"{source_name} sample identifiers")
    _reject_duplicate_keys(sample_keys, source_name, "sample IDs")

    raw_columns = pd.Series(list(df.columns), dtype=object)
    gene_keys = raw_columns.astype(str).str.strip()
    if raw_columns.isna().any() or (gene_keys == "").any():
        raise ValueError(f"{source_name} gene identifiers must be non-empty")
    if (raw_columns.astype(str) != gene_keys).any():
        raise ValueError(
            f"{source_name} gene identifiers must not contain leading or trailing whitespace"
        )
    _reject_duplicate_keys(gene_keys, source_name, "gene columns")

    result = df.copy(deep=False)
    result.index = pd.Index(sample_keys.to_numpy(dtype=str), name=df.index.name)
    result.columns = pd.Index(gene_keys.to_numpy(dtype=str), name=df.columns.name)
    return result
