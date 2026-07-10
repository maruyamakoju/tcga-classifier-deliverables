"""Shared CLI input validation for tcga_rnaseq.

These were previously duplicated (or cross-imported from one leaf CLI script,
calibrate_threshold.py, into five others) across the scoring CLI suite.
"""
import numpy as np


def validate_threshold(value, name="threshold"):
    """Validate a 0-1 probability/fraction argument shared by the scoring CLIs.

    Used for --threshold, --max-invalid-cell-fraction, --min-model-gene-match-rate,
    and similar 0-1 range arguments across the CLI suite.
    """
    value = float(value)
    if not np.isfinite(value) or not 0 <= value <= 1:
        raise ValueError(f"{name} must be between 0 and 1")
    return value


def normalize_label(value):
    """Map a tumor/normal label value (string or 0/1) to a binary int."""
    text = str(value).strip().lower()
    if text in {"1", "tumor", "tumour", "primary tumor", "cancer", "positive", "pos", "true"}:
        return 1
    if text in {"0", "normal", "solid tissue normal", "healthy", "negative", "neg", "false"}:
        return 0
    raise ValueError(f"Unrecognized label: {value!r}")


def sample_key(series):
    """Trimmed, string-typed sample identifiers; rejects any empty/NaN entry."""
    keys = series.astype(str).str.strip()
    if series.isna().any() or (keys == "").any():
        raise ValueError("sample identifiers must be non-empty")
    return keys


def require_unique_samples(df, sample_col, source_name):
    """sample_key(df[sample_col]), raising ValueError on any duplicate."""
    keys = sample_key(df[sample_col])
    duplicated = sorted(keys[keys.duplicated()].unique())
    if duplicated:
        preview = ", ".join(duplicated[:5])
        raise ValueError(f"{source_name} contains duplicate sample IDs: {preview}")
    return keys
