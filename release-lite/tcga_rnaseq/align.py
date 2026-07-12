"""Gene-column alignment for tcga_rnaseq."""
from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd


def _preview(values: Any, limit: int = 5) -> str:
    text = [str(value) for value in values[:limit]]
    suffix = ", ..." if len(values) > limit else ""
    return ", ".join(text) + suffix


def strip_version(gene_id: Any) -> str:
    """Drop a trailing numeric version suffix only.

    ``ENSG00000005.6`` becomes ``ENSG00000005`` while identifiers containing
    meaningful dots, such as ``HLA.DRA`` or ``GENE.A.2``, retain everything
    except an actual final numeric suffix.
    """
    text = str(gene_id)
    base, separator, suffix = text.rpartition(".")
    if separator and base and suffix.isdecimal():
        return base
    return text


def build_gene_column_lookups(columns: Any) -> tuple[dict[str, Any], dict[str, Any]]:
    """Build exact and version-stripped lookups for input gene columns.

    The model accepts Ensembl IDs with or without ``.version`` suffixes, but
    duplicate or version-colliding columns are ambiguous enough to reject.
    """
    exact = {}
    duplicate_exact = []
    for column in columns:
        key = str(column)
        if key in exact:
            duplicate_exact.append(key)
        else:
            exact[key] = column
    if duplicate_exact:
        raise ValueError(
            "Duplicate gene columns are not allowed: " + _preview(sorted(set(duplicate_exact)))
        )

    stripped: dict[str, Any] = {}
    collisions: dict[str, list[str]] = {}
    for column in columns:
        base = strip_version(column)
        if base in stripped:
            collisions.setdefault(base, [str(stripped[base])]).append(str(column))
        else:
            stripped[base] = column
    if collisions:
        examples = []
        for base in sorted(collisions)[:5]:
            examples.append(f"{base} -> {', '.join(collisions[base])}")
        suffix = ", ..." if len(collisions) > 5 else ""
        raise ValueError(
            "Ambiguous gene columns after removing Ensembl version suffix: "
            + "; ".join(examples)
            + suffix
        )
    return exact, stripped


def _fraction(numerator: float, denominator: float) -> float:
    return float(numerator) / float(denominator) if denominator else 0.0


def align_to_genes_with_report(X: pd.DataFrame, genes: Any, impute_mean: np.ndarray | None = None) -> tuple[np.ndarray, dict[str, Any]]:
    """Reindex an expression DataFrame to the model's gene order.

    Matches Ensembl IDs with OR without the ``.version`` suffix (so a user CSV
    using unversioned IDs still aligns to a versioned model, and vice versa).
    Genes present in the input are coerced to numeric; genes missing from the
    input are filled with ``impute_mean`` (a neutral, standardized-zero value).

    X            samples x genes DataFrame
    genes        (g,) model gene order (may be versioned)
    impute_mean  (g,) per-gene training mean, or None to fill missing with NaN

    Returns ``(values, report)``. The report describes matched/missing genes and
    any matched cells that were non-numeric, NaN, or infinite before imputation.
    """
    genes = [str(g) for g in genes]
    exact, stripped = build_gene_column_lookups(X.columns)

    if impute_mean is None:
        means = np.full(len(genes), np.nan)
    else:
        means = np.asarray(impute_mean, dtype=float)
        if means.shape != (len(genes),):
            raise ValueError(
                f"impute_mean length/shape {means.shape} does not match gene count {len(genes)}"
            )
        if not np.all(np.isfinite(means)):
            raise ValueError("impute_mean must contain only finite values")

    out = np.empty((X.shape[0], len(genes)), dtype=float)
    missing = []
    invalid_cells = 0
    invalid_gene_reports = []
    all_invalid_genes = []
    sample_invalid_counts = np.zeros(X.shape[0], dtype=int)

    for j, g in enumerate(genes):
        src = exact.get(g)
        if src is None:
            src = stripped.get(strip_version(g))
        if src is None:
            out[:, j] = means[j]
            missing.append(g)
        else:
            values = pd.to_numeric(X[src], errors="coerce").to_numpy(dtype=float)
            finite = np.isfinite(values)
            invalid_mask = ~finite
            n_invalid = int(invalid_mask.sum())
            if n_invalid:
                invalid_cells += n_invalid
                sample_invalid_counts += invalid_mask.astype(int)
                invalid_gene_reports.append({
                    "gene": g,
                    "source_column": str(src),
                    "invalid_cells": n_invalid,
                    "total_cells": int(values.shape[0]),
                    "invalid_fraction": _fraction(n_invalid, values.shape[0]),
                })
                if values.shape[0] and int(finite.sum()) == 0:
                    all_invalid_genes.append(g)
            out[:, j] = np.where(finite, values, means[j])

    n_samples = int(X.shape[0])
    n_model_genes = len(genes)
    n_matched = n_model_genes - len(missing)
    matched_cells = n_samples * n_matched
    invalid_sample_indices = np.flatnonzero(sample_invalid_counts)
    all_invalid_sample_indices = (
        np.flatnonzero(sample_invalid_counts == n_matched) if n_matched else np.array([], dtype=int)
    )
    first_invalid_samples = []
    for i in invalid_sample_indices[:20]:
        n_invalid = int(sample_invalid_counts[i])
        first_invalid_samples.append({
            "sample": str(X.index[i]),
            "invalid_cells": n_invalid,
            "matched_genes": int(n_matched),
            "invalid_fraction": _fraction(n_invalid, n_matched),
        })
    max_sample_fraction = (
        _fraction(int(sample_invalid_counts.max()), n_matched)
        if n_samples and n_matched else 0.0
    )

    report = {
        "n_samples": n_samples,
        "n_model_genes": int(n_model_genes),
        "n_matched_genes": int(n_matched),
        "n_missing_genes": int(len(missing)),
        "missing_genes": [str(g) for g in missing],
        "matched_cells": int(matched_cells),
        "invalid_matched_cells": int(invalid_cells),
        "invalid_matched_fraction": _fraction(invalid_cells, matched_cells),
        "n_genes_with_invalid_values": int(len(invalid_gene_reports)),
        "n_genes_with_all_invalid_values": int(len(all_invalid_genes)),
        "first_genes_with_invalid_values": invalid_gene_reports[:20],
        "first_genes_with_all_invalid_values": [str(g) for g in all_invalid_genes[:20]],
        "n_samples_with_invalid_values": int(len(invalid_sample_indices)),
        "n_samples_with_all_invalid_values": int(len(all_invalid_sample_indices)),
        "max_invalid_matched_cell_fraction_per_sample": max_sample_fraction,
        "first_samples_with_invalid_values": first_invalid_samples,
        "first_samples_with_all_invalid_values": [
            str(X.index[i]) for i in all_invalid_sample_indices[:20]
        ],
    }
    return out, report


def validate_alignment_report(report: dict[str, Any], max_invalid_cell_fraction: float = 0.0) -> list[str]:
    """Return blocking issues for invalid values in matched model-gene cells."""
    max_invalid_cell_fraction = float(max_invalid_cell_fraction)
    issues: list[str] = []
    invalid_cells = int(report.get("invalid_matched_cells", 0))
    if invalid_cells <= 0:
        return issues

    all_invalid_genes = int(report.get("n_genes_with_all_invalid_values", 0))
    all_invalid_samples = int(report.get("n_samples_with_all_invalid_values", 0))
    invalid_fraction = float(report.get("invalid_matched_fraction", 0.0))
    max_sample_fraction = float(
        report.get("max_invalid_matched_cell_fraction_per_sample", 0.0)
    )
    if all_invalid_genes:
        examples = ", ".join(report.get("first_genes_with_all_invalid_values", [])[:5])
        suffix = f" Examples: {examples}." if examples else ""
        issues.append(
            f"{all_invalid_genes} matched model genes have no finite values.{suffix}"
        )
    if all_invalid_samples:
        examples = ", ".join(report.get("first_samples_with_all_invalid_values", [])[:5])
        suffix = f" Examples: {examples}." if examples else ""
        issues.append(
            f"{all_invalid_samples} samples have no finite matched model-gene values.{suffix}"
        )
    if invalid_fraction > max_invalid_cell_fraction:
        issues.append(
            "Invalid matched-value fraction "
            f"{invalid_fraction:.3%} exceeds --max-invalid-cell-fraction "
            f"{max_invalid_cell_fraction:.3%}."
        )
    if max_sample_fraction > max_invalid_cell_fraction:
        issues.append(
            "Worst-sample invalid matched-value fraction "
            f"{max_sample_fraction:.3%} exceeds --max-invalid-cell-fraction "
            f"{max_invalid_cell_fraction:.3%}."
        )
    return issues


def validate_gene_match_report(report: dict[str, Any], min_match_rate: float = 0.5) -> list[str]:
    """Return blocking issues for low model-gene coverage."""
    min_match_rate = float(min_match_rate)
    n_model_genes = int(report.get("n_model_genes", 0))
    n_matched_genes = int(report.get("n_matched_genes", 0))
    match_rate = _fraction(n_matched_genes, n_model_genes)
    if match_rate >= min_match_rate:
        return []

    if n_matched_genes == 0:
        return [
            "No model genes matched the input columns; check gene IDs and "
            "row/column orientation."
        ]
    return [
        f"Only {match_rate:.1%} of model genes matched "
        f"({n_matched_genes}/{n_model_genes}); required at least "
        f"{min_match_rate:.1%}. Check gene IDs and row/column orientation."
    ]


def format_gene_match_issues(report: dict[str, Any], min_match_rate: float = 0.5) -> str:
    """Return a single ValueError-ready message for low model-gene coverage."""
    issues = validate_gene_match_report(report, min_match_rate=min_match_rate)
    if not issues:
        return ""
    return "low model-gene coverage: " + " ".join(issues)


def format_alignment_issues(report: dict[str, Any], max_invalid_cell_fraction: float = 0.0) -> str:
    """Return a single ValueError-ready message for invalid matched values."""
    issues = validate_alignment_report(
        report,
        max_invalid_cell_fraction=max_invalid_cell_fraction,
    )
    if not issues:
        return ""
    return "invalid matched values: " + " ".join(issues)


def print_invalid_alignment_summary(report: dict[str, Any], stream: Any, prefix: str = "[score]") -> None:
    """Print a concise invalid matched-value summary to ``stream``."""
    invalid_cells = int(report.get("invalid_matched_cells", 0))
    if invalid_cells <= 0:
        return
    matched_cells = int(report.get("matched_cells", 0))
    invalid_fraction = float(report.get("invalid_matched_fraction", 0.0))
    print(
        f"{prefix} invalid matched values: "
        f"{invalid_cells}/{matched_cells} ({invalid_fraction:.3%}); "
        f"{report.get('n_genes_with_invalid_values', 0)} genes, "
        f"{report.get('n_samples_with_invalid_values', 0)} samples",
        file=stream,
    )
    gene_examples = report.get("first_genes_with_invalid_values", [])[:3]
    if gene_examples:
        text = ", ".join(
            f"{item['gene']}:{item['invalid_cells']}/{item['total_cells']}"
            for item in gene_examples
        )
        print(f"{prefix} invalid gene examples: {text}", file=stream)
    sample_examples = report.get("first_samples_with_invalid_values", [])[:3]
    if sample_examples:
        text = ", ".join(
            f"{item['sample']}:{item['invalid_cells']}/{item['matched_genes']}"
            for item in sample_examples
        )
        print(f"{prefix} invalid sample examples: {text}", file=stream)


def align_to_genes(
    X: pd.DataFrame,
    genes: Any,
    impute_mean: np.ndarray | None = None,
    max_invalid_cell_fraction: float = 0.0,
    allow_invalid_values: bool = False,
) -> tuple[np.ndarray, int, list[str]]:
    """Reindex an expression DataFrame to the model's gene order.

    Returns ``(values ndarray (n_samples, g), n_matched int, missing list)``.
    Missing model genes are imputed at ``impute_mean``. Invalid values in
    matched model-gene cells raise ``ValueError`` by default because this
    legacy return shape cannot carry the invalid-value report; pass
    ``allow_invalid_values=True`` only after reviewing mean imputation.
    """
    out, report = align_to_genes_with_report(X, genes, impute_mean=impute_mean)
    if not allow_invalid_values:
        message = format_alignment_issues(
            report,
            max_invalid_cell_fraction=max_invalid_cell_fraction,
        )
        if message:
            raise ValueError(message)
    missing = report["missing_genes"]
    return out, report["n_matched_genes"], missing
