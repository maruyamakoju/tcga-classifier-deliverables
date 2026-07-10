#!/usr/bin/env python3
"""
cohort_adapt_score.py -- cross-platform (domain-adapted) scoring for the
TCGA/GDC tumor-vs-normal logistic-regression release.

The deployed model is calibrated to GDC STAR-Counts log2(TPM+1). On foreign
RNA-seq pipelines (e.g. UCSC Xena Toil/RSEM, GTEx/Toil) discrimination (AUC)
transfers but the 0.5 decision threshold does not: probabilities saturate and
almost every sample is called tumor.

This tool applies a label-free, no-retraining domain adaptation: it standardizes
each gene using the *input cohort's own* per-gene mean and standard deviation,
then applies the frozen logistic-regression coefficients. This realigns a
foreign cohort's per-gene marginal distribution onto the training marginal and
restores the default 0.5 threshold, provided the cohort contains an internal
tumor/normal contrast (see the important limitation below).

Adaptation modes (from tcga_rnaseq.score.standardize):
  none            deployed scoring: z = (x - train_mean)/train_scale
  cohort_zscore   z = (x - cohort_mean)/cohort_std               (recommended)
  cohort_center   z = (x - cohort_mean)/train_scale              (location-only)

IMPORTANT LIMITATION: cohort standardization assumes the input cohort has an
internal mix of tumor and normal. A near-single-class cohort (e.g. an all-normal
QC panel) has no internal contrast to anchor the recentering and is only
partially corrected -- prefer an explicit labeled-anchor recalibration
(calibrate_threshold.py) for such cohorts.

Requires only numpy and pandas (via the tcga_rnaseq shared core).
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from calibrate_threshold import normalize_label, validate_threshold  # noqa: E402
from tcga_rnaseq import (  # noqa: E402
    load_lr_model,
    predict_proba,
    print_invalid_alignment_summary,
    read_matrix,
    score_binary_dataframe,
    validate_alignment_report,
    validate_gene_match_report,
)
from tcga_rnaseq import metrics as M  # noqa: E402
from tcga_rnaseq.score import ADAPT_MODES  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def load_label_vector(labels_path, sample_index, sample_col="sample", label_col="label"):
    labels = pd.read_csv(labels_path)
    if sample_col not in labels.columns:
        raise ValueError(f"labels CSV must contain {sample_col!r}")
    if label_col not in labels.columns:
        raise ValueError(f"labels CSV must contain {label_col!r}")

    labels = labels.copy()
    sample_values = labels[sample_col]
    labels["_sample_key"] = sample_values.astype(str).str.strip()
    if sample_values.isna().any() or (labels["_sample_key"] == "").any():
        raise ValueError("labels CSV sample identifiers must be non-empty")
    duplicated = sorted(labels.loc[labels["_sample_key"].duplicated(), "_sample_key"].unique())
    if duplicated:
        raise ValueError("labels CSV contains duplicate sample IDs: " + ", ".join(duplicated[:5]))

    labels["label_binary"] = labels[label_col].map(normalize_label)
    sample_keys = pd.Index(sample_index.astype(str), name="_sample_key")
    label_series = labels.set_index("_sample_key")["label_binary"]
    aligned = label_series.reindex(sample_keys)
    matched = aligned.notna().to_numpy()
    extra_labels = int((~label_series.index.isin(sample_keys)).sum())
    stats = {
        "n_labels": int(len(labels)),
        "n_labeled": int(matched.sum()),
        "n_unmatched_samples": int((~matched).sum()),
        "n_extra_labels": extra_labels,
    }
    return aligned.to_numpy(dtype=float), matched, stats


def main(argv=None):
    ap = argparse.ArgumentParser(description="Cross-platform domain-adapted tumor-vs-normal scoring.")
    ap.add_argument("input_csv", help="expression matrix CSV: rows=samples, cols=Ensembl gene IDs, values=log2(TPM+1)")
    ap.add_argument("--adapt", default="cohort_zscore", choices=list(ADAPT_MODES),
                    help="domain-adaptation mode (default: cohort_zscore)")
    ap.add_argument("--threshold", type=float, default=0.5)
    ap.add_argument("--labels", default=None,
                    help="optional CSV with columns sample,label (label in {0,1} or {tumor,normal}) to report metrics")
    ap.add_argument("--sample-column", default="sample")
    ap.add_argument("--label-column", default="label")
    ap.add_argument("--weights", default=os.path.join(HERE, "deployable_lr_weights.npz"))
    ap.add_argument("--out", default=None)
    ap.add_argument("--min-samples", type=int, default=20)
    ap.add_argument("--max-invalid-cell-fraction", type=float, default=0.0,
                    help=("maximum allowed missing, non-numeric, NaN, or infinite cells "
                          "among matched model genes before failing (default 0)"))
    ap.add_argument("--allow-invalid-values", action="store_true",
                    help=("warn instead of failing when matched model-gene cells are "
                          "missing, non-numeric, NaN, or infinite"))
    ap.add_argument("--min-model-gene-match-rate", type=float, default=0.5,
                    help=("minimum fraction of model genes that must match input columns "
                          "before scoring (default 0.5)"))
    ap.add_argument("--allow-low-gene-coverage", action="store_true",
                    help=("warn instead of failing when too few model genes match; use "
                          "only after reviewing gene IDs and imputation"))
    args = ap.parse_args(argv)

    try:
        validate_threshold(args.threshold, "--threshold")
        validate_threshold(args.max_invalid_cell_fraction, "--max-invalid-cell-fraction")
        validate_threshold(args.min_model_gene_match_rate, "--min-model-gene-match-rate")
    except ValueError as exc:
        ap.error(str(exc))

    model = load_lr_model(args.weights)
    try:
        X = read_matrix(args.input_csv)
    except ValueError as exc:
        ap.error(str(exc))
    n = X.shape[0]
    n_genes = len(model["genes"])
    out_df, n_matched, missing, alignment_report = score_binary_dataframe(
        model,
        X,
        threshold=args.threshold,
        adapt=args.adapt,
        allow_invalid_values=True,
        return_alignment_report=True,
    )

    warnings = []
    print(
        f"[adapt] {n} samples; matched {n_matched}/{n_genes} model genes "
        f"({len(missing)} filled with training mean)",
        file=sys.stderr,
    )
    print_invalid_alignment_summary(alignment_report, sys.stderr, prefix="[adapt]")
    gene_match_issues = validate_gene_match_report(
        alignment_report,
        min_match_rate=args.min_model_gene_match_rate,
    )
    if gene_match_issues and not args.allow_low_gene_coverage:
        for issue in gene_match_issues:
            print(f"[adapt] ERROR: {issue}", file=sys.stderr)
        print(
            "[adapt] Refusing to write adapted scores with low model-gene coverage; "
            "fix the gene IDs/orientation or pass --allow-low-gene-coverage after "
            "reviewing the imputation.",
            file=sys.stderr,
        )
        return 1
    if gene_match_issues:
        for issue in gene_match_issues:
            warnings.append(issue)
            print(f"[adapt] WARNING: {issue}", file=sys.stderr)
    alignment_issues = validate_alignment_report(
        alignment_report,
        max_invalid_cell_fraction=args.max_invalid_cell_fraction,
    )
    if alignment_issues and not args.allow_invalid_values:
        for issue in alignment_issues:
            print(f"[adapt] ERROR: {issue}", file=sys.stderr)
        print(
            "[adapt] Refusing to write adapted scores with invalid matched expression "
            "values; fix the input or pass --allow-invalid-values after reviewing the "
            "imputation.",
            file=sys.stderr,
        )
        return 1
    if alignment_issues:
        for issue in alignment_issues:
            warnings.append(issue)
            print(f"[adapt] WARNING: {issue}", file=sys.stderr)
    if args.adapt != "none" and n < args.min_samples:
        warnings.append(f"cohort has only {n} samples (< {args.min_samples}); "
                        f"cohort statistics may be unreliable")
    if args.adapt != "none":
        warnings.append("cohort standardization assumes an internal tumor/normal mix; "
                        "a near-single-class cohort is only partially corrected")

    # Use the raw (unrounded) probability for metrics/summary stats, not the
    # rounded tumor_probability CSV column: score_binary_dataframe derives
    # `call` from the unrounded probability before rounding it for display,
    # so re-thresholding the rounded column here could disagree with `call`
    # for a sample whose raw probability rounds across the threshold.
    p = predict_proba(model, X, adapt=args.adapt, allow_invalid_values=True)

    out_path = args.out or (os.path.splitext(args.input_csv)[0] + ".adapted_scores.csv")
    out_df.to_csv(out_path, index=False)

    metrics = None
    if args.labels:
        try:
            y, m, label_stats = load_label_vector(
                args.labels, X.index, args.sample_column, args.label_column
            )
        except ValueError as exc:
            label_stats = None
            warnings.append(f"labels CSV could not be used: {exc}; metrics were not computed")
            print(f"[adapt] WARNING: labels CSV could not be used: {exc}", file=sys.stderr)
        if label_stats is not None and label_stats["n_labeled"] == 0:
            warnings.append("labels CSV did not match any input samples; metrics were not computed")
            label_stats = None
        if label_stats is not None:
            if label_stats["n_unmatched_samples"]:
                warnings.append(
                    f"labels matched {label_stats['n_labeled']}/{n} input samples; "
                    "metrics use matched samples only"
                )
            if label_stats["n_extra_labels"]:
                warnings.append(
                    f"labels CSV contains {label_stats['n_extra_labels']} rows not present in input"
                )
            if len(set(y[m])) > 1:
                cm = M.classification_metrics(y[m].astype(int), p[m], args.threshold)
                metrics = {**label_stats, "auc": round(cm["auc"], 4),
                           "accuracy": round(cm["accuracy"], 4),
                           "balanced_accuracy": round(cm["balanced_accuracy"], 4),
                           "sensitivity": round(cm["sensitivity"], 4),
                           "specificity": round(cm["specificity"], 4)}
            else:
                warnings.append("matched labels contain only one class; metrics were not computed")

    report = {"n_samples": int(n), "adapt": args.adapt, "threshold": args.threshold,
              "matched_model_genes": int(alignment_report["n_matched_genes"]),
              "missing_model_genes": int(alignment_report["n_missing_genes"]),
              "invalid_matched_cells": int(alignment_report["invalid_matched_cells"]),
              "invalid_matched_fraction": float(alignment_report["invalid_matched_fraction"]),
              "tumor_calls": int((out_df["call"] == "tumor").sum()),
              "normal_calls": int((out_df["call"] == "normal").sum()),
              "median_tumor_probability": float(np.median(p)),
              "scores_csv": out_path, "warnings": warnings, "metrics": metrics}
    print(__import__("json").dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
