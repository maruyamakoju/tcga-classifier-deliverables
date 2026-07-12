#!/usr/bin/env python3
"""
cohort_adapt_score.py -- cross-platform (domain-adapted) scoring for the
TCGA/GDC tumor-vs-normal logistic-regression release.

The deployed model is calibrated to GDC STAR-Counts log2(TPM+1). On foreign
RNA-seq pipelines (e.g. UCSC Xena Toil/RSEM, GTEx/Toil) discrimination (AUC)
transfers but the 0.5 decision threshold does not: probabilities saturate and
almost every sample is called tumor.

This tool can apply an experimental, label-free, no-retraining transform: it
standardizes each gene using the *input cohort's own* per-gene mean and standard
deviation, then applies the frozen logistic-regression coefficients. Historical
benchmarks suggest this can reduce some platform shift, but it does not restore
calibration or guarantee that the default 0.5 threshold is valid.

Adaptation modes (from tcga_rnaseq.score.standardize):
  none            deployed scoring: z = (x - train_mean)/train_scale
  cohort_zscore   z = (x - cohort_mean)/cohort_std               (experimental)
  cohort_center   z = (x - cohort_mean)/train_scale              (location-only)

IMPORTANT LIMITATION: cohort standardization is transductive and composition
dependent: adding, removing, or regrouping samples changes every score. It
assumes an internal tumor/normal mix and is especially unsuitable for a
near-single-class cohort. Treat adapted values as experimental model scores,
not clinical risks or independently calibrated probabilities.

Requires only numpy and pandas (via the tcga_rnaseq shared core).
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tcga_rnaseq import (  # noqa: E402
    ensure_distinct_paths,
    load_lr_model,
    normalize_label,
    print_invalid_alignment_summary,
    read_csv_table,
    read_matrix,
    require_unique_samples,
    sample_key,
    score_binary_dataframe,
    validate_alignment_report,
    validate_gene_match_report,
    validate_threshold,
    write_dataframe_csv,
)
from tcga_rnaseq import metrics as M  # noqa: E402
from tcga_rnaseq.score import ADAPT_MODES  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def load_label_vector(labels_path, sample_index, sample_col="sample", label_col="label"):
    labels = read_csv_table(labels_path, string_columns=(sample_col,))
    if sample_col not in labels.columns:
        raise ValueError(f"labels CSV must contain {sample_col!r}")
    if label_col not in labels.columns:
        raise ValueError(f"labels CSV must contain {label_col!r}")

    labels = labels.copy()
    labels["_sample_key"] = require_unique_samples(labels, sample_col, "labels CSV")
    labels["label_binary"] = labels[label_col].map(normalize_label)
    sample_values = sample_key(sample_index, "expression matrix sample identifiers")
    if sample_values.duplicated().any():
        duplicates = ", ".join(sorted(sample_values[sample_values.duplicated()].unique())[:5])
        raise ValueError(f"expression matrix contains duplicate sample IDs: {duplicates}")
    sample_keys = pd.Index(sample_values.to_numpy(dtype=str), name="_sample_key")
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
    ap.add_argument("--adapt", default="none", choices=list(ADAPT_MODES),
                    help="domain-adaptation mode; explicit opt-in (default: none)")
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
    if args.min_samples < 2:
        ap.error("--min-samples must be >= 2")

    out_path = args.out or (os.path.splitext(args.input_csv)[0] + ".adapted_scores.csv")
    try:
        ensure_distinct_paths(
            {"adapted scores output": out_path},
            {
                "expression input": args.input_csv,
                "LR weights": args.weights,
                "labels input": args.labels,
            },
        )
        model = load_lr_model(args.weights)
        X = read_matrix(args.input_csv)
        n = X.shape[0]
        if args.adapt != "none" and n < args.min_samples:
            raise ValueError(
                f"adapted mode requires at least --min-samples={args.min_samples} "
                f"samples; received {n}. Use --adapt none or provide a larger cohort"
            )
        n_genes = len(model["genes"])
        out_df, n_matched, missing, alignment_report = score_binary_dataframe(
            model,
            X,
            threshold=args.threshold,
            adapt=args.adapt,
            allow_invalid_values=True,
            allow_low_gene_coverage=True,
            return_alignment_report=True,
        )
    except ValueError as exc:
        ap.error(str(exc))

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
    if args.adapt != "none":
        adaptation_warnings = [
            "cohort standardization assumes an internal tumor/normal mix; cohort "
            "composition can materially change every sample's model score",
            "adapted scoring is transductive: a sample's score depends on the other "
            "samples in the same batch, so scores from different batches are not "
            "directly comparable",
        ]
        warnings.extend(adaptation_warnings)
        for warning in adaptation_warnings:
            print(f"[adapt] WARNING: {warning}", file=sys.stderr)

    # The public score frame now retains full float precision, so this is the
    # exact probability used to derive `call`; do not recompute it a second time.
    p = out_df["tumor_probability"].to_numpy(dtype=float)
    metrics = None
    if args.labels:
        try:
            y, m, label_stats = load_label_vector(
                args.labels, X.index, args.sample_column, args.label_column
            )
        except ValueError as exc:
            ap.error(str(exc))
        if label_stats["n_labeled"] == 0:
            ap.error("labels CSV did not match any input samples")
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

    try:
        write_dataframe_csv(out_df, out_path, index=False)
    except ValueError as exc:
        ap.error(str(exc))

    report = {"n_samples": int(n), "adapt": args.adapt, "threshold": args.threshold,
              "score_interpretation": "model score, not calibrated clinical risk",
              "matched_model_genes": int(alignment_report["n_matched_genes"]),
              "missing_model_genes": int(alignment_report["n_missing_genes"]),
              "invalid_matched_cells": int(alignment_report["invalid_matched_cells"]),
              "invalid_matched_fraction": float(alignment_report["invalid_matched_fraction"]),
              "tumor_calls": int((out_df["call"] == "tumor").sum()),
              "normal_calls": int((out_df["call"] == "normal").sum()),
              "median_tumor_probability": float(np.median(p)),
              "scores_csv": out_path, "warnings": warnings, "metrics": metrics}
    print(json.dumps(report, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
