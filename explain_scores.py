#!/usr/bin/env python3
"""Explain LR tumor-vs-normal scores by per-gene logit contributions."""
import argparse
import os
import sys

import numpy as np
import pandas as pd

from calibrate_threshold import validate_threshold
from score_tumor_normal import load_lr_weights
from tcga_rnaseq import (
    align_to_genes_with_report,
    print_invalid_alignment_summary,
    read_matrix,
    sigmoid,
    strip_version,
    validate_alignment_report,
    validate_gene_match_report,
)


EXPLANATION_COLUMNS = [
    "sample", "tumor_probability", "logit", "direction", "rank", "gene_id",
    "gene_name", "contribution_logit", "expression_log2_tpm1", "training_mean",
    "scaled_value", "lr_coef",
]


def load_gene_metadata(path):
    if not path or not os.path.exists(path):
        return {}
    df = pd.read_csv(path)
    out = {}
    for row in df.itertuples(index=False):
        name = getattr(row, "gene_name", "")
        out[str(row.gene_id)] = name if isinstance(name, str) else ""
        out[str(row.gene_id_base)] = out[str(row.gene_id)]
    return out


def explain_dataframe(df, weights, top_n, gene_names=None, return_alignment_report=False):
    genes = weights["selected_genes"]
    mean = weights["scaler_mean"]
    scale = weights["scaler_scale"]
    coef = weights["coef"]
    intercept = weights["intercept"]
    gene_names = gene_names or {}

    X, alignment_report = align_to_genes_with_report(df, genes, mean)
    n_matched = alignment_report["n_matched_genes"]
    missing = alignment_report["missing_genes"]
    X_scaled = (X - mean) / scale
    contributions = X_scaled * coef
    logits = contributions.sum(axis=1) + intercept
    probabilities = sigmoid(logits)

    rows = []
    for i, sample in enumerate(df.index):
        contrib = contributions[i]
        top_tumor = np.argsort(contrib)[::-1][:top_n]
        top_normal = np.argsort(contrib)[:top_n]
        for direction, indices in [("tumor", top_tumor), ("normal", top_normal)]:
            for rank, j in enumerate(indices, start=1):
                gene = genes[j]
                rows.append({
                    "sample": sample,
                    "tumor_probability": round(float(probabilities[i]), 6),
                    "logit": float(logits[i]),
                    "direction": direction,
                    "rank": rank,
                    "gene_id": gene,
                    "gene_name": gene_names.get(gene, gene_names.get(strip_version(gene), "")),
                    "contribution_logit": float(contrib[j]),
                    "expression_log2_tpm1": float(X[i, j]),
                    "training_mean": float(mean[j]),
                    "scaled_value": float(X_scaled[i, j]),
                    "lr_coef": float(coef[j]),
                })
    result = pd.DataFrame(rows, columns=EXPLANATION_COLUMNS)
    if return_alignment_report:
        return result, n_matched, missing, alignment_report
    return result, n_matched, missing


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="expression matrix (samples x genes)")
    parser.add_argument("-o", "--output", help="output CSV (default: <input>.explanations.csv)")
    parser.add_argument("--lr-weights", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "deployable_lr_weights.npz"))
    parser.add_argument("--gene-metadata", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "model_gene_metadata.csv"))
    parser.add_argument("--top-n", type=int, default=10)
    parser.add_argument("--max-invalid-cell-fraction", type=float, default=0.0,
                        help=("maximum allowed missing, non-numeric, NaN, or infinite cells "
                              "among matched model genes before failing (default 0)"))
    parser.add_argument("--allow-invalid-values", action="store_true",
                        help=("warn instead of failing when matched model-gene cells are "
                              "missing, non-numeric, NaN, or infinite"))
    parser.add_argument("--min-model-gene-match-rate", type=float, default=0.5,
                        help=("minimum fraction of model genes that must match input columns "
                              "before writing explanations (default 0.5)"))
    parser.add_argument("--allow-low-gene-coverage", action="store_true",
                        help=("warn instead of failing when too few model genes match; use "
                              "only after reviewing gene IDs and imputation"))
    parser.add_argument("--transpose", action="store_true")
    args = parser.parse_args(argv)

    if args.top_n < 1:
        parser.error("--top-n must be >= 1")
    try:
        validate_threshold(args.max_invalid_cell_fraction, "--max-invalid-cell-fraction")
        validate_threshold(args.min_model_gene_match_rate, "--min-model-gene-match-rate")
    except ValueError as exc:
        parser.error(str(exc))

    weights = load_lr_weights(args.lr_weights)
    try:
        df = read_matrix(args.input, transpose=args.transpose)
    except ValueError as exc:
        parser.error(str(exc))
    gene_names = load_gene_metadata(args.gene_metadata)
    explanations, n_matched, missing, alignment_report = explain_dataframe(
        df,
        weights,
        args.top_n,
        gene_names,
        return_alignment_report=True,
    )

    print(f"[explain] {df.shape[0]} samples; matched {n_matched}/{len(weights['selected_genes'])} "
          f"model genes ({len(missing)} filled with training mean)", file=sys.stderr)
    print_invalid_alignment_summary(alignment_report, sys.stderr, prefix="[explain]")
    gene_match_issues = validate_gene_match_report(
        alignment_report,
        min_match_rate=args.min_model_gene_match_rate,
    )
    if gene_match_issues and not args.allow_low_gene_coverage:
        for issue in gene_match_issues:
            print(f"[explain] ERROR: {issue}", file=sys.stderr)
        print(
            "[explain] Refusing to write explanations with low model-gene coverage; "
            "fix the gene IDs/orientation or pass --allow-low-gene-coverage after "
            "reviewing the imputation.",
            file=sys.stderr,
        )
        return 1
    if gene_match_issues:
        for issue in gene_match_issues:
            print(f"[explain] WARNING: {issue}", file=sys.stderr)
    alignment_issues = validate_alignment_report(
        alignment_report,
        max_invalid_cell_fraction=args.max_invalid_cell_fraction,
    )
    if alignment_issues and not args.allow_invalid_values:
        for issue in alignment_issues:
            print(f"[explain] ERROR: {issue}", file=sys.stderr)
        print(
            "[explain] Refusing to write explanations with invalid matched expression values; "
            "fix the input or pass --allow-invalid-values after reviewing the imputation.",
            file=sys.stderr,
        )
        return 1
    if alignment_issues:
        for issue in alignment_issues:
            print(f"[explain] WARNING: {issue}", file=sys.stderr)

    out = args.output or (os.path.splitext(args.input)[0] + ".explanations.csv")
    explanations.to_csv(out, index=False)
    print(f"[explain] wrote {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
