#!/usr/bin/env python3
"""
predict_cancer_type.py -- score expression samples for TCGA cancer type
(tissue of origin) with the deployable multinomial logistic-regression model.

Pipeline (via the tcga_rnaseq shared core, pure numpy):
  z = (x_selected - scaler_mean) / scaler_scale
  logits = z @ coef.T + intercept        # coef: (17, k), one row per cancer type
  p = softmax(logits)                    # per-type probability
  call = argmax type

Input CSV: rows = samples, columns = Ensembl gene IDs, values = GDC STAR-Counts
log2(TPM+1). Ensembl version suffixes are matched with or without ".N"; missing
model genes are imputed at the training mean. Requires only numpy and pandas.
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))  # deliverables root, for tcga_rnaseq
from tcga_rnaseq import (  # noqa: E402
    load_lr_model,
    print_invalid_alignment_summary,
    read_matrix,
    validate_alignment_report,
    validate_gene_match_report,
    validate_threshold,
)
from tcga_rnaseq.align import align_to_genes_with_report  # noqa: E402
from tcga_rnaseq.score import predict_proba_from_aligned  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description="Predict TCGA cancer type (tissue of origin).")
    ap.add_argument("input_csv", help="rows=samples, cols=Ensembl gene IDs, values=log2(TPM+1)")
    ap.add_argument("--weights", default=os.path.join(HERE, "cancer_type_lr_weights.npz"))
    ap.add_argument("--out", default=None)
    ap.add_argument("--topk", type=int, default=3, help="report top-k types per sample")
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

    if not os.path.exists(args.weights):
        ap.error(f"weights file not found: {args.weights}")
    try:
        model = load_lr_model(args.weights)
    except ValueError as exc:
        ap.error(str(exc))
    if model["kind"] != "multiclass":
        ap.error("weights file is not a multi-class cancer-type model")
    if args.topk < 1:
        ap.error("--topk must be >= 1")
    if args.topk > len(model["classes"]):
        ap.error(f"--topk must be <= number of classes ({len(model['classes'])})")
    try:
        validate_threshold(args.max_invalid_cell_fraction, "--max-invalid-cell-fraction")
        validate_threshold(args.min_model_gene_match_rate, "--min-model-gene-match-rate")
    except ValueError as exc:
        ap.error(str(exc))
    try:
        X = read_matrix(args.input_csv)
    except ValueError as exc:
        ap.error(str(exc))
    values, alignment_report = align_to_genes_with_report(
        X, model["genes"], impute_mean=model["mean"]
    )
    print(
        f"[cancer-type] {X.shape[0]} samples; matched "
        f"{alignment_report['n_matched_genes']}/{alignment_report['n_model_genes']} "
        f"model genes ({alignment_report['n_missing_genes']} filled with training mean)",
        file=sys.stderr,
    )
    print_invalid_alignment_summary(alignment_report, sys.stderr, prefix="[cancer-type]")
    gene_match_issues = validate_gene_match_report(
        alignment_report,
        min_match_rate=args.min_model_gene_match_rate,
    )
    if gene_match_issues and not args.allow_low_gene_coverage:
        for issue in gene_match_issues:
            print(f"[cancer-type] ERROR: {issue}", file=sys.stderr)
        print(
            "[cancer-type] Refusing to write predictions with low model-gene coverage; "
            "fix the gene IDs/orientation or pass --allow-low-gene-coverage after "
            "reviewing the imputation.",
            file=sys.stderr,
        )
        return 1
    if gene_match_issues:
        for issue in gene_match_issues:
            print(f"[cancer-type] WARNING: {issue}", file=sys.stderr)
    alignment_issues = validate_alignment_report(
        alignment_report,
        max_invalid_cell_fraction=args.max_invalid_cell_fraction,
    )
    if alignment_issues and not args.allow_invalid_values:
        for issue in alignment_issues:
            print(f"[cancer-type] ERROR: {issue}", file=sys.stderr)
        print(
            "[cancer-type] Refusing to write predictions with invalid matched "
            "expression values; fix the input or pass --allow-invalid-values after "
            "reviewing the imputation.",
            file=sys.stderr,
        )
        return 1
    if alignment_issues:
        for issue in alignment_issues:
            print(f"[cancer-type] WARNING: {issue}", file=sys.stderr)
    P = predict_proba_from_aligned(model, values)
    classes = model["classes"]
    order = np.argsort(-P, axis=1)

    rows = []
    import pandas as pd
    for i, s in enumerate(X.index.astype(str)):
        top = order[i, :args.topk]
        row = {"sample": s, "predicted_type": classes[top[0]],
               "probability": round(float(P[i, top[0]]), 4)}
        for r, ci in enumerate(top, 1):
            row[f"top{r}"] = f"{classes[ci]}:{P[i, ci]:.3f}"
        rows.append(row)
    out = pd.DataFrame(rows)
    out_path = args.out or (os.path.splitext(args.input_csv)[0] + ".cancer_type_pred.csv")
    out.to_csv(out_path, index=False)
    print(json.dumps({"n_samples": int(X.shape[0]), "n_types": len(classes),
                      "matched_model_genes": int(alignment_report["n_matched_genes"]),
                      "missing_model_genes": int(alignment_report["n_missing_genes"]),
                      "invalid_matched_cells": int(alignment_report["invalid_matched_cells"]),
                      "invalid_matched_fraction": float(alignment_report["invalid_matched_fraction"]),
                      "scores_csv": out_path,
                      "call_distribution": out["predicted_type"].value_counts().to_dict()}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
