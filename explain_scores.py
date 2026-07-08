#!/usr/bin/env python3
"""Explain LR tumor-vs-normal scores by per-gene logit contributions."""
import argparse
import os
import sys

import numpy as np
import pandas as pd

from score_tumor_normal import load_lr_weights
from tcga_rnaseq import align_to_genes, read_matrix, sigmoid


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


def explain_dataframe(df, weights, top_n, gene_names=None):
    genes = weights["selected_genes"]
    mean = weights["scaler_mean"]
    scale = weights["scaler_scale"]
    coef = weights["coef"]
    intercept = weights["intercept"]
    gene_names = gene_names or {}

    X, n_matched, missing = align_to_genes(df, genes, mean)
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
                    "gene_name": gene_names.get(gene, gene_names.get(gene.split(".")[0], "")),
                    "contribution_logit": float(contrib[j]),
                    "expression_log2_tpm1": float(X[i, j]),
                    "training_mean": float(mean[j]),
                    "scaled_value": float(X_scaled[i, j]),
                    "lr_coef": float(coef[j]),
                })
    return pd.DataFrame(rows, columns=EXPLANATION_COLUMNS), n_matched, missing


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="expression matrix (samples x genes)")
    parser.add_argument("-o", "--output", help="output CSV (default: <input>.explanations.csv)")
    parser.add_argument("--lr-weights", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "deployable_lr_weights.npz"))
    parser.add_argument("--gene-metadata", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "model_gene_metadata.csv"))
    parser.add_argument("--top-n", type=int, default=10)
    parser.add_argument("--transpose", action="store_true")
    args = parser.parse_args(argv)

    if args.top_n < 1:
        parser.error("--top-n must be >= 1")

    weights = load_lr_weights(args.lr_weights)
    df = read_matrix(args.input, transpose=args.transpose)
    gene_names = load_gene_metadata(args.gene_metadata)
    explanations, n_matched, missing = explain_dataframe(df, weights, args.top_n, gene_names)

    out = args.output or (os.path.splitext(args.input)[0] + ".explanations.csv")
    explanations.to_csv(out, index=False)
    print(f"[explain] {df.shape[0]} samples; matched {n_matched}/{len(weights['selected_genes'])} "
          f"model genes ({len(missing)} filled with training mean)", file=sys.stderr)
    print(f"[explain] wrote {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
