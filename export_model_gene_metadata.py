#!/usr/bin/env python3
"""Create a transparent gene metadata table for the lightweight LR scorer."""
import argparse
import os

import numpy as np
import pandas as pd


def load_known_gene_names(paths):
    names = {}
    for path in paths:
        if not os.path.exists(path):
            continue
        df = pd.read_csv(path)
        if "gene_id" not in df.columns or "gene_name" not in df.columns:
            continue
        for gene_id, gene_name in zip(df["gene_id"], df["gene_name"]):
            if pd.notna(gene_name):
                names[str(gene_id)] = str(gene_name)
                names[str(gene_id).split(".")[0]] = str(gene_name)
    return names


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--weights", default="deployable_lr_weights.npz")
    parser.add_argument("--output", default="model_gene_metadata.csv")
    parser.add_argument("--name-source", action="append",
                        default=["top_genes_logreg.csv", "top_genes_xgboost.csv"],
                        help="CSV with gene_id,gene_name columns; may be repeated")
    args = parser.parse_args(argv)

    weights = np.load(args.weights, allow_pickle=False)
    genes = weights["selected_genes"].astype(str)
    coef = weights["coef"].astype(float)
    names = load_known_gene_names(args.name_source)

    df = pd.DataFrame({
        "gene_id": genes,
        "gene_id_base": [gene.split(".")[0] for gene in genes],
        "gene_name": [names.get(gene, names.get(gene.split(".")[0], "")) for gene in genes],
        "lr_coef": coef,
        "abs_lr_coef": np.abs(coef),
        "direction_if_high": np.where(coef >= 0, "tumor", "normal"),
        "scaler_mean": weights["scaler_mean"].astype(float),
        "scaler_scale": weights["scaler_scale"].astype(float),
    })
    df = df.sort_values("abs_lr_coef", ascending=False).reset_index(drop=True)
    df.insert(0, "rank_abs_lr_coef", np.arange(1, len(df) + 1))
    df.to_csv(args.output, index=False)
    print(f"[metadata] wrote {args.output} ({len(df)} genes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
