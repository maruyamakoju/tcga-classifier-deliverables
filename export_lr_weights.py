#!/usr/bin/env python3
"""Export the deployed logistic-regression scorer to a pure NumPy weight file."""
import argparse
import os
import sys

import numpy as np

from tcga_rnaseq import load_pipeline


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--pipeline", default="deployable_pipeline.pkl")
    parser.add_argument("--output", default="deployable_lr_weights.npz")
    args = parser.parse_args(argv)

    pipe = load_pipeline(args.pipeline)
    lr = pipe["logistic_regression_model"]
    scaler = pipe["scaler"]

    if list(getattr(lr, "classes_", [])) != [0, 1]:
        raise ValueError(f"Unexpected LR classes: {getattr(lr, 'classes_', None)}")

    selected_genes = np.array(pipe["selected_genes"], dtype=str)
    scaler_mean = np.asarray(scaler.mean_, dtype=np.float64)
    scaler_scale = np.asarray(scaler.scale_, dtype=np.float64)
    coef = np.asarray(lr.coef_[0], dtype=np.float64)
    intercept = np.asarray(lr.intercept_[0], dtype=np.float64)

    n_genes = len(selected_genes)
    if scaler_mean.shape != (n_genes,):
        raise ValueError(
            f"scaler.mean_ shape {scaler_mean.shape} does not match "
            f"selected_genes count {n_genes}"
        )
    if scaler_scale.shape != (n_genes,):
        raise ValueError(
            f"scaler.scale_ shape {scaler_scale.shape} does not match "
            f"selected_genes count {n_genes}"
        )
    if coef.shape != (n_genes,):
        raise ValueError(
            f"lr.coef_[0] shape {coef.shape} does not match selected_genes count {n_genes}"
        )

    np.savez_compressed(
        args.output,
        selected_genes=selected_genes,
        scaler_mean=scaler_mean,
        scaler_scale=scaler_scale,
        coef=coef,
        intercept=intercept,
        class_order=np.array(lr.classes_, dtype=np.int64),
        source_pipeline=os.path.basename(args.pipeline),
        notes=np.array(
            "Pure NumPy logistic-regression export. Input must be GDC STAR-Counts "
            "log2(TPM+1), rows=samples, columns=Ensembl genes."
        ),
    )
    print(f"[export] wrote {args.output} ({len(pipe['selected_genes'])} genes)",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
