#!/usr/bin/env python3
"""Export the deployed logistic-regression scorer to a pure NumPy weight file."""
import argparse
import os
import sys
import tempfile
from pathlib import Path

import numpy as np

from tcga_rnaseq import (
    ensure_distinct_paths,
    load_lr_model,
    load_pipeline,
    validate_lr_model,
)


def _atomic_write_validated_npz(path, arrays):
    """Write a model archive through a validated sibling temporary file."""
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{destination.name}.", suffix=".tmp", dir=destination.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            np.savez_compressed(handle, **arrays)
            handle.flush()
            os.fsync(handle.fileno())
        load_lr_model(temporary)
        os.replace(temporary, destination)
    except BaseException:
        try:
            os.close(descriptor)
        except OSError:
            pass
        temporary.unlink(missing_ok=True)
        raise


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--pipeline", default="deployable_pipeline.pkl")
    parser.add_argument("--output", default="deployable_lr_weights.npz")
    parser.add_argument(
        "--trusted-pipeline",
        action="store_true",
        help="acknowledge that the pickle source was independently verified as trusted",
    )
    args = parser.parse_args(argv)
    if not args.trusted_pipeline:
        parser.error(
            "--trusted-pipeline is required because loading a pickle can execute code"
        )
    try:
        ensure_distinct_paths(
            {"weights output": args.output}, {"trusted pipeline input": args.pipeline}
        )
        pipe = load_pipeline(args.pipeline, trusted=True)
        lr = pipe["logistic_regression_model"]
        scaler = pipe["scaler"]

        if list(getattr(lr, "classes_", [])) != [0, 1]:
            raise ValueError(f"Unexpected LR classes: {getattr(lr, 'classes_', None)}")

        selected_genes = np.array(pipe["selected_genes"], dtype=str)
        scaler_mean = np.asarray(scaler.mean_, dtype=np.float64)
        scaler_scale = np.asarray(scaler.scale_, dtype=np.float64)
        coef_values = np.asarray(lr.coef_, dtype=np.float64)
        intercept_values = np.asarray(lr.intercept_, dtype=np.float64)
        if coef_values.ndim != 2 or coef_values.shape[0] != 1:
            raise ValueError(f"Unexpected LR coefficient shape: {coef_values.shape}")
        if intercept_values.shape != (1,):
            raise ValueError(f"Unexpected LR intercept shape: {intercept_values.shape}")
        coef = coef_values[0]
        intercept = intercept_values[0]

        canonical = validate_lr_model(
            {
                "genes": selected_genes,
                "mean": scaler_mean,
                "scale": scaler_scale,
                "coef": coef,
                "intercept": intercept,
                "classes": np.asarray(lr.classes_),
            }
        )
        arrays = {
            "selected_genes": canonical["genes"],
            "scaler_mean": canonical["mean"],
            "scaler_scale": canonical["scale"],
            "coef": canonical["coef"],
            "intercept": np.asarray(canonical["intercept"], dtype=np.float64),
            "class_order": np.asarray(canonical["classes"], dtype=np.int64),
            "source_pipeline": np.array(os.path.basename(args.pipeline)),
            "notes": np.array(
                "Pure NumPy logistic-regression export. Input must be GDC STAR-Counts "
                "log2(TPM+1), rows=samples, columns=Ensembl genes."
            ),
        }
        _atomic_write_validated_npz(args.output, arrays)
    except (ValueError, TypeError, KeyError, AttributeError, IndexError, OSError) as exc:
        parser.error(str(exc))
    print(f"[export] wrote {args.output} ({len(pipe['selected_genes'])} genes)",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
