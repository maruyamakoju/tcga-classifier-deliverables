#!/usr/bin/env python3
"""Score new RNA-seq samples as tumor vs normal with the TCGA pan-cancer classifier.

Uses `deployable_lr_weights.npz` by default for pure NumPy logistic-regression
scoring (2,000 genes; test AUC 0.997, leave-one-cancer-out AUC 0.994). The
legacy `deployable_pipeline.pkl` path remains available for RF scoring.

INPUT
  A table of expression values, rows = samples, columns = genes (Ensembl gene IDs,
  e.g. ENSG00000000005 or ENSG00000000005.6). Values must be **log2(TPM+1)** on the
  GDC STAR-Counts scale, the same as training. Accepted: .csv .tsv .parquet .pkl
  (a pickled pandas DataFrame). Use --transpose if your genes are rows.

USAGE
  python score_tumor_normal.py expr.csv                     # -> expr.scored.csv
  python score_tumor_normal.py expr.csv -o out.csv --threshold 0.5
  python score_tumor_normal.py expr.pkl --model rf          # random forest instead
  python score_tumor_normal.py expr.csv --use-pickle-lr     # legacy sklearn LR path
  python score_tumor_normal.py --self-test                  # verify bundled example

THRESHOLD NOTE
  Ranking (AUC) transfers across cancer types, but the 0.5 decision threshold does
  not always transfer to a tissue the model never trained on (prostate/liver were
  under-called at 0.5). If you are scoring a genuinely new tissue and have a few
  labeled samples, pick a threshold on them and pass it with --threshold. For a
  cohort from a different RNA-seq pipeline, see cohort_adapt_score.py.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tcga_rnaseq import (load_lr_model, load_pipeline, read_matrix,  # noqa: E402
                         align_to_genes, score_binary_dataframe)

HERE = os.path.dirname(os.path.abspath(__file__))


def load_lr_weights(path):
    """Legacy dict shape (selected_genes/scaler_mean/scaler_scale/coef/intercept)
    for callers predating tcga_rnaseq.load_lr_model. Prefer load_lr_model."""
    m = load_lr_model(path)
    if m["kind"] != "binary":
        raise ValueError("tumor-vs-normal scoring requires a binary LR weights file")
    return {"selected_genes": m["genes"].astype(str).tolist(),
            "scaler_mean": m["mean"], "scaler_scale": m["scale"],
            "coef": m["coef"], "intercept": m["intercept"]}


def _as_model(model):
    """Accept either a tcga_rnaseq model dict or the legacy load_lr_weights dict."""
    if "genes" in model:
        if model.get("kind") != "binary":
            raise ValueError("score_dataframe_lr_weights requires a binary model")
        return model
    coef = np.asarray(model["coef"], dtype=float)
    if coef.ndim != 1:
        raise ValueError("score_dataframe_lr_weights requires a 1-D binary coefficient vector")
    return {"genes": np.asarray(model["selected_genes"], dtype=str),
            "mean": np.asarray(model["scaler_mean"], dtype=float),
            "scale": np.asarray(model["scaler_scale"], dtype=float),
            "coef": coef,
            "intercept": float(model["intercept"]),
            "classes": np.array([0, 1]),
            "kind": "binary"}


def score_dataframe_lr_weights(df, model, threshold=0.5):
    """Score with the pure-NumPy npz logistic-regression model.
    `model` may be a tcga_rnaseq model dict or a legacy load_lr_weights dict."""
    return score_binary_dataframe(_as_model(model), df, threshold=threshold)


def score_dataframe(df, pipe, model_name="lr", threshold=0.5):
    """Score with the legacy pickled sklearn pipeline (LR scaled, RF unscaled)."""
    selected_genes = pipe["selected_genes"]
    scaler = pipe["scaler"]
    model = (pipe["logistic_regression_model"]
             if model_name == "lr" else pipe["random_forest_model"])
    X, n_matched, missing = align_to_genes(df, selected_genes, impute_mean=scaler.mean_)
    if model_name == "lr":
        proba = model.predict_proba(scaler.transform(X))[:, 1]
    else:
        proba = model.predict_proba(X)[:, 1]  # RF trained on unscaled selected features
    call = np.where(proba >= threshold, "tumor", "normal")
    res = pd.DataFrame({"sample": df.index, "tumor_probability": proba.round(6), "call": call})
    return res, n_matched, missing


def run_self_test(pipeline_path, lr_weights_path, use_pickle_lr=False):
    example_in = os.path.join(HERE, "example_input.csv")
    expected_out = os.path.join(HERE, "example_output.csv")

    if not use_pickle_lr and os.path.exists(lr_weights_path):
        model = load_lr_model(lr_weights_path)
        observed, n_matched, missing = score_dataframe_lr_weights(read_matrix(example_in), model)
        n_genes = len(model["genes"])
    else:
        pipe = load_pipeline(pipeline_path)
        observed, n_matched, missing = score_dataframe(read_matrix(example_in), pipe)
        n_genes = len(pipe["selected_genes"])
    expected = pd.read_csv(expected_out)

    same_samples = observed["sample"].tolist() == expected["sample"].tolist()
    same_calls = observed["call"].tolist() == expected["call"].tolist()
    same_length = len(observed) == len(expected)
    if same_length:
        max_delta = float(np.max(np.abs(observed["tumor_probability"].to_numpy()
                                        - expected["tumor_probability"].to_numpy())))
    else:
        max_delta = float("inf")
    ok = same_length and same_samples and same_calls and max_delta <= 1e-6 and not missing

    print(f"[self-test] matched {n_matched}/{n_genes} model genes", file=sys.stderr)
    print(f"[self-test] max probability delta vs expected: {max_delta:.6g}", file=sys.stderr)
    if ok:
        print("[self-test] PASS: bundled example reproduces expected calls", file=sys.stderr)
        return 0
    print("[self-test] FAIL: bundled example does not match expected output", file=sys.stderr)
    if not same_samples:
        print("[self-test] sample order differs", file=sys.stderr)
    if not same_calls:
        print("[self-test] calls differ", file=sys.stderr)
    if not same_length:
        print(
            f"[self-test] row count differs: observed={len(observed)} expected={len(expected)}",
            file=sys.stderr,
        )
    if missing:
        print(f"[self-test] {len(missing)} model genes missing from example input", file=sys.stderr)
    return 1


def main(argv=None):
    ap = argparse.ArgumentParser(description="Score samples as tumor vs normal.")
    ap.add_argument("input", nargs="?", help="expression matrix (samples x genes)")
    ap.add_argument("-o", "--output", help="output CSV (default: <input>.scored.csv)")
    ap.add_argument("-p", "--pipeline", default=os.path.join(HERE, "deployable_pipeline.pkl"),
                    help="path to deployable_pipeline.pkl (needed for --model rf or --use-pickle-lr)")
    ap.add_argument("--lr-weights", default=os.path.join(HERE, "deployable_lr_weights.npz"),
                    help="path to pure NumPy LR weights (default: deployable_lr_weights.npz)")
    ap.add_argument("--model", choices=["lr", "rf"], default="lr",
                    help="lr = logistic regression (default, best); rf = random forest")
    ap.add_argument("--threshold", type=float, default=0.5,
                    help="probability cutoff for the tumor call (default 0.5)")
    ap.add_argument("--transpose", action="store_true",
                    help="input has genes as rows, samples as columns")
    ap.add_argument("--self-test", action="store_true",
                    help="score bundled example_input.csv and compare with example_output.csv")
    ap.add_argument("--use-pickle-lr", action="store_true",
                    help="use legacy sklearn LR object from deployable_pipeline.pkl")
    args = ap.parse_args(argv)

    if not 0 <= args.threshold <= 1:
        ap.error("--threshold must be between 0 and 1")
    if args.self_test:
        return run_self_test(args.pipeline, args.lr_weights, args.use_pickle_lr)
    if not args.input:
        ap.error("input is required unless --self-test is used")

    df = read_matrix(args.input, transpose=args.transpose)
    if args.model == "lr" and not args.use_pickle_lr and os.path.exists(args.lr_weights):
        model = load_lr_model(args.lr_weights)
        res, n_matched, missing = score_dataframe_lr_weights(df, model, args.threshold)
        n_genes = len(model["genes"])
        scorer = "lr-numpy"
    else:
        if args.model == "lr" and not args.use_pickle_lr:
            print(f"[score] WARNING: {args.lr_weights} not found; falling back to pickle LR",
                  file=sys.stderr)
        pipe = load_pipeline(args.pipeline)
        res, n_matched, missing = score_dataframe(df, pipe, args.model, args.threshold)
        n_genes = len(pipe["selected_genes"])
        scorer = args.model

    print(f"[score] scorer={scorer}; {df.shape[0]} samples; matched {n_matched}/{n_genes} "
          f"model genes ({len(missing)} filled with training mean)", file=sys.stderr)
    if n_matched < 0.5 * n_genes:
        print("[score] WARNING: <50% of model genes found - check gene IDs and that "
              "values are log2(TPM+1). Results may be unreliable.", file=sys.stderr)

    out = args.output or (os.path.splitext(args.input)[0] + ".scored.csv")
    res.to_csv(out, index=False)
    call = res["call"].to_numpy()
    print(f"[score] wrote {out}  ({(call=='tumor').sum()} tumor / "
          f"{(call=='normal').sum()} normal at threshold {args.threshold})", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
