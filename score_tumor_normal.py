#!/usr/bin/env python3
"""Score new RNA-seq samples as tumor vs normal with the TCGA pan-cancer classifier.

Uses `deployable_lr_weights.npz` by default for pure NumPy logistic-regression
scoring (2,000 genes; test AUC 0.997, leave-one-cancer-out AUC 0.994).

INPUT
  A table of expression values, rows = samples, columns = genes (Ensembl gene IDs,
  e.g. ENSG00000000005 or ENSG00000000005.6). Values must be **log2(TPM+1)** on the
  GDC STAR-Counts scale, the same as training. Accepted: .csv .tsv .txt .parquet.
  Use --transpose if your genes are rows. Very low model-gene coverage is refused
  by default before scores are written.

USAGE
  python score_tumor_normal.py expr.csv                     # -> expr.scored.csv
  python score_tumor_normal.py expr.csv -o out.csv --threshold 0.5
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
from tcga_rnaseq import (  # noqa: E402
    load_lr_model,
    print_invalid_alignment_summary,
    read_matrix,
    score_binary_dataframe,
    validate_alignment_report,
    validate_gene_match_report,
    validate_threshold,
)

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


def score_dataframe_lr_weights(
    df,
    model,
    threshold=0.5,
    max_invalid_cell_fraction=0.0,
    allow_invalid_values=False,
    return_alignment_report=False,
):
    """Score with the pure-NumPy npz logistic-regression model.
    `model` may be a tcga_rnaseq model dict or a legacy load_lr_weights dict."""
    return score_binary_dataframe(
        _as_model(model),
        df,
        threshold=threshold,
        max_invalid_cell_fraction=max_invalid_cell_fraction,
        allow_invalid_values=allow_invalid_values,
        return_alignment_report=return_alignment_report,
    )


def run_self_test(lr_weights_path):
    example_in = os.path.join(HERE, "example_input.csv")
    expected_out = os.path.join(HERE, "example_output.csv")

    model = load_lr_model(lr_weights_path)
    observed, n_matched, missing = score_dataframe_lr_weights(read_matrix(example_in), model)
    n_genes = len(model["genes"])
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
    ap.add_argument("--lr-weights", default=os.path.join(HERE, "deployable_lr_weights.npz"),
                    help="path to pure NumPy LR weights (default: deployable_lr_weights.npz)")
    ap.add_argument("--model", choices=["lr"], default="lr",
                    help="logistic regression scorer (default; only public-bundle model)")
    ap.add_argument("--threshold", type=float, default=0.5,
                    help="probability cutoff for the tumor call (default 0.5)")
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
    ap.add_argument("--transpose", action="store_true",
                    help="input has genes as rows, samples as columns")
    ap.add_argument("--self-test", action="store_true",
                    help="score bundled example_input.csv and compare with example_output.csv")
    ap.add_argument("-p", "--pipeline", help=argparse.SUPPRESS)
    ap.add_argument("--use-pickle-lr", action="store_true", help=argparse.SUPPRESS)
    args = ap.parse_args(argv)

    try:
        validate_threshold(args.threshold, "--threshold")
        validate_threshold(args.max_invalid_cell_fraction, "--max-invalid-cell-fraction")
        validate_threshold(args.min_model_gene_match_rate, "--min-model-gene-match-rate")
    except ValueError as exc:
        ap.error(str(exc))
    if args.pipeline or args.use_pickle_lr:
        ap.error(
            "legacy pickle/RF scoring is not available in the public lightweight release; "
            "use deployable_lr_weights.npz with the default NumPy LR scorer"
        )
    if args.self_test:
        return run_self_test(args.lr_weights)
    if not args.input:
        ap.error("input is required unless --self-test is used")

    try:
        df = read_matrix(args.input, transpose=args.transpose)
    except ValueError as exc:
        ap.error(str(exc))
    if not os.path.exists(args.lr_weights):
        ap.error(f"LR weights file not found: {args.lr_weights}")
    model = load_lr_model(args.lr_weights)
    res, n_matched, missing, alignment_report = score_dataframe_lr_weights(
        df,
        model,
        args.threshold,
        allow_invalid_values=True,
        return_alignment_report=True,
    )
    n_genes = len(model["genes"])
    # --model is kept for CLI compatibility (TROUBLESHOOTING.md documents that
    # --model rf is rejected); choices=["lr"] means this can only ever be "lr".
    scorer = f"{args.model}-numpy"

    print(f"[score] scorer={scorer}; {df.shape[0]} samples; matched {n_matched}/{n_genes} "
          f"model genes ({len(missing)} filled with training mean)", file=sys.stderr)
    print_invalid_alignment_summary(alignment_report, sys.stderr)
    gene_match_issues = validate_gene_match_report(
        alignment_report,
        min_match_rate=args.min_model_gene_match_rate,
    )
    if gene_match_issues and not args.allow_low_gene_coverage:
        for issue in gene_match_issues:
            print(f"[score] ERROR: {issue}", file=sys.stderr)
        print(
            "[score] Refusing to write scores with low model-gene coverage; fix the "
            "gene IDs/orientation or pass --allow-low-gene-coverage after reviewing "
            "the imputation.",
            file=sys.stderr,
        )
        return 1
    if gene_match_issues:
        for issue in gene_match_issues:
            print(f"[score] WARNING: {issue}", file=sys.stderr)
    alignment_issues = validate_alignment_report(
        alignment_report,
        max_invalid_cell_fraction=args.max_invalid_cell_fraction,
    )
    if alignment_issues and not args.allow_invalid_values:
        for issue in alignment_issues:
            print(f"[score] ERROR: {issue}", file=sys.stderr)
        print(
            "[score] Refusing to write scores with invalid matched expression values; "
            "fix the input or pass --allow-invalid-values after reviewing the imputation.",
            file=sys.stderr,
        )
        return 1
    if alignment_issues:
        for issue in alignment_issues:
            print(f"[score] WARNING: {issue}", file=sys.stderr)

    out = args.output or (os.path.splitext(args.input)[0] + ".scored.csv")
    res.to_csv(out, index=False)
    call = res["call"].to_numpy()
    print(f"[score] wrote {out}  ({(call=='tumor').sum()} tumor / "
          f"{(call=='normal').sum()} normal at threshold {args.threshold})", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
