#!/usr/bin/env python3
"""Inspect an expression matrix before tumor-vs-normal scoring.

This is a lightweight guardrail for the bundled GDC STAR-Counts-scale model. It
checks gene coverage, expression value range, standardized distribution shift,
and score saturation using only NumPy/pandas plus deployable_lr_weights.npz.
"""
import argparse
import json
import os
import sys

import numpy as np

from score_tumor_normal import load_lr_weights
from tcga_rnaseq import read_matrix, validate_threshold, write_json
from tcga_rnaseq.align import align_to_genes_with_report
from tcga_rnaseq.score import predict_proba_from_aligned, standardize


DEFAULT_RULES = {
    "min_match_rate_fail": 0.50,
    "min_match_rate_warn": 0.95,
    "negative_fraction_warn": 0.001,
    "value_p99_warn_above": 12.0,
    "value_p99_fail_above": 20.0,
    "value_max_warn_above": 30.0,
    "value_max_fail_above": 100.0,
    "nonfinite_fraction_warn": 0.0,
    "abs_z_gt_6_fraction_warn": 0.005,
    "abs_z_gt_10_fraction_warn": 0.001,
    "gene_mean_abs_z_p99_warn": 3.0,
    "gene_mean_abs_z_max_warn": 8.0,
    "expected_class_warn_fraction": 0.20,
}


def finite_or_none(value):
    value = float(value)
    return value if np.isfinite(value) else None


def quantiles(values, probs):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return {f"p{int(p * 100):02d}": None for p in probs}
    out = {}
    for p in probs:
        key = f"p{int(p * 100):02d}"
        if p == 0:
            key = "min"
        elif p == 1:
            key = "max"
        out[key] = finite_or_none(np.quantile(values, p))
    return out


def load_reference(path):
    if not path or not os.path.exists(path):
        return {"rules": DEFAULT_RULES, "reference_source": "built-in heuristic rules"}
    with open(path, "r", encoding="utf-8") as handle:
        ref = json.load(handle)
    rules = DEFAULT_RULES.copy()
    rules.update(ref.get("rules", {}))
    ref["rules"] = rules
    return ref


def add_message(messages, level, code, text):
    messages.append({"level": level, "code": code, "message": text})


def inspect_dataframe(df, weights, threshold=0.5, expected_class="unknown",
                      reference=None):
    selected_genes = weights["selected_genes"]
    mean = np.asarray(weights["scaler_mean"], dtype=float)
    scale = np.asarray(weights["scaler_scale"], dtype=float)
    coef = weights["coef"]
    intercept = weights["intercept"]
    model = {"genes": selected_genes, "mean": mean, "scale": scale,
             "coef": coef, "intercept": intercept, "kind": "binary"}
    reference = reference or {"rules": DEFAULT_RULES}
    rules = reference.get("rules", DEFAULT_RULES)
    messages = []

    X_raw, align_report = align_to_genes_with_report(df, selected_genes, impute_mean=None)
    n_model_genes = align_report["n_model_genes"]
    n_matched = align_report["n_matched_genes"]
    missing = align_report["missing_genes"]
    match_rate = n_matched / n_model_genes if n_model_genes else 0.0

    if n_matched == 0:
        add_message(messages, "ERROR", "no_model_genes_matched",
                    "No model genes matched the input columns.")
    elif match_rate < rules["min_match_rate_fail"]:
        add_message(messages, "ERROR", "very_low_model_gene_match",
                    f"Only {match_rate:.1%} of model genes matched.")
    elif match_rate < rules["min_match_rate_warn"]:
        add_message(messages, "WARNING", "low_model_gene_match",
                    f"{match_rate:.1%} of model genes matched; missing genes are filled with training means.")

    finite_mask = np.isfinite(X_raw)
    finite_values = X_raw[finite_mask]
    total_cells = X_raw.size
    nonfinite_fraction = 1.0 - (float(finite_mask.sum()) / total_cells if total_cells else 0.0)
    if nonfinite_fraction > rules["nonfinite_fraction_warn"]:
        add_message(messages, "WARNING", "nonfinite_or_missing_values",
                    f"{nonfinite_fraction:.3%} of model-gene cells are missing, non-numeric, or non-finite.")

    if finite_values.size:
        value_q = quantiles(finite_values, [0, 0.01, 0.50, 0.99, 1])
        negative_fraction = float(np.mean(finite_values < -1e-9))
        if negative_fraction > rules["negative_fraction_warn"]:
            add_message(messages, "WARNING", "negative_expression_values",
                        f"{negative_fraction:.3%} of finite values are negative; expected log2(TPM+1) is non-negative.")
        p99 = value_q["p99"]
        vmax = value_q["max"]
        if p99 is not None and p99 > rules["value_p99_fail_above"]:
            add_message(messages, "ERROR", "expression_values_too_large",
                        f"Expression p99={p99:.3g}, which is too high for expected log2(TPM+1) values.")
        elif p99 is not None and p99 > rules["value_p99_warn_above"]:
            add_message(messages, "WARNING", "expression_values_high",
                        f"Expression p99={p99:.3g}; check for raw TPM/counts or incompatible normalization.")
        if vmax is not None and vmax > rules["value_max_fail_above"]:
            add_message(messages, "ERROR", "expression_max_too_large",
                        f"Expression max={vmax:.3g}, which strongly suggests unlogged counts/TPM.")
        elif vmax is not None and vmax > rules["value_max_warn_above"]:
            add_message(messages, "WARNING", "expression_max_high",
                        f"Expression max={vmax:.3g}; check that values are log2(TPM+1).")
    else:
        value_q = {"min": None, "p01": None, "p50": None, "p99": None, "max": None}
        negative_fraction = None

    X = np.where(finite_mask, X_raw, mean)
    X_scaled = standardize(X, model, adapt="none")
    abs_z = np.abs(X_scaled)
    sample_median_abs_z = np.median(abs_z, axis=1) if df.shape[0] else np.array([])
    sample_p95_abs_z = np.percentile(abs_z, 95, axis=1) if df.shape[0] else np.array([])
    gene_mean_abs_z = np.abs(np.mean(X_scaled, axis=0)) if df.shape[0] else np.array([])
    abs_z_gt_6_fraction = float(np.mean(abs_z > 6)) if abs_z.size else None
    abs_z_gt_10_fraction = float(np.mean(abs_z > 10)) if abs_z.size else None

    if abs_z_gt_6_fraction is not None and abs_z_gt_6_fraction > rules["abs_z_gt_6_fraction_warn"]:
        add_message(messages, "WARNING", "many_extreme_standardized_values",
                    f"{abs_z_gt_6_fraction:.3%} of standardized model-gene values have |z| > 6.")
    if abs_z_gt_10_fraction is not None and abs_z_gt_10_fraction > rules["abs_z_gt_10_fraction_warn"]:
        add_message(messages, "WARNING", "very_extreme_standardized_values",
                    f"{abs_z_gt_10_fraction:.3%} of standardized model-gene values have |z| > 10.")

    gene_mean_q = quantiles(gene_mean_abs_z, [0.50, 0.90, 0.95, 0.99, 1])
    if df.shape[0] >= 20:
        p99_shift = gene_mean_q["p99"]
        max_shift = gene_mean_q["max"]
        if p99_shift is not None and p99_shift > rules["gene_mean_abs_z_p99_warn"]:
            add_message(messages, "WARNING", "cohort_distribution_shift",
                        f"Cohort gene-mean |z| p99={p99_shift:.3g}; this can indicate platform or tissue shift.")
        if max_shift is not None and max_shift > rules["gene_mean_abs_z_max_warn"]:
            add_message(messages, "WARNING", "large_cohort_gene_shift",
                        f"Maximum cohort gene-mean |z|={max_shift:.3g}; check normalization and gene annotation.")

    probabilities = predict_proba_from_aligned(model, X, adapt="none")
    calls = probabilities >= threshold
    tumor_fraction = float(np.mean(calls)) if calls.size else None
    normal_fraction = None if tumor_fraction is None else 1.0 - tumor_fraction
    expected_warn_fraction = rules["expected_class_warn_fraction"]
    if expected_class == "normal" and tumor_fraction is not None and tumor_fraction > expected_warn_fraction:
        add_message(messages, "WARNING", "unexpected_tumor_calls",
                    f"{tumor_fraction:.1%} of samples are called tumor in a cohort expected to be normal.")
    if expected_class == "tumor" and normal_fraction is not None and normal_fraction > expected_warn_fraction:
        add_message(messages, "WARNING", "unexpected_normal_calls",
                    f"{normal_fraction:.1%} of samples are called normal in a cohort expected to be tumor.")

    levels = {m["level"] for m in messages}
    status = "FAIL" if "ERROR" in levels else "WARN" if "WARNING" in levels else "PASS"

    return {
        "status": status,
        "threshold": threshold,
        "expected_class": expected_class,
        "reference_source": reference.get("reference_source", "built-in heuristic rules"),
        "shape": {
            "samples": int(df.shape[0]),
            "input_genes": int(df.shape[1]),
            "model_genes": int(n_model_genes),
        },
        "gene_match": {
            "matched_model_genes": int(n_matched),
            "missing_model_genes": int(len(missing)),
            "match_rate": finite_or_none(match_rate),
            "first_missing_model_genes": [str(g) for g in missing[:20]],
        },
        "value_summary": {
            "finite_model_gene_values": int(finite_values.size),
            "nonfinite_or_missing_fraction": finite_or_none(nonfinite_fraction),
            "negative_fraction": None if negative_fraction is None else finite_or_none(negative_fraction),
            **value_q,
        },
        "distribution_summary": {
            "sample_median_abs_z": quantiles(sample_median_abs_z, [0.50, 0.90, 0.95, 0.99, 1]),
            "sample_p95_abs_z": quantiles(sample_p95_abs_z, [0.50, 0.90, 0.95, 0.99, 1]),
            "cohort_gene_mean_abs_z": gene_mean_q,
            "abs_z_gt_6_fraction": None if abs_z_gt_6_fraction is None else finite_or_none(abs_z_gt_6_fraction),
            "abs_z_gt_10_fraction": None if abs_z_gt_10_fraction is None else finite_or_none(abs_z_gt_10_fraction),
        },
        "score_summary": {
            "tumor_calls": int(calls.sum()),
            "normal_calls": int((~calls).sum()),
            "tumor_call_fraction": None if tumor_fraction is None else finite_or_none(tumor_fraction),
            "tumor_probability": quantiles(probabilities, [0, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 1]),
            "probability_ge_0_99_fraction": finite_or_none(np.mean(probabilities >= 0.99)) if probabilities.size else None,
            "probability_le_0_01_fraction": finite_or_none(np.mean(probabilities <= 0.01)) if probabilities.size else None,
        },
        "messages": messages,
    }


def print_human_summary(report, output_path):
    gene = report["gene_match"]
    values = report["value_summary"]
    dist = report["distribution_summary"]
    score = report["score_summary"]
    print(
        f"[qc] status={report['status']}; matched "
        f"{gene['matched_model_genes']}/{report['shape']['model_genes']} model genes "
        f"({gene['match_rate']:.1%})",
        file=sys.stderr,
    )
    print(
        f"[qc] values: min={values['min']}, median={values['p50']}, "
        f"p99={values['p99']}, max={values['max']}",
        file=sys.stderr,
    )
    print(
        f"[qc] distribution: |z|>6 fraction={dist['abs_z_gt_6_fraction']}; "
        f"cohort gene-mean |z| p99={dist['cohort_gene_mean_abs_z']['p99']}",
        file=sys.stderr,
    )
    print(
        f"[qc] scores: {score['tumor_calls']} tumor / {score['normal_calls']} normal "
        f"at threshold {report['threshold']}; median probability="
        f"{score['tumor_probability']['p50']}",
        file=sys.stderr,
    )
    for msg in report["messages"]:
        print(f"[qc] {msg['level']}: {msg['message']}", file=sys.stderr)
    if output_path:
        print(f"[qc] wrote {output_path}", file=sys.stderr)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Inspect expression input compatibility before scoring.")
    parser.add_argument("input", help="expression matrix (samples x genes)")
    parser.add_argument("-o", "--output", help="output QC JSON (default: <input>.qc.json)")
    parser.add_argument("--lr-weights", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "deployable_lr_weights.npz"))
    parser.add_argument("--qc-reference", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "model_qc_reference.json"))
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--transpose", action="store_true")
    parser.add_argument("--expected-class", choices=["unknown", "mixed", "tumor", "normal"],
                        default="unknown",
                        help="optional cohort expectation used only for warning about surprising calls")
    parser.add_argument("--strict", action="store_true",
                        help="return non-zero on warnings as well as errors")
    args = parser.parse_args(argv)

    try:
        validate_threshold(args.threshold, "--threshold")
    except ValueError as exc:
        parser.error(str(exc))

    weights = load_lr_weights(args.lr_weights)
    reference = load_reference(args.qc_reference)
    try:
        df = read_matrix(args.input, transpose=args.transpose)
    except ValueError as exc:
        parser.error(str(exc))
    report = inspect_dataframe(df, weights, args.threshold, args.expected_class, reference)

    out = args.output or (os.path.splitext(args.input)[0] + ".qc.json")
    write_json(report, out, sort_keys=True)
    print_human_summary(report, out)

    if report["status"] == "FAIL" or (args.strict and report["status"] == "WARN"):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
