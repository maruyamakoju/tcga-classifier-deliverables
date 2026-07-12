#!/usr/bin/env python3
"""Choose a tumor-probability threshold from labeled scored samples."""
import argparse
import os
import sys

import numpy as np
import pandas as pd

from tcga_rnaseq import metrics as M
from tcga_rnaseq import (
    ensure_distinct_paths,
    normalize_label,
    read_csv_table,
    require_unique_samples,
    validate_threshold,
    write_dataframe_csv,
    write_json,
)


def _validated_probabilities(series):
    values = pd.to_numeric(series, errors="coerce")
    arr = values.to_numpy(dtype=float)
    if not np.all(np.isfinite(arr)):
        raise ValueError("tumor_probability must contain only finite numeric values")
    if not np.all((arr >= 0) & (arr <= 1)):
        raise ValueError("tumor_probability values must be between 0 and 1")
    return values


def load_scores_and_labels(scores_path, labels_path, sample_col, label_col,
                           min_match_fraction=1.0):
    min_match_fraction = float(min_match_fraction)
    if not np.isfinite(min_match_fraction) or not 0 < min_match_fraction <= 1:
        raise ValueError("min_match_fraction must be finite and in (0, 1]")
    if sample_col == label_col:
        raise ValueError("sample and label columns must be different")
    if label_col == "tumor_probability":
        raise ValueError("label column must not be tumor_probability")
    scores = read_csv_table(scores_path, string_columns=(sample_col,))
    if "tumor_probability" not in scores.columns:
        raise ValueError("scores CSV must contain tumor_probability")
    if sample_col not in scores.columns:
        raise ValueError(f"scores CSV must contain {sample_col!r}")
    scores = scores.copy()
    scores["_sample_key"] = require_unique_samples(scores, sample_col, "scores CSV")
    scores["tumor_probability"] = _validated_probabilities(scores["tumor_probability"])

    if labels_path:
        labels = read_csv_table(labels_path, string_columns=(sample_col,))
        if sample_col not in labels.columns:
            raise ValueError(f"labels CSV must contain {sample_col!r}")
        if label_col not in labels.columns:
            raise ValueError(f"labels CSV must contain {label_col!r}")
        labels = labels.copy()
        labels["_sample_key"] = require_unique_samples(labels, sample_col, "labels CSV")
        label_value_column = "_calibration_label_value"
        merged = scores[[sample_col, "_sample_key", "tumor_probability"]].merge(
            labels[["_sample_key", label_col]].rename(
                columns={label_col: label_value_column}
            ),
            on="_sample_key", how="inner",
            validate="one_to_one"
        )
        match_fraction = len(merged) / len(scores) if len(scores) else 0.0
        if len(merged) != len(scores):
            message = f"matched {len(merged)}/{len(scores)} scored samples"
            if match_fraction < min_match_fraction:
                raise ValueError(
                    f"{message}; below --min-match-fraction {min_match_fraction:g}"
                )
            print(f"[calibrate] WARNING: {message}", file=sys.stderr)
        extra_labels = int((~labels["_sample_key"].isin(scores["_sample_key"])).sum())
        if extra_labels:
            print(
                f"[calibrate] WARNING: labels CSV contains {extra_labels} rows "
                "not present in scores CSV",
                file=sys.stderr,
            )
    else:
        if label_col not in scores.columns:
            raise ValueError(f"scores CSV must contain {label_col!r} when labels CSV is omitted")
        merged = scores.copy()
        label_value_column = label_col

    merged["label_binary"] = merged[label_value_column].map(normalize_label)
    if merged["label_binary"].nunique() < 2:
        raise ValueError("Need at least one tumor and one normal labeled sample")
    return merged


def confusion_at_threshold(y_true, scores, threshold):
    return M.confusion_at(y_true, scores, threshold)


def metrics_at_threshold(y_true, scores, threshold, name):
    threshold = validate_threshold(threshold)
    row = M.classification_metrics(y_true, scores, threshold)
    return {
        "threshold_name": name,
        "threshold": float(threshold),
        "accuracy": row["accuracy"],
        "f1": row["f1"],
        "precision": row["precision"],
        "recall": row["recall"],
        "specificity": row["specificity"],
        "tn": row["tn"],
        "fp": row["fp"],
        "fn": row["fn"],
        "tp": row["tp"],
    }


def rank_auc(y_true, scores):
    return M.roc_auc(y_true, scores)


def choose_youden_threshold(y_true, scores):
    scores = np.asarray(scores, dtype=float)
    if not np.all(np.isfinite(scores)) or not np.all((scores >= 0) & (scores <= 1)):
        raise ValueError("scores must be finite probabilities in [0, 1]")
    if scores.size == 0:
        raise ValueError("Need at least one scored sample")
    best = M.youden_threshold(y_true, scores)
    return {
        "threshold_name": "youden_j",
        "threshold": best["threshold"],
        "accuracy": best["accuracy"],
        "f1": best["f1"],
        "precision": best["precision"],
        "recall": best["recall"],
        "specificity": best["specificity"],
        "tn": best["tn"],
        "fp": best["fp"],
        "fn": best["fn"],
        "tp": best["tp"],
        "youden_j": best["sensitivity"] + best["specificity"] - 1.0,
    }


def build_calibration_summary(y_true, scores, best, small_class_threshold=10):
    """Build the shared calibration summary with honest evaluation labeling."""
    y_true = np.asarray(y_true, dtype=int)
    n_tumor = int((y_true == 1).sum())
    n_normal = int((y_true == 0).sum())
    warnings = []
    if min(n_tumor, n_normal) < small_class_threshold:
        warnings.append(
            "Calibration metrics are unstable because at least one class has fewer "
            f"than {small_class_threshold} labeled samples "
            f"({n_tumor} tumor, {n_normal} normal)."
        )
    return {
        "n": int(len(y_true)),
        "n_tumor": n_tumor,
        "n_normal": n_normal,
        "auc": float(rank_auc(y_true, scores)),
        "recommended_threshold": float(best["threshold"]),
        "recommended_metric": "youden_j",
        "recommended_accuracy": float(best["accuracy"]),
        "recommended_recall": float(best["recall"]),
        "recommended_specificity": float(best["specificity"]),
        "evaluation_type": "apparent_resubstitution",
        "evaluation_note": (
            "Threshold and metrics were estimated and evaluated on the same labeled "
            "samples; they are apparent/resubstitution performance, not independent "
            "validation."
        ),
        "warnings": warnings,
    }


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("scores", help="CSV from score_tumor_normal.py")
    parser.add_argument("labels", nargs="?", help="CSV with sample + label columns")
    parser.add_argument("-o", "--output", help="threshold metrics CSV")
    parser.add_argument("--json-output", help="write compact calibration JSON")
    parser.add_argument("--sample-column", default="sample")
    parser.add_argument("--label-column", default="label")
    parser.add_argument("--default-threshold", type=float, default=0.5)
    parser.add_argument("--min-match-fraction", type=float, default=1.0,
                        help="minimum fraction of scored samples that must have labels (default 1.0)")
    parser.add_argument("--extra-threshold", type=float, action="append", default=[],
                        help="additional threshold to evaluate; may be repeated")
    args = parser.parse_args(argv)

    try:
        validate_threshold(args.default_threshold, "--default-threshold")
        for threshold in args.extra_threshold:
            validate_threshold(threshold, "--extra-threshold")
    except ValueError as exc:
        parser.error(str(exc))

    output = args.output or os.path.splitext(args.scores)[0] + ".thresholds.csv"
    try:
        ensure_distinct_paths(
            {"threshold metrics output": output, "calibration JSON output": args.json_output},
            {"scores input": args.scores, "labels input": args.labels},
        )
        data = load_scores_and_labels(
            args.scores, args.labels, args.sample_column, args.label_column,
            args.min_match_fraction,
        )
        y_true = data["label_binary"].to_numpy(dtype=int)
        scores = data["tumor_probability"].to_numpy(dtype=float)

        best = choose_youden_threshold(y_true, scores)
        rows = [
            metrics_at_threshold(y_true, scores, args.default_threshold, "default"),
            best,
        ]
        for threshold in args.extra_threshold:
            rows.append(
                metrics_at_threshold(
                    y_true, scores, threshold, f"threshold_{threshold:g}"
                )
            )
        metrics = pd.DataFrame(rows)
        summary = build_calibration_summary(y_true, scores, best)
        write_dataframe_csv(metrics, output, index=False)
        if args.json_output:
            write_json(summary, args.json_output, sort_keys=True)
    except ValueError as exc:
        parser.error(str(exc))

    print(f"[calibrate] n={summary['n']} tumor={summary['n_tumor']} normal={summary['n_normal']}")
    print(f"[calibrate] AUC={summary['auc']:.4f}")
    print(f"[calibrate] recommended threshold={summary['recommended_threshold']:.6f} "
          f"(accuracy={summary['recommended_accuracy']:.4f}, "
          f"recall={summary['recommended_recall']:.4f}, "
          f"specificity={summary['recommended_specificity']:.4f})")
    print(f"[calibrate] wrote {output}")
    print(
        "[calibrate] NOTE: reported calibration metrics are apparent/resubstitution "
        "estimates from the same labeled samples",
        file=sys.stderr,
    )
    for warning in summary["warnings"]:
        print(f"[calibrate] WARNING: {warning}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
