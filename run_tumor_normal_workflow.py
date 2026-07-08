#!/usr/bin/env python3
"""Run the complete lightweight tumor-vs-normal scoring workflow."""
import argparse
import json
import os
import sys
from pathlib import Path

import pandas as pd

from calibrate_threshold import (
    choose_youden_threshold,
    load_scores_and_labels,
    metrics_at_threshold,
    rank_auc,
    validate_threshold,
)
from explain_scores import explain_dataframe, load_gene_metadata
from inspect_expression_input import inspect_dataframe, load_reference
from score_tumor_normal import (
    load_lr_weights,
    print_invalid_alignment_summary,
    read_matrix,
    score_dataframe_lr_weights,
    validate_alignment_report,
)


def write_json(path, data):
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, sort_keys=True)
        handle.write("\n")


def read_json(path):
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def calibration_from_labels(scores_path, labels_path, sample_column, label_column,
                            default_threshold, extra_thresholds,
                            min_match_fraction):
    data = load_scores_and_labels(str(scores_path), str(labels_path), sample_column,
                                  label_column, min_match_fraction)
    y_true = data["label_binary"].to_numpy(dtype=int)
    scores = data["tumor_probability"].to_numpy(dtype=float)

    best = choose_youden_threshold(y_true, scores)
    rows = [
        metrics_at_threshold(y_true, scores, default_threshold, "default"),
        best,
    ]
    for threshold in extra_thresholds:
        rows.append(metrics_at_threshold(y_true, scores, threshold, f"threshold_{threshold:g}"))
    metrics = pd.DataFrame(rows)
    summary = {
        "n": int(len(data)),
        "n_tumor": int((y_true == 1).sum()),
        "n_normal": int((y_true == 0).sum()),
        "auc": float(rank_auc(y_true, scores)),
        "recommended_threshold": float(best["threshold"]),
        "recommended_metric": "youden_j",
        "recommended_accuracy": float(best["accuracy"]),
        "recommended_recall": float(best["recall"]),
        "recommended_specificity": float(best["specificity"]),
    }
    return metrics, summary


def markdown_table(rows):
    if not rows:
        return ""
    headers = list(rows[0].keys())
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(row.get(h, "")) for h in headers) + " |")
    return "\n".join(lines)


def format_float(value, digits=4):
    if value is None:
        return "NA"
    return f"{float(value):.{digits}f}"


def build_report(input_path, qc, scores, out_files, calibration_summary=None,
                 explanations_rows=None, stop_reason="QC"):
    gene = qc["gene_match"]
    value = qc["value_summary"]
    dist = qc["distribution_summary"]
    score_q = qc["score_summary"]["tumor_probability"]
    call_counts = scores["call"].value_counts().to_dict() if scores is not None else {}

    lines = [
        "# Tumor-vs-normal workflow report",
        "",
        f"Input: `{input_path}`",
        "",
        "## QC",
        "",
        f"- Status: **{qc['status']}**",
        f"- Model genes matched: {gene['matched_model_genes']}/{qc['shape']['model_genes']} "
        f"({gene['match_rate']:.1%})",
        f"- Expression median / p99 / max: {format_float(value['p50'])} / "
        f"{format_float(value['p99'])} / {format_float(value['max'])}",
        f"- |z| > 6 fraction: {format_float(dist['abs_z_gt_6_fraction'], 6)}",
        f"- Cohort gene-mean |z| p99: "
        f"{format_float(dist['cohort_gene_mean_abs_z']['p99'])}",
        "",
    ]
    if qc["messages"]:
        lines.append("### QC messages")
        lines.append("")
        for msg in qc["messages"]:
            lines.append(f"- {msg['level']}: {msg['message']}")
        lines.append("")
    else:
        lines.append("No QC warnings or errors.")
        lines.append("")

    lines.extend(["## Scores", ""])
    if scores is None:
        lines.extend([
            f"- Scoring was not run because the workflow stopped after {stop_reason}.",
            f"- Input samples inspected: {qc['shape']['samples']}",
            "",
        ])
    else:
        lines.extend([
            f"- Samples: {len(scores)}",
            f"- Tumor calls: {int(call_counts.get('tumor', 0))}",
            f"- Normal calls: {int(call_counts.get('normal', 0))}",
            f"- Tumor probability median / p90 / max: {format_float(score_q['p50'])} / "
            f"{format_float(score_q['p90'])} / {format_float(score_q['max'])}",
            "",
        ])

    if calibration_summary:
        lines.extend([
            "## Calibration",
            "",
            f"- Labeled samples: {calibration_summary['n']} "
            f"({calibration_summary['n_tumor']} tumor / "
            f"{calibration_summary['n_normal']} normal)",
            f"- AUC: {format_float(calibration_summary['auc'])}",
            f"- Recommended threshold: "
            f"{calibration_summary['recommended_threshold']:.6f} "
            f"({calibration_summary['recommended_metric']})",
            f"- Recommended accuracy / recall / specificity: "
            f"{format_float(calibration_summary['recommended_accuracy'])} / "
            f"{format_float(calibration_summary['recommended_recall'])} / "
            f"{format_float(calibration_summary['recommended_specificity'])}",
            "",
        ])

    if explanations_rows is not None:
        lines.extend([
            "## Explanations",
            "",
            f"- Explanation rows: {explanations_rows}",
            "- Positive rows push the LR logit toward tumor; negative rows push it toward normal.",
            "",
        ])

    file_rows = [{"file": name, "path": path} for name, path in out_files.items()]
    lines.extend([
        "## Output files",
        "",
        markdown_table(file_rows),
        "",
    ])
    return "\n".join(lines)


def default_output_dir(input_path):
    path = Path(input_path)
    return path.with_suffix("").name + "_tumor_normal_workflow"


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Run QC, scoring, optional calibration, explanations, and a report."
    )
    parser.add_argument("input", help="expression matrix (samples x genes)")
    parser.add_argument("-o", "--output-dir", default=None,
                        help="directory for workflow outputs")
    parser.add_argument("--labels", help="optional CSV with sample + label columns")
    parser.add_argument("--sample-column", default="sample")
    parser.add_argument("--label-column", default="label")
    parser.add_argument("--min-match-fraction", type=float, default=1.0,
                        help="minimum fraction of scored samples that must have labels (default 1.0)")
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--max-invalid-cell-fraction", type=float, default=0.0,
                        help=("maximum allowed missing, non-numeric, NaN, or infinite cells "
                              "among matched model genes before failing (default 0)"))
    parser.add_argument("--allow-invalid-values", action="store_true",
                        help=("warn instead of failing when matched model-gene cells are "
                              "missing, non-numeric, NaN, or infinite"))
    parser.add_argument("--extra-threshold", type=float, action="append", default=[])
    parser.add_argument("--top-n", type=int, default=10)
    parser.add_argument("--skip-explanations", action="store_true")
    parser.add_argument("--expected-class", choices=["unknown", "mixed", "tumor", "normal"],
                        default="unknown")
    parser.add_argument("--transpose", action="store_true")
    parser.add_argument("--strict-qc", action="store_true",
                        help="exit non-zero on QC WARN or FAIL")
    parser.add_argument("--allow-qc-fail", action="store_true",
                        help="continue even if QC status is FAIL")
    parser.add_argument("--lr-weights", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "deployable_lr_weights.npz"))
    parser.add_argument("--qc-reference", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "model_qc_reference.json"))
    parser.add_argument("--gene-metadata", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "model_gene_metadata.csv"))
    args = parser.parse_args(argv)

    if not 0 <= args.threshold <= 1:
        parser.error("--threshold must be between 0 and 1")
    if not 0 <= args.max_invalid_cell_fraction <= 1:
        parser.error("--max-invalid-cell-fraction must be between 0 and 1")
    for threshold in args.extra_threshold:
        try:
            validate_threshold(threshold, "--extra-threshold")
        except ValueError as exc:
            parser.error(str(exc))
    if args.top_n < 1 and not args.skip_explanations:
        parser.error("--top-n must be >= 1 unless --skip-explanations is used")

    output_dir = Path(args.output_dir or default_output_dir(args.input))
    output_dir.mkdir(parents=True, exist_ok=True)

    paths = {
        "qc_json": output_dir / "qc.json",
        "scores_csv": output_dir / "scores.csv",
        "report_md": output_dir / "workflow_report.md",
        "manifest_json": output_dir / "manifest.json",
    }
    if args.labels:
        paths["thresholds_csv"] = output_dir / "thresholds.csv"
        paths["calibration_json"] = output_dir / "calibration.json"
    if not args.skip_explanations:
        paths["explanations_csv"] = output_dir / "explanations.csv"

    weights = load_lr_weights(args.lr_weights)
    reference = load_reference(args.qc_reference)
    df = read_matrix(args.input, transpose=args.transpose)

    qc = inspect_dataframe(df, weights, threshold=args.threshold,
                           expected_class=args.expected_class, reference=reference)
    write_json(paths["qc_json"], qc)

    if qc["status"] == "FAIL" and not args.allow_qc_fail:
        manifest = {
            "status": "stopped_after_qc_fail",
            "input": args.input,
            "outputs": {"qc_json": paths["qc_json"].name},
            "qc_status": qc["status"],
        }
        write_json(paths["manifest_json"], manifest)
        report = build_report(args.input, qc, None,
                              {"qc_json": paths["qc_json"].name,
                               "manifest_json": paths["manifest_json"].name})
        paths["report_md"].write_text(report, encoding="utf-8")
        print(f"[workflow] QC status=FAIL; wrote {paths['qc_json']}")
        print("[workflow] stopped before scoring; pass --allow-qc-fail to continue")
        return 1

    scores, n_matched, missing, alignment_report = score_dataframe_lr_weights(
        df, weights, args.threshold, return_alignment_report=True
    )
    print_invalid_alignment_summary(alignment_report, sys.stderr)
    alignment_issues = validate_alignment_report(
        alignment_report,
        max_invalid_cell_fraction=args.max_invalid_cell_fraction,
    )
    if alignment_issues and not args.allow_invalid_values:
        manifest = {
            "status": "stopped_after_invalid_input",
            "input": args.input,
            "outputs": {"qc_json": paths["qc_json"].name},
            "qc_status": qc["status"],
            "alignment": {
                "invalid_matched_cells": int(alignment_report["invalid_matched_cells"]),
                "matched_cells": int(alignment_report["matched_cells"]),
                "invalid_matched_fraction": float(alignment_report["invalid_matched_fraction"]),
                "n_genes_with_invalid_values": int(
                    alignment_report["n_genes_with_invalid_values"]
                ),
                "n_samples_with_invalid_values": int(
                    alignment_report["n_samples_with_invalid_values"]
                ),
                "first_genes_with_invalid_values": alignment_report[
                    "first_genes_with_invalid_values"
                ],
                "first_samples_with_invalid_values": alignment_report[
                    "first_samples_with_invalid_values"
                ],
            },
        }
        write_json(paths["manifest_json"], manifest)
        report = build_report(
            args.input,
            qc,
            None,
            {"qc_json": paths["qc_json"].name,
             "manifest_json": paths["manifest_json"].name},
            stop_reason="invalid matched expression values",
        )
        paths["report_md"].write_text(report, encoding="utf-8")
        for issue in alignment_issues:
            print(f"[workflow] ERROR: {issue}", file=sys.stderr)
        print(
            "[workflow] stopped before writing scores; fix invalid values or "
            "pass --allow-invalid-values",
            file=sys.stderr,
        )
        return 1
    if alignment_issues:
        for issue in alignment_issues:
            print(f"[workflow] WARNING: {issue}", file=sys.stderr)
    scores.to_csv(paths["scores_csv"], index=False)

    calibration_summary = None
    if args.labels:
        threshold_metrics, calibration_summary = calibration_from_labels(
            paths["scores_csv"], args.labels, args.sample_column, args.label_column,
            args.threshold, args.extra_threshold, args.min_match_fraction
        )
        threshold_metrics.to_csv(paths["thresholds_csv"], index=False)
        write_json(paths["calibration_json"], calibration_summary)

    explanations_rows = None
    if not args.skip_explanations:
        gene_names = load_gene_metadata(args.gene_metadata)
        explanations, _, _ = explain_dataframe(df, weights, args.top_n, gene_names)
        explanations.to_csv(paths["explanations_csv"], index=False)
        explanations_rows = int(len(explanations))

    out_files = {name: path.name for name, path in paths.items()}
    report = build_report(args.input, qc, scores, out_files, calibration_summary,
                          explanations_rows)
    paths["report_md"].write_text(report, encoding="utf-8")

    manifest = {
        "status": "complete",
        "input": args.input,
        "threshold": args.threshold,
        "samples": int(df.shape[0]),
        "input_genes": int(df.shape[1]),
        "matched_model_genes": int(n_matched),
        "missing_model_genes": int(len(missing)),
        "qc_status": qc["status"],
        "tumor_calls": int((scores["call"] == "tumor").sum()),
        "normal_calls": int((scores["call"] == "normal").sum()),
        "labels": args.labels,
        "calibration": calibration_summary,
        "outputs": out_files,
    }
    write_json(paths["manifest_json"], manifest)

    print(f"[workflow] QC status={qc['status']}; matched {n_matched}/{len(weights['selected_genes'])} model genes")
    print(f"[workflow] wrote {paths['scores_csv']} ({manifest['tumor_calls']} tumor / {manifest['normal_calls']} normal)")
    if calibration_summary:
        print(f"[workflow] recommended threshold={calibration_summary['recommended_threshold']:.6f}")
    if explanations_rows is not None:
        print(f"[workflow] wrote {paths['explanations_csv']} ({explanations_rows} rows)")
    print(f"[workflow] wrote {paths['report_md']}")

    if args.strict_qc and qc["status"] != "PASS":
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
