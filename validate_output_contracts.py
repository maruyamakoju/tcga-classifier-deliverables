#!/usr/bin/env python3
"""Validate stable CSV/JSON output contracts shipped in the release."""
import argparse
import json
import sys
from collections.abc import Mapping
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent

SCORE_COLUMNS = ["sample", "tumor_probability", "call"]
THRESHOLD_COLUMNS = [
    "threshold_name", "threshold", "accuracy", "f1", "precision", "recall",
    "specificity", "tn", "fp", "fn", "tp", "youden_j",
]
EXPLANATION_COLUMNS = [
    "sample", "tumor_probability", "logit", "direction", "rank", "gene_id",
    "gene_name", "contribution_logit", "expression_log2_tpm1", "training_mean",
    "scaled_value", "lr_coef",
]
GENE_METADATA_COLUMNS = [
    "rank_abs_lr_coef", "gene_id", "gene_id_base", "gene_name", "lr_coef",
    "abs_lr_coef", "direction_if_high", "scaler_mean", "scaler_scale",
]

QC_TOP_LEVEL_KEYS = {
    "distribution_summary",
    "expected_class",
    "gene_match",
    "messages",
    "reference_source",
    "score_summary",
    "shape",
    "status",
    "threshold",
    "value_summary",
}
MANIFEST_TOP_LEVEL_KEYS = {
    "input",
    "input_genes",
    "matched_model_genes",
    "missing_model_genes",
    "normal_calls",
    "outputs",
    "qc_status",
    "samples",
    "status",
    "threshold",
    "tumor_calls",
}
CALIBRATION_KEYS = {
    "auc",
    "n",
    "n_normal",
    "n_tumor",
    "recommended_accuracy",
    "recommended_metric",
    "recommended_recall",
    "recommended_specificity",
    "recommended_threshold",
}


def add_message(messages, level, code, message, path=None):
    item = {"level": level, "code": code, "message": message}
    if path is not None:
        item["path"] = str(path)
    messages.append(item)


def require_file(rel, messages):
    path = ROOT / rel
    if not path.exists():
        add_message(messages, "ERROR", "required_file_missing",
                    f"Required contract file is missing: {rel}", path)
        return None
    return path


def read_csv(rel, messages):
    path = require_file(rel, messages)
    if path is None:
        return None
    try:
        return pd.read_csv(path)
    except Exception as exc:
        add_message(messages, "ERROR", "csv_read_failed",
                    f"Could not read {rel}: {exc}", path)
        return None


def read_json(rel, messages):
    path = require_file(rel, messages)
    if path is None:
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        add_message(messages, "ERROR", "json_read_failed",
                    f"Could not read {rel}: {exc}", path)
        return None


def require_json_object(value, rel, messages):
    if not isinstance(value, Mapping):
        add_message(messages, "ERROR", "json_top_level_not_object",
                    f"{rel} top-level value must be a JSON object.", ROOT / rel)
        return None
    return value


def check_numeric_range(value, rel, key, messages, minimum=0, maximum=1):
    if isinstance(value, bool):
        numeric = None
    else:
        try:
            numeric = float(value)
        except (TypeError, ValueError):
            numeric = None
    if numeric is None or not (minimum <= numeric <= maximum):
        add_message(messages, "ERROR", "json_metric_out_of_range",
                    f"{rel}:{key} must be numeric and in [{minimum}, {maximum}].",
                    ROOT / rel)


def check_columns(df, expected, rel, messages):
    actual = list(df.columns)
    if actual != expected:
        add_message(messages, "ERROR", "column_contract_mismatch",
                    f"{rel} columns changed. expected={expected} actual={actual}", ROOT / rel)


def check_probability_series(series, rel, column, messages):
    values = pd.to_numeric(series, errors="coerce")
    if values.isna().any():
        add_message(messages, "ERROR", "non_numeric_probability",
                    f"{rel}:{column} contains non-numeric values.", ROOT / rel)
    elif not ((values >= 0) & (values <= 1)).all():
        add_message(messages, "ERROR", "probability_out_of_range",
                    f"{rel}:{column} must be in [0, 1].", ROOT / rel)


def check_score_csv(rel, messages):
    df = read_csv(rel, messages)
    if df is None:
        return None
    check_columns(df, SCORE_COLUMNS, rel, messages)
    if "sample" in df:
        if df["sample"].isna().any() or (df["sample"].astype(str).str.len() == 0).any():
            add_message(messages, "ERROR", "empty_sample_id",
                        f"{rel}:sample contains empty values.", ROOT / rel)
        if df["sample"].duplicated().any():
            add_message(messages, "ERROR", "duplicate_sample_id",
                        f"{rel}:sample contains duplicate IDs.", ROOT / rel)
    if "tumor_probability" in df:
        check_probability_series(df["tumor_probability"], rel, "tumor_probability", messages)
    if "call" in df and not set(df["call"]).issubset({"tumor", "normal"}):
        add_message(messages, "ERROR", "invalid_call_value",
                    f"{rel}:call must contain only tumor/normal.", ROOT / rel)
    return df


def check_labels(messages):
    df = read_csv("example_labels.csv", messages)
    if df is None:
        return
    check_columns(df, ["sample", "label"], "example_labels.csv", messages)
    if "label" in df and not set(df["label"]).issubset({"tumor", "normal", 1, 0, "1", "0"}):
        add_message(messages, "ERROR", "invalid_label_value",
                    "example_labels.csv:label has unsupported values.", ROOT / "example_labels.csv")


def check_thresholds(messages):
    df = read_csv("example_workflow_output/thresholds.csv", messages)
    if df is None:
        return
    check_columns(df, THRESHOLD_COLUMNS, "example_workflow_output/thresholds.csv", messages)
    if "threshold_name" in df:
        names = set(df["threshold_name"])
        if not {"default", "youden_j"}.issubset(names):
            add_message(messages, "ERROR", "threshold_rows_missing",
                        "thresholds.csv must contain default and youden_j rows.",
                        ROOT / "example_workflow_output/thresholds.csv")
    for column in ["threshold", "accuracy", "f1", "precision", "recall", "specificity", "youden_j"]:
        if column in df:
            values = pd.to_numeric(df[column], errors="coerce").dropna()
            if not ((values >= 0) & (values <= 1)).all():
                add_message(messages, "ERROR", "metric_out_of_range",
                            f"thresholds.csv:{column} must be in [0, 1] when present.",
                            ROOT / "example_workflow_output/thresholds.csv")
    for column in ["tn", "fp", "fn", "tp"]:
        if column in df:
            values = pd.to_numeric(df[column], errors="coerce")
            if values.isna().any():
                add_message(messages, "ERROR", "non_numeric_count",
                            f"thresholds.csv:{column} must be numeric.",
                            ROOT / "example_workflow_output/thresholds.csv")
            elif (values < 0).any():
                add_message(messages, "ERROR", "negative_count",
                            f"thresholds.csv:{column} must be non-negative.",
                            ROOT / "example_workflow_output/thresholds.csv")


def check_explanations(messages):
    df = read_csv("example_workflow_output/explanations.csv", messages)
    if df is None:
        return
    check_columns(df, EXPLANATION_COLUMNS, "example_workflow_output/explanations.csv", messages)
    if len(df) != 30:
        add_message(messages, "ERROR", "example_explanation_row_count",
                    f"Example explanations should have 30 rows; found {len(df)}.",
                    ROOT / "example_workflow_output/explanations.csv")
    if "direction" in df and not set(df["direction"]).issubset({"tumor", "normal"}):
        add_message(messages, "ERROR", "invalid_explanation_direction",
                    "explanations.csv:direction must contain only tumor/normal.",
                    ROOT / "example_workflow_output/explanations.csv")
    if "rank" in df:
        ranks = pd.to_numeric(df["rank"], errors="coerce")
        if ranks.isna().any() or (ranks < 1).any():
            add_message(messages, "ERROR", "invalid_explanation_rank",
                        "explanations.csv:rank must be positive integers.",
                        ROOT / "example_workflow_output/explanations.csv")
    for column in ["tumor_probability"]:
        if column in df:
            check_probability_series(df[column], "example_workflow_output/explanations.csv",
                                     column, messages)


def check_json_contracts(messages):
    qc = read_json("example_workflow_output/qc.json", messages)
    if qc is not None:
        qc = require_json_object(qc, "example_workflow_output/qc.json", messages)
    if qc is not None:
        missing = QC_TOP_LEVEL_KEYS - set(qc)
        if missing:
            add_message(messages, "ERROR", "qc_keys_missing",
                        f"qc.json missing keys: {sorted(missing)}",
                        ROOT / "example_workflow_output/qc.json")
        if qc.get("status") not in {"PASS", "WARN", "FAIL"}:
            add_message(messages, "ERROR", "invalid_qc_status",
                        "qc.json:status must be PASS/WARN/FAIL.",
                        ROOT / "example_workflow_output/qc.json")
        gene_match = qc.get("gene_match", {})
        if gene_match.get("matched_model_genes") != 2000:
            add_message(messages, "ERROR", "example_qc_gene_match_changed",
                        "Example QC should match 2000 model genes.",
                        ROOT / "example_workflow_output/qc.json")

    manifest = read_json("example_workflow_output/manifest.json", messages)
    if manifest is not None:
        manifest = require_json_object(
            manifest,
            "example_workflow_output/manifest.json",
            messages,
        )
    if manifest is not None:
        missing = MANIFEST_TOP_LEVEL_KEYS - set(manifest)
        if missing:
            add_message(messages, "ERROR", "manifest_keys_missing",
                        f"manifest.json missing keys: {sorted(missing)}",
                        ROOT / "example_workflow_output/manifest.json")
        if manifest.get("status") != "complete" or manifest.get("qc_status") != "PASS":
            add_message(messages, "ERROR", "example_manifest_status_changed",
                        "Example workflow manifest should be complete with QC PASS.",
                        ROOT / "example_workflow_output/manifest.json")
        outputs = manifest.get("outputs", {})
        if not isinstance(outputs, Mapping):
            add_message(messages, "ERROR", "manifest_outputs_not_object",
                        "manifest.json:outputs must be a JSON object.",
                        ROOT / "example_workflow_output/manifest.json")
        else:
            for key, rel in outputs.items():
                if not isinstance(rel, str):
                    add_message(messages, "ERROR", "manifest_output_path_not_string",
                                f"manifest output path must be a string: {key}={rel!r}",
                                ROOT / "example_workflow_output/manifest.json")
                    continue
                if not (ROOT / "example_workflow_output" / rel).exists():
                    add_message(messages, "ERROR", "manifest_output_missing",
                                f"manifest output is missing: {key}={rel}",
                                ROOT / "example_workflow_output" / rel)

    calibration = read_json("example_workflow_output/calibration.json", messages)
    if calibration is not None:
        calibration = require_json_object(
            calibration,
            "example_workflow_output/calibration.json",
            messages,
        )
    if calibration is not None:
        missing = CALIBRATION_KEYS - set(calibration)
        if missing:
            add_message(messages, "ERROR", "calibration_keys_missing",
                        f"calibration.json missing keys: {sorted(missing)}",
                        ROOT / "example_workflow_output/calibration.json")
        for key in ["auc", "recommended_accuracy", "recommended_recall",
                    "recommended_specificity", "recommended_threshold"]:
            check_numeric_range(
                calibration.get(key),
                "example_workflow_output/calibration.json",
                key,
                messages,
            )


def check_gene_metadata(messages):
    df = read_csv("model_gene_metadata.csv", messages)
    if df is None:
        return
    check_columns(df, GENE_METADATA_COLUMNS, "model_gene_metadata.csv", messages)
    if len(df) != 2000:
        add_message(messages, "ERROR", "model_gene_count_changed",
                    f"model_gene_metadata.csv should have 2000 rows; found {len(df)}.",
                    ROOT / "model_gene_metadata.csv")
    if "gene_id" in df and df["gene_id"].duplicated().any():
        add_message(messages, "ERROR", "duplicate_model_gene_id",
                    "model_gene_metadata.csv:gene_id contains duplicates.",
                    ROOT / "model_gene_metadata.csv")
    if "direction_if_high" in df and not set(df["direction_if_high"]).issubset({"tumor", "normal"}):
        add_message(messages, "ERROR", "invalid_model_gene_direction",
                    "model_gene_metadata.csv:direction_if_high must be tumor/normal.",
                    ROOT / "model_gene_metadata.csv")
    if "scaler_scale" in df and (pd.to_numeric(df["scaler_scale"], errors="coerce") <= 0).any():
        add_message(messages, "ERROR", "invalid_scaler_scale",
                    "model_gene_metadata.csv:scaler_scale must be positive.",
                    ROOT / "model_gene_metadata.csv")


def check_score_consistency(messages):
    expected = check_score_csv("example_output.csv", messages)
    workflow = check_score_csv("example_workflow_output/scores.csv", messages)
    if expected is None or workflow is None:
        return
    if expected["sample"].tolist() != workflow["sample"].tolist():
        add_message(messages, "ERROR", "example_score_sample_order_changed",
                    "example_output.csv and workflow scores.csv sample order differs.",
                    ROOT / "example_workflow_output/scores.csv")
    delta = (expected["tumor_probability"] - workflow["tumor_probability"]).abs().max()
    if float(delta) > 1e-6:
        add_message(messages, "ERROR", "example_score_probability_changed",
                    f"example_output.csv and workflow scores.csv differ by {delta}.",
                    ROOT / "example_workflow_output/scores.csv")
    if not (expected["call"] == workflow["call"]).all():
        add_message(messages, "ERROR", "example_score_call_changed",
                    "example_output.csv and workflow scores.csv calls differ.",
                    ROOT / "example_workflow_output/scores.csv")


def check_qc_reference(messages):
    report = read_json("model_qc_reference.json", messages)
    if report is None:
        return
    report = require_json_object(report, "model_qc_reference.json", messages)
    if report is None:
        return
    for key in ["intended_input", "reference_reports", "rules"]:
        if key not in report:
            add_message(messages, "ERROR", "qc_reference_key_missing",
                        f"model_qc_reference.json missing key: {key}",
                        ROOT / "model_qc_reference.json")


def build_report():
    messages = []
    check_score_consistency(messages)
    check_labels(messages)
    check_thresholds(messages)
    check_explanations(messages)
    check_json_contracts(messages)
    check_gene_metadata(messages)
    check_qc_reference(messages)
    levels = {item["level"] for item in messages}
    status = "FAIL" if "ERROR" in levels else "WARN" if "WARNING" in levels else "PASS"
    return {
        "schema_version": "1.0",
        "status": status,
        "root": str(ROOT),
        "messages": messages,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description="Validate bundled output schemas and contracts.")
    parser.add_argument("-o", "--output", help="write JSON report")
    parser.add_argument("--strict", action="store_true",
                        help="return non-zero on warnings as well as errors")
    args = parser.parse_args(argv)

    report = build_report()
    for message in report["messages"]:
        stream = sys.stderr if message["level"] in {"ERROR", "WARNING"} else sys.stdout
        print(f"[contracts] {message['level']}: {message['message']}", file=stream)
    print(f"[contracts] status={report['status']}")
    if args.output:
        out_path = Path(args.output)
        if not out_path.is_absolute():
            out_path = ROOT / out_path
        out_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                            encoding="utf-8")
        print(f"[contracts] wrote {out_path}")
    if report["status"] == "FAIL" or (args.strict and report["status"] == "WARN"):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
