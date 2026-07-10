#!/usr/bin/env python3
"""Run negative-path guardrail tests for the lightweight release."""
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

from release_tools.common import append_timeout_message, subprocess_output_text


ROOT = Path(__file__).resolve().parent


def run(cmd, timeout_seconds=300):
    print("[safety]", " ".join(str(x) for x in cmd))
    try:
        return subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True,
                              encoding="utf-8", errors="replace", timeout=timeout_seconds)
    except subprocess.TimeoutExpired as exc:
        return subprocess.CompletedProcess(
            cmd, 124,
            stdout=subprocess_output_text(exc.stdout),
            stderr=append_timeout_message(subprocess_output_text(exc.stderr), timeout_seconds),
        )


def require(condition, message):
    if not condition:
        raise SystemExit(f"[safety] FAIL: {message}")


def require_ok(result, label):
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr, file=sys.stderr)
        raise SystemExit(f"[safety] FAIL: {label} returned {result.returncode}")


def require_fail(result, label):
    if result.returncode == 0:
        print(result.stdout)
        print(result.stderr, file=sys.stderr)
        raise SystemExit(f"[safety] FAIL: {label} unexpectedly succeeded")


def message_codes(qc_report):
    return {message["code"] for message in qc_report.get("messages", [])}


def write_no_gene_match_input(example, out_path):
    bad = pd.DataFrame(
        {
            "NOT_A_MODEL_GENE_1": [1.0] * len(example),
            "NOT_A_MODEL_GENE_2": [2.0] * len(example),
        },
        index=example.index,
    )
    bad.to_csv(out_path)


def write_raw_count_like_input(example, out_path):
    raw_like = example.copy()
    raw_like.iloc[:, :] = 1000.0
    raw_like.to_csv(out_path)


def write_invalid_matched_input(example, out_path):
    invalid = example.copy().astype(object)
    invalid.iloc[:, 0] = "not_numeric"
    invalid.to_csv(out_path)


def main():
    example = pd.read_csv(ROOT / "example_input.csv", index_col=0)
    temp_root = Path(tempfile.mkdtemp(prefix="tcga_safety_", dir=ROOT))
    try:
        no_match_input = temp_root / "no_model_genes.csv"
        no_match_qc = temp_root / "no_model_genes.qc.json"
        no_match_scores = temp_root / "no_model_genes.scored.csv"
        no_match_explanations = temp_root / "no_model_genes.explanations.csv"
        no_match_adapted_scores = temp_root / "no_model_genes.adapted_scores.csv"
        no_match_workflow = temp_root / "no_model_genes_workflow"
        raw_like_input = temp_root / "raw_counts_like.csv"
        raw_like_qc = temp_root / "raw_counts_like.qc.json"
        expected_normal_qc = temp_root / "expected_normal.qc.json"
        invalid_input = temp_root / "invalid_matched_values.csv"
        invalid_workflow = temp_root / "invalid_matched_workflow"
        bad_labels = temp_root / "bad_labels.csv"
        bad_labels_workflow = temp_root / "bad_labels_workflow"
        missing_input = temp_root / "missing_expression.csv"
        missing_output = temp_root / "missing_expression.scored.csv"
        pickle_input = temp_root / "untrusted_expression.pkl"
        pickle_output = temp_root / "untrusted_expression.scored.csv"

        write_no_gene_match_input(example, no_match_input)
        write_raw_count_like_input(example, raw_like_input)
        write_invalid_matched_input(example, invalid_input)
        bad_labels.write_text(
            f"sample,label\n{example.index[0]},tumor\n",
            encoding="utf-8",
        )
        example.to_pickle(pickle_input)

        result = run([sys.executable, "score_tumor_normal.py", "example_input.csv",
                      "--threshold", "1.5"])
        require_fail(result, "out-of-range threshold")
        require("threshold must be between 0 and 1" in result.stderr,
                "threshold validation message missing")

        result = run([sys.executable, "score_tumor_normal.py", "example_input.csv",
                      "--use-pickle-lr"])
        require_fail(result, "legacy pickle LR")
        require("legacy pickle/RF scoring is not available" in result.stderr,
                "legacy scorer rejection message missing")

        result = run([sys.executable, "score_tumor_normal.py", "example_input.csv",
                      "--model", "rf"])
        require_fail(result, "legacy RF scorer")
        require("invalid choice" in result.stderr and "'rf'" in result.stderr,
                "legacy RF rejection message missing")

        result = run([sys.executable, "score_tumor_normal.py", str(pickle_input),
                      "-o", str(pickle_output)])
        require_fail(result, "pickle expression input")
        require("Pickle expression inputs are disabled by default" in result.stderr,
                "pickle expression input rejection message missing")
        require(not pickle_output.exists(),
                "score_tumor_normal.py wrote scores after pickle input rejection")

        result = run([sys.executable, "score_tumor_normal.py", str(missing_input),
                      "-o", str(missing_output)])
        require_fail(result, "missing expression input")
        require("expression matrix file not found" in result.stderr,
                "missing expression input message missing")
        require("Traceback" not in result.stderr,
                "missing expression input produced a traceback")
        require(not missing_output.exists(),
                "score_tumor_normal.py wrote scores after missing input rejection")

        result = run([sys.executable, "score_tumor_normal.py", str(no_match_input),
                      "-o", str(no_match_scores)])
        require_fail(result, "low gene coverage scorer")
        require("low model-gene coverage" in result.stderr,
                "low gene coverage scorer message missing")
        require("Refusing to write scores" in result.stderr,
                "low gene coverage scorer refusal missing")
        require(not no_match_scores.exists(),
                "score_tumor_normal.py wrote scores after low gene coverage")

        result = run([sys.executable, "score_tumor_normal.py", str(no_match_input),
                      "-o", str(no_match_scores), "--allow-low-gene-coverage"])
        require_ok(result, "low gene coverage scorer explicit allow")

        result = run([sys.executable, "score_tumor_normal.py", str(invalid_input)])
        require_fail(result, "invalid matched values scorer")
        require("invalid matched values" in result.stderr,
                "invalid matched-value summary missing")
        require("Refusing to write scores" in result.stderr,
                "invalid matched-value refusal missing")

        result = run([sys.executable, "score_tumor_normal.py", str(invalid_input),
                      "--allow-invalid-values"])
        require_ok(result, "invalid matched values explicit allow")

        invalid_explanations = invalid_input.with_suffix(".explanations.csv")
        result = run([sys.executable, "explain_scores.py", str(invalid_input)])
        require_fail(result, "invalid matched values explanations")
        require("invalid matched values" in result.stderr,
                "invalid matched-value explanation summary missing")
        require("Refusing to write explanations" in result.stderr,
                "invalid matched-value explanation refusal missing")
        require(not invalid_explanations.exists(),
                "explain_scores.py wrote explanations after invalid matched values")

        result = run([sys.executable, "explain_scores.py", str(invalid_input),
                      "--allow-invalid-values"])
        require_ok(result, "invalid matched values explanations explicit allow")

        result = run([sys.executable, "explain_scores.py", str(no_match_input),
                      "-o", str(no_match_explanations)])
        require_fail(result, "low gene coverage explanations")
        require("low model-gene coverage" in result.stderr,
                "low gene coverage explanation message missing")
        require("Refusing to write explanations" in result.stderr,
                "low gene coverage explanation refusal missing")
        require(not no_match_explanations.exists(),
                "explain_scores.py wrote explanations after low gene coverage")

        invalid_adapted_scores = invalid_input.with_suffix(".adapted_scores.csv")
        result = run([sys.executable, "cohort_adapt_score.py", str(invalid_input),
                      "--adapt", "none"])
        require_fail(result, "invalid matched values cohort adaptation")
        require("invalid matched values" in result.stderr,
                "invalid matched-value adaptation summary missing")
        require("Refusing to write adapted scores" in result.stderr,
                "invalid matched-value adaptation refusal missing")
        require(not invalid_adapted_scores.exists(),
                "cohort_adapt_score.py wrote scores after invalid matched values")

        result = run([sys.executable, "cohort_adapt_score.py", str(invalid_input),
                      "--adapt", "none", "--allow-invalid-values"])
        require_ok(result, "invalid matched values cohort adaptation explicit allow")

        result = run([sys.executable, "cohort_adapt_score.py", str(no_match_input),
                      "--adapt", "none", "--out", str(no_match_adapted_scores)])
        require_fail(result, "low gene coverage cohort adaptation")
        require("low model-gene coverage" in result.stderr,
                "low gene coverage adaptation message missing")
        require("Refusing to write adapted scores" in result.stderr,
                "low gene coverage adaptation refusal missing")
        require(not no_match_adapted_scores.exists(),
                "cohort_adapt_score.py wrote scores after low gene coverage")

        result = run([sys.executable, "explain_scores.py", "example_input.csv",
                      "--top-n", "0"])
        require_fail(result, "invalid explanation top-n")
        require("--top-n must be >= 1" in result.stderr,
                "top-n validation message missing")

        result = run([sys.executable, "inspect_expression_input.py",
                      str(no_match_input), "-o", str(no_match_qc)])
        require_fail(result, "QC no model genes")
        no_match_report = json.loads(no_match_qc.read_text(encoding="utf-8"))
        require(no_match_report["status"] == "FAIL", "no-match QC did not FAIL")
        require("no_model_genes_matched" in message_codes(no_match_report),
                "no-match QC code missing")

        result = run([sys.executable, "run_tumor_normal_workflow.py",
                      str(no_match_input), "--output-dir", str(no_match_workflow)])
        require_fail(result, "workflow no model genes")
        manifest = json.loads((no_match_workflow / "manifest.json").read_text(encoding="utf-8"))
        require(manifest["status"] == "stopped_after_qc_fail",
                "workflow did not stop after QC FAIL")
        require(not (no_match_workflow / "scores.csv").exists(),
                "workflow wrote scores.csv after QC FAIL")

        result = run([sys.executable, "run_tumor_normal_workflow.py",
                      str(invalid_input), "--output-dir", str(invalid_workflow)])
        require_fail(result, "workflow invalid matched values")
        manifest = json.loads((invalid_workflow / "manifest.json").read_text(encoding="utf-8"))
        require(manifest["status"] == "stopped_after_invalid_input",
                "workflow did not stop after invalid matched values")
        require(not (invalid_workflow / "scores.csv").exists(),
                "workflow wrote scores.csv after invalid matched values")

        result = run([sys.executable, "run_tumor_normal_workflow.py",
                      "example_input.csv", "--labels", str(bad_labels),
                      "--output-dir", str(bad_labels_workflow)])
        require_fail(result, "workflow bad calibration labels")
        require("calibration failed" in result.stderr,
                "workflow calibration failure message missing")
        manifest = json.loads((bad_labels_workflow / "manifest.json").read_text(encoding="utf-8"))
        require(manifest["status"] == "stopped_after_calibration_error",
                "workflow did not record calibration failure status")
        require("calibration_error" in manifest,
                "workflow calibration failure did not record error")
        require((bad_labels_workflow / "scores.csv").exists(),
                "workflow did not preserve scores after calibration failure")
        require(not (bad_labels_workflow / "thresholds.csv").exists(),
                "workflow wrote thresholds.csv after calibration failure")
        require(not (bad_labels_workflow / "calibration.json").exists(),
                "workflow wrote calibration.json after calibration failure")
        require(not (bad_labels_workflow / "explanations.csv").exists(),
                "workflow wrote explanations.csv after calibration failure")

        result = run([sys.executable, "inspect_expression_input.py",
                      str(raw_like_input), "-o", str(raw_like_qc)])
        require_fail(result, "QC raw-count-like values")
        raw_report = json.loads(raw_like_qc.read_text(encoding="utf-8"))
        require(raw_report["status"] == "FAIL", "raw-count-like QC did not FAIL")
        raw_codes = message_codes(raw_report)
        require("expression_values_too_large" in raw_codes, "raw p99 error missing")
        require("expression_max_too_large" in raw_codes, "raw max error missing")

        result = run([sys.executable, "inspect_expression_input.py",
                      "example_input.csv", "-o", str(expected_normal_qc),
                      "--expected-class", "normal", "--strict"])
        require_fail(result, "strict expected-normal warning")
        expected_report = json.loads(expected_normal_qc.read_text(encoding="utf-8"))
        require(expected_report["status"] == "WARN", "expected-normal QC did not WARN")
        require("unexpected_tumor_calls" in message_codes(expected_report),
                "expected-normal warning code missing")

        release_dir = ROOT / "release-lite"
        zip_path = ROOT / "tcga-tumor-normal-release-lite.zip"
        if release_dir.exists() and zip_path.exists():
            result = run([sys.executable, "validate_release_lite.py",
                          "--release-dir", str(release_dir), "--zip", str(zip_path)])
            require_ok(result, "release validator")

        print("[safety] PASS")
        return 0
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
