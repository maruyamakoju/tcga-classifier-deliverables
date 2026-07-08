#!/usr/bin/env python3
"""Run negative-path guardrail tests for the lightweight release."""
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent


def run(cmd):
    print("[safety]", " ".join(str(x) for x in cmd))
    return subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True)


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


def main():
    example = pd.read_csv(ROOT / "example_input.csv", index_col=0)
    temp_root = Path(tempfile.mkdtemp(prefix="tcga_safety_", dir=ROOT))
    try:
        no_match_input = temp_root / "no_model_genes.csv"
        no_match_qc = temp_root / "no_model_genes.qc.json"
        no_match_workflow = temp_root / "no_model_genes_workflow"
        raw_like_input = temp_root / "raw_counts_like.csv"
        raw_like_qc = temp_root / "raw_counts_like.qc.json"
        expected_normal_qc = temp_root / "expected_normal.qc.json"

        write_no_gene_match_input(example, no_match_input)
        write_raw_count_like_input(example, raw_like_input)

        result = run([sys.executable, "score_tumor_normal.py", "example_input.csv",
                      "--threshold", "1.5"])
        require_fail(result, "out-of-range threshold")
        require("threshold must be between 0 and 1" in result.stderr,
                "threshold validation message missing")

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
