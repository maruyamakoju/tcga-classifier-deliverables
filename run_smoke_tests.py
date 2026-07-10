#!/usr/bin/env python3
"""Run lightweight release smoke tests."""
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent


def run(cmd, timeout_seconds=300):
    print("[smoke]", " ".join(cmd))
    return subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True,
                          encoding="utf-8", errors="replace", timeout=timeout_seconds)


def require_ok(result):
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr, file=sys.stderr)
        raise SystemExit(result.returncode)
    if result.stdout:
        print(result.stdout.rstrip())
    if result.stderr:
        print(result.stderr.rstrip(), file=sys.stderr)


def compare_outputs(observed_path, expected_path):
    observed = pd.read_csv(observed_path)
    expected = pd.read_csv(expected_path)
    if observed["sample"].tolist() != expected["sample"].tolist():
        raise SystemExit("[smoke] sample order mismatch")
    max_delta = float((observed["tumor_probability"]
                       - expected["tumor_probability"]).abs().max())
    if max_delta > 1e-6:
        raise SystemExit(f"[smoke] probability mismatch: max_delta={max_delta}")
    if not (observed["call"] == expected["call"]).all():
        raise SystemExit("[smoke] call mismatch")
    print(f"[smoke] scored output matches expected; max_delta={max_delta:.6g}")


def main():
    out_path = ROOT / "_smoke_scored.csv"
    thresholds_path = ROOT / "_smoke_thresholds.csv"
    explanations_path = ROOT / "_smoke_explanations.csv"
    qc_path = ROOT / "_smoke_qc.json"
    env_path = ROOT / "_smoke_environment.json"
    workflow_dir = ROOT / "_smoke_workflow"
    try:
        require_ok(run([sys.executable, "check_environment.py", "--self-test",
                        "-o", str(env_path)]))
        with open(env_path, "r", encoding="utf-8") as handle:
            env = json.load(handle)
        if env["status"] == "FAIL":
            raise SystemExit("[smoke] environment check failed")
        require_ok(run([sys.executable, "score_tumor_normal.py", "--self-test"]))
        require_ok(run([sys.executable, "inspect_expression_input.py", "example_input.csv",
                        "-o", str(qc_path)]))
        with open(qc_path, "r", encoding="utf-8") as handle:
            qc = json.load(handle)
        if qc["status"] != "PASS":
            raise SystemExit(f"[smoke] input QC did not pass: {qc['status']}")
        require_ok(run([sys.executable, "score_tumor_normal.py", "example_input.csv",
                        "-o", str(out_path)]))
        compare_outputs(out_path, ROOT / "example_output.csv")
        require_ok(run([sys.executable, "calibrate_threshold.py", str(out_path),
                        "example_labels.csv", "-o", str(thresholds_path)]))
        thresholds = pd.read_csv(thresholds_path)
        if "youden_j" not in thresholds["threshold_name"].tolist():
            raise SystemExit("[smoke] missing youden_j threshold row")
        require_ok(run([sys.executable, "explain_scores.py", "example_input.csv",
                        "-o", str(explanations_path), "--top-n", "3"]))
        explanations = pd.read_csv(explanations_path)
        expected_rows = 5 * 2 * 3
        if len(explanations) != expected_rows:
            raise SystemExit(f"[smoke] explanation row count mismatch: {len(explanations)}")
        print(f"[smoke] explanations generated: {len(explanations)} rows")
        require_ok(run([sys.executable, "run_tumor_normal_workflow.py", "example_input.csv",
                        "--labels", "example_labels.csv", "--output-dir",
                        str(workflow_dir), "--top-n", "3"]))
        compare_outputs(workflow_dir / "scores.csv", ROOT / "example_output.csv")
        with open(workflow_dir / "manifest.json", "r", encoding="utf-8") as handle:
            manifest = json.load(handle)
        if manifest["qc_status"] != "PASS":
            raise SystemExit(f"[smoke] workflow QC did not pass: {manifest['qc_status']}")
        workflow_explanations = pd.read_csv(workflow_dir / "explanations.csv")
        if len(workflow_explanations) != expected_rows:
            raise SystemExit("[smoke] workflow explanation row count mismatch")
        print("[smoke] workflow generated expected outputs")
        print("[smoke] PASS")
        return 0
    finally:
        for path in [out_path, thresholds_path, explanations_path, qc_path, env_path]:
            if path.exists():
                path.unlink()
        if workflow_dir.exists():
            for child in workflow_dir.iterdir():
                child.unlink()
            workflow_dir.rmdir()


if __name__ == "__main__":
    raise SystemExit(main())
