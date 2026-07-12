"""Subprocess-level failure contracts for the public lightweight CLIs."""
import hashlib
import os
import shutil
import subprocess
import sys

import numpy as np
import pytest


def run_cli(root, *arguments):
    return subprocess.run(
        [sys.executable, *map(str, arguments)],
        cwd=root,
        text=True,
        capture_output=True,
    )


def assert_clean_failure(result):
    assert result.returncode != 0
    assert "Traceback" not in result.stderr


@pytest.mark.parametrize(
    "arguments",
    [
        ("score_tumor_normal.py", "example_input.csv", "--lr-weights", "missing.npz"),
        ("inspect_expression_input.py", "example_input.csv", "--lr-weights", "missing.npz"),
        ("explain_scores.py", "example_input.csv", "--lr-weights", "missing.npz"),
        ("cohort_adapt_score.py", "example_input.csv", "--weights", "missing.npz"),
        ("run_tumor_normal_workflow.py", "example_input.csv", "--lr-weights", "missing.npz"),
        ("calibrate_threshold.py", "missing_scores.csv"),
    ],
)
def test_missing_public_inputs_never_print_tracebacks(root, arguments):
    result = run_cli(root, *arguments)
    assert_clean_failure(result)


def test_trusted_weight_export_missing_pipeline_is_a_clean_failure(tmp_path, root):
    output = tmp_path / "weights.npz"
    result = run_cli(
        root,
        "export_lr_weights.py",
        "--trusted-pipeline",
        "--pipeline",
        tmp_path / "missing.pkl",
        "--output",
        output,
    )

    assert_clean_failure(result)
    assert "trusted legacy pipeline" in result.stderr
    assert not output.exists()


@pytest.mark.parametrize(
    ("script", "weight_flag", "extra"),
    [
        ("score_tumor_normal.py", "--lr-weights", ()),
        ("inspect_expression_input.py", "--lr-weights", ()),
        ("explain_scores.py", "--lr-weights", ()),
        ("cohort_adapt_score.py", "--weights", ("--adapt", "none")),
        ("run_tumor_normal_workflow.py", "--lr-weights", ()),
    ],
)
def test_missing_key_model_never_prints_traceback(
    tmp_path, root, script, weight_flag, extra
):
    malformed = tmp_path / "missing_coef.npz"
    np.savez(
        malformed,
        selected_genes=np.array(["g1"]),
        scaler_mean=np.array([0.0]),
        scaler_scale=np.array([1.0]),
        intercept=np.array(0.0),
    )
    result = run_cli(
        root, script, os.path.join(root, "example_input.csv"), weight_flag, malformed, *extra
    )

    assert_clean_failure(result)
    assert "missing required arrays: coef" in result.stderr


@pytest.mark.parametrize(
    ("script", "output_flag", "extra"),
    [
        ("score_tumor_normal.py", "-o", ()),
        ("inspect_expression_input.py", "-o", ()),
        ("explain_scores.py", "-o", ()),
        ("cohort_adapt_score.py", "--out", ("--adapt", "none")),
    ],
)
def test_cli_refuses_to_overwrite_expression_input(
    tmp_path, root, script, output_flag, extra
):
    expression = tmp_path / "expression.csv"
    shutil.copyfile(os.path.join(root, "example_input.csv"), expression)
    before = hashlib.sha256(expression.read_bytes()).hexdigest()

    result = run_cli(root, script, expression, output_flag, expression, *extra)

    assert_clean_failure(result)
    assert "refusing to overwrite expression input" in result.stderr
    assert hashlib.sha256(expression.read_bytes()).hexdigest() == before


def test_calibration_refuses_to_overwrite_scores_input(tmp_path, root):
    scores = tmp_path / "scores.csv"
    shutil.copyfile(os.path.join(root, "example_output.csv"), scores)
    labels = os.path.join(root, "example_labels.csv")
    before = hashlib.sha256(scores.read_bytes()).hexdigest()

    result = run_cli(root, "calibrate_threshold.py", scores, labels, "-o", scores)

    assert_clean_failure(result)
    assert "refusing to overwrite scores input" in result.stderr
    assert hashlib.sha256(scores.read_bytes()).hexdigest() == before


def test_calibration_accepts_scores_with_stale_label_when_labels_file_is_supplied(
    tmp_path, root
):
    scores = tmp_path / "scores.csv"
    labels = tmp_path / "labels.csv"
    output = tmp_path / "thresholds.csv"
    scores.write_text(
        "sample,tumor_probability,label\ns1,0.1,stale\ns2,0.9,stale\n",
        encoding="utf-8",
    )
    labels.write_text("sample,label\ns1,normal\ns2,tumor\n", encoding="utf-8")

    result = run_cli(root, "calibrate_threshold.py", scores, labels, "-o", output)

    assert result.returncode == 0, result.stderr
    assert "Traceback" not in result.stderr
    assert output.is_file()


def test_workflow_refuses_managed_output_collision_with_input(tmp_path, root):
    output_dir = tmp_path / "workflow"
    output_dir.mkdir()
    expression = output_dir / "scores.csv"
    shutil.copyfile(os.path.join(root, "example_input.csv"), expression)
    before = hashlib.sha256(expression.read_bytes()).hexdigest()

    result = run_cli(
        root,
        "run_tumor_normal_workflow.py",
        expression,
        "--output-dir",
        output_dir,
    )

    assert_clean_failure(result)
    assert "refusing to overwrite expression input" in result.stderr
    assert hashlib.sha256(expression.read_bytes()).hexdigest() == before


def test_duplicate_gene_header_is_cleanly_rejected(tmp_path, root):
    expression = tmp_path / "duplicate_genes.csv"
    expression.write_text(
        "sample,ENSG00000000001.1,ENSG00000000001.1\ns1,1,2\n",
        encoding="utf-8",
    )
    output = tmp_path / "scores.csv"

    result = run_cli(root, "score_tumor_normal.py", expression, "-o", output)

    assert_clean_failure(result)
    assert "duplicate gene columns" in result.stderr
    assert not output.exists()


def test_malformed_qc_reference_and_gene_metadata_are_clean_failures(tmp_path, root):
    bad_reference = tmp_path / "bad_reference.json"
    bad_reference.write_text(
        '{"rules": {"min_match_rate_fail": 2}}', encoding="utf-8"
    )
    qc_output = tmp_path / "qc.json"
    qc_result = run_cli(
        root,
        "inspect_expression_input.py",
        "example_input.csv",
        "--qc-reference",
        bad_reference,
        "-o",
        qc_output,
    )
    assert_clean_failure(qc_result)
    assert "must be between 0 and 1" in qc_result.stderr
    assert not qc_output.exists()

    bad_metadata = tmp_path / "bad_metadata.csv"
    bad_metadata.write_text("wrong,column\na,b\n", encoding="utf-8")
    explanation_output = tmp_path / "explanations.csv"
    explanation_result = run_cli(
        root,
        "explain_scores.py",
        "example_input.csv",
        "--gene-metadata",
        bad_metadata,
        "-o",
        explanation_output,
    )
    assert_clean_failure(explanation_result)
    assert "missing required columns" in explanation_result.stderr
    assert not explanation_output.exists()


def test_unknown_qc_rule_is_a_clean_failure(tmp_path, root):
    bad_reference = tmp_path / "typo_reference.json"
    bad_reference.write_text(
        '{"rules": {"min_match_rate_fali": 0.99}}', encoding="utf-8"
    )
    output = tmp_path / "qc.json"

    result = run_cli(
        root,
        "inspect_expression_input.py",
        "example_input.csv",
        "--qc-reference",
        bad_reference,
        "-o",
        output,
    )

    assert_clean_failure(result)
    assert "unknown rule keys: min_match_rate_fali" in result.stderr
    assert not output.exists()


def test_empty_expression_input_is_cleanly_rejected(tmp_path, root):
    empty = tmp_path / "empty.csv"
    empty.write_text("sample,ENSG00000000001.1\n", encoding="utf-8")
    output = tmp_path / "scores.csv"

    result = run_cli(root, "score_tumor_normal.py", empty, "-o", output)

    assert_clean_failure(result)
    assert "at least one sample" in result.stderr
    assert not output.exists()


def test_atomic_writer_creates_nested_output_parent(tmp_path, root):
    output = tmp_path / "new" / "nested" / "scores.csv"
    result = run_cli(root, "score_tumor_normal.py", "example_input.csv", "-o", output)

    assert result.returncode == 0, result.stderr
    assert output.is_file()
    assert not list(output.parent.glob(f".{output.name}.*.tmp"))
