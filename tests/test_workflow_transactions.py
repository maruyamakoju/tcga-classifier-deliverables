"""Managed-artifact transaction and calibration disclosure tests."""
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

import run_tumor_normal_workflow as workflow


OPTIONAL_OUTPUTS = {"thresholds.csv", "calibration.json", "explanations.csv"}


def run_workflow(root, *arguments):
    return subprocess.run(
        [sys.executable, "run_tumor_normal_workflow.py", *map(str, arguments)],
        cwd=root,
        text=True,
        capture_output=True,
    )


def sha256_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def managed_hashes(output_dir):
    output_dir = Path(output_dir)
    return {
        name: sha256_file(output_dir / name)
        for name in workflow.MANAGED_OUTPUT_NAMES
        if (output_dir / name).is_file()
    }


def actual_managed_names(output_dir):
    output_dir = Path(output_dir)
    return {
        name for name in workflow.MANAGED_OUTPUT_NAMES if (output_dir / name).is_file()
    }


def direct_arguments(root, output_dir, *, labels=False):
    arguments = [
        os.path.join(root, "example_input.csv"),
        "--output-dir",
        str(output_dir),
        "--skip-explanations",
    ]
    if labels:
        arguments.extend(["--labels", os.path.join(root, "example_labels.csv")])
    return arguments


def assert_manifest_matches_managed_files(output_dir, expected_status):
    manifest = json.loads((Path(output_dir) / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["status"] == expected_status
    assert set(manifest["outputs"].values()) == actual_managed_names(output_dir)
    return manifest


def test_workflow_rerun_removes_stale_optional_outputs_and_failed_scores(tmp_path, root):
    output_dir = tmp_path / "workflow"
    first = run_workflow(
        root,
        os.path.join(root, "example_input.csv"),
        "--labels",
        os.path.join(root, "example_labels.csv"),
        "--output-dir",
        output_dir,
        "--top-n",
        "1",
    )
    assert first.returncode == 0, first.stderr
    assert all((output_dir / name).is_file() for name in OPTIONAL_OUTPUTS)
    first_manifest = assert_manifest_matches_managed_files(output_dir, "complete")
    assert first_manifest["input_sha256"] == sha256_file(
        os.path.join(root, "example_input.csv")
    )
    assert first_manifest["model_sha256"] == sha256_file(
        os.path.join(root, "deployable_lr_weights.npz")
    )
    assert set(first_manifest["input_artifacts"]) == {
        "expression_input",
        "gene_metadata",
        "labels",
        "lr_weights",
        "qc_reference",
    }
    for record in first_manifest["input_artifacts"].values():
        assert record["bytes"] > 0
        assert len(record["sha256"]) == 64

    second = run_workflow(
        root,
        os.path.join(root, "example_input.csv"),
        "--output-dir",
        output_dir,
        "--skip-explanations",
    )
    assert second.returncode == 0, second.stderr
    assert not any((output_dir / name).exists() for name in OPTIONAL_OUTPUTS)
    assert (output_dir / "scores.csv").is_file()
    second_manifest = json.loads((output_dir / "manifest.json").read_text(encoding="utf-8"))
    assert second_manifest["status"] == "complete"
    assert second_manifest["labels"] is None
    assert set(second_manifest["outputs"].values()) == actual_managed_names(output_dir)

    no_match = tmp_path / "no_match.csv"
    pd.DataFrame({"NOT_A_MODEL_GENE": [1.0]}, index=["s1"]).to_csv(
        no_match, index_label="sample"
    )
    third = run_workflow(root, no_match, "--output-dir", output_dir)
    assert third.returncode != 0
    third_manifest = json.loads((output_dir / "manifest.json").read_text(encoding="utf-8"))
    assert third_manifest["status"] == "stopped_after_qc_fail"
    assert not (output_dir / "scores.csv").exists()
    assert not any((output_dir / name).exists() for name in OPTIONAL_OUTPUTS)
    assert set(third_manifest["outputs"].values()) == {
        "manifest.json",
        "qc.json",
        "workflow_report.md",
    }
    assert set(third_manifest["outputs"].values()) == actual_managed_names(output_dir)


def test_workflow_labels_apparent_metrics_and_small_class_warning(tmp_path, root):
    output_dir = tmp_path / "workflow"
    result = run_workflow(
        root,
        os.path.join(root, "example_input.csv"),
        "--labels",
        os.path.join(root, "example_labels.csv"),
        "--output-dir",
        output_dir,
        "--skip-explanations",
    )

    assert result.returncode == 0, result.stderr
    calibration = json.loads((output_dir / "calibration.json").read_text(encoding="utf-8"))
    report = (output_dir / "workflow_report.md").read_text(encoding="utf-8")
    assert calibration["evaluation_type"] == "apparent_resubstitution"
    assert "not independent validation" in calibration["evaluation_note"]
    assert calibration["warnings"]
    assert "apparent/resubstitution" in report
    assert "fewer than 10" in report
    assert "apparent/resubstitution" in result.stderr


def test_staging_write_failure_preserves_existing_generation_byte_for_byte(
    tmp_path, root, monkeypatch
):
    output_dir = tmp_path / "workflow"
    arguments = direct_arguments(root, output_dir)
    assert workflow.main(arguments) == 0
    before = managed_hashes(output_dir)
    real_write_text = workflow.write_text

    def fail_report_write(path, text, encoding="utf-8"):
        if Path(path).name == "workflow_report.md":
            raise ValueError("injected staging write failure")
        return real_write_text(path, text, encoding=encoding)

    monkeypatch.setattr(workflow, "write_text", fail_report_write)
    with pytest.raises(SystemExit) as exc_info:
        workflow.main(arguments)

    assert exc_info.value.code == 2
    assert managed_hashes(output_dir) == before
    assert_manifest_matches_managed_files(output_dir, "complete")


def test_fresh_staging_failure_never_publishes_a_complete_manifest(
    tmp_path, root, monkeypatch
):
    output_dir = tmp_path / "workflow"
    arguments = direct_arguments(root, output_dir)

    def fail_report_write(path, text, encoding="utf-8"):
        raise ValueError("injected fresh staging write failure")

    monkeypatch.setattr(workflow, "write_text", fail_report_write)
    with pytest.raises(SystemExit) as exc_info:
        workflow.main(arguments)

    assert exc_info.value.code == 2
    assert not (output_dir / "manifest.json").exists()
    assert not actual_managed_names(output_dir)


def test_commit_failure_after_data_promotions_rolls_back_complete_generation(
    tmp_path, root, monkeypatch
):
    output_dir = tmp_path / "workflow"
    arguments = direct_arguments(root, output_dir, labels=True)
    assert workflow.main(arguments) == 0
    before = managed_hashes(output_dir)
    real_replace = workflow._replace_file
    promoted = []

    def fail_during_promotion(source, destination):
        source = Path(source)
        destination = Path(destination)
        is_stage_promotion = (
            ".workflow-stage." in source.parent.name
            and destination.parent == output_dir.resolve()
        )
        if is_stage_promotion:
            promoted.append(source.name)
            if source.name == "qc.json":
                raise OSError("injected commit promotion failure")
        return real_replace(source, destination)

    monkeypatch.setattr(workflow, "_replace_file", fail_during_promotion)
    with pytest.raises(SystemExit) as exc_info:
        workflow.main(arguments)

    assert exc_info.value.code == 2
    assert promoted == ["calibration.json", "qc.json"]
    assert managed_hashes(output_dir) == before
    assert_manifest_matches_managed_files(output_dir, "complete")
    assert not list(tmp_path.glob(".workflow-backup.*"))


def test_manifest_is_the_last_file_promoted(tmp_path, root, monkeypatch):
    output_dir = tmp_path / "workflow"
    arguments = direct_arguments(root, output_dir)
    real_replace = workflow._replace_file
    promotions = []

    def record_promotions(source, destination):
        source = Path(source)
        destination = Path(destination)
        if (
            ".workflow-stage." in source.parent.name
            and destination.parent == output_dir.resolve()
        ):
            promotions.append(source.name)
        return real_replace(source, destination)

    monkeypatch.setattr(workflow, "_replace_file", record_promotions)
    assert workflow.main(arguments) == 0

    assert promotions
    assert promotions[-1] == "manifest.json"
    assert "manifest.json" not in promotions[:-1]
    assert_manifest_matches_managed_files(output_dir, "complete")


def test_same_path_corrupt_rerun_preserves_detectable_prior_generation(tmp_path, root):
    expression = tmp_path / "expression.csv"
    expression.write_bytes(Path(root, "example_input.csv").read_bytes())
    output_dir = tmp_path / "workflow"
    arguments = [
        str(expression),
        "--output-dir",
        str(output_dir),
        "--skip-explanations",
    ]
    assert workflow.main(arguments) == 0
    before = managed_hashes(output_dir)
    old_manifest = json.loads(
        (output_dir / "manifest.json").read_text(encoding="utf-8")
    )

    expression.write_text("not,a,valid,expression\n", encoding="utf-8")
    with pytest.raises(SystemExit) as exc_info:
        workflow.main(arguments)

    assert exc_info.value.code == 2
    assert managed_hashes(output_dir) == before
    assert old_manifest["input_sha256"] != sha256_file(expression)


def test_designed_terminal_failures_publish_exact_output_sets(tmp_path, root):
    invalid_expression = tmp_path / "invalid_expression.csv"
    pd.DataFrame(
        {"ENSG00000000005.6": ["not-a-number"]}, index=["sample-1"]
    ).to_csv(invalid_expression, index_label="sample")
    invalid_output = tmp_path / "invalid_output"
    invalid_result = run_workflow(
        root,
        invalid_expression,
        "--output-dir",
        invalid_output,
        "--allow-qc-fail",
        "--skip-explanations",
    )
    assert invalid_result.returncode == 1, invalid_result.stderr
    invalid_manifest = assert_manifest_matches_managed_files(
        invalid_output, "stopped_after_invalid_input"
    )
    assert set(invalid_manifest["outputs"].values()) == {
        "manifest.json",
        "qc.json",
        "workflow_report.md",
    }

    labels = tmp_path / "unmatched_labels.csv"
    labels.write_text("sample,label\nnot-a-scored-sample,tumor\n", encoding="utf-8")
    calibration_output = tmp_path / "calibration_output"
    calibration_result = run_workflow(
        root,
        os.path.join(root, "example_input.csv"),
        "--labels",
        labels,
        "--output-dir",
        calibration_output,
        "--skip-explanations",
    )
    assert calibration_result.returncode == 1, calibration_result.stderr
    calibration_manifest = assert_manifest_matches_managed_files(
        calibration_output, "stopped_after_calibration_error"
    )
    assert set(calibration_manifest["outputs"].values()) == {
        "manifest.json",
        "qc.json",
        "scores.csv",
        "workflow_report.md",
    }
