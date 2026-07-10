"""Unit tests for output contract validation edge cases."""
import json

import validate_output_contracts as contracts


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value) + "\n", encoding="utf-8")


def write_thresholds(path, metric_value="1.0", youden_default=""):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join([
            "threshold_name,threshold,accuracy,f1,precision,recall,specificity,tn,fp,fn,tp,youden_j",
            f"default,0.5,{metric_value},1.0,1.0,1.0,1.0,2,0,0,3,{youden_default}",
            "youden_j,0.989957,1.0,1.0,1.0,1.0,1.0,2,0,0,3,1.0",
            "",
        ]),
        encoding="utf-8",
    )


def test_json_contracts_report_non_object_top_levels_without_crashing(tmp_path, monkeypatch):
    monkeypatch.setattr(contracts, "ROOT", tmp_path)
    out = tmp_path / "example_workflow_output"
    write_json(out / "qc.json", [])
    write_json(out / "manifest.json", [])
    write_json(out / "calibration.json", [])

    messages = []
    contracts.check_json_contracts(messages)

    assert [message["code"] for message in messages].count("json_top_level_not_object") == 3


def test_json_contracts_report_bad_calibration_numbers_without_crashing(tmp_path, monkeypatch):
    monkeypatch.setattr(contracts, "ROOT", tmp_path)
    out = tmp_path / "example_workflow_output"
    write_json(out / "qc.json", {"status": "PASS", "gene_match": {"matched_model_genes": 2000}})
    write_json(out / "manifest.json", {"status": "complete", "qc_status": "PASS", "outputs": {}})
    write_json(
        out / "calibration.json",
        {
            "auc": "not-a-number",
            "n": 5,
            "n_normal": 2,
            "n_tumor": 3,
            "recommended_accuracy": True,
            "recommended_metric": "youden_j",
            "recommended_recall": 0.9,
            "recommended_specificity": 0.8,
            "recommended_threshold": 1.5,
        },
    )

    messages = []
    contracts.check_json_contracts(messages)

    codes = [message["code"] for message in messages]
    assert codes.count("json_metric_out_of_range") == 3


def test_json_contracts_validate_manifest_outputs_shape(tmp_path, monkeypatch):
    monkeypatch.setattr(contracts, "ROOT", tmp_path)
    out = tmp_path / "example_workflow_output"
    write_json(out / "qc.json", {"status": "PASS", "gene_match": {"matched_model_genes": 2000}})
    write_json(out / "calibration.json", {})
    write_json(
        out / "manifest.json",
        {"status": "complete", "qc_status": "PASS", "outputs": {"scores": 123}},
    )

    messages = []
    contracts.check_json_contracts(messages)

    assert "manifest_output_path_not_string" in {message["code"] for message in messages}


def test_threshold_contracts_reject_non_numeric_metrics(tmp_path, monkeypatch):
    monkeypatch.setattr(contracts, "ROOT", tmp_path)
    write_thresholds(
        tmp_path / "example_workflow_output" / "thresholds.csv",
        metric_value="not-a-number",
    )

    messages = []
    contracts.check_thresholds(messages)

    assert "non_numeric_metric" in {message["code"] for message in messages}


def test_threshold_contracts_allow_blank_optional_youden_j(tmp_path, monkeypatch):
    monkeypatch.setattr(contracts, "ROOT", tmp_path)
    write_thresholds(tmp_path / "example_workflow_output" / "thresholds.csv")

    messages = []
    contracts.check_thresholds(messages)

    assert messages == []
