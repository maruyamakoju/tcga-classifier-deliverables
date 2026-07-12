"""Offline end-to-end checks for locked external-validation cohorts."""
import os
import sys

import numpy as np
import pandas as pd
import pytest


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXTERNAL_VALIDATION = os.path.join(ROOT, "external-validation")
sys.path.insert(0, EXTERNAL_VALIDATION)

import provenance as pv  # noqa: E402
import validate_cptac_gdc as cg  # noqa: E402
import validate_gtex_xena as gx  # noqa: E402
import validate_tcga_toil_xena as tx  # noqa: E402


GDC_NORMAL_ID = "10000000-0000-4000-8000-000000000001"
GDC_TUMOR_ID = "10000000-0000-4000-8000-000000000002"


def _write_model(path):
    np.savez_compressed(
        path,
        selected_genes=np.array(["ENSG1.1", "ENSG2.1"]),
        scaler_mean=np.array([0.0, 0.0]),
        scaler_scale=np.array([1.0, 1.0]),
        coef=np.array([1.0, -0.5]),
        intercept=np.array(0.1),
        class_order=np.array([0, 1]),
    )


def _never_network():
    raise AssertionError("valid locked-cohort caches must not access the network")


def test_gtex_locked_cohort_offline_run_writes_provenance(tmp_path, monkeypatch):
    model_path = tmp_path / "model.npz"
    _write_model(model_path)
    locked_path = tmp_path / "locked_gtex.csv"
    samples = ["GTEX-A-0001", "GTEX-B-0001"]
    pd.DataFrame({
        "sample": samples,
        "detailed_category": ["Tissue A", "Tissue B"],
        "primary disease or tissue": ["Tissue A", "Tissue B"],
        "_primary_site": ["Site A", "Site B"],
        "_sample_type": ["Normal Tissue", "Normal Tissue"],
        "_gender": ["Female", "Male"],
        "_study": ["GTEX", "GTEX"],
    }).to_csv(locked_path, index=False)
    out_dir = tmp_path / "gtex-output"
    out_dir.mkdir()
    matrix_path = out_dir / "gtex_selected_genes_model_scale.parquet"
    genes = ["ENSG1.1", "ENSG2.1"]
    source = {
        "dataset_id": gx.GTEX_DATASET_ID,
        "url": gx.GTEX_TPM_URL,
        "revision": "fixture-v1",
    }
    fingerprint = gx.cache_fingerprint(
        samples, genes, gx.GTEX_TPM_URL, pv.sha256_file(model_path), source
    )
    pv.write_dataframe_cache(
        matrix_path,
        pd.DataFrame([[0.2, 0.1], [0.8, 0.4]], index=samples, columns=genes),
        fingerprint=fingerprint,
        fingerprint_inputs=gx._xena_cache_inputs(
            samples, genes, gx.GTEX_TPM_URL, pv.sha256_file(model_path), source
        ),
        cache_kind="test_xena_matrix",
    )
    monkeypatch.setattr(gx, "_http_client", _never_network)
    monkeypatch.setattr(sys, "argv", [
        "validate_gtex_xena.py", "--sample-manifest", str(locked_path),
        "--source-revision", "fixture-v1", "--offline", "--out-dir", str(out_dir),
        "--weights", str(model_path),
    ])

    assert gx.main() == 0

    predictions = pd.read_csv(out_dir / "gtex_predictions.csv")
    expected = 1.0 / (1.0 + np.exp(-np.array([0.25, 0.7])))
    np.testing.assert_allclose(predictions["tumor_probability"], expected, rtol=0, atol=1e-15)
    manifest = __import__("json").loads(
        (out_dir / "run_manifest.json").read_text(encoding="utf-8")
    )
    assert manifest["parameters"]["cohort_mode"] == "locked_manifest"
    assert manifest["model"]["sha256"] == pv.sha256_file(model_path)
    assert manifest["alignment"]["n_missing_genes"] == 0
    assert manifest["source_code"]["validator"]["sha256"] == pv.sha256_file(
        gx.__file__
    )
    for item in manifest["outputs"].values():
        assert item["sha256"] == pv.sha256_file(item["path"])


def test_gtex_derived_output_failure_preserves_previous_complete_generation(
    tmp_path, monkeypatch
):
    model_path = tmp_path / "model.npz"
    _write_model(model_path)
    locked_path = tmp_path / "locked_gtex.csv"
    samples = ["GTEX-A-0001", "GTEX-B-0001"]
    pd.DataFrame({
        "sample": samples,
        "primary disease or tissue": ["Tissue A", "Tissue B"],
        "_primary_site": ["Site A", "Site B"],
        "_sample_type": ["Normal Tissue", "Normal Tissue"],
        "_study": ["GTEX", "GTEX"],
    }).to_csv(locked_path, index=False)
    out_dir = tmp_path / "gtex-output"
    out_dir.mkdir()
    genes = ["ENSG1.1", "ENSG2.1"]
    source = {
        "dataset_id": gx.GTEX_DATASET_ID,
        "url": gx.GTEX_TPM_URL,
        "revision": "fixture-v1",
    }
    inputs = gx._xena_cache_inputs(
        samples, genes, gx.GTEX_TPM_URL, pv.sha256_file(model_path), source
    )
    pv.write_dataframe_cache(
        out_dir / "gtex_selected_genes_model_scale.parquet",
        pd.DataFrame([[0.2, 0.1], [0.8, 0.4]], index=samples, columns=genes),
        fingerprint=pv.semantic_fingerprint(inputs),
        fingerprint_inputs=inputs,
        cache_kind="test_xena_matrix",
    )
    old_files = {
        "gtex_predictions.csv": "old-predictions\n",
        "gtex_summary.csv": "old-summary\n",
        "run_manifest.json": '{"old": true}\n',
    }
    for name, text in old_files.items():
        (out_dir / name).write_text(text, encoding="utf-8")
    monkeypatch.setattr(gx, "_http_client", _never_network)
    monkeypatch.setattr(
        gx, "write_report",
        lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError("injected report failure")),
    )
    monkeypatch.setattr(sys, "argv", [
        "validate_gtex_xena.py", "--sample-manifest", str(locked_path),
        "--source-revision", "fixture-v1", "--offline", "--out-dir", str(out_dir),
        "--weights", str(model_path),
    ])

    with pytest.raises(RuntimeError, match="injected report failure"):
        gx.main()

    for name, text in old_files.items():
        assert (out_dir / name).read_text(encoding="utf-8") == text
    assert not (out_dir / "sampled_gtex_manifest.csv").exists()


def test_tcga_locked_cohort_offline_run_uses_dependency_light_metrics(tmp_path, monkeypatch):
    model_path = tmp_path / "model.npz"
    _write_model(model_path)
    locked_path = tmp_path / "locked_tcga.csv"
    samples = ["TCGA-AA-0001-11", "TCGA-BB-0002-01"]
    pd.DataFrame({
        "sample": samples,
        "detailed_category": ["Normal", "Tumor"],
        "primary disease or tissue": ["Kidney", "Kidney"],
        "_primary_site": ["Kidney", "Kidney"],
        "_sample_type": ["Solid Tissue Normal", "Primary Tumor"],
        "_gender": ["Female", "Male"],
        "_study": ["TCGA", "TCGA"],
        "label": [0, 1],
    }).to_csv(locked_path, index=False)
    out_dir = tmp_path / "tcga-output"
    out_dir.mkdir()
    matrix_path = out_dir / "tcga_toil_selected_genes_model_scale.parquet"
    genes = ["ENSG1.1", "ENSG2.1"]
    source = {
        "dataset_id": tx.TCGA_DATASET_ID,
        "url": tx.TCGA_TPM_URL,
        "revision": "fixture-v1",
    }
    fingerprint = gx.cache_fingerprint(
        samples, genes, tx.TCGA_TPM_URL, pv.sha256_file(model_path), source
    )
    pv.write_dataframe_cache(
        matrix_path,
        pd.DataFrame([[0.0, 1.0], [2.0, 0.0]], index=samples, columns=genes),
        fingerprint=fingerprint,
        fingerprint_inputs=gx._xena_cache_inputs(
            samples, genes, tx.TCGA_TPM_URL, pv.sha256_file(model_path), source
        ),
        cache_kind="test_xena_matrix",
    )
    monkeypatch.setattr(gx, "_http_client", _never_network)
    monkeypatch.setattr(sys, "argv", [
        "validate_tcga_toil_xena.py", "--sample-manifest", str(locked_path),
        "--source-revision", "fixture-v1", "--offline", "--out-dir", str(out_dir),
        "--weights", str(model_path),
    ])

    assert tx.main() == 0

    summary = pd.read_csv(out_dir / "tcga_toil_summary.csv").iloc[0]
    assert summary["auc"] == 1.0
    assert summary["average_precision"] == 1.0
    assert (out_dir / "run_manifest.json").is_file()


def test_cptac_locked_cohort_offline_run_reuses_expression_cache(tmp_path, monkeypatch):
    model_path = tmp_path / "model.npz"
    _write_model(model_path)
    locked_path = tmp_path / "locked_cptac.csv"
    locked = pd.DataFrame({
        "file_id": [GDC_NORMAL_ID, GDC_TUMOR_ID],
        "file_name": ["normal.tsv", "tumor.tsv"],
        "file_size": [100, 200],
        "md5sum": ["0" * 32, "1" * 32],
        "project": ["CPTAC-3", "CPTAC-3"],
        "case_submitter_id": ["case-normal", "case-tumor"],
        "sample_submitter_id": ["sample-normal", "sample-tumor"],
        "sample_type": ["Solid Tissue Normal", "Primary Tumor"],
        "n_sample_types": [1, 1],
        "label": [0, 1],
    })
    locked.to_csv(locked_path, index=False)
    sampled = cg.validate_locked_sample_manifest(locked_path)
    out_dir = tmp_path / "cptac-output"
    out_dir.mkdir()
    matrix_path = out_dir / "expression_selected_genes.parquet"
    genes = ["ENSG1.1", "ENSG2.1"]
    inputs = cg._expression_cache_inputs(
        sampled, genes, pv.sha256_file(model_path), "fixture-v1"
    )
    pv.write_dataframe_cache(
        matrix_path,
        pd.DataFrame(
            [[0.0, 1.0], [2.0, 0.0]],
            index=[GDC_NORMAL_ID, GDC_TUMOR_ID],
            columns=genes,
        ),
        fingerprint=pv.semantic_fingerprint(inputs),
        fingerprint_inputs=inputs,
        cache_kind="test_cptac_matrix",
    )
    monkeypatch.setattr(cg, "_http_client", _never_network)
    monkeypatch.setattr(sys, "argv", [
        "validate_cptac_gdc.py", "--sample-manifest", str(locked_path),
        "--source-revision", "fixture-v1", "--offline", "--out-dir", str(out_dir),
        "--weights", str(model_path),
    ])

    assert cg.main() == 0

    summary = pd.read_csv(out_dir / "cptac_summary.csv").iloc[0]
    assert summary["auc"] == 1.0
    assert summary["average_precision"] == 1.0
    assert (out_dir / "run_manifest.json").is_file()
