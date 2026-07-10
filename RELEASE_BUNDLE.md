# Release-lite bundle

This folder is the lightweight deployment bundle for the TCGA/GDC tumor-vs-normal classifier.

Release: `v1.1.19-gdc-starcounts` (`2026-07-10`)

## Contents

- `score_tumor_normal.py`: default pure NumPy LR scorer.
- `deployable_lr_weights.npz`: model genes, scaler, coefficients, and intercept.
- `EXECUTIVE_SUMMARY.md`: short handoff summary.
- `USER_GUIDE.md`: practical input preparation and QC interpretation guide.
- `DATA_DICTIONARY.md`: stable input/output column and JSON contract reference.
- `TROUBLESHOOTING.md`: install, input-QC, threshold, and release-integrity fixes.
- `VERSION` and `RELEASE_METADATA.json`: versioned release metadata.
- `LICENSE`, `NOTICE.md`, `CITATION.cff`, `.zenodo.json`, and
  `codemeta.json`: license, third-party data notice, and citation/software
  metadata.
- `run_tumor_normal_workflow.py`: one-command QC, scoring, calibration, explanations, and report.
- `check_environment.py`: runtime/package/file diagnostic and optional self-test.
- `audit_lightweight_dependencies.py`: release import and minimal dependency audit.
- `audit_cli_entrypoints.py`: release CLI `--help` and shebang audit.
- `audit_release_docs.py`: documentation and release-bundle reference audit.
- `validate_output_contracts.py`: bundled CSV/JSON output contract validator.
- `inspect_expression_input.py`: pre-scoring input compatibility QC.
- `model_qc_reference.json`: QC heuristic rules and compact reference summaries.
- `model_gene_metadata.csv`: all 2,000 model genes with coefficients and scaling metadata.
- `calibrate_threshold.py`: threshold calibration from labeled scored samples.
- `explain_scores.py`: top per-sample LR logit contribution report.
- `run_smoke_tests.py`: verifies scoring, QC, calibration, explanations, and workflow output.
- `run_safety_tests.py`: verifies invalid-input guardrails and QC-fail workflow stop behavior.
- `run_release_acceptance.py`: runs environment, smoke, safety, and release-integrity checks.
- `validate_release_lite.py`: verifies checksums, manifest, zip parity, and absence of large training artifacts.
- `validate_zip_bundle.py`: extracts the zip into a clean temp directory and runs acceptance.
- `release_manifest.json`: generated file list with byte sizes and SHA256 hashes.
- `example_input.csv`, `example_output.csv`, `example_labels.csv`: runnable demo.
- `example_workflow_output/`: reference output from the bundled workflow demo.
- `templates/`: minimal CSV sketches for expression input and labels.
- `MODEL_CARD.md`, `REPORT.md`, `REPRODUCIBILITY.md`: documentation.
- `external-validation/`: compact validation reports and summary CSVs.

The large training matrices and pickle model checkpoints are intentionally not included. Use a private full-artifact checkout if you need retraining or legacy pickle parity checks.

In the full deliverables folder, regenerate this bundle with:

```bash
python build_release_lite.py --smoke
python validate_zip_bundle.py tcga-tumor-normal-release-lite.zip
```

## Quick start

```bash
pip install -r requirements-light.txt
python check_environment.py --self-test
python audit_lightweight_dependencies.py
python audit_cli_entrypoints.py
python audit_release_docs.py
python validate_output_contracts.py
python run_release_acceptance.py
python validate_zip_bundle.py ../tcga-tumor-normal-release-lite.zip
python run_smoke_tests.py
python run_safety_tests.py
python run_tumor_normal_workflow.py example_input.csv --labels example_labels.csv
python inspect_expression_input.py example_input.csv
python score_tumor_normal.py example_input.csv
python explain_scores.py example_input.csv --top-n 5
```

## Boundary

Use this bundle for GDC STAR-Counts-style log2(TPM+1) matrices. Do not use it for direct hard calls on Toil/RSEM/GTEx/GEO-style matrices without refitting or pipeline-specific threshold calibration.
