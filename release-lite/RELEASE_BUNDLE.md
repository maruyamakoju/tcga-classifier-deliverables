# Release-lite bundle

This folder is the lightweight deployment bundle for the TCGA/GDC tumor-vs-normal classifier.
Lightweight scoring supports Python 3.11 or newer (CI exercises 3.11 and 3.13).
Exact model refitting is a full-tree workflow pinned to Python 3.11 and is not
part of the lite bundle.

Release: `v2.3.0-gdc-starcounts` (`2026-07-12`; public scoring-library API `3.0.0`)

## Contents

- `score_tumor_normal.py`: default pure NumPy LR scorer.
- `cohort_adapt_score.py`: scoring with adaptation disabled by default and
  explicit experimental adapted modes.
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
python validate_zip_bundle.py tcga-tumor-normal-release-lite.zip --expected-sha256 <trusted-published-sha256>
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
python validate_zip_bundle.py ../tcga-tumor-normal-release-lite.zip --expected-sha256 <trusted-published-sha256>
python run_smoke_tests.py
python run_safety_tests.py
python run_tumor_normal_workflow.py example_input.csv --labels example_labels.csv
python inspect_expression_input.py example_input.csv
python score_tumor_normal.py example_input.csv
python explain_scores.py example_input.csv --top-n 5
```

## Boundary

Use this bundle for GDC STAR-Counts-style log2(TPM+1) matrices. Do not use it
for direct hard calls on Toil/RSEM/GTEx/GEO-style matrices without refitting or
pipeline-specific threshold calibration. `tumor_probability` is a model
logistic score, not clinical risk or a calibrated diagnostic probability.

Adaptation defaults to `none`. Adapted modes are explicit, transductive and
composition-dependent opt-ins, require at least 20 samples by default, and
produce scores that cannot be compared across separately adapted batches.
Calibration metrics are same-sample apparent/resubstitution estimates. The
committed external metrics are historical; v2.0.0 adds locked cohorts and
cache/run provenance, and v2.2.0 live-reconfirmed CPTAC-3 (GDC Data Release 45.0,
AUC 0.9886 reproduced) while TCGA-Toil and GTEx remain historical snapshots.

Never extract or execute a downloaded ZIP based only on its filename. Supply a
trusted published digest with `--expected-sha256`; if no trusted digest is
available, `--skip-acceptance` performs structural inspection without
extracting or running archive content.
