# Executive summary

## What this is

This package is a lightweight deployment release for a TCGA/GDC pan-cancer
bulk RNA-seq classifier that scores samples as **tumor** vs **normal**. The
default deployable model is logistic regression over 2,000 selected genes,
exported to `deployable_lr_weights.npz` so ordinary scoring uses only NumPy and
pandas rather than scikit-learn pickle loading.

Release version: `v1.1.5-gdc-starcounts`
Release date: `2026-07-09`

## Validated use

Use this release for expression matrices with:

- rows = samples, columns = Ensembl gene IDs
- values = **log2(TPM+1)**
- source scale compatible with **GDC STAR-Counts**
- tumor-vs-solid-tissue-normal / adjacent-normal research contrasts

The model is strongest inside the GDC STAR-Counts ecosystem. Internal TCGA
validation and an external CPTAC-3/GDC smoke validation support this boundary.

## Key results

| Setting | Result |
|---|---:|
| TCGA patient-held-out test | AUC 0.997, accuracy 0.979 |
| TCGA grouped 5-fold CV | AUC 0.997 +/- 0.003 |
| TCGA leave-one-cancer-type-out | macro AUC 0.994 |
| CPTAC-3/GDC STAR-Counts smoke test | AUC 0.989, accuracy 0.955 |
| TCGA Toil/RSEM pipeline check | AUC 0.992, but default-threshold accuracy 0.515 |
| GTEx/Toil normal-tissue check | 538/540 normals called tumor at threshold 0.5 |

## Important boundary

Do **not** use this release for direct hard calls on Toil/RSEM, GTEx, GEO, raw
counts, FPKM, microarray, single-cell, or spatial data without pipeline-specific
refitting or calibration. Cross-platform checks showed that ranking can remain
strong while probabilities and thresholds shift severely.

This is not a clinical diagnostic model and should not be used for patient
management.

## Fast path for users

From the lightweight release folder:

```bash
pip install -r requirements-light.txt
python check_environment.py --self-test
python audit_lightweight_dependencies.py
python audit_cli_entrypoints.py
python audit_release_docs.py
python validate_output_contracts.py
python run_release_acceptance.py
python validate_zip_bundle.py tcga-tumor-normal-release-lite.zip
python run_smoke_tests.py
python run_safety_tests.py
python run_tumor_normal_workflow.py example_input.csv --labels example_labels.csv
```

For a new matrix:

```bash
python run_tumor_normal_workflow.py input.csv --labels labels.csv
```

The workflow writes:

- `qc.json`
- `scores.csv`
- `thresholds.csv` and `calibration.json` when labels are supplied
- `explanations.csv`
- `manifest.json`
- `workflow_report.md`

Review `workflow_report.md` first. Treat any QC `WARN` or `FAIL` as a reason
to inspect normalization, gene IDs, and platform compatibility before using
hard calls. Invalid matched expression values now stop scoring and explanation
outputs before files are written unless explicitly allowed. Use `TROUBLESHOOTING.md`
for common install, input-QC, threshold, and release-integrity failures.

## Release integrity

Build and validate the release from the full deliverables folder with:

```bash
python build_release_lite.py --smoke
python validate_release_lite.py --release-dir release-lite --zip tcga-tumor-normal-release-lite.zip
python validate_zip_bundle.py tcga-tumor-normal-release-lite.zip
```

The build writes `release-lite/SHA256SUMS.txt`, `release-lite/release_manifest.json`,
and the sidecar `RELEASE_ARTIFACTS.json` containing the zip size and SHA256.

## Main files

- `release-lite/` and `tcga-tumor-normal-release-lite.zip`: deployment bundle
- `run_tumor_normal_workflow.py`: one-command QC, scoring, calibration, explanations, report
- `score_tumor_normal.py`: scoring CLI
- `inspect_expression_input.py`: input compatibility QC
- `calibrate_threshold.py`: threshold calibration from labels
- `explain_scores.py`: per-sample LR contribution report
- `MODEL_CARD.md`: intended use, limitations, and validation boundary
- `REPORT.md`: methods and results
- `REPRODUCIBILITY.md`: environment, tests, and release build procedure
- `USER_GUIDE.md`: practical input preparation and QC interpretation guide
- `DATA_DICTIONARY.md`: stable input/output columns and JSON contract reference
- `TROUBLESHOOTING.md`: common failure modes and fixes
- `check_environment.py`: package, required-file, and self-test diagnostic
- `audit_lightweight_dependencies.py`: lightweight runtime dependency audit
- `audit_cli_entrypoints.py`: release CLI `--help` and shebang audit
- `audit_release_docs.py`: documentation and release-bundle reference audit
- `validate_output_contracts.py`: bundled CSV/JSON output contract validator
- `run_release_acceptance.py`: end-to-end environment, smoke, safety, and release-integrity checks
- `validate_zip_bundle.py`: clean zip extraction and acceptance check
- `RELEASE_METADATA.json`, `VERSION`, `release_manifest.json`, `SHA256SUMS.txt`: release metadata and integrity files
