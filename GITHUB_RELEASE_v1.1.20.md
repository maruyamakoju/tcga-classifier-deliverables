# v1.1.20-gdc-starcounts

Workflow calibration failure handling release for the GDC STAR-Counts
tumor-vs-normal lightweight bundle. Model weights, training data, and headline
validation metrics are unchanged from v1.1.19.

## What changed

- `run_tumor_normal_workflow.py` now records
  `stopped_after_calibration_error` in `manifest.json` when label calibration
  fails after scoring.
- Valid `qc.json` and `scores.csv` are preserved on calibration failure.
- Downstream `thresholds.csv`, `calibration.json`, and `explanations.csv` are
  not written after calibration failure.
- `workflow_report.md` records the calibration failure message.
- Expanded release safety tests for bad workflow calibration labels.

## Validation

- `python run_safety_tests.py`
- `python -m pytest tests/test_core_units.py -q`
- `python build_release_lite.py --smoke --timeout-seconds 300`
- `python -m pytest -q -rs`
- `python audit_release_docs.py`
- `python audit_publication_readiness.py`
- `python validate_release_lite.py --release-dir release-lite --zip tcga-tumor-normal-release-lite.zip --source-root . --artifacts RELEASE_ARTIFACTS.json`
- `python validate_zip_bundle.py tcga-tumor-normal-release-lite.zip`
- `python run_release_acceptance.py --timeout-seconds 300`
- `git diff --check`

## Release asset

- Asset: `tcga-tumor-normal-release-lite.zip`
- SHA256: `c86802ee092d7a265862107ec49d57d69907d71eee76c0b15908b214452afecf`
- Size: `317053` bytes
- Zip entries: `72`
