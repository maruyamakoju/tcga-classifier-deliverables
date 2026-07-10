# v1.1.19-gdc-starcounts

Cohort label sample validation release for the GDC STAR-Counts tumor-vs-normal
lightweight bundle. Model weights, training data, and headline validation
metrics are unchanged from v1.1.18.

## What changed

- `cohort_adapt_score.py` now rejects missing or blank label sample IDs before
  metric alignment.
- Numeric string labels and unmatched input-sample reporting remain unchanged.
- Added unit coverage for missing label sample identifiers in cohort adaptation.

## Validation

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
- SHA256: `08f1a93135bc19eeb7e3968f65345491482cac75650c9eeb6b7133629f0e10e7`
- Size: `316011` bytes
- Zip entries: `72`
