# v1.1.13-gdc-starcounts

Output contract JSON hardening release for the GDC STAR-Counts tumor-vs-normal
lightweight bundle. Model weights, training data, and headline validation
metrics are unchanged from v1.1.12.

## What changed

- `validate_output_contracts.py` now reports non-object JSON top-level values
  as contract errors for bundled QC, manifest, calibration, and QC-reference
  files.
- Calibration metric validation now rejects non-numeric, boolean, and
  out-of-range values with explicit diagnostics.
- Workflow manifest `outputs` validation now rejects non-object output maps and
  non-string output paths.
- Added unit coverage for malformed JSON output contracts.

## Validation

- `python audit_release_docs.py`
- `python build_release_lite.py --smoke --timeout-seconds 300`
- `python audit_publication_readiness.py`
- `python validate_release_lite.py --release-dir release-lite --zip tcga-tumor-normal-release-lite.zip --source-root . --artifacts RELEASE_ARTIFACTS.json`
- `python validate_zip_bundle.py tcga-tumor-normal-release-lite.zip`
- `python run_release_acceptance.py --timeout-seconds 300`
- `python -m pytest -q -rs`
- `git diff --check`

## Release asset

- Asset: `tcga-tumor-normal-release-lite.zip`
- SHA256: `e32c3b350eab9072c84ab6196ea287b7924a3ca8c10ba2f768ec8d72d58b3976`
- Size: `311551` bytes
- Zip entries: `72`
