# v1.1.14-gdc-starcounts

Threshold contract metric hardening release for the GDC STAR-Counts
tumor-vs-normal lightweight bundle. Model weights, training data, and headline
validation metrics are unchanged from v1.1.13.

## What changed

- `validate_output_contracts.py` now rejects non-empty threshold metric fields
  that cannot be parsed as numbers instead of silently ignoring them as missing.
- Blank optional `youden_j` values remain valid for rows where that metric is
  intentionally absent.
- Added unit coverage for non-numeric threshold metrics and blank optional
  `youden_j` handling.

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
- SHA256: `24ffae7af5db4a8c04ef47d527969f11a8208d55df2223bf42f193a71d1c4faf`
- Size: `311884` bytes
- Zip entries: `72`
