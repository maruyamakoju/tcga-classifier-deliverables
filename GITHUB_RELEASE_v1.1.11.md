# v1.1.11-gdc-starcounts

Timeout diagnostics hardening release for the GDC STAR-Counts tumor-vs-normal
lightweight bundle. Model weights, training data, and headline validation
metrics are unchanged from v1.1.10.

## What changed

- `run_release_acceptance.py` and `validate_zip_bundle.py` now normalize
  timeout stdout/stderr safely when Python supplies captured timeout output as
  bytes.
- Timeout stderr reports now append the timeout summary without inserting an
  extra blank line after existing stderr output.
- Added unit coverage for timeout reporting helpers in release acceptance and
  zip-bundle validation.

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
- SHA256: `5c1495e06e19e1fca97734b728a2f6eb657496be3006886b42c3b200eb87c50d`
- Size: `310555` bytes
- Zip entries: `72`
