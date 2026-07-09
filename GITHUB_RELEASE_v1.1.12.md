# v1.1.12-gdc-starcounts

Duplicate zip member validation release for the GDC STAR-Counts tumor-vs-normal
lightweight bundle. Model weights, training data, and headline validation
metrics are unchanged from v1.1.11.

## What changed

- `validate_release_lite.py` now rejects duplicate file entries in the release
  zip instead of silently collapsing names into a set during parity checks.
- `validate_zip_bundle.py` now rejects duplicate member paths before extraction,
  preventing overwrite ambiguity in malformed archives.
- Added unit coverage for duplicate zip member detection in both validation
  paths.

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
- SHA256: `bb7bac96a6b020d3398ae20be797d3f36ca02d38341758b54c99da5c79f1ce97`
- Size: `310907` bytes
- Zip entries: `72`
