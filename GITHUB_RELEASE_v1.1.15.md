# v1.1.15-gdc-starcounts

Release validation malformed-input hardening release for the GDC STAR-Counts
tumor-vs-normal lightweight bundle. Model weights, training data, and headline
validation metrics are unchanged from v1.1.14.

## What changed

- `validate_release_lite.py` now reports malformed `release_manifest.json`
  content as validation errors instead of crashing during release directory or
  source-parity validation.
- `validate_release_lite.py` now reports corrupt or mislabeled zip files as
  explicit validation failures.
- `validate_zip_bundle.py` now reports corrupt or mislabeled zip files before
  extraction instead of raising raw zip exceptions.
- Added unit coverage for malformed release manifests and bad zip archives in
  both release validation paths.

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
- SHA256: `ed4696fa5c6a3728f2f50cf74db6058788a8a9853d5216213e28b0b8f5b19e2f`
- Size: `312251` bytes
- Zip entries: `72`
