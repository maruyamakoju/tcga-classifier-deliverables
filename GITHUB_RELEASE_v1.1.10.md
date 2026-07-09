# v1.1.10-gdc-starcounts

Release manifest metadata audit release for the GDC STAR-Counts tumor-vs-normal
lightweight bundle. Model weights, training data, and headline validation
metrics are unchanged from v1.1.9.

## What changed

- `validate_release_lite.py` now checks `release_manifest.json` top-level
  metadata against `VERSION`, `RELEASE_METADATA.json`, the manifest file count,
  and the validator's forbidden-artifact deny-list.
- Malformed manifest shapes, including a non-object top level or non-list
  `files` value, now report validation errors instead of risking a validator
  crash.
- Added unit coverage for valid manifest metadata, stale manifest metadata,
  non-integer manifest file counts, and malformed manifest handling.

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
- SHA256: `1faf5ce29a5e53c2dd13ae231fbb79f683b522a71370039634ea7477435dd7cf`
- Size: `310019` bytes
- Zip entries: `72`
