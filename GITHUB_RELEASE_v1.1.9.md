# v1.1.9-gdc-starcounts

Quiet documentation audit release for the GDC STAR-Counts tumor-vs-normal
lightweight bundle. Model weights, training data, and headline validation
metrics are unchanged from v1.1.8.

## What changed

- `audit_release_docs.py` now hides informational messages during normal CLI
  runs and prints a compact hidden-info count instead.
- `audit_release_docs.py --show-info` restores the full informational listing
  when debugging full-development-tree-only or release-sidecar references.
- Added unit coverage for audit output formatting so errors remain visible and
  INFO messages remain available on demand.

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
- SHA256: `82ce63eb4c2b4ba3587f1cc07f79c4390e8e5db54336d642a3bcde8e0643b460`
- Size: `308891` bytes
- Zip entries: `72`
