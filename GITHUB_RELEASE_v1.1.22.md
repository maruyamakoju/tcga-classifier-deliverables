# v1.1.22-gdc-starcounts

Expression input read error handling release for the GDC STAR-Counts
tumor-vs-normal lightweight bundle. Model weights, training data, and headline
validation metrics are unchanged from v1.1.21.

## What changed

- `tcga_rnaseq.read_matrix()` now reports missing expression files and directory
  paths as user-facing errors.
- Common pandas/OS read failures are normalized as expression matrix read errors
  so bundled CLIs can report them without Python tracebacks.
- The scorer safety suite now verifies that missing expression inputs do not
  write partial score outputs.
- Added unit coverage for missing expression input handling.

## Validation

- `python -m pytest tests/test_core_units.py -q`
- `python run_safety_tests.py`
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
- SHA256: `86ccea066043fadb6932d98e00274a2cdd945ed2ca8f0c409db52ffa09bc5128`
- Size: `318271` bytes
- Zip entries: `72`
