# v1.1.21-gdc-starcounts

Output sample ID contract hardening release for the GDC STAR-Counts
tumor-vs-normal lightweight bundle. Model weights, training data, and headline
validation metrics are unchanged from v1.1.20.

## What changed

- `validate_output_contracts.py` now validates scored-output and label sample
  identifiers with a shared trimmed, non-empty uniqueness contract.
- Output contract validation now reports sample IDs with leading/trailing
  whitespace.
- Output contract validation now reports duplicate sample IDs after trimming
  whitespace.
- Expanded output contract unit tests for whitespace-only, padded, and
  trim-colliding sample identifiers.

## Validation

- `python -m pytest tests/test_output_contracts.py -q`
- `python validate_output_contracts.py`
- `python build_release_lite.py --smoke --timeout-seconds 300`
- `python run_safety_tests.py`
- `python -m pytest -q -rs`
- `python audit_release_docs.py`
- `python audit_publication_readiness.py`
- `python validate_release_lite.py --release-dir release-lite --zip tcga-tumor-normal-release-lite.zip --source-root . --artifacts RELEASE_ARTIFACTS.json`
- `python validate_zip_bundle.py tcga-tumor-normal-release-lite.zip`
- `python run_release_acceptance.py --timeout-seconds 300`
- `git diff --check`

## Release asset

- Asset: `tcga-tumor-normal-release-lite.zip`
- SHA256: `54a7d583a895a729d8663c8ffe1365468884821b799d83b1c90b82a7d7a1bf86`
- Size: `317624` bytes
- Zip entries: `72`
