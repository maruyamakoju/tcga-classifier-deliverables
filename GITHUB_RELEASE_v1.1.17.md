# v1.1.17-gdc-starcounts

Score consistency contract hardening release for the GDC STAR-Counts
tumor-vs-normal lightweight bundle. Model weights, training data, and headline
validation metrics are unchanged from v1.1.16.

## What changed

- `validate_output_contracts.py` now validates score CSV consistency
  preconditions before comparing example probabilities or calls.
- Missing required score columns now report explicit contract errors instead of
  surfacing `KeyError` during consistency validation.
- Non-numeric or out-of-range `tumor_probability` values now fail with targeted
  diagnostics before probability deltas are computed.
- Row-count and sample-order mismatches between `example_output.csv` and
  `example_workflow_output/scores.csv` now stop comparison with clear contract
  errors, avoiding misleading downstream probability or call mismatch reports.
- Added unit coverage for malformed score consistency inputs.

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
- SHA256: `648839500c4a79ecb76c9222f6d079d5f777ddf99bfb517235d4c3a042c12352`
- Size: `313819` bytes
- Zip entries: `72`
