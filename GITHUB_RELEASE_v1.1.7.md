# v1.1.7-gdc-starcounts

Core guardrail release for the GDC STAR-Counts tumor-vs-normal lightweight bundle.
Model weights, training data, and headline validation metrics are unchanged from
v1.1.6.

## What changed

- `tcga_rnaseq.predict_proba()` and `score_binary_dataframe()` now reject invalid
  matched model-gene values by default unless explicitly allowed.
- The legacy low-level `align_to_genes()` helper is strict by default because it
  cannot return invalid-value reports.
- `predict_proba_from_aligned()` now validates aligned matrix shape and finite
  values before scoring.
- Invalid alignment helpers now live in the shared `tcga_rnaseq` core.
- External validation and cross-platform benchmark scripts use the deployable
  `.npz` model and explicit invalid-value handling.
- `audit_release_docs.py` now checks `INDEX.md` version drift and broken local
  Markdown links.
- Release artifact metadata now advertises the strict source-parity validation
  command.

## Validation

- `python -m pytest -q -rs`
- `python run_safety_tests.py`
- `python audit_cli_entrypoints.py`
- `python audit_release_docs.py`
- `python validate_output_contracts.py`
- `python build_release_lite.py --smoke --timeout-seconds 300`
- `python audit_publication_readiness.py`
- `python validate_release_lite.py --release-dir release-lite --zip tcga-tumor-normal-release-lite.zip --source-root . --artifacts RELEASE_ARTIFACTS.json`
- `python run_release_acceptance.py --timeout-seconds 300`
- `python validate_zip_bundle.py tcga-tumor-normal-release-lite.zip`

## Release asset

- Asset: `tcga-tumor-normal-release-lite.zip`
- SHA256: `c5effeb237bbcb84e7163bb51ff504ab7205a414afd060439041159edcc32c33`
- Size: `307291` bytes
- Zip entries: `72`
