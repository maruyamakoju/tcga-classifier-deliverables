# v1.1.18-gdc-starcounts

Low gene coverage scoring guardrail release for the GDC STAR-Counts
tumor-vs-normal lightweight bundle. Model weights, training data, and headline
validation metrics are unchanged from v1.1.17.

## What changed

- `score_tumor_normal.py`, `explain_scores.py`, and `cohort_adapt_score.py`
  now refuse to write outputs when too few model genes match the input matrix by
  default.
- `cancer-type-classifier/predict_cancer_type.py` uses the same low
  model-gene-coverage guardrail for full-repository cancer-type predictions.
- Added shared gene-match validation helpers in `tcga_rnaseq.align`.
- Added `--min-model-gene-match-rate` and `--allow-low-gene-coverage` controls
  for reviewed override cases where missing-gene mean imputation is intentional.
- Expanded safety tests so no-model-gene inputs cannot silently produce score,
  explanation, or adapted-score outputs.

## Validation

- `python -m pytest tests/test_core_units.py -q`
- `python run_safety_tests.py`
- `python -m pytest -q -rs`
- `python audit_release_docs.py`
- `python build_release_lite.py --smoke --timeout-seconds 300`
- `python audit_publication_readiness.py`
- `python validate_release_lite.py --release-dir release-lite --zip tcga-tumor-normal-release-lite.zip --source-root . --artifacts RELEASE_ARTIFACTS.json`
- `python validate_zip_bundle.py tcga-tumor-normal-release-lite.zip`
- `python run_release_acceptance.py --timeout-seconds 300`
- `git diff --check`

## Release asset

- Asset: `tcga-tumor-normal-release-lite.zip`
- SHA256: `6f38249bb5f677444768fa598f430a1fa6c0b2666899d1f17e28d908e04db88a`
- Size: `315821` bytes
- Zip entries: `72`
