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
- Fixed a real `IndentationError` in `cross-platform-adaptation/run_adaptation_benchmark.py`
  (an internal, non-CI-tested analysis script) introduced by an earlier commit;
  it now runs and reproduces its previously committed benchmark numbers exactly.
- `cohort_adapt_score.py`, `calibrate_threshold.py`, and `inspect_expression_input.py`
  now call the shared `tcga_rnaseq` core (alignment, standardization, scoring,
  Youden threshold) instead of partially reimplementing it; fixed an output-ordering
  bug where a malformed `--labels` file could discard an already-successful
  `cohort_adapt_score.py` scoring run.
- `export_lr_weights.py` and `export_model_gene_metadata.py` now validate/load
  through the same shape-checked path as the shipped scorer.
- Consolidated the audit/validate/build script family's duplicated boilerplate
  (issue collection, JSON report writing, subprocess step running) into a new
  `release_tools/common.py`, shipped alongside the scripts that use it.
- Added a `ruff` lint step to CI and 41 new unit tests (`tcga_rnaseq` edge
  cases, `release_tools.common`).
- Fixed several docs-vs-code inconsistencies, including an undocumented CLI
  and a rounding error in a reported metric (0.9963 shown as "0.997" in
  README.md/MODEL_CARD.md/EXECUTIVE_SUMMARY.md/INDEX.md/REPORT.md/REPRODUCIBILITY.md).
- Independent review pass caught and fixed: `cohort_adapt_score.py` computed
  its reported metrics from the *rounded* `tumor_probability` column instead
  of the raw probability, which could disagree with the CSV's own `call`
  column for a sample whose probability rounds across the threshold; added
  regression coverage. `run_safety_tests.py`/`run_smoke_tests.py` gained a
  subprocess timeout without catching it, so a hung step would crash with a
  raw traceback instead of failing cleanly. `audit_publication_readiness.py`
  could crash instead of reporting a clean error if `VERSION` were missing.
  `explain_scores.py` still hand-rolled the standardization formula the rest
  of this pass moved into `tcga_rnaseq`. `cohort_adapt_score.py` was missing
  from the release-lite required-file lists despite being shipped/documented.
- Moved `validate_threshold`/`normalize_label`/sample-key helpers (previously
  cross-imported from `calibrate_threshold.py` by five other CLIs, or
  duplicated between it and `cohort_adapt_score.py`) into a new
  `tcga_rnaseq/validation.py`.
- Fixed a real crash in `cancer-type-classifier/predict_cancer_type.py`: a
  missing or malformed `--weights` file raised a raw traceback instead of a
  clean error, unlike its sibling `score_tumor_normal.py`.
- Fixed a real data-integrity bug in `external-validation/`: the GDC/Xena
  matrix cache was keyed only on a fixed output path with no fingerprint of
  sampling parameters, so a rerun with different arguments could silently
  reuse a cache built for a different cohort. Also fixed a silently-cached
  empty download, silent gene-ID collisions, non-atomic cache writes, and
  added post-merge integrity checks. Added 16 tests (mocked, no live network
  calls -- `external-validation/` is not part of the lightweight release).
- Added 61 total new unit tests this pass (86 -> 150 passing).

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
- SHA256: `dc7e3acf32f8c11b0de262094fd9d16b35b5bcb48a8f37c3988cc4f21c2f782b`
- Size: `322431` bytes
- Zip entries: `75`
