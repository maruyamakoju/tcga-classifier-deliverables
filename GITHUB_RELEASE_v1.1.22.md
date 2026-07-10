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
  and a rounding error in a reported metric (0.9963 shown as "0.997").

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
- SHA256: `ad0a7a66a5e6a91c4e454f9c56bef42f7905fa14f6aa57b16093199bf0bbfe89`
- Size: `321465` bytes
- Zip entries: `74`
