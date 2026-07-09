# Release notes

## v1.1.8-gdc-starcounts — 2026-07-09

Lite documentation audit release. Model weights, training data, and headline
validation metrics are unchanged from v1.1.7; this release tightens the
documentation integrity checks and makes the lightweight bundle's file map
clearer when read outside the full development repository.

### Added

- `audit_release_docs.py` now validates code-spanned local file and directory
  references, not only Markdown links and Python command references.
- Central allowlists for generated workflow outputs, release sidecars, and
  full-development-tree-only references, so unexpected stale paths fail the
  audit.

### Changed

- `README.md` now separates lightweight-bundle files from full development
  tree-only maintenance and historical-analysis assets.
- `REPORT.md` now marks LOCO, cross-platform adaptation benchmark,
  cancer-type classifier, and maintenance/regeneration artifacts as full
  development tree-only when they are not part of the lightweight zip.

## v1.1.7-gdc-starcounts — 2026-07-09

Core guardrail release. Model weights, training data, and headline validation
metrics are unchanged from v1.1.6; this release moves invalid matched-value
protection from only the public CLIs into the shared scoring API and associated
reproducibility scripts.

### Added

- Strict default invalid matched-value checks in `tcga_rnaseq.predict_proba()`
  and `tcga_rnaseq.score_binary_dataframe()`.
- Strict default invalid matched-value checks in the legacy low-level
  `align_to_genes()` helper; callers that intentionally reproduce historical
  mean-imputed benchmarks must pass `allow_invalid_values=True`.
- Shape and finite-value validation in `predict_proba_from_aligned()`.
- Markdown local-link validation in `audit_release_docs.py`, including
  `INDEX.md` version coverage.

### Changed

- Invalid alignment summary and validation helpers now live in the shared
  `tcga_rnaseq` core instead of only `score_tumor_normal.py`.
- External validation and cross-platform benchmark scripts use
  `deployable_lr_weights.npz` and explicit invalid-value handling.
- Release artifact metadata now advertises the strict source-parity validation
  command with `--source-root . --artifacts RELEASE_ARTIFACTS.json`.
- `INDEX.md` is release-lite-safe and no longer links to full-repo-only paths.

### Tests

- Unit coverage for strict core scoring, aligned-matrix shape/finiteness checks,
  and default rejection of invalid matched values.
- Reproducibility tests now explicitly mark historical Toil/GTEx external
  benchmarks as mean-imputed via `allow_invalid_values=True`.

## v1.1.6-gdc-starcounts — 2026-07-09

Remaining CLI guardrail update. The fitted model weights and headline
validation metrics are unchanged from v1.1.5; this release extends invalid
matched expression value checks to the remaining prediction/adaptation entry
points.

### Added

- `--max-invalid-cell-fraction` and `--allow-invalid-values` on
  `cohort_adapt_score.py`.
- Invalid matched-value reporting on
  `cancer-type-classifier/predict_cancer_type.py`.
- Adaptation output metadata for matched/missing model genes and invalid
  matched-value counts.
- Safety coverage showing `cohort_adapt_score.py` refuses invalid matched
  values unless explicitly allowed.
- Unit coverage showing the cancer-type CLI refuses invalid matched values
  unless explicitly allowed.

### Changed

- `cohort_adapt_score.py` now aligns once, validates matched values, and then
  scores from the validated aligned matrix.
- `cancer-type-classifier/predict_cancer_type.py` now stops by default before
  writing predictions if any matched model-gene values are invalid.

## v1.1.5-gdc-starcounts — 2026-07-09

Explanation guardrail update. The fitted model weights and headline validation
metrics are unchanged from v1.1.4; this release applies the invalid matched
expression value policy consistently to explanation generation as well as
scoring.

### Added

- `--max-invalid-cell-fraction` and `--allow-invalid-values` on
  `explain_scores.py`.
- Unit coverage for explanation alignment reports.
- Safety coverage showing `explain_scores.py` refuses to write
  `explanations.csv` when invalid matched expression values are present unless
  explicitly allowed.
- Hosted GitHub repository audit checks for administrator-enforced branch
  protection, linear history, conversation resolution, protected `v*` release
  tag rulesets, and Dependabot vulnerability alerts.

### Changed

- `explain_scores.py` now stops by default before writing explanations if any
  matched model-gene values are invalid. Users must fix the input or explicitly
  opt into reviewed mean imputation, matching the scorer/workflow behavior.

## v1.1.4-gdc-starcounts — 2026-07-08

Invalid-input guardrail update. The fitted model weights and headline
validation metrics are unchanged from v1.1.3; this release prevents malformed
matched expression values from being silently mean-imputed into deployable
scores.

### Added

- Alignment diagnostics for matched model-gene cells that are missing,
  non-numeric, `NaN`, or infinite before training-mean imputation.
- `--max-invalid-cell-fraction` and `--allow-invalid-values` on
  `score_tumor_normal.py` and `run_tumor_normal_workflow.py`.
- Safety coverage showing the scorer and workflow stop before writing
  `scores.csv` when invalid matched expression values are present.

### Changed

- `score_tumor_normal.py` and the one-command workflow now stop by default
  before writing scores if any matched model-gene values are invalid. Users
  must fix the input or explicitly opt into reviewed mean imputation.

## v1.1.3-gdc-starcounts — 2026-07-08

Quality hardening update. The fitted model weights and headline validation
metrics are unchanged from v1.1.2; this release tightens public CLI behavior,
input validation, calibration safety, and hosted repository auditing.

### Added

- `audit_github_repository.py`, a hosted-repository audit for public visibility,
  branch protection, required CI contexts, release asset digest/size, topics,
  and stale pip Dependabot PRs.
- Unit coverage for hosted repository audit helpers, duplicate/colliding gene
  columns, cohort-adaptation label joins, and accidental calibration subsets.
- Safety coverage for unsupported legacy pickle/RF scorer options.

### Changed

- The public `score_tumor_normal.py` CLI now exposes only the pure NumPy
  logistic-regression scorer from `deployable_lr_weights.npz`; legacy pickle/RF
  options fail clearly because those artifacts are not in the lightweight
  public release.
- `calibrate_threshold.py` and the one-command workflow now require all scored
  samples to have labels by default. Use `--min-match-fraction` only when a
  partial calibration subset is intentional.

### Fixed

- Duplicate input gene columns and Ensembl-version collisions are rejected
  before scoring/QC instead of silently choosing one column.
- `cohort_adapt_score.py` now normalizes labels with the shared label parser,
  preserves missing label matches, reports label join counts, and no longer
  maps string labels such as `"1"`/`"0"` incorrectly.
- `cancer-type-classifier/predict_cancer_type.py` now validates `--topk` bounds.

## v1.1.2-gdc-starcounts — 2026-07-08

Publication metadata refresh. The fitted model weights, scoring behavior, and
validated headline metrics are unchanged from v1.1.1; this release improves
public repository metadata, citation readiness, and hosted CI coverage.

### Added

- `.zenodo.json` and `codemeta.json` for machine-readable software metadata.
- `audit_publication_readiness.py`, a public-release audit for secret-like
  strings, oversized tracked/history blobs, line endings, and release asset
  metadata consistency.
- CI coverage for the publication readiness audit.
- README badges for CI, release, license, and citation metadata.

### Changed

- Updated hosted CI actions to current pinned releases.
- Clarified that legacy pickle/RF training artifacts are intentionally excluded
  from the public Git history and are not needed for the lightweight release.
- Added Zenodo/CodeMeta metadata to the lightweight release bundle.

## v1.1.1-gdc-starcounts — 2026-07-08

Release-quality hardening update. The fitted model weights and validated
headline metrics are unchanged from v1.1.0; this release tightens input
validation, release reproducibility, and regression coverage around the
deployment bundle.

### Added

- Source-parity validation for `release-lite/`: `validate_release_lite.py
  --source-root .` now fails when copied payload files drift from the full
  deliverables tree.
- `RELEASE_ARTIFACTS.json` validation for zip byte size, entry count, SHA256,
  release file count, and total release bytes.
- Per-step subprocess timeouts for release build, release acceptance, zip
  validation, and optional smoke validation.
- Always-on pytest assertions for shipped headline metrics, patient-disjoint
  train/test metadata, calibration validation, Youden-J tie-breaking, empty
  explanation output schema, and multiclass-weight rejection.

### Fixed

- Non-finite expression values (`NaN`, `inf`, `-inf`) are now imputed at the
  training mean during scoring rather than leaking into model probabilities.
- Tumor-vs-normal scoring, QC, and explanation paths now reject multiclass
  `.npz` weights with a clear error.
- `calibrate_threshold.py` now rejects duplicate sample IDs, non-finite or
  out-of-range probabilities, and out-of-range thresholds before computing
  metrics.
- Youden-J threshold tie-breaking is now consistent with the shared metrics
  core.
- QC-failed workflows now report that scoring was not run instead of implying
  that zero samples were scored.
- Empty explanation outputs preserve the stable CSV column contract.

### Release engineering

- Zip generation is now byte-reproducible for identical release contents.
- `validate_zip_bundle.py` can resolve a sibling/parent release zip when run
  from inside `release-lite/`.
- The rebuilt zip SHA256 is recorded in `RELEASE_ARTIFACTS.json`.

## v1.1.0-gdc-starcounts — 2026-07-06

Engineering-quality release. The tumor-vs-normal model and its validated results
are unchanged from v1.0.0; this release hardens and unifies the code around it.

### Added

- Shared core library `tcga_rnaseq/` (`io`, `align`, `score`, `metrics`). Every
  scoring CLI now shares one hardened, unit-tested implementation instead of
  re-implementing model loading, gene alignment, and scoring.
- `cohort_adapt_score.py`: label-free cross-platform domain-adaptation scorer
  (cohort standardization) is now part of the bundle.
- `INDEX.md`: a single guided reading path across the base model, LOCO, external
  validation, cross-platform adaptation, and the cancer-type classifier.
- A `pytest` regression suite (dev tree) that locks in the verified metrics
  (held-out AUC 0.997, LOCO macro 0.994, Toil-adapted acc 0.935, cancer-type
  patient-held-out acc 0.930) and the numerically stable sigmoid / gene-alignment.

### Fixed

- Silent gene-alignment failure: Ensembl IDs now match with or without the
  `.version` suffix in every scorer (previously `cohort_adapt_score.py` and the
  cancer-type scorer could NaN-impute every gene and score at the bare intercept).
- Numerically stable `sigmoid` used everywhere (no overflow on large logits).
- `train_classifier.py`: cross-validation now re-fits feature selection inside
  each fold (no selection leakage), guards the optional xgboost import, and runs
  under a `__main__` guard with anchored paths.
- Removed the stale self-nested `release-lite/release-lite` copy.

## v1.0.0-gdc-starcounts — 2026-07-03

Lightweight scoring release for GDC STAR-Counts-style tumor-vs-adjacent-normal
RNA-seq matrices.

Artifact byte size and SHA256 are written to `RELEASE_ARTIFACTS.json` during
`python build_release_lite.py --smoke`.

This release separates the deployment path from the full training artifacts.

### Added

- Pure NumPy LR scorer weights: `deployable_lr_weights.npz`.
- `EXECUTIVE_SUMMARY.md`, `VERSION`, and `RELEASE_METADATA.json` for handoff
  and versioned release metadata.
- `USER_GUIDE.md` and `templates/` for practical input preparation and QC
  interpretation.
- `TROUBLESHOOTING.md` for install, input-QC, threshold, and release-integrity
  fixes.
- `check_environment.py` for runtime/package/file diagnostics and optional
  bundled self-test.
- `audit_lightweight_dependencies.py` for guarding the lightweight runtime
  against accidental heavy dependencies.
- `audit_cli_entrypoints.py` for guarding release CLI `--help` entry points.
- `audit_release_docs.py` for documentation and release-bundle reference
  checks.
- `DATA_DICTIONARY.md` and `validate_output_contracts.py` for stable
  CSV/JSON output schema documentation and validation.
- `example_workflow_output/` as a checked-in reference output for the bundled
  one-command workflow.
- Default scorer path in `score_tumor_normal.py` no longer needs scikit-learn.
- `requirements-light.txt` for minimal scoring installs.
- `inspect_expression_input.py` and `model_qc_reference.json` for pre-scoring
  compatibility QC.
- `run_tumor_normal_workflow.py` for one-command QC, scoring, optional
  calibration, explanations, manifest, and Markdown report.
- `build_release_lite.py` to regenerate the lightweight bundle, checksums, zip,
  and release smoke test from one command.
- `validate_release_lite.py` and `release_manifest.json` for independent
  release integrity checks and forbidden-artifact detection.
- `run_safety_tests.py` for negative-path checks around QC failures, invalid
  CLI arguments, and workflow stop behavior.
- `run_release_acceptance.py` for one-command environment, smoke, safety, and
  release-integrity acceptance checks.
- `validate_zip_bundle.py` for clean zip extraction and acceptance checks.
- `calibrate_threshold.py` for labeled threshold calibration.
- `explain_scores.py` for per-sample LR logit contribution reports.
- `model_gene_metadata.csv` with all 2,000 model genes, coefficients, and
  scaling metadata.
- `MODEL_CARD.md` with intended use, limitations, and validation boundaries.
- `run_smoke_tests.py` for release sanity checks.
- `external-validation/` reports for CPTAC/GDC, GTEx/Toil, and TCGA/Toil.

### Current validated boundary

- Strong within GDC STAR-Counts-scale TCGA/CPTAC-style data.
- Not safe for direct hard calls on Toil/RSEM/GTEx/GEO-style matrices without
  refitting or calibration.
- Input QC now flags the Toil/RSEM and GTEx boundary checks as WARN while the
  bundled TCGA example and CPTAC/GDC validation matrix pass.

### Recommended smoke test

```bash
python check_environment.py --self-test
python audit_lightweight_dependencies.py
python audit_cli_entrypoints.py
python audit_release_docs.py
python validate_output_contracts.py
python run_release_acceptance.py
python validate_zip_bundle.py tcga-tumor-normal-release-lite.zip
python run_smoke_tests.py
```
