# Release notes

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
