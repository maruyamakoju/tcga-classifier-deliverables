# INDEX — a guided reading path

Release: `v1.1.19-gdc-starcounts`

This lightweight release focuses on the deployable tumor-vs-normal classifier,
its validation boundary, and the runnable scoring workflow. Read it in this
order for the public bundle; the full development tree also contains historical
training, leave-one-cancer-out, cross-platform benchmark, and cancer-type
classifier sources that are intentionally outside the lite bundle.

## 1. Base tumor-vs-normal model

- [`REPORT.md`](REPORT.md) — full methods and results for the pan-cancer
  tumor-vs-normal logistic-regression classifier (held-out test AUC 0.997,
  accuracy 0.979; grouped 5-fold CV AUC 0.997±0.003).
- [`MODEL_CARD.md`](MODEL_CARD.md) — one-page fact sheet: intended use, limits,
  and validation summary.

## 2. External validation and limits

- [`external-validation/cptac_gdc/CPTAC_EXTERNAL_VALIDATION.md`](external-validation/cptac_gdc/CPTAC_EXTERNAL_VALIDATION.md)
  — external non-TCGA CPTAC-3 cohort on the same GDC STAR-Counts pipeline
  (AUC 0.989).
- [`external-validation/tcga_toil_xena/TCGA_TOIL_PIPELINE_CHECK.md`](external-validation/tcga_toil_xena/TCGA_TOIL_PIPELINE_CHECK.md)
  — cross-platform UCSC Xena Toil/RSEM check: ranking holds (AUC 0.992) but the
  default 0.5 threshold breaks.
- [`external-validation/gtex_xena/GTEX_NORMAL_VALIDATION.md`](external-validation/gtex_xena/GTEX_NORMAL_VALIDATION.md)
  — 540 GTEx normals across 27 primary sites; false-positive rate 0.996 at the
  0.5 threshold (cross-platform deployment boundary).

## 3. Reproducibility and bundle contents

- [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md) — exact commands and expected checks
  for rebuilding and validating the lightweight release.
- [`RELEASE_BUNDLE.md`](RELEASE_BUNDLE.md) — file-by-file contents of the public
  bundle.
- [`DATA_DICTIONARY.md`](DATA_DICTIONARY.md) — schemas for inputs, outputs,
  manifests, and validation reports.

## Code & tests

- [`tcga_rnaseq/`](tcga_rnaseq/) — shared core library (I/O, gene alignment,
  scoring, metrics) reused across the scoring entry points.
- [`run_smoke_tests.py`](run_smoke_tests.py) — lightweight end-to-end smoke tests.
- [`run_safety_tests.py`](run_safety_tests.py) — public safety tests for malformed
  inputs and invalid matched expression values.

Scoring entry points: `score_tumor_normal.py` (tumor-vs-normal),
`cohort_adapt_score.py` (cross-platform adaptation + scoring), and
`run_tumor_normal_workflow.py` (QC, scoring, optional calibration, explanations,
and report generation).
