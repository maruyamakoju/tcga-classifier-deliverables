# INDEX — a guided reading path

Release: `v1.1.2-gdc-starcounts`

This deliverable grew in layers, from a base tumor-vs-normal model to external
validation, cross-platform adaptation, and a separate cancer-type classifier.
Read it in this order for the complete story; each stop below has a one-line
description and a link.

## 1. Base tumor-vs-normal model

- [`REPORT.md`](REPORT.md) — full methods and results for the pan-cancer
  tumor-vs-normal logistic-regression classifier (held-out test AUC 0.997,
  accuracy 0.979; grouped 5-fold CV AUC 0.997±0.003).
- [`MODEL_CARD.md`](MODEL_CARD.md) — one-page fact sheet: intended use, limits,
  and validation summary.

## 2. Cross-cancer generalization (LOCO)

- [`cross-cancer-holdout/LOCO_REPORT.md`](cross-cancer-holdout/LOCO_REPORT.md) —
  leave-one-cancer-type-out test: pooled AUC 0.988, macro-mean AUC 0.994, worst
  held-out type PRAD 0.950. Shows the signature transfers to cancer types never
  seen in training.

## 3. External validation

- [`external-validation/cptac_gdc/CPTAC_EXTERNAL_VALIDATION.md`](external-validation/cptac_gdc/CPTAC_EXTERNAL_VALIDATION.md)
  — external non-TCGA CPTAC-3 cohort on the same GDC STAR-Counts pipeline
  (AUC 0.989).
- [`external-validation/tcga_toil_xena/TCGA_TOIL_PIPELINE_CHECK.md`](external-validation/tcga_toil_xena/TCGA_TOIL_PIPELINE_CHECK.md)
  — cross-platform UCSC Xena Toil/RSEM check: ranking holds (AUC 0.992) but the
  default 0.5 threshold breaks.
- [`external-validation/gtex_xena/GTEX_NORMAL_VALIDATION.md`](external-validation/gtex_xena/GTEX_NORMAL_VALIDATION.md)
  — 540 GTEx normals across 27 primary sites; false-positive rate 0.996 at the
  0.5 threshold (cross-platform deployment boundary).

## 4. Cross-platform adaptation

- [`cross-platform-adaptation/CROSS_PLATFORM_ADAPTATION.md`](cross-platform-adaptation/CROSS_PLATFORM_ADAPTATION.md)
  — label-free cohort standardization that restores Toil accuracy from 0.515 to
  0.935 without refitting the model.

## 5. Cancer-type classifier

- [`cancer-type-classifier/CANCER_TYPE_CLASSIFIER.md`](cancer-type-classifier/CANCER_TYPE_CLASSIFIER.md)
  — a separate 17-class tissue-of-origin classifier (patient-held-out accuracy
  0.930, balanced accuracy 0.878, macro-F1 0.877).

## Code & tests

- [`tcga_rnaseq/`](tcga_rnaseq/) — shared core library (I/O, gene alignment,
  scoring, metrics) reused across the scoring entry points.
- [`tests/`](tests/) — pytest suite covering core units and numerical
  reproducibility.

Scoring entry points: `score_tumor_normal.py` (tumor-vs-normal),
`cohort_adapt_score.py` (cross-platform adaptation + scoring), and
`cancer-type-classifier/predict_cancer_type.py` (tissue-of-origin).
