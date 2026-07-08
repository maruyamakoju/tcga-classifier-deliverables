# Cross-cancer-type generalization test (Leave-One-Cancer-type-Out)

This follow-up directly addresses the main caveat of the original report:
performance was measured only on the *same* 17 cancer types seen during training,
so it did not show whether the classifier generalizes to a tumor type it was
**never trained on**. Here each of the 17 cancer types is held out in turn: the
model trains on the other 16 and is tested on the held-out one.

## Method

- **Split:** Leave-One-Cancer-type-Out (LOCO). For each cancer type C, train on all
  samples of the other 16 types, test on all samples of C. 17 train/test rounds.
- **Pipeline:** identical to the original — SelectKBest ANOVA F-test (top 2,000
  genes) fit on the training types only, features standardized, then the model.
  Feature selection and scaling never see the held-out type (no leakage).
- **Model:** **Logistic regression only** (the original report's best model,
  test AUC 0.997). Random forest and XGBoost were dropped *from this test only*
  because of a process-level instability in the Claude Science conda env —
  `import xgboost` segfaults deterministically and RandomForest's parallel
  `.fit()` segfaults/bus-errors intermittently. Logistic regression's lbfgs fit
  is single-threaded and ran cleanly (a few folds needed an automatic retry after
  an unrelated segfault; each succeeded on retry). This is an environment issue,
  not a modeling one, and does not affect the conclusion — LR was already the top
  model and the two tree models were within a fraction of a percent of it.

## Result: the classifier generalizes almost perfectly to unseen cancer types

Pooled over all held-out predictions (every sample scored by a model that never
saw its cancer type):

| Metric | Within-distribution (original) | LOCO (unseen cancer type) |
|---|---|---|
| AUC (pooled) | 0.997 | **0.988** |
| AUC (macro-mean over 17 types) | ~1.000 | **0.994** |
| AUC (worst type) | 0.965 (PRAD) | **0.950 (PRAD)** |
| Accuracy (pooled, 0.5 threshold) | 0.979 | 0.947 |
| Average precision (pooled) | — | 0.994 |

The AUC barely moves. On 12 of 17 held-out types LOCO AUC is ≥ 0.998, and on 7 it
is exactly 1.000. The largest drop versus the within-distribution number is only
**0.015** (PRAD). This is strong evidence that the model relies on a **shared,
tissue-agnostic tumor-vs-normal signal** (the stromal / extracellular-matrix
remodeling signature identified in the original report) rather than memorizing
cancer-type-specific patterns.

## The one real caveat: discrimination transfers, the decision threshold does not

AUC (ranking quality) generalizes; the fixed 0.5 operating threshold — calibrated
on the *training* cancer types — does not always transfer to an unseen tissue.
Two held-out types rank almost perfectly but are mis-thresholded:

- **TCGA-PRAD** (prostate): AUC 0.950 but accuracy 0.750, recall 0.654 — many
  prostate tumors score just below 0.5 and get called "normal."
- **TCGA-LIHC** (liver): AUC 0.999 but accuracy 0.847, recall 0.770 — same effect.

So for a genuinely new cancer type you should expect near-perfect *ranking* but
should **re-calibrate the probability threshold** on a small labeled sample of
that tissue before trusting hard tumor/normal calls. Types with high AUC *and*
high accuracy (BRCA, the lung/colorectal/kidney types) need no adjustment.

## Bottom line

The original headline — a pan-cancer tumor-vs-normal classifier at ~1.0 AUC — is
**not an artifact of testing on seen cancer types**. It holds up when the entire
test cancer type is withheld from training (macro-mean AUC 0.994). The only
practical qualification is threshold calibration on a couple of tissues, not
discrimination ability.

## Files

- `loco_per_cancer_metrics.csv` — every held-out type × model, full metrics
- `loco_pooled_summary.csv` — pooled + macro AUC/accuracy/F1
- `loco_vs_within_comparison.csv` — per-type within-distribution vs LOCO, with AUC drop
- `loco_predictions.pkl` — pooled held-out labels + scores (for re-plotting / re-thresholding)
- `loco_report.html` — interactive figure + table (open in a browser)
