# Leave-one-cancer-type-out (LOCO) sensitivity analysis

This follow-up directly addresses the main caveat of the original report:
performance was measured only on the *same* 17 cancer types seen during training,
so it did not test project-level transport. Here each of the 17 TCGA cancer
types is held out in turn: the model trains on the other 16 and is evaluated on
the held-out project. This is a retrospective internal sensitivity analysis,
not evidence for arbitrary unseen tissues or external clinical populations.

## Method

- **Split:** Leave-One-Cancer-type-Out (LOCO). For each cancer type C, train on all
  samples of the other 16 types, test on all samples of C. 17 train/test rounds.
- **Pipeline:** identical to the original — SelectKBest ANOVA F-test (top 2,000
  genes) fit on the training types only, features standardized, then the model.
  Feature selection and scaling never see the held-out type (no leakage).
- **Model:** **Logistic regression only** (the original report's best model,
  test AUC 0.997). Random forest and XGBoost were dropped *from this test only*
  to keep the analysis aligned with the deployable model. The canonical
  `run_loco.py` reproduction uses the pinned scientific environment and fails on
  convergence rather than retrying or substituting a fit.

## Result: high retrospective LOCO discrimination within the GDC cohort

Pooled over all held-out predictions (every sample scored by a model that never
saw its cancer type):

| Metric | Within-distribution (original) | LOCO (unseen cancer type) |
|---|---|---|
| AUC (pooled) | 0.997 | **0.988** |
| AUC (macro-mean over 17 types) | ~1.000 | **0.994** |
| AUC (worst type) | 0.965 (PRAD) | **0.950 (PRAD)** |
| Accuracy (pooled, 0.5 threshold) | 0.979 | 0.947 |
| Average precision (pooled) | — | 0.994 |

The pooled values are **cross-fitted aggregates across 17 different
project-specific models**, not a global metric from one frozen model. Pooling
assumes that scores from those separately fitted models have a comparable
scale; their intercepts, selected features, and calibration can differ. Treat
the per-type AUCs and their macro summary as the primary ranking evidence, and
the pooled AUC as a secondary descriptive aggregate.

On 12 of 17 held-out types LOCO AUC is ≥ 0.998, and on 7 it is exactly 1.000 in
this dataset. The largest drop versus the within-distribution number is 0.015
(PRAD). These results are consistent with a shared signal, but LOCO cannot
separate tumor biology from GDC-wide processing, adjacent-normal sampling,
procurement, center, purity, or project-specific batch effects. It therefore
does not establish a tissue-agnostic causal mechanism.

## Threshold behavior and broader validity limits

Within this LOCO analysis, ranking is stronger than fixed-threshold accuracy.
The 0.5 operating threshold, selected in the original setting, does not always
transfer to a held-out TCGA project. Two examples are:

- **TCGA-PRAD** (prostate): AUC 0.950 but accuracy 0.750, recall 0.654 — many
  prostate tumors score just below 0.5 and get called "normal."
- **TCGA-LIHC** (liver): AUC 0.999 but accuracy 0.847, recall 0.770 — same effect.

Do not extrapolate these values into an expectation of near-perfect ranking for
a genuinely new cancer type. Any new tissue, center, protocol, or population
needs independent validation. A threshold estimated on a small labeled sample
is exploratory apparent/resubstitution calibration until confirmed on separate
data; no listed type should be treated as requiring “no adjustment” by default.

## Bottom line

Withholding one TCGA project at a time produced macro-mean AUC 0.994; the
cross-fitted pooled aggregate was 0.988. The macro/per-type results are the more
directly interpretable LOCO ranking summary because the pooled value combines
scores from separately fitted models. This reduces one specific concern about
direct project-label reuse, but does not rule out shared GDC/procurement
confounding or demonstrate external clinical transport. Threshold calibration
is not the only qualification.

## Files

- `loco_per_cancer_metrics.csv` — every held-out type × model, full metrics
- `loco_pooled_summary.csv` — pooled + macro AUC/accuracy/F1
- `loco_vs_within_comparison.csv` — per-type within-distribution vs LOCO, with AUC drop
- `run_loco.py` — canonical float64 reproduction and strict artifact verifier
- `loco_oof_predictions.csv` — optional regenerated sample-level scores (not a
  committed headline artifact)
- `loco_report.html` — interactive figure + table (open in a browser)
