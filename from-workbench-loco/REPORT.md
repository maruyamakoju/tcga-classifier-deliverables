> Superseded: this is an early workbench snapshot. The canonical write-up is ../REPORT.md and ../cross-cancer-holdout/LOCO_REPORT.md.

# Pan-cancer tumor-vs-normal classifier from TCGA RNA-seq

## Data

RNA-seq gene expression (STAR gene counts, GENCODE v36) was pulled directly
from the NCI Genomic Data Commons (GDC) API for 17 TCGA cancer types that
have matched "primary tumor" and "solid tissue normal" samples available:
BRCA, KIRC, LUAD, LUSC, THCA, HNSC, LIHC, COAD, STAD, PRAD, KIRP, BLCA, UCEC,
KICH, READ, ESCA, CHOL.

- 720 solid-tissue-normal samples (all available in these 17 projects) and
  1,440 primary-tumor samples (capped at 2:1 tumor:normal per cancer type to
  limit class imbalance while keeping the tumor class in the majority, as is
  typical of TCGA cohort composition) — **2,160 samples total**.
- Expression quantified as TPM by the GDC STAR-Counts pipeline; values were
  log2(TPM+1) transformed.
- Restricted to protein-coding genes (GENCODE biotype) expressed (TPM>1) in
  >20% of samples: 14,850 genes retained out of 19,962 protein-coding genes.

## Train/test split

Samples were split with `StratifiedGroupKFold` (5-fold, one fold held out as
test), **grouped by patient** (`submitter_id`) so that no patient's tumor and
matched-normal pair could appear in both train and test, and **stratified**
by cancer type × label to preserve class balance and per-cancer-type
representation in both splits. 179 patients had a matched tumor/normal pair;
none of these are split across train/test.

- Train: 1,727 samples (576 normal / 1,151 tumor)
- Test: 433 samples (144 normal / 289 tumor), held out entirely by patient

## Feature selection & modeling

Genes were ranked by ANOVA F-test on the training set only; the top 2,000
genes were retained as features. Three classifiers were trained on this
2,000-gene training set:

- Logistic regression (L2, C=0.1, class-balanced, standardized features)
- Random forest (500 trees, max depth 8, class-balanced)
- XGBoost (300 trees, max depth 4, learning rate 0.05, scale_pos_weight
  balancing)

## Results (held-out test set, 433 samples from patients unseen in training)

| Model | AUC | Accuracy | F1 | Precision | Recall |
|---|---|---|---|---|---|
| Logistic regression | 0.997 | 0.979 | 0.984 | 0.986 | 0.983 |
| Random forest | 0.992 | 0.975 | 0.981 | 0.983 | 0.979 |
| XGBoost | 0.994 | 0.982 | 0.986 | 0.990 | 0.983 |

Logistic regression edges out on AUC; XGBoost on accuracy/F1/precision — all
three are within a fraction of a percent of each other. 5-fold grouped
cross-validation on the training set alone confirms this is not a lucky split:
logistic regression 0.997 ± 0.003 AUC, random forest 0.996 ± 0.002 AUC across
folds.

Per-cancer-type breakdown (logistic regression) shows AUC = 1.000 on 15 of 17
cancer types in the test set, with the two exceptions — TCGA-PRAD (0.965) and
TCGA-THCA (0.972) — still well above 0.96. Confusion matrix: 4 false
positives, 5 false negatives out of 433 test samples.

## Top genes driving the classification

By XGBoost importance: **PDGFRA, MMP11, ASPA, ANGPTL1, RNF150, ABCA8, HPSE2,
CDH19, LYVE1, DCN** — several of these (MMP11, a matrix metalloproteinase;
DCN/decorin; LYVE1, a lymphatic marker) are established markers of
tumor-associated stromal remodeling, consistent with prior literature on
epithelial-cancer transcriptomic signatures. By |logistic-regression
coefficient|: BMPER, GABRD, COL10A1, ACAN, MMRN1, DNASE1L3, ELF5, COMP,
HAGHL, ESM1 — again dominated by extracellular-matrix and vasculature genes
(COL10A1, COMP, MMRN1), a signature broadly reported for pan-cancer stromal
activation.

## Leave-one-cancer-type-out (LOCO) generalization

To test generalization to a cancer type never seen during training (the gap
noted in the original caveats below), the pipeline was retrained 17 times,
each time holding out one entire cancer type from training (feature
selection, scaling, and logistic regression all refit on the remaining 16
cancer types) and testing exclusively on the excluded type.

- **Mean AUC across all 17 held-out cancer types: 0.994** (range 0.950–1.000).
  Every held-out type reaches AUC ≥ 0.95 despite the model never seeing a
  single sample of that cancer type during training.
- Raw accuracy at the default 0.5 threshold is more variable (e.g., LIHC
  0.847, PRAD 0.750, UCEC 0.886) — but re-thresholding per cancer type at its
  own Youden's-J-optimal cutoff recovers accuracy to within a few points of
  the AUC (LIHC 0.847→0.993, STAD 0.907→1.000, PRAD 0.750→0.936). This shows
  the drop in raw accuracy for a few cancer types is a **threshold
  miscalibration** issue (the decision boundary learned on 16 other cancer
  types doesn't perfectly transfer), not a failure to rank tumor vs. normal
  correctly for the unseen type.
- Prostate (PRAD) is the hardest case to generalize to (AUC 0.950),
  consistent with prostate adenocarcinoma's comparatively subtle
  transcriptomic tumor/normal contrast relative to other solid tumors.

See `loco_generalization.png`, `loco_generalization.csv`, and
`loco_generalization_with_optimal_threshold.csv`.

## Caveats

- This is a **pan-cancer** tumor-vs-normal classifier, not a cancer-type
  classifier — it distinguishes malignant from non-malignant tissue across
  cancer types, largely by picking up a shared stromal/ECM remodeling
  signature. The LOCO analysis above confirms this signature generalizes to
  entirely unseen cancer types (mean AUC 0.994), addressing the original
  concern that held-out performance was measured only on cancer types also
  seen in training.
- Solid-tissue-normal samples in TCGA are usually resected margin tissue,
  not tissue from cancer-free individuals — this is the correct
  tumor-vs-adjacent-normal contrast but is a different (easier) comparison
  than tumor-vs-healthy-population.
- No batch-effect correction across cancer types or sequencing centers was
  applied beyond TPM normalization; the near-perfect per-cancer AUCs partly
  reflect the strength of a shared cancer-vs-normal signal, which is common
  and expected for this task on TCGA data (this exact contrast is
  near-saturated in the field).

## Files

- `model_performance.png` — ROC curves, confusion matrix, per-cancer-type AUC
- `feature_importance.png` — top gene importances and PDGFRA expression by group
- `loco_generalization.png` — leave-one-cancer-type-out AUC and threshold-calibration figure
- `loco_generalization.csv`, `loco_generalization_with_optimal_threshold.csv` — per-cancer-type LOCO metrics
- `test_metrics.csv` — held-out test metrics for all 3 models
- `per_cancer_type_performance.csv` — per-cancer-type AUC/accuracy breakdown
- `top_genes_xgboost.csv`, `top_genes_logreg.csv` — top 30 genes by each model
- `deployable_pipeline.pkl` — fitted selector + scaler + all 3 trained models
  + the exact 2,000-gene feature list, ready to score new TCGA-style
  log2(TPM+1) expression vectors
- `selected_files.csv` — full manifest of the 2,160 GDC file IDs / case IDs /
  cancer types / labels used
