# TCGA cancer-type (tissue-of-origin) classifier

**Date:** 2026-07-06
**Task:** given a tumor's bulk RNA-seq profile, predict which of 17 TCGA cancer
types (tissue of origin) it is.
**Status:** new multi-class extension of the tumor-vs-normal release.

## Data

- 1,440 TCGA/GDC tumor samples across 17 cancer types (1,438 unique patients).
- Features: 14,850 genes, GDC STAR-Counts `log2(TPM+1)` (`X_full_filtered.pkl`).
- Labels: `project` from `selected_files.csv`; patient = `case_id`.
- Class sizes range from BRCA (226) to CHOL (18), READ (20), ESCA (26).

## Model

`StandardScaler -> SelectKBest(f_classif, k=1000) -> multinomial LogisticRegression(C=2)`.

Multinomial logistic regression was chosen over gradient boosting: on the same
grouped out-of-fold protocol a `HistGradientBoostingClassifier` scored essentially
the same accuracy (0.932) but **lower balanced accuracy (0.865 vs 0.878)** and ran
~50x slower (102 s vs 2 s). Logistic regression is also directly interpretable as
per-type gene markers.

## Evaluation (patient-grouped out-of-fold)

5-fold `StratifiedGroupKFold` grouped by `case_id` (no patient in both train and
test), with out-of-fold predictions pooled. These estimate internal performance
for this fixed TCGA cohort while preventing the same `case_id` from crossing a
fold. They are not an external or prospective generalization estimate; model
and hyperparameter choices were developed in the same overall research cohort.

| Metric | Value |
|---|---:|
| Accuracy | **0.930** |
| Balanced accuracy | 0.878 |
| Macro F1 | 0.877 |
| Weighted F1 | 0.928 |

### Per-type F1 (patient-held-out)

| Cancer type | n | Precision | Recall | F1 |
|---|---:|---:|---:|---:|
| THCA (thyroid) | 118 | 1.00 | 1.00 | **1.00** |
| PRAD (prostate) | 104 | 1.00 | 1.00 | **1.00** |
| BRCA (breast) | 226 | 0.99 | 0.98 | 0.98 |
| UCEC (endometrium) | 70 | 0.97 | 0.97 | 0.97 |
| LUAD (lung adeno) | 118 | 0.94 | 0.96 | 0.95 |
| LIHC (liver) | 100 | 0.94 | 0.95 | 0.95 |
| KIRC (kidney clear-cell) | 144 | 0.96 | 0.92 | 0.94 |
| KIRP (kidney papillary) | 64 | 0.91 | 0.95 | 0.93 |
| HNSC (head & neck) | 88 | 0.90 | 0.91 | 0.90 |
| KICH (kidney chromophobe) | 50 | 0.88 | 0.92 | 0.90 |
| STAD (stomach) | 72 | 0.87 | 0.92 | 0.89 |
| BLCA (bladder) | 38 | 0.89 | 0.87 | 0.88 |
| LUSC (lung squamous) | 102 | 0.88 | 0.86 | 0.87 |
| COAD (colon) | 82 | 0.83 | 0.89 | 0.86 |
| CHOL (bile duct) | 18 | 0.75 | 0.83 | 0.79 |
| ESCA (esophagus) | 26 | 0.79 | 0.73 | 0.76 |
| READ (rectum) | 20 | 0.45 | 0.25 | **0.32** |

## Error-pattern hypotheses

Many misclassifications are between anatomically or developmentally related
tissues (dominant off-diagonal confusions, true -> predicted). This pattern is
biologically plausible but does not prove why the model made an error:

- **READ -> COAD: 15 of 20.** Rectal and colon adenocarcinoma are essentially one
  disease (colorectal); TCGA labels them separately but they are transcriptionally
  near-identical, so READ is the one genuinely hard class.
- **ESCA <-> STAD** (6 / 5): esophageal and gastric, the gastro-esophageal junction.
- **KIRC <-> KIRP <-> KICH** (5 / 5 / ...): the three kidney subtypes.
- **LUAD <-> LUSC**, **LUSC <-> HNSC** (7 / 4 / 6): lung and squamous histologies.
- **CHOL -> LIHC** (2 of 18): bile-duct vs hepatocellular, both hepatobiliary.

Tissues with unique, strongly expressed markers (thyroid, prostate) are classified
perfectly.

## Marker-gene consistency

Several top positive per-type coefficients agree with canonical tissue markers
(see `cancer_type_top_genes.csv`). This is a qualitative consistency check, not
causal interpretation or proof that predictions are free of technical confounding:

| Type | Top markers |
|---|---|
| THCA | TG, TPO, IYD, FOXE1 (thyroid hormone synthesis) |
| PRAD | ACP3, NKX3-1, KLK4, OR51E2 (prostate) |
| LIHC | SHBG, AMBP, GC, LBP (liver-secreted) |
| KIRC | FXYD2, CLCNKB, BHMT (kidney transport) |
| LUSC | SFTPB, SFTPA1/A2 (lung surfactant) |
| BRCA | PRLR, TRPS1, AZGP1 (breast) |

## Limitations

- Trained only on GDC STAR-Counts `log2(TPM+1)`; the same cross-platform threshold/
  scale caveats as the tumor-vs-normal release apply (see `../cross-platform-adaptation/`).
- 17 TCGA types only; a sample from a tissue outside this set is forced into one of
  the 17 (no out-of-distribution / "unknown" option).
- Small classes (READ, ESCA, CHOL) are the least reliable; READ is not separable
  from COAD at the transcriptome level.
- Research use only; not a clinical diagnostic.

## Usage

```bash
# score a matrix (rows=samples, cols=Ensembl gene IDs, values=log2(TPM+1))
python predict_cancer_type.py input.csv --topk 3 --out predictions.csv
```

Output columns: `sample`, `predicted_type`, `probability`, `top1..topk`
(each `TYPE:prob`). The pure-numpy model reproduces the scikit-learn pipeline
exactly (argmax agreement 1.0, max |Δp| = 1.5e-8). Requires only numpy + pandas.

## Reproduce

```bash
# 1) one-time trusted-pickle migration in the exact, separate converter
#    environment from training_tools/requirements-feature-export.txt
python cancer-type-classifier/export_features_npy.py \
  --source X_full_filtered.pkl \
  --output-dir cancer-type-classifier \
  --dtype both \
  --trusted-source-pickle

# 2) switch to the exact requirements-training.txt environment, then write a
#    fresh candidate generation without replacing the committed artifacts
python cancer-type-classifier/train_cancer_type_classifier.py \
  --features cancer-type-classifier/X_full.npy \
  --metadata selected_files.csv \
  --output-dir <fresh-output-dir> \
  --gene-symbols cancer-type-classifier/gene_id_to_name.csv \
  --verify-shipped cancer-type-classifier/cancer_type_lr_weights.npz
```

The float32 input is intentional for exact reproduction of this historical
17-class pipeline. The committed path-neutral feature-export lock verifies the
source/output hashes, shapes, dtypes, and exact converter versions without
committing the large arrays. See the top-level `REPRODUCIBILITY.md` for the
isolated converter setup, manifest, overwrite, and numerical-drift rules.

## Files

- `train_cancer_type_classifier.py` - training, patient-grouped OOF evaluation, exports
- `predict_cancer_type.py` - pure-numpy scoring CLI
- `export_features_npy.py` - locked, version-neutral trusted-pickle migration
- `cancer_type_lr_weights.npz` - deployable model (scaler + 17x1000 coef + intercept)
- `cancer_type_per_class_metrics.csv`, `cancer_type_confusion_matrix.csv`
- `cancer_type_oof_predictions.csv`, `cancer_type_top_genes.csv`, `cancer_type_summary.json`
- `cancer_type_classifier.html` - visual summary
- `gene_id_to_name.csv` - Ensembl ID -> symbol map (version-neutral)
