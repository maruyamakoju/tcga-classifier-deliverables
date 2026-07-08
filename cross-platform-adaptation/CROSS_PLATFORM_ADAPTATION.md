# Cross-platform calibration by label-free cohort standardization

**Date:** 2026-07-05
**Model:** frozen release `v1.0.0-gdc-starcounts` (`deployable_lr_weights.npz`)
**Status:** research add-on to the released tumor-vs-normal classifier

## Problem

The released model is a logistic regression on 2,000 genes, calibrated to **GDC
STAR-Counts `log2(TPM+1)`**. External validation established a sharp boundary:

- Discrimination (AUC) transfers across pipelines: TCGA-Toil/RSEM **AUC 0.992**,
  CPTAC-3/GDC **AUC 0.989**.
- The fixed **0.5 decision threshold does not** transfer. On the Toil/RSEM
  pipeline the deployed model calls 197/200 samples tumor (accuracy **0.515**,
  specificity **0.03**); on 540 GTEx/Toil normal tissues it calls 538 tumor
  (false-positive rate **0.996**, median probability 0.9999).

The release therefore tells users *not* to make hard calls on non-GDC pipelines
without recalibration, but shipped no recalibration recipe that does not require
labels. This document provides and quantifies one.

## Method

The failure is a per-gene **location/scale shift** between pipelines: after
`log2(TPM+1)`, Toil/RSEM gene distributions sit at systematically different
means and spreads than GDC STAR-Counts, so the training-fitted standardizer
`z = (x - train_mean)/train_scale` pushes foreign samples far into "tumor"
territory.

**Cohort standardization** removes that shift without labels and without
retraining: standardize each gene using the *input cohort's own* per-gene mean
and standard deviation before applying the frozen coefficients.

```
z_i = (x_i - cohort_mean) / cohort_std        # cohort_zscore  (recommended)
z_i = (x_i - cohort_mean) / train_scale       # cohort_center  (location only)
p   = sigmoid(coef . z_i + intercept)
```

This realigns the cohort's per-gene marginal onto the training marginal (a
location-scale / ComBat-lite alignment). It is applied at inference, needs only
the expression matrix, and preserves the ranking (AUC) while restoring the
threshold.

## Results

Frozen model, three external cohorts, three scoring modes. The numpy
reimplementation reproduces the deployed pipeline exactly
(`max |p - reported| = 4.9e-7`).

| Cohort | Mode | AUC | acc@0.5 | balanced acc | FPR@0.5 |
|---|---|---:|---:|---:|---:|
| TCGA-Toil/RSEM (foreign) | deployed | 0.992 | 0.515 | 0.515 | 0.97 |
| TCGA-Toil/RSEM (foreign) | cohort center | 0.992 | 0.930 | 0.930 | 0.12 |
| **TCGA-Toil/RSEM (foreign)** | **cohort z-score** | **0.994** | **0.935** | **0.935** | **0.11** |
| CPTAC-3/GDC (native) | deployed | 0.989 | 0.955 | 0.955 | 0.08 |
| CPTAC-3/GDC (native) | cohort z-score | 0.990 | 0.930 | 0.930 | 0.14 |
| GTEx/Toil normals (all-normal) | deployed | - | - | - | 0.996 |
| GTEx/Toil normals (all-normal) | cohort center | - | - | - | 0.907 |
| GTEx/Toil normals (all-normal) | cohort z-score | - | - | - | 0.915 |

**Headline:** on the foreign Toil/RSEM pipeline, cohort standardization lifts
default-threshold accuracy from **0.515 to 0.935** (specificity 0.03 -> 0.89),
label-free and without retraining, while AUC is unchanged (0.992 -> 0.994).

**Native pipeline:** on CPTAC (already GDC STAR-Counts) it costs a little
(0.955 -> 0.930). Use adaptation only for non-GDC pipelines; keep `--adapt none`
for GDC STAR-Counts input.

## Important limitation: cohorts need internal contrast

Cohort standardization recenters on the cohort mean, which is only the training
mean if the cohort has a comparable tumor/normal mix. A near-single-class cohort
has no internal contrast to anchor the recentering. On the **all-normal** GTEx
panel it only partially helps (FPR 0.996 -> 0.915 at 0.5; 0.965 -> 0.41 at 0.99).

Subsampling the Toil cohort to different tumor fractions makes the dependence
explicit (balanced accuracy, mean of 200 resamples, n=80):

| Cohort tumor fraction | deployed | cohort z-score | AUC (adapted) |
|---:|---:|---:|---:|
| 0.10 | 0.514 | 0.607 | 0.996 |
| 0.25 | 0.515 | 0.757 | 0.994 |
| 0.50 | 0.516 | 0.920 | 0.994 |
| 0.75 | 0.516 | 0.974 | 0.996 |
| 0.90 | 0.513 | 0.956 | 0.997 |

AUC is preserved at every composition; threshold recovery is strong for mixed
cohorts (>= ~0.9 balanced accuracy for 50-90% tumor) and weak for near-pure
cohorts. For a near-single-class cohort, prefer an explicit labeled-anchor
recalibration (`calibrate_threshold.py`) over cohort standardization.

## Usage

```bash
# foreign pipeline (e.g. Toil/RSEM), mixed cohort -> restore the 0.5 threshold
python cohort_adapt_score.py input.csv --adapt cohort_zscore --out scores.csv

# with labels, also print AUC / accuracy / balanced accuracy
python cohort_adapt_score.py input.csv --adapt cohort_zscore --labels labels.csv

# native GDC STAR-Counts input -> no adaptation
python cohort_adapt_score.py input.csv --adapt none
```

`input.csv`: rows = samples, columns = Ensembl gene IDs, values = `log2(TPM+1)`.
Missing model genes are imputed at the training mean (neutral). Requires only
numpy and pandas.

## Reproduce

```bash
cd cross-platform-adaptation
python run_adaptation_benchmark.py     # writes adaptation_benchmark.csv, adaptation_imbalance.csv
```

Reads the sibling `external-validation/` cohort matrices and the parent
`deployable_lr_weights.npz`.

## Files

- `../cohort_adapt_score.py` - domain-adapted scoring library + CLI
- `run_adaptation_benchmark.py` - reproduces the tables from released data
- `adaptation_benchmark.csv` - per-cohort, per-mode metrics
- `adaptation_imbalance.csv` - balanced accuracy vs cohort tumor fraction
- `cross_platform_adaptation.html` - visual summary
