# Experimental cross-platform cohort standardization

**Date:** 2026-07-05
**Model:** frozen `deployable_lr_weights.npz` (weights unchanged in
`v2.0.0-gdc-starcounts`)
**Status:** retrospective research benchmark; not a calibrated deployment mode

## Problem

The released model is a logistic regression on 2,000 genes fitted to **GDC
STAR-Counts `log2(TPM+1)`**. Historical external benchmarks showed a sharp
boundary:

- Discrimination (AUC) transfers across pipelines: TCGA-Toil/RSEM **AUC 0.992**,
  CPTAC-3/GDC **AUC 0.989**.
- The fixed **0.5 decision threshold does not** transfer. On the Toil/RSEM
  pipeline the deployed model calls 197/200 samples tumor (accuracy **0.515**,
  specificity **0.03**); on 540 GTEx/Toil normal tissues it calls 538 tumor
  (false-positive rate **0.996**, median probability 0.9999).

The release therefore tells users *not* to make hard calls on non-GDC pipelines
without separate validation. This document records an experimental, label-free
transform evaluated on the frozen historical cohorts. It is not calibration in
the statistical sense and no post-cache-fix live-network rerun has been made.

## Method

One plausible contributor is per-gene **location/scale shift** between
pipelines: after
`log2(TPM+1)`, Toil/RSEM gene distributions sit at systematically different
means and spreads than GDC STAR-Counts, so the training-fitted standardizer
`z = (x - train_mean)/train_scale` pushes foreign samples far into "tumor"
territory.

**Cohort standardization** attempts to reduce that shift without labels or
retraining by standardizing each gene using the *input cohort's own* per-gene
mean and standard deviation before applying the frozen coefficients.

```
z_i = (x_i - cohort_mean) / cohort_std        # cohort_zscore  (experimental)
z_i = (x_i - cohort_mean) / train_scale       # cohort_center  (location only)
p   = sigmoid(coef . z_i + intercept)
```

This is a transductive, composition-dependent transform: changing the samples
in the batch changes every score. It does not make the foreign distribution
identical to training data, preserve ranking by mathematical guarantee, restore
probability calibration, or establish that threshold 0.5 is valid.

## Results

Frozen model, three historical external cohort snapshots, three scoring modes.
The NumPy reimplementation reproduced the recorded deployed scores
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

**Historical observation:** in this fixed, balanced TCGA-Toil/RSEM cohort,
cohort standardization changed default-threshold accuracy from **0.515 to
0.935** (specificity 0.03 -> 0.89), without labels or retraining. This is a
retrospective result on the same cohort used to characterize the transform, not
an independent estimate of performance on a new batch.

**Native pipeline:** on CPTAC (already GDC STAR-Counts) it costs a little
(0.955 -> 0.930). The safe default is `--adapt none`, including for GDC
STAR-Counts input. Foreign pipelines require separate validation; adaptation is
an explicit experimental option, not a general remedy.

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

In these particular resamples AUC stayed similar, while threshold behavior
changed sharply with composition. Neither result is guaranteed for another
cohort. For a near-single-class cohort, do not use cohort standardization. If
labels are available, `calibrate_threshold.py` can estimate a cohort-specific
threshold, but its reported metrics are apparent/resubstitution estimates until
confirmed on independent data.

## Usage

```bash
# Explicit experimental opt-in for a reviewed foreign-pipeline mixed cohort
python cohort_adapt_score.py input.csv --adapt cohort_zscore --out scores.csv

# with labels, also print AUC / accuracy / balanced accuracy
python cohort_adapt_score.py input.csv --adapt cohort_zscore --labels labels.csv

# Default and recommended starting point: no adaptation
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

Reads the trusted local historical cohort matrices under sibling
`external-validation/` and the parent `deployable_lr_weights.npz`. The command
recomputes the retrospective tables; it does not redownload or independently
revalidate the cohorts.

## Files

- `../cohort_adapt_score.py` - domain-adapted scoring library + CLI
- `run_adaptation_benchmark.py` - reproduces the tables from released data
- `adaptation_benchmark.csv` - per-cohort, per-mode metrics
- `adaptation_imbalance.csv` - balanced accuracy vs cohort tumor fraction
- `cross_platform_adaptation.html` - visual summary
