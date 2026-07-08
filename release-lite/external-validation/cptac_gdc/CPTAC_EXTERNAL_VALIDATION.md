# CPTAC external smoke validation

## Data source

- Source: NCI GDC API, CPTAC program
- Files queried: 2879
- Project scored: CPTAC-3
- Sampled files: 200 (100 primary tumor, 100 solid tissue normal)
- Workflow: GDC STAR-Counts, `tpm_unstranded` converted to log2(TPM+1)
- Model: bundled logistic regression from `deployable_pipeline.pkl`

## Result at threshold 0.5

| Metric | Value |
|---|---:|
| AUC | 0.9886 |
| Average precision | 0.9851 |
| Accuracy | 0.9550 |
| F1 | 0.9565 |
| Precision | 0.9252 |
| Recall | 0.9900 |
| True normal / false tumor | 92 / 8 |
| False normal / true tumor | 1 / 99 |

## Probability by label

| Sample type | n | Mean tumor probability | Tumor calls |
|---|---:|---:|---:|
| Primary Tumor | 100 | 0.9834 | 99 |
| Solid Tissue Normal | 100 | 0.1132 | 8 |

## Threshold sensitivity

| Threshold | Cutoff | Accuracy | Precision | Recall | TN / FP | FN / TP |
|---|---:|---:|---:|---:|---:|---:|
| default_0.5 | 0.500000 | 0.9550 | 0.9252 | 0.9900 | 92 / 8 | 1 / 99 |
| youden_j | 0.492041 | 0.9600 | 0.9259 | 1.0000 | 92 / 8 | 0 / 100 |
| high_specificity_0.75 | 0.750000 | 0.9500 | 0.9245 | 0.9800 | 92 / 8 | 2 / 98 |
| high_specificity_0.9 | 0.900000 | 0.9500 | 0.9412 | 0.9600 | 94 / 6 | 4 / 96 |
| high_specificity_0.95 | 0.950000 | 0.9500 | 0.9688 | 0.9300 | 97 / 3 | 7 / 93 |
| very_high_specificity_0.99 | 0.990000 | 0.9250 | 0.9885 | 0.8600 | 99 / 1 | 14 / 86 |

## Score quantiles

| Sample type | min | p10 | p25 | median | p75 | p90 | max |
|---|---:|---:|---:|---:|---:|---:|---:|
| Primary Tumor | 0.4920 | 0.9799 | 0.9983 | 0.9999 | 1.0000 | 1.0000 | 1.0000 |
| Solid Tissue Normal | 0.0001 | 0.0003 | 0.0009 | 0.0044 | 0.0476 | 0.3979 | 1.0000 |

## Interpretation

This is an external smoke test, not a replacement for a full independent
benchmark. It is outside TCGA but still uses GDC harmonized STAR-Counts, so it
tests cohort transfer more than cross-platform transfer. The stricter remaining
gap is non-GDC RNA-seq, where normalization and gene annotation differences can
dominate.
