# TCGA Toil/RSEM pipeline-transfer sanity check

## Data source

- Source: UCSC Xena Toil RNA-seq recompute compendium
- Dataset: `tcga_RSEM_gene_tpm`
- Phenotype table: `TcgaTargetGTEX_phenotype`
- Samples scored: 200 TCGA samples (100 primary tumor, 100 solid tissue normal)
- Input transform: Xena log2(TPM+0.001) -> TPM -> log2(TPM+1)
- Model: bundled logistic regression from `deployable_pipeline.pkl`

## Result at threshold 0.5

| Metric | Value |
|---|---:|
| AUC | 0.9923 |
| Average precision | 0.9936 |
| Accuracy | 0.5150 |
| F1 | 0.6734 |
| Precision | 0.5076 |
| Recall | 1.0000 |
| True normal / false tumor | 3 / 97 |
| False normal / true tumor | 0 / 100 |

## Probability by label

| Sample type | n | Mean tumor probability | Median tumor probability | Tumor calls |
|---|---:|---:|---:|---:|
| Primary Tumor | 100 | 1.0000 | 1.0000 | 100 |
| Solid Tissue Normal | 100 | 0.9536 | 0.9905 | 97 |

## Threshold sensitivity

| Threshold | Cutoff | Accuracy | Precision | Recall | TN / FP | FN / TP |
|---|---:|---:|---:|---:|---:|---:|
| default_0.5 | 0.500000 | 0.5150 | 0.5076 | 1.0000 | 3 / 97 | 0 / 100 |
| youden_j | 0.999975 | 0.9700 | 0.9796 | 0.9600 | 98 / 2 | 4 / 96 |
| high_0.99 | 0.990000 | 0.7500 | 0.6667 | 1.0000 | 50 / 50 | 0 / 100 |
| high_0.999 | 0.999000 | 0.8450 | 0.7674 | 0.9900 | 70 / 30 | 1 / 99 |
| high_0.9999 | 0.999900 | 0.9550 | 0.9417 | 0.9700 | 94 / 6 | 3 / 97 |
| high_0.99999 | 0.999990 | 0.9600 | 0.9894 | 0.9300 | 99 / 1 | 7 / 93 |

## Interpretation

This is not a new biological validation cohort because TCGA overlaps the model's
training source. It is a pipeline-transfer sanity check. Poor hard-call behavior
here means the deployed GDC STAR-Counts model should not be applied directly to
Toil/RSEM matrices without refitting or recalibrating on that expression
pipeline.
