# GTEx normal-tissue cross-platform check

## Data source

- Source: UCSC Xena Toil RNA-seq recompute compendium
- Dataset: `gtex_RSEM_gene_tpm`
- Phenotype table: `TcgaTargetGTEX_phenotype`
- Samples scored: 540 GTEx normal samples, stratified by primary site
- Input transform: Xena log2(TPM+0.001) -> TPM -> log2(TPM+1)
- Model: bundled logistic regression from `deployable_pipeline.pkl`

## Result at threshold 0.5

| Metric | Value |
|---|---:|
| Samples | 540 |
| Normal calls | 2 |
| Tumor calls / false positives | 538 |
| False positive rate | 0.9963 |
| Median tumor probability | 0.9999 |
| p90 tumor probability | 1.0000 |
| p95 tumor probability | 1.0000 |
| Max tumor probability | 1.0000 |

## Threshold sensitivity

| Threshold | Tumor calls | False positive rate |
|---:|---:|---:|
| 0.500000 | 538 | 0.996 |
| 0.750000 | 538 | 0.996 |
| 0.900000 | 537 | 0.994 |
| 0.950000 | 536 | 0.993 |
| 0.990000 | 521 | 0.965 |
| 0.999000 | 451 | 0.835 |
| 0.999900 | 305 | 0.565 |
| 0.999975 | 216 | 0.400 |
| 0.999990 | 161 | 0.298 |
| 0.999999 | 101 | 0.187 |

## Per-site summary

| Primary site | n | Tumor calls | False positive rate | Median tumor probability | Max tumor probability |
|---|---:|---:|---:|---:|---:|
| Adrenal Gland | 20 | 20 | 1.000 | 1.0000 | 1.0000 |
| Blood Vessel | 20 | 20 | 1.000 | 1.0000 | 1.0000 |
| Heart | 20 | 20 | 1.000 | 0.9999 | 1.0000 |
| Liver | 20 | 20 | 1.000 | 0.9998 | 1.0000 |
| Muscle | 20 | 20 | 1.000 | 1.0000 | 1.0000 |
| Pituitary | 20 | 20 | 1.000 | 1.0000 | 1.0000 |
| Prostate | 20 | 20 | 1.000 | 1.0000 | 1.0000 |
| Salivary Gland | 20 | 20 | 1.000 | 1.0000 | 1.0000 |
| Small Intestine | 20 | 20 | 1.000 | 0.9999 | 1.0000 |
| Spleen | 20 | 20 | 1.000 | 1.0000 | 1.0000 |
| Thyroid | 20 | 20 | 1.000 | 1.0000 | 1.0000 |
| Lung | 20 | 20 | 1.000 | 0.9998 | 1.0000 |
| Nerve | 20 | 20 | 1.000 | 0.9998 | 1.0000 |
| Ovary | 20 | 20 | 1.000 | 1.0000 | 1.0000 |
| Esophagus | 20 | 20 | 1.000 | 0.9998 | 1.0000 |
| Pancreas | 20 | 20 | 1.000 | 1.0000 | 1.0000 |
| Uterus | 20 | 20 | 1.000 | 0.9998 | 1.0000 |
| Stomach | 20 | 20 | 1.000 | 0.9993 | 1.0000 |
| Vagina | 20 | 20 | 1.000 | 0.9998 | 1.0000 |
| Skin | 20 | 20 | 1.000 | 0.9999 | 1.0000 |
| Testis | 20 | 20 | 1.000 | 0.9999 | 1.0000 |
| Breast | 20 | 20 | 1.000 | 0.9985 | 1.0000 |
| Colon | 20 | 20 | 1.000 | 0.9997 | 1.0000 |
| Adipose Tissue | 20 | 20 | 1.000 | 0.9976 | 1.0000 |
| Blood | 20 | 20 | 1.000 | 0.9994 | 1.0000 |
| Kidney | 20 | 19 | 0.950 | 1.0000 | 1.0000 |
| Brain | 20 | 19 | 0.950 | 0.9986 | 1.0000 |

## Highest-scoring GTEx normals

| Sample | Primary site | Tissue | Tumor probability | Call |
|---|---|---|---:|---|
| GTEX-11EMC-0526-SM-5EGJN | Adrenal Gland | Adrenal Gland | 1.000000 | tumor |
| GTEX-11PRG-2226-SM-5GU5R | Kidney | Kidney - Cortex | 1.000000 | tumor |
| GTEX-ZYWO-0126-SM-5GZWQ | Adrenal Gland | Adrenal Gland | 1.000000 | tumor |
| GTEX-ZE9C-1426-SM-4WKGM | Kidney | Kidney - Cortex | 1.000000 | tumor |
| GTEX-13S86-0126-SM-5S2PI | Spleen | Spleen | 1.000000 | tumor |
| GTEX-13O21-0126-SM-5IJE8 | Adrenal Gland | Adrenal Gland | 1.000000 | tumor |
| GTEX-12584-3126-SM-5EGKR | Pituitary | Pituitary | 1.000000 | tumor |
| GTEX-13PVR-0226-SM-5RQJI | Adrenal Gland | Adrenal Gland | 1.000000 | tumor |
| GTEX-WQUQ-0126-SM-4OOSS | Adrenal Gland | Adrenal Gland | 1.000000 | tumor |
| GTEX-1477Z-2826-SM-5SI9J | Pituitary | Pituitary | 1.000000 | tumor |
| GTEX-12WSB-0726-SM-5N9GD | Thyroid | Thyroid | 1.000000 | tumor |
| GTEX-XOTO-2826-SM-4B65I | Prostate | Prostate | 1.000000 | tumor |
| GTEX-12WSN-0126-SM-5DUX5 | Spleen | Spleen | 1.000000 | tumor |
| GTEX-11TUW-0226-SM-5LU8X | Thyroid | Thyroid | 1.000000 | tumor |
| GTEX-ZAB5-2326-SM-5IJFR | Salivary Gland | Minor Salivary Gland | 1.000000 | tumor |
| GTEX-13OW8-0426-SM-5J2NR | Prostate | Prostate | 1.000000 | tumor |
| GTEX-14E6E-1226-SM-5S2R5 | Spleen | Spleen | 1.000000 | tumor |
| GTEX-1211K-0126-SM-59HJE | Adrenal Gland | Adrenal Gland | 1.000000 | tumor |
| GTEX-11ZVC-3226-SM-5FQV1 | Pituitary | Pituitary | 1.000000 | tumor |
| GTEX-14BMV-2326-SM-5RQJ4 | Prostate | Prostate | 1.000000 | tumor |

## Interpretation

This is a stricter platform-transfer check than the CPTAC/GDC validation because
it uses GTEx/Toil normal tissues rather than GDC STAR-Counts. It only measures
normal-sample false positives, not tumor-vs-normal AUC. A high false positive
rate would indicate that the TCGA/GDC-trained threshold does not transfer cleanly
to GTEx/Toil. The companion TCGA Toil/RSEM check shows that even TCGA samples
require an extreme threshold shift on this pipeline, so this result should be
treated as a pipeline/domain-transfer failure of the deployed GDC STAR-Counts
model rather than evidence that GTEx tissues are biologically tumor-like.
