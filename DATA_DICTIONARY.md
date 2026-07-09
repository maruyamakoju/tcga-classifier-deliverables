# Data dictionary

Release: `v1.1.11-gdc-starcounts` (`2026-07-09`)

This file defines the stable input and output fields used by the lightweight
release. Run `python validate_output_contracts.py` to verify the bundled
examples still match these contracts.

## Input Matrix

Expression input files use samples as rows and Ensembl gene IDs as columns.
The first column is the sample identifier / row index. Values must be
`log2(TPM+1)` on a GDC STAR-Counts-style scale. Ensembl version suffixes are
accepted.

## Labels CSV

`example_labels.csv` and calibration label files use:

| Column | Type | Meaning |
|---|---|---|
| `sample` | string | Sample identifier matching the scored matrix |
| `label` | string/int | `tumor`, `normal`, `1`, or `0` |

## Scores CSV

`score_tumor_normal.py` and the workflow write:

| Column | Type | Meaning |
|---|---|---|
| `sample` | string | Sample identifier |
| `tumor_probability` | float | Logistic-regression tumor probability in `[0, 1]` |
| `call` | string | `tumor` if probability is at or above the threshold, else `normal` |

## Thresholds CSV

`calibrate_threshold.py` writes:

| Column | Type | Meaning |
|---|---|---|
| `threshold_name` | string | `default` or `youden_j` |
| `threshold` | float | Probability cutoff in `[0, 1]` |
| `accuracy` | float | Accuracy at the cutoff |
| `f1` | float | Tumor-class F1 |
| `precision` | float | Tumor-class precision |
| `recall` | float | Tumor-class recall |
| `specificity` | float | Normal-class recall |
| `tn`, `fp`, `fn`, `tp` | int | Confusion-matrix counts |
| `youden_j` | float/blank | Recall + specificity - 1 |

## Explanations CSV

`explain_scores.py` writes per-sample LR logit contributions:

| Column | Type | Meaning |
|---|---|---|
| `sample` | string | Sample identifier |
| `tumor_probability` | float | Tumor probability in `[0, 1]` |
| `logit` | float | Logistic-regression logit before sigmoid |
| `direction` | string | `tumor` for positive contributions or `normal` for negative contributions |
| `rank` | int | Rank within sample and direction |
| `gene_id` | string | Ensembl gene ID used by the model |
| `gene_name` | string/blank | Gene symbol when available |
| `contribution_logit` | float | Gene contribution to the LR logit |
| `expression_log2_tpm1` | float | Input expression value |
| `training_mean` | float | Training-set scaler mean |
| `scaled_value` | float | Standardized input value |
| `lr_coef` | float | Logistic-regression coefficient |

## QC JSON

`inspect_expression_input.py` writes a JSON report with these top-level keys:

| Key | Meaning |
|---|---|
| `status` | `PASS`, `WARN`, or `FAIL` |
| `shape` | sample count, input gene count, model gene count |
| `gene_match` | matched/missing model-gene counts and match rate |
| `value_summary` | finite/non-finite, negative, and expression range summaries |
| `distribution_summary` | standardized distribution-shift summaries |
| `score_summary` | score/call distribution at the selected threshold |
| `threshold` | probability threshold used for calls |
| `expected_class` | expected cohort label used for optional strict checks |
| `messages` | warning/error messages with machine-readable codes |
| `reference_source` | provenance note for QC rules |

## Workflow Manifest

`run_tumor_normal_workflow.py` writes `manifest.json` with input paths, sample
and gene counts, QC status, call counts, threshold, output filenames, and
optional calibration summary.

When the workflow stops before scoring, `status` can be
`stopped_after_qc_fail` or `stopped_after_invalid_input`. The invalid-input
case includes an `alignment` object summarizing invalid matched model-gene
cells and example affected genes/samples.

## Model Gene Metadata

`model_gene_metadata.csv` contains one row per model gene:

| Column | Meaning |
|---|---|
| `rank_abs_lr_coef` | Rank by absolute LR coefficient |
| `gene_id` | Versioned Ensembl gene ID |
| `gene_id_base` | Ensembl gene ID without version suffix |
| `gene_name` | Gene symbol when available |
| `lr_coef` | Logistic-regression coefficient |
| `abs_lr_coef` | Absolute coefficient |
| `direction_if_high` | Direction implied by high expression, `tumor` or `normal` |
| `scaler_mean` | Training-set scaler mean |
| `scaler_scale` | Training-set scaler scale |
