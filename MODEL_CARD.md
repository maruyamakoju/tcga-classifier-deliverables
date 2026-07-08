# Model card: TCGA/GDC tumor-vs-normal RNA-seq classifier

Release: `v1.1.3-gdc-starcounts` (`2026-07-08`)

## Intended use

This model scores bulk RNA-seq samples as tumor vs normal for **GDC
STAR-Counts-style log2(TPM+1)** expression matrices. It is best treated as a
research classifier for TCGA/CPTAC-like tumor-vs-adjacent-normal contrasts.

Default deployment uses logistic regression exported to
`deployable_lr_weights.npz`, which avoids scikit-learn pickle compatibility
issues for ordinary scoring.

## Not intended for

- Clinical diagnosis or patient management.
- Cancer-type classification by THIS tumor-vs-normal model (a separate 17-class
  classifier ships in `cancer-type-classifier/`).
- Direct hard calls on arbitrary TPM/RSEM/GTEx/GEO matrices.
- Single-cell RNA-seq, spatial transcriptomics, microarray, raw counts, FPKM, or
  unlogged TPM without conversion to the expected scale.
- Healthy-population normal-vs-tumor claims without additional validation.

## Training data

- Source: NCI GDC TCGA STAR-Counts RNA-seq.
- Samples: 2,160 total from 17 TCGA projects.
- Labels: 1,440 primary tumor, 720 solid tissue normal.
- Split: patient-grouped train/test split; matched tumor/normal patients were
  kept within the same split.
- Features: top 2,000 genes selected by ANOVA F-test on training data only.

## Model

- Default model: logistic regression, L2 regularized, class-balanced,
  standardized features.
- Input genes: Ensembl IDs, with or without version suffixes.
- Missing model genes: filled with the training mean.
- Default threshold: 0.5, but threshold calibration is recommended for new
  tissues or pipelines.
- Input QC: `inspect_expression_input.py` checks model-gene coverage,
  expression range, standardized distribution shift, and score saturation before
  hard calls.
- Explanations: per-gene contributions from `explain_scores.py` are additive
  terms in the logistic-regression logit. They are useful for debugging model
  behavior, but they are not causal biological explanations.

## Validation summary

| Setting | Key result |
|---|---:|
| TCGA patient-held-out test | AUC 0.997, accuracy 0.979 |
| TCGA grouped 5-fold CV | AUC 0.997 +/- 0.003 |
| TCGA leave-one-cancer-type-out | macro AUC 0.994 |
| CPTAC-3/GDC STAR-Counts smoke test | AUC 0.989, accuracy 0.955 |
| TCGA Toil/RSEM pipeline check | AUC 0.992, but default-threshold accuracy 0.515 |
| GTEx/Toil normal tissue check | 538/540 normals called tumor at threshold 0.5 |

## Key limitations

The model is **GDC STAR-Counts-scale specific**. It transfers well to CPTAC-3
when the GDC STAR-Counts pipeline is held constant, but direct use on
UCSC Xena Toil/RSEM matrices is not safe for hard calls. On Toil/RSEM, ranking
can remain strong while probabilities and thresholds shift dramatically.

TCGA solid tissue normals are tumor-adjacent tissues, not healthy-donor normals.
Therefore the validated target is tumor-vs-adjacent-normal, not a general
healthy-vs-cancer diagnostic setting.

## Recommended deployment workflow

1. Prepare expression matrix as rows=samples, columns=Ensembl genes, values
   log2(TPM+1) from GDC STAR-Counts.
2. Prefer the bundled workflow:
   `python run_tumor_normal_workflow.py input.csv --labels labels.csv`.
3. Review `workflow_report.md`, especially the QC status and messages.
4. Treat any WARN/FAIL QC result as a reason to check normalization, gene IDs,
   and whether the matrix is truly GDC STAR-Counts-style.
5. If labeled calibration samples are available, use the reported
   Youden's-J threshold instead of assuming the default 0.5 cutoff transfers.
6. Use `scores.csv` for calls and `explanations.csv` only for model debugging.

The individual tools (`inspect_expression_input.py`, `score_tumor_normal.py`,
`calibrate_threshold.py`, and `explain_scores.py`) remain available for custom
pipelines.

## Files to cite

- `REPORT.md`: full methods and results.
- `REPRODUCIBILITY.md`: runtime, smoke tests, and validation status.
- `external-validation/cptac_gdc/CPTAC_EXTERNAL_VALIDATION.md`: CPTAC-3 check.
- `external-validation/gtex_xena/GTEX_NORMAL_VALIDATION.md`: GTEx boundary check.
- `external-validation/tcga_toil_xena/TCGA_TOIL_PIPELINE_CHECK.md`: Toil/RSEM
  pipeline-transfer check.
