# Pan-cancer tumor-vs-normal classifier from TCGA RNA-seq

Release: `v1.1.5-gdc-starcounts` (`2026-07-09`)

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

See the root-level LOCO artifacts in `cross-cancer-holdout/`:
`LOCO_REPORT.md`, `loco_report.html`, `loco_per_cancer_metrics.csv`,
`loco_pooled_summary.csv`, `loco_vs_within_comparison.csv`, and
`loco_predictions.pkl`. (The raw `loco_generalization*` figure/CSVs from the
original workbench run live in `from-workbench-loco/`.)

## External CPTAC-3 smoke validation

As a first external check outside TCGA, the deployed logistic-regression model
was scored on a stratified sample of 200 CPTAC-3 RNA-seq files pulled from the
GDC API: 100 primary tumors and 100 solid tissue normals. These files are not
from TCGA, but they are still processed by the GDC harmonized STAR-Counts
workflow, so this test primarily measures transfer to an independent cohort
while holding the expression pipeline constant.

- **AUC: 0.989**, average precision 0.985.
- At the default 0.5 threshold: accuracy 0.955, precision 0.925, recall 0.990.
  Confusion matrix: 92 true normals, 8 false positives, 1 false negative,
  99 true tumors.
- The model matched all 2,000 selected genes in the CPTAC files. Tumor scores
  remained highly separated from normal scores (tumor median 0.9999; normal
  median 0.0044), but a small tail of CPTAC normal samples scored high,
  reinforcing that hard calls depend on cohort- and tissue-specific threshold
  calibration.

See `external-validation/cptac_gdc/CPTAC_EXTERNAL_VALIDATION.md`,
`external-validation/cptac_gdc/cptac_predictions.csv`, and
`external-validation/cptac_gdc/cptac_threshold_sweep.csv`.

## Cross-platform Toil/GTEx boundary check

The stricter question is whether the deployed GDC STAR-Counts model can be
applied directly to non-GDC RNA-seq matrices. To test this, GTEx normal tissues
were sampled from the UCSC Xena Toil RNA-seq recompute compendium
(`gtex_RSEM_gene_tpm`). Values were converted from Xena's log2(TPM+0.001) scale
back to TPM and then to log2(TPM+1), matching the nominal input transform used
by the deployed model.

Result: the default deployed model is **not cross-platform safe**.

- GTEx normal-tissue check: 540 GTEx normal samples across 27 primary sites.
  At threshold 0.5, 538/540 were called tumor (false-positive rate 0.996).
- A companion TCGA Toil/RSEM sanity check showed that this is mainly a
  pipeline/threshold-transfer problem, not simply a GTEx biology claim. On 100
  TCGA primary tumors and 100 TCGA solid tissue normals from Toil, ranking
  remained strong (AUC 0.992), but the 0.5 threshold miscalled 97/100 normals
  as tumor. Re-thresholding within Toil at Youden's-J cutoff 0.999975 recovered
  accuracy to 0.970.
- Applying that high Toil-derived threshold to GTEx still left a high normal
  false-positive rate (216/540 = 0.400), so GTEx/Toil should not be used for
  hard calls without a dedicated refit/recalibration strategy.

This sharply defines the current model's deployment boundary: it is strong for
GDC STAR-Counts-scale TCGA/CPTAC-style tumor-vs-adjacent-normal contrasts, but
not a plug-in classifier for arbitrary TPM/RSEM/GTEx/GEO matrices.

See `external-validation/gtex_xena/GTEX_NORMAL_VALIDATION.md` and
`external-validation/tcga_toil_xena/TCGA_TOIL_PIPELINE_CHECK.md`.

## Input compatibility QC

To reduce accidental misuse, the lightweight release includes
`inspect_expression_input.py`. It checks the input matrix before hard calls:

- 2,000 model-gene match rate, including Ensembl IDs with or without version
  suffixes.
- Missing, non-numeric, non-finite, negative, or unusually large expression
  values.
- Standardized distribution-shift summaries using the logistic-regression
  training scaler.
- Score distribution and optional expected-cohort warnings, e.g. normal-only
  input with most samples called tumor.

The bundled TCGA example and CPTAC/GDC validation matrix pass this QC. The
TCGA Toil/RSEM and GTEx/Toil boundary checks are flagged as WARN because they
show excess standardized outliers and cohort-level gene-mean shifts, matching
the observed threshold-transfer failure. A PASS result is not a substitute for
external validation, but WARN/FAIL is a clear stop signal before using hard
calls.

For deployment, `run_tumor_normal_workflow.py` wraps the release tools into a
single command. It writes QC, scores, optional labeled-threshold calibration,
per-sample explanations, a machine-readable manifest, and a Markdown workflow
report in one output directory.

Matched model-gene cells that are missing, non-numeric, `NaN`, or infinite stop
the scorer, workflow, and explainer before output files are written unless a
reviewed tolerance or explicit override is supplied. This prevents a malformed
matched input column or sample from being silently neutral-imputed.

## Cross-platform adaptation (restoring Toil accuracy)

The Toil/GTEx boundary above is a threshold- and scale-transfer failure, not a
loss of biological signal. A label-free cohort standardization step aligns an
incoming non-GDC matrix to the model's training distribution before scoring,
which restores Toil accuracy from 0.515 to 0.935 without refitting the model.
This is packaged as `cohort_adapt_score.py`. See
`cross-platform-adaptation/CROSS_PLATFORM_ADAPTATION.md` for the method and
benchmark.

## Cancer-type classifier (tissue of origin)

This tumor-vs-normal model deliberately does not identify which cancer type a
sample is. A separate 17-class tissue-of-origin classifier ships alongside it,
reaching patient-held-out accuracy 0.930 (balanced accuracy 0.878, macro-F1
0.877). See `cancer-type-classifier/CANCER_TYPE_CLASSIFIER.md` and the
`cancer-type-classifier/predict_cancer_type.py` entry point.

## Caveats

- This is a **pan-cancer** tumor-vs-normal classifier, not a cancer-type
  classifier — it distinguishes malignant from non-malignant tissue across
  cancer types, largely by picking up a shared stromal/ECM remodeling
  signature. The LOCO analysis above confirms this signature generalizes to
  entirely unseen cancer types (mean AUC 0.994), addressing the original
  concern that held-out performance was measured only on cancer types also
  seen in training.
- CPTAC-3 provides an external non-TCGA cohort check (AUC 0.989), but it is
  still GDC harmonized STAR-Counts. Non-GDC / cross-platform RNA-seq was tested
  with UCSC Xena Toil/GTEx and failed as a direct hard-call deployment target.
  The model should be refit or recalibrated for those pipelines.
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
- `cross-cancer-holdout/LOCO_REPORT.md`, `cross-cancer-holdout/loco_report.html` — leave-one-cancer-type-out write-up and interactive figure/table
- `cross-cancer-holdout/loco_per_cancer_metrics.csv`, `cross-cancer-holdout/loco_pooled_summary.csv`, `cross-cancer-holdout/loco_vs_within_comparison.csv`, `cross-cancer-holdout/loco_predictions.pkl` — per-cancer-type LOCO metrics, pooled summary, comparison, and pooled predictions
- `from-workbench-loco/loco_generalization.png`, `from-workbench-loco/loco_generalization.csv`, `from-workbench-loco/loco_generalization_with_optimal_threshold.csv` — raw LOCO figure/CSVs from the original workbench run
- `external-validation/cptac_gdc/CPTAC_EXTERNAL_VALIDATION.md` — external CPTAC-3 smoke validation
- `external-validation/validate_cptac_gdc.py` — reproducible GDC/CPTAC validation script
- `external-validation/gtex_xena/GTEX_NORMAL_VALIDATION.md` — GTEx normal-tissue cross-platform check
- `external-validation/tcga_toil_xena/TCGA_TOIL_PIPELINE_CHECK.md` — Toil/RSEM pipeline-transfer check
- `test_metrics.csv` — held-out test metrics for all 3 models
- `per_cancer_type_performance.csv` — per-cancer-type AUC/accuracy breakdown
- `top_genes_xgboost.csv`, `top_genes_logreg.csv` — top 30 genes by each model
- `EXECUTIVE_SUMMARY.md`, `VERSION`, `RELEASE_METADATA.json` — short handoff and release metadata
- `USER_GUIDE.md`, `templates/` — practical input-preparation guide and CSV sketches
- `DATA_DICTIONARY.md` — stable input/output columns and JSON contract reference
- `TROUBLESHOOTING.md` — common install, input-QC, threshold, and release-integrity fixes
- `example_workflow_output/` — reference output from the one-command example workflow
- `deployable_pipeline.pkl` — fitted selector + scaler + all 3 trained models
  + the exact 2,000-gene feature list, ready to score new TCGA-style
  log2(TPM+1) expression vectors
- `deployable_lr_weights.npz` — pure NumPy export of the default LR scorer
- `run_tumor_normal_workflow.py` — one-command QC, scoring, calibration, explanations, and report
- `check_environment.py` — runtime/package/file diagnostic and bundled self-test wrapper
- `audit_lightweight_dependencies.py` — release import and minimal dependency audit
- `audit_cli_entrypoints.py` — release CLI `--help` and shebang audit
- `audit_release_docs.py` — documentation and release-bundle reference audit
- `validate_output_contracts.py` — bundled CSV/JSON output contract validator
- `inspect_expression_input.py`, `model_qc_reference.json` — pre-scoring input compatibility QC
- `model_gene_metadata.csv` — model gene coefficients, scaling metadata, and known gene names
- `calibrate_threshold.py` — threshold calibration utility for labeled scored samples
- `explain_scores.py` — per-sample LR logit contribution explanations
- `export_lr_weights.py` — regenerates the pure NumPy LR weights from `deployable_pipeline.pkl`
- `export_qc_reference.py` — regenerates the QC reference JSON from validation matrices
- `export_model_gene_metadata.py` — regenerates the model gene metadata table
- `build_release_lite.py` — regenerates the lightweight release folder, checksums, and zip
- `validate_release_lite.py` — validates release folder, manifest, checksums, zip, and forbidden artifacts
- `validate_zip_bundle.py` — extracts the zip into a clean temp directory and runs acceptance
- `run_safety_tests.py` — verifies guardrails for invalid inputs and QC-fail workflow behavior
- `run_release_acceptance.py` — one-command environment, smoke, safety, and release-integrity checks
- `RELEASE_ARTIFACTS.json` — generated sidecar with zip size and SHA256
- `selected_files.csv` — full manifest of the 2,160 GDC file IDs / case IDs /
  cancer types / labels used
