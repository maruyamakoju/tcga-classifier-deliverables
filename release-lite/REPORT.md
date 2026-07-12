# Pan-cancer tumor-vs-normal classifier from TCGA RNA-seq

Release: `v2.2.0-gdc-starcounts` (`2026-07-12`; public scoring-library API `3.0.0`)

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
cross-validation within the training cohort gives a similar estimate:
logistic regression 0.996 ± 0.003 AUC, random forest 0.996 ± 0.002 AUC across
folds. This consistency reduces dependence on the single held-out split, but
does not constitute external or prospective validation.

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

## Leave-one-cancer-type-out (LOCO) sensitivity analysis

To measure performance when a TCGA cancer project is omitted from fitting, the
pipeline was retrained 17 times,
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
  that a different threshold can fit the observed held-out labels, not that
  threshold shift is the sole cause of the default-threshold errors. Because
  the same held-out labels select and evaluate the cutoff, those improved
  accuracies are apparent/resubstitution results. The AUC supports ranking
  within the observed project but does not identify the source of that ranking.
- Prostate (PRAD) has the lowest observed LOCO AUC (0.950). A comparatively
  subtle tumor/normal transcriptional contrast is one possible explanation,
  not a conclusion established by this analysis.

Detailed LOCO artifacts are retained in the full development repository, not
inside the lightweight zip, under `cross-cancer-holdout/`:
`LOCO_REPORT.md`, `loco_report.html`, `loco_per_cancer_metrics.csv`,
`loco_pooled_summary.csv`, `loco_vs_within_comparison.csv`, and
`loco_predictions.pkl`. (The raw `loco_generalization*` figure/CSVs from the
original workbench run live in `from-workbench-loco/`.)

LOCO is a project-level holdout, not a causal transport experiment. Cancer type
is entangled with GDC project, tissue procurement, center, and batch. The design
does not remove those confounders and must not be described as proving
biological generalization to every unseen cancer type or deployment setting.

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
`external-validation/cptac_gdc/cptac_summary.csv`, and
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

The CPTAC, GTEx, and Toil metrics in this repository are committed historical
snapshots. Version 2.0.0 adds locked cohort manifests, semantic cache
fingerprints, content hashes, atomic cache publication, and run provenance, but
no post-fix live-network rerun was performed. These values therefore document
the earlier runs; they are not evidence that the corrected fetch/cache paths
have already reproduced the metrics live.

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
the scorer, workflow, explainer, adaptation scorer, and cancer-type predictor
before output files are written unless a reviewed tolerance or explicit
override is supplied. This prevents a malformed matched input column or sample
from being silently neutral-imputed.

## Experimental cross-platform cohort adaptation

The historical Toil benchmark found accuracy 0.935 after label-free cohort
standardization, compared with 0.515 for unadapted 0.5-threshold calls. This is
an exploratory, transductive result and does not prove that pipeline shift was
the only failure or that biological signal was preserved. Adaptation is
disabled by default (`--adapt none`). The `cohort_zscore` and `cohort_center`
modes are explicit opt-ins, depend on the other samples and class composition
in the batch, require at least `--min-samples` (default 20), and make scores
from separately adapted batches non-comparable. They are unsuitable for
single-sample inference and do not validate arbitrary non-GDC inputs.

The implementation is packaged as `cohort_adapt_score.py`; the historical
benchmark write-up lives under
`cross-platform-adaptation/CROSS_PLATFORM_ADAPTATION.md` in the full tree.

## Cancer-type classifier (tissue of origin)

This tumor-vs-normal model deliberately does not identify which cancer type a
sample is. A separate 17-class tissue-of-origin classifier ships alongside it,
reaching patient-held-out accuracy 0.930 (balanced accuracy 0.878, macro-F1
0.877). Its write-up, weights, and prediction entry point live in the full
development repository under `cancer-type-classifier/`; they are not part of
the lightweight tumor-vs-normal zip.

## Caveats

- This is a **pan-cancer** tumor-vs-normal classifier, not a cancer-type
  classifier — it distinguishes malignant from non-malignant tissue across
  cancer types. The LOCO result (mean AUC 0.994) is consistent with a shared
  signal in these data, but does not remove project/procurement/batch
  confounding or prove the proposed stromal/ECM mechanism.
- CPTAC-3 provides an external non-TCGA cohort check (AUC 0.989), but it is
  still GDC harmonized STAR-Counts. Non-GDC / cross-platform RNA-seq was tested
  with UCSC Xena Toil/GTEx and failed as a direct hard-call deployment target.
  The model should be refit or recalibrated for those pipelines.
- Solid-tissue-normal samples in TCGA are usually resected margin tissue,
  not tissue from cancer-free individuals — this is the correct
  tumor-vs-adjacent-normal contrast but is a different (easier) comparison
  than tumor-vs-healthy-population.
- No batch-effect correction across cancer types or sequencing centers was
  applied beyond TPM normalization. The high per-project AUCs may reflect a
  shared cancer-vs-adjacent-normal signal, technical/project structure, or
  both; this study does not separate those contributions.
- `tumor_probability` is the logistic model score, not clinical risk or a
  calibrated diagnostic probability. If a threshold is selected with
  `calibrate_threshold.py`, its reported metrics are same-sample
  apparent/resubstitution estimates; either class below 10 samples is warned.
- Literature agreement of highly weighted genes provides qualitative context,
  not causal mechanism proof or independent validation.

## Files

The lightweight release contains the deployable tumor-vs-normal scorer,
documentation, examples, validation checks, and summary external-validation
artifacts. Full development tree-only artifacts are listed separately below so
the extracted zip does not imply those paths are present locally.

### Lightweight bundle

- `external-validation/cptac_gdc/CPTAC_EXTERNAL_VALIDATION.md` — external CPTAC-3 smoke validation
- `external-validation/gtex_xena/GTEX_NORMAL_VALIDATION.md` — GTEx normal-tissue cross-platform check
- `external-validation/tcga_toil_xena/TCGA_TOIL_PIPELINE_CHECK.md` — Toil/RSEM pipeline-transfer check
- `test_metrics.csv` — held-out test metrics for all 3 models
- `per_cancer_type_performance.csv` — per-cancer-type AUC/accuracy breakdown
- `top_genes_xgboost.csv`, `top_genes_logreg.csv` — top 30 genes by each model
- `EXECUTIVE_SUMMARY.md`, `VERSION`, `RELEASE_METADATA.json`,
  `release_manifest.json`, and `SHA256SUMS.txt` — short handoff and release
  metadata/integrity files
- `USER_GUIDE.md`, `templates/` — practical input-preparation guide and CSV sketches
- `DATA_DICTIONARY.md` — stable input/output columns and JSON contract reference
- `TROUBLESHOOTING.md` — common install, input-QC, threshold, and release-integrity fixes
- `example_workflow_output/` — reference output from the one-command example workflow
- `deployable_lr_weights.npz` — pure NumPy export of the default LR scorer
- `run_tumor_normal_workflow.py` — one-command QC, scoring, calibration, explanations, and report
- `score_tumor_normal.py` and `cohort_adapt_score.py` — public
  scoring/adaptation entry points
- `tcga_rnaseq/` — shared dependency-light scoring core
- `check_environment.py` — runtime/package/file diagnostic and bundled self-test wrapper
- `audit_lightweight_dependencies.py` — release import and minimal dependency audit
- `audit_cli_entrypoints.py` — release CLI `--help` and shebang audit
- `audit_release_docs.py` — documentation and release-bundle reference audit
- `validate_output_contracts.py` — bundled CSV/JSON output contract validator
- `inspect_expression_input.py`, `model_qc_reference.json` — pre-scoring input compatibility QC
- `model_gene_metadata.csv` — model gene coefficients, scaling metadata, and known gene names
- `calibrate_threshold.py` — threshold calibration utility for labeled scored samples
- `explain_scores.py` — per-sample LR logit contribution explanations
- `validate_release_lite.py` — validates release folder, manifest, checksums, zip, and forbidden artifacts
- `validate_zip_bundle.py` — extracts the zip into a clean temp directory and runs acceptance
- `run_safety_tests.py` — verifies guardrails for invalid inputs and QC-fail workflow behavior
- `run_release_acceptance.py` — one-command environment, smoke, safety, and release-integrity checks

### Full development tree only

- `model_performance.png` and `feature_importance.png` — original figures.
- `cross-cancer-holdout/LOCO_REPORT.md`, `cross-cancer-holdout/loco_report.html`,
  `cross-cancer-holdout/loco_per_cancer_metrics.csv`,
  `cross-cancer-holdout/loco_pooled_summary.csv`,
  `cross-cancer-holdout/loco_vs_within_comparison.csv`, and
  `cross-cancer-holdout/loco_predictions.pkl` — detailed LOCO artifacts.
- `from-workbench-loco/loco_generalization.png`,
  `from-workbench-loco/loco_generalization.csv`, and
  `from-workbench-loco/loco_generalization_with_optimal_threshold.csv` — raw
  workbench LOCO figure/CSVs.
- `external-validation/validate_cptac_gdc.py`,
  `external-validation/validate_gtex_xena.py`, and
  `external-validation/validate_tcga_toil_xena.py` — data-fetching validation
  scripts; the lightweight zip ships their summary outputs.
- `cross-platform-adaptation/` — expanded adaptation benchmark write-up and CSVs.
- `cancer-type-classifier/` — separate tissue-of-origin classifier artifacts.
- `build_release_lite.py`, `export_lr_weights.py`, `export_qc_reference.py`,
  and `export_model_gene_metadata.py` — maintenance/regeneration helpers.
- `selected_files.csv` — full manifest of the 2,160 GDC file IDs / case IDs /
  cancer types / labels used.
- `RELEASE_ARTIFACTS.json` — generated sidecar next to the release zip with
  zip size and SHA256.
- `deployable_pipeline.pkl` and other training/checkpoint pickle artifacts are
  intentionally excluded from the public Git history and are not required by
  the lightweight scorer.
