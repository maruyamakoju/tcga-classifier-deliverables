# TCGA tumor-vs-normal classifier — deliverables

[![CI](https://github.com/maruyamakoju/tcga-classifier-deliverables/actions/workflows/ci.yml/badge.svg)](https://github.com/maruyamakoju/tcga-classifier-deliverables/actions/workflows/ci.yml)
[![Release](https://img.shields.io/github/v/release/maruyamakoju/tcga-classifier-deliverables?display_name=tag)](https://github.com/maruyamakoju/tcga-classifier-deliverables/releases/tag/v1.1.6-gdc-starcounts)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Citation](https://img.shields.io/badge/citation-CITATION.cff-blue.svg)](CITATION.cff)

Release: `v1.1.6-gdc-starcounts` (`2026-07-09`). For a single guided path
through the whole deliverable (base model → generalization → external
validation → cross-platform adaptation → cancer-type classifier), start with
`INDEX.md`. Otherwise start with `EXECUTIVE_SUMMARY.md` if you need a short
handoff/readout, or `USER_GUIDE.md` if you are preparing a new input matrix.

A pan-cancer RNA-seq classifier that calls a sample **tumor** vs **normal**, trained on
2,160 TCGA samples across 17 cancer types. Best model: logistic regression on 2,000
genes. Validated three ways: patient-held-out AUC 0.997, 5-fold grouped CV 0.997±0.003,
leave-one-**cancer-type**-out macro-mean AUC 0.994 (pooled 0.988) (generalizes to cancer types never trained on),
and an external CPTAC-3/GDC STAR-Counts smoke test with AUC 0.989.

Important boundary: the deployable model is **GDC STAR-Counts-scale specific**. A
cross-platform UCSC Xena Toil/RSEM check showed strong ranking on sampled TCGA Toil
samples (AUC 0.992) but severe threshold shift at the default 0.5 cutoff, and GTEx
normal tissues were not safely callable without refitting/recalibration. Do not apply the
bundled model directly to Toil/RSEM, GTEx, GEO, or other non-GDC pipelines as hard calls.

## Score new samples — `score_tumor_normal.py`

```bash
python run_tumor_normal_workflow.py example_input.csv --labels example_labels.csv
python check_environment.py --self-test          # verify runtime, required files, and bundled example
python audit_lightweight_dependencies.py         # verify lightweight runtime imports stay minimal
python audit_cli_entrypoints.py                  # verify release CLI --help entry points
python audit_release_docs.py                     # check docs and command references
python validate_output_contracts.py              # check bundled output schemas
python run_release_acceptance.py                 # run environment, smoke, safety, and release checks
python validate_zip_bundle.py tcga-tumor-normal-release-lite.zip  # clean zip extraction acceptance
python score_tumor_normal.py example_input.csv          # -> example_input.scored.csv
python score_tumor_normal.py --self-test                # verify bundled example, no sklearn needed
python inspect_expression_input.py example_input.csv     # QC gene coverage, scale, and shift
python score_tumor_normal.py expr.csv -o calls.csv --threshold 0.5
python calibrate_threshold.py calls.csv labels.csv       # choose a threshold from labeled samples
python explain_scores.py expr.csv --top-n 10             # per-sample LR gene contributions
```

- **Input:** rows = samples, columns = genes (Ensembl IDs, with or without the `.version`
  suffix), values = **log2(TPM+1)** on the GDC STAR-Counts scale. Formats: `.csv .tsv
  .parquet .pkl`. Add `--transpose` if genes are rows. Genes missing from the input are
  filled with the training mean (neutral after standardization) and reported.
- **Output CSV:** `sample, tumor_probability, call`.
- **One-command workflow:** `run_tumor_normal_workflow.py` writes `qc.json`,
  `scores.csv`, optional `thresholds.csv` / `calibration.json`,
  `explanations.csv`, `manifest.json`, and `workflow_report.md` into one output
  directory.
- **Input QC:** `inspect_expression_input.py` writes a JSON report with gene match rate,
  expression range, standardized distribution-shift metrics, and score summary. Run it
  before scoring when the matrix came from a new pipeline or collaborator.
- **Invalid matched values:** matched model-gene cells that are missing, non-numeric,
  `NaN`, or infinite now stop scoring, workflow, explanation, adaptation, and
  cancer-type prediction CLIs by default before outputs are written.
  Fix the input, set a reviewed tolerance with `--max-invalid-cell-fraction`, or pass
  `--allow-invalid-values` only when mean imputation is intentional.
- **Runnable example:** `example_input.csv` (5 real samples) → `example_output.csv`
  (first 3 tumor at p>0.99, last 2 normal at p<0.06 — matches their true labels).
- **Default scorer:** pure NumPy logistic regression from `deployable_lr_weights.npz`
  (small, no scikit-learn pickle warning). The public lightweight CLI does not expose
  legacy pickle/RF scoring because those artifacts are intentionally excluded from the
  public Git history.

### Threshold calibration (important for a new tissue)
Ranking (AUC) transfers across cancer types, but the fixed **0.5 threshold does not**
always transfer to a tissue the model never trained on — prostate and liver tumors were
under-called at 0.5 despite AUC ≈ 0.95–1.0. If you score a genuinely new tissue and have
a few labeled samples, choose a cutoff on them and pass `--threshold`. See
`cross-cancer-holdout/` for the per-type calibration analysis. Use
`calibrate_threshold.py` with a scored CSV and a `sample,label` CSV to compute a
Youden's-J threshold.

### Explanations
`explain_scores.py` reports the top positive and negative per-gene contributions to the
logistic-regression logit for each sample. Use it for model debugging and sanity checks,
not as a causal biological explanation. `model_gene_metadata.csv` lists all 2,000 model
genes, coefficients, training means/scales, and the direction implied by high expression.

### Running environment
For default LR scoring, use `requirements-light.txt` (NumPy + pandas only, with pyarrow
for parquet input). Use `requirements.txt` or `environment.yml` only when retraining,
running external validation scripts, or doing full-artifact maintenance outside the
public lightweight bundle.

Run `python check_environment.py --self-test` after installation. If it or the
workflow QC reports WARN/FAIL, start with `TROUBLESHOOTING.md`.

This model is validated for TCGA/GDC-style log2(TPM+1) tumor-vs-adjacent-normal RNA-seq
contrasts and now has an external CPTAC-3 smoke validation within the same GDC
harmonized STAR-Counts ecosystem. Non-GDC / cross-platform RNA-seq has been tested via
UCSC Xena Toil/GTEx and should be treated as **not directly compatible** without
pipeline-specific refitting or threshold calibration.

## Files
- `INDEX.md` — guided single reading path through the whole deliverable
- `release-lite/` / `tcga-tumor-normal-release-lite.zip` — lightweight deployment bundle
- `EXECUTIVE_SUMMARY.md` — short handoff summary for reviewers/users
- `USER_GUIDE.md` — practical input preparation, QC interpretation, and workflow guide
- `DATA_DICTIONARY.md` — stable input/output columns and JSON contract reference
- `TROUBLESHOOTING.md` — fixes for install, input, QC, threshold, and release-integrity issues
- `VERSION`, `RELEASE_METADATA.json` — release version and metadata included in the bundle
- `RELEASE_ARTIFACTS.json` — generated sidecar next to the zip with artifact size and SHA256
- `run_tumor_normal_workflow.py` — one-command QC, scoring, calibration, explanations, and report
- `score_tumor_normal.py` — scoring CLI (this tool)
- `deployable_lr_weights.npz` — default pure NumPy LR scorer weights
- `inspect_expression_input.py` / `model_qc_reference.json` — input compatibility QC
- `model_gene_metadata.csv` — 2,000 model genes with LR coefficients and scaling metadata
- `check_environment.py` — runtime/package/file/self-test diagnostic
- `audit_lightweight_dependencies.py` — release import and `requirements-light.txt` dependency audit
- `audit_cli_entrypoints.py` — release CLI `--help` and shebang audit
- `audit_release_docs.py` — documentation and release-bundle reference audit
- `audit_publication_readiness.py` — public-release audit for secrets, large blobs,
  line endings, and release metadata consistency
- `audit_github_repository.py` — hosted GitHub settings and release-asset audit
- `validate_output_contracts.py` — bundled CSV/JSON output contract validator
- `calibrate_threshold.py` — choose a cutoff from labeled scored samples
- `explain_scores.py` — per-sample LR contribution explanations
- `export_lr_weights.py` — regenerate `deployable_lr_weights.npz` from the pickle pipeline
- `export_qc_reference.py` — regenerate `model_qc_reference.json` from validation matrices
- `export_model_gene_metadata.py` — regenerate `model_gene_metadata.csv`
- `run_smoke_tests.py` — local release sanity check including the one-command workflow
- `run_safety_tests.py` — negative-path guardrail tests for QC and workflow failure behavior
- `run_release_acceptance.py` — one-command environment, smoke, safety, and release-integrity check
- `validate_zip_bundle.py` — extract zip into a clean temp directory and run acceptance
- `build_release_lite.py` — rebuild `release-lite/`, `SHA256SUMS.txt`, and the zip archive
- `validate_release_lite.py` — verify checksums, manifest, zip contents, and absence of large training artifacts
- `release_manifest.json` — generated inside `release-lite/` with file hashes and sizes
- `MODEL_CARD.md`, `RELEASE_NOTES.md` — deployment boundary and release summary
- `requirements-light.txt`, `requirements.txt`, `environment.yml`, `REPRODUCIBILITY.md` — reproducible setup
- `example_input.csv` / `example_output.csv` / `example_labels.csv` — runnable demo + expected result + labels
- `example_workflow_output/` — reference output from the one-command demo workflow
- `templates/` — minimal input and label CSV format templates
- `REPORT.md` — full methods/results (now includes the leave-one-cancer-type-out section)
- `model_performance.png`, `feature_importance.png` — original figures
- `test_metrics.csv`, `per_cancer_type_performance.csv`, `top_genes_*.csv`, `selected_files.csv`
- `cross-cancer-holdout/` — generalization test: `LOCO_REPORT.md`, `loco_report.html`
  (interactive figure + table), per-cancer metrics CSVs, pooled predictions
- `external-validation/` — CPTAC-3 external smoke validation script, manifest, cached
  selected-gene matrix, predictions, threshold sweep, and report; GTEx/Toil and
  TCGA/Toil cross-platform checks
- `from-workbench-loco/` — the Claude Science workbench's own LOCO run
  (`loco_generalization.png` + optimal-threshold CSVs); agrees with `cross-cancer-holdout/`
- `cross-platform-adaptation/` — label-free cohort standardization that restores
  Toil accuracy 0.515→0.935; `CROSS_PLATFORM_ADAPTATION.md` + benchmark CSVs
- `cancer-type-classifier/` — separate 17-class tissue-of-origin classifier
  (`CANCER_TYPE_CLASSIFIER.md`, `predict_cancer_type.py`, weights, and metrics)
- `tcga_rnaseq/` — shared core library (I/O, gene alignment, scoring, metrics)
  used by the scoring entry points
- `tests/` — pytest suite for core units and numerical reproducibility
- `.zenodo.json`, `codemeta.json`, `CITATION.cff` — machine-readable citation and
  software metadata
- full training/checkpoint artifacts such as `model_lr.pkl`, `model_rf.pkl`,
  `feature_selection.pkl`, `X_full_filtered.pkl`, `y_full.pkl`, and
  `deployable_pipeline.pkl` are intentionally excluded from the public Git history; the
  lightweight release does not require them

## Citation and license

Use `CITATION.cff`, `.zenodo.json`, and `codemeta.json` for software citation
metadata. Repository code and project-authored documentation are MIT licensed;
see `LICENSE`. Third-party source datasets remain subject to their original
provider terms; see `NOTICE.md`.
