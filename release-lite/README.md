# TCGA tumor-vs-normal classifier — deliverables

[![CI](https://github.com/maruyamakoju/tcga-classifier-deliverables/actions/workflows/ci.yml/badge.svg)](https://github.com/maruyamakoju/tcga-classifier-deliverables/actions/workflows/ci.yml)
[![Release](https://img.shields.io/github/v/release/maruyamakoju/tcga-classifier-deliverables?display_name=tag)](https://github.com/maruyamakoju/tcga-classifier-deliverables/releases/tag/v2.1.0-gdc-starcounts)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Citation](https://img.shields.io/badge/citation-CITATION.cff-blue.svg)](CITATION.cff)

Release: `v2.1.0-gdc-starcounts` (`2026-07-12`; public scoring-library API
`3.0.0`). For a single guided path
through the public lightweight bundle, start with `INDEX.md`. Otherwise start
with `EXECUTIVE_SUMMARY.md` if you need a short
handoff/readout, or `USER_GUIDE.md` if you are preparing a new input matrix.

A pan-cancer RNA-seq classifier that calls a sample **tumor** vs **normal**, trained on
2,160 TCGA samples across 17 cancer types. Best model: logistic regression on 2,000
genes. Validated three ways: patient-held-out AUC 0.997, 5-fold grouped CV 0.996±0.003,
leave-one-**cancer-type**-out macro-mean AUC 0.994 (pooled 0.988),
and an external CPTAC-3/GDC STAR-Counts smoke test with AUC 0.989.
LOCO measures held-out project discrimination, but it does not remove
GDC/project/procurement/batch confounding and is not proof of causal biological
generalization.

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
python validate_zip_bundle.py tcga-tumor-normal-release-lite.zip --expected-sha256 <trusted-published-sha256>
python score_tumor_normal.py example_input.csv          # -> example_input.scored.csv
python score_tumor_normal.py --self-test                # verify bundled example, no sklearn needed
python inspect_expression_input.py example_input.csv     # QC gene coverage, scale, and shift
python score_tumor_normal.py expr.csv -o calls.csv --threshold 0.5
python calibrate_threshold.py calls.csv labels.csv       # choose a threshold from labeled samples
python explain_scores.py expr.csv --top-n 10             # per-sample LR gene contributions
python cohort_adapt_score.py expr.csv --adapt cohort_zscore  # explicit experimental cohort adaptation
```

- **Input:** rows = samples, columns = genes (Ensembl IDs, with or without the `.version`
  suffix), values = **log2(TPM+1)** on the GDC STAR-Counts scale. Formats:
  `.csv .tsv .txt .parquet`. Pickled expression matrices are intentionally
  rejected by the public CLIs. Add `--transpose` if genes are rows (accepted by
  `score_tumor_normal.py`, `explain_scores.py`, `inspect_expression_input.py`, and
  `run_tumor_normal_workflow.py`; `cohort_adapt_score.py` currently expects
  samples as rows). Genes missing from the input are filled with the training
  mean (neutral after standardization) and reported.
- **Cross-platform cohorts:** adaptation defaults to `none`. The
  `--adapt cohort_zscore` and `cohort_center` modes are explicit experimental
  opt-ins: they are transductive, depend on the other samples and class mixture
  in the batch, require at least `--min-samples` (default 20), and produce scores
  that are not comparable across separately adapted batches. They do not
  validate arbitrary foreign pipelines; see
  `cross-platform-adaptation/CROSS_PLATFORM_ADAPTATION.md`.
- **Output CSV:** `sample, tumor_probability, call`. `tumor_probability` is the
  model's logistic score, not clinical risk or a calibrated diagnostic
  probability.
- **One-command workflow:** `run_tumor_normal_workflow.py` writes `qc.json`,
  `scores.csv`, optional `thresholds.csv` / `calibration.json`,
  `explanations.csv`, `manifest.json`, and `workflow_report.md` into one output
  directory.
- **Input QC:** `inspect_expression_input.py` writes a JSON report with gene match rate,
  expression range, standardized distribution-shift metrics, and score summary. Run it
  before scoring when the matrix came from a new pipeline or collaborator.
- **Low model-gene coverage:** direct scoring, explanation, adaptation, and
  cancer-type prediction CLIs now refuse to write outputs when fewer than 50%
  of model genes match by default. Fix gene IDs/orientation first; use
  `--allow-low-gene-coverage` only after reviewing mean imputation.
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
Youden's-J threshold. Reported calibration metrics are apparent/resubstitution
estimates computed on the same samples used to choose the threshold, not
independent performance estimates. The CLI warns when either class has fewer
than 10 samples.

### Explanations
`explain_scores.py` reports the top positive and negative per-gene contributions to the
logistic-regression logit for each sample. Use it for model debugging and sanity checks,
not as a causal biological explanation. `model_gene_metadata.csv` lists all 2,000 model
genes, coefficients, training means/scales, and the direction implied by high expression.

Literature consistency of selected genes is qualitative context only. It does
not establish a causal mechanism for the classifier or its predictions.

### Running environment
Python 3.11 or newer is supported; CI exercises Python 3.11 and 3.13. For
default LR scoring, use `requirements-light.txt` (NumPy + pandas only, with pyarrow
for parquet input). Use `requirements.txt` or `environment.yml` only when retraining,
running external validation scripts, or doing full-artifact maintenance outside the
public lightweight bundle.

That broad runtime statement applies to lightweight scoring, not exact model
refitting. Strict shipped-model reproduction uses `requirements-training.txt`
on Python 3.11 with NumPy 1.26.4, pandas 2.3.3, SciPy 1.15.3, and
scikit-learn 1.8.0. Python 3.13 is covered only for lightweight scoring
acceptance; an investigated Python 3.13/scikit-learn 1.9 refit showed
coefficient and out-of-fold drift.
Treat such drift as an environment mismatch; do not relax committed golden or
parity tolerances to make a different stack pass.

Run `python check_environment.py --self-test` after installation. If it or the
workflow QC reports WARN/FAIL, start with `TROUBLESHOOTING.md`.

This model was evaluated internally on TCGA/GDC-style log2(TPM+1)
tumor-vs-adjacent-normal RNA-seq contrasts and has a historical external CPTAC-3
smoke check within the same GDC harmonized STAR-Counts ecosystem. Non-GDC /
cross-platform RNA-seq was explored via UCSC Xena Toil/GTEx and should be treated
as **not directly compatible** without pipeline-specific refitting and independent
threshold evaluation.

The committed CPTAC/GTEx/Toil metrics are historical snapshots. Version 2.0.0
adds locked cohort manifests, cache fingerprints, hashes, atomic cache writes,
and run manifests, but no post-fix live-network validation rerun was performed;
do not describe those committed numbers as a newly reproduced live run.
External validators now require a concrete provider revision for live access,
offer explicit `--offline` / `--cache-only` operation, strictly bind locked
cohort semantics, and publish derived outputs as a staged set with the manifest
last. The historical CPTAC locked manifest has no provider MD5 and is therefore
eligible only for offline reuse with an already-valid historical cache; a new
download requires reviewed GDC metadata refresh and provider MD5 verification.

### What changed in 2.0.0

The public `tcga_rnaseq` API is now 3.0.0 and intentionally breaks unsafe
legacy behavior: model/sample contracts fail closed, public output paths cannot
collide with inputs or each other, and each output file is written atomically.
The workflow is not an all-or-nothing transaction: it removes stale downstream
files, records documented stop states, and writes the manifest last. The
development tree also adds separated lightweight/external/dev dependency
profiles, Ruff and full pytest CI, cross-platform release acceptance, canonical
ZIP metadata, and a non-mutating deterministic release drift check.

For a downloaded release ZIP, first obtain its SHA-256 from a trusted release
channel and pass it with `--expected-sha256`. If no trusted digest is available,
use `--skip-acceptance` only for non-executing structural inspection; the
validator will not extract or execute the archive.

## Files

### Lightweight bundle

- `INDEX.md` — guided single reading path through the public bundle.
- `EXECUTIVE_SUMMARY.md`, `USER_GUIDE.md`, `DATA_DICTIONARY.md`,
  `TROUBLESHOOTING.md`, `MODEL_CARD.md`, `REPORT.md`, `REPRODUCIBILITY.md`,
  and `RELEASE_NOTES.md` — user-facing documentation.
- `VERSION`, `RELEASE_METADATA.json`, `release_manifest.json`, and
  `SHA256SUMS.txt` — release metadata and integrity files.
- `run_tumor_normal_workflow.py`, `score_tumor_normal.py`,
  `inspect_expression_input.py`, `calibrate_threshold.py`, and
  `explain_scores.py` — deployable command-line tools.
- `deployable_lr_weights.npz`, `model_qc_reference.json`, and
  `model_gene_metadata.csv` — model payload and QC/gene metadata.
- `tcga_rnaseq/` — shared core library used by the scoring entry points.
- `check_environment.py`, `run_smoke_tests.py`, `run_safety_tests.py`,
  `run_release_acceptance.py`, `validate_release_lite.py`,
  `validate_zip_bundle.py`, `validate_output_contracts.py`,
  `audit_lightweight_dependencies.py`, `audit_cli_entrypoints.py`, and
  `audit_release_docs.py` — bundled validation and acceptance checks.
- `requirements-light.txt`, `example_input.csv`, `example_output.csv`,
  `example_labels.csv`, `example_workflow_output/`, and `templates/` —
  runnable examples and input templates.
- `external-validation/` — bundled summary reports and CSVs for CPTAC/GDC and
  Toil/GTEx boundary checks.
- `LICENSE`, `NOTICE.md`, `.zenodo.json`, `codemeta.json`, and `CITATION.cff`
  — license, attribution, and citation metadata.

### Full development tree only

The full repository, not the lightweight zip, also contains maintenance and
historical-analysis assets such as `build_release_lite.py`,
`audit_publication_readiness.py`, `audit_github_repository.py`,
`requirements.txt`, `requirements-training.txt`, `environment.yml`, `tests/`, `cross-cancer-holdout/`,
`from-workbench-loco/`, `cross-platform-adaptation/`, and
`cancer-type-classifier/`. Full training/checkpoint artifacts such as
`model_lr.pkl`, `model_rf.pkl`, `feature_selection.pkl`,
`X_full_filtered.pkl`, `y_full.pkl`, and `deployable_pipeline.pkl` are
intentionally excluded from the public Git history; the lightweight release
does not require them. `RELEASE_ARTIFACTS.json` is a sidecar next to the zip in
the source tree and GitHub release metadata, not a file inside the extracted
bundle.

## Citation and license

Use `CITATION.cff`, `.zenodo.json`, and `codemeta.json` for software citation
metadata. Repository code and project-authored documentation are MIT licensed;
see `LICENSE`. Third-party source datasets remain subject to their original
provider terms; see `NOTICE.md`.
