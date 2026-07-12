# User guide

Release: `v2.2.0-gdc-starcounts` (`2026-07-12`; public scoring-library API `3.0.0`)

This guide is for people who want to run the lightweight classifier on a new
expression matrix. For a short project-level handoff, read
`EXECUTIVE_SUMMARY.md` first.

## 1. Install

From the `release-lite/` folder:

Lightweight scoring supports Python 3.11 or newer. The hosted scoring
compatibility matrix exercises Python 3.11 and 3.13; exact model refitting is
separately pinned to the canonical Python 3.11 stack in `REPRODUCIBILITY.md`.

```bash
pip install -r requirements-light.txt
```

Then run the bundled tests:

```bash
python check_environment.py --self-test
python audit_lightweight_dependencies.py
python audit_cli_entrypoints.py
python audit_release_docs.py
python validate_output_contracts.py
python run_release_acceptance.py
python validate_zip_bundle.py ../tcga-tumor-normal-release-lite.zip --expected-sha256 <trusted-published-sha256>
python run_smoke_tests.py
python run_safety_tests.py
```

All should pass before scoring new data.
If installation or QC fails, see `TROUBLESHOOTING.md`.

The ZIP command requires a SHA-256 obtained from a trusted release channel
before it extracts or executes archive content. With no trusted digest, use
`--skip-acceptance` only for non-executing structural inspection.

## 2. Prepare input

The expression matrix must be:

- rows = samples
- columns = Ensembl gene IDs
- first column = sample identifier / row index
- values = `log2(TPM+1)`
- source scale = GDC STAR-Counts-style TPM

Accepted file types:

- `.csv`
- `.tsv` / `.txt`
- `.parquet`

Pickled pandas DataFrames are intentionally rejected by the public CLIs because
unpickling user-controlled files can execute code. Convert collaborator-provided
matrices to CSV, TSV, or Parquet before scoring.

Ensembl version suffixes are accepted. For example, both `ENSG00000123456` and
`ENSG00000123456.7` can match model genes.

Use `templates/input_matrix_template.csv` only as a format sketch. A real input
should contain most or all of the 2,000 model genes listed in
`model_gene_metadata.csv`.

## 3. Run the workflow

For unlabeled samples:

```bash
python run_tumor_normal_workflow.py input.csv
```

For a calibration or evaluation set with labels:

```bash
python run_tumor_normal_workflow.py input.csv --labels labels.csv
```

`labels.csv` should contain:

- `sample`
- `label`

Accepted labels include `tumor`, `normal`, `1`, and `0`. See
`templates/labels_template.csv`.

## 4. Read the output

The workflow creates a folder named `<input>_tumor_normal_workflow/` unless
`--output-dir` is supplied.

Start with:

- `workflow_report.md`: human-readable summary
- `qc.json`: detailed input compatibility checks
- `scores.csv`: `sample,tumor_probability,call`

`tumor_probability` is the model's unrounded logistic score. It is not clinical
risk and must not be presented as a calibrated diagnostic probability.

When labels are supplied:

- `thresholds.csv`: default and Youden's-J threshold metrics
- `calibration.json`: compact recommended-threshold summary

These threshold metrics are apparent/resubstitution estimates: the same
labeled samples both choose and evaluate the threshold. They are not an
independent validation result. The CLI warns when either class has fewer than
10 samples.

For model debugging:

- `explanations.csv`: top positive and negative per-gene LR logit contributions

The explanation file is not a causal biological interpretation.

## 5. Interpret QC

QC status can be:

- `PASS`: input is compatible with the current heuristic checks.
- `WARN`: do not trust hard calls until the warning is understood.
- `FAIL`: workflow stops before scoring unless `--allow-qc-fail` is used.

Common messages:

| Code | Meaning | Typical fix |
|---|---|---|
| `low_model_gene_match` | Too many model genes are absent | Check gene IDs, row/column orientation, Ensembl version handling |
| `no_model_genes_matched` | No model genes matched | Verify that genes are columns or use `--transpose` |
| `expression_values_high` | Values look too large for log2(TPM+1) | Check that data are not raw TPM/counts |
| `expression_values_too_large` | Values are far too large | Convert from TPM to log2(TPM+1), or stop if counts are raw |
| `nonfinite_or_missing_values` | Some model-gene cells are blank, non-numeric, `NaN`, or infinite | Fix the matrix before scoring, or explicitly review imputation |
| `cohort_distribution_shift` | Cohort differs from model scaler distribution | Check platform, normalization, tissue mix, and batch |
| `unexpected_tumor_calls` | Normal-expected cohort has many tumor calls | Check domain compatibility before hard calls |

A `PASS` result does not validate a new RNA-seq pipeline. It only says the
input did not trip the current guardrails.

The direct scoring, explanation, adaptation, and cancer-type prediction CLIs
also enforce the severe gene-coverage boundary by default. If fewer than 50% of
model genes match, they stop before writing outputs. Fix gene IDs, transpose
orientation, or matrix export first; use `--allow-low-gene-coverage` or a lower
`--min-model-gene-match-rate` only after reviewing the missing-gene imputation.

Matched model-gene cells that are missing, non-numeric, `NaN`, or infinite are
not silently accepted by the scorer, workflow, explainer, or adaptation scorer.
By default they stop before output files are written. Fix the input first; use
`--max-invalid-cell-fraction` or `--allow-invalid-values` only when you have
reviewed and accepted training mean imputation for those cells.

## 6. Thresholds

The default threshold is `0.5`. It is validated for GDC STAR-Counts-style
TCGA/CPTAC-like data, but it may not transfer to new tissues or pipelines.

If you have labeled samples, use the workflow with `--labels` and review the
recommended Youden's-J threshold. If the recommended threshold is extreme or QC
is `WARN`, treat that as evidence of domain or platform shift.

If label calibration fails after scoring, the workflow keeps `qc.json` and
`scores.csv`, writes `manifest.json` with
`status: stopped_after_calibration_error`, and records the label error in
`workflow_report.md`. Fix the labels and rerun before using threshold metrics.

## 7. Strong do-not-use cases

Do not use this release for direct hard calls on:

- Toil/RSEM matrices
- GTEx
- GEO or other arbitrary TPM matrices
- raw counts
- FPKM
- microarray
- single-cell RNA-seq
- spatial transcriptomics
- clinical diagnosis or patient management

These require additional validation, recalibration, or refitting. For a
foreign-pipeline cohort where the model's ranking (AUC) still holds but the
default 0.5 threshold does not, adaptation is still disabled by default
(`--adapt none`). `cohort_adapt_score.py --adapt cohort_zscore` is an explicit
experimental opt-in that re-centers the batch on its own mean/variance before
scoring; see
`cross-platform-adaptation/CROSS_PLATFORM_ADAPTATION.md` for when this
applies and its limits (it needs an internal tumor/normal contrast in the
cohort, so it does not help an all-normal or all-tumor batch). Adaptation is
transductive and composition-dependent, requires at least `--min-samples`
(default 20), changes a sample's score when its batch changes, and makes scores
from separately adapted batches non-comparable.

## 8. Quick reference

```bash
python check_environment.py --self-test
python audit_lightweight_dependencies.py
python audit_cli_entrypoints.py
python audit_release_docs.py
python validate_output_contracts.py
python run_release_acceptance.py
python validate_zip_bundle.py ../tcga-tumor-normal-release-lite.zip --expected-sha256 <trusted-published-sha256>
python run_tumor_normal_workflow.py input.csv --labels labels.csv
python inspect_expression_input.py input.csv -o input.qc.json
python score_tumor_normal.py input.csv -o scores.csv --threshold 0.5
python calibrate_threshold.py scores.csv labels.csv
python explain_scores.py input.csv --top-n 10
python cohort_adapt_score.py input.csv --adapt cohort_zscore
```
