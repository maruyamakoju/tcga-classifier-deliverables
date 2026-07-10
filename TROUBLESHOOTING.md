# Troubleshooting

Release: `v1.1.16-gdc-starcounts` (`2026-07-10`)

Start with:

```bash
python check_environment.py --self-test
python audit_lightweight_dependencies.py
python audit_cli_entrypoints.py
python run_smoke_tests.py
python run_safety_tests.py
python validate_zip_bundle.py tcga-tumor-normal-release-lite.zip
```

If these pass, most issues are input-specific rather than installation-specific.

## Installation issues

### `ModuleNotFoundError: No module named 'numpy'` or `pandas`

Install the lightweight runtime:

```bash
pip install -r requirements-light.txt
```

### Parquet input fails

Install `pyarrow`:

```bash
pip install pyarrow
```

CSV and TSV inputs do not require parquet support.

### Unsupported legacy scorer options

Default LR scoring does not need scikit-learn or pickle loading. Use:

```bash
python score_tumor_normal.py input.csv
```

`--use-pickle-lr` and `--model rf` are rejected in the public lightweight
release because the pickle/RF artifacts are intentionally not included.

## Input orientation and gene IDs

### QC says `no_model_genes_matched`

Likely causes:

- genes are rows, not columns
- columns are gene symbols instead of Ensembl IDs
- the first column was not parsed as sample IDs

Try:

```bash
python inspect_expression_input.py input.csv --transpose
```

If that fixes gene matching, run the workflow with `--transpose`.

### QC says `low_model_gene_match`

Check that gene IDs are Ensembl IDs and that the matrix contains most of the
2,000 model genes in `model_gene_metadata.csv`. Version suffixes are fine.

## Expression-scale problems

### QC says `expression_values_high` or `expression_values_too_large`

The model expects `log2(TPM+1)`. Raw TPM, raw counts, or integer count matrices
can be much too large.

If your data are TPM, convert with:

```text
log2(TPM + 1)
```

Do not apply this conversion directly to raw counts and assume the result is
GDC STAR-Counts-compatible.

### QC says `cohort_distribution_shift`

The input distribution differs from the model scaler. Common causes:

- non-GDC RNA-seq pipeline
- RSEM/Toil/GTEx/GEO source
- tissue mix far outside TCGA adjacent-normal contrast
- batch or normalization shift

Do not use hard calls until this is understood. Use labeled calibration data or
refit/recalibrate for that pipeline.

## Workflow behavior

### Workflow stops before writing `scores.csv`

This is expected when QC status is `FAIL` or when matched model-gene cells are
missing, non-numeric, `NaN`, or infinite. Review:

- `qc.json`
- `workflow_report.md`
- `manifest.json`

You can force scoring with `--allow-qc-fail`, but that is only for debugging.
For invalid matched expression values, fix the matrix first. If reviewed mean
imputation is intentional, use `--max-invalid-cell-fraction` to set an explicit
tolerance or `--allow-invalid-values` to downgrade the scorer stop to warnings.

### Scorer, explainer, or adaptation CLI says `invalid matched values`

The input has at least one matched model-gene cell that cannot be used as a
finite number. Common causes are blank cells, `NA` strings, spreadsheet export
markers, or infinite values from an upstream transform. The CLI reports example
genes and samples before refusing to write outputs.

### Many normal samples are called tumor

This is a known failure mode for cross-platform normal datasets such as
GTEx/Toil. Check:

```bash
python inspect_expression_input.py input.csv --expected-class normal --strict
```

If this warns, do not use hard calls without additional validation.

## Thresholds

### The recommended threshold is far from 0.5

This can happen with new tissues or new pipelines. It means the default 0.5
probability cutoff may not transfer. Use a labeled calibration set and report
the calibrated threshold with the results.

### AUC is high but accuracy is poor

Ranking can transfer while probability calibration fails. This was observed in
the TCGA Toil/RSEM check. Treat this as a threshold/pipeline-transfer problem,
not as validation of default hard calls.

## Release integrity

Validate a built release with:

```bash
python validate_release_lite.py --release-dir release-lite --zip tcga-tumor-normal-release-lite.zip
python validate_zip_bundle.py tcga-tumor-normal-release-lite.zip
```

This checks:

- manifest and SHA256 consistency
- zip/folder parity
- absence of large training artifacts
- absence of transient files such as `__pycache__`
