# Troubleshooting

Release: `v2.3.0-gdc-starcounts` (`2026-07-12`; public scoring-library API `3.0.0`)

Start with:

```bash
python check_environment.py --self-test
python audit_lightweight_dependencies.py
python audit_cli_entrypoints.py
python run_smoke_tests.py
python run_safety_tests.py
python validate_zip_bundle.py tcga-tumor-normal-release-lite.zip --expected-sha256 <trusted-published-sha256>
```

If these pass, most issues are input-specific rather than installation-specific.
The lightweight scoring runtime supports Python 3.11 or newer; CI exercises
3.11 and 3.13. Exact model refitting uses the canonical Python 3.11 training
stack documented in `REPRODUCIBILITY.md`.

### CLI says `expression matrix file not found`

The input path does not point to a readable expression matrix file. Check the
path, working directory, and file extension before rerunning. The public CLIs
expect `.csv`, `.tsv`, `.txt`, or `.parquet` expression matrices by default;
pickled expression inputs are blocked unless an internal trusted caller opts in.

### CLI rejects a `.pkl` expression matrix

This is intentional. Pickle loading can execute code. Convert untrusted or
collaborator-provided data to CSV, TSV, or Parquet. The low-level
`read_matrix(..., allow_pickle=True)` path is only for an already trusted local
artifact; it is not a way to make an unknown pickle safe.

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

`--use-pickle-lr` and `--model rf` are both rejected in the public lightweight
release because the pickle/RF artifacts are intentionally not included, but
with different error text: `--use-pickle-lr` reaches
`score_tumor_normal.py` and fails with an explanatory
`legacy pickle/RF scoring is not available` message; `--model` only accepts
`lr` in this release, so `--model rf` is rejected by argument parsing itself
with a generic `invalid choice: 'rf'` message before the script runs.

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

Direct scoring, explanation, adaptation, and cancer-type prediction CLIs refuse
to write outputs when fewer than 50% of model genes match by default. This is
usually an orientation or identifier problem, not a condition to override. Use
`--allow-low-gene-coverage` or lower `--min-model-gene-match-rate` only after
reviewing the missing-gene imputation.

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

### Workflow stops with `stopped_after_calibration_error`

Scoring finished, but the label CSV could not be used for threshold calibration.
The workflow preserves `qc.json` and `scores.csv`, records the calibration error
in `manifest.json` and `workflow_report.md`, and does not write threshold or
calibration outputs. Check label sample IDs, duplicates, supported label values,
and `--min-match-fraction`, then rerun.

### Output contract audit reports sample ID errors

`validate_output_contracts.py` treats sample IDs in scored outputs and
`example_labels.csv` as trimmed identifiers. Empty IDs, leading/trailing
whitespace, and duplicate IDs after trimming are contract failures. Fix the
source CSV so sample identifiers are exact and unique before rebuilding outputs.

### Scorer, explainer, or adaptation CLI says `invalid matched values`

The input has at least one matched model-gene cell that cannot be used as a
finite number. Common causes are blank cells, `NA` strings, spreadsheet export
markers, or infinite values from an upstream transform. The CLI reports example
genes and samples before refusing to write outputs.

### Scorer, explainer, or adaptation CLI says `low model-gene coverage`

The input matched too few model genes for direct output writing. Check that
genes are columns, use `--transpose` where supported, and verify that identifiers
are Ensembl gene IDs rather than gene symbols. Review `inspect_expression_input.py`
before using any override.

### CLI reports an invalid model artifact

The v3 API validates every model array, shape, gene identifier, scale, and
finite value before scoring. Do not patch around this error or partially load
the artifact. Restore a verified `deployable_lr_weights.npz` and rerun the
release checks.

### CLI rejects output paths or says paths collide

Inputs, model files, labels, and outputs must not resolve to the same path, and
multiple public outputs must remain distinct. Choose a fresh output file or
directory. Each file is published atomically, so an error should not leave a
partially written file. The workflow as a whole is not all-or-nothing: it may
retain documented valid earlier outputs, removes stale downstream artifacts at
startup, records stop states, and writes its manifest last.

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

Threshold-table and JSON metrics are apparent/resubstitution estimates because
the same labeled samples select and evaluate the cutoff. Do not report them as
independent test performance. A warning that either class has fewer than 10
samples is a precision warning, not a successful validation signal.

### AUC is high but accuracy is poor

Ranking can transfer while probability calibration fails. This was observed in
the TCGA Toil/RSEM check. Treat this as a threshold/pipeline-transfer problem,
not as validation of default hard calls.

### Adaptation fails because the cohort is too small

Adaptation defaults to `none`. Experimental `cohort_zscore` and
`cohort_center` modes require at least `--min-samples` (default 20). Do not
lower this just to score a single sample or tiny batch: adaptation is
transductive and composition-dependent, a sample's result changes with its
batch, and scores from separately adapted batches are not comparable.

### `tumor_probability` looks like a clinical probability

Despite the historical column name, it is the logistic model score. It is not
clinical risk and is not a calibrated diagnostic probability. Interpret hard
calls only within the validated GDC STAR-Counts research boundary.

## External-validation maintenance

### Validator rejects `unversioned` or tries no network in offline mode

Live external validation requires a concrete provider snapshot in
`--source-revision`. `--offline` / `--cache-only` is deliberately fail-closed:
if the Parquet cache fingerprint, hash, or axes do not match, it stops instead
of downloading replacement data.

### CPTAC locked manifest has no provider MD5

The committed CPTAC cohort manifest under
`external-validation/cptac_gdc/sampled_manifest.csv` is a historical pre-MD5
snapshot. It may be
used only with `--offline` and a matching already-trusted cache. For a new live
download, pass `--refresh-manifest` with the locked cohort; the validator keeps
the exact file IDs, checks that provider identity fields did not change, and
requires the refreshed GDC MD5. Do not invent or copy checksums from unrelated
files.

### Protected input collides with a managed output or cache

Use a fresh live `--out-dir` outside the model and locked-manifest paths. The
validators reject aliases, symlinks/hardlinks, and protected files inside a
mutable cache directory before network access or writes. Derived results are
staged and the run manifest is published last, so an ordinary failed run keeps
the previous complete result generation unchanged.

## Release integrity

Validate a built release with:

```bash
python validate_release_lite.py --release-dir release-lite --zip tcga-tumor-normal-release-lite.zip
python validate_zip_bundle.py tcga-tumor-normal-release-lite.zip --expected-sha256 <trusted-published-sha256>
```

This checks:

- manifest and SHA256 consistency
- zip/folder parity
- absence of large training artifacts
- absence of transient files such as `__pycache__`

`validate_zip_bundle.py` will not extract or execute a downloaded archive for
acceptance without `--expected-sha256` matching a trusted published digest.
If you do not yet have a trusted digest, use `--skip-acceptance` only for
non-executing structural inspection; obtain the digest from a trusted release
channel before running bundled code.
