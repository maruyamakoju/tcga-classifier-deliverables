# Reproducibility notes

Release: `v1.1.9-gdc-starcounts` (`2026-07-09`)

## Recommended scoring environment

Default logistic-regression scoring uses `deployable_lr_weights.npz`, a pure
NumPy export of the fitted model. It does not need scikit-learn:

```bash
pip install -r requirements-light.txt
python check_environment.py --self-test
python audit_lightweight_dependencies.py
python audit_cli_entrypoints.py
python audit_release_docs.py
python validate_output_contracts.py
python validate_zip_bundle.py tcga-tumor-normal-release-lite.zip
```

For retraining, external-validation maintenance, or regenerating lightweight
weights from a full local training-artifact checkout, use one of:

```bash
pip install -r requirements.txt
```

or:

```bash
conda env create -f environment.yml
conda activate tcga-tumor-normal
```

Full pickle artifacts are intentionally excluded from the public Git history
and are not needed for the lightweight release. The public scoring CLI uses
only the pure NumPy logistic-regression weights in `deployable_lr_weights.npz`;
legacy pickle/RF CLI modes are not exposed in the lightweight bundle.

## Smoke test

Run the bundled self-test after installing dependencies:

```bash
python check_environment.py --self-test
python audit_lightweight_dependencies.py
python audit_cli_entrypoints.py
python audit_release_docs.py
python validate_output_contracts.py
python score_tumor_normal.py --self-test
```

Expected behavior:

- 5 samples scored from `example_input.csv`
- 2000/2000 model genes matched
- maximum probability delta vs `example_output.csv` is 0
- final line: `PASS: bundled example reproduces expected calls`

For the full lightweight release sanity check, run:

```bash
python run_release_acceptance.py
python run_smoke_tests.py
```

`run_release_acceptance.py` wraps the environment check, lightweight dependency
audit, CLI audit, documentation audit, output-contract validation, smoke test,
safety test, and release-integrity validator when a built release is present.
For a focused smoke-only pass, run `run_smoke_tests.py`. This checks scoring, input QC,
threshold calibration, per-sample explanations, and the one-command workflow.
To test negative-path guardrails alone, run:

```bash
python run_safety_tests.py
```

This verifies that invalid thresholds/top-N values fail, unsupported legacy
pickle/RF scorer options fail clearly, invalid matched expression values stop
score, explanation, and adaptation outputs before files are written unless
explicitly allowed, QC rejects inputs with no model genes or raw-count-like
values, and the workflow stops before scoring when QC status is FAIL.

For the development regression suite, run pytest through the active Python
interpreter:

```bash
python -m pytest -q -rs
```

The default suite checks core unit behavior, stable shipped summary artifacts,
patient-disjoint split metadata, lightweight external-validation artifacts, and
golden-number reproduction where the required matrices are bundled. The
full-data reproduction tests skip unless `TCGA_FEATURES` points to an exported
`X_full.npy` with sibling `X_genes.npy` and `X_samples.npy`.

For common installation, input-QC, threshold, and release-integrity failures,
see `TROUBLESHOOTING.md`.

## Release build

To rebuild the lightweight deployment bundle from the full deliverables folder:

```bash
python build_release_lite.py --smoke
```

This recreates `release-lite/`, writes `release-lite/release_manifest.json`,
writes `release-lite/SHA256SUMS.txt`, rebuilds
`tcga-tumor-normal-release-lite.zip`, validates the folder/manifest/checksums
and zip, and runs the bundled smoke test inside `release-lite/`.
The `--smoke` build also runs `run_safety_tests.py` inside `release-lite/`.
It also writes the sidecar `RELEASE_ARTIFACTS.json` with the zip byte size and
SHA256.

To verify the zip as a standalone artifact, extract it into a temporary clean
directory and run acceptance there:

```bash
python validate_zip_bundle.py tcga-tumor-normal-release-lite.zip
```

To validate an already-built bundle without rebuilding:

```bash
python validate_release_lite.py --release-dir release-lite --zip tcga-tumor-normal-release-lite.zip
python validate_zip_bundle.py tcga-tumor-normal-release-lite.zip
```

The default self-test uses `deployable_lr_weights.npz` and should not emit a
scikit-learn pickle warning. The current local Windows Python environment used
during cleanup was:

- Python 3.11.9
- pandas 2.3.3
- scikit-learn 1.8.0
- numpy 1.26.4
- scipy 1.15.3

The lightweight release does not load scikit-learn pickle artifacts during
scoring. For full-artifact maintenance such as regenerating `deployable_lr_weights.npz`,
prefer the pinned `requirements.txt` / `environment.yml`.

## Input contract

`score_tumor_normal.py` expects:

- rows = samples, columns = genes
- Ensembl gene IDs, with or without version suffixes
- values = log2(TPM+1) on the GDC STAR-Counts scale
- `.csv`, `.tsv`, `.txt`, `.parquet`, or pickled pandas DataFrame input

Missing model genes are filled with the training mean, which is neutral after
standardization for logistic regression. If fewer than half of the 2,000 model
genes match, the CLI warns that gene identifiers or expression scale are likely
wrong.

Before scoring a new matrix, run:

```bash
python run_tumor_normal_workflow.py expr.csv --labels labels.csv
python inspect_expression_input.py expr.csv -o expr.qc.json
```

The workflow writes a complete result folder with `qc.json`, `scores.csv`,
optional calibration files, `explanations.csv`, `manifest.json`, and
`workflow_report.md`.

See `example_workflow_output/` for a small checked-in reference output from the
bundled example input.

The QC report checks model-gene match rate, non-finite values, value range,
standardized distribution shift against the model scaler, and score summary.
`model_qc_reference.json` stores the heuristic rules and compact reference
summaries. A PASS result does not prove that a new pipeline is validated, but a
WARN/FAIL result is a strong reason to stop and inspect normalization or
calibrate/refit before making hard calls.

## Threshold calibration

If you have labeled samples from a new tissue or pipeline, score them and choose
a cutoff:

```bash
python score_tumor_normal.py expr.csv -o calls.csv
python calibrate_threshold.py calls.csv labels.csv -o calibration_thresholds.csv
```

`labels.csv` should contain `sample,label`, where label is `tumor` or `normal`
(1/0 is also accepted). The recommended threshold is the Youden's-J cutoff. Pass
that value back into scoring with `--threshold`.

## Model explanations

For per-sample model debugging:

```bash
python explain_scores.py expr.csv -o explanations.csv --top-n 10
```

The output lists the largest positive and negative contributions to the LR
logit. Positive contributions push the probability toward tumor; negative
contributions push it toward normal. The companion `model_gene_metadata.csv`
contains all model genes, LR coefficients, and scaling parameters.

## Validation status

Internal validation is complete:

- patient-held-out test AUC: 0.997
- grouped 5-fold CV AUC: 0.997 +/- 0.003
- leave-one-cancer-type-out macro-mean AUC: 0.994 (pooled 0.988)

External non-TCGA smoke validation is also available:

- CPTAC-3/GDC STAR-Counts, 100 primary tumor + 100 solid tissue normal files
- AUC: 0.989
- Accuracy at threshold 0.5: 0.955
- Files: `external-validation/cptac_gdc/`

Cross-platform checks are also available:

- GTEx/Toil normal tissues, 540 healthy normal samples across 27 primary sites:
  538/540 called tumor at the default 0.5 threshold.
- TCGA/Toil pipeline sanity check, 100 primary tumor + 100 solid tissue normal:
  AUC 0.992 but default-threshold accuracy only 0.515; Toil-specific
  re-thresholding at 0.999975 recovered accuracy to 0.970.
- Files: `external-validation/gtex_xena/` and
  `external-validation/tcga_toil_xena/`

The model should currently be described as a **GDC STAR-Counts-style
tumor-vs-adjacent-normal classifier**, not a clinical diagnostic model and not
a plug-in classifier for arbitrary TPM/RSEM/GTEx/GEO matrices.
