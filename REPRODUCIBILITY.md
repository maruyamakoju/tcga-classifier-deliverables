# Reproducibility notes

Release: `v2.0.0-gdc-starcounts` (`2026-07-12`; public scoring-library API `3.0.0`)

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
python validate_zip_bundle.py tcga-tumor-normal-release-lite.zip --expected-sha256 <trusted-published-sha256>
```

Dependency profiles are intentionally separated by role:

- `requirements-light.txt`: public scoring/runtime bundle.
- `requirements-external-validation.txt`: network validation maintenance.
- `requirements-dev.txt`: lint, tests, audits, and release tooling.
- `requirements-training.txt`: canonical exact scientific reproduction stack.
- `requirements.txt`: complete maintenance/retraining profile composed from
  the external-validation and canonical training profiles.

For retraining or full-artifact maintenance use:

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
- maximum probability delta vs the rounded `example_output.csv` is at most
  `1e-6` (the scorer writes full-precision probabilities)
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
`X_full_float64.npy` with sibling `X_genes.npy` and `X_samples.npy`.

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

The builder stages into a temporary directory and publishes only after
validation. ZIP members are sorted, use fixed release-date timestamps,
`create_system=3`, and mode `0644`, so identical source produces identical
bytes across operating systems. Check for drift without modifying committed
artifacts:

```bash
python build_release_lite.py --check --smoke --timeout-seconds 300
```

CI installs separated role-specific dependency profiles plus the exact-pinned
canonical training profile, runs `pip check`, Ruff, and the full pytest suite,
then runs deterministic release acceptance on
Windows, Linux, and macOS. Python 3.11 is the scientific-reproduction job;
Python 3.13 is a lightweight scoring-acceptance job only. GitHub Actions are
pinned by commit.

To verify the zip as a standalone artifact, extract it into a temporary clean
directory and run acceptance there:

```bash
python validate_zip_bundle.py tcga-tumor-normal-release-lite.zip --expected-sha256 <trusted-published-sha256>
```

To validate an already-built bundle without rebuilding:

```bash
python validate_release_lite.py --release-dir release-lite --zip tcga-tumor-normal-release-lite.zip
python validate_zip_bundle.py tcga-tumor-normal-release-lite.zip --expected-sha256 <trusted-published-sha256>
```

The default self-test uses `deployable_lr_weights.npz` and should not emit a
scikit-learn pickle warning. The canonical stack that exactly reproduced the
shipped model is pinned by `requirements-training.txt`:

- CPython 3.11.x (the reference run used 3.11.9)
- pandas 2.3.3
- scikit-learn 1.8.0
- numpy 1.26.4
- scipy 1.15.3

Python 3.13 with scikit-learn 1.9 produced coefficient/weight and out-of-fold
drift during refitting, so it is supported only for lightweight scoring
acceptance, not exact scientific reproduction. The lightweight release does
not load scikit-learn pickle artifacts during scoring. For regenerating
`deployable_lr_weights.npz` or verifying training/LOCO artifacts, use the
canonical Python 3.11 training stack.

Do not loosen `tests/golden_numbers.json`, shipped-weight parity, LOCO, or other
committed numerical tolerances to make a noncanonical environment pass. A drift
is a failed reproduction or environment mismatch that must be investigated and
reported.

## Exact full-data reproduction

The historical pandas feature object is trusted local input. Pickle loading can
execute code: use this step only after independently verifying the source. The
one-time converter has a separate exact environment because model fitting must
remain on the canonical Python 3.11 stack:

```bash
python3.13 -m venv .venv-feature-export
source .venv-feature-export/bin/activate
# Windows PowerShell instead: .\.venv-feature-export\Scripts\Activate.ps1
python -m pip install -r training_tools/requirements-feature-export.txt
python cancer-type-classifier/export_features_npy.py \
  --source X_full_filtered.pkl \
  --output-dir cancer-type-classifier \
  --dtype both \
  --trusted-source-pickle
```

This pickle-to-NPY conversion is a separate trusted-artifact migration step;
do not use its pandas 3 environment for model refitting. After export, switch
to the canonical Python 3.11 `requirements-training.txt` environment for the
commands below. `training_tools/feature_export_lock.json` commits the canonical
source and output SHA-256 values, sizes, shapes, dtypes, and converter versions
without committing the large pickle or arrays. The generated schema-v3
`cancer-type-classifier/X_full.export.json` contains only a source basename, not a machine-specific
absolute path. A mismatch fails before publication; the explicit
`--allow-noncanonical-export` escape is for development only, and canonical
trainers reject such a manifest unless `--allow-unverified-features` is also
deliberately supplied. The NPY inputs avoid pandas-pickle runtime coupling.

Use the float32 export only for the historical 17-class cancer-type path:

```bash
python cancer-type-classifier/train_cancer_type_classifier.py \
  --features cancer-type-classifier/X_full.npy \
  --metadata selected_files.csv \
  --output-dir <fresh-output-dir> \
  --gene-symbols cancer-type-classifier/gene_id_to_name.csv \
  --verify-shipped cancer-type-classifier/cancer_type_lr_weights.npz
```

Use the float64 export for exact binary and LOCO reproduction:

```bash
python train_classifier.py \
  --features cancer-type-classifier/X_full_float64.npy \
  --metadata selected_files.csv \
  --train-index train_idx.npy \
  --test-index test_idx.npy \
  --verify-shipped deployable_lr_weights.npz \
  --output-dir <fresh-output-dir>

python cross-cancer-holdout/run_loco.py \
  --features cancer-type-classifier/X_full_float64.npy \
  --metadata selected_files.csv \
  --verify-existing cross-cancer-holdout \
  --tolerance 1e-10
```

Using float32 for binary/LOCO changes the numerical contract; using float64
for the historical cancer-type pipeline changes its recorded reproduction
path.

## Input contract

`score_tumor_normal.py` expects:

- rows = samples, columns = genes
- Ensembl gene IDs, with or without version suffixes
- values = log2(TPM+1) on the GDC STAR-Counts scale
- `.csv`, `.tsv`, `.txt`, or `.parquet` input; pickled expression matrices are
  intentionally rejected by the public CLIs

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

The reported threshold metrics are apparent/resubstitution metrics: the same
labeled samples select and evaluate the cutoff. They are not an independent
generalization estimate. A warning is emitted when either class contains fewer
than 10 samples.

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
- grouped 5-fold CV AUC: 0.996 +/- 0.003
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

These committed metrics are a historical snapshot. Version 2.0.0 fixed cache
identity, empty/colliding data checks, atomic cache writes, locked-cohort
handling, content hashes, and generated run-manifest provenance, but no post-fix
live-network rerun was performed. Therefore the historical numbers above must
not be represented as newly reproduced by the fixed scripts.

The validators distinguish two modes:

- Live mode requires a concrete, non-`unversioned` `--source-revision`. Locked
  study/sample-type/project/label fields are validated before access. New GDC
  downloads additionally require a valid provider MD5.
- `--offline` (alias `--cache-only`) never accesses the network and fails if a
  semantic fingerprint, content hash, or cache axis does not match. It may be
  used only with an already-valid cache carrying the requested source revision.

The committed CPTAC `external-validation/cptac_gdc/sampled_manifest.csv`
predates provider-MD5 capture and
does not contain `md5sum`. This is disclosed historical metadata, not a value to
infer. It can identify an offline historical cohort when a matching private
cache exists, but cannot authorize a new download. For a future live run,
`--refresh-manifest` re-queries metadata for the exact locked file IDs, refuses
cohort identity changes, and supplies the provider MD5 before expression fetch.

Derived CSV/Markdown outputs are staged and published as one rollback-protected
set, with the run manifest published last. Individually reusable caches are
still updated atomically and may survive a later analysis failure. Safe hashed
cache keys replace raw provider IDs, so legacy raw-ID and pickle caches are not
read as current cache generations.

For a future, separately reported live rerun, use the committed locked cohorts,
record an actual provider revision, refresh the caches, and write into fresh
directories:

```bash
python external-validation/validate_gtex_xena.py --sample-manifest external-validation/gtex_xena/sampled_gtex_manifest.csv --source-revision <provider-revision> --refresh --out-dir <fresh-dir>
python external-validation/validate_tcga_toil_xena.py --sample-manifest external-validation/tcga_toil_xena/sampled_tcga_toil_manifest.csv --source-revision <provider-revision> --refresh --out-dir <fresh-dir>
python external-validation/validate_cptac_gdc.py --sample-manifest external-validation/cptac_gdc/sampled_manifest.csv --source-revision <gdc-release> --refresh-manifest --refresh-expression-cache --out-dir <fresh-dir>
```

Archive each resulting run manifest with the outputs. Do not overwrite
the committed snapshot until the live run, hashes, cohort identity, and metric
diff have been reviewed.

For a cache-only rerun in a private/full checkout, point `--out-dir` to the
directory that already contains the matching Parquet cache and add `--offline`.
Do not use a downloaded cache unless its provenance and content hashes have
been independently trusted.

The model should currently be described as a **GDC STAR-Counts-style
tumor-vs-adjacent-normal classifier**, not a clinical diagnostic model and not
a plug-in classifier for arbitrary TPM/RSEM/GTEx/GEO matrices.

In all scoring outputs, `tumor_probability` is a logistic model score, not
clinical risk or a calibrated diagnostic probability. Cohort adaptation is
disabled by default; adapted modes are explicit transductive,
composition-dependent opt-ins requiring at least 20 samples by default, and
scores from separately adapted batches are not comparable. LOCO does not
remove project/procurement/center/batch confounding, and literature consistency
does not prove a causal mechanism.
