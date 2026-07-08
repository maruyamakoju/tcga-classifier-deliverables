# v1.1.3-gdc-starcounts

Quality hardening release for the TCGA/GDC STAR-Counts tumor-vs-normal
RNA-seq classifier. The fitted model weights and headline validation metrics
are unchanged from v1.1.2.

## Highlights

- Public scoring CLI now exposes only the supported pure NumPy logistic
  regression path from `deployable_lr_weights.npz`.
- Unsupported legacy pickle/RF scorer options now fail clearly instead of
  reaching missing artifacts in the lightweight bundle.
- Duplicate gene columns and Ensembl-version collisions are rejected before
  scoring/QC.
- Cohort-adaptation label handling now correctly supports numeric-string
  labels, missing label matches, and extra label rows.
- Threshold calibration now requires all scored samples to have labels by
  default; use `--min-match-fraction` only for intentional subsets.
- Added hosted GitHub repository audit coverage for branch protection, release
  asset digest, repository topics, and stale pip Dependabot PRs.

## Validation

```text
python -m pytest -q -rs
31 passed, 2 skipped

python build_release_lite.py --smoke --timeout-seconds 300
PASS

python run_release_acceptance.py --timeout-seconds 300
PASS
```

## Release Asset

```text
tcga-tumor-normal-release-lite.zip
SHA256: 8d809178bc8931cf486d911b55ebcfcdb8f4a82283a94a330d02f03001549e06
Bytes: 299259
Entries: 72
```

## Intended Use

Use this release for research scoring of GDC STAR-Counts-style `log2(TPM+1)`
bulk RNA-seq matrices with rows as samples and Ensembl genes as columns.

Do not use this model for clinical diagnosis or patient management. Do not use
it for direct hard calls on Toil/RSEM, GTEx, GEO, raw counts, FPKM, microarray,
single-cell, or spatial data without pipeline-specific validation,
recalibration, or refitting.
