# v1.1.1-gdc-starcounts

Release-quality hardening update for the TCGA/GDC STAR-Counts tumor-vs-normal
RNA-seq classifier. The fitted model weights and headline validation metrics are
unchanged from v1.1.0; this release improves input validation, release
reproducibility, regression coverage, and GitHub-ready project metadata.

## Highlights

- Hardened scoring against non-finite expression values.
- Rejected multiclass `.npz` weights in tumor-vs-normal scoring paths with clear
  errors.
- Validated calibration inputs for duplicate samples, non-finite probabilities,
  and out-of-range thresholds.
- Aligned Youden-J threshold tie-breaking with the shared metrics core.
- Added source-parity validation so stale `release-lite/` bundles fail
  acceptance.
- Validated `RELEASE_ARTIFACTS.json` against the zip and release directory.
- Made release zip generation byte-reproducible.
- Added CI, citation metadata, license, notice, contributing guide, and security
  policy.

## Validation

```text
python -m pytest
25 passed, 2 skipped

python run_release_acceptance.py --timeout-seconds 300
PASS

python validate_release_lite.py --release-dir release-lite --zip tcga-tumor-normal-release-lite.zip --source-root . --artifacts RELEASE_ARTIFACTS.json
PASS
```

## Release Asset

```text
tcga-tumor-normal-release-lite.zip
SHA256: 96c7e5e24a08ca5852b44201d43489798ea2badb474c43757dc975d4ea8d3ad2
Bytes: 295293
Entries: 70
```

## Intended Use

Use this release for research scoring of GDC STAR-Counts-style `log2(TPM+1)`
bulk RNA-seq matrices with rows as samples and Ensembl genes as columns.

Do not use this model for clinical diagnosis or patient management. Do not use
it for direct hard calls on Toil/RSEM, GTEx, GEO, raw counts, FPKM, microarray,
single-cell, or spatial data without pipeline-specific validation,
recalibration, or refitting.
