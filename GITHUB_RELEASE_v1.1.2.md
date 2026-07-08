# v1.1.2-gdc-starcounts

Publication metadata refresh for the TCGA/GDC STAR-Counts tumor-vs-normal
RNA-seq classifier. The fitted model weights, scoring behavior, and headline
validation metrics are unchanged from v1.1.1.

## Highlights

- Added `.zenodo.json` and `codemeta.json` for machine-readable software
  metadata.
- Added README badges for CI, release, license, and citation metadata.
- Included Zenodo/CodeMeta metadata in the lightweight release bundle.
- Added hosted CI coverage for `audit_publication_readiness.py`.
- Clarified that legacy pickle/RF training artifacts are intentionally excluded
  from the public Git history and are not needed for the lightweight release.

## Validation

```text
python -m pytest -q -rs
25 passed, 2 skipped

python build_release_lite.py --smoke --timeout-seconds 300
PASS

python run_release_acceptance.py --timeout-seconds 300
PASS
```

## Release Asset

```text
tcga-tumor-normal-release-lite.zip
SHA256: 85e6d228c03f9f14f34c54031a30d0e6e9e9d6c20c6902b0cb217724c7ff3193
Bytes: 297655
Entries: 72
```

## Intended Use

Use this release for research scoring of GDC STAR-Counts-style `log2(TPM+1)`
bulk RNA-seq matrices with rows as samples and Ensembl genes as columns.

Do not use this model for clinical diagnosis or patient management. Do not use
it for direct hard calls on Toil/RSEM, GTEx, GEO, raw counts, FPKM, microarray,
single-cell, or spatial data without pipeline-specific validation,
recalibration, or refitting.
