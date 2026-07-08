# v1.1.4-gdc-starcounts

Invalid-input guardrail release for the TCGA/GDC STAR-Counts tumor-vs-normal
RNA-seq classifier. The fitted model weights and headline validation metrics
are unchanged from v1.1.3.

## Highlights

- Matched model-gene cells that are missing, non-numeric, `NaN`, or infinite
  are now counted before training-mean imputation.
- `score_tumor_normal.py` and `run_tumor_normal_workflow.py` now stop before
  writing scores when invalid matched values are present unless explicitly
  allowed.
- New `--max-invalid-cell-fraction` and `--allow-invalid-values` flags make any
  reviewed tolerance/override explicit.
- Safety tests now cover invalid matched-value scorer refusal, explicit allow
  behavior, and workflow stop behavior.

## Validation

```text
python -m pytest -q -rs
34 passed, 2 skipped

python build_release_lite.py --smoke --timeout-seconds 300
PASS

python audit_publication_readiness.py
PASS
python audit_cli_entrypoints.py
PASS
python audit_release_docs.py
PASS
python validate_output_contracts.py
PASS
```

## Release Asset

```text
tcga-tumor-normal-release-lite.zip
SHA256: d8b2e5dcef61802291f4a0acb3e251183908723388e0fcd268ea4a8a6c8c5169
Bytes: 303149
Entries: 72
```

## Intended Use

Use this release for research scoring of GDC STAR-Counts-style `log2(TPM+1)`
bulk RNA-seq matrices with rows as samples and Ensembl genes as columns.

Do not use this model for clinical diagnosis or patient management. Do not use
it for direct hard calls on Toil/RSEM, GTEx, GEO, raw counts, FPKM, microarray,
single-cell, or spatial data without pipeline-specific validation,
recalibration, or refitting.
