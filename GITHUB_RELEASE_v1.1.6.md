# v1.1.6-gdc-starcounts

Remaining CLI guardrail release for the TCGA/GDC STAR-Counts tumor-vs-normal
RNA-seq classifier. The fitted model weights and headline validation metrics
are unchanged from v1.1.5.

## Highlights

- `cohort_adapt_score.py` now rejects invalid matched model-gene values before
  writing adapted scores unless explicitly allowed.
- `cancer-type-classifier/predict_cancer_type.py` now rejects invalid matched
  model-gene values before writing predictions unless explicitly allowed.
- Adaptation and cancer-type outputs now report matched/missing model-gene
  counts and invalid matched-value counts.
- Safety/unit tests now cover the remaining invalid matched-value CLI paths.

## Validation

```text
python -m pytest -q -rs
36 passed, 2 skipped

python build_release_lite.py --smoke --timeout-seconds 300
PASS

python audit_release_docs.py
PASS
python validate_output_contracts.py
PASS
python audit_cli_entrypoints.py
PASS
```

## Release Asset

```text
tcga-tumor-normal-release-lite.zip
SHA256: feb800b657d199ef6ab48d0ed376013a301c405f9200393054d2881b48d9b685
Bytes: 305259
Entries: 72
```

## Intended Use

Use this release for research scoring of GDC STAR-Counts-style `log2(TPM+1)`
bulk RNA-seq matrices with rows as samples and Ensembl genes as columns.

Do not use this model for clinical diagnosis or patient management. Do not use
it for direct hard calls on Toil/RSEM, GTEx, GEO, raw counts, FPKM, microarray,
single-cell, or spatial data without pipeline-specific validation,
recalibration, or refitting.
