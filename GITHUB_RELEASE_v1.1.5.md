# v1.1.5-gdc-starcounts

Explanation guardrail release for the TCGA/GDC STAR-Counts tumor-vs-normal
RNA-seq classifier. The fitted model weights and headline validation metrics
are unchanged from v1.1.4.

## Highlights

- `explain_scores.py` now rejects invalid matched model-gene values before
  writing explanations, matching the scorer and workflow behavior.
- New `--max-invalid-cell-fraction` and `--allow-invalid-values` flags are
  available on explanation generation.
- Safety tests now cover invalid matched-value explanation refusal and explicit
  allow behavior.
- Hosted repository audit coverage now includes administrator-enforced branch
  protection, linear history, conversation resolution, protected `v*` release
  tag rulesets, and Dependabot vulnerability alerts.

## Validation

```text
python -m pytest -q -rs
35 passed, 2 skipped

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
SHA256: dfe952a183a318d3a9ee3e58612982dffaf7a68f39da727325bcfc3a5a77422b
Bytes: 304175
Entries: 72
```

## Intended Use

Use this release for research scoring of GDC STAR-Counts-style `log2(TPM+1)`
bulk RNA-seq matrices with rows as samples and Ensembl genes as columns.

Do not use this model for clinical diagnosis or patient management. Do not use
it for direct hard calls on Toil/RSEM, GTEx, GEO, raw counts, FPKM, microarray,
single-cell, or spatial data without pipeline-specific validation,
recalibration, or refitting.
