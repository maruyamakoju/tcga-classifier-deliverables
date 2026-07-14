# v2.4.0-gdc-starcounts

Live external re-validation of **all three cohorts**, with archived provenance. No
shipped scoring code, deployed weights, bundled payload logic, or headline metric
changed.

## Highlight: all three external cohorts reproduce live

A post-fix live re-fetch of each locked cohort from its current provider source,
scored by the current code, reproduced every committed metric:

| Cohort | Source | Metric | Live | Committed |
| --- | --- | --- | --- | --- |
| CPTAC-3 | GDC Data Release 45.0 | AUC | **0.9886** | 0.9886 |
| TCGA-Toil | UCSC Xena | AUC / acc | **0.9923 / 0.515** | 0.9923 / 0.515 |
| GTEx normals | UCSC Xena | FPR@0.5 (n=540) | **0.9963** | 0.9963 |

Each run's full provenance (source revision, cache fingerprints, per-file content
hashes, code and environment hashes) plus the reproduced summary and threshold
sweep are archived under `external-validation/<cohort>/revalidation/`, and
`external-validation/REVALIDATION.md` is the canonical record. This closes the
prior "committed external metrics are a historical snapshot / no post-fix live
rerun" gap for all three cohorts (the committed CSVs remain the frozen snapshot;
the live runs confirm the current code reproduces them).

## Fix that unblocked the Xena cohorts

- **Xena decode tolerance** (`validate_gtex_xena.py`): a truly zero-expression
  gene stored as rounded `log2(TPM+0.001)` round-trips to a tiny negative TPM
  (~1e-8). The decode's `1e-9` rejection tolerance was too tight and rejected
  legitimate live Xena data; widened to `1e-6` (the value is clamped to zero).
  This is the third over-strict fail-closed check the live runs surfaced (after
  the CPTAC multi-biospecimen validator and the Windows long-path cache check in
  v2.2.0). External-validation tooling only; the shipped GDC scoring path is
  unaffected.

## Also

- Updated the external-validation provenance disclosure across the docs to point
  at `external-validation/REVALIDATION.md`.

## Validation

- `python -m pytest -q`: 417 passed.
- `python -m ruff check .` and `python -m mypy`: passed.
- Deterministic build, release acceptance, and cross-platform CI: green.

## Security boundary for the ZIP

Do not extract or execute the downloaded asset until its SHA-256 matches the
trusted value published below. Acceptance requires the digest explicitly:

```bash
python validate_zip_bundle.py tcga-tumor-normal-release-lite.zip --expected-sha256 2355392b19e6ce9a2f08697f247c905edeeb06f280110ff83a9a053fa8febd9e
```

Without a trusted digest, `--skip-acceptance` is limited to non-executing
structural inspection.

## Release asset

- Asset: `tcga-tumor-normal-release-lite.zip`
- SHA256: `2355392b19e6ce9a2f08697f247c905edeeb06f280110ff83a9a053fa8febd9e`
- Size: `949729` bytes
- Zip entries: `76`

These asset values are taken from the generated `RELEASE_ARTIFACTS.json` after the
deterministic builder and validators passed.
