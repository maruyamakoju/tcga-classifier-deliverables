# v2.3.0-gdc-starcounts

Documentation release. The external-validation provenance disclosure across the
project docs now reflects the **v2.2.0 live CPTAC-3 re-validation**. No shipped
scoring code, deployed weights, bundled payload logic, or headline metric changed.

## What changed

The "committed external metrics are a historical snapshot / no post-fix
live-network rerun" disclosure was updated in `README.md`, `MODEL_CARD.md`,
`EXECUTIVE_SUMMARY.md`, `REPORT.md`, `REPRODUCIBILITY.md`, `INDEX.md`,
`MAINTENANCE.md`, `RELEASE_BUNDLE.md`,
`cross-platform-adaptation/CROSS_PLATFORM_ADAPTATION.md`, `.zenodo.json`, and
`model_qc_reference.json` to state accurately that:

- In v2.2.0 a post-fix live re-fetch of the locked CPTAC-3 cohort from **GDC Data
  Release 45.0** (2025-12-04), scored by the current code, **reproduced its
  committed AUC (0.9886)** (accuracy 0.96 vs 0.955, one threshold-boundary
  sample).
- The **TCGA-Toil and GTEx** cross-platform summaries have **not** been re-run
  live and **remain historical snapshots**.

Historical changelog entries (which describe earlier versions' state) are left
unchanged. The publication checklist's scientific-framing item was updated to
match.

## Validation

- `python -m pytest -q`: 416 passed.
- `python -m ruff check .` and `python -m mypy`: passed.
- Deterministic build, release acceptance, and cross-platform CI: green.

## Security boundary for the ZIP

Do not extract or execute the downloaded asset until its SHA-256 matches the
trusted value published below. Acceptance requires the digest explicitly:

```bash
python validate_zip_bundle.py tcga-tumor-normal-release-lite.zip --expected-sha256 011184ceb16285758e672cf12470fc1bc2f929b3a97d3d150cb3ef4b1e261bd9
```

Without a trusted digest, `--skip-acceptance` is limited to non-executing
structural inspection.

## Release asset

- Asset: `tcga-tumor-normal-release-lite.zip`
- SHA256: `011184ceb16285758e672cf12470fc1bc2f929b3a97d3d150cb3ef4b1e261bd9`
- Size: `947519` bytes
- Zip entries: `76`

These asset values are taken from the generated `RELEASE_ARTIFACTS.json` after the
deterministic builder and validators passed.
