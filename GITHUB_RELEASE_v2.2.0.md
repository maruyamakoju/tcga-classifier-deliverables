# v2.2.0-gdc-starcounts

External-validation robustness release. Two real bugs fixed in the
(development-only) external-validation tooling, and a post-fix **live**
re-validation performed. The shipped scoring code, deployed weights, bundled
payload logic, and headline metrics are unchanged.

## Highlight: live CPTAC-3 re-validation reproduces the headline metric

A fresh fetch of the locked 200-sample CPTAC-3 cohort from **GDC Data Release 45.0
(2025-12-04)**, scored by the current code, reproduced the committed result:

| Metric | Live (GDC 45.0) | Committed snapshot |
| --- | --- | --- |
| AUC | **0.9886** | 0.9886 |
| Average precision | 0.98505 | 0.98506 |
| Accuracy @ 0.5 | 0.960 | 0.955 |
| n (tumor / normal) | 200 (100 / 100) | 200 (100 / 100) |
| Model genes matched | 2000 / 2000 | — |

The AUC reproduces exactly; accuracy differs by a single threshold-boundary
sample. The cross-platform TCGA-Toil and GTEx summaries remain historical
snapshots (not re-run live in this release).

## Fixes that made the live run possible

- **CPTAC/GDC manifest validator** no longer rejects a file mapped to several
  biospecimen submitter IDs of one case and one `sample_type`: the tumor/normal
  label is unambiguous and already guaranteed by the `n_sample_types==1` and
  multi-*case* checks. Only a genuine multi-*case* mapping (different patients) is
  rejected. This unblocks `--refresh-manifest` on current GDC releases without
  ever guessing a label.
- **Windows long-path robustness** in `provenance.contained_cache_path`: it now
  verifies the hashed cache key is a bare basename instead of re-resolving the
  child path, which `Path.resolve()` could prefix with a `\\?\` extended-length
  marker on only one side past `MAX_PATH`, spuriously rejecting a legitimate long
  cache directory.

## Also

- Finalized the publication checklist; added an **Archival & DOI** section
  tracking the one outstanding step (mint a Zenodo DOI and backfill it into
  `CITATION.cff`, `.zenodo.json`, and a `README` badge).

## Validation

- `python -m pytest -q`: 416 passed (+3: long-path cache containment, CPTAC
  same-case multi-biospecimen accepted, multi-case still rejected).
- `python -m ruff check .` and `python -m mypy`: passed.
- Deterministic build, release acceptance, and cross-platform CI: green.

## Security boundary for the ZIP

Do not extract or execute the downloaded asset until its SHA-256 matches the
trusted value published below. Acceptance requires the digest explicitly:

```bash
python validate_zip_bundle.py tcga-tumor-normal-release-lite.zip --expected-sha256 b4332c69acd081ede3a98d74fb060d225406c13e1e969ef1fb533e1d2f3955d0
```

Without a trusted digest, `--skip-acceptance` is limited to non-executing
structural inspection.

## Release asset

- Asset: `tcga-tumor-normal-release-lite.zip`
- SHA256: `b4332c69acd081ede3a98d74fb060d225406c13e1e969ef1fb533e1d2f3955d0`
- Size: `945560` bytes
- Zip entries: `76`

These asset values are taken from the generated `RELEASE_ARTIFACTS.json` after the
deterministic builder and validators passed.
