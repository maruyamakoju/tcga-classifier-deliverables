# v2.0.1-gdc-starcounts

Test-hygiene patch on top of `v2.0.0-gdc-starcounts`. No shipped code, deployed
weights, bundled artifact payload, or headline metric changed from v2.0.0; this
release only makes continuous integration fully green on a clean checkout.

## What changed

- Fixed `tests/test_reproducibility.py::test_train_test_split_is_patient_disjoint`
  so it skips (via the `features_npy` fixture) when the gitignored full-data
  feature matrix is unavailable, instead of hard-failing with `FileNotFoundError`
  on CI. The other full-data reproduction tests already skip this way; this test
  was missing the same guard, so it passed locally (artifact present) but failed
  the non-required `full-unit` CI job.

The lightweight release bundle is regenerated so its committed VERSION and docs
read `v2.0.1-gdc-starcounts`; the scoring payload is otherwise identical to
v2.0.0.

## Validation

- `python -m pytest -q`: 413 passed (0 failed on a checkout without the local
  feature matrix, where this test now skips).
- `python -m ruff check .`: passed.
- Deterministic build, release acceptance, and cross-platform CI: green.

## Security boundary for the ZIP

Do not extract or execute the downloaded asset until its SHA-256 matches the
trusted value published below. Acceptance requires the digest explicitly:

```bash
python validate_zip_bundle.py tcga-tumor-normal-release-lite.zip --expected-sha256 c66f207bf7572c443c02f8aaa658ea70344d656054179b1576ebb7ff6fc39040
```

Without a trusted digest, `--skip-acceptance` is limited to non-executing
structural inspection.

## Release asset

- Asset: `tcga-tumor-normal-release-lite.zip`
- SHA256: `c66f207bf7572c443c02f8aaa658ea70344d656054179b1576ebb7ff6fc39040`
- Size: `937564` bytes
- Zip entries: `75`

These asset values are taken from the generated `RELEASE_ARTIFACTS.json` after the
deterministic builder and validators passed.
