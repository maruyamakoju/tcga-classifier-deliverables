# v2.0.2-gdc-starcounts

Maintenance patch on top of `v2.0.1-gdc-starcounts`. No shipped scoring code,
deployed weights, bundled payload logic, or headline metric changed from
v2.0.0/v2.0.1; the release-lite payload differs only in the committed VERSION and
doc strings.

## What changed

- Configured dependabot to ignore the version-locked scientific reproduction
  stack (`numpy`, `pandas`, `scipy`, `scikit-learn`). These are pinned exactly in
  `requirements-training.txt` so the shipped weights and golden metrics stay
  reproducible, and are bumped manually with a golden re-verification. This stops
  the recurring automated PRs that fight the pins (e.g. the pandas 3.x bump, whose
  CI is red because it breaks the suite). Development tooling, `pyarrow`,
  `requests`, and GitHub Actions still receive automated updates.
- Raised development-tooling lower bounds: `pytest>=9.1.1`, `ruff>=0.15.21`,
  `requests>=2.34.2`. CI already resolved to these ranges, so this is a
  lower-bound documentation change with no behavioral effect.

## Validation

- `python -m pytest -q`: 413 passed.
- `python -m ruff check .`: passed.
- Deterministic build, release acceptance, and cross-platform CI: green.

## Security boundary for the ZIP

Do not extract or execute the downloaded asset until its SHA-256 matches the
trusted value published below. Acceptance requires the digest explicitly:

```bash
python validate_zip_bundle.py tcga-tumor-normal-release-lite.zip --expected-sha256 a521a49c0c38eb704e69a3505a1c285d70d72d185f6f4663be0c2336dce514fe
```

Without a trusted digest, `--skip-acceptance` is limited to non-executing
structural inspection.

## Release asset

- Asset: `tcga-tumor-normal-release-lite.zip`
- SHA256: `a521a49c0c38eb704e69a3505a1c285d70d72d185f6f4663be0c2336dce514fe`
- Size: `938822` bytes
- Zip entries: `75`

These asset values are taken from the generated `RELEASE_ARTIFACTS.json` after the
deterministic builder and validators passed.
