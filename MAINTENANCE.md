# Maintenance guide

This repository is public and research-only. Keep operational changes small,
auditable, and backed by the release validators.

## Routine checks

Run before merging changes to `main`:

```bash
python -m pytest -q -rs
python audit_publication_readiness.py
python run_release_acceptance.py --timeout-seconds 300
```

Run before publishing a release ZIP:

```bash
python build_release_lite.py --smoke --timeout-seconds 300
python validate_release_lite.py --release-dir release-lite --zip tcga-tumor-normal-release-lite.zip --source-root . --artifacts RELEASE_ARTIFACTS.json
python validate_zip_bundle.py tcga-tumor-normal-release-lite.zip
```

## Release rules

- Bump `VERSION`, `RELEASE_METADATA.json`, `CITATION.cff`, `.zenodo.json`,
  `codemeta.json`, and versioned documentation together.
- Rebuild `release-lite/` and `tcga-tumor-normal-release-lite.zip`; never edit
  generated manifest/checksum files by hand.
- Update or create `GITHUB_RELEASE_<version>.md` with the final ZIP SHA256,
  byte size, and entry count from `RELEASE_ARTIFACTS.json`.
- Keep previous tags immutable after public release unless correcting a broken
  upload before users have consumed it.
- Do not add full training matrices, private manifests, credentials, PHI,
  private sample identifiers, or unpublished raw data to Git.

## Dependency updates

Dependabot opens dependency PRs for GitHub Actions only. Python dependency
updates are manual because `requirements-light.txt` and `requirements.txt` are
part of the release payload; changing them requires rebuilding `release-lite/`,
`SHA256SUMS.txt`, `release_manifest.json`, and the release ZIP in the same
commit. Keep the lightweight runtime dependency contract minimal: default
scoring should require only NumPy, pandas, and pyarrow.

## Scope guardrails

The release remains validated for GDC STAR-Counts-style `log2(TPM+1)` bulk
RNA-seq tumor-vs-adjacent-normal contrasts. Do not broaden claims to Toil/RSEM,
GTEx, GEO, raw counts, FPKM, microarray, single-cell, spatial data, cancer
screening, diagnosis, or patient management without separate validation.
