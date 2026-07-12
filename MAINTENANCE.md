# Maintenance guide

Current release: `v2.0.0-gdc-starcounts` (`2026-07-12`); public scoring-library
API `3.0.0`.

This repository is public and research-only. Keep operational changes small,
auditable, and backed by the release validators.

## Routine checks

Run before merging changes to `main`:

```bash
python -m pytest -q -rs
python -m ruff check .
python -m pip check
python audit_publication_readiness.py
python run_release_acceptance.py --timeout-seconds 300
```

Run before publishing a release ZIP:

```bash
python build_release_lite.py --smoke --timeout-seconds 300
python validate_release_lite.py --release-dir release-lite --zip tcga-tumor-normal-release-lite.zip --source-root . --artifacts RELEASE_ARTIFACTS.json
python validate_zip_bundle.py tcga-tumor-normal-release-lite.zip --expected-sha256 <trusted-published-sha256>
```

Run after publishing or changing hosted repository settings:

```bash
python audit_github_repository.py
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
- Never rewrite a historical GitHub release body with later work. Add new
  findings only to the next release note.
- Do not add full training matrices, private manifests, credentials, PHI,
  private sample identifiers, or unpublished raw data to Git.

## Dependency updates

Dependabot covers GitHub Actions and pip. Keep role-specific dependency
profiles synchronized: `requirements-light.txt` for public scoring,
`requirements-external-validation.txt` for network validation,
`requirements-dev.txt` for tests/lint/release tooling, and `requirements.txt`
for the complete development environment. Exact shipped-model reproduction is
separately pinned in `requirements-training.txt` to Python 3.11, NumPy 1.26.4,
pandas 2.3.3, SciPy 1.15.3, and scikit-learn 1.8.0. Python 3.13 is a lightweight
scoring-acceptance target only. Any payload dependency change
requires rebuilding the release in the same commit. Keep the lightweight
runtime contract minimal: default scoring should require only NumPy, pandas,
and pyarrow.

Never broaden golden, shipped-weight, OOF, or LOCO tolerances to conceal drift
from a noncanonical training stack. Investigate and report the environment
delta; update a golden value only through a separately reviewed scientific
change with regenerated provenance.

## Deterministic release maintenance

The builder must stage, validate, and atomically publish canonical release
outputs. Run `python build_release_lite.py --check --smoke --timeout-seconds
300` to detect drift without mutating committed artifacts. Never edit
`release_manifest.json`, `SHA256SUMS.txt`, the ZIP, or
`RELEASE_ARTIFACTS.json` by hand.

Treat a downloaded ZIP as untrusted until its SHA-256 matches a trusted
published digest. `validate_zip_bundle.py` requires `--expected-sha256` before
it extracts or executes archive content. Without a trusted digest, use only
`--skip-acceptance` for non-executing structural inspection.

## External-validation maintenance

Use the committed sample manifests as locked cohorts, provide a real
`--source-revision`, refresh caches explicitly, write to a fresh output
directory, and archive the generated run manifest. Review semantic cache
fingerprints, content hashes, cohort identity, and metric diffs before updating
committed results. Cache/provenance fixes in v2.0.0 have not yet been followed
by a live-network rerun, so current external metrics remain historical.

Live validation must not use `unversioned`. CPTAC downloads require provider
MD5; the committed historical locked manifest lacks that field, so use
`--refresh-manifest` to bind the same file IDs to reviewed live GDC metadata
before any new download. `--offline`/`--cache-only` forbids network fallback and
requires an already-valid cache. Use a fresh live output directory so protected
locked manifests/models cannot alias managed cache or result paths. Derived
outputs publish from staging with the run manifest last; cache files remain
individually atomic and may be retained after a later failed step.

## Scope guardrails

The release remains validated for GDC STAR-Counts-style `log2(TPM+1)` bulk
RNA-seq tumor-vs-adjacent-normal contrasts. Do not broaden claims to Toil/RSEM,
GTEx, GEO, raw counts, FPKM, microarray, single-cell, spatial data, cancer
screening, diagnosis, or patient management without separate validation.
