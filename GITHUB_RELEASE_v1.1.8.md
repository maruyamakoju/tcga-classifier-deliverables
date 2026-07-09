# v1.1.8-gdc-starcounts

Lite documentation audit release for the GDC STAR-Counts tumor-vs-normal
lightweight bundle. Model weights, training data, and headline validation
metrics are unchanged from v1.1.7.

## What changed

- `audit_release_docs.py` now validates code-spanned local file and directory
  references, not only Markdown links and Python command references.
- Generated workflow outputs, release sidecars, and full-development-tree-only
  references are explicitly categorized so unexpected stale paths fail the
  documentation audit.
- `README.md` now separates lightweight-bundle contents from full development
  tree-only maintenance and historical-analysis assets.
- `REPORT.md` now marks LOCO, cross-platform adaptation benchmark,
  cancer-type classifier, and maintenance/regeneration artifacts as full
  development tree-only when they are not part of the lightweight zip.

## Validation

- `python audit_release_docs.py`
- `python build_release_lite.py --smoke --timeout-seconds 300`
- `python audit_publication_readiness.py`
- `python validate_release_lite.py --release-dir release-lite --zip tcga-tumor-normal-release-lite.zip --source-root . --artifacts RELEASE_ARTIFACTS.json`
- `python run_release_acceptance.py --timeout-seconds 300`
- `python validate_zip_bundle.py tcga-tumor-normal-release-lite.zip`
- `python -m pytest -q -rs`
- `git diff --check`

## Release asset

- Asset: `tcga-tumor-normal-release-lite.zip`
- SHA256: `0334150ee068e5dd6597f95b5a686509682a04f594d7c0d74259a9cb8522f6fa`
- Size: `308456` bytes
- Zip entries: `72`
