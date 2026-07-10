# v1.1.16-gdc-starcounts

Pickle expression input rejection release for the GDC STAR-Counts tumor-vs-normal
lightweight bundle. Model weights, training data, and headline validation
metrics are unchanged from v1.1.15.

## What changed

- `tcga_rnaseq.read_matrix()` now rejects `.pkl` expression matrices by default
  because unpickling user-controlled files can execute code.
- Public CLIs surface the `.pkl` rejection through normal argument errors instead
  of loading the file.
- Trusted internal cross-platform benchmark code that reads local pickle caches
  now opts in explicitly with `allow_pickle=True`.
- Updated public input-format documentation to require CSV, TSV, TXT, or
  Parquet expression matrices.
- Added core unit and release safety coverage for pickle input rejection.

## Validation

- `python audit_release_docs.py`
- `python build_release_lite.py --smoke --timeout-seconds 300`
- `python audit_publication_readiness.py`
- `python validate_release_lite.py --release-dir release-lite --zip tcga-tumor-normal-release-lite.zip --source-root . --artifacts RELEASE_ARTIFACTS.json`
- `python validate_zip_bundle.py tcga-tumor-normal-release-lite.zip`
- `python run_release_acceptance.py --timeout-seconds 300`
- `python -m pytest -q -rs`
- `git diff --check`

## Release asset

- Asset: `tcga-tumor-normal-release-lite.zip`
- SHA256: `c854cf435d6e626c96e695618a46ca642ed7107e1c8a0bcf3f591d20d1036d2d`
- Size: `313223` bytes
- Zip entries: `72`
