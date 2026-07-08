# Contributing

This repository is maintained as a release-quality research artifact. Changes
should preserve the public scoring contracts unless the release version and
documentation are updated together.

## Development Checks

Run these before proposing changes:

```bash
python -m pytest -q -rs
python run_release_acceptance.py --timeout-seconds 300
python validate_release_lite.py --release-dir release-lite --zip tcga-tumor-normal-release-lite.zip --source-root . --artifacts RELEASE_ARTIFACTS.json
```

When changing release payload files, rebuild the bundle:

```bash
python build_release_lite.py --smoke --timeout-seconds 300
```

## Compatibility Rules

- Keep `score_tumor_normal.py` output columns as `sample,tumor_probability,call`.
- Keep workflow output filenames and JSON top-level keys stable unless this is
  a planned breaking release.
- Do not add heavy runtime dependencies to the lightweight path.
- Do not commit full training matrices, large pickle checkpoints, or external
  validation cache files.
- Preserve the documented deployment boundary: GDC STAR-Counts-style
  `log2(TPM+1)` tumor-vs-adjacent-normal research use.

## Release Updates

For a new release, update:

- `VERSION`
- `RELEASE_METADATA.json`
- `RELEASE_NOTES.md`
- version headers in the release documentation
- `release-lite/`, `tcga-tumor-normal-release-lite.zip`, and
  `RELEASE_ARTIFACTS.json` via `build_release_lite.py`
