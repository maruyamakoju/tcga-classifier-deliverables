# v2.1.0-gdc-starcounts

Type-safety release. The public `tcga_rnaseq` scoring library is now fully
type-annotated and ships a PEP 561 `py.typed` marker, and a mypy gate runs in
continuous integration. No scoring code, deployed weights, bundled payload logic,
or headline metric changed from v2.0.x — the runtime API and behavior are
identical, and the library API version stays `3.0.0` because the annotations are
additive.

## What changed

- **Typed public API** — every function across `io`, `align`, `score`, `metrics`,
  and `validation` now carries type annotations. The scoring layer stays
  numpy + pandas only; numpy ships `py.typed`, pandas is treated as untyped rather
  than pulling in a stub package.
- **`py.typed` marker** — the bundled `tcga_rnaseq` package now advertises its
  inline types (PEP 561), so a downstream project type-checking against the
  library picks up the real signatures instead of `Any`.
- **mypy CI gate** — a scoped `[tool.mypy]` config (`disallow_untyped_defs`,
  strict-equality, no implicit optional) checks the public package on every push
  and pull request, alongside Ruff and the full pytest suite.

Annotating the code surfaced and fixed a handful of latent type inconsistencies
(an intercept re-typed across the binary/multiclass branches, a set re-bound to a
list, and a couple of unannotated containers) — all behavior-preserving.

## Validation

- `python -m mypy`: success, no issues in 6 source files.
- `python -m pytest -q`: 413 passed.
- `python -m ruff check .`: passed.
- Deterministic build, release acceptance, and cross-platform CI: green.

## Security boundary for the ZIP

Do not extract or execute the downloaded asset until its SHA-256 matches the
trusted value published below. Acceptance requires the digest explicitly:

```bash
python validate_zip_bundle.py tcga-tumor-normal-release-lite.zip --expected-sha256 2b8ee4456564aa73de5d664752369174e762894e2cddbc487d997c769c7a8e15
```

Without a trusted digest, `--skip-acceptance` is limited to non-executing
structural inspection.

## Release asset

- Asset: `tcga-tumor-normal-release-lite.zip`
- SHA256: `2b8ee4456564aa73de5d664752369174e762894e2cddbc487d997c769c7a8e15`
- Size: `943291` bytes
- Zip entries: `76`

These asset values are taken from the generated `RELEASE_ARTIFACTS.json` after the
deterministic builder and validators passed.
