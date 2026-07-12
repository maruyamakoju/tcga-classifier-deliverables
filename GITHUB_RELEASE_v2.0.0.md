# v2.0.0-gdc-starcounts

Safety, provenance, and deterministic reproducibility release for the
research-only GDC STAR-Counts tumor-vs-adjacent-normal classifier. The shipped
weights and committed headline metrics are unchanged. This is a major release
because the public `tcga_rnaseq` API is now 3.0.0 and intentionally fails closed
on unsafe or ambiguous contracts.

## What changed

- Fully validates model arrays, gene IDs, sample IDs, matched values, thresholds,
  and public output paths. Each output file is written atomically; the workflow
  records documented stop states and writes its manifest last rather than
  claiming an all-or-nothing directory transaction.
- Preserves full-precision logistic scores and clarifies that
  `tumor_probability` is a model score, not clinical risk or a calibrated
  diagnostic probability.
- Makes cohort adaptation default to `none`. Adapted modes are explicit,
  experimental, transductive and composition-dependent opt-ins; they require at
  least `--min-samples` (default 20), and scores from separately adapted batches
  are not comparable.
- Labels calibration metrics as apparent/resubstitution estimates computed on
  the same samples used to select the threshold, with a warning when either
  class contains fewer than 10 samples.
- Adds exact float32 cancer-type and float64 binary/LOCO reproduction paths.
- Adds locked external-validation cohorts, semantic cache fingerprints, content
  hashes, atomic caches, and run manifests. No post-fix live-network rerun was
  performed; committed CPTAC/GTEx/Toil metrics remain a historical snapshot.
- Documents that LOCO does not remove project/procurement/center/batch
  confounding and literature consistency is not causal-mechanism proof.
- Supports Python 3.11+ for lightweight scoring and tests 3.11/3.13 scoring in
  CI. Exact shipped-model reproduction is pinned to Python 3.11, NumPy 1.26.4,
  pandas 2.3.3, SciPy 1.15.3, and scikit-learn 1.8.0; the Python
  3.13/scikit-learn 1.9 refit drifted, so golden tolerances are not relaxed to
  make that noncanonical stack pass.
- Separates role-specific dependency profiles, adds the exact-pinned canonical
  training profile, runs Ruff/full pytest, and validates deterministic release
  behavior across Windows, Linux, and macOS.
- Builds a canonical, byte-reproducible ZIP through staged atomic publication
  and a non-mutating `--check` drift mode.

## Security boundary for the ZIP

Do not extract or execute the downloaded asset until its SHA-256 matches the
trusted value published below. Acceptance requires the digest explicitly:

```bash
python validate_zip_bundle.py tcga-tumor-normal-release-lite.zip --expected-sha256 a23301d6d194f4a91f7a2d59ed427749d92eca6a26da0b84e17400897fe86b6a
```

Without a trusted digest, `--skip-acceptance` is limited to non-executing
structural inspection.

## Validation

- `python -m pytest -q`: 412 passed.
- `python -m ruff check .`: passed.
- `python run_smoke_tests.py`: passed.
- Deterministic build and local artifact checks passed; hosted cross-platform CI
  runs on the publication gate.

## Release asset

- Asset: `tcga-tumor-normal-release-lite.zip`
- SHA256: `a23301d6d194f4a91f7a2d59ed427749d92eca6a26da0b84e17400897fe86b6a`
- Size: `937046` bytes
- Zip entries: `75`

These asset values are taken from the generated `RELEASE_ARTIFACTS.json` after the
deterministic builder and validators passed.
