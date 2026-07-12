# Publication Checklist

Use this checklist before making the repository or a release public.

Status for `v2.0.0-gdc-starcounts`: release candidate on `2026-07-12`.
Unchecked publication and asset items must be completed after the final
deterministic build.

## Repository

- [x] Confirm the GitHub repository visibility is intentional.
- [ ] Confirm CI passes on the hosted repository for the release commit.
- [x] Confirm `.gitignore` excludes large training artifacts and caches.
- [x] Confirm no private credentials, tokens, or unpublished data are present.
- [x] Confirm `LICENSE`, `NOTICE.md`, `CITATION.cff`, `.zenodo.json`, and
      `codemeta.json` are present.
- [ ] Run `python audit_publication_readiness.py` after staging the new release
      note and final deterministic artifacts.

## Release

- [x] Confirm `VERSION` matches `RELEASE_METADATA.json`.
- [ ] Run `python build_release_lite.py --smoke --timeout-seconds 300`.
- [ ] Run `python build_release_lite.py --check --smoke --timeout-seconds 300`
      and confirm the check is non-mutating.
- [ ] Run `python run_release_acceptance.py --timeout-seconds 300`.
- [ ] Run `python validate_release_lite.py --release-dir release-lite --zip
      tcga-tumor-normal-release-lite.zip --source-root . --artifacts
      RELEASE_ARTIFACTS.json`.
- [ ] Confirm `RELEASE_ARTIFACTS.json` records the final release ZIP SHA256,
      byte size, entry count, and trusted-digest acceptance command.
- [ ] Run `python validate_zip_bundle.py tcga-tumor-normal-release-lite.zip
      --expected-sha256 <trusted-final-sha256>`.
- [ ] Replace every explicit TBD field in `GITHUB_RELEASE_v2.0.0.md` from the
      final generated sidecar.
- [ ] Upload `tcga-tumor-normal-release-lite.zip` as the GitHub Release asset.
- [ ] Paste `GITHUB_RELEASE_v2.0.0.md` as the GitHub Release body.

## Scientific Framing

- [x] Describe the model as research-only.
- [x] State the validated input boundary: GDC STAR-Counts-style `log2(TPM+1)`.
- [x] State that direct hard calls on Toil/RSEM, GTEx, GEO, raw counts, FPKM,
      microarray, single-cell, or spatial data are not validated.
- [x] Include the main internal, LOCO, CPTAC/GDC, and cross-platform boundary
      results.
- [x] State that `tumor_probability` is a logistic model score, not clinical
      risk or a calibrated diagnostic probability.
- [x] State that adaptation defaults to `none`; adapted modes are explicit,
      transductive, composition-dependent opt-ins with at least 20 samples by
      default, and separately adapted scores are not comparable.
- [x] Label threshold-calibration metrics as same-sample
      apparent/resubstitution estimates and disclose the per-class `<10`
      warning.
- [x] State that LOCO does not remove project/procurement/center/batch
      confounding and literature consistency is not causal-mechanism proof.
- [x] State that external cache/provenance fixes have no post-fix live rerun
      and committed external metrics are a historical snapshot.
- [x] Include TCGA/CPTAC/GTEx/GDC/UCSC Xena citation requirements in any paper
      or formal manuscript.
