# Publication Checklist

Use this checklist before making the repository or a release public.

Status for `v1.1.8-gdc-starcounts`: completed for the public GitHub release on
2026-07-09. Re-run this checklist before any future public release.

## Repository

- [x] Confirm the GitHub repository visibility is intentional.
- [x] Confirm CI passes on the hosted repository.
- [x] Confirm `.gitignore` excludes large training artifacts and caches.
- [x] Confirm no private credentials, tokens, or unpublished data are present.
- [x] Confirm `LICENSE`, `NOTICE.md`, `CITATION.cff`, `.zenodo.json`, and
      `codemeta.json` are present.
- [x] Run `python audit_publication_readiness.py`.

## Release

- [x] Confirm `VERSION` matches `RELEASE_METADATA.json`.
- [x] Run `python build_release_lite.py --smoke --timeout-seconds 300`.
- [x] Run `python run_release_acceptance.py --timeout-seconds 300`.
- [x] Run `python validate_release_lite.py --release-dir release-lite --zip
      tcga-tumor-normal-release-lite.zip --source-root . --artifacts
      RELEASE_ARTIFACTS.json`.
- [x] Confirm `RELEASE_ARTIFACTS.json` records the release zip SHA256.
- [x] Upload `tcga-tumor-normal-release-lite.zip` as the GitHub Release asset.
- [x] Paste `GITHUB_RELEASE_v1.1.8.md` as the GitHub Release body.

## Scientific Framing

- [x] Describe the model as research-only.
- [x] State the validated input boundary: GDC STAR-Counts-style `log2(TPM+1)`.
- [x] State that direct hard calls on Toil/RSEM, GTEx, GEO, raw counts, FPKM,
      microarray, single-cell, or spatial data are not validated.
- [x] Include the main internal, LOCO, CPTAC/GDC, and cross-platform boundary
      results.
- [x] Include TCGA/CPTAC/GTEx/GDC/UCSC Xena citation requirements in any paper
      or formal manuscript.
