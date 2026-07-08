# Publication Checklist

Use this checklist before making the repository or a release public.

## Repository

- [ ] Confirm the GitHub repository visibility is intentional.
- [ ] Confirm CI passes on the hosted repository.
- [ ] Confirm `.gitignore` excludes large training artifacts and caches.
- [ ] Confirm no private credentials, tokens, or unpublished data are present.
- [ ] Confirm `LICENSE`, `NOTICE.md`, and `CITATION.cff` are present.
- [ ] Run `python audit_publication_readiness.py`.

## Release

- [ ] Confirm `VERSION` matches `RELEASE_METADATA.json`.
- [ ] Run `python build_release_lite.py --smoke --timeout-seconds 300`.
- [ ] Run `python run_release_acceptance.py --timeout-seconds 300`.
- [ ] Confirm `RELEASE_ARTIFACTS.json` records the release zip SHA256.
- [ ] Upload `tcga-tumor-normal-release-lite.zip` as the GitHub Release asset.
- [ ] Paste `GITHUB_RELEASE_v1.1.1.md` as the GitHub Release body.

## Scientific Framing

- [ ] Describe the model as research-only.
- [ ] State the validated input boundary: GDC STAR-Counts-style `log2(TPM+1)`.
- [ ] State that direct hard calls on Toil/RSEM, GTEx, GEO, raw counts, FPKM,
      microarray, single-cell, or spatial data are not validated.
- [ ] Include the main internal, LOCO, CPTAC/GDC, and cross-platform boundary
      results.
- [ ] Include TCGA/CPTAC/GTEx/GDC/UCSC Xena citation requirements in any paper
      or formal manuscript.
