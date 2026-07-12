# Publication Checklist

Use this checklist before making the repository or a release public.

Status for `v2.2.0-gdc-starcounts`: released and public. Every item below has been
executed and verified: the Repository and Scientific Framing items are durably
true of the tree, and the Release items run automatically through the
deterministic builder, `run_release_acceptance.py`, and hosted CI on each version
tag (v2.0.0 through the current release). The one outstanding publication step is
**archival: mint a Zenodo DOI and backfill it into `CITATION.cff`, `.zenodo.json`,
and the `README` badge** — this requires enabling the GitHub-Zenodo integration,
which then archives the release and issues the DOI.

## Repository

- [x] Confirm the GitHub repository visibility is intentional.
- [x] Confirm CI passes on the hosted repository for the release commit.
- [x] Confirm `.gitignore` excludes large training artifacts and caches.
- [x] Confirm no private credentials, tokens, or unpublished data are present.
- [x] Confirm `LICENSE`, `NOTICE.md`, `CITATION.cff`, `.zenodo.json`, and
      `codemeta.json` are present.
- [x] Run `python audit_publication_readiness.py` after staging the new release
      note and final deterministic artifacts.

## Release

- [x] Confirm `VERSION` matches `RELEASE_METADATA.json`.
- [x] Run `python build_release_lite.py --smoke --timeout-seconds 300`.
- [x] Run `python build_release_lite.py --check --smoke --timeout-seconds 300`
      and confirm the check is non-mutating.
- [x] Run `python run_release_acceptance.py --timeout-seconds 300`.
- [x] Run `python validate_release_lite.py --release-dir release-lite --zip
      tcga-tumor-normal-release-lite.zip --source-root . --artifacts
      RELEASE_ARTIFACTS.json`.
- [x] Confirm `RELEASE_ARTIFACTS.json` records the final release ZIP SHA256,
      byte size, entry count, and trusted-digest acceptance command.
- [x] Run `python validate_zip_bundle.py tcga-tumor-normal-release-lite.zip
      --expected-sha256 <trusted-final-sha256>`.
- [x] Replace every explicit TBD field in `GITHUB_RELEASE_v2.2.0.md` from the
      final generated sidecar.
- [x] Upload `tcga-tumor-normal-release-lite.zip` as the GitHub Release asset.
- [x] Paste `GITHUB_RELEASE_v2.2.0.md` as the GitHub Release body.

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

## Archival & DOI

- [x] Confirm `.zenodo.json`, `CITATION.cff`, and `codemeta.json` describe the
      current release version and metadata.
- [ ] Enable the GitHub-Zenodo integration for the repository so a published
      release is archived automatically.
- [ ] Obtain the minted Zenodo DOI for the release (and the version-independent
      concept DOI).
- [ ] Backfill the DOI into `CITATION.cff` (`doi:` / `identifiers:`),
      `.zenodo.json`, and a `README` DOI badge, then cut a follow-up release so
      the archived record and the repository agree.
