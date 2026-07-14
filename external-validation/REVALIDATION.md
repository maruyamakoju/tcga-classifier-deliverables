# Live external re-validation record

This is the canonical, provenance-backed record of post-fix **live** re-validation
runs: fresh fetches of each locked external cohort from the current provider
sources, scored by the current code, compared to the committed historical
metrics. The committed per-cohort artifacts remain the frozen historical
snapshot; this document records that the current code reproduces them live.

## Result — all three cohorts reproduced (July 2026)

| Cohort | Source revision | Headline metric | Live | Committed | Match |
| --- | --- | --- | --- | --- | --- |
| CPTAC-3 / GDC | GDC Data Release 45.0 (2025-12-04) | AUC (n=200) | **0.9886** | 0.9886 | exact |
| CPTAC-3 / GDC | " | accuracy @ 0.5 | 0.960 | 0.955 | +1 boundary sample |
| TCGA-Toil / UCSC Xena | Xena Toil RSEM (accessed 2026-07-14) | AUC (n=200) | **0.9923** | 0.9923 | exact |
| TCGA-Toil / UCSC Xena | " | accuracy @ 0.5 | 0.515 | 0.515 | exact |
| GTEx normals / UCSC Xena | Xena GTEx RSEM (accessed 2026-07-14) | FPR @ 0.5 (n=540) | **0.9963** | 0.9963 | exact |

All three headline metrics reproduce. The single CPTAC accuracy difference is one
threshold-boundary sample (the AUC, a ranking metric, is identical). The GTEx and
TCGA-Toil "poor" numbers (high false-positive rate; shifted 0.5 cutoff) are the
**expected documented cross-platform boundary** — the model is GDC STAR-Counts
specific and is not a hard-call classifier on Toil/RSEM or GTEx.

## Provenance

Each run's full provenance — source revision, semantic cache fingerprints,
per-file content hashes, code source hashes, and environment — is archived
alongside the reproduced summary and threshold sweep:

- `cptac_gdc/revalidation/run_manifest.json`
- `tcga_toil_xena/revalidation/run_manifest.json`
- `gtex_xena/revalidation/run_manifest.json`

The large downloaded expression matrices are not committed (they are re-fetchable
from the providers at the recorded revisions).

## What the live runs surfaced

Reaching a clean live run required fixing three real, over-strict fail-closed
checks that rejected legitimate current-source data (none affect the shipped GDC
scoring path):

- CPTAC manifest validator rejected files mapped to several same-case,
  same-`sample_type` biospecimen IDs, though the tumor/normal label is
  unambiguous (fixed in v2.2.0).
- `provenance.contained_cache_path` false-rejected legitimate long Windows cache
  directories at the `MAX_PATH` boundary (fixed in v2.2.0).
- The Xena `log2(TPM+0.001)` decode rejected the tiny negative TPM that a
  truly-zero-expression gene produces from storage rounding (~1e-8); the
  tolerance was widened from `1e-9` to `1e-6` and the value clamped to zero
  (fixed in v2.4.0).

## Reproducing this

See the commands and prerequisites in `REPRODUCIBILITY.md`. In short, run each
`validate_*.py` with the committed `--sample-manifest`, a concrete
`--source-revision`, `--refresh`/`--refresh-expression-cache`, and a fresh
`--out-dir`, then compare the fresh `*_summary.csv` to the committed values within
a tolerance of 0.01.
