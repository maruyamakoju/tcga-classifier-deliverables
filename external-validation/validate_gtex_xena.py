#!/usr/bin/env python3
"""Cross-platform normal-tissue check on GTEx via UCSC Xena Toil.

This validation asks a narrower question than the CPTAC test: if healthy GTEx
normal tissues are scored by the TCGA/GDC-trained classifier, how often are they
called tumor? It is not an AUC benchmark because GTEx contributes normals only.

The Toil matrix stores log2(TPM + 0.001). Before scoring, values are converted
back to TPM and then to the model's expected log2(TPM + 1) scale.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import os
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from tcga_rnaseq import load_lr_model, score_binary_dataframe  # noqa: E402
from tcga_rnaseq.align import strip_version  # noqa: E402

PHENOTYPE_URL = "https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz"
GTEX_TPM_URL = "https://toil.xenahubs.net/download/gtex_RSEM_gene_tpm.gz"


def read_gzip_tsv_from_url(url: str) -> pd.DataFrame:
    response = requests.get(url, timeout=120)
    response.raise_for_status()
    payload = gzip.decompress(response.content)
    return pd.read_csv(io.BytesIO(payload), sep="\t", encoding="latin-1")


def atomic_write_bytes(path: Path, write_fn) -> None:
    """Call write_fn(tmp_path) then atomically replace path with tmp_path.

    Protects a concurrently-reading process (or a rerun after this one is
    killed mid-write) from ever observing a truncated cache file.
    """
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    write_fn(tmp_path)
    os.replace(tmp_path, path)


def load_or_download_phenotype(path: Path, refresh: bool) -> pd.DataFrame:
    if path.exists() and not refresh:
        return pd.read_csv(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    phenotype = read_gzip_tsv_from_url(PHENOTYPE_URL)
    atomic_write_bytes(path, lambda tmp: phenotype.to_csv(tmp, index=False))
    return phenotype


def choose_gtex_samples(phenotype: pd.DataFrame, n_per_site: int,
                        min_site_n: int, seed: int) -> pd.DataFrame:
    gtex = phenotype[phenotype["_study"] == "GTEX"].copy()
    gtex = gtex[gtex["_sample_type"].fillna("").str.contains("Normal", case=False)]
    sampled = []
    for site, part in sorted(gtex.groupby("_primary_site")):
        if len(part) < min_site_n:
            continue
        n = min(n_per_site, len(part))
        sampled.append(part.sample(n=n, random_state=seed).copy())
    if not sampled:
        raise ValueError("No GTEx sites met sampling criteria")
    out = pd.concat(sampled, ignore_index=True)
    return out.sample(frac=1, random_state=seed).reset_index(drop=True)


def xena_log2_tpm001_to_model_scale(values: np.ndarray) -> np.ndarray:
    tpm = np.exp2(values) - 0.001
    tpm = np.clip(tpm, 0, None)
    return np.log2(tpm + 1.0)


def cache_fingerprint(sample_ids: list[str], selected_genes: list[str], matrix_url: str) -> str:
    """Fingerprint of everything that determines extract_matrix_from_xena's output.

    A cache keyed only on a fixed output path (the prior behavior) silently
    returns stale data -- for the wrong sample_ids, or extracted against an
    older model's selected_genes -- if the script is rerun with different
    sampling parameters or a different --weights model without --refresh.
    """
    digest = hashlib.sha256()
    digest.update(matrix_url.encode("utf-8"))
    for sample_id in sample_ids:
        digest.update(b"\x00" + str(sample_id).encode("utf-8"))
    for gene in selected_genes:
        digest.update(b"\x01" + str(gene).encode("utf-8"))
    return digest.hexdigest()


def _cache_meta_path(cache_path: Path) -> Path:
    return cache_path.with_suffix(cache_path.suffix + ".meta.json")


def _load_cached_matrix(cache_path: Path, fingerprint: str) -> pd.DataFrame | None:
    meta_path = _cache_meta_path(cache_path)
    if not cache_path.exists() or not meta_path.exists():
        return None
    try:
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if meta.get("fingerprint") != fingerprint:
        print(f"[xena] cache at {cache_path} is stale (sample_ids/selected_genes/url changed "
              "since it was written); re-extracting", file=sys.stderr)
        return None
    try:
        return pd.read_pickle(cache_path)
    except (OSError, EOFError, pickle.UnpicklingError):
        print(f"[xena] cache at {cache_path} could not be read; re-extracting", file=sys.stderr)
        return None


def extract_matrix_from_xena(sample_ids: list[str], selected_genes: list[str],
                             cache_path: Path, refresh: bool,
                             matrix_url: str = GTEX_TPM_URL) -> pd.DataFrame:
    fingerprint = cache_fingerprint(sample_ids, selected_genes, matrix_url)
    if not refresh:
        cached = _load_cached_matrix(cache_path, fingerprint)
        if cached is not None:
            return cached

    wanted_bases = {strip_version(gene): gene for gene in selected_genes}
    selected_set = set(selected_genes)
    rows = {}
    collisions = []
    header = None
    sample_indices = None

    print(f"[xena] streaming {matrix_url}", file=sys.stderr)
    with requests.get(matrix_url, stream=True, timeout=120) as response:
        response.raise_for_status()
        with gzip.GzipFile(fileobj=response.raw) as gz:
            for line_no, raw in enumerate(gz):
                line = raw.decode("utf-8", errors="replace").rstrip("\n")
                parts = line.split("\t")
                if line_no == 0:
                    header = parts
                    sample_to_idx = {sample: i for i, sample in enumerate(header)}
                    missing_samples = [s for s in sample_ids if s not in sample_to_idx]
                    if missing_samples:
                        raise ValueError(f"{len(missing_samples)} sampled GTEx IDs missing from matrix")
                    sample_indices = [sample_to_idx[s] for s in sample_ids]
                    continue

                gene_id = parts[0]
                target = gene_id if gene_id in selected_set else wanted_bases.get(strip_version(gene_id))
                if target is None:
                    continue
                if target in rows:
                    # A second raw row maps to the same selected gene after
                    # version-stripping (duplicate/ambiguous annotation row):
                    # keep the first value seen instead of silently
                    # overwriting it, and surface the collision.
                    collisions.append(gene_id)
                    continue
                vals = np.array([float(parts[i]) for i in sample_indices], dtype=float)
                rows[target] = xena_log2_tpm001_to_model_scale(vals)
                if len(rows) % 250 == 0:
                    print(f"[gtex] extracted {len(rows)}/{len(selected_genes)} genes",
                          file=sys.stderr)
                if len(rows) == len(selected_genes):
                    break

    if not rows:
        raise ValueError(
            f"No selected genes matched any row in {matrix_url}; the download may be "
            "empty, truncated, or in an unexpected format."
        )
    if collisions:
        print(f"[xena] WARNING: {len(collisions)} gene rows collided with an already-matched "
              f"selected gene after version-stripping and were skipped: {collisions[:5]}",
              file=sys.stderr)

    matrix = pd.DataFrame(rows, index=sample_ids)
    matrix = matrix.reindex(columns=selected_genes)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    atomic_write_bytes(cache_path, matrix.to_pickle)
    _cache_meta_path(cache_path).write_text(
        json.dumps({"fingerprint": fingerprint}, indent=2) + "\n", encoding="utf-8"
    )
    return matrix


def require_complete_merge(predictions: pd.DataFrame, expected_n: int, score_col: str,
                           context: str) -> None:
    """Raise if a label/score merge silently dropped, duplicated, or left rows unscored.

    A left-merge that doesn't match every row (e.g. a stale/misaligned cache,
    per cache_fingerprint's docstring) produces NaN scores that downstream
    aggregate metrics (means, rates) compute over silently, without ever
    raising -- a 0% or 100% "false positive rate" from an all-NaN column is a
    real, previously-possible failure mode this guards against.
    """
    if len(predictions) != expected_n:
        raise ValueError(
            f"{context}: merge produced {len(predictions)} rows, expected {expected_n} "
            "(duplicate or missing join keys)"
        )
    n_missing = int(predictions[score_col].isna().sum())
    if n_missing:
        raise ValueError(
            f"{context}: {n_missing}/{expected_n} samples have no score after the merge "
            "(sample ID mismatch between the extracted matrix and the sampled manifest)"
        )


def summarize(predictions: pd.DataFrame, threshold: float) -> dict:
    tumor_calls = int((predictions["call"] == "tumor").sum())
    n = len(predictions)
    return {
        "n": n,
        "threshold": threshold,
        "tumor_calls": tumor_calls,
        "normal_calls": int((predictions["call"] == "normal").sum()),
        "false_positive_rate": tumor_calls / n,
        "median_tumor_probability": float(predictions["tumor_probability"].median()),
        "p90_tumor_probability": float(predictions["tumor_probability"].quantile(0.90)),
        "p95_tumor_probability": float(predictions["tumor_probability"].quantile(0.95)),
        "max_tumor_probability": float(predictions["tumor_probability"].max()),
    }


def threshold_sweep(predictions: pd.DataFrame) -> pd.DataFrame:
    rows = []
    scores = predictions["tumor_probability"].to_numpy()
    for threshold in [0.5, 0.75, 0.9, 0.95, 0.99, 0.999, 0.9999,
                      0.999975, 0.99999, 0.999999]:
        tumor_calls = int((scores >= threshold).sum())
        rows.append({
            "threshold": threshold,
            "tumor_calls": tumor_calls,
            "normal_calls": int(len(scores) - tumor_calls),
            "false_positive_rate": tumor_calls / len(scores),
        })
    return pd.DataFrame(rows)


def write_report(path: Path, summary: dict, predictions: pd.DataFrame,
                 by_site: pd.DataFrame, sweep: pd.DataFrame) -> None:
    site_rows = "\n".join(
        f"| {row['_primary_site']} | {int(row['n'])} | {int(row['tumor_calls'])} | "
        f"{row['false_positive_rate']:.3f} | {row['median_tumor_probability']:.4f} | "
        f"{row['max_tumor_probability']:.4f} |"
        for _, row in by_site.iterrows()
    )
    sweep_rows = "\n".join(
        f"| {row.threshold:.6f} | {int(row.tumor_calls)} | "
        f"{row.false_positive_rate:.3f} |"
        for row in sweep.itertuples(index=False)
    )
    top = predictions.sort_values("tumor_probability", ascending=False).head(20)
    top_rows = "\n".join(
        f"| {row['sample']} | {row['_primary_site']} | {row['primary disease or tissue']} | "
        f"{row['tumor_probability']:.6f} | {row['call']} |"
        for _, row in top.iterrows()
    )
    text = f"""# GTEx normal-tissue cross-platform check

## Data source

- Source: UCSC Xena Toil RNA-seq recompute compendium
- Dataset: `gtex_RSEM_gene_tpm`
- Phenotype table: `TcgaTargetGTEX_phenotype`
- Samples scored: {summary['n']} GTEx normal samples, stratified by primary site
- Input transform: Xena log2(TPM+0.001) -> TPM -> log2(TPM+1)
- Model: bundled logistic regression from `deployable_lr_weights.npz`

## Result at threshold {summary['threshold']}

| Metric | Value |
|---|---:|
| Samples | {summary['n']} |
| Normal calls | {summary['normal_calls']} |
| Tumor calls / false positives | {summary['tumor_calls']} |
| False positive rate | {summary['false_positive_rate']:.4f} |
| Median tumor probability | {summary['median_tumor_probability']:.4f} |
| p90 tumor probability | {summary['p90_tumor_probability']:.4f} |
| p95 tumor probability | {summary['p95_tumor_probability']:.4f} |
| Max tumor probability | {summary['max_tumor_probability']:.4f} |

## Threshold sensitivity

| Threshold | Tumor calls | False positive rate |
|---:|---:|---:|
{sweep_rows}

## Per-site summary

| Primary site | n | Tumor calls | False positive rate | Median tumor probability | Max tumor probability |
|---|---:|---:|---:|---:|---:|
{site_rows}

## Highest-scoring GTEx normals

| Sample | Primary site | Tissue | Tumor probability | Call |
|---|---|---|---:|---|
{top_rows}

## Interpretation

This is a stricter platform-transfer check than the CPTAC/GDC validation because
it uses GTEx/Toil normal tissues rather than GDC STAR-Counts. It only measures
normal-sample false positives, not tumor-vs-normal AUC. A high false positive
rate would indicate that the TCGA/GDC-trained threshold does not transfer cleanly
to GTEx/Toil. The companion TCGA Toil/RSEM check shows that even TCGA samples
require an extreme threshold shift on this pipeline, so this result should be
treated as a pipeline/domain-transfer failure of the deployed GDC STAR-Counts
model rather than evidence that GTEx tissues are biologically tumor-like.
"""
    path.write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-per-site", type=int, default=20)
    parser.add_argument("--min-site-n", type=int, default=20)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument("--out-dir", default=str(ROOT / "external-validation" / "gtex_xena"))
    parser.add_argument("--weights", default=str(ROOT / "deployable_lr_weights.npz"))
    parser.add_argument("--max-invalid-cell-fraction", type=float, default=0.0,
                        help=("maximum allowed missing, non-numeric, NaN, or infinite cells "
                              "among matched model genes before failing (default 0)"))
    parser.add_argument("--allow-invalid-values", action="store_true",
                        help=("warn instead of failing when matched model-gene cells are "
                              "missing, non-numeric, NaN, or infinite"))
    args = parser.parse_args()
    if not 0 <= args.max_invalid_cell_fraction <= 1:
        parser.error("--max-invalid-cell-fraction must be between 0 and 1")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    phenotype_path = out_dir / "TcgaTargetGTEX_phenotype.csv"
    sampled_path = out_dir / "sampled_gtex_manifest.csv"
    matrix_path = out_dir / "gtex_selected_genes_model_scale.pkl"

    phenotype = load_or_download_phenotype(phenotype_path, args.refresh)
    sampled = choose_gtex_samples(phenotype, args.n_per_site, args.min_site_n, args.seed)
    sampled.to_csv(sampled_path, index=False)
    print(f"[gtex] sampled {len(sampled)} GTEx normals across "
          f"{sampled['_primary_site'].nunique()} primary sites", file=sys.stderr)

    model = load_lr_model(args.weights)
    selected_genes = model["genes"].astype(str).tolist()
    matrix = extract_matrix_from_xena(sampled["sample"].tolist(), selected_genes,
                                      matrix_path, args.refresh)
    n_matched = int(matrix.notna().any(axis=0).sum())
    if n_matched < len(selected_genes):
        print(f"[gtex] WARNING: matched {n_matched}/{len(selected_genes)} genes",
              file=sys.stderr)

    scored, _, _, alignment_report = score_binary_dataframe(
        model,
        matrix,
        threshold=args.threshold,
        max_invalid_cell_fraction=args.max_invalid_cell_fraction,
        allow_invalid_values=args.allow_invalid_values,
        return_alignment_report=True,
    )
    if alignment_report["invalid_matched_cells"]:
        print(
            f"[gtex] WARNING: imputed {alignment_report['invalid_matched_cells']} "
            "invalid matched cells",
            file=sys.stderr,
        )
    predictions = sampled.merge(scored, left_on="sample", right_on="sample", how="left")
    require_complete_merge(predictions, len(sampled), "tumor_probability", "[gtex]")
    predictions_path = out_dir / "gtex_predictions.csv"
    predictions.to_csv(predictions_path, index=False)

    summary = summarize(predictions, args.threshold)
    summary_path = out_dir / "gtex_summary.csv"
    pd.DataFrame([summary]).to_csv(summary_path, index=False)

    sweep = threshold_sweep(predictions)
    sweep_path = out_dir / "gtex_threshold_sweep.csv"
    sweep.to_csv(sweep_path, index=False)

    by_site = predictions.groupby("_primary_site").agg(
        n=("sample", "size"),
        tumor_calls=("call", lambda s: int((s == "tumor").sum())),
        false_positive_rate=("call", lambda s: float((s == "tumor").mean())),
        median_tumor_probability=("tumor_probability", "median"),
        max_tumor_probability=("tumor_probability", "max"),
    ).reset_index().sort_values(["false_positive_rate", "max_tumor_probability"],
                                ascending=False)
    by_site_path = out_dir / "gtex_per_site_summary.csv"
    by_site.to_csv(by_site_path, index=False)

    report_path = out_dir / "GTEX_NORMAL_VALIDATION.md"
    write_report(report_path, summary, predictions, by_site, sweep)

    print(f"[gtex] summary -> {summary_path}", file=sys.stderr)
    print(f"[gtex] per-site -> {by_site_path}", file=sys.stderr)
    print(f"[gtex] predictions -> {predictions_path}", file=sys.stderr)
    print(f"[gtex] report -> {report_path}", file=sys.stderr)
    print(pd.DataFrame([summary]).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
