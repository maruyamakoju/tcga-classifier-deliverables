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
import sys
from pathlib import Path

import numpy as np
import pandas as pd
try:
    import requests
except ImportError:  # Keep pure-logic imports usable in the lightweight environment.
    requests = None

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from tcga_rnaseq import load_lr_model, predict_proba, validate_threshold  # noqa: E402
from tcga_rnaseq.align import strip_version  # noqa: E402
from provenance import (  # noqa: E402
    atomic_write_csv,
    atomic_write_path,
    atomic_write_text,
    cache_meta_path,
    group_audit,
    load_csv_cache,
    load_dataframe_cache,
    read_cache_metadata,
    publish_staged_files,
    scored_dataframe,
    semantic_fingerprint,
    sha256_file,
    staged_output_directory,
    utc_now,
    validate_identifier_column,
    validate_managed_paths,
    validate_source_revision,
    validate_unique_strings,
    write_csv_cache,
    write_dataframe_cache,
    write_run_manifest,
)

PHENOTYPE_URL = "https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz"
GTEX_TPM_URL = "https://toil.xenahubs.net/download/gtex_RSEM_gene_tpm.gz"
PHENOTYPE_DATASET_ID = "TcgaTargetGTEX_phenotype"
GTEX_DATASET_ID = "gtex_RSEM_gene_tpm"
PHENOTYPE_PARSER_VERSION = "phenotype-gzip-tsv-v2"
XENA_MATRIX_PARSER_VERSION = "xena-gzip-matrix-selected-genes-v3"
XENA_TRANSFORM_VERSION = "log2-tpm-plus-0.001_to_log2-tpm-plus-1-v1"


def _http_client():
    if requests is None:
        raise RuntimeError(
            "live external validation requires the optional 'requests' package"
        )
    return requests


def read_gzip_tsv_from_url(url: str) -> tuple[pd.DataFrame, dict[str, str | None]]:
    response = _http_client().get(url, timeout=120)
    response.raise_for_status()
    payload = gzip.decompress(response.content)
    source = {
        "url": url,
        "etag": response.headers.get("ETag") if hasattr(response, "headers") else None,
        "last_modified": (
            response.headers.get("Last-Modified") if hasattr(response, "headers") else None
        ),
        "content_length": (
            response.headers.get("Content-Length") if hasattr(response, "headers") else None
        ),
        "download_sha256": hashlib.sha256(response.content).hexdigest(),
    }
    return pd.read_csv(io.BytesIO(payload), sep="\t", encoding="latin-1"), source


def atomic_write_bytes(path: Path, write_fn) -> None:
    """Compatibility wrapper around the shared unique-temp atomic writer."""
    atomic_write_path(path, write_fn)


def _phenotype_cache_inputs(url: str, source_revision: str) -> dict:
    return {
        "cache_schema": "phenotype-csv-v2",
        "parser_version": PHENOTYPE_PARSER_VERSION,
        "source": {
            "dataset_id": PHENOTYPE_DATASET_ID,
            "url": url,
            "revision": source_revision,
        },
    }


def load_or_download_phenotype(
    path: Path,
    refresh: bool,
    source_revision: str = "unversioned",
    offline: bool = False,
) -> pd.DataFrame:
    inputs = _phenotype_cache_inputs(PHENOTYPE_URL, source_revision)
    fingerprint = semantic_fingerprint(inputs)
    required = ["sample", "_study", "_sample_type", "_primary_site"]
    if not refresh:
        cached = load_csv_cache(path, fingerprint=fingerprint, required_columns=required)
        if cached is not None:
            return cached
    if offline:
        raise ValueError(
            f"offline/cache-only mode requires a valid phenotype cache: {path}"
        )
    validate_source_revision(source_revision, live=True)
    path.parent.mkdir(parents=True, exist_ok=True)
    phenotype, source_metadata = read_gzip_tsv_from_url(PHENOTYPE_URL)
    missing = sorted(set(required) - set(phenotype.columns))
    if missing:
        raise ValueError(f"Xena phenotype is missing required columns: {missing}")
    write_csv_cache(
        path,
        phenotype,
        fingerprint=fingerprint,
        fingerprint_inputs=inputs,
        cache_kind="xena_phenotype",
        extra_metadata={"source_response": source_metadata},
    )
    return phenotype


def donor_id_from_sample(sample_id: str) -> str:
    value = str(sample_id).strip()
    if value.startswith("GTEX-"):
        return "-".join(value.split("-")[:2])
    if value.startswith("TCGA-"):
        return "-".join(value.split("-")[:3])
    return value


def _validate_sampling_count(value: int, name: str) -> int:
    if value <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return value


def choose_gtex_samples(
    phenotype: pd.DataFrame,
    n_per_site: int,
    min_site_n: int,
    seed: int,
    allow_donor_overlap: bool = False,
) -> pd.DataFrame:
    _validate_sampling_count(n_per_site, "n_per_site")
    _validate_sampling_count(min_site_n, "min_site_n")
    validate_identifier_column(phenotype, "sample", "Xena phenotype")
    gtex = phenotype[phenotype["_study"] == "GTEX"].copy()
    gtex = gtex[gtex["_sample_type"].fillna("").str.contains("Normal", case=False)]
    gtex["donor_id"] = gtex["sample"].map(donor_id_from_sample)
    sampled = []
    used_donors: set[str] = set()
    site_parts = []
    for site, part in gtex.groupby("_primary_site", sort=True):
        part = part.sort_values("sample").drop_duplicates("donor_id", keep="first")
        if len(part) < min_site_n:
            continue
        site_parts.append((site, part))
    # Constrain the rarest sites first so common tissues cannot consume their
    # donors before a small site is sampled.
    site_parts.sort(key=lambda item: (len(item[1]), str(item[0])))
    for site, part in site_parts:
        if not allow_donor_overlap:
            part = part[~part["donor_id"].isin(used_donors)]
        n = min(n_per_site, len(part))
        if n < min(n_per_site, min_site_n):
            raise ValueError(
                f"Only {n} unused GTEx donors remain for site {site!r}; pass "
                "--allow-donor-overlap to reproduce a sample-level panel"
            )
        sampled.append(part.sample(n=n, random_state=seed).copy())
        used_donors.update(sampled[-1]["donor_id"].astype(str))
    if not sampled:
        raise ValueError("No GTEx sites met sampling criteria")
    out = pd.concat(sampled, ignore_index=True)
    return out.sample(frac=1, random_state=seed).reset_index(drop=True)


def load_locked_sample_manifest(path: Path, study: str) -> pd.DataFrame:
    manifest = pd.read_csv(path, dtype={"sample": str}, keep_default_na=False)
    manifest = manifest.copy()
    required = {"_study", "_sample_type", "_primary_site"}
    if study == "GTEX":
        required.add("primary disease or tissue")
    missing = sorted(required - set(manifest.columns))
    if missing:
        raise ValueError(f"locked {study} sample manifest is missing columns: {missing}")
    manifest["sample"] = validate_identifier_column(
        manifest, "sample", f"locked {study} sample manifest"
    )
    studies = manifest["_study"].astype(str).str.strip()
    if studies.eq("").any() or not studies.eq(study).all():
        found = sorted(set(studies))
        raise ValueError(
            f"locked sample manifest must contain only non-empty _study={study!r}; "
            f"found {found}"
        )
    manifest["_study"] = studies
    sites = manifest["_primary_site"].astype(str).str.strip()
    if sites.eq("").any():
        raise ValueError(f"locked {study} sample manifest contains blank _primary_site")
    manifest["_primary_site"] = sites
    sample_types = manifest["_sample_type"].astype(str).str.strip()
    allowed_types = {
        "GTEX": {"Normal Tissue"},
        "TCGA": {"Solid Tissue Normal", "Primary Tumor"},
    }.get(study)
    if allowed_types is None:
        raise ValueError(f"unsupported locked Xena study: {study!r}")
    if not sample_types.isin(allowed_types).all():
        bad = sorted(set(sample_types[~sample_types.isin(allowed_types)]))
        raise ValueError(
            f"locked {study} sample manifest has unsupported _sample_type values: {bad}"
        )
    manifest["_sample_type"] = sample_types
    expected_prefix = f"{study}-"
    if not manifest["sample"].str.startswith(expected_prefix).all():
        raise ValueError(
            f"locked {study} sample IDs must start with {expected_prefix!r}"
        )
    if study == "GTEX":
        tissues = manifest["primary disease or tissue"].astype(str).str.strip()
        if tissues.eq("").any():
            raise ValueError("locked GTEX sample manifest contains blank tissue values")
        manifest["primary disease or tissue"] = tissues
    manifest["donor_id"] = manifest["sample"].map(donor_id_from_sample)
    return manifest


def xena_log2_tpm001_to_model_scale(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("Xena expression values must be finite")
    with np.errstate(over="ignore", invalid="ignore"):
        tpm = np.exp2(values) - 0.001
    if not np.isfinite(tpm).all():
        raise ValueError("Xena expression values overflow while converting to TPM")
    # Xena stores log2(TPM+0.001) rounded, so a truly zero-expression gene decodes
    # to a tiny negative TPM (observed ~-1e-8) purely from storage/float round-trip.
    # Reject only genuinely negative TPM (data corruption is orders of magnitude
    # larger); clamp the sub-tolerance round-trip noise to exactly zero below.
    tolerance = 1e-6
    if bool((tpm < -tolerance).any()):
        minimum = float(tpm.min())
        raise ValueError(
            "Xena log2(TPM+0.001) values decode below zero TPM "
            f"(minimum {minimum:.6g})"
        )
    tpm = np.maximum(tpm, 0.0)
    converted = np.log2(tpm + 1.0)
    if not np.isfinite(converted).all():
        raise ValueError("converted Xena expression values must be finite")
    return converted


def _xena_cache_inputs(
    sample_ids: list[str],
    selected_genes: list[str],
    matrix_url: str,
    model_sha256: str,
    source_identity: dict | None,
) -> dict:
    sample_ids = validate_unique_strings(sample_ids, "Xena sample IDs")
    selected_genes = validate_unique_strings(selected_genes, "selected model genes")
    bases = [strip_version(gene) for gene in selected_genes]
    if len(set(bases)) != len(bases):
        raise ValueError("selected model genes collide after removing Ensembl version suffixes")
    return {
        "cache_schema": "xena-selected-gene-matrix-parquet-v4",
        "parser_version": XENA_MATRIX_PARSER_VERSION,
        "transform_version": XENA_TRANSFORM_VERSION,
        "source": source_identity or {"url": matrix_url, "revision": "unversioned"},
        "matrix_url": matrix_url,
        "model_sha256": model_sha256,
        "sample_ids": [str(value) for value in sample_ids],
        "selected_genes": [str(value) for value in selected_genes],
    }


def cache_fingerprint(
    sample_ids: list[str],
    selected_genes: list[str],
    matrix_url: str,
    model_sha256: str = "unspecified",
    source_identity: dict | None = None,
) -> str:
    """Fingerprint of everything that determines extract_matrix_from_xena's output.

    A cache keyed only on a fixed output path (the prior behavior) silently
    returns stale data -- for the wrong sample_ids, or extracted against an
    older model's selected_genes -- if the script is rerun with different
    sampling parameters or a different --weights model without --refresh.
    """
    return semantic_fingerprint(
        _xena_cache_inputs(
            sample_ids, selected_genes, matrix_url, model_sha256, source_identity
        )
    )


def _cache_meta_path(cache_path: Path) -> Path:
    return cache_meta_path(cache_path)


def _load_cached_matrix(
    cache_path: Path,
    fingerprint: str,
    sample_ids: list[str] | None = None,
    selected_genes: list[str] | None = None,
) -> pd.DataFrame | None:
    cached = load_dataframe_cache(
        cache_path,
        fingerprint=fingerprint,
        expected_index=sample_ids,
        expected_columns=selected_genes,
        columns_may_be_subset=selected_genes is not None,
        reject_all_missing_columns=True,
    )
    if cached is None and cache_path.exists():
        print(
            f"[xena] cache at {cache_path} failed fingerprint/content/axis validation; "
            "re-extracting",
            file=sys.stderr,
        )
    return cached


def extract_matrix_from_xena(sample_ids: list[str], selected_genes: list[str],
                             cache_path: Path, refresh: bool,
                             matrix_url: str = GTEX_TPM_URL,
                             model_sha256: str = "unspecified",
                             source_identity: dict | None = None,
                             offline: bool = False) -> pd.DataFrame:
    inputs = _xena_cache_inputs(
        sample_ids, selected_genes, matrix_url, model_sha256, source_identity
    )
    fingerprint = semantic_fingerprint(inputs)
    if not refresh:
        cached = _load_cached_matrix(cache_path, fingerprint, sample_ids, selected_genes)
        if cached is not None:
            return cached
    if offline:
        raise ValueError(
            f"offline/cache-only mode requires a valid Xena expression cache: {cache_path}"
        )
    revision = (source_identity or {}).get("revision", "unversioned")
    validate_source_revision(revision, live=True)

    wanted_bases = {strip_version(gene): gene for gene in selected_genes}
    selected_set = set(selected_genes)
    rows = {}
    collisions = []
    sample_indices = None
    extract_digest = hashlib.sha256()
    response_metadata = {}

    print(f"[xena] streaming {matrix_url}", file=sys.stderr)
    with _http_client().get(matrix_url, stream=True, timeout=120) as response:
        response.raise_for_status()
        headers = getattr(response, "headers", {})
        response_metadata = {
            "url": matrix_url,
            "etag": headers.get("ETag"),
            "last_modified": headers.get("Last-Modified"),
            "content_length": headers.get("Content-Length"),
        }
        with gzip.GzipFile(fileobj=response.raw) as gz:
            for line_no, raw in enumerate(gz):
                line = raw.decode("utf-8", errors="replace").rstrip("\r\n")
                parts = line.split("\t")
                if line_no == 0:
                    extract_digest.update(raw)
                    if len(parts) < 2:
                        raise ValueError("Xena matrix header contains no sample columns")
                    if len(set(parts)) != len(parts):
                        duplicates = sorted({value for value in parts if parts.count(value) > 1})
                        raise ValueError(
                            f"Xena matrix header contains duplicate sample IDs: {duplicates[:5]}"
                        )
                    sample_to_idx = {sample: i for i, sample in enumerate(parts)}
                    missing_samples = [s for s in sample_ids if s not in sample_to_idx]
                    if missing_samples:
                        raise ValueError(f"{len(missing_samples)} sampled GTEx IDs missing from matrix")
                    sample_indices = [sample_to_idx[s] for s in sample_ids]
                    continue

                gene_id = parts[0]
                if not gene_id:
                    raise ValueError(f"Xena matrix row {line_no + 1} has a blank gene ID")
                target = gene_id if gene_id in selected_set else wanted_bases.get(strip_version(gene_id))
                if target is None:
                    continue
                if target in rows:
                    # A second raw row maps to the same selected gene after
                    # version-stripping (duplicate/ambiguous annotation row):
                    collisions.append(gene_id)
                    extract_digest.update(raw)
                    continue
                if sample_indices and max(sample_indices) >= len(parts):
                    raise ValueError(
                        f"Xena matrix row {line_no + 1} is shorter than its header"
                    )
                try:
                    vals = np.array([float(parts[i]) for i in sample_indices], dtype=float)
                    rows[target] = xena_log2_tpm001_to_model_scale(vals)
                except (TypeError, ValueError) as exc:
                    raise ValueError(
                        f"Xena matrix has invalid expression values for {gene_id!r}"
                    ) from exc
                extract_digest.update(raw)
                if len(rows) % 250 == 0:
                    print(f"[gtex] extracted {len(rows)}/{len(selected_genes)} genes",
                          file=sys.stderr)

    if not rows:
        raise ValueError(
            f"No selected genes matched any row in {matrix_url}; the download may be "
            "empty, truncated, or in an unexpected format."
        )
    if collisions:
        raise ValueError(
            f"Xena matrix has {len(collisions)} gene rows colliding with an "
            "already-matched selected gene after version stripping: "
            f"{collisions[:5]}"
        )

    matrix = pd.DataFrame(rows, index=sample_ids)
    matrix = matrix[[gene for gene in selected_genes if gene in matrix.columns]]
    # A source gene absent from the Xena matrix is a missing model gene, not a
    # matched gene full of NaNs.  Omitting all-missing columns lets the shared
    # scorer apply its documented training-mean imputation contract while still
    # rejecting isolated invalid cells in genes that really are present.
    all_missing = matrix.columns[matrix.isna().all(axis=0)].tolist()
    if all_missing:
        print(
            f"[xena] WARNING: dropping {len(all_missing)} all-missing source genes so they "
            f"are handled as missing model genes: {all_missing[:5]}",
            file=sys.stderr,
        )
        matrix = matrix.drop(columns=all_missing)
    if matrix.shape[1] == 0:
        raise ValueError("No usable selected genes remained after Xena extraction")
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    write_dataframe_cache(
        cache_path,
        matrix,
        fingerprint=fingerprint,
        fingerprint_inputs=inputs,
        cache_kind="xena_selected_gene_matrix",
        extra_metadata={
            "source_response": response_metadata,
            "selected_source_lines_sha256": extract_digest.hexdigest(),
            "matched_genes": int(matrix.shape[1]),
            "missing_genes": int(len(selected_genes) - matrix.shape[1]),
            "collisions": int(len(collisions)),
            "dropped_all_missing_genes": all_missing,
        },
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
                 by_site: pd.DataFrame, sweep: pd.DataFrame,
                 model_description: str) -> None:
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
- Model: {model_description}

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
    atomic_write_text(path, text)


def main() -> int:
    started_at = utc_now()
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-per-site", type=int, default=20)
    parser.add_argument("--min-site-n", type=int, default=20)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument(
        "--offline", "--cache-only", dest="offline", action="store_true",
        help="forbid network access and require semantically valid local caches",
    )
    parser.add_argument(
        "--sample-manifest",
        help="locked cohort CSV; use these exact sample IDs and do not resample phenotype data",
    )
    parser.add_argument(
        "--allow-donor-overlap",
        action="store_true",
        help="allow multiple tissues from one donor when constructing a new sampled panel",
    )
    parser.add_argument(
        "--source-revision",
        default="unversioned",
        help="provider snapshot/revision identifier included in cache provenance",
    )
    parser.add_argument("--out-dir", default=str(ROOT / "external-validation" / "gtex_xena"))
    parser.add_argument("--weights", default=str(ROOT / "deployable_lr_weights.npz"))
    parser.add_argument("--max-invalid-cell-fraction", type=float, default=0.0,
                        help=("maximum allowed missing, non-numeric, NaN, or infinite cells "
                              "among matched model genes before failing (default 0)"))
    parser.add_argument("--allow-invalid-values", action="store_true",
                        help=("warn instead of failing when matched model-gene cells are "
                              "missing, non-numeric, NaN, or infinite"))
    args = parser.parse_args()
    try:
        validate_threshold(args.threshold, "--threshold")
        validate_threshold(args.max_invalid_cell_fraction, "--max-invalid-cell-fraction")
        _validate_sampling_count(args.n_per_site, "--n-per-site")
        _validate_sampling_count(args.min_site_n, "--min-site-n")
        if args.offline and args.refresh:
            raise ValueError("--offline/--cache-only cannot be combined with --refresh")
        args.source_revision = validate_source_revision(
            args.source_revision, live=not args.offline
        )
    except ValueError as exc:
        parser.error(str(exc))

    out_dir = Path(args.out_dir).resolve(strict=False)
    phenotype_path = out_dir / "TcgaTargetGTEX_phenotype.csv"
    sampled_path = out_dir / "sampled_gtex_manifest.csv"
    matrix_path = out_dir / "gtex_selected_genes_model_scale.parquet"
    locked_manifest_path = Path(args.sample_manifest).resolve() if args.sample_manifest else None
    final_outputs = {
        "sampled_manifest": sampled_path,
        "predictions": out_dir / "gtex_predictions.csv",
        "summary": out_dir / "gtex_summary.csv",
        "threshold_sweep": out_dir / "gtex_threshold_sweep.csv",
        "per_site_summary": out_dir / "gtex_per_site_summary.csv",
        "report": out_dir / "GTEX_NORMAL_VALIDATION.md",
        "run_manifest": out_dir / "run_manifest.json",
    }
    protected_inputs = {"model": Path(args.weights)}
    if locked_manifest_path:
        protected_inputs["locked_sample_manifest"] = locked_manifest_path
    validate_managed_paths(
        protected_inputs=protected_inputs,
        managed_files={
            **final_outputs,
            "phenotype_cache": phenotype_path,
            "phenotype_cache_metadata": cache_meta_path(phenotype_path),
            "expression_matrix_cache": matrix_path,
            "expression_matrix_cache_metadata": cache_meta_path(matrix_path),
        },
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    if locked_manifest_path:
        sampled = load_locked_sample_manifest(locked_manifest_path, "GTEX")
        cohort_mode = "locked_manifest"
    else:
        phenotype = load_or_download_phenotype(
            phenotype_path, args.refresh, args.source_revision, args.offline
        )
        sampled = choose_gtex_samples(
            phenotype,
            args.n_per_site,
            args.min_site_n,
            args.seed,
            allow_donor_overlap=args.allow_donor_overlap,
        )
        cohort_mode = "sampled_from_phenotype"
    print(f"[gtex] sampled {len(sampled)} GTEx normals across "
          f"{sampled['_primary_site'].nunique()} primary sites", file=sys.stderr)

    model = load_lr_model(args.weights)
    model_sha256 = sha256_file(args.weights)
    selected_genes = model["genes"].astype(str).tolist()
    source_identity = {
        "dataset_id": GTEX_DATASET_ID,
        "url": GTEX_TPM_URL,
        "revision": args.source_revision,
    }
    matrix = extract_matrix_from_xena(sampled["sample"].tolist(), selected_genes,
                                      matrix_path, args.refresh,
                                      model_sha256=model_sha256,
                                      source_identity=source_identity,
                                      offline=args.offline)
    n_matched = int(matrix.shape[1])
    if n_matched < len(selected_genes):
        print(f"[gtex] WARNING: matched {n_matched}/{len(selected_genes)} genes",
              file=sys.stderr)

    probabilities, alignment_report = predict_proba(
        model,
        matrix,
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
    scored = scored_dataframe(matrix.index, probabilities, args.threshold)
    predictions = sampled.merge(
        scored, left_on="sample", right_on="sample", how="left", validate="one_to_one"
    )
    require_complete_merge(predictions, len(sampled), "tumor_probability", "[gtex]")

    summary = summarize(predictions, args.threshold)

    sweep = threshold_sweep(predictions)

    by_site = predictions.groupby("_primary_site").agg(
        n=("sample", "size"),
        tumor_calls=("call", lambda s: int((s == "tumor").sum())),
        false_positive_rate=("call", lambda s: float((s == "tumor").mean())),
        median_tumor_probability=("tumor_probability", "median"),
        max_tumor_probability=("tumor_probability", "max"),
    ).reset_index().sort_values(["false_positive_rate", "max_tumor_probability"],
                                ascending=False)
    model_description = f"`{Path(args.weights).name}` (SHA256 `{model_sha256[:12]}...`)"

    cohort_details = group_audit(sampled, "sample", "donor_id")
    cohort_details["cohort_mode"] = cohort_mode
    matrix_metadata = read_cache_metadata(matrix_path) or {}
    input_paths = {"expression_matrix_cache": matrix_path}
    if locked_manifest_path:
        input_paths["locked_sample_manifest"] = locked_manifest_path
    else:
        input_paths["phenotype_cache"] = phenotype_path
    with staged_output_directory(out_dir) as stage_dir:
        staged_outputs = {
            name: stage_dir / path.name for name, path in final_outputs.items()
        }
        atomic_write_csv(staged_outputs["sampled_manifest"], sampled)
        atomic_write_csv(staged_outputs["predictions"], predictions)
        atomic_write_csv(staged_outputs["summary"], pd.DataFrame([summary]))
        atomic_write_csv(staged_outputs["threshold_sweep"], sweep)
        atomic_write_csv(staged_outputs["per_site_summary"], by_site)
        write_report(
            staged_outputs["report"], summary, predictions, by_site, sweep,
            model_description,
        )
        manifest_outputs = {
            name: path for name, path in staged_outputs.items()
            if name != "run_manifest"
        }
        manifest_final_outputs = {
            name: path for name, path in final_outputs.items()
            if name != "run_manifest"
        }
        write_run_manifest(
            staged_outputs["run_manifest"],
            root=ROOT,
            run_kind="gtex_xena_normal_validation",
            started_at_utc=started_at,
            argv=sys.argv,
            parameters={**vars(args), "cohort_mode": cohort_mode},
            model_path=args.weights,
            sources={"phenotype": PHENOTYPE_URL, "expression": source_identity},
            inputs=input_paths,
            outputs=manifest_outputs,
            output_display_paths=manifest_final_outputs,
            alignment=alignment_report,
            cohort_audit=cohort_details,
            cache_details={"expression_matrix": matrix_metadata},
            source_code={
                "validator": Path(__file__).resolve(),
                "score_core": ROOT / "tcga_rnaseq" / "score.py",
                "alignment_core": ROOT / "tcga_rnaseq" / "align.py",
            },
        )
        publish_staged_files(staged_outputs, final_outputs)

    print(f"[gtex] summary -> {final_outputs['summary']}", file=sys.stderr)
    print(f"[gtex] per-site -> {final_outputs['per_site_summary']}", file=sys.stderr)
    print(f"[gtex] predictions -> {final_outputs['predictions']}", file=sys.stderr)
    print(f"[gtex] report -> {final_outputs['report']}", file=sys.stderr)
    print(f"[gtex] provenance -> {final_outputs['run_manifest']}", file=sys.stderr)
    print(pd.DataFrame([summary]).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
