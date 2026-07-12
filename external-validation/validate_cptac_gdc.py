#!/usr/bin/env python3
"""External smoke validation on CPTAC STAR-Counts files from the GDC.

This is intentionally a lightweight external check: CPTAC is not TCGA, but the
files are harmonized by the same GDC STAR-Counts pipeline, so input scale and
gene identifiers match the deployed model. The script downloads only the
selected 2,000 genes per file and caches those extracted vectors.
"""
from __future__ import annotations

import argparse
import hashlib
import math
import sys
import time
import uuid
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
try:
    import requests
except ImportError:  # Keep pure-logic imports usable in the lightweight environment.
    requests = None


class SourceIntegrityError(ValueError):
    """Downloaded source bytes did not match their provider checksum."""


REQUEST_EXCEPTIONS = (requests.RequestException,) if requests is not None else ()
RETRYABLE_EXCEPTIONS = REQUEST_EXCEPTIONS + (SourceIntegrityError,)

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from tcga_rnaseq import load_lr_model, predict_proba, validate_threshold  # noqa: E402
from tcga_rnaseq import metrics as M  # noqa: E402
from tcga_rnaseq.align import strip_version  # noqa: E402
from provenance import (  # noqa: E402
    atomic_write_csv,
    atomic_write_text,
    cache_meta_path,
    contained_cache_path,
    group_audit,
    load_csv_cache,
    load_dataframe_cache,
    load_series_cache,
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
    write_series_cache,
)

GDC_FILES_ENDPOINT = "https://api.gdc.cancer.gov/files"
GDC_DATA_ENDPOINT = "https://api.gdc.cancer.gov/data"
LABEL_MAP = {"Solid Tissue Normal": 0, "Primary Tumor": 1}
GDC_MANIFEST_PARSER_VERSION = "gdc-cptac-manifest-v3"
GDC_STAR_COUNTS_PARSER_VERSION = "gdc-star-counts-selected-genes-v4"
GDC_TPM_TRANSFORM_VERSION = "tpm_unstranded_to_log2-tpm-plus-1-v1"


def _http_client():
    if requests is None:
        raise RuntimeError(
            "live external validation requires the optional 'requests' package"
        )
    return requests


def _iter_hashed_response_lines(response, digest):
    """Yield decoded lines while hashing the exact downloaded response bytes."""
    pending = b""
    for chunk in response.iter_content(chunk_size=1024 * 1024):
        if not chunk:
            continue
        digest.update(chunk)
        pending += chunk
        while b"\n" in pending:
            line, pending = pending.split(b"\n", 1)
            yield line.rstrip(b"\r").decode("utf-8", errors="replace")
    if pending:
        yield pending.rstrip(b"\r").decode("utf-8", errors="replace")


def _new_md5_digest():
    """Create a non-security MD5 digest, including on FIPS-aware Python builds."""
    try:
        return hashlib.md5(usedforsecurity=False)
    except TypeError:  # Older Python/OpenSSL builds do not expose this keyword.
        return hashlib.md5()


def query_cptac_manifest() -> pd.DataFrame:
    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.program.name",
                                       "value": ["CPTAC"]}},
            {"op": "in", "content": {"field": "files.data_type",
                                       "value": ["Gene Expression Quantification"]}},
            {"op": "in", "content": {"field": "files.analysis.workflow_type",
                                       "value": ["STAR - Counts"]}},
            {"op": "in", "content": {"field": "files.data_format",
                                       "value": ["TSV"]}},
            {"op": "in", "content": {"field": "files.experimental_strategy",
                                       "value": ["RNA-Seq"]}},
            {"op": "in", "content": {"field": "cases.samples.sample_type",
                                       "value": list(LABEL_MAP)}},
        ],
    }
    fields = ",".join([
        "file_id",
        "file_name",
        "file_size",
        "md5sum",
        "data_format",
        "experimental_strategy",
        "created_datetime",
        "updated_datetime",
        "analysis.workflow_type",
        "cases.case_id",
        "cases.submitter_id",
        "cases.project.project_id",
        "cases.samples.sample_type",
        "cases.samples.submitter_id",
    ])
    payload = {"filters": filters, "fields": fields, "format": "json", "size": 5000}
    response = _http_client().post(GDC_FILES_ENDPOINT, json=payload, timeout=60)
    response.raise_for_status()

    rows = []
    for hit in response.json()["data"]["hits"]:
        sample_types = set()
        sample_ids = []
        projects = set()
        case_ids = []
        for case in hit.get("cases", []):
            if case.get("submitter_id"):
                case_ids.append(case["submitter_id"])
            if case.get("project"):
                projects.add(case["project"].get("project_id"))
            for sample in case.get("samples", []):
                if sample.get("sample_type"):
                    sample_types.add(sample["sample_type"])
                if sample.get("submitter_id"):
                    sample_ids.append(sample["submitter_id"])
        rows.append({
            "file_id": hit["file_id"],
            "file_name": hit.get("file_name"),
            "file_size": hit.get("file_size"),
            "md5sum": hit.get("md5sum"),
            "data_format": hit.get("data_format"),
            "experimental_strategy": hit.get("experimental_strategy"),
            "workflow_type": (hit.get("analysis") or {}).get("workflow_type"),
            "created_datetime": hit.get("created_datetime"),
            "updated_datetime": hit.get("updated_datetime"),
            "project": ";".join(sorted(projects)),
            "case_submitter_id": ";".join(sorted(set(case_ids))),
            "sample_submitter_id": ";".join(sorted(set(sample_ids))),
            "sample_type": ";".join(sorted(sample_types)),
            "n_sample_types": len(sample_types),
        })
    manifest = pd.DataFrame(rows).drop_duplicates("file_id")
    manifest = manifest[manifest["n_sample_types"] == 1].copy()
    manifest = manifest[manifest["sample_type"].isin(LABEL_MAP)].copy()
    return manifest


def _manifest_cache_inputs(source_revision: str) -> dict:
    return {
        "cache_schema": "gdc-cptac-manifest-csv-v2",
        "parser_version": GDC_MANIFEST_PARSER_VERSION,
        "source": {
            "endpoint": GDC_FILES_ENDPOINT,
            "program": "CPTAC",
            "workflow": "STAR - Counts",
            "revision": source_revision,
        },
    }


def load_or_query_manifest(
    path: Path,
    refresh: bool,
    source_revision: str = "unversioned",
    offline: bool = False,
) -> pd.DataFrame:
    inputs = _manifest_cache_inputs(source_revision)
    fingerprint = semantic_fingerprint(inputs)
    required = [
        "file_id", "file_name", "file_size", "md5sum", "project",
        "case_submitter_id", "sample_submitter_id", "sample_type", "n_sample_types",
    ]
    if not refresh:
        cached = load_csv_cache(path, fingerprint=fingerprint, required_columns=required)
        if cached is not None:
            return cached
    if offline:
        raise ValueError(
            f"offline/cache-only mode requires a valid GDC manifest cache: {path}"
        )
    validate_source_revision(source_revision, live=True)
    manifest = query_cptac_manifest()
    path.parent.mkdir(parents=True, exist_ok=True)
    write_csv_cache(
        path,
        manifest,
        fingerprint=fingerprint,
        fingerprint_inputs=inputs,
        cache_kind="gdc_cptac_manifest",
    )
    return manifest


def refresh_locked_sample_metadata(
    sampled: pd.DataFrame,
    *,
    expected_project: str,
) -> pd.DataFrame:
    """Refresh provider metadata for fixed file IDs without resampling the cohort."""
    provider = query_cptac_manifest()
    if "file_id" not in provider.columns:
        raise ValueError("refreshed GDC CPTAC manifest is missing file_id")
    provider["file_id"] = validate_identifier_column(
        provider, "file_id", "refreshed GDC CPTAC manifest"
    ).map(validate_gdc_file_id)
    wanted = sampled["file_id"].tolist()
    provider = provider[
        provider["project"].eq(expected_project) & provider["file_id"].isin(wanted)
    ].copy()
    provider = _validate_cptac_sample_frame(
        provider,
        expected_project=expected_project,
        require_provider_md5=True,
        context="refreshed GDC CPTAC manifest",
    )
    provider_by_id = provider.set_index("file_id", verify_integrity=True)
    missing = [file_id for file_id in wanted if file_id not in provider_by_id.index]
    if missing:
        raise ValueError(
            f"GDC metadata refresh did not return {len(missing)} locked file IDs: {missing[:5]}"
        )
    refreshed = provider_by_id.loc[wanted].reset_index()
    for column in ("project", "case_submitter_id", "sample_submitter_id", "sample_type"):
        old = sampled[column].astype(str).tolist()
        new = refreshed[column].astype(str).tolist()
        if old != new:
            raise ValueError(
                f"GDC metadata refresh changed locked cohort field {column!r}; "
                "refusing to rebind the cohort"
            )
    for column in ("file_name", "file_size", "md5sum"):
        if column not in sampled.columns:
            continue
        old = sampled[column].astype(str).str.lower().tolist()
        new = refreshed[column].astype(str).str.lower().tolist()
        if old != new:
            raise ValueError(
                f"GDC metadata refresh changed locked provider identity field {column!r}; "
                "refusing to rebind the cohort"
            )
    return _validate_cptac_sample_frame(
        refreshed,
        expected_project=expected_project,
        require_provider_md5=True,
        context="refreshed locked CPTAC sample manifest",
    )


def _normalized_nonempty(frame: pd.DataFrame, column: str, context: str) -> pd.Series:
    if column not in frame.columns:
        raise ValueError(f"{context} must contain {column!r}")
    if frame[column].isna().any():
        raise ValueError(f"{context} contains missing {column} values")
    raw = frame[column].astype(str)
    values = raw.str.strip()
    if values.eq("").any():
        raise ValueError(f"{context} contains blank {column} values")
    if not raw.eq(values).all():
        raise ValueError(f"{context} contains padded {column} values")
    return values


def validate_gdc_file_id(value: object) -> str:
    text = str(value).strip().lower()
    try:
        parsed = uuid.UUID(text)
    except (ValueError, AttributeError) as exc:
        raise ValueError(f"invalid canonical GDC file UUID: {value!r}") from exc
    canonical = str(parsed)
    if text != canonical:
        raise ValueError(f"invalid canonical GDC file UUID: {value!r}")
    return canonical


def _validate_provider_md5(values: pd.Series, context: str) -> pd.Series:
    normalized = values.astype(str).str.strip().str.lower()
    valid = normalized.str.fullmatch(r"[0-9a-f]{32}", na=False)
    if not valid.all():
        bad = normalized[~valid].tolist()[:5]
        raise ValueError(f"{context} contains invalid provider md5sum values: {bad}")
    return normalized


def _validate_cptac_sample_frame(
    sampled: pd.DataFrame,
    *,
    expected_project: str | None,
    require_provider_md5: bool,
    context: str,
) -> pd.DataFrame:
    sampled = sampled.copy()
    sampled["file_id"] = validate_identifier_column(
        sampled, "file_id", context
    )
    sampled["file_id"] = sampled["file_id"].map(validate_gdc_file_id)
    sampled["sample_submitter_id"] = validate_identifier_column(
        sampled, "sample_submitter_id", context
    )
    sampled["case_submitter_id"] = _normalized_nonempty(
        sampled, "case_submitter_id", context
    )
    sampled["project"] = _normalized_nonempty(
        sampled, "project", context
    )
    sampled["sample_type"] = _normalized_nonempty(
        sampled, "sample_type", context
    )
    if expected_project is not None and not sampled["project"].eq(expected_project).all():
        found = sorted(set(sampled["project"]))
        raise ValueError(
            f"{context} must contain only project={expected_project!r}; found {found}"
        )
    if require_provider_md5 and "md5sum" not in sampled.columns:
        raise ValueError(
            f"{context} has no provider md5sum; historical manifests may be used only "
            "with --offline and valid caches, or refreshed from GDC metadata"
        )
    if (
        require_provider_md5
        and sampled["sample_submitter_id"].str.contains(";", regex=False).any()
    ):
        raise ValueError(f"{context} contains ambiguous multi-biospecimen identifiers")
    if require_provider_md5 and sampled["case_submitter_id"].str.contains(";", regex=False).any():
        raise ValueError(f"{context} contains ambiguous multi-case identifiers")
    if not sampled["sample_type"].isin(LABEL_MAP).all():
        bad = sorted(set(sampled.loc[~sampled["sample_type"].isin(LABEL_MAP), "sample_type"]))
        raise ValueError(f"{context} has unsupported sample_type values: {bad}")
    expected_labels = sampled["sample_type"].map(LABEL_MAP).astype(int)
    if "label" in sampled.columns:
        labels = pd.to_numeric(sampled["label"], errors="coerce")
        if labels.isna().any() or not np.array_equal(labels.astype(int), expected_labels):
            raise ValueError(f"{context} labels do not agree with sample_type")
    sampled["label"] = expected_labels
    if sampled["label"].nunique() != 2:
        raise ValueError(f"{context} must contain both tumor and normal")
    if "n_sample_types" in sampled.columns:
        counts = pd.to_numeric(sampled["n_sample_types"], errors="coerce")
        if counts.isna().any() or not counts.eq(1).all():
            raise ValueError(f"{context} requires n_sample_types=1 for every row")
    elif require_provider_md5:
        raise ValueError(f"{context} must contain n_sample_types")
    for column in ("file_name", "file_size"):
        if require_provider_md5 and column not in sampled.columns:
            raise ValueError(f"{context} must contain {column!r}")
    if "file_name" in sampled.columns:
        sampled["file_name"] = _normalized_nonempty(sampled, "file_name", context)
    if "file_size" in sampled.columns:
        sizes = pd.to_numeric(sampled["file_size"], errors="coerce")
        if sizes.isna().any() or not np.equal(sizes, np.floor(sizes)).all() or not sizes.gt(0).all():
            raise ValueError(f"{context} contains invalid file_size values")
        sampled["file_size"] = sizes.astype(np.int64)
    if "md5sum" in sampled.columns:
        sampled["md5sum"] = _validate_provider_md5(sampled["md5sum"], context)
    return sampled


def validate_locked_sample_manifest(
    path: Path,
    *,
    expected_project: str | None = None,
    require_provider_md5: bool = False,
) -> pd.DataFrame:
    sampled = pd.read_csv(path, dtype={"file_id": str}, keep_default_na=False)
    return _validate_cptac_sample_frame(
        sampled,
        expected_project=expected_project,
        require_provider_md5=require_provider_md5,
        context="locked CPTAC sample manifest",
    )


def choose_files(
    manifest: pd.DataFrame,
    project: str,
    n_per_class: int,
    seed: int,
    allow_case_overlap: bool = False,
) -> pd.DataFrame:
    if n_per_class <= 0:
        raise ValueError("n_per_class must be a positive integer")
    manifest = manifest.copy()
    manifest["file_id"] = validate_identifier_column(
        manifest, "file_id", "CPTAC manifest"
    ).map(validate_gdc_file_id)
    _normalized_nonempty(manifest, "sample_submitter_id", "CPTAC manifest")
    _normalized_nonempty(manifest, "case_submitter_id", "CPTAC manifest")
    _normalized_nonempty(manifest, "project", "CPTAC manifest")
    _normalized_nonempty(manifest, "sample_type", "CPTAC manifest")
    for column in ("file_name", "file_size", "md5sum", "n_sample_types"):
        if column not in manifest.columns:
            raise ValueError(f"CPTAC manifest must contain {column!r}")
    manifest["md5sum"] = _validate_provider_md5(manifest["md5sum"], "CPTAC manifest")
    sizes = pd.to_numeric(manifest["file_size"], errors="coerce")
    if sizes.isna().any() or not sizes.gt(0).all():
        raise ValueError("CPTAC manifest contains invalid file_size values")
    counts = pd.to_numeric(manifest["n_sample_types"], errors="coerce")
    if counts.isna().any() or not counts.eq(1).all():
        raise ValueError("CPTAC manifest requires n_sample_types=1 for every row")
    df = manifest[manifest["project"] == project].copy()
    if df.empty:
        raise ValueError(f"No files found for project {project}")

    sampled = []
    used_cases: set[str] = set()
    for sample_type in LABEL_MAP:
        part = df[df["sample_type"] == sample_type].copy()
        # A GDC query can expose multiple STAR-Counts files for one biospecimen.
        # Pick one deterministically before random case-level sampling.
        part = part.sort_values("file_id").drop_duplicates("sample_submitter_id", keep="first")
        part = part.drop_duplicates("case_submitter_id", keep="first")
        if not allow_case_overlap:
            part = part[~part["case_submitter_id"].astype(str).isin(used_cases)]
        if part.empty:
            raise ValueError(f"No {sample_type} files found for {project}")
        n = min(n_per_class, len(part))
        sampled.append(part.sample(n=n, random_state=seed).copy())
        used_cases.update(sampled[-1]["case_submitter_id"].astype(str))
    out = pd.concat(sampled, ignore_index=True)
    out["label"] = out["sample_type"].map(LABEL_MAP)
    return out.sample(frac=1, random_state=seed).reset_index(drop=True)


def _gene_cache_inputs(
    file_id: str,
    selected_genes: list[str],
    model_sha256: str,
    source_identity: dict | None,
) -> dict:
    file_id = validate_gdc_file_id(validate_unique_strings([file_id], "GDC file ID")[0])
    selected_genes = validate_unique_strings(selected_genes, "selected model genes")
    bases = [strip_version(gene) for gene in selected_genes]
    if len(set(bases)) != len(bases):
        raise ValueError("selected model genes collide after removing Ensembl version suffixes")
    return {
        "cache_schema": "gdc-selected-gene-series-parquet-v4",
        "parser_version": GDC_STAR_COUNTS_PARSER_VERSION,
        "transform_version": GDC_TPM_TRANSFORM_VERSION,
        "source": source_identity or {
            "endpoint": GDC_DATA_ENDPOINT,
            "file_id": file_id,
            "revision": "unversioned",
        },
        "model_sha256": model_sha256,
        "selected_genes": [str(gene) for gene in selected_genes],
    }


def gene_cache_path(cache_dir: Path, file_id: str) -> Path:
    file_id = validate_gdc_file_id(file_id)
    return contained_cache_path(
        cache_dir,
        file_id,
        namespace="gdc-file-id",
        suffix=".parquet",
    )


def extract_selected_genes(
    file_id: str,
    selected_genes: list[str],
    cache_dir: Path,
    retries: int = 4,
    refresh: bool = False,
    model_sha256: str = "unspecified",
    source_identity: dict | None = None,
    offline: bool = False,
) -> pd.Series:
    if retries <= 0:
        raise ValueError("retries must be a positive integer")
    file_id = validate_gdc_file_id(file_id)
    cache_path = gene_cache_path(cache_dir, file_id)
    inputs = _gene_cache_inputs(file_id, selected_genes, model_sha256, source_identity)
    fingerprint = semantic_fingerprint(inputs)
    if not refresh:
        cached = load_series_cache(
            cache_path,
            fingerprint=fingerprint,
            expected_name=file_id,
            allowed_index=selected_genes,
        )
        if cached is not None:
            expected_order = [gene for gene in selected_genes if gene in set(cached.index.astype(str))]
            if cached.index.astype(str).tolist() != expected_order:
                cached = None
        if cached is not None:
            return cached
        if cache_path.exists():
            print(
                f"[cptac] cache for {file_id} failed fingerprint/content/index validation; "
                "re-downloading",
                file=sys.stderr,
            )
    if offline:
        raise ValueError(
            f"offline/cache-only mode requires a valid gene cache for GDC file {file_id}"
        )
    validate_source_revision(
        (source_identity or {}).get("revision", "unversioned"), live=True
    )

    wanted = set(selected_genes)
    wanted_by_base = {strip_version(gene): gene for gene in selected_genes}
    values = {}
    url = f"{GDC_DATA_ENDPOINT}/{file_id}"
    last_error = None
    selected_source_digest = hashlib.sha256()
    downloaded_md5 = None
    response_metadata = {}
    expected_md5 = str((source_identity or {}).get("md5sum", "")).strip().lower()
    if (
        len(expected_md5) != 32
        or any(character not in "0123456789abcdef" for character in expected_md5)
    ):
        raise ValueError(
            f"new GDC downloads require a valid provider MD5 for file {file_id}; "
            f"found {expected_md5!r}"
        )
    cache_dir.mkdir(parents=True, exist_ok=True)
    for attempt in range(1, retries + 1):
        values = {}
        selected_source_digest = hashlib.sha256()
        downloaded_md5 = None
        try:
            with _http_client().get(url, stream=True, timeout=120) as response:
                response.raise_for_status()
                headers = getattr(response, "headers", {})
                response_metadata = {
                    "url": url,
                    "etag": headers.get("ETag"),
                    "last_modified": headers.get("Last-Modified"),
                    "content_length": headers.get("Content-Length"),
                }
                header = None
                gene_idx = tpm_idx = None
                download_digest = _new_md5_digest()
                line_iterator = _iter_hashed_response_lines(response, download_digest)
                for raw_line in line_iterator:
                    if isinstance(raw_line, bytes):
                        raw_line = raw_line.decode("utf-8", errors="replace")
                    if not raw_line or raw_line.startswith("#"):
                        continue
                    parts = raw_line.split("\t")
                    if header is None:
                        header = parts
                        gene_idx = header.index("gene_id")
                        tpm_idx = header.index("tpm_unstranded")
                        selected_source_digest.update(raw_line.encode("utf-8"))
                        continue

                    gene_id = parts[gene_idx]
                    target = gene_id if gene_id in wanted else wanted_by_base.get(strip_version(gene_id))
                    if target is None:
                        continue
                    if target in values:
                        raise ValueError(
                            f"GDC file {file_id} has multiple rows mapping to selected gene "
                            f"{target!r} after Ensembl version stripping (latest row {gene_id!r})"
                        )
                    try:
                        tpm = float(parts[tpm_idx])
                    except ValueError:
                        raise ValueError(
                            f"GDC file {file_id} has non-numeric tpm_unstranded for {gene_id!r}"
                        )
                    if not math.isfinite(tpm) or tpm < 0:
                        raise ValueError(
                            f"GDC file {file_id} has invalid tpm_unstranded={tpm!r} for "
                            f"{gene_id!r}; expected a finite non-negative value"
                        )
                    values[target] = math.log2(tpm + 1.0)
                    selected_source_digest.update(raw_line.encode("utf-8"))
                downloaded_md5 = download_digest.hexdigest().lower()
                if downloaded_md5 != expected_md5:
                    raise SourceIntegrityError(
                        f"GDC file {file_id} MD5 mismatch: expected {expected_md5}, "
                        f"downloaded {downloaded_md5}"
                    )
            break
        except RETRYABLE_EXCEPTIONS as exc:
            last_error = exc
            if attempt == retries:
                raise
            time.sleep(2 ** (attempt - 1))
            print(f"[cptac] retry {attempt}/{retries} for {file_id}: {exc}",
                  file=sys.stderr)
    if last_error is not None and not values:
        raise last_error
    if not values:
        raise ValueError(
            f"No selected genes matched GDC file {file_id}; the download may be "
            "empty, corrupt, or in an unexpected format."
        )
    ordered_genes = [gene for gene in selected_genes if gene in values]
    series = pd.Series([values[gene] for gene in ordered_genes], index=ordered_genes,
                       name=file_id, dtype=float)
    write_series_cache(
        cache_path,
        series,
        fingerprint=fingerprint,
        fingerprint_inputs=inputs,
        cache_kind="gdc_selected_gene_series",
        extra_metadata={
            "source_response": response_metadata,
            "provider_md5": expected_md5,
            "downloaded_md5": downloaded_md5,
            "selected_source_lines_sha256": selected_source_digest.hexdigest(),
            "matched_genes": int(len(series)),
            "missing_genes": int(len(selected_genes) - len(series)),
            "collisions": 0,
        },
    )
    return series


def _source_identity_from_row(row, source_revision: str) -> dict:
    identity = {
        "endpoint": GDC_DATA_ENDPOINT,
        "file_id": str(row.file_id),
        "revision": source_revision,
    }
    for field in ("file_name", "file_size", "md5sum", "updated_datetime"):
        if hasattr(row, field):
            value = getattr(row, field)
            if pd.notna(value):
                identity[field] = value.item() if isinstance(value, np.generic) else value
    return identity


def _expression_cache_inputs(
    files: pd.DataFrame,
    selected_genes: list[str],
    model_sha256: str,
    source_revision: str,
) -> dict:
    selected_genes = validate_unique_strings(selected_genes, "selected model genes")
    bases = [strip_version(gene) for gene in selected_genes]
    if len(set(bases)) != len(bases):
        raise ValueError("selected model genes collide after removing Ensembl version suffixes")
    identities = [
        _source_identity_from_row(row, source_revision)
        for row in files.itertuples(index=False)
    ]
    return {
        "cache_schema": "gdc-cptac-expression-matrix-parquet-v3",
        "parser_version": GDC_STAR_COUNTS_PARSER_VERSION,
        "transform_version": GDC_TPM_TRANSFORM_VERSION,
        "model_sha256": model_sha256,
        "selected_genes": [str(gene) for gene in selected_genes],
        "files": identities,
    }


def build_expression_matrix(
    files: pd.DataFrame,
    selected_genes: list[str],
    cache_dir: Path,
    workers: int,
    retries: int,
    *,
    matrix_cache_path: Path | None = None,
    refresh: bool = False,
    model_sha256: str = "unspecified",
    source_revision: str = "unversioned",
    offline: bool = False,
) -> pd.DataFrame:
    if workers <= 0:
        raise ValueError("workers must be a positive integer")
    if retries <= 0:
        raise ValueError("retries must be a positive integer")
    file_ids = [validate_gdc_file_id(value) for value in files["file_id"].tolist()]
    inputs = _expression_cache_inputs(files, selected_genes, model_sha256, source_revision)
    fingerprint = semantic_fingerprint(inputs)
    if matrix_cache_path is not None and not refresh:
        cached = load_dataframe_cache(
            matrix_cache_path,
            fingerprint=fingerprint,
            expected_index=file_ids,
            expected_columns=selected_genes,
            columns_may_be_subset=True,
            reject_all_missing_columns=True,
        )
        if cached is not None:
            return cached
    rows = {}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        future_to_file = {
            pool.submit(
                extract_selected_genes,
                row.file_id,
                selected_genes,
                cache_dir,
                retries,
                refresh,
                model_sha256,
                _source_identity_from_row(row, source_revision),
                offline,
            ): row.file_id
            for row in files.itertuples(index=False)
        }
        done = 0
        total = len(future_to_file)
        for future in as_completed(future_to_file):
            file_id = future_to_file[future]
            rows[file_id] = future.result()
            done += 1
            if done == total or done % 10 == 0:
                print(f"[cptac] extracted {done}/{total} files", file=sys.stderr)

    matrix = pd.DataFrame.from_dict(rows, orient="index")
    matrix = matrix.reindex(file_ids)
    matrix = matrix[[gene for gene in selected_genes if gene in matrix.columns]]
    all_missing = matrix.columns[matrix.isna().all(axis=0)].tolist()
    if all_missing:
        matrix = matrix.drop(columns=all_missing)
    if matrix.shape[1] == 0:
        raise ValueError("No usable selected genes remained in the CPTAC expression matrix")
    if matrix_cache_path is not None:
        write_dataframe_cache(
            matrix_cache_path,
            matrix,
            fingerprint=fingerprint,
            fingerprint_inputs=inputs,
            cache_kind="gdc_cptac_expression_matrix",
            extra_metadata={
                "matched_genes": int(matrix.shape[1]),
                "missing_genes": int(len(selected_genes) - matrix.shape[1]),
                "dropped_all_missing_genes": all_missing,
            },
        )
    return matrix


def summarize(y_true: np.ndarray, scores: np.ndarray, threshold: float) -> dict:
    metrics = M.classification_metrics(y_true, scores, threshold)
    return {
        "n": int(len(y_true)),
        "n_tumor": int(y_true.sum()),
        "n_normal": int((y_true == 0).sum()),
        "auc": metrics["auc"],
        "average_precision": M.average_precision(y_true, scores),
        "accuracy": metrics["accuracy"],
        "f1": metrics["f1"],
        "precision": metrics["precision"],
        "recall": metrics["recall"],
        "tn": metrics["tn"],
        "fp": metrics["fp"],
        "fn": metrics["fn"],
        "tp": metrics["tp"],
        "threshold": threshold,
    }


def threshold_sweep(y_true: np.ndarray, scores: np.ndarray) -> pd.DataFrame:
    youden_threshold = M.youden_threshold(y_true, scores)["threshold"]
    threshold_specs = [
        ("default_0.5", 0.5),
        ("youden_j", float(youden_threshold)),
        ("high_specificity_0.75", 0.75),
        ("high_specificity_0.9", 0.9),
        ("high_specificity_0.95", 0.95),
        ("very_high_specificity_0.99", 0.99),
    ]
    rows = []
    for name, threshold in threshold_specs:
        metrics = M.classification_metrics(y_true, scores, threshold)
        rows.append({
            "threshold_name": name,
            "threshold": threshold,
            "accuracy": metrics["accuracy"],
            "f1": metrics["f1"],
            "precision": metrics["precision"],
            "recall": metrics["recall"],
            "tn": metrics["tn"],
            "fp": metrics["fp"],
            "fn": metrics["fn"],
            "tp": metrics["tp"],
        })
    return pd.DataFrame(rows)


def write_report(path: Path, summary: dict, sampled: pd.DataFrame,
                 predictions: pd.DataFrame, n_manifest: int | None,
                 sweep: pd.DataFrame, model_description: str) -> None:
    by_type = predictions.groupby("sample_type").agg(
        n=("label", "size"),
        mean_tumor_probability=("tumor_probability", "mean"),
        tumor_calls=("call", lambda s: int((s == "tumor").sum())),
    )
    by_type_rows = "\n".join(
        f"| {idx} | {int(row.n)} | {row.mean_tumor_probability:.4f} | {int(row.tumor_calls)} |"
        for idx, row in by_type.iterrows()
    )
    sweep_rows = "\n".join(
        f"| {row.threshold_name} | {row.threshold:.6f} | {row.accuracy:.4f} | "
        f"{row.precision:.4f} | {row.recall:.4f} | {int(row.tn)} / {int(row.fp)} | "
        f"{int(row.fn)} / {int(row.tp)} |"
        for row in sweep.itertuples(index=False)
    )
    quantiles = predictions.groupby("sample_type")["tumor_probability"].quantile(
        [0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0]
    ).unstack()
    quantile_rows = "\n".join(
        f"| {idx} | " + " | ".join(f"{val:.4f}" for val in row) + " |"
        for idx, row in quantiles.iterrows()
    )
    manifest_description = str(n_manifest) if n_manifest is not None else "locked cohort manifest"
    text = f"""# CPTAC external smoke validation

## Data source

- Source: NCI GDC API, CPTAC program
- Files queried: {manifest_description}
- Project scored: {sampled['project'].iloc[0]}
- Sampled files: {len(sampled)} ({summary['n_tumor']} primary tumor, {summary['n_normal']} solid tissue normal)
- Workflow: GDC STAR-Counts, `tpm_unstranded` converted to log2(TPM+1)
- Model: {model_description}

## Result at threshold {summary['threshold']}

| Metric | Value |
|---|---:|
| AUC | {summary['auc']:.4f} |
| Average precision | {summary['average_precision']:.4f} |
| Accuracy | {summary['accuracy']:.4f} |
| F1 | {summary['f1']:.4f} |
| Precision | {summary['precision']:.4f} |
| Recall | {summary['recall']:.4f} |
| True normal / false tumor | {summary['tn']} / {summary['fp']} |
| False normal / true tumor | {summary['fn']} / {summary['tp']} |

## Probability by label

| Sample type | n | Mean tumor probability | Tumor calls |
|---|---:|---:|---:|
{by_type_rows}

## Threshold sensitivity

| Threshold | Cutoff | Accuracy | Precision | Recall | TN / FP | FN / TP |
|---|---:|---:|---:|---:|---:|---:|
{sweep_rows}

## Score quantiles

| Sample type | min | p10 | p25 | median | p75 | p90 | max |
|---|---:|---:|---:|---:|---:|---:|---:|
{quantile_rows}

## Interpretation

This is an external smoke test, not a replacement for a full independent
benchmark. It is outside TCGA but still uses GDC harmonized STAR-Counts, so it
tests cohort transfer more than cross-platform transfer. The stricter remaining
gap is non-GDC RNA-seq, where normalization and gene annotation differences can
dominate.
"""
    atomic_write_text(path, text)


def main() -> int:
    started_at = utc_now()
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", default="CPTAC-3")
    parser.add_argument("--n-per-class", type=int, default=50)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--retries", type=int, default=4)
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--refresh-manifest", action="store_true")
    parser.add_argument(
        "--refresh-expression-cache",
        action="store_true",
        help="ignore full-matrix and per-file expression caches and download fresh data",
    )
    parser.add_argument(
        "--offline", "--cache-only", dest="offline", action="store_true",
        help="forbid network access and require semantically valid local caches",
    )
    parser.add_argument(
        "--sample-manifest",
        help="locked sampled cohort CSV; use these exact file IDs and do not query/resample",
    )
    parser.add_argument(
        "--allow-case-overlap",
        action="store_true",
        help="allow more than one sampled row from the same case in a newly sampled cohort",
    )
    parser.add_argument(
        "--source-revision",
        default="unversioned",
        help="GDC data release/snapshot identifier included in cache provenance",
    )
    parser.add_argument("--out-dir", default=str(ROOT / "external-validation" / "cptac_gdc"))
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
        if args.n_per_class <= 0:
            raise ValueError("--n-per-class must be a positive integer")
        if args.workers <= 0:
            raise ValueError("--workers must be a positive integer")
        if args.retries <= 0:
            raise ValueError("--retries must be a positive integer")
        if not args.project or args.project != args.project.strip():
            raise ValueError("--project must be non-empty and unpadded")
        if args.offline and (args.refresh_manifest or args.refresh_expression_cache):
            raise ValueError(
                "--offline/--cache-only cannot be combined with refresh options"
            )
        args.source_revision = validate_source_revision(
            args.source_revision, live=not args.offline
        )
    except ValueError as exc:
        parser.error(str(exc))

    out_dir = Path(args.out_dir).resolve(strict=False)
    manifest_path = out_dir / "cptac_gdc_manifest.csv"
    cache_dir = out_dir / "gene_cache"
    sampled_path = out_dir / "sampled_manifest.csv"
    locked_manifest_path = Path(args.sample_manifest).resolve() if args.sample_manifest else None
    matrix_path = out_dir / "expression_selected_genes.parquet"
    final_outputs = {
        "sampled_manifest": sampled_path,
        "predictions": out_dir / "cptac_predictions.csv",
        "summary": out_dir / "cptac_summary.csv",
        "threshold_sweep": out_dir / "cptac_threshold_sweep.csv",
        "report": out_dir / "CPTAC_EXTERNAL_VALIDATION.md",
        "run_manifest": out_dir / "run_manifest.json",
    }
    protected_inputs = {"model": Path(args.weights)}
    if locked_manifest_path:
        protected_inputs["locked_sample_manifest"] = locked_manifest_path
    fixed_managed_files = {
        **final_outputs,
        "gdc_manifest_cache": manifest_path,
        "gdc_manifest_cache_metadata": cache_meta_path(manifest_path),
        "expression_matrix_cache": matrix_path,
        "expression_matrix_cache_metadata": cache_meta_path(matrix_path),
    }
    validate_managed_paths(
        protected_inputs=protected_inputs,
        managed_files=fixed_managed_files,
        managed_directories={"gene_cache": cache_dir},
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    if locked_manifest_path:
        sampled = validate_locked_sample_manifest(
            locked_manifest_path,
            expected_project=args.project,
            require_provider_md5=not (args.offline or args.refresh_manifest),
        )
        if args.refresh_manifest:
            sampled = refresh_locked_sample_metadata(
                sampled, expected_project=args.project
            )
        manifest = None
        n_manifest = None
        cohort_mode = "locked_manifest"
    else:
        manifest = load_or_query_manifest(
            manifest_path, args.refresh_manifest, args.source_revision, args.offline
        )
        sampled = choose_files(
            manifest,
            args.project,
            args.n_per_class,
            args.seed,
            allow_case_overlap=args.allow_case_overlap,
        )
        n_manifest = len(manifest)
        cohort_mode = "sampled_from_gdc_manifest"
    print(f"[cptac] sampled {len(sampled)} files -> {sampled_path}", file=sys.stderr)

    gene_cache_files = {}
    for file_id in sampled["file_id"].tolist():
        path = gene_cache_path(cache_dir, file_id)
        gene_cache_files[f"gene_cache:{file_id}"] = path
        gene_cache_files[f"gene_cache_metadata:{file_id}"] = cache_meta_path(path)
    validate_managed_paths(
        protected_inputs=protected_inputs,
        managed_files={**fixed_managed_files, **gene_cache_files},
        managed_directories={"gene_cache": cache_dir},
    )

    model = load_lr_model(args.weights)
    model_sha256 = sha256_file(args.weights)
    selected_genes = model["genes"].astype(str).tolist()
    matrix = build_expression_matrix(
        sampled,
        selected_genes,
        cache_dir,
        args.workers,
        args.retries,
        matrix_cache_path=matrix_path,
        refresh=args.refresh_expression_cache,
        model_sha256=model_sha256,
        source_revision=args.source_revision,
        offline=args.offline,
    )

    probabilities, alignment_report = predict_proba(
        model,
        matrix,
        max_invalid_cell_fraction=args.max_invalid_cell_fraction,
        allow_invalid_values=args.allow_invalid_values,
        return_alignment_report=True,
    )
    n_matched = int(alignment_report["n_matched_genes"])
    missing = alignment_report["missing_genes"]
    if missing:
        print(f"[cptac] WARNING: {len(missing)} missing model genes", file=sys.stderr)
    if alignment_report["invalid_matched_cells"]:
        print(
            f"[cptac] WARNING: imputed {alignment_report['invalid_matched_cells']} "
            "invalid matched cells",
            file=sys.stderr,
        )
    print(f"[cptac] matched {n_matched}/{len(selected_genes)} model genes", file=sys.stderr)

    scored = scored_dataframe(matrix.index, probabilities, args.threshold)

    predictions = sampled[["file_id", "project", "case_submitter_id",
                           "sample_submitter_id", "sample_type", "label"]].merge(
        scored, left_on="file_id", right_on="sample", how="left", validate="one_to_one"
    ).drop(columns=["sample"])
    if len(predictions) != len(sampled):
        raise ValueError(
            f"[cptac] merge produced {len(predictions)} rows, expected {len(sampled)} "
            "(duplicate or missing file_id join keys)"
        )
    n_missing = int(predictions["tumor_probability"].isna().sum())
    if n_missing:
        raise ValueError(
            f"[cptac] {n_missing}/{len(sampled)} samples have no score after the merge "
            "(file_id mismatch between the extracted matrix and the sampled manifest)"
        )

    summary = summarize(predictions["label"].to_numpy(),
                        predictions["tumor_probability"].to_numpy(),
                        args.threshold)
    sweep = threshold_sweep(predictions["label"].to_numpy(),
                            predictions["tumor_probability"].to_numpy())
    model_description = f"`{Path(args.weights).name}` (SHA256 `{model_sha256[:12]}...`)"

    cohort_details = group_audit(sampled, "file_id", "case_submitter_id")
    cohort_details.update({
        "cohort_mode": cohort_mode,
        "n_unique_biospecimens": int(sampled["sample_submitter_id"].astype(str).nunique()),
        "n_ambiguous_biospecimen_ids": int(
            sampled["sample_submitter_id"].astype(str).str.contains(";", regex=False).sum()
        ),
    })
    input_paths = {"expression_matrix_cache": matrix_path}
    if locked_manifest_path:
        input_paths["locked_sample_manifest"] = locked_manifest_path
    else:
        input_paths["gdc_manifest_cache"] = manifest_path
    gene_meta_files = list(cache_dir.glob("*.parquet.meta.json"))
    with staged_output_directory(out_dir) as stage_dir:
        staged_outputs = {
            name: stage_dir / path.name for name, path in final_outputs.items()
        }
        atomic_write_csv(staged_outputs["sampled_manifest"], sampled)
        atomic_write_csv(staged_outputs["predictions"], predictions)
        atomic_write_csv(staged_outputs["summary"], pd.DataFrame([summary]))
        atomic_write_csv(staged_outputs["threshold_sweep"], sweep)
        write_report(
            staged_outputs["report"],
            summary,
            sampled,
            predictions,
            n_manifest,
            sweep,
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
            run_kind="cptac_gdc_external_validation",
            started_at_utc=started_at,
            argv=sys.argv,
            parameters={**vars(args), "cohort_mode": cohort_mode},
            model_path=args.weights,
            sources={
                "files_endpoint": GDC_FILES_ENDPOINT,
                "data_endpoint": GDC_DATA_ENDPOINT,
                "revision": args.source_revision,
            },
            inputs=input_paths,
            outputs=manifest_outputs,
            output_display_paths=manifest_final_outputs,
            alignment=alignment_report,
            cohort_audit=cohort_details,
            cache_details={
                "expression_matrix": read_cache_metadata(matrix_path) or {},
                "gene_cache": {
                    "directory": str(cache_dir.resolve()),
                    "metadata_files": int(len(gene_meta_files)),
                },
            },
            source_code={
                "validator": Path(__file__).resolve(),
                "metrics_core": ROOT / "tcga_rnaseq" / "metrics.py",
                "score_core": ROOT / "tcga_rnaseq" / "score.py",
                "alignment_core": ROOT / "tcga_rnaseq" / "align.py",
            },
        )
        publish_staged_files(staged_outputs, final_outputs)

    print(f"[cptac] summary -> {final_outputs['summary']}", file=sys.stderr)
    print(f"[cptac] threshold sweep -> {final_outputs['threshold_sweep']}", file=sys.stderr)
    print(f"[cptac] predictions -> {final_outputs['predictions']}", file=sys.stderr)
    print(f"[cptac] report -> {final_outputs['report']}", file=sys.stderr)
    print(f"[cptac] provenance -> {final_outputs['run_manifest']}", file=sys.stderr)
    print(pd.DataFrame([summary]).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
