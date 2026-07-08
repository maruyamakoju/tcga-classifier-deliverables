#!/usr/bin/env python3
"""External smoke validation on CPTAC STAR-Counts files from the GDC.

This is intentionally a lightweight external check: CPTAC is not TCGA, but the
files are harmonized by the same GDC STAR-Counts pipeline, so input scale and
gene identifiers match the deployed model. The script downloads only the
selected 2,000 genes per file and caches those extracted vectors.
"""
from __future__ import annotations

import argparse
import math
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from sklearn.metrics import (accuracy_score, average_precision_score,
                             confusion_matrix, f1_score, precision_score,
                             recall_score, roc_auc_score, roc_curve)

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from score_tumor_normal import load_pipeline, score_dataframe  # noqa: E402

GDC_FILES_ENDPOINT = "https://api.gdc.cancer.gov/files"
GDC_DATA_ENDPOINT = "https://api.gdc.cancer.gov/data"
LABEL_MAP = {"Solid Tissue Normal": 0, "Primary Tumor": 1}


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
            {"op": "in", "content": {"field": "cases.samples.sample_type",
                                       "value": list(LABEL_MAP)}},
        ],
    }
    fields = ",".join([
        "file_id",
        "file_name",
        "file_size",
        "cases.case_id",
        "cases.submitter_id",
        "cases.project.project_id",
        "cases.samples.sample_type",
        "cases.samples.submitter_id",
    ])
    payload = {"filters": filters, "fields": fields, "format": "json", "size": 5000}
    response = requests.post(GDC_FILES_ENDPOINT, json=payload, timeout=60)
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


def load_or_query_manifest(path: Path, refresh: bool) -> pd.DataFrame:
    if path.exists() and not refresh:
        return pd.read_csv(path)
    manifest = query_cptac_manifest()
    path.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(path, index=False)
    return manifest


def choose_files(manifest: pd.DataFrame, project: str, n_per_class: int,
                 seed: int) -> pd.DataFrame:
    df = manifest[manifest["project"] == project].copy()
    if df.empty:
        raise ValueError(f"No files found for project {project}")

    sampled = []
    for sample_type in LABEL_MAP:
        part = df[df["sample_type"] == sample_type].copy()
        if part.empty:
            raise ValueError(f"No {sample_type} files found for {project}")
        n = min(n_per_class, len(part))
        sampled.append(part.sample(n=n, random_state=seed).copy())
    out = pd.concat(sampled, ignore_index=True)
    out["label"] = out["sample_type"].map(LABEL_MAP)
    return out.sample(frac=1, random_state=seed).reset_index(drop=True)


def extract_selected_genes(file_id: str, selected_genes: list[str],
                           cache_dir: Path, retries: int = 4) -> pd.Series:
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / f"{file_id}.pkl"
    if cache_path.exists():
        return pd.read_pickle(cache_path)

    wanted = set(selected_genes)
    wanted_by_base = {gene.split(".")[0]: gene for gene in selected_genes}
    values = {}
    url = f"{GDC_DATA_ENDPOINT}/{file_id}"
    last_error = None
    for attempt in range(1, retries + 1):
        values = {}
        try:
            with requests.get(url, stream=True, timeout=120) as response:
                response.raise_for_status()
                header = None
                gene_idx = tpm_idx = None
                for raw_line in response.iter_lines(decode_unicode=True):
                    if isinstance(raw_line, bytes):
                        raw_line = raw_line.decode("utf-8", errors="replace")
                    if not raw_line or raw_line.startswith("#"):
                        continue
                    parts = raw_line.split("\t")
                    if header is None:
                        header = parts
                        gene_idx = header.index("gene_id")
                        tpm_idx = header.index("tpm_unstranded")
                        continue

                    gene_id = parts[gene_idx]
                    target = gene_id if gene_id in wanted else wanted_by_base.get(gene_id.split(".")[0])
                    if target is None:
                        continue
                    try:
                        tpm = float(parts[tpm_idx])
                    except ValueError:
                        continue
                    values[target] = math.log2(tpm + 1.0)
                    if len(values) == len(selected_genes):
                        break
            break
        except requests.RequestException as exc:
            last_error = exc
            if attempt == retries:
                raise
            time.sleep(2 ** (attempt - 1))
            print(f"[cptac] retry {attempt}/{retries} for {file_id}: {exc}",
                  file=sys.stderr)
    if last_error is not None and not values:
        raise last_error

    series = pd.Series(values, name=file_id, dtype=float)
    pd.to_pickle(series, cache_path)
    return series


def build_expression_matrix(files: pd.DataFrame, selected_genes: list[str],
                            cache_dir: Path, workers: int, retries: int) -> pd.DataFrame:
    rows = {}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        future_to_file = {
            pool.submit(extract_selected_genes, row.file_id, selected_genes, cache_dir,
                        retries): row.file_id
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
    return matrix.reindex(files["file_id"])[selected_genes]


def summarize(y_true: np.ndarray, scores: np.ndarray, threshold: float) -> dict:
    pred = (scores >= threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, pred).ravel()
    return {
        "n": int(len(y_true)),
        "n_tumor": int(y_true.sum()),
        "n_normal": int((y_true == 0).sum()),
        "auc": roc_auc_score(y_true, scores),
        "average_precision": average_precision_score(y_true, scores),
        "accuracy": accuracy_score(y_true, pred),
        "f1": f1_score(y_true, pred),
        "precision": precision_score(y_true, pred),
        "recall": recall_score(y_true, pred),
        "tn": int(tn),
        "fp": int(fp),
        "fn": int(fn),
        "tp": int(tp),
        "threshold": threshold,
    }


def threshold_sweep(y_true: np.ndarray, scores: np.ndarray) -> pd.DataFrame:
    fpr, tpr, thresholds = roc_curve(y_true, scores)
    best_idx = int(np.argmax(tpr - fpr))
    threshold_specs = [
        ("default_0.5", 0.5),
        ("youden_j", float(thresholds[best_idx])),
        ("high_specificity_0.75", 0.75),
        ("high_specificity_0.9", 0.9),
        ("high_specificity_0.95", 0.95),
        ("very_high_specificity_0.99", 0.99),
    ]
    rows = []
    for name, threshold in threshold_specs:
        pred = (scores >= threshold).astype(int)
        tn, fp, fn, tp = confusion_matrix(y_true, pred).ravel()
        rows.append({
            "threshold_name": name,
            "threshold": threshold,
            "accuracy": accuracy_score(y_true, pred),
            "f1": f1_score(y_true, pred),
            "precision": precision_score(y_true, pred, zero_division=0),
            "recall": recall_score(y_true, pred),
            "tn": int(tn),
            "fp": int(fp),
            "fn": int(fn),
            "tp": int(tp),
        })
    return pd.DataFrame(rows)


def write_report(path: Path, summary: dict, sampled: pd.DataFrame,
                 predictions: pd.DataFrame, n_manifest: int,
                 sweep: pd.DataFrame) -> None:
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
    text = f"""# CPTAC external smoke validation

## Data source

- Source: NCI GDC API, CPTAC program
- Files queried: {n_manifest}
- Project scored: {sampled['project'].iloc[0]}
- Sampled files: {len(sampled)} ({summary['n_tumor']} primary tumor, {summary['n_normal']} solid tissue normal)
- Workflow: GDC STAR-Counts, `tpm_unstranded` converted to log2(TPM+1)
- Model: bundled logistic regression from `deployable_pipeline.pkl`

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
    path.write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", default="CPTAC-3")
    parser.add_argument("--n-per-class", type=int, default=50)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--retries", type=int, default=4)
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--refresh-manifest", action="store_true")
    parser.add_argument("--out-dir", default=str(ROOT / "external-validation" / "cptac_gdc"))
    parser.add_argument("--pipeline", default=str(ROOT / "deployable_pipeline.pkl"))
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = out_dir / "cptac_gdc_manifest.csv"
    cache_dir = out_dir / "gene_cache"

    manifest = load_or_query_manifest(manifest_path, args.refresh_manifest)
    sampled = choose_files(manifest, args.project, args.n_per_class, args.seed)
    sampled_path = out_dir / "sampled_manifest.csv"
    sampled.to_csv(sampled_path, index=False)
    print(f"[cptac] sampled {len(sampled)} files -> {sampled_path}", file=sys.stderr)

    pipe = load_pipeline(args.pipeline)
    selected_genes = pipe["selected_genes"]
    matrix = build_expression_matrix(sampled, selected_genes, cache_dir, args.workers,
                                     args.retries)
    matrix_path = out_dir / "expression_selected_genes.pkl"
    matrix.to_pickle(matrix_path)

    scored, n_matched, missing = score_dataframe(matrix, pipe, "lr", args.threshold)
    if missing:
        print(f"[cptac] WARNING: {len(missing)} missing model genes", file=sys.stderr)
    print(f"[cptac] matched {n_matched}/{len(selected_genes)} model genes", file=sys.stderr)

    predictions = sampled[["file_id", "project", "case_submitter_id",
                           "sample_submitter_id", "sample_type", "label"]].merge(
        scored, left_on="file_id", right_on="sample", how="left"
    ).drop(columns=["sample"])
    predictions_path = out_dir / "cptac_predictions.csv"
    predictions.to_csv(predictions_path, index=False)

    summary = summarize(predictions["label"].to_numpy(),
                        predictions["tumor_probability"].to_numpy(),
                        args.threshold)
    summary_path = out_dir / "cptac_summary.csv"
    pd.DataFrame([summary]).to_csv(summary_path, index=False)
    sweep = threshold_sweep(predictions["label"].to_numpy(),
                            predictions["tumor_probability"].to_numpy())
    sweep_path = out_dir / "cptac_threshold_sweep.csv"
    sweep.to_csv(sweep_path, index=False)
    report_path = out_dir / "CPTAC_EXTERNAL_VALIDATION.md"
    write_report(report_path, summary, sampled, predictions, len(manifest), sweep)

    print(f"[cptac] summary -> {summary_path}", file=sys.stderr)
    print(f"[cptac] threshold sweep -> {sweep_path}", file=sys.stderr)
    print(f"[cptac] predictions -> {predictions_path}", file=sys.stderr)
    print(f"[cptac] report -> {report_path}", file=sys.stderr)
    print(pd.DataFrame([summary]).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
