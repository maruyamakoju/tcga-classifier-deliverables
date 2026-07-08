#!/usr/bin/env python3
"""Sanity check the deployed model on TCGA samples from UCSC Xena Toil.

The purpose is to separate two effects:

1. biology/source cohort transfer (GTEx healthy normals), and
2. expression-pipeline transfer (GDC STAR-Counts training -> Toil/RSEM values).

If TCGA samples processed by Toil/RSEM are also mis-thresholded, then the GTEx
failure is primarily a pipeline/normalization incompatibility rather than a
statement about healthy tissue biology.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
from sklearn.metrics import (accuracy_score, average_precision_score,
                             confusion_matrix, f1_score, precision_score,
                             recall_score, roc_auc_score, roc_curve)

ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(SCRIPT_DIR))

from score_tumor_normal import load_pipeline, score_dataframe  # noqa: E402
from validate_gtex_xena import (PHENOTYPE_URL, extract_matrix_from_xena,  # noqa: E402
                                load_or_download_phenotype)

TCGA_TPM_URL = "https://toil.xenahubs.net/download/tcga_RSEM_gene_tpm.gz"
LABEL_MAP = {"Solid Tissue Normal": 0, "Primary Tumor": 1}


def choose_tcga_samples(phenotype: pd.DataFrame, n_per_class: int,
                        seed: int) -> pd.DataFrame:
    tcga = phenotype[phenotype["_study"] == "TCGA"].copy()
    tcga = tcga[tcga["_sample_type"].isin(LABEL_MAP)].copy()
    sampled = []
    for sample_type in LABEL_MAP:
        part = tcga[tcga["_sample_type"] == sample_type].copy()
        if part.empty:
            raise ValueError(f"No {sample_type} samples found in Toil phenotype")
        sampled.append(part.sample(n=min(n_per_class, len(part)), random_state=seed))
    out = pd.concat(sampled, ignore_index=True)
    out["label"] = out["_sample_type"].map(LABEL_MAP)
    return out.sample(frac=1, random_state=seed).reset_index(drop=True)


def summarize(predictions: pd.DataFrame, threshold: float) -> dict:
    y_true = predictions["label"].to_numpy()
    scores = predictions["tumor_probability"].to_numpy()
    pred = (scores >= threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, pred).ravel()
    return {
        "n": int(len(predictions)),
        "n_tumor": int(y_true.sum()),
        "n_normal": int((y_true == 0).sum()),
        "auc": roc_auc_score(y_true, scores),
        "average_precision": average_precision_score(y_true, scores),
        "accuracy": accuracy_score(y_true, pred),
        "f1": f1_score(y_true, pred),
        "precision": precision_score(y_true, pred, zero_division=0),
        "recall": recall_score(y_true, pred),
        "tn": int(tn),
        "fp": int(fp),
        "fn": int(fn),
        "tp": int(tp),
        "threshold": threshold,
    }


def threshold_sweep(predictions: pd.DataFrame) -> pd.DataFrame:
    y_true = predictions["label"].to_numpy()
    scores = predictions["tumor_probability"].to_numpy()
    fpr, tpr, thresholds = roc_curve(y_true, scores)
    best_idx = int((tpr - fpr).argmax())
    rows = []
    for name, threshold in [
        ("default_0.5", 0.5),
        ("youden_j", float(thresholds[best_idx])),
        ("high_0.99", 0.99),
        ("high_0.999", 0.999),
        ("high_0.9999", 0.9999),
        ("high_0.99999", 0.99999),
    ]:
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


def write_report(path: Path, summary: dict, predictions: pd.DataFrame,
                 sweep: pd.DataFrame) -> None:
    by_label = predictions.groupby("_sample_type").agg(
        n=("sample", "size"),
        mean_tumor_probability=("tumor_probability", "mean"),
        median_tumor_probability=("tumor_probability", "median"),
        tumor_calls=("call", lambda s: int((s == "tumor").sum())),
    )
    by_label_rows = "\n".join(
        f"| {idx} | {int(row.n)} | {row.mean_tumor_probability:.4f} | "
        f"{row.median_tumor_probability:.4f} | {int(row.tumor_calls)} |"
        for idx, row in by_label.iterrows()
    )
    sweep_rows = "\n".join(
        f"| {row.threshold_name} | {row.threshold:.6f} | {row.accuracy:.4f} | "
        f"{row.precision:.4f} | {row.recall:.4f} | {int(row.tn)} / {int(row.fp)} | "
        f"{int(row.fn)} / {int(row.tp)} |"
        for row in sweep.itertuples(index=False)
    )
    text = f"""# TCGA Toil/RSEM pipeline-transfer sanity check

## Data source

- Source: UCSC Xena Toil RNA-seq recompute compendium
- Dataset: `tcga_RSEM_gene_tpm`
- Phenotype table: `TcgaTargetGTEX_phenotype`
- Samples scored: {summary['n']} TCGA samples ({summary['n_tumor']} primary tumor, {summary['n_normal']} solid tissue normal)
- Input transform: Xena log2(TPM+0.001) -> TPM -> log2(TPM+1)
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

| Sample type | n | Mean tumor probability | Median tumor probability | Tumor calls |
|---|---:|---:|---:|---:|
{by_label_rows}

## Threshold sensitivity

| Threshold | Cutoff | Accuracy | Precision | Recall | TN / FP | FN / TP |
|---|---:|---:|---:|---:|---:|---:|
{sweep_rows}

## Interpretation

This is not a new biological validation cohort because TCGA overlaps the model's
training source. It is a pipeline-transfer sanity check. Poor hard-call behavior
here means the deployed GDC STAR-Counts model should not be applied directly to
Toil/RSEM matrices without refitting or recalibrating on that expression
pipeline.
"""
    path.write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-per-class", type=int, default=100)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument("--out-dir", default=str(ROOT / "external-validation" / "tcga_toil_xena"))
    parser.add_argument("--pipeline", default=str(ROOT / "deployable_pipeline.pkl"))
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    phenotype_path = out_dir / "TcgaTargetGTEX_phenotype.csv"
    matrix_path = out_dir / "tcga_toil_selected_genes_model_scale.pkl"
    sampled_path = out_dir / "sampled_tcga_toil_manifest.csv"

    phenotype = load_or_download_phenotype(phenotype_path, args.refresh)
    sampled = choose_tcga_samples(phenotype, args.n_per_class, args.seed)
    sampled.to_csv(sampled_path, index=False)
    print(f"[tcga-toil] sampled {len(sampled)} TCGA Toil samples", file=sys.stderr)

    pipe = load_pipeline(args.pipeline)
    matrix = extract_matrix_from_xena(sampled["sample"].tolist(), pipe["selected_genes"],
                                      matrix_path, args.refresh, TCGA_TPM_URL)
    matrix = matrix.fillna(pd.Series(pipe["scaler"].mean_, index=pipe["selected_genes"]))

    scored, _, _ = score_dataframe(matrix, pipe, "lr", args.threshold)
    predictions = sampled.merge(scored, on="sample", how="left")
    predictions_path = out_dir / "tcga_toil_predictions.csv"
    predictions.to_csv(predictions_path, index=False)

    summary = summarize(predictions, args.threshold)
    summary_path = out_dir / "tcga_toil_summary.csv"
    pd.DataFrame([summary]).to_csv(summary_path, index=False)
    sweep = threshold_sweep(predictions)
    sweep_path = out_dir / "tcga_toil_threshold_sweep.csv"
    sweep.to_csv(sweep_path, index=False)

    report_path = out_dir / "TCGA_TOIL_PIPELINE_CHECK.md"
    write_report(report_path, summary, predictions, sweep)

    print(f"[tcga-toil] summary -> {summary_path}", file=sys.stderr)
    print(f"[tcga-toil] threshold sweep -> {sweep_path}", file=sys.stderr)
    print(f"[tcga-toil] predictions -> {predictions_path}", file=sys.stderr)
    print(f"[tcga-toil] report -> {report_path}", file=sys.stderr)
    print(pd.DataFrame([summary]).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
