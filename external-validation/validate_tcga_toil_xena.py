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

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(SCRIPT_DIR))

from tcga_rnaseq import load_lr_model, predict_proba, validate_threshold  # noqa: E402
from tcga_rnaseq import metrics as M  # noqa: E402
from provenance import (  # noqa: E402
    atomic_write_csv,
    atomic_write_text,
    cache_meta_path,
    group_audit,
    publish_staged_files,
    read_cache_metadata,
    scored_dataframe,
    sha256_file,
    staged_output_directory,
    utc_now,
    validate_identifier_column,
    validate_managed_paths,
    validate_source_revision,
    write_run_manifest,
)
from validate_gtex_xena import (donor_id_from_sample,  # noqa: E402
                                PHENOTYPE_URL,
                                extract_matrix_from_xena,
                                load_locked_sample_manifest,
                                load_or_download_phenotype,
                                require_complete_merge)

TCGA_TPM_URL = "https://toil.xenahubs.net/download/tcga_RSEM_gene_tpm.gz"
TCGA_DATASET_ID = "tcga_RSEM_gene_tpm"
LABEL_MAP = {"Solid Tissue Normal": 0, "Primary Tumor": 1}


def validate_locked_tcga_sample_manifest(path: Path) -> pd.DataFrame:
    sampled = load_locked_sample_manifest(path, "TCGA")
    expected_labels = sampled["_sample_type"].map(LABEL_MAP).astype(int)
    if "label" in sampled.columns:
        labels = pd.to_numeric(sampled["label"], errors="coerce")
        if labels.isna().any() or not np.array_equal(
            labels.astype(int).to_numpy(), expected_labels.to_numpy()
        ):
            raise ValueError("locked TCGA labels do not agree with _sample_type")
    sampled["label"] = expected_labels
    if sampled["label"].nunique() != 2:
        raise ValueError("locked TCGA sample manifest must contain both tumor and normal")
    return sampled


def choose_tcga_samples(
    phenotype: pd.DataFrame,
    n_per_class: int,
    seed: int,
    allow_donor_overlap: bool = False,
) -> pd.DataFrame:
    if n_per_class <= 0:
        raise ValueError("n_per_class must be a positive integer")
    validate_identifier_column(phenotype, "sample", "Xena phenotype")
    tcga = phenotype[phenotype["_study"] == "TCGA"].copy()
    tcga = tcga[tcga["_sample_type"].isin(LABEL_MAP)].copy()
    tcga["donor_id"] = tcga["sample"].map(donor_id_from_sample)
    sampled = []
    used_donors: set[str] = set()
    for sample_type in LABEL_MAP:
        part = tcga[tcga["_sample_type"] == sample_type].copy()
        part = part.sort_values("sample").drop_duplicates("donor_id", keep="first")
        if not allow_donor_overlap:
            part = part[~part["donor_id"].isin(used_donors)]
        if part.empty:
            raise ValueError(f"No {sample_type} samples found in Toil phenotype")
        n = min(n_per_class, len(part))
        sampled.append(part.sample(n=n, random_state=seed))
        used_donors.update(sampled[-1]["donor_id"].astype(str))
    out = pd.concat(sampled, ignore_index=True)
    out["label"] = out["_sample_type"].map(LABEL_MAP)
    return out.sample(frac=1, random_state=seed).reset_index(drop=True)


def summarize(predictions: pd.DataFrame, threshold: float) -> dict:
    y_true = predictions["label"].to_numpy()
    scores = predictions["tumor_probability"].to_numpy()
    metrics = M.classification_metrics(y_true, scores, threshold)
    return {
        "n": int(len(predictions)),
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


def threshold_sweep(predictions: pd.DataFrame) -> pd.DataFrame:
    y_true = predictions["label"].to_numpy()
    scores = predictions["tumor_probability"].to_numpy()
    youden_threshold = M.youden_threshold(y_true, scores)["threshold"]
    rows = []
    for name, threshold in [
        ("default_0.5", 0.5),
        ("youden_j", float(youden_threshold)),
        ("high_0.99", 0.99),
        ("high_0.999", 0.999),
        ("high_0.9999", 0.9999),
        ("high_0.99999", 0.99999),
    ]:
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


def write_report(path: Path, summary: dict, predictions: pd.DataFrame,
                 sweep: pd.DataFrame, model_description: str) -> None:
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
    atomic_write_text(path, text)


def main() -> int:
    started_at = utc_now()
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-per-class", type=int, default=100)
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
        help="allow tumor and normal rows from the same TCGA donor in a newly sampled panel",
    )
    parser.add_argument(
        "--source-revision",
        default="unversioned",
        help="provider snapshot/revision identifier included in cache provenance",
    )
    parser.add_argument("--out-dir", default=str(ROOT / "external-validation" / "tcga_toil_xena"))
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
        if args.offline and args.refresh:
            raise ValueError("--offline/--cache-only cannot be combined with --refresh")
        args.source_revision = validate_source_revision(
            args.source_revision, live=not args.offline
        )
    except ValueError as exc:
        parser.error(str(exc))

    out_dir = Path(args.out_dir).resolve(strict=False)
    phenotype_path = out_dir / "TcgaTargetGTEX_phenotype.csv"
    matrix_path = out_dir / "tcga_toil_selected_genes_model_scale.parquet"
    sampled_path = out_dir / "sampled_tcga_toil_manifest.csv"

    locked_manifest_path = Path(args.sample_manifest).resolve() if args.sample_manifest else None
    final_outputs = {
        "sampled_manifest": sampled_path,
        "predictions": out_dir / "tcga_toil_predictions.csv",
        "summary": out_dir / "tcga_toil_summary.csv",
        "threshold_sweep": out_dir / "tcga_toil_threshold_sweep.csv",
        "report": out_dir / "TCGA_TOIL_PIPELINE_CHECK.md",
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
        sampled = validate_locked_tcga_sample_manifest(locked_manifest_path)
        cohort_mode = "locked_manifest"
    else:
        phenotype = load_or_download_phenotype(
            phenotype_path, args.refresh, args.source_revision, args.offline
        )
        sampled = choose_tcga_samples(
            phenotype,
            args.n_per_class,
            args.seed,
            allow_donor_overlap=args.allow_donor_overlap,
        )
        cohort_mode = "sampled_from_phenotype"
    print(f"[tcga-toil] sampled {len(sampled)} TCGA Toil samples", file=sys.stderr)

    model = load_lr_model(args.weights)
    model_sha256 = sha256_file(args.weights)
    selected_genes = model["genes"].astype(str).tolist()
    source_identity = {
        "dataset_id": TCGA_DATASET_ID,
        "url": TCGA_TPM_URL,
        "revision": args.source_revision,
    }
    matrix = extract_matrix_from_xena(sampled["sample"].tolist(), selected_genes,
                                      matrix_path, args.refresh, TCGA_TPM_URL,
                                      model_sha256=model_sha256,
                                      source_identity=source_identity,
                                      offline=args.offline)

    probabilities, alignment_report = predict_proba(
        model,
        matrix,
        max_invalid_cell_fraction=args.max_invalid_cell_fraction,
        allow_invalid_values=args.allow_invalid_values,
        return_alignment_report=True,
    )
    if alignment_report["invalid_matched_cells"]:
        print(
            f"[tcga-toil] WARNING: imputed {alignment_report['invalid_matched_cells']} "
            "invalid matched cells",
            file=sys.stderr,
        )
    scored = scored_dataframe(matrix.index, probabilities, args.threshold)
    predictions = sampled.merge(scored, on="sample", how="left", validate="one_to_one")
    require_complete_merge(predictions, len(sampled), "tumor_probability", "[tcga-toil]")

    summary = summarize(predictions, args.threshold)
    sweep = threshold_sweep(predictions)

    model_description = f"`{Path(args.weights).name}` (SHA256 `{model_sha256[:12]}...`)"

    cohort_details = group_audit(sampled, "sample", "donor_id")
    cohort_details["cohort_mode"] = cohort_mode
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
        write_report(
            staged_outputs["report"], summary, predictions, sweep, model_description
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
            run_kind="tcga_toil_pipeline_validation",
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
            cache_details={"expression_matrix": read_cache_metadata(matrix_path) or {}},
            source_code={
                "validator": Path(__file__).resolve(),
                "xena_extractor": SCRIPT_DIR / "validate_gtex_xena.py",
                "metrics_core": ROOT / "tcga_rnaseq" / "metrics.py",
                "score_core": ROOT / "tcga_rnaseq" / "score.py",
            },
        )
        publish_staged_files(staged_outputs, final_outputs)

    print(f"[tcga-toil] summary -> {final_outputs['summary']}", file=sys.stderr)
    print(f"[tcga-toil] threshold sweep -> {final_outputs['threshold_sweep']}", file=sys.stderr)
    print(f"[tcga-toil] predictions -> {final_outputs['predictions']}", file=sys.stderr)
    print(f"[tcga-toil] report -> {final_outputs['report']}", file=sys.stderr)
    print(f"[tcga-toil] provenance -> {final_outputs['run_manifest']}", file=sys.stderr)
    print(pd.DataFrame([summary]).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
