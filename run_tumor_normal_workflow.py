#!/usr/bin/env python3
"""Run the complete lightweight tumor-vs-normal scoring workflow."""
import argparse
import contextlib
import hashlib
import os
import shutil
import sys
import tempfile
from pathlib import Path

import pandas as pd

from calibrate_threshold import (
    build_calibration_summary,
    choose_youden_threshold,
    load_scores_and_labels,
    metrics_at_threshold,
)
from explain_scores import explain_dataframe, load_gene_metadata
from inspect_expression_input import inspect_dataframe, load_reference
from score_tumor_normal import (
    load_lr_weights,
    score_dataframe_lr_weights,
)
from tcga_rnaseq import (
    ensure_distinct_paths,
    print_invalid_alignment_summary,
    read_matrix,
    validate_alignment_report,
    validate_threshold,
    write_dataframe_csv,
    write_json,
    write_text,
)


MANAGED_OUTPUT_NAMES = {
    "qc.json",
    "scores.csv",
    "workflow_report.md",
    "manifest.json",
    "thresholds.csv",
    "calibration.json",
    "explanations.csv",
}


def sha256_file(path):
    digest = hashlib.sha256()
    try:
        with Path(path).open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
    except OSError as exc:
        raise ValueError(f"could not hash workflow input {path}: {exc}") from exc
    return digest.hexdigest()


def snapshot_inputs(paths):
    """Hash every file consumed by the workflow before computation."""
    records = {}
    for name, raw_path in paths.items():
        if not raw_path:
            continue
        path = Path(raw_path).resolve()
        if not path.is_file():
            raise ValueError(f"workflow input is not a regular file ({name}): {path}")
        try:
            size = int(path.stat().st_size)
        except OSError as exc:
            raise ValueError(f"could not inspect workflow input ({name}) {path}: {exc}") from exc
        try:
            display_path = path.relative_to(Path.cwd().resolve()).as_posix()
        except ValueError:
            display_path = str(path)
        records[name] = {
            "path": display_path,
            "bytes": size,
            "sha256": sha256_file(path),
            "_resolved_path": str(path),
        }
    return records


def verify_input_snapshot(records):
    """Refuse publication if an input changed while outputs were computed."""
    for name, expected in records.items():
        path = Path(expected["_resolved_path"])
        if not path.is_file():
            raise ValueError(f"workflow input changed during the run ({name}): {path}")
        try:
            size = int(path.stat().st_size)
        except OSError as exc:
            raise ValueError(
                f"could not verify workflow input snapshot ({name}) {path}: {exc}"
            ) from exc
        if size != expected["bytes"] or sha256_file(path) != expected["sha256"]:
            raise ValueError(f"workflow input changed during the run ({name}): {path}")


def manifest_provenance(records):
    """Stable top-level hashes plus the complete consumed-input snapshot."""
    public_records = {
        name: {key: value for key, value in record.items() if not key.startswith("_")}
        for name, record in records.items()
    }
    return {
        "input_sha256": records["expression_input"]["sha256"],
        "model_sha256": records["lr_weights"]["sha256"],
        "input_artifacts": public_records,
    }


def _remove_managed_path(path):
    path = Path(path)
    if path.is_symlink() or path.is_file():
        path.unlink(missing_ok=True)
    elif path.exists():
        shutil.rmtree(path)


def _replace_file(source, destination):
    """Small indirection used by rollback fault-injection tests."""
    os.replace(source, destination)


@contextlib.contextmanager
def staged_generation(output_dir):
    """Yield a same-filesystem sibling stage and always clean it afterward."""
    output_dir = Path(output_dir).resolve()
    if output_dir.exists() and not output_dir.is_dir():
        raise ValueError(f"workflow output path is not a directory: {output_dir}")
    try:
        output_dir.parent.mkdir(parents=True, exist_ok=True)
        stage = Path(
            tempfile.mkdtemp(
                prefix=f".{output_dir.name}.workflow-stage.",
                dir=output_dir.parent,
            )
        )
    except OSError as exc:
        raise ValueError(f"could not create workflow staging directory: {exc}") from exc
    try:
        yield stage
    finally:
        if stage.exists():
            shutil.rmtree(stage, ignore_errors=True)


def commit_generation(stage_dir, output_dir, staged_names):
    """Rollback-capable managed-file commit with the manifest installed last."""
    stage_dir = Path(stage_dir).resolve()
    output_dir = Path(output_dir).resolve()
    staged_names = list(staged_names)
    if not staged_names or staged_names[-1] != "manifest.json":
        raise ValueError("workflow manifest must be staged and committed last")
    if len(set(staged_names)) != len(staged_names):
        raise ValueError("workflow staged output names must be unique")
    if not set(staged_names) <= MANAGED_OUTPUT_NAMES:
        raise ValueError("workflow stage contains an unmanaged output name")
    for name in staged_names:
        if Path(name).name != name or not (stage_dir / name).is_file():
            raise ValueError(f"invalid or missing staged workflow output: {name}")

    try:
        output_dir.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise ValueError(f"could not create workflow output directory: {exc}") from exc
    for name in MANAGED_OUTPUT_NAMES:
        path = output_dir / name
        if path.exists() and path.is_dir() and not path.is_symlink():
            raise ValueError(f"managed workflow output path is a directory: {path}")

    try:
        backup = Path(
            tempfile.mkdtemp(prefix=".workflow-backup.", dir=output_dir.parent)
        )
    except OSError as exc:
        raise ValueError(f"could not create workflow rollback directory: {exc}") from exc
    backed_up = []
    installed = []
    completed = False
    rolled_back = False
    try:
        backup_order = ["manifest.json"] + sorted(
            MANAGED_OUTPUT_NAMES - {"manifest.json"}
        )
        for name in backup_order:
            source = output_dir / name
            if source.exists() or source.is_symlink():
                _replace_file(source, backup / name)
                backed_up.append(name)
        for name in staged_names:
            _replace_file(stage_dir / name, output_dir / name)
            installed.append(name)
        completed = True
    except BaseException as exc:
        rollback_errors = []
        for name in reversed(installed):
            try:
                _remove_managed_path(output_dir / name)
            except BaseException as rollback_exc:
                rollback_errors.append(str(rollback_exc))
        for name in backed_up:
            previous = backup / name
            if not previous.exists() and not previous.is_symlink():
                continue
            try:
                _replace_file(previous, output_dir / name)
            except BaseException as rollback_exc:
                rollback_errors.append(str(rollback_exc))
        rolled_back = not rollback_errors
        if rollback_errors:
            raise RuntimeError(
                "workflow commit failed and rollback was incomplete: "
                + "; ".join(rollback_errors)
                + f"; recovery files remain in {backup}"
            ) from exc
        if isinstance(exc, ValueError):
            raise
        if isinstance(exc, Exception):
            raise ValueError(f"could not commit workflow outputs: {exc}") from exc
        raise
    finally:
        if backup.exists() and (completed or rolled_back):
            shutil.rmtree(backup, ignore_errors=True)


def publish_generation(stage_dir, output_dir, manifest, outputs, input_records):
    """Validate and commit one coherent terminal workflow generation."""
    stage_dir = Path(stage_dir)
    output_names = list(outputs.values())
    if len(set(output_names)) != len(output_names):
        raise ValueError("workflow manifest output paths must be unique")
    if set(output_names) - MANAGED_OUTPUT_NAMES:
        raise ValueError("workflow manifest references an unmanaged output")
    if "manifest.json" not in output_names:
        raise ValueError("workflow generation must publish manifest.json")
    for name in output_names:
        if name == "manifest.json":
            continue
        if not (stage_dir / name).is_file():
            raise ValueError(f"workflow manifest references a missing staged output: {name}")

    manifest = {
        **manifest,
        **manifest_provenance(input_records),
        "outputs": outputs,
    }
    write_json(manifest, stage_dir / "manifest.json", sort_keys=True)
    verify_input_snapshot(input_records)
    staged_names = sorted(name for name in output_names if name != "manifest.json")
    staged_names.append("manifest.json")
    commit_generation(stage_dir, output_dir, staged_names)
    return manifest


def calibration_from_labels(scores_path, labels_path, sample_column, label_column,
                            default_threshold, extra_thresholds,
                            min_match_fraction):
    data = load_scores_and_labels(str(scores_path), str(labels_path), sample_column,
                                  label_column, min_match_fraction)
    y_true = data["label_binary"].to_numpy(dtype=int)
    scores = data["tumor_probability"].to_numpy(dtype=float)

    best = choose_youden_threshold(y_true, scores)
    rows = [
        metrics_at_threshold(y_true, scores, default_threshold, "default"),
        best,
    ]
    for threshold in extra_thresholds:
        rows.append(metrics_at_threshold(y_true, scores, threshold, f"threshold_{threshold:g}"))
    metrics = pd.DataFrame(rows)
    summary = build_calibration_summary(y_true, scores, best)
    return metrics, summary


def markdown_table(rows):
    if not rows:
        return ""
    headers = list(rows[0].keys())
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(row.get(h, "")) for h in headers) + " |")
    return "\n".join(lines)


def format_float(value, digits=4):
    if value is None:
        return "NA"
    return f"{float(value):.{digits}f}"


def build_report(input_path, qc, scores, out_files, calibration_summary=None,
                 explanations_rows=None, stop_reason="QC", calibration_error=None):
    gene = qc["gene_match"]
    value = qc["value_summary"]
    dist = qc["distribution_summary"]
    score_q = qc["score_summary"]["tumor_probability"]
    call_counts = scores["call"].value_counts().to_dict() if scores is not None else {}

    lines = [
        "# Tumor-vs-normal workflow report",
        "",
        f"Input: `{input_path}`",
        "",
        "## QC",
        "",
        f"- Status: **{qc['status']}**",
        f"- Model genes matched: {gene['matched_model_genes']}/{qc['shape']['model_genes']} "
        f"({gene['match_rate']:.1%})",
        f"- Expression median / p99 / max: {format_float(value['p50'])} / "
        f"{format_float(value['p99'])} / {format_float(value['max'])}",
        f"- |z| > 6 fraction: {format_float(dist['abs_z_gt_6_fraction'], 6)}",
        f"- Cohort gene-mean |z| p99: "
        f"{format_float(dist['cohort_gene_mean_abs_z']['p99'])}",
        "",
    ]
    if qc["messages"]:
        lines.append("### QC messages")
        lines.append("")
        for msg in qc["messages"]:
            lines.append(f"- {msg['level']}: {msg['message']}")
        lines.append("")
    else:
        lines.append("No QC warnings or errors.")
        lines.append("")

    lines.extend(["## Scores", ""])
    if scores is None:
        lines.extend([
            f"- Scoring was not run because the workflow stopped after {stop_reason}.",
            f"- Input samples inspected: {qc['shape']['samples']}",
            "",
        ])
    else:
        lines.extend([
            f"- Samples: {len(scores)}",
            f"- Tumor calls: {int(call_counts.get('tumor', 0))}",
            f"- Normal calls: {int(call_counts.get('normal', 0))}",
            f"- Tumor probability median / p90 / max: {format_float(score_q['p50'])} / "
            f"{format_float(score_q['p90'])} / {format_float(score_q['max'])}",
            "",
        ])

    if calibration_summary:
        lines.extend([
            "## Calibration",
            "",
            f"- Labeled samples: {calibration_summary['n']} "
            f"({calibration_summary['n_tumor']} tumor / "
            f"{calibration_summary['n_normal']} normal)",
            f"- AUC: {format_float(calibration_summary['auc'])}",
            "- Evaluation: **apparent/resubstitution** (threshold and metrics were "
            "estimated on these same labeled samples; this is not independent validation)",
            f"- Recommended threshold: "
            f"{calibration_summary['recommended_threshold']:.6f} "
            f"({calibration_summary['recommended_metric']})",
            f"- Recommended accuracy / recall / specificity: "
            f"{format_float(calibration_summary['recommended_accuracy'])} / "
            f"{format_float(calibration_summary['recommended_recall'])} / "
            f"{format_float(calibration_summary['recommended_specificity'])}",
            "",
        ])
        for warning in calibration_summary.get("warnings", []):
            lines.append(f"- WARNING: {warning}")
        if calibration_summary.get("warnings"):
            lines.append("")
    elif calibration_error:
        lines.extend([
            "## Calibration",
            "",
            f"- Calibration was not completed: {calibration_error}",
            "",
        ])

    if explanations_rows is not None:
        lines.extend([
            "## Explanations",
            "",
            f"- Explanation rows: {explanations_rows}",
            "- Positive rows push the LR logit toward tumor; negative rows push it toward normal.",
            "",
        ])

    file_rows = [{"file": name, "path": path} for name, path in out_files.items()]
    lines.extend([
        "## Output files",
        "",
        markdown_table(file_rows),
        "",
    ])
    return "\n".join(lines).rstrip("\n") + "\n"


def default_output_dir(input_path):
    path = Path(input_path)
    return path.with_suffix("").name + "_tumor_normal_workflow"


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Run QC, scoring, optional calibration, explanations, and a report."
    )
    parser.add_argument("input", help="expression matrix (samples x genes)")
    parser.add_argument("-o", "--output-dir", default=None,
                        help="directory for workflow outputs")
    parser.add_argument("--labels", help="optional CSV with sample + label columns")
    parser.add_argument("--sample-column", default="sample")
    parser.add_argument("--label-column", default="label")
    parser.add_argument("--min-match-fraction", type=float, default=1.0,
                        help="minimum fraction of scored samples that must have labels (default 1.0)")
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--max-invalid-cell-fraction", type=float, default=0.0,
                        help=("maximum allowed missing, non-numeric, NaN, or infinite cells "
                              "among matched model genes before failing (default 0)"))
    parser.add_argument("--allow-invalid-values", action="store_true",
                        help=("warn instead of failing when matched model-gene cells are "
                              "missing, non-numeric, NaN, or infinite"))
    parser.add_argument("--extra-threshold", type=float, action="append", default=[])
    parser.add_argument("--top-n", type=int, default=10)
    parser.add_argument("--skip-explanations", action="store_true")
    parser.add_argument("--expected-class", choices=["unknown", "mixed", "tumor", "normal"],
                        default="unknown")
    parser.add_argument("--transpose", action="store_true")
    parser.add_argument("--strict-qc", action="store_true",
                        help="exit non-zero on QC WARN or FAIL")
    parser.add_argument("--allow-qc-fail", action="store_true",
                        help="continue even if QC status is FAIL")
    parser.add_argument("--lr-weights", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "deployable_lr_weights.npz"))
    parser.add_argument("--qc-reference", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "model_qc_reference.json"))
    parser.add_argument("--gene-metadata", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "model_gene_metadata.csv"))
    args = parser.parse_args(argv)

    try:
        validate_threshold(args.threshold, "--threshold")
        validate_threshold(args.max_invalid_cell_fraction, "--max-invalid-cell-fraction")
        for threshold in args.extra_threshold:
            validate_threshold(threshold, "--extra-threshold")
    except ValueError as exc:
        parser.error(str(exc))
    if args.top_n < 1 and not args.skip_explanations:
        parser.error("--top-n must be >= 1 unless --skip-explanations is used")

    output_dir = Path(args.output_dir or default_output_dir(args.input))
    paths = {
        "qc_json": output_dir / "qc.json",
        "scores_csv": output_dir / "scores.csv",
        "report_md": output_dir / "workflow_report.md",
        "manifest_json": output_dir / "manifest.json",
    }
    if args.labels:
        paths["thresholds_csv"] = output_dir / "thresholds.csv"
        paths["calibration_json"] = output_dir / "calibration.json"
    if not args.skip_explanations:
        paths["explanations_csv"] = output_dir / "explanations.csv"

    consumed_inputs = {
        "expression_input": args.input,
        "lr_weights": args.lr_weights,
        "qc_reference": args.qc_reference,
    }
    input_roles = {
        "expression input": args.input,
        "LR weights": args.lr_weights,
        "QC reference": args.qc_reference,
    }
    if args.labels:
        consumed_inputs["labels"] = args.labels
        input_roles["labels input"] = args.labels
    if not args.skip_explanations:
        consumed_inputs["gene_metadata"] = args.gene_metadata
        input_roles["gene metadata"] = args.gene_metadata

    try:
        all_managed_paths = {
            f"workflow output {name}": output_dir / name
            for name in MANAGED_OUTPUT_NAMES
        }
        ensure_distinct_paths(all_managed_paths, input_roles)
        input_records = snapshot_inputs(consumed_inputs)
        weights = load_lr_weights(args.lr_weights)
        reference = load_reference(args.qc_reference)
        df = read_matrix(args.input, transpose=args.transpose)
        qc = inspect_dataframe(
            df,
            weights,
            threshold=args.threshold,
            expected_class=args.expected_class,
            reference=reference,
        )
    except ValueError as exc:
        parser.error(str(exc))

    if qc["status"] == "FAIL" and not args.allow_qc_fail:
        outputs = {
            "qc_json": paths["qc_json"].name,
            "report_md": paths["report_md"].name,
            "manifest_json": paths["manifest_json"].name,
        }
        manifest = {
            "status": "stopped_after_qc_fail",
            "input": args.input,
            "qc_status": qc["status"],
        }
        report = build_report(args.input, qc, None, outputs)
        try:
            with staged_generation(output_dir) as stage:
                write_json(qc, stage / outputs["qc_json"], sort_keys=True)
                write_text(stage / outputs["report_md"], report)
                publish_generation(
                    stage, output_dir, manifest, outputs, input_records
                )
        except (ValueError, RuntimeError) as exc:
            parser.error(str(exc))
        print(f"[workflow] QC status=FAIL; wrote {paths['qc_json']}")
        print("[workflow] stopped before scoring; pass --allow-qc-fail to continue")
        return 1

    try:
        scores, n_matched, missing, alignment_report = score_dataframe_lr_weights(
            df,
            weights,
            args.threshold,
            allow_invalid_values=True,
            allow_low_gene_coverage=True,
            return_alignment_report=True,
        )
    except ValueError as exc:
        parser.error(str(exc))
    print_invalid_alignment_summary(alignment_report, sys.stderr)
    alignment_issues = validate_alignment_report(
        alignment_report,
        max_invalid_cell_fraction=args.max_invalid_cell_fraction,
    )
    if alignment_issues and not args.allow_invalid_values:
        outputs = {
            "qc_json": paths["qc_json"].name,
            "report_md": paths["report_md"].name,
            "manifest_json": paths["manifest_json"].name,
        }
        manifest = {
            "status": "stopped_after_invalid_input",
            "input": args.input,
            "qc_status": qc["status"],
            "alignment": {
                "invalid_matched_cells": int(alignment_report["invalid_matched_cells"]),
                "matched_cells": int(alignment_report["matched_cells"]),
                "invalid_matched_fraction": float(alignment_report["invalid_matched_fraction"]),
                "n_genes_with_invalid_values": int(
                    alignment_report["n_genes_with_invalid_values"]
                ),
                "n_samples_with_invalid_values": int(
                    alignment_report["n_samples_with_invalid_values"]
                ),
                "first_genes_with_invalid_values": alignment_report[
                    "first_genes_with_invalid_values"
                ],
                "first_samples_with_invalid_values": alignment_report[
                    "first_samples_with_invalid_values"
                ],
            },
        }
        report = build_report(
            args.input,
            qc,
            None,
            outputs,
            stop_reason="invalid matched expression values",
        )
        try:
            with staged_generation(output_dir) as stage:
                write_json(qc, stage / outputs["qc_json"], sort_keys=True)
                write_text(stage / outputs["report_md"], report)
                publish_generation(
                    stage, output_dir, manifest, outputs, input_records
                )
        except (ValueError, RuntimeError) as exc:
            parser.error(str(exc))
        for issue in alignment_issues:
            print(f"[workflow] ERROR: {issue}", file=sys.stderr)
        print(
            "[workflow] stopped before writing scores; fix invalid values or "
            "pass --allow-invalid-values",
            file=sys.stderr,
        )
        return 1
    if alignment_issues:
        for issue in alignment_issues:
            print(f"[workflow] WARNING: {issue}", file=sys.stderr)
    explanations = None
    explanations_rows = None
    if not args.skip_explanations:
        try:
            gene_names = load_gene_metadata(args.gene_metadata)
            explanations, _, _ = explain_dataframe(df, weights, args.top_n, gene_names)
            explanations_rows = int(len(explanations))
        except ValueError as exc:
            parser.error(str(exc))

    calibration_summary = None
    calibration_error = None
    try:
        with staged_generation(output_dir) as stage:
            write_json(qc, stage / paths["qc_json"].name, sort_keys=True)
            write_dataframe_csv(scores, stage / paths["scores_csv"].name, index=False)

            if args.labels:
                try:
                    threshold_metrics, calibration_summary = calibration_from_labels(
                        stage / paths["scores_csv"].name,
                        args.labels,
                        args.sample_column,
                        args.label_column,
                        args.threshold,
                        args.extra_threshold,
                        args.min_match_fraction,
                    )
                except ValueError as exc:
                    calibration_error = str(exc)

            if calibration_error is not None:
                outputs = {
                    "qc_json": paths["qc_json"].name,
                    "scores_csv": paths["scores_csv"].name,
                    "report_md": paths["report_md"].name,
                    "manifest_json": paths["manifest_json"].name,
                }
                manifest = {
                    "status": "stopped_after_calibration_error",
                    "input": args.input,
                    "labels": args.labels,
                    "qc_status": qc["status"],
                    "calibration_error": calibration_error,
                }
                report = build_report(
                    args.input,
                    qc,
                    scores,
                    outputs,
                    calibration_error=calibration_error,
                )
                write_text(stage / outputs["report_md"], report)
                manifest = publish_generation(
                    stage, output_dir, manifest, outputs, input_records
                )
            else:
                if args.labels:
                    write_dataframe_csv(
                        threshold_metrics,
                        stage / paths["thresholds_csv"].name,
                        index=False,
                    )
                    write_json(
                        calibration_summary,
                        stage / paths["calibration_json"].name,
                        sort_keys=True,
                    )

                if explanations is not None:
                    write_dataframe_csv(
                        explanations,
                        stage / paths["explanations_csv"].name,
                        index=False,
                    )

                out_files = {name: path.name for name, path in paths.items()}
                report = build_report(
                    args.input,
                    qc,
                    scores,
                    out_files,
                    calibration_summary,
                    explanations_rows,
                )
                manifest = {
                    "status": "complete",
                    "input": args.input,
                    "threshold": args.threshold,
                    "samples": int(df.shape[0]),
                    "input_genes": int(df.shape[1]),
                    "matched_model_genes": int(n_matched),
                    "missing_model_genes": int(len(missing)),
                    "qc_status": qc["status"],
                    "tumor_calls": int((scores["call"] == "tumor").sum()),
                    "normal_calls": int((scores["call"] == "normal").sum()),
                    "labels": args.labels,
                    "calibration": calibration_summary,
                }
                write_text(stage / paths["report_md"].name, report)
                manifest = publish_generation(
                    stage, output_dir, manifest, out_files, input_records
                )
    except (ValueError, RuntimeError) as exc:
        parser.error(str(exc))

    if calibration_error is not None:
        print(
            f"[workflow] ERROR: calibration failed: {calibration_error}",
            file=sys.stderr,
        )
        print(
            "[workflow] stopped after writing scores; fix labels and rerun",
            file=sys.stderr,
        )
        return 1

    print(f"[workflow] QC status={qc['status']}; matched {n_matched}/{len(weights['selected_genes'])} model genes")
    print(f"[workflow] wrote {paths['scores_csv']} ({manifest['tumor_calls']} tumor / {manifest['normal_calls']} normal)")
    if calibration_summary:
        print(f"[workflow] recommended threshold={calibration_summary['recommended_threshold']:.6f}")
        print(
            "[workflow] NOTE: calibration metrics are apparent/resubstitution estimates",
            file=sys.stderr,
        )
        for warning in calibration_summary.get("warnings", []):
            print(f"[workflow] WARNING: {warning}", file=sys.stderr)
    if explanations_rows is not None:
        print(f"[workflow] wrote {paths['explanations_csv']} ({explanations_rows} rows)")
    print(f"[workflow] wrote {paths['report_md']}")

    if args.strict_qc and qc["status"] != "PASS":
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
