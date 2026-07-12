#!/usr/bin/env python3
"""Regenerate model_qc_reference.json from reviewed historical matrices."""
import argparse
from pathlib import Path

import pandas as pd

from inspect_expression_input import DEFAULT_RULES, inspect_dataframe
from score_tumor_normal import load_lr_weights, read_matrix
from tcga_rnaseq import ensure_distinct_paths, write_json


ROOT = Path(__file__).resolve().parent


def load_matrix(path, *, trusted_historical_pickles=False):
    path = Path(path)
    if not path.is_absolute():
        path = ROOT / path
    if not path.exists():
        return None
    if path.suffix == ".pkl":
        if not trusted_historical_pickles:
            raise ValueError(
                "historical pickle inputs require --trusted-historical-pickles; "
                "never use this flag for collaborator-controlled files"
            )
        return read_matrix(str(path), allow_pickle=True)
    return read_matrix(str(path))


def compact_metrics(report):
    return {
        "shape": report["shape"],
        "gene_match": {
            "matched_model_genes": report["gene_match"]["matched_model_genes"],
            "match_rate": report["gene_match"]["match_rate"],
        },
        "value_summary": report["value_summary"],
        "distribution_summary": report["distribution_summary"],
        "score_summary": report["score_summary"],
        "status_under_current_rules": report["status"],
        "messages_under_current_rules": report["messages"],
    }


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Regenerate the QC reference from the fixed reviewed historical matrices."
        )
    )
    parser.add_argument("--weights", type=Path, default=ROOT / "deployable_lr_weights.npz")
    parser.add_argument("--output", type=Path, default=ROOT / "model_qc_reference.json")
    parser.add_argument(
        "--trusted-historical-pickles",
        action="store_true",
        help="acknowledge that the fixed repository-local pickle inputs were reviewed",
    )
    parser.add_argument(
        "--allow-missing-sources",
        action="store_true",
        help="development only: generate an explicitly incomplete reference",
    )
    parser.add_argument(
        "--dry-run", action="store_true", help="validate and build in memory without writing"
    )
    args = parser.parse_args(argv)

    sources = {
        "example_tcga_gdc": "example_input.csv",
        "cptac3_gdc_star_counts": "external-validation/cptac_gdc/expression_selected_genes.pkl",
        "tcga_toil_rsem_pipeline_check": "external-validation/tcga_toil_xena/tcga_toil_selected_genes_model_scale.pkl",
        "gtex_toil_normals_check": "external-validation/gtex_xena/gtex_selected_genes_model_scale.pkl",
    }
    source_paths = {name: ROOT / path for name, path in sources.items()}
    try:
        ensure_distinct_paths(
            {"QC reference output": args.output},
            {"weights input": args.weights, **source_paths},
        )
        weights = load_lr_weights(args.weights)
        reports = {}
        compatible_frames = []
        missing_sources = []
        for name, source_path in source_paths.items():
            df = load_matrix(
                source_path,
                trusted_historical_pickles=args.trusted_historical_pickles,
            )
            if df is None:
                missing_sources.append(source_path)
                continue
            report = inspect_dataframe(df, weights, reference={"rules": DEFAULT_RULES})
            reports[name] = compact_metrics(report)
            if name in {"example_tcga_gdc", "cptac3_gdc_star_counts"}:
                compatible_frames.append(df)
        if missing_sources and not args.allow_missing_sources:
            raise ValueError(
                "required historical QC sources are missing: "
                + ", ".join(str(path) for path in missing_sources)
                + "; use --allow-missing-sources only for an intentionally incomplete "
                "development reference"
            )

        if compatible_frames:
            combined = pd.concat(compatible_frames, axis=0)
            reports["compatible_gdc_combined"] = compact_metrics(
                inspect_dataframe(combined, weights, reference={"rules": DEFAULT_RULES})
            )

        reference = {
            "schema_version": "1.0",
            "reference_source": (
                "Rules use the model training scaler plus empirical summaries from the "
                "reviewed historical TCGA example and CPTAC-3/GDC STAR-Counts smoke "
                "validation. Toil/RSEM/GTEx summaries are historical cross-platform "
                "boundary checks, not a post-fix live rerun."
            ),
            "intended_input": (
                "GDC STAR-Counts-style log2(TPM+1), rows=samples, "
                "columns=Ensembl genes."
            ),
            "rules": DEFAULT_RULES,
            "reference_reports": reports,
        }
        if not args.dry_run:
            write_json(reference, args.output, sort_keys=True)
    except (ValueError, OSError, TypeError) as exc:
        parser.error(str(exc))

    action = "validated without writing" if args.dry_run else f"wrote {args.output}"
    print(f"[qc-reference] {action}; {len(reports)} reference reports")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
