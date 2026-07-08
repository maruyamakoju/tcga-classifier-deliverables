#!/usr/bin/env python3
"""Regenerate model_qc_reference.json from bundled validation matrices."""
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from inspect_expression_input import DEFAULT_RULES, inspect_dataframe
from score_tumor_normal import load_lr_weights, read_matrix


ROOT = Path(__file__).resolve().parent


def load_matrix(path):
    path = ROOT / path
    if not path.exists():
        return None
    if path.suffix == ".pkl":
        return pd.read_pickle(path)
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


def main():
    weights = load_lr_weights(ROOT / "deployable_lr_weights.npz")
    sources = {
        "example_tcga_gdc": "example_input.csv",
        "cptac3_gdc_star_counts": "external-validation/cptac_gdc/expression_selected_genes.pkl",
        "tcga_toil_rsem_pipeline_check": "external-validation/tcga_toil_xena/tcga_toil_selected_genes_model_scale.pkl",
        "gtex_toil_normals_check": "external-validation/gtex_xena/gtex_selected_genes_model_scale.pkl",
    }
    reports = {}
    compatible_frames = []
    for name, rel_path in sources.items():
        df = load_matrix(rel_path)
        if df is None:
            continue
        report = inspect_dataframe(df, weights, reference={"rules": DEFAULT_RULES})
        reports[name] = compact_metrics(report)
        if name in {"example_tcga_gdc", "cptac3_gdc_star_counts"}:
            compatible_frames.append(df)

    if compatible_frames:
        combined = pd.concat(compatible_frames, axis=0)
        reports["compatible_gdc_combined"] = compact_metrics(
            inspect_dataframe(combined, weights, reference={"rules": DEFAULT_RULES})
        )

    reference = {
        "schema_version": "1.0",
        "reference_source": (
            "Rules use the model training scaler plus empirical summaries from the bundled "
            "TCGA example and CPTAC-3/GDC STAR-Counts smoke validation. Toil/RSEM/GTEx "
            "summaries are included as known cross-platform boundary checks."
        ),
        "intended_input": "GDC STAR-Counts-style log2(TPM+1), rows=samples, columns=Ensembl genes.",
        "rules": DEFAULT_RULES,
        "reference_reports": reports,
    }
    out = ROOT / "model_qc_reference.json"
    with open(out, "w", encoding="utf-8") as handle:
        json.dump(reference, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
