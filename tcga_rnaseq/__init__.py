"""
tcga_rnaseq -- shared core for the TCGA/GDC RNA-seq classifier deliverables.

A single, dependency-light (numpy + pandas) implementation of the logic that was
previously copy-pasted across the scoring CLIs: model I/O, gene alignment,
standardization (incl. cross-platform domain adaptation), logistic scoring
(binary tumor-vs-normal and multi-class cancer-type), and metrics.

Public API:
    from tcga_rnaseq import load_lr_model, read_expression_csv, predict_proba
    from tcga_rnaseq.metrics import roc_auc, accuracy, balanced_accuracy
"""
from .io import (
    ensure_distinct_paths,
    load_lr_model,
    load_pipeline,
    read_csv_table,
    read_expression_csv,
    read_matrix,
    validate_lr_model,
    write_dataframe_csv,
    write_json,
    write_text,
)
from .align import (
    align_to_genes,
    align_to_genes_with_report,
    format_alignment_issues,
    format_gene_match_issues,
    print_invalid_alignment_summary,
    strip_version,
    validate_alignment_report,
    validate_gene_match_report,
)
from .score import (
    predict_proba,
    predict_proba_from_aligned,
    score_binary_dataframe,
    predict,
    sigmoid,
    softmax,
    standardize,
    validate_tumor_binary_model,
    ADAPT_MODES,
)
from .validation import (
    normalize_label,
    require_unique_samples,
    sample_key,
    validate_expression_matrix,
    validate_threshold,
)

__all__ = [
    "load_lr_model", "validate_lr_model", "load_pipeline", "read_matrix",
    "read_csv_table", "read_expression_csv", "ensure_distinct_paths",
    "write_dataframe_csv", "write_json", "write_text",
    "align_to_genes", "align_to_genes_with_report", "strip_version",
    "validate_alignment_report", "format_alignment_issues",
    "validate_gene_match_report", "format_gene_match_issues",
    "print_invalid_alignment_summary",
    "predict_proba", "predict_proba_from_aligned",
    "score_binary_dataframe", "predict", "sigmoid", "softmax",
    "validate_tumor_binary_model",
    "standardize", "ADAPT_MODES",
    "validate_threshold", "normalize_label", "sample_key", "require_unique_samples",
    "validate_expression_matrix",
]

# tcga_rnaseq's own library version, independent of the release VERSION file
# (currently v2.0.0-gdc-starcounts) -- this bumps only on changes to this
# package's public API/behavior, not on every release.
__version__ = "3.0.0"
