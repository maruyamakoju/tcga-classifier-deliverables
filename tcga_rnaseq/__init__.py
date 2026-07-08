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
from .io import (load_lr_model, load_pipeline, read_matrix,
                 read_expression_csv, write_json)
from .align import align_to_genes, align_to_genes_with_report, strip_version
from .score import (
    predict_proba,
    predict_proba_from_aligned,
    score_binary_dataframe,
    predict,
    sigmoid,
    softmax,
    standardize,
    ADAPT_MODES,
)

__all__ = [
    "load_lr_model", "load_pipeline", "read_matrix",
    "read_expression_csv", "write_json",
    "align_to_genes", "align_to_genes_with_report", "strip_version",
    "predict_proba", "predict_proba_from_aligned",
    "score_binary_dataframe", "predict", "sigmoid", "softmax",
    "standardize", "ADAPT_MODES",
]

__version__ = "2.0.0"
