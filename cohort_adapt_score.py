#!/usr/bin/env python3
"""
cohort_adapt_score.py -- cross-platform (domain-adapted) scoring for the
TCGA/GDC tumor-vs-normal logistic-regression release.

The deployed model is calibrated to GDC STAR-Counts log2(TPM+1). On foreign
RNA-seq pipelines (e.g. UCSC Xena Toil/RSEM, GTEx/Toil) discrimination (AUC)
transfers but the 0.5 decision threshold does not: probabilities saturate and
almost every sample is called tumor.

This tool applies a label-free, no-retraining domain adaptation: it standardizes
each gene using the *input cohort's own* per-gene mean and standard deviation,
then applies the frozen logistic-regression coefficients. This realigns a
foreign cohort's per-gene marginal distribution onto the training marginal and
restores the default 0.5 threshold, provided the cohort contains an internal
tumor/normal contrast (see the important limitation below).

Adaptation modes (from tcga_rnaseq.score.standardize):
  none            deployed scoring: z = (x - train_mean)/train_scale
  cohort_zscore   z = (x - cohort_mean)/cohort_std               (recommended)
  cohort_center   z = (x - cohort_mean)/train_scale              (location-only)

IMPORTANT LIMITATION: cohort standardization assumes the input cohort has an
internal mix of tumor and normal. A near-single-class cohort (e.g. an all-normal
QC panel) has no internal contrast to anchor the recentering and is only
partially corrected -- prefer an explicit labeled-anchor recalibration
(calibrate_threshold.py) for such cohorts.

Requires only numpy and pandas (via the tcga_rnaseq shared core).
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tcga_rnaseq import load_lr_model, read_matrix, predict_proba, write_json  # noqa: E402
from tcga_rnaseq import metrics as M  # noqa: E402
from tcga_rnaseq.score import ADAPT_MODES  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def main(argv=None):
    ap = argparse.ArgumentParser(description="Cross-platform domain-adapted tumor-vs-normal scoring.")
    ap.add_argument("input_csv", help="expression matrix CSV: rows=samples, cols=Ensembl gene IDs, values=log2(TPM+1)")
    ap.add_argument("--adapt", default="cohort_zscore", choices=list(ADAPT_MODES),
                    help="domain-adaptation mode (default: cohort_zscore)")
    ap.add_argument("--threshold", type=float, default=0.5)
    ap.add_argument("--labels", default=None,
                    help="optional CSV with columns sample,label (label in {0,1} or {tumor,normal}) to report metrics")
    ap.add_argument("--weights", default=os.path.join(HERE, "deployable_lr_weights.npz"))
    ap.add_argument("--out", default=None)
    ap.add_argument("--min-samples", type=int, default=20)
    args = ap.parse_args(argv)

    if not 0.0 <= args.threshold <= 1.0:
        ap.error("--threshold must be between 0 and 1")

    model = load_lr_model(args.weights)
    X = read_matrix(args.input_csv)
    n = X.shape[0]

    warnings = []
    if args.adapt != "none" and n < args.min_samples:
        warnings.append(f"cohort has only {n} samples (< {args.min_samples}); "
                        f"cohort statistics may be unreliable")
    if args.adapt != "none":
        warnings.append("cohort standardization assumes an internal tumor/normal mix; "
                        "a near-single-class cohort is only partially corrected")

    p = predict_proba(model, X, adapt=args.adapt)
    calls = np.where(p >= args.threshold, "tumor", "normal")
    out_df = pd.DataFrame({"sample": X.index.astype(str),
                           "tumor_probability": p.round(6), "call": calls})

    metrics = None
    if args.labels:
        lab = pd.read_csv(args.labels)
        key = lab.columns[0]
        s = lab.set_index(lab[key].astype(str))["label"]
        y = s.reindex(X.index.astype(str))
        if y.dtype == object:
            y = (y.astype(str).str.lower() == "tumor").astype(float)
        y = y.values
        m = ~np.isnan(y)
        if m.sum() > 0 and len(set(y[m])) > 1:
            cm = M.classification_metrics(y[m].astype(int), p[m], args.threshold)
            metrics = {"n_labeled": int(m.sum()), "auc": round(cm["auc"], 4),
                       "accuracy": round(cm["accuracy"], 4),
                       "balanced_accuracy": round(cm["balanced_accuracy"], 4),
                       "sensitivity": round(cm["sensitivity"], 4),
                       "specificity": round(cm["specificity"], 4)}

    out_path = args.out or (os.path.splitext(args.input_csv)[0] + ".adapted_scores.csv")
    out_df.to_csv(out_path, index=False)

    report = {"n_samples": int(n), "adapt": args.adapt, "threshold": args.threshold,
              "tumor_calls": int((p >= args.threshold).sum()),
              "normal_calls": int((p < args.threshold).sum()),
              "median_tumor_probability": float(np.median(p)),
              "scores_csv": out_path, "warnings": warnings, "metrics": metrics}
    print(__import__("json").dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
