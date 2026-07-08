#!/usr/bin/env python3
"""
predict_cancer_type.py -- score expression samples for TCGA cancer type
(tissue of origin) with the deployable multinomial logistic-regression model.

Pipeline (via the tcga_rnaseq shared core, pure numpy):
  z = (x_selected - scaler_mean) / scaler_scale
  logits = z @ coef.T + intercept        # coef: (17, k), one row per cancer type
  p = softmax(logits)                    # per-type probability
  call = argmax type

Input CSV: rows = samples, columns = Ensembl gene IDs, values = GDC STAR-Counts
log2(TPM+1). Ensembl version suffixes are matched with or without ".N"; missing
model genes are imputed at the training mean. Requires only numpy and pandas.
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))  # deliverables root, for tcga_rnaseq
from tcga_rnaseq import load_lr_model, read_matrix, predict_proba  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description="Predict TCGA cancer type (tissue of origin).")
    ap.add_argument("input_csv", help="rows=samples, cols=Ensembl gene IDs, values=log2(TPM+1)")
    ap.add_argument("--weights", default=os.path.join(HERE, "cancer_type_lr_weights.npz"))
    ap.add_argument("--out", default=None)
    ap.add_argument("--topk", type=int, default=3, help="report top-k types per sample")
    args = ap.parse_args(argv)

    model = load_lr_model(args.weights)
    if model["kind"] != "multiclass":
        ap.error("weights file is not a multi-class cancer-type model")
    if args.topk < 1:
        ap.error("--topk must be >= 1")
    X = read_matrix(args.input_csv)
    P = predict_proba(model, X)
    classes = model["classes"]
    if args.topk > len(classes):
        ap.error(f"--topk must be <= number of classes ({len(classes)})")
    order = np.argsort(-P, axis=1)

    rows = []
    import pandas as pd
    for i, s in enumerate(X.index.astype(str)):
        top = order[i, :args.topk]
        row = {"sample": s, "predicted_type": classes[top[0]],
               "probability": round(float(P[i, top[0]]), 4)}
        for r, ci in enumerate(top, 1):
            row[f"top{r}"] = f"{classes[ci]}:{P[i, ci]:.3f}"
        rows.append(row)
    out = pd.DataFrame(rows)
    out_path = args.out or (os.path.splitext(args.input_csv)[0] + ".cancer_type_pred.csv")
    out.to_csv(out_path, index=False)
    print(json.dumps({"n_samples": int(X.shape[0]), "n_types": len(classes),
                      "scores_csv": out_path,
                      "call_distribution": out["predicted_type"].value_counts().to_dict()}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
