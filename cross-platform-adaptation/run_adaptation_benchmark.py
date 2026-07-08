#!/usr/bin/env python
"""
Reproduce the cross-platform domain-adaptation benchmark.

Scores the three external validation cohorts (GTEx/Toil normals, TCGA-Toil/RSEM,
CPTAC-3/GDC) with the frozen release model under three scoring modes and writes:
  adaptation_benchmark.csv   per cohort x mode: AUC / acc@0.5 / balanced-acc / FPR
  adaptation_imbalance.csv   Toil: balanced-acc vs cohort tumor fraction

Uses the tcga_rnaseq shared core. Run from this folder; it reads the sibling
external-validation/ pkls and the parent deployable_lr_weights.npz.
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from tcga_rnaseq import load_lr_model, read_matrix, predict_proba  # noqa: E402
from tcga_rnaseq import metrics as M  # noqa: E402
from tcga_rnaseq.align import align_to_genes  # noqa: E402

EV = os.path.join(ROOT, "external-validation")
MODES = ["none", "cohort_center", "cohort_zscore"]
MODE_LABEL = {"none": "deployed (no adaptation)",
              "cohort_center": "cohort mean-centering",
              "cohort_zscore": "cohort standardization (z-score)"}
COHORTS = [
    ("TCGA-Toil/RSEM (foreign pipeline)",
     "tcga_toil_xena/tcga_toil_selected_genes_model_scale.pkl",
     "tcga_toil_xena/tcga_toil_predictions.csv", "sample", False),
    ("CPTAC-3/GDC STAR-Counts (native pipeline)",
     "cptac_gdc/expression_selected_genes.pkl",
     "cptac_gdc/cptac_predictions.csv", "file_id", False),
    ("GTEx/Toil normals (foreign, all-normal)",
     "gtex_xena/gtex_selected_genes_model_scale.pkl", None, None, True),
]


def _metrics(y, p):
    cm = M.classification_metrics(y, p, 0.5)
    auc = M.roc_auc(y, p)
    return {"auc": auc, "acc": cm["accuracy"], "bacc": cm["balanced_accuracy"],
            "fpr": (cm["fp"] / (cm["fp"] + cm["tn"])) if (cm["fp"] + cm["tn"]) else float("nan")}


def _labels(model, X, lp, key, all_normal):
    if all_normal:
        return np.zeros(len(X), dtype=int)
    lab = pd.read_csv(os.path.join(EV, lp))
    return lab.set_index(lab[key].astype(str))["label"].reindex(X.index.astype(str)).astype(int).values


def benchmark(model):
    rows = []
    for name, mp, lp, key, all_normal in COHORTS:
        X = read_matrix(os.path.join(EV, mp))
        y = _labels(model, X, lp, key, all_normal)
        for mode in MODES:
            p = predict_proba(model, X, adapt=mode)
            m = _metrics(y, p)
            has_two = len(set(y)) > 1
            rows.append(dict(cohort=name, n=len(X), mode=MODE_LABEL[mode],
                             AUC=round(m["auc"], 4) if has_two else "",
                             acc_at_0p5=round(m["acc"], 4),
                             balanced_acc=round(m["bacc"], 4) if has_two else "",
                             FPR_at_0p5=round(m["fpr"], 4) if m["fpr"] == m["fpr"] else ""))
            print(f"{name[:38]:38s} | {MODE_LABEL[mode]:34s} | AUC {rows[-1]['AUC']!s:>6} "
                  f"acc {rows[-1]['acc_at_0p5']:.3f} bacc {rows[-1]['balanced_acc']!s:>6} "
                  f"FPR {rows[-1]['FPR_at_0p5']!s:>6}")
    return pd.DataFrame(rows)


def imbalance_curve(model):
    X = read_matrix(os.path.join(EV, "tcga_toil_xena/tcga_toil_selected_genes_model_scale.pkl"))
    lab = pd.read_csv(os.path.join(EV, "tcga_toil_xena/tcga_toil_predictions.csv"))
    y = lab.set_index(lab["sample"].astype(str))["label"].reindex(X.index.astype(str)).astype(int).values
    V, _, _ = align_to_genes(X, model["genes"], impute_mean=model["mean"])
    Vt, Vn = V[y == 1], V[y == 0]
    rng = np.random.default_rng(0)
    irows = []
    for frac in [0.10, 0.25, 0.50, 0.75, 0.90]:
        b_base, b_da, auc_da = [], [], []
        for _ in range(200):
            n = 80
            nt = max(1, min(n - 1, int(round(frac * n))))
            it = rng.choice(len(Vt), nt, replace=True)
            ino = rng.choice(len(Vn), n - nt, replace=True)
            Xs = pd.DataFrame(np.vstack([Vt[it], Vn[ino]]), columns=model["genes"])
            ys = np.r_[np.ones(nt), np.zeros(n - nt)].astype(int)
            pb = predict_proba(model, Xs, adapt="none")
            pdz = predict_proba(model, Xs, adapt="cohort_zscore")
            b_base.append(_metrics(ys, pb)["bacc"])
            b_da.append(_metrics(ys, pdz)["bacc"])
            auc_da.append(M.roc_auc(ys, pdz))
        irows.append(dict(tumor_fraction=frac, cohort_n=80,
                          balanced_acc_deployed=round(np.mean(b_base), 4),
                          balanced_acc_adapted=round(np.mean(b_da), 4),
                          AUC_adapted=round(np.mean(auc_da), 4)))
        print(f"frac {frac:.2f}  bacc deployed {irows[-1]['balanced_acc_deployed']:.3f} "
              f"-> adapted {irows[-1]['balanced_acc_adapted']:.3f}  (AUC {irows[-1]['AUC_adapted']:.3f})")
    return pd.DataFrame(irows)


def main():
    model = load_lr_model(os.path.join(ROOT, "deployable_lr_weights.npz"))
    benchmark(model).to_csv(os.path.join(HERE, "adaptation_benchmark.csv"), index=False)
    imbalance_curve(model).to_csv(os.path.join(HERE, "adaptation_imbalance.csv"), index=False)
    print("\nwrote adaptation_benchmark.csv and adaptation_imbalance.csv")


if __name__ == "__main__":
    main()
