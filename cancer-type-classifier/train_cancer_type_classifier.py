#!/usr/bin/env python
"""
Train and patient-held-out-evaluate a 17-class TCGA cancer-type (tissue-of-origin)
classifier on tumor samples, then export a pure-numpy deployable model.

Model: StandardScaler -> SelectKBest(f_classif, k=1000) -> multinomial LogisticRegression(C=2).
Evaluation: 5-fold StratifiedGroupKFold grouped by patient (case_id); out-of-fold
predictions pooled for the confusion matrix and per-class metrics (no patient
appears in both train and test).

Inputs (all version-neutral):
  --features  X_full.npy  (n_samples x 14850 float32, GDC STAR-Counts log2(TPM+1))
              regenerate from X_full_filtered.pkl with export_features_npy.py
  X_genes.npy, X_samples.npy   (sibling of --features)
  ../selected_files.csv         (file_id -> project [cancer type], label, case_id)

Outputs (this folder):
  cancer_type_oof_predictions.csv   per-sample true/pred/confidence (patient-held-out)
  cancer_type_per_class_metrics.csv precision/recall/F1/support per cancer type
  cancer_type_confusion_matrix.csv  rows=true, cols=pred
  cancer_type_summary.json          headline metrics + config
  cancer_type_top_genes.csv         top marker genes per type (final model)
  cancer_type_lr_weights.npz        pure-numpy deployable model
"""
import argparse, json, os, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import (accuracy_score, balanced_accuracy_score, f1_score,
                             confusion_matrix, precision_recall_fscore_support)

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
K, C, NFOLDS, SEED = 1000, 2.0, 5, 0


def load(features):
    X = np.load(features)
    d = os.path.dirname(features)
    genes = np.load(os.path.join(d, "X_genes.npy"), allow_pickle=True).astype(str)
    samp = np.load(os.path.join(d, "X_samples.npy"), allow_pickle=True).astype(str)
    sf = pd.read_csv(os.path.join(ROOT, "selected_files.csv")).set_index("file_id")
    meta = sf.reindex(samp)
    ctype = meta["project"].str.replace("TCGA-", "", regex=False).values
    grp = meta["case_id"].astype(str).values
    tumor = meta["label"].values == "tumor"
    return X[tumor], ctype[tumor], grp[tumor], samp[tumor], genes


def make_pipe():
    return Pipeline([("sc", StandardScaler()),
                     ("sel", SelectKBest(f_classif, k=K)),
                     ("lr", LogisticRegression(C=C, max_iter=3000))])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--features", default=os.path.join(HERE, "X_full.npy"))
    args = ap.parse_args()

    X, y, g, samp, genes = load(args.features)
    classes = sorted(set(y))
    print(f"tumors={len(y)} classes={len(classes)} patients={len(set(g))} genes={X.shape[1]}")

    # ---- patient-held-out out-of-fold predictions ----
    cv = StratifiedGroupKFold(n_splits=NFOLDS, shuffle=True, random_state=SEED)
    oof = np.empty(len(y), dtype=object)
    oof_conf = np.zeros(len(y))
    for tr, te in cv.split(X, y, g):
        p = make_pipe().fit(X[tr], y[tr])
        proba = p.predict_proba(X[te])
        oof[te] = p.classes_[proba.argmax(1)]
        oof_conf[te] = proba.max(1)

    acc = accuracy_score(y, oof)
    bacc = balanced_accuracy_score(y, oof)
    mf1 = f1_score(y, oof, average="macro")
    wf1 = f1_score(y, oof, average="weighted")
    print(f"acc={acc:.4f} balanced_acc={bacc:.4f} macroF1={mf1:.4f} weightedF1={wf1:.4f}")

    pd.DataFrame({"file_id": samp, "true": y, "pred": oof,
                  "correct": (oof == y).astype(int), "confidence": oof_conf.round(4)}
                 ).to_csv(os.path.join(HERE, "cancer_type_oof_predictions.csv"), index=False)

    pr, rc, f1s, sup = precision_recall_fscore_support(y, oof, labels=classes, zero_division=0)
    pd.DataFrame({"cancer_type": classes, "support": sup,
                  "precision": pr.round(4), "recall": rc.round(4), "f1": f1s.round(4)}
                 ).sort_values("f1", ascending=False
                 ).to_csv(os.path.join(HERE, "cancer_type_per_class_metrics.csv"), index=False)

    cm = confusion_matrix(y, oof, labels=classes)
    pd.DataFrame(cm, index=classes, columns=classes).to_csv(
        os.path.join(HERE, "cancer_type_confusion_matrix.csv"))

    # ---- final deployable model on all tumors ----
    sc = StandardScaler().fit(X)
    Xs = sc.transform(X)
    sel = SelectKBest(f_classif, k=K).fit(Xs, y)
    mask = sel.get_support()
    lr = LogisticRegression(C=C, max_iter=5000).fit(Xs[:, mask], y)

    sel_idx = np.where(mask)[0]
    np.savez(os.path.join(HERE, "cancer_type_lr_weights.npz"),
             selected_genes=genes[sel_idx].astype(str),
             selected_gene_index=sel_idx.astype(np.int32),
             scaler_mean=sc.mean_[sel_idx], scaler_scale=sc.scale_[sel_idx],
             coef=lr.coef_, intercept=lr.intercept_, classes=np.array(lr.classes_, dtype=str),
             notes="Pure-numpy multinomial LR cancer-type classifier. Input: GDC "
                   "STAR-Counts log2(TPM+1), columns=Ensembl gene IDs. "
                   "score: softmax(coef @ ((x_sel-scaler_mean)/scaler_scale) + intercept).")

    # ---- top marker genes per type (largest positive standardized coefficient) ----
    sg = genes[sel_idx]
    sym = {}
    smap = os.path.join(HERE, "gene_id_to_name.csv")
    if os.path.exists(smap):
        sym = pd.read_csv(smap, index_col=0)["symbol"].astype(str).to_dict()
    rows = []
    for ci, cname in enumerate(lr.classes_):
        order = np.argsort(lr.coef_[ci])[::-1][:15]
        for rank, gi in enumerate(order, 1):
            rows.append({"cancer_type": cname, "rank": rank, "gene": sg[gi],
                         "symbol": sym.get(sg[gi], ""), "coef": round(float(lr.coef_[ci][gi]), 4)})
    pd.DataFrame(rows).to_csv(os.path.join(HERE, "cancer_type_top_genes.csv"), index=False)

    json.dump({"n_tumors": int(len(y)), "n_classes": len(classes), "n_patients": int(len(set(g))),
               "n_input_genes": int(X.shape[1]), "selected_genes": K,
               "model": f"StandardScaler+SelectKBest(f_classif,k={K})+LogisticRegression(C={C})",
               "evaluation": f"{NFOLDS}-fold StratifiedGroupKFold by case_id (patient-held-out)",
               "accuracy": round(acc, 4), "balanced_accuracy": round(bacc, 4),
               "macro_f1": round(mf1, 4), "weighted_f1": round(wf1, 4)},
              open(os.path.join(HERE, "cancer_type_summary.json"), "w"), indent=2)
    print("wrote outputs to", HERE)


if __name__ == "__main__":
    main()
