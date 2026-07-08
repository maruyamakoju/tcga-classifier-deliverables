#!/usr/bin/env python3
"""Train and evaluate the TCGA tumor-vs-normal classifiers (provenance script).

Reproduces the shipped models: univariate ANOVA feature selection (top 2,000
genes, fit on TRAIN only), then Logistic Regression / Random Forest / XGBoost,
evaluated on a patient-disjoint held-out test set plus grouped 5-fold CV.

Requires the training matrices (X_full_filtered.pkl etc.), which were pickled
with numpy>=2 / pandas>=3; run in that environment. XGBoost is optional (it
segfaults in some conda envs) and is skipped with a note if unavailable.

    python train_classifier.py
"""
import os
import pickle
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))


def _p(name):
    return os.path.join(HERE, name)


def main():
    t0 = time.time()
    X_full = pd.read_pickle(_p("X_full_filtered.pkl"))
    y_full = pd.read_pickle(_p("y_full.pkl"))
    groups_full = pd.read_pickle(_p("groups_full.pkl"))
    projects_full = pd.read_pickle(_p("projects_full.pkl"))
    train_idx = np.load(_p("train_idx.npy"))
    test_idx = np.load(_p("test_idx.npy"))

    X_train, X_test = X_full.iloc[train_idx], X_full.iloc[test_idx]
    y_train, y_test = y_full.iloc[train_idx], y_full.iloc[test_idx]
    proj_test = projects_full.iloc[test_idx]
    print(f"[{time.time()-t0:.0f}s] loaded: train={X_train.shape} test={X_test.shape}", flush=True)

    from sklearn.feature_selection import SelectKBest, f_classif
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import LogisticRegression
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.pipeline import Pipeline
    from sklearn.model_selection import StratifiedGroupKFold
    from sklearn.metrics import (roc_auc_score, accuracy_score, f1_score, precision_score,
                                 recall_score, confusion_matrix, roc_curve)
    try:
        import xgboost as xgb
        have_xgb = True
    except Exception as exc:  # xgboost import segfaults / is absent in some envs
        have_xgb = False
        print(f"[warn] xgboost unavailable ({exc}); skipping XGB model", flush=True)

    # --- Feature selection: univariate ANOVA F-test on TRAIN only, top 2000 genes ---
    selector = SelectKBest(f_classif, k=2000)
    X_train_sel = selector.fit_transform(X_train, y_train)
    X_test_sel = selector.transform(X_test)
    selected_genes = X_train.columns[selector.get_support()]
    print(f"[{time.time()-t0:.0f}s] feature selection done: {X_train_sel.shape}", flush=True)

    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train_sel)
    X_test_scaled = scaler.transform(X_test_sel)

    results = {}
    lr = LogisticRegression(max_iter=5000, C=0.1, class_weight="balanced", random_state=42)
    lr.fit(X_train_scaled, y_train)
    results["logistic_regression"] = lr
    print(f"[{time.time()-t0:.0f}s] LR trained", flush=True)

    rf = RandomForestClassifier(n_estimators=500, max_depth=8, min_samples_leaf=3,
                                class_weight="balanced", n_jobs=-1, random_state=42)
    rf.fit(X_train_sel, y_train)
    results["random_forest"] = rf
    print(f"[{time.time()-t0:.0f}s] RF trained", flush=True)

    if have_xgb:
        spw = (y_train == 0).sum() / (y_train == 1).sum()
        xgb_clf = xgb.XGBClassifier(n_estimators=300, max_depth=4, learning_rate=0.05,
                                    subsample=0.8, colsample_bytree=0.8, scale_pos_weight=spw,
                                    eval_metric="logloss", random_state=42, n_jobs=-1)
        xgb_clf.fit(X_train_sel, y_train)
        results["xgboost"] = xgb_clf
        print(f"[{time.time()-t0:.0f}s] XGB trained", flush=True)

    # --- Evaluate on held-out test set ---
    metrics_rows, roc_data = [], {}
    for name, model in results.items():
        Xte = X_test_scaled if name == "logistic_regression" else X_test_sel
        proba = model.predict_proba(Xte)[:, 1]
        pred = model.predict(Xte)
        cm = confusion_matrix(y_test, pred)
        fpr, tpr, _ = roc_curve(y_test, proba)
        roc_data[name] = (fpr.tolist(), tpr.tolist(), roc_auc_score(y_test, proba))
        metrics_rows.append({"model": name, "test_auc": roc_auc_score(y_test, proba),
                             "test_accuracy": accuracy_score(y_test, pred),
                             "test_f1": f1_score(y_test, pred),
                             "test_precision": precision_score(y_test, pred),
                             "test_recall": recall_score(y_test, pred),
                             "tn": int(cm[0, 0]), "fp": int(cm[0, 1]),
                             "fn": int(cm[1, 0]), "tp": int(cm[1, 1])})
        print(f"{name}: AUC={roc_auc_score(y_test, proba):.4f} "
              f"ACC={accuracy_score(y_test, pred):.4f}", flush=True)
    metrics_df = pd.DataFrame(metrics_rows)
    metrics_df.to_csv(_p("test_metrics.csv"), index=False)

    # --- Grouped 5-fold CV on training set, feature selection re-fit INSIDE each fold ---
    # (Selecting inside the fold avoids leaking validation-fold labels into gene
    #  selection; the effect here is negligible, but keeps the CV number honest.)
    sgkf = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=1)
    strat = projects_full.iloc[train_idx].astype(str) + "_" + y_train.astype(str)
    Xtr_df = X_train.reset_index(drop=True)
    ytr = y_train.reset_index(drop=True)
    gtr = groups_full.iloc[train_idx].reset_index(drop=True)
    cv_aucs = {}
    for name, ctor in [
        ("logistic_regression", lambda: Pipeline([("sel", SelectKBest(f_classif, k=2000)),
                                                   ("scale", StandardScaler()),
                                                   ("clf", LogisticRegression(max_iter=5000, C=0.1,
                                                    class_weight="balanced", random_state=42))])),
        ("random_forest", lambda: Pipeline([("sel", SelectKBest(f_classif, k=2000)),
                                            ("clf", RandomForestClassifier(n_estimators=300, max_depth=8,
                                             min_samples_leaf=3, class_weight="balanced",
                                             n_jobs=-1, random_state=42))])),
    ]:
        aucs = []
        for tr_i, va_i in sgkf.split(Xtr_df, strat.values, groups=gtr):
            m = ctor()
            m.fit(Xtr_df.iloc[tr_i], ytr.iloc[tr_i])
            p = m.predict_proba(Xtr_df.iloc[va_i])[:, 1]
            aucs.append(roc_auc_score(ytr.iloc[va_i], p))
        cv_aucs[name] = aucs
        print(f"CV {name}: mean AUC={np.mean(aucs):.4f} +/- {np.std(aucs):.4f}", flush=True)

    # --- Per-cancer-type test breakdown (best model) ---
    best_name = metrics_df.sort_values("test_auc", ascending=False).iloc[0]["model"]
    best_model = results[best_name]
    Xte_best = X_test_scaled if best_name == "logistic_regression" else X_test_sel
    proba_best = best_model.predict_proba(Xte_best)[:, 1]
    pred_best = (proba_best >= 0.5).astype(int)
    per_cancer = []
    for proj in sorted(proj_test.unique()):
        mask = (proj_test == proj).values
        if mask.sum() < 3:
            continue
        yt = y_test.values[mask]
        auc_c = roc_auc_score(yt, proba_best[mask]) if len(set(yt)) == 2 else None
        per_cancer.append({"project": proj, "n": int(mask.sum()), "auc": auc_c,
                           "accuracy": accuracy_score(yt, pred_best[mask])})
    pd.DataFrame(per_cancer).to_csv(_p("per_cancer_type_performance.csv"), index=False)

    # --- Feature importance (top genes) ---
    gene_name_map = pd.read_pickle(_p("gene_id_to_name.pkl"))
    importances = (best_model.feature_importances_ if best_name in ("xgboost", "random_forest")
                   else np.abs(best_model.coef_[0]))
    top_idx = np.argsort(importances)[::-1][:30]
    pd.DataFrame({"gene_id": selected_genes[top_idx],
                  "gene_name": [gene_name_map.get(g, g) for g in selected_genes[top_idx]],
                  "importance": importances[top_idx]}).to_csv(_p("top_important_genes.csv"), index=False)

    with open(_p("model_results.pkl"), "wb") as f:
        pickle.dump({"results": results, "metrics_df": metrics_df, "roc_data": roc_data,
                     "cv_aucs": cv_aucs, "best_name": best_name, "selected_genes": selected_genes,
                     "y_test": y_test, "proba_best": proba_best, "pred_best": pred_best,
                     "scaler": scaler, "selector": selector}, f)
    print(f"[{time.time()-t0:.0f}s] ALL DONE. best_model={best_name}", flush=True)


if __name__ == "__main__":
    main()
