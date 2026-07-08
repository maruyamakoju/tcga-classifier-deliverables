"""Regression tests: current code must reproduce the verified golden numbers."""
import json
import os

import numpy as np
import pandas as pd
import pytest

from tcga_rnaseq import predict_proba
from tcga_rnaseq import metrics as M
from conftest import load_external


def _acc(p, y, thr=0.5):
    return float(np.mean((p >= thr).astype(int) == y))


def test_binary_reproduces_example_output(binary_model, example_input, root, golden):
    p = predict_proba(binary_model, example_input)
    exp = pd.read_csv(os.path.join(root, "example_output.csv"))
    col = [c for c in exp.columns if "prob" in c.lower()][0]
    tol = golden["binary_tumor_normal"]["reproduce_example_output_max_abs_delta"]
    assert np.abs(p - exp[col].values).max() < tol


def test_shipped_summary_artifacts_match_golden(root, golden):
    metrics = pd.read_csv(os.path.join(root, "test_metrics.csv")).set_index("model")
    lr = metrics.loc["logistic_regression"]
    assert lr["test_auc"] == pytest.approx(
        golden["binary_tumor_normal"]["heldout_test_auc"], abs=1e-4
    )
    assert lr["test_accuracy"] == pytest.approx(
        golden["binary_tumor_normal"]["heldout_test_accuracy"], abs=1e-4
    )

    loco = pd.read_csv(
        os.path.join(root, "cross-cancer-holdout", "loco_pooled_summary.csv")
    ).set_index("model")
    assert loco.loc["logistic_regression", "macro_mean_auc"] == pytest.approx(
        golden["binary_tumor_normal"]["loco_macro_mean_auc"], abs=5e-4
    )

    with open(os.path.join(root, "cancer-type-classifier", "cancer_type_summary.json")) as f:
        cancer_type = json.load(f)
    expected = golden["cancer_type"]
    assert cancer_type["n_tumors"] == expected["n_tumors"]
    assert cancer_type["n_classes"] == expected["n_classes"]
    assert cancer_type["accuracy"] == pytest.approx(
        expected["patient_heldout_accuracy"], abs=1e-4
    )
    assert cancer_type["balanced_accuracy"] == pytest.approx(
        expected["patient_heldout_balanced_accuracy"], abs=1e-4
    )
    assert cancer_type["macro_f1"] == pytest.approx(
        expected["patient_heldout_macro_f1"], abs=1e-4
    )


def test_train_test_split_is_patient_disjoint(root):
    train_idx = np.load(os.path.join(root, "train_idx.npy"))
    test_idx = np.load(os.path.join(root, "test_idx.npy"))
    meta = pd.read_csv(os.path.join(root, "selected_files.csv"))
    assert set(train_idx).isdisjoint(set(test_idx))
    for column in ["case_id", "submitter_id"]:
        groups = meta[column].astype(str).to_numpy()
        assert set(groups[train_idx]).isdisjoint(set(groups[test_idx]))


def test_toil_baseline_and_adaptation(binary_model, golden):
    g = golden["external_validation"]["toil_rsem"]
    tol = golden["tolerance"]
    X, y = load_external("toil")
    # The cached Toil/RSEM matrix has four selected genes absent from the
    # source extract. These regression tests reproduce the historical
    # mean-imputed benchmark explicitly; public scoring APIs reject this by
    # default.
    p0 = predict_proba(binary_model, X, adapt="none", allow_invalid_values=True)
    assert M.roc_auc(y, p0) == pytest.approx(g["baseline_auc"], abs=tol)
    assert _acc(p0, y) == pytest.approx(g["baseline_acc_at_0p5"], abs=tol)
    assert int(((p0 >= 0.5).astype(int) == y).sum()) == round(
        g["baseline_acc_at_0p5"] * len(y)
    )
    pa = predict_proba(binary_model, X, adapt="cohort_zscore", allow_invalid_values=True)
    assert _acc(pa, y) == pytest.approx(g["cohort_zscore_acc_at_0p5"], abs=tol)
    assert M.roc_auc(y, pa) == pytest.approx(g["cohort_zscore_auc"], abs=tol)
    assert int(((pa >= 0.5).astype(int) == y).sum()) == round(
        g["cohort_zscore_acc_at_0p5"] * len(y)
    )


def test_cptac_baseline(binary_model, golden):
    g = golden["external_validation"]["cptac_gdc"]
    tol = golden["tolerance"]
    X, y = load_external("cptac")
    p0 = predict_proba(binary_model, X, adapt="none")
    assert M.roc_auc(y, p0) == pytest.approx(g["baseline_auc"], abs=tol)
    assert _acc(p0, y) == pytest.approx(g["baseline_acc_at_0p5"], abs=tol)
    assert int(((p0 >= 0.5).astype(int) == y).sum()) == round(
        g["baseline_acc_at_0p5"] * len(y)
    )


def test_gtex_normals_fpr(binary_model, golden):
    g = golden["external_validation"]["gtex_normals"]
    tol = golden["tolerance"]
    X, _ = load_external("gtex")
    # The cached GTEx/Toil matrix has the same four absent selected genes as the
    # TCGA-Toil extract; keep the historical benchmark explicit.
    p0 = predict_proba(binary_model, X, adapt="none", allow_invalid_values=True)
    assert float(np.mean(p0 >= 0.5)) == pytest.approx(g["baseline_fpr_at_0p5"], abs=tol)
    assert int((p0 >= 0.5).sum()) == round(g["baseline_fpr_at_0p5"] * len(X))
    pa = predict_proba(binary_model, X, adapt="cohort_zscore", allow_invalid_values=True)
    assert float(np.mean(pa >= 0.5)) == pytest.approx(g["cohort_zscore_fpr_at_0p5"], abs=0.02)
    assert int((pa >= 0.5).sum()) == round(g["cohort_zscore_fpr_at_0p5"] * len(X))


def test_cancer_type_npz_reproduces_sklearn(cancer_type_model, features_npy, root, golden):
    """The exported multiclass npz must match a freshly refit sklearn pipeline."""
    sk_pipe = pytest.importorskip("sklearn.pipeline")
    from sklearn.preprocessing import StandardScaler
    from sklearn.feature_selection import SelectKBest, f_classif
    from sklearn.linear_model import LogisticRegression
    d = os.path.dirname(features_npy)
    X = np.load(features_npy)
    genes = np.load(os.path.join(d, "X_genes.npy"), allow_pickle=True).astype(str)
    samp = np.load(os.path.join(d, "X_samples.npy"), allow_pickle=True).astype(str)
    sf = pd.read_csv(os.path.join(root, "selected_files.csv")).set_index("file_id")
    meta = sf.reindex(samp)
    ct = meta["project"].str.replace("TCGA-", "", regex=False).values
    tum = meta["label"].values == "tumor"
    Xt, yt = X[tum], ct[tum]
    pipe = sk_pipe.Pipeline([("sc", StandardScaler()),
                             ("sel", SelectKBest(f_classif, k=1000)),
                             ("lr", LogisticRegression(C=2.0, max_iter=5000))]).fit(Xt, yt)
    proba_npz = predict_proba(cancer_type_model, pd.DataFrame(Xt, columns=genes))
    proba_sklearn = pipe.predict_proba(Xt)
    assert list(pipe.classes_) == list(cancer_type_model["classes"])
    assert np.abs(proba_npz - proba_sklearn).max() <= golden[
        "cancer_type"
    ]["npz_reproduces_sklearn_max_abs_delta"]
    agree = np.mean(pipe.predict(Xt) == cancer_type_model["classes"][proba_npz.argmax(1)])
    assert agree == pytest.approx(1.0)


@pytest.mark.slow
def test_binary_heldout_auc_full_data(binary_model, features_npy, root, golden):
    """Full-data reproduction of the headline held-out test AUC (patient-disjoint split)."""
    from sklearn.feature_selection import SelectKBest, f_classif
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    d = os.path.dirname(features_npy)
    X = np.load(features_npy)
    samp = np.load(os.path.join(d, "X_samples.npy"), allow_pickle=True).astype(str)
    sf = pd.read_csv(os.path.join(root, "selected_files.csv")).set_index("file_id")
    y = (sf.reindex(samp)["label"].values == "tumor").astype(int)
    tr = np.load(os.path.join(root, "train_idx.npy"))
    te = np.load(os.path.join(root, "test_idx.npy"))
    sel = SelectKBest(f_classif, k=2000).fit(X[tr], y[tr])
    sc = StandardScaler().fit(sel.transform(X[tr]))
    lr = LogisticRegression(max_iter=5000, C=0.1, class_weight="balanced",
                            random_state=42).fit(sc.transform(sel.transform(X[tr])), y[tr])
    p = lr.predict_proba(sc.transform(sel.transform(X[te])))[:, 1]
    exp = golden["binary_tumor_normal"]["heldout_test_auc"]
    assert roc_auc_score(y[te], p) == pytest.approx(exp, abs=golden["tolerance"])
