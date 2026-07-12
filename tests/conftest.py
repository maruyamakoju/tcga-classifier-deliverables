"""Shared pytest fixtures for the tcga_rnaseq test suite."""
import json
import os
import sys

import numpy as np
import pandas as pd
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)


@pytest.fixture(scope="session")
def root():
    return ROOT


@pytest.fixture(scope="session")
def golden():
    with open(os.path.join(os.path.dirname(__file__), "golden_numbers.json")) as f:
        return json.load(f)


@pytest.fixture(scope="session")
def binary_model():
    from tcga_rnaseq import load_lr_model
    return load_lr_model(os.path.join(ROOT, "deployable_lr_weights.npz"))


@pytest.fixture(scope="session")
def cancer_type_model():
    from tcga_rnaseq import load_lr_model
    return load_lr_model(os.path.join(ROOT, "cancer-type-classifier", "cancer_type_lr_weights.npz"))


@pytest.fixture(scope="session")
def example_input():
    from tcga_rnaseq import read_expression_csv
    return read_expression_csv(os.path.join(ROOT, "example_input.csv"))


def load_external(name):
    """Load a shipped external-validation selected-gene matrix + 0/1 labels.

    These per-cohort matrices are large and may not be bundled in a lightweight
    release; the test skips when they are absent.
    """
    ev = os.path.join(ROOT, "external-validation")
    spec = {
        "toil": ("tcga_toil_xena/tcga_toil_selected_genes_model_scale.pkl",
                 "tcga_toil_xena/tcga_toil_predictions.csv", "sample"),
        "cptac": ("cptac_gdc/expression_selected_genes.pkl",
                  "cptac_gdc/cptac_predictions.csv", "file_id"),
        "gtex": ("gtex_xena/gtex_selected_genes_model_scale.pkl", None, None),
    }
    mp, lp, key = spec[name]
    if not os.path.exists(os.path.join(ev, mp)):
        pytest.skip(f"external cohort matrix not bundled: {mp}")
    X = pd.read_pickle(os.path.join(ev, mp))
    if lp is None:
        y = np.zeros(len(X), dtype=int)  # GTEx is all-normal
    else:
        lab = pd.read_csv(os.path.join(ev, lp))
        y = lab.set_index(lab[key].astype(str))["label"].reindex(X.index.astype(str)).astype(int).values
    return X, y


@pytest.fixture(scope="session")
def features_npy():
    """Path to the full 14850-gene feature matrix, if available.

    Not shipped (128 MB). Set TCGA_FEATURES to an X_full.npy (with sibling
    X_genes.npy / X_samples.npy) to enable the full-data reproduction tests;
    otherwise those tests skip. Regenerate with the path-neutral canonical lock
    and exact converter environment documented in REPRODUCIBILITY.md.
    """
    for cand in [os.environ.get("TCGA_FEATURES"),
                 os.path.join(ROOT, "cancer-type-classifier", "X_full.npy")]:
        if cand and os.path.exists(cand):
            values = np.load(cand, mmap_mode="r", allow_pickle=False)
            if values.dtype != np.float32:
                pytest.fail(f"exact cancer-type feature matrix must be float32: {cand}")
            return cand
    pytest.skip("full feature matrix not available (set TCGA_FEATURES)")


@pytest.fixture(scope="session")
def features_float64_npy():
    """Exact float64 matrix used by the binary and LOCO analyses, if available."""
    candidates = [
        os.environ.get("TCGA_FEATURES_FLOAT64"),
        os.path.join(ROOT, "cancer-type-classifier", "X_full_float64.npy"),
    ]
    for candidate in candidates:
        if candidate and os.path.exists(candidate):
            values = np.load(candidate, mmap_mode="r", allow_pickle=False)
            if values.dtype != np.float64:
                pytest.fail(f"exact LOCO feature matrix must be float64: {candidate}")
            return candidate
    pytest.skip(
        "exact float64 feature matrix not available (set TCGA_FEATURES_FLOAT64)"
    )
