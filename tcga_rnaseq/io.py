"""Model and data I/O for tcga_rnaseq. numpy + pandas only."""
import json
import os
import pickle

import numpy as np
import pandas as pd


def load_lr_model(path):
    """Load a deployable logistic-regression model from a .npz file.

    Handles both the binary tumor-vs-normal export and the multi-class
    cancer-type export, normalizing their slightly different key names into one
    dict:

        genes      (g,)   str      selected gene IDs, model column order
        mean       (g,)   float    per-gene standardizer mean (training)
        scale      (g,)   float    per-gene standardizer scale (training)
        coef       binary: (g,)    multiclass: (n_classes, g)
        intercept  binary: float   multiclass: (n_classes,)
        classes    (n_classes,)    class labels
        kind       "binary" | "multiclass"
        meta       dict of any extra arrays/notes in the file

    The npz holds only arrays/strings, so it loads with allow_pickle=False
    (no code-execution surface).
    """
    z = np.load(path, allow_pickle=False)
    keys = set(z.files)
    genes = z["selected_genes"].astype(str)
    mean = z["scaler_mean"].astype(float)
    scale = z["scaler_scale"].astype(float)
    coef = z["coef"].astype(float)
    kind = "multiclass" if coef.ndim == 2 else "binary"
    if "classes" in keys:
        classes = z["classes"].astype(str)
    elif "class_order" in keys:
        classes = np.asarray(z["class_order"])
    else:
        classes = np.array([0, 1])
    intercept = float(z["intercept"]) if kind == "binary" else z["intercept"].astype(float)
    n_genes = len(genes)
    if mean.shape != (n_genes,) or scale.shape != (n_genes,):
        raise ValueError(
            "Model scaler arrays must have one value per selected gene: "
            f"genes={n_genes}, mean={mean.shape}, scale={scale.shape}"
        )
    if not np.all(np.isfinite(mean)) or not np.all(np.isfinite(scale)):
        raise ValueError("Model scaler mean/scale arrays must be finite")
    if np.any(scale <= 0):
        raise ValueError("Model scaler scale values must be positive")
    if kind == "binary":
        if coef.shape != (n_genes,):
            raise ValueError(
                f"Binary coefficient shape {coef.shape} does not match genes={n_genes}"
            )
    else:
        if coef.shape[1] != n_genes:
            raise ValueError(
                f"Multiclass coefficient shape {coef.shape} does not match genes={n_genes}"
            )
        if np.asarray(intercept).shape != (coef.shape[0],):
            raise ValueError(
                "Multiclass intercept length must match number of coefficient rows"
            )
        if len(classes) != coef.shape[0]:
            raise ValueError("Class label count must match number of coefficient rows")
    reserved = {"selected_genes", "scaler_mean", "scaler_scale", "coef",
                "intercept", "classes", "class_order"}
    meta = {k: (z[k].item() if z[k].ndim == 0 else z[k]) for k in keys - reserved}
    return {"genes": genes, "mean": mean, "scale": scale, "coef": coef,
            "intercept": intercept, "classes": classes, "kind": kind, "meta": meta}


class _XgbStub:
    """Absorbs a pickled xgboost model so unpickling never triggers
    `import xgboost` (which segfaults in the project conda env). Never scored."""
    def __init__(self, *a, **k):
        pass

    def __setstate__(self, state):
        self._state = state


class _SafeUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if module.startswith("xgboost"):
            return _XgbStub
        return super().find_class(module, name)


def load_pipeline(path):
    """Load the legacy pickled sklearn pipeline, stubbing xgboost safely."""
    with open(path, "rb") as f:
        return _SafeUnpickler(f).load()


def read_matrix(path, transpose=False):
    """Read an expression matrix by extension (.csv/.tsv/.txt/.parquet/.pkl)."""
    ext = os.path.splitext(path)[1].lower()
    if ext == ".pkl":
        df = pd.read_pickle(path)
    elif ext == ".parquet":
        df = pd.read_parquet(path)
    elif ext in (".tsv", ".txt"):
        df = pd.read_csv(path, sep="\t", index_col=0)
    else:
        df = pd.read_csv(path, index_col=0)
    return df.T if transpose else df


def read_expression_csv(path, index_col=0):
    """Read an expression matrix CSV (rows=samples, cols=Ensembl gene IDs)."""
    return pd.read_csv(path, index_col=index_col)


def write_json(obj, path, indent=2, sort_keys=False):
    """Write a JSON report (trailing newline), creating parent dirs as needed."""
    parent = os.path.dirname(os.path.abspath(path))
    os.makedirs(parent, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=indent, sort_keys=sort_keys)
        f.write("\n")
    return path
