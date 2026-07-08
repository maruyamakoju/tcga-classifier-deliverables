#!/usr/bin/env python
"""
Re-export the training feature matrix to a version-neutral .npy so the rest of
the cancer-type pipeline runs on any numpy/pandas.

X_full_filtered.pkl was pickled with numpy 2.x + pandas 3.0 (its gene-ID column
index uses the pandas 3.0 string dtype), so it will NOT unpickle on older
numpy/pandas. Run THIS script once in an environment with numpy>=2 and
pandas>=3 (e.g. a fresh venv: `pip install "numpy>=2.1" "pandas>=3.0"`).

Writes, next to this script:
  X_full.npy      (n_samples x n_genes float32)
  X_genes.npy     (n_genes,)  Ensembl gene IDs
  X_samples.npy   (n_samples,) GDC file_id (matches ../selected_files.csv)
"""
import os, pickle
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def main():
    with open(os.path.join(ROOT, "X_full_filtered.pkl"), "rb") as f:
        X = pickle.load(f)
    np.save(os.path.join(HERE, "X_full.npy"), np.asarray(X.values, dtype=np.float32))
    np.save(os.path.join(HERE, "X_genes.npy"), np.asarray([str(c) for c in X.columns]))
    np.save(os.path.join(HERE, "X_samples.npy"), np.asarray([str(i) for i in X.index]))
    print(f"exported X_full.npy {X.shape} + X_genes.npy + X_samples.npy to {HERE}")


if __name__ == "__main__":
    main()
