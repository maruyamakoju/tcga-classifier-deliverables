"""Gene-column alignment for tcga_rnaseq."""
import numpy as np
import pandas as pd


def strip_version(gene_id):
    """Drop the Ensembl version suffix: 'ENSG00000005.6' -> 'ENSG00000005'."""
    return str(gene_id).split(".")[0]


def align_to_genes(X, genes, impute_mean=None):
    """Reindex an expression DataFrame to the model's gene order.

    Matches Ensembl IDs with OR without the ``.version`` suffix (so a user CSV
    using unversioned IDs still aligns to a versioned model, and vice versa).
    Genes present in the input are coerced to numeric; genes missing from the
    input are filled with ``impute_mean`` (a neutral, standardized-zero value).

    X            samples x genes DataFrame
    genes        (g,) model gene order (may be versioned)
    impute_mean  (g,) per-gene training mean, or None to fill missing with NaN

    Returns (values ndarray (n_samples, g), n_matched int, missing list).
    """
    genes = [str(g) for g in genes]
    exact = {str(c): c for c in X.columns}
    stripped = {}
    for c in X.columns:  # first column wins for a given base id
        stripped.setdefault(strip_version(c), c)

    if impute_mean is None:
        means = np.full(len(genes), np.nan)
    else:
        means = np.asarray(impute_mean, dtype=float)
        if means.shape[0] != len(genes):
            raise ValueError(
                f"impute_mean length {means.shape[0]} does not match gene count {len(genes)}"
            )

    out = np.empty((X.shape[0], len(genes)), dtype=float)
    missing = []
    for j, g in enumerate(genes):
        src = exact.get(g)
        if src is None:
            src = stripped.get(strip_version(g))
        if src is None:
            out[:, j] = means[j]
            missing.append(g)
        else:
            values = pd.to_numeric(X[src], errors="coerce").to_numpy(dtype=float)
            out[:, j] = np.where(np.isfinite(values), values, means[j])
    return out, len(genes) - len(missing), missing
