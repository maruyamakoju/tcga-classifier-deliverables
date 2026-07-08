"""Gene-column alignment for tcga_rnaseq."""
import numpy as np
import pandas as pd


def _preview(values, limit=5):
    text = [str(value) for value in values[:limit]]
    suffix = ", ..." if len(values) > limit else ""
    return ", ".join(text) + suffix


def strip_version(gene_id):
    """Drop the Ensembl version suffix: 'ENSG00000005.6' -> 'ENSG00000005'."""
    return str(gene_id).split(".")[0]


def build_gene_column_lookups(columns):
    """Build exact and version-stripped lookups for input gene columns.

    The model accepts Ensembl IDs with or without ``.version`` suffixes, but
    duplicate or version-colliding columns are ambiguous enough to reject.
    """
    exact = {}
    duplicate_exact = []
    for column in columns:
        key = str(column)
        if key in exact:
            duplicate_exact.append(key)
        else:
            exact[key] = column
    if duplicate_exact:
        raise ValueError(
            "Duplicate gene columns are not allowed: " + _preview(sorted(set(duplicate_exact)))
        )

    stripped = {}
    collisions = {}
    for column in columns:
        base = strip_version(column)
        if base in stripped:
            collisions.setdefault(base, [str(stripped[base])]).append(str(column))
        else:
            stripped[base] = column
    if collisions:
        examples = []
        for base in sorted(collisions)[:5]:
            examples.append(f"{base} -> {', '.join(collisions[base])}")
        suffix = ", ..." if len(collisions) > 5 else ""
        raise ValueError(
            "Ambiguous gene columns after removing Ensembl version suffix: "
            + "; ".join(examples)
            + suffix
        )
    return exact, stripped


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
    exact, stripped = build_gene_column_lookups(X.columns)

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
