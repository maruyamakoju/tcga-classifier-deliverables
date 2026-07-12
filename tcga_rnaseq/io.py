"""Model and data I/O for tcga_rnaseq. numpy + pandas only."""
import json
import os
import pickle
import tempfile
import csv
import numbers
from collections.abc import Mapping
from pathlib import Path

import numpy as np
import pandas as pd

from .validation import validate_expression_matrix


_MODEL_REQUIRED_KEYS = {
    "selected_genes", "scaler_mean", "scaler_scale", "coef", "intercept",
}


def _model_text(value, field_name):
    if isinstance(value, (bytes, np.bytes_)):
        try:
            return bytes(value).decode("utf-8")
        except UnicodeDecodeError as exc:
            raise ValueError(f"{field_name} byte strings must be valid UTF-8") from exc
    return str(value)


def _validate_class_labels(classes, expected_count):
    classes = np.asarray(classes)
    if classes.ndim != 1 or classes.shape != (expected_count,):
        raise ValueError(
            f"Class label count must match number of coefficient rows ({expected_count})"
        )
    numeric_labels = [value for value in classes if isinstance(value, numbers.Real)]
    if numeric_labels and not all(np.isfinite(float(value)) for value in numeric_labels):
        raise ValueError("Class labels must be finite")
    raw_text = np.asarray(
        [_model_text(value, "Class labels") for value in classes], dtype=str
    )
    text = np.char.strip(raw_text)
    missing = np.asarray([
        value is None or (isinstance(value, (float, np.floating)) and np.isnan(value))
        for value in classes
    ])
    if missing.any() or np.any(text == ""):
        raise ValueError("Class labels must be non-empty")
    if np.any(raw_text != text):
        raise ValueError("Class labels must not contain leading or trailing whitespace")
    if len(set(text.tolist())) != len(text):
        raise ValueError("Class labels must be unique")
    if all(isinstance(value, (str, bytes, np.str_, np.bytes_)) for value in classes):
        return text
    return classes


def validate_lr_model(model):
    """Return a canonical, fully validated logistic-regression model dict."""
    if not isinstance(model, dict):
        raise ValueError("Model must be a dictionary")
    missing = sorted({"genes", "mean", "scale", "coef", "intercept"} - set(model))
    if missing:
        raise ValueError(f"Model is missing required fields: {', '.join(missing)}")

    genes_array = np.asarray(model["genes"])
    if genes_array.ndim != 1 or genes_array.size == 0:
        raise ValueError("Model selected_genes must be a non-empty 1-D array")
    raw_gene_text = np.asarray(
        [_model_text(value, "Model selected_genes") for value in genes_array],
        dtype=str,
    )
    genes = np.char.strip(raw_gene_text)
    missing_gene = np.asarray([
        value is None or (isinstance(value, (float, np.floating)) and np.isnan(value))
        for value in genes_array
    ])
    if missing_gene.any() or np.any(genes == ""):
        raise ValueError("Model selected_genes must be non-empty strings")
    if np.any(raw_gene_text != genes):
        raise ValueError(
            "Model selected_genes must not contain leading or trailing whitespace"
        )
    if len(set(genes.tolist())) != len(genes):
        raise ValueError("Model selected_genes must be unique")
    from .align import strip_version
    bases = [strip_version(gene) for gene in genes]
    if len(set(bases)) != len(bases):
        raise ValueError(
            "Model selected_genes contain ambiguous version-colliding identifiers"
        )

    n_genes = len(genes)
    try:
        mean = np.asarray(model["mean"], dtype=float)
        scale = np.asarray(model["scale"], dtype=float)
        coef = np.asarray(model["coef"], dtype=float)
        intercept_array = np.asarray(model["intercept"], dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Model numeric arrays could not be converted to float: {exc}") from exc

    if mean.shape != (n_genes,) or scale.shape != (n_genes,):
        raise ValueError(
            "Model scaler arrays must have one value per selected gene: "
            f"genes={n_genes}, mean={mean.shape}, scale={scale.shape}"
        )
    if not np.all(np.isfinite(mean)) or not np.all(np.isfinite(scale)):
        raise ValueError("Model scaler mean/scale arrays must be finite")
    if np.any(scale <= 0):
        raise ValueError("Model scaler scale values must be positive")
    if coef.ndim == 1:
        kind = "binary"
        if coef.shape != (n_genes,):
            raise ValueError(
                f"Binary coefficient shape {coef.shape} does not match genes={n_genes}"
            )
        if intercept_array.shape != ():
            raise ValueError("Binary intercept must be a scalar")
        intercept = float(intercept_array)
        classes = _validate_class_labels(model.get("classes", np.array([0, 1])), 2)
    elif coef.ndim == 2:
        kind = "multiclass"
        if coef.shape[0] < 2:
            raise ValueError("Multiclass model must contain at least two classes")
        if coef.shape[1] != n_genes:
            raise ValueError(
                f"Multiclass coefficient shape {coef.shape} does not match genes={n_genes}"
            )
        if intercept_array.shape != (coef.shape[0],):
            raise ValueError(
                "Multiclass intercept length must match number of coefficient rows"
            )
        intercept = intercept_array
        classes = _validate_class_labels(
            model.get("classes", np.arange(coef.shape[0])), coef.shape[0]
        )
    else:
        raise ValueError("Model coef must be a 1-D binary or 2-D multiclass array")
    if not np.all(np.isfinite(coef)) or not np.all(np.isfinite(intercept_array)):
        raise ValueError("Model coefficients and intercept must be finite")
    declared_kind = model.get("kind")
    if declared_kind is not None and declared_kind != kind:
        raise ValueError(
            f"Model kind {declared_kind!r} does not match coefficient shape ({kind})"
        )
    metadata = model.get("meta", {})
    if not isinstance(metadata, Mapping):
        raise ValueError("Model meta must be a mapping")
    return {
        "genes": genes,
        "mean": mean,
        "scale": scale,
        "coef": coef,
        "intercept": intercept,
        "classes": classes,
        "kind": kind,
        "meta": dict(metadata),
    }


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
    path = os.fspath(path)
    if not os.path.exists(path):
        raise ValueError(f"model artifact file not found: {path}")
    if os.path.isdir(path):
        raise ValueError(f"model artifact path is a directory: {path}")
    z = None
    try:
        z = np.load(path, allow_pickle=False)
        if not hasattr(z, "files"):
            raise ValueError("model artifact must be an NPZ archive")
        keys = set(z.files)
        missing = sorted(_MODEL_REQUIRED_KEYS - keys)
        if missing:
            raise ValueError(
                "model artifact is missing required arrays: " + ", ".join(missing)
            )
        if "classes" in keys and "class_order" in keys:
            declared_classes = np.asarray(z["classes"])
            declared_order = np.asarray(z["class_order"])
            classes_text = [_model_text(value, "classes") for value in declared_classes.flat]
            order_text = [_model_text(value, "class_order") for value in declared_order.flat]
            if declared_classes.shape != declared_order.shape or classes_text != order_text:
                raise ValueError("classes and class_order arrays disagree")
        if "classes" in keys:
            classes = np.asarray(z["classes"])
        elif "class_order" in keys:
            classes = np.asarray(z["class_order"])
        else:
            classes = np.array([0, 1])
        reserved = {
            "selected_genes", "scaler_mean", "scaler_scale", "coef",
            "intercept", "classes", "class_order",
        }
        meta = {
            key: (z[key].item() if z[key].ndim == 0 else np.asarray(z[key]))
            for key in keys - reserved
        }
        return validate_lr_model({
            "genes": np.asarray(z["selected_genes"]),
            "mean": np.asarray(z["scaler_mean"]),
            "scale": np.asarray(z["scaler_scale"]),
            "coef": np.asarray(z["coef"]),
            "intercept": np.asarray(z["intercept"]),
            "classes": classes,
            "meta": meta,
        })
    except ValueError as exc:
        raise ValueError(f"invalid model artifact {path}: {exc}") from exc
    except Exception as exc:
        raise ValueError(f"could not load model artifact {path}: {exc}") from exc
    finally:
        if z is not None and hasattr(z, "close"):
            z.close()


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


def load_pipeline(path, *, trusted=False):
    """Load a *trusted* legacy pickle, stubbing xgboost imports.

    Pickle can execute arbitrary code.  The xgboost stub prevents an unwanted
    import but does not make a pickle safe; callers must explicitly acknowledge
    that the file is trusted.
    """
    if not trusted:
        raise ValueError(
            "legacy pipeline pickle loading requires trusted=True; never load an "
            "untrusted pickle"
        )
    try:
        with open(path, "rb") as f:
            return _SafeUnpickler(f).load()
    except Exception as exc:
        raise ValueError(f"could not load trusted legacy pipeline {path}: {exc}") from exc


def _looks_like_ensembl_gene_id(value):
    text = str(value).strip()
    base = text.rsplit(".", 1)[0] if text.rsplit(".", 1)[-1].isdigit() else text
    return base.startswith("ENSG") and base[4:].isdigit()


def _read_delimited_header(path, delimiter):
    try:
        with open(path, "r", encoding="utf-8-sig", newline="") as handle:
            header = next(csv.reader(handle, delimiter=delimiter), None)
    except (OSError, UnicodeError, csv.Error) as exc:
        raise ValueError(f"could not read CSV header {path}: {exc}") from exc
    if not header:
        raise ValueError(f"delimited file {path} is empty")
    return header


def _validate_delimited_header(path, delimiter):
    """Detect headers pandas would otherwise silently mangle (``x`` -> ``x.1``)."""
    header = _read_delimited_header(path, delimiter)
    genes = [str(value).strip() for value in header[1:]]
    if any(value == "" for value in genes):
        raise ValueError("expression matrix gene identifiers must be non-empty")
    seen = set()
    duplicates = set()
    for value in genes:
        if value in seen:
            duplicates.add(value)
        seen.add(value)
    duplicates = sorted(duplicates)
    if duplicates:
        raise ValueError(
            "expression matrix contains duplicate gene columns: "
            + ", ".join(duplicates[:5])
        )


def _promote_parquet_row_identifier(df):
    """Restore the documented first-column row identifier for indexless Parquet.

    Pandas restores an index that was written into Parquet metadata.  With
    ``to_parquet(index=False)``, however, it returns a fresh default RangeIndex
    and leaves the documented sample/gene identifier as the first data column.
    Treat that representation like CSV/TSV instead of silently emitting
    ``0..n-1`` as sample IDs.  A deliberate numeric RangeIndex must be exported
    as an explicit column; the two encodings are otherwise indistinguishable
    after loading.
    """
    if not isinstance(df, pd.DataFrame):
        return df
    index = df.index
    is_default_range = (
        isinstance(index, pd.RangeIndex)
        and index.name is None
        and index.start == 0
        and index.stop == len(df)
        and index.step == 1
    )
    if not is_default_range:
        return df
    if df.shape[1] == 0:
        return df
    if isinstance(df.columns, pd.MultiIndex):
        return df  # the shared matrix validator emits the canonical error
    first_column = df.columns[0]
    if list(df.columns).count(first_column) != 1:
        raise ValueError(
            "indexless Parquet expression matrix has an ambiguous duplicate "
            "first-column row identifier"
        )
    return df.set_index(first_column, drop=True)


def read_matrix(path, transpose=False, allow_pickle=False):
    """Read an expression matrix by extension.

    Supported public input formats are .csv, .tsv, .txt, and .parquet. Pickled
    pandas DataFrames are only allowed for trusted internal artifacts when
    allow_pickle=True because unpickling user-controlled files can execute code.
    """
    path = os.fspath(path)
    if not os.path.exists(path):
        raise ValueError(f"expression matrix file not found: {path}")
    if os.path.isdir(path):
        raise ValueError(f"expression matrix path is a directory: {path}")
    ext = os.path.splitext(path)[1].lower()
    if ext == ".pkl" and not allow_pickle:
        raise ValueError(
            "Pickle expression inputs are disabled by default because unpickling "
            "user-controlled files can execute code; convert the matrix to CSV, "
            "TSV, or Parquet instead."
        )
    if ext not in {".pkl", ".parquet", ".tsv", ".txt", ".csv"}:
        raise ValueError(
            "unsupported expression matrix format; use .csv, .tsv, .txt, "
            ".parquet, or an explicitly trusted .pkl"
        )
    if ext in (".tsv", ".txt"):
        _validate_delimited_header(path, "\t")
    elif ext == ".csv":
        _validate_delimited_header(path, ",")

    try:
        if ext == ".pkl":
            df = pd.read_pickle(path)
        elif ext == ".parquet":
            df = pd.read_parquet(path)
        elif ext in (".tsv", ".txt"):
            df = pd.read_csv(
                path, sep="\t", index_col=0, converters={0: str}, keep_default_na=False
            )
        else:  # .csv
            df = pd.read_csv(
                path, index_col=0, converters={0: str}, keep_default_na=False
            )
    except Exception as exc:
        # Parquet engines expose backend-specific corruption exceptions (for
        # example pyarrow.lib.ArrowInvalid). Normalize all ordinary reader
        # failures into the public ValueError contract; never leak a backend
        # traceback for malformed user input.
        raise ValueError(f"could not read expression matrix {path}: {exc}") from exc
    if ext == ".parquet":
        df = _promote_parquet_row_identifier(df)
    if not transpose and _looks_like_ensembl_gene_id(df.index.name):
        raise ValueError(
            "expression matrix first column appears to be an Ensembl gene, not a "
            "sample identifier column; add a sample ID column or use --transpose"
        )
    result = df.T if transpose else df
    return validate_expression_matrix(result, f"expression matrix {path}")


def read_expression_csv(path, index_col=0):
    """Read an expression matrix CSV (rows=samples, cols=Ensembl gene IDs)."""
    if index_col == 0:
        return read_matrix(path)
    try:
        _validate_delimited_header(path, ",")
        header = pd.read_csv(path, nrows=0, keep_default_na=False)
        if isinstance(index_col, int):
            if not -len(header.columns) <= index_col < len(header.columns):
                raise ValueError(f"index_col {index_col} is outside the CSV columns")
            index_name = header.columns[index_col]
        else:
            index_name = index_col
        if index_name not in header.columns:
            raise ValueError(f"sample identifier column {index_name!r} was not found")
        frame = pd.read_csv(
            path,
            dtype={index_name: str},
            keep_default_na=False,
        ).set_index(index_name)
    except ValueError:
        raise
    except (OSError, pd.errors.ParserError, UnicodeError) as exc:
        raise ValueError(f"could not read expression matrix {path}: {exc}") from exc
    return validate_expression_matrix(frame, f"expression matrix {path}")


def read_csv_table(path, string_columns=()):
    """Read a public CSV while preserving selected identifier columns."""
    path = os.fspath(path)
    if not os.path.exists(path):
        raise ValueError(f"CSV file not found: {path}")
    if os.path.isdir(path):
        raise ValueError(f"CSV path is a directory: {path}")
    try:
        raw_header = _read_delimited_header(path, ",")
        normalized_header = [str(value).strip() for value in raw_header]
        if any(not value for value in normalized_header):
            raise ValueError(f"CSV {path} contains a blank column name")
        if raw_header != normalized_header:
            raise ValueError(
                f"CSV {path} column names must not contain leading or trailing whitespace"
            )
        if len(set(normalized_header)) != len(normalized_header):
            raise ValueError(f"CSV {path} contains duplicate column names")
        header = pd.read_csv(path, nrows=0, keep_default_na=False)
        dtype = {column: str for column in string_columns if column in header.columns}
        return pd.read_csv(path, dtype=dtype or None, keep_default_na=False)
    except (OSError, pd.errors.ParserError, UnicodeError) as exc:
        raise ValueError(f"could not read CSV {path}: {exc}") from exc


def _canonical_path(path):
    return os.path.normcase(os.path.realpath(os.path.abspath(os.fspath(path))))


def ensure_distinct_paths(outputs, inputs=None):
    """Reject output/output and output/input path collisions.

    ``outputs`` and ``inputs`` are mappings from human-readable names to paths;
    false/None paths are ignored.
    """
    if not isinstance(outputs, Mapping):
        raise ValueError("outputs must be a mapping of names to paths")
    if inputs is None:
        inputs = {}
    if not isinstance(inputs, Mapping):
        raise ValueError("inputs must be a mapping of names to paths")
    output_items = [(name, path) for name, path in outputs.items() if path]
    input_items = [(name, path) for name, path in inputs.items() if path]
    seen = {}
    for name, path in output_items:
        canonical = _canonical_path(path)
        if canonical in seen:
            raise ValueError(
                f"output paths for {seen[canonical]} and {name} refer to the same file: {path}"
            )
        for previous_name, previous_path in output_items:
            if previous_name == name:
                break
            if os.path.exists(path) and os.path.exists(previous_path):
                try:
                    if os.path.samefile(path, previous_path):
                        raise ValueError(
                            f"output paths for {previous_name} and {name} refer to "
                            f"the same file: {path}"
                        )
                except OSError:
                    pass
        seen[canonical] = name
        for input_name, input_path in input_items:
            same = canonical == _canonical_path(input_path)
            if not same and os.path.exists(path) and os.path.exists(input_path):
                try:
                    same = os.path.samefile(path, input_path)
                except OSError:
                    same = False
            if same:
                raise ValueError(
                    f"refusing to overwrite {input_name} with {name}: {path}"
                )


def _atomic_text_write(path, text, encoding="utf-8"):
    target = Path(path)
    try:
        target.parent.mkdir(parents=True, exist_ok=True)
        fd, temporary = tempfile.mkstemp(
            prefix=f".{target.name}.", suffix=".tmp", dir=str(target.parent)
        )
        try:
            with os.fdopen(fd, "w", encoding=encoding, newline="") as handle:
                handle.write(text)
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(temporary, target)
        except Exception:
            try:
                os.unlink(temporary)
            except OSError:
                pass
            raise
    except Exception as exc:
        raise ValueError(f"could not write {target}: {exc}") from exc
    return path


def write_dataframe_csv(df, path, index=False, **kwargs):
    """Atomically write a DataFrame CSV, creating its parent directory."""
    target = Path(path)
    kwargs.setdefault("lineterminator", "\n")
    try:
        target.parent.mkdir(parents=True, exist_ok=True)
        fd, temporary = tempfile.mkstemp(
            prefix=f".{target.name}.", suffix=".tmp", dir=str(target.parent)
        )
        os.close(fd)
        try:
            df.to_csv(temporary, index=index, **kwargs)
            # Windows' CRT rejects fsync on a read-only descriptor; r+b keeps
            # the already-written bytes unchanged while providing a durable
            # writable handle on every supported OS.
            with open(temporary, "r+b") as handle:
                os.fsync(handle.fileno())
            os.replace(temporary, target)
        except Exception:
            try:
                os.unlink(temporary)
            except OSError:
                pass
            raise
    except Exception as exc:
        raise ValueError(f"could not write CSV {target}: {exc}") from exc
    return path


def write_text(path, text, encoding="utf-8"):
    """Atomically write text, creating its parent directory."""
    return _atomic_text_write(path, text, encoding=encoding)


def write_json(obj, path, indent=2, sort_keys=False):
    """Write a JSON report (trailing newline), creating parent dirs as needed."""
    try:
        content = json.dumps(
            obj, indent=indent, sort_keys=sort_keys, allow_nan=False
        ) + "\n"
    except (TypeError, ValueError) as exc:
        raise ValueError(f"could not serialize JSON for {path}: {exc}") from exc
    return _atomic_text_write(path, content)
