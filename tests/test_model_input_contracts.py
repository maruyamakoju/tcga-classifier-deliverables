"""Regression tests for model artifacts and expression identity contracts."""
import numpy as np
import pandas as pd
import pytest

from export_lr_weights import _atomic_write_validated_npz
from tcga_rnaseq import (
    ensure_distinct_paths,
    load_lr_model,
    load_pipeline,
    read_csv_table,
    read_matrix,
    validate_lr_model,
    write_dataframe_csv,
    write_json,
)
from tcga_rnaseq.align import strip_version


def valid_model_arrays():
    return {
        "selected_genes": np.array(["ENSG00000000001.1", "ENSG00000000002.2"]),
        "scaler_mean": np.array([0.0, 0.0]),
        "scaler_scale": np.array([1.0, 1.0]),
        "coef": np.array([1.0, -1.0]),
        "intercept": np.array(0.0),
        "class_order": np.array([0, 1]),
    }


def test_weight_export_writer_validates_before_atomic_replace(tmp_path, monkeypatch):
    output = tmp_path / "weights.npz"
    output.write_bytes(b"previous-valid-generation")

    def fail_after_partial_write(handle, **_arrays):
        handle.write(b"partial")
        raise OSError("injected write failure")

    monkeypatch.setattr(np, "savez_compressed", fail_after_partial_write)
    with pytest.raises(OSError, match="injected write failure"):
        _atomic_write_validated_npz(output, valid_model_arrays())

    assert output.read_bytes() == b"previous-valid-generation"
    assert not list(tmp_path.glob(f".{output.name}.*.tmp"))


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda values: values.pop("coef"), "missing required arrays: coef"),
        (
            lambda values: values.update(coef=np.array([np.nan, 1.0])),
            "coefficients and intercept must be finite",
        ),
        (
            lambda values: values.update(intercept=np.array(np.inf)),
            "coefficients and intercept must be finite",
        ),
        (
            lambda values: values.update(
                selected_genes=np.array(["ENSG00000000001.1", "ENSG00000000001.2"])
            ),
            "version-colliding",
        ),
        (
            lambda values: values.update(selected_genes=np.array(["g1", "g1"])),
            "must be unique",
        ),
        (
            lambda values: values.update(selected_genes=np.array([" g1", "g2"])),
            "leading or trailing whitespace",
        ),
        (
            lambda values: values.update(intercept=np.array([0.0])),
            "intercept must be a scalar",
        ),
    ],
)
def test_model_artifact_rejects_invalid_schema(tmp_path, mutation, message):
    values = valid_model_arrays()
    mutation(values)
    path = tmp_path / "bad.npz"
    np.savez(path, **values)

    with pytest.raises(ValueError, match=message):
        load_lr_model(path)


def test_corrupt_model_file_is_normalized_to_valueerror(tmp_path):
    path = tmp_path / "not-an-npz.npz"
    path.write_text("not an archive", encoding="utf-8")

    with pytest.raises(ValueError, match="model artifact"):
        load_lr_model(path)


def test_model_loader_decodes_utf8_byte_identifiers(tmp_path):
    values = valid_model_arrays()
    values["selected_genes"] = np.array([b"g1", b"g2"])
    values["class_order"] = np.array([b"normal", b"tumor"])
    path = tmp_path / "bytes.npz"
    np.savez(path, **values)

    model = load_lr_model(path)

    assert model["genes"].tolist() == ["g1", "g2"]
    assert model["classes"].tolist() == ["normal", "tumor"]


def test_model_loader_rejects_conflicting_class_declarations(tmp_path):
    values = valid_model_arrays()
    values["classes"] = np.array(["tumor", "normal"])
    path = tmp_path / "conflicting-classes.npz"
    np.savez(path, **values)

    with pytest.raises(ValueError, match="classes and class_order arrays disagree"):
        load_lr_model(path)


def test_direct_model_validation_rejects_nonmapping_metadata():
    model = {
        "genes": np.array(["g1"]),
        "mean": np.array([0.0]),
        "scale": np.array([1.0]),
        "coef": np.array([1.0]),
        "intercept": 0.0,
        "meta": None,
    }
    with pytest.raises(ValueError, match="meta must be a mapping"):
        validate_lr_model(model)


def test_pickle_loader_requires_explicit_trust(tmp_path):
    path = tmp_path / "plain.pkl"
    path.write_bytes(b"N.")  # pickle protocol 0: None

    with pytest.raises(ValueError, match="trusted=True"):
        load_pipeline(path)
    assert load_pipeline(path, trusted=True) is None


def test_csv_sample_ids_preserve_leading_zero_and_literal_na(tmp_path):
    path = tmp_path / "ids.csv"
    path.write_text(
        "sample,ENSG00000000001.1\n001,1.0\nNA,2.0\n",
        encoding="utf-8",
    )

    frame = read_matrix(path)

    assert frame.index.tolist() == ["001", "NA"]
    assert all(isinstance(sample_id, str) for sample_id in frame.index.tolist())


@pytest.mark.parametrize(
    "rows",
    [
        ",1.0\ns2,2.0\n",
        " s1 ,1.0\ns2,2.0\n",
        "s1,1.0\ns1,2.0\n",
    ],
)
def test_csv_rejects_invalid_sample_identifiers(tmp_path, rows):
    path = tmp_path / "bad_ids.csv"
    path.write_text("sample,ENSG00000000001.1\n" + rows, encoding="utf-8")

    with pytest.raises(ValueError, match="sample identifiers|duplicate sample"):
        read_matrix(path)


def test_csv_rejects_empty_matrix_and_missing_sample_column(tmp_path):
    empty = tmp_path / "empty.csv"
    empty.write_text("sample,ENSG00000000001.1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="at least one sample"):
        read_matrix(empty)

    no_sample = tmp_path / "no_sample.csv"
    no_sample.write_text(
        "ENSG00000000001.1,ENSG00000000002.2\n1.0,2.0\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="first column appears to be an Ensembl gene"):
        read_matrix(no_sample)

    padded_gene = tmp_path / "padded_first_gene.csv"
    padded_gene.write_text(
        " ENSG00000000001.1 ,ENSG00000000002.2\n1.0,2.0\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="first column appears to be an Ensembl gene"):
        read_matrix(padded_gene)


def test_numeric_only_gene_version_stripping():
    assert strip_version("ENSG00000000001.12") == "ENSG00000000001"
    assert strip_version("HLA.DRA") == "HLA.DRA"
    assert strip_version("GENE.A.2") == "GENE.A"
    assert strip_version("GENE.2A") == "GENE.2A"


def test_parquet_sample_contract_is_applied(tmp_path):
    pytest.importorskip("pyarrow")
    path = tmp_path / "bad.parquet"
    pd.DataFrame({"g1": [1.0, 2.0]}, index=["s1", "s1"]).to_parquet(path)

    with pytest.raises(ValueError, match="duplicate sample"):
        read_matrix(path)


def test_indexless_parquet_promotes_explicit_sample_column(tmp_path):
    pytest.importorskip("pyarrow")
    path = tmp_path / "indexless.parquet"
    pd.DataFrame(
        {
            "sample": ["001", "NA"],
            "ENSG00000000001.1": [1.0, 2.0],
        }
    ).to_parquet(path, index=False)

    observed = read_matrix(path)

    assert observed.index.tolist() == ["001", "NA"]
    assert observed.columns.tolist() == ["ENSG00000000001.1"]


def test_indexless_parquet_without_row_identifier_fails_closed(tmp_path):
    pytest.importorskip("pyarrow")
    path = tmp_path / "genes_only.parquet"
    pd.DataFrame(
        {
            "ENSG00000000001.1": [1.0],
            "ENSG00000000002.2": [2.0],
        }
    ).to_parquet(path, index=False)

    with pytest.raises(ValueError, match="first column appears to be an Ensembl gene"):
        read_matrix(path)


def test_parquet_preserves_explicit_nondefault_index(tmp_path):
    pytest.importorskip("pyarrow")
    path = tmp_path / "indexed.parquet"
    expected = pd.DataFrame(
        {"ENSG00000000001.1": [1.0, 2.0]},
        index=pd.Index(["sample-a", "sample-b"], name="specimen"),
    )
    expected.to_parquet(path)

    pd.testing.assert_frame_equal(read_matrix(path), expected)


def test_corrupt_parquet_is_normalized_to_public_valueerror(tmp_path):
    pytest.importorskip("pyarrow")
    path = tmp_path / "corrupt.parquet"
    path.write_bytes(b"not-a-parquet-file")

    with pytest.raises(ValueError, match="could not read expression matrix"):
        read_matrix(path)


def test_expression_reader_rejects_unsupported_extension(tmp_path):
    path = tmp_path / "matrix.xlsx"
    path.write_text("sample,g1\ns1,1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="unsupported expression matrix format"):
        read_matrix(path)


def test_public_csv_table_rejects_duplicate_headers_before_pandas_mangling(tmp_path):
    path = tmp_path / "labels.csv"
    path.write_text("sample,sample,label\ns1,s1,tumor\n", encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate column names"):
        read_csv_table(path, string_columns=("sample",))


def test_json_writer_rejects_nonstandard_nan_without_partial_file(tmp_path):
    path = tmp_path / "report.json"
    with pytest.raises(ValueError, match="serialize JSON"):
        write_json({"value": np.nan}, path)
    assert not path.exists()


def test_atomic_csv_writer_uses_cross_platform_lf_newlines(tmp_path):
    path = tmp_path / "output.csv"
    write_dataframe_csv(pd.DataFrame({"sample": ["s1"], "value": [1]}), path)
    assert path.read_bytes() == b"sample,value\ns1,1\n"


def test_distinct_path_api_allows_omitted_inputs_and_validates_mappings(tmp_path):
    ensure_distinct_paths({"output": tmp_path / "out.csv"})
    with pytest.raises(ValueError, match="outputs must be a mapping"):
        ensure_distinct_paths([("output", tmp_path / "out.csv")])
    with pytest.raises(ValueError, match="inputs must be a mapping"):
        ensure_distinct_paths({"output": tmp_path / "out.csv"}, [])
