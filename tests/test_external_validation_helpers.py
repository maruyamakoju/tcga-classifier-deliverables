"""Unit tests for external-validation cache/merge-integrity helpers.

Pure-logic tests only -- no real network calls. external-validation/ is not
part of the lightweight release (not in release_tools.common.RELEASE_FILES),
so this intentionally does not attempt to exercise the live GDC/Xena download
paths; it covers the cache-fingerprinting and merge-integrity logic added to
fix a real bug (a cache keyed only on a fixed path, with no fingerprint of
sampling parameters/selected genes, could silently return stale data for the
wrong samples after a rerun with different arguments).
"""
import gzip
import hashlib
import io
import json
import os
import pickle
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXTERNAL_VALIDATION = os.path.join(ROOT, "external-validation")
sys.path.insert(0, EXTERNAL_VALIDATION)

import validate_gtex_xena as gx  # noqa: E402
import validate_cptac_gdc as cg  # noqa: E402
import validate_tcga_toil_xena as tx  # noqa: E402
import provenance as pv  # noqa: E402

GDC_FILE_IDS = [f"00000000-0000-4000-8000-{index:012d}" for index in range(1, 20)]


def _gdc_source(file_id, payload, *, revision="fixture-v1", md5sum=None):
    return {
        "file_id": file_id,
        "revision": revision,
        "md5sum": md5sum or hashlib.md5(payload).hexdigest(),
    }


class _FakeGdcResponse:
    headers = {}

    def __init__(self, payload):
        self.payload = payload

    def raise_for_status(self):
        pass

    def iter_content(self, chunk_size):
        yield self.payload

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def _write_tiny_model(path):
    np.savez_compressed(
        path,
        selected_genes=np.array(["ENSG1.1", "ENSG2.1"]),
        scaler_mean=np.array([0.0, 0.0]),
        scaler_scale=np.array([1.0, 1.0]),
        coef=np.array([1.0, -0.5]),
        intercept=np.array(0.1),
        class_order=np.array([0, 1]),
    )


def _write_matrix_cache(path, matrix, fingerprint):
    pv.write_dataframe_cache(
        path,
        matrix,
        fingerprint=fingerprint,
        fingerprint_inputs={"test": fingerprint},
        cache_kind="test_matrix",
    )


def test_cache_fingerprint_changes_with_sample_ids():
    a = gx.cache_fingerprint(["s1", "s2"], ["ENSG1", "ENSG2"], "http://example/matrix")
    b = gx.cache_fingerprint(["s1", "s3"], ["ENSG1", "ENSG2"], "http://example/matrix")
    assert a != b


def test_cache_fingerprint_changes_with_selected_genes():
    a = gx.cache_fingerprint(["s1", "s2"], ["ENSG1", "ENSG2"], "http://example/matrix")
    b = gx.cache_fingerprint(["s1", "s2"], ["ENSG1", "ENSG3"], "http://example/matrix")
    assert a != b


def test_cache_fingerprint_changes_with_url():
    a = gx.cache_fingerprint(["s1"], ["ENSG1"], "http://example/a")
    b = gx.cache_fingerprint(["s1"], ["ENSG1"], "http://example/b")
    assert a != b


def test_cache_fingerprint_changes_with_model_source_and_parser(monkeypatch):
    base = gx.cache_fingerprint(
        ["s1"], ["ENSG1"], "http://example/a", "model-a",
        {"dataset_id": "dataset", "revision": "one"},
    )
    assert base != gx.cache_fingerprint(
        ["s1"], ["ENSG1"], "http://example/a", "model-b",
        {"dataset_id": "dataset", "revision": "one"},
    )
    assert base != gx.cache_fingerprint(
        ["s1"], ["ENSG1"], "http://example/a", "model-a",
        {"dataset_id": "dataset", "revision": "two"},
    )
    monkeypatch.setattr(gx, "XENA_MATRIX_PARSER_VERSION", "future-parser")
    assert base != gx.cache_fingerprint(
        ["s1"], ["ENSG1"], "http://example/a", "model-a",
        {"dataset_id": "dataset", "revision": "one"},
    )


def test_cache_fingerprint_stable_for_identical_inputs():
    a = gx.cache_fingerprint(["s1", "s2"], ["ENSG1", "ENSG2"], "http://example/matrix")
    b = gx.cache_fingerprint(["s1", "s2"], ["ENSG1", "ENSG2"], "http://example/matrix")
    assert a == b


def test_load_cached_matrix_rejects_stale_fingerprint(tmp_path):
    cache_path = tmp_path / "matrix.parquet"
    matrix = pd.DataFrame({"ENSG1": [1.0]}, index=["s1"])
    _write_matrix_cache(cache_path, matrix, "old-fingerprint")

    result = gx._load_cached_matrix(cache_path, "new-fingerprint")

    assert result is None


def test_load_cached_matrix_accepts_matching_fingerprint(tmp_path):
    cache_path = tmp_path / "matrix.parquet"
    matrix = pd.DataFrame({"ENSG1": [1.0]}, index=["s1"])
    _write_matrix_cache(cache_path, matrix, "abc")

    result = gx._load_cached_matrix(cache_path, "abc", ["s1"], ["ENSG1"])

    assert result is not None
    assert cache_path.read_bytes()[:4] == b"PAR1"
    pd.testing.assert_frame_equal(result, matrix)


def test_load_cached_matrix_missing_meta_treated_as_stale(tmp_path):
    cache_path = tmp_path / "matrix.pkl"
    cache_path.write_bytes(b"legacy pickle cache must never be read")

    result = gx._load_cached_matrix(cache_path, "anything")

    assert result is None


def test_legacy_pickle_cache_is_never_unpickled(tmp_path):
    marker = tmp_path / "executed.txt"

    class Malicious:
        def __reduce__(self):
            expression = f"open({str(marker)!r}, 'w').write('executed')"
            return eval, (expression,)

    cache_path = tmp_path / "legacy.pkl"
    cache_path.write_bytes(pickle.dumps(Malicious()))
    gx._cache_meta_path(cache_path).write_text(
        json.dumps({"fingerprint": "abc"}), encoding="utf-8"
    )

    assert gx._load_cached_matrix(cache_path, "abc") is None
    assert not marker.exists()


def test_new_cache_writer_rejects_pickle_extension(tmp_path):
    with pytest.raises(ValueError, match="expected .parquet"):
        _write_matrix_cache(
            tmp_path / "unsafe.pkl",
            pd.DataFrame({"ENSG1": [1.0]}, index=["s1"]),
            "abc",
        )
    assert not (tmp_path / "unsafe.pkl").exists()


def test_load_cached_matrix_handles_corrupt_parquet(tmp_path):
    cache_path = tmp_path / "matrix.parquet"
    _write_matrix_cache(cache_path, pd.DataFrame({"ENSG1": [1.0]}, index=["s1"]), "abc")
    cache_path.write_bytes(b"not a pickle")

    result = gx._load_cached_matrix(cache_path, "abc")

    assert result is None


def test_load_cached_matrix_rejects_axis_tampering_even_with_updated_byte_hash(tmp_path):
    cache_path = tmp_path / "matrix.parquet"
    _write_matrix_cache(cache_path, pd.DataFrame({"ENSG1": [1.0]}, index=["s1"]), "abc")
    pd.DataFrame({"ENSG1": [1.0]}, index=["wrong-sample"]).to_parquet(cache_path)
    metadata = json.loads(gx._cache_meta_path(cache_path).read_text(encoding="utf-8"))
    metadata["content_sha256"] = pv.sha256_file(cache_path)
    gx._cache_meta_path(cache_path).write_text(json.dumps(metadata), encoding="utf-8")

    assert gx._load_cached_matrix(cache_path, "abc", ["s1"], ["ENSG1"]) is None


def test_load_cached_matrix_accepts_ordered_selected_gene_subset(tmp_path):
    cache_path = tmp_path / "matrix.parquet"
    matrix = pd.DataFrame({"ENSG1": [1.0], "ENSG3": [3.0]}, index=["s1"])
    _write_matrix_cache(cache_path, matrix, "abc")

    result = gx._load_cached_matrix(
        cache_path, "abc", ["s1"], ["ENSG1", "ENSG2", "ENSG3"]
    )

    pd.testing.assert_frame_equal(result, matrix)


def test_atomic_write_bytes_leaves_no_tmp_file_behind(tmp_path):
    path = tmp_path / "out.csv"
    gx.atomic_write_bytes(path, lambda p: p.write_text("data", encoding="utf-8"))

    assert path.read_text(encoding="utf-8") == "data"
    assert not path.with_suffix(path.suffix + ".tmp").exists()
    assert not list(tmp_path.glob(".out.csv.*.tmp"))


def test_atomic_csv_is_utf8_with_lf_newlines(tmp_path):
    path = tmp_path / "table.csv"
    pv.atomic_write_csv(path, pd.DataFrame({"value": ["α", "β"]}))
    payload = path.read_bytes()
    assert b"\r\n" not in payload
    assert payload.count(b"\n") == 3


def test_require_complete_merge_accepts_complete_merge():
    predictions = pd.DataFrame({"sample": ["a", "b"], "tumor_probability": [0.1, 0.9]})
    gx.require_complete_merge(predictions, 2, "tumor_probability", "[test]")


def test_require_complete_merge_rejects_row_count_mismatch():
    predictions = pd.DataFrame({"sample": ["a", "b", "b"], "tumor_probability": [0.1, 0.9, 0.9]})
    with pytest.raises(ValueError, match="merge produced 3 rows, expected 2"):
        gx.require_complete_merge(predictions, 2, "tumor_probability", "[test]")


def test_require_complete_merge_rejects_unmatched_scores():
    predictions = pd.DataFrame({"sample": ["a", "b"], "tumor_probability": [0.1, np.nan]})
    with pytest.raises(ValueError, match="1/2 samples have no score"):
        gx.require_complete_merge(predictions, 2, "tumor_probability", "[test]")


def test_extract_matrix_from_xena_uses_cache_when_fingerprint_matches(tmp_path, monkeypatch):
    cache_path = tmp_path / "matrix.parquet"
    sample_ids = ["s1", "s2"]
    selected_genes = ["ENSG1", "ENSG2"]
    matrix = pd.DataFrame({"ENSG1": [1.0, 2.0], "ENSG2": [3.0, 4.0]}, index=sample_ids)
    fingerprint = gx.cache_fingerprint(sample_ids, selected_genes, gx.GTEX_TPM_URL)
    _write_matrix_cache(cache_path, matrix, fingerprint)

    def fail_if_called(*args, **kwargs):
        raise AssertionError("should not hit the network when cache is valid")

    monkeypatch.setattr(gx.requests, "get", fail_if_called)

    result = gx.extract_matrix_from_xena(sample_ids, selected_genes, cache_path, refresh=False)

    pd.testing.assert_frame_equal(result, matrix)


def test_extract_matrix_from_xena_ignores_cache_for_different_samples(tmp_path, monkeypatch):
    """The bug this test guards against: rerunning with a different sample
    set (e.g. a different --seed) must not silently reuse a cache built for
    a different cohort."""
    cache_path = tmp_path / "matrix.parquet"
    old_sample_ids = ["s1", "s2"]
    new_sample_ids = ["s3", "s4"]
    selected_genes = ["ENSG1"]
    old_matrix = pd.DataFrame({"ENSG1": [1.0, 2.0]}, index=old_sample_ids)
    stale_fingerprint = gx.cache_fingerprint(old_sample_ids, selected_genes, gx.GTEX_TPM_URL)
    _write_matrix_cache(cache_path, old_matrix, stale_fingerprint)

    called = {}

    def fake_get(url, stream=True, timeout=120):
        called["hit"] = True
        raise RuntimeError("network stub: extraction attempted, as expected")

    monkeypatch.setattr(gx.requests, "get", fake_get)

    with pytest.raises(RuntimeError, match="network stub"):
        gx.extract_matrix_from_xena(
            new_sample_ids, selected_genes, cache_path, refresh=False,
            source_identity={"dataset_id": "fixture", "revision": "fixture-v1"},
        )

    assert called.get("hit"), "expected a re-extraction attempt for the new sample set"


class _FakeXenaResponse:
    def __init__(self, text):
        self.raw = io.BytesIO(gzip.compress(text.encode("utf-8")))
        self.headers = {"ETag": '"fixture-etag"', "Content-Length": str(len(text))}

    def raise_for_status(self):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def test_xena_all_nan_selected_gene_is_rejected_before_cache_write(
    tmp_path, monkeypatch
):
    payload = (
        "sample\ts1\ts2\n"
        "ENSG_PRESENT.1\t0.0\t1.0\n"
        "ENSG_PRESENT_2.1\t1.0\t2.0\n"
        "ENSG_PRESENT_3.1\t2.0\t3.0\n"
        "ENSG_ALL_NAN.1\tnan\tnan\n"
    )
    monkeypatch.setattr(
        gx.requests, "get", lambda *args, **kwargs: _FakeXenaResponse(payload)
    )
    selected = [
        "ENSG_PRESENT.1", "ENSG_PRESENT_2.1", "ENSG_PRESENT_3.1",
        "ENSG_ABSENT.1", "ENSG_ALL_NAN.1",
    ]
    with pytest.raises(ValueError, match="invalid expression values"):
        gx.extract_matrix_from_xena(
            ["s1", "s2"],
            selected,
            tmp_path / "matrix.parquet",
            refresh=True,
            model_sha256="model-sha",
            source_identity={"dataset_id": "fixture", "revision": "1"},
        )
    assert not (tmp_path / "matrix.parquet").exists()


def test_scored_dataframe_preserves_raw_probability_precision():
    probability = np.nextafter(0.5, 1.0)
    scored = pv.scored_dataframe(["s1"], [probability], 0.5)

    assert scored.loc[0, "tumor_probability"] == probability
    assert scored.loc[0, "call"] == "tumor"


@pytest.mark.parametrize("sample_ids", [[""], [" padded"], ["padded "], ["dup", "dup"], [None]])
def test_scored_dataframe_rejects_invalid_sample_ids(sample_ids):
    with pytest.raises(ValueError, match="sample ID"):
        pv.scored_dataframe(sample_ids, [0.5] * len(sample_ids), 0.5)


@pytest.mark.parametrize("threshold", [np.nan, np.inf, -0.01, 1.01, "invalid"])
def test_scored_dataframe_rejects_invalid_threshold(threshold):
    with pytest.raises(ValueError, match="threshold"):
        pv.scored_dataframe(["s1"], [0.5], threshold)


@pytest.mark.parametrize("probability", [np.nan, np.inf, -0.01, 1.01])
def test_scored_dataframe_rejects_invalid_probability(probability):
    with pytest.raises(ValueError, match="probabilities"):
        pv.scored_dataframe(["s1"], [probability], 0.5)


@pytest.mark.parametrize("value", [np.nan, np.inf, -np.inf, np.float64(np.nan)])
def test_json_provenance_rejects_nonfinite_values(tmp_path, value):
    with pytest.raises(ValueError, match="must not contain"):
        pv.canonical_json_bytes({"bad": value})
    path = tmp_path / "metadata.json"
    with pytest.raises(ValueError, match="must not contain"):
        pv.atomic_write_json(path, {"bad": value})
    assert not path.exists()


def test_cache_metadata_reader_rejects_nonfinite_json(tmp_path):
    cache_path = tmp_path / "matrix.parquet"
    cache_path.write_bytes(b"placeholder")
    pv.cache_meta_path(cache_path).write_text(
        '{"fingerprint": NaN}', encoding="utf-8"
    )
    assert pv.read_cache_metadata(cache_path) is None


def test_extract_selected_genes_raises_on_empty_successful_download(tmp_path, monkeypatch):
    """A GDC response that parses cleanly but matches zero selected genes
    (e.g. a corrupted/withdrawn file) must not be silently cached as valid."""

    payload = b"gene_id\ttpm_unstranded\nENSG_NOT_SELECTED.1\t5.0\n"
    file_id = GDC_FILE_IDS[0]
    monkeypatch.setattr(
        cg.requests, "get", lambda *a, **k: _FakeGdcResponse(payload)
    )

    with pytest.raises(ValueError, match="No selected genes matched"):
        cg.extract_selected_genes(
            file_id, ["ENSG_SELECTED"], tmp_path, retries=1,
            source_identity=_gdc_source(file_id, payload),
        )

    assert not cg.gene_cache_path(tmp_path, file_id).exists()


def test_extract_selected_genes_rejects_version_strip_collision(tmp_path, monkeypatch):
    payload = (
        b"gene_id\ttpm_unstranded\nENSG_SELECTED.1\t5.0\n"
        b"ENSG_SELECTED.2\t9.0\n"
    )
    file_id = GDC_FILE_IDS[1]
    monkeypatch.setattr(
        cg.requests, "get", lambda *a, **k: _FakeGdcResponse(payload)
    )

    with pytest.raises(ValueError, match="multiple rows mapping to selected gene"):
        cg.extract_selected_genes(
            file_id, ["ENSG_SELECTED"], tmp_path, retries=1,
            source_identity=_gdc_source(file_id, payload),
        )
    assert not cg.gene_cache_path(tmp_path, file_id).exists()


@pytest.mark.parametrize("tpm", ["nan", "inf", "-0.1", "not-a-number"])
def test_cptac_parser_rejects_nonfinite_negative_or_nonnumeric_tpm(
    tmp_path, monkeypatch, tpm
):
    payload = f"gene_id\ttpm_unstranded\nENSG1.1\t{tpm}\n".encode()
    file_id = GDC_FILE_IDS[2]
    monkeypatch.setattr(
        cg.requests, "get", lambda *args, **kwargs: _FakeGdcResponse(payload)
    )
    with pytest.raises(ValueError, match="tpm_unstranded"):
        cg.extract_selected_genes(
            file_id, ["ENSG1.1"], tmp_path, retries=1,
            source_identity=_gdc_source(file_id, payload),
        )
    assert not cg.gene_cache_path(tmp_path, file_id).exists()


def test_md5_helper_falls_back_when_usedforsecurity_keyword_is_unavailable(monkeypatch):
    real_md5 = hashlib.md5
    calls = []

    def compatibility_md5(*args, **kwargs):
        calls.append(kwargs)
        if "usedforsecurity" in kwargs:
            raise TypeError("keyword unsupported")
        return real_md5(*args)

    monkeypatch.setattr(cg.hashlib, "md5", compatibility_md5)
    digest = cg._new_md5_digest()
    digest.update(b"fixture")

    assert digest.hexdigest() == real_md5(b"fixture").hexdigest()
    assert calls == [{"usedforsecurity": False}, {}]


def test_cptac_gene_cache_fingerprints_genes_model_and_source(tmp_path, monkeypatch):
    calls = []
    payload = b"gene_id\ttpm_unstranded\nENSG1.1\t3.0\nENSG2.1\t7.0\n"
    file_id = GDC_FILE_IDS[3]

    class FakeResponse(_FakeGdcResponse):
        headers = {"ETag": '"gdc-fixture"'}

        def __enter__(self):
            calls.append("network")
            return self

    monkeypatch.setattr(
        cg.requests, "get", lambda *args, **kwargs: FakeResponse(payload)
    )
    source = {
        **_gdc_source(file_id, payload, revision="1"),
        "provider_token": "abc",
    }
    first = cg.extract_selected_genes(
        file_id, ["ENSG1.1"], tmp_path, retries=1,
        model_sha256="model-a", source_identity=source,
    )
    second = cg.extract_selected_genes(
        file_id, ["ENSG1.1"], tmp_path, retries=1,
        model_sha256="model-a", source_identity=source,
    )
    assert len(calls) == 1
    pd.testing.assert_series_equal(first, second)

    changed = cg.extract_selected_genes(
        file_id, ["ENSG2.1"], tmp_path, retries=1,
        model_sha256="model-b",
        source_identity={**source, "revision": "2"},
    )
    assert len(calls) == 2
    assert changed.index.tolist() == ["ENSG2.1"]

    cg.extract_selected_genes(
        file_id, ["ENSG2.1"], tmp_path, retries=1, refresh=True,
        model_sha256="model-b", source_identity={**source, "revision": "2"},
    )
    assert len(calls) == 3


def test_cptac_gene_cache_rejects_content_tampering(tmp_path, monkeypatch):
    payload = b"gene_id\ttpm_unstranded\nENSG1.1\t3.0\n"
    file_id = GDC_FILE_IDS[4]
    calls = {"n": 0}

    def get(*args, **kwargs):
        calls["n"] += 1
        return _FakeGdcResponse(payload)

    monkeypatch.setattr(cg.requests, "get", get)
    source = _gdc_source(file_id, payload)
    cg.extract_selected_genes(
        file_id, ["ENSG1.1"], tmp_path, retries=1, source_identity=source
    )
    cg.gene_cache_path(tmp_path, file_id).write_bytes(b"tampered")
    cg.extract_selected_genes(
        file_id, ["ENSG1.1"], tmp_path, retries=1, source_identity=source
    )
    assert calls["n"] == 2


def test_cptac_download_verifies_provider_md5(tmp_path, monkeypatch):
    payload = b"gene_id\ttpm_unstranded\nENSG1.1\t3.0\n"
    file_id = GDC_FILE_IDS[5]

    class FakeResponse:
        headers = {}

        def raise_for_status(self):
            pass

        def iter_content(self, chunk_size):
            yield payload[:10]
            yield payload[10:]

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    monkeypatch.setattr(cg.requests, "get", lambda *args, **kwargs: FakeResponse())
    expected_md5 = hashlib.md5(payload).hexdigest()
    series = cg.extract_selected_genes(
        file_id, ["ENSG1.1"], tmp_path, retries=1,
        source_identity={
            "file_id": file_id, "md5sum": expected_md5, "revision": "fixture-v1"
        },
    )

    assert series["ENSG1.1"] == pytest.approx(2.0)
    metadata = pv.read_cache_metadata(cg.gene_cache_path(tmp_path, file_id))
    assert metadata["extra"]["downloaded_md5"] == expected_md5


def test_cptac_download_rejects_provider_md5_mismatch(tmp_path, monkeypatch):
    payload = b"gene_id\ttpm_unstranded\nENSG1.1\t3.0\n"
    file_id = GDC_FILE_IDS[6]

    class FakeResponse:
        headers = {}

        def raise_for_status(self):
            pass

        def iter_content(self, chunk_size):
            yield payload

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    monkeypatch.setattr(cg.requests, "get", lambda *args, **kwargs: FakeResponse())
    with pytest.raises(cg.SourceIntegrityError, match="MD5 mismatch"):
        cg.extract_selected_genes(
            file_id, ["ENSG1.1"], tmp_path, retries=1,
            source_identity={
                "file_id": file_id, "md5sum": "0" * 32,
                "revision": "fixture-v1",
            },
        )
    assert not cg.gene_cache_path(tmp_path, file_id).exists()


def test_cptac_download_rejects_malformed_provider_md5_before_network(tmp_path, monkeypatch):
    file_id = GDC_FILE_IDS[7]
    monkeypatch.setattr(
        cg.requests,
        "get",
        lambda *args, **kwargs: pytest.fail("network must not be called"),
    )
    with pytest.raises(ValueError, match="valid provider MD5"):
        cg.extract_selected_genes(
            file_id, ["ENSG1.1"], tmp_path, retries=1,
            source_identity={
                "file_id": file_id, "md5sum": "not-an-md5",
                "revision": "fixture-v1",
            },
        )


def test_cptac_expression_cache_reuses_valid_matrix_and_refreshes(tmp_path, monkeypatch):
    files = pd.DataFrame({
        "file_id": GDC_FILE_IDS[8:10],
        "file_name": ["one.tsv", "two.tsv"],
        "file_size": [10, 20],
        "md5sum": ["a" * 32, "b" * 32],
    })
    calls = []

    def fake_extract(file_id, selected_genes, cache_dir, retries, refresh,
                     model_sha256, source_identity, offline):
        calls.append((file_id, refresh))
        return pd.Series({"ENSG1": 1.0, "ENSG2": 2.0}, name=file_id)

    monkeypatch.setattr(cg, "extract_selected_genes", fake_extract)
    matrix_path = tmp_path / "expression.parquet"
    first = cg.build_expression_matrix(
        files, ["ENSG1", "ENSG2"], tmp_path / "genes", 2, 1,
        matrix_cache_path=matrix_path, model_sha256="model", source_revision="fixture-v1",
    )
    second = cg.build_expression_matrix(
        files, ["ENSG1", "ENSG2"], tmp_path / "genes", 2, 1,
        matrix_cache_path=matrix_path, model_sha256="model", source_revision="fixture-v1",
    )
    assert len(calls) == 2
    pd.testing.assert_frame_equal(first, second)

    cg.build_expression_matrix(
        files, ["ENSG1", "ENSG2"], tmp_path / "genes", 2, 1,
        matrix_cache_path=matrix_path, refresh=True, model_sha256="model",
        source_revision="fixture-v1",
    )
    assert len(calls) == 4
    assert all(refresh for _, refresh in calls[-2:])


def test_cptac_sampling_deduplicates_biospecimens_and_cases():
    rows = [
        (GDC_FILE_IDS[10], "normal-1", "case-1", "Solid Tissue Normal"),
        (GDC_FILE_IDS[11], "normal-1", "case-1", "Solid Tissue Normal"),
        (GDC_FILE_IDS[12], "normal-2", "case-2", "Solid Tissue Normal"),
        (GDC_FILE_IDS[13], "tumor-1", "case-1", "Primary Tumor"),
        (GDC_FILE_IDS[14], "tumor-2", "case-3", "Primary Tumor"),
        (GDC_FILE_IDS[15], "tumor-3", "case-4", "Primary Tumor"),
    ]
    manifest = pd.DataFrame(rows, columns=[
        "file_id", "sample_submitter_id", "case_submitter_id", "sample_type"
    ])
    manifest["project"] = "CPTAC-3"
    manifest["file_name"] = [f"{index}.tsv" for index in range(len(manifest))]
    manifest["file_size"] = 100
    manifest["md5sum"] = [f"{index:x}" * 32 for index in range(len(manifest))]
    manifest["n_sample_types"] = 1

    sampled = cg.choose_files(manifest, "CPTAC-3", 2, 42)

    assert len(sampled) == 4
    assert sampled["sample_submitter_id"].is_unique
    assert sampled["case_submitter_id"].is_unique
    assert sampled["label"].value_counts().to_dict() == {0: 2, 1: 2}


def test_tcga_sampling_avoids_same_donor_across_classes():
    phenotype = pd.DataFrame({
        "sample": ["TCGA-AA-0001-11", "TCGA-BB-0002-11", "TCGA-AA-0001-01",
                   "TCGA-CC-0003-01", "TCGA-DD-0004-01"],
        "_study": ["TCGA"] * 5,
        "_sample_type": ["Solid Tissue Normal", "Solid Tissue Normal",
                         "Primary Tumor", "Primary Tumor", "Primary Tumor"],
        "_primary_site": ["A", "B", "A", "C", "D"],
    })

    sampled = tx.choose_tcga_samples(phenotype, 2, 42)

    assert len(sampled) == 4
    assert sampled["donor_id"].is_unique


def test_gtex_sampling_can_enforce_unique_donors_across_sites():
    phenotype = pd.DataFrame({
        "sample": ["GTEX-A-0001", "GTEX-B-0001", "GTEX-A-0002", "GTEX-C-0001"],
        "_study": ["GTEX"] * 4,
        "_sample_type": ["Normal Tissue"] * 4,
        "_primary_site": ["Site1", "Site1", "Site2", "Site2"],
    })

    sampled = gx.choose_gtex_samples(phenotype, 1, 1, 42)

    assert len(sampled) == 2
    assert sampled["donor_id"].is_unique


def test_locked_cptac_manifest_rejects_duplicate_biospecimens(tmp_path):
    path = tmp_path / "locked.csv"
    pd.DataFrame({
        "file_id": GDC_FILE_IDS[16:18],
        "project": ["CPTAC-3", "CPTAC-3"],
        "case_submitter_id": ["c1", "c2"],
        "sample_submitter_id": ["sample", "sample"],
        "sample_type": ["Primary Tumor", "Primary Tumor"],
        "label": [1, 1],
    }).to_csv(path, index=False)

    with pytest.raises(ValueError, match="duplicate sample_submitter_id"):
        cg.validate_locked_sample_manifest(path)


def test_locked_xena_manifest_preserves_exact_order_and_audits_donors(tmp_path):
    path = tmp_path / "locked.csv"
    pd.DataFrame({
        "sample": ["GTEX-A-0001", "GTEX-A-0002"],
        "_study": ["GTEX", "GTEX"],
        "_sample_type": ["Normal Tissue", "Normal Tissue"],
        "_primary_site": ["A", "B"],
        "primary disease or tissue": ["A tissue", "B tissue"],
    }).to_csv(path, index=False)

    locked = gx.load_locked_sample_manifest(path, "GTEX")
    audit = pv.group_audit(locked, "sample", "donor_id")

    assert locked["sample"].tolist() == ["GTEX-A-0001", "GTEX-A-0002"]
    assert audit["n_repeated_groups"] == 1


def test_run_manifest_records_hashes_parameters_git_model_and_alignment(tmp_path):
    model = tmp_path / "model.npz"
    source = tmp_path / "source.csv"
    output = tmp_path / "output.csv"
    model.write_bytes(b"model")
    source.write_text("source", encoding="utf-8")
    output.write_text("output", encoding="utf-8")

    manifest = pv.write_run_manifest(
        tmp_path / "run_manifest.json",
        root=ROOT,
        run_kind="unit_test",
        started_at_utc="2026-07-12T00:00:00Z",
        argv=["validator.py", "--seed", "7"],
        parameters={"seed": 7},
        model_path=model,
        sources={"dataset": "fixture"},
        inputs={"source": source},
        outputs={"result": output},
        alignment={"n_matched_genes": 3, "missing_genes": ["g4"]},
        cohort_audit={"n_rows": 2},
    )

    assert manifest["parameters"]["seed"] == 7
    assert manifest["model"]["sha256"] == pv.sha256_file(model)
    assert manifest["inputs"]["source"]["sha256"] == pv.sha256_file(source)
    assert manifest["outputs"]["result"]["sha256"] == pv.sha256_file(output)
    assert "commit" in manifest["git"]
    assert manifest["alignment"]["missing_genes"] == ["g4"]
    assert manifest["source_code"]["provenance"]["sha256"] == pv.sha256_file(
        pv.PROVENANCE_SOURCE_PATH
    )


def test_git_state_fails_closed_on_timeout(monkeypatch):
    def timeout(*args, **kwargs):
        raise subprocess.TimeoutExpired(cmd="git", timeout=5)

    monkeypatch.setattr(pv.subprocess, "run", timeout)
    assert pv.git_state(ROOT) == {"commit": None, "dirty": None}


def test_external_validation_modules_import_without_sklearn_or_requests():
    code = f"""
import importlib.abc, sys
class Block(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if (fullname == 'sklearn' or fullname.startswith('sklearn.') or
                fullname == 'requests' or fullname.startswith('requests.')):
            raise ImportError('optional package intentionally unavailable')
        return None
sys.meta_path.insert(0, Block())
sys.path.insert(0, {EXTERNAL_VALIDATION!r})
sys.path.insert(0, {ROOT!r})
import validate_gtex_xena, validate_tcga_toil_xena, validate_cptac_gdc
"""
    result = subprocess.run(
        [sys.executable, "-c", code], cwd=ROOT, capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr


def test_gdc_gene_cache_uses_hashed_contained_key_and_rejects_path_traversal(tmp_path):
    file_id = GDC_FILE_IDS[0]
    path = cg.gene_cache_path(tmp_path / "gene_cache", file_id)

    assert path.parent == (tmp_path / "gene_cache").resolve()
    assert path.name.endswith(".parquet")
    assert file_id not in path.name
    assert len(path.stem) == 64
    with pytest.raises(ValueError, match="canonical GDC file UUID"):
        cg.gene_cache_path(tmp_path / "gene_cache", "../escape")


@pytest.mark.parametrize(
    ("module", "script", "predictions_name"),
    [
        (gx, "validate_gtex_xena.py", "gtex_predictions.csv"),
        (tx, "validate_tcga_toil_xena.py", "tcga_toil_predictions.csv"),
        (cg, "validate_cptac_gdc.py", "cptac_predictions.csv"),
    ],
)
def test_validator_rejects_locked_input_output_collision_before_other_work(
    tmp_path, monkeypatch, module, script, predictions_name
):
    out_dir = tmp_path / module.__name__
    out_dir.mkdir()
    collision = out_dir / predictions_name
    collision.write_text("protected-input", encoding="utf-8")
    model = tmp_path / f"{module.__name__}.npz"
    _write_tiny_model(model)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            script,
            "--sample-manifest", str(collision),
            "--source-revision", "fixture-v1",
            "--out-dir", str(out_dir),
            "--weights", str(model),
        ],
    )

    with pytest.raises(ValueError, match="protected input .* collides"):
        module.main()

    assert collision.read_text(encoding="utf-8") == "protected-input"
    assert not any(out_dir.glob("*.parquet"))


@pytest.mark.parametrize(
    "bad_column,bad_values,match",
    [
        ("_study", ["GTEX", ""], "non-empty _study"),
        ("_sample_type", ["Normal Tissue", "Primary Tumor"], "unsupported _sample_type"),
        ("_primary_site", ["Site", ""], "blank _primary_site"),
    ],
)
def test_locked_gtex_manifest_rejects_provider_semantic_mismatch(
    tmp_path, bad_column, bad_values, match
):
    frame = pd.DataFrame({
        "sample": ["GTEX-A-0001", "GTEX-B-0001"],
        "_study": ["GTEX", "GTEX"],
        "_sample_type": ["Normal Tissue", "Normal Tissue"],
        "_primary_site": ["Site A", "Site B"],
        "primary disease or tissue": ["Tissue A", "Tissue B"],
    })
    frame[bad_column] = bad_values
    path = tmp_path / "locked.csv"
    frame.to_csv(path, index=False)

    with pytest.raises(ValueError, match=match):
        gx.load_locked_sample_manifest(path, "GTEX")


def test_locked_tcga_manifest_rejects_label_sample_type_mismatch(tmp_path):
    path = tmp_path / "locked.csv"
    pd.DataFrame({
        "sample": ["TCGA-AA-0001-11", "TCGA-BB-0002-01"],
        "_study": ["TCGA", "TCGA"],
        "_sample_type": ["Solid Tissue Normal", "Primary Tumor"],
        "_primary_site": ["Kidney", "Kidney"],
        "label": [1, 0],
    }).to_csv(path, index=False)

    with pytest.raises(ValueError, match="labels do not agree"):
        tx.validate_locked_tcga_sample_manifest(path)


def _locked_cptac_frame(include_md5=True):
    frame = pd.DataFrame({
        "file_id": GDC_FILE_IDS[:2],
        "file_name": ["normal.tsv", "tumor.tsv"],
        "file_size": [100, 200],
        "project": ["CPTAC-3", "CPTAC-3"],
        "case_submitter_id": ["case-normal", "case-tumor"],
        "sample_submitter_id": ["sample-normal", "sample-tumor"],
        "sample_type": ["Solid Tissue Normal", "Primary Tumor"],
        "n_sample_types": [1, 1],
        "label": [0, 1],
    })
    if include_md5:
        frame["md5sum"] = ["0" * 32, "1" * 32]
    return frame


def test_locked_cptac_manifest_requires_project_and_md5_for_live_use(tmp_path):
    frame = _locked_cptac_frame(include_md5=False)
    path = tmp_path / "locked.csv"
    frame.to_csv(path, index=False)

    historical = cg.validate_locked_sample_manifest(
        path, expected_project="CPTAC-3", require_provider_md5=False
    )
    assert "md5sum" not in historical
    with pytest.raises(ValueError, match="historical manifests may be used only"):
        cg.validate_locked_sample_manifest(
            path, expected_project="CPTAC-3", require_provider_md5=True
        )

    frame["project"] = ["TCGA-X", "CPTAC-3"]
    frame.to_csv(path, index=False)
    with pytest.raises(ValueError, match="must contain only project"):
        cg.validate_locked_sample_manifest(
            path, expected_project="CPTAC-3", require_provider_md5=False
        )


def test_committed_cptac_locked_manifest_is_historical_offline_only():
    path = Path(ROOT) / "external-validation" / "cptac_gdc" / "sampled_manifest.csv"

    historical = cg.validate_locked_sample_manifest(
        path, expected_project="CPTAC-3", require_provider_md5=False
    )
    assert "md5sum" not in historical.columns
    with pytest.raises(ValueError, match="historical manifests may be used only"):
        cg.validate_locked_sample_manifest(
            path, expected_project="CPTAC-3", require_provider_md5=True
        )


def test_refresh_locked_cptac_metadata_binds_exact_cohort_and_adds_md5(monkeypatch):
    locked = _locked_cptac_frame(include_md5=False)
    provider = _locked_cptac_frame(include_md5=True)
    monkeypatch.setattr(cg, "query_cptac_manifest", lambda: provider)

    refreshed = cg.refresh_locked_sample_metadata(
        locked, expected_project="CPTAC-3"
    )

    assert refreshed["file_id"].tolist() == locked["file_id"].tolist()
    assert refreshed["md5sum"].tolist() == ["0" * 32, "1" * 32]


def test_new_cptac_download_rejects_missing_md5_before_network(tmp_path, monkeypatch):
    monkeypatch.setattr(
        cg.requests, "get", lambda *args, **kwargs: pytest.fail("network must not be called")
    )
    with pytest.raises(ValueError, match="require a valid provider MD5"):
        cg.extract_selected_genes(
            GDC_FILE_IDS[2], ["ENSG1.1"], tmp_path, retries=1,
            source_identity={"file_id": GDC_FILE_IDS[2], "revision": "fixture-v1"},
        )


def test_xena_version_strip_collision_after_all_genes_is_fatal(tmp_path, monkeypatch):
    payload = (
        "sample\ts1\n"
        "ENSG1.1\t1.0\n"
        "ENSG2.1\t2.0\n"
        "ENSG1.2\t3.0\n"
    )
    monkeypatch.setattr(
        gx.requests, "get", lambda *args, **kwargs: _FakeXenaResponse(payload)
    )
    cache = tmp_path / "matrix.parquet"
    with pytest.raises(ValueError, match="colliding"):
        gx.extract_matrix_from_xena(
            ["s1"], ["ENSG1.9", "ENSG2.9"], cache, refresh=True,
            source_identity={"dataset_id": "fixture", "revision": "fixture-v1"},
        )
    assert not cache.exists()


@pytest.mark.parametrize("value", ["nan", "inf", "-100", "not-a-number"])
def test_xena_rejects_nonfinite_or_out_of_domain_source_values(
    tmp_path, monkeypatch, value
):
    payload = f"sample\ts1\nENSG1.1\t{value}\n"
    monkeypatch.setattr(
        gx.requests, "get", lambda *args, **kwargs: _FakeXenaResponse(payload)
    )
    with pytest.raises(ValueError, match="invalid expression values"):
        gx.extract_matrix_from_xena(
            ["s1"], ["ENSG1.1"], tmp_path / "matrix.parquet", refresh=True,
            source_identity={"dataset_id": "fixture", "revision": "fixture-v1"},
        )


@pytest.mark.parametrize(
    "payload,match",
    [
        ("sample\ts1\ts1\nENSG1.1\t1\t1\n", "duplicate sample IDs"),
        ("sample\ts1\ts2\nENSG1.1\t1\n", "shorter than its header"),
    ],
)
def test_xena_rejects_ambiguous_header_and_short_rows(
    tmp_path, monkeypatch, payload, match
):
    monkeypatch.setattr(
        gx.requests, "get", lambda *args, **kwargs: _FakeXenaResponse(payload)
    )
    samples = ["s1"] if "s1\ts1" in payload else ["s1", "s2"]
    with pytest.raises(ValueError, match=match):
        gx.extract_matrix_from_xena(
            samples, ["ENSG1.1"], tmp_path / "matrix.parquet", refresh=True,
            source_identity={"dataset_id": "fixture", "revision": "fixture-v1"},
        )


def test_offline_cache_miss_never_attempts_network(tmp_path, monkeypatch):
    monkeypatch.setattr(
        gx, "_http_client", lambda: pytest.fail("network must not be called")
    )
    with pytest.raises(ValueError, match="offline/cache-only"):
        gx.extract_matrix_from_xena(
            ["s1"], ["ENSG1.1"], tmp_path / "matrix.parquet", refresh=False,
            source_identity={"dataset_id": "fixture", "revision": "unversioned"},
            offline=True,
        )


def test_live_fetch_requires_concrete_source_revision_before_network(tmp_path, monkeypatch):
    monkeypatch.setattr(
        gx, "_http_client", lambda: pytest.fail("network must not be called")
    )
    with pytest.raises(ValueError, match="provider snapshot/revision"):
        gx.extract_matrix_from_xena(
            ["s1"], ["ENSG1.1"], tmp_path / "matrix.parquet", refresh=True,
        )


def test_publish_staged_files_rolls_back_entire_output_set(tmp_path, monkeypatch):
    stage = tmp_path / "stage"
    final = tmp_path / "final"
    stage.mkdir()
    final.mkdir()
    staged = {
        "a": stage / "a.txt",
        "b": stage / "b.txt",
        "run_manifest": stage / "run_manifest.json",
    }
    destinations = {
        "a": final / "a.txt",
        "b": final / "b.txt",
        "run_manifest": final / "run_manifest.json",
    }
    for name, path in staged.items():
        path.write_text(f"new-{name}", encoding="utf-8")
    for name, path in destinations.items():
        path.write_text(f"old-{name}", encoding="utf-8")
    real_replace = pv.os.replace

    def fail_second_staged_move(source, destination):
        if Path(source) == staged["b"]:
            raise OSError("injected publish failure")
        return real_replace(source, destination)

    monkeypatch.setattr(pv.os, "replace", fail_second_staged_move)
    with pytest.raises(OSError, match="injected publish failure"):
        pv.publish_staged_files(staged, destinations)

    for name, path in destinations.items():
        assert path.read_text(encoding="utf-8") == f"old-{name}"


def test_publish_staged_files_moves_manifest_last(tmp_path, monkeypatch):
    stage = tmp_path / "stage"
    final = tmp_path / "final"
    stage.mkdir()
    final.mkdir()
    staged = {
        "a": stage / "a.txt",
        "run_manifest": stage / "run_manifest.json",
    }
    destinations = {
        "a": final / "a.txt",
        "run_manifest": final / "run_manifest.json",
    }
    for path in staged.values():
        path.write_text("new", encoding="utf-8")
    moved_to_final = []
    real_replace = pv.os.replace

    def record_replace(source, destination):
        destination = Path(destination)
        if destination.parent == final:
            moved_to_final.append(destination.name)
        return real_replace(source, destination)

    monkeypatch.setattr(pv.os, "replace", record_replace)
    pv.publish_staged_files(staged, destinations)

    assert moved_to_final == ["a.txt", "run_manifest.json"]
