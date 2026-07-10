"""Unit tests for external-validation cache/merge-integrity helpers.

Pure-logic tests only -- no real network calls. external-validation/ is not
part of the lightweight release (not in release_tools.common.RELEASE_FILES),
so this intentionally does not attempt to exercise the live GDC/Xena download
paths; it covers the cache-fingerprinting and merge-integrity logic added to
fix a real bug (a cache keyed only on a fixed path, with no fingerprint of
sampling parameters/selected genes, could silently return stale data for the
wrong samples after a rerun with different arguments).
"""
import json
import math
import os
import sys

import numpy as np
import pandas as pd
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXTERNAL_VALIDATION = os.path.join(ROOT, "external-validation")
sys.path.insert(0, EXTERNAL_VALIDATION)

import validate_gtex_xena as gx  # noqa: E402
import validate_cptac_gdc as cg  # noqa: E402


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


def test_cache_fingerprint_stable_for_identical_inputs():
    a = gx.cache_fingerprint(["s1", "s2"], ["ENSG1", "ENSG2"], "http://example/matrix")
    b = gx.cache_fingerprint(["s1", "s2"], ["ENSG1", "ENSG2"], "http://example/matrix")
    assert a == b


def test_load_cached_matrix_rejects_stale_fingerprint(tmp_path):
    cache_path = tmp_path / "matrix.pkl"
    matrix = pd.DataFrame({"ENSG1": [1.0]}, index=["s1"])
    matrix.to_pickle(cache_path)
    gx._cache_meta_path(cache_path).write_text(
        json.dumps({"fingerprint": "old-fingerprint"}), encoding="utf-8"
    )

    result = gx._load_cached_matrix(cache_path, "new-fingerprint")

    assert result is None


def test_load_cached_matrix_accepts_matching_fingerprint(tmp_path):
    cache_path = tmp_path / "matrix.pkl"
    matrix = pd.DataFrame({"ENSG1": [1.0]}, index=["s1"])
    matrix.to_pickle(cache_path)
    gx._cache_meta_path(cache_path).write_text(
        json.dumps({"fingerprint": "abc"}), encoding="utf-8"
    )

    result = gx._load_cached_matrix(cache_path, "abc")

    assert result is not None
    pd.testing.assert_frame_equal(result, matrix)


def test_load_cached_matrix_missing_meta_treated_as_stale(tmp_path):
    cache_path = tmp_path / "matrix.pkl"
    pd.DataFrame({"ENSG1": [1.0]}, index=["s1"]).to_pickle(cache_path)
    # No .meta.json sidecar written (e.g. a cache from before fingerprinting existed).

    result = gx._load_cached_matrix(cache_path, "anything")

    assert result is None


def test_load_cached_matrix_handles_corrupt_pickle(tmp_path):
    cache_path = tmp_path / "matrix.pkl"
    cache_path.write_bytes(b"not a pickle")
    gx._cache_meta_path(cache_path).write_text(
        json.dumps({"fingerprint": "abc"}), encoding="utf-8"
    )

    result = gx._load_cached_matrix(cache_path, "abc")

    assert result is None


def test_atomic_write_bytes_leaves_no_tmp_file_behind(tmp_path):
    path = tmp_path / "out.csv"
    gx.atomic_write_bytes(path, lambda p: p.write_text("data", encoding="utf-8"))

    assert path.read_text(encoding="utf-8") == "data"
    assert not path.with_suffix(path.suffix + ".tmp").exists()


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
    cache_path = tmp_path / "matrix.pkl"
    sample_ids = ["s1", "s2"]
    selected_genes = ["ENSG1", "ENSG2"]
    matrix = pd.DataFrame({"ENSG1": [1.0, 2.0], "ENSG2": [3.0, 4.0]}, index=sample_ids)
    matrix.to_pickle(cache_path)
    fingerprint = gx.cache_fingerprint(sample_ids, selected_genes, gx.GTEX_TPM_URL)
    gx._cache_meta_path(cache_path).write_text(
        json.dumps({"fingerprint": fingerprint}), encoding="utf-8"
    )

    def fail_if_called(*args, **kwargs):
        raise AssertionError("should not hit the network when cache is valid")

    monkeypatch.setattr(gx.requests, "get", fail_if_called)

    result = gx.extract_matrix_from_xena(sample_ids, selected_genes, cache_path, refresh=False)

    pd.testing.assert_frame_equal(result, matrix)


def test_extract_matrix_from_xena_ignores_cache_for_different_samples(tmp_path, monkeypatch):
    """The bug this test guards against: rerunning with a different sample
    set (e.g. a different --seed) must not silently reuse a cache built for
    a different cohort."""
    cache_path = tmp_path / "matrix.pkl"
    old_sample_ids = ["s1", "s2"]
    new_sample_ids = ["s3", "s4"]
    selected_genes = ["ENSG1"]
    old_matrix = pd.DataFrame({"ENSG1": [1.0, 2.0]}, index=old_sample_ids)
    old_matrix.to_pickle(cache_path)
    stale_fingerprint = gx.cache_fingerprint(old_sample_ids, selected_genes, gx.GTEX_TPM_URL)
    gx._cache_meta_path(cache_path).write_text(
        json.dumps({"fingerprint": stale_fingerprint}), encoding="utf-8"
    )

    called = {}

    def fake_get(url, stream=True, timeout=120):
        called["hit"] = True
        raise RuntimeError("network stub: extraction attempted, as expected")

    monkeypatch.setattr(gx.requests, "get", fake_get)

    with pytest.raises(RuntimeError, match="network stub"):
        gx.extract_matrix_from_xena(new_sample_ids, selected_genes, cache_path, refresh=False)

    assert called.get("hit"), "expected a re-extraction attempt for the new sample set"


def test_extract_selected_genes_raises_on_empty_successful_download(tmp_path, monkeypatch):
    """A GDC response that parses cleanly but matches zero selected genes
    (e.g. a corrupted/withdrawn file) must not be silently cached as valid."""

    class FakeResponse:
        def raise_for_status(self):
            pass

        def iter_lines(self, decode_unicode=True):
            yield "gene_id\ttpm_unstranded"
            yield "ENSG_NOT_SELECTED.1\t5.0"

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    monkeypatch.setattr(cg.requests, "get", lambda *a, **k: FakeResponse())

    with pytest.raises(ValueError, match="No selected genes matched"):
        cg.extract_selected_genes("fake-file-id", ["ENSG_SELECTED"], tmp_path, retries=1)

    assert not (tmp_path / "fake-file-id.pkl").exists()


def test_extract_selected_genes_keeps_first_value_on_collision(tmp_path, monkeypatch):
    class FakeResponse:
        def raise_for_status(self):
            pass

        def iter_lines(self, decode_unicode=True):
            yield "gene_id\ttpm_unstranded"
            yield "ENSG_SELECTED.1\t5.0"
            yield "ENSG_SELECTED.2\t9.0"  # collides with .1 after version-stripping

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    monkeypatch.setattr(cg.requests, "get", lambda *a, **k: FakeResponse())

    series = cg.extract_selected_genes("fake-file-id-2", ["ENSG_SELECTED"], tmp_path, retries=1)

    # First row (.1, tpm=5.0) wins; the colliding second row (.2, tpm=9.0) is skipped.
    assert series["ENSG_SELECTED"] == pytest.approx(math.log2(5.0 + 1.0))
