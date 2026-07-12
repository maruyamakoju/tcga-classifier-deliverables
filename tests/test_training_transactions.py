"""Generation-level integrity tests for training artifacts."""

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import training_tools as TT


ROOT = Path(__file__).resolve().parents[1]
EXPORT_SCRIPT = ROOT / "cancer-type-classifier" / "export_features_npy.py"
SPEC = importlib.util.spec_from_file_location("export_features_npy", EXPORT_SCRIPT)
EXPORT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(EXPORT)


def write_source(path, value=1.0):
    frame = pd.DataFrame(
        [[value, value + 1], [value + 2, value + 3]],
        index=["sample-1", "sample-2"],
        columns=["gene-1", "gene-2"],
        dtype=np.float64,
    )
    frame.to_pickle(path)


def test_committed_feature_export_lock_is_path_neutral_and_complete():
    lock = TT.load_feature_export_lock()
    assert lock["converter_environment"] == {
        "implementation": "CPython",
        "python": "3.13",
        "numpy": "2.4.1",
        "pandas": "3.0.0",
    }
    assert set(lock["outputs"]) == {
        "features_float32",
        "features_float64",
        "genes",
        "samples",
    }
    records = [lock["source"], *lock["outputs"].values()]
    for record in records:
        assert Path(record["path"]).name == record["path"]
        assert not Path(record["path"]).is_absolute()
        assert record["bytes"] > 0
        assert len(record["sha256"]) == 64


def test_exporter_hash_binds_one_fresh_generation(tmp_path, monkeypatch):
    source = tmp_path / "source.pkl"
    output = tmp_path / "feature-generation"
    write_source(source)
    monkeypatch.setattr(EXPORT.pd, "__version__", "3.0.0")

    manifest = EXPORT.export(
        source,
        output,
        "both",
        trusted_source_pickle=True,
        allow_noncanonical_export=True,
    )

    assert manifest["schema_version"] == "3.0"
    assert manifest["canonical_lock_verified"] is False
    assert manifest["source"]["path"] == source.name
    assert set(manifest["outputs"]) == {
        "features_float32",
        "features_float64",
        "genes",
        "samples",
    }
    TT.validate_feature_manifest(
        output / "X_full.npy",
        output / "X_genes.npy",
        output / "X_samples.npy",
        expected_dtype=np.float32,
        allow_unverified=False,
        require_canonical_lock=False,
    )
    TT.validate_feature_manifest(
        output / "X_full_float64.npy",
        output / "X_genes.npy",
        output / "X_samples.npy",
        expected_dtype=np.float64,
        allow_unverified=False,
        require_canonical_lock=False,
    )
    with pytest.raises(ValueError, match="canonical lock"):
        TT.validate_feature_manifest(
            output / "X_full.npy",
            output / "X_genes.npy",
            output / "X_samples.npy",
            expected_dtype=np.float32,
            allow_unverified=False,
        )
    with pytest.raises(ValueError, match="--force"):
        EXPORT.export(
            source,
            output,
            "both",
            trusted_source_pickle=True,
            allow_noncanonical_export=True,
        )
    with (output / "X_full.npy").open("ab") as handle:
        handle.write(b"tamper")
    with pytest.raises(ValueError, match="hash mismatch"):
        TT.validate_feature_manifest(
            output / "X_full.npy",
            output / "X_genes.npy",
            output / "X_samples.npy",
            expected_dtype=np.float32,
            allow_unverified=False,
            require_canonical_lock=False,
        )


def test_exporter_failure_before_commit_preserves_previous_generation(
    tmp_path, monkeypatch
):
    source = tmp_path / "source.pkl"
    output = tmp_path / "feature-generation"
    write_source(source, value=1.0)
    monkeypatch.setattr(EXPORT.pd, "__version__", "3.0.0")
    EXPORT.export(
        source,
        output,
        "both",
        trusted_source_pickle=True,
        allow_noncanonical_export=True,
    )
    before = {
        name: TT.sha256_file(output / name)
        for name in EXPORT.MANAGED_EXPORT_NAMES
    }
    write_source(source, value=100.0)

    def injected_failure(*args, **kwargs):
        raise OSError("injected commit failure")

    monkeypatch.setattr(EXPORT, "commit_file_generation", injected_failure)
    with pytest.raises(OSError, match="injected"):
        EXPORT.export(
            source,
            output,
            "both",
            force=True,
            trusted_source_pickle=True,
            allow_noncanonical_export=True,
        )
    after = {
        name: TT.sha256_file(output / name)
        for name in EXPORT.MANAGED_EXPORT_NAMES
    }
    assert after == before


def test_file_generation_rolls_back_failure_after_data_move(tmp_path, monkeypatch):
    output = tmp_path / "output"
    stage = tmp_path / "stage"
    output.mkdir()
    stage.mkdir()
    (output / "data.bin").write_bytes(b"old-data")
    (output / "manifest.json").write_text("old-manifest", encoding="utf-8")
    (stage / "data.bin").write_bytes(b"new-data")
    (stage / "manifest.json").write_text("new-manifest", encoding="utf-8")
    real_replace = TT.os.replace

    def fail_manifest_promotion(source, destination):
        if Path(source).parent == stage and Path(source).name == "manifest.json":
            raise OSError("injected manifest promotion failure")
        return real_replace(source, destination)

    monkeypatch.setattr(TT.os, "replace", fail_manifest_promotion)
    with pytest.raises(OSError, match="injected"):
        TT.commit_file_generation(
            stage,
            output,
            staged_names=["data.bin", "manifest.json"],
            managed_names=["data.bin", "manifest.json"],
            manifest_name="manifest.json",
            force=True,
        )
    assert (output / "data.bin").read_bytes() == b"old-data"
    assert (output / "manifest.json").read_text(encoding="utf-8") == "old-manifest"


def test_staged_output_failure_and_force_have_generation_semantics(tmp_path):
    output = tmp_path / "run"
    output.mkdir()
    (output / "old.txt").write_text("old", encoding="utf-8")

    with pytest.raises(RuntimeError, match="injected"):
        with TT.staged_output_directory(
            output, force=True, protected_paths=[]
        ) as stage:
            (stage / "new.txt").write_text("new", encoding="utf-8")
            raise RuntimeError("injected")
    assert (output / "old.txt").read_text(encoding="utf-8") == "old"

    with TT.staged_output_directory(output, force=True, protected_paths=[]) as stage:
        (stage / "new.txt").write_text("new", encoding="utf-8")
    assert not (output / "old.txt").exists()
    assert (output / "new.txt").read_text(encoding="utf-8") == "new"


def test_input_snapshot_and_output_isolation_fail_closed(tmp_path):
    source = tmp_path / "input.bin"
    source.write_bytes(b"before")
    snapshot = TT.snapshot_inputs({"source": source})
    source.write_bytes(b"after")
    with pytest.raises(ValueError, match="changed during the run"):
        TT.verify_input_snapshot(snapshot)
    with pytest.raises(ValueError, match="must not contain an input"):
        with TT.staged_output_directory(
            tmp_path,
            force=True,
            protected_paths=[source],
        ):
            pass


def test_staged_output_rejects_existing_regular_file_even_with_force(tmp_path):
    output = tmp_path / "important.txt"
    output.write_text("must survive", encoding="utf-8")

    with pytest.raises(ValueError, match="not a directory"):
        with TT.staged_output_directory(output, force=True, protected_paths=[]):
            pass
    assert output.read_text(encoding="utf-8") == "must survive"


def test_staged_output_rejects_symlink_or_reparse_point(tmp_path):
    target = tmp_path / "target"
    target.mkdir()
    link = tmp_path / "linked-output"
    try:
        link.symlink_to(target, target_is_directory=True)
    except OSError as exc:
        pytest.skip(f"directory symlinks unavailable: {exc}")

    with pytest.raises(ValueError, match="symlink or reparse"):
        with TT.staged_output_directory(link, force=True, protected_paths=[]):
            pass


def test_single_writer_lock_rejects_overlapping_generation(tmp_path):
    output = tmp_path / "generation"
    with TT.exclusive_output_lock(output):
        with pytest.raises(ValueError, match="another process"):
            with TT.staged_output_directory(
                output, force=False, protected_paths=[]
            ):
                pass


def test_exporter_requires_explicit_pickle_trust(tmp_path):
    source = tmp_path / "source.pkl"
    write_source(source)
    with pytest.raises(ValueError, match="loading a pickle can execute code"):
        EXPORT.export(
            source,
            tmp_path / "output",
            "both",
            allow_noncanonical_export=True,
        )


def test_exporter_rejects_source_managed_destination_before_unpickle(
    tmp_path, monkeypatch
):
    output = tmp_path / "output"
    output.mkdir()
    source = output / "X_full.npy"
    write_source(source)
    monkeypatch.setattr(EXPORT.pd, "__version__", "3.0.0")

    def forbidden_load(*args, **kwargs):
        raise AssertionError("collision must fail before pickle.load")

    monkeypatch.setattr(EXPORT.pickle, "load", forbidden_load)
    with pytest.raises(ValueError, match="collide"):
        EXPORT.export(
            source,
            output,
            "float32",
            force=True,
            trusted_source_pickle=True,
            allow_noncanonical_export=True,
        )
    assert source.is_file()
