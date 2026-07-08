#!/usr/bin/env python3
"""Build the lightweight release bundle and zip archive."""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import zipfile
from datetime import datetime
from pathlib import Path

from validate_release_lite import FORBIDDEN_NAMES


ROOT = Path(__file__).resolve().parent
RELEASE_DIR = ROOT / "release-lite"
ZIP_PATH = ROOT / "tcga-tumor-normal-release-lite.zip"
ARTIFACTS_PATH = ROOT / "RELEASE_ARTIFACTS.json"

RELEASE_FILES = [
    "audit_cli_entrypoints.py",
    "audit_lightweight_dependencies.py",
    "audit_release_docs.py",
    "calibrate_threshold.py",
    "check_environment.py",
    "cohort_adapt_score.py",
    "DATA_DICTIONARY.md",
    "deployable_lr_weights.npz",
    "EXECUTIVE_SUMMARY.md",
    "INDEX.md",
    "LICENSE",
    "NOTICE.md",
    "CITATION.cff",
    ".zenodo.json",
    "codemeta.json",
    "tcga_rnaseq/__init__.py",
    "tcga_rnaseq/io.py",
    "tcga_rnaseq/align.py",
    "tcga_rnaseq/score.py",
    "tcga_rnaseq/metrics.py",
    "example_input.csv",
    "example_labels.csv",
    "example_output.csv",
    "explain_scores.py",
    "inspect_expression_input.py",
    "run_tumor_normal_workflow.py",
    "LITERATURE_CHECK.md",
    "MODEL_CARD.md",
    "model_gene_metadata.csv",
    "model_qc_reference.json",
    "per_cancer_type_performance.csv",
    "README.md",
    "RELEASE_BUNDLE.md",
    "RELEASE_METADATA.json",
    "RELEASE_NOTES.md",
    "REPORT.md",
    "REPRODUCIBILITY.md",
    "requirements-light.txt",
    "run_release_acceptance.py",
    "run_safety_tests.py",
    "run_smoke_tests.py",
    "score_tumor_normal.py",
    "test_metrics.csv",
    "top_genes_logreg.csv",
    "top_genes_xgboost.csv",
    "TROUBLESHOOTING.md",
    "USER_GUIDE.md",
    "validate_output_contracts.py",
    "validate_release_lite.py",
    "validate_zip_bundle.py",
    "VERSION",
    "templates/input_matrix_template.csv",
    "templates/labels_template.csv",
    "example_workflow_output/README.md",
    "example_workflow_output/calibration.json",
    "example_workflow_output/explanations.csv",
    "example_workflow_output/manifest.json",
    "example_workflow_output/qc.json",
    "example_workflow_output/scores.csv",
    "example_workflow_output/thresholds.csv",
    "example_workflow_output/workflow_report.md",
    "external-validation/cptac_gdc/CPTAC_EXTERNAL_VALIDATION.md",
    "external-validation/cptac_gdc/cptac_summary.csv",
    "external-validation/cptac_gdc/cptac_threshold_sweep.csv",
    "external-validation/gtex_xena/GTEX_NORMAL_VALIDATION.md",
    "external-validation/gtex_xena/gtex_summary.csv",
    "external-validation/gtex_xena/gtex_threshold_sweep.csv",
    "external-validation/tcga_toil_xena/TCGA_TOIL_PIPELINE_CHECK.md",
    "external-validation/tcga_toil_xena/tcga_toil_summary.csv",
    "external-validation/tcga_toil_xena/tcga_toil_threshold_sweep.csv",
]

BINARY_SUFFIXES = {".npy", ".npz", ".pkl", ".png", ".zip"}
TEXT_NAMES = {"LICENSE", "VERSION"}


def assert_inside_root(path):
    resolved = path.resolve()
    root = ROOT.resolve()
    if resolved != root and root not in resolved.parents:
        raise RuntimeError(f"Refusing to operate outside project root: {resolved}")
    return resolved


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_text_lf(path, text, encoding="utf-8"):
    path.write_bytes(text.encode(encoding))


def clean_release_dir():
    target = assert_inside_root(RELEASE_DIR)
    if target.exists():
        if target.name != "release-lite":
            raise RuntimeError(f"Unexpected release directory name: {target}")
        shutil.rmtree(target)
    target.mkdir(parents=True)


def is_text_release_file(rel):
    path = Path(rel)
    if path.name in TEXT_NAMES:
        return True
    return path.suffix.lower() not in BINARY_SUFFIXES


def sorted_release_paths():
    return sorted(
        (path for path in RELEASE_DIR.rglob("*") if path.is_file()),
        key=lambda path: path.relative_to(RELEASE_DIR).as_posix(),
    )


def write_release_file(src, dst, rel):
    if is_text_release_file(rel):
        data = src.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
        dst.write_bytes(data)
    else:
        shutil.copy2(src, dst)


def copy_release_files():
    missing = []
    for rel in RELEASE_FILES:
        src = ROOT / rel
        if not src.exists():
            missing.append(rel)
            continue
        dst = RELEASE_DIR / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        write_release_file(src, dst, rel)
    if missing:
        raise FileNotFoundError("Missing release source files:\n" + "\n".join(missing))


def manifest_file_records():
    records = []
    for path in sorted_release_paths():
        if path.name in {"SHA256SUMS.txt", "release_manifest.json"}:
            continue
        records.append({
            "path": path.relative_to(RELEASE_DIR).as_posix(),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        })
    return records


def write_manifest():
    records = manifest_file_records()
    metadata_path = RELEASE_DIR / "RELEASE_METADATA.json"
    release_metadata = {}
    if metadata_path.exists():
        release_metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    manifest = {
        "schema_version": "1.0",
        "bundle_name": "tcga-tumor-normal-release-lite",
        "version": release_metadata.get("version"),
        "release_date": release_metadata.get("release_date"),
        "intended_input": "GDC STAR-Counts-style log2(TPM+1), rows=samples, columns=Ensembl genes.",
        "builder": "build_release_lite.py",
        "validation_command": (
            "python validate_release_lite.py --release-dir release-lite "
            "--zip tcga-tumor-normal-release-lite.zip"
        ),
        "file_count_excluding_manifest_and_checksums": len(records),
        "forbidden_artifact_names": sorted(FORBIDDEN_NAMES),
        "files": records,
    }
    write_text_lf(
        RELEASE_DIR / "release_manifest.json",
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
    )
    return len(records)


def write_checksums():
    rows = []
    for path in sorted_release_paths():
        if path.name == "SHA256SUMS.txt":
            continue
        rel = path.relative_to(RELEASE_DIR).as_posix()
        rows.append(f"{sha256(path)}  {rel}")
    checksum_path = RELEASE_DIR / "SHA256SUMS.txt"
    write_text_lf(checksum_path, "\n".join(rows) + "\n", encoding="ascii")
    return len(rows)


def release_zip_datetime():
    metadata_path = ROOT / "RELEASE_METADATA.json"
    if metadata_path.exists():
        try:
            release_date = json.loads(metadata_path.read_text(encoding="utf-8")).get(
                "release_date"
            )
            if release_date:
                dt = datetime.strptime(release_date, "%Y-%m-%d")
                return (dt.year, dt.month, dt.day, 0, 0, 0)
        except (ValueError, json.JSONDecodeError):
            pass
    return (1980, 1, 1, 0, 0, 0)


def write_zip():
    if ZIP_PATH.exists():
        ZIP_PATH.unlink()
    fixed_datetime = release_zip_datetime()
    with zipfile.ZipFile(ZIP_PATH, "w", compression=zipfile.ZIP_DEFLATED,
                         compresslevel=9) as zf:
        for path in sorted_release_paths():
            rel = path.relative_to(RELEASE_DIR).as_posix()
            info = zipfile.ZipInfo(rel, date_time=fixed_datetime)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            zf.writestr(info, path.read_bytes(), compress_type=zipfile.ZIP_DEFLATED,
                        compresslevel=9)
    with zipfile.ZipFile(ZIP_PATH) as zf:
        bad = zf.testzip()
        if bad is not None:
            raise RuntimeError(f"Zip validation failed at {bad}")
        return len(zf.infolist())


def write_artifact_metadata(zip_entries):
    release_metadata_path = ROOT / "RELEASE_METADATA.json"
    release_metadata = {}
    if release_metadata_path.exists():
        release_metadata = json.loads(release_metadata_path.read_text(encoding="utf-8"))
    files = sorted_release_paths()
    artifact = {
        "schema_version": "1.0",
        "version": release_metadata.get("version"),
        "release_date": release_metadata.get("release_date"),
        "release_dir": RELEASE_DIR.relative_to(ROOT).as_posix(),
        "release_file_count": len(files),
        "release_total_bytes": sum(path.stat().st_size for path in files),
        "zip_path": ZIP_PATH.relative_to(ROOT).as_posix(),
        "zip_entries": zip_entries,
        "zip_bytes": ZIP_PATH.stat().st_size,
        "zip_sha256": sha256(ZIP_PATH),
        "validation_command": (
            "python validate_release_lite.py --release-dir release-lite "
            "--zip tcga-tumor-normal-release-lite.zip"
        ),
        "zip_acceptance_command": (
            "python validate_zip_bundle.py tcga-tumor-normal-release-lite.zip"
        ),
    }
    write_text_lf(ARTIFACTS_PATH, json.dumps(artifact, indent=2, sort_keys=True) + "\n")


def cleanup_transient_files():
    patterns = ["_smoke_*", "_acceptance_*", "*.pyc"]
    for pattern in patterns:
        for path in RELEASE_DIR.rglob(pattern):
            if path.is_file():
                path.unlink()
    for path in sorted(RELEASE_DIR.rglob("__pycache__"), reverse=True):
        if path.is_dir():
            shutil.rmtree(path)


def _run_checked(cmd, cwd, timeout_seconds):
    env = dict(os.environ)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    try:
        result = subprocess.run(cmd, cwd=cwd, text=True, env=env, timeout=timeout_seconds)
    except subprocess.TimeoutExpired:
        raise SystemExit(f"Command timed out after {timeout_seconds}s: {' '.join(cmd)}")
    if result.returncode != 0:
        raise SystemExit(result.returncode)


def run_smoke_test(timeout_seconds):
    for script in ["run_smoke_tests.py", "run_safety_tests.py"]:
        _run_checked([sys.executable, script], RELEASE_DIR, timeout_seconds)


def run_release_validation(timeout_seconds):
    _run_checked(
        [
            sys.executable,
            "validate_release_lite.py",
            "--release-dir",
            str(RELEASE_DIR),
            "--zip",
            str(ZIP_PATH),
            "--source-root",
            str(ROOT),
            "--artifacts",
            str(ARTIFACTS_PATH),
        ],
        ROOT,
        timeout_seconds,
    )


def build(run_smoke=False, timeout_seconds=300):
    clean_release_dir()
    copy_release_files()
    if run_smoke:
        run_smoke_test(timeout_seconds)
    cleanup_transient_files()
    manifest_count = write_manifest()
    checksum_count = write_checksums()
    zip_entries = write_zip()
    write_artifact_metadata(zip_entries)
    run_release_validation(timeout_seconds)
    print(f"[release] files in manifest: {manifest_count}")
    print(f"[release] files with checksums: {checksum_count}")
    print(f"[release] zip entries: {zip_entries}")
    print(f"[release] zip bytes: {ZIP_PATH.stat().st_size}")
    print(f"[release] zip sha256: {sha256(ZIP_PATH)}")
    print(f"[release] wrote {RELEASE_DIR}")
    print(f"[release] wrote {ZIP_PATH}")
    print(f"[release] wrote {ARTIFACTS_PATH}")


def main(argv=None):
    parser = argparse.ArgumentParser(description="Build release-lite and its zip archive.")
    parser.add_argument("--smoke", action="store_true",
                        help="run release-lite smoke and safety tests after copying files")
    parser.add_argument("--timeout-seconds", type=int, default=300,
                        help="per-step subprocess timeout (default: 300)")
    args = parser.parse_args(argv)
    build(run_smoke=args.smoke, timeout_seconds=args.timeout_seconds)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
