"""Shared helpers for the release audit/validate/build script family.

Consolidates patterns that were previously copy-pasted (with small,
drifting variations) across audit_cli_entrypoints.py,
audit_lightweight_dependencies.py, audit_release_docs.py,
audit_publication_readiness.py, validate_output_contracts.py,
validate_release_lite.py, validate_zip_bundle.py, run_release_acceptance.py,
run_safety_tests.py, run_smoke_tests.py, and build_release_lite.py.

Several of the scripts above are shipped inside release-lite/ so this module
ships alongside them (see build_release_lite.RELEASE_FILES below) -- keep it
stdlib-only, same constraint as tcga_rnaseq.
"""
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

# The single "what does the lightweight release ship" list. build_release_lite.py
# copies these paths (relative to the full-deliverables root) into release-lite/
# and the zip; validate_release_lite.py's REQUIRED_FALLBACK and
# audit_release_docs.py's CORE_FILES are intentionally smaller, curated subsets
# and assert (at import time) that they stay subsets of this list.
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
    "release_tools/__init__.py",
    "release_tools/common.py",
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

# Never-ship training/full-artifact names inside the release-lite payload.
# Distinct in scope from audit_publication_readiness.FORBIDDEN_TRACKED_NAMES,
# which guards the full tracked repo (e.g. it forbids operon_watchdog.log,
# a log file that would never be a release-lite candidate; this list forbids
# training-provenance files like train_classifier.py that are legitimately
# tracked in the full repo but must never end up in the lightweight bundle).
FORBIDDEN_NAMES = {
    "deployable_pipeline.pkl",
    "feature_selection.pkl",
    "final_model_results.pkl",
    "gene_id_to_name.pkl",
    "groups_full.pkl",
    "model_lr.pkl",
    "model_rf.pkl",
    "model_xgb.pkl",
    "projects_full.pkl",
    "sample_metadata.pkl",
    "selected_files.csv",
    "train_classifier.py",
    "train_idx.npy",
    "test_idx.npy",
    "X_full_filtered.pkl",
    "y_full.pkl",
}

RELEASE_ZIP_NAME = "tcga-tumor-normal-release-lite.zip"
RELEASE_BUNDLE_NAME = "tcga-tumor-normal-release-lite"


def add_message(messages, level, code, message, path=None):
    """Append a structured {level, code, message[, path]} audit message."""
    item = {"level": level, "code": code, "message": message}
    if path is not None:
        item["path"] = str(path)
    messages.append(item)


def status_from_levels(messages):
    """FAIL if any ERROR message, WARN if any WARNING, else PASS."""
    levels = {item["level"] for item in messages}
    if "ERROR" in levels:
        return "FAIL"
    if "WARNING" in levels:
        return "WARN"
    return "PASS"


def exit_code_for_status(status, strict=False):
    """1 on FAIL, or on WARN when --strict was requested; else 0."""
    if status == "FAIL" or (strict and status == "WARN"):
        return 1
    return 0


def sha256_file(path):
    """Streaming SHA-256 of a file (safe for large release artifacts)."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json_report(path, report, root=None, prefix="report"):
    """Write an indented, sorted-keys JSON report with a trailing newline.

    Always creates parent directories (some call sites previously omitted
    this and would crash on a non-existent -o subdirectory). If root is
    given and path is relative, it is resolved against root.
    """
    out_path = Path(path)
    if root is not None and not out_path.is_absolute():
        out_path = Path(root) / out_path
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")
    print(f"[{prefix}] wrote {out_path}")
    return out_path


def release_target_root(root):
    """Prefer root/release-lite if it looks like a built release, else root."""
    root = Path(root)
    if (root / "release_manifest.json").exists() and (root / "SHA256SUMS.txt").exists():
        return root
    nested = root / "release-lite"
    if (nested / "release_manifest.json").exists() and (nested / "SHA256SUMS.txt").exists():
        return nested
    return root


def load_manifest_paths(target_root):
    """Relative paths listed in target_root/release_manifest.json, or None."""
    manifest_path = Path(target_root) / "release_manifest.json"
    if not manifest_path.exists():
        return None
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None
    return [item["path"] for item in manifest.get("files", [])]


def release_python_files(target_root):
    """The .py files the release manifest lists, or a bare glob fallback."""
    target_root = Path(target_root)
    manifest_paths = load_manifest_paths(target_root)
    if manifest_paths is not None:
        paths = [target_root / rel for rel in manifest_paths if rel.endswith(".py")]
    else:
        paths = sorted(target_root.glob("*.py"))
    return [path for path in paths if path.exists() and path.is_file()]


def subprocess_output_text(value):
    """Decode subprocess output that may be None, bytes, or already str."""
    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return str(value)


def append_timeout_message(stderr, timeout_seconds):
    message = f"Timed out after {timeout_seconds}s"
    if stderr:
        return stderr.rstrip("\n") + "\n" + message
    return message


def run_subprocess_step(label, cmd, cwd, timeout_seconds=300, required=True, prefix="step"):
    """Run a subprocess step, capturing output and timing, as a report dict.

    Unifies the run_step()/subprocess_output_text()/append_timeout_message()
    trio previously duplicated (with small drifting variations, including an
    inconsistent Windows text-decoding codec) across run_release_acceptance.py
    and validate_zip_bundle.py.
    """
    print(f"[{prefix}] {label}: {' '.join(str(x) for x in cmd)}")
    started = time.perf_counter()
    env = dict(os.environ)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    try:
        result = subprocess.run(
            cmd, cwd=cwd, text=True, capture_output=True,
            encoding="utf-8", errors="replace",
            timeout=timeout_seconds, env=env,
        )
    except subprocess.TimeoutExpired as exc:
        duration = time.perf_counter() - started
        stdout = subprocess_output_text(exc.stdout)
        stderr = subprocess_output_text(exc.stderr)
        print(f"[{prefix}] {label}: FAIL timeout after {duration:.1f}s", file=sys.stderr)
        return {
            "label": label,
            "command": [str(x) for x in cmd],
            "cwd": str(cwd),
            "required": required,
            "returncode": 124,
            "status": "FAIL",
            "duration_seconds": round(duration, 3),
            "stdout": stdout,
            "stderr": append_timeout_message(stderr, timeout_seconds),
        }
    duration = time.perf_counter() - started
    if result.stdout:
        print(result.stdout.rstrip())
    if result.stderr:
        print(result.stderr.rstrip(), file=sys.stderr)
    status = "PASS" if result.returncode == 0 else "FAIL"
    if not required and result.returncode != 0:
        status = "WARN"
    print(f"[{prefix}] {label}: {status} ({duration:.1f}s)")
    return {
        "label": label,
        "command": [str(x) for x in cmd],
        "cwd": str(cwd),
        "required": required,
        "returncode": result.returncode,
        "status": status,
        "duration_seconds": round(duration, 3),
        "stdout": result.stdout,
        "stderr": result.stderr,
    }
