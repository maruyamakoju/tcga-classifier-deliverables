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
import stat
import subprocess
import sys
import time
import unicodedata
import zipfile
from datetime import datetime
from pathlib import Path, PurePosixPath

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
    "tcga_rnaseq/validation.py",
    "tcga_rnaseq/py.typed",
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
RELEASE_SCHEMA_VERSION = "1.0"

# The archive format is deliberately boring. Stored members avoid zlib-version
# drift, and every header field that ZipInfo otherwise derives from the host OS
# is fixed explicitly. The resulting bytes are reproducible on Windows, Linux,
# and macOS for identical payload bytes.
CANONICAL_ZIP_COMPRESSION = zipfile.ZIP_STORED
CANONICAL_ZIP_CREATE_SYSTEM = 3  # Unix, regardless of the build host.
CANONICAL_ZIP_CREATE_VERSION = 20
CANONICAL_ZIP_EXTRACT_VERSION = 20
CANONICAL_ZIP_FILE_MODE = stat.S_IFREG | 0o644

# Hard safety ceilings used before CRC checks or extraction. These are much
# larger than the current ~1 MiB bundle while remaining small enough to reject
# accidental or malicious decompression bombs early.
DEFAULT_MAX_ZIP_ENTRIES = 1_024
DEFAULT_MAX_ZIP_ARCHIVE_BYTES = 100_000_000
DEFAULT_MAX_ZIP_MEMBER_BYTES = 10_000_000
DEFAULT_MAX_ZIP_TOTAL_BYTES = 100_000_000
DEFAULT_MAX_ZIP_COMPRESSION_RATIO = 200.0

RELEASE_VALIDATION_COMMAND = (
    "python validate_release_lite.py --release-dir release-lite "
    "--zip tcga-tumor-normal-release-lite.zip --source-root . "
    "--artifacts RELEASE_ARTIFACTS.json"
)
ZIP_ACCEPTANCE_COMMAND = (
    "python validate_zip_bundle.py tcga-tumor-normal-release-lite.zip"
)

_WINDOWS_RESERVED_COMPONENTS = {
    "AUX",
    "CON",
    "CONIN$",
    "CONOUT$",
    "NUL",
    "PRN",
    *(f"COM{number}" for number in range(1, 10)),
    *(f"LPT{number}" for number in range(1, 10)),
}


def normalize_release_path(rel):
    """Return a canonical POSIX release path or raise ValueError.

    PurePosixPath normalizes ``./`` and repeated slashes, so compare its output
    to the original string as well as rejecting traversal and absolute paths.
    This prevents two spellings from referring to the same extracted member.
    """
    if not isinstance(rel, str) or not rel.strip():
        raise ValueError(f"Invalid empty release path: {rel!r}")
    if any(ord(character) < 32 or ord(character) == 127 for character in rel):
        raise ValueError(f"Release path contains a control character: {rel!r}")
    if "\\" in rel:
        raise ValueError(f"Release paths must use forward slashes: {rel!r}")
    if any(character in '<>"|?*' for character in rel):
        raise ValueError(f"Release path contains a Windows-forbidden character: {rel!r}")
    if ":" in rel:
        raise ValueError(
            f"Release path contains a Windows drive/stream separator: {rel!r}"
        )
    if unicodedata.normalize("NFC", rel) != rel:
        raise ValueError(f"Release path is not Unicode NFC-normalized: {rel!r}")
    path = PurePosixPath(rel)
    normalized = path.as_posix()
    if path.is_absolute() or rel.startswith("/"):
        raise ValueError(f"Release path must be relative: {rel!r}")
    if any(part in {"", ".", ".."} for part in path.parts):
        raise ValueError(f"Release path contains unsafe component: {rel!r}")
    for component in path.parts:
        if component.endswith((" ", ".")):
            raise ValueError(
                f"Release path contains a Windows-ambiguous component: {rel!r}"
            )
        reserved_stem = component.split(".", 1)[0].upper()
        if reserved_stem in _WINDOWS_RESERVED_COMPONENTS:
            raise ValueError(
                f"Release path contains a Windows reserved name: {rel!r}"
            )
    if normalized != rel:
        raise ValueError(f"Release path is not canonical: {rel!r}")
    return normalized


def release_path_collision_key(rel):
    """Return the cross-platform identity key for a normalized path.

    Windows and the default macOS filesystem are case-insensitive. Treating
    case-fold-equivalent names as the same member keeps a bundle unambiguous
    when it moves between supported platforms.
    """
    return unicodedata.normalize("NFC", normalize_release_path(rel)).casefold()


def canonical_zip_datetime(release_date):
    """Convert an ISO release date to the one allowed archive timestamp."""
    if not isinstance(release_date, str):
        raise ValueError("release_date must be a YYYY-MM-DD string")
    try:
        parsed = datetime.strptime(release_date, "%Y-%m-%d")
    except ValueError as exc:
        raise ValueError(
            f"release_date must be a valid YYYY-MM-DD date: {release_date!r}"
        ) from exc
    if parsed.strftime("%Y-%m-%d") != release_date:
        raise ValueError(f"release_date is not canonical: {release_date!r}")
    if parsed.year < 1980:
        raise ValueError("release_date must be 1980-01-01 or later for ZIP")
    return (parsed.year, parsed.month, parsed.day, 0, 0, 0)


def canonical_zip_info(filename, release_datetime):
    """Create a fully specified, cross-platform deterministic file header."""
    normalized = normalize_release_path(filename)
    info = zipfile.ZipInfo(normalized, date_time=release_datetime)
    info.compress_type = CANONICAL_ZIP_COMPRESSION
    info.create_system = CANONICAL_ZIP_CREATE_SYSTEM
    info.create_version = CANONICAL_ZIP_CREATE_VERSION
    info.extract_version = CANONICAL_ZIP_EXTRACT_VERSION
    info.flag_bits = 0
    info.internal_attr = 0
    info.external_attr = CANONICAL_ZIP_FILE_MODE << 16
    info.extra = b""
    info.comment = b""
    return info


def zip_ratio(info):
    """Return an archive member's uncompressed/compressed size ratio."""
    if info.file_size == 0:
        return 1.0
    if info.compress_size == 0:
        return float("inf")
    return info.file_size / info.compress_size


def zip_safety_errors(
    infos,
    *,
    max_entries=DEFAULT_MAX_ZIP_ENTRIES,
    max_member_bytes=DEFAULT_MAX_ZIP_MEMBER_BYTES,
    max_total_bytes=DEFAULT_MAX_ZIP_TOTAL_BYTES,
    max_compression_ratio=DEFAULT_MAX_ZIP_COMPRESSION_RATIO,
):
    """Validate resource limits and path safety before reading ZIP contents."""
    errors = []
    if len(infos) > max_entries:
        errors.append(f"Zip has {len(infos)} entries; limit is {max_entries}")
    total = 0
    seen = set()
    seen_collision_keys = {}
    for info in infos:
        if getattr(info, "orig_filename", info.filename) != info.filename:
            errors.append(f"Zip member filename contains a NUL byte: {info.orig_filename!r}")
        try:
            rel = normalize_release_path(info.filename)
        except ValueError as exc:
            errors.append(f"Zip member {info.filename!r}: {exc}")
            continue
        if rel in seen:
            errors.append(f"Zip contains duplicate member path: {info.filename}")
        else:
            collision_key = release_path_collision_key(rel)
            previous = seen_collision_keys.get(collision_key)
            if previous is not None:
                errors.append(
                    "Zip contains case-insensitive member collision: "
                    f"{previous} and {info.filename}"
                )
            else:
                seen_collision_keys[collision_key] = info.filename
        seen.add(rel)
        if info.is_dir() or info.filename.endswith("/"):
            errors.append(f"Zip contains a directory entry: {info.filename}")
            continue
        if info.flag_bits & 0x1:
            errors.append(f"Zip contains encrypted member: {info.filename}")
        member_mode = (info.external_attr >> 16) & 0xFFFF
        member_type = stat.S_IFMT(member_mode)
        if member_type not in {0, stat.S_IFREG}:
            errors.append(f"Zip contains special-file member: {info.filename}")
        if info.file_size > max_member_bytes:
            errors.append(
                f"Zip member exceeds {max_member_bytes} uncompressed bytes: "
                f"{info.filename} ({info.file_size})"
            )
        total += info.file_size
        ratio = zip_ratio(info)
        if ratio > max_compression_ratio:
            errors.append(
                f"Zip member compression ratio exceeds {max_compression_ratio:g}: "
                f"{info.filename} ({ratio:.1f})"
            )
    if total > max_total_bytes:
        errors.append(f"Zip expands to {total} bytes; limit is {max_total_bytes}")
    return errors


def canonical_zip_errors(zf, infos, expected_datetime):
    """Validate deterministic archive metadata after basic safety checks."""
    errors = []
    if zf.comment != b"":
        errors.append("Zip archive comment must be empty")
    names = [info.filename for info in infos]
    if names != sorted(names):
        errors.append("Zip member order is not canonical (sorted POSIX paths required)")
    for info in infos:
        if info.date_time != expected_datetime:
            errors.append(f"Zip member timestamp is not canonical: {info.filename}")
        if info.compress_type != CANONICAL_ZIP_COMPRESSION:
            errors.append(f"Zip member compression is not canonical: {info.filename}")
        if info.create_system != CANONICAL_ZIP_CREATE_SYSTEM:
            errors.append(f"Zip member create_system is not canonical: {info.filename}")
        if info.create_version != CANONICAL_ZIP_CREATE_VERSION:
            errors.append(f"Zip member create_version is not canonical: {info.filename}")
        if info.extract_version != CANONICAL_ZIP_EXTRACT_VERSION:
            errors.append(f"Zip member extract_version is not canonical: {info.filename}")
        if info.external_attr != CANONICAL_ZIP_FILE_MODE << 16:
            errors.append(f"Zip member file mode is not canonical: {info.filename}")
        if info.internal_attr != 0:
            errors.append(f"Zip member internal attributes are not canonical: {info.filename}")
        if info.extra:
            errors.append(f"Zip member extra field must be empty: {info.filename}")
        if info.comment:
            errors.append(f"Zip member comment must be empty: {info.filename}")
        # Only the UTF-8 filename flag may be set, and all current payload names
        # are ASCII so canonical builds set no flags at all.
        if info.flag_bits != 0:
            errors.append(f"Zip member flags are not canonical: {info.filename}")
    return errors


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
