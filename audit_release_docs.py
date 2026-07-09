#!/usr/bin/env python3
"""Audit release documentation for stale scripts and bundle references."""
import argparse
import json
import re
import sys
from urllib.parse import unquote
from pathlib import Path


ROOT = Path(__file__).resolve().parent

VERSIONED_DOCS = [
    "DATA_DICTIONARY.md",
    "INDEX.md",
    "README.md",
    "USER_GUIDE.md",
    "RELEASE_BUNDLE.md",
    "RELEASE_NOTES.md",
    "EXECUTIVE_SUMMARY.md",
    "REPRODUCIBILITY.md",
    "REPORT.md",
    "TROUBLESHOOTING.md",
]

COMMAND_DOCS = VERSIONED_DOCS + ["MODEL_CARD.md"]

OPTIONAL_FULL_REPO_DOCS = [
    "PUBLICATION_CHECKLIST.md",
]

CORE_FILES = [
    "audit_cli_entrypoints.py",
    "audit_lightweight_dependencies.py",
    "audit_release_docs.py",
    "calibrate_threshold.py",
    "check_environment.py",
    "deployable_lr_weights.npz",
    "DATA_DICTIONARY.md",
    "example_input.csv",
    "example_labels.csv",
    "example_output.csv",
    "LICENSE",
    "NOTICE.md",
    "CITATION.cff",
    ".zenodo.json",
    "codemeta.json",
    "explain_scores.py",
    "inspect_expression_input.py",
    "INDEX.md",
    "model_gene_metadata.csv",
    "model_qc_reference.json",
    "README.md",
    "RELEASE_BUNDLE.md",
    "RELEASE_METADATA.json",
    "requirements-light.txt",
    "run_release_acceptance.py",
    "run_safety_tests.py",
    "run_smoke_tests.py",
    "run_tumor_normal_workflow.py",
    "score_tumor_normal.py",
    "TROUBLESHOOTING.md",
    "USER_GUIDE.md",
    "validate_output_contracts.py",
    "validate_release_lite.py",
    "validate_zip_bundle.py",
    "VERSION",
]

FULL_DELIVERABLE_COMMANDS = {"build_release_lite.py"}

PATH_LIKE_SUFFIXES = {
    ".cff",
    ".csv",
    ".html",
    ".json",
    ".md",
    ".npz",
    ".pkl",
    ".png",
    ".py",
    ".txt",
    ".yaml",
    ".yml",
    ".zip",
}

FULL_DELIVERABLE_ONLY_REFS = {
    "CANCER_TYPE_CLASSIFIER.md",
    "audit_github_repository.py",
    "audit_publication_readiness.py",
    "build_release_lite.py",
    "CROSS_PLATFORM_ADAPTATION.md",
    "deployable_pipeline.pkl",
    "environment.yml",
    "feature_importance.png",
    "export_lr_weights.py",
    "export_model_gene_metadata.py",
    "export_qc_reference.py",
    "feature_selection.pkl",
    "final_model_results.pkl",
    "model_performance.png",
    "loco_generalization.png",
    "loco_per_cancer_metrics.csv",
    "loco_pooled_summary.csv",
    "loco_predictions.pkl",
    "loco_report.html",
    "LOCO_REPORT.md",
    "loco_vs_within_comparison.csv",
    "model_lr.pkl",
    "model_rf.pkl",
    "model_xgb.pkl",
    "predict_cancer_type.py",
    "requirements.txt",
    "selected_files.csv",
    "tcga-tumor-normal-release-lite.zip",
    "train_classifier.py",
    "X_full_filtered.pkl",
    "y_full.pkl",
    "external-validation/cptac_gdc/cptac_gdc_manifest.csv",
    "external-validation/cptac_gdc/cptac_predictions.csv",
    "external-validation/cptac_gdc/sampled_manifest.csv",
    "external-validation/gtex_xena/TcgaTargetGTEX_phenotype.csv",
    "external-validation/gtex_xena/gtex_per_site_summary.csv",
    "external-validation/gtex_xena/gtex_predictions.csv",
    "external-validation/gtex_xena/sampled_gtex_manifest.csv",
    "external-validation/tcga_toil_xena/TcgaTargetGTEX_phenotype.csv",
    "external-validation/tcga_toil_xena/sampled_tcga_toil_manifest.csv",
    "external-validation/tcga_toil_xena/tcga_toil_predictions.csv",
    "external-validation/validate_cptac_gdc.py",
    "external-validation/validate_gtex_xena.py",
    "external-validation/validate_tcga_toil_xena.py",
}

FULL_DELIVERABLE_ONLY_PREFIXES = (
    "cancer-type-classifier/",
    "cross-cancer-holdout/",
    "cross-platform-adaptation/",
    "from-workbench-loco/",
    "tests/",
)

RELEASE_SIDECAR_REFS = {
    "RELEASE_ARTIFACTS.json",
    "release_manifest.json",
    "SHA256SUMS.txt",
}

GENERATED_OUTPUT_REFS = {
    "calibration.json",
    "explanations.csv",
    "labels.csv",
    "manifest.json",
    "qc.json",
    "scores.csv",
    "thresholds.csv",
    "workflow_report.md",
}

HISTORICAL_RELEASE_NOTE_REFS = {
    "release-lite/release-lite",
}

CODE_SPAN_RE = re.compile(r"`([^`\n]+)`")
MARKDOWN_LINK_RE = re.compile(r"!?\[[^\]\n]*\]\(([^)\n]+)\)")
PYTHON_SCRIPT_COMMAND_RE = re.compile(
    r"(?m)^\s*python\s+([A-Za-z0-9_./\\-]+\.py)\b"
)


def add_message(messages, level, code, message, path=None):
    item = {"level": level, "code": code, "message": message}
    if path is not None:
        item["path"] = str(path)
    messages.append(item)


def read_text(path, messages):
    if not path.exists():
        add_message(messages, "ERROR", "doc_missing", f"Document is missing: {path.name}", path)
        return None
    try:
        return path.read_text(encoding="utf-8")
    except UnicodeDecodeError as exc:
        add_message(messages, "ERROR", "doc_not_utf8", f"Document is not UTF-8: {exc}", path)
        return None


def release_target_root():
    if (ROOT / "release_manifest.json").exists() and (ROOT / "SHA256SUMS.txt").exists():
        return ROOT
    nested = ROOT / "release-lite"
    if (nested / "release_manifest.json").exists() and (nested / "SHA256SUMS.txt").exists():
        return nested
    return ROOT


def check_core_files(messages):
    for rel in CORE_FILES:
        path = ROOT / rel
        if not path.exists():
            add_message(messages, "ERROR", "core_file_missing",
                        f"Core release file is missing: {rel}", path)


def check_version_consistency(messages):
    version_path = ROOT / "VERSION"
    metadata_path = ROOT / "RELEASE_METADATA.json"
    if not version_path.exists() or not metadata_path.exists():
        return
    version = version_path.read_text(encoding="utf-8").strip()
    try:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        add_message(messages, "ERROR", "release_metadata_invalid_json",
                    f"RELEASE_METADATA.json is invalid JSON: {exc}", metadata_path)
        return
    if metadata.get("version") != version:
        add_message(messages, "ERROR", "release_version_mismatch",
                    "VERSION and RELEASE_METADATA.json disagree.", metadata_path)
    for rel in VERSIONED_DOCS:
        path = ROOT / rel
        text = read_text(path, messages)
        if text is not None and version not in text:
            add_message(messages, "ERROR", "doc_version_missing",
                        f"{rel} does not mention release version {version}.", path)
    for rel in OPTIONAL_FULL_REPO_DOCS:
        path = ROOT / rel
        if not path.exists():
            continue
        text = read_text(path, messages)
        if text is not None and version not in text:
            add_message(messages, "ERROR", "doc_version_missing",
                        f"{rel} does not mention release version {version}.", path)


def check_release_bundle_contents(messages):
    target_root = release_target_root()
    bundle_path = target_root / "RELEASE_BUNDLE.md"
    text = read_text(bundle_path, messages)
    if text is None:
        return
    lines = text.splitlines()
    in_contents = False
    found_refs = 0
    for line in lines:
        if line.strip() == "## Contents":
            in_contents = True
            continue
        if in_contents and line.startswith("## "):
            break
        if not in_contents or not line.lstrip().startswith("- "):
            continue
        for ref in CODE_SPAN_RE.findall(line):
            if " " in ref or ref.startswith("-"):
                continue
            found_refs += 1
            normalized = ref.replace("\\", "/").rstrip("/")
            path = target_root / normalized
            if ref.endswith("/"):
                if not path.is_dir():
                    add_message(messages, "ERROR", "bundle_directory_missing",
                                f"RELEASE_BUNDLE.md references missing directory: {ref}", path)
            elif not path.exists():
                add_message(messages, "ERROR", "bundle_file_missing",
                            f"RELEASE_BUNDLE.md references missing file: {ref}", path)
    if found_refs == 0:
        add_message(messages, "ERROR", "bundle_contents_empty",
                    "No code-spanned file references found in RELEASE_BUNDLE.md Contents.",
                    bundle_path)


def check_python_commands(messages):
    for rel in COMMAND_DOCS:
        path = ROOT / rel
        text = read_text(path, messages)
        if text is None:
            continue
        for match in PYTHON_SCRIPT_COMMAND_RE.finditer(text):
            script = match.group(1).replace("\\", "/")
            script_path = ROOT / script
            if script_path.exists():
                continue
            if script in FULL_DELIVERABLE_COMMANDS:
                add_message(messages, "INFO", "full_deliverable_command",
                            f"{rel} references full-deliverables-only command: {script}", path)
                continue
            add_message(messages, "ERROR", "command_script_missing",
                        f"{rel} references missing Python script: {script}", script_path)


def normalize_inline_path(raw_target):
    target = raw_target.strip()
    if not target or target.startswith(("-", "$")):
        return None
    if "<" in target or ">" in target:
        return None
    if any(ch.isspace() for ch in target):
        return None
    target = target.replace("\\", "/")
    target = target.split("#", 1)[0].split("?", 1)[0]
    target = target.rstrip(".,;:")
    if not target or target in {".", ".."}:
        return None
    if target.startswith(("/", "\\")) or re.match(r"^[A-Za-z]:", target):
        return target
    if "/" in target:
        return target
    suffix = Path(target).suffix.lower()
    if suffix in PATH_LIKE_SUFFIXES:
        return target
    if target in CORE_FILES or target in RELEASE_SIDECAR_REFS:
        return target
    return None


def is_intentionally_external_to_lite(target):
    normalized = target.replace("\\", "/").rstrip("/")
    if (
        normalized in FULL_DELIVERABLE_ONLY_REFS
        or normalized in RELEASE_SIDECAR_REFS
        or normalized in HISTORICAL_RELEASE_NOTE_REFS
    ):
        return True
    return any(
        normalized == prefix.rstrip("/") or normalized.startswith(prefix)
        for prefix in FULL_DELIVERABLE_ONLY_PREFIXES
    )


def is_release_root():
    return (ROOT / "release_manifest.json").exists() and (ROOT / "SHA256SUMS.txt").exists()


def candidate_exists(base_path, target):
    if "*" in target or "?" in target:
        return any(base_path.parent.glob(target))
    target_path = (base_path.parent / target).resolve()
    return target_path.is_dir() if target.endswith("/") else target_path.exists()


def release_lite_prefix_exists(base_path, target):
    if not is_release_root():
        return False
    normalized = target.replace("\\", "/")
    stripped = normalized.rstrip("/")
    if stripped == "release-lite":
        mapped_target = "."
    elif normalized.startswith("release-lite/"):
        mapped_target = normalized[len("release-lite/"):]
    else:
        return False
    if not mapped_target:
        mapped_target = "."
    return candidate_exists(base_path, mapped_target)


def check_code_spanned_paths(messages):
    docs = list(dict.fromkeys(COMMAND_DOCS + OPTIONAL_FULL_REPO_DOCS))
    for rel in docs:
        path = ROOT / rel
        if not path.exists():
            continue
        text = read_text(path, messages)
        if text is None:
            continue
        seen = set()
        for match in CODE_SPAN_RE.finditer(text):
            target = normalize_inline_path(match.group(1))
            if target is None or target in seen:
                continue
            seen.add(target)
            if target.startswith(("/", "\\")) or re.match(r"^[A-Za-z]:", target):
                add_message(messages, "ERROR", "code_span_absolute_local_path",
                            f"{rel} uses an absolute local code-spanned path: {target}",
                            path)
                continue
            target_path = (path.parent / target).resolve()
            try:
                target_path.relative_to(ROOT.resolve())
            except ValueError:
                add_message(messages, "ERROR", "code_span_path_outside_root",
                            f"{rel} references a code-spanned path outside the root: {target}",
                            path)
                continue
            exists = candidate_exists(path, target) or release_lite_prefix_exists(path, target)
            if exists:
                continue
            if target.replace("\\", "/").rstrip("/") in GENERATED_OUTPUT_REFS:
                continue
            if is_intentionally_external_to_lite(target):
                add_message(messages, "INFO", "full_deliverable_reference",
                            f"{rel} references a full-deliverable or sidecar path: {target}",
                            path)
                continue
            add_message(messages, "ERROR", "code_span_path_missing",
                        f"{rel} references missing code-spanned local path: {target}",
                        target_path)


def is_external_link(target):
    lower = target.lower()
    return (
        lower.startswith(("http://", "https://", "mailto:"))
        or lower.startswith("#")
    )


def normalize_markdown_target(raw_target):
    target = raw_target.strip()
    if not target or is_external_link(target):
        return None
    if target.startswith("<") and target.endswith(">"):
        target = target[1:-1].strip()
    target = target.split("#", 1)[0].split("?", 1)[0]
    if not target:
        return None
    return unquote(target)


def check_markdown_links(messages):
    docs = list(dict.fromkeys(COMMAND_DOCS + OPTIONAL_FULL_REPO_DOCS))
    for rel in docs:
        path = ROOT / rel
        if not path.exists():
            continue
        text = read_text(path, messages)
        if text is None:
            continue
        for match in MARKDOWN_LINK_RE.finditer(text):
            target = normalize_markdown_target(match.group(1))
            if target is None:
                continue
            if target.startswith(("/", "\\")) or re.match(r"^[A-Za-z]:", target):
                add_message(messages, "ERROR", "markdown_absolute_local_link",
                            f"{rel} uses an absolute local Markdown link: {target}", path)
                continue
            link_path = (path.parent / target).resolve()
            try:
                link_path.relative_to(ROOT.resolve())
            except ValueError:
                add_message(messages, "ERROR", "markdown_link_outside_root",
                            f"{rel} links outside the release root: {target}", path)
                continue
            if not link_path.exists():
                add_message(messages, "ERROR", "markdown_link_missing",
                            f"{rel} links to missing local path: {target}", link_path)


def build_report():
    messages = []
    check_core_files(messages)
    check_version_consistency(messages)
    check_release_bundle_contents(messages)
    check_python_commands(messages)
    check_code_spanned_paths(messages)
    check_markdown_links(messages)
    levels = {item["level"] for item in messages}
    status = "FAIL" if "ERROR" in levels else "WARN" if "WARNING" in levels else "PASS"
    return {
        "schema_version": "1.0",
        "status": status,
        "root": str(ROOT),
        "release_target_root": str(release_target_root()),
        "messages": messages,
    }


def print_report(report, show_info=False):
    hidden_info = 0
    for message in report["messages"]:
        if message["level"] == "INFO" and not show_info:
            hidden_info += 1
            continue
        stream = sys.stderr if message["level"] in {"ERROR", "WARNING"} else sys.stdout
        print(f"[docs-audit] {message['level']}: {message['message']}", file=stream)
    if hidden_info:
        print(f"[docs-audit] info_messages={hidden_info} hidden; use --show-info to display")
    print(f"[docs-audit] status={report['status']}")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Audit release documentation for stale local references."
    )
    parser.add_argument("-o", "--output", help="write JSON report")
    parser.add_argument("--show-info", action="store_true",
                        help="print informational messages; JSON output always includes them")
    parser.add_argument("--strict", action="store_true",
                        help="return non-zero on warnings as well as errors")
    args = parser.parse_args(argv)

    report = build_report()
    print_report(report, show_info=args.show_info)
    if args.output:
        out_path = Path(args.output)
        if not out_path.is_absolute():
            out_path = ROOT / out_path
        out_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                            encoding="utf-8")
        print(f"[docs-audit] wrote {out_path}")
    if report["status"] == "FAIL" or (args.strict and report["status"] == "WARN"):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
