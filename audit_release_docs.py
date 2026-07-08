#!/usr/bin/env python3
"""Audit release documentation for stale scripts and bundle references."""
import argparse
import json
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent

VERSIONED_DOCS = [
    "DATA_DICTIONARY.md",
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

CODE_SPAN_RE = re.compile(r"`([^`\n]+)`")
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


def build_report():
    messages = []
    check_core_files(messages)
    check_version_consistency(messages)
    check_release_bundle_contents(messages)
    check_python_commands(messages)
    levels = {item["level"] for item in messages}
    status = "FAIL" if "ERROR" in levels else "WARN" if "WARNING" in levels else "PASS"
    return {
        "schema_version": "1.0",
        "status": status,
        "root": str(ROOT),
        "release_target_root": str(release_target_root()),
        "messages": messages,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Audit release documentation for stale local references."
    )
    parser.add_argument("-o", "--output", help="write JSON report")
    parser.add_argument("--strict", action="store_true",
                        help="return non-zero on warnings as well as errors")
    args = parser.parse_args(argv)

    report = build_report()
    for message in report["messages"]:
        stream = sys.stderr if message["level"] in {"ERROR", "WARNING"} else sys.stdout
        print(f"[docs-audit] {message['level']}: {message['message']}", file=stream)
    print(f"[docs-audit] status={report['status']}")
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
