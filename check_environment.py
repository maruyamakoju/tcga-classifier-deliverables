#!/usr/bin/env python3
"""Check runtime environment and required release files."""
import argparse
import importlib
import importlib.metadata
import subprocess
import sys
from pathlib import Path

from tcga_rnaseq import write_json

ROOT = Path(__file__).resolve().parent

REQUIRED_FILES = [
    "audit_release_docs.py",
    "audit_lightweight_dependencies.py",
    "audit_cli_entrypoints.py",
    "check_environment.py",
    "validate_output_contracts.py",
    "validate_zip_bundle.py",
    "deployable_lr_weights.npz",
    "score_tumor_normal.py",
    "inspect_expression_input.py",
    "run_tumor_normal_workflow.py",
    "run_release_acceptance.py",
    "run_smoke_tests.py",
    "run_safety_tests.py",
    "model_qc_reference.json",
    "model_gene_metadata.csv",
    "example_input.csv",
    "example_output.csv",
    "RELEASE_BUNDLE.md",
    "DATA_DICTIONARY.md",
    "TROUBLESHOOTING.md",
    "USER_GUIDE.md",
    "tcga_rnaseq/__init__.py",
]

PACKAGE_RULES = {
    "numpy": {"required": True, "min": (1, 26), "max_major": 3},
    "pandas": {"required": True, "min": (2, 3), "max_major": 4},
    "pyarrow": {"required": False, "min": (16,), "max_major": None},
    "tcga_rnaseq": {"required": True, "min": None, "max_major": None},
}


def parse_version(version):
    parts = []
    for token in version.replace("-", ".").split("."):
        digits = "".join(ch for ch in token if ch.isdigit())
        if digits == "":
            break
        parts.append(int(digits))
    return tuple(parts) if parts else (0,)


def add_message(messages, level, code, text):
    messages.append({"level": level, "code": code, "message": text})


def check_python(messages):
    version = tuple(sys.version_info[:3])
    if version < (3, 9):
        add_message(messages, "ERROR", "python_too_old",
                    f"Python {version[0]}.{version[1]} detected; use Python >=3.9.")
    return {
        "executable": sys.executable,
        "version": ".".join(str(x) for x in version),
    }


def check_package(name, rule, messages):
    try:
        module = importlib.import_module(name)
    except Exception as exc:
        level = "ERROR" if rule["required"] else "WARNING"
        add_message(messages, level, f"{name}_missing",
                    f"{name} is not importable: {exc}")
        return {"installed": False, "version": None}

    try:
        version = importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        version = getattr(module, "__version__", "unknown")

    parsed = parse_version(version)
    min_version = rule.get("min")
    max_major = rule.get("max_major")
    if min_version and parsed < min_version:
        level = "ERROR" if rule["required"] else "WARNING"
        add_message(messages, level, f"{name}_too_old",
                    f"{name} {version} detected; expected >= {'.'.join(map(str, min_version))}.")
    if max_major is not None and parsed and parsed[0] >= max_major:
        level = "ERROR" if rule["required"] else "WARNING"
        add_message(messages, level, f"{name}_too_new",
                    f"{name} {version} detected; expected major version < {max_major}.")

    return {"installed": True, "version": version}


def check_files(messages):
    files = {}
    for rel in REQUIRED_FILES:
        path = ROOT / rel
        exists = path.exists()
        files[rel] = {"exists": exists, "bytes": path.stat().st_size if exists else None}
        if not exists:
            add_message(messages, "ERROR", "required_file_missing",
                        f"Required file is missing: {rel}")
    return files


def run_self_test(messages):
    result = subprocess.run(
        [sys.executable, "score_tumor_normal.py", "--self-test"],
        cwd=ROOT,
        text=True,
        capture_output=True,
    )
    if result.returncode != 0:
        add_message(messages, "ERROR", "self_test_failed",
                    "score_tumor_normal.py --self-test failed.")
    return {
        "returncode": result.returncode,
        "stdout": result.stdout.strip(),
        "stderr": result.stderr.strip(),
    }


def build_report(run_test=False):
    messages = []
    report = {
        "python": check_python(messages),
        "packages": {},
        "required_files": check_files(messages),
        "self_test": None,
    }
    for name, rule in PACKAGE_RULES.items():
        report["packages"][name] = check_package(name, rule, messages)
    if run_test:
        report["self_test"] = run_self_test(messages)

    levels = {message["level"] for message in messages}
    report["status"] = "FAIL" if "ERROR" in levels else "WARN" if "WARNING" in levels else "PASS"
    report["messages"] = messages
    return report


def print_summary(report, output_path):
    print(f"[env] status={report['status']}", file=sys.stderr)
    print(f"[env] python={report['python']['version']} ({report['python']['executable']})",
          file=sys.stderr)
    for name, info in report["packages"].items():
        if info["installed"]:
            print(f"[env] {name}={info['version']}", file=sys.stderr)
        else:
            print(f"[env] {name}=missing", file=sys.stderr)
    for message in report["messages"]:
        print(f"[env] {message['level']}: {message['message']}", file=sys.stderr)
    if output_path:
        print(f"[env] wrote {output_path}", file=sys.stderr)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Check runtime environment and release files.")
    parser.add_argument("-o", "--output", help="write JSON report")
    parser.add_argument("--self-test", action="store_true",
                        help="also run score_tumor_normal.py --self-test")
    parser.add_argument("--strict", action="store_true",
                        help="return non-zero on warnings as well as errors")
    args = parser.parse_args(argv)

    report = build_report(run_test=args.self_test)
    if args.output:
        write_json(report, args.output, sort_keys=True)
    print_summary(report, args.output)

    if report["status"] == "FAIL" or (args.strict and report["status"] == "WARN"):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
