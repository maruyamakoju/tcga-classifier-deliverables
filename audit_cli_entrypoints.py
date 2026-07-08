#!/usr/bin/env python3
"""Audit release CLI entry points for usable --help output."""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent

HELP_COMMANDS = [
    "score_tumor_normal.py",
    "cohort_adapt_score.py",
    "inspect_expression_input.py",
    "run_tumor_normal_workflow.py",
    "calibrate_threshold.py",
    "explain_scores.py",
    "check_environment.py",
    "audit_lightweight_dependencies.py",
    "audit_release_docs.py",
    "audit_cli_entrypoints.py",
    "validate_output_contracts.py",
    "run_release_acceptance.py",
    "validate_release_lite.py",
    "validate_zip_bundle.py",
]


def add_message(messages, level, code, message, path=None):
    item = {"level": level, "code": code, "message": message}
    if path is not None:
        item["path"] = str(path)
    messages.append(item)


def release_target_root():
    if (ROOT / "release_manifest.json").exists() and (ROOT / "SHA256SUMS.txt").exists():
        return ROOT
    nested = ROOT / "release-lite"
    if (nested / "release_manifest.json").exists() and (nested / "SHA256SUMS.txt").exists():
        return nested
    return ROOT


def load_manifest_paths(target_root):
    manifest_path = target_root / "release_manifest.json"
    if not manifest_path.exists():
        return None
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None
    return [item["path"] for item in manifest.get("files", [])]


def release_python_files(target_root):
    manifest_paths = load_manifest_paths(target_root)
    if manifest_paths is not None:
        paths = [target_root / rel for rel in manifest_paths if rel.endswith(".py")]
    else:
        paths = sorted(target_root.glob("*.py"))
    return [path for path in paths if path.exists() and path.is_file()]


def check_shebangs(target_root, messages):
    for path in release_python_files(target_root):
        # Package modules (a directory holding an __init__.py) are importable
        # libraries, not executable CLIs, so they do not need a shebang.
        if (path.parent / "__init__.py").exists():
            continue
        rel = path.relative_to(target_root).as_posix()
        first_line = path.read_text(encoding="utf-8").splitlines()[0:1]
        if not first_line or first_line[0].strip() != "#!/usr/bin/env python3":
            add_message(messages, "ERROR", "python_shebang_missing",
                        f"{rel} should start with #!/usr/bin/env python3.", path)


def run_help(target_root, script, timeout_seconds, messages):
    path = target_root / script
    if not path.exists():
        add_message(messages, "ERROR", "help_script_missing",
                    f"Expected CLI script is missing: {script}", path)
        return None
    env = dict(os.environ)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    try:
        result = subprocess.run(
            [sys.executable, script, "--help"],
            cwd=target_root,
            text=True,
            capture_output=True,
            timeout=timeout_seconds,
            env=env,
        )
    except subprocess.TimeoutExpired:
        add_message(messages, "ERROR", "help_timeout",
                    f"{script} --help timed out after {timeout_seconds}s.", path)
        return None

    combined = (result.stdout or "") + "\n" + (result.stderr or "")
    if result.returncode != 0:
        add_message(messages, "ERROR", "help_nonzero_exit",
                    f"{script} --help returned {result.returncode}.", path)
    if "traceback" in combined.lower():
        add_message(messages, "ERROR", "help_traceback",
                    f"{script} --help printed a traceback.", path)
    if "usage:" not in combined.lower():
        add_message(messages, "ERROR", "help_usage_missing",
                    f"{script} --help did not include a usage line.", path)
    return {
        "script": script,
        "returncode": result.returncode,
        "stdout_bytes": len(result.stdout or ""),
        "stderr_bytes": len(result.stderr or ""),
    }


def build_report(timeout_seconds):
    messages = []
    target_root = release_target_root()
    check_shebangs(target_root, messages)
    help_results = []
    for script in HELP_COMMANDS:
        result = run_help(target_root, script, timeout_seconds, messages)
        if result is not None:
            help_results.append(result)
    levels = {item["level"] for item in messages}
    status = "FAIL" if "ERROR" in levels else "WARN" if "WARNING" in levels else "PASS"
    return {
        "schema_version": "1.0",
        "status": status,
        "root": str(ROOT),
        "target_root": str(target_root),
        "help_commands": HELP_COMMANDS,
        "help_results": help_results,
        "messages": messages,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Audit release CLI scripts for usable --help output."
    )
    parser.add_argument("-o", "--output", help="write JSON report")
    parser.add_argument("--timeout-seconds", type=float, default=20.0)
    parser.add_argument("--strict", action="store_true",
                        help="return non-zero on warnings as well as errors")
    args = parser.parse_args(argv)

    report = build_report(args.timeout_seconds)
    for message in report["messages"]:
        stream = sys.stderr if message["level"] in {"ERROR", "WARNING"} else sys.stdout
        print(f"[cli-audit] {message['level']}: {message['message']}", file=stream)
    print(f"[cli-audit] checked {len(report['help_results'])} help commands")
    print(f"[cli-audit] status={report['status']}")
    if args.output:
        out_path = Path(args.output)
        if not out_path.is_absolute():
            out_path = ROOT / out_path
        out_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                            encoding="utf-8")
        print(f"[cli-audit] wrote {out_path}")
    if report["status"] == "FAIL" or (args.strict and report["status"] == "WARN"):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
