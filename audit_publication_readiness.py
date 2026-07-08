#!/usr/bin/env python3
"""Audit repository state before making a release repository public."""
import argparse
import hashlib
import json
import re
import subprocess
import sys
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parent
ZIP_NAME = "tcga-tumor-normal-release-lite.zip"

MAX_TRACKED_BYTES = 25 * 1024 * 1024
MAX_HISTORY_BLOB_BYTES = 25 * 1024 * 1024

BINARY_SUFFIXES = {".npy", ".npz", ".pkl", ".png", ".zip"}
TEXT_NAMES = {"LICENSE", "VERSION"}

REQUIRED_FILES = {
    ".github/workflows/ci.yml",
    ".gitattributes",
    ".gitignore",
    ".zenodo.json",
    "CITATION.cff",
    "CONTRIBUTING.md",
    "GITHUB_RELEASE_v1.1.2.md",
    "LICENSE",
    "NOTICE.md",
    "PUBLICATION_CHECKLIST.md",
    "README.md",
    "RELEASE_ARTIFACTS.json",
    "SECURITY.md",
    "VERSION",
    "codemeta.json",
    ZIP_NAME,
}

REQUIRED_GITIGNORE_MARKERS = {
    "*.pkl",
    "X_full_filtered.pkl",
    "deployable_pipeline.pkl",
    "external-validation/**/gene_cache/",
    "*.log",
}

FORBIDDEN_TRACKED_NAMES = {
    "deployable_pipeline.pkl",
    "feature_selection.pkl",
    "final_model_results.pkl",
    "gene_id_to_name.pkl",
    "groups_full.pkl",
    "model_lr.pkl",
    "model_rf.pkl",
    "model_xgb.pkl",
    "operon_watchdog.log",
    "projects_full.pkl",
    "sample_metadata.pkl",
    "X_full_filtered.pkl",
    "y_full.pkl",
}

SECRET_PATTERNS = [
    ("private_key", re.compile(rb"-----BEGIN (?:RSA |DSA |EC |OPENSSH )?PRIVATE KEY-----")),
    ("github_token", re.compile(rb"\b(?:ghp|gho|ghu|ghs|ghr)_[A-Za-z0-9_]{20,}\b")),
    ("github_pat", re.compile(rb"\bgithub_pat_[A-Za-z0-9_]{20,}\b")),
    ("openai_key", re.compile(rb"\bsk-[A-Za-z0-9]{20,}\b")),
    ("aws_access_key", re.compile(rb"\bAKIA[0-9A-Z]{16}\b")),
    (
        "assigned_secret",
        re.compile(
            rb"(?i)\b(api[_-]?key|access[_-]?token|refresh[_-]?token|auth[_-]?token|"
            rb"client[_-]?secret|private[_-]?key|password|passwd)\b\s*[:=]\s*['\"]?[^'\"\s]+"
        ),
    ),
]


def add_message(messages, level, code, message, path=None):
    item = {"level": level, "code": code, "message": message}
    if path is not None:
        item["path"] = str(path)
    messages.append(item)


def run_git(args):
    result = subprocess.run(
        ["git"] + list(args),
        cwd=ROOT,
        text=True,
        capture_output=True,
    )
    if result.returncode != 0:
        raise RuntimeError(result.stderr.strip() or "git command failed")
    return result.stdout


def tracked_files():
    data = run_git(["ls-files", "-z"])
    return [item for item in data.split("\0") if item]


def is_text_path(rel):
    path = Path(rel)
    if path.name in TEXT_NAMES:
        return True
    return path.suffix.lower() not in BINARY_SUFFIXES


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def check_required_files(files, messages):
    file_set = set(files)
    for rel in sorted(REQUIRED_FILES):
        if rel not in file_set or not (ROOT / rel).exists():
            add_message(messages, "ERROR", "required_file_missing",
                        f"Required publication file is missing: {rel}", ROOT / rel)


def check_gitignore(messages):
    path = ROOT / ".gitignore"
    if not path.exists():
        add_message(messages, "ERROR", "gitignore_missing", ".gitignore is missing.", path)
        return
    text = path.read_text(encoding="utf-8")
    for marker in sorted(REQUIRED_GITIGNORE_MARKERS):
        if marker not in text:
            add_message(messages, "ERROR", "gitignore_marker_missing",
                        f".gitignore does not include expected marker: {marker}", path)


def check_tracked_file_sizes(files, messages, max_bytes):
    for rel in files:
        path = ROOT / rel
        if not path.is_file():
            continue
        name = Path(rel).name
        if name in FORBIDDEN_TRACKED_NAMES:
            add_message(messages, "ERROR", "forbidden_tracked_artifact",
                        f"Forbidden large/private artifact is tracked: {rel}", path)
        size = path.stat().st_size
        if size > max_bytes:
            add_message(messages, "ERROR", "tracked_file_too_large",
                        f"Tracked file exceeds {max_bytes} bytes: {rel} ({size})", path)


def check_line_endings(files, messages):
    for rel in files:
        if not is_text_path(rel):
            continue
        path = ROOT / rel
        if path.is_file() and b"\r" in path.read_bytes():
            add_message(messages, "ERROR", "cr_newline",
                        f"Tracked text file contains CR newline bytes: {rel}", path)


def check_secret_patterns(files, messages):
    for rel in files:
        if not is_text_path(rel):
            continue
        path = ROOT / rel
        if not path.is_file():
            continue
        data = path.read_bytes()
        for code, pattern in SECRET_PATTERNS:
            if pattern.search(data):
                add_message(messages, "ERROR", "secret_pattern",
                            f"Potential secret pattern '{code}' found in {rel}", path)


def check_history_blob_sizes(messages, max_bytes):
    rev_list = run_git(["rev-list", "--objects", "--all"]).splitlines()
    for line in rev_list:
        parts = line.split(" ", 1)
        sha = parts[0]
        rel = parts[1] if len(parts) > 1 else ""
        object_type = run_git(["cat-file", "-t", sha]).strip()
        if object_type != "blob":
            continue
        size = int(run_git(["cat-file", "-s", sha]).strip())
        if size > max_bytes:
            add_message(messages, "ERROR", "history_blob_too_large",
                        f"Git history blob exceeds {max_bytes} bytes: {rel or sha} ({size})")


def check_release_artifacts(messages):
    artifacts_path = ROOT / "RELEASE_ARTIFACTS.json"
    zip_path = ROOT / ZIP_NAME
    version = (ROOT / "VERSION").read_text(encoding="utf-8").strip()
    release_prefix = version.split("-", 1)[0]
    note_path = ROOT / f"GITHUB_RELEASE_{release_prefix}.md"
    if not artifacts_path.exists() or not zip_path.exists():
        return
    try:
        artifacts = json.loads(artifacts_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        add_message(messages, "ERROR", "release_artifacts_invalid_json",
                    f"RELEASE_ARTIFACTS.json is invalid JSON: {exc}", artifacts_path)
        return
    actual_sha = sha256_file(zip_path)
    actual_bytes = zip_path.stat().st_size
    with zipfile.ZipFile(zip_path) as zf:
        actual_entries = len(zf.infolist())
    expected = {
        "zip_sha256": actual_sha,
        "zip_bytes": actual_bytes,
        "zip_entries": actual_entries,
    }
    for key, actual in expected.items():
        if artifacts.get(key) != actual:
            add_message(messages, "ERROR", "release_artifact_mismatch",
                        f"RELEASE_ARTIFACTS.json {key} does not match zip.", artifacts_path)
    if note_path.exists():
        note = note_path.read_text(encoding="utf-8")
        if actual_sha not in note or str(actual_bytes) not in note:
            add_message(messages, "ERROR", "github_release_note_stale",
                        "GitHub release notes do not contain current zip SHA and byte size.",
                        note_path)


def check_json_metadata(messages):
    checks = {
        ".zenodo.json": {
            "access_right",
            "creators",
            "description",
            "license",
            "title",
            "upload_type",
            "version",
        },
        "codemeta.json": {
            "@context",
            "@type",
            "author",
            "codeRepository",
            "description",
            "license",
            "name",
            "version",
        },
    }
    for rel, required_keys in checks.items():
        path = ROOT / rel
        if not path.exists():
            continue
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            add_message(messages, "ERROR", "metadata_json_invalid",
                        f"{rel} is invalid JSON: {exc}", path)
            continue
        missing = sorted(key for key in required_keys if key not in data)
        if missing:
            add_message(messages, "ERROR", "metadata_key_missing",
                        f"{rel} is missing required keys: {', '.join(missing)}", path)
        if data.get("version") != (ROOT / "VERSION").read_text(encoding="utf-8").strip():
            add_message(messages, "ERROR", "metadata_version_mismatch",
                        f"{rel} version does not match VERSION.", path)


def check_workflow(messages):
    path = ROOT / ".github/workflows/ci.yml"
    if not path.exists():
        return
    text = path.read_text(encoding="utf-8")
    for marker in [
        "audit_publication_readiness.py",
        "python -m pytest",
        "run_release_acceptance.py",
        "validate_release_lite.py",
    ]:
        if marker not in text:
            add_message(messages, "ERROR", "ci_marker_missing",
                        f"CI workflow does not run expected check: {marker}", path)


def write_json_report(path, report):
    out_path = Path(path)
    if not out_path.is_absolute():
        out_path = ROOT / out_path
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")
    print(f"[publication-audit] wrote {out_path}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-tracked-mb", type=int, default=25)
    parser.add_argument("--max-history-blob-mb", type=int, default=25)
    parser.add_argument("-o", "--output", help="optional JSON report path")
    args = parser.parse_args(argv)

    messages = []
    files = tracked_files()
    check_required_files(files, messages)
    check_gitignore(messages)
    check_tracked_file_sizes(files, messages, args.max_tracked_mb * 1024 * 1024)
    check_line_endings(files, messages)
    check_secret_patterns(files, messages)
    check_history_blob_sizes(messages, args.max_history_blob_mb * 1024 * 1024)
    check_release_artifacts(messages)
    check_json_metadata(messages)
    check_workflow(messages)

    errors = [item for item in messages if item["level"] == "ERROR"]
    report = {
        "schema_version": "1.0",
        "status": "FAIL" if errors else "PASS",
        "tracked_file_count": len(files),
        "messages": messages,
    }
    for item in messages:
        path = f" ({item['path']})" if "path" in item else ""
        print(f"[publication-audit] {item['level']}: {item['code']}: {item['message']}{path}")
    print(f"[publication-audit] tracked files: {len(files)}")
    print(f"[publication-audit] status={report['status']}")
    if args.output:
        write_json_report(args.output, report)
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
