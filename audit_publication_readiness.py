#!/usr/bin/env python3
"""Audit repository state before making a release repository public."""
import argparse
import json
import re
import subprocess
import zipfile
from pathlib import Path

from release_tools.common import (
    RELEASE_SCHEMA_VERSION,
    RELEASE_VALIDATION_COMMAND,
    RELEASE_ZIP_NAME,
    ZIP_ACCEPTANCE_COMMAND,
    add_message,
    sha256_file,
    write_json_report,
)


ROOT = Path(__file__).resolve().parent
ZIP_NAME = RELEASE_ZIP_NAME

MAX_TRACKED_BYTES = 25 * 1024 * 1024
MAX_HISTORY_BLOB_BYTES = 25 * 1024 * 1024

BINARY_SUFFIXES = {".npy", ".npz", ".pkl", ".png", ".zip"}
TEXT_NAMES = {"LICENSE", "VERSION"}

REQUIRED_FILES = {
    ".github/workflows/ci.yml",
    ".github/dependabot.yml",
    ".github/ISSUE_TEMPLATE/bug_report.yml",
    ".github/ISSUE_TEMPLATE/config.yml",
    ".github/ISSUE_TEMPLATE/question.yml",
    ".github/PULL_REQUEST_TEMPLATE.md",
    ".gitattributes",
    ".gitignore",
    ".zenodo.json",
    "CITATION.cff",
    "CONTRIBUTING.md",
    "environment.yml",
    "LICENSE",
    "MAINTENANCE.md",
    "NOTICE.md",
    "PUBLICATION_CHECKLIST.md",
    "README.md",
    "RELEASE_ARTIFACTS.json",
    "requirements-dev.txt",
    "requirements-external-validation.txt",
    "requirements-light.txt",
    "requirements-training.txt",
    "requirements.txt",
    "SECURITY.md",
    "VERSION",
    "audit_github_repository.py",
    "audit_release_identity.py",
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


def current_github_release_note_name():
    """GITHUB_RELEASE_<prefix>.md for the current VERSION (e.g. GITHUB_RELEASE_v1.1.22.md).

    Historical release notes for older versions are never deleted, so this
    must be derived from VERSION rather than hardcoded -- a hardcoded name
    only fails once the file it names has actually been removed, silently
    skipping the real "does the current release note exist" check.

    Returns None if VERSION itself is missing; VERSION is already in
    REQUIRED_FILES, so that gets reported as a normal missing-file message
    instead of this function raising before that check runs.
    """
    version_path = ROOT / "VERSION"
    if not version_path.exists():
        return None
    version = version_path.read_text(encoding="utf-8").strip()
    release_prefix = version.split("-", 1)[0]
    return f"GITHUB_RELEASE_{release_prefix}.md"


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


def check_required_files(files, messages):
    file_set = set(files)
    required = set(REQUIRED_FILES)
    note_name = current_github_release_note_name()
    if note_name is not None:
        required.add(note_name)
    for rel in sorted(required):
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
    if not artifacts_path.exists() or not zip_path.exists():
        return
    note_name = current_github_release_note_name()
    note_path = ROOT / note_name if note_name is not None else None
    try:
        artifacts = json.loads(artifacts_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        add_message(messages, "ERROR", "release_artifacts_invalid_json",
                    f"RELEASE_ARTIFACTS.json is invalid JSON: {exc}", artifacts_path)
        return
    if not isinstance(artifacts, dict):
        add_message(messages, "ERROR", "release_artifacts_not_object",
                    "RELEASE_ARTIFACTS.json top-level value must be an object.",
                    artifacts_path)
        return
    actual_sha = sha256_file(zip_path)
    actual_bytes = zip_path.stat().st_size
    try:
        with zipfile.ZipFile(zip_path) as zf:
            actual_entries = len(zf.infolist())
    except (OSError, zipfile.BadZipFile) as exc:
        add_message(messages, "ERROR", "release_zip_invalid",
                    f"Could not read release ZIP: {exc}", zip_path)
        return
    version_path = ROOT / "VERSION"
    metadata_path = ROOT / "RELEASE_METADATA.json"
    version = None
    release_date = None
    if version_path.is_file():
        version = version_path.read_text(encoding="utf-8").strip()
    if metadata_path.is_file():
        try:
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            if isinstance(metadata, dict):
                release_date = metadata.get("release_date")
        except (OSError, UnicodeError, json.JSONDecodeError):
            pass
    expected = {
        "schema_version": RELEASE_SCHEMA_VERSION,
        "version": version,
        "release_date": release_date,
        "release_dir": "release-lite",
        "zip_sha256": actual_sha,
        "zip_bytes": actual_bytes,
        "zip_entries": actual_entries,
        "zip_path": ZIP_NAME,
        "validation_command": RELEASE_VALIDATION_COMMAND,
        "zip_acceptance_command": (
            f"{ZIP_ACCEPTANCE_COMMAND} --expected-sha256 {actual_sha}"
        ),
    }
    for key, actual in expected.items():
        if artifacts.get(key) != actual:
            add_message(messages, "ERROR", "release_artifact_mismatch",
                        f"RELEASE_ARTIFACTS.json {key} does not match zip.", artifacts_path)
    if note_path is not None and note_path.exists():
        note = note_path.read_text(encoding="utf-8")
        if actual_sha not in note or str(actual_bytes) not in note:
            add_message(messages, "ERROR", "github_release_note_stale",
                        "GitHub release notes do not contain current zip SHA and byte size.",
                        note_path)


def check_json_metadata(messages):
    version_path = ROOT / "VERSION"
    if not version_path.is_file():
        # check_required_files already emits the primary diagnostic. Avoid a
        # secondary FileNotFoundError while collecting the remaining findings.
        return
    try:
        version = version_path.read_text(encoding="utf-8").strip()
    except (OSError, UnicodeError) as exc:
        add_message(messages, "ERROR", "version_unreadable",
                    f"VERSION could not be read: {exc}", version_path)
        return
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
        except (OSError, UnicodeError, json.JSONDecodeError) as exc:
            add_message(messages, "ERROR", "metadata_json_invalid",
                         f"{rel} is invalid JSON: {exc}", path)
            continue
        if not isinstance(data, dict):
            add_message(messages, "ERROR", "metadata_not_object",
                        f"{rel} top-level value must be an object.", path)
            continue
        missing = sorted(key for key in required_keys if key not in data)
        if missing:
            add_message(messages, "ERROR", "metadata_key_missing",
                        f"{rel} is missing required keys: {', '.join(missing)}", path)
        if data.get("version") != version:
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
        "build_release_lite.py",
        "audit_release_identity.py",
        "requirements-dev.txt",
        "requirements-training.txt",
    ]:
        if marker not in text:
            add_message(messages, "ERROR", "ci_marker_missing",
                        f"CI workflow does not run expected check: {marker}", path)


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
        write_json_report(args.output, report, root=ROOT, prefix="publication-audit")
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
