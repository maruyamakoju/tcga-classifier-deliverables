#!/usr/bin/env python3
"""Prevent reuse of an existing release VERSION for a different commit."""
import argparse
import re
import subprocess
import sys
from pathlib import Path

from release_tools.common import (
    add_message,
    exit_code_for_status,
    status_from_levels,
    write_json_report,
)


ROOT = Path(__file__).resolve().parent
VERSION_RE = re.compile(r"^v\d+\.\d+\.\d+(?:-[0-9A-Za-z][0-9A-Za-z.-]*)?$")


def run_git(args, root=ROOT, check=True):
    try:
        result = subprocess.run(
            ["git", *args],
            cwd=root,
            text=True,
            capture_output=True,
            timeout=30,
        )
    except (FileNotFoundError, subprocess.TimeoutExpired) as exc:
        raise RuntimeError(f"Could not run git {' '.join(args)}: {exc}") from exc
    if check and result.returncode != 0:
        detail = (result.stderr or result.stdout).strip()
        raise RuntimeError(detail or f"git {' '.join(args)} failed")
    return result


def evaluate_release_identity(version, head_sha, tag_sha, tag_version=None):
    messages = []
    if not version:
        add_message(messages, "ERROR", "version_empty", "VERSION must be non-empty.")
        return messages
    if not VERSION_RE.fullmatch(version):
        add_message(
            messages,
            "ERROR",
            "version_format_invalid",
            f"VERSION is not a supported release identifier: {version!r}.",
        )
    if tag_sha is None:
        add_message(
            messages,
            "INFO",
            "release_tag_not_created",
            f"Tag {version!r} does not exist yet; release preparation is allowed.",
        )
        return messages
    if tag_version is not None and tag_version != version:
        add_message(
            messages,
            "ERROR",
            "release_tag_version_mismatch",
            f"Tag {version!r} contains VERSION {tag_version!r}.",
        )
    if tag_sha != head_sha:
        add_message(
            messages,
            "ERROR",
            "released_version_reused",
            f"VERSION {version!r} already tags {tag_sha}, but HEAD is {head_sha}; "
            "bump VERSION instead of changing a published release.",
        )
    return messages


def build_report(root=ROOT):
    root = Path(root)
    messages = []
    version_path = root / "VERSION"
    try:
        version = version_path.read_text(encoding="utf-8").strip()
    except (OSError, UnicodeError) as exc:
        add_message(
            messages,
            "ERROR",
            "version_unreadable",
            f"Could not read VERSION: {exc}",
            version_path,
        )
        version = ""
    head_sha = None
    tag_sha = None
    tag_version = None
    try:
        head_sha = run_git(["rev-parse", "HEAD"], root=root).stdout.strip()
        if version:
            tag_result = run_git(
                ["rev-parse", "--verify", f"refs/tags/{version}^{{commit}}"],
                root=root,
                check=False,
            )
            if tag_result.returncode == 0:
                tag_sha = tag_result.stdout.strip()
                tag_version = run_git(
                    ["show", f"{tag_sha}:VERSION"], root=root
                ).stdout.strip()
    except RuntimeError as exc:
        add_message(messages, "ERROR", "git_query_failed", str(exc))
    if head_sha is not None:
        messages.extend(evaluate_release_identity(version, head_sha, tag_sha, tag_version))
    return {
        "schema_version": "1.0",
        "status": status_from_levels(messages),
        "root": str(root.resolve()),
        "version": version,
        "head_sha": head_sha,
        "tag_sha": tag_sha,
        "tag_version": tag_version,
        "messages": messages,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o", "--output", help="optional JSON report path")
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args(argv)
    report = build_report()
    for message in report["messages"]:
        stream = sys.stderr if message["level"] in {"ERROR", "WARNING"} else sys.stdout
        print(f"[release-identity] {message['level']}: {message['message']}", file=stream)
    print(f"[release-identity] status={report['status']}")
    if args.output:
        write_json_report(args.output, report, root=ROOT, prefix="release-identity")
    return exit_code_for_status(report["status"], strict=args.strict)


if __name__ == "__main__":
    raise SystemExit(main())
