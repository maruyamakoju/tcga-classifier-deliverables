#!/usr/bin/env python3
"""Audit hosted GitHub repository settings and release asset consistency."""
import argparse
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
ZIP_NAME = "tcga-tumor-normal-release-lite.zip"

EXPECTED_REQUIRED_CONTEXTS = {
    "windows-latest / py3.11",
    "ubuntu-latest / py3.11",
    "macos-latest / py3.11",
}

EXPECTED_TOPICS = {
    "bioinformatics",
    "cancer-genomics",
    "machine-learning",
    "python",
    "rnaseq",
    "tcga",
}


def add_message(messages, level, code, message, path=None):
    item = {"level": level, "code": code, "message": message}
    if path is not None:
        item["path"] = str(path)
    messages.append(item)


def run_command(args):
    try:
        result = subprocess.run(
            args,
            cwd=ROOT,
            text=True,
            capture_output=True,
        )
    except FileNotFoundError as exc:
        raise RuntimeError(f"Command not found: {args[0]}") from exc
    if result.returncode != 0:
        detail = (result.stderr or result.stdout or "").strip()
        raise RuntimeError(detail or f"{args[0]} exited with {result.returncode}")
    return result.stdout


def infer_repo_from_remote(remote):
    remote = remote.strip()
    patterns = [
        r"^https://github\.com/([^/]+/[^/.]+?)(?:\.git)?/?$",
        r"^git@github\.com:([^/]+/[^/.]+?)(?:\.git)?$",
        r"^ssh://git@github\.com/([^/]+/[^/.]+?)(?:\.git)?/?$",
    ]
    for pattern in patterns:
        match = re.match(pattern, remote)
        if match:
            return match.group(1)
    return None


def default_repo():
    try:
        remote = run_command(["git", "remote", "get-url", "origin"])
    except RuntimeError:
        return None
    return infer_repo_from_remote(remote)


def gh_json(endpoint):
    data = run_command(["gh", "api", endpoint])
    return json.loads(data)


def load_json_file(path, messages):
    if not path.exists():
        add_message(messages, "ERROR", "json_file_missing", f"Missing JSON file: {path.name}", path)
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        add_message(messages, "ERROR", "json_invalid", f"{path.name} is invalid JSON: {exc}", path)
        return None


def check_repository_metadata(repo_data, expected_repo, messages):
    if repo_data.get("full_name") != expected_repo:
        add_message(messages, "ERROR", "repo_name_mismatch",
                    f"GitHub returned {repo_data.get('full_name')} for {expected_repo}.")
    if repo_data.get("private") is not False:
        add_message(messages, "ERROR", "repo_not_public", "Repository is not public.")
    if repo_data.get("default_branch") != "main":
        add_message(messages, "ERROR", "default_branch_not_main",
                    f"Default branch is {repo_data.get('default_branch')!r}, expected 'main'.")
    if repo_data.get("archived"):
        add_message(messages, "ERROR", "repo_archived", "Repository is archived.")
    if not repo_data.get("description"):
        add_message(messages, "WARNING", "repo_description_missing",
                    "Repository description is empty.")
    license_data = repo_data.get("license") or {}
    if license_data.get("spdx_id") != "MIT":
        add_message(messages, "WARNING", "repo_license_unexpected",
                    f"Repository license is {license_data.get('spdx_id')!r}, expected 'MIT'.")


def check_topics(topics_data, messages):
    topics = set(topics_data.get("names") or [])
    missing = sorted(EXPECTED_TOPICS - topics)
    if missing:
        add_message(messages, "WARNING", "repo_topics_missing",
                    "Repository topics are missing: " + ", ".join(missing))


def check_branch_protection(protection, messages):
    checks = protection.get("required_status_checks") or {}
    contexts = set(checks.get("contexts") or [])
    missing = sorted(EXPECTED_REQUIRED_CONTEXTS - contexts)
    if checks.get("strict") is not True:
        add_message(messages, "ERROR", "required_checks_not_strict",
                    "Branch protection does not require up-to-date status checks.")
    if missing:
        add_message(messages, "ERROR", "required_context_missing",
                    "Branch protection is missing required contexts: " + ", ".join(missing))
    if (protection.get("allow_force_pushes") or {}).get("enabled") is not False:
        add_message(messages, "ERROR", "force_pushes_allowed",
                    "Branch protection allows force pushes.")
    if (protection.get("allow_deletions") or {}).get("enabled") is not False:
        add_message(messages, "ERROR", "branch_deletions_allowed",
                    "Branch protection allows branch deletion.")


def find_release_asset(release_data, name):
    for asset in release_data.get("assets") or []:
        if asset.get("name") == name:
            return asset
    return None


def check_release(release_data, artifacts, version, messages):
    if release_data.get("tag_name") != version:
        add_message(messages, "ERROR", "release_tag_mismatch",
                    f"Latest audited release tag is {release_data.get('tag_name')!r}, expected {version!r}.")
    if release_data.get("draft"):
        add_message(messages, "ERROR", "release_is_draft", "Audited release is still a draft.")
    if release_data.get("prerelease"):
        add_message(messages, "WARNING", "release_is_prerelease", "Audited release is marked prerelease.")

    asset = find_release_asset(release_data, ZIP_NAME)
    if asset is None:
        add_message(messages, "ERROR", "release_asset_missing",
                    f"Release asset is missing: {ZIP_NAME}")
        return

    expected_size = artifacts.get("zip_bytes")
    if asset.get("size") != expected_size:
        add_message(messages, "ERROR", "release_asset_size_mismatch",
                    f"Release asset size is {asset.get('size')}, expected {expected_size}.")

    expected_digest = f"sha256:{artifacts.get('zip_sha256')}"
    actual_digest = asset.get("digest")
    if actual_digest is None:
        add_message(messages, "WARNING", "release_asset_digest_missing",
                    "GitHub API did not return an asset digest; checked size only.")
    elif actual_digest.lower() != expected_digest.lower():
        add_message(messages, "ERROR", "release_asset_digest_mismatch",
                    f"Release asset digest is {actual_digest}, expected {expected_digest}.")


def check_open_pull_requests(pulls, messages):
    for pull in pulls:
        user = pull.get("user") or {}
        head = pull.get("head") or {}
        if user.get("login") == "dependabot[bot]" and head.get("ref", "").startswith("dependabot/pip"):
            add_message(messages, "ERROR", "open_pip_dependabot_pr",
                        f"Open pip Dependabot PR remains: #{pull.get('number')} {pull.get('title')}")


def build_report(repo):
    messages = []
    version = (ROOT / "VERSION").read_text(encoding="utf-8").strip()
    artifacts = load_json_file(ROOT / "RELEASE_ARTIFACTS.json", messages)
    if artifacts is None:
        artifacts = {}

    if not repo:
        add_message(messages, "ERROR", "repo_not_inferred",
                    "Could not infer GitHub repository from origin remote; pass --repo owner/name.")
    elif shutil.which("gh") is None:
        add_message(messages, "ERROR", "gh_missing",
                    "GitHub CLI 'gh' is required for hosted repository audit.")
    else:
        endpoints = {
            "repo": f"repos/{repo}",
            "topics": f"repos/{repo}/topics",
            "protection": f"repos/{repo}/branches/main/protection",
            "release": f"repos/{repo}/releases/tags/{version}",
            "pulls": f"repos/{repo}/pulls?state=open&per_page=100",
        }
        loaded = {}
        for name, endpoint in endpoints.items():
            try:
                loaded[name] = gh_json(endpoint)
            except Exception as exc:
                add_message(messages, "ERROR", f"github_{name}_query_failed",
                            f"Could not query GitHub {name}: {exc}")

        if "repo" in loaded:
            check_repository_metadata(loaded["repo"], repo, messages)
        if "topics" in loaded:
            check_topics(loaded["topics"], messages)
        if "protection" in loaded:
            check_branch_protection(loaded["protection"], messages)
        if "release" in loaded:
            check_release(loaded["release"], artifacts, version, messages)
        if "pulls" in loaded:
            check_open_pull_requests(loaded["pulls"], messages)

    levels = {item["level"] for item in messages}
    status = "FAIL" if "ERROR" in levels else "WARN" if "WARNING" in levels else "PASS"
    return {
        "schema_version": "1.0",
        "status": status,
        "repo": repo,
        "version": version,
        "messages": messages,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Audit hosted GitHub repository settings and release asset consistency."
    )
    parser.add_argument("--repo", help="GitHub repository as owner/name; defaults to origin remote")
    parser.add_argument("-o", "--output", help="write JSON report")
    parser.add_argument("--strict", action="store_true",
                        help="return non-zero on warnings as well as errors")
    args = parser.parse_args(argv)

    report = build_report(args.repo or default_repo())
    for message in report["messages"]:
        stream = sys.stderr if message["level"] in {"ERROR", "WARNING"} else sys.stdout
        print(f"[github-audit] {message['level']}: {message['message']}", file=stream)
    print(f"[github-audit] repo={report['repo']}")
    print(f"[github-audit] version={report['version']}")
    print(f"[github-audit] status={report['status']}")
    if args.output:
        out_path = Path(args.output)
        if not out_path.is_absolute():
            out_path = ROOT / out_path
        out_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                            encoding="utf-8")
        print(f"[github-audit] wrote {out_path}")
    if report["status"] == "FAIL" or (args.strict and report["status"] == "WARN"):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
