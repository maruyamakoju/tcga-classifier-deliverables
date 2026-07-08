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

EXPECTED_TAG_RULESET_REF = "refs/tags/v*"
EXPECTED_TAG_RULESET_RULES = {"deletion", "non_fast_forward", "update"}


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


def gh_status(endpoint):
    data = run_command(["gh", "api", "-i", endpoint])
    first_line = data.splitlines()[0] if data.splitlines() else ""
    match = re.match(r"^HTTP/\S+\s+(\d+)", first_line)
    if not match:
        raise RuntimeError(f"Could not parse GitHub API status line: {first_line!r}")
    return int(match.group(1))


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
    if (protection.get("enforce_admins") or {}).get("enabled") is not True:
        add_message(messages, "ERROR", "admins_not_enforced",
                    "Branch protection does not apply to administrators.")
    if (protection.get("required_linear_history") or {}).get("enabled") is not True:
        add_message(messages, "ERROR", "linear_history_not_required",
                    "Branch protection does not require linear history.")
    if (protection.get("required_conversation_resolution") or {}).get("enabled") is not True:
        add_message(messages, "ERROR", "conversation_resolution_not_required",
                    "Branch protection does not require conversation resolution before merge.")
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


def check_release_tag_rulesets(rulesets, messages):
    candidates = []
    for ruleset in rulesets:
        conditions = ruleset.get("conditions") or {}
        ref_name = conditions.get("ref_name") or {}
        includes = set(ref_name.get("include") or [])
        if ruleset.get("target") == "tag" and EXPECTED_TAG_RULESET_REF in includes:
            candidates.append(ruleset)
    if not candidates:
        add_message(messages, "ERROR", "release_tag_ruleset_missing",
                    f"No tag ruleset protects {EXPECTED_TAG_RULESET_REF}.")
        return

    protected = False
    for ruleset in candidates:
        rule_types = {rule.get("type") for rule in ruleset.get("rules") or []}
        missing_rules = sorted(EXPECTED_TAG_RULESET_RULES - rule_types)
        if ruleset.get("enforcement") != "active":
            add_message(messages, "ERROR", "release_tag_ruleset_inactive",
                        f"Tag ruleset {ruleset.get('name')!r} is not active.")
            continue
        if missing_rules:
            add_message(messages, "ERROR", "release_tag_rules_missing",
                        f"Tag ruleset {ruleset.get('name')!r} is missing rules: "
                        + ", ".join(missing_rules))
            continue
        if ruleset.get("current_user_can_bypass") != "never":
            add_message(messages, "ERROR", "release_tag_ruleset_bypass_allowed",
                        f"Current user can bypass tag ruleset {ruleset.get('name')!r}.")
            continue
        if ruleset.get("bypass_actors"):
            add_message(messages, "WARNING", "release_tag_ruleset_bypass_actors",
                        f"Tag ruleset {ruleset.get('name')!r} has bypass actors configured.")
            continue
        protected = True
    if not protected:
        add_message(messages, "ERROR", "release_tag_ruleset_not_protective",
                    f"No active tag ruleset fully protects {EXPECTED_TAG_RULESET_REF}.")


def check_vulnerability_alerts(status_code, messages):
    if status_code != 204:
        add_message(messages, "ERROR", "vulnerability_alerts_disabled",
                    f"Dependabot vulnerability alerts endpoint returned HTTP {status_code}, expected 204.")


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
            "rulesets": f"repos/{repo}/rulesets",
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
        if "rulesets" in loaded:
            details = []
            for ruleset in loaded["rulesets"]:
                ruleset_id = ruleset.get("id")
                if ruleset_id is None:
                    continue
                try:
                    details.append(gh_json(f"repos/{repo}/rulesets/{ruleset_id}"))
                except Exception as exc:
                    add_message(messages, "ERROR", "github_ruleset_detail_query_failed",
                                f"Could not query GitHub ruleset {ruleset_id}: {exc}")
            loaded["ruleset_details"] = details
        try:
            loaded["vulnerability_alerts_status"] = gh_status(
                f"repos/{repo}/vulnerability-alerts"
            )
        except Exception as exc:
            add_message(messages, "ERROR", "github_vulnerability_alerts_query_failed",
                        f"Could not query GitHub vulnerability alerts: {exc}")

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
        if "ruleset_details" in loaded:
            check_release_tag_rulesets(loaded["ruleset_details"], messages)
        if "vulnerability_alerts_status" in loaded:
            check_vulnerability_alerts(loaded["vulnerability_alerts_status"], messages)

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
