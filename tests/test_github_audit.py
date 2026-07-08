"""Unit tests for hosted GitHub repository audit helpers."""
from audit_github_repository import (
    check_branch_protection,
    check_release_tag_rulesets,
    check_release,
    check_vulnerability_alerts,
    infer_repo_from_remote,
)


def test_infer_repo_from_common_github_remotes():
    assert (
        infer_repo_from_remote("https://github.com/owner/repo.git")
        == "owner/repo"
    )
    assert infer_repo_from_remote("git@github.com:owner/repo.git") == "owner/repo"
    assert (
        infer_repo_from_remote("ssh://git@github.com/owner/repo.git")
        == "owner/repo"
    )
    assert infer_repo_from_remote("https://example.com/owner/repo.git") is None


def test_branch_protection_requires_all_expected_contexts():
    messages = []
    check_branch_protection(
        {
            "required_status_checks": {
                "strict": False,
                "contexts": ["ubuntu-latest / py3.11"],
            },
            "allow_force_pushes": {"enabled": True},
            "allow_deletions": {"enabled": False},
        },
        messages,
    )
    codes = {message["code"] for message in messages}
    assert "required_checks_not_strict" in codes
    assert "required_context_missing" in codes
    assert "admins_not_enforced" in codes
    assert "linear_history_not_required" in codes
    assert "conversation_resolution_not_required" in codes
    assert "force_pushes_allowed" in codes


def test_release_audit_detects_asset_mismatch():
    messages = []
    check_release(
        {
            "tag_name": "v1",
            "draft": False,
            "prerelease": False,
            "assets": [
                {
                    "name": "tcga-tumor-normal-release-lite.zip",
                    "size": 9,
                    "digest": "sha256:bad",
                }
            ],
        },
        {"zip_bytes": 10, "zip_sha256": "abc"},
        "v1",
        messages,
    )
    codes = {message["code"] for message in messages}
    assert "release_asset_size_mismatch" in codes
    assert "release_asset_digest_mismatch" in codes


def test_release_tag_ruleset_requires_active_v_tag_protection():
    messages = []
    check_release_tag_rulesets(
        [
            {
                "name": "Protect release tags",
                "target": "tag",
                "enforcement": "active",
                "conditions": {"ref_name": {"include": ["refs/tags/v*"]}},
                "rules": [
                    {"type": "deletion"},
                    {"type": "non_fast_forward"},
                    {"type": "update"},
                ],
                "current_user_can_bypass": "never",
                "bypass_actors": [],
            }
        ],
        messages,
    )
    assert messages == []

    messages = []
    check_release_tag_rulesets(
        [
            {
                "name": "Weak tag rules",
                "target": "tag",
                "enforcement": "disabled",
                "conditions": {"ref_name": {"include": ["refs/tags/v*"]}},
                "rules": [{"type": "deletion"}],
                "current_user_can_bypass": "pull_request",
            }
        ],
        messages,
    )
    codes = {message["code"] for message in messages}
    assert "release_tag_ruleset_inactive" in codes
    assert "release_tag_ruleset_not_protective" in codes


def test_vulnerability_alerts_require_204_status():
    messages = []
    check_vulnerability_alerts(204, messages)
    assert messages == []
    check_vulnerability_alerts(404, messages)
    assert {message["code"] for message in messages} == {"vulnerability_alerts_disabled"}
