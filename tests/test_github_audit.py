"""Unit tests for hosted GitHub repository audit helpers."""
from audit_github_repository import (
    check_branch_protection,
    check_release,
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
