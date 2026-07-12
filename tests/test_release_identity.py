"""Tests for immutable VERSION/tag identity enforcement."""

from audit_release_identity import evaluate_release_identity


def codes(messages):
    return {message["code"] for message in messages}


def test_release_preparation_allows_version_without_existing_tag():
    messages = evaluate_release_identity("v1.2.0-gdc-starcounts", "head", None)

    assert codes(messages) == {"release_tag_not_created"}
    assert messages[0]["level"] == "INFO"


def test_existing_version_must_point_to_head():
    messages = evaluate_release_identity(
        "v1.1.22-gdc-starcounts",
        "new-head",
        "published-commit",
        "v1.1.22-gdc-starcounts",
    )

    assert "released_version_reused" in codes(messages)


def test_existing_version_at_head_with_matching_version_passes():
    assert evaluate_release_identity("v1.2.0-gdc-starcounts", "same", "same", "v1.2.0-gdc-starcounts") == []


def test_release_identity_rejects_bad_version_and_tag_contents():
    messages = evaluate_release_identity("not a version", "same", "same", "v0.0.0")

    assert codes(messages) == {"version_format_invalid", "release_tag_version_mismatch"}
