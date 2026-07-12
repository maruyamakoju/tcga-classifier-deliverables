"""Unit tests for release documentation auditing."""
import audit_release_docs
from audit_release_docs import check_publication_release_note_reference, print_report


def sample_report():
    return {
        "status": "FAIL",
        "messages": [
            {"level": "INFO", "message": "intentional full-tree reference"},
            {"level": "ERROR", "message": "missing local path"},
        ],
    }


def test_print_report_hides_info_by_default(capsys):
    print_report(sample_report())

    captured = capsys.readouterr()
    assert "intentional full-tree reference" not in captured.out
    assert "info_messages=1 hidden" in captured.out
    assert "status=FAIL" in captured.out
    assert "missing local path" in captured.err


def test_print_report_can_show_info(capsys):
    print_report(sample_report(), show_info=True)

    captured = capsys.readouterr()
    assert "intentional full-tree reference" in captured.out
    assert "info_messages=" not in captured.out
    assert "missing local path" in captured.err


def test_publication_checklist_must_reference_current_release_note(tmp_path, monkeypatch):
    monkeypatch.setattr(audit_release_docs, "ROOT", tmp_path)
    (tmp_path / "VERSION").write_text("v1.2.0-gdc-starcounts\n", encoding="utf-8")
    checklist = tmp_path / "PUBLICATION_CHECKLIST.md"
    checklist.write_text("Paste `GITHUB_RELEASE_v1.1.20.md`.\n", encoding="utf-8")
    messages = []

    check_publication_release_note_reference(messages)

    assert {message["code"] for message in messages} == {"publication_release_note_stale"}
    checklist.write_text("Paste `GITHUB_RELEASE_v1.2.0.md`.\n", encoding="utf-8")
    messages = []
    check_publication_release_note_reference(messages)
    assert messages == []
