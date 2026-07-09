"""Unit tests for release documentation audit output formatting."""
from audit_release_docs import print_report


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
