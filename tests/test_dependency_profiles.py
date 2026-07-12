"""Keep CI/test imports aligned with the declared dependency profiles."""
from pathlib import Path

import check_environment


ROOT = Path(__file__).resolve().parents[1]


def requirement_names(path, seen=None):
    seen = set() if seen is None else seen
    path = path.resolve()
    if path in seen:
        return set()
    seen.add(path)
    names = set()
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.split("#", 1)[0].strip()
        if line.startswith(("-r ", "--requirement ")):
            included = line.split(maxsplit=1)[1]
            names.update(requirement_names(path.parent / included, seen))
            continue
        if not line or line.startswith("-"):
            continue
        name = line
        for marker in ["<", ">", "=", "!", "~", "[", ";"]:
            name = name.split(marker, 1)[0]
        names.add(name.strip().lower().replace("-", "_"))
    return names


def test_external_validation_profile_declares_imported_packages():
    names = requirement_names(ROOT / "requirements-external-validation.txt")

    assert {"numpy", "pandas", "pyarrow", "requests"} <= names
    assert "scikit_learn" not in names


def test_development_profile_declares_test_tools_and_external_profile():
    text = (ROOT / "requirements-dev.txt").read_text(encoding="utf-8")
    names = requirement_names(ROOT / "requirements-dev.txt")

    assert "-r requirements-external-validation.txt" in text
    assert "-r requirements-training.txt" in text
    assert {"pytest", "ruff", "scipy", "scikit_learn"} <= names


def test_training_profile_pins_canonical_scientific_stack():
    lines = {
        line.strip()
        for line in (ROOT / "requirements-training.txt").read_text(
            encoding="utf-8"
        ).splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    }

    assert lines == {
        "numpy==1.26.4",
        "pandas==2.3.3",
        "scipy==1.15.3",
        "scikit-learn==1.8.0",
    }


def test_feature_export_profile_is_exact_and_separate_from_model_fitting():
    lines = {
        line.strip()
        for line in (
            ROOT / "training_tools" / "requirements-feature-export.txt"
        ).read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    }

    assert lines == {"numpy==2.4.1", "pandas==3.0.0"}


def test_full_runtime_composes_external_and_canonical_training_profiles():
    text = (ROOT / "requirements.txt").read_text(encoding="utf-8")

    assert "-r requirements-external-validation.txt" in text
    assert "-r requirements-training.txt" in text


def test_conda_environment_matches_canonical_training_stack():
    text = (ROOT / "environment.yml").read_text(encoding="utf-8")

    for requirement in [
        "python=3.11",
        "numpy=1.26.4",
        "pandas=2.3.3",
        "scipy=1.15.3",
        "scikit-learn=1.8.0",
    ]:
        assert requirement in text


def test_ci_separates_full_unit_and_lightweight_bundle_profiles():
    workflow = (ROOT / ".github/workflows/ci.yml").read_text(encoding="utf-8")

    assert "full-unit:" in workflow
    assert "test-release:" in workflow
    assert "pip install -r requirements-dev.txt" in workflow
    assert "pip install -r requirements-light.txt" in workflow
    assert "build_release_lite.py --check" in workflow
    assert "name: full-unit / py3.11" in workflow
    assert 'python-version: ["3.11", "3.13"]' in workflow


def test_environment_contract_requires_supported_python(monkeypatch):
    assert check_environment.MIN_PYTHON == (3, 11)
    messages = []
    monkeypatch.setattr(check_environment.sys, "version_info", (3, 10, 14))

    check_environment.check_python(messages)

    assert messages == [
        {
            "level": "ERROR",
            "code": "python_too_old",
            "message": "Python 3.10 detected; use Python >=3.11.",
        }
    ]


def test_pytest_configuration_has_no_undeclared_plugin_options():
    config = (ROOT / "pytest.ini").read_text(encoding="utf-8")

    assert "asyncio_" not in config
    assert "--strict-config" in config
    assert "--strict-markers" in config
