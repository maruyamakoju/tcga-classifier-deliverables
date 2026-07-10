#!/usr/bin/env python3
"""Audit release Python imports against the lightweight runtime contract."""
import argparse
import ast
import re
import sys
from pathlib import Path

from release_tools.common import (
    add_message,
    exit_code_for_status,
    release_python_files as python_files,
    release_target_root as _release_target_root,
    status_from_levels,
    write_json_report,
)


ROOT = Path(__file__).resolve().parent


def release_target_root():
    return _release_target_root(ROOT)

ALLOWED_EXTERNAL = {"numpy", "pandas", "pyarrow"}
BANNED_IMPORT_ROOTS = {
    "joblib",
    "lightgbm",
    "matplotlib",
    "scipy",
    "sklearn",
    "tensorflow",
    "torch",
    "xgboost",
}
IGNORED_MODULE_ROOTS = {"__future__"}
MINIMAL_REQUIREMENT_ROOTS = {"numpy", "pandas", "pyarrow"}


def stdlib_roots():
    names = set(getattr(sys, "stdlib_module_names", set()))
    names.update(sys.builtin_module_names)
    # Common implementation/private roots that can appear on some platforms.
    names.update({"nt", "posix"})
    return names


def local_module_roots(target_root, files):
    roots = {path.stem for path in files}
    for package_init in target_root.rglob("__init__.py"):
        roots.add(package_init.parent.name)
    return roots


def import_root(name):
    return name.split(".", 1)[0]


def literal_string(node):
    return node.value if isinstance(node, ast.Constant) and isinstance(node.value, str) else None


def collect_imports(path, messages):
    try:
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    except SyntaxError as exc:
        add_message(messages, "ERROR", "python_parse_failed",
                    f"Could not parse {path.name}: {exc}", path)
        return []
    imports = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                imports.append((alias.name, node.lineno, "import"))
        elif isinstance(node, ast.ImportFrom):
            if node.level == 0 and node.module:
                imports.append((node.module, node.lineno, "from"))
        elif isinstance(node, ast.Call):
            func = node.func
            name = None
            if isinstance(func, ast.Name) and func.id == "__import__" and node.args:
                name = literal_string(node.args[0])
            elif (
                isinstance(func, ast.Attribute)
                and func.attr == "import_module"
                and isinstance(func.value, ast.Name)
                and func.value.id == "importlib"
                and node.args
            ):
                name = literal_string(node.args[0])
            if name:
                imports.append((name, node.lineno, "dynamic"))
    return imports


def check_python_imports(target_root, messages):
    files = python_files(target_root)
    local_roots = local_module_roots(target_root, files)
    stdlib = stdlib_roots()
    rows = []
    for path in files:
        for module, lineno, kind in collect_imports(path, messages):
            root = import_root(module)
            rel = path.relative_to(target_root).as_posix()
            rows.append({"file": rel, "line": lineno, "kind": kind, "module": module})
            if root in IGNORED_MODULE_ROOTS:
                continue
            if root in BANNED_IMPORT_ROOTS:
                add_message(messages, "ERROR", "banned_dependency_import",
                            f"{rel}:{lineno} imports banned dependency root '{root}'.", path)
                continue
            if root in stdlib or root in local_roots or root in ALLOWED_EXTERNAL:
                continue
            add_message(messages, "ERROR", "undeclared_external_import",
                        f"{rel}:{lineno} imports '{module}', which is not in the lightweight contract.",
                        path)
    return files, rows


def requirement_root(line):
    line = line.split("#", 1)[0].strip()
    if not line or line.startswith("-"):
        return None
    match = re.match(r"([A-Za-z0-9_.-]+)", line)
    if not match:
        return None
    return match.group(1).lower().replace("-", "_")


def check_requirements_light(target_root, messages):
    path = target_root / "requirements-light.txt"
    if not path.exists():
        add_message(messages, "ERROR", "requirements_light_missing",
                    "requirements-light.txt is missing.", path)
        return []
    roots = []
    for line in path.read_text(encoding="utf-8").splitlines():
        root = requirement_root(line)
        if root:
            roots.append(root)
    missing = sorted(MINIMAL_REQUIREMENT_ROOTS - set(roots))
    extra_banned = sorted(set(roots) & BANNED_IMPORT_ROOTS)
    if missing:
        add_message(messages, "ERROR", "minimal_requirement_missing",
                    f"requirements-light.txt missing expected packages: {missing}", path)
    if extra_banned:
        add_message(messages, "ERROR", "banned_light_requirement",
                    f"requirements-light.txt includes banned packages: {extra_banned}", path)
    return roots


def build_report():
    messages = []
    target_root = release_target_root()
    files, imports = check_python_imports(target_root, messages)
    requirements = check_requirements_light(target_root, messages)
    return {
        "schema_version": "1.0",
        "status": status_from_levels(messages),
        "root": str(ROOT),
        "target_root": str(target_root),
        "allowed_external": sorted(ALLOWED_EXTERNAL),
        "banned_import_roots": sorted(BANNED_IMPORT_ROOTS),
        "python_files_checked": [path.relative_to(target_root).as_posix() for path in files],
        "requirements_light_roots": requirements,
        "imports": imports,
        "messages": messages,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Audit release imports against the lightweight dependency contract."
    )
    parser.add_argument("-o", "--output", help="write JSON report")
    parser.add_argument("--strict", action="store_true",
                        help="return non-zero on warnings as well as errors")
    args = parser.parse_args(argv)

    report = build_report()
    for message in report["messages"]:
        stream = sys.stderr if message["level"] in {"ERROR", "WARNING"} else sys.stdout
        print(f"[deps-audit] {message['level']}: {message['message']}", file=stream)
    print(f"[deps-audit] checked {len(report['python_files_checked'])} Python files")
    print(f"[deps-audit] status={report['status']}")
    if args.output:
        write_json_report(args.output, report, root=ROOT, prefix="deps-audit")
    return exit_code_for_status(report["status"], strict=args.strict)


if __name__ == "__main__":
    raise SystemExit(main())
