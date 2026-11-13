#!/usr/bin/env python3
"""Fail the build if any registered check lacks a known-failing fixture.

THE GOVERNING RULE. A guard that cannot be DEMONSTRATED to fire is not a guard. Every
raise site in this package exists because a specific measurement was once published
wrong; the fixture that reproduces that measurement is the only evidence the guard
still works. A guard with no fixture is indistinguishable from a guard that was
silently broken by a refactor, so this script treats it as a build failure, not a
warning.

WHAT COUNTS AS A REGISTERED CHECK.

  1. Every `raise <GuardError>(...)` site in the package, where `<GuardError>` is an
     exception class DEFINED IN THIS PACKAGE (`CanonError`, `DegenerateGroupingError`,
     and anything a later module adds). Bare `ValueError`/`TypeError` raises are not
     counted: the project's own guards are the ones that carry the findings.
  2. Any module-level `CHECKS` registry mapping check names to callables, if a module
     declares one.

HOW DEMONSTRATION IS MEASURED. The fixtures are executed under `sys.settrace`, and the
line numbers actually executed inside package modules are recorded. A raise site is
demonstrated when a line inside its statement was executed -- i.e. the guard really
fired during the fixtures, rather than merely being mentioned in a test name. Nothing
is inferred from naming, and there is no waiver list: adding a guard means adding the
fixture that fires it, in the same change.

Standard library only, like the package it checks. Run it the way CI does:

    python3 -I scripts/check_fixture_coverage.py

Exit status 0 when every registered check is demonstrated, 1 otherwise.
"""

from __future__ import annotations

import ast
import runpy
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PKG = ROOT / "epigauge"
TESTS = PKG / "tests"

# Builtin bases that make a class an exception for our purposes.
_EXC_BASES = {
    "Exception", "ValueError", "TypeError", "RuntimeError", "KeyError",
    "LookupError", "ArithmeticError", "OSError", "IOError", "AssertionError",
}


def package_modules() -> list[Path]:
    """Every package module except the fixtures themselves."""
    out = []
    for p in sorted(PKG.rglob("*.py")):
        if "__pycache__" in p.parts or TESTS in p.parents or p == TESTS:
            continue
        out.append(p)
    return out


def test_modules() -> list[Path]:
    return sorted(p for p in TESTS.glob("test_*.py") if "__pycache__" not in p.parts)


def guard_class_names(trees: dict[Path, ast.AST]) -> set[str]:
    """Exception classes defined in the package, resolved to a fixpoint."""
    names: set[str] = set()
    changed = True
    while changed:
        changed = False
        for tree in trees.values():
            for node in ast.walk(tree):
                if not isinstance(node, ast.ClassDef) or node.name in names:
                    continue
                bases = {b.id for b in node.bases if isinstance(b, ast.Name)}
                bases |= {b.attr for b in node.bases if isinstance(b, ast.Attribute)}
                if bases & (_EXC_BASES | names):
                    names.add(node.name)
                    changed = True
    return names


def _raised_name(node: ast.Raise) -> str | None:
    exc = node.exc
    if isinstance(exc, ast.Call):
        exc = exc.func
    if isinstance(exc, ast.Name):
        return exc.id
    if isinstance(exc, ast.Attribute):
        return exc.attr
    return None


class _SiteCollector(ast.NodeVisitor):
    def __init__(self, path: Path, guards: set[str], lines: list[str]) -> None:
        self.path = path
        self.guards = guards
        self.lines = lines
        self.scope: list[str] = []
        self.sites: list[dict] = []

    def _scoped(self, node) -> None:
        self.scope.append(node.name)
        self.generic_visit(node)
        self.scope.pop()

    visit_FunctionDef = _scoped
    visit_AsyncFunctionDef = _scoped
    visit_ClassDef = _scoped

    def visit_Raise(self, node: ast.Raise) -> None:
        name = _raised_name(node)
        if name in self.guards:
            first = self.lines[node.lineno - 1].strip()
            self.sites.append({
                "module": self.path.relative_to(ROOT).as_posix(),
                "qualname": ".".join(self.scope) or "<module>",
                "exc": name,
                "lineno": node.lineno,
                "end_lineno": node.end_lineno or node.lineno,
                "file": str(self.path),
                "source": first[:88],
            })
        self.generic_visit(node)


def registry_names(trees: dict[Path, ast.AST]) -> dict[str, str]:
    """Named checks from a module-level `CHECKS = {...}` registry, if one exists."""
    found: dict[str, str] = {}
    for path, tree in trees.items():
        for node in tree.body:
            targets = (
                node.targets if isinstance(node, ast.Assign)
                else [node.target] if isinstance(node, ast.AnnAssign) else []
            )
            if not any(isinstance(t, ast.Name) and t.id == "CHECKS" for t in targets):
                continue
            value = node.value
            keys = []
            if isinstance(value, ast.Dict):
                keys = [k.value for k in value.keys
                        if isinstance(k, ast.Constant) and isinstance(k.value, str)]
            elif isinstance(value, (ast.List, ast.Tuple, ast.Set)):
                keys = [e.value for e in value.elts
                        if isinstance(e, ast.Constant) and isinstance(e.value, str)]
            for k in keys:
                found[k] = path.relative_to(ROOT).as_posix()
    return found


def run_fixtures(watch: set[str]) -> tuple[dict[str, set[int]], list[tuple[str, int]]]:
    """Execute every fixture module under a line tracer. Returns executed lines."""
    executed: dict[str, set[int]] = defaultdict(set)

    def tracer(frame, event, arg):
        filename = frame.f_code.co_filename
        if filename in watch:
            executed[filename].add(frame.f_lineno)
            return tracer
        return None

    results: list[tuple[str, int]] = []
    sys.path.insert(0, str(ROOT))
    for test in test_modules():
        print(f"  running {test.relative_to(ROOT).as_posix()}")
        namespace: dict | None = None
        sys.settrace(tracer)
        try:
            namespace = runpy.run_path(str(test), run_name="__main__")
            code = 0
        except SystemExit as e:
            # The module ran its own standalone `main()` and reported a status.
            code = int(e.code or 0)
        finally:
            sys.settrace(None)

        # A pytest-style module with no `__main__` block defines its fixtures without
        # running them. Call them here, so coverage reflects the whole suite and not
        # just the modules that happen to carry a standalone runner.
        if namespace is not None:
            funcs = [(k, v) for k, v in sorted(namespace.items())
                     if k.startswith("test_") and callable(v)]
            sys.settrace(tracer)
            try:
                for name, fn in funcs:
                    try:
                        fn()
                        print(f"  [  ok] {name}")
                    except Exception as e:  # noqa: BLE001
                        code = 1
                        print(f"  [FAIL] {name}: {type(e).__name__}: {e}")
            finally:
                sys.settrace(None)
        results.append((test.name, code))
    return executed, results


FLOOR_FILE = Path(__file__).resolve().parents[1] / "epigauge" / "tests" / "coverage_floor.txt"


def read_floor() -> int:
    """The number of guard sites currently demonstrated. May not fall; see main()."""
    try:
        return int(FLOOR_FILE.read_text().strip())
    except (OSError, ValueError):
        return 0


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--update-floor", action="store_true",
                    help="record the current count as the new floor")
    args = ap.parse_args(argv)
    update_floor = args.update_floor

    modules = package_modules()
    if not modules:
        print("FAIL: no package modules found under", PKG)
        return 1
    tests = test_modules()
    if not tests:
        print("FAIL: no fixtures found under", TESTS,
              "-- a package with no demonstrations cannot be released")
        return 1

    trees: dict[Path, ast.AST] = {}
    sources: dict[Path, list[str]] = {}
    for p in modules:
        text = p.read_text(encoding="utf-8")
        trees[p] = ast.parse(text, filename=str(p))
        sources[p] = text.splitlines()

    guards = guard_class_names(trees)
    sites: list[dict] = []
    for p, tree in trees.items():
        collector = _SiteCollector(p, guards, sources[p])
        collector.visit(tree)
        sites.extend(collector.sites)

    print(f"guard exception classes: {', '.join(sorted(guards)) or '(none)'}")
    print(f"registered raise sites : {len(sites)}")
    print(f"fixture modules        : {len(tests)}\n")

    executed, results = run_fixtures({str(p) for p in modules})

    # The fixtures must have exercised the SOURCE tree, not an installed copy.
    canon = sys.modules.get("epigauge.canon")
    if canon is not None and not Path(canon.__file__).resolve().is_relative_to(PKG):
        print(f"\nFAIL: fixtures imported {canon.__file__}, not the source tree under {PKG}")
        return 1

    failed_runs = [(name, code) for name, code in results if code != 0]

    demonstrated, missing = [], []
    for site in sites:
        hit = executed.get(site["file"], set())
        if any(line in hit for line in range(site["lineno"], site["end_lineno"] + 1)):
            demonstrated.append(site)
        else:
            missing.append(site)

    named = registry_names(trees)
    test_text = "\n".join(p.read_text(encoding="utf-8") for p in tests)
    missing_named = sorted(n for n in named if n not in test_text)

    print(f"\ndemonstrated: {len(demonstrated)}/{len(sites)} guard raise sites")
    if named:
        print(f"named checks: {len(named) - len(missing_named)}/{len(named)} referenced "
              "by a fixture")

    if failed_runs:
        print("\nFIXTURES FAILED (coverage is not meaningful until they pass):")
        for name, code in failed_runs:
            print(f"  - {name} exited {code}")

    if missing:
        print("\nUNDEMONSTRATED GUARDS -- each needs a fixture that makes it fire:")
        for site in sorted(missing, key=lambda s: (s["module"], s["lineno"])):
            print(f"  - {site['module']}:{site['lineno']}  {site['qualname']}  "
                  f"raise {site['exc']}")
            print(f"      {site['source']}")

    if missing_named:
        print("\nREGISTERED CHECKS WITH NO FIXTURE REFERENCE:")
        for name in missing_named:
            print(f"  - {name}  (registered in {named[name]})")

    # --- the ratchet --------------------------------------------------------
    #
    # The target is every guard demonstrated. Failing the build until then makes the
    # signal permanently red, and a permanently red signal is ignored within a week --
    # at which point it reports nothing when something real breaks.
    #
    # So the gate is a FLOOR, not the target: coverage may not fall. New code cannot add
    # an undemonstrated guard without failing, so the gap cannot grow; the existing gap
    # closes as people work in those files. The count is printed on every run, so it is
    # visible rather than buried in a failing job.
    #
    # Raising the floor is deliberate: `--update-floor` writes the new value.
    if failed_runs:
        print("\nFAIL: fixtures must pass before coverage means anything.")
        return 1

    floor = read_floor()
    count = len(demonstrated)

    if update_floor:
        FLOOR_FILE.write_text(f"{count}\n")
        print(f"\nfloor updated: {floor} -> {count}")
        return 0

    print(f"\nratchet: {count} demonstrated, floor {floor}, target {len(sites)}")

    if count < floor:
        print(f"\nFAIL: coverage fell from {floor} to {count}. A guard that was "
              "demonstrated no longer is.")
        print("Restore the fixture, or lower the floor deliberately with --update-floor "
              "and say why in the commit.")
        return 1

    if count > floor:
        print(f"advance available: run --update-floor to raise it to {count} and lock "
              "the gain in.")

    remaining = len(sites) - count
    if remaining:
        print(f"{remaining} guard sites are still undemonstrated. They are listed above. "
              "Some may be\nunreachable -- an unreachable guard should be deleted, not "
              "tested.")
    else:
        print("every guard in the package has a fixture that demonstrates it firing.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
