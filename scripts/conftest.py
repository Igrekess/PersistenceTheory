"""
Persistence Theory v7 -- pytest integration.

Collects all proof_* and test_* scripts from chapter subdirectories
and runs them as subprocess tests. Each script is an independent process,
ensuring no shared state between tests.

Usage:
    pytest                          # run all scripts
    pytest -k "ch01"                # run ch01 only
    pytest --co                     # list all discovered tests
"""

import os
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT_DIR = Path(__file__).resolve().parent

_SKIP_PREFIXES = ("pt_", "_", "__", "conftest")
_SKIP_FILES = {"run_all.py", "conftest.py"}
_SKIP_DIR_PREFIXES = (".", "__", "archive", "lib", "reports")


def _discover_domains():
    """
    Auto-discover all subdirectories containing runnable .py scripts.

    Replaces the previous manually-maintained DOMAINS list so that
    new chapter directories (e.g. Class B Extensions ch20[bcdefg]_*)
    are picked up automatically without code changes.
    """
    for entry in sorted(os.listdir(SCRIPT_DIR)):
        path = SCRIPT_DIR / entry
        if not path.is_dir():
            continue
        if entry.startswith(_SKIP_DIR_PREFIXES):
            continue
        # Confirm at least one runnable script exists in this directory
        has_script = any(
            f.endswith(".py")
            and not f.startswith("._")
            and not any(f.startswith(p) for p in _SKIP_PREFIXES)
            and f not in _SKIP_FILES
            for f in os.listdir(path) if os.path.isfile(path / f)
        )
        if has_script:
            yield entry


def _discover_scripts():
    """Yield (domain, script_path) for all runnable scripts."""
    for domain in _discover_domains():
        domain_dir = SCRIPT_DIR / domain
        for f in sorted(os.listdir(domain_dir)):
            # Skip macOS AppleDouble metadata files (._foo.py etc.)
            if f.startswith("._"):
                continue
            if not f.endswith(".py"):
                continue
            if any(f.startswith(p) for p in _SKIP_PREFIXES):
                continue
            if f in _SKIP_FILES:
                continue
            yield domain, domain_dir / f


_ALL_SCRIPTS = list(_discover_scripts())


def pytest_ignore_collect(collection_path, config):
    """
    Prevent pytest from importing chapter test_*.py and proof_*.py scripts.

    These scripts are designed to run as standalone processes (they call
    `ck.summary()` -> `sys.exit()` at module level via the Checker framework
    in lib/pt_check.py). Importing them via pytest's normal collection
    mechanism triggers an INTERNALERROR.

    Instead, every script is exposed to pytest as a parametrised case of
    `test_pt_script` (defined in test_runner.py), which invokes the script
    via `subprocess.run`. The only files pytest should collect from this
    directory are conftest.py and test_runner.py.
    """
    p = collection_path if hasattr(collection_path, "name") else None
    if p is None:
        return False
    # Allow pytest to collect conftest.py and test_runner.py at top level
    if p.parent == SCRIPT_DIR and p.name in ("conftest.py", "test_runner.py"):
        return False
    # Ignore any other Python file in chapter subdirectories
    if p.suffix == ".py":
        try:
            p.relative_to(SCRIPT_DIR)
            return True
        except ValueError:
            return False
    return False


def _make_id(domain, path):
    return f"{domain}/{path.name}"


@pytest.fixture(
    params=_ALL_SCRIPTS,
    ids=[_make_id(d, p) for d, p in _ALL_SCRIPTS]
)
def pt_script(request):
    return request.param
