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


def _make_id(domain, path):
    return f"{domain}/{path.name}"


@pytest.fixture(
    params=_ALL_SCRIPTS,
    ids=[_make_id(d, p) for d, p in _ALL_SCRIPTS]
)
def pt_script(request):
    return request.param


def test_pt_script(pt_script):
    """Run a PT test script as subprocess; PASS iff exit code == 0."""
    domain, script_path = pt_script
    result = subprocess.run(
        [sys.executable, str(script_path)],
        capture_output=True,
        text=True,
        timeout=300,
        cwd=str(SCRIPT_DIR),
        env={**os.environ, "PYTHONPATH": str(SCRIPT_DIR)},
        encoding="utf-8",
        errors="replace",
    )
    if result.returncode != 0:
        output = (result.stdout + result.stderr)[-500:]
        pytest.fail(
            f"{script_path.name} exited with code {result.returncode}\n{output}"
        )
