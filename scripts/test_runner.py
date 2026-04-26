"""
test_runner.py -- pytest entry point for the PT companion suite.

Each chapter script (proof_*.py, test_*.py in chapter subdirectories) is
executed as a standalone subprocess; PASS iff the script returns exit
code 0. The list of scripts is auto-discovered in conftest.py.

Run:
    pytest                   # all scripts
    pytest -k "ch10"         # one chapter
    pytest --co              # list discovered cases (no run)
"""

import os
import subprocess
import sys

import pytest


def test_pt_script(pt_script):
    """Run a PT test script as subprocess; PASS iff exit code == 0."""
    domain, script_path = pt_script
    result = subprocess.run(
        [sys.executable, str(script_path)],
        capture_output=True,
        text=True,
        timeout=300,
        cwd=str(script_path.parent.parent),
        env={**os.environ, "PYTHONPATH": str(script_path.parent.parent)},
        encoding="utf-8",
        errors="replace",
    )
    if result.returncode != 0:
        output = (result.stdout + result.stderr)[-500:]
        pytest.fail(
            f"{script_path.name} exited with code {result.returncode}\n{output}"
        )
