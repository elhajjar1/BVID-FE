"""Tests for BVID-FE CLI entry point."""

import json
import subprocess
import sys


def _run_cli(*args):
    # Invoke the CLI via ``python -m bvidfe.cli`` so the same command works
    # whether the package is installed in .venv/ (local dev) or on the system
    # PATH (CI runners). Hardcoding "./.venv/bin/bvidfe" breaks in CI.
    return subprocess.run(
        [sys.executable, "-m", "bvidfe.cli", *args],
        capture_output=True,
        text=True,
        check=False,
    )


def test_cli_runs_empirical_impact():
    res = _run_cli(
        "--tier",
        "empirical",
        "--material",
        "IM7/8552",
        "--layup",
        "0,45,-45,90,0,45,-45,90",
        "--thickness",
        "0.152",
        "--panel",
        "150x100",
        "--loading",
        "compression",
        "--energy",
        "25",
    )
    assert res.returncode == 0, res.stderr
    data = json.loads(res.stdout)
    assert "knockdown" in data
    assert 0 < data["knockdown"] <= 1.0


def test_cli_runs_tension():
    res = _run_cli(
        "--tier",
        "empirical",
        "--material",
        "IM7/8552",
        "--layup",
        "0,45,-45,90,0,45,-45,90",
        "--thickness",
        "0.152",
        "--panel",
        "150x100",
        "--loading",
        "tension",
        "--energy",
        "25",
    )
    assert res.returncode == 0, res.stderr
    data = json.loads(res.stdout)
    assert data["tier_used"] == "empirical"


def test_cli_help_works():
    res = _run_cli("--help")
    assert res.returncode == 0
    assert "BVID" in res.stdout or "bvid" in res.stdout.lower()


def test_cli_rejects_bad_panel_format():
    res = _run_cli(
        "--tier",
        "empirical",
        "--material",
        "IM7/8552",
        "--layup",
        "0,90,0,90",
        "--thickness",
        "0.2",
        "--panel",
        "notaxspec",
        "--loading",
        "compression",
        "--energy",
        "10",
    )
    assert res.returncode != 0
