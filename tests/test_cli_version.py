"""Tests for CLI version + quick flags."""

import subprocess


def test_cli_version_flag_prints_version():
    res = subprocess.run(
        ["./.venv/bin/bvidfe", "--version"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert res.returncode == 0
    assert res.stdout.strip().startswith("bvidfe ")
    # Version is consumed from bvidfe.__version__
    import bvidfe

    assert bvidfe.__version__ in res.stdout


def test_cli_quick_flag_prints_only_knockdown():
    """--quick prints just the knockdown as a scalar, no JSON."""
    res = subprocess.run(
        [
            "./.venv/bin/bvidfe",
            "--material",
            "IM7/8552",
            "--layup",
            "0,45,-45,90,90,-45,45,0",
            "--thickness",
            "0.152",
            "--panel",
            "150x100",
            "--loading",
            "compression",
            "--tier",
            "empirical",
            "--energy",
            "20",
            "--quick",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert res.returncode == 0
    # stdout should be a single float (to 6 decimals), no JSON braces
    stdout = res.stdout.strip()
    assert "{" not in stdout
    assert "}" not in stdout
    # Must parse as a float in (0, 1]
    kd = float(stdout)
    assert 0.0 < kd <= 1.0
