"""Test CLI --version flag exists and reports the package version."""

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
