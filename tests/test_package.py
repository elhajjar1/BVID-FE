"""Package import and version-consistency checks.

The version string lives in four places: ``pyproject.toml`` (what PyPI and
``pip show`` report), ``bvidfe.__version__`` (what ``bvidfe --version``
prints; kept as a literal because the Streamlit Cloud deployment runs from
a source checkout without installing the package, so package metadata is
unavailable there), ``CITATION.cff`` and the README BibTeX entry. v0.2.0
shipped with all four disagreeing; this test keeps them in lockstep.
"""

import re
from pathlib import Path

import bvidfe

_ROOT = Path(__file__).resolve().parents[1]


def _read(name: str) -> str:
    return (_ROOT / name).read_text(encoding="utf-8")


def _pyproject_version() -> str:
    match = re.search(r'^version\s*=\s*"([^"]+)"', _read("pyproject.toml"), re.MULTILINE)
    assert match, "no [project] version line in pyproject.toml"
    return match.group(1)


def test_bvidfe_importable():
    assert isinstance(bvidfe.__version__, str) and bvidfe.__version__


def test_dunder_version_matches_pyproject():
    assert bvidfe.__version__ == _pyproject_version()


def test_citation_cff_version_matches_pyproject():
    match = re.search(r"^version:\s*(\S+)", _read("CITATION.cff"), re.MULTILINE)
    assert match, "no version field in CITATION.cff"
    assert match.group(1).strip("\"'") == _pyproject_version()


def test_readme_bibtex_version_matches_pyproject():
    match = re.search(r"^\s*version\s*=\s*\{([^}]+)\}", _read("README.md"), re.MULTILINE)
    assert match, "no version field in the README BibTeX entry"
    assert match.group(1) == _pyproject_version()
