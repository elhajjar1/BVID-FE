"""The citation every animation frame carries (bvidfe.citation) matches
CITATION.cff: title, author, year, repository and version."""

import re
from pathlib import Path

from bvidfe import __version__, citation

_CFF = (Path(__file__).resolve().parents[1] / "CITATION.cff").read_text(encoding="utf-8")


def _field(name: str, text: str = _CFF) -> str:
    match = re.search(rf"^\s*-?\s*{re.escape(name)}:\s*(.+?)\s*$", text, re.MULTILINE)
    assert match, f"no {name} in CITATION.cff"
    return match.group(1).strip("\"'")


def _first_author() -> str:
    block = _CFF.split("authors:", 1)[1].split("\nkeywords:", 1)[0]
    return block


def test_title_matches():
    assert citation.TITLE == _field("title")


def test_author_matches_the_first_author():
    authors = _first_author()
    assert citation.AUTHOR_FAMILY == _field("family-names", authors)
    assert citation.AUTHOR_GIVEN == _field("given-names", authors)


def test_repository_matches():
    assert citation.REPOSITORY == _field("repository-code")


def test_version_matches():
    assert citation.VERSION == _field("version") == __version__


def test_year_is_the_release_year():
    assert str(citation.YEAR) == _field("date-released")[:4]


def test_frame_credit_carries_every_part():
    credit = citation.frame_credit()
    for part in (
        f"v{__version__}",
        "Elhajjar, R.",
        f"({citation.YEAR})",
        citation.TITLE,
        citation.REPOSITORY,
        "please cite",
    ):
        assert part in credit


def test_short_credit_carries_author_year_repository():
    short = citation.short_credit()
    for part in (
        f"v{__version__}",
        "Elhajjar, R.",
        f"({citation.YEAR})",
        "github.com/elhajjar1/BVID-FE",
        "please cite",
    ):
        assert part in short
