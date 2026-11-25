"""How to cite BVID-FE.

The one place in the package that holds the citation: every animation
frame carries it in its footer (``bvidfe.viz.anim``), so a video that is
shared still names the work. ``tests/test_citation.py`` checks the strings
against ``CITATION.cff`` (title, author, year, repository and version), so
the two cannot drift; update both together.

>>> from bvidfe.citation import citation_text
>>> citation_text().startswith("Elhajjar, R. (2026). BVID-FE: ")
True
"""

from __future__ import annotations

from bvidfe import __version__

TITLE = "BVID-FE: Barely Visible Impact Damage residual-strength analysis for composite laminates"
"""Software title, as ``title`` in ``CITATION.cff``."""

AUTHOR_FAMILY = "Elhajjar"
AUTHOR_GIVEN = "Rani"
YEAR = 2026
"""Year of the release (``date-released`` in ``CITATION.cff``)."""

REPOSITORY = "https://github.com/elhajjar1/BVID-FE"
"""``repository-code`` in ``CITATION.cff``."""

VERSION = __version__


def citation_text() -> str:
    """The citation in one line: author, year, title and repository."""
    return f"{AUTHOR_FAMILY}, {AUTHOR_GIVEN[0]}. ({YEAR}). {TITLE}. {REPOSITORY}"


def frame_credit() -> str:
    """The footer line burned into every animation frame."""
    return f"BVID-FE v{VERSION} — {citation_text()} — please cite"


def short_credit() -> str:
    """The footer of a small frame (a GIF preview), where the full line
    would be too small to read: version, author, year and repository."""
    repo = REPOSITORY.removeprefix("https://")
    return (
        f"BVID-FE v{VERSION} · {AUTHOR_FAMILY}, {AUTHOR_GIVEN[0]}. ({YEAR}) "
        f"· {repo} · please cite"
    )
