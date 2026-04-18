"""BVID-FE desktop app entry point."""

from __future__ import annotations

import logging
import sys
import traceback

from PyQt6.QtWidgets import QApplication

from bvidfe.gui.main_window import BvidMainWindow


def _install_excepthook() -> None:
    """Install a top-level exception hook that logs + suppresses abort.

    Background: PyQt6 on Qt 6.5+ will call ``qFatal()`` (which aborts the
    process) on any unhandled Python exception raised inside a Qt slot.
    That behavior is surprising for users — a minor matplotlib or plotting
    hiccup takes down the whole app with a macOS crash report.

    Our fix: replace ``sys.excepthook`` so the traceback goes to stderr
    (and the bvidfe.gui logger) but the process keeps running. Individual
    tab-update slots in ``_on_analysis_ready`` are already wrapped in
    try/except for the same reason; this excepthook is the belt-and-
    suspenders backstop for anything we didn't wrap.
    """
    log = logging.getLogger("bvidfe.gui")

    def _hook(exc_type, exc_value, exc_tb):
        # Route KeyboardInterrupt normally
        if issubclass(exc_type, KeyboardInterrupt):
            sys.__excepthook__(exc_type, exc_value, exc_tb)
            return
        log.error(
            "Unhandled exception in GUI: %s: %s",
            exc_type.__name__,
            exc_value,
        )
        traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stderr)

    sys.excepthook = _hook


def main() -> int:
    _install_excepthook()
    app = QApplication(sys.argv)
    window = BvidMainWindow()
    window.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
