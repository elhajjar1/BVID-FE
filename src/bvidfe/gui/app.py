"""BVID-FE desktop app entry point."""

from __future__ import annotations

import sys

from PyQt6.QtWidgets import QApplication

from bvidfe.gui.main_window import BvidMainWindow


def main() -> int:
    app = QApplication(sys.argv)
    window = BvidMainWindow()
    window.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
