"""Summary tab: text-only display of AnalysisResults."""

from __future__ import annotations

from PyQt6.QtWidgets import QPlainTextEdit, QVBoxLayout, QWidget

from bvidfe.analysis import AnalysisResults


class SummaryTab(QWidget):
    """Text-only summary of an AnalysisResults."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.text_area = QPlainTextEdit(self)
        self.text_area.setReadOnly(True)
        self.text_area.setPlaceholderText("Run an analysis to see its summary here.")
        lay = QVBoxLayout(self)
        lay.addWidget(self.text_area)

    def update(self, results: AnalysisResults) -> None:  # type: ignore[override]
        self.text_area.setPlainText(results.summary())
