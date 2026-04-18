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
        text = results.summary()
        # Append a user-facing note when the fe3d tier is used explaining why
        # knockdown may not track impact energy the way empirical does.
        if results.tier_used == "fe3d" and results.knockdown > 0:
            text += (
                "\n\nNote: the fe3d tier's residual-strength prediction is "
                "controlled by stress concentration at the healthy/damaged "
                "element boundary and in this v0.2.0-dev release does not "
                "strongly scale with impact energy once damage exceeds the "
                "Olsson threshold. For energy-dependent knockdown curves "
                "switch the tier to 'empirical' (Soutis scales with DPA) or "
                "'semi_analytical' (Rayleigh-Ritz sublaminate buckling "
                "scales with the largest ellipse)."
            )
        self.text_area.setPlainText(text)
