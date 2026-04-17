"""Input mode selector panel (impact-driven vs. damage-driven)."""

from __future__ import annotations

from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import QButtonGroup, QHBoxLayout, QRadioButton, QWidget


class InputModePanel(QWidget):
    """Radio button panel to choose between impact-driven and damage-driven workflows."""

    configChanged = pyqtSignal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.impact_radio = QRadioButton("Impact-driven")
        self.damage_radio = QRadioButton("Damage-driven (inspection)")
        self.impact_radio.setChecked(True)

        self.group = QButtonGroup(self)
        self.group.addButton(self.impact_radio)
        self.group.addButton(self.damage_radio)

        for r in (self.impact_radio, self.damage_radio):
            r.toggled.connect(lambda _: self.configChanged.emit())

        lay = QHBoxLayout(self)
        lay.addWidget(self.impact_radio)
        lay.addWidget(self.damage_radio)

    def current_mode(self) -> str:
        """Return 'impact' or 'damage'."""
        return "impact" if self.impact_radio.isChecked() else "damage"
