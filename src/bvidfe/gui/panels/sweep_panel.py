"""Parametric energy sweep panel."""

from __future__ import annotations

from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import (
    QFormLayout,
    QLineEdit,
    QPushButton,
    QVBoxLayout,
    QWidget,
)


class SweepPanel(QWidget):
    """Input panel for the parametric energy sweep."""

    configChanged = pyqtSignal()
    sweepRequested = pyqtSignal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)

        self.energies_edit = QLineEdit("5, 10, 20, 30, 40")
        self.csv_edit = QLineEdit("")
        self.csv_edit.setPlaceholderText("optional CSV output path")

        self.sweep_button = QPushButton("Run energy sweep")

        self.energies_edit.editingFinished.connect(lambda: self.configChanged.emit())
        self.csv_edit.editingFinished.connect(lambda: self.configChanged.emit())
        self.sweep_button.clicked.connect(self.sweepRequested.emit)

        v = QVBoxLayout(self)
        form = QFormLayout()
        form.addRow("Energies (J):", self.energies_edit)
        form.addRow("CSV output:", self.csv_edit)
        v.addLayout(form)
        v.addWidget(self.sweep_button)

    def get_energies_J(self) -> list[float]:
        """Parse the comma-separated energies field into a list of floats."""
        try:
            return [float(x.strip()) for x in self.energies_edit.text().split(",") if x.strip()]
        except ValueError:
            return []

    def get_csv_path(self) -> str:
        """Return the optional CSV output path (empty string if not set)."""
        return self.csv_edit.text().strip()
