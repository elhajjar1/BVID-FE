"""Panel geometry input panel (Lx, Ly, boundary condition)."""

from __future__ import annotations

from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import QComboBox, QDoubleSpinBox, QFormLayout, QWidget


class PanelPanel(QWidget):
    """Input panel for PanelGeometry (dimensions + boundary condition)."""

    configChanged = pyqtSignal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.lx_spin = QDoubleSpinBox()
        self.lx_spin.setRange(1.0, 10000.0)
        self.lx_spin.setValue(150.0)
        self.lx_spin.setDecimals(1)

        self.ly_spin = QDoubleSpinBox()
        self.ly_spin.setRange(1.0, 10000.0)
        self.ly_spin.setValue(100.0)
        self.ly_spin.setDecimals(1)

        self.boundary_combo = QComboBox()
        self.boundary_combo.addItems(["simply_supported", "clamped", "free"])

        for w in (self.lx_spin, self.ly_spin):
            w.valueChanged.connect(self._on_changed)
        self.boundary_combo.currentIndexChanged.connect(self._on_changed)

        form = QFormLayout(self)
        form.addRow("Lx (mm):", self.lx_spin)
        form.addRow("Ly (mm):", self.ly_spin)
        form.addRow("Boundary:", self.boundary_combo)

    def _on_changed(self, *_: object) -> None:
        self.configChanged.emit()

    def get_Lx_mm(self) -> float:
        return float(self.lx_spin.value())

    def get_Ly_mm(self) -> float:
        return float(self.ly_spin.value())

    def get_boundary(self) -> str:
        return self.boundary_combo.currentText()
