"""Material selection + layup + ply-thickness input panel."""

from __future__ import annotations

from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import (
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QLineEdit,
    QWidget,
)

from bvidfe.core.material import MATERIAL_LIBRARY


class MaterialPanel(QWidget):
    """Select material preset + enter layup + ply thickness."""

    configChanged = pyqtSignal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.material_combo = QComboBox()
        # Add IM7/8552 first so it sits at index 0 (default selection);
        # then append the remaining presets in library order.
        preferred = "IM7/8552"
        names = list(MATERIAL_LIBRARY.keys())
        ordered = [preferred] + [n for n in names if n != preferred]
        for name in ordered:
            self.material_combo.addItem(name)
        self.material_combo.setCurrentIndex(0)  # IM7/8552
        self.material_combo.currentIndexChanged.connect(self._on_changed)

        self.layup_edit = QLineEdit("0, 45, -45, 90, 90, -45, 45, 0")
        self.layup_edit.editingFinished.connect(self._on_changed)

        self.thickness_spin = QDoubleSpinBox()
        self.thickness_spin.setRange(0.001, 10.0)
        self.thickness_spin.setDecimals(4)
        self.thickness_spin.setSingleStep(0.001)
        self.thickness_spin.setValue(0.152)
        self.thickness_spin.valueChanged.connect(self._on_changed)

        form = QFormLayout(self)
        form.addRow("Material preset:", self.material_combo)
        form.addRow("Layup (deg):", self.layup_edit)
        form.addRow("Ply thickness (mm):", self.thickness_spin)

    def _on_changed(self, *_: object) -> None:
        self.configChanged.emit()

    # --- API for the main window ---

    def get_material_name(self) -> str:
        return self.material_combo.currentText()

    def get_layup_deg(self) -> list[float]:
        try:
            return [float(x.strip()) for x in self.layup_edit.text().split(",")]
        except ValueError:
            return []

    def get_ply_thickness_mm(self) -> float:
        return float(self.thickness_spin.value())
