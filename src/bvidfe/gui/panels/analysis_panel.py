"""Analysis configuration panel (tier, loading, mesh params, Run button)."""

from __future__ import annotations

from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import (
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QPushButton,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)


class AnalysisPanel(QWidget):
    """Input panel for analysis tier, loading mode, mesh parameters, and run trigger."""

    configChanged = pyqtSignal()
    runRequested = pyqtSignal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)

        self.tier_combo = QComboBox()
        self.tier_combo.addItems(["empirical", "semi_analytical", "fe3d"])

        self.loading_combo = QComboBox()
        self.loading_combo.addItems(["compression", "tension"])

        # Mesh parameters group (only relevant for fe3d tier)
        self.mesh_group = QGroupBox("Mesh parameters (fe3d only)")
        mg = QFormLayout(self.mesh_group)

        self.elements_per_ply = QSpinBox()
        self.elements_per_ply.setRange(1, 20)
        self.elements_per_ply.setValue(4)

        self.in_plane_size = QDoubleSpinBox()
        self.in_plane_size.setRange(0.1, 100.0)
        self.in_plane_size.setValue(1.0)
        self.in_plane_size.setDecimals(2)

        mg.addRow("Elements per ply:", self.elements_per_ply)
        mg.addRow("In-plane size (mm):", self.in_plane_size)

        self.run_button = QPushButton("Run analysis")

        # Signal connections
        self.tier_combo.currentIndexChanged.connect(self._on_tier_changed)
        self.loading_combo.currentIndexChanged.connect(lambda _: self.configChanged.emit())
        self.elements_per_ply.valueChanged.connect(lambda _: self.configChanged.emit())
        self.in_plane_size.valueChanged.connect(lambda _: self.configChanged.emit())
        self.run_button.clicked.connect(self.runRequested.emit)

        v = QVBoxLayout(self)
        form = QFormLayout()
        form.addRow("Tier:", self.tier_combo)
        form.addRow("Loading:", self.loading_combo)
        v.addLayout(form)
        v.addWidget(self.mesh_group)
        v.addWidget(self.run_button)

        self._on_tier_changed()

    def _on_tier_changed(self, *_: object) -> None:
        self.mesh_group.setVisible(self.tier_combo.currentText() == "fe3d")
        self.configChanged.emit()

    def get_tier(self) -> str:
        return self.tier_combo.currentText()

    def get_loading(self) -> str:
        return self.loading_combo.currentText()

    def get_elements_per_ply(self) -> int:
        return int(self.elements_per_ply.value())

    def get_in_plane_size_mm(self) -> float:
        return float(self.in_plane_size.value())
