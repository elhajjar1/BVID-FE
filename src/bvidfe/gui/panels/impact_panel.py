"""Impact event input panel (energy, impactor, mass, location, live E_onset display)."""

from __future__ import annotations

from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import (
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QLabel,
    QWidget,
)

from bvidfe.core.geometry import ImpactorGeometry
from bvidfe.impact.mapping import ImpactEvent


class ImpactPanel(QWidget):
    """Input panel for a single ImpactEvent."""

    configChanged = pyqtSignal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)

        self.energy_spin = QDoubleSpinBox()
        self.energy_spin.setRange(0.01, 1000.0)
        self.energy_spin.setValue(30.0)
        self.energy_spin.setDecimals(2)

        self.diameter_spin = QDoubleSpinBox()
        self.diameter_spin.setRange(1.0, 100.0)
        self.diameter_spin.setValue(16.0)
        self.diameter_spin.setDecimals(2)

        self.shape_combo = QComboBox()
        self.shape_combo.addItems(["hemispherical", "flat", "conical"])

        self.mass_spin = QDoubleSpinBox()
        self.mass_spin.setRange(0.01, 1000.0)
        self.mass_spin.setValue(5.5)
        self.mass_spin.setDecimals(3)

        self.location_x = QDoubleSpinBox()
        self.location_x.setRange(-1e4, 1e4)
        self.location_x.setValue(0.0)

        self.location_y = QDoubleSpinBox()
        self.location_y.setRange(-1e4, 1e4)
        self.location_y.setValue(0.0)

        self.onset_label = QLabel("E_onset: \u2014 J")

        for w in (
            self.energy_spin,
            self.diameter_spin,
            self.mass_spin,
            self.location_x,
            self.location_y,
        ):
            w.valueChanged.connect(self._on_changed)
        self.shape_combo.currentIndexChanged.connect(self._on_changed)

        form = QFormLayout(self)
        form.addRow("Energy (J):", self.energy_spin)
        form.addRow("Impactor diameter (mm):", self.diameter_spin)
        form.addRow("Impactor shape:", self.shape_combo)
        form.addRow("Impactor mass (kg):", self.mass_spin)
        form.addRow("Location X (mm):", self.location_x)
        form.addRow("Location Y (mm):", self.location_y)
        form.addRow("", self.onset_label)

    def _on_changed(self, *_: object) -> None:
        self.configChanged.emit()

    def get_impact_event(self) -> ImpactEvent:
        """Return the current ImpactEvent from widget values."""
        return ImpactEvent(
            energy_J=float(self.energy_spin.value()),
            impactor=ImpactorGeometry(
                diameter_mm=float(self.diameter_spin.value()),
                shape=self.shape_combo.currentText(),
            ),
            mass_kg=float(self.mass_spin.value()),
            location_xy_mm=(
                float(self.location_x.value()),
                float(self.location_y.value()),
            ),
        )

    def set_onset_energy(self, E_onset_J: float) -> None:
        """Update the live E_onset display label."""
        self.onset_label.setText(f"E_onset: {E_onset_J:.2f} J")
