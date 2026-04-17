"""Publication-quality matplotlib style constants for BVID-FE plots."""

from __future__ import annotations

FIGSIZE_LANDSCAPE = (8.0, 5.0)
FIGSIZE_SQUARE = (6.0, 6.0)

DPI_PUBLICATION = 300
DPI_SCREEN = 100

FONT = {
    "family": "sans-serif",
    "size": 11,
}

COLORS = {
    "empirical": "#1f77b4",
    "semi_analytical": "#ff7f0e",
    "fe3d": "#2ca02c",
    "pristine": "#7f7f7f",
    "knockdown": "#d62728",
}

ELLIPSE_CMAP = "viridis"  # colour by interface_index
