"""Style constants shared by the BVID-FE plots and animations.

The light-surface tokens and the categorical order are the validated palette
of the predicted-against-measured plots (``validation/plot_validation.py``
imports them from here), so a figure, a validation plot and an animation
frame use the same inks and series colours.
"""

from __future__ import annotations

from typing import Tuple, Union

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

# ---- Light-surface tokens -------------------------------------------------
SURFACE = "#fcfcfb"  # chart and frame background
TEXT = "#0b0b0b"  # primary ink: titles, values
TEXT_2 = "#52514e"  # secondary ink: axis labels, ticks, captions
MUTED = "#898781"  # muted ink: hairline annotations
GRID = "#e4e3df"  # hairline gridlines
BAND = "#efeeea"  # shaded bands (the validation +-15% band)
BASELINE = "#c3c2b7"  # axis baselines and reference rules

# ---- Categorical series colours --------------------------------------------
# The eight-slot categorical order validated for colour-vision deficiency on
# the light surface (adjacent pairs). With every pair on screen at once
# (scatter, overlapping shapes) only the first three slots separate, so a
# fourth series and beyond carries a second cue: a line style (see
# ``series_linestyle``), or a hatch on a fill. Never a generated ninth hue.
CATEGORICAL: Tuple[str, ...] = (
    "#2a78d6",  # 1 blue
    "#eb6834",  # 2 orange
    "#1baf7a",  # 3 aqua
    "#eda100",  # 4 yellow
    "#e87ba4",  # 5 magenta
    "#008300",  # 6 green
    "#6250d6",  # 7 violet
    "#e34948",  # 8 red
)

# Status colours, reserved for state (always with an icon and a label).
STATUS = {
    "good": "#0ca30c",
    "warning": "#fab219",
    "serious": "#ec835a",
    "critical": "#d03b3b",
}

LineStyle = Union[str, Tuple[float, Tuple[float, ...]]]

# Dash patterns (in line widths) for the series past the third. Five
# patterns against eight colours: (colour, pattern) repeats only every 40
# series.
_DASHES: Tuple[Tuple[float, ...], ...] = (
    (5.0, 2.0),
    (1.2, 1.6),
    (6.0, 2.0, 1.2, 2.0),
    (2.8, 1.6),
    (6.0, 1.6, 1.2, 1.6, 1.2, 1.6),
)


def series_color(k: int) -> str:
    """Colour of categorical series ``k`` (0-based): slot ``k mod 8``."""
    return CATEGORICAL[k % len(CATEGORICAL)]


def series_linestyle(k: int) -> LineStyle:
    """Line style of categorical series ``k``: solid for the first three,
    then a dash pattern, so no two of the first 40 series share both colour
    and line style."""
    if k < 3:
        return "solid"
    return (0.0, _DASHES[(k - 3) % len(_DASHES)])


# ---- Sequential ramp ----------------------------------------------------------
# One hue, light to dark, for a magnitude (the load ratio to local failure
# of the fe3d load ramp): the blue ramp of the validated palette, steps 100
# to 700. The lightest step means "near zero" and recedes toward the surface.
SEQUENTIAL: Tuple[str, ...] = (
    "#cde2fb",  # 100
    "#b7d3f6",  # 150
    "#9ec5f4",  # 200
    "#86b6ef",  # 250
    "#6da7ec",  # 300
    "#5598e7",  # 350
    "#3987e5",  # 400
    "#2a78d6",  # 450
    "#256abf",  # 500
    "#1c5cab",  # 550
    "#184f95",  # 600
    "#104281",  # 650
    "#0d366b",  # 700
)


def sequential_cmap(name: str = "bvidfe_sequential"):
    """The sequential ramp as a matplotlib colormap (light at 0, dark at 1)."""
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list(name, list(SEQUENTIAL))
