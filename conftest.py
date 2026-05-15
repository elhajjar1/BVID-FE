"""Top-level pytest config.

Forces a headless matplotlib backend so plot-rendering tests work on CI
runners without a display server, and so importing modules that touch
``matplotlib.pyplot`` (e.g. ``bvidfe.viz.plots_2d``) doesn't try to open
a Tk window.
"""

import matplotlib

matplotlib.use("Agg")
