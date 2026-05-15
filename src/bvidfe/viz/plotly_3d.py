"""Plotly 3D viz helpers for the BVID-FE Streamlit app.

Replaces the desktop-only PyVista plotters with pure-numpy + Plotly
equivalents that render in the browser. Functions are pure: they take
the ``FeMesh`` from ``bvidfe.analysis.fe_mesh`` plus parameters and
return a ``plotly.graph_objects.Figure``.

Two plot families are exposed:

- ``mesh_damage_figure`` — boundary surface of the hex mesh, coloured by
  per-element damage factor (1 - oop_damage). Used for the "3D Damage"
  tab.
- ``stress_field_figure`` — same boundary surface coloured by an
  arbitrary per-element scalar (e.g. damage severity, von Mises). Used
  for the "Stress Field" tab when the FE tier produces field data.

For the through-thickness damage projection rendered in the old desktop
``Mesh3DTab``, see ``mesh_3d_orthographic_figure`` in :mod:`plots_2d`-
style matplotlib — that view is more informative for BVID inspectors
than a rotating hex mesh and is the default in the Streamlit app.
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import plotly.graph_objects as go

from bvidfe.analysis.fe_mesh import FeMesh

# Local node indices for the 6 faces of a linear hex element, using the
# VTK_HEXAHEDRON / CGNS convention that ``bvidfe.analysis.fe_mesh`` emits.
HEX_FACES: np.ndarray = np.array(
    [
        [0, 1, 2, 3],  # z-low
        [4, 5, 6, 7],  # z-high
        [0, 1, 5, 4],  # y-low
        [3, 2, 6, 7],  # y-high
        [0, 3, 7, 4],  # x-low
        [1, 2, 6, 5],  # x-high
    ],
    dtype=np.int64,
)


def _boundary_faces(elements: np.ndarray) -> np.ndarray:
    """Return boundary quad faces with parent element id.

    Output shape is ``(n_boundary_faces, 5)``: cols 0-3 are global node
    indices in the element's local face ordering, col 4 is the parent
    element id (used to look up cell scalars).
    """
    face_owners: dict[tuple[int, ...], list[tuple[int, int]]] = {}
    for ei in range(elements.shape[0]):
        conn = elements[ei]
        for fi, face in enumerate(HEX_FACES):
            key = tuple(sorted(int(v) for v in conn[face]))
            face_owners.setdefault(key, []).append((ei, fi))

    out: list[list[int]] = []
    for owners in face_owners.values():
        if len(owners) == 1:
            ei, fi = owners[0]
            face = HEX_FACES[fi]
            out.append([int(elements[ei][n]) for n in face] + [ei])
    return np.asarray(out, dtype=np.int64)


def _quads_to_triangles(quad_faces: np.ndarray) -> np.ndarray:
    """Split each quad face into two triangles, preserving parent elem id."""
    n = quad_faces.shape[0]
    tri = np.empty((n * 2, 4), dtype=np.int64)
    tri[0::2, 0] = quad_faces[:, 0]
    tri[0::2, 1] = quad_faces[:, 1]
    tri[0::2, 2] = quad_faces[:, 2]
    tri[0::2, 3] = quad_faces[:, 4]
    tri[1::2, 0] = quad_faces[:, 0]
    tri[1::2, 1] = quad_faces[:, 2]
    tri[1::2, 2] = quad_faces[:, 3]
    tri[1::2, 3] = quad_faces[:, 4]
    return tri


def _scene_layout() -> dict:
    return dict(
        xaxis_title="x [mm]",
        yaxis_title="y [mm]",
        zaxis_title="z [mm]",
        aspectmode="data",
    )


def _mesh3d_figure(
    vertices: np.ndarray,
    elements: np.ndarray,
    *,
    cell_scalar: np.ndarray | None = None,
    colorscale: str = "Viridis",
    colorbar_title: str = "",
    title: str = "",
    cmin: Optional[float] = None,
    cmax: Optional[float] = None,
    height: int = 480,
) -> go.Figure:
    bf = _boundary_faces(elements)
    tri = _quads_to_triangles(bf)

    kwargs: dict = dict(
        x=vertices[:, 0],
        y=vertices[:, 1],
        z=vertices[:, 2],
        i=tri[:, 0],
        j=tri[:, 1],
        k=tri[:, 2],
        flatshading=True,
    )
    if cell_scalar is not None:
        intensity = np.asarray(cell_scalar)[tri[:, 3]]
        kwargs.update(
            intensity=intensity,
            intensitymode="cell",
            colorscale=colorscale,
            colorbar=dict(title=colorbar_title, len=0.7),
            showscale=True,
        )
        if cmin is not None:
            kwargs["cmin"] = cmin
        if cmax is not None:
            kwargs["cmax"] = cmax
    else:
        kwargs.update(color="#9ec5fe", showscale=False)

    fig = go.Figure(go.Mesh3d(**kwargs))
    fig.update_layout(
        title=title,
        scene=_scene_layout(),
        margin=dict(l=0, r=0, t=40, b=0),
        height=height,
    )
    return fig


def mesh_damage_figure(mesh: FeMesh, *, title: str = "Damage factor") -> go.Figure:
    """Render the hex mesh boundary surface coloured by damage factor.

    ``damage_factor`` is the per-element out-of-plane stiffness multiplier
    (1.0 for pristine, ``DAMAGE_OOP_FACTOR`` ≈ 0.05 inside delaminations
    and the fiber-break core). Lower values = more severe damage; we
    invert the scale for the colorbar so red = damaged.
    """
    return _mesh3d_figure(
        mesh.node_coords,
        mesh.element_connectivity,
        cell_scalar=1.0 - mesh.damage_factors,
        colorscale="Reds",
        colorbar_title="1 − damage_factor",
        title=title,
        cmin=0.0,
        cmax=1.0,
    )


def stress_field_figure(
    mesh: FeMesh,
    element_scalar: np.ndarray,
    *,
    label: str = "stress",
    title: Optional[str] = None,
    colorscale: str = "Plasma",
) -> go.Figure:
    """Render the hex mesh boundary surface coloured by an element scalar.

    ``element_scalar`` must have length ``mesh.n_elements``. Use this for
    von Mises, max failure index, or any other per-element field produced
    by the fe3d tier.
    """
    if element_scalar.size != mesh.n_elements:
        raise ValueError(
            f"element_scalar must have length n_elements ({mesh.n_elements}); "
            f"got {element_scalar.size}"
        )
    return _mesh3d_figure(
        mesh.node_coords,
        mesh.element_connectivity,
        cell_scalar=np.asarray(element_scalar),
        colorscale=colorscale,
        colorbar_title=label,
        title=title or label,
    )
