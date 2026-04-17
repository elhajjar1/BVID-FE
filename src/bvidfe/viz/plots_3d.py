"""3D PyVista visualizations for BVID-FE.

Exposes:
- mesh_to_pyvista: convert FeMesh to pyvista.UnstructuredGrid with cell data
- plot_mesh_with_damage: highlight damaged elements by stiffness factor
- plot_mode_shape: color nodes by a scalar mode-shape field
- plot_stress_field: color elements by a scalar stress measure
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import pyvista as pv

from bvidfe.analysis.fe_mesh import FeMesh

VTK_HEXAHEDRON = 12  # VTK cell type code for linear hex


def mesh_to_pyvista(mesh: FeMesh) -> pv.UnstructuredGrid:
    """Convert an FeMesh to a pyvista.UnstructuredGrid with ply + damage cell data."""
    n_elem = mesh.n_elements
    conn = mesh.element_connectivity
    # pyvista expects a flat "cells" array: [n_pts_per_cell, p0, p1, ..., n_pts_per_cell, ...]
    cells = np.empty((n_elem, 9), dtype=int)
    cells[:, 0] = 8
    cells[:, 1:] = conn
    cell_types = np.full(n_elem, VTK_HEXAHEDRON, dtype=np.uint8)
    grid = pv.UnstructuredGrid(cells.flatten(), cell_types, mesh.node_coords.copy())
    grid.cell_data["damage_factor"] = mesh.damage_factors.copy()
    grid.cell_data["ply_index"] = mesh.ply_indices.copy()
    return grid


def plot_mesh_with_damage(
    mesh: FeMesh,
    title: Optional[str] = None,
    show_edges: bool = True,
) -> pv.Plotter:
    """Return a pyvista Plotter rendering the mesh colored by damage factor."""
    grid = mesh_to_pyvista(mesh)
    p = pv.Plotter(off_screen=True)
    p.add_mesh(
        grid,
        scalars="damage_factor",
        cmap="RdYlGn",
        clim=(0.0, 1.0),
        show_edges=show_edges,
        scalar_bar_args={"title": "Damage factor"},
    )
    p.add_title(title or "BVID mesh with damage factor")
    return p


def plot_mode_shape(
    mesh: FeMesh,
    mode_vector: np.ndarray,
    title: Optional[str] = None,
) -> pv.Plotter:
    """Plot a nodal scalar field (e.g. eigenvector component) on the mesh."""
    grid = mesh_to_pyvista(mesh)
    if mode_vector.shape[0] == mesh.n_dof:
        # Interpret as (n_nodes, 3) displacement; take the x-component magnitude
        disp = mode_vector.reshape(mesh.n_nodes, 3)
        scalars = np.linalg.norm(disp, axis=1)
    else:
        scalars = np.asarray(mode_vector).flatten()
        if scalars.size != mesh.n_nodes:
            raise ValueError(f"mode_vector must have length n_dof or n_nodes; got {scalars.size}")
    grid.point_data["mode_scalar"] = scalars
    p = pv.Plotter(off_screen=True)
    p.add_mesh(
        grid,
        scalars="mode_scalar",
        cmap="viridis",
        show_edges=False,
        scalar_bar_args={"title": "|mode|"},
    )
    p.add_title(title or "Mode shape")
    return p


def plot_stress_field(
    mesh: FeMesh,
    element_stress_scalar: np.ndarray,
    title: Optional[str] = None,
    label: str = "stress [MPa]",
) -> pv.Plotter:
    """Plot an element-valued scalar (e.g., von Mises, max Tsai-Wu index) on the mesh."""
    grid = mesh_to_pyvista(mesh)
    if element_stress_scalar.size != mesh.n_elements:
        raise ValueError(
            f"element_stress_scalar must have length n_elements ({mesh.n_elements});"
            f" got {element_stress_scalar.size}"
        )
    grid.cell_data[label] = np.asarray(element_stress_scalar).copy()
    p = pv.Plotter(off_screen=True)
    p.add_mesh(
        grid,
        scalars=label,
        cmap="plasma",
        show_edges=False,
        scalar_bar_args={"title": label},
    )
    p.add_title(title or label)
    return p
