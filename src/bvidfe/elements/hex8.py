"""8-node isoparametric hexahedral element (C3D8) for 3D composite analysis.

Standard trilinear brick with 2x2x2 Gauss quadrature.

Node ordering (Abaqus/VTK convention):
    Bottom face (zeta = -1): 0, 1, 2, 3 CCW from +z
    Top face    (zeta = +1): 4, 5, 6, 7 CCW
Natural coordinates xi, eta, zeta in [-1, 1].
"""

from __future__ import annotations

import numpy as np

from bvidfe.core.material import OrthotropicMaterial
from bvidfe.elements.gauss import gauss_points_hex

# Node natural coordinates (xi, eta, zeta)
_NODE_COORDS = np.array(
    [
        [-1, -1, -1],
        [+1, -1, -1],
        [+1, +1, -1],
        [-1, +1, -1],
        [-1, -1, +1],
        [+1, -1, +1],
        [+1, +1, +1],
        [-1, +1, +1],
    ],
    dtype=float,
)


def _T_sigma_z(theta_rad: float) -> np.ndarray:
    """Voigt stress transformation matrix for rotation about the z-axis."""
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    T = np.array(
        [
            [c * c, s * s, 0, 0, 0, 2 * s * c],
            [s * s, c * c, 0, 0, 0, -2 * s * c],
            [0, 0, 1, 0, 0, 0],
            [0, 0, 0, c, -s, 0],
            [0, 0, 0, s, c, 0],
            [-s * c, s * c, 0, 0, 0, c * c - s * s],
        ],
        dtype=float,
    )
    return T


class Hex8Element:
    """8-node isoparametric hex element with orthotropic material and optional ply rotation."""

    def __init__(
        self,
        node_coords: np.ndarray,
        material: OrthotropicMaterial,
        ply_angle_deg: float = 0.0,
    ) -> None:
        node_coords = np.asarray(node_coords, dtype=float)
        if node_coords.shape != (8, 3):
            raise ValueError(f"node_coords must be (8,3), got {node_coords.shape}")
        self.node_coords = node_coords
        self.material = material
        self.ply_angle_deg = ply_angle_deg
        self._C_global = self._compute_global_stiffness()

    def _compute_global_stiffness(self) -> np.ndarray:
        C_mat = self.material.get_stiffness_matrix()
        theta = np.radians(self.ply_angle_deg)
        if abs(theta) < 1e-14:
            return C_mat
        T = _T_sigma_z(theta)
        return T @ C_mat @ T.T

    # --- Shape functions and derivatives ---

    def shape_functions(self, xi: float, eta: float, zeta: float) -> np.ndarray:
        """Eight trilinear shape functions at (xi, eta, zeta)."""
        N = np.empty(8)
        for i in range(8):
            xi_i, eta_i, zeta_i = _NODE_COORDS[i]
            N[i] = 0.125 * (1 + xi * xi_i) * (1 + eta * eta_i) * (1 + zeta * zeta_i)
        return N

    def shape_derivatives(self, xi: float, eta: float, zeta: float) -> np.ndarray:
        """Shape function derivatives d N_i / d {xi, eta, zeta}. Returns (3, 8)."""
        dN = np.empty((3, 8))
        for i in range(8):
            xi_i, eta_i, zeta_i = _NODE_COORDS[i]
            dN[0, i] = 0.125 * xi_i * (1 + eta * eta_i) * (1 + zeta * zeta_i)
            dN[1, i] = 0.125 * (1 + xi * xi_i) * eta_i * (1 + zeta * zeta_i)
            dN[2, i] = 0.125 * (1 + xi * xi_i) * (1 + eta * eta_i) * zeta_i
        return dN

    def jacobian(self, xi: float, eta: float, zeta: float) -> np.ndarray:
        """3x3 Jacobian matrix: J_ij = sum_k dN_k/d(xi_i) * x_k_j."""
        dN = self.shape_derivatives(xi, eta, zeta)
        return dN @ self.node_coords  # (3,8) @ (8,3) = (3,3)

    def B_matrix(self, xi: float, eta: float, zeta: float) -> tuple[np.ndarray, float]:
        """Strain-displacement matrix B (6, 24) and det(J) at (xi, eta, zeta).

        Voigt strain = [e_xx, e_yy, e_zz, 2*e_yz, 2*e_xz, 2*e_xy] (engineering shear).
        """
        dN_nat = self.shape_derivatives(xi, eta, zeta)  # (3, 8)
        J = self.jacobian(xi, eta, zeta)  # (3, 3)
        detJ = np.linalg.det(J)
        J_inv = np.linalg.inv(J)
        dN_phys = J_inv @ dN_nat  # (3, 8) — d N_k / d x, d y, d z
        B = np.zeros((6, 24))
        for k in range(8):
            Nx, Ny, Nz = dN_phys[0, k], dN_phys[1, k], dN_phys[2, k]
            col = 3 * k
            # e_xx
            B[0, col + 0] = Nx
            # e_yy
            B[1, col + 1] = Ny
            # e_zz
            B[2, col + 2] = Nz
            # 2*e_yz
            B[3, col + 1] = Nz
            B[3, col + 2] = Ny
            # 2*e_xz
            B[4, col + 0] = Nz
            B[4, col + 2] = Nx
            # 2*e_xy
            B[5, col + 0] = Ny
            B[5, col + 1] = Nx
        return B, detJ

    def stiffness_matrix(self) -> np.ndarray:
        """Element stiffness 24x24 via 2x2x2 Gauss quadrature."""
        gp, wt = gauss_points_hex(order=2)
        K = np.zeros((24, 24))
        C = self._C_global
        for ig in range(gp.shape[0]):
            xi, eta, zeta = gp[ig]
            B, detJ = self.B_matrix(xi, eta, zeta)
            K += np.dot(B.T, np.dot(C, B)) * (detJ * wt[ig])
        return K

    def stress_at_gauss_points(self, u_elem: np.ndarray) -> np.ndarray:
        """Recover Voigt stress (n_gp, 6) at Gauss points from element DOF vector (24,)."""
        gp, _ = gauss_points_hex(order=2)
        out = np.empty((gp.shape[0], 6))
        C = self._C_global
        for ig in range(gp.shape[0]):
            xi, eta, zeta = gp[ig]
            B, _ = self.B_matrix(xi, eta, zeta)
            eps = B @ u_elem
            out[ig] = C @ eps
        return out

    def strain_at_gauss_points(self, u_elem: np.ndarray) -> np.ndarray:
        """Recover Voigt strain (n_gp, 6) at Gauss points from element DOF vector (24,)."""
        gp, _ = gauss_points_hex(order=2)
        out = np.empty((gp.shape[0], 6))
        for ig in range(gp.shape[0]):
            xi, eta, zeta = gp[ig]
            B, _ = self.B_matrix(xi, eta, zeta)
            out[ig] = B @ u_elem
        return out
