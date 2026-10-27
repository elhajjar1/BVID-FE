"""Orthotropic composite material model + built-in material library."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np


@dataclass
class OrthotropicMaterial:
    """Orthotropic composite ply material with strengths, toughnesses, and
    impact-mapping calibration constants.

    Units: moduli and strengths in MPa (N/mm^2), fracture toughness in N/mm,
    density in kg/mm^3 (e.g. 1.58e-6 for CFRP, matching 1.58 g/cm^3).

    Through-thickness strengths (Zt, Zc) and the 1-3 shear strength (S13)
    are optional. For unidirectional CFRP plies under the transverse-isotropy
    assumption, leaving these as ``None`` falls back to the in-plane
    transverse values via the ``Zt_resolved`` / ``Zc_resolved`` /
    ``S13_resolved`` properties (Zt → Yt, Zc → Yc, S13 → S12). Override only
    when material-specific through-thickness data are available — measured
    Zt is typically 20-40% lower than Yt for void-prone CFRP (Makeev et al.
    2012, *Composites Part A*).

    Attributes
    ----------
    name : str
        Label (``MATERIAL_LIBRARY`` key for the presets).
    E11, E22 : float
        Fibre-direction and transverse Young's moduli (MPa); E33 = E22.
    nu12 : float
        Major Poisson's ratio, ``-1 < nu12 < 0.5``.
    G12, G13, G23 : float
        Shear moduli (MPa).
    Xt, Xc, Yt, Yc : float
        Fibre-direction and transverse tensile / compressive strengths
        (MPa, positive magnitudes).
    S12, S23 : float
        In-plane and transverse (2-3) shear strengths (MPa).
    G_Ic, G_IIc : float
        Mode I and mode II interlaminar fracture toughness (N/mm). G_IIc
        sets the impact onset load and the damage area; both set the
        delamination-growth toughness (Benzeggagh-Kenane) of the
        semi-analytical and fe3d CAI; G_Ic sets the dent depth.
    rho : float
        Density (kg/mm^3), ``1e-7 <= rho <= 1e-5``; used for the
        impactor-to-plate mass ratio.
    olsson_alpha : float
        Fraction of the impact energy above onset that delaminates:
        ``DPA = olsson_alpha * (E - E_onset) / G_IIc``. Default 0.0865
        (fitted to Lovejoy & Scotti's IM7/8552 C-scan areas); 0.0648 for
        T800S/3900-2B.
    dent_beta, dent_gamma : float
        Dent depth ``d = h * dent_beta * ((E - E_onset) / (G_Ic h^2)) **
        dent_gamma``, capped at ``h / 2`` (h the laminate thickness).
        Defaults 0.05 and 0.5. The dent does not enter any strength.
    fiber_break_eta, fiber_break_E_threshold : float
        Fibre-break core radius ``r = fiber_break_eta * sqrt(E -
        fiber_break_E_threshold)`` (mm, E in J), used only by fe3d. Off by
        default (``eta = 0``, threshold 1e9 J) and in every preset.
    soutis_k_s, soutis_m : float
        Soutis CAI knockdown ``1 / (1 + k_s (DPA / A_ref) ** m)``, see
        :func:`bvidfe.failure.soutis_openhole.soutis_cai`. Defaults 4.8
        (calibrated on the NCAMP AS4/8552 and Lovejoy & Scotti IM7/8552 CAI
        data) and 0.5.
    wn_d0_mm : float
        Whitney-Nuismer characteristic distance (mm) of the TAI model;
        default 1.0.
    tai_notch_fraction : float
        TAI equivalent-hole diameter as a fraction of the delamination
        footprint width across the load, ``0 < f <= 1``, see
        :func:`bvidfe.failure.soutis_openhole.fibre_damage_tai`. Default
        0.05, fitted (with ``wn_d0_mm`` = 1.0) to one measured TAI
        condition (Körbelin et al. 2020, M21/T800S).
    nu23 : float, optional
        Transverse (2-3) Poisson's ratio of the 3D compliance; ``None``
        resolves (``nu23_resolved``) to ``E22 / (2 G23) - 1``, the value G23
        implies for a transversely isotropic ply.
    Zt, Zc, S13 : float, optional
        Through-thickness strengths and 1-3 shear strength (MPa), see above;
        used by Tsai-Wu only.
    puck_p_nt_plus, puck_p_nt_minus : float
        Puck inclination parameters (``>= 0``), defaults 0.30 / 0.25.

    Raises ``ValueError`` for a non-positive modulus, strength, toughness or
    density, ``nu12`` outside ``(-1, 0.5)``, ``rho`` outside
    ``[1e-7, 1e-5]``, a non-positive ``Zt`` / ``Zc`` / ``S13``, a
    negative Puck parameter, or ``tai_notch_fraction`` outside ``(0, 1]``.
    """

    name: str
    E11: float
    E22: float
    nu12: float
    G12: float
    G13: float
    G23: float
    Xt: float
    Xc: float
    Yt: float
    Yc: float
    S12: float
    S23: float
    G_Ic: float
    G_IIc: float
    rho: float
    # Impact-mapping calibration (defaults chosen for typical CFRP; per-material override).
    # olsson_alpha: fraction of the impact energy above onset that goes into
    # mode-II delamination, DPA = alpha (E - E_onset) / G_IIc. Fitted to the
    # C-scan areas of Lovejoy & Scotti's eight IM7/8552 CAI coupons (median),
    # with the corrected Olsson threshold load (0.078 while the onset energy
    # was ~10x too low).
    olsson_alpha: float = 0.0865
    dent_beta: float = 0.05
    dent_gamma: float = 0.5
    fiber_break_eta: float = 0.0
    fiber_break_E_threshold: float = 1e9  # disabled unless explicitly overridden
    # Calibrated on the measured CAI datasets in validation/datasets (NCAMP
    # AS4/8552, Lovejoy & Scotti IM7/8552), with the damage area scaled by the
    # fixed ASTM D7137 coupon area and predicted by the energy balance: 4.8 is
    # the pooled optimum of both tiers (4.7 with the onset energy ~10x too
    # low). Fitted on NCAMP alone (4.8) it predicts Lovejoy & Scotti at 8.7%;
    # on Lovejoy & Scotti alone (5.2), NCAMP at 4.9%.
    soutis_k_s: float = 4.8
    soutis_m: float = 0.5
    wn_d0_mm: float = 1.0
    # TAI equivalent hole: diameter = tai_notch_fraction x the delamination
    # footprint width across the load. Delaminations cut no fibres, so the
    # hole stands for the fibre damage under the impactor, a small part of
    # the footprint. Fitted with wn_d0_mm = 1.0 to Korbelin, Dreiner &
    # Fiedler (Compos Part C 3, 2020, 100074): five M21/T800S [45/0/-45/90]2s
    # coupons impacted at 25 J keep 0.889 of their tension strength (0.046
    # fits the mean; rounded to 0.05). One condition, in-sample: see
    # validation/datasets/korbelin_m21_t800s_tai.json.
    tai_notch_fraction: float = 0.05
    # Through-thickness strengths and 1-3 shear (None = transverse-isotropy
    # fallback to Yt / Yc / S12 via the *_resolved properties below).
    Zt: Optional[float] = None
    Zc: Optional[float] = None
    S13: Optional[float] = None
    # Puck (1998) action-plane inclination parameters. Tension (``puck_p_nt_plus``)
    # and compression (``puck_p_nt_minus``) slopes of the master fracture
    # envelope at sigma_n = 0; conservative CFRP/epoxy defaults from
    # Puck-Schürmann 1998 (Composites Sci. Tech. 58, 1045-1067).
    puck_p_nt_plus: float = 0.30
    puck_p_nt_minus: float = 0.25
    # Transverse (2-3 plane) Poisson ratio of the 3D compliance. None = the
    # transversely isotropic value E22 / (2 G23) - 1 implied by G23 (see
    # ``nu23_resolved``); set it when nu23 is measured independently of G23.
    nu23: Optional[float] = None

    def __post_init__(self) -> None:
        positive_fields = (
            "E11",
            "E22",
            "G12",
            "G13",
            "G23",
            "Xt",
            "Xc",
            "Yt",
            "Yc",
            "S12",
            "S23",
            "G_Ic",
            "G_IIc",
            "rho",
        )
        for k in positive_fields:
            v = getattr(self, k)
            if v <= 0:
                raise ValueError(f"{k} must be > 0 (got {v})")
        for k in ("Zt", "Zc", "S13"):
            v = getattr(self, k)
            if v is not None and v <= 0:
                raise ValueError(f"{k} must be > 0 if provided (got {v})")
        for k in ("puck_p_nt_plus", "puck_p_nt_minus"):
            v = getattr(self, k)
            if v < 0:
                raise ValueError(f"{k} must be >= 0 (got {v})")
        if not 0.0 < self.tai_notch_fraction <= 1.0:
            raise ValueError(
                f"tai_notch_fraction must be in (0, 1] (got {self.tai_notch_fraction})"
            )
        if not -1.0 < self.nu12 < 0.5:
            raise ValueError(f"nu12 out of physical range (got {self.nu12})")
        nu23 = self.nu23_resolved
        source = "nu23" if self.nu23 is not None else "nu23 = E22 / (2 G23) - 1"
        if not -1.0 < nu23 < 1.0:
            raise ValueError(
                f"{source} = {nu23:.3g} is outside (-1, 1); "
                "set nu23 explicitly or check E22 and G23"
            )
        # Positive-definite 3D compliance (Lempriere 1968): with E33 = E22 and
        # nu13 = nu12, the normal block is positive definite iff
        # (1 + nu23) (1 - nu23 - 2 nu12 nu21) > 0.
        if 1.0 - nu23 - 2.0 * self.nu12 * self.nu21 <= 0.0:
            raise ValueError(
                f"{source} = {nu23:.3g} with nu12 = {self.nu12} gives a 3D compliance "
                "that is not positive definite (need nu23 < 1 - 2 nu12 nu21)"
            )
        # rho is kg/mm^3; realistic engineering materials span ~[1e-7, 1e-5]
        # (0.1 .. 10 g/cm^3). A value outside this range almost always means
        # the user followed the old (wrong) "t/mm^3" docstring and is 1000x
        # off, which silently corrupts the impact-mapping mass ratio.
        if not 1e-7 <= self.rho <= 1e-5:
            raise ValueError(
                f"rho={self.rho} kg/mm^3 is outside the realistic range "
                f"[1e-7, 1e-5] (0.1-10 g/cm^3). Note rho is kg/mm^3, "
                f"e.g. 1.58e-6 for CFRP."
            )

    @property
    def nu21(self) -> float:
        return self.nu12 * self.E22 / self.E11

    @property
    def nu23_resolved(self) -> float:
        """Transverse Poisson ratio nu23 of the 3D compliance. Falls back to
        the transversely isotropic value ``E22 / (2 G23) - 1`` when not set."""
        if self.nu23 is not None:
            return self.nu23
        return self.E22 / (2.0 * self.G23) - 1.0

    @property
    def Zt_resolved(self) -> float:
        """Through-thickness tensile strength; falls back to Yt under
        transverse isotropy when not explicitly set."""
        return self.Yt if self.Zt is None else self.Zt

    @property
    def Zc_resolved(self) -> float:
        """Through-thickness compressive strength; falls back to Yc under
        transverse isotropy when not explicitly set."""
        return self.Yc if self.Zc is None else self.Zc

    @property
    def S13_resolved(self) -> float:
        """1-3 plane shear strength; falls back to S12 under transverse
        isotropy when not explicitly set."""
        return self.S12 if self.S13 is None else self.S13

    def get_compliance_matrix(self) -> np.ndarray:
        """Return the 6x6 orthotropic compliance matrix in Voigt notation
        (engineering shear strains), material frame.

        Transversely isotropic in the 2-3 plane for the moduli: E33 = E22,
        nu13 = nu12. The 2-3 coupling uses ``nu23_resolved`` (about 0.3-0.7
        for CFRP), not nu21 = nu12 E22 / E11 (about 0.02), which is the
        minor 1-2 ratio. Plane-stress CLT does not use this matrix.
        """
        S = np.zeros((6, 6))
        S[0, 0] = 1.0 / self.E11
        S[1, 1] = 1.0 / self.E22
        S[2, 2] = 1.0 / self.E22
        S[0, 1] = S[1, 0] = -self.nu12 / self.E11
        S[0, 2] = S[2, 0] = -self.nu12 / self.E11
        S[1, 2] = S[2, 1] = -self.nu23_resolved / self.E22
        S[3, 3] = 1.0 / self.G23
        S[4, 4] = 1.0 / self.G13
        S[5, 5] = 1.0 / self.G12
        return S

    def get_stiffness_matrix(self) -> np.ndarray:
        """Return the 6x6 orthotropic stiffness matrix (inverse of compliance)."""
        return np.linalg.inv(self.get_compliance_matrix())


#: Built-in presets by name: ``"AS4/3501-6"``, ``"AS4/8552"``, ``"IM7/8552"``,
#: ``"T700/2510"``, ``"T800S/3900-2B"`` and ``"T800/epoxy"``. Pass a name as
#: ``AnalysisConfig.material``; to change a value, copy a preset with
#: ``dataclasses.replace(MATERIAL_LIBRARY["IM7/8552"], name="...", ...)``
#: rather than mutating the shared instance. AS4/8552, IM7/8552 and
#: T800S/3900-2B cite their sources below; the other three are generic values.
MATERIAL_LIBRARY: dict[str, OrthotropicMaterial] = {
    "AS4/3501-6": OrthotropicMaterial(
        name="AS4/3501-6",
        E11=138000,
        E22=9000,
        nu12=0.30,
        G12=6900,
        G13=6900,
        G23=3450,
        Xt=2280,
        Xc=1440,
        Yt=57,
        Yc=228,
        S12=71,
        S23=50,
        G_Ic=0.26,
        G_IIc=1.0,
        rho=1.58e-6,
    ),
    # Hexcel 8552 AS4 unidirectional tape. E11, E22, nu12, G12, strengths and
    # density: NCAMP CAM-RP-2010-002 Rev A (2011), Table 2-1, RTD means
    # (fiber-dominated values normalized to CPT = 0.0074 in = 0.188 mm); S12 is
    # F12s at 5% strain. Xt is NCAMP's F1tu derived from the [0/90] cross-ply
    # tension tests with the back-out factor (271.47 ksi), the 0 deg ply
    # strength as realised inside a laminate; the unidirectional-coupon value
    # (289.47 ksi = 1996 MPa) overpredicted NCAMP's unnotched multidirectional
    # tension strengths. G23 = E22 / (2 (1 + nu23)) with nu23 = 0.45 measured on
    # AS4/8552 by Gonzalez et al. (Compos Sci Technol 71, 2011). S23 is the
    # LaRC03 transverse shear strength from Yc at a 53 deg fracture plane
    # (Davila, Camanho & Rose, J Compos Mater 39, 2005). G_Ic and G_IIc are the
    # IM7/8552 values (same 8552 matrix), as Gonzalez et al. (Compos Struct 94,
    # 2012) use for AS4/8552.
    "AS4/8552": OrthotropicMaterial(
        name="AS4/8552",
        E11=127300,
        E22=9240,
        nu12=0.30,
        G12=4830,
        G13=4830,
        G23=3186,
        Xt=1872,
        Xc=1398,
        Yt=64,
        Yc=268,
        S12=92,
        S23=101,
        G_Ic=0.28,
        G_IIc=0.79,
        rho=1.59e-6,
    ),
    # Xt is NCAMP's F1tu derived from the [0/90] cross-ply tension tests with
    # the back-out factor: CAM-RP-2009-015 Rev B (2019), Table 2-1, RTD,
    # normalized (324.62 ksi), as for AS4/8552. The other constants are generic
    # IM7/8552 values, not NCAMP's.
    "IM7/8552": OrthotropicMaterial(
        name="IM7/8552",
        E11=165000,
        E22=8400,
        nu12=0.34,
        G12=5600,
        G13=5600,
        G23=2800,
        Xt=2238,
        Xc=1590,
        Yt=73,
        Yc=185,
        S12=90,
        S23=55,
        G_Ic=0.28,
        G_IIc=0.79,
        rho=1.57e-6,
    ),
    "T700/2510": OrthotropicMaterial(
        name="T700/2510",
        E11=127000,
        E22=8400,
        nu12=0.31,
        G12=4200,
        G13=4200,
        G23=2500,
        Xt=2100,
        Xc=1200,
        Yt=58,
        Yc=175,
        S12=66,
        S23=45,
        G_Ic=0.22,
        G_IIc=0.60,
        rho=1.55e-6,
    ),
    # Toray T800S/3900-2B unidirectional tape (interlayer-toughened). E11, E22,
    # nu12, G12, G_Ic and G_IIc: Hasebe et al., Data in Brief 60 (2025) 111509,
    # Table 1. Strengths and density: Toray 3900 prepreg data sheet
    # (Rev. 7.30.2020), P2362W-19L T800S-24K UD, RTA. Hasebe's Table 1 gives
    # Xc = 719.4 MPa, under half the data sheet's ASTM D6641 1779 MPa, with no
    # test method; the data sheet value is used. G23 and S23 as for AS4/8552
    # (nu23 = 0.45; LaRC03 from Yc). olsson_alpha is fitted to the projected
    # delamination areas Hasebe et al. measured on their quasi-isotropic
    # laminates struck by the 12.7 mm hemisphere (median of 38 specimens);
    # the default 0.0865 overpredicts them about 1.3x (interlayer-toughened).
    # Over the 329 specimens in validation/datasets with a measured area it
    # predicts the area at a median 1.01x (0.78-1.13x per layup).
    "T800S/3900-2B": OrthotropicMaterial(
        name="T800S/3900-2B",
        E11=153000,
        E22=8000,
        nu12=0.34,
        G12=4030,
        G13=4030,
        G23=2759,
        Xt=3006,
        Xc=1779,
        Yt=60.8,
        Yc=216,
        S12=68.9,
        S23=81,
        G_Ic=0.54,
        G_IIc=1.64,
        rho=1.54e-6,
        olsson_alpha=0.0648,
    ),
    "T800/epoxy": OrthotropicMaterial(
        name="T800/epoxy",
        E11=155000,
        E22=8500,
        nu12=0.33,
        G12=5000,
        G13=5000,
        G23=2600,
        Xt=2700,
        Xc=1680,
        Yt=68,
        Yc=200,
        S12=85,
        S23=52,
        G_Ic=0.25,
        G_IIc=0.70,
        rho=1.58e-6,
    ),
}
