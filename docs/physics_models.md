# Physics Models

BVID-FE exposes three residual-strength modeling tiers that share a common
data contract (`DamageState` → tier engine → `AnalysisResults`) but differ in
fidelity, runtime, and the failure mechanisms they capture. This page
summarises what each tier solves, and the impact-to-damage mapping that
feeds them; see the [Python API reference](python_api.md) for usage.

Units throughout: lengths in mm, stresses and moduli in MPa, toughness in
N/mm, impact energy in J (multiplied by 10³ to get N·mm where it meets a
toughness).

## Impact-to-damage mapping

When an `AnalysisConfig` carries an `ImpactEvent` rather than a measured
`DamageState`, `impact_to_damage()` (`src/bvidfe/impact/`) builds the damage
state in five steps. It uses `AnalysisConfig.impact_panel` when the plate
was impacted on a different panel or window than the one tested (for
example the 125 × 75 mm ASTM D7136 support window and a 150 × 100 mm CAI
coupon, or a 60 × 60 mm clamped window and an 80 × 50 mm coupon); the
damage is then moved to the centre of the test panel. The impacted panel
also bounds the delamination footprint (step 5).

1. **Onset energy.** The threshold load in N is the quasi-static,
   large-mass load at which one delamination grows unstably in a centrally
   loaded plate (Davies & Robinson 1992; Davies & Zhang 1995; Suemasu &
   Majima 1996), with Olsson's effective stiffness for orthotropic plates
   (Olsson 1992, 2001):

    ```text
    P_c = π √(32 D* G_IIc / 3)
    D*  = √(D11 D22 (A + 1) / 2),   A = (D12 + 2 D66) / √(D11 D22)
    ```

   (`impact/olsson.py:threshold_load`), with `G_IIc` in N/mm and the CLT
   bending stiffnesses in N·mm. It is the `n` = 1 case of
   `π √(32 D* G_IIc / (n + 2))` for `n` delaminations (González et al.
   2011, Eq. 8; Salvetti et al. 2016, Eq. 4), and for an isotropic plate
   it is Davies' `P_c² = 8 π² E h³ G_IIc / (9 (1 − ν²))`. `D*` is the
   harmonic mean over directions of the directional bending stiffness,
   which sets the static point-load deflection (exact for `D11 = D22`):
   `A` = 1.0–1.2 and `D* ≈ √(D11 D22)` on quasi-isotropic laminates, and
   `D*` = 0.77 √(D11 D22) on a [0/90]2s cross-ply (`A` = 0.17). No
   dynamic factor is applied: the 1.213× of Olsson, Donadon & Falzon
   (2006) is for small-mass, inertia-governed impacts, and ASTM D7136 drop
   weights are large-mass. Up to 0.2.0 the code used
   `π √(8 G_IIc √(D11 D22) / 9)`, Davies' form with `D` written for
   `E h³ / (1 − ν²)` (which is `12 D`), so `P_c` was √12 = 3.5× and the
   onset energy 8–12× too low.

   Against measured threshold loads (the first large load drop of
   instrumented drop-weight tests), predicted / measured:

   | Laminate | Measured | 0.2.0 | Now |
   |---|---|---|---|
   | González et al. 2011, AS4/8552 [(45/0/−45/90)4]s, 5.8 mm, `G_IIc` 0.79 | 9.94 kN | 2.45 kN (0.25×) | 8.58 kN (0.86×) |
   | Salvetti et al. 2016, IM7/MTM45-1 QI, 32 plies, 4.13 mm, `G_IIc` 0.8 assumed | about 5.7 kN | 1.64 kN (0.29×) | 5.67 kN (0.99×) |
   | NIAR (CAM-RP-2013-020), IM7/8552 QI, 4.45 mm (24 plies assumed) | about 6.7 kN | 1.83 kN (0.27×) | 6.45 kN (0.96×) |

   The González row uses the AS4/8552 preset; the other two the paper's
   ply properties (Salvetti) or the IM7/8552 preset (NIAR). Clustered plies
   lower the measured load (González's 2- and 4-ply clusters: 7.6 and
   5.5 kN, predicted 1.11× and 1.51×), which the single-delamination
   formula does not see. `P_c` does not depend on the panel or the
   impactor.

   The onset energy is the energy the plate and the contact store when
   the load reaches `P_c` (`impact/olsson.py:onset_energy`):

    ```text
    P = k_b w + k_m w³        plate: deflection w under the impactor
    P = K_c δⁿ                contact: indentation δ
    E_onset = k_b w_c² / 2 + k_m w_c⁴ / 4 + P_c δ_c / (n + 1)
    ```

   where `w_c` is the real root of `k_b w + k_m w³ = P_c` and `δ_c` the
   indentation at `P_c`. Both energies are the exact integrals `∫ P dw`
   and `∫ P dδ`.
    - `k_b` is the point-load bending stiffness of the simply supported
      panel at the impact site (Navier series, 11 × 11 modes, using D11,
      D12, D22 and D66 only: D16, D26 and the coupling matrix B are
      ignored), times 1.0 (`simply_supported`), 2.5 (`clamped`) or 0.4
      (`free`, with a `UserWarning`). The 2.5 and 0.4 are uncalibrated
      round numbers, not solutions for those edges.
    - `k_m` is the membrane (large-deflection) stiffness: the single-term
      solution of the von Kármán plate equations for the panel deflected
      in an assumed mode, with the stress function solved exactly for
      that mode and the equilibrium equation satisfied in the Galerkin
      sense (Levy 1942 for the isotropic simply supported plate; Chia 1980
      for orthotropic laminates). Shivakumar, Elber & Illg (1985) and
      Olsson use the same `k_b w + k_m w³` law for impacted plates. With
      `a = A⁻¹` the membrane compliance, `α = π/Lx` and `β = π/Ly`:
      - `simply_supported`: mode `sin(αx) sin(βy)`, edges free to move
        in-plane but kept straight (an ASTM D7136 coupon lies on its
        support window under four toggle clamps):
        `k_m = (Lx Ly / 64) (α⁴ / a11 + β⁴ / a22)`. For an isotropic
        square plate of side `L` this is `π⁴ E h / (32 L²)`.
      - `clamped`: mode `sin²(αx) sin²(βy)`, edges held in-plane (gripped
        between plates): a sum over the seven Fourier terms of the
        mode's stress function plus the edge forces that cancel the edge
        shortening, `(9 Lx Ly / 512) (A11 α⁴ + 2 A12 α² β² + A22 β⁴)`
        (the formula is in `olsson.py`).
      - `free`: `k_m` = 0.

      An independent Ritz minimisation of the membrane energy over the
      in-plane displacements gives the same coefficients
      (`tests/impact/test_olsson_membrane.py`). Approximations: one
      assumed mode (a point load deflects the plate more locally, and the
      shape changes as membrane action grows); the membrane part uses the
      mode while the bending part keeps the Navier `k_b`; A16, A26 and B
      are ignored; the central-load value is used at every impact point
      (one mode cannot represent a deflection concentrated near an edge);
      and `P_c` is the small-deflection threshold load.
    - The contact law is Hertz sphere-on-flat,
      `P = (4/3) E* √R δ^1.5`, for `hemispherical` (`R` = diameter / 2),
      a flat punch, `P = 2 E* R δ`, for `flat`, and Love's cone,
      `P = (2 E* tan α / π) δ²` with `α` = 30° from the axis, for
      `conical` (Johnson 1985). `E*` combines a steel impactor (200 GPa,
      ν = 0.3) with the ply's `E22` (ν = 0.3). The contact energy
      `P_c δ_c / (n + 1)` is `2 P_c δ_c / 5` for Hertz, `P_c δ_c / 2` for
      the punch and `P_c δ_c / 3` for the cone.

    Below `E_onset` the damage state is empty. Until the membrane term was
    added the plate was linear, `E_onset = P_c² / (2 k_cb)` with
    `k_cb = 1 / (1/k_b + 1/k_c)` and a contact stiffness linearised at
    `P_c` (the Hertz secant, which stores `P_c δ / 2` instead of
    `2 P_c δ / 5`, and the cone's tangent, `P_c δ / 4` instead of
    `P_c δ / 3`). A thin or large plate deflects many times its thickness
    before `P_c` on the linear law, where membrane stretching carries most
    of the load, so the linear onset grew with the panel area. IM7/8552,
    16 mm hemisphere, simply supported unless noted (`w_c / h` is the
    deflection at `P_c` over the thickness, linear law in brackets):

    | Laminate, panel | `w_c / h` | `E_onset` linear | Now |
    |---|---|---|---|
    | 8-ply QI [0/45/−45/90]s, 150 × 100 mm | 2.6 (13.9) | 7.2 J | 0.84 J |
    | 8-ply [0/90]2s, 150 × 100 mm | 2.4 (15.8) | 8.5 J | 0.79 J |
    | 16-ply QI [0/45/−45/90]4, 150 × 100 mm | 1.2 (2.0) | 6.6 J | 3.3 J |
    | 24-ply QI [0/45/−45/90]3s, 150 × 100 mm | 0.64 (0.77) | 7.4 J | 5.7 J |
    | 16-ply QI, 150 × 100 mm, clamped | 0.66 (0.80) | 2.8 J | 2.1 J |
    | 8-ply QI, 300 × 200 mm | 4.4 (56) | 29 J | 1.2 J |
    | 16-ply QI, 500 × 500 mm | 4.4 (38) | 122 J | 8.1 J |
    | NCAMP AS4/8552 coupons, 125 × 75 mm window (median) | | 4.65 J | 4.27 J |
    | Lovejoy & Scotti IM7/8552 coupons, same window (median) | | 7.2 J | 6.5 J |
    | Hasebe et al. T800S/3900-2B, 60 × 60 mm clamped window (median) | | 2.5 J | 2.3 J |

    The onset now grows with thickness (`P_c² ∝ h³` while the membrane
    stiffness grows only as `h`) and only slowly with the panel size. The
    validation coupons are thick for their windows (`w_c / h` = 0.28 on
    NCAMP's, 0.10–0.12 on Lovejoy & Scotti's, a median 0.16 on Hasebe
    et al.'s), so their onset falls by 5–10%, mostly from the exact Hertz
    energy: of the median drop, the membrane term is 0.21 of 0.38 J
    (NCAMP), 0.03 of 0.67 J (Lovejoy & Scotti) and 0.02 of 0.23 J (Hasebe
    et al.). Up to 0.2.0, with `P_c` 3.5× too low, the onset was
    0.6–0.9 J on 150 × 100 mm panels.

    The cone's contact is the softest, and with the exact contact energy
    its onset is now further above the hemisphere's: 7.0 J against 5.7 J
    before on the 24-ply Hasebe et al. cross-ply specimens, where the
    elastic cone indents 1.7 mm at `P_c` and the contact stores 5.2 J.
    Three of those specimens struck by the 60° cone at 1.1 J/mm (5.0 J)
    show 62–72 mm² of damage where the model predicts none, as it did
    before (the membrane term is negligible there, `w_c / h` = 0.09). The
    data have no specimen near the predicted onset to check it against:
    Hasebe et al.'s lowest level is 1.1 J/mm (24-ply only), and their 14
    undamaged specimens (11 struck by the flat-ended impactor, 3 by the
    25.4 mm hemisphere) were struck at 2.2–5.9× the predicted onset (the
    model predicts damage for all of them, as it did with the linear
    plate).

2. **Damage area.** The target projected damage area is

    ```text
    DPA = α (E − E_onset) · 10³ / G_IIc × (5.5 / m)^0.1 × f_shape × (16 / D)^0.2
    ```

    in mm². It is an energy balance: a fraction `α = olsson_alpha` of the
    impact energy above onset goes into mode-II delamination. The default
    0.0855 is fitted to the median of Lovejoy & Scotti's eight IM7/8552
    C-scan areas. T800S/3900-2B (interlayer-toughened) has 0.0605, fitted
    to Hasebe et al.'s quasi-isotropic laminates struck by the 12.7 mm
    hemisphere; over the 329 of their 343 specimens with a measured area
    the median predicted / measured area is 1.03× (per-layup medians
    0.77–1.26×, the highest on the thinnest cross-ply, c8). Both were
    refitted when the threshold load of step 1 was corrected (they were
    0.078 and 0.047 with an onset energy about 10× too low; most test
    energies are 5–20× the onset, so `α` absorbed most of the error), and
    again when the membrane term and the exact contact energy were added
    to the onset (0.0865 and 0.0648 with the linear plate). The fit and the check share data: only the 291 specimens
    outside the 38 used for `α` are out of sample, and the same specimens
    set the impactor factors below.

    The area used to be divided by the laminate thickness `h` as well. That
    gave mm rather than mm² and left the area flat in thickness at constant
    energy per thickness, where measured areas grow with thickness:
    Hasebe et al.'s 8- to 24-ply laminates at fixed J/mm, and Lovejoy &
    Scotti's 7.6 mm coupons, which had 1.21× the C-scan area of their
    6.1 mm coupons at the same 8.9 J/mm. Thickness now enters only through
    `E_onset`.

    `m` is the impactor mass (kg), `D` its diameter (mm; the diameter
    factor is 1 for a cone, whose contact is set by its angle) and
    `f_shape` = 1.0 (hemispherical), 0.2 (flat) or 1.9 (conical). `f_shape`
    and the exponent 0.2 are fitted to Hasebe et al. (2025), the only open
    data that vary the impactor on one material: over same-layup,
    same-energy groups a 15.9 mm flat-ended impactor left 0.20× and a 60°
    cone 1.20× the C-scan area of a 12.7 mm hemisphere, and 25.4 and
    31.8 mm hemispheres 0.88× and 0.87×. The tip also sets the onset
    energy through the contact energy (the cone's is the highest), so the
    factors are fitted to make the model's predicted / measured area for
    each tip match the 12.7 mm hemisphere's in the same groups: 0.22 → 0.2
    (flat), 1.86 → 1.9 (cone; 1.54 → 1.5 with the linearised contact,
    whose cone stored a quarter of `P_c δ` rather than a third); the
    diameter ratios moved by under 5% with the corrected onset (1.03 and
    1.26 for the 25.4 and 31.8 mm hemispheres now), so the exponent stays
    0.2. The per-impactor median predicted / measured area is then
    0.94–1.07×. At the ASTM D7136 reference
    (5.5 kg, 16 mm, hemispherical) the three correction factors are 1. DPA
    is capped at 80% of the panel area (`DPACapClipWarning`). On a
    150 × 100 mm IM7/8552 panel the area grows by about 110 mm² per joule
    above onset (on 8- to 24-ply quasi-isotropic laminates 0–450 mm² at
    5 J, 470–990 mm² at 10 J and 2,640–3,160 mm² at 30 J, the thinnest
    the largest). The footprint confinement of step 5 then caps it below
    the 12,000 mm² limit, at a value set by the layup: about 11,500 mm²
    from 107–112 J on symmetric quasi-isotropic laminates, 10,410 mm²
    from about 99 J on [0/45/−45/90]4 and 6,490 mm² from about 61 J on
    [0/90]2s. The
    mass factor uses the impacted panel's mass (`rho · Lx · Ly · h`), and
    a `SmallMassQuasiStaticWarning` is raised when the impactor is lighter
    than it.

3. **Dent depth**: `d = β (E − E_onset)^γ · S_shape` (mm, energies in J),
   zero at or below onset and capped at `h / 2`, with `β = dent_beta`
   (0.044 mm/J^0.57), `γ = dent_gamma` (0.57) and `S_shape` 1 for a
   hemispherical impactor, 8.5 for the 60° cone and 0.03 for a flat-ended
   one. A power law of the impact energy is the usual empirical
   indentation law; Caprino & Lopresto (Compos. Sci. Technol. 60, 2000)
   found theirs scarcely affected by fibre type, orientation or matrix,
   and the measured dents here agree: one `(β, γ)` serves the untoughened
   IM7/8552 and the interlayer-toughened T800S/3900-2B, and neither
   thickness nor `G_Ic` is in the law. On a 16-ply (2.4 mm) IM7/8552
   panel at 20 J (3.3 J onset) the dent is 0.22 mm; on Lovejoy & Scotti's
   6.1 and 7.6 mm coupons at 54 and 68 J, 0.40 and 0.46 mm (measured
   0.36–0.41 and 0.43–0.48 mm).
   The constants are a log least-squares fit to the positive measured
   dents of the hemispherical cases (8 Lovejoy & Scotti coupons, 210
   Hasebe et al. specimens, each dataset weighted equally), with `E_onset`
   computed by the code (`validation/README.md`, *Dent depth*); they were
   fitted with the linear-plate onset energy and kept with the
   membrane-stiffened one, which moves the fit to 0.039 mm/J^0.60 but the
   medians below by under 3%. Median
   predicted / measured is 1.01 on Lovejoy & Scotti (0.96–1.10) and 0.94
   on the 308 Hasebe et al. specimens with a positive dent (10–90%:
   0.52–1.44; hemispherical 0.99, cone 0.91). Fitted on Hasebe et al.
   alone, the law predicts Lovejoy & Scotti at 1.38×; fitted on Lovejoy &
   Scotti (with Hasebe's exponent), Hasebe at 0.71×. The scatter within
   Hasebe et al. is a factor of about 1.7 (log-sd 0.54), and the law does
   not carry the impactor diameter: the 12.7, 25.4 and 31.8 mm
   hemispheres are predicted at 0.78, 1.12 and 1.44× their dents. The
   cap binds on seven of the measured cases (the cone on the 1.5 mm 8-ply
   laminates at 5.2 J, up to 1.04× the cap uncapped); 18 of those cone
   dents exceed `h / 2` (up to 1.2 h), so on thin laminates struck by a
   sharp impactor the cap under-predicts. The dent is the
   as-measured depth: dents relax with time, moisture and temperature
   (by up to 3× in some studies), and certification criteria for "barely
   visible" (about 0.3–0.5 mm, often after relaxation) refer to the
   relaxed dent. Before this fit the law was
   `d = h · 0.05 · ((E − E_onset) · 10³ / (G_Ic h²))^0.5` (dimensionally
   consistent, but never calibrated): about 3 mm × √(E − E_onset) for
   8552, on its `h / 2` cap for 345 of the 351 measured cases and
   7.5–8.8× the measured dents (uncapped, a median 33× on Hasebe et al.
   and 51× on Lovejoy & Scotti).

4. **Fibre-break radius**: `r_fb = η √max(0, E − E_fb)` with
   `η = fiber_break_eta` and `E_fb = fiber_break_E_threshold`. `η` is 0 in
   every preset, so the fibre-break core is off unless a custom material
   or a C-scan file sets it.

5. **Distribution over interfaces** ("peanut" template). Ply 0 of
   `layup_deg` (at `z = −h/2`) is the impacted face, and interface `i`
   (`0 … n_plies − 2`) lies between ply `i` (impact side) and ply `i + 1`
   (back-face side). Fibre directions are axial, so every angle is taken
   modulo 180°: `−45` and `135`, or `90` and `−90`, give the same
   footprint. Each interface whose two plies differ in direction gets one
   ellipse centred on the impact point. Its orientation is the fibre
   angle of the back-face-side ply `i + 1`: a peanut-shaped delamination
   elongates along the fibres of the ply below the interface (Liu 1988;
   Abrate 1998). Its aspect ratio is `1 + 0.025 Δθ`, where `Δθ` ∈ [0°, 90°]
   is the acute angle between the two fibre directions (1 to 3.25; a
   −45/90 interface is 45°, AR 2.125). Its size weight
   `0.3 + 0.7 (i + 1) / n_interfaces` grows from the impact face to the
   back face. Interfaces between plies of the same direction (the doubled
   mid-plane ply of a symmetric laminate) carry no delamination: there is
   no bending-stiffness mismatch to drive one (Liu 1988), and
   delaminations are observed only between plies of different orientation
   (Abrate 1991). A layup whose plies all share one direction keeps a
   circle at every interface. One scale factor is solved so the union of
   the footprints equals the target DPA within 1%. The weights take ply 0
   (the first entry of `layup_deg`, at `z = −h/2` in the CLT stack and at
   `z = 0` in the fe3d mesh) as the impact face, so interface 0 gets the
   smallest delamination; list the layup from the impacted face.

   Every ellipse is kept inside the impacted panel (`impact_panel`, else
   `panel`): an ellipse whose bounding box would reach past the panel edge
   nearest the impact point is shrunk, keeping its aspect ratio, and the
   scale factor is re-solved for the same union area. If the panel cannot
   hold the target area, the largest footprint that fits is used with a
   `UserWarning`. The union of differently oriented ellipses spreads
   beyond the circle of the same area: unconfined, the template footprint
   is 67 × 62 mm on the NCAMP quasi-isotropic coupon (same-area circle
   61 mm), 89 × 85 mm on Lovejoy & Scotti's quasi-isotropic coupon
   (circle 82 mm) and 125 × 122 mm on their hard laminate (circle 92 mm),
   wider than the 102 mm coupon. (Before the angles were taken modulo
   180°, the −45/90 interfaces of the quasi-isotropic coupons were
   ellipses of aspect ratio 4 at 22.5°, and the footprints were 99 × 52
   and 134 × 70 mm; the hard laminate's was 108 × 108 mm.) For an ASTM
   D7136 coupon the impacted panel is the 125 × 75 mm support window,
   which is where the damage stays.

The dent depth is carried in `DamageState` and reported, but no tier uses
it in the residual-strength calculation.

## Empirical tier

Both models are closed-form and run in milliseconds. CAI sees the damage
only through its projected area `DPA` (the union of all ellipse
footprints), TAI only through the width of that footprint across the
load.

**CAI** (after Soutis & Curtis 1996):

```text
σ_CAI = σ₀ / (1 + k_s · (DPA / A_ref)^m)
```

`σ₀` is the pristine strength (see
[Knockdown definition](#knockdown-definition-and-cross-tier-comparability)),
`A_ref` = 15,000 mm² is the area of the 150 × 100 mm ASTM D7137 coupon
(`SOUTIS_REFERENCE_AREA_MM2` in `bvidfe.failure.soutis_openhole`), on which
`k_s = soutis_k_s` (default 4.8, calibrated on measured CAI data, see
`validation/README.md`) is calibrated, and `m = soutis_m` (default 0.5).
`A_ref` does not change with the analysed panel. It used to be the panel
area `Lx · Ly`, which sent the knockdown to 1 as the panel grew, although
impact damage is a local notch, and over-penalised Hasebe et al.'s
80 × 50 mm coupons. The panel size still enters the onset energy, the 80%
area cap, the semi-analytical growth check and fe3d.

This is an empirical power-law curve, not Soutis & Curtis's model. They
treat the impact damage as an equivalent open hole and predict its
compressive strength with a fracture-mechanics (damage-zone) model; the
code borrows only the idea that the damage size sets the knockdown. The
form `1 / (1 + k_s (DPA/A_ref)^m)` and both constants are BVID-FE's
own, fitted to the measured datasets.

**TAI**: Whitney-Nuismer point-stress criterion on an equivalent circular
hole that stands for the fibre damage, with a finite-width correction
(`fibre_damage_tai` in `bvidfe.failure.soutis_openhole`):

```text
2R = f · w
σ_TAI / σ₀ = 2 / (2 + ξ² + 3ξ⁴ − (K_T∞ − 3)(5ξ⁶ − 7ξ⁸)) · (K_T∞ / K_T)(2R / W),   ξ = R / (R + d₀)
```

Tension strength after impact is controlled by the broken fibres, which
sit in a small region under the impactor; the delaminations around it cut
no fibres, and the delaminated plies still carry tension. The hole is
therefore a fraction of the damage, not the whole footprint: `w` is the
width of the delamination footprint across the load (the y extent of the
union of the ellipses, at most the panel width `W = Ly`) and
`f = tai_notch_fraction` (default 0.05). A damage state with a fibre-break
core (`fiber_break_radius_mm > 0`) sets a minimum `R`.

`d₀ = wn_d0_mm` (default 1.0 mm) is the characteristic distance, and
`K_T∞ = 1 + √(2(√(Ex/Ey) − ν_xy) + Ex/Gxy)` is Lekhnitskii's
infinite-plate stress concentration for the orthotropic laminate (3.0 for a
quasi-isotropic one). `Ex`, `Ey`, `Gxy` and `ν_xy` come from the
average-stress compliance `a* = h A⁻¹`. When `K_T∞ ≳ 9.2` the polynomial
would give a residual above `σ₀` (or diverge for `K_T∞ ≳ 20.3`), so the
ratio is capped at 1 with a `UserWarning`. The built-in presets stay below
`K_T∞ ≈ 7.8` (the highest is an all-0° T800S/3900-2B laminate; 5.2 for
IM7/8552 [0/90]2s).

`K_T∞ / K_T` is Tan's (1988) finite-width correction for a central hole
in an orthotropic strip, with `x = 2R / W`:

```text
K_T∞ / K_T = 3(1 − x) / (2 + (1 − x)³) + ½ (xM)⁶ (K_T∞ − 3)(1 − (xM)²)
M² = (√(1 − 8(3(1 − x) / (2 + (1 − x)³) − 1)) − 1) / (2x²)
```

For a quasi-isotropic laminate it is the isotropic factor
`3(1 − x) / (2 + (1 − x)³)`, the inverse of `[2 + (1 − 2R/W)³] / [3(1 − 2R/W)]`.
It tends to 1 for a small hole and to 0 as the hole fills the panel. With
`f` = 0.05 the hole stays below 5% of the panel width, so the correction
is below 0.5% unless `f` is raised.

`f` is fitted, with `d₀` = 1 mm, to the one measured TAI condition in
`validation/datasets/`: Körbelin, Dreiner & Fiedler (2020), five
M21/T800S [45/0/−45/90]2s coupons 60 mm wide impacted at 25 J on a
50 × 30 mm window, which keep 0.889 of their tension strength
(0.85–0.97). The predicted footprint fills the 30 mm window, so the
required hole is about 1.4 mm (`f` = 0.046 fits the mean; 0.05 is
adopted, knockdown 0.875). The fit is in-sample, on one energy and one
layup, and only `f · w / d₀` is constrained: the TAI model is not
validated, and `AnalysisResults.notes` says so on every closed-form TAI
run.

On 150 × 100 mm IM7/8552 panels (0.152 mm plies, 16 mm hemisphere,
5.5 kg; below the onset energy, 0.8 J on the 8-ply and 3.2 J on the
16-ply laminates here, both knockdowns are 1):

| Laminate | Loading | 5 J | 10 J | 15 J | 30 J |
|---|---|---|---|---|---|
| 8-ply QI [0/45/−45/90]s | TAI | 0.892 | 0.819 | 0.771 | 0.684 |
| | CAI | 0.546 | 0.448 | 0.395 | 0.312 |
| 16-ply QI [0/45/−45/90]2s | TAI | 0.932 | 0.829 | 0.768 | 0.669 |
| | CAI | 0.646 | 0.485 | 0.417 | 0.321 |
| 8-ply cross-ply [0/90]2s | TAI | 0.886 | 0.835 | 0.802 | 0.733 |
| | CAI | 0.544 | 0.447 | 0.394 | 0.312 |

TAI now falls with energy and stays above CAI, the usual test trend. The
earlier model sized the hole from the whole footprint (`R = √(DPA / π)`,
infinite plate; still available as `whitney_nuismer_tai`). With `R ≫ d₀`
its knockdown sat near the large-hole limit `1 / K_T∞` at every BVID
energy (QI 0.41 at 10 J, 0.36 at 30 J; cross-ply 0.40 and 0.24), below
the CAI knockdown, and it scored 57% against the Körbelin coupons. The
point-stress form also makes the higher-`K_T∞` cross-ply laminate keep
more strength than the quasi-isotropic one from 10 J, at these small
holes: its stress concentration decays faster from the hole edge.

## Semi-analytical tier

For TAI this tier calls the same fibre-damage open-hole model as the
empirical tier, so its tension results are identical. For CAI it adds a
delamination growth bound.

At each delaminated interface the thinner sublaminate (by thickness; on
a tie, plies 0 … i, the impact side) over the largest delamination at
that interface is treated as a plate over the delamination's enclosing
rectangle: the rotated ellipse's bounding box in the panel frame,
`a = 2 · ex` along x (the load direction) and `b = 2 · ey` along y, with
the half-extents `ex = hypot(major_mm · cos θ, minor_mm · sin θ)` and
`ey = hypot(major_mm · sin θ, minor_mm · cos θ)` for θ = `orientation_deg`
(`DelaminationEllipse.half_extents_mm`). These are the extents the
[damage template](#impact-to-damage-mapping) confines and the TAI damage
width reads, so a delamination the template elongates along a ±45° or
90° ply is as wide across the load here as it is for the open-hole
model (it used to be taken as `2 · major_mm` along the load whatever its
orientation).
The sublaminate's `Dij` are taken about its own mid-plane and its
coupling matrix B is ignored, although a sublaminate cut from a
symmetric laminate is in general unsymmetric. The Navier (Rayleigh-Ritz) closed form for a
simply supported orthotropic rectangle under uniaxial compression gives its
buckling load per unit width

```text
N_cr = min over m, n ∈ 1..5 of  (π² / a²) · [D11 m⁴ + 2(D12 + 2D66)(m a/b)² n² + D22 (n a/b)⁴] / m²
```

with the sublaminate's CLT bending stiffnesses `Dij` (D16 and D26
ignored), times a boundary factor of 1.0 (`simply_supported`), 1.9
(`clamped`) or 0.5 (`free`) taken from the **panel's** `boundary`. The
factor transfers the test panel's edge condition to the sublaminate,
whose own edges are the delamination front; it is a ratio of
square-plate buckling coefficients (Timoshenko & Gere 1961), not a
derived sublaminate boundary condition. An aspect ratio `a/b` above 50
is clipped to 50 (`DegenerateThinSublaminateWarning`). The film buckling
stress is `σ_c = N_cr / h`, where `h` is the sublaminate thickness.

Buckling alone does not fail the laminate: the sublaminate keeps
carrying load after it buckles. The CAI strength is the far-field stress
at which the buckled sublaminate grows its delamination. The single-ply
sublaminate over the largest (back-face) delamination does buckle at a
few MPa (3–21 MPa on the panels below), but it is not the one that
grows first: a thin film needs a high film stress to release `G_c`. At the governing interface the film
buckling stress is about 60–210 MPa on 150 × 100 mm IM7/8552 panels at
10–30 J (8-ply [0/45/−45/90]s 60–125 MPa, 16-ply [0/45/−45/90]4
105–210 MPa, 8-ply [0/90]2s 90–190 MPa), 85–87 MPa on the NCAMP coupons
and 110–125 MPa on Lovejoy & Scotti's (the value is in
`AnalysisResults.notes` for fe3d). The thin-film energy release rate
(Chai, Babcock & Knauss 1981; Hutchinson & Suo 1992, straight-sided
blister)

```text
G = h / (2 E_f) · (σ − σ_c)(σ + 3σ_c)
```

is set equal to the mixed-mode toughness, giving the film stress at growth

```text
σ_g = −σ_c + √(4σ_c² + 2 E_f G_c / h)
```

- `σ` is the sublaminate (film) stress along x.
- `G_c` is the Benzeggagh-Kenane mix `G_Ic + (G_IIc − G_Ic)(sin²ψ)^η` at
  the thin-film phase angle ψ = 52.1° with η = 2 (0.478 N/mm for IM7/8552
  and AS4/8552, which share the 8552 toughness values).
- `E_f` is the film's x-stress per unit laminate x-strain, with the
  laminate's own Poisson contraction. Film and laminate share the far-field
  strain, so the laminate stress at growth is `Ex · σ_g / E_f`.
- The weakest interface governs, and the result is capped by the Soutis
  empirical strength: `σ_CAI = min(Soutis, σ_growth)`.

The energy release rate above is Hutchinson & Suo's result for a
straight-sided film with **clamped** edges, where `σ_c` is the clamped
film's own Euler stress. The code instead inserts the `σ_c` of a simply
supported rectangle times the panel boundary factor, so the two halves
of the model assume different edge conditions; and it uses the laminate
`E_f` defined above rather than the film's plane-strain modulus.

The governing interface is `AnalysisResults.critical_sublaminate` and its
buckling load `N_cr` (N/mm) is in `AnalysisResults.buckling_eigenvalues`.
With an empty damage state the tier returns `σ₀`.

With the calibrated `soutis_k_s`, the Soutis bound governs on many
configurations. On 150 × 100 mm IM7/8552 panels (8-ply QI
[0/45/−45/90]s, 16-ply QI [0/45/−45/90]4 as in the synthetic dataset,
and 8-ply cross-ply [0/90]2s; 0.152 mm plies, 10–30 J; the onset energy
is 0.8 J on the 8-ply and 3.3 J on the 16-ply panels) the growth
knockdown is 0.35–0.50 and the Soutis knockdown 0.31–0.49. Growth
governs only on the cross-ply panel at 10–15 J (0.389 and 0.394 against
0.447 and 0.394); Soutis governs elsewhere. Closer to onset, at 5 J,
growth governs all three (8-ply QI 0.469 against 0.546, 16-ply QI 0.557
against 0.651, cross-ply 0.424 against 0.544; the symmetric 16-ply
[0/45/−45/90]2s 0.520 against 0.650).

**The growth channel is effectively unvalidated.** On the NCAMP and
Lovejoy & Scotti validation cases (`validation/datasets/`), the two
gated measured datasets, Soutis governs every coupon (growth stress
205–206 MPa against a 176–177 MPa Soutis bound on NCAMP, 182
against 156 MPa on the Lovejoy & Scotti quasi-isotropic coupons, and 227–228
against 224 MPa on their 40-ply hard laminate), so the two tiers give
the same CAI strengths there and the gate never exercises the growth
stress. Growth governs on 256 of Hasebe et al.'s 343 coupons (154 of the
157 cross-ply coupons, every hard-laminate r0 coupon, and 75 of the 132
quasi-isotropic ones; no soft r45 coupon), and there it makes the prediction worse: on
those 256 the MAE is 25.6% against 23.0% for the empirical tier (on the
other 87 both are 17.2%), so the semi-analytical tier scores 23.4%
against 21.5% over the whole set.

## 3D FE tier

A structured hexahedral mesh of the whole panel is built for the damaged
laminate (`MeshParams`: `in_plane_size_mm`, default 5 mm, and
`elements_per_ply`, default 1). Ply interfaces lie on element faces.
Each ply's 3D stiffness is orthotropic with the 2-3 plane transversely
isotropic in its moduli: `E33 = E22`, `ν13 = ν12`, and the 2-3 Poisson
ratio `ν23` is `OrthotropicMaterial.nu23` or, unset, the value
`E22 / (2 G23) − 1` the preset's `G23` implies (0.30 to 0.68 over the
presets; 0.45 for AS4/8552, the measured value its `G23` was derived
from). Until this release the compliance used `ν21 = ν12 E22 / E11`
(about 0.02) for the 2-3 term. The plane-stress tiers do not use `ν23`.
Delaminated interfaces are approximated by a **component-wise
stiffness-reduction model** (true cohesive surfaces deferred to a future
release):

- Elements that touch a delaminated interface plane (the element layers
  on both sides) and whose centroid lies inside that delamination's ellipse carry
  an out-of-plane factor `DAMAGE_OOP_FACTOR = 0.05`. It scales the
  through-thickness and transverse-shear stiffness and the in-plane /
  out-of-plane Poisson coupling.
- The damage zone as a whole (every ply under the union of the
  delamination footprints, where matrix cracks and fibre damage accompany
  the delaminations; confined to the impacted panel, see
  [step 5](#impact-to-damage-mapping)) is a soft inclusion with in-plane
  stiffness ×`DAMAGE_ZONE_INPLANE_FACTOR` (0.15, fitted to three measured
  CAI coupons, see [Limitations](#limitations); Soutis & Curtis 1996).
  Without it, delaminations left the in-plane stress field unchanged and
  first-ply failure could not see the impact. An element straddling the
  zone's edge takes the factor weighted by the fraction of its area inside
  the zone (sampled on a 6 × 6 grid), so the zone does not move in
  element-sized steps as the mesh changes.
- The zone's in-plane softening depends on the sign of the fibre load.
  In compression (CAI) the factor scales the whole in-plane stiffness:
  delaminated sublaminates buckle, and matrix cracks take away the
  lateral support that keeps the fibres from kinking; the factor is
  fitted to that case. In tension (TAI) the zone's fibres are intact,
  since delaminations and matrix cracks cut none, so only its
  matrix-dominated moduli `E22` and `G12` take the factor and `E11` keeps
  its pristine value. The damage acts on the ply compliance, as in
  continuum damage models with a fibre-tension damage variable set only
  by fibre failure (Matzenmiller, Lubliner & Taylor 1995; Maimí et al.
  2007): in each ply's material frame `S22` and `S66` are divided by the
  in-plane factor, `S33`, `S44` and `S55` by the out-of-plane factor, and
  the Poisson compliances are kept. This is the ply-discount picture
  (`E11` kept, `E22` and `G12` discounted) that the unnotched pristine
  strength already uses. Until this release the zone was softened in
  every in-plane direction in tension too, so it failed as if its fibres
  were cut (see [Limitations](#limitations)).
- The optional fibre-break core (elements whose centre is within
  `fiber_break_radius_mm` of any delamination centroid; off in every
  preset) gets both factors: out-of-plane ×0.05 and in-plane
  ×`DAMAGE_FIBER_BREAK_INPLANE_FACTOR` (0.30). Its broken fibres are not a
  third, fibre-direction softening in tension: they are the hole below,
  whose radius is at least the core's (issue #164). The impact mapping
  sets the core radius to `fiber_break_eta · √(E − fiber_break_E_threshold)`,
  zero for every preset, so only a damage state that gives
  `fiber_break_radius_mm` (a C-scan file, or a custom material card) has
  one. A core that holds no element centre (at 5 mm, one under 3.5 mm:
  the impact point is a node, 3.54 mm from the nearest centre) softens
  nothing and raises `FiberBreakCoreUnmeshedWarning`, recorded in
  `AnalysisResults.warnings`.
- In tension the broken fibres under the impactor are the closed-form
  tiers' equivalent hole (radius `tai_hole_radius_mm`, see the empirical
  TAI model above: the footprint's hole, or the fibre-break core when
  that is larger), with its
  Whitney-Nuismer point-stress knockdown `k_h` (`wn_d0_mm` = 1 mm, Tan's
  finite-width correction). The hole is about 1–4 mm across at BVID
  energies, below the mesh (2–5 mm elements and the 2.5 mm averaging
  disk), so it is not meshed: every element whose footprint overlaps it
  (`fe_mesh.fibre_damage_hole_elements`, all plies) has its disk-averaged
  stress multiplied by `1 / k_h` before the failure check (Tsai-Wu index
  `a c / k_h + b c² / k_h²`; `fe_tier._tai_fibre_damage_notch`). The
  hole thus sits in the strain field of the matrix-cracked zone around it,
  and with the zone removed fe3d returns the closed-form knockdown. The
  two tiers share one fibre-damage size and notch law, and no constant is
  added. Compression has no hole.

The panel is loaded by a prescribed x displacement between rigid,
frictionless platens: `u_x = 0` on every node of the `x = 0` face and
`u_x = ε · Lx` on every node of the `x = Lx` face, with `u_y` and `u_z`
free there. The remaining rigid-body modes are removed by `u_z = 0` on
the mid-surface node layer and `u_y = 0` at one node (on the `x = 0` face,
mid-surface, nearest `y = Ly / 2`). The plies can therefore expand freely
by Poisson in y and through the thickness: an undamaged unidirectional
coupon fails at exactly its ply strength (`Xc / E11`, `Yc / E22`, and in
tension `Xt / E11`, `Yt / E22`), and an undamaged laminate at its CLT
first-ply failure load less the free-edge stresses at `y = 0` and
`y = Ly` (8-ply QI IM7/8552, 150 × 100 mm, 5 mm mesh: 599.3 against
601.3 MPa in compression, 494.5 against 502.5 MPa in tension). For a
laminate and damage symmetric about the mid-plane (the soft zone covers
every ply), the exact solution has `u_z = 0` on the mid-plane, so that
restraint is a symmetry condition and carries no load; an unsymmetric
layup is held flat instead of bending (with an odd number of element
layers the restrained layer is the one half an element below the
mid-plane). The solve is a membrane, pre-buckling state: bending is left
to the closed-form sublaminate growth and panel buckling, which use the
panel `boundary`, and `boundary` does not change these restraints, since
out-of-plane edge supports carry no load in it. Supporting only the
mid-surface lines of the edges and leaving the panel free to bend
changed the knockdowns below by under 1% and made the solves 2–4 times
slower (the ILU-preconditioned Krylov solve no longer converges on a
40-ply coupon), and the bending of fully integrated bricks this thin
would be locked. Until this release `u_z = 0` was held on the whole
`z = 0` face, on the whole loaded faces (`simply_supported`; all four
edge faces for `clamped`) and `u_y = 0` on the whole `y = 0` face; the
edge faces blocked the plies' through-thickness Poisson expansion
(`ν23` = 0.3–0.7), which lowered the undamaged 8-ply QI tension strength
12.6% (439.0 MPa) and made a `[90]_4` coupon fail 3–4% early.
One linear static solve at a reference strain of 5% gives the stress field;
because the response is linear, the load multiplier at which the failure
index first reaches 1 follows analytically (Tsai-Wu: root of a quadratic;
LaRC05: `1/√index`; Puck: `1/index`). If nothing fails up to 5% strain, the
5% state is used.

First-ply failure is evaluated on stress rotated into each ply's material
frame, using LaRC05 (CAI) and Tsai-Wu (TAI). The failure index is averaged
over each element's Gauss points (2 × 2 × 2) and then over a disk of
radius `FPF_AVERAGING_RADIUS_MM` = 2.5 mm in the same through-thickness
layer, each element weighted by the fraction of its area inside the disk
(16 × 16 sample points per element; after Whitney & Nuismer's 1974
average-stress criterion, though here the failure index is averaged,
not the stress), clipped at the panel edges; the averaged index is linear in
the coefficients of the load multiplier, so the analytic scaling above
still holds. A check at single Gauss points has no length scale: the stress
at the re-entrant corners of the zone's footprint is singular, and the
pointwise failure load kept moving as the mesh was refined. Counting
whole elements by their centres instead of weighting them by overlap, the
disk covered only the element and its four neighbours at 2–2.5 mm meshes
and shrank as the mesh was refined. 2.5 mm is within the 1-4 mm
characteristic distances of CFRP notched strength; at the default 5 mm
mesh the disk holds only the element itself, and a damage zone only a few
elements across is under-resolved (see [Limitations](#limitations)). The solver also implements Puck; `BvidAnalysis` does not
expose a criterion choice, and an unknown criterion name raises
`ValueError`.

- Damage-zone elements fail at the strain intact material would (strain
  equivalence; Lemaitre 1992). The criterion sees their effective stress:
  the in-plane stress the pristine ply would carry at the same strain, with
  the actual (damaged) out-of-plane tractions. Their strength therefore
  drops with their stiffness; with full strengths, a zone spanning most of
  the panel lost its stress concentration and made the panel look stronger.
  In tension the effective stress is the continuum-damage one that goes
  with the compliance above, each ply-frame component divided by its
  factor (`σ₂₂ / f`, `τ₁₂ / f`, `σ₁₁` by the fibre factor, 1 outside the
  fibre-break core; the out-of-plane tractions as they are). The zone's
  fibres therefore fail at the intact `Xt`: a `[0]₄` coupon inside the
  zone fails in tension at the load of the intact one, where it carried
  0.15 of it before. Only the matrix-dominated strengths drop. In
  compression the effective stress is still the pristine in-plane
  stiffness times the strain, unchanged.
- The failure stress is the reaction force on the loaded edge at failure
  over the gross section `Ly · h`, so a softened panel is not credited with
  load it cannot carry.

For CAI, the lower of first-ply failure and the semi-analytical delamination
growth stress governs. Buckling onset is reported in
`AnalysisResults.buckling_eigenvalues` (MPa) but does not set the strength.
It is the lower of whole-panel buckling (the Navier closed form above
applied to the full laminate over `Lx × Ly`) and sublaminate buckling at
the interface with the largest delamination area × distance to the farther
laminate surface (`find_critical_interface`), which need not be
`critical_sublaminate`. The closed form replaced a 3D `K_g` eigensolve,
which locked too aggressively on thin laminates at affordable mesh sizes
(issue #129). Whole-panel buckling is left out of the knockdown because CAI
fixtures carry anti-buckling guides; its stress is recorded in
`AnalysisResults.notes`.

The linear solve eliminates the constrained DOFs and runs conjugate
gradients preconditioned by an incomplete LU factorisation to a relative
residual of 1e-10, retrying with GMRES and then a direct sparse LU if it
does not converge (`bvidfe.solver.static`). The fe3d run is refused with
`FE3DSizeError` above 150,000 degrees of freedom (`BVIDFE_FE3D_MAX_DOF`),
which admits an 8-ply 150 × 100 mm panel down to a ~1.7 mm mesh.

### Failure criteria

All three take the ply-frame Voigt stress
`[σ1, σ2, σ3, τ23, τ13, τ12]` (`src/bvidfe/failure/`). Failure is at
index 1.

- **"LaRC05"** (`larc05.py`) is not LaRC05. It is a quadratic,
  Hashin-type index:
  `max((σ1 / X)², (σ2 / Y)² + (τ12 / S12)² + (τ23 / S23)²)`, with
  `X` = Xt or Xc and `Y` = Yt or Yc by the sign of σ1 and σ2. Fibre
  compression is therefore a maximum-stress check on Xc, with no
  fibre-kinking or misalignment analysis; matrix compression uses Yc
  directly, with no fracture-plane search and no friction (normal-stress)
  term; σ3 and τ13 are dropped; and in-situ strengths are not used. It is
  homogeneous of degree 2 in stress.
- **Tsai-Wu** (`tsai_wu.py`, Tsai & Wu 1971) is the full 3D form,
  `Fᵢ σᵢ + Fᵢⱼ σᵢ σⱼ`, with `F11 = 1/(Xt Xc)`, `F1 = 1/Xt − 1/Xc` (and
  likewise for 2 and 3), `F44 = 1/S23²`, `F55 = 1/S13²`, `F66 = 1/S12²`
  and the interaction terms `Fij = −½ √(Fii Fjj)`. Zt, Zc and S13 fall
  back to Yt, Yc and S12 when the material card leaves them unset (every
  preset does).
- **Puck** (`puck.py`) is a simplified Puck & Schürmann (1998)
  inter-fibre-failure check: the action-plane angle is swept in 5° steps
  over ±90°, one pair of inclination parameters (`p+` = 0.30,
  `p−` = 0.25) is used on every plane instead of Puck's interpolation
  between the ⊥∥ and ⊥⊥ values, and `R⊥⊥A = Yc / (2(1 + p−))`. Planes
  with σn ≥ 0 use Puck's mode-A envelope, and every plane with σn < 0
  uses one compressive envelope, with no separate mode-C branch. Fibre
  failure is maximum stress on Xt / Xc. It is homogeneous of degree 1.
  `BvidAnalysis` never selects it.

## Knockdown definition and cross-tier comparability

`AnalysisResults.knockdown` is computed in exactly one place —
`BvidAnalysis.run()` (`src/bvidfe/analysis/bvid.py`):

```python
knockdown = residual_strength_MPa / pristine_strength_MPa
```

**The denominator is identical across all three tiers.**
`_pristine_strength()` is the unnotched laminate strength from classical
lamination theory, by ply-discount last-ply failure
(`bvidfe.failure.laminate_strength`; Daniel & Ishai 2006):

- The laminate is loaded by a uniaxial average stress σₓ.
- A ply whose matrix fails (max stress: σ₂ against Yt/Yc, |τ₁₂| against
  S12) keeps 1% of its E22 and G12, and the load it sheds is redistributed
  at the same applied stress.
- The laminate fails when any ply's fibres fail (σ₁ against Xt/Xc), or when
  every ply has cracked (a [90]ₙ laminate fails at Yt/Yc).

Plies share the laminate strain, so the 0° plies carry most of an axial
load. The thickness-weighted average of the ply strengths used before,
`Σ tᵢ (X·cos²θᵢ + Y·sin²θᵢ) / Σ tᵢ`, credited the off-axis plies with
strength they never reach: it was 1.17–2.85× the NCAMP unnotched strengths
of AS4/8552 and IM7/8552 laminates. Ply-discount matches them within
0.83–1.14 (mean 0.97; `tests/validation/test_unnotched_strength_ncamp.py`):

| Material | Layup | Tension: NCAMP / model (MPa) | Compression: NCAMP / model (MPa) |
| --- | --- | --- | --- |
| AS4/8552 | 25/50/25 | 611 / 699 | 561 / 543 |
| AS4/8552 | 10/80/10 | 439 / 392 | 430 / 357 |
| AS4/8552 | 50/40/10 | 1050 / 1037 | 904 / 830 |
| IM7/8552 | 25/50/25 | 722 / 753 | 600 / 601 |
| IM7/8552 | 10/80/10 | 462 / 415 | 458 / 388 |
| IM7/8552 | 50/40/10 | 1211 / 1237 | 833 / 933 |

Xt for both presets is NCAMP's 0° tension strength derived from the
[0/90] cross-ply tests with the back-out factor, the strength the 0° plies
reach inside a laminate. The unidirectional-coupon value is 6% (AS4/8552)
and 12% (IM7/8552, against the previous generic 2560 MPa) higher, and with
it the 0° plies outlived the ±45 cracking and broke 5–19% above the
measured tension strengths. In the 50/40/10 laminates the 0° plies now
break at Xt after the ±45s crack (0.99–1.02). In the 25/50/25 and 10/80/10
laminates they break as soon as the ±45 plies crack in shear, so S12 sets
the strength: 25/50/25 is 4% (IM7/8552) and 14% (AS4/8552) high, 10/80/10
10% low in tension and 15–17% low in compression. One shear strength
cannot move these in opposite directions.

**The numerator (residual strength) is what differs between tiers:**

| Tier | CAI residual stress | TAI residual stress |
| --- | --- | --- |
| `empirical` | Soutis: `σ₀ / (1 + k_s·(DPA/A_ref)^m)`, `A_ref` = 15,000 mm² | Whitney-Nuismer point stress on a hole of diameter `tai_notch_fraction` × footprint width, finite-width corrected |
| `semi_analytical` | `min(Soutis, σ_growth)` | Delegates to the empirical TAI model (identical to `empirical`) |
| `fe3d` | `σ₀ × min(1, D/U)` with `D`, `U` = `min(σ_growth, FPF_LaRC05)` on the damaged / undamaged panel | `σ₀ × min(1, D/U)` with `D`, `U` = `FPF_Tsai-Wu` damaged (fibre-damage hole as a notch at the impact site) / undamaged |

fe3d's own failure stresses are not the same kind of quantity as σ₀, so
fe3d is normalised by an undamaged fe3d run of the same panel (a second FE
solve; with no delaminations the undamaged growth stress is infinite, so
`U` is first-ply failure). The raw stresses are recorded in
`AnalysisResults.notes`.

**What this means for users:**

- All three tiers report knockdown on the **same scale** (ratio relative to
  the same pristine baseline), so values are *qualitatively comparable*.
- They are **not** numerically interchangeable: each tier captures
  different failure mechanisms.
    - For **TAI**, `empirical` and `semi_analytical` are mathematically
      identical; `fe3d` uses the same fibre-damage hole and adds the
      matrix-cracked zone around it, so it falls with energy like the
      closed form and sits 2–8% below it (8-ply QI 150 × 100 mm IM7/8552
      at 5 / 15 / 30 J: 0.834 / 0.714 / 0.629 against 0.892 / 0.771 /
      0.684; see [Limitations](#limitations)).
    - For **CAI**, `semi_analytical ≤ empirical` always (the growth stress
      only lowers the residual). `fe3d` divides by its own undamaged
      strength (first-ply failure, typically below σ₀), so for the same
      growth stress its knockdown is higher than `semi_analytical`'s.
- For **energy-scaling studies**, prefer `empirical` (Soutis scales with
  DPA). `fe3d` CAI also falls with damage area, through its softened
  damage zone (16-ply QI at 5 / 15 / 30 J: 0.430 / 0.302 / 0.263, see
  [Limitations](#limitations)), but its damage-zone stiffness is fitted to
  only three CAI coupons; fe3d TAI follows the closed-form TAI. Where the delamination growth stress
  governs `semi_analytical`, that tier is flatter in energy: the growth
  stress changes little with delamination size once the sublaminate
  buckles well below it.

## Limitations

- `soutis_k_s` (4.8) is calibrated on two measured CAI datasets (AS4/8552 and
  IM7/8552, both on the 8552 matrix; see `validation/README.md`). Fitted on
  either dataset alone it predicts the other within 5.0–8.8%, so expect
  errors of that order on a new material or layup. The default
  `olsson_alpha` (0.0855) is fitted to Lovejoy & Scotti's C-scan areas, the
  T800S/3900-2B value (0.0605) and the impactor shape and diameter factors
  to Hasebe et al.'s, and `dent_beta` and `dent_gamma` to both datasets'
  measured dents. The other calibration constants (`soutis_m`,
  `wn_d0_mm`, the impactor mass factor, the boundary factors, the mixed-mode phase angle and B-K exponent, and the
  fe3d out-of-plane and fibre-break factors) are uncalibrated defaults
  chosen for typical CFRP systems. Precise values need calibration against
  material-specific coupon test data before use in certification.
- **Material presets.** AS4/8552 and T800S/3900-2B cite their sources in
  `src/bvidfe/core/material.py`; IM7/8552 cites NCAMP for Xt only and
  calls the rest generic values. AS4/3501-6, T700/2510 and T800/epoxy
  have no cited source.
- Validation covers CAI, apart from one TAI condition, and fe3d is not
  run on the measured datasets in CI (it gives 141.0 MPa against the
  measured 177.7 MPa, 20.6% low, on the first NCAMP coupon, one of the
  three coupons its damage-zone factor is fitted to, `validation/README.md`;
  on the measured TAI dataset it is 3.6% off, partly in-sample, see the
  fe3d TAI item below).
- **The TAI model is not validated.** Its one constant,
  `tai_notch_fraction` (0.05), is fitted in-sample to a single measured
  condition: five M21/T800S quasi-isotropic coupons at 25 J (Körbelin et
  al. 2020), run with the T800S/3900-2B card, whose footprint fills the
  impact window. Nothing tests its energy, layup or thickness
  dependence, or the cross-ply result that a higher `K_T∞` keeps more
  strength. The hole stands for fibre damage that BVID-FE does not
  predict: a fixed fraction of the footprint width cannot follow the
  fibre failure that impacts at elevated temperature, sharper impactors
  or energies near penetration cause (Körbelin et al.'s coupons impacted
  at 80 °C lost 0.19 instead of 0.11). With `d₀` = 1 mm the 1–4 mm holes
  are within a few `d₀` of the edge, where the point-stress strength is
  sensitive to `d₀`. See [Empirical tier](#empirical-tier).
- **The CAI delamination-growth channel is effectively unvalidated.**
  Soutis governs every coupon of the two gated measured datasets, and on Hasebe
  et al.'s set, where growth governs 256 of 343 coupons, it raises the
  error (25.6% against 23.0% on those coupons); see
  [Semi-analytical tier](#semi-analytical-tier).
- **The dent depth is an empirical fit with a factor-of-1.7 scatter.**
  `dent_beta` and `dent_gamma` are fitted to two datasets (IM7/8552 and
  T800S/3900-2B; median predicted / measured 1.01 and 0.94), the shape
  factors to Hasebe et al. alone, and the impactor diameter is not in the
  law (step 3 of the [impact mapping](#impact-to-damage-mapping)). It is
  the as-measured dent, before relaxation; no tier uses it.
- **"LaRC05" is a quadratic Hashin-type index**, not LaRC05: no
  fibre-kinking analysis (fibre compression is maximum stress on Xc), no
  matrix fracture-plane search or friction term, no in-situ strengths,
  and σ₃ and τ₁₃ do not enter (see [Failure criteria](#failure-criteria)).
  It is used for fe3d CAI first-ply failure only.
- The `fe3d` tier uses component-wise stiffness reduction (out-of-plane at
  delaminated interfaces, in-plane over the damage zone) instead of true
  cohesive surfaces with bilinear traction-separation laws and a continuum
  damage model. Cohesive surfaces are deferred to a future release. The
  damage-zone in-plane factor is uniform over the footprint. Its value,
  0.15, is fitted in-sample to three measured CAI coupons, with the
  delaminations kept inside the ASTM D7136 impact window: it puts the first
  Lovejoy & Scotti IM7/8552 quasi-isotropic coupon and their first 40-ply
  hard laminate 1.2% and 0.3% high and the first NCAMP AS4/8552 coupon
  20.6% low, the lowest mean error (7.4%) of the factors tried (0.13 /
  0.17 / 0.20 give 13.8% / 10.5% / 12.4%). The NCAMP coupon wants about
  0.22 and the Lovejoy & Scotti coupons about 0.15, so
  no single factor fits both families within 5%; fe3d strengths are
  qualitative. (With the template before it took ply angles modulo 180°,
  step 5, and the ply stiffness before it used `ν23`, 0.17 fitted all
  three within 5%.) Before the footprint was
  confined, the hard laminate's template footprint was 108 × 108 mm on a
  102 mm wide coupon: the soft zone spanned the whole section, the 0°
  fibres failed at the coupon's side edge, and no single factor fitted
  both families (0.20 left the hard laminate 47% low). A damage area at
  the 80%-of-panel cap would put the soft zone across the loaded edges.
- **`fe3d` TAI takes its fibre damage from the closed form, and is
  validated on one condition only.** Tension strength after impact is
  set by the broken fibres under the impactor (Cantwell & Morton 1991);
  delaminations and matrix splits cut no fibres, and in notched tension
  they relieve the stress concentration rather than weaken the laminate
  (Hallett et al. 2009). fe3d's damage zone therefore keeps its fibres
  in tension, and the broken fibres are the closed-form equivalent hole,
  applied as a sub-mesh notch at the impact site (see *3D FE*). Without
  the hole the knockdown was only the strain concentration of the
  matrix-cracked zone in the 90° plies: nearly flat in energy, since that
  depends on the zone's shape rather than its size, and above the closed
  form from about 15 J. On 150 × 100 mm IM7/8552 panels (5 mm mesh) the
  fe3d tension knockdown at 5 / 15 / 30 J is:

    | Laminate | Zone cut in tension | Zone only | Now (zone + hole) | Closed form |
    |---|---|---|---|---|
    | 8-ply QI [0/45/−45/90]s | 0.319 / 0.242 / 0.220 | 0.928 / 0.912 / 0.894 | 0.834 / 0.714 / 0.629 | 0.892 / 0.771 / 0.684 |
    | 16-ply QI [0/45/−45/90]4 | 0.422 / 0.282 / 0.253 | 0.940 / 0.915 / 0.892 | 0.895 / 0.759 / 0.665 | 0.952 / 0.813 / 0.718 |
    | 8-ply cross-ply [0/90]2s | 0.288 / 0.156 / 0.129 | 0.907 / 0.917 / 0.938 | 0.860 / 0.782 / 0.715 | 0.886 / 0.802 / 0.733 |

    fe3d now falls with energy like the closed form, 2–8% below it: at
    the zone's centre the matrix-cracked plies are strained a little more
    than intact ones, and the hole adds its notch on top. The hole
    (0.8–4.3 mm across here, 1.5 mm on Körbelin et al.'s coupons) is far
    below the mesh. Meshed, it did not work: the fibre-break core, which
    takes elements by their centroid, caught none of it at 5 mm (the
    impact point is a mesh node, 3.5 mm from the nearest element centres)
    and the four central elements at 2.5 mm once its radius passed
    1.8 mm (8-ply QI 30 J: 0.894 at 5 mm, 0.603 at 2.5 mm). Weighted by
    area fraction as a near-empty core (stiffness ×0.01) it was smeared
    into the elements around it and did not converge (8-ply QI 15 / 30 J:
    0.898 / 0.864 at 5 mm, 0.852 / 0.737 at 2.5 mm; ×0.30: 0.907 / 0.882
    and 0.881 / 0.807). As a notch on the elements under it, where the
    resolved field is smooth, it converges (see the mesh item below). On
    the one measured TAI dataset (`korbelin_m21_t800s_tai`, M21/T800S at
    25 J) fe3d gives 0.870 against the measured 0.889, an MAE of 3.6%
    (closed-form tiers 3.5%; 0.965 and 9.0% with the zone only, 0.335
    and 62.2% with the zone cut in tension). This is partly in-sample:
    the hole size uses `tai_notch_fraction`, fitted to the same coupons,
    and only the zone around it is fitted elsewhere (on CAI), so fe3d TAI
    is no better validated than the closed-form TAI. The notch scales
    every stress component of the elements under the hole by the
    laminate-level factor, an approximation of the ply stresses at a hole
    edge. A damage state with its own fibre-break core
    (`DamageState.fiber_break_radius_mm`, a C-scan) has it as the hole's
    minimum radius, as in the closed-form tiers. Until issue #164 the
    core switched the hole off in favour of the meshed core, which held
    no element at 5 mm below a 3.5 mm radius, so a measured 1.5 mm core
    made the 8-ply QI 15 J knockdown 0.912 against 0.714 without one
    (closed form 0.727 / 0.771). Now the hole covers the core and the
    meshed core no longer softens the fibres as well (the same broken
    fibres would count twice): 0.673 for the 1.5 mm core, 0.525 / 0.484
    for 3 / 4 mm cores (closed form 0.567 / 0.513; the 4 mm core is
    meshed, four columns, and the hole's notch sits on those elements).
    Use the closed-form tiers for TAI design values.
- **The `fe3d` first-ply failure length scale is fixed, not fitted.** The
  failure index is averaged over a 2.5 mm disk
  (`FPF_AVERAGING_RADIUS_MM`, within the 1-4 mm characteristic distances
  of CFRP notched strength), each element weighted by its overlap with the
  disk. Apart from that length, the stress
  concentration at the edge of the soft damage zone depends on the zone's
  shape and on the panel width, not on its absolute size, so fe3d's
  notch-size effect is weaker than that of the empirical tier's Soutis and
  Whitney-Nuismer models. At the default 5 mm mesh the disk holds only
  the element itself. Cohesive surfaces and proper load-introduction BCs
  are left to a future release.
- **At the default 5 mm mesh a small `fe3d` damage zone is
  under-resolved.** On 8-ply IM7/8552 150 × 100 mm panels at 15 J (damage
  area about 1,533 mm² QI, 1,538 mm² cross-ply) the knockdown is:

    | Case | 5 mm | 2.5 mm | 2.0 mm |
    |---|---|---|---|
    | QI [0/45/−45/90]s, tension | 0.714 | 0.720 | not run |
    | QI, compression | 0.256 | 0.274 | 0.268 |
    | Cross-ply [0/90/0/90]s, tension | 0.782 | 0.783 | not run |
    | Cross-ply, compression | 0.167 | 0.168 | 0.163 |

    In compression 2.5 to 2.0 mm is within about 2–4%, and the 5 mm
    default within about 3–10% of the 2.0 mm value (mostly below it). In
    tension the fibre-damage hole sets the knockdown at the zone's smooth
    centre, and 5 mm is within 1% of 2.5 mm; with the zone only the
    tension rows were 0.912 / 0.930 (QI) and 0.917 / 0.954 (cross-ply),
    5 mm 2–4% below 2.5 mm. Before the damage zone kept its fibres in
    tension they were 0.242 / 0.273 / 0.269 (QI) and 0.156 / 0.168 / 0.162
    (cross-ply). Smaller zones, only a few elements across at 5 mm, are less well
    resolved: with earlier code, at damage areas of about 710–855 mm²,
    the 5 mm cross-ply compression knockdown was about 25% above its
    converged value. The
    calibration coupons' damage areas (2,800–6,500 mm²) are well resolved
    at 5 mm. Use 2.5 mm or finer when the damage area is below about
    1,000 mm² (about 35 s for an 8-ply 150 × 100 mm panel);
    `FE3D_MAX_DOF` (150,000) admits such a panel down to about 1.67 mm.
    Before the averaging disk was weighted by element overlap it counted
    whole elements by their centres and shrank with the mesh, and the
    cross-ply compression knockdown kept falling (0.186 / 0.165 at 2.5 /
    2.0 mm).

- **The delamination growth stress is a 1D thin-film estimate.** It uses
  the straight-sided blister energy release rate over the delamination's
  enclosing rectangle, a fixed mixed-mode phase angle (52.1°; at a buckled
  strip's edge the phase angle moves from about 38° at buckling onset
  towards mode II as the stress rises), and treats the rest of the laminate
  as rigid. Its buckling stress comes from a simply supported rectangle
  scaled by the panel's boundary factor, while the energy release rate
  assumes a clamped film, and the sublaminate's extension-bending
  coupling (B) is ignored although most sublaminates are unsymmetric. It
  is least accurate when the weakest sublaminate is near the
  mid-plane, which is common: the minimum over interfaces often lands on a
  3-5 ply sublaminate. Once a sublaminate buckles at less than about 0.3×
  its steady-state growth stress `√(2E_f·G_c/h)`, the growth stress stays
  between 0.87× and 1× that value, so larger delaminations at the same
  interface lower it little and can raise it slightly. The growth
  knockdown therefore follows impact energy mainly when a different
  interface becomes the weakest. On 150 × 100 mm IM7/8552 panels from 10
  to 30 J it goes from 0.57 to 0.47 on 16-ply QI [0/45/−45/90]4 (0.60 at
  15 J; the weakest interface moves from 0 to 11), from 0.51 to 0.46 on
  8-ply QI [0/45/−45/90]s and from 0.56 to 0.39 on [0/90]2s (0.39 from
  15 J). In `fe3d` the first-ply failure of the softened damage zone falls
  with energy as well: at 5 / 15 / 30 J (5 mm mesh, every energy above the
  onset) the fe3d CAI knockdown is 0.430 / 0.302 / 0.263 on 16-ply QI
  [0/45/−45/90]4 (damage area 187 / 1,270 / 2,893 mm²; empirical 0.651 /
  0.417 / 0.322), 0.333 / 0.256 / 0.225 on 8-ply QI (450 / 1,533 /
  3,156 mm²; empirical 0.546 / 0.395 / 0.312) and 0.307 / 0.167 / 0.134
  on [0/90]2s (456 / 1,538 / 3,162 mm²; empirical 0.544 / 0.394 /
  0.312).
- `fe3d` CAI ignores whole-panel buckling (CAI fixtures carry anti-buckling
  guides). A slender panel without guides can buckle well below the
  reported residual.

## References

- Davies, G.A.O. & Zhang, X. (1995). Impact damage prediction in carbon
  composite structures. *International Journal of Impact Engineering*,
  16(1), 149-170.
- Suemasu, H. & Majima, O. (1996). Multiple delaminations and their
  severity in circular axisymmetric plates subjected to transverse loading.
  *Journal of Composite Materials*, 30(4), 441.
- Olsson, R. (1992). Impact response of orthotropic composite plates
  predicted from a one-parameter differential equation. *AIAA Journal*,
  30(6), 1587-1596.
- Olsson, R. (2001). Analytical prediction of large mass impact damage in
  composite laminates. *Composites Part A*, 32(9), 1207-1215.
- Olsson, R., Donadon, M.V. & Falzon, B.G. (2006). Delamination threshold
  load for dynamic impact on plates. *International Journal of Solids and
  Structures*, 43(10), 3124-3141.
- González, E.V., Maimí, P., Camanho, P.P., Lopes, C.S. & Blanco, N.
  (2011). Effects of ply clustering in laminated composite plates under
  low-velocity impact loading. *Composites Science and Technology*, 71(6),
  805-817.
- Salvetti, M., Gilioli, A., Sbarufatti, C., Dragan, K., Chalimoniuk, M.,
  Manes, A. & Giglio, M. (2016). Analytical model to describe damage in
  CFRP specimen when subjected to low velocity impacts. *Procedia
  Engineering*, 167, 2-9.
- Levy, S. (1942). Bending of rectangular plates with large deflections.
  NACA Report 737.
- Chia, C.Y. (1980). *Nonlinear Analysis of Plates*. McGraw-Hill.
- Shivakumar, K.N., Elber, W. & Illg, W. (1985). Prediction of impact
  force and duration due to low-velocity impact on circular composite
  laminates. *Journal of Applied Mechanics*, 52(3), 674-680.
- Soutis, C. (1996). Compressive strength of unidirectional composites:
  measurement and prediction. *ASTM STP*, 1242, 168-176.
- Whitney, J.M. & Nuismer, R.J. (1974). Stress fracture criteria for
  laminated composites containing stress concentrations.
  *Journal of Composite Materials*, 8(3), 253-265.
- Tan, S.C. (1988). Finite-width correction factors for anisotropic plate
  containing a central opening. *Journal of Composite Materials*, 22(11),
  1080-1097.
- Körbelin, J., Dreiner, C. & Fiedler, B. (2020). Impact of temperature on
  LVI-damage and tensile and compressive residual strength of CFRP.
  *Composites Part C: Open Access*, 3, 100074 (CC BY 4.0).
- Tsai, S.W. & Wu, E.M. (1971). A general theory of strength for anisotropic
  materials. *Journal of Composite Materials*, 5(1), 58-80.
- Davila, C.G., Camanho, P.P., & Rose, C.A. (2005). Failure criteria for FRP
  laminates. NASA/TM-2005-213530 (LaRC05).
- Puck, A. & Schürmann, H. (1998). Failure analysis of FRP laminates by
  means of physically based phenomenological models. *Composites Science
  and Technology*, 58(7), 1045-1067.
- Daniel, I.M. & Ishai, O. (2006). *Engineering Mechanics of Composite
  Materials* (2nd ed.). Oxford University Press.
- Timoshenko, S.P. & Gere, J.M. (1961). *Theory of Elastic Stability*
  (2nd ed.). McGraw-Hill.
- Johnson, K.L. (1985). *Contact Mechanics*. Cambridge University Press
  (Hertz, flat-punch and conical contact).
- Chai, H., Babcock, C.D., & Knauss, W.G. (1981). One dimensional modelling
  of failure in laminated plates by delamination buckling. *International
  Journal of Solids and Structures*, 17(11), 1069-1083.
- Hutchinson, J.W. & Suo, Z. (1992). Mixed mode cracking in layered
  materials. *Advances in Applied Mechanics*, 29, 63-191.
- Lemaitre, J. (1992). *A Course on Damage Mechanics*. Springer.
- Matzenmiller, A., Lubliner, J. & Taylor, R.L. (1995). A constitutive
  model for anisotropic damage in fiber-composites. *Mechanics of
  Materials*, 20(2), 125-152.
- Maimí, P., Camanho, P.P., Mayugo, J.A. & Dávila, C.G. (2007). A
  continuum damage model for composite laminates: Part I – Constitutive
  model. *Mechanics of Materials*, 39(10), 897-908.
- Cantwell, W.J. & Morton, J. (1991). The impact resistance of composite
  materials – a review. *Composites*, 22(5), 347-362.
- Hallett, S.R., Green, B.G., Jiang, W.G. & Wisnom, M.R. (2009). An
  experimental and numerical investigation into the damage mechanisms in
  notched composites. *Composites Part A*, 40(5), 613-624.
- Liu, D. (1988). Impact-induced delamination — a view of bending stiffness
  mismatching. *Journal of Composite Materials*, 22(7), 674-692.
- Abrate, S. (1991). Impact on laminated composite materials. *Applied
  Mechanics Reviews*, 44(4), 155-190.
- Abrate, S. (1998). *Impact on Composite Structures*. Cambridge University
  Press.
- Soutis, C. & Curtis, P.T. (1996). Prediction of the post-impact
  compressive strength of CFRP laminated composites. *Composites Science
  and Technology*, 56(6), 677-684.
- Benzeggagh, M.L. & Kenane, M. (1996). Measurement of mixed-mode
  delamination fracture toughness of unidirectional glass/epoxy composites
  with mixed-mode bending apparatus. *Composites Science and Technology*,
  56(4), 439-449.
