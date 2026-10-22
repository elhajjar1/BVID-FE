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

1. **Onset energy** (Olsson 2001). The threshold load is
   `P_c = π √(8 G_IIc D_eff / 9)`, with `D_eff = √(D11 D22)` from the
   laminate's CLT bending stiffness. The onset energy is
   `E_onset = P_c² / (2 k_cb)`, where `k_cb = 1 / (1/k_b + 1/k_c)`:
    - `k_b` is the point-load bending stiffness of the simply supported
      panel at the impact site (Navier series, 11 × 11 modes), times 1.0
      (`simply_supported`), 2.5 (`clamped`) or 0.4 (`free`, with a
      `UserWarning`).
    - `k_c` is the contact stiffness linearised at `P_c`: Hertz
      sphere-on-flat for `hemispherical`, `2 E* R` for `flat`, and Love's
      cone (30° half-angle) for `conical`. `E*` combines a steel impactor
      with the ply's `E22`.

    Below `E_onset` the damage state is empty.

2. **Damage area.** The target projected damage area is

    ```text
    DPA = α (E − E_onset) · 10³ / G_IIc × (5.5 / m)^0.1 × f_shape × (16 / D)^0.2
    ```

    in mm². It is an energy balance: a fraction `α = olsson_alpha` of the
    impact energy above onset goes into mode-II delamination. The default
    0.078 is fitted to the median of Lovejoy & Scotti's eight IM7/8552
    C-scan areas. T800S/3900-2B (interlayer-toughened) has 0.047, fitted to
    Hasebe et al.'s quasi-isotropic laminates struck by the 12.7 mm
    hemisphere; over all 343 of their specimens the median predicted /
    measured area is 1.03× (0.81–1.31× per layup).

    The area used to be divided by the laminate thickness `h` as well. That
    gave mm rather than mm² and left the area flat in thickness at constant
    energy per thickness, where measured areas grow with thickness:
    Hasebe et al.'s 8- to 24-ply laminates at fixed J/mm, and Lovejoy &
    Scotti's 7.6 mm coupons, which had 1.21× the C-scan area of their
    6.1 mm coupons at the same 8.9 J/mm. Thickness now enters only through
    `E_onset`.

    `m` is the impactor mass (kg), `D` its diameter (mm; the diameter
    factor is 1 for a cone, whose contact is set by its angle) and
    `f_shape` = 1.0 (hemispherical), 0.2 (flat) or 1.3 (conical). `f_shape`
    and the exponent 0.2 are fitted to Hasebe et al. (2025), the only open
    data that vary the impactor on one material: over same-layup,
    same-energy groups a 15.9 mm flat-ended impactor left 0.20× and a 60°
    cone 1.20× the C-scan area of a 12.7 mm hemisphere, and 25.4 and
    31.8 mm hemispheres 0.88× and 0.87×. At the ASTM D7136 reference
    (5.5 kg, 16 mm, hemispherical) the three correction factors are 1. DPA
    is capped at 80% of the panel area (`DPACapClipWarning`). On a
    150 × 100 mm IM7/8552 panel the area grows by about 99 mm² per joule;
    the footprint confinement of step 5 then holds it at about 11,200 mm²
    from about 114 J, before the 12,000 mm² cap. A
    `SmallMassQuasiStaticWarning` is raised when the impactor is lighter
    than the panel.

3. **Dent depth**: `d = h · β · ((E − E_onset) · 10³ / (G_Ic h²))^γ`,
   capped at `h / 2`, with `β = dent_beta` (0.05) and `γ = dent_gamma`
   (0.5).

4. **Fibre-break radius**: `r_fb = η √max(0, E − E_fb)` with
   `η = fiber_break_eta` and `E_fb = fiber_break_E_threshold`. `η` is 0 in
   every preset, so the fibre-break core is off unless a custom material
   or a C-scan file sets it.

5. **Distribution over interfaces** ("peanut" template). Every ply
   interface `i` (`0 … n_plies − 2`) gets one ellipse centred on the impact
   point. Its orientation is the mean of the two adjacent ply angles, its
   aspect ratio is `1 + 0.025 |Δθ|` clipped to [1, 4], and its size weight
   `0.3 + 0.7 (i + 1) / n_interfaces` grows from the impact face to the
   back face. One scale factor is solved so the union of the footprints
   equals the target DPA within 1%.

   Every ellipse is kept inside the impacted panel (`impact_panel`, else
   `panel`): an ellipse whose bounding box would reach past the panel edge
   nearest the impact point is shrunk, keeping its aspect ratio, and the
   scale factor is re-solved for the same union area. If the panel cannot
   hold the target area, the largest footprint that fits is used with a
   `UserWarning`. The union of differently oriented ellipses spreads well
   beyond the circle of the same area: unconfined, the template footprint
   was 99 × 52 mm on the NCAMP quasi-isotropic coupon (same-area circle
   61 mm), 134 × 70 mm on Lovejoy & Scotti's quasi-isotropic coupon
   (circle 82 mm) and 108 × 108 mm on their hard laminate (circle 92 mm),
   wider than the 102 mm coupon. For an ASTM D7136 coupon the impacted
   panel is the 125 × 75 mm support window, which is where the damage
   stays.

The dent depth is carried in `DamageState` and reported, but no tier uses
it in the residual-strength calculation.

## Empirical tier

Both models are closed-form and run in milliseconds. They see the damage
only through its projected area `DPA` (the union of all ellipse
footprints).

**CAI** (after Soutis & Curtis 1996):

```text
σ_CAI = σ₀ / (1 + k_s · (DPA / A_ref)^m)
```

`σ₀` is the pristine strength (see
[Knockdown definition](#knockdown-definition-and-cross-tier-comparability)),
`A_ref` = 15,000 mm² is the area of the 150 × 100 mm ASTM D7137 coupon
(`SOUTIS_REFERENCE_AREA_MM2` in `bvidfe.failure.soutis_openhole`), on which
`k_s = soutis_k_s` (default 4.7, calibrated on measured CAI data, see
`validation/README.md`) is calibrated, and `m = soutis_m` (default 0.5).
`A_ref` does not change with the analysed panel. It used to be the panel
area `Lx · Ly`, which sent the knockdown to 1 as the panel grew, although
impact damage is a local notch, and over-penalised Hasebe et al.'s
80 × 50 mm coupons. The panel size still enters the onset energy, the 80%
area cap, the semi-analytical growth check and fe3d.

**TAI**: Whitney-Nuismer point-stress criterion on an equivalent circular
hole of radius `R = √(DPA / π)`:

```text
σ_TAI / σ₀ = 2 / (2 + ξ² + 3ξ⁴ − (K_T∞ − 3)(5ξ⁶ − 7ξ⁸)),   ξ = R / (R + d₀)
```

`d₀ = wn_d0_mm` (default 1.0 mm) is the characteristic distance, and
`K_T∞ = 1 + √(2(√(Ex/Ey) − ν_xy) + Ex/Gxy)` is Lekhnitskii's
infinite-plate stress concentration for the orthotropic laminate (3.0 for a
quasi-isotropic one). `Ex`, `Ey`, `Gxy` and `ν_xy` come from the
average-stress compliance `a* = h A⁻¹`. When `K_T∞ ≳ 9.2` the polynomial
would give a residual above `σ₀` (or diverge for `K_T∞ ≳ 20.3`), so the
ratio is capped at 1 with a `UserWarning`. The built-in presets stay below
`K_T∞ ≈ 7.3`.

## Semi-analytical tier

For TAI this tier calls the same Whitney-Nuismer model as the empirical
tier, so its tension results are identical. For CAI it adds a
delamination growth bound.

At each delaminated interface the thinner sublaminate (by thickness) over
the largest delamination at that interface is treated as a plate over the
delamination's enclosing rectangle, `a = 2 · major_mm` along x (the load
direction) and `b = 2 · minor_mm` along y; the ellipse's `orientation_deg`
is not used. The Navier (Rayleigh-Ritz) closed form for a
simply supported orthotropic rectangle under uniaxial compression gives its
buckling load per unit width

```text
N_cr = min over m, n ∈ 1..5 of  (π² / a²) · [D11 m⁴ + 2(D12 + 2D66)(m a/b)² n² + D22 (n a/b)⁴] / m²
```

with the sublaminate's CLT bending stiffnesses `Dij`, times a boundary
factor of 1.0 (`simply_supported`), 1.9 (`clamped`) or 0.5 (`free`). The
aspect ratio `a/b` is clipped to 50 (`DegenerateThinSublaminateWarning`).
The film buckling stress is `σ_c = N_cr / h`, where `h` is the
sublaminate thickness.

Buckling alone does not fail the laminate: a thin sublaminate over a
BVID-sized delamination buckles at a few MPa and keeps carrying load. The
CAI strength is the far-field stress at which the buckled sublaminate grows
its delamination. The thin-film energy release rate (Chai, Babcock & Knauss
1981; Hutchinson & Suo 1992, straight-sided blister)

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

The governing interface is `AnalysisResults.critical_sublaminate` and its
buckling load `N_cr` (N/mm) is in `AnalysisResults.buckling_eigenvalues`.
With an empty damage state the tier returns `σ₀`.

With the calibrated `soutis_k_s`, the Soutis bound governs on many
configurations. On 150 × 100 mm IM7/8552 panels (8-ply QI, 16-ply QI and
8-ply cross-ply, 0.152 mm plies, 5–30 J) the growth knockdown is
0.38–0.60 and the Soutis knockdown 0.33–0.57. Growth governs on the 8-ply
QI panel at 5 J (0.500 against 0.556) and on the cross-ply panel at
5–15 J (0.39–0.40 against 0.41–0.57); Soutis governs elsewhere. On the
NCAMP and Lovejoy & Scotti validation cases (`validation/datasets/`)
Soutis governs throughout (on the 40-ply hard laminate the growth stress
is 252 MPa against a 226 MPa Soutis bound), so the two tiers give the
same CAI strengths there. Growth governs on 225 of Hasebe et al.'s 343
coupons (every cross-ply and hard-laminate coupon, and some
quasi-isotropic ones), where the semi-analytical MAE is 23.3% against
21.3% for the empirical tier.

## 3D FE tier

A structured hexahedral mesh of the whole panel is built for the damaged
laminate (`MeshParams`: `in_plane_size_mm`, default 5 mm, and
`elements_per_ply`, default 1). Ply interfaces lie on element faces.
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
  stiffness ×`DAMAGE_ZONE_INPLANE_FACTOR` (0.17, fitted to three measured
  CAI coupons, see [Limitations](#limitations); Soutis & Curtis 1996).
  Without it, delaminations left the in-plane stress field unchanged and
  first-ply failure could not see the impact. An element straddling the
  zone's edge takes the factor weighted by the fraction of its area inside
  the zone (sampled on a 6 × 6 grid), so the zone does not move in
  element-sized steps as the mesh changes.
- The optional fibre-break core (elements within `fiber_break_radius_mm`
  of any delamination centroid; off in every preset) gets both factors:
  out-of-plane ×0.05 and in-plane ×`DAMAGE_FIBER_BREAK_INPLANE_FACTOR`
  (0.30).

The panel is loaded by a prescribed x displacement: `u_x = 0` on the
`x = 0` face and `u_x = ε · Lx` on the `x = Lx` face, with symmetry
restraints `u_y = 0` on `y = 0` and `u_z = 0` on `z = 0`. The panel
`boundary` adds a `u_z = 0` restraint on the two loaded edges
(`simply_supported`), on all four edges (`clamped`) or nowhere (`free`).
One linear static solve at a reference strain of 5% gives the stress field;
because the response is linear, the load multiplier at which the failure
index first reaches 1 follows analytically (Tsai-Wu: root of a quadratic;
LaRC05: `1/√index`; Puck: `1/index`). If nothing fails up to 5% strain, the
5% state is used.

First-ply failure is evaluated on stress rotated into each ply's material
frame, using LaRC05 (CAI) and Tsai-Wu (TAI). The failure index is averaged
over each element's Gauss points (2 × 2 × 2) and then over the elements of
the same through-thickness layer whose centres lie within
`FPF_AVERAGING_RADIUS_MM` = 2.5 mm (an average-stress criterion; Whitney &
Nuismer 1974), clipped at the panel edges; the averaged index is linear in
the coefficients of the load multiplier, so the analytic scaling above
still holds. A check at single Gauss points has no length scale: the stress
at the re-entrant corners of the zone's footprint is singular, and the
pointwise failure load kept moving as the mesh was refined. 2.5 mm is
within the 1-4 mm characteristic distances of CFRP notched strength; at
the default 5 mm mesh the disk holds only the element itself. The solver also implements Puck; `BvidAnalysis` does not
expose a criterion choice, and an unknown criterion name raises
`ValueError`.

- Damage-zone elements fail at the strain intact material would (strain
  equivalence; Lemaitre 1992). The criterion sees their effective stress:
  the in-plane stress the pristine ply would carry at the same strain, with
  the actual (damaged) out-of-plane tractions. Their strength therefore
  drops with their stiffness; with full strengths, a zone spanning most of
  the panel lost its stress concentration and made the panel look stronger.
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

The fe3d run is refused with `FE3DSizeError` above 500,000 degrees of
freedom (`BVIDFE_FE3D_MAX_DOF`).

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
| `empirical` | Soutis: `σ₀ / (1 + k_s·(DPA/A_ref)^m)`, `A_ref` = 15,000 mm² | Whitney-Nuismer point-stress on equivalent hole |
| `semi_analytical` | `min(Soutis, σ_growth)` | Delegates to Whitney-Nuismer (identical to `empirical`) |
| `fe3d` | `σ₀ × min(1, D/U)` with `D`, `U` = `min(σ_growth, FPF_LaRC05)` on the damaged / undamaged panel | `σ₀ × min(1, D/U)` with `D`, `U` = `FPF_Tsai-Wu` damaged / undamaged |

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
      identical; `fe3d` differs.
    - For **CAI**, `semi_analytical ≤ empirical` always (the growth stress
      only lowers the residual). `fe3d` divides by its own undamaged
      strength (first-ply failure, typically below σ₀), so for the same
      growth stress its knockdown is higher than `semi_analytical`'s.
- For **energy-scaling studies**, prefer `empirical` (Soutis scales with
  DPA). `fe3d` also falls with damage area, through its softened damage
  zone (16-ply QI tension: 0.40/0.31/0.31 at 5/15/30 J, against empirical
  0.40/0.37/0.36), but its damage-zone stiffness is fitted to only three
  CAI coupons. Where the delamination growth stress
  governs `semi_analytical`, that tier is flatter in energy: the growth
  stress changes little with delamination size once the sublaminate
  buckles well below it.

## Limitations

- `soutis_k_s` (4.7) is calibrated on two measured CAI datasets (AS4/8552 and
  IM7/8552, both on the 8552 matrix; see `validation/README.md`). Fitted on
  either dataset alone it predicts the other within 6.2–9.0%, so expect
  errors of that order on a new material or layup. The default
  `olsson_alpha` (0.078) is fitted to Lovejoy & Scotti's C-scan areas, the
  T800S/3900-2B value (0.047) and the impactor shape and diameter factors
  to Hasebe et al.'s. The other calibration constants (`soutis_m`,
  `wn_d0_mm`, `dent_beta`, `dent_gamma`, the impactor mass factor, the
  boundary factors, and the fe3d out-of-plane and fibre-break factors)
  are uncalibrated defaults
  chosen for typical CFRP systems. Precise values need calibration against
  material-specific coupon test data before use in certification.
- Validation covers CAI only. TAI has no measured dataset, and fe3d is not
  run on the measured datasets in CI (it gives 171 MPa against the measured
  178 MPa on the first NCAMP coupon, one of the three coupons its
  damage-zone factor is fitted to, `validation/README.md`).
- LaRC05 is implemented as a minimal Hashin-3D reduction: σ₃ and τ₁₃ do
  not enter the index. Full plane-search fiber-kinking is deferred to a
  future release.
- The `fe3d` tier uses component-wise stiffness reduction (out-of-plane at
  delaminated interfaces, in-plane over the damage zone) instead of true
  cohesive surfaces with bilinear traction-separation laws and a continuum
  damage model. Cohesive surfaces are deferred to a future release. The
  damage-zone in-plane factor is uniform over the footprint. Its value,
  0.17, is fitted to three measured CAI coupons, with the delaminations
  kept inside the ASTM D7136 impact window: the first NCAMP AS4/8552 and
  the first Lovejoy & Scotti IM7/8552 quasi-isotropic coupons (3.2% and
  1.2% low) and the first Lovejoy & Scotti 40-ply hard laminate (4.7%
  high), so both laminate families fit within 5%. Before the footprint was
  confined, the hard laminate's template footprint was 108 × 108 mm on a
  102 mm wide coupon: the soft zone spanned the whole section, the 0°
  fibres failed at the coupon's side edge, and no single factor fitted
  both families (0.20 left the hard laminate 47% low). On
  150 × 100 mm IM7/8552 panels in tension at 30 J the fe3d knockdown is
  below the empirical tier's (16-ply QI 0.33 against 0.36, 24-ply QI 0.31
  against 0.36); the damage area
  there (about 2,900 mm²) is far below its 80%-of-panel cap. A damage area
  at the cap would put the soft zone across the loaded edges.
- **The `fe3d` first-ply failure length scale is fixed, not fitted.** The
  failure index is averaged over a 2.5 mm disk
  (`FPF_AVERAGING_RADIUS_MM`, within the 1-4 mm characteristic distances
  of CFRP notched strength). Apart from that length, the stress
  concentration at the edge of the soft damage zone depends on the zone's
  shape and on the panel width, not on its absolute size, so fe3d's
  notch-size effect is weaker than that of the empirical tier's Soutis and
  Whitney-Nuismer models. With the averaging and the area-weighted zone
  the result is mesh-converged: going from 5 mm to 2.5 mm elements changes
  the 15 J knockdown on 8-ply IM7/8552 150 × 100 mm panels by 0.8% (QI
  tension, 0.4209 to 0.4243), 3.1% (cross-ply tension, 0.2296 to 0.2368),
  0.2% (cross-ply compression, 0.2488 to 0.2493) and 0.9% (QI
  compression, 0.4260 to 0.4300). A centroid-assigned zone checked at
  single Gauss points changed by 5.3%, 24.2%, 27.6% and 3.3%, and its
  undamaged QI first-ply failure stress fell from 594.8 to 506.6 MPa (now
  598.2 and 597.2 MPa). At the default 5 mm mesh the disk holds only the
  element itself. Cohesive surfaces and proper load-introduction BCs are
  left to a future release.
- **The delamination growth stress is a 1D thin-film estimate.** It uses
  the straight-sided blister energy release rate over the delamination's
  enclosing rectangle, a fixed mixed-mode phase angle (52.1°; at a buckled
  strip's edge the phase angle moves from about 38° at buckling onset
  towards mode II as the stress rises), and treats the rest of the laminate
  as rigid. It is least accurate when the weakest sublaminate is near the
  mid-plane, which is common: the minimum over interfaces often lands on a
  3-5 ply sublaminate. Once a sublaminate buckles at less than about 0.3×
  its steady-state growth stress `√(2E_f·G_c/h)`, the growth stress stays
  between 0.87× and 1× that value, so larger delaminations at the same
  interface lower it little and can raise it slightly. The growth
  knockdown therefore follows impact energy mainly when a different
  interface becomes the weakest. On 150 × 100 mm IM7/8552 panels from 5 to
  30 J it goes from 0.57 to 0.46 on 16-ply QI (0.60 at 10 J; the weakest
  interface moves from 0 to 11), from 0.50 to 0.46 on 8-ply QI and from
  0.40 to 0.38 on 8-ply cross-ply.
- `fe3d` CAI ignores whole-panel buckling (CAI fixtures carry anti-buckling
  guides). A slender panel without guides can buckle well below the
  reported residual.

## References

- Olsson, R. (2001). Analytical prediction of large mass impact damage in
  composite laminates. *Composites Part A*, 32(9), 1207-1215.
- Soutis, C. (1996). Compressive strength of unidirectional composites:
  measurement and prediction. *ASTM STP*, 1242, 168-176.
- Whitney, J.M. & Nuismer, R.J. (1974). Stress fracture criteria for
  laminated composites containing stress concentrations.
  *Journal of Composite Materials*, 8(3), 253-265.
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
- Chai, H., Babcock, C.D., & Knauss, W.G. (1981). One dimensional modelling
  of failure in laminated plates by delamination buckling. *International
  Journal of Solids and Structures*, 17(11), 1069-1083.
- Hutchinson, J.W. & Suo, Z. (1992). Mixed mode cracking in layered
  materials. *Advances in Applied Mechanics*, 29, 63-191.
- Lemaitre, J. (1992). *A Course on Damage Mechanics*. Springer.
- Soutis, C. & Curtis, P.T. (1996). Prediction of the post-impact
  compressive strength of CFRP laminated composites. *Composites Science
  and Technology*, 56(6), 677-684.
- Benzeggagh, M.L. & Kenane, M. (1996). Measurement of mixed-mode
  delamination fracture toughness of unidirectional glass/epoxy composites
  with mixed-mode bending apparatus. *Composites Science and Technology*,
  56(4), 439-449.
