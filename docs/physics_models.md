# Physics Models

BVID-FE exposes three residual-strength modeling tiers that share a common
data contract (`DamageState` → tier engine → `AnalysisResults`) but differ in
fidelity, runtime, and the failure mechanisms they capture. This page
summarises what each tier solves; see the
[Python API reference](python_api.md) for usage.

## Empirical tier

Soutis open-hole-equivalent CAI model relates the delamination-affected stress
concentration around the projected damage area to a net-section failure.
Whitney-Nuismer point-stress and average-stress criteria are used for TAI.
Both models are closed-form and run in milliseconds.

## Semi-analytical tier

At each delaminated interface the thinner sublaminate over the largest
delamination is treated as a plate over the delamination's enclosing
rectangle, and a Rayleigh-Ritz closed form gives its buckling load `N_cr`
(film buckling stress `σ_c = N_cr / h`).

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

- `G_c` is the Benzeggagh-Kenane mix `G_Ic + (G_IIc − G_Ic)(sin²ψ)^η` at
  the thin-film phase angle ψ = 52.1° with η = 2 (0.48 N/mm for IM7/8552).
- `E_f` is the film's x-stress per unit laminate x-strain, with the
  laminate's own Poisson contraction. Film and laminate share the far-field
  strain, so the laminate stress at growth is `Ex · σ_g / E_f`.
- The weakest interface governs, and the result is capped by the Soutis
  empirical strength: `σ_CAI = min(Soutis, σ_growth)`.

Whitney-Nuismer is retained for TAI. The governing interface is
`AnalysisResults.critical_sublaminate` and its buckling load (N/mm) is in
`AnalysisResults.buckling_eigenvalues`.

## 3D FE tier

A structured hexahedral mesh is built for the damaged laminate. Delaminated
interfaces are approximated by a **component-wise stiffness-reduction model**
(true cohesive surfaces deferred to a future release): each damaged element
carries an out-of-plane factor (`DAMAGE_OOP_FACTOR ≈ 0.05`) that scales the
through-thickness and transverse-shear stiffness, while in-plane stiffness is
preserved (the plies themselves remain intact). Inside the fiber-break core
under the impact site, in-plane stiffness is also reduced
(`DAMAGE_FIBER_BREAK_INPLANE_FACTOR ≈ 0.30`) to represent fiber bundle
fracture. First-ply-failure is evaluated at all Gauss points, on stress
rotated into each ply's material frame, using LaRC05 (CAI) and Tsai-Wu (TAI).
For CAI, the lower of first-ply failure and the semi-analytical delamination
growth stress governs. Buckling onset (whole panel or sublaminate, from the
Rayleigh-Ritz closed form, issue #129 — the 3D K_g eigensolve previously
used here was retired because 3D Hex on thin laminates locks too
aggressively at affordable mesh sizes) is reported in
`AnalysisResults.buckling_eigenvalues` but does not set the strength.
Whole-panel buckling is left out of the knockdown because CAI fixtures carry
anti-buckling guides; its stress is recorded in `AnalysisResults.notes`.

## Knockdown definition and cross-tier comparability

`AnalysisResults.knockdown` is computed in exactly one place —
`BvidAnalysis.run()` (`src/bvidfe/analysis/bvid.py`):

```python
knockdown = residual_strength_MPa / pristine_strength_MPa
```

**The denominator is identical across all three tiers.**
`_pristine_strength()` is a thickness-weighted ply-average of the
lamina-level strengths from the material card:

- CAI: `Σ tᵢ (Xc·cos²θᵢ + Yc·sin²θᵢ) / Σ tᵢ`
- TAI: `Σ tᵢ (Xt·cos²θᵢ + Yt·sin²θᵢ) / Σ tᵢ`

**The numerator (residual strength) is what differs between tiers:**

| Tier | CAI residual stress | TAI residual stress |
| --- | --- | --- |
| `empirical` | Soutis: `σ₀ / (1 + k_s·(DPA/A_panel)^m)` | Whitney-Nuismer point-stress on equivalent hole |
| `semi_analytical` | `min(Soutis, σ_growth)` | Delegates to Whitney-Nuismer (identical to `empirical`) |
| `fe3d` | `σ₀ × min(σ_growth, FPF_LaRC05)` damaged ÷ the same quantity undamaged | `σ₀ × FPF_Tsai-Wu` damaged ÷ undamaged |

fe3d's own failure stresses are not the same kind of quantity as σ₀, so
fe3d is normalised by an undamaged fe3d run of the same panel; the raw
stresses are recorded in `AnalysisResults.notes`.

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
      growth stress its knockdown is higher than `semi_analytical`'s — see
      "Limitations" below.
- For **energy-scaling studies**, prefer `empirical` (Soutis scales with
  DPA). The delamination growth stress in `semi_analytical` and `fe3d`
  changes little with delamination size once the sublaminate buckles well
  below it, so those tiers are flatter in energy. `fe3d` is intended for
  stress-field context and through-thickness damage visualization, not
  energy-dependent knockdown curves.

## Limitations

- Material calibration constants (`olsson_alpha`, `soutis_k_s`, `dent_beta`,
  and related parameters) are reasonable defaults for typical CFRP systems.
  Precise values need calibration against material-specific coupon test data
  before use in certification.
- LaRC05 is implemented as a minimal Hashin-3D reduction. Full plane-search
  fiber-kinking is deferred to a future release.
- The `fe3d` tier uses component-wise stiffness reduction at delaminated
  interfaces (in-plane preserved, out-of-plane reduced) instead of true
  cohesive surfaces with bilinear traction-separation laws. Cohesive surfaces
  are deferred to a future release.
- **The `fe3d` tier's knockdown is partially insensitive to impact energy**
  above the Olsson threshold. The delamination growth stress changes
  little with delamination size (see below), and the FPF strain is
  controlled by stress concentration at the healthy/damaged boundary
  rather than damage magnitude. For energy-dependent knockdown curves
  prefer `tier="empirical"`. Full energy-monotonicity (cohesive surfaces +
  proper load-introduction BCs) is v0.3.0 scope.
- **The delamination growth stress is a 1D thin-film estimate.** It uses
  the straight-sided blister energy release rate over the delamination's
  enclosing rectangle, a fixed mixed-mode phase angle (52.1°; at a buckled
  strip's edge the phase angle moves from about 38° at buckling onset
  towards mode II as the stress rises), and treats the rest of the laminate
  as rigid. It is least accurate when the weakest sublaminate is near the
  mid-plane, which is common: the minimum over interfaces often lands on a
  3-5 ply sublaminate. Once a sublaminate buckles at less than about 0.3×
  its steady-state growth stress `√(2E_f·G_c/h)`, the growth stress stays
  between 0.87× and 1× that value, so larger delaminations lower it little
  and can raise it slightly; knockdowns from this channel are therefore
  nearly flat with impact energy.
- `fe3d` CAI ignores whole-panel buckling (CAI fixtures carry anti-buckling
  guides). A slender panel without guides can buckle well below the
  reported residual.
- No validated datasets included; comparison against published Soutis,
  Caprino, and NASA datasets is on the roadmap.

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
- Chai, H., Babcock, C.D., & Knauss, W.G. (1981). One dimensional modelling
  of failure in laminated plates by delamination buckling. *International
  Journal of Solids and Structures*, 17(11), 1069-1083.
- Hutchinson, J.W. & Suo, Z. (1992). Mixed mode cracking in layered
  materials. *Advances in Applied Mechanics*, 29, 63-191.
- Benzeggagh, M.L. & Kenane, M. (1996). Measurement of mixed-mode
  delamination fracture toughness of unidirectional glass/epoxy composites
  with mixed-mode bending apparatus. *Composites Science and Technology*,
  56(4), 439-449.
