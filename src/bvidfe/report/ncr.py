"""Nonconformance Report (NCR) with a proposed MRB disposition.

Decision support for the Streamlit app's Export NCR tab, kept out of
``app.py`` so it can be tested: :func:`recommend_disposition` bands an
:meth:`AnalysisResults.to_dict` result into a proposed Material Review
Board (MRB) disposition path, and :func:`build_ncr` assembles the NCR as a
Markdown document and a machine-readable dict. Both are pure: dicts and
an :class:`~bvidfe.analysis.AnalysisConfig` in, text and dicts out.

Disposition bands
-----------------
The knockdown bands below (``KD_USE_AS_IS`` = 0.97, ``KD_SUBSTANTIATE`` =
0.85, ``KD_REPAIR`` = 0.70) are **placeholder engineering-judgement
values, not taken from a cited standard or program document.** Nothing in
CMH-17, AS9100 or an FAA advisory circular sets them; they only order the
recommendation. Likewise ``BVID_DENT_THRESHOLD_MM`` (0.30 mm) is a common
design damage-tolerance value for the barely-visible dent depth, not a
requirement. A program must replace both with its own damage-tolerance
substantiation, and a qualified MRB issues the disposition: the output is
always marked as proposed and requiring review.

A biaxial analysis (``AnalysisConfig.transverse_stress_ratio`` != 0,
experimental and not validated, issue #115) is recorded as such: the NCR
gives the ratio and the sigma_y at failure, and the recommendation carries
a flag that its knockdown is not validated.
"""

from __future__ import annotations

from datetime import datetime, timezone
from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:  # pragma: no cover
    from bvidfe.analysis.config import AnalysisConfig

__all__ = [
    "BVID_DENT_THRESHOLD_MM",
    "DISPOSITION_STATUS",
    "KD_REPAIR",
    "KD_SUBSTANTIATE",
    "KD_USE_AS_IS",
    "build_ncr",
    "recommend_disposition",
]

#: Residual surface-dent depth [mm] above which impact damage is treated as
#: visible (VID) rather than barely visible (BVID). 0.30 mm (~0.012 in) is a
#: common design damage-tolerance value; placeholder, see the module docstring.
BVID_DENT_THRESHOLD_MM = 0.30

#: Knockdown at or above which, with no delamination and a sub-threshold dent
#: and no fibre breakage, Use-As-Is is proposed. Placeholder band.
KD_USE_AS_IS = 0.97

#: Knockdown at or above which (below ``KD_USE_AS_IS``) Use-As-Is with
#: engineering substantiation, or Repair, is proposed. Placeholder band.
KD_SUBSTANTIATE = 0.85

#: Knockdown at or above which (below ``KD_SUBSTANTIATE``) Repair is
#: proposed; below it, structural Repair or Scrap. Placeholder band.
KD_REPAIR = 0.70

#: Status carried by every recommendation.
DISPOSITION_STATUS = "PROPOSED — REQUIRES MRB REVIEW AND APPROVAL"


def _format_ply_thickness(t) -> str:
    if isinstance(t, (list, tuple)):
        return ", ".join(f"{x:g}" for x in t) + " mm (per ply)"
    return f"{t:g} mm"


def recommend_disposition(
    result_dict: dict, dent_threshold_mm: float = BVID_DENT_THRESHOLD_MM
) -> dict:
    """Heuristic MRB disposition *recommendation* from BVID-FE outputs.

    This is structured decision support only: it bands the predicted
    strength retention (``KD_*`` placeholder bands, see the module
    docstring), raises the governing flags, and cites the criteria a board
    must check. It does NOT issue an approved disposition — a qualified
    Material Review Board must review, modify, and sign.

    Parameters
    ----------
    result_dict : dict
        :meth:`AnalysisResults.to_dict` output: ``knockdown``, ``dpa_mm2``,
        ``residual_strength_MPa``, ``pristine_strength_MPa``, ``tier_used``
        and ``damage`` (``dent_depth_mm``, ``fiber_break_radius_mm``,
        ``delaminations``).
    dent_threshold_mm : float
        Detectability threshold: a dent at or above it is flagged visible.

    Returns
    -------
    dict
        ``recommended_path``, ``rationale``, ``criteria`` (list),
        ``flags`` (list), ``knockdown`` and ``status``.
    """
    kd = float(result_dict["knockdown"])
    dent = float(result_dict["damage"]["dent_depth_mm"])
    fbr = float(result_dict["damage"]["fiber_break_radius_mm"])
    n_delam = len(result_dict["damage"]["delaminations"])
    dpa = float(result_dict["dpa_mm2"])

    flags: list[str] = []
    if dent >= dent_threshold_mm:
        flags.append(
            f"Residual dent {dent:.3f} mm ≥ detectability threshold "
            f"{dent_threshold_mm:.3f} mm — damage is visible (VID, not "
            "BVID); disposition under the visible-damage damage-tolerance "
            "branch."
        )
    if fbr > 0.0:
        flags.append(
            f"Fiber breakage present (r ≈ {fbr:.2f} mm) — primary load "
            "path severed; Use-As-Is is not appropriate without structural "
            "substantiation."
        )
    if n_delam > 0:
        flags.append(
            f"{n_delam} delamination(s) modelled — assess sublaminate "
            "stability and growth under spectrum/limit load."
        )
    if result_dict.get("biaxial_interaction_factor") is not None:
        ratio = result_dict.get("config_snapshot", {}).get("transverse_stress_ratio")
        ratio_txt = "" if ratio is None else f" σy/σx = {float(ratio):+.2f}"
        flags.append(
            f"Biaxial load{ratio_txt} — EXPERIMENTAL, NOT VALIDATED: no biaxial "
            "residual-strength data back this knockdown; substantiate it by test "
            "or a validated analysis before relying on it."
        )

    if fbr > 0.0:
        path = "Repair"
        band = "fiber breakage present"
    elif kd >= KD_USE_AS_IS and dent < dent_threshold_mm and n_delam == 0:
        path = "Use-As-Is (UAI)"
        band = f"knockdown {kd:.3f} ≥ {KD_USE_AS_IS:.2f}, no delamination, sub-threshold dent"
    elif kd >= KD_SUBSTANTIATE:
        path = "Use-As-Is (UAI) with engineering substantiation, or Repair"
        band = f"{KD_SUBSTANTIATE:.2f} ≤ knockdown {kd:.3f} < {KD_USE_AS_IS:.2f}"
    elif kd >= KD_REPAIR:
        path = "Repair"
        band = f"{KD_REPAIR:.2f} ≤ knockdown {kd:.3f} < {KD_SUBSTANTIATE:.2f}"
    else:
        path = "Repair (structural) or Scrap"
        band = f"knockdown {kd:.3f} < {KD_REPAIR:.2f}"

    rationale = (
        f"Predicted residual strength retains {kd * 100:.1f}% of pristine "
        f"({result_dict['residual_strength_MPa']:.1f} MPa vs "
        f"{result_dict['pristine_strength_MPa']:.1f} MPa) at projected "
        f"damage area {dpa:.0f} mm² (analysis tier: "
        f"{result_dict['tier_used']}). Banding: {band}. The board must "
        "confirm a positive Margin of Safety against design ULTIMATE load "
        "(and LIMIT load with no detrimental deformation) for the as-found "
        "damage, including no-growth / arrested-growth substantiation over "
        "the service inspection interval."
    )

    criteria = [
        "Margin of Safety ≥ 0 vs design ultimate load for the as-found "
        "damage state (governing design load case).",
        "Damage-tolerance / no-growth substantiation per CMH-17 Vol. 3 "
        "guidance and the program structural substantiation document.",
        "BVID detectability basis (residual dent depth + NDI/tap-test) "
        f"against the {dent_threshold_mm:.3f} mm design threshold.",
        "Approved Structural Repair Manual (SRM) limits if a Repair path "
        "is selected; engineering disposition required if outside SRM.",
        "Control of nonconforming output per AS9100 §8.7 and the site MRB procedure.",
    ]

    return {
        "recommended_path": path,
        "rationale": rationale,
        "criteria": criteria,
        "flags": flags,
        "knockdown": kd,
        "status": DISPOSITION_STATUS,
    }


def build_ncr(
    form: dict,
    cfg: "AnalysisConfig",
    result_dict: Optional[dict],
    rec: Optional[dict],
    *,
    generated_utc: Optional[str] = None,
) -> tuple[str, dict]:
    """Assemble a structured NCR as ``(Markdown document, machine-readable dict)``.

    Parameters
    ----------
    form : dict
        The NCR fields: ``ncr_number``, ``revision``, ``originator``,
        ``date``, ``program``, ``part_number``, ``part_name``,
        ``serial_lot``, ``work_order``, ``quantity``, ``station``,
        ``found_during``, ``nc_type``, ``dimensions`` and ``description``.
    cfg : AnalysisConfig
        Material, layup, ply thickness and panel of the nonconforming part.
    result_dict : dict or None
        :meth:`AnalysisResults.to_dict` output, or ``None`` when no analysis
        is attached (the analysis section is then marked pending).
    rec : dict or None
        :func:`recommend_disposition` output, or ``None``.
    generated_utc : str, optional
        Generation timestamp; defaults to now (UTC, ISO 8601). Pass one for
        a reproducible document.
    """
    if generated_utc is None:
        generated_utc = datetime.now(timezone.utc).isoformat(timespec="seconds")
    layup_str = "[" + "/".join(f"{a:g}" for a in cfg.layup_deg) + "]"

    analysis: Optional[dict] = None
    if result_dict is not None:
        analysis = {
            "tool": "BVID-FE",
            "tier_used": result_dict["tier_used"],
            "pristine_strength_MPa": result_dict["pristine_strength_MPa"],
            "residual_strength_MPa": result_dict["residual_strength_MPa"],
            "knockdown": result_dict["knockdown"],
            "projected_damage_area_mm2": result_dict["dpa_mm2"],
            "dent_depth_mm": result_dict["damage"]["dent_depth_mm"],
            "fiber_break_radius_mm": result_dict["damage"]["fiber_break_radius_mm"],
            "n_delaminations": len(result_dict["damage"]["delaminations"]),
            "pristine_transverse_stress_MPa": result_dict.get("pristine_transverse_stress_MPa"),
            "residual_transverse_stress_MPa": result_dict.get("residual_transverse_stress_MPa"),
            "biaxial_interaction_factor": result_dict.get("biaxial_interaction_factor"),
            "notes": list(result_dict.get("notes", [])),
            "warnings": list(result_dict.get("warnings", [])),
        }

    ncr = {
        "ncr_number": form["ncr_number"],
        "revision": form["revision"],
        "originated_by": form["originator"],
        "date": form["date"],
        "program": form["program"],
        "part_number": form["part_number"],
        "part_name": form["part_name"],
        "serial_or_lot": form["serial_lot"],
        "work_order": form["work_order"],
        "quantity_affected": form["quantity"],
        "found_at_station": form["station"],
        "found_during": form["found_during"],
        "nonconformance_type": form["nc_type"],
        "nonconformance_description": form["description"],
        "as_found_dimensions": form["dimensions"],
        "material_configuration": {
            "material": cfg.material,
            "layup_deg": list(cfg.layup_deg),
            "layup_string": layup_str,
            "ply_thickness_mm": cfg.ply_thickness_mm,
            "panel_Lx_mm": cfg.panel.Lx_mm,
            "panel_Ly_mm": cfg.panel.Ly_mm,
            "panel_boundary": cfg.panel.boundary,
            "loading": cfg.loading,
            "transverse_stress_ratio": cfg.transverse_stress_ratio,
        },
        "engineering_analysis": analysis,
        "proposed_disposition": rec,
        "generated_utc": generated_utc,
    }

    lines: list[str] = []
    lines.append(f"# Nonconformance Report — {form['ncr_number']}")
    lines.append("")
    lines.append(
        "> **Decision-support document.** The disposition below is a "
        "structured *recommendation* generated from BVID-FE analysis. It "
        "is NOT an approved disposition. A qualified Material Review Board "
        "must review, modify as required, and approve."
    )
    lines.append("")
    lines.append("## 1. Identification")
    lines.append("")
    lines.append(f"- **NCR number:** {form['ncr_number']}  (rev {form['revision']})")
    lines.append(f"- **Originated by:** {form['originator']}")
    lines.append(f"- **Date:** {form['date']}")
    lines.append(f"- **Program / project:** {form['program']}")
    lines.append(f"- **Part number:** {form['part_number']}")
    lines.append(f"- **Part name:** {form['part_name']}")
    lines.append(f"- **Serial / lot:** {form['serial_lot']}")
    lines.append(f"- **Work order:** {form['work_order']}")
    lines.append(f"- **Quantity affected:** {form['quantity']}")
    lines.append(f"- **Found at / station:** {form['station']}")
    lines.append(f"- **Found during:** {form['found_during']}")
    lines.append("")
    lines.append("## 2. Nonconformance")
    lines.append("")
    lines.append(f"- **Type:** {form['nc_type']}")
    lines.append(f"- **As-found dimensions:** {form['dimensions']}")
    lines.append("")
    lines.append(form["description"] or "_No description provided._")
    lines.append("")
    lines.append("## 3. Material / configuration")
    lines.append("")
    lines.append(f"- **Material:** {cfg.material}")
    lines.append(f"- **Layup:** {layup_str} ({len(cfg.layup_deg)} plies)")
    lines.append(f"- **Ply thickness:** {_format_ply_thickness(cfg.ply_thickness_mm)}")
    lines.append(f"- **Panel:** {cfg.panel.Lx_mm:g} × {cfg.panel.Ly_mm:g} mm, {cfg.panel.boundary}")
    if cfg.transverse_stress_ratio != 0:
        lines.append(
            f"- **Biaxial load:** σy/σx = {cfg.transverse_stress_ratio:+.2f} "
            f"({cfg.loading}) — experimental, not validated"
        )
    lines.append("")
    lines.append("## 4. Engineering analysis (BVID-FE)")
    lines.append("")
    if analysis is None:
        lines.append(
            "_Engineering analysis pending — run a BVID-FE analysis and "
            "regenerate this NCR to attach residual-strength results._"
        )
    else:
        lines.append(f"- **Analysis tier:** {analysis['tier_used']}")
        lines.append(f"- **Pristine strength:** {analysis['pristine_strength_MPa']:.1f} MPa")
        lines.append(f"- **Residual strength:** {analysis['residual_strength_MPa']:.1f} MPa")
        lines.append(f"- **Knockdown:** {analysis['knockdown']:.3f}")
        lines.append(
            f"- **Projected damage area:** {analysis['projected_damage_area_mm2']:.0f} mm²"
        )
        lines.append(f"- **Dent depth:** {analysis['dent_depth_mm']:.3f} mm")
        lines.append(f"- **Fiber-break radius:** {analysis['fiber_break_radius_mm']:.2f} mm")
        lines.append(f"- **Delaminations modelled:** {analysis['n_delaminations']}")
        if analysis["biaxial_interaction_factor"] is not None:
            lines.append(
                f"- **σy at failure (experimental):** "
                f"{analysis['pristine_transverse_stress_MPa']:.1f} MPa pristine, "
                f"{analysis['residual_transverse_stress_MPa']:.1f} MPa residual"
            )
            lines.append(
                f"- **Biaxial interaction factor (experimental):** "
                f"{analysis['biaxial_interaction_factor']:.3f}"
            )
        if analysis["notes"]:
            lines.append("- **Runtime notes:**")
            for n in analysis["notes"]:
                lines.append(f"  - {n}")
        if analysis["warnings"]:
            lines.append("- **Warnings:**")
            for w in analysis["warnings"]:
                lines.append(f"  - {w}")
    lines.append("")
    lines.append("## 5. Proposed disposition (requires MRB review)")
    lines.append("")
    if rec is None:
        lines.append("_No disposition recommendation — attach a BVID-FE analysis first._")
    else:
        lines.append(f"- **Status:** {rec['status']}")
        lines.append(f"- **Recommended path:** {rec['recommended_path']}")
        lines.append("")
        lines.append(f"**Rationale.** {rec['rationale']}")
        lines.append("")
        if rec["flags"]:
            lines.append("**Flags raised:**")
            for fl in rec["flags"]:
                lines.append(f"- {fl}")
            lines.append("")
        lines.append("**Criteria the board must verify:**")
        for c in rec["criteria"]:
            lines.append(f"- {c}")
        lines.append("")
    lines.append("## 6. MRB review & approval")
    lines.append("")
    lines.append("| Role | Name | Disposition concurrence | Signature | Date |")
    lines.append("|------|------|-------------------------|-----------|------|")
    lines.append("| Engineering (stress) |  |  |  |  |")
    lines.append("| Quality |  |  |  |  |")
    lines.append("| Materials & Processes |  |  |  |  |")
    lines.append("| MRB chair |  |  |  |  |")
    lines.append("")
    lines.append(f"_Generated {ncr['generated_utc']} UTC by BVID-FE._")
    lines.append("")

    return "\n".join(lines), ncr
