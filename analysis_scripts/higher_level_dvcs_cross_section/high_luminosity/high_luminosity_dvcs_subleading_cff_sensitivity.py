#!/usr/bin/env python3
"""
High-luminosity CLAS12 DVCS: local subleading-CFF separation study.

Purpose
-------
Use the same matched Pass-2 XS+BSA pseudo-data inputs as the D-term study, but
ask a different question: how well can an unpolarized proton target separate
H, Htilde, and E as luminosity increases?

No dispersion relation and no D-term constraint are used here.  In each
(xB,Q2,t) cell the selected real and imaginary CFF components are varied
independently around the KM15 pseudo-truth and propagated through Gepard to
the measured phi-dependent unpolarized cross section and beam-spin asymmetry.

Fit modes
---------
  H              : ReH, ImH
  H_Ht           : ReH, ImH, ReHt, ImHt
  H_Ht_E         : ReH, ImH, ReHt, ImHt, ReE, ImE

The comparison H_Ht -> H_Ht_E directly quantifies how much apparent Htilde
precision is lost when E is no longer held fixed.

Experimental treatment
----------------------
  * same Stage-2 Pass-2 kinematics and pseudo-data uncertainties
  * XS and BSA used jointly in every fit
  * one shared XS normalization nuisance
  * one shared BSA beam-polarization nuisance
  * statistical uncertainties scale with luminosity in the Stage-2 inputs
  * point-to-point and scale systematics remain as supplied by Stage 2

Outputs
-------
  output/stage4_cff_separation/
      tables/
      figures/

Important
---------
Relative errors can become meaningless near a KM15 zero crossing.  Absolute
uncertainties are always saved.  Relative uncertainties are saved only when
|CFF| exceeds --relative-floor.

This is a local sensitivity/separation projection, not a model-independent
global CFF extraction.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import contextlib
import math
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


LUMI_FACTORS = (1, 2, 5, 10)

FIT_MODES = {
    "H": ("ReH", "ImH"),
    "H_Ht": ("ReH", "ImH", "ReHt", "ImHt"),
    "H_Ht_E": ("ReH", "ImH", "ReHt", "ImHt", "ReE", "ImE"),
}

ALL_CFFS = ("ReH", "ImH", "ReHt", "ImHt", "ReE", "ImE")

LABELS = {
    "ReH": r"$\mathrm{Re}\,\mathcal{H}$",
    "ImH": r"$\mathrm{Im}\,\mathcal{H}$",
    "ReHt": r"$\mathrm{Re}\,\widetilde{\mathcal{H}}$",
    "ImHt": r"$\mathrm{Im}\,\widetilde{\mathcal{H}}$",
    "ReE": r"$\mathrm{Re}\,\mathcal{E}$",
    "ImE": r"$\mathrm{Im}\,\mathcal{E}$",
}

DEFAULT_MIN_NPHI = 3
DEFAULT_REL_STEP = 0.02
DEFAULT_ABS_STEP = 0.05
DEFAULT_RELATIVE_FLOOR = 0.10
SVD_RTOL = 1.0e-11

# RGC A_UL precision model calibrated to Samy Polcher Rafael's preliminary
# Summer-2022 target-spin asymmetry, thesis Fig. 5.10.  The values below are
# digitized directly from the vector error bars in that figure.  The thesis says
# those bars contain DVCS-yield counting statistics plus the statistical part of
# the pi0 subtraction uncertainty; several final systematics were not yet included.
#
# Projection convention requested for the workshop:
#   Su22 thesis sample = S
#   already collected RGC data = 3 S (conservative allowance for Fa22 FTOFF)
#   remaining approved beam time at nominal performance = another 3 S
# Therefore the future remaining-time scenarios 0.5x, 1x, 2x give total
# statistics 4.5 S, 6 S, 9 S, respectively.
DEFAULT_CURRENT_EQUIV_SU22 = 3.0
DEFAULT_AUL_SCALE_FRAC = 0.05
DEFAULT_AUT_SCALE_FRAC = 0.05
RGC_REMAINING_FACTORS = (0.5, 1.0, 2.0)

SAMY_SU22_AUL_ERRORS = [
    (3.25, 0.24, 0.24, 7.61, 0.08008),
    (3.25, 0.24, 0.24, 22.25, 0.08341),
    (3.25, 0.24, 0.24, 36.98, 0.10208),
    (3.25, 0.24, 0.24, 53.52, 0.12695),
    (3.25, 0.24, 0.24, 87.43, 0.11771),
    (3.25, 0.24, 0.24, 193.87, 0.11655),
    (3.25, 0.24, 0.24, 305.36, 0.13115),
    (3.25, 0.24, 0.24, 325.67, 0.10923),
    (3.25, 0.24, 0.24, 347.10, 0.06266),
    (3.61, 0.31, 0.42, 7.45, 0.09003),
    (3.61, 0.31, 0.42, 22.42, 0.09114),
    (3.61, 0.31, 0.42, 37.12, 0.11304),
    (3.61, 0.31, 0.42, 59.08, 0.13470),
    (3.61, 0.31, 0.42, 94.51, 0.14624),
    (3.61, 0.31, 0.42, 140.55, 0.14577),
    (3.61, 0.31, 0.42, 209.27, 0.14317),
    (3.61, 0.31, 0.42, 285.09, 0.14624),
    (3.61, 0.31, 0.42, 320.41, 0.14583),
    (3.61, 0.31, 0.42, 335.24, 0.10530),
    (3.61, 0.31, 0.42, 351.52, 0.08206),
    (4.01, 0.40, 1.05, 7.35, 0.09926),
    (4.01, 0.40, 1.05, 22.31, 0.11677),
    (4.01, 0.40, 1.05, 39.09, 0.16151),
    (4.01, 0.40, 1.05, 67.57, 0.17079),
    (4.01, 0.40, 1.05, 102.24, 0.20688),
    (4.01, 0.40, 1.05, 135.08, 0.18231),
    (4.01, 0.40, 1.05, 174.46, 0.19173),
    (4.01, 0.40, 1.05, 222.65, 0.20515),
    (4.01, 0.40, 1.05, 266.20, 0.23513),
    (4.01, 0.40, 1.05, 301.09, 0.17832),
    (4.01, 0.40, 1.05, 324.59, 0.15486),
    (4.01, 0.40, 1.05, 347.40, 0.07766),
    (2.57, 0.17, 0.22, 7.52, 0.10494),
    (2.57, 0.17, 0.22, 22.68, 0.09016),
    (2.57, 0.17, 0.22, 37.21, 0.08875),
    (2.57, 0.17, 0.22, 51.62, 0.11146),
    (2.57, 0.17, 0.22, 71.15, 0.13341),
    (2.57, 0.17, 0.22, 199.22, 0.11598),
    (2.57, 0.17, 0.22, 301.79, 0.13649),
    (2.57, 0.17, 0.22, 316.23, 0.09362),
    (2.57, 0.17, 0.22, 330.76, 0.09213),
    (2.57, 0.17, 0.22, 348.65, 0.08667),
    (2.71, 0.19, 0.40, 7.86, 0.11820),
    (2.71, 0.19, 0.40, 23.31, 0.09997),
    (2.71, 0.19, 0.40, 37.48, 0.09604),
    (2.71, 0.19, 0.40, 51.01, 0.14097),
    (2.71, 0.19, 0.40, 114.13, 0.12890),
    (2.71, 0.19, 0.40, 242.75, 0.11942),
    (2.71, 0.19, 0.40, 309.08, 0.13570),
    (2.71, 0.19, 0.40, 323.38, 0.09301),
    (2.71, 0.19, 0.40, 343.45, 0.07664),
    (2.90, 0.21, 0.85, 8.11, 0.13511),
    (2.90, 0.21, 0.85, 23.04, 0.09521),
    (2.90, 0.21, 0.85, 36.45, 0.10682),
    (2.90, 0.21, 0.85, 84.86, 0.14699),
    (2.90, 0.21, 0.85, 190.97, 0.15889),
    (2.90, 0.21, 0.85, 283.94, 0.14984),
    (2.90, 0.21, 0.85, 319.22, 0.12172),
    (2.90, 0.21, 0.85, 333.20, 0.08476),
    (2.90, 0.21, 0.85, 348.60, 0.10004),
    (1.96, 0.17, 0.20, 7.87, 0.08792),
    (1.96, 0.17, 0.20, 22.17, 0.08825),
    (1.96, 0.17, 0.20, 36.49, 0.10972),
    (1.96, 0.17, 0.20, 52.48, 0.15289),
    (1.96, 0.17, 0.20, 91.12, 0.14054),
    (1.96, 0.17, 0.20, 187.31, 0.13188),
    (1.96, 0.17, 0.20, 304.56, 0.17063),
    (1.96, 0.17, 0.20, 322.42, 0.12384),
    (1.96, 0.17, 0.20, 336.97, 0.09912),
    (1.96, 0.17, 0.20, 352.19, 0.09169),
    (1.95, 0.20, 0.33, 7.86, 0.09480),
    (1.95, 0.20, 0.33, 21.98, 0.09441),
    (1.95, 0.20, 0.33, 36.09, 0.13678),
    (1.95, 0.20, 0.33, 69.26, 0.17604),
    (1.95, 0.20, 0.33, 128.66, 0.16172),
    (1.95, 0.20, 0.33, 201.92, 0.15759),
    (1.95, 0.20, 0.33, 299.31, 0.17113),
    (1.95, 0.20, 0.33, 328.70, 0.12915),
    (1.95, 0.20, 0.33, 347.63, 0.07397),
    (1.92, 0.22, 0.78, 7.72, 0.10620),
    (1.92, 0.22, 0.78, 21.80, 0.11549),
    (1.92, 0.22, 0.78, 45.31, 0.16601),
    (1.92, 0.22, 0.78, 103.39, 0.17649),
    (1.92, 0.22, 0.78, 164.88, 0.17447),
    (1.92, 0.22, 0.78, 241.71, 0.20191),
    (1.92, 0.22, 0.78, 291.36, 0.18799),
    (1.92, 0.22, 0.78, 321.02, 0.14555),
    (1.92, 0.22, 0.78, 345.62, 0.06871),
    (1.64, 0.11, 0.19, 10.54, 0.12171),
    (1.64, 0.11, 0.19, 27.60, 0.10130),
    (1.64, 0.11, 0.19, 42.33, 0.08750),
    (1.64, 0.11, 0.19, 56.63, 0.09607),
    (1.64, 0.11, 0.19, 71.33, 0.11542),
    (1.64, 0.11, 0.19, 105.15, 0.12415),
    (1.64, 0.11, 0.19, 247.05, 0.12274),
    (1.64, 0.11, 0.19, 291.06, 0.11709),
    (1.64, 0.11, 0.19, 305.59, 0.08419),
    (1.64, 0.11, 0.19, 319.54, 0.08766),
    (1.64, 0.11, 0.19, 340.37, 0.08580),
    (1.69, 0.12, 0.32, 15.91, 0.14339),
    (1.69, 0.12, 0.32, 33.53, 0.09154),
    (1.69, 0.12, 0.32, 46.97, 0.08426),
    (1.69, 0.12, 0.32, 61.58, 0.12184),
    (1.69, 0.12, 0.32, 112.27, 0.13310),
    (1.69, 0.12, 0.32, 243.03, 0.12962),
    (1.69, 0.12, 0.32, 299.75, 0.12393),
    (1.69, 0.12, 0.32, 314.20, 0.08059),
    (1.69, 0.12, 0.32, 332.51, 0.07933),
    (1.71, 0.12, 0.68, 9.21, 0.17690),
    (1.71, 0.12, 0.68, 26.41, 0.09579),
    (1.71, 0.12, 0.68, 38.64, 0.08503),
    (1.71, 0.12, 0.68, 56.38, 0.14516),
    (1.71, 0.12, 0.68, 159.59, 0.15565),
    (1.71, 0.12, 0.68, 280.68, 0.14461),
    (1.71, 0.12, 0.68, 315.85, 0.11232),
    (1.71, 0.12, 0.68, 329.18, 0.08057),
    (1.71, 0.12, 0.68, 344.58, 0.11960),
]


def samy_su22_sigma_aul(Q2, xB, t_abs, phi_deg):
    """Nearest-cell/nearest-phi empirical Su22 A_UL uncertainty from thesis Fig. 5.10."""
    a = np.asarray(SAMY_SU22_AUL_ERRORS, dtype=float)
    # First choose the nearest of Samy's 12 (Q2,xB,|t|) cells.  Scales roughly
    # reflect the spacing of the published cells so no one coordinate dominates.
    cells = np.unique(a[:, :3], axis=0)
    d2 = ((cells[:,0]-Q2)/1.0)**2 + ((cells[:,1]-xB)/0.10)**2 + ((cells[:,2]-t_abs)/0.30)**2
    q, x, t = cells[np.argmin(d2)]
    m = (np.isclose(a[:,0],q) & np.isclose(a[:,1],x) & np.isclose(a[:,2],t))
    b = a[m]
    # Use nearest measured phi rather than interpolating across large phi gaps;
    # this is deliberately conservative and keeps the input tied to a real error bar.
    dphi = np.abs((b[:,3] - phi_deg + 180.0) % 360.0 - 180.0)
    return float(b[np.argmin(dphi),4]), (float(q),float(x),float(t),float(b[np.argmin(dphi),3]))

def rgc_scenarios(current_equiv_su22=DEFAULT_CURRENT_EQUIV_SU22):
    out = [("Su22 thesis", 1.0), ("current collected", current_equiv_su22)]
    for f in RGC_REMAINING_FACTORS:
        out.append((f"remaining {f:g}x", current_equiv_su22 * (1.0 + f)))
    return out


def make_point(g, row, helicity: int = 0):
    phi_trento = math.pi - math.radians(float(row.phi_deg))
    pt = g.DataPoint(
        xB=float(row.xB),
        t=-abs(float(row.t_abs)),
        Q2=float(row.Q2),
        phi=float(phi_trento),
        observable="XS",
        frame="trento",
        process="ep2epgamma",
        exptype="fixed target",
        in1energy=float(row.ebeam),
        in1charge=-1,
        in1polarization=int(helicity),
        in2particle="p",
    )
    pt.prepare()
    return pt


def make_aul_point(g, row):
    """Longitudinal-target A_UL point at the same RGA kinematics."""
    phi_trento = math.pi - math.radians(float(row.phi_deg))
    pt = g.DataPoint(
        xB=float(row.xB),
        t=-abs(float(row.t_abs)),
        Q2=float(row.Q2),
        phi=float(phi_trento),
        observable="TSA",
        frame="trento",
        process="ep2epgamma",
        exptype="fixed target",
        in1energy=float(row.ebeam),
        in1charge=-1,
        in1polarization=0,
        in2particle="p",
        in2polarization=1,
        in2polarizationvector="L",
    )
    pt.prepare()
    return pt


def aul_value(th, pt) -> float:
    """Gepard longitudinal target-spin asymmetry."""
    if hasattr(th, "AUL"):
        return float(th.AUL(pt))
    return float(th.predict(pt))


def make_aut_point(g, row):
    """
    Transverse-target A_UT point at the same RGA kinematics.

    Gepard uses in2polarizationvector="T" for A_UT.  We explicitly set
    varFTn=-1, selecting Gepard's dominant sine(varphi) transverse-spin
    component while retaining phi, so the observable is still evaluated at
    each measured DVCS azimuth.  In Gepard's BMK implementation this sets
    varphi=pi/2 inside the transverse-target cross section.  This is the
    E-sensitive transverse-target handle used for the RGH projection.
    """
    phi_trento = math.pi - math.radians(float(row.phi_deg))
    pt = g.DataPoint(
        xB=float(row.xB),
        t=-abs(float(row.t_abs)),
        Q2=float(row.Q2),
        phi=float(phi_trento),
        observable="AUT",
        frame="trento",
        process="ep2epgamma",
        exptype="fixed target",
        in1energy=float(row.ebeam),
        in1charge=-1,
        in1polarization=0,
        in2particle="p",
        in2polarization=1,
        in2polarizationvector="T",
        # Gepard transverse-target convention: varFTn=-1 selects the
        # dominant sine(varphi) component.  With phi retained, AUT remains
        # differential in the DVCS azimuth phi; internally Gepard evaluates
        # the transverse-spin term at varphi=pi/2 (BMK convention).
        varFTn=-1,
    )
    pt.prepare()
    return pt


def aut_value(th, pt) -> float:
    """Gepard transverse target-spin asymmetry A_UT."""
    if hasattr(th, "AUT"):
        return float(th.AUT(pt))
    return float(th.predict(pt))


def cff_owner(th, name: str):
    candidates = [th.m] if hasattr(th, "m") else []
    candidates.append(th)
    seen = set()
    for obj in candidates:
        if obj is None or id(obj) in seen:
            continue
        seen.add(id(obj))
        if hasattr(obj, name) and callable(getattr(obj, name)):
            return obj
    raise AttributeError(f"Could not find callable Gepard CFF method {name}")


@contextlib.contextmanager
def shifted_cff(th, name: str, shift: float):
    owner = cff_owner(th, name)
    old = getattr(owner, name)
    setattr(
        owner,
        name,
        lambda pt, _old=old, _shift=float(shift): float(_old(pt)) + _shift,
    )
    try:
        yield
    finally:
        setattr(owner, name, old)


def pred(th, pt) -> float:
    return float(th.predict(pt))


def cff_value(th, name: str, pt) -> float:
    return float(getattr(cff_owner(th, name), name)(pt))


def bsa_from_points(th, plus_pt, minus_pt) -> float:
    sp = pred(th, plus_pt)
    sm = pred(th, minus_pt)
    den = sp + sm
    if not np.isfinite(den) or abs(den) < 1e-30:
        return np.nan
    return (sp - sm) / den


def finite_derivative(
    th,
    points: Sequence,
    observable: str,
    name: str,
    nominal_cff: float,
    rel_step: float,
    abs_step: float,
) -> Tuple[np.ndarray, float]:
    step = max(abs_step, rel_step * max(abs(nominal_cff), 1.0))

    def evaluate():
        if observable == "xs":
            return np.asarray([pred(th, p) for p in points], dtype=float)
        if observable == "aul":
            return np.asarray([aul_value(th, p) for p in points], dtype=float)
        if observable == "aut":
            return np.asarray([aut_value(th, p) for p in points], dtype=float)
        return np.asarray(
            [bsa_from_points(th, pair[0], pair[1]) for pair in points],
            dtype=float,
        )

    with shifted_cff(th, name, +step):
        plus = evaluate()
    with shifted_cff(th, name, -step):
        minus = evaluate()

    return (plus - minus) / (2.0 * step), step


def load_inputs(indir: Path) -> Dict[int, pd.DataFrame]:
    out = {}
    for factor in LUMI_FACTORS:
        path = indir / f"joint_fit_input_km15_{factor}x.csv"
        if not path.exists():
            raise FileNotFoundError(path)
        d = pd.read_csv(path)
        d["bin"] = d["bin"].astype(int)
        out[factor] = d.sort_values(["bin", "phi_deg"]).reset_index(drop=True)

    ref = out[1][["point_id", "bin"]].astype(str)
    for factor in LUMI_FACTORS[1:]:
        if not ref.equals(out[factor][["point_id", "bin"]].astype(str)):
            raise RuntimeError("Luminosity input kinematics/order mismatch.")
    return out


def build_derivatives(df, th, g, rel_step, abs_step, min_nphi):
    point_rows, cell_rows = [], []
    groups = list(df.groupby("bin", sort=True))
    print(f"[CFF separation] calculating derivatives for {len(groups)} cells")

    for icell, (bid, cell) in enumerate(groups, 1):
        cell = cell.sort_values("phi_deg")
        if len(cell) < min_nphi:
            print(f"[skip] bin {bid}: only {len(cell)} phi points")
            continue

        xs_pts = [make_point(g, r, 0) for r in cell.itertuples(index=False)]
        bsa_pts = [
            (make_point(g, r, +1), make_point(g, r, -1))
            for r in cell.itertuples(index=False)
        ]
        aul_pts = [make_aul_point(g, r) for r in cell.itertuples(index=False)]
        aut_pts = [make_aut_point(g, r) for r in cell.itertuples(index=False)]
        pt0 = xs_pts[0]

        cff0 = {name: cff_value(th, name, pt0) for name in ALL_CFFS}
        xs0 = np.asarray([pred(th, p) for p in xs_pts], dtype=float)
        a0 = np.asarray([bsa_from_points(th, p, m) for p, m in bsa_pts], dtype=float)
        aul0 = np.asarray([aul_value(th, p) for p in aul_pts], dtype=float)
        aut0 = np.asarray([aut_value(th, p) for p in aut_pts], dtype=float)

        dx, da, dul, dut, steps = {}, {}, {}, {}, {}
        for name in ALL_CFFS:
            dx[name], steps[name] = finite_derivative(
                th, xs_pts, "xs", name, cff0[name], rel_step, abs_step
            )
            da[name], _ = finite_derivative(
                th, bsa_pts, "bsa", name, cff0[name], rel_step, abs_step
            )
            dul[name], _ = finite_derivative(
                th, aul_pts, "aul", name, cff0[name], rel_step, abs_step
            )
            dut[name], _ = finite_derivative(
                th, aut_pts, "aut", name, cff0[name], rel_step, abs_step
            )

        # A transparent observable-level diagnostic: response to a +10% change
        # in Im(Htilde).  This is not a fit result.
        delta_imht_10 = 0.10 * cff0["ImHt"]
        xs_frac_response_10 = np.full(len(cell), np.nan)
        good_xs = np.isfinite(xs0) & (np.abs(xs0) > 1e-30)
        xs_frac_response_10[good_xs] = (
            dx["ImHt"][good_xs] * delta_imht_10 / xs0[good_xs]
        )
        bsa_abs_response_10 = da["ImHt"] * delta_imht_10

        for j, row in enumerate(cell.itertuples(index=False)):
            rec = {
                "point_id": row.point_id,
                "bin": int(bid),
                "phi_deg": float(row.phi_deg),
                "xs_km15": float(xs0[j]),
                "bsa_km15": float(a0[j]),
                "aul_km15": float(aul0[j]),
                "aut_km15": float(aut0[j]),
                "xs_frac_response_to_10pct_ImHt": float(xs_frac_response_10[j]),
                "bsa_abs_response_to_10pct_ImHt": float(bsa_abs_response_10[j]),
            }
            for name in ALL_CFFS:
                rec[f"d_xs_d_{name}"] = float(dx[name][j])
                rec[f"d_bsa_d_{name}"] = float(da[name][j])
                rec[f"d_aul_d_{name}"] = float(dul[name][j])
                rec[f"d_aut_d_{name}"] = float(dut[name][j])
            point_rows.append(rec)

        xi = float(np.median(cell["xB"])) / (2.0 - float(np.median(cell["xB"])))
        meta = {
            "bin": int(bid),
            "xB": float(np.median(cell["xB"])),
            "xi": xi,
            "Q2": float(np.median(cell["Q2"])),
            "t_abs": float(np.median(cell["t_abs"])),
            "t_over_Q2": float(np.median(cell["t_over_Q2"])),
            "n_phi": len(cell),
            "median_abs_xs_frac_response_to_10pct_ImHt": float(
                np.nanmedian(np.abs(xs_frac_response_10))
            ),
            "median_abs_bsa_response_to_10pct_ImHt": float(
                np.nanmedian(np.abs(bsa_abs_response_10))
            ),
        }
        for name in ALL_CFFS:
            meta[f"{name}_KM15"] = cff0[name]
            meta[f"finite_step_{name}"] = steps[name]
        cell_rows.append(meta)

        if icell % 10 == 0 or icell == len(groups):
            print(f"[CFF separation] {icell:3d}/{len(groups)} cells")

    return pd.DataFrame(point_rows), pd.DataFrame(cell_rows)


def _derivative_worker(payload):
    """Process one chunk of cells with an independent Gepard model instance."""
    chunk, rel_step, abs_step, min_nphi = payload
    import gepard as g
    from gepard.fits import th_KM15
    return build_derivatives(chunk, th_KM15, g, rel_step, abs_step, min_nphi)


def build_derivatives_parallel(df, rel_step, abs_step, min_nphi, workers):
    """Split complete kinematic cells across independent processes."""
    workers = max(1, min(int(workers), 4))
    if workers == 1:
        import gepard as g
        from gepard.fits import th_KM15
        return build_derivatives(df, th_KM15, g, rel_step, abs_step, min_nphi)

    bins = sorted(df["bin"].astype(int).unique())
    chunks = []
    for ib in np.array_split(np.asarray(bins, dtype=int), workers):
        if len(ib):
            chunks.append(df[df["bin"].isin(ib.tolist())].copy())
    print(f"[CFF separation] distributing {len(bins)} cells across {len(chunks)} workers")
    point_parts, cell_parts = [], []
    payloads = [(c, rel_step, abs_step, min_nphi) for c in chunks]
    with ProcessPoolExecutor(max_workers=len(chunks)) as ex:
        futs = [ex.submit(_derivative_worker, x) for x in payloads]
        for i, fut in enumerate(as_completed(futs), 1):
            a, b = fut.result()
            point_parts.append(a); cell_parts.append(b)
            print(f"[CFF separation] worker chunk {i}/{len(futs)} complete")
    points = pd.concat(point_parts, ignore_index=True).sort_values(["bin","phi_deg"]).reset_index(drop=True)
    cells = pd.concat(cell_parts, ignore_index=True).sort_values("bin").reset_index(drop=True)
    return points, cells


def layout_for(cell_meta, mode, include_aul=False, include_aut=False):
    cffs = FIT_MODES[mode]
    pars = []
    for b in cell_meta["bin"].astype(int):
        for name in cffs:
            pars.append((b, name))
    pars += [(-1, "xs_scale_beta"), (-1, "bsa_scale_beta")]
    if include_aul:
        pars += [(-1, "aul_scale_beta")]
    if include_aut:
        pars += [(-1, "aut_scale_beta")]
    return pars


def covariance(
    df,
    deriv,
    cell_meta,
    mode,
    include_aul=False,
    aul_stats_multiplier=1.0,
    aul_scenario="none",
    aul_scale_frac=DEFAULT_AUL_SCALE_FRAC,
    include_aut=False,
    aut_stats_multiplier=1.0,
    aut_scenario="none",
    aut_scale_frac=DEFAULT_AUT_SCALE_FRAC,
):
    use = df.merge(deriv, on=["point_id", "bin", "phi_deg"], how="inner")
    layout = layout_for(cell_meta, mode, include_aul=include_aul, include_aut=include_aut)
    idx = {k: i for i, k in enumerate(layout)}
    npar = len(layout)
    cffs = FIT_MODES[mode]

    rows, sigmas = [], []

    for r in use.itertuples(index=False):
        v = np.zeros(npar)
        b = int(r.bin)
        for name in cffs:
            v[idx[(b, name)]] = getattr(r, f"d_xs_d_{name}")
        v[idx[(-1, "xs_scale_beta")]] = r.xs_scale_frac * r.xs_km15
        srow = math.hypot(r.xs_stat_pseudo_abs, r.xs_ptp_sys_pseudo_abs)
        if np.isfinite(srow) and srow > 0:
            rows.append(v)
            sigmas.append(srow)

    for r in use.itertuples(index=False):
        v = np.zeros(npar)
        b = int(r.bin)
        for name in cffs:
            v[idx[(b, name)]] = getattr(r, f"d_bsa_d_{name}")
        v[idx[(-1, "bsa_scale_beta")]] = r.bsa_scale_frac * r.bsa_km15
        srow = math.hypot(r.bsa_stat_pseudo_abs, r.bsa_ptp_sys_pseudo_abs)
        if np.isfinite(srow) and srow > 0:
            rows.append(v)
            sigmas.append(srow)

    if include_aul:
        # Empirical Su22 precision from Samy's thesis Fig. 5.10.  Statistical
        # precision scales as 1/sqrt(N); it is independent of the RGA luminosity
        # factor used for the unpolarized XS+BSA projection.
        for r in use.itertuples(index=False):
            v = np.zeros(npar)
            b = int(r.bin)
            for name in cffs:
                v[idx[(b, name)]] = getattr(r, f"d_aul_d_{name}")
            v[idx[(-1, "aul_scale_beta")]] = aul_scale_frac * r.aul_km15
            su22_sigma, _ = samy_su22_sigma_aul(r.Q2, r.xB, r.t_abs, r.phi_deg)
            srow = su22_sigma / math.sqrt(aul_stats_multiplier)
            if np.isfinite(srow) and srow > 0:
                rows.append(v)
                sigmas.append(srow)

    if include_aut:
        # RGH has not run yet.  For this sensitivity study assume its A_UT
        # statistical precision at matched kinematics is the same as the
        # corresponding RGC A_UL scenario.  We therefore reuse the empirical
        # Samy-Su22 error map only as an absolute precision template.
        for r in use.itertuples(index=False):
            v = np.zeros(npar)
            b = int(r.bin)
            for name in cffs:
                v[idx[(b, name)]] = getattr(r, f"d_aut_d_{name}")
            v[idx[(-1, "aut_scale_beta")]] = aut_scale_frac * r.aut_km15
            su22_sigma, _ = samy_su22_sigma_aul(r.Q2, r.xB, r.t_abs, r.phi_deg)
            srow = su22_sigma / math.sqrt(aut_stats_multiplier)
            if np.isfinite(srow) and srow > 0:
                rows.append(v)
                sigmas.append(srow)

    A = np.asarray(rows, dtype=float)
    sigmas = np.asarray(sigmas, dtype=float)
    B = A / sigmas[:, None]

    prior_keys = [(-1, "xs_scale_beta"), (-1, "bsa_scale_beta")]
    if include_aul:
        prior_keys.append((-1, "aul_scale_beta"))
    if include_aut:
        prior_keys.append((-1, "aut_scale_beta"))
    priors = []
    for key in prior_keys:
        p = np.zeros(npar)
        p[idx[key]] = 1.0
        priors.append(p)
    B = np.vstack([B] + priors)

    U, sv, Vt = np.linalg.svd(B, full_matrices=False)
    tol = max(sv[0], 1.0) * SVD_RTOL
    keep = sv > tol
    rank = int(np.sum(keep))
    rank_deficit = npar - rank

    inv_s2 = np.zeros_like(sv)
    inv_s2[keep] = 1.0 / (sv[keep] ** 2)
    cov = (Vt.T * inv_s2) @ Vt

    unresolved_fraction = np.zeros(npar)
    if rank_deficit > 0:
        null = Vt[~keep, :]
        unresolved_fraction = np.sum(null * null, axis=0)

    cond = float(sv[keep][0] / sv[keep][-1]) if rank > 1 else np.inf
    diag = {
        "fit_mode": mode,
        "include_aul": bool(include_aul),
        "aul_stats_multiplier_vs_su22": float(aul_stats_multiplier) if include_aul else np.nan,
        "aul_scenario": str(aul_scenario) if include_aul else "none",
        "aul_scale_frac": float(aul_scale_frac) if include_aul else np.nan,
        "include_aut": bool(include_aut),
        "aut_stats_multiplier_vs_su22": float(aut_stats_multiplier) if include_aut else np.nan,
        "aut_scenario": str(aut_scenario) if include_aut else "none",
        "aut_scale_frac": float(aut_scale_frac) if include_aut else np.nan,
        "n_measurement_rows": len(A),
        "n_parameters": npar,
        "rank": rank,
        "rank_deficit": rank_deficit,
        "condition_number_effective": cond,
        "smallest_kept_singular_value": float(sv[keep][-1]) if rank else np.nan,
        "largest_singular_value": float(sv[0]) if len(sv) else np.nan,
    }
    return cov, layout, unresolved_fraction, diag


def summarize(cell_meta, cov, layout, unresolved, mode, factor, relative_floor):
    idx = {k: i for i, k in enumerate(layout)}
    cffs = FIT_MODES[mode]
    rows = []

    for c in cell_meta.itertuples(index=False):
        b = int(c.bin)
        row = {
            "fit_mode": mode,
            "luminosity_factor": factor,
            "bin": b,
            "xB": c.xB,
            "xi": c.xi,
            "Q2": c.Q2,
            "t_abs": c.t_abs,
            "t_over_Q2": c.t_over_Q2,
            "n_phi": c.n_phi,
        }

        local_indices = [idx[(b, name)] for name in cffs]

        for name in ALL_CFFS:
            truth = getattr(c, f"{name}_KM15")
            row[f"{name}_KM15"] = truth
            row[f"sigma_{name}"] = np.nan
            row[f"relative_sigma_{name}"] = np.nan
            row[f"unresolved_fraction_{name}"] = np.nan

        for name in cffs:
            j = idx[(b, name)]
            sig = math.sqrt(max(cov[j, j], 0.0))
            truth = getattr(c, f"{name}_KM15")
            row[f"sigma_{name}"] = sig
            row[f"unresolved_fraction_{name}"] = unresolved[j]
            if abs(truth) >= relative_floor and unresolved[j] < 1e-8:
                row[f"relative_sigma_{name}"] = sig / abs(truth)

        # Pairwise local correlations among all floated CFFs.
        for i, n1 in enumerate(cffs):
            j1 = idx[(b, n1)]
            for n2 in cffs[i + 1 :]:
                j2 = idx[(b, n2)]
                den = math.sqrt(max(cov[j1, j1] * cov[j2, j2], 0.0))
                row[f"corr_{n1}_{n2}"] = cov[j1, j2] / den if den > 0 else np.nan

        # Maximum absolute correlation of each CFF with another floated CFF.
        for name in cffs:
            vals = []
            j1 = idx[(b, name)]
            for other in cffs:
                if other == name:
                    continue
                j2 = idx[(b, other)]
                den = math.sqrt(max(cov[j1, j1] * cov[j2, j2], 0.0))
                if den > 0:
                    vals.append(abs(cov[j1, j2] / den))
            row[f"max_abs_corr_{name}"] = max(vals) if vals else np.nan

        rows.append(row)

    return pd.DataFrame(rows)


def savefig(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=250)
    plt.close(fig)


def finite_rel(d, col):
    return d[np.isfinite(d[col])].copy()


def make_plots(results, cell_meta, deriv, figures):
    # 1: Im Htilde relative precision vs xB, H+Htilde fit.
    fig, ax = plt.subplots(figsize=(8.8, 6.0))
    for factor in LUMI_FACTORS:
        d = results[
            (results.fit_mode == "H_Ht") &
            (results.luminosity_factor == factor)
        ]
        d = finite_rel(d, "relative_sigma_ImHt")
        ax.scatter(d.xB, 100 * d.relative_sigma_ImHt, s=30, alpha=.70, label=f"{factor}x")
    ax.set_xlabel(r"$x_B$")
    ax.set_ylabel(r"Relative uncertainty on $\mathrm{Im}\,\widetilde{\mathcal{H}}$ (%)")
    ax.set_title(r"Local $H+\widetilde H$ fit: luminosity dependence")
    ax.grid(alpha=.2); ax.legend(title="Luminosity")
    savefig(fig, figures / "01_ImHt_relative_uncertainty_vs_xB.png")

    # 2: same versus xi.
    fig, ax = plt.subplots(figsize=(8.8, 6.0))
    for factor in LUMI_FACTORS:
        d = results[
            (results.fit_mode == "H_Ht") &
            (results.luminosity_factor == factor)
        ]
        d = finite_rel(d, "relative_sigma_ImHt")
        ax.scatter(d.xi, 100 * d.relative_sigma_ImHt, s=30, alpha=.70, label=f"{factor}x")
    ax.set_xlabel(r"$\xi \simeq x_B/(2-x_B)$")
    ax.set_ylabel(r"Relative uncertainty on $\mathrm{Im}\,\widetilde{\mathcal{H}}$ (%)")
    ax.set_title(r"Does $\widetilde H$ sensitivity emerge at high $\xi$?")
    ax.grid(alpha=.2); ax.legend(title="Luminosity")
    savefig(fig, figures / "02_ImHt_relative_uncertainty_vs_xi.png")

    # 3: freeing E penalty at 10x.
    fig, ax = plt.subplots(figsize=(8.8, 6.0))
    for mode, label in [
        ("H_Ht", r"$H+\widetilde H$"),
        ("H_Ht_E", r"$H+\widetilde H+E$"),
    ]:
        d = results[(results.fit_mode == mode) & (results.luminosity_factor == 10)]
        d = finite_rel(d, "relative_sigma_ImHt")
        ax.scatter(d.xB, 100 * d.relative_sigma_ImHt, s=32, alpha=.72, label=label)
    ax.set_xlabel(r"$x_B$")
    ax.set_ylabel(r"Relative uncertainty on $\mathrm{Im}\,\widetilde{\mathcal{H}}$ (%)")
    ax.set_title(r"CFF-separation penalty when $\mathcal{E}$ is also free at 10x")
    ax.grid(alpha=.2); ax.legend()
    savefig(fig, figures / "03_ImHt_HHt_vs_HHtE_10x.png")

    # 4: absolute uncertainty, robust through zero crossings.
    fig, ax = plt.subplots(figsize=(8.8, 6.0))
    for factor in LUMI_FACTORS:
        d = results[
            (results.fit_mode == "H_Ht_E") &
            (results.luminosity_factor == factor)
        ]
        ax.scatter(d.xB, d.sigma_ImHt, s=30, alpha=.70, label=f"{factor}x")
    ax.set_xlabel(r"$x_B$")
    ax.set_ylabel(r"Absolute uncertainty on $\mathrm{Im}\,\widetilde{\mathcal{H}}$")
    ax.set_title(r"Full local $H+\widetilde H+E$ fit")
    ax.grid(alpha=.2); ax.legend(title="Luminosity")
    savefig(fig, figures / "04_ImHt_absolute_uncertainty_full_fit.png")

    # 5: 10x kinematic map.  Marker size inversely tracks absolute uncertainty.
    d = results[(results.fit_mode == "H_Ht_E") & (results.luminosity_factor == 10)].copy()
    finite = np.isfinite(d.sigma_ImHt) & (d.sigma_ImHt > 0)
    d = d[finite]
    inv = 1.0 / d.sigma_ImHt
    sizes = 25 + 175 * (inv - inv.min()) / max(inv.max() - inv.min(), 1e-30)
    fig, ax = plt.subplots(figsize=(8.8, 6.2))
    sc = ax.scatter(d.xB, d.t_abs, c=d.sigma_ImHt, s=sizes, alpha=.78)
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label(r"$\sigma[\mathrm{Im}\,\widetilde{\mathcal{H}}]$")
    ax.set_xlabel(r"$x_B$")
    ax.set_ylabel(r"$|t|$ (GeV$^2$)")
    ax.set_title(r"Where does XS+BSA constrain $\widetilde H$?  Full fit, 10x")
    ax.grid(alpha=.2)
    savefig(fig, figures / "05_ImHt_sensitivity_map_10x.png")

    # 6: maximum local CFF correlation for ImHt.
    fig, ax = plt.subplots(figsize=(8.8, 6.0))
    for mode, label in [
        ("H_Ht", r"$H+\widetilde H$"),
        ("H_Ht_E", r"$H+\widetilde H+E$"),
    ]:
        d = results[(results.fit_mode == mode) & (results.luminosity_factor == 10)]
        col = "max_abs_corr_ImHt"
        ax.scatter(d.xB, d[col], s=30, alpha=.72, label=label)
    ax.set_xlabel(r"$x_B$")
    ax.set_ylabel(r"Maximum $|\rho|$ involving $\mathrm{Im}\,\widetilde{\mathcal{H}}$")
    ax.set_ylim(0, 1.03)
    ax.set_title(r"Is $\widetilde H$ precision limited by CFF separation?")
    ax.grid(alpha=.2); ax.legend()
    savefig(fig, figures / "06_ImHt_max_correlation_10x.png")

    # 7: observable-level response to a +10% ImHt change.
    fig, ax = plt.subplots(figsize=(8.8, 6.0))
    ax.scatter(
        cell_meta.xB,
        100 * cell_meta.median_abs_xs_frac_response_to_10pct_ImHt,
        s=34, alpha=.72, label="XS: median fractional response"
    )
    ax.set_xlabel(r"$x_B$")
    ax.set_ylabel(r"Median $|\Delta\sigma/\sigma|$ for a 10% change in $\mathrm{Im}\,\widetilde{\mathcal{H}}$ (%)")
    ax.set_title(r"Direct observable sensitivity to $\mathrm{Im}\,\widetilde{\mathcal{H}}$")
    ax.grid(alpha=.2)
    savefig(fig, figures / "07_XS_response_to_10pct_ImHt.png")

    fig, ax = plt.subplots(figsize=(8.8, 6.0))
    ax.scatter(
        cell_meta.xB,
        cell_meta.median_abs_bsa_response_to_10pct_ImHt,
        s=34, alpha=.72
    )
    ax.set_xlabel(r"$x_B$")
    ax.set_ylabel(r"Median $|\Delta A_{LU}|$ for a 10% change in $\mathrm{Im}\,\widetilde{\mathcal{H}}$")
    ax.set_title(r"Direct BSA sensitivity to $\mathrm{Im}\,\widetilde{\mathcal{H}}$")
    ax.grid(alpha=.2)
    savefig(fig, figures / "08_BSA_response_to_10pct_ImHt.png")

    # 9: H degradation as progressively more CFFs are released.
    fig, ax = plt.subplots(figsize=(8.8, 6.0))
    for mode, label in [
        ("H", r"$H$ only"),
        ("H_Ht", r"$H+\widetilde H$"),
        ("H_Ht_E", r"$H+\widetilde H+E$"),
    ]:
        d = results[(results.fit_mode == mode) & (results.luminosity_factor == 10)]
        d = finite_rel(d, "relative_sigma_ImH")
        ax.scatter(d.xB, 100 * d.relative_sigma_ImH, s=30, alpha=.70, label=label)
    ax.set_xlabel(r"$x_B$")
    ax.set_ylabel(r"Relative uncertainty on $\mathrm{Im}\,\mathcal{H}$ (%)")
    ax.set_title(r"Cost of progressively relaxing the $H$-dominance assumption at 10x")
    ax.grid(alpha=.2); ax.legend()
    savefig(fig, figures / "09_ImH_degradation_when_subleading_CFFs_float.png")



def make_presentation_cff_band_figure(results, figures):
    """Presentation figure: constraint collapse from complementary observables.

    This deliberately does *not* connect independent local CFF fits into a
    pseudo-continuous CFF band.  Each panel summarizes the distribution over
    all controlled kinematic cells.  Points are medians and vertical bars are
    the central 68% of the cell-by-cell projected uncertainties.

    ImH and ImHt use sigma/|truth| where the KM15 truth is safely nonzero.
    ImE uses absolute sigma because KM15 has ImE=0 in these cells.
    Polarized scenarios use nominal completion (remaining 1x = 6 Su22-equivalent
    statistics).  For each observable set, 1x and 10x RGA are shown side by side.
    """
    specs = [
        ("XS+BSA", "none", "none", r"$\sigma+A_{LU}$"),
        ("XS+BSA+AUL", "remaining 1x", "none", r"$+A_{UL}$"),
        ("XS+BSA+AUL+AUT", "remaining 1x", "RGH same as remaining 1x", r"$+A_{UT}$"),
    ]
    panels = [
        ("sigma_ImH", "ImH_KM15", True, r"$\mathrm{Im}\,\mathcal{H}$", "Relative uncertainty (%)"),
        ("sigma_ImHt", "ImHt_KM15", True, r"$\mathrm{Im}\,\widetilde{\mathcal{H}}$", "Relative uncertainty (%)"),
        ("sigma_ImE", "ImE_KM15", False, r"$\mathrm{Im}\,\mathcal{E}$", "Absolute uncertainty"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(14.8, 4.9))
    xpos = np.arange(len(specs), dtype=float)
    offsets = {1: -0.11, 10: +0.11}
    markers = {1: "o", 10: "s"}

    for ax, (sig_col, truth_col, relative, title, ylabel) in zip(axes, panels):
        for L in (1, 10):
            meds=[]; lo=[]; hi=[]
            for obs, aul, aut, _ in specs:
                d = results[(results.fit_mode == "H_Ht_E") &
                            (results.observable_set == obs) &
                            (results.luminosity_factor == L) &
                            (results.aul_scenario == aul) &
                            (results.aut_scenario == aut)].copy()
                vals = d[sig_col].to_numpy(float)
                if relative:
                    truth = np.abs(d[truth_col].to_numpy(float))
                    good = np.isfinite(vals) & np.isfinite(truth) & (truth >= RELATIVE_TRUTH_FLOOR)
                    vals = 100.0 * vals[good] / truth[good]
                else:
                    vals = vals[np.isfinite(vals)]
                if len(vals):
                    q16,q50,q84=np.percentile(vals,[16,50,84])
                else:
                    q16=q50=q84=np.nan
                meds.append(q50); lo.append(q50-q16); hi.append(q84-q50)
            x=xpos+offsets[L]
            ax.errorbar(x, meds, yerr=np.vstack([lo,hi]), fmt=markers[L], ms=6,
                        capsize=4, lw=1.5, label=f"RGA {L}x")
        ax.set_xticks(xpos)
        ax.set_xticklabels([q[3] for q in specs])
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", alpha=.18)
        ax.legend(frameon=False, fontsize=9)
    fig.suptitle(r"Complementary observables break CFF degeneracies; luminosity then improves precision")
    fig.text(.5,.012,
             r"Full local $H+\widetilde H+E$ fit. Points: median across controlled cells; bars: central 68% of cell-by-cell projected uncertainties. "
             r"RGC/RGH use nominal completed polarized-target precision (6$\times$ Su22-equivalent).",
             ha="center", fontsize=8.5)
    fig.tight_layout(rect=(0,.06,1,.94))
    fig.savefig(figures / "presentation_CFF_constraint_collapse_1x_vs_10x.png", dpi=300)
    plt.close(fig)


def make_representative_low_t_cff_band_diagnostic(results, figures):
    """Optional diagnostic only: same old envelope idea in 0.20<|t|<0.40.

    Kept because this is the high-precision region requested for inspection,
    but it is explicitly labelled as a diagnostic: the cells differ in Q2 and t
    and therefore do not define a continuous CFF-vs-xi function.
    """
    base = results[(results.fit_mode == "H_Ht_E") &
                   (results.t_abs >= 0.20) & (results.t_abs <= 0.40)].copy()
    specs = [
        ("XS+BSA", "none", "none", r"$\sigma+A_{LU}$"),
        ("XS+BSA+AUL", "remaining 1x", "none", r"$+A_{UL}$ (RGC)"),
        ("XS+BSA+AUL+AUT", "remaining 1x", "RGH same as remaining 1x", r"$+A_{UT}$ (RGH)"),
    ]
    cffs = [("ImH_KM15","sigma_ImH",r"$\mathrm{Im}\,\mathcal{H}$"),
            ("ImHt_KM15","sigma_ImHt",r"$\mathrm{Im}\,\widetilde{\mathcal{H}}$"),
            ("ImE_KM15","sigma_ImE",r"$\mathrm{Im}\,\mathcal{E}$")]
    fig,axes=plt.subplots(1,3,figsize=(15,4.9),sharex=True)
    for ax,(truth_col,sig_col,title) in zip(axes,cffs):
        for obs,aul,aut,label in specs:
            for L,alpha,ls in [(1,.10,"--"),(10,.25,"-")]:
                d=base[(base.observable_set==obs)&(base.luminosity_factor==L)&
                       (base.aul_scenario==aul)&(base.aut_scenario==aut)].sort_values("xi")
                d=d[np.isfinite(d[truth_col])&np.isfinite(d[sig_col])]
                if d.empty: continue
                x=d.xi.to_numpy(float); y=d[truth_col].to_numpy(float); e=d[sig_col].to_numpy(float)
                ax.fill_between(x,y-e,y+e,alpha=alpha,linewidth=0,label=(f"{label}, {L}x" if L==10 else None))
                ax.plot(x,y-e,ls=ls,lw=.7,alpha=.55); ax.plot(x,y+e,ls=ls,lw=.7,alpha=.55)
        d0=base[(base.observable_set=="XS+BSA+AUL+AUT")&(base.luminosity_factor==10)&
                (base.aul_scenario=="remaining 1x")].sort_values("xi")
        ax.plot(d0.xi,d0[truth_col],"o-",ms=3,lw=1.2,label="KM15 pseudo-truth")
        ax.set_title(title); ax.set_xlabel(r"$\xi \simeq x_B/(2-x_B)$"); ax.grid(alpha=.15)
    axes[0].set_ylabel("CFF value with projected 68% interval")
    axes[-1].legend(fontsize=8)
    fig.suptitle(r"Diagnostic only: local CFF intervals in the high-precision $0.20<|t|<0.40$ GeV$^2$ region")
    fig.text(.5,.01,r"Cells differ in $Q^2$ and $t$; connected envelopes are visual guides only, not a continuous CFF fit.",ha="center",fontsize=9)
    fig.tight_layout(rect=(0,.045,1,.94))
    fig.savefig(figures / "diagnostic_CFF_bands_t0p2_0p4_1x_vs_10x.png",dpi=300)
    plt.close(fig)

def print_summary(results, diagnostics):
    print("\n" + "=" * 108)
    print("LOCAL H / Htilde / E CFF-SEPARATION PROJECTION")
    print("=" * 108)

    rows = []
    for mode in FIT_MODES:
        for factor in LUMI_FACTORS:
            d = results[
                (results.fit_mode == mode) &
                (results.luminosity_factor == factor)
            ]
            diag = diagnostics[
                (diagnostics.fit_mode == mode) &
                (diagnostics.luminosity_factor == factor)
            ].iloc[0]
            rec = {
                "mode": mode,
                "L": f"{factor}x",
                "rank_def": int(diag.rank_deficit),
                "cond": diag.condition_number_effective,
                "med_rel_ImH_%": 100 * np.nanmedian(d.relative_sigma_ImH),
                "med_rel_ImHt_%": np.nan,
                "med_abs_ImHt": np.nan,
                "med_maxcorr_ImHt": np.nan,
            }
            if "ImHt" in FIT_MODES[mode]:
                rec["med_rel_ImHt_%"] = 100 * np.nanmedian(d.relative_sigma_ImHt)
                rec["med_abs_ImHt"] = np.nanmedian(d.sigma_ImHt)
                rec["med_maxcorr_ImHt"] = np.nanmedian(d.max_abs_corr_ImHt)
            rows.append(rec)

    print(pd.DataFrame(rows).to_string(index=False, float_format=lambda x: f"{x:.4g}"))

    # Explicitly quantify the penalty from freeing E.
    print("\n10x Htilde separation penalty from freeing E:")
    a = results[(results.fit_mode == "H_Ht") & (results.luminosity_factor == 10)]
    b = results[(results.fit_mode == "H_Ht_E") & (results.luminosity_factor == 10)]
    m = a[["bin", "sigma_ImHt"]].merge(
        b[["bin", "sigma_ImHt"]], on="bin", suffixes=("_HHt", "_HHtE")
    )
    good = np.isfinite(m.sigma_ImHt_HHt) & np.isfinite(m.sigma_ImHt_HHtE) & (m.sigma_ImHt_HHt > 0)
    ratio = m.loc[good, "sigma_ImHt_HHtE"] / m.loc[good, "sigma_ImHt_HHt"]
    if len(ratio):
        print(
            f"  median sigma(ImHt) ratio [H+Ht+E]/[H+Ht] = {np.median(ratio):.3g}"
        )
        print(
            f"  central 68% range of ratio = "
            f"{np.percentile(ratio,16):.3g} -- {np.percentile(ratio,84):.3g}"
        )


def main(argv: Optional[List[str]] = None) -> int:
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser()
    p.add_argument(
        "--input-dir",
        default=str(here / "output" / "stage2" / "tables"),
        help="Directory containing joint_fit_input_km15_{1,2,5,10}x.csv",
    )
    p.add_argument(
        "--outdir",
        default=str(here / "output" / "stage4_cff_separation"),
    )
    p.add_argument("--min-nphi", type=int, default=DEFAULT_MIN_NPHI)
    p.add_argument("--workers", type=int, default=1, choices=(1,2,3,4), help="Processes for Gepard derivative calculation (max 4).")
    p.add_argument("--force-derivatives", action="store_true", help="Recompute Gepard derivatives even if cached tables exist.")
    p.add_argument("--relative-step", type=float, default=DEFAULT_REL_STEP)
    p.add_argument("--absolute-step", type=float, default=DEFAULT_ABS_STEP)
    p.add_argument(
        "--current-rgc-equiv-su22",
        type=float,
        default=DEFAULT_CURRENT_EQUIV_SU22,
        help="Effective already-collected RGC statistics in units of Samy's Su22 sample; default 3.",
    )
    p.add_argument(
        "--aul-scale-frac",
        type=float,
        default=DEFAULT_AUL_SCALE_FRAC,
        help="Common fractional A_UL scale uncertainty; default 0.05.",
    )
    p.add_argument(
        "--aut-scale-frac",
        type=float,
        default=DEFAULT_AUT_SCALE_FRAC,
        help="Assumed common fractional A_UT scale uncertainty for RGH; default 0.05.",
    )
    p.add_argument(
        "--relative-floor",
        type=float,
        default=DEFAULT_RELATIVE_FLOOR,
        help="Do not quote relative CFF error when |CFF_truth| is below this.",
    )
    args = p.parse_args(argv)

    indir = Path(args.input_dir).resolve()
    outdir = Path(args.outdir).resolve()
    tables = outdir / "tables"
    figures = outdir / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    inputs = load_inputs(indir)
    one = inputs[1]

    try:
        import gepard as g
        from gepard.fits import th_KM15
    except Exception as exc:
        raise RuntimeError("Could not import Gepard/KM15.") from exc
    th = th_KM15

    # Preflight verifies the exact Gepard interface and demonstrates observable
    # response before the expensive derivative pass.
    row = one.iloc[len(one) // 2]
    xs0 = make_point(g, row, 0)
    bp = make_point(g, row, +1)
    bm = make_point(g, row, -1)
    ul = make_aul_point(g, row)
    ut = make_aut_point(g, row)

    print("[preflight] representative KM15 point")
    print("[preflight] AUT convention: transverse target, varFTn=-1 (sine-varphi); "
          "Gepard BMK sets varphi=pi/2 while retaining differential phi")
    print(
        f"[preflight] xB={row.xB:.4g}, Q2={row.Q2:.4g}, "
        f"|t|={row.t_abs:.4g}, phi={row.phi_deg:.4g}"
    )
    print(
        f"[preflight] XS={pred(th,xs0):.8g}, "
        f"BSA={bsa_from_points(th,bp,bm):+.8g}, "
        f"AUL={aul_value(th,ul):+.8g}, "
        f"AUT={aut_value(th,ut):+.8g}"
    )
    for name in ALL_CFFS:
        val = cff_value(th, name, xs0)
        step = max(args.absolute_step, args.relative_step * max(abs(val), 1.0))
        base_xs = pred(th, xs0)
        base_a = bsa_from_points(th, bp, bm)
        base_ul = aul_value(th, ul)
        base_ut = aut_value(th, ut)
        with shifted_cff(th, name, step):
            xs1 = pred(th, xs0)
            a1 = bsa_from_points(th, bp, bm)
            ul1 = aul_value(th, ul)
            ut1 = aut_value(th, ut)
        print(
            f"[preflight] {name:4s}: value={val:+.6g}, step={step:.4g}, "
            f"XS response={100*(xs1/base_xs-1):+.4g}%, "
            f"BSA change={a1-base_a:+.4g}, "
            f"AUL change={ul1-base_ul:+.4g}, "
            f"AUT change={ut1-base_ut:+.4g}"
        )

    deriv_path = tables / "point_xs_bsa_allcff_derivatives.csv"
    cell_path = tables / "cell_km15_allcffs.csv"
    required_deriv_cols = {"aul_km15", "d_aul_d_ReH", "d_aul_d_ImHt", "d_aul_d_ImE",
                           "aut_km15", "d_aut_d_ReH", "d_aut_d_ImHt", "d_aut_d_ImE"}
    can_reuse = deriv_path.exists() and cell_path.exists() and not args.force_derivatives
    if can_reuse:
        deriv = pd.read_csv(deriv_path)
        cell_meta = pd.read_csv(cell_path)
        if not required_deriv_cols.issubset(deriv.columns):
            can_reuse = False
            print("[cache] derivative table predates AUL/AUT support; recomputing")
    if can_reuse:
        print(f"[cache] reusing {deriv_path} ({len(deriv)} points); skipping Gepard derivative pass")
    else:
        deriv, cell_meta = build_derivatives_parallel(
            one, args.relative_step, args.absolute_step, args.min_nphi, args.workers
        )
        deriv.to_csv(deriv_path, index=False)
        cell_meta.to_csv(cell_path, index=False)
        print(f"[cache] saved derivatives for future reruns: {deriv_path}")

    result_frames = []
    diag_rows = []

    for mode in FIT_MODES:
        for factor in LUMI_FACTORS:
            print(f"[fit] mode={mode:7s} luminosity={factor}x  XS+BSA")
            cov, layout, unresolved, diag = covariance(
                inputs[factor], deriv, cell_meta, mode, include_aul=False
            )
            diag["luminosity_factor"] = factor
            diag_rows.append(diag)
            rr = summarize(
                cell_meta, cov, layout, unresolved, mode, factor, args.relative_floor
            )
            rr["observable_set"] = "XS+BSA"
            rr["aul_scenario"] = "none"
            rr["aul_stats_multiplier_vs_su22"] = np.nan
            rr["aut_scenario"] = "none"
            rr["aut_stats_multiplier_vs_su22"] = np.nan
            result_frames.append(rr)

            # A_UL scan is most relevant once Htilde is released.
            if "ImHt" in FIT_MODES[mode]:
                for scenario, n_su22 in rgc_scenarios(args.current_rgc_equiv_su22):
                    print(
                        f"[fit] mode={mode:7s} luminosity={factor}x  "
                        f"XS+BSA+AUL  {scenario} ({n_su22:g}xSu22 stats)"
                    )
                    cov, layout, unresolved, diag = covariance(
                        inputs[factor],
                        deriv,
                        cell_meta,
                        mode,
                        include_aul=True,
                        aul_stats_multiplier=n_su22,
                        aul_scenario=scenario,
                        aul_scale_frac=args.aul_scale_frac,
                    )
                    diag["luminosity_factor"] = factor
                    diag_rows.append(diag)
                    rr = summarize(
                        cell_meta, cov, layout, unresolved, mode, factor,
                        args.relative_floor
                    )
                    rr["observable_set"] = "XS+BSA+AUL"
                    rr["aul_scenario"] = scenario
                    rr["aul_stats_multiplier_vs_su22"] = n_su22
                    rr["aut_scenario"] = "none"
                    rr["aut_stats_multiplier_vs_su22"] = np.nan
                    result_frames.append(rr)

                    # Add the future RGH transverse-target measurement.  The
                    # requested working assumption is that RGH achieves the same
                    # matched-bin statistical precision as RGC A_UL for this
                    # scenario.
                    print(
                        f"[fit] mode={mode:7s} luminosity={factor}x  "
                        f"XS+BSA+AUL+AUT  AUL={scenario}; "
                        f"RGH AUT same precision ({n_su22:g}xSu22-equivalent)"
                    )
                    cov, layout, unresolved, diag = covariance(
                        inputs[factor], deriv, cell_meta, mode,
                        include_aul=True,
                        aul_stats_multiplier=n_su22,
                        aul_scenario=scenario,
                        aul_scale_frac=args.aul_scale_frac,
                        include_aut=True,
                        aut_stats_multiplier=n_su22,
                        aut_scenario=f"RGH same as {scenario}",
                        aut_scale_frac=args.aut_scale_frac,
                    )
                    diag["luminosity_factor"] = factor
                    diag_rows.append(diag)
                    rr = summarize(
                        cell_meta, cov, layout, unresolved, mode, factor,
                        args.relative_floor
                    )
                    rr["observable_set"] = "XS+BSA+AUL+AUT"
                    rr["aul_scenario"] = scenario
                    rr["aul_stats_multiplier_vs_su22"] = n_su22
                    rr["aut_scenario"] = f"RGH same as {scenario}"
                    rr["aut_stats_multiplier_vs_su22"] = n_su22
                    result_frames.append(rr)

    results = pd.concat(result_frames, ignore_index=True)
    diagnostics = pd.DataFrame(diag_rows)

    results.to_csv(tables / "cff_separation_cell_uncertainties.csv", index=False)
    diagnostics.to_csv(tables / "cff_separation_fit_diagnostics.csv", index=False)

    # Save a compact direct-sensitivity table for easy inspection.
    sens_cols = [
        "bin", "xB", "xi", "Q2", "t_abs", "t_over_Q2", "n_phi",
        "ReH_KM15", "ImH_KM15", "ReHt_KM15", "ImHt_KM15",
        "ReE_KM15", "ImE_KM15",
        "median_abs_xs_frac_response_to_10pct_ImHt",
        "median_abs_bsa_response_to_10pct_ImHt",
    ]
    cell_meta[sens_cols].to_csv(
        tables / "ImHt_direct_observable_sensitivity.csv", index=False
    )

    # Focused A_UL figures using the empirical Su22 precision calibration.
    scenarios = rgc_scenarios(args.current_rgc_equiv_su22)
    for mode in ("H_Ht", "H_Ht_E"):
        fig, ax = plt.subplots(figsize=(8.8, 6.0))
        base = results[(results.fit_mode == mode) & (results.observable_set == "XS+BSA")]
        med = base.groupby("luminosity_factor").sigma_ImHt.median()
        ax.plot(med.index, med.values, marker="o", label="XS+BSA only")
        for scenario, n_su22 in scenarios:
            d = results[(results.fit_mode == mode) &
                        (results.observable_set == "XS+BSA+AUL") &
                        (results.aul_scenario == scenario)]
            med = d.groupby("luminosity_factor").sigma_ImHt.median()
            ax.plot(med.index, med.values, marker="o", label=f"+ AUL: {scenario}")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xticks(LUMI_FACTORS); ax.set_xticklabels([f"{x}x" for x in LUMI_FACTORS])
        ax.set_xlabel("Unpolarized RGA luminosity")
        ax.set_ylabel(r"Median absolute uncertainty on $\mathrm{Im}\,\widetilde{\mathcal{H}}$")
        ax.set_title(f"{mode}: Samy-Su22-calibrated polarized-target constraint")
        ax.grid(alpha=.2); ax.legend(fontsize=8)
        savefig(fig, figures / f"AUL_SamySu22_scan_{mode}_median_sigma_ImHt.png")

        # At 10x RGA, isolate the polarized-target running leverage.
        fig, ax = plt.subplots(figsize=(8.8, 6.0))
        labels=[]; vals=[]
        for scenario, n_su22 in scenarios:
            d = results[(results.fit_mode == mode) &
                        (results.observable_set == "XS+BSA+AUL") &
                        (results.luminosity_factor == 10) &
                        (results.aul_scenario == scenario)]
            labels.append(scenario); vals.append(np.nanmedian(d.sigma_ImHt))
        ax.plot(range(len(vals)), vals, marker="o")
        ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, rotation=20, ha="right")
        ax.set_ylabel(r"Median absolute uncertainty on $\mathrm{Im}\,\widetilde{\mathcal{H}}$")
        ax.set_title(f"{mode}: polarized-target leverage at 10x RGA")
        ax.grid(alpha=.2)
        savefig(fig, figures / f"AUL_SamySu22_{mode}_polarized_running_at_10xRGA.png")

    # Direct comparison of A_UL alone with A_UL+A_UT.  This is the key RGH
    # diagnostic: does transverse-target information release E and recover
    # useful Htilde precision in the six-CFF local fit?
    for mode in ("H_Ht", "H_Ht_E"):
        fig, ax = plt.subplots(figsize=(8.8, 6.0))
        for obs, lslabel in (("XS+BSA+AUL", "+ AUL"),
                             ("XS+BSA+AUL+AUT", "+ AUL + AUT (RGH)")):
            d = results[(results.fit_mode == mode) &
                        (results.observable_set == obs) &
                        (results.aul_scenario == "remaining 1x")]
            med = d.groupby("luminosity_factor").sigma_ImHt.median()
            ax.plot(med.index, med.values, marker="o", label=lslabel)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xticks(LUMI_FACTORS); ax.set_xticklabels([f"{x}x" for x in LUMI_FACTORS])
        ax.set_xlabel("Unpolarized RGA luminosity")
        ax.set_ylabel(r"Median absolute uncertainty on $\mathrm{Im}\,\widetilde{\mathcal{H}}$")
        ax.set_title(f"{mode}: impact of future transverse-target A_UT")
        ax.grid(alpha=.2); ax.legend()
        savefig(fig, figures / f"AUT_RGH_{mode}_impact_remaining1x.png")

    # At 10x RGA show the full polarized-running scan with and without RGH.
    mode = "H_Ht_E"
    fig, ax = plt.subplots(figsize=(8.8, 6.0))
    labels = [x[0] for x in scenarios]
    for obs, lab in (("XS+BSA+AUL", "+ AUL"),
                     ("XS+BSA+AUL+AUT", "+ AUL + AUT (RGH)")):
        vals=[]
        for scenario, _ in scenarios:
            d = results[(results.fit_mode == mode) &
                        (results.observable_set == obs) &
                        (results.luminosity_factor == 10) &
                        (results.aul_scenario == scenario)]
            vals.append(np.nanmedian(d.sigma_ImHt))
        ax.plot(range(len(vals)), vals, marker="o", label=lab)
    ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_ylabel(r"Median absolute uncertainty on $\mathrm{Im}\,\widetilde{\mathcal{H}}$")
    ax.set_title("H+Htilde+E: RGH transverse-target leverage at 10x RGA")
    ax.grid(alpha=.2); ax.legend()
    savefig(fig, figures / "AUT_RGH_H_Ht_E_polarized_running_at_10xRGA.png")

    make_presentation_cff_band_figure(results, figures)
    make_representative_low_t_cff_band_diagnostic(results, figures)

    print("\n" + "=" * 124)
    print("A_UL CONSTRAINT: SAMY SU22 PRECISION + CONSERVATIVE RGC RUNNING PROJECTION")
    print("=" * 124)
    print(f"Already-collected RGC effective statistics = {args.current_rgc_equiv_su22:g} x Samy Su22.")
    print("Remaining approved beam time is taken equal to the already-collected half at nominal performance.")
    print("Remaining-time factors 0.5x, 1x, 2x therefore add 0.5, 1, 2 times the current effective statistics.")
    print(f"AUL common scale nuisance = {100*args.aul_scale_frac:.2f}%.")
    print("Su22 point errors are digitized from thesis Fig. 5.10; they contain yield statistics plus statistical pi0-subtraction uncertainty.")
    rows = []
    for mode in ("H_Ht", "H_Ht_E"):
        for factor in LUMI_FACTORS:
            for scenario, n_su22 in scenarios:
                d = results[(results.fit_mode == mode) &
                            (results.observable_set == "XS+BSA+AUL") &
                            (results.luminosity_factor == factor) &
                            (results.aul_scenario == scenario)]
                dg = diagnostics[(diagnostics.fit_mode == mode) &
                                 (diagnostics.luminosity_factor == factor) &
                                 (diagnostics.include_aul == True) &
                                 (diagnostics.aul_scenario == scenario)]
                rows.append({
                    "mode": mode, "RGA_L": f"{factor}x", "RGC_scenario": scenario,
                    "N/Su22": n_su22,
                    "AUL_err/Su22": 1/math.sqrt(n_su22),
                    "med_sigma_ImHt": np.nanmedian(d.sigma_ImHt),
                    "med_rel_ImHt_%": 100*np.nanmedian(d.relative_sigma_ImHt),
                    "med_maxcorr_ImHt": np.nanmedian(d.max_abs_corr_ImHt),
                    "rank_def": int(dg.iloc[0].rank_deficit),
                    "cond": dg.iloc[0].condition_number_effective,
                })
    print(pd.DataFrame(rows).to_string(index=False, float_format=lambda x: f"{x:.4g}"))

    print("\n" + "=" * 124)
    print("RGH A_UT IMPACT: SAME MATCHED-BIN PRECISION AS THE CORRESPONDING RGC A_UL SCENARIO")
    print("=" * 124)
    print(f"AUT common scale nuisance = {100*args.aut_scale_frac:.2f}%.")
    print("RGH has not run: its A_UT precision is an explicit projection assumption, not measured input.")
    rows = []
    for mode in ("H_Ht", "H_Ht_E"):
        for factor in LUMI_FACTORS:
            for scenario, n_su22 in scenarios:
                a = results[(results.fit_mode == mode) &
                            (results.observable_set == "XS+BSA+AUL") &
                            (results.luminosity_factor == factor) &
                            (results.aul_scenario == scenario)]
                u = results[(results.fit_mode == mode) &
                            (results.observable_set == "XS+BSA+AUL+AUT") &
                            (results.luminosity_factor == factor) &
                            (results.aul_scenario == scenario)]
                dg = diagnostics[(diagnostics.fit_mode == mode) &
                                 (diagnostics.luminosity_factor == factor) &
                                 (diagnostics.include_aul == True) &
                                 (diagnostics.include_aut == True) &
                                 (diagnostics.aul_scenario == scenario)]
                sa = np.nanmedian(a.sigma_ImHt); su = np.nanmedian(u.sigma_ImHt)
                rows.append({
                    "mode": mode, "RGA_L": f"{factor}x", "pol_scenario": scenario,
                    "N/Su22_each": n_su22,
                    "sigma_ImHt_AUL": sa,
                    "sigma_ImHt_AUL_AUT": su,
                    "AUT_gain": sa/su if su > 0 else np.nan,
                    "rel_ImHt_AUL_AUT_%": 100*np.nanmedian(u.relative_sigma_ImHt),
                    "maxcorr_ImHt_AUL_AUT": np.nanmedian(u.max_abs_corr_ImHt),
                    "cond_AUL_AUT": dg.iloc[0].condition_number_effective,
                })
    print(pd.DataFrame(rows).to_string(index=False, float_format=lambda x: f"{x:.4g}"))

    print("[interpretation]")
    print("  H       : optimistic H-dominance reference.")
    print("  H_Ht    : asks whether XS+BSA can separate Htilde once H is free.")
    print("  H_Ht_E  : asks how much that conclusion survives when E is also free.")
    print("  Absolute errors remain meaningful near CFF zero crossings.")
    print("  Relative errors are suppressed below the configured truth floor.")
    print("  AUL precision is calibrated to Samy's Su22 thesis Fig. 5.10, not to RGA BSA errors.")
    print("  Current collected RGC is conservatively treated as 3x Su22 effective statistics.")
    print("  Remaining 0.5x/1x/2x scenarios give final 4.5x/6x/9x Su22 statistics.")
    print("  Polarized-target precision is independent of the RGA unpolarized luminosity factor.")
    print("  RGH AUT is projected with the same matched-bin statistical precision as AUL in each scenario.")
    print("  Gepard AUT uses a transverse target and its default dominant sine-varphi harmonic.")
    print("  Any nonzero rank deficit is a warning that the corresponding fit")
    print("  contains an exactly unresolved CFF combination and must not be")
    print("  interpreted from pseudoinverse errors alone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
