#!/usr/bin/env python3
"""
clas6_cross_check.py

Compare current CLAS12 / Hall B pass-2 DVCS cross sections against:
  - previous pass-1 CLAS12 / Hall B Fa18 result
  - CLAS6 unpolarized cross sections from Gepard

This script uses the CLAS6 comparison points listed in the previous pass-1
analysis note and makes one PNG per comparison point.

Important conventions
---------------------

1. CLAS6 data are loaded from Gepard, by default dataset id 98:

     [98] CLAS 2640 XUU 1504.02009 CLAS data base E145M1

   corresponding to H.S. Jo et al., Phys. Rev. Lett. 115, 212003 (2015).

2. Gepard zero-placeholder points are removed.

   Dataset 98 contains entries with val=0 and uncertainties=0 in some phi bins.
   These are not physical cross-section measurements and are not plotted.

3. The selected CLAS12 pass-2 bin is found automatically.

   For each CLAS6 comparison target, the script:
     - selects the matching CLAS6 kinematic group from Gepard,
     - scans all pass-2 CLAS12 bins,
     - computes weighted-average CLAS12 kinematics for each bin,
     - chooses the closest pass-2 bin in scaled (xB,Q2,|t|) distance.

4. CLAS6 is corrected to the selected CLAS12 pass-2 kinematics using KM15.

   For each CLAS6 phi point, the correction is:

     C(phi) =
       sigma_KM15(E=10.604, xB=<CLAS12>, Q2=<CLAS12>, |t|=<CLAS12>, phi)
       /
       sigma_KM15(E=5.75,  xB=<CLAS6>,  Q2=<CLAS6>,  |t|=<CLAS6>,  phi)

   and the plotted CLAS6 point is:

     sigma_CLAS6_corrected = C(phi) * sigma_CLAS6_raw

   with CLAS6 uncertainties scaled by the same factor.

   This puts CLAS6 onto the selected CLAS12/pass-2 kinematics and beam energy.

5. Ratio plots are pass-2 / corrected-CLAS6 only.

   The corrected CLAS6 points are linearly interpolated only inside the actual
   CLAS6 phi coverage. No extrapolation is performed below the lowest CLAS6 phi
   point or above the highest CLAS6 phi point.

Plot layout
-----------

Default layout:

  top-left:     pass-2 10.6 GeV, pass-1 Fa18, corrected CLAS6
  top-right:    pass-2 10.6 GeV / corrected CLAS6
  bottom-left:  individual pass-2 run periods and corrected CLAS6
  bottom-right: individual pass-2 run periods / corrected CLAS6

Simple layout, enabled with --simple:

  left:   pass-2 10.6 GeV and corrected CLAS6
  right:  pass-2 10.6 GeV / corrected CLAS6

In --simple mode:
  - pass-1 is not loaded or plotted.
  - individual pass-2 run periods are not extracted or plotted.
  - only the combined 10.6 GeV pass-2 cross-section column is required.

Uncertainty treatment
---------------------

CLAS6:
  total uncertainty is sqrt(errstat^2 + errsyst^2) when available. If no
  stat/sys decomposition exists, the Gepard err value is used. The KM15
  kinematic-transfer factor scales both central values and uncertainties.

pass-2 CLAS12:
  CSV tuple is read as:

    (value, stat, sys_csv)

  By default, an additional 15% bin-to-bin systematic is added:

    sys_est = 0.15 * sigma

  and the plotted uncertainty is:

    err_total = sqrt(stat^2 + sys_csv^2 + sys_est^2)

pass-1 CLAS12:
  CSV columns are read as:

    cross sections, ep->epg, exp
    cross sections, ep->epg, exp, stat. unc.
    cross sections, ep->epg, exp, syst. unc. (up)
    cross sections, ep->epg, exp, syst. unc. (down)

  A 31% normalization uncertainty is added in quadrature:

    norm = 0.31 * sigma

  so:

    err_up   = sqrt(stat^2 + sys_up^2 + norm^2)
    err_down = sqrt(stat^2 + sys_down^2 + norm^2)

Output
------

  output/clas6_cross_check/clas6_cross_check_<label>_p2bin<bin>.png
"""

import argparse
import copy
import math
import os
import re
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# -----------------------------------------------------------------------------
# CLAS6 comparison targets from the pass-1 analysis note.
# -----------------------------------------------------------------------------

CLAS6_EBEAM = 5.75
CLAS12_EBEAM = 10.604

CLAS6_TARGETS = [
    # label, xB, Q2, |t|
    ("26",  0.185, 1.630, 0.200),
    ("27",  0.185, 1.630, 0.260),
    ("28",  0.185, 1.630, 0.340),
    ("29",  0.185, 1.630, 0.450),
    ("52",  0.245, 2.120, 0.340),
    ("53",  0.245, 2.120, 0.450),
    ("64",  0.275, 2.350, 0.340),
    ("86",  0.335, 2.780, 0.340),
    ("108", 0.451, 3.630, 0.340),
]


# -----------------------------------------------------------------------------
# Default file locations.
# -----------------------------------------------------------------------------

DEFAULT_PASS1_CSV = (
    "/u/home/thayward/clas12_analysis_software/analysis_scripts/"
    "dvcs_cross_section/imports/all_bin_v3.csv"
)


# -----------------------------------------------------------------------------
# Defaults.
# -----------------------------------------------------------------------------

DEFAULT_PASS2_SCALE = 1.0
DEFAULT_PASS1_SCALE = 1.0
DEFAULT_CLAS6_SCALE = 1.0

DEFAULT_PASS2_BIN_TO_BIN_SYS_FRAC = 0.15
DEFAULT_PASS1_NORM_SYS_FRAC = 0.31

DEFAULT_CLAS6_DATASET_ID = 98
DEFAULT_OUTPUT_DIR = "output/clas6_cross_check"


# -----------------------------------------------------------------------------
# Plot labels and styles.
# -----------------------------------------------------------------------------

PASS2_COMBINED_DISPLAY_LABEL = "pass-2 10.6 GeV"
PASS2_COMBINED_CSV_PERIOD = "10.6 GeV"

PASS1_DISPLAY_LABEL = "pass-1 Fa18"
CLAS6_DISPLAY_LABEL = "CLAS6 corrected to pass-2"

PASS2_PERIOD_DISPLAY_TO_CSV_PERIOD = {
    "Sp18 Inb": "Sp18 Inb",
    "Sp18 Out": "Sp18 Out",
    "Fa18 Inb": "Fa18 Inb",
    "Fa18 Out": "Fa18 Out",
}

TOP_CROSS_SECTION_SERIES = [
    PASS2_COMBINED_DISPLAY_LABEL,
    PASS1_DISPLAY_LABEL,
]

TOP_RATIO_SERIES = [
    PASS2_COMBINED_DISPLAY_LABEL,
]

BOTTOM_CROSS_SECTION_SERIES = [
    "Sp18 Inb",
    "Sp18 Out",
    "Fa18 Inb",
    "Fa18 Out",
]

BOTTOM_RATIO_SERIES = [
    "Sp18 Inb",
    "Sp18 Out",
    "Fa18 Inb",
    "Fa18 Out",
]

ALL_PASS2_CSV_PERIODS = [
    PASS2_COMBINED_CSV_PERIOD,
    "Sp18 Inb",
    "Sp18 Out",
    "Fa18 Inb",
    "Fa18 Out",
]

SIMPLE_PASS2_CSV_PERIODS = [
    PASS2_COMBINED_CSV_PERIOD,
]

SERIES_STYLE = {
    CLAS6_DISPLAY_LABEL: dict(marker="s", linestyle="None", color="black"),
    PASS2_COMBINED_DISPLAY_LABEL: dict(marker="o", linestyle="None", color="tab:red"),
    PASS1_DISPLAY_LABEL: dict(marker="o", linestyle="None", color="tab:cyan"),
    "Sp18 Inb": dict(marker="D", linestyle="None", color="tab:green"),
    "Sp18 Out": dict(marker="P", linestyle="None", color="tab:purple"),
    "Fa18 Inb": dict(marker="^", linestyle="None", color="tab:blue"),
    "Fa18 Out": dict(marker="v", linestyle="None", color="tab:orange"),
}

FLOAT_PATTERN = re.compile(
    r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
)


@dataclass
class DataPoint:
    phi: float
    sigma: float
    err_low: float
    err_high: float
    stat: float = 0.0
    sys: float = 0.0
    sys_csv: float = 0.0
    sys_est: float = 0.0
    norm: float = 0.0
    model_corr: float = 1.0


@dataclass
class KinematicPoint:
    xB: float
    Q2: float
    t_abs: float


@dataclass
class KinematicAverages:
    xB: float
    Q2: float
    t_abs: float
    weight_period: str
    n_points: int


@dataclass
class TargetSpec:
    label: str
    xB: float
    Q2: float
    t_abs: float


@dataclass
class BinChoice:
    bin_name: int
    xB: float
    Q2: float
    t_abs: float
    distance: float
    n_rows: int


# -----------------------------------------------------------------------------
# Generic helpers.
# -----------------------------------------------------------------------------

def parse_tuple3(value) -> Tuple[float, float, float]:
    if value is None:
        return (math.nan, math.nan, math.nan)
    # endif

    if isinstance(value, float) and math.isnan(value):
        return (math.nan, math.nan, math.nan)
    # endif

    if isinstance(value, (int, float, np.integer, np.floating)):
        return (float(value), 0.0, 0.0)
    # endif

    text = str(value).strip()

    if text == "" or text.lower() in {"nan", "none", "null"}:
        return (math.nan, math.nan, math.nan)
    # endif

    numbers = FLOAT_PATTERN.findall(text)

    if len(numbers) == 0:
        return (math.nan, math.nan, math.nan)
    # endif

    parsed = [float(x) for x in numbers]

    while len(parsed) < 3:
        parsed.append(0.0)
    # endwhile

    return (parsed[0], parsed[1], parsed[2])


def parse_scalar_from_cell(value) -> float:
    if value is None:
        return math.nan
    # endif

    if isinstance(value, float) and math.isnan(value):
        return math.nan
    # endif

    if isinstance(value, (int, float, np.integer, np.floating)):
        return float(value)
    # endif

    text = str(value).strip()

    if text == "" or text.lower() in {"nan", "none", "null"}:
        return math.nan
    # endif

    numbers = FLOAT_PATTERN.findall(text)

    if len(numbers) == 0:
        return math.nan
    # endif

    return float(numbers[0])


def finite_or_zero(value: float) -> float:
    if np.isfinite(value):
        return float(value)
    # endif

    return 0.0


def is_zero_placeholder(*values: float, tolerance: float = 1.0e-15) -> bool:
    finite_values = []

    for value in values:
        try:
            value_float = float(value)
        except Exception:
            continue
        # endtry

        if np.isfinite(value_float):
            finite_values.append(abs(value_float))
        # endif
    # endfor

    if len(finite_values) == 0:
        return True
    # endif

    return all(value <= tolerance for value in finite_values)


def symmetric_error(point: DataPoint) -> float:
    return 0.5 * (point.err_low + point.err_high)


def points_to_arrays(points: List[DataPoint]) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    phi = np.array([p.phi for p in points], dtype=float)
    sigma = np.array([p.sigma for p in points], dtype=float)
    err_low = np.array([p.err_low for p in points], dtype=float)
    err_high = np.array([p.err_high for p in points], dtype=float)

    return phi, sigma, err_low, err_high


def require_columns(df: pd.DataFrame, columns: Iterable[str], context: str) -> None:
    missing = [c for c in columns if c not in df.columns]

    if missing:
        raise KeyError(
            f"The {context} CSV is missing required columns:\n  "
            + "\n  ".join(missing)
        )
    # endif


# -----------------------------------------------------------------------------
# CSV column helpers.
# -----------------------------------------------------------------------------

def pass2_cross_section_column(csv_period: str) -> str:
    return f"normed cross sections, ep->epg, exp, {csv_period}, unpol"


def average_phi_column(csv_period: str) -> str:
    return f"phiavg, {csv_period}"


def avg_column(quantity: str, csv_period: str) -> str:
    return f"{quantity}, {csv_period}"


# -----------------------------------------------------------------------------
# Gepard / KM15 helpers.
# -----------------------------------------------------------------------------

def import_gepard():
    try:
        import gepard as g
    except ImportError as exc:
        raise RuntimeError(
            "Could not import the Gepard package.\n\n"
            "Install it into the Python interpreter used to run this script.\n"
            "On ifarm, prefer:\n"
            "  python -m pip install --user gepard\n"
        ) from exc
    # endtry

    return g


def import_km15():
    try:
        from gepard.fits import th_KM15
    except ImportError as exc:
        raise RuntimeError("Could not import th_KM15 from gepard.fits.") from exc
    # endtry

    return th_KM15


def get_attr_any(obj, names: List[str], default=math.nan):
    for name in names:
        if hasattr(obj, name):
            return getattr(obj, name)
        # endif
    # endfor

    return default


def make_gepard_xs_point(
    g,
    xB: float,
    Q2: float,
    t_abs: float,
    phi_deg_trento: float,
    ebeam: float,
):
    pt = g.DataPoint(
        xB=float(xB),
        Q2=float(Q2),
        t=-abs(float(t_abs)),
        phi=math.radians(float(phi_deg_trento)),
        frame="Trento",
        process="ep2epgamma",
        exptype="fixed target",
        in1energy=float(ebeam),
        in1charge=-1,
        in1polarization=0,
        observable="XS",
    )

    if hasattr(pt, "to_conventions"):
        pt.to_conventions()
    # endif

    return pt


def km15_xs(
    g,
    th_km15,
    xB: float,
    Q2: float,
    t_abs: float,
    phi_deg_trento: float,
    ebeam: float,
) -> float:
    pt = make_gepard_xs_point(
        g=g,
        xB=xB,
        Q2=Q2,
        t_abs=t_abs,
        phi_deg_trento=phi_deg_trento,
        ebeam=ebeam,
    )

    return float(th_km15.predict(pt))


def km15_transfer_factor_clas6_to_clas12(
    g,
    th_km15,
    phi_deg_trento: float,
    clas6_kin: KinematicPoint,
    clas6_ebeam: float,
    clas12_kin: KinematicPoint,
    clas12_ebeam: float,
) -> float:
    """
    Transfer CLAS6 data to CLAS12/pass-2 kinematics:

      C(phi) =
        KM15(CLAS12 kinematics, 10.604 GeV)
        /
        KM15(CLAS6 kinematics, 5.75 GeV)
    """

    numerator = km15_xs(
        g=g,
        th_km15=th_km15,
        xB=clas12_kin.xB,
        Q2=clas12_kin.Q2,
        t_abs=clas12_kin.t_abs,
        phi_deg_trento=phi_deg_trento,
        ebeam=clas12_ebeam,
    )

    denominator = km15_xs(
        g=g,
        th_km15=th_km15,
        xB=clas6_kin.xB,
        Q2=clas6_kin.Q2,
        t_abs=clas6_kin.t_abs,
        phi_deg_trento=phi_deg_trento,
        ebeam=clas6_ebeam,
    )

    if not np.isfinite(numerator) or not np.isfinite(denominator) or denominator == 0.0:
        return math.nan
    # endif

    return numerator / denominator


def apply_model_correction_to_clas6(
    raw_points: List[DataPoint],
    g,
    th_km15,
    clas6_kin: KinematicPoint,
    clas6_ebeam: float,
    clas12_kin: KinematicPoint,
    clas12_ebeam: float,
) -> List[DataPoint]:
    corrected_points: List[DataPoint] = []

    for point in raw_points:
        corr = km15_transfer_factor_clas6_to_clas12(
            g=g,
            th_km15=th_km15,
            phi_deg_trento=point.phi,
            clas6_kin=clas6_kin,
            clas6_ebeam=clas6_ebeam,
            clas12_kin=clas12_kin,
            clas12_ebeam=clas12_ebeam,
        )

        if not np.isfinite(corr):
            continue
        # endif

        corrected = copy.copy(point)

        corrected.sigma *= corr
        corrected.err_low *= abs(corr)
        corrected.err_high *= abs(corr)
        corrected.stat *= abs(corr)
        corrected.sys *= abs(corr)
        corrected.sys_csv *= abs(corr)
        corrected.sys_est *= abs(corr)
        corrected.norm *= abs(corr)
        corrected.model_corr = corr

        corrected_points.append(corrected)
    # endfor

    corrected_points.sort(key=lambda p: p.phi)

    return corrected_points


# -----------------------------------------------------------------------------
# Gepard / CLAS6 data.
# -----------------------------------------------------------------------------

def convert_gepard_phi_to_degrees(pt, phi_convention: str) -> float:
    """
    Convert Gepard point phi to degrees.

    Gepard stores available datasets internally in BMK convention. For datasets
    originally given in Trento, converting back to the published-style angle is:

      phi_Trento = pi - phi_BMK

    modulo 2pi.
    """

    phi_raw = get_attr_any(pt, ["phi"])
    phi = float(phi_raw)

    is_radians = abs(phi) <= 2.0 * math.pi * 1.5
    frame = str(get_attr_any(pt, ["frame"], default="")).lower()

    if phi_convention == "trento" and is_radians and frame == "trento":
        phi = (math.pi - phi) % (2.0 * math.pi)
    # endif

    if is_radians:
        phi_deg = math.degrees(phi)
    else:
        phi_deg = phi
    # endif

    return phi_deg % 360.0


def clas6_point_from_gepard(
    pt,
    clas6_scale: float,
    phi_convention: str,
) -> Optional[Tuple[KinematicPoint, DataPoint]]:
    xB = float(get_attr_any(pt, ["xB", "xb"]))
    Q2 = float(get_attr_any(pt, ["Q2", "q2"]))
    t = float(get_attr_any(pt, ["t"]))
    val = float(get_attr_any(pt, ["val", "value"]))

    if not (np.isfinite(xB) and np.isfinite(Q2) and np.isfinite(t) and np.isfinite(val)):
        return None
    # endif

    stat = get_attr_any(pt, ["errstat", "err_stat", "stat"], default=math.nan)
    sys = get_attr_any(pt, ["errsyst", "errsys", "err_syst", "sys"], default=math.nan)
    err = get_attr_any(pt, ["err"], default=math.nan)

    stat = float(stat) if np.isfinite(stat) else math.nan
    sys = float(sys) if np.isfinite(sys) else math.nan
    err = float(err) if np.isfinite(err) else math.nan

    if is_zero_placeholder(val, stat, sys, err):
        return None
    # endif

    if np.isfinite(stat) or np.isfinite(sys):
        stat = finite_or_zero(stat)
        sys = finite_or_zero(sys)
        total = math.sqrt(stat * stat + sys * sys)
    elif np.isfinite(err):
        stat = float(err)
        sys = 0.0
        total = float(err)
    else:
        stat = 0.0
        sys = 0.0
        total = 0.0
    # endif

    val *= clas6_scale
    stat *= clas6_scale
    sys *= clas6_scale
    total *= clas6_scale

    phi_deg = convert_gepard_phi_to_degrees(pt, phi_convention=phi_convention)

    kin = KinematicPoint(
        xB=xB,
        Q2=Q2,
        t_abs=abs(t),
    )

    data = DataPoint(
        phi=phi_deg,
        sigma=val,
        err_low=total,
        err_high=total,
        stat=stat,
        sys=sys,
        sys_csv=sys,
        sys_est=0.0,
        norm=0.0,
        model_corr=1.0,
    )

    return kin, data


def load_all_clas6_points(
    g,
    dataset_id: int,
    clas6_scale: float,
    phi_convention: str,
) -> List[Tuple[KinematicPoint, DataPoint]]:
    if dataset_id not in g.dset:
        raise RuntimeError(f"Gepard dataset id {dataset_id} is not available in g.dset.")
    # endif

    dataset = g.dset[dataset_id]
    points: List[Tuple[KinematicPoint, DataPoint]] = []

    for pt in dataset:
        converted = clas6_point_from_gepard(
            pt=pt,
            clas6_scale=clas6_scale,
            phi_convention=phi_convention,
        )

        if converted is None:
            continue
        # endif

        points.append(converted)
    # endfor

    return points


def group_clas6_points(
    clas6_points: List[Tuple[KinematicPoint, DataPoint]],
    round_digits: int,
) -> Dict[Tuple[float, float, float], List[Tuple[KinematicPoint, DataPoint]]]:
    groups: Dict[Tuple[float, float, float], List[Tuple[KinematicPoint, DataPoint]]] = {}

    for kin, point in clas6_points:
        key = (
            round(kin.xB, round_digits),
            round(kin.Q2, round_digits),
            round(kin.t_abs, round_digits),
        )

        if key not in groups:
            groups[key] = []
        # endif

        groups[key].append((kin, point))
    # endfor

    return groups


def mean_kinematics_for_group(group: List[Tuple[KinematicPoint, DataPoint]]) -> KinematicPoint:
    return KinematicPoint(
        xB=float(np.mean([item[0].xB for item in group])),
        Q2=float(np.mean([item[0].Q2 for item in group])),
        t_abs=float(np.mean([item[0].t_abs for item in group])),
    )


def scaled_kinematic_distance(
    a: KinematicPoint,
    b: KinematicPoint,
    xB_scale: float,
    Q2_scale: float,
    t_scale: float,
) -> float:
    return math.sqrt(
        ((a.xB - b.xB) / xB_scale) ** 2
        + ((a.Q2 - b.Q2) / Q2_scale) ** 2
        + ((a.t_abs - b.t_abs) / t_scale) ** 2
    )


def select_nearest_clas6_group(
    groups: Dict[Tuple[float, float, float], List[Tuple[KinematicPoint, DataPoint]]],
    target: TargetSpec,
    xB_scale: float,
    Q2_scale: float,
    t_scale: float,
) -> Tuple[KinematicPoint, List[DataPoint], float]:
    target_kin = KinematicPoint(xB=target.xB, Q2=target.Q2, t_abs=target.t_abs)

    best_key = None
    best_distance = math.inf
    best_kin = None

    for key, group in groups.items():
        kin = mean_kinematics_for_group(group)
        distance = scaled_kinematic_distance(
            a=kin,
            b=target_kin,
            xB_scale=xB_scale,
            Q2_scale=Q2_scale,
            t_scale=t_scale,
        )

        if distance < best_distance:
            best_key = key
            best_distance = distance
            best_kin = kin
        # endif
    # endfor

    if best_key is None or best_kin is None:
        raise RuntimeError("No CLAS6 groups were available after loading Gepard data.")
    # endif

    points = [item[1] for item in groups[best_key]]
    points.sort(key=lambda p: p.phi)

    return best_kin, points, best_distance


# -----------------------------------------------------------------------------
# CLAS12 pass-2 bin matching.
# -----------------------------------------------------------------------------

def compute_weighted_kinematic_averages(
    selected: pd.DataFrame,
    csv_period: str = "10.6 GeV",
    fallback_periods: Optional[List[str]] = None,
) -> KinematicAverages:
    if fallback_periods is None:
        fallback_periods = [
            "10.6 GeV",
            "Fa18 Inb",
            "Fa18 Out",
            "Sp18 Inb",
            "Sp18 Out",
        ]
    # endif

    candidate_periods = [csv_period]

    for p in fallback_periods:
        if p not in candidate_periods:
            candidate_periods.append(p)
        # endif
    # endfor

    chosen_period = csv_period

    for p in candidate_periods:
        required_avg_cols = [
            avg_column("xBavg", p),
            avg_column("Q2avg", p),
            avg_column("t_abs_avg", p),
        ]

        if all(col in selected.columns for col in required_avg_cols):
            chosen_period = p
            break
        # endif
    # endfor

    xb_col = avg_column("xBavg", chosen_period)
    q2_col = avg_column("Q2avg", chosen_period)
    t_col = avg_column("t_abs_avg", chosen_period)
    xs_col = pass2_cross_section_column(chosen_period)

    def fallback_midpoint(min_col: str, max_col: str) -> float:
        lo = pd.to_numeric(selected[min_col], errors="coerce")
        hi = pd.to_numeric(selected[max_col], errors="coerce")

        if len(lo) > 0 and np.isfinite(lo.iloc[0]) and np.isfinite(hi.iloc[0]):
            return 0.5 * (float(lo.iloc[0]) + float(hi.iloc[0]))
        # endif

        return math.nan

    if xb_col not in selected.columns:
        return KinematicAverages(
            xB=fallback_midpoint("xBmin", "xBmax"),
            Q2=fallback_midpoint("Q2min", "Q2max"),
            t_abs=fallback_midpoint("t_abs_min", "t_abs_max"),
            weight_period="bin midpoint fallback",
            n_points=0,
        )
    # endif

    weighted_sums = {"xB": 0.0, "Q2": 0.0, "t_abs": 0.0}
    unweighted_values = {"xB": [], "Q2": [], "t_abs": []}

    weight_sum = 0.0
    n_points = 0

    for _, row in selected.iterrows():
        xB = parse_scalar_from_cell(row[xb_col])
        Q2 = parse_scalar_from_cell(row[q2_col])
        t_abs = parse_scalar_from_cell(row[t_col])

        if not (np.isfinite(xB) and np.isfinite(Q2) and np.isfinite(t_abs)):
            continue
        # endif

        unweighted_values["xB"].append(xB)
        unweighted_values["Q2"].append(Q2)
        unweighted_values["t_abs"].append(t_abs)

        weight = 1.0

        if xs_col in selected.columns:
            _, stat, _ = parse_tuple3(row[xs_col])

            if np.isfinite(stat) and stat > 0.0:
                weight = 1.0 / (stat * stat)
            else:
                weight = 1.0
            # endif
        # endif

        weighted_sums["xB"] += weight * xB
        weighted_sums["Q2"] += weight * Q2
        weighted_sums["t_abs"] += weight * t_abs
        weight_sum += weight
        n_points += 1
    # endfor

    if weight_sum > 0.0:
        return KinematicAverages(
            xB=weighted_sums["xB"] / weight_sum,
            Q2=weighted_sums["Q2"] / weight_sum,
            t_abs=weighted_sums["t_abs"] / weight_sum,
            weight_period=chosen_period,
            n_points=n_points,
        )
    # endif

    if len(unweighted_values["xB"]) > 0:
        return KinematicAverages(
            xB=float(np.mean(unweighted_values["xB"])),
            Q2=float(np.mean(unweighted_values["Q2"])),
            t_abs=float(np.mean(unweighted_values["t_abs"])),
            weight_period=f"{chosen_period}, unweighted fallback",
            n_points=len(unweighted_values["xB"]),
        )
    # endif

    return KinematicAverages(
        xB=fallback_midpoint("xBmin", "xBmax"),
        Q2=fallback_midpoint("Q2min", "Q2max"),
        t_abs=fallback_midpoint("t_abs_min", "t_abs_max"),
        weight_period="bin midpoint fallback",
        n_points=0,
    )


def select_pass2_bin_nearest_to_kinematics(
    pass2_df: pd.DataFrame,
    target_kin: KinematicPoint,
    xB_scale: float,
    Q2_scale: float,
    t_scale: float,
    avg_period: str,
) -> Tuple[pd.DataFrame, BinChoice]:
    best_choice: Optional[BinChoice] = None
    best_rows: Optional[pd.DataFrame] = None

    bin_names = sorted(
        int(x)
        for x in pd.to_numeric(pass2_df["Bin Name"], errors="coerce").dropna().unique()
    )

    bin_numeric_all = pd.to_numeric(pass2_df["Bin Name"], errors="coerce")

    for bin_name in bin_names:
        rows = pass2_df.loc[bin_numeric_all == bin_name].copy()

        if len(rows) == 0:
            continue
        # endif

        avg = compute_weighted_kinematic_averages(
            selected=rows,
            csv_period=avg_period,
        )

        if not (np.isfinite(avg.xB) and np.isfinite(avg.Q2) and np.isfinite(avg.t_abs)):
            continue
        # endif

        kin = KinematicPoint(xB=avg.xB, Q2=avg.Q2, t_abs=avg.t_abs)

        distance = scaled_kinematic_distance(
            a=kin,
            b=target_kin,
            xB_scale=xB_scale,
            Q2_scale=Q2_scale,
            t_scale=t_scale,
        )

        if best_choice is None or distance < best_choice.distance:
            best_choice = BinChoice(
                bin_name=bin_name,
                xB=avg.xB,
                Q2=avg.Q2,
                t_abs=avg.t_abs,
                distance=distance,
                n_rows=len(rows),
            )
            best_rows = rows.sort_values(["phimin", "phimax"]).copy()
        # endif
    # endfor

    if best_choice is None or best_rows is None:
        raise RuntimeError("Could not select a nearest pass-2 CLAS12 bin.")
    # endif

    return best_rows, best_choice


def select_bin_by_name(df: pd.DataFrame, bin_name: int) -> pd.DataFrame:
    selected = df.loc[pd.to_numeric(df["Bin Name"], errors="coerce") == int(bin_name)].copy()
    selected = selected.sort_values(["phimin", "phimax"]).copy()
    return selected


# -----------------------------------------------------------------------------
# Pass-2 and pass-1 extraction.
# -----------------------------------------------------------------------------

def pass2_points_for_period(
    selected: pd.DataFrame,
    csv_period: str,
    pass2_scale: float,
    include_pass2_estimated_sys: bool,
    pass2_bin_to_bin_sys_frac: float,
) -> List[DataPoint]:
    xs_col = pass2_cross_section_column(csv_period)
    phi_avg_col = average_phi_column(csv_period)

    if xs_col not in selected.columns:
        return []
    # endif

    points: List[DataPoint] = []

    for _, row in selected.sort_values("phimin").iterrows():
        sigma, stat, sys_csv = parse_tuple3(row[xs_col])

        if not np.isfinite(sigma):
            continue
        # endif

        stat = finite_or_zero(stat)
        sys_csv = finite_or_zero(sys_csv)

        if is_zero_placeholder(sigma, stat, sys_csv):
            continue
        # endif

        sigma *= pass2_scale
        stat *= pass2_scale
        sys_csv *= pass2_scale

        if include_pass2_estimated_sys:
            sys_est = pass2_bin_to_bin_sys_frac * abs(sigma)
        else:
            sys_est = 0.0
        # endif

        sys_total = math.sqrt(sys_csv * sys_csv + sys_est * sys_est)

        if phi_avg_col in selected.columns:
            phi = parse_scalar_from_cell(row[phi_avg_col])
        else:
            phi = math.nan
        # endif

        if not np.isfinite(phi):
            phi_min = parse_scalar_from_cell(row["phimin"])
            phi_max = parse_scalar_from_cell(row["phimax"])
            phi = 0.5 * (phi_min + phi_max)
        # endif

        err_total = math.sqrt(stat * stat + sys_total * sys_total)

        points.append(
            DataPoint(
                phi=phi,
                sigma=sigma,
                err_low=err_total,
                err_high=err_total,
                stat=stat,
                sys=sys_total,
                sys_csv=sys_csv,
                sys_est=sys_est,
                norm=0.0,
                model_corr=1.0,
            )
        )
    # endfor

    points.sort(key=lambda p: p.phi)
    return points


def pass1_points(
    selected: pd.DataFrame,
    pass1_scale: float,
    pass1_norm_sys_frac: float,
) -> List[DataPoint]:
    required = [
        "cross sections, ep->epg, exp",
        "cross sections, ep->epg, exp, stat. unc.",
        "cross sections, ep->epg, exp, syst. unc. (up)",
        "cross sections, ep->epg, exp, syst. unc. (down)",
        "phiavg",
        "phimin",
        "phimax",
    ]

    for col in required:
        if col not in selected.columns:
            return []
        # endif
    # endfor

    points: List[DataPoint] = []

    for _, row in selected.sort_values("phimin").iterrows():
        sigma = parse_scalar_from_cell(row["cross sections, ep->epg, exp"])
        stat = parse_scalar_from_cell(row["cross sections, ep->epg, exp, stat. unc."])
        sys_up = parse_scalar_from_cell(row["cross sections, ep->epg, exp, syst. unc. (up)"])
        sys_down = parse_scalar_from_cell(row["cross sections, ep->epg, exp, syst. unc. (down)"])

        if not np.isfinite(sigma):
            continue
        # endif

        stat = finite_or_zero(stat)
        sys_up = finite_or_zero(sys_up)
        sys_down = finite_or_zero(sys_down)

        if is_zero_placeholder(sigma, stat, sys_up, sys_down):
            continue
        # endif

        sigma *= pass1_scale
        stat *= pass1_scale
        sys_up *= pass1_scale
        sys_down *= pass1_scale

        norm = pass1_norm_sys_frac * abs(sigma)

        err_high = math.sqrt(stat * stat + sys_up * sys_up + norm * norm)
        err_low = math.sqrt(stat * stat + sys_down * sys_down + norm * norm)

        phi = parse_scalar_from_cell(row["phiavg"])

        if not np.isfinite(phi):
            phi_min = parse_scalar_from_cell(row["phimin"])
            phi_max = parse_scalar_from_cell(row["phimax"])
            phi = 0.5 * (phi_min + phi_max)
        # endif

        points.append(
            DataPoint(
                phi=phi,
                sigma=sigma,
                err_low=err_low,
                err_high=err_high,
                stat=stat,
                sys=0.5 * (abs(sys_up) + abs(sys_down)),
                sys_csv=0.5 * (abs(sys_up) + abs(sys_down)),
                sys_est=0.0,
                norm=norm,
                model_corr=1.0,
            )
        )
    # endfor

    points.sort(key=lambda p: p.phi)
    return points


# -----------------------------------------------------------------------------
# Non-extrapolating interpolation, ratios, chi2.
# -----------------------------------------------------------------------------

def bounded_linear_interp(
    phi_query: float,
    phi_values: np.ndarray,
    y_values: np.ndarray,
) -> float:
    """
    Linear interpolation only inside the actual reference phi range.

    This intentionally does not use periodic extension. It prevents ratio points
    from being made where CLAS6 has no actual nearby data and the old script
    would have extrapolated across the 0/360 boundary.
    """

    phi_base = np.array(phi_values, dtype=float)
    y_base = np.array(y_values, dtype=float)

    if len(phi_base) < 2:
        return math.nan
    # endif

    order = np.argsort(phi_base)
    phi_base = phi_base[order]
    y_base = y_base[order]

    phi_min = float(np.nanmin(phi_base))
    phi_max = float(np.nanmax(phi_base))

    if phi_query < phi_min or phi_query > phi_max:
        return math.nan
    # endif

    return float(np.interp(phi_query, phi_base, y_base))


def interpolate_reference_point_bounded(
    phi_query: float,
    reference_points: List[DataPoint],
) -> DataPoint:
    phi, sigma, err_low, err_high = points_to_arrays(reference_points)

    return DataPoint(
        phi=phi_query,
        sigma=bounded_linear_interp(phi_query, phi, sigma),
        err_low=bounded_linear_interp(phi_query, phi, err_low),
        err_high=bounded_linear_interp(phi_query, phi, err_high),
    )


def ratio_to_reference(
    numerator_points: List[DataPoint],
    reference_points: List[DataPoint],
) -> List[DataPoint]:
    ratios: List[DataPoint] = []

    for num in numerator_points:
        ref = interpolate_reference_point_bounded(
            phi_query=num.phi,
            reference_points=reference_points,
        )

        if (
            not np.isfinite(ref.sigma)
            or not np.isfinite(ref.err_low)
            or not np.isfinite(ref.err_high)
            or ref.sigma == 0.0
        ):
            continue
        # endif

        ratio = num.sigma / ref.sigma

        ratio_err_high = math.sqrt(
            (num.err_high / ref.sigma) ** 2
            + (num.sigma * ref.err_low / (ref.sigma * ref.sigma)) ** 2
        )

        ratio_err_low = math.sqrt(
            (num.err_low / ref.sigma) ** 2
            + (num.sigma * ref.err_high / (ref.sigma * ref.sigma)) ** 2
        )

        ratios.append(
            DataPoint(
                phi=num.phi,
                sigma=ratio,
                err_low=ratio_err_low,
                err_high=ratio_err_high,
            )
        )
    # endfor

    ratios.sort(key=lambda p: p.phi)
    return ratios


def chi2_ndf_to_reference(
    comparison_points: List[DataPoint],
    reference_points: List[DataPoint],
) -> Tuple[float, int, float]:
    chi2 = 0.0
    ndf = 0

    for comp in comparison_points:
        ref = interpolate_reference_point_bounded(
            phi_query=comp.phi,
            reference_points=reference_points,
        )

        if not np.isfinite(ref.sigma):
            continue
        # endif

        err_comp = symmetric_error(comp)
        err_ref = symmetric_error(ref)

        variance = err_comp * err_comp + err_ref * err_ref

        if variance <= 0.0 or not np.isfinite(variance):
            continue
        # endif

        residual = comp.sigma - ref.sigma
        chi2 += residual * residual / variance
        ndf += 1
    # endfor

    if ndf <= 0:
        return (math.nan, 0, math.nan)
    # endif

    return (chi2, ndf, chi2 / ndf)


def format_label_with_chi2(label: str, chi2_info: Tuple[float, int, float]) -> str:
    chi2, ndf, chi2ndf = chi2_info

    if ndf <= 0 or not np.isfinite(chi2ndf):
        return f"{label} ($\\chi^2$/ndf=N/A)"
    # endif

    return f"{label} ($\\chi^2$/ndf={chi2ndf:.2f})"


# -----------------------------------------------------------------------------
# Plotting.
# -----------------------------------------------------------------------------

def plot_dataset(
    ax,
    points: List[DataPoint],
    label: str,
    legend_label: Optional[str] = None,
) -> None:
    if len(points) == 0:
        return
    # endif

    phi, sigma, err_low, err_high = points_to_arrays(points)
    style = SERIES_STYLE.get(label, dict(marker="o", linestyle="None"))

    if legend_label is None:
        legend_label = label
    # endif

    if label == PASS1_DISPLAY_LABEL:
        ax.errorbar(
            phi,
            sigma,
            yerr=np.vstack([err_low, err_high]),
            label=legend_label,
            markersize=6.0,
            linewidth=0.0,
            elinewidth=1.0,
            capsize=2.5,
            marker=style["marker"],
            linestyle=style["linestyle"],
            color=style["color"],
            markerfacecolor="none",
            markeredgecolor=style["color"],
            markeredgewidth=1.3,
        )
    else:
        ax.errorbar(
            phi,
            sigma,
            yerr=np.vstack([err_low, err_high]),
            label=legend_label,
            markersize=5.5,
            linewidth=0.0,
            elinewidth=1.0,
            capsize=2.5,
            **style,
        )
    # endif


def auto_ratio_ylim(ax, points_by_label: Dict[str, List[DataPoint]]) -> None:
    yvals: List[float] = []
    yerr_lows: List[float] = []
    yerr_highs: List[float] = []

    for points in points_by_label.values():
        for p in points:
            if np.isfinite(p.sigma):
                yvals.append(p.sigma)
                yerr_lows.append(p.err_low if np.isfinite(p.err_low) else 0.0)
                yerr_highs.append(p.err_high if np.isfinite(p.err_high) else 0.0)
            # endif
        # endfor
    # endfor

    if len(yvals) == 0:
        ax.set_ylim(0.0, 2.0)
        return
    # endif

    lows = np.array(yvals) - np.array(yerr_lows)
    highs = np.array(yvals) + np.array(yerr_highs)

    ymin = float(np.nanmin(lows))
    ymax = float(np.nanmax(highs))

    span_low = abs(1.0 - ymin)
    span_high = abs(ymax - 1.0)
    span = max(span_low, span_high, 0.35)

    ax.set_ylim(1.0 - 1.10 * span, 1.0 + 1.10 * span)


def figure_suptitle_text(
    target: TargetSpec,
    raw_clas6_kin: KinematicPoint,
    pass2_choice: BinChoice,
) -> str:
    return (
        f"CLAS6 / CLAS12 DVCS cross-section cross-check: note label {target.label}\n"
        rf"CLAS6 raw: $\langle x_B\rangle={raw_clas6_kin.xB:.3f}$, "
        rf"$\langle Q^2\rangle={raw_clas6_kin.Q2:.3f}~{{\rm GeV}}^2$, "
        rf"$\langle |t|\rangle={raw_clas6_kin.t_abs:.3f}~{{\rm GeV}}^2$, "
        rf"$E={CLAS6_EBEAM:.2f}~{{\rm GeV}}$"
        "\n"
        rf"pass-2 bin {pass2_choice.bin_name}: "
        rf"$\langle x_B\rangle={pass2_choice.xB:.3f}$, "
        rf"$\langle Q^2\rangle={pass2_choice.Q2:.3f}~{{\rm GeV}}^2$, "
        rf"$\langle |t|\rangle={pass2_choice.t_abs:.3f}~{{\rm GeV}}^2$, "
        rf"$E={CLAS12_EBEAM:.3f}~{{\rm GeV}}$"
        "\n"
        r"CLAS6 points corrected to selected pass-2 kinematics with KM15"
    )


def output_path_for_target(
    output_dir: str,
    target: TargetSpec,
    pass2_choice: BinChoice,
    simple: bool,
) -> str:
    if simple:
        filename = f"clas6_cross_check_{target.label}_p2bin{pass2_choice.bin_name}_simple.png"
    else:
        filename = f"clas6_cross_check_{target.label}_p2bin{pass2_choice.bin_name}.png"
    # endif

    return os.path.join(output_dir, filename)


def make_simple_plot(
    target: TargetSpec,
    raw_clas6_kin: KinematicPoint,
    corrected_clas6_points: List[DataPoint],
    pass2_choice: BinChoice,
    points_by_label: Dict[str, List[DataPoint]],
    ratios_by_label: Dict[str, List[DataPoint]],
    chi2_by_label: Dict[str, Tuple[float, int, float]],
    output_dir: str,
    y_units: str,
    use_log_cross_section: bool,
) -> None:
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(15.0, 5.6),
        constrained_layout=True,
        sharex=True,
    )

    fig.suptitle(
        figure_suptitle_text(
            target=target,
            raw_clas6_kin=raw_clas6_kin,
            pass2_choice=pass2_choice,
        ),
        fontsize=12,
    )

    left = axes[0]
    right = axes[1]

    plot_dataset(
        ax=left,
        points=corrected_clas6_points,
        label=CLAS6_DISPLAY_LABEL,
        legend_label=CLAS6_DISPLAY_LABEL,
    )

    plot_dataset(
        ax=left,
        points=points_by_label.get(PASS2_COMBINED_DISPLAY_LABEL, []),
        label=PASS2_COMBINED_DISPLAY_LABEL,
        legend_label=format_label_with_chi2(
            PASS2_COMBINED_DISPLAY_LABEL,
            chi2_by_label.get(PASS2_COMBINED_DISPLAY_LABEL, (math.nan, 0, math.nan)),
        ),
    )

    left.set_xlabel(r"$\phi$ [deg]")
    left.set_ylabel(rf"$d^4\sigma/(dx_B\,dQ^2\,d|t|\,d\phi)$ [{y_units}]")
    left.set_title("Cross sections: pass-2 combined and corrected CLAS6")
    left.grid(True, alpha=0.25)
    left.legend(fontsize=8, frameon=True)

    if use_log_cross_section:
        left.set_yscale("log")
    # endif

    ratios = ratios_by_label.get(PASS2_COMBINED_DISPLAY_LABEL, [])

    plot_dataset(
        ax=right,
        points=ratios,
        label=PASS2_COMBINED_DISPLAY_LABEL,
        legend_label=format_label_with_chi2(
            PASS2_COMBINED_DISPLAY_LABEL,
            chi2_by_label.get(PASS2_COMBINED_DISPLAY_LABEL, (math.nan, 0, math.nan)),
        ),
    )

    right.axhline(1.0, color="0.35", linewidth=1.0, linestyle="--", zorder=0)
    right.set_xlabel(r"$\phi$ [deg]")
    right.set_ylabel(r"pass-2 CLAS12 / corrected CLAS6")
    right.set_title("Ratio to corrected CLAS6: pass-2 combined")
    right.grid(True, alpha=0.25)
    right.legend(fontsize=8, frameon=True)
    auto_ratio_ylim(right, {PASS2_COMBINED_DISPLAY_LABEL: ratios})

    os.makedirs(output_dir, exist_ok=True)

    output_path = output_path_for_target(
        output_dir=output_dir,
        target=target,
        pass2_choice=pass2_choice,
        simple=True,
    )

    fig.savefig(output_path, dpi=200)
    plt.close(fig)

    print(f"Wrote: {output_path}")


def make_full_plot(
    target: TargetSpec,
    raw_clas6_kin: KinematicPoint,
    corrected_clas6_points: List[DataPoint],
    pass2_choice: BinChoice,
    points_by_label: Dict[str, List[DataPoint]],
    ratios_by_label: Dict[str, List[DataPoint]],
    chi2_by_label: Dict[str, Tuple[float, int, float]],
    output_dir: str,
    y_units: str,
    use_log_cross_section: bool,
) -> None:
    fig, axes = plt.subplots(
        2,
        2,
        figsize=(15.0, 10.0),
        constrained_layout=True,
        sharex=True,
    )

    fig.suptitle(
        figure_suptitle_text(
            target=target,
            raw_clas6_kin=raw_clas6_kin,
            pass2_choice=pass2_choice,
        ),
        fontsize=13,
    )

    top_left = axes[0, 0]
    top_right = axes[0, 1]
    bottom_left = axes[1, 0]
    bottom_right = axes[1, 1]

    plot_dataset(
        ax=top_left,
        points=corrected_clas6_points,
        label=CLAS6_DISPLAY_LABEL,
        legend_label=CLAS6_DISPLAY_LABEL,
    )

    for label in TOP_CROSS_SECTION_SERIES:
        if label not in points_by_label:
            continue
        # endif

        plot_dataset(
            ax=top_left,
            points=points_by_label[label],
            label=label,
            legend_label=format_label_with_chi2(
                label,
                chi2_by_label.get(label, (math.nan, 0, math.nan)),
            ),
        )
    # endfor

    top_left.set_ylabel(rf"$d^4\sigma/(dx_B\,dQ^2\,d|t|\,d\phi)$ [{y_units}]")
    top_left.set_title("Cross sections: pass-2, pass-1, and corrected CLAS6")
    top_left.grid(True, alpha=0.25)
    top_left.legend(fontsize=8, frameon=True)

    if use_log_cross_section:
        top_left.set_yscale("log")
    # endif

    top_ratios = {}

    for label in TOP_RATIO_SERIES:
        if label not in ratios_by_label:
            continue
        # endif

        ratios = ratios_by_label[label]
        top_ratios[label] = ratios

        plot_dataset(
            ax=top_right,
            points=ratios,
            label=label,
            legend_label=format_label_with_chi2(
                label,
                chi2_by_label.get(label, (math.nan, 0, math.nan)),
            ),
        )
    # endfor

    top_right.axhline(1.0, color="0.35", linewidth=1.0, linestyle="--", zorder=0)
    top_right.set_ylabel(r"pass-2 CLAS12 / corrected CLAS6")
    top_right.set_title("Ratio to corrected CLAS6: pass-2 combined")
    top_right.grid(True, alpha=0.25)
    top_right.legend(fontsize=8, frameon=True)
    auto_ratio_ylim(top_right, top_ratios)

    plot_dataset(
        ax=bottom_left,
        points=corrected_clas6_points,
        label=CLAS6_DISPLAY_LABEL,
        legend_label=CLAS6_DISPLAY_LABEL,
    )

    for label in BOTTOM_CROSS_SECTION_SERIES:
        plot_dataset(
            ax=bottom_left,
            points=points_by_label.get(label, []),
            label=label,
            legend_label=format_label_with_chi2(
                label,
                chi2_by_label.get(label, (math.nan, 0, math.nan)),
            ),
        )
    # endfor

    bottom_left.set_xlabel(r"$\phi$ [deg]")
    bottom_left.set_ylabel(rf"$d^4\sigma/(dx_B\,dQ^2\,d|t|\,d\phi)$ [{y_units}]")
    bottom_left.set_title("Cross sections: individual pass-2 periods and corrected CLAS6")
    bottom_left.grid(True, alpha=0.25)
    bottom_left.legend(fontsize=8, frameon=True)

    if use_log_cross_section:
        bottom_left.set_yscale("log")
    # endif

    bottom_ratios = {}

    for label in BOTTOM_RATIO_SERIES:
        ratios = ratios_by_label.get(label, [])
        bottom_ratios[label] = ratios

        plot_dataset(
            ax=bottom_right,
            points=ratios,
            label=label,
            legend_label=format_label_with_chi2(
                label,
                chi2_by_label.get(label, (math.nan, 0, math.nan)),
            ),
        )
    # endfor

    bottom_right.axhline(1.0, color="0.35", linewidth=1.0, linestyle="--", zorder=0)
    bottom_right.set_xlabel(r"$\phi$ [deg]")
    bottom_right.set_ylabel(r"pass-2 CLAS12 / corrected CLAS6")
    bottom_right.set_title("Ratio to corrected CLAS6: individual pass-2 periods")
    bottom_right.grid(True, alpha=0.25)
    bottom_right.legend(fontsize=8, frameon=True)
    auto_ratio_ylim(bottom_right, bottom_ratios)

    os.makedirs(output_dir, exist_ok=True)

    output_path = output_path_for_target(
        output_dir=output_dir,
        target=target,
        pass2_choice=pass2_choice,
        simple=False,
    )

    fig.savefig(output_path, dpi=200)
    plt.close(fig)

    print(f"Wrote: {output_path}")


def make_one_plot(
    target: TargetSpec,
    raw_clas6_kin: KinematicPoint,
    corrected_clas6_points: List[DataPoint],
    pass2_choice: BinChoice,
    points_by_label: Dict[str, List[DataPoint]],
    ratios_by_label: Dict[str, List[DataPoint]],
    chi2_by_label: Dict[str, Tuple[float, int, float]],
    output_dir: str,
    y_units: str,
    use_log_cross_section: bool,
    simple: bool,
) -> None:
    if simple:
        make_simple_plot(
            target=target,
            raw_clas6_kin=raw_clas6_kin,
            corrected_clas6_points=corrected_clas6_points,
            pass2_choice=pass2_choice,
            points_by_label=points_by_label,
            ratios_by_label=ratios_by_label,
            chi2_by_label=chi2_by_label,
            output_dir=output_dir,
            y_units=y_units,
            use_log_cross_section=use_log_cross_section,
        )
    else:
        make_full_plot(
            target=target,
            raw_clas6_kin=raw_clas6_kin,
            corrected_clas6_points=corrected_clas6_points,
            pass2_choice=pass2_choice,
            points_by_label=points_by_label,
            ratios_by_label=ratios_by_label,
            chi2_by_label=chi2_by_label,
            output_dir=output_dir,
            y_units=y_units,
            use_log_cross_section=use_log_cross_section,
        )
    # endif


# -----------------------------------------------------------------------------
# Diagnostics.
# -----------------------------------------------------------------------------

def print_point_summary(label: str, points: List[DataPoint]) -> None:
    total_values = [symmetric_error(p) for p in points if np.isfinite(symmetric_error(p))]
    stat_values = [p.stat for p in points if np.isfinite(p.stat)]
    sys_values = [p.sys for p in points if np.isfinite(p.sys)]
    corr_values = [p.model_corr for p in points if np.isfinite(p.model_corr)]

    mean_total = float(np.mean(total_values)) if len(total_values) > 0 else math.nan
    mean_stat = float(np.mean(stat_values)) if len(stat_values) > 0 else math.nan
    mean_sys = float(np.mean(sys_values)) if len(sys_values) > 0 else math.nan
    mean_corr = float(np.mean(corr_values)) if len(corr_values) > 0 else math.nan
    min_corr = float(np.min(corr_values)) if len(corr_values) > 0 else math.nan
    max_corr = float(np.max(corr_values)) if len(corr_values) > 0 else math.nan

    print(f"  {label}")
    print(f"    points          = {len(points)}")
    print(f"    mean stat err   = {mean_stat:.6g}")
    print(f"    mean sys err    = {mean_sys:.6g}")
    print(f"    mean total err  = {mean_total:.6g}")
    print(f"    KM15 corr mean  = {mean_corr:.6g}")
    print(f"    KM15 corr range = [{min_corr:.6g}, {max_corr:.6g}]")

    if len(points) > 0:
        first = points[0]
        print(
            f"    first point     = phi {first.phi:.3f} deg, "
            f"sigma {first.sigma:.6g}, "
            f"stat {first.stat:.6g}, "
            f"sys {first.sys:.6g}, "
            f"total {symmetric_error(first):.6g}, "
            f"corr {first.model_corr:.6g}"
        )
    # endif


# -----------------------------------------------------------------------------
# Main workflow.
# -----------------------------------------------------------------------------

def parse_targets(selection: str) -> List[TargetSpec]:
    all_targets = [
        TargetSpec(label=label, xB=xB, Q2=Q2, t_abs=t_abs)
        for label, xB, Q2, t_abs in CLAS6_TARGETS
    ]

    if selection.strip().lower() in {"all", ""}:
        return all_targets
    # endif

    wanted = {x.strip() for x in selection.split(",") if x.strip() != ""}
    selected = [t for t in all_targets if t.label in wanted]

    if len(selected) == 0:
        raise ValueError(
            f"No valid targets selected from --targets={selection!r}. "
            f"Available labels: {', '.join(t.label for t in all_targets)}"
        )
    # endif

    return selected


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare pass-2 CLAS12 DVCS cross sections to pass-1 and CLAS6."
    )

    parser.add_argument(
        "pass2_csv",
        help="Input final pass-2 DVCS analysis CSV.",
    )

    parser.add_argument(
        "--simple",
        action="store_true",
        help=(
            "Make a compact 1x2 plot containing only corrected CLAS6, "
            "pass-2 combined 10.6 GeV, and pass-2/CLAS6 ratio. "
            "This mode automatically disables pass-1 and individual run-period panels."
        ),
    )

    parser.add_argument(
        "--pass1-csv",
        default=DEFAULT_PASS1_CSV,
        help=f"Input pass-1 CSV. Default: {DEFAULT_PASS1_CSV}",
    )

    parser.add_argument(
        "--no-pass1",
        action="store_true",
        help="Disable pass-1 Fa18 overlay.",
    )

    parser.add_argument(
        "--clas6-dataset-id",
        type=int,
        default=DEFAULT_CLAS6_DATASET_ID,
        help=f"Gepard dataset id for CLAS6 unpolarized cross sections. Default: {DEFAULT_CLAS6_DATASET_ID}.",
    )

    parser.add_argument(
        "--targets",
        default="all",
        help="Comma-separated target labels, e.g. '26,27,52', or 'all'. Default: all.",
    )

    parser.add_argument(
        "--output-dir",
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory. Default: {DEFAULT_OUTPUT_DIR}",
    )

    parser.add_argument(
        "--pass2-scale",
        type=float,
        default=DEFAULT_PASS2_SCALE,
        help=f"Scale factor for pass-2 CLAS12 values. Default: {DEFAULT_PASS2_SCALE:g}.",
    )

    parser.add_argument(
        "--pass1-scale",
        type=float,
        default=DEFAULT_PASS1_SCALE,
        help=f"Scale factor for pass-1 CLAS12 values. Default: {DEFAULT_PASS1_SCALE:g}.",
    )

    parser.add_argument(
        "--clas6-scale",
        type=float,
        default=DEFAULT_CLAS6_SCALE,
        help=f"Scale factor for CLAS6 Gepard values. Default: {DEFAULT_CLAS6_SCALE:g}.",
    )

    parser.add_argument(
        "--y-units",
        default=r"nb/GeV$^4$",
        help=r"Y-axis unit label. Default: nb/GeV$^4$.",
    )

    parser.add_argument(
        "--avg-period",
        default="10.6 GeV",
        help="Pass-2 period used for weighted average kinematics. Default: '10.6 GeV'.",
    )

    parser.add_argument(
        "--no-clas12-estimated-sys",
        action="store_true",
        help="Disable the extra estimated pass-2 CLAS12 bin-to-bin systematic.",
    )

    parser.add_argument(
        "--clas12-bin-to-bin-sys-frac",
        type=float,
        default=DEFAULT_PASS2_BIN_TO_BIN_SYS_FRAC,
        help=f"Estimated fractional pass-2 systematic. Default: {DEFAULT_PASS2_BIN_TO_BIN_SYS_FRAC:.2f}.",
    )

    parser.add_argument(
        "--pass1-norm-sys-frac",
        type=float,
        default=DEFAULT_PASS1_NORM_SYS_FRAC,
        help=f"Additional pass-1 normalization systematic. Default: {DEFAULT_PASS1_NORM_SYS_FRAC:.2f}.",
    )

    parser.add_argument(
        "--clas6-round-digits",
        type=int,
        default=6,
        help="Rounding digits used to group CLAS6 Gepard points by xB,Q2,|t|. Default: 6.",
    )

    parser.add_argument(
        "--clas6-phi-convention",
        choices=["trento", "gepard"],
        default="trento",
        help=(
            "Use 'trento' to convert Gepard BMK phi back to published-style Trento phi "
            "for Trento-framed data; use 'gepard' to leave phi untouched. Default: trento."
        ),
    )

    parser.add_argument(
        "--match-xB-scale",
        type=float,
        default=0.05,
        help="Scale for xB in nearest-bin matching distance. Default: 0.05.",
    )

    parser.add_argument(
        "--match-Q2-scale",
        type=float,
        default=0.50,
        help="Scale for Q2 in nearest-bin matching distance. Default: 0.50.",
    )

    parser.add_argument(
        "--match-t-scale",
        type=float,
        default=0.10,
        help="Scale for |t| in nearest-bin matching distance. Default: 0.10.",
    )

    parser.add_argument(
        "--linear-cross-section",
        action="store_true",
        help="Use linear y-scale for cross-section panels. Default is log scale.",
    )

    parser.add_argument(
        "--no-km15-correction",
        action="store_true",
        help="Disable KM15 correction and plot raw CLAS6 points.",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.simple:
        args.no_pass1 = True
    # endif

    if args.pass2_scale <= 0.0:
        raise ValueError("--pass2-scale must be positive.")
    # endif

    if args.pass1_scale <= 0.0:
        raise ValueError("--pass1-scale must be positive.")
    # endif

    if args.clas6_scale <= 0.0:
        raise ValueError("--clas6-scale must be positive.")
    # endif

    if args.clas12_bin_to_bin_sys_frac < 0.0:
        raise ValueError("--clas12-bin-to-bin-sys-frac must be non-negative.")
    # endif

    if args.pass1_norm_sys_frac < 0.0:
        raise ValueError("--pass1-norm-sys-frac must be non-negative.")
    # endif

    targets = parse_targets(args.targets)

    g = import_gepard()
    th_km15 = import_km15()

    pass2_required = [
        "Bin Name",
        "xBmin",
        "xBmax",
        "Q2min",
        "Q2max",
        "t_abs_min",
        "t_abs_max",
        "phimin",
        "phimax",
    ]

    required_pass2_periods = SIMPLE_PASS2_CSV_PERIODS if args.simple else ALL_PASS2_CSV_PERIODS

    pass2_required += [
        pass2_cross_section_column(csv_period)
        for csv_period in required_pass2_periods
    ]

    pass2_df = pd.read_csv(args.pass2_csv)
    require_columns(pass2_df, pass2_required, context="pass-2")

    pass1_df: Optional[pd.DataFrame] = None

    if not args.no_pass1:
        if not os.path.exists(args.pass1_csv):
            print()
            print("WARNING: pass-1 CSV was requested but not found:")
            print(f"  {args.pass1_csv}")
            print("  Continuing without pass-1 Fa18.")
        else:
            pass1_required = [
                "Bin Name",
                "phimin",
                "phimax",
                "phiavg",
                "cross sections, ep->epg, exp",
                "cross sections, ep->epg, exp, stat. unc.",
                "cross sections, ep->epg, exp, syst. unc. (up)",
                "cross sections, ep->epg, exp, syst. unc. (down)",
            ]

            pass1_df = pd.read_csv(args.pass1_csv)
            require_columns(pass1_df, pass1_required, context="pass-1")
        # endif
    # endif

    print("Loading CLAS6 data from Gepard...")
    print(f"  dataset id = {args.clas6_dataset_id}")

    if args.simple:
        print("  simple mode = enabled")
        print("  pass-1 overlay = disabled")
        print("  individual pass-2 run-period panels = disabled")
    # endif

    clas6_all_points = load_all_clas6_points(
        g=g,
        dataset_id=args.clas6_dataset_id,
        clas6_scale=args.clas6_scale,
        phi_convention=args.clas6_phi_convention,
    )

    clas6_groups = group_clas6_points(
        clas6_points=clas6_all_points,
        round_digits=args.clas6_round_digits,
    )

    print(f"  loaded nonzero CLAS6 points = {len(clas6_all_points)}")
    print(f"  grouped CLAS6 kinematic bins = {len(clas6_groups)}")

    for target in targets:
        print()
        print("=" * 80)
        print(f"Processing CLAS6 comparison target {target.label}")
        print(
            f"  requested CLAS6 target: "
            f"xB={target.xB:.6f}, Q2={target.Q2:.6f}, |t|={target.t_abs:.6f}"
        )

        raw_clas6_kin, raw_clas6_points, clas6_distance = select_nearest_clas6_group(
            groups=clas6_groups,
            target=target,
            xB_scale=args.match_xB_scale,
            Q2_scale=args.match_Q2_scale,
            t_scale=args.match_t_scale,
        )

        print("  selected CLAS6 group:")
        print(
            f"    xB={raw_clas6_kin.xB:.6f}, "
            f"Q2={raw_clas6_kin.Q2:.6f}, "
            f"|t|={raw_clas6_kin.t_abs:.6f}"
        )
        print(f"    nonzero points={len(raw_clas6_points)}, scaled distance={clas6_distance:.6f}")

        pass2_selected, pass2_choice = select_pass2_bin_nearest_to_kinematics(
            pass2_df=pass2_df,
            target_kin=raw_clas6_kin,
            xB_scale=args.match_xB_scale,
            Q2_scale=args.match_Q2_scale,
            t_scale=args.match_t_scale,
            avg_period=args.avg_period,
        )

        pass2_kin = KinematicPoint(
            xB=pass2_choice.xB,
            Q2=pass2_choice.Q2,
            t_abs=pass2_choice.t_abs,
        )

        print("  selected pass-2 CLAS12 bin:")
        print(f"    Bin Name={pass2_choice.bin_name}, rows={pass2_choice.n_rows}")
        print(
            f"    xB={pass2_choice.xB:.6f}, "
            f"Q2={pass2_choice.Q2:.6f}, "
            f"|t|={pass2_choice.t_abs:.6f}"
        )
        print(f"    scaled distance to CLAS6={pass2_choice.distance:.6f}")

        if args.no_km15_correction:
            corrected_clas6_points = raw_clas6_points
        else:
            corrected_clas6_points = apply_model_correction_to_clas6(
                raw_points=raw_clas6_points,
                g=g,
                th_km15=th_km15,
                clas6_kin=raw_clas6_kin,
                clas6_ebeam=CLAS6_EBEAM,
                clas12_kin=pass2_kin,
                clas12_ebeam=CLAS12_EBEAM,
            )
        # endif

        points_by_label: Dict[str, List[DataPoint]] = {}

        points_by_label[PASS2_COMBINED_DISPLAY_LABEL] = pass2_points_for_period(
            selected=pass2_selected,
            csv_period=PASS2_COMBINED_CSV_PERIOD,
            pass2_scale=args.pass2_scale,
            include_pass2_estimated_sys=not args.no_clas12_estimated_sys,
            pass2_bin_to_bin_sys_frac=args.clas12_bin_to_bin_sys_frac,
        )

        if not args.simple:
            for display_label, csv_period in PASS2_PERIOD_DISPLAY_TO_CSV_PERIOD.items():
                points_by_label[display_label] = pass2_points_for_period(
                    selected=pass2_selected,
                    csv_period=csv_period,
                    pass2_scale=args.pass2_scale,
                    include_pass2_estimated_sys=not args.no_clas12_estimated_sys,
                    pass2_bin_to_bin_sys_frac=args.clas12_bin_to_bin_sys_frac,
                )
            # endfor
        # endif

        if pass1_df is not None:
            pass1_selected = select_bin_by_name(
                df=pass1_df,
                bin_name=pass2_choice.bin_name,
            )

            if len(pass1_selected) > 0:
                points_by_label[PASS1_DISPLAY_LABEL] = pass1_points(
                    selected=pass1_selected,
                    pass1_scale=args.pass1_scale,
                    pass1_norm_sys_frac=args.pass1_norm_sys_frac,
                )

                print("  selected pass-1 CLAS12 bin:")
                print(f"    Bin Name={pass2_choice.bin_name}, rows={len(pass1_selected)}")
            else:
                print("  pass-1 CLAS12 bin not found for selected pass-2 Bin Name.")
            # endif
        # endif

        chi2_by_label = {
            label: chi2_ndf_to_reference(
                comparison_points=points,
                reference_points=corrected_clas6_points,
            )
            for label, points in points_by_label.items()
        }

        ratios_by_label = {
            label: ratio_to_reference(
                numerator_points=points,
                reference_points=corrected_clas6_points,
            )
            for label, points in points_by_label.items()
            if label != PASS1_DISPLAY_LABEL
        }

        print("  dataset summaries:")
        print_point_summary(CLAS6_DISPLAY_LABEL, corrected_clas6_points)

        if args.simple:
            summary_labels = [
                PASS2_COMBINED_DISPLAY_LABEL,
            ]
        else:
            summary_labels = TOP_CROSS_SECTION_SERIES + BOTTOM_CROSS_SECTION_SERIES
        # endif

        for label in summary_labels:
            if label in points_by_label:
                print_point_summary(label, points_by_label[label])
            # endif
        # endfor

        print("  chi2/ndf against corrected CLAS6:")
        for label, chi2_info in chi2_by_label.items():
            chi2, ndf, chi2ndf = chi2_info

            if ndf > 0 and np.isfinite(chi2ndf):
                print(f"    {label:18s}: chi2={chi2:.4f}, ndf={ndf:d}, chi2/ndf={chi2ndf:.4f}")
            else:
                print(f"    {label:18s}: chi2/ndf=N/A")
            # endif
        # endfor

        make_one_plot(
            target=target,
            raw_clas6_kin=raw_clas6_kin,
            corrected_clas6_points=corrected_clas6_points,
            pass2_choice=pass2_choice,
            points_by_label=points_by_label,
            ratios_by_label=ratios_by_label,
            chi2_by_label=chi2_by_label,
            output_dir=args.output_dir,
            y_units=args.y_units,
            use_log_cross_section=not args.linear_cross_section,
            simple=args.simple,
        )
    # endfor


if __name__ == "__main__":
    main()
# endif