# BUILD_MARKER: EXTERNAL_F2_RUNTIME_HARDENED_20260904
#!/usr/bin/env python3
"""
CLAS12 reproduction of the analysis in

  Moradi et al., "Determination of proton electromagnetic form factors
  from DVCS measurements", arXiv:2512.06554.

This version intentionally removes the Hayward et al. (2018) P/IP/CF Sachs
fit families and follows the fitting strategy of arXiv:2512.06554 directly.

Workflow
--------
1. Read the preliminary CLAS12 10.6-GeV unpolarized ep -> e'p'gamma
   four-fold cross sections.
2. Evaluate the KM15 full EP, BH, DVCS, and interference contributions at
   each measured CLAS12 bin center.
3. Construct the five nested BH-dominated samples used in the paper:
      |1 - sigma_BH/sigma_EP| <= 1%, 2%, 3%, 4%, 5%.
4. Add the paper's corresponding 1%, 2%, 3%, 4%, or 5% uncorrelated
   BH-selection uncertainty to each retained cross-section point.
5. Reproduce Fits 1-5 using
      F1(t) = (1 - t/aE)^(-2)
      F2(t) = kappa (1 - t/aM)^(-2).
6. On Set 5 reproduce the P-pole studies:
      Fit 6: aE = aM = a, with PE and PM free.
      Fit 7: PE = PM = P, with aE and aM free.
7. Reproduce Fit 8:
      F2(t) fixed to Kelly,
      F1(t) = (1 - t/aE)^(-PE).
8. Convert F1,F2 -> GE,GM and extract rE,rM from their t=0 slopes.
9. Produce paper-style FF/radius plots plus CLAS12 observable-space plots
   showing the measured cross-section points together with the fitted BH
   predictions. No pseudo "CLAS12 GE/GM data points" are constructed.

10. For the final KM15/VGG99/GK16 purity-model study, retain Moradi's
    observable-level BH criterion.  A point is legitimately BH-like when the
    NET non-BH contribution sigma_DVCS + sigma_INT is small, including through
    cancellation.  Model/data agreement is reported only as a diagnostic and
    never feeds back into purity-model selection, avoiding circular logic.
11. Before any alternative-model closure study or radius fit, require the
    PARTONS BH subprocess to reproduce the same point-by-point BH angular
    structure as Gepard within deliberately loose sanity bounds.  This
    BH-vs-BH preflight is independent of measured data and prevents a phi or
    other kinematic-convention error from masquerading as GPD-model dependence.

Fit implementation
------------------
The paper uses iMinuit, the standard Hessian errors, and Delta chi^2 = 1
for 68% confidence intervals. This script does the same.

By default the CLAS12 fit uncertainty treatment is

  delta_i^2 =
      delta_stat,i^2
    + delta_point-to-point,i^2
    + (c_BH * sigma_i)^2,

where c_BH is 0.01 ... 0.05 for Fits 1 ... 5, following the additional
uncorrelated uncertainty prescription of arXiv:2512.06554.

The CLAS12 correlated scale uncertainties are represented by two
independent Gaussian nuisance parameters:

  beta_global ~ N(0,1)
  beta_comb   ~ N(0,1)

with the per-point multiplicative scale factor

  S_i = (1 + beta_global * 0.0476)
        (1 + beta_comb * combination_sys_i).

Thus the 4.76%% target-thickness/absolute-charge uncertainty moves every
point by the same fraction, while the combination systematic remains fully
correlated but has the correct bin-dependent magnitude.

The fit chi2 is

  chi2 =
      sum_i [(S_i sigma_BH,i - sigma_i)/delta_i]^2
    + beta_global^2 + beta_comb^2.

The point-to-point systematic and both correlated scale nuisances are
included by default.

Efficiency
----------
KM15 is evaluated only once per measured point and cached. The BH cross
section is then represented exactly as

  sigma_BH = A F1^2 + B F1 F2 + C F2^2,

so the Minuit fits use only fast vectorized NumPy operations.

Default paths assume this script is run from

  /u/home/thayward/clas12_analysis_software/analysis_scripts/
      dvcs_cross_section/external_scripts/

Input:
  ../output/csvs/dvcs_pass2_analysis.csv

Output:
  output/emff_from_bh_paper_method/

Run modes
---------
The preliminary publication-facing default uses only the two CLAS
measurements currently retained for the combined extraction:

  * CLAS6 Jo et al. 2015 (Gepard dataset 98)
  * CLAS12 Lee 2026 (authoritative CLAS Physics Database E214M1 release)

CLAS6 Saylor 2018 and Hall A Defurne 2015 are excluded by default and can be
enabled explicitly for diagnostic studies.  They do not enter the preliminary
Jo+Lee model-selected radius prescription.

CLAS12 pass 2 is deliberately excluded from the default publication-facing
workflow and is retained only as an explicit option.

The script also contains:
  * KM15 BH-purity threshold scans from 1--10%;
  * separate threshold scans with published experimental errors only and
    with the Moradi threshold-dependent extra uncertainty;
  * direct elastic comparisons (Kelly, AMT, A1/Bernauer);
  * low-Q2 GE/GD, GM/(mu_p GD), and mu_p GE/GM diagnostics;
  * local BH F1/F2 logarithmic sensitivity diagnostics;
  * a radius-extrapolation bias/variance pseudodata framework using
    polynomial, inverse-polynomial, and continued-fraction Sachs forms.

to bypass the CLAS12 CSV and reproduce the Moradi-et-al. analysis using
Gepard's built-in CLAS6 2015 XUU dataset (ID 98). Results are written to

  output/emff_from_bh_paper_method/clas6_validation/
"""

from __future__ import annotations

import argparse
import ast
import copy
import json
import math
import os
import re
import itertools
import sys
import time
import warnings
import urllib.request
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    from iminuit import Minuit
except ImportError:
    Minuit = None
#endtry
from matplotlib.backends.backend_pdf import PdfPages


# -----------------------------------------------------------------------------
# Constants
# -----------------------------------------------------------------------------
MP = 0.9382720813          # GeV
MP2 = MP * MP
MU_P = 2.79284734463       # proton magnetic moment
KAPPA_P = MU_P - 1.0

# Neutron elastic constants.  The neutron electric "radius" is a mean-square
# charge radius because GE_n(0)=0; it is therefore reported as <r_E,n^2> in
# fm^2 rather than as a positive square-root radius.
MU_N = -1.91304273
KAPPA_N = MU_N
DEFAULT_NDVCS_UNPOLARIZED = "import/ndvcs_clas12_preliminary_unpolarized.txt"
DEFAULT_NDVCS_EBEAM = 10.45


A1_BERNAUER_Q2_MIN = 3.0e-3
PRAD_Q2_MIN = 2.0e-4
GLOBAL_SCALE_FRAC = 0.0476

# Gepard CLAS:2015uuo (H.S. Jo et al., arXiv:1504.02009)
# dataset 98 = 2640-point beam-helicity-independent XUU cross section.
CLAS6_GEPARD_DATASET_ID = 98

# Hall A E00-110 / Defurne et al. 2015, arXiv:1504.05453.
# Original phi-dependent helicity-independent XUU datasets in Gepard:
# 107=Kin2, 108=Kin3, 112=KinX2, 113=KinX3.
HALLA_DEFURNE_GEPARD_DATASET_IDS = (107, 108, 112, 113)
HALLA_DEFURNE_GLOBAL_SCALE_FRAC = 0.028

# Hall A E07-007 / Defurne et al., Nature Communications 8, 1408 (2017),
# arXiv:1703.09442.  Gepard IDs 129--134 are the 450 beam-spin-sum
# (helicity-independent / XUU-equivalent) phi-dependent cross-section points.
# The paper quotes 3.2% point-to-point systematics and correlated components
# 1.0% HRS acceptance, 2.0% luminosity+DAQ dead time, 0.5% trigger efficiency.
HALLA_DEFURNE2017_GEPARD_DATASET_IDS = (129, 130, 131, 132, 133, 134)
HALLA_DEFURNE2017_POINT_TO_POINT_FRAC = 0.032
HALLA_DEFURNE2017_GLOBAL_SCALE_FRAC = math.sqrt(0.010**2 + 0.020**2 + 0.005**2)

JO_GLOBAL_SCALE_FRAC = 0.05

SAYLOR_EBEAM = 5.88
SAYLOR_GLOBAL_SCALE_FRAC = 0.04
DEFAULT_SAYLOR_FILE = "import/saylor_CLAS6.txt"
DEFAULT_PRAD_FILE = "import/prad_normalized_ge.csv"
DEFAULT_GEORGES_FILE = "import/E12-06-114.xlsx"

# Georges et al., PRL 128, 252002 (2022), Hall A E12-06-114.
# Nominal (xB,Q2,Ebeam) settings from Table I.  The supplemental spreadsheet
# contains the actual average xB and Q2 in each t bin, so every row is mapped
# to the nearest nominal setting in the (xB,Q2) plane.
GEORGES_KINEMATIC_SETTINGS = (
    ("Kin-36-1", 0.36, 3.20, 7.38),
    ("Kin-36-2", 0.36, 3.60, 8.52),
    ("Kin-36-3", 0.36, 4.47, 10.59),
    ("Kin-48-1", 0.48, 2.70, 4.49),
    ("Kin-48-2", 0.48, 4.37, 8.85),
    ("Kin-48-3", 0.48, 5.33, 8.85),
    ("Kin-48-4", 0.48, 6.90, 10.99),
    ("Kin-60-1", 0.60, 5.54, 8.52),
    ("Kin-60-3", 0.60, 8.40, 10.59),
)
SAYLOR_SUPPLEMENT_URL = (
    "https://uknowledge.uky.edu/cgi/viewcontent.cgi?"
    "article=1625&context=physastron_facpub&filename=1&type=additional"
)

CLAS6_BEAM_ENERGY = 5.75
DEFAULT_PASS1_CSV = "../imports/clasdb_E214M1.txt"
PASS1_GLOBAL_SCALE_FRAC = 0.31

# Georges 2022 supplemental tables provide statistical and TOTAL systematic
# uncertainties, while the publication states that the latter contain both
# point-to-point and correlated components.  The public table does not provide
# a decomposition that can be reconstructed unambiguously here.  Therefore the
# production fit keeps the quoted total systematic pointwise and does NOT invent
# a Gaussian normalization prior.  Dedicated normalization-tension diagnostics
# below additionally allow the Georges scale to float freely (diagnostic only)
# so its preferred relative normalization can be measured without pretending
# that an authoritative prior width is known.
GEORGES_GLOBAL_SCALE_FRAC = 0.0

# Sachs-fit optimizer controls.  Hayward & Griffioen emphasize that the chosen
# fit must converge robustly over the full replica study.  Start all candidate
# families at a conventional proton-radius scale rather than at an arbitrary
# dipole coefficient, and impose only broad physical guards on the slope.
SACHS_INITIAL_RADIUS_FM = 0.84
SACHS_MIN_RADIUS_FM = 0.20
SACHS_MAX_RADIUS_FM = 1.50

PAPER_SELECTION_COUNTS = {0.01: 98, 0.02: 212, 0.03: 292, 0.04: 404, 0.05: 473}
PAPER_CHI2_NDOF = {1: 0.82, 2: 0.75, 3: 0.73, 4: 0.68, 5: 0.69, 6: 0.67, 7: 0.67, 8: 0.68}
PAPER_DIPOLE_PARAMS = {
    1: (0.966, 0.075, 0.585, 0.096),
    2: (0.982, 0.058, 0.531, 0.073),
    3: (1.033, 0.048, 0.476, 0.054),
    4: (1.005, 0.043, 0.516, 0.053),
    5: (0.958, 0.040, 0.562, 0.055),
}
HBARC = 0.1973269804       # GeV fm

DEFAULT_CSV = "../output/csvs/dvcs_pass2_analysis.csv"
DEFAULT_OUTDIR = "output/emff_from_bh_paper_method"

COL_XB = "xBavg, 10.6 GeV"
COL_Q2 = "Q2avg, 10.6 GeV"
COL_T = "t_abs_avg, 10.6 GeV"
COL_PHI = "phiavg, 10.6 GeV"
COL_XS = "normed cross sections, ep->epg, exp, 10.6 GeV, unpol"
COL_PTP_SYS = "Syst. err (point-to-point total)"
COL_COMB_SYS = "normed cross sections, ep->epg, exp, 10.6 GeV, unpol, combination sys"
COL_SCALE_SYS = "normed cross sections, ep->epg, exp, 10.6 GeV, unpol, total scale sys"


# -----------------------------------------------------------------------------
# Data containers
# -----------------------------------------------------------------------------
@dataclass
class FitResult:
    name: str
    category: str
    bh_cut: float
    npts: int
    npar: int
    chi2: float
    ndof: int
    chi2_ndof: float
    success: bool
    message: str
    params: np.ndarray
    param_names: List[str]
    covariance: np.ndarray
    rE_fm: float
    rE_err_fm: float
    rM_fm: float
    rM_err_fm: float
    norm: float
    norm_err: float
    model_kind: str
    meta: Dict[str, object]


# -----------------------------------------------------------------------------
# Parsing / I/O
# -----------------------------------------------------------------------------
def parse_tuple_cell(value) -> Tuple[float, float, float]:
    """Parse CSV cells such as '(2.01, 0.10, 0)' robustly."""
    if pd.isna(value):
        return np.nan, np.nan, np.nan
    #endif
    if isinstance(value, (tuple, list, np.ndarray)):
        vals = list(value)
    else:
        text = str(value).strip()
        try:
            vals = list(ast.literal_eval(text))
        except Exception:
            text = text.strip("()[]")
            vals = [float(x.strip()) for x in text.split(",") if x.strip()]
        #endtry
    #endif
    vals = vals + [np.nan] * (3 - len(vals))
    return float(vals[0]), float(vals[1]), float(vals[2])
#enddef


def load_clas12_csv(path: Path) -> pd.DataFrame:
    """Load only the columns needed by this analysis and unpack the XS tuple."""
    header = pd.read_csv(path, nrows=0).columns
    required = [COL_XB, COL_Q2, COL_T, COL_PHI, COL_XS]
    optional = [COL_PTP_SYS, COL_COMB_SYS, COL_SCALE_SYS, "bin index", "Bin Name",
                "xBmin", "xBmax", "Q2min", "Q2max",
                "t_abs_min", "t_abs_max", "phimin", "phimax"]
    missing = [c for c in required if c not in header]
    if missing:
        raise KeyError(f"Missing required CSV columns: {missing}")
    #endif

    usecols = required + [c for c in optional if c in header]
    df = pd.read_csv(path, usecols=usecols, low_memory=False)

    parsed = df[COL_XS].apply(parse_tuple_cell)
    df["xs"] = parsed.str[0]
    df["xs_stat"] = parsed.str[1]
    df["xs_aux"] = parsed.str[2]
    df["xB"] = pd.to_numeric(df[COL_XB], errors="coerce")
    df["Q2"] = pd.to_numeric(df[COL_Q2], errors="coerce")
    df["t_abs"] = pd.to_numeric(df[COL_T], errors="coerce")
    df["phi_deg"] = pd.to_numeric(df[COL_PHI], errors="coerce")

    if COL_PTP_SYS in df:
        # Absolute point-to-point systematic uncertainty, in the same units as
        # the measured cross section.
        df["ptp_sys_abs"] = pd.to_numeric(df[COL_PTP_SYS], errors="coerce").fillna(0.0)
    else:
        df["ptp_sys_abs"] = 0.0
    #endif

    if COL_COMB_SYS in df:
        # Fractional run-period-combination uncertainty. Retained for
        # diagnostics; it is already folded into the total-scale column.
        df["comb_sys_frac"] = pd.to_numeric(df[COL_COMB_SYS], errors="coerce").fillna(0.0)
    else:
        df["comb_sys_frac"] = 0.0
    #endif

    if COL_SCALE_SYS in df:
        # Fractional correlated scale uncertainty.
        df["scale_sys_frac"] = pd.to_numeric(df[COL_SCALE_SYS], errors="coerce")
    else:
        df["scale_sys_frac"] = np.nan
    #endif

    good = (
        np.isfinite(df["xB"]) & np.isfinite(df["Q2"]) &
        np.isfinite(df["t_abs"]) & np.isfinite(df["phi_deg"]) &
        np.isfinite(df["xs"]) & np.isfinite(df["xs_stat"]) &
        (df["xB"] > 0) & (df["xB"] < 1) &
        (df["Q2"] > 0) & (df["t_abs"] > 0) &
        (df["xs"] > 0) & (df["xs_stat"] > 0)
    )
    df = df.loc[good].copy().reset_index(drop=True)
    return df
#enddef



def load_clas12_pass1_csv(csv_path: Path) -> pd.DataFrame:
    """
    Load the authoritative CLAS Physics Database release of the CLAS12
    RG-A DVCS cross sections (measurement E214M1).

    The production default is the tab-delimited CLAS database file
    ``clasdb_E214M1.txt``.  For backward compatibility, a legacy preliminary
    ``all_bin_v3.csv`` can still be supplied explicitly with --pass1-csv.

    Authoritative CLAS-database columns:
      bin, x, Q^2, t_average, phi,
      d4sigma/dQ2 dxb dt dphi, Stat. error, Syst. error

    The published Lee et al. result states that the database systematic error
    is the TOTAL systematic uncertainty, containing the 31% overall
    normalization uncertainty in quadrature with the bin-dependent
    systematic uncertainty.  The fit treats that 31% normalization as a
    single correlated nuisance parameter, so for the authoritative database
    table the pointwise component is reconstructed as

        sqrt(max(total_syst^2 - (0.31 * sigma)^2, 0)).

    This avoids double counting the published normalization uncertainty.
    """
    csv_path = Path(csv_path)

    if csv_path.suffix.lower() in {".txt", ".dat"}:
        rows = []
        with csv_path.open("r", encoding="utf-8") as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 8:
                    continue
                #endif
                try:
                    row = [
                        int(parts[0]),
                        float(parts[1]),
                        float(parts[2]),
                        float(parts[3]),
                        float(parts[4]),
                        float(parts[5]),
                        float(parts[6]),
                        float(parts[7]),
                    ]
                except (TypeError, ValueError):
                    continue
                #endtry
                rows.append(row)
            #endfor
        #endwith

        if not rows:
            raise ValueError(
                f"No CLAS E214M1 data rows could be parsed from {csv_path}"
            )
        #endif

        raw = pd.DataFrame(
            rows,
            columns=[
                "bin",
                "xB",
                "Q2",
                "t_abs",
                "phi_deg",
                "xs",
                "xs_stat",
                "xs_sys_total",
            ],
        )

        outdf = raw.copy()
        outdf["_row"] = np.arange(len(outdf), dtype=int)

        total_sys = outdf["xs_sys_total"].to_numpy(float)
        xs = outdf["xs"].to_numpy(float)
        norm_abs = PASS1_GLOBAL_SCALE_FRAC * xs
        residual2 = total_sys * total_sys - norm_abs * norm_abs

        # Consistency audit: if the released systematic column is truly the
        # quadrature sum of a 31% correlated scale uncertainty and nonnegative
        # point-to-point terms, it cannot be smaller than 0.31*xs.  The current
        # E214M1 table contains a small number of exceptions.  Do not silently
        # explain these as rounding: flag them prominently for the Lee-analysis
        # authors while clipping the negative residual variance to zero so the
        # present diagnostic extraction can continue reproducibly.
        bad_norm = np.isfinite(residual2) & (residual2 < 0.0)
        outdf["released_syst_fraction"] = total_sys / xs
        outdf["released_syst_below_nominal_scale"] = bad_norm
        outdf["released_syst_minus_nominal_scale_abs"] = total_sys - norm_abs
        if np.any(bad_norm):
            bad = outdf.loc[
                bad_norm,
                [
                    "bin", "xB", "Q2", "t_abs", "phi_deg", "xs",
                    "xs_stat", "xs_sys_total", "released_syst_fraction",
                    "released_syst_minus_nominal_scale_abs",
                ],
            ].copy()
            detail = "\n".join(
                "  "
                + ", ".join([
                    f"bin={int(r.bin)}",
                    f"xB={r.xB:.6g}",
                    f"Q2={r.Q2:.6g}",
                    f"|t|={r.t_abs:.6g}",
                    f"phi={r.phi_deg:.6g}",
                    f"xs={r.xs:.7g}",
                    f"syst={r.xs_sys_total:.7g}",
                    f"syst/xs={100.0*r.released_syst_fraction:.4f}%",
                ])
                for r in bad.itertuples(index=False)
            )
            warnings.warn(
                "LEE2026 released-systematic consistency warning: "
                f"{int(np.sum(bad_norm))} point(s) have total systematic "
                f"below the nominal {100*PASS1_GLOBAL_SCALE_FRAC:.1f}% "
                "correlated normalization contribution. This cannot follow "
                "from a nonnegative quadrature decomposition at exactly 31%. "
                "The pointwise residual is clipped to zero for the current "
                "analysis; inspect the release/source calculation before a "
                "final publication result.\n" + detail
            )
        #endif
        outdf["ptp_sys_abs"] = np.sqrt(np.maximum(residual2, 0.0))
        outdf["pass1_sys_up"] = outdf["xs_sys_total"]
        outdf["pass1_sys_down"] = outdf["xs_sys_total"]

        outdf["stat"] = outdf["xs_stat"]
        outdf["ptp_sys"] = outdf["ptp_sys_abs"]
        outdf["comb_sys_frac"] = 0.0
        outdf["scale_sys_frac"] = PASS1_GLOBAL_SCALE_FRAC
        outdf["ebeam"] = 10.604

        finite = (
            np.isfinite(outdf["xB"])
            & np.isfinite(outdf["Q2"])
            & np.isfinite(outdf["t_abs"])
            & np.isfinite(outdf["phi_deg"])
            & np.isfinite(outdf["xs"])
            & np.isfinite(outdf["xs_stat"])
            & np.isfinite(outdf["ptp_sys_abs"])
            & (outdf["xs"] > 0.0)
            & (outdf["xs_stat"] > 0.0)
            & (outdf["t_abs"] > 0.0)
        )
        outdf = outdf.loc[finite].copy().reset_index(drop=True)

        print(
            f"[LEE2026] Loaded {len(outdf)} authoritative CLAS Physics "
            f"Database E214M1 points from {csv_path}"
        )
        rel_stat = outdf["xs_stat"] / outdf["xs"]
        rel_total_sys = outdf["xs_sys_total"] / outdf["xs"]
        rel_ptp = outdf["ptp_sys_abs"] / outdf["xs"]
        print(
            f"[LEE2026] median stat/xs={100*np.nanmedian(rel_stat):.1f}%, "
            f"median released total syst/xs="
            f"{100*np.nanmedian(rel_total_sys):.1f}%, "
            f"inferred pointwise syst/xs={100*np.nanmedian(rel_ptp):.1f}%, "
            f"correlated scale={100*PASS1_GLOBAL_SCALE_FRAC:.1f}%"
        )
        return outdf
    #endif

    # Legacy preliminary CSV compatibility.  This branch is retained only so
    # old validation runs can still be reproduced explicitly.
    raw = pd.read_csv(csv_path, low_memory=False)

    required = [
        "xBavg",
        "Q2avg",
        "t_abs_avg",
        "phiavg",
        "cross sections, ep->epg, exp",
        "cross sections, ep->epg, exp, stat. unc.",
        "cross sections, ep->epg, exp, syst. unc. (up)",
        "cross sections, ep->epg, exp, syst. unc. (down)",
        "valid bin",
    ]
    missing = [c for c in required if c not in raw.columns]
    if missing:
        raise KeyError(
            "Legacy CLAS12 CSV is missing required columns: "
            + ", ".join(missing)
        )
    #endif

    valid = pd.to_numeric(raw["valid bin"], errors="coerce") == 1
    work = raw.loc[valid].copy()

    outdf = pd.DataFrame({
        "_row": work.index.to_numpy(int),
        "xB": pd.to_numeric(work["xBavg"], errors="coerce"),
        "Q2": pd.to_numeric(work["Q2avg"], errors="coerce"),
        "t_abs": pd.to_numeric(work["t_abs_avg"], errors="coerce"),
        "phi_deg": pd.to_numeric(work["phiavg"], errors="coerce"),
        "xs": pd.to_numeric(
            work["cross sections, ep->epg, exp"], errors="coerce"
        ),
        "xs_stat": pd.to_numeric(
            work["cross sections, ep->epg, exp, stat. unc."],
            errors="coerce",
        ),
        "pass1_sys_up": pd.to_numeric(
            work["cross sections, ep->epg, exp, syst. unc. (up)"],
            errors="coerce",
        ),
        "pass1_sys_down": pd.to_numeric(
            work["cross sections, ep->epg, exp, syst. unc. (down)"],
            errors="coerce",
        ),
    })

    outdf["ptp_sys_abs"] = 0.5 * (
        np.abs(outdf["pass1_sys_up"].to_numpy(float))
        + np.abs(outdf["pass1_sys_down"].to_numpy(float))
    )
    outdf["stat"] = outdf["xs_stat"]
    outdf["ptp_sys"] = outdf["ptp_sys_abs"]
    outdf["comb_sys_frac"] = 0.0
    outdf["scale_sys_frac"] = PASS1_GLOBAL_SCALE_FRAC
    outdf["ebeam"] = 10.604

    passthrough = [
        "Bin Name",
        "xBmin", "xBmax",
        "Q2min", "Q2max",
        "t_abs_min", "t_abs_max",
        "phimin", "phimax",
        "cross sections, ep->epg, sim(KM15)",
        "cross sections, ep->epg, sim(pureBH)",
    ]
    for col in passthrough:
        if col in work.columns:
            outdf[col] = work[col].to_numpy()
        #endif
    #endfor

    finite = (
        np.isfinite(outdf["xB"])
        & np.isfinite(outdf["Q2"])
        & np.isfinite(outdf["t_abs"])
        & np.isfinite(outdf["phi_deg"])
        & np.isfinite(outdf["xs"])
        & np.isfinite(outdf["xs_stat"])
        & np.isfinite(outdf["ptp_sys_abs"])
        & (outdf["xs"] > 0.0)
        & (outdf["xs_stat"] > 0.0)
        & (outdf["t_abs"] > 0.0)
    )
    outdf = outdf.loc[finite].copy().reset_index(drop=True)

    print(
        f"[LEE2026 legacy] Loaded {len(outdf)} preliminary CLAS12 "
        f"cross-section points from {csv_path}"
    )
    rel_stat = outdf["xs_stat"] / outdf["xs"]
    rel_sys = outdf["ptp_sys_abs"] / outdf["xs"]
    print(
        f"[LEE2026 legacy] median stat/xs="
        f"{100*np.nanmedian(rel_stat):.1f}%, "
        f"median preliminary pointwise syst/xs="
        f"{100*np.nanmedian(rel_sys):.1f}%, "
        f"global scale={100*PASS1_GLOBAL_SCALE_FRAC:.1f}%"
    )
    return outdf
#enddef




def resolve_script_relative_path(path_value: str) -> Path:
    """
    Resolve a user path robustly.

    Absolute paths are used as-is.  Relative paths are interpreted first
    relative to the current working directory for backward compatibility; if
    they do not exist there, they are interpreted relative to this script.
    """
    p = Path(path_value).expanduser()
    if p.is_absolute():
        return p
    #endif

    cwd_candidate = p.resolve()
    if cwd_candidate.exists():
        return cwd_candidate
    #endif

    return (Path(__file__).resolve().parent / p).resolve()
#enddef


def print_run_plan(args) -> None:
    """Print the exact analysis branch and expected output products."""
    if args.only_pass2:
        mode_name = "ONLY PASS-2"
    elif args.only_pass1:
        mode_name = "ONLY CLAS12 LEE 2026"
    elif args.only_clas6:
        mode_name = "ONLY CLAS6 JO 2015"
    elif args.only_saylor:
        mode_name = "ONLY CLAS6 SAYLOR 2018"
    elif args.only_halla_defurne:
        mode_name = "ONLY HALL A DEFURNE 2015"
    elif args.only_clas6_pass1_combined:
        mode_name = "ONLY JO 2015 + CLAS12 LEE 2026 (LEGACY DEBUG MODE)"
    elif args.all_including_pass2:
        mode_name = "FULL THREE-MEASUREMENT STUDY + PASS-2"
    else:
        mode_name = "FULL THREE-MEASUREMENT PUBLISHED STUDY"
    #endif

    print("\n" + "=" * 78)
    print(f"[RUN MODE] {mode_name}")
    print(f"[cwd]      {Path.cwd()}")
    print(f"[script]   {Path(__file__).resolve()}")
    print(f"[outdir]   {Path(args.outdir).expanduser().resolve()}")
    print(f"[Saylor]   {resolve_script_relative_path(args.saylor_file)}")

    if (
        not args.only_pass2
        and not args.only_pass1
        and not args.only_clas6
        and not args.only_saylor
        and not args.only_halla_defurne
        and not args.only_clas6_pass1_combined
    ):
        print("[expected published-study products]")
        print("  singles:")
        print("    CLAS6 Jo 2015")
        print("    CLAS6 Saylor 2018")
        print("    Hall A Defurne 2015")
        print("    CLAS12 Lee 2026")
        print("  pairwise:")
        print("    Jo 2015 + Saylor 2018")
        print("    Jo 2015 + Lee 2026")
        print("    Saylor 2018 + Lee 2026")
        print("  combined:")
        print("    Jo 2015 + Saylor 2018 + Lee 2026")
        print("  summaries:")
        print("    BH-cut plateau")
        print("    separate-vs-combined comparison")
    #endif
    print("=" * 78)
#enddef



def maybe_download_saylor_supplement(path: Path, force: bool = False) -> Path:
    """
    Try to download the public Saylor et al. supplemental TXT.

    Automated downloads from institutional repositories can occasionally
    reject non-browser clients.  If this fails, download Supplemental-material.txt
    manually from the PRC/UKnowledge paper page and save it at ``path``.
    """
    path = Path(path).expanduser().resolve()
    if path.exists() and not force:
        return path
    #endif

    path.parent.mkdir(parents=True, exist_ok=True)
    print(f"[SAYLOR] Downloading supplemental table -> {path}")
    request = urllib.request.Request(
        SAYLOR_SUPPLEMENT_URL,
        headers={"User-Agent": "Mozilla/5.0"},
    )
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            payload = response.read()
        #endwith
        if len(payload) < 1000:
            raise RuntimeError(
                f"downloaded payload is unexpectedly small ({len(payload)} bytes)"
            )
        #endif
        path.write_bytes(payload)
        print(f"[SAYLOR] Downloaded {len(payload):,} bytes.")
    except Exception as exc:
        raise RuntimeError(
            "Automatic Saylor supplemental download failed.  Download the "
            "'Supplemental-material.txt' attachment for Phys. Rev. C 98, "
            "045203 (2018) manually and save it as:\n"
            f"  {path}\n"
            "Then rerun without --download-saylor.\n"
            f"Underlying error: {exc}"
        ) from exc
    #endtry
    return path
#enddef


def _normalize_saylor_header_token(token: str) -> str:
    token = token.strip().lower()
    token = token.replace("<", "").replace(">", "")
    token = token.replace("^", "").replace("{", "").replace("}", "")
    token = token.replace("(", "").replace(")", "")
    token = token.replace("[", "").replace("]", "")
    token = token.replace("-", "")
    token = token.replace("_", "")
    token = token.replace("/", "")
    return token
#enddef


def _saylor_column_role(token: str) -> Optional[str]:
    """Map common supplemental-table header spellings to canonical roles."""
    t = _normalize_saylor_header_token(token)

    if t in {"xb", "x_b", "xbavg", "xmean"} or "xb" == t:
        return "xB"
    #endif
    if t in {"q2", "q2avg", "q2mean"} or t.startswith("q2"):
        return "Q2"
    #endif
    if t in {"t", "tavg", "tmean", "abst", "minust"}:
        return "t_abs"
    #endif
    if t.startswith("phi"):
        return "phi_deg"
    #endif

    # Unpolarized cross section and its errors.  The supplement also contains
    # polarized quantities; require an unpolarized-like marker where possible.
    if any(k in t for k in ["sigmaunp", "sigunp", "xsecunp", "crosssectionunp"]):
        if "stat" in t:
            return "stat"
        #endif
        if "sys" in t or "syst" in t:
            return "sys_total"
        #endif
        return "xs"
    #endif
    if t in {"sig", "sigma", "xsec", "crosssection"}:
        return "xs"
    #endif
    if "stat" in t and not any(k in t for k in ["pol", "asym"]):
        return "stat"
    #endif
    if ("sys" in t or "syst" in t) and not any(k in t for k in ["pol", "asym"]):
        return "sys_total"
    #endif
    return None
#enddef


def load_saylor_supplement(path: str) -> pd.DataFrame:
    """
    Load the Hirlinger Saylor et al. 2018 unpolarized cross-section table.

    The published supplemental TXT contains TWO tables:
      1. unpolarized cross section:
           xB  Q2  -t  phi  sig  stat-err  sys-err
      2. helicity-dependent difference:
           xB  Q2  -t  phi  Dx10^l-sig  stat-err  sys-err

    Only the FIRST table belongs in the unpolarized BH/form-factor analysis.
    The parser therefore stops when the second header is reached.

    The published 4% elastic-normalization uncertainty is removed in quadrature
    from the tabulated total systematic uncertainty to define ``ptp_sys`` and
    is subsequently treated as a correlated normalization nuisance.

    Malformed/non-numeric rows are reported and dropped rather than guessed.
    """
    path = Path(path).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(
            f"Saylor supplemental file not found: {path}\n"
            "Place saylor_CLAS6.txt there or run with --download-saylor."
        )
    #endif

    lines = path.read_text(errors="replace").splitlines()

    # Locate the published unpolarized and helicity-difference table headers.
    unpol_header_idx = None
    helicity_header_idx = None
    for iline, line in enumerate(lines):
        tokens = re.split(r"\s+", line.strip())
        if tokens[:7] == [
            "xB", "Q2", "-t", "phi", "sig", "stat-err", "sys-err"
        ]:
            unpol_header_idx = iline
        #endif
        if len(tokens) >= 5 and tokens[:4] == ["xB", "Q2", "-t", "phi"]:
            if "Dx10^l-sig" in tokens[4]:
                helicity_header_idx = iline
                if unpol_header_idx is not None:
                    break
                #endif
            #endif
        #endif
    #endfor

    if unpol_header_idx is None:
        raise RuntimeError(
            "Could not find the published Saylor unpolarized header "
            "'xB Q2 -t phi sig stat-err sys-err'."
        )
    #endif

    stop_idx = helicity_header_idx if helicity_header_idx is not None else len(lines)
    data_lines = lines[unpol_header_idx + 1:stop_idx]

    rows = []
    malformed = []
    for iline, line in enumerate(
            data_lines, start=unpol_header_idx + 2):
        stripped = line.strip()
        if not stripped:
            continue
        #endif

        tokens = re.split(r"\s+", stripped)
        if len(tokens) < 7:
            malformed.append((iline, stripped, "fewer than 7 columns"))
            continue
        #endif

        try:
            vals = [float(tok) for tok in tokens[:7]]
        except ValueError:
            malformed.append((iline, stripped, "non-numeric token"))
            continue
        #endtry

        rows.append(vals)
    #endfor

    df = pd.DataFrame(
        rows,
        columns=["xB", "Q2", "t_abs", "phi_deg", "xs", "stat", "sys_total"],
    )

    print(
        f"[SAYLOR] Parsed unpolarized table only: {len(df)} numeric rows "
        f"(lines {unpol_header_idx + 2}--{stop_idx})."
    )
    if helicity_header_idx is not None:
        print(
            f"[SAYLOR] Ignoring second supplemental table beginning at line "
            f"{helicity_header_idx + 1}: helicity-dependent Dx10^l-sig."
        )
    #endif

    if malformed:
        print(
            f"[SAYLOR] WARNING: dropped {len(malformed)} malformed "
            f"unpolarized-table row(s):"
        )
        for iline, raw, reason in malformed[:10]:
            print(f"  line {iline}: {reason}: {raw}")
        #endfor
        if len(malformed) > 10:
            print(f"  ... plus {len(malformed) - 10} additional malformed rows")
        #endif
    #endif

    # Canonicalize.
    for col in ["xB", "Q2", "t_abs", "phi_deg", "xs", "stat", "sys_total"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    #endfor

    df["t_abs"] = np.abs(df["t_abs"])
    df["stat"] = np.abs(df["stat"])
    df["sys_total"] = np.abs(df["sys_total"])

    finite = (
        np.isfinite(df["xB"])
        & np.isfinite(df["Q2"])
        & np.isfinite(df["t_abs"])
        & np.isfinite(df["phi_deg"])
        & np.isfinite(df["xs"])
        & np.isfinite(df["stat"])
        & np.isfinite(df["sys_total"])
        & (df["xs"] > 0.0)
        & (df["stat"] > 0.0)
        & (df["sys_total"] >= 0.0)
    )

    n_before = len(df)
    df = df.loc[finite].copy().reset_index(drop=True)
    if len(df) != n_before:
        print(
            f"[SAYLOR] Dropped {n_before - len(df)} additional row(s) failing "
            f"finite/positive cross-section uncertainty requirements."
        )
    #endif

    # Separate the final-paper 5% correlated normalization contribution.
    global_abs = SAYLOR_GLOBAL_SCALE_FRAC * np.abs(df["xs"])
    df["ptp_sys"] = np.sqrt(
        np.maximum(df["sys_total"]**2 - global_abs**2, 0.0)
    )

    df["phi_deg"] = np.mod(df["phi_deg"], 360.0)
    df["source"] = "CLAS6 Saylor 2018"
    df["ebeam"] = SAYLOR_EBEAM

    print(
        f"[SAYLOR] Loaded {len(df)} valid unpolarized points; "
        f"median stat/xs={100*np.nanmedian(df['stat']/df['xs']):.1f}%, "
        f"median total sys/xs="
        f"{100*np.nanmedian(df['sys_total']/df['xs']):.1f}%, "
        f"median point-to-point sys/xs="
        f"{100*np.nanmedian(df['ptp_sys']/df['xs']):.1f}%, "
        f"global scale={100*SAYLOR_GLOBAL_SCALE_FRAC:.1f}%."
    )

    rel_sys_total = df["sys_total"].to_numpy(float) / df["xs"].to_numpy(float)
    n_gt_100 = int(np.sum(rel_sys_total > 1.0))
    n_gt_1000 = int(np.sum(rel_sys_total > 10.0))
    if n_gt_100 > 0:
        print(
            f"[SAYLOR] NOTE: {n_gt_100} unpolarized row(s) have "
            f"sys-err > 100% of the cross section; {n_gt_1000} exceed "
            f"1000%. These numerically valid rows are retained exactly as "
            f"published and therefore receive very small statistical weight."
        )
    #endif

    return df
#enddef


def _nearest_georges_setting(xb: float, q2: float):
    """Return the unique nearest nominal E12-06-114 setting."""
    # Scale the two coordinates by their natural separation so Q2 does not
    # numerically dominate the xB distance.
    ranked = []
    for name, xb0, q20, ebeam in GEORGES_KINEMATIC_SETTINGS:
        d2 = ((float(xb) - xb0) / 0.12)**2 + ((float(q2) - q20) / 1.0)**2
        ranked.append((d2, name, xb0, q20, ebeam))
    #endfor
    ranked.sort(key=lambda item: item[0])
    return ranked[0]
#enddef


def load_georges_supplement(path: str) -> pd.DataFrame:
    """
    Load the Hall A E12-06-114 supplemental spreadsheet (Georges et al. 2022).

    Only the helicity-INDEPENDENT cross section is used here.  The spreadsheet
    also contains the helicity-dependent cross section in columns 9--12; those
    values are deliberately ignored for the present unpolarized BH study.

    IMPORTANT: although the spreadsheet header labels column 3 as t, the
    E12-06-114 analysis bins this quantity as (t - t_min).  Therefore the
    physical momentum transfer used by the BH/DVCS calculation is reconstructed
    row-by-row as t = t_min(xB,Q2) + (t - t_min).  The internal convention is
    then t_abs=|t|.  The two asymmetric systematic columns are retained and a symmetric
    pointwise diagnostic uncertainty is defined as their average magnitude.
    No separate correlated normalization nuisance is invented here because the
    spreadsheet does not decompose its systematic uncertainty that way.
    """
    path = Path(path).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"Georges E12-06-114 spreadsheet not found: {path}")
    #endif

    # The farm Python installation does not necessarily provide pandas' optional
    # Excel engines (openpyxl/xlrd).  For the committed .xlsx supplement, read
    # the first worksheet directly with the Python standard library so Georges
    # has no optional package dependency.  Keep pandas Excel support as a
    # fallback for a user-supplied legacy .xls file.
    if path.suffix.lower() == ".xlsx":
        import re
        import zipfile
        import xml.etree.ElementTree as ET

        ns = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}

        def _xlsx_col_index(cell_ref: str) -> int:
            letters = re.match(r"[A-Za-z]+", str(cell_ref)).group(0).upper()
            idx = 0
            for ch in letters:
                idx = 26 * idx + (ord(ch) - ord("A") + 1)
            #endfor
            return idx - 1
        #enddef

        with zipfile.ZipFile(path, "r") as zf:
            shared = []
            if "xl/sharedStrings.xml" in zf.namelist():
                root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
                for si in root.findall("m:si", ns):
                    shared.append("".join(t.text or "" for t in si.iterfind(".//m:t", ns)))
                #endfor
            #endif

            # The official E12-06-114 workbook stores the published table in
            # its first worksheet.  Resolve that sheet through workbook rels
            # rather than assuming a particular sheet filename.
            wb = ET.fromstring(zf.read("xl/workbook.xml"))
            first_sheet = wb.find("m:sheets/m:sheet", ns)
            if first_sheet is None:
                raise RuntimeError(f"No worksheet found in Georges workbook: {path}")
            #endif
            rel_id = first_sheet.attrib.get("{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id")
            rel_ns = {"r": "http://schemas.openxmlformats.org/package/2006/relationships"}
            rels = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
            target = None
            for rel in rels.findall("r:Relationship", rel_ns):
                if rel.attrib.get("Id") == rel_id:
                    target = rel.attrib.get("Target")
                    break
                #endif
            #endfor
            if target is None:
                raise RuntimeError(f"Could not resolve first worksheet in Georges workbook: {path}")
            #endif
            target = target.lstrip("/")
            if not target.startswith("xl/"):
                target = "xl/" + target
            #endif

            sheet = ET.fromstring(zf.read(target))
            rows = []
            for row in sheet.findall(".//m:sheetData/m:row", ns):
                values = {}
                for cell in row.findall("m:c", ns):
                    ref = cell.attrib.get("r", "A1")
                    col = _xlsx_col_index(ref)
                    typ = cell.attrib.get("t")
                    if typ == "inlineStr":
                        node = cell.find("m:is/m:t", ns)
                        value = "" if node is None else (node.text or "")
                    else:
                        node = cell.find("m:v", ns)
                        value = "" if node is None else (node.text or "")
                        if typ == "s" and value != "":
                            value = shared[int(value)]
                        #endif
                    #endif
                    values[col] = value
                #endfor
                if values:
                    width = max(values) + 1
                    rows.append([values.get(i, "") for i in range(width)])
                #endif
            #endfor
        #endwith

        if len(rows) < 2:
            raise RuntimeError(f"Unexpected empty Georges workbook: {path}")
        #endif
        ncol = max(len(r) for r in rows)
        rows = [r + [""] * (ncol - len(r)) for r in rows]
        raw = pd.DataFrame(rows[1:], columns=rows[0])
    else:
        try:
            raw = pd.read_excel(path, sheet_name=0)
        except ImportError as exc:
            raise ImportError(
                f"Reading Georges input {path.name} requires a pandas Excel engine. "
                "Use the repository default import/E12-06-114.xlsx, which is read "
                "without openpyxl/xlrd, or install the engine required for this file."
            ) from exc
        #endtry
    #endif

    if raw.shape[1] < 8:
        raise RuntimeError(
            f"Unexpected E12-06-114 spreadsheet layout: only {raw.shape[1]} columns"
        )
    #endif

    # The official supplement has fixed semantic column order.  Use positions
    # rather than fragile embedded-newline Excel header strings.
    out = pd.DataFrame({
        "xB": pd.to_numeric(raw.iloc[:, 0], errors="coerce"),
        "Q2": pd.to_numeric(raw.iloc[:, 1], errors="coerce"),
        "t_raw": pd.to_numeric(raw.iloc[:, 2], errors="coerce"),
        "phi_deg": pd.to_numeric(raw.iloc[:, 3], errors="coerce"),
        "xs": pd.to_numeric(raw.iloc[:, 4], errors="coerce"),
        "stat": pd.to_numeric(raw.iloc[:, 5], errors="coerce"),
        "sys_upper": pd.to_numeric(raw.iloc[:, 6], errors="coerce"),
        "sys_lower": pd.to_numeric(raw.iloc[:, 7], errors="coerce"),
    })
    finite = np.all(np.isfinite(out[["xB", "Q2", "t_raw", "phi_deg", "xs", "stat"]]), axis=1)
    out = out.loc[finite].copy().reset_index(drop=True)

    # E12-06-114 publishes its t binning in delta_t = t - t_min (negative).
    # The spreadsheet header simply says t, which is easy to misread.  Rebuild
    # the physical t before passing the points to Gepard/PARTONS.  Using the
    # spreadsheet value directly would put many points at |t| < |t_min| and
    # hence outside the physical DVCS phase space.
    xb_arr = out["xB"].to_numpy(float)
    q2_arr = out["Q2"].to_numpy(float)
    eps2 = 4.0 * MP * MP * xb_arr * xb_arr / q2_arr
    tmin = -q2_arr * (
        2.0 * (1.0 - xb_arr) * (1.0 - np.sqrt(1.0 + eps2)) + eps2
    ) / (4.0 * xb_arr * (1.0 - xb_arr) + eps2)
    out["t_minus_tmin"] = out["t_raw"].to_numpy(float)
    out["t_min"] = tmin
    out["t_physical"] = tmin + out["t_minus_tmin"].to_numpy(float)
    out["t_abs"] = np.abs(out["t_physical"].to_numpy(float))
    out["phi_deg"] = np.mod(out["phi_deg"].to_numpy(float), 360.0)
    out["stat"] = np.abs(out["stat"].to_numpy(float))
    out["sys_upper"] = np.abs(out["sys_upper"].to_numpy(float))
    out["sys_lower"] = np.abs(out["sys_lower"].to_numpy(float))
    out["ptp_sys_abs"] = 0.5 * (out["sys_upper"] + out["sys_lower"])

    setting_names = []
    ebeams = []
    nominal_xb = []
    nominal_q2 = []
    for xb, q2 in zip(out["xB"].to_numpy(float), out["Q2"].to_numpy(float)):
        _, name, xb0, q20, ebeam = _nearest_georges_setting(xb, q2)
        setting_names.append(name)
        nominal_xb.append(xb0)
        nominal_q2.append(q20)
        ebeams.append(ebeam)
    #endfor
    out["setting"] = setting_names
    out["ebeam"] = np.asarray(ebeams, dtype=float)
    out["nominal_xB"] = np.asarray(nominal_xb, dtype=float)
    out["nominal_Q2"] = np.asarray(nominal_q2, dtype=float)
    out["source"] = "Hall A Georges 2022"
    out["_georges_row"] = np.arange(len(out), dtype=int)

    counts = out.groupby("setting", sort=False).size()
    print(f"[GEORGES] Loaded {len(out)} helicity-independent points from {path.name}.")
    print("[GEORGES] Setting counts: " + ", ".join(f"{k}={int(v)}" for k, v in counts.items()))
    print(
        f"[GEORGES] Reconstructed physical |t| from published (t-t_min): "
        f"{out['t_abs'].min():.4f}--{out['t_abs'].max():.4f} GeV^2; "
        f"Q2 range={out['Q2'].min():.2f}--{out['Q2'].max():.2f} GeV^2."
    )
    return out
#enddef


def evaluate_km15_georges_dataframe(
        df: pd.DataFrame,
        workers: int,
        cache_dir: Path,
        force: bool = False) -> pd.DataFrame:
    """KM15/BH decomposition for Georges rows with setting-dependent Ebeam."""
    pieces = []
    for setting_name, group in df.groupby("setting", sort=False):
        ebeams = np.unique(group["ebeam"].to_numpy(float))
        if len(ebeams) != 1:
            raise RuntimeError(f"{setting_name}: non-unique Ebeam values {ebeams}")
        #endif
        tmp = group.sort_values("_georges_row").reset_index(drop=True)
        cache_path = cache_dir / f"km15_bh_decomposition_georges2022_{setting_name}.csv"
        print(
            f"[GEORGES KM15] {setting_name}: {len(tmp)} points at "
            f"Ebeam={ebeams[0]:.2f} GeV"
        )
        evaluated = evaluate_km15_dataframe(
            tmp,
            ebeam=float(ebeams[0]),
            workers=workers,
            cache_path=cache_path,
            force=force,
        )
        pieces.append(evaluated)
    #endfor
    out = pd.concat(pieces, ignore_index=True)
    out = out.sort_values("_georges_row").reset_index(drop=True)
    return out
#enddef


def make_georges_diagnostic_bundle(args) -> Dict[str, object]:
    """Build the Georges bundle for diagnostics and optional KM15 ensemble fits."""
    outdir = Path(args.outdir).expanduser().resolve() / "halla_georges2022"
    outdir.mkdir(parents=True, exist_ok=True)
    df = load_georges_supplement(str(resolve_script_relative_path(args.georges_file)))
    df = evaluate_km15_georges_dataframe(
        df,
        workers=args.workers,
        cache_dir=outdir,
        force=args.force_km15,
    )
    # Canonical aliases used by the shared selection/export machinery.
    df["rbh"] = df["R_BH"].to_numpy(float)
    df["delta_bh"] = df["bh_delta"].to_numpy(float)
    set5 = df.loc[df["delta_bh"] <= 0.05].copy().reset_index(drop=True)
    print(f"[GEORGES] KM15 5% BH-purity diagnostic: {len(set5)}/{len(df)} points")
    return {
        "key": "halla_georges2022",
        "kind": "georges2022",
        "label": "Hall A Georges 2022",
        "all_data": df.reset_index(drop=True),
        "set5": set5,
        "norm_frac": GEORGES_GLOBAL_SCALE_FRAC,
        "georges_systematic_policy": "published total systematic used pointwise; no extra normalization nuisance",
        "diagnostic_only": False,
    }
#enddef


def evaluate_km15_saylor_dataframe(
        df: pd.DataFrame,
        workers: int,
        cache_path: Path,
        force: bool = False) -> pd.DataFrame:
    """KM15/BH decomposition for the 5.88-GeV Saylor bin centers."""
    tmp = df.copy()
    print(
        f"[SAYLOR KM15] Preparing {len(tmp)} unpolarized Saylor points "
        f"at Ebeam={SAYLOR_EBEAM:.2f} GeV."
    )
    return evaluate_km15_dataframe(
        tmp,
        ebeam=SAYLOR_EBEAM,
        workers=workers,
        cache_path=cache_path,
        force=force,
    )
#enddef


def load_clas6_gepard_dataset(dataset_id: int = CLAS6_GEPARD_DATASET_ID) -> pd.DataFrame:
    """
    Load the published CLAS6 unpolarized four-fold cross sections directly
    from Gepard.

    Gepard dataset 98 is:
      CLAS:2015uuo, H.S. Jo et al., arXiv:1504.02009
      2640 points, observable XUU, beam energy 5.75 GeV.

    Gepard's shipped datasets are already transformed internally to the BMK
    convention. We therefore retain the exact stored BMK phi (radians) for
    the validation calculation rather than passing through the CLAS12
    Trento/BMK conversion helper.
    """
    try:
        import gepard as g
    except Exception as exc:
        raise RuntimeError(
            "Could not import gepard. The CLAS6 validation requires the "
            "same Gepard installation used for the KM15 calculation."
        ) from exc
    #endtry

    if dataset_id not in g.dset:
        raise RuntimeError(
            f"Gepard dataset {dataset_id} is not available in this installation."
        )
    #endif

    dset = g.dset[dataset_id]
    rows = []

    for i, pt in enumerate(dset):
        observable = str(getattr(pt, "observable", ""))
        if observable.upper() != "XUU":
            continue
        #endif

        err_total = float(getattr(pt, "err", np.nan))
        err_stat = float(getattr(pt, "errstat", np.nan))
        err_syst = float(getattr(pt, "errsyst", np.nan))

        # The paper fits the published cross section with its published total
        # uncertainty, then adds the 1--5% BH-selection uncertainty.
        if not np.isfinite(err_total) or err_total <= 0.0:
            if np.isfinite(err_stat) and np.isfinite(err_syst):
                err_total = math.hypot(err_stat, err_syst)
            elif np.isfinite(err_stat):
                err_total = err_stat
            #endif
        #endif

        beam_energy = float(
            getattr(pt, "in1energy", CLAS6_BEAM_ENERGY)
        )

        rows.append({
            "_row": i,
            "xB": float(pt.xB),
            "Q2": float(pt.Q2),
            "t_abs": abs(float(pt.t)),
            "phi_bmk_rad": float(pt.phi),
            "phi_deg": math.degrees(float(pt.phi)),
            "xs": float(pt.val),
            # fit_sigma_errors() uses xs_stat as its base error. In validation
            # mode this is intentionally the published TOTAL CLAS6 error.
            "xs_stat": err_total,
            "clas6_err_total": err_total,
            "clas6_err_stat": err_stat,
            "clas6_err_syst": err_syst,
            # Jo et al. quote a 5% uncertainty on the global elastic
            # renormalization factor. Keep the original published total error
            # untouched for the Moradi validation, but provide a second
            # pointwise-only total for combined fits so the 5% contribution is
            # not counted both diagonally and as a correlated nuisance.
            "clas6_err_pointwise": math.sqrt(max(
                err_total**2 - (JO_GLOBAL_SCALE_FRAC * abs(float(pt.val)))**2,
                err_stat**2 if np.isfinite(err_stat) and err_stat > 0.0 else 0.0,
            )),
            "ptp_sys_abs": 0.0,
            "comb_sys_frac": 0.0,
            "scale_sys_frac": 0.0,
            "ebeam": beam_energy,
            "gepard_dataset_id": int(dataset_id),
            "observable": observable,
        })
    #endfor

    df = pd.DataFrame(rows)
    good = (
        np.isfinite(df["xB"])
        & np.isfinite(df["Q2"])
        & np.isfinite(df["t_abs"])
        & np.isfinite(df["phi_bmk_rad"])
        & np.isfinite(df["xs"])
        & np.isfinite(df["xs_stat"])
        & (df["xs_stat"] > 0.0)
        & (df["t_abs"] > 0.0)
    )
    df = df.loc[good].copy().reset_index(drop=True)

    print(
        f"[CLAS6] Loaded Gepard dataset {dataset_id}: "
        f"{len(df)} XUU points"
    )
    return df
#enddef



def load_halla_defurne_gepard_datasets(
        dataset_ids: Sequence[int] = HALLA_DEFURNE_GEPARD_DATASET_IDS
        ) -> pd.DataFrame:
    """Load Hall A Defurne 2015 original phi-dependent XUU data from Gepard."""
    try:
        import gepard as g
    except Exception as exc:
        raise RuntimeError(
            "Could not import gepard. Hall A Defurne requires the same Gepard "
            "installation used by the Jo/KM15 analysis."
        ) from exc
    #endtry

    rows = []
    for dataset_id in dataset_ids:
        if dataset_id not in g.dset:
            raise RuntimeError(f"Gepard dataset {dataset_id} is unavailable.")
        #endif

        for i, pt in enumerate(g.dset[dataset_id]):
            observable = str(getattr(pt, "observable", ""))
            if observable.upper() != "XUU":
                continue
            #endif

            err_total = float(getattr(pt, "err", np.nan))
            err_stat = float(getattr(pt, "errstat", np.nan))
            err_syst = float(getattr(pt, "errsyst", np.nan))
            if not np.isfinite(err_syst) or err_syst < 0.0:
                err_syst = 0.0
            #endif
            if not np.isfinite(err_total) or err_total <= 0.0:
                if np.isfinite(err_stat):
                    err_total = math.hypot(err_stat, err_syst)
                #endif
            #endif

            beam_energy = float(getattr(pt, "in1energy", np.nan))
            if not np.isfinite(beam_energy) or beam_energy <= 0.0:
                raise ValueError(
                    f"Hall A dataset {dataset_id}, point {i}: missing beam energy."
                )
            #endif

            rows.append({
                "_row": len(rows),
                "xB": float(pt.xB),
                "Q2": float(pt.Q2),
                "t_abs": abs(float(pt.t)),
                "phi_bmk_rad": float(pt.phi),
                "phi_deg": math.degrees(float(pt.phi)),
                "xs": float(pt.val),
                "xs_stat": err_stat,
                "halla_err_stat": err_stat,
                "halla_err_syst": err_syst,
                "halla_err_total": err_total,
                "ptp_sys_abs": err_syst,
                "comb_sys_frac": 0.0,
                "scale_sys_frac": HALLA_DEFURNE_GLOBAL_SCALE_FRAC,
                "ebeam": beam_energy,
                "gepard_dataset_id": int(dataset_id),
                "observable": observable,
                "source": "Hall A Defurne 2015",
            })
        #endfor
    #endfor

    df = pd.DataFrame(rows)
    good = (
        np.isfinite(df["xB"])
        & np.isfinite(df["Q2"])
        & np.isfinite(df["t_abs"])
        & np.isfinite(df["phi_bmk_rad"])
        & np.isfinite(df["xs"])
        & np.isfinite(df["halla_err_stat"])
        & np.isfinite(df["halla_err_syst"])
        & (df["xs"] > 0.0)
        & (df["halla_err_stat"] > 0.0)
        & (df["halla_err_syst"] >= 0.0)
        & (df["t_abs"] > 0.0)
    )
    df = df.loc[good].copy().reset_index(drop=True)
    df["_row"] = np.arange(len(df), dtype=int)

    counts = df.groupby("gepard_dataset_id").size().to_dict()
    print(
        f"[HALLA] Loaded {len(df)} XUU points from Gepard datasets "
        f"{tuple(dataset_ids)}; per-dataset counts={counts}; "
        f"global scale={100*HALLA_DEFURNE_GLOBAL_SCALE_FRAC:.1f}%."
    )
    print(
        f"[HALLA] median stat/xs="
        f"{100*np.nanmedian(df['halla_err_stat']/df['xs']):.1f}%, "
        f"median pointwise syst/xs="
        f"{100*np.nanmedian(df['halla_err_syst']/df['xs']):.1f}%."
    )
    return df
#enddef



def load_halla_defurne2017_gepard_datasets(
        dataset_ids: Sequence[int] = HALLA_DEFURNE2017_GEPARD_DATASET_IDS
        ) -> pd.DataFrame:
    """
    Load Hall A E07-007 / Defurne et al. 2017 unpolarized cross sections.

    Gepard labels the original beam-spin-sum datasets 129--134 as ``BSS``.
    The Gepard dataset documentation identifies these six sets as the 450 XUU
    measurements at Ebeam = 3.355, 4.455, and 5.55 GeV.

    Error treatment follows the publication rather than double-counting its
    total systematic uncertainty:
      * statistical uncertainty: Gepard errstat;
      * point-to-point systematic: 3.2% of the measured cross section;
      * correlated normalization nuisance:
            sqrt(1.0%^2 + 2.0%^2 + 0.5%^2) = 2.291...%.
    The latter combines HRS acceptance, luminosity/dead-time, and trigger.
    """
    try:
        import gepard as g
    except Exception as exc:
        raise RuntimeError(
            "Could not import gepard. Hall A Defurne 2017 requires the same "
            "Gepard installation used by the Jo/KM15 analysis."
        ) from exc
    #endtry

    rows = []
    for dataset_id in dataset_ids:
        if dataset_id not in g.dset:
            raise RuntimeError(
                f"Gepard Defurne-2017 dataset {dataset_id} is unavailable. "
                "Update/check the Gepard installation; expected IDs are "
                f"{tuple(dataset_ids)}."
            )
        #endif

        for i, pt in enumerate(g.dset[dataset_id]):
            observable = str(getattr(pt, "observable", "")).upper()
            # The shipped E07-007 originals are BSS.  Accept XUU as a
            # defensive compatibility alias in case a Gepard release renames
            # the observable while retaining the same dataset IDs.
            if observable not in ("BSS", "XUU"):
                continue
            #endif

            xs = float(getattr(pt, "val", np.nan))
            err_stat = float(getattr(pt, "errstat", np.nan))
            err_total_gepard = float(getattr(pt, "err", np.nan))
            err_syst_gepard = float(getattr(pt, "errsyst", np.nan))

            beam_energy = float(getattr(pt, "in1energy", np.nan))
            if not np.isfinite(beam_energy) or beam_energy <= 0.0:
                raise ValueError(
                    f"Hall A Defurne 2017 dataset {dataset_id}, point {i}: "
                    "missing beam energy."
                )
            #endif

            point_sys = (
                HALLA_DEFURNE2017_POINT_TO_POINT_FRAC * abs(xs)
                if np.isfinite(xs) else np.nan
            )

            rows.append({
                "_row": len(rows),
                "xB": float(pt.xB),
                "Q2": float(pt.Q2),
                "t_abs": abs(float(pt.t)),
                "phi_bmk_rad": float(pt.phi),
                "phi_deg": math.degrees(float(pt.phi)),
                "xs": xs,
                "xs_stat": err_stat,
                "halla_err_stat": err_stat,
                "halla_err_syst": point_sys,
                "halla_err_total": math.hypot(err_stat, point_sys),
                "gepard_err_total_original": err_total_gepard,
                "gepard_err_syst_original": err_syst_gepard,
                "ptp_sys_abs": point_sys,
                "comb_sys_frac": 0.0,
                "scale_sys_frac": HALLA_DEFURNE2017_GLOBAL_SCALE_FRAC,
                "ebeam": beam_energy,
                "gepard_dataset_id": int(dataset_id),
                "observable": observable,
                "source": "Hall A Defurne 2017 (E07-007)",
            })
        #endfor
    #endfor

    df = pd.DataFrame(rows)
    if len(df) == 0:
        raise RuntimeError(
            "No Defurne 2017 BSS/XUU points were loaded from Gepard IDs "
            f"{tuple(dataset_ids)}."
        )
    #endif

    good = (
        np.isfinite(df["xB"])
        & np.isfinite(df["Q2"])
        & np.isfinite(df["t_abs"])
        & np.isfinite(df["phi_bmk_rad"])
        & np.isfinite(df["xs"])
        & np.isfinite(df["halla_err_stat"])
        & (df["xs"] > 0.0)
        & (df["halla_err_stat"] > 0.0)
        & (df["t_abs"] > 0.0)
    )
    df = df.loc[good].copy().reset_index(drop=True)
    df["_row"] = np.arange(len(df), dtype=int)

    counts = df.groupby("gepard_dataset_id").size().to_dict()
    print(
        f"[HALLA 2017] Loaded {len(df)} unpolarized E07-007 points from "
        f"Gepard datasets {tuple(dataset_ids)}; per-dataset counts={counts}."
    )
    print(
        f"[HALLA 2017] uncertainty prescription: "
        f"{100*HALLA_DEFURNE2017_POINT_TO_POINT_FRAC:.1f}% point-to-point + "
        f"{100*HALLA_DEFURNE2017_GLOBAL_SCALE_FRAC:.3f}% correlated scale "
        "(1.0% HRS acceptance, 2.0% luminosity/dead time, 0.5% trigger)."
    )
    print(
        f"[HALLA 2017] median stat/xs="
        f"{100*np.nanmedian(df['halla_err_stat']/df['xs']):.1f}%."
    )
    return df
#enddef


def evaluate_km15_halla_dataframe(
        df: pd.DataFrame,
        workers: int,
        cache_path: Path,
        force: bool = False) -> pd.DataFrame:
    """Use exact Gepard-stored BMK phi and per-point beam energy."""
    print(
        f"[HALLA KM15] Evaluating {len(df)} Defurne XUU points with "
        f"{workers} worker(s)..."
    )
    return evaluate_km15_clas6_dataframe(
        df,
        workers=workers,
        cache_path=cache_path,
        force=force,
    )
#enddef



# -----------------------------------------------------------------------------
# Gepard / KM15 calculations
# -----------------------------------------------------------------------------
def make_gepard_point(g, xB: float, Q2: float, t_abs: float,
                      phi_deg: float, ebeam: float):
    phi_rad = math.radians(phi_deg)
    phi_trento = math.pi - phi_rad
    pt = g.DataPoint(
        xB=float(xB),
        t=-abs(float(t_abs)),
        Q2=float(Q2),
        phi=float(phi_trento),
        observable="XS",
        frame="trento",
        process="ep2epgamma",
        exptype="fixed target",
        in1energy=float(ebeam),
        in1charge=-1,
        in1polarization=0,
        in2particle="p",
    )
    pt.prepare()
    return pt
#enddef


def _set_eff_basis(model, f1_value: float, f2_value: float):
    """
    Temporarily replace F1/F2 on the model instance.

    Instance-assigned callables are not descriptors, so the lambda receives
    only the DataPoint argument, which matches self.m.F1(pt) usage in Gepard.
    """
    old_f1 = model.F1
    old_f2 = model.F2
    model.F1 = lambda pt: float(f1_value)
    model.F2 = lambda pt: float(f2_value)
    return old_f1, old_f2
#enddef


def _restore_eff(model, old_f1, old_f2) -> None:
    model.F1 = old_f1
    model.F2 = old_f2
#enddef


def evaluate_km15_point(args: Tuple[int, float, float, float, float, float, str]) -> Dict[str, float]:
    """
    Worker-safe KM15 evaluation for one point.

    Returns full EP, BH, DVCS^2, INT, R_BH, and the quadratic BH coefficients
    A/B/C such that sigma_BH=A*F1^2+B*F1*F2+C*F2^2.
    """
    idx, xB, Q2, t_abs, phi_deg, ebeam, phi_convention = args
    try:
        import gepard as g
        from gepard.fits import th_KM15
    except Exception as exc:
        raise RuntimeError(
            "Could not import gepard/KM15. Activate the Python environment used "
            "by the existing km15_cli.py before running this script."
        ) from exc
    #endtry

    # Each worker process evaluates tasks sequentially. Reuse the imported KM15
    # object rather than deep-copying the full CFF/GPD model for every bin.
    # F1/F2 are restored immediately after the three BH basis evaluations.
    th = th_KM15
    pt = make_gepard_point(g, xB, Q2, t_abs, phi_deg, ebeam)

    pref = float(th.PreFacSigma(pt))
    bh_amp = float(th.TBH2unp(pt))
    int_amp = float(th.TINTunp(pt))
    dvcs_amp = float(th.TDVCS2unp(pt))

    sigma_bh = pref * bh_amp
    sigma_int = pref * int_amp
    sigma_dvcs = pref * dvcs_amp
    sigma_ep = sigma_bh + sigma_int + sigma_dvcs

    sigma_predict = float(th.predict(pt))
    ep_relerr = abs(sigma_ep - sigma_predict) / max(abs(sigma_predict), 1e-30)

    # Build exact quadratic coefficients for the BH-only fit.
    # Gepard's BH term is quadratic in F1,F2, so three basis evaluations suffice.
    model = th.m
    try:
        old_f1, old_f2 = _set_eff_basis(model, 1.0, 0.0)
        A = pref * float(th.TBH2unp(pt))
        model.F1 = lambda p: 0.0
        model.F2 = lambda p: 1.0
        C = pref * float(th.TBH2unp(pt))
        model.F1 = lambda p: 1.0
        model.F2 = lambda p: 1.0
        one_one = pref * float(th.TBH2unp(pt))
        B = one_one - A - C
    finally:
        _restore_eff(model, old_f1, old_f2)
    #endtry

    # Validate the quadratic representation using the model's nominal FFs.
    f1_nom = float(model.F1(pt))
    f2_nom = float(model.F2(pt))
    sigma_bh_reco = A * f1_nom**2 + B * f1_nom * f2_nom + C * f2_nom**2
    denom = max(abs(sigma_bh), 1e-30)
    quad_relerr = abs(sigma_bh_reco - sigma_bh) / denom

    rbh = sigma_bh / sigma_ep if sigma_ep != 0 else np.nan
    return {
        "_row": idx,
        "km15_ep": sigma_ep,
        "km15_bh": sigma_bh,
        "km15_dvcs": sigma_dvcs,
        "km15_int": sigma_int,
        "R_BH": rbh,
        "bh_delta": abs(1.0 - rbh) if np.isfinite(rbh) else np.nan,
        "bh_A": A,
        "bh_B": B,
        "bh_C": C,
        "bh_quad_relerr": quad_relerr,
        "km15_F1": f1_nom,
        "km15_F2": f2_nom,
        "km15_ep_predict": sigma_predict,
        "km15_ep_decomp_relerr": ep_relerr,
        "phi_model_deg": transform_neutron_phi_deg(phi_deg, phi_convention),
        "phi_convention": phi_convention,
    }
#enddef



def evaluate_km15_point_bmk(args: Tuple[int, float, float, float, float, float]) -> Dict[str, float]:
    """
    Evaluate one Gepard CLAS6 point using its exact stored BMK azimuth.

    This avoids any round-trip through the CLAS12 phi convention and is the
    cleanest possible validation against the paper.
    """
    idx, xB, Q2, t_abs, phi_bmk_rad, ebeam = args

    try:
        import gepard as g
        from gepard.fits import th_KM15
    except Exception as exc:
        raise RuntimeError("Could not import gepard/KM15.") from exc
    #endtry

    pt = g.DataPoint(
        xB=float(xB),
        t=-abs(float(t_abs)),
        Q2=float(Q2),
        phi=float(phi_bmk_rad),
        observable="XUU",
        frame="BMK",
        process="ep2epgamma",
        exptype="fixed target",
        in1energy=float(ebeam),
        in1charge=-1,
        in1polarization=0,
        in2particle="p",
    )
    pt.prepare()

    th = th_KM15
    pref = float(th.PreFacSigma(pt))
    sigma_bh = pref * float(th.TBH2unp(pt))
    sigma_int = pref * float(th.TINTunp(pt))
    sigma_dvcs = pref * float(th.TDVCS2unp(pt))
    sigma_ep = sigma_bh + sigma_int + sigma_dvcs

    sigma_predict = float(th.predict(pt))
    ep_relerr = abs(sigma_ep - sigma_predict) / max(abs(sigma_predict), 1.0e-30)

    model = th.m
    try:
        old_f1, old_f2 = _set_eff_basis(model, 1.0, 0.0)
        A = pref * float(th.TBH2unp(pt))

        model.F1 = lambda p: 0.0
        model.F2 = lambda p: 1.0
        C = pref * float(th.TBH2unp(pt))

        model.F1 = lambda p: 1.0
        model.F2 = lambda p: 1.0
        one_one = pref * float(th.TBH2unp(pt))
        B = one_one - A - C
    finally:
        _restore_eff(model, old_f1, old_f2)
    #endtry

    f1_nom = float(model.F1(pt))
    f2_nom = float(model.F2(pt))
    sigma_bh_reco = A * f1_nom**2 + B * f1_nom * f2_nom + C * f2_nom**2
    quad_relerr = abs(sigma_bh_reco - sigma_bh) / max(abs(sigma_bh), 1.0e-30)

    rbh = sigma_bh / sigma_ep if sigma_ep != 0.0 else np.nan

    return {
        "_row": idx,
        "km15_ep": sigma_ep,
        "km15_bh": sigma_bh,
        "km15_dvcs": sigma_dvcs,
        "km15_int": sigma_int,
        "R_BH": rbh,
        "bh_delta": abs(1.0 - rbh) if np.isfinite(rbh) else np.nan,
        "bh_A": A,
        "bh_B": B,
        "bh_C": C,
        "bh_quad_relerr": quad_relerr,
        "km15_F1": f1_nom,
        "km15_F2": f2_nom,
        "km15_ep_predict": sigma_predict,
        "km15_ep_decomp_relerr": ep_relerr,
        "phi_model_deg": float(np.degrees(phi_bmk_rad) % 360.0),
        "phi_convention": "Gepard BMK direct",
    }
#enddef


def evaluate_km15_clas6_dataframe(df: pd.DataFrame,
                                  workers: int,
                                  cache_path: Path,
                                  force: bool = False) -> pd.DataFrame:
    """Evaluate/cache KM15 at the exact Gepard CLAS6 BMK data points."""
    core_cache_cols = [
        "_row", "km15_ep", "km15_bh", "km15_dvcs", "km15_int",
        "R_BH", "bh_delta", "bh_A", "bh_B", "bh_C",
        "bh_quad_relerr", "km15_F1", "km15_F2",
        "km15_ep_predict", "km15_ep_decomp_relerr",
    ]
    metadata_cache_cols = ["phi_model_deg", "phi_convention"]
    cache_cols = core_cache_cols + metadata_cache_cols

    if cache_path.exists() and not force:
        cache = pd.read_csv(cache_path)
        core_cache_valid = (
            len(cache) == len(df)
            and all(c in cache.columns for c in core_cache_cols)
        )
        if core_cache_valid:
            cache = cache.copy()
            augmented_cache = False

            if "phi_model_deg" not in cache.columns:
                cache["phi_model_deg"] = (
                    np.degrees(df["phi_bmk_rad"].to_numpy(float)) % 360.0
                )
                augmented_cache = True
            #endif

            if "phi_convention" not in cache.columns:
                cache["phi_convention"] = "Gepard BMK direct"
                augmented_cache = True
            #endif

            if augmented_cache:
                cache.to_csv(cache_path, index=False)
                print(
                    "[CLAS6 KM15] Loaded legacy physics cache and "
                    f"backfilled phi metadata without reevaluation: {cache_path}"
                )
            else:
                print(f"[CLAS6 KM15] Loaded cache: {cache_path}")
            #endif

            outdf = df.copy()
            for c in cache_cols:
                if c != "_row":
                    outdf[c] = cache[c].to_numpy()
                #endif
            #endfor
            return outdf
        #endif
    #endif

    tasks = [
        (
            int(i),
            float(r.xB),
            float(r.Q2),
            float(r.t_abs),
            float(r.phi_bmk_rad),
            float(r.ebeam),
        )
        for i, r in df[
            ["xB", "Q2", "t_abs", "phi_bmk_rad", "ebeam"]
        ].iterrows()
    ]

    print(
        f"[CLAS6 KM15] Evaluating {len(tasks)} Gepard points "
        f"with {workers} worker(s)..."
    )
    t0 = time.time()

    if workers <= 1:
        rows = []
        for j, task in enumerate(tasks, start=1):
            rows.append(evaluate_km15_point_bmk(task))
            if j % 200 == 0 or j == len(tasks):
                print(
                    f"[CLAS6 KM15] {j:5d}/{len(tasks)} "
                    f"elapsed={time.time()-t0:8.1f} s"
                )
            #endif
        #endfor
    else:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            rows = list(ex.map(evaluate_km15_point_bmk, tasks, chunksize=12))
        #endwith
    #endif

    km = pd.DataFrame(rows).sort_values("_row").reset_index(drop=True)

    missing = [c for c in cache_cols if c not in km.columns]
    if missing:
        raise RuntimeError(
            "Internal CLAS6 KM15 evaluator/cache schema mismatch; "
            f"missing columns: {missing}"
        )
    #endif

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    km.to_csv(cache_path, index=False)

    print(
        f"[CLAS6 KM15] Finished in {time.time()-t0:.1f} s; "
        f"cache -> {cache_path}"
    )

    outdf = df.copy()
    for c in cache_cols:
        if c != "_row":
            outdf[c] = km[c].to_numpy()
        #endif
    #endfor
    return outdf
#enddef


def evaluate_km15_dataframe(df: pd.DataFrame, ebeam: float, workers: int,
                            cache_path: Path, force: bool = False) -> pd.DataFrame:
    """Evaluate or load cached KM15/BH decomposition."""
    cache_cols = ["_row", "km15_ep", "km15_bh", "km15_dvcs", "km15_int",
                  "R_BH", "bh_delta", "bh_A", "bh_B", "bh_C",
                  "bh_quad_relerr", "km15_F1", "km15_F2",
                  "km15_ep_predict", "km15_ep_decomp_relerr"]
    if cache_path.exists() and not force:
        cache = pd.read_csv(cache_path)
        if len(cache) == len(df) and all(c in cache.columns for c in cache_cols):
            out = df.copy()
            for c in cache_cols:
                if c != "_row":
                    out[c] = cache[c].to_numpy()
                #endif
            #endfor
            print(f"[KM15] Loaded cache: {cache_path}")
            return out
        #endif
    #endif

    # Proton CLAS12 inputs use the native/identity phi convention.
    # evaluate_km15_point() also serves phi-convention-aware bookkeeping and
    # therefore expects the convention as the seventh task element.
    tasks = [
        (
            int(i),
            float(r.xB),
            float(r.Q2),
            float(r.t_abs),
            float(r.phi_deg),
            float(ebeam),
            "identity",
        )
        for i, r in df[["xB", "Q2", "t_abs", "phi_deg"]].iterrows()
    ]

    print(f"[KM15] Evaluating {len(tasks)} CLAS12 points with {workers} worker(s)...")
    t0 = time.time()
    if workers <= 1:
        rows = []
        for j, task in enumerate(tasks, start=1):
            rows.append(evaluate_km15_point(task))
            if j % 100 == 0 or j == len(tasks):
                print(f"[KM15] {j:5d}/{len(tasks)}  elapsed={time.time()-t0:8.1f} s")
            #endif
        #endfor
    else:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            rows = list(ex.map(evaluate_km15_point, tasks, chunksize=8))
        #endwith
    #endif

    km = pd.DataFrame(rows).sort_values("_row").reset_index(drop=True)
    if np.nanmax(km["bh_quad_relerr"].to_numpy()) > 1e-8:
        raise RuntimeError(
            "The quadratic BH coefficient reconstruction failed validation. "
            f"max relative error={np.nanmax(km['bh_quad_relerr']):.3e}."
        )
    #endif
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    km.to_csv(cache_path, index=False)
    print(f"[KM15] Finished in {time.time()-t0:.1f} s; cache -> {cache_path}")

    out = df.copy()
    for c in cache_cols:
        if c != "_row":
            out[c] = km[c].to_numpy()
        #endif
    #endfor
    return out
#enddef


# -----------------------------------------------------------------------------
# Elastic FF utilities
# -----------------------------------------------------------------------------
def kelly_f1_f2(q: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Kelly (2004) proton elastic form factors, evaluated at q=|t|.

    G_E and G_M are written in tau=q/(4 M_p^2), then converted to F1/F2.
    """
    q = np.asarray(q, dtype=float)
    tau = q / (4.0 * MP2)

    ge = (
        1.0 - 0.24 * tau
    ) / (
        1.0 + 10.98 * tau + 12.82 * tau**2 + 21.97 * tau**3
    )

    gm = MU_P * (
        1.0 + 0.12 * tau
    ) / (
        1.0 + 10.97 * tau + 18.86 * tau**2 + 6.55 * tau**3
    )

    return sachs_to_f1f2(q, ge, gm)
#enddef



def standard_dipole_sachs(q: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Standard dipole Sachs form factors."""
    q = np.asarray(q, dtype=float)
    gd = (1.0 + q / 0.71) ** -2
    return gd, MU_P * gd
#enddef


def kelly_sachs(q: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Kelly (2004) proton Sachs form factors."""
    f1, f2 = kelly_f1_f2(q)
    return f1f2_to_sachs(np.asarray(q, dtype=float), f1, f2)
#enddef


def amt2007_sachs(q: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Arrington-Melnitchouk-Tjon (2007) TPE-corrected global fit.

    Rational form in tau = Q2/(4 Mp^2):
      G / G(0) = (1 + a1*tau + a2*tau^2 + a3*tau^3) /
                 (1 + b1*tau + ... + b5*tau^5)

    Coefficients are from Table I of Phys. Rev. C 76, 035205 (2007).
    """
    q = np.asarray(q, dtype=float)
    tau = q / (4.0 * MP2)

    def rational(a, b):
        num = 1.0
        for i, ai in enumerate(a, start=1):
            num = num + ai * tau**i
        #endfor
        den = 1.0
        for i, bi in enumerate(b, start=1):
            den = den + bi * tau**i
        #endfor
        return num / den
    #enddef

    ge = rational(
        [3.439, -1.602, 0.068],
        [15.055, 48.061, 99.304, 0.012, 8.650],
    )
    gm_over_mu = rational(
        [-1.465, 1.260, 0.262],
        [9.627, 0.000, 0.000, 11.179, 13.245],
    )
    return ge, MU_P * gm_over_mu
#enddef


def bernauer_polyxdipole_sachs(q: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    A1/Bernauer polynomial-times-standard-dipole parameterization.

    Coefficients are from Appendix J.1 of J.C. Bernauer's Mainz thesis
    for the n=8 polynomial x dipole fit.  GM is returned with its physical
    normalization mu_p.
    """
    q = np.asarray(q, dtype=float)
    gd = (1.0 + q / 0.71) ** -2

    ae = np.array([
        -0.4980, +5.4592, -34.7281, +124.3173,
        -262.9808, +329.1395, -227.3306, +66.6980,
    ])
    am = np.array([
        +0.2472, -4.9123, +29.7509, -84.0430,
        +129.3256, -111.1068, +49.9753, -9.1659,
    ])

    pe = np.ones_like(q)
    pm = np.ones_like(q)
    for i, coeff in enumerate(ae, start=1):
        pe += coeff * q**i
    #endfor
    for i, coeff in enumerate(am, start=1):
        pm += coeff * q**i
    #endfor

    return gd * pe, MU_P * gd * pm
#enddef



# YAHL18: Ye, Arrington, Hill, Lee, Phys. Lett. B 777 (2018) 8--15.
# Central proton Sachs form factors in the published bounded z expansion:
# G(Q^2)=sum_{k=0}^{12} a_k z^k, t_cut=4 m_pi^2, t0=-0.7 GeV^2.
YAHL18_MPI_GEV = 0.13957039
YAHL18_TCUT_GEV2 = 4.0 * YAHL18_MPI_GEV**2
YAHL18_T0_GEV2 = -0.7

YAHL18_GE_COEFFS = np.array([
    +0.23916, -1.10986, +1.44438, +0.47957, -2.28689,
    +1.12663, +1.25062, -3.63102, +4.08222, +0.50410,
    -5.08512, +3.96774, -0.98153,
], dtype=float)

YAHL18_GM_OVER_MU_COEFFS = np.array([
    +0.26414, -1.09531, +1.21855, +0.66114, -1.40568,
    -1.35642, +1.44703, +4.23567, -5.33405, -2.91630,
    +8.70740, -5.70700, +1.28081,
], dtype=float)


def yahl18_z(q: np.ndarray) -> np.ndarray:
    """YAHL18 conformal variable for spacelike Q^2=q>=0."""
    q = np.asarray(q, dtype=float)
    a = np.sqrt(YAHL18_TCUT_GEV2 + q)
    b = math.sqrt(YAHL18_TCUT_GEV2 - YAHL18_T0_GEV2)
    return (a - b) / (a + b)
#enddef


def _yahl18_eval(coeffs: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Evaluate a YAHL18 central-value z polynomial via Horner's method."""
    z = yahl18_z(q)
    out = np.zeros_like(z, dtype=float)
    for coeff in coeffs[::-1]:
        out = out * z + float(coeff)
    #endfor
    return out
#enddef


def yahl18_sachs(q: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """YAHL18 proton GE and physical GM (GM(0)=mu_p)."""
    q = np.asarray(q, dtype=float)
    ge = _yahl18_eval(YAHL18_GE_COEFFS, q)
    gm_over_mu = _yahl18_eval(YAHL18_GM_OVER_MU_COEFFS, q)
    return ge, MU_P * gm_over_mu
#enddef


def yahl18_f1_f2(q: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """YAHL18 proton Dirac and Pauli form factors."""
    ge, gm = yahl18_sachs(q)
    return sachs_to_f1f2(np.asarray(q, dtype=float), ge, gm)
#enddef


def validate_yahl18_implementation() -> Dict[str, float]:
    """Validate normalization and the radii imposed in the YAHL18 central fit."""
    ge0, gm0 = yahl18_sachs(np.array([0.0], dtype=float))
    rE = radius_from_shape(
        lambda qq: yahl18_sachs(np.asarray(qq, dtype=float))[0], 1.0
    )
    rM = radius_from_shape(
        lambda qq: yahl18_sachs(np.asarray(qq, dtype=float))[1], MU_P
    )
    if abs(float(ge0[0]) - 1.0) > 2.0e-5:
        raise RuntimeError(f"YAHL18 GE(0) check failed: {ge0[0]}")
    #endif
    if abs(float(gm0[0]) / MU_P - 1.0) > 2.0e-5:
        raise RuntimeError(
            f"YAHL18 GM(0)/mu_p check failed: {gm0[0]/MU_P}"
        )
    #endif
    if abs(float(rE) - 0.879) > 3.0e-4:
        raise RuntimeError(f"YAHL18 rE check failed: {rE:.6f} fm")
    #endif
    if abs(float(rM) - 0.851) > 3.0e-4:
        raise RuntimeError(f"YAHL18 rM check failed: {rM:.6f} fm")
    #endif
    return {
        "GE0": float(ge0[0]),
        "GM0_over_mu": float(gm0[0] / MU_P),
        "rE_fm": float(rE),
        "rM_fm": float(rM),
    }
#enddef


def bernauer_rosenbluth_data() -> pd.DataFrame:
    """
    Free-GM Rosenbluth-separated A1/Bernauer proton form factors.

    Values are from Table K.3 of the Bernauer Mainz thesis.  The four
    constrained-GM points listed separately in the source are intentionally
    omitted.  The quoted errors are pointwise 68.3% errors; GE-GM covariance
    is not available in this simple tabulation, so pulls made with these
    values are diagnostic rather than a covariance-complete chi2.
    """
    raw = """
0.0152 0.9433 0.0071 3.2859 0.3955
0.0162 0.9430 0.0042 3.1148 0.2332
0.0172 0.9346 0.0045 3.2880 0.2134
0.0182 0.9434 0.0037 2.6790 0.2057
0.0192 0.9328 0.0053 3.1029 0.2263
0.0202 0.9354 0.0041 2.7368 0.1861
0.0213 0.9375 0.0046 2.4972 0.2096
0.0236 0.9349 0.0039 2.3928 0.1712
0.0255 0.9250 0.0029 2.4870 0.1138
0.0272 0.9170 0.0027 2.6058 0.0894
0.0304 0.9119 0.0030 2.4760 0.0973
0.0322 0.9053 0.0026 2.5483 0.0768
0.0340 0.9021 0.0017 2.4843 0.0480
0.0357 0.8921 0.0021 2.6470 0.0460
0.0376 0.8921 0.0024 2.5045 0.0477
0.0394 0.8888 0.0016 2.4698 0.0386
0.0413 0.8800 0.0014 2.5439 0.0259
0.0433 0.8767 0.0010 2.4878 0.0201
0.0473 0.8657 0.0010 2.4788 0.0164
0.0494 0.8622 0.0012 2.4500 0.0174
0.0515 0.8579 0.0013 2.4175 0.0135
0.0537 0.8534 0.0009 2.4014 0.0111
0.0557 0.8474 0.0013 2.3997 0.0123
0.0580 0.8419 0.0009 2.3938 0.0121
0.0602 0.8364 0.0012 2.3790 0.0088
0.0625 0.8326 0.0008 2.3574 0.0076
0.0648 0.8274 0.0007 2.3404 0.0051
0.0672 0.8213 0.0010 2.3275 0.0082
0.0695 0.8162 0.0008 2.3157 0.0077
0.0711 0.8115 0.0007 2.3072 0.0055
0.0742 0.8033 0.0007 2.3028 0.0058
0.0773 0.7968 0.0007 2.2861 0.0040
0.0805 0.7909 0.0005 2.2611 0.0031
0.0837 0.7858 0.0005 2.2327 0.0039
0.0867 0.7794 0.0038 2.2087 0.0416
0.0908 0.7707 0.0024 2.1924 0.0257
0.0949 0.7625 0.0031 2.1873 0.0319
0.1033 0.7427 0.0029 2.1669 0.0260
0.1076 0.7379 0.0026 2.1125 0.0213
0.1120 0.7272 0.0011 2.1069 0.0087
0.1163 0.7212 0.0017 2.0652 0.0100
0.1255 0.7055 0.0016 2.0243 0.0110
0.1300 0.6990 0.0013 1.9862 0.0071
0.1348 0.6923 0.0012 1.9522 0.0062
0.1395 0.6857 0.0015 1.9296 0.0065
0.1441 0.6779 0.0012 1.9096 0.0059
0.1490 0.6676 0.0009 1.8941 0.0034
0.1538 0.6642 0.0016 1.8522 0.0081
0.1588 0.6551 0.0012 1.8419 0.0070
0.1637 0.6483 0.0013 1.8164 0.0053
0.1687 0.6362 0.0013 1.8162 0.0064
0.1739 0.6305 0.0009 1.7866 0.0028
0.1789 0.6247 0.0011 1.7642 0.0043
0.1840 0.6154 0.0010 1.7447 0.0036
0.1893 0.6111 0.0018 1.7222 0.0031
0.1943 0.6029 0.0014 1.7027 0.0035
0.1997 0.5931 0.0008 1.6891 0.0023
0.2051 0.5831 0.0009 1.6746 0.0030
0.2103 0.5785 0.0010 1.6523 0.0028
0.2156 0.5730 0.0008 1.6293 0.0026
0.2210 0.5687 0.0021 1.6015 0.0068
0.2318 0.5575 0.0021 1.5609 0.0067
0.2427 0.5370 0.0019 1.5429 0.0056
0.2481 0.5355 0.0022 1.5166 0.0059
0.2592 0.5257 0.0023 1.4813 0.0054
0.2758 0.5067 0.0022 1.4351 0.0040
0.2814 0.4998 0.0026 1.4222 0.0063
0.3069 0.4723 0.0012 1.3619 0.0024
0.3354 0.4447 0.0012 1.2929 0.0020
0.3638 0.4192 0.0010 1.2310 0.0015
0.3916 0.4031 0.0019 1.1633 0.0029
0.4201 0.3824 0.0020 1.1105 0.0027
0.4471 0.3609 0.0027 1.0697 0.0034
0.4748 0.3448 0.0026 1.0214 0.0029
0.5011 0.3301 0.0018 0.9789 0.0017
0.5274 0.3105 0.0021 0.9428 0.0018
0.5524 0.2972 0.0021 0.9108 0.0016
"""
    rows = []
    for line in raw.strip().splitlines():
        q2, ge, dge, gm, dgm = map(float, line.split())
        rows.append({
            "Q2": q2,
            "GE": ge,
            "GE_err": dge,
            "GM": gm,
            "GM_err": dgm,
        })
    #endfor
    return pd.DataFrame(rows)
#enddef


def load_prad_normalized_ge(
        path: str | Path = DEFAULT_PRAD_FILE) -> pd.DataFrame:
    """
    Load the published PRad normalized proton electric form-factor points.

    The authoritative numerical values are from the Jefferson Lab PRad public
    data release, readme.pdf pp. 4--5:
      https://www.jlab.org/prad/data/readme.pdf

    These are G_E^p(Q2)=f(Q2)/n, i.e. with the fitted beam-energy-dependent
    floating normalizations removed.  The public release quotes
      n = 1.0002 for 1.101 GeV,
      n = 0.9983 for 2.143 GeV.

    Statistical and systematic uncertainties are retained separately.  The
    plotting comparison uses their quadrature sum for a compact visual
    uncertainty; the source columns remain available in the import CSV.
    """
    p = Path(path).expanduser()
    if not p.exists():
        raise FileNotFoundError(
            f"PRad normalized-GE table not found: {p}. "
            "Place prad_normalized_ge.csv in external_scripts/import/."
        )
    #endif

    df = pd.read_csv(p)
    required = {
        "beam_energy_MeV", "theta_deg", "Q2_GeV2",
        "GE", "GE_stat", "GE_syst",
    }
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(
            "PRad table is missing required columns: "
            + ", ".join(sorted(missing))
        )
    #endif

    for col in required:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    #endfor
    df = df.loc[
        np.isfinite(df["Q2_GeV2"])
        & np.isfinite(df["GE"])
        & np.isfinite(df["GE_stat"])
        & np.isfinite(df["GE_syst"])
    ].copy()
    df["GE_total"] = np.hypot(
        df["GE_stat"].to_numpy(float),
        df["GE_syst"].to_numpy(float),
    )
    df["beam_energy_GeV"] = df["beam_energy_MeV"] / 1000.0
    return df.sort_values(
        ["beam_energy_MeV", "Q2_GeV2"]
    ).reset_index(drop=True)
#enddef


def elastic_reference_curves(q: np.ndarray) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    """Named direct-elastic reference curves used consistently in all plots."""
    return {
        "Kelly 2004": kelly_sachs(q),
        "AMT 2007": amt2007_sachs(q),
        "YAHL18": yahl18_sachs(q),
        "A1/Bernauer order-8 poly×dipole": bernauer_polyxdipole_sachs(q),
    }
#enddef


def f1f2_to_sachs(q: np.ndarray, f1: np.ndarray, f2: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    tau = np.asarray(q, dtype=float) / (4.0 * MP2)
    ge = f1 - tau * f2
    gm = f1 + f2
    return ge, gm
#enddef


def sachs_to_f1f2(q: np.ndarray, ge: np.ndarray, gm: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    tau = np.asarray(q, dtype=float) / (4.0 * MP2)
    f1 = (ge + tau * gm) / (1.0 + tau)
    f2 = (gm - ge) / (1.0 + tau)
    return f1, f2
#enddef


def bh_from_f1f2(A: np.ndarray, B: np.ndarray, C: np.ndarray,
                  f1: np.ndarray, f2: np.ndarray) -> np.ndarray:
    return A * f1**2 + B * f1 * f2 + C * f2**2
#enddef



# -----------------------------------------------------------------------------
# Paper-method fit helpers
# -----------------------------------------------------------------------------
def numerical_gradient(fun: Callable[[np.ndarray], float],
                       p: np.ndarray) -> np.ndarray:
    """Central finite-difference gradient used only for propagated quantities."""
    p = np.asarray(p, dtype=float)
    grad = np.zeros_like(p)
    for i in range(len(p)):
        step = 1.0e-5 * max(1.0, abs(p[i]))
        pp = p.copy()
        pm = p.copy()
        pp[i] += step
        pm[i] -= step
        grad[i] = (fun(pp) - fun(pm)) / (2.0 * step)
    #endfor
    return grad
#enddef


def propagate_scalar(fun: Callable[[np.ndarray], float],
                     p: np.ndarray,
                     cov: np.ndarray) -> Tuple[float, float]:
    value = float(fun(p))
    if cov.size == 0 or not np.all(np.isfinite(cov)):
        return value, np.nan
    #endif
    grad = numerical_gradient(fun, p)
    var = float(grad @ cov @ grad)
    return value, math.sqrt(max(var, 0.0)) if np.isfinite(var) else np.nan
#enddef


def radius_from_shape(shape_fun: Callable[[np.ndarray], np.ndarray],
                      g0: float) -> float:
    """
    RMS radius in fm from a form factor written as a function of q=|t|=-t.

      r^2 = -6/G(0) dG/dq |_(q=0)

    which is equivalent to the paper's r^2 = +6/G(0) dG/dt at t=0.
    """
    h = 1.0e-6
    g_at_0 = float(shape_fun(np.array([0.0]))[0])
    g_at_h = float(shape_fun(np.array([h]))[0])
    slope_q = (g_at_h - g_at_0) / h
    r2_gev2 = -6.0 * slope_q / g0
    if not np.isfinite(r2_gev2) or r2_gev2 < 0.0:
        return np.nan
    #endif
    return math.sqrt(r2_gev2) * HBARC
#enddef


def fit_sigma_errors(data: pd.DataFrame,
                     bh_cut: float,
                     include_clas12_ptp_sys: bool,
                     include_bh_sys: bool) -> np.ndarray:
    """
    Pointwise uncertainty for the paper-method fit.

    Nominal CLAS12 treatment:
      sqrt(stat^2 + ptp_total^2 + (BH_cut * sigma)^2)
    """
    stat = data["xs_stat"].to_numpy(float)
    var = stat**2

    if include_clas12_ptp_sys:
        ptp = data["ptp_sys_abs"].to_numpy(float)
        var = var + ptp**2
    #endif

    if include_bh_sys:
        var = var + (bh_cut * data["xs"].to_numpy(float))**2
    #endif

    err = np.sqrt(var)
    bad = ~np.isfinite(err) | (err <= 0.0)
    if np.any(bad):
        err[bad] = stat[bad]
    #endif
    return err
#enddef


# -----------------------------------------------------------------------------
# Moradi et al. F1/F2 models: Fits 1-8
# -----------------------------------------------------------------------------
def paper_model_f1f2(kind: str,
                     q: np.ndarray,
                     pars: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Evaluate the exact fit forms used in arXiv:2512.06554.

    q = |t| = -t, hence (1 - t/a) = (1 + q/a).
    """
    q = np.asarray(q, dtype=float)
    kind = kind.lower()

    if kind == "dipole":
        aE, aM = pars
        f1 = (1.0 + q / aE)**(-2.0)
        f2 = KAPPA_P * (1.0 + q / aM)**(-2.0)
        return f1, f2
    #endif

    if kind == "fit6_same_a":
        a, pE, pM = pars
        f1 = (1.0 + q / a)**(-pE)
        f2 = KAPPA_P * (1.0 + q / a)**(-pM)
        return f1, f2
    #endif

    if kind == "fit7_same_p":
        aE, aM, p = pars
        f1 = (1.0 + q / aE)**(-p)
        f2 = KAPPA_P * (1.0 + q / aM)**(-p)
        return f1, f2
    #endif

    if kind == "fit8_f2_kelly":
        aE, pE = pars
        f1 = (1.0 + q / aE)**(-pE)
        _, f2 = kelly_f1_f2(q)
        return f1, f2
    #endif

    raise ValueError(f"Unknown paper fit kind: {kind}")
#enddef


def paper_model_setup(kind: str) -> Tuple[List[str], np.ndarray]:
    """Initial values are seeded near the paper's reported solutions."""
    kind = kind.lower()

    if kind == "dipole":
        return ["aE", "aM"], np.array([0.95, 0.56], dtype=float)
    #endif

    if kind == "fit6_same_a":
        return ["a", "PE", "PM"], np.array([0.53, 1.2, 1.96], dtype=float)
    #endif

    if kind == "fit7_same_p":
        return ["aE", "aM", "P"], np.array([0.65, 0.35, 1.42], dtype=float)
    #endif

    if kind == "fit8_f2_kelly":
        return ["aE", "PE"], np.array([0.565, 1.318], dtype=float)
    #endif

    raise ValueError(f"Unknown paper fit kind: {kind}")
#enddef


def fit_paper_model(data: pd.DataFrame,
                    kind: str,
                    fit_name: str,
                    bh_cut: float,
                    include_clas12_ptp_sys: bool,
                    include_bh_sys: bool,
                    include_combination_norm_sys: bool,
                    include_global_norm_sys: bool,
                    global_scale_frac: float = GLOBAL_SCALE_FRAC) -> FitResult:
    """
    Paper-method chi2 fit with two correlated CLAS12 scale nuisances.

      beta_global : common 4.76% scale uncertainty
      beta_comb   : common nuisance with per-bin coefficient comb_sys_frac_i

      pred_i = sigma_BH,i
               * (1 + beta_global * global_scale_frac)
               * (1 + beta_comb * comb_sys_frac_i)

      chi2 += beta_global^2 + beta_comb^2
    """
    q = data["t_abs"].to_numpy(float)
    y = data["xs"].to_numpy(float)
    err = fit_sigma_errors(
        data,
        bh_cut,
        include_clas12_ptp_sys,
        include_bh_sys,
    )
    A = data["bh_A"].to_numpy(float)
    B = data["bh_B"].to_numpy(float)
    C = data["bh_C"].to_numpy(float)

    comb = data["comb_sys_frac"].to_numpy(float)
    comb = np.where(np.isfinite(comb) & (comb >= 0.0), comb, 0.0)

    names, p0 = paper_model_setup(kind)
    fit_names = list(names)
    fit_p0 = list(p0)

    if include_global_norm_sys:
        fit_names.append("beta_global")
        fit_p0.append(0.0)
    #endif

    use_comb = include_combination_norm_sys and np.any(comb > 0.0)
    if use_comb:
        fit_names.append("beta_comb")
        fit_p0.append(0.0)
    #endif

    fit_p0 = np.asarray(fit_p0, dtype=float)

    def split_params(p: np.ndarray):
        p = np.asarray(p, dtype=float)
        ff_pars = p[:len(names)]
        cursor = len(names)

        beta_global = 0.0
        beta_comb = 0.0

        if include_global_norm_sys:
            beta_global = float(p[cursor])
            cursor += 1
        #endif

        if use_comb:
            beta_comb = float(p[cursor])
        #endif

        return ff_pars, beta_global, beta_comb
    #enddef

    def scale_vector(beta_global: float,
                     beta_comb: float) -> np.ndarray:
        global_factor = (
            1.0 + beta_global * global_scale_frac
            if include_global_norm_sys
            else 1.0
        )
        comb_factor = (
            1.0 + beta_comb * comb
            if use_comb
            else np.ones_like(comb)
        )
        return global_factor * comb_factor
    #enddef

    def chi2_array(p: np.ndarray) -> float:
        ff_pars, beta_global, beta_comb = split_params(p)
        f1, f2 = paper_model_f1f2(kind, q, ff_pars)
        bare_bh = bh_from_f1f2(A, B, C, f1, f2)
        pred = scale_vector(beta_global, beta_comb) * bare_bh
        pull = (pred - y) / err
        chi2 = float(np.dot(pull, pull))

        if include_global_norm_sys:
            chi2 += beta_global**2
        #endif
        if use_comb:
            chi2 += beta_comb**2
        #endif
        return chi2
    #enddef

    def minuit_chi2(*values):
        return chi2_array(np.asarray(values, dtype=float))
    #enddef

    m = Minuit(minuit_chi2, *fit_p0, name=tuple(fit_names))
    m.errordef = Minuit.LEAST_SQUARES

    for name in names:
        m.limits[name] = (1.0e-6, None)
    #endfor
    if "beta_global" in fit_names:
        m.limits["beta_global"] = (-10.0, 10.0)
    #endif
    if "beta_comb" in fit_names:
        m.limits["beta_comb"] = (-10.0, 10.0)
    #endif

    m.migrad()
    m.hesse()

    all_params = np.array(
        [float(m.values[name]) for name in fit_names],
        dtype=float,
    )

    full_cov = np.full((len(fit_names), len(fit_names)), np.nan, dtype=float)
    if m.covariance is not None:
        for i, ni in enumerate(fit_names):
            for j, nj in enumerate(fit_names):
                full_cov[i, j] = float(m.covariance[ni, nj])
            #endfor
        #endfor
    #endif

    params = all_params[:len(names)]
    cov = full_cov[:len(names), :len(names)]

    beta_global = (
        float(all_params[fit_names.index("beta_global")])
        if "beta_global" in fit_names
        else 0.0
    )
    beta_comb = (
        float(all_params[fit_names.index("beta_comb")])
        if "beta_comb" in fit_names
        else 0.0
    )

    def par_error(name: str) -> float:
        if name not in fit_names or not np.all(np.isfinite(full_cov)):
            return np.nan
        #endif
        idx = fit_names.index(name)
        return math.sqrt(max(full_cov[idx, idx], 0.0))
    #enddef

    beta_global_err = par_error("beta_global")
    beta_comb_err = par_error("beta_comb")

    global_norm = 1.0 + beta_global * global_scale_frac
    global_norm_err = (
        global_scale_frac * beta_global_err
        if np.isfinite(beta_global_err)
        else np.nan
    )

    chi2 = float(m.fval)
    npar = len(all_params)
    ndof = max(len(y) - npar, 1)

    def re_of_p(p):
        return radius_from_shape(
            lambda qq: f1f2_to_sachs(
                qq, *paper_model_f1f2(kind, qq, p)
            )[0],
            1.0,
        )
    #enddef

    def rm_of_p(p):
        return radius_from_shape(
            lambda qq: f1f2_to_sachs(
                qq, *paper_model_f1f2(kind, qq, p)
            )[1],
            MU_P,
        )
    #enddef

    rE, drE = propagate_scalar(re_of_p, params, cov)
    rM, drM = propagate_scalar(rm_of_p, params, cov)

    message = (
        f"valid={m.valid}; accurate={m.accurate}; "
        f"edm={float(m.fmin.edm):.3e}; "
        f"has_posdef_covar={bool(m.fmin.has_posdef_covar)}"
    )

    return FitResult(
        name=fit_name,
        category="paper_F1F2",
        bh_cut=bh_cut,
        npts=len(y),
        npar=npar,
        chi2=chi2,
        ndof=ndof,
        chi2_ndof=chi2 / ndof,
        success=bool(m.valid),
        message=message,
        params=params,
        param_names=names,
        covariance=cov,
        rE_fm=rE,
        rE_err_fm=drE,
        rM_fm=rM,
        rM_err_fm=drM,
        norm=global_norm,
        norm_err=global_norm_err,
        model_kind=kind,
        meta={
            "minuit_valid": bool(m.valid),
            "minuit_accurate": bool(m.accurate),
            "minuit_edm": float(m.fmin.edm),
            "has_posdef_covar": bool(m.fmin.has_posdef_covar),
            "include_clas12_ptp_sys": bool(include_clas12_ptp_sys),
            "include_bh_selection_sys": bool(include_bh_sys),
            "include_combination_norm_sys": bool(use_comb),
            "include_global_norm_sys": bool(include_global_norm_sys),
            "global_scale_frac": float(global_scale_frac),
            "beta_global": float(beta_global),
            "beta_global_err": float(beta_global_err),
            "beta_comb": float(beta_comb),
            "beta_comb_err": float(beta_comb_err),
            "combination_sys_median": float(np.nanmedian(comb)),
            "combination_sys_min": float(np.nanmin(comb)),
            "combination_sys_max": float(np.nanmax(comb)),
        },
    )
#enddef




def dataset_statistical_errors(
        data: pd.DataFrame,
        dataset_kind: str) -> np.ndarray:
    """
    Return the statistical uncertainty for one published measurement.

    The three loaders intentionally preserve their native column names, so all
    cross-dataset code must come through this helper rather than assuming a
    universal ``stat`` or ``err`` column.

      Jo 2015:       clas6_err_stat when available; published total error as
                     a defensive fallback because some Gepard points do not
                     carry a separate finite statistical component.
      Saylor 2018:   stat
      CLAS12 Lee:    xs_stat
    """
    if dataset_kind == "jo2015":
        if "clas6_err_stat" in data.columns:
            arr = data["clas6_err_stat"].to_numpy(float)
            fallback = (
                data["clas6_err_total"].to_numpy(float)
                if "clas6_err_total" in data.columns
                else data["xs_stat"].to_numpy(float)
            )
            bad = ~np.isfinite(arr) | (arr <= 0.0)
            arr = arr.copy()
            arr[bad] = fallback[bad]
            return arr
        elif "clas6_err_total" in data.columns:
            return data["clas6_err_total"].to_numpy(float)
        elif "xs_stat" in data.columns:
            return data["xs_stat"].to_numpy(float)
        elif "err" in data.columns:
            return data["err"].to_numpy(float)
        else:
            raise KeyError(
                "Jo 2015 dataframe has no usable statistical/error column; "
                "expected clas6_err_stat, clas6_err_total, or xs_stat."
            )
        #endif
    elif dataset_kind == "saylor2018":
        if "stat" not in data.columns:
            raise KeyError(
                "Saylor 2018 dataframe is missing canonical column 'stat'."
            )
        #endif
        return data["stat"].to_numpy(float)
    elif dataset_kind == "georges2022":
        if "stat" not in data.columns:
            raise KeyError("Hall A Georges dataframe is missing 'stat'.")
        #endif
        return data["stat"].to_numpy(float)
    elif dataset_kind in ("halla_defurne2015", "halla_defurne2017"):
        if "halla_err_stat" in data.columns:
            return data["halla_err_stat"].to_numpy(float)
        elif "xs_stat" in data.columns:
            return data["xs_stat"].to_numpy(float)
        else:
            raise KeyError(
                "Hall A Defurne dataframe is missing 'halla_err_stat'."
            )
        #endif
    elif dataset_kind in ("pass1", "pass2"):
        if "xs_stat" in data.columns:
            return data["xs_stat"].to_numpy(float)
        elif "stat" in data.columns:
            return data["stat"].to_numpy(float)
        else:
            label = "CLAS12 pass-1/Lee 2026" if dataset_kind == "pass1" else "CLAS12 pass-2"
            raise KeyError(
                f"{label} dataframe is missing statistical uncertainty "
                "column 'xs_stat'."
            )
        #endif
    else:
        raise ValueError(f"unknown dataset kind: {dataset_kind}")
    #endif
#enddef


def dataset_point_errors(
        data: pd.DataFrame,
        dataset_kind: str,
        bh_cut: float,
        add_moradi_bh_systematic: bool) -> np.ndarray:
    """Canonical uncorrelated fit errors for each published measurement."""
    y = data["xs"].to_numpy(float)

    if dataset_kind == "jo2015":
        # Combined/global fits treat the quoted 5% elastic-renormalization
        # uncertainty as correlated.  Therefore use the Jo published total
        # error with that contribution removed in quadrature.  The full
        # ``clas6_err_total`` remains intact for the standalone Moradi
        # reproduction, where the published total error is used directly.
        if "clas6_err_pointwise" in data.columns:
            base = data["clas6_err_pointwise"].to_numpy(float)
        elif "clas6_err_total" in data.columns:
            base = data["clas6_err_total"].to_numpy(float)
        elif "xs_stat" in data.columns:
            base = data["xs_stat"].to_numpy(float)
        elif "err" in data.columns:
            # Defensive compatibility with any older externally prepared
            # Jo dataframe.
            base = data["err"].to_numpy(float)
        else:
            raise KeyError(
                "Jo 2015 dataframe has no published-total-error column; "
                "expected clas6_err_total or xs_stat."
            )
        #endif
    elif dataset_kind == "saylor2018":
        stat = dataset_statistical_errors(data, dataset_kind)
        if "ptp_sys" not in data.columns:
            raise KeyError(
                "Saylor 2018 dataframe is missing canonical column 'ptp_sys'."
            )
        #endif
        ptp = data["ptp_sys"].to_numpy(float)
        base = np.sqrt(stat**2 + ptp**2)
    elif dataset_kind == "georges2022":
        stat = dataset_statistical_errors(data, dataset_kind)
        if "ptp_sys_abs" not in data.columns:
            raise KeyError("Hall A Georges dataframe is missing 'ptp_sys_abs'.")
        #endif
        ptp = data["ptp_sys_abs"].to_numpy(float)
        base = np.sqrt(stat**2 + ptp**2)
    elif dataset_kind in ("halla_defurne2015", "halla_defurne2017"):
        stat = dataset_statistical_errors(data, dataset_kind)
        if "halla_err_syst" in data.columns:
            ptp = data["halla_err_syst"].to_numpy(float)
        elif "ptp_sys_abs" in data.columns:
            ptp = data["ptp_sys_abs"].to_numpy(float)
        else:
            raise KeyError(
                "Hall A Defurne dataframe is missing 'halla_err_syst'."
            )
        #endif
        base = np.sqrt(stat**2 + ptp**2)
    elif dataset_kind in ("pass1", "pass2"):
        stat = dataset_statistical_errors(data, dataset_kind)
        if "ptp_sys_abs" in data.columns:
            ptp = data["ptp_sys_abs"].to_numpy(float)
        elif "ptp_sys" in data.columns:
            ptp = data["ptp_sys"].to_numpy(float)
        else:
            label = "CLAS12 Lee 2026" if dataset_kind == "pass1" else "CLAS12 pass-2"
            raise KeyError(
                f"{label} dataframe is missing point-to-point "
                "systematic column 'ptp_sys_abs'."
            )
        #endif
        base = np.sqrt(stat**2 + ptp**2)
    else:
        raise ValueError(f"unknown dataset kind: {dataset_kind}")
    #endif

    if add_moradi_bh_systematic:
        base = np.sqrt(base**2 + (float(bh_cut) * y)**2)
    #endif
    return base
#enddef


def fit_multi_measurements(
        datasets: Sequence[Dict[str, object]],
        kind: str,
        fit_name: str,
        bh_cut: float = 0.05,
        add_moradi_bh_systematic: bool = True,
        unconstrained_norm_keys: Optional[Sequence[str]] = None) -> FitResult:
    """
    Simultaneous common-form-factor fit to any subset of Jo 2015,
    Saylor 2018, and CLAS12 Lee 2026.

    Each measurement has its own cross-section data and uncertainty model.
    Correlated normalization nuisances are independent between experiments:
      Jo 2015: 5% global-renormalization nuisance in combined fits;
      Saylor 2018: 5% global-normalization nuisance;
      Hall A Defurne 2015: 2.8% global-normalization nuisance;
      Hall A Defurne 2017: 2.291% correlated normalization nuisance;
      CLAS12 Lee 2026: 31% conservative global normalization nuisance.
    """
    names, p0 = paper_model_setup(kind)

    unconstrained_norm_keys = set(
        str(k) for k in (unconstrained_norm_keys or [])
    )
    # Measurement specs can carry the production normalization prescription
    # directly. Explicit function arguments are unioned with those flags.
    unconstrained_norm_keys.update(
        str(spec["key"]) for spec in datasets
        if bool(spec.get("unconstrained_norm", False))
    )
    nuisance_names = []
    nuisance_fracs = {}
    nuisance_is_free = {}
    for spec in datasets:
        key = str(spec["key"])
        free_norm = key in unconstrained_norm_keys
        frac = 1.0 if free_norm else float(spec.get("norm_frac", 0.0))
        if frac > 0.0:
            nname = "beta_" + re.sub(r"[^A-Za-z0-9]+", "_", key)
            nuisance_names.append(nname)
            nuisance_fracs[nname] = frac
            nuisance_is_free[nname] = bool(free_norm)
        #endif
    #endfor

    fit_names = list(names) + nuisance_names
    fit_p0 = np.concatenate([p0, np.zeros(len(nuisance_names), dtype=float)])

    prepared = []
    for spec in datasets:
        d = spec["data"]

        required_common = {"t_abs", "xs", "bh_A", "bh_B", "bh_C"}
        missing_common = sorted(required_common.difference(d.columns))
        if missing_common:
            raise KeyError(
                f"{spec['label']} is missing required combined-fit columns: "
                + ", ".join(missing_common)
            )
        #endif

        point_errors = dataset_point_errors(
            d, str(spec["kind"]), bh_cut, add_moradi_bh_systematic
        )
        if (
            len(point_errors) != len(d)
            or not np.all(np.isfinite(point_errors))
            or np.any(point_errors <= 0.0)
        ):
            raise ValueError(
                f"{spec['label']} produced invalid combined-fit point errors: "
                f"Ndata={len(d)}, Nerr={len(point_errors)}, "
                f"Nnonfinite={int(np.sum(~np.isfinite(point_errors)))}, "
                f"Nnonpositive={int(np.sum(point_errors <= 0.0))}."
            )
        #endif

        key = str(spec["key"])
        effective_norm_frac = (
            1.0 if key in unconstrained_norm_keys
            else float(spec.get("norm_frac", 0.0))
        )
        prepared.append({
            "key": key,
            "kind": spec["kind"],
            "norm_frac": effective_norm_frac,
            "q": d["t_abs"].to_numpy(float),
            "y": d["xs"].to_numpy(float),
            "e": point_errors,
            "A": d["bh_A"].to_numpy(float),
            "B": d["bh_B"].to_numpy(float),
            "C": d["bh_C"].to_numpy(float),
            "N": len(d),
        })
    #endfor

    nuisance_index = {
        name: len(names) + i for i, name in enumerate(nuisance_names)
    }

    def chi2_minuit(*values):
        p = np.asarray(values, dtype=float)
        ff = p[:len(names)]
        total = 0.0

        for item in prepared:
            f1, f2 = paper_model_f1f2(kind, item["q"], ff)
            pred = bh_from_f1f2(item["A"], item["B"], item["C"], f1, f2)

            frac = item["norm_frac"]
            if frac > 0.0:
                nname = "beta_" + re.sub(
                    r"[^A-Za-z0-9]+", "_", str(item["key"])
                )
                beta = p[nuisance_index[nname]]
                pred = pred * (1.0 + frac * beta)
            #endif

            pulls = (pred - item["y"]) / item["e"]
            total += float(np.dot(pulls, pulls))
        #endfor

        for nname in nuisance_names:
            if not nuisance_is_free.get(nname, False):
                beta = p[nuisance_index[nname]]
                total += float(beta**2)
            #endif
        #endfor
        return total
    #enddef

    m = Minuit(chi2_minuit, *fit_p0, name=tuple(fit_names))
    m.errordef = Minuit.LEAST_SQUARES
    for name in names:
        m.limits[name] = (1.0e-6, None)
    #endfor
    for nname in nuisance_names:
        if nuisance_is_free.get(nname, False):
            # For an unconstrained diagnostic nuisance frac=1, so beta is
            # directly the fractional normalization shift. Keep scale positive.
            m.limits[nname] = (-0.50, 0.50)
        else:
            m.limits[nname] = (-10.0, 10.0)
        #endif
    #endfor

    m.migrad()
    m.hesse()

    full = np.array([float(m.values[n]) for n in fit_names])
    params = full[:len(names)]

    full_cov = np.full((len(fit_names), len(fit_names)), np.nan)
    if m.covariance is not None:
        for i, ni in enumerate(fit_names):
            for j, nj in enumerate(fit_names):
                full_cov[i, j] = float(m.covariance[ni, nj])
            #endfor
        #endfor
    #endif
    cov = full_cov[:len(names), :len(names)]

    def re_of_p(p):
        return radius_from_shape(
            lambda qq: f1f2_to_sachs(
                qq, *paper_model_f1f2(kind, qq, p)
            )[0],
            1.0,
        )
    #enddef

    def rm_of_p(p):
        return radius_from_shape(
            lambda qq: f1f2_to_sachs(
                qq, *paper_model_f1f2(kind, qq, p)
            )[1],
            MU_P,
        )
    #enddef

    rE, drE = propagate_scalar(re_of_p, params, cov)
    rM, drM = propagate_scalar(rm_of_p, params, cov)

    nuisance_meta = {}
    for nname in nuisance_names:
        i = nuisance_index[nname]
        beta = float(full[i])
        beta_err = (
            math.sqrt(max(full_cov[i, i], 0.0))
            if np.isfinite(full_cov[i, i])
            else np.nan
        )
        nuisance_meta[nname] = beta
        nuisance_meta[nname + "_err"] = beta_err
        nuisance_meta[nname + "_scale_factor"] = (
            1.0 + nuisance_fracs[nname] * beta
        )
    #endfor

    npts = sum(item["N"] for item in prepared)
    npar = len(fit_names)
    ndof = max(npts - npar, 1)

    labels = [str(spec["label"]) for spec in datasets]
    return FitResult(
        name=fit_name,
        category="multi_measurement",
        bh_cut=bh_cut,
        npts=npts,
        npar=npar,
        chi2=float(m.fval),
        ndof=ndof,
        chi2_ndof=float(m.fval) / ndof,
        success=bool(m.valid),
        message=(
            f"valid={m.valid}; accurate={m.accurate}; "
            f"edm={float(m.fmin.edm):.3e}"
        ),
        params=params,
        param_names=names,
        covariance=cov,
        rE_fm=rE,
        rE_err_fm=drE,
        rM_fm=rM,
        rM_err_fm=drM,
        norm=1.0,
        norm_err=0.0,
        model_kind=kind,
        meta={
            "dataset": " + ".join(labels),
            "dataset_labels": labels,
            "dataset_keys": [str(spec["key"]) for spec in datasets],
            "add_moradi_bh_systematic": bool(add_moradi_bh_systematic),
            **nuisance_meta,
        },
    )
#enddef


def measurement_spec(
        key: str,
        label: str,
        kind: str,
        data: pd.DataFrame,
        norm_frac: float = 0.0,
        unconstrained_norm: bool = False) -> Dict[str, object]:
    """
    Canonical measurement specification used by combined/global fits.

    ``unconstrained_norm`` is optional for backward compatibility.  When True,
    the corresponding fit machinery treats the dataset scale as a free
    normalization rather than assigning a Gaussian normalization penalty.
    """
    return {
        "key": key,
        "label": label,
        "kind": kind,
        "data": data,
        "norm_frac": float(norm_frac),
        "unconstrained_norm": bool(unconstrained_norm),
    }
#enddef


def fit_combined_clas6_pass1(
        clas6: pd.DataFrame,
        pass1: pd.DataFrame,
        kind: str,
        fit_name: str,
        bh_cut: float = 0.05) -> FitResult:
    """
    Simultaneous fit with common F1/F2 shape parameters.

    CLAS6:
      published total error + BH-selection systematic; no normalization nuisance.

    CLAS12 Lee 2026:
      stat + point-to-point systematic + BH-selection systematic, plus its
      independent 31% correlated global normalization nuisance.

    The form-factor parameters are common to both datasets.
    """
    names, p0 = paper_model_setup(kind)
    fit_names = list(names) + ["beta_pass1"]
    fit_p0 = np.concatenate([p0, [0.0]])

    q6 = clas6["t_abs"].to_numpy(float)
    y6 = clas6["xs"].to_numpy(float)
    e6 = fit_sigma_errors(
        clas6, bh_cut,
        include_clas12_ptp_sys=False,
        include_bh_sys=True,
    )
    A6 = clas6["bh_A"].to_numpy(float)
    B6 = clas6["bh_B"].to_numpy(float)
    C6 = clas6["bh_C"].to_numpy(float)

    q1 = pass1["t_abs"].to_numpy(float)
    y1 = pass1["xs"].to_numpy(float)
    e1 = fit_sigma_errors(
        pass1, bh_cut,
        include_clas12_ptp_sys=True,
        include_bh_sys=True,
    )
    A1 = pass1["bh_A"].to_numpy(float)
    B1 = pass1["bh_B"].to_numpy(float)
    C1 = pass1["bh_C"].to_numpy(float)

    def chi2_minuit(*values):
        p = np.asarray(values, dtype=float)
        ff = p[:len(names)]
        beta = float(p[-1])

        f16, f26 = paper_model_f1f2(kind, q6, ff)
        pred6 = bh_from_f1f2(A6, B6, C6, f16, f26)

        f11, f21 = paper_model_f1f2(kind, q1, ff)
        pred1 = bh_from_f1f2(A1, B1, C1, f11, f21)
        pred1 *= 1.0 + PASS1_GLOBAL_SCALE_FRAC * beta

        pull6 = (pred6 - y6) / e6
        pull1 = (pred1 - y1) / e1
        return float(
            np.dot(pull6, pull6)
            + np.dot(pull1, pull1)
            + beta**2
        )
    #enddef

    m = Minuit(chi2_minuit, *fit_p0, name=tuple(fit_names))
    m.errordef = Minuit.LEAST_SQUARES
    for name in names:
        m.limits[name] = (1.0e-6, None)
    #endfor
    m.limits["beta_pass1"] = (-10.0, 10.0)
    m.migrad()
    m.hesse()

    full = np.array([float(m.values[n]) for n in fit_names])
    params = full[:len(names)]
    beta = float(full[-1])

    full_cov = np.full((len(fit_names), len(fit_names)), np.nan)
    if m.covariance is not None:
        for i, ni in enumerate(fit_names):
            for j, nj in enumerate(fit_names):
                full_cov[i, j] = float(m.covariance[ni, nj])
            #endfor
        #endfor
    #endif
    cov = full_cov[:len(names), :len(names)]

    def re_of_p(p):
        return radius_from_shape(
            lambda qq: f1f2_to_sachs(
                qq, *paper_model_f1f2(kind, qq, p)
            )[0],
            1.0,
        )
    #enddef
    def rm_of_p(p):
        return radius_from_shape(
            lambda qq: f1f2_to_sachs(
                qq, *paper_model_f1f2(kind, qq, p)
            )[1],
            MU_P,
        )
    #enddef

    rE, drE = propagate_scalar(re_of_p, params, cov)
    rM, drM = propagate_scalar(rm_of_p, params, cov)

    beta_err = (
        math.sqrt(max(full_cov[-1, -1], 0.0))
        if np.all(np.isfinite(full_cov))
        else np.nan
    )

    npts = len(clas6) + len(pass1)
    npar = len(fit_names)
    ndof = max(npts - npar, 1)

    return FitResult(
        name=fit_name,
        category="combined_CLAS6_pass1",
        bh_cut=bh_cut,
        npts=npts,
        npar=npar,
        chi2=float(m.fval),
        ndof=ndof,
        chi2_ndof=float(m.fval) / ndof,
        success=bool(m.valid),
        message=(
            f"valid={m.valid}; accurate={m.accurate}; "
            f"edm={float(m.fmin.edm):.3e}"
        ),
        params=params,
        param_names=names,
        covariance=cov,
        rE_fm=rE,
        rE_err_fm=drE,
        rM_fm=rM,
        rM_err_fm=drM,
        norm=1.0 + PASS1_GLOBAL_SCALE_FRAC * beta,
        norm_err=PASS1_GLOBAL_SCALE_FRAC * beta_err,
        model_kind=kind,
        meta={
            "dataset": "CLAS6 + CLAS12 Lee 2026",
            "N_CLAS6": len(clas6),
            "N_pass1": len(pass1),
            "beta_pass1": beta,
            "beta_pass1_err": beta_err,
            "global_scale_frac": PASS1_GLOBAL_SCALE_FRAC,
            "include_global_norm_sys": False,
            "include_combination_norm_sys": False,
        },
    )
#enddef


def evaluate_fit_form_factors(fr: FitResult,
                              q: np.ndarray) -> Tuple[np.ndarray, np.ndarray,
                                                      np.ndarray, np.ndarray]:
    """Return F1, F2, GE, GM for a paper-method FitResult."""
    f1, f2 = paper_model_f1f2(fr.model_kind, q, fr.params)
    ge, gm = f1f2_to_sachs(q, f1, f2)
    return f1, f2, ge, gm
#enddef


def evaluate_fit_cross_section(fr: FitResult,
                               data: pd.DataFrame,
                               apply_nuisances: bool = True) -> np.ndarray:
    """Evaluate fitted BH at the measured CLAS12 points."""
    q = data["t_abs"].to_numpy(float)
    f1, f2 = paper_model_f1f2(fr.model_kind, q, fr.params)
    bare = bh_from_f1f2(
        data["bh_A"].to_numpy(float),
        data["bh_B"].to_numpy(float),
        data["bh_C"].to_numpy(float),
        f1,
        f2,
    )

    if not apply_nuisances:
        return bare
    #endif

    beta_global = float(fr.meta.get("beta_global", 0.0))
    beta_comb = float(fr.meta.get("beta_comb", 0.0))
    use_global = bool(fr.meta.get("include_global_norm_sys", False))
    use_comb = bool(fr.meta.get("include_combination_norm_sys", False))

    global_factor = (
        1.0 + beta_global * float(fr.meta.get('global_scale_frac', GLOBAL_SCALE_FRAC))
        if use_global
        else 1.0
    )

    comb = data["comb_sys_frac"].to_numpy(float)
    comb = np.where(np.isfinite(comb) & (comb >= 0.0), comb, 0.0)
    combination_factor = (
        1.0 + beta_comb * comb
        if use_comb
        else np.ones_like(comb)
    )

    return bare * global_factor * combination_factor
#enddef


# -----------------------------------------------------------------------------
# Fit bookkeeping and Hessian uncertainty bands
# -----------------------------------------------------------------------------
def fitresult_to_record(fr: FitResult) -> Dict[str, object]:
    rec = {
        "name": fr.name,
        "category": fr.category,
        "bh_cut": fr.bh_cut,
        "npts": fr.npts,
        "npar": fr.npar,
        "chi2": fr.chi2,
        "ndof": fr.ndof,
        "chi2_ndof": fr.chi2_ndof,
        "success": fr.success,
        "message": fr.message,
        "rE_fm": fr.rE_fm,
        "rE_err_fm": fr.rE_err_fm,
        "rM_fm": fr.rM_fm,
        "rM_err_fm": fr.rM_err_fm,
        "norm": fr.norm,
        "norm_err": fr.norm_err,
        "global_norm_factor": fr.norm,
        "global_norm_factor_err": fr.norm_err,
        "beta_global": fr.meta.get("beta_global", np.nan),
        "beta_global_err": fr.meta.get("beta_global_err", np.nan),
        "beta_comb": fr.meta.get("beta_comb", np.nan),
        "beta_comb_err": fr.meta.get("beta_comb_err", np.nan),
        "combination_sys_median": fr.meta.get("combination_sys_median", np.nan),
        "combination_sys_min": fr.meta.get("combination_sys_min", np.nan),
        "combination_sys_max": fr.meta.get("combination_sys_max", np.nan),
    }
    for name, value in zip(fr.param_names, fr.params):
        rec[name] = value
    #endfor
    return rec
#enddef


def form_factor_band(fr: FitResult,
                     q: np.ndarray,
                     which: str) -> Tuple[np.ndarray, np.ndarray]:
    """Central value and 1-sigma Hessian band."""
    index = {"F1": 0, "F2": 1, "GE": 2, "GM": 3}[which]
    central = evaluate_fit_form_factors(fr, q)[index]

    if not np.all(np.isfinite(fr.covariance)):
        return central, np.full_like(central, np.nan)
    #endif

    sigma = np.zeros_like(q, dtype=float)
    for iq, q_value in enumerate(q):
        def scalar(p):
            tmp = FitResult(**{**fr.__dict__, "params": np.asarray(p, dtype=float)})
            return float(
                evaluate_fit_form_factors(
                    tmp, np.array([q_value], dtype=float)
                )[index][0]
            )
        #enddef

        grad = numerical_gradient(scalar, fr.params)
        var = float(grad @ fr.covariance @ grad)
        sigma[iq] = math.sqrt(max(var, 0.0)) if np.isfinite(var) else np.nan
    #endfor

    return central, sigma
#enddef



def derived_form_factor_band(
        fr: FitResult,
        q: np.ndarray,
        quantity: str) -> Tuple[np.ndarray, np.ndarray]:
    """Propagated 1-sigma Hessian band for common normalized FF ratios."""
    q = np.asarray(q, dtype=float)

    def evaluate_with_params(params: np.ndarray) -> np.ndarray:
        tmp = FitResult(**{**fr.__dict__, "params": np.asarray(params, dtype=float)})
        _, _, ge, gm = evaluate_fit_form_factors(tmp, q)
        gd = (1.0 + q / 0.71) ** -2

        if quantity == "GE_over_GD":
            return ge / gd
        #endif
        if quantity == "GM_over_muGD":
            return gm / (MU_P * gd)
        #endif
        if quantity == "muGE_over_GM":
            return MU_P * ge / gm
        #endif
        raise ValueError(quantity)
    #enddef

    central = evaluate_with_params(fr.params)
    if not np.all(np.isfinite(fr.covariance)):
        return central, np.full_like(central, np.nan)
    #endif

    sigma = np.zeros_like(q)
    for iq in range(len(q)):
        def scalar(p):
            vals = evaluate_with_params(p)
            return float(vals[iq])
        #enddef
        grad = numerical_gradient(scalar, fr.params)
        var = float(grad @ fr.covariance @ grad)
        sigma[iq] = math.sqrt(max(var, 0.0)) if np.isfinite(var) else np.nan
    #endfor
    return central, sigma
#enddef


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------
def draw_t_reference_lines(ax, data: pd.DataFrame) -> None:
    """Draw only the fitted |t|/Q2 limits; ultra-low-Q2 reference lines are omitted."""
    if len(data) == 0:
        return
    #endif

    ax.axvline(
        float(data["t_abs"].min()),
        linewidth=0.8,
        linestyle="-",
        label="fitted |t| limits",
    )
    ax.axvline(
        float(data["t_abs"].max()),
        linewidth=0.8,
        linestyle="-",
    )
#enddef


def save_bh_selection_plots(df: pd.DataFrame,
                            bh_cuts: Sequence[float],
                            outdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    vals = df["R_BH"].to_numpy(float)
    vals = vals[np.isfinite(vals)]
    ax.hist(vals, bins=80, histtype="step", linewidth=1.5)

    for cut in bh_cuts:
        ax.axvline(
            1.0 - cut,
            linestyle="--",
            linewidth=1.0,
            label=rf"$\Delta\sigma/\sigma\leq{100*cut:.0f}\%$",
        )
        ax.axvline(1.0 + cut, linestyle="--", linewidth=1.0)
    #endfor

    ax.set_xlabel(r"$R_{\rm BH}=\sigma_{\rm BH}/\sigma_{\rm EP}^{\rm KM15}$")
    ax.set_ylabel("CLAS12 cross-section points")
    ax.set_title("KM15 BH-dominance selection")
    ax.legend(fontsize=7.8, ncol=2)
    fig.tight_layout()
    fig.savefig(outdir / "01_bh_fraction_distribution.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.0, 5.4))
    sc = ax.scatter(
        df["phi_deg"],
        df["xB"],
        c=df["R_BH"],
        s=15,
        vmin=0.8,
        vmax=1.05,
    )
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label(r"$R_{\rm BH}$")
    ax.set_xlabel(r"$\phi$ (deg)")
    ax.set_ylabel(r"$x_B$")
    ax.set_title("Preliminary CLAS12 kinematics colored by KM15 BH fraction")
    fig.tight_layout()
    fig.savefig(outdir / "02_bh_fraction_phi_xB.png")
    plt.close(fig)
#enddef


def save_selected_cross_section_pages(data: pd.DataFrame,
                                      fit5: FitResult,
                                      fit8: FitResult,
                                      outdir: Path) -> None:
    """
    Save one PNG per measured (xB,Q2,|t|) group.

    Upper panel:
      - filled points: measured preliminary CLAS12 cross sections
      - open markers: fitted BH predictions evaluated at those exact measured
        CLAS12 points only
      - no connecting theory lines are drawn

    Lower panel:
      - point-by-point pulls using exactly the uncertainty entering chi2

    The continuous fitted quantities are F1(t), F2(t), GE(t), and GM(t);
    those are shown as smooth analytic curves in the form-factor figures.
    """
    plotdir = outdir / "03_cross_section_fits"
    plotdir.mkdir(parents=True, exist_ok=True)

    work = data.copy()
    work["fit5_bh"] = evaluate_fit_cross_section(fit5, work, True)
    work["fit8_bh"] = evaluate_fit_cross_section(fit8, work, True)
    work["fit5_bh_bare"] = evaluate_fit_cross_section(fit5, work, False)
    work["fit8_bh_bare"] = evaluate_fit_cross_section(fit8, work, False)

    work["fit5_pull"] = (
        (work["fit5_bh"] - work["xs"]) / work["fit_sigma_default"]
    )
    work["fit8_pull"] = (
        (work["fit8_bh"] - work["xs"]) / work["fit_sigma_default"]
    )
    work["fit5_frac_residual"] = (
        (work["fit5_bh"] - work["xs"]) / work["xs"]
    )
    work["fit8_frac_residual"] = (
        (work["fit8_bh"] - work["xs"]) / work["xs"]
    )

    group_cols = [
        c for c in
        ["xBmin", "xBmax", "Q2min", "Q2max", "t_abs_min", "t_abs_max"]
        if c in work.columns
    ]
    if len(group_cols) < 6:
        work["_gx"] = work["xB"].round(3)
        work["_gq"] = work["Q2"].round(2)
        work["_gt"] = work["t_abs"].round(3)
        group_cols = ["_gx", "_gq", "_gt"]
    #endif

    groups = [
        (key, group.sort_values("phi_deg"))
        for key, group in work.groupby(group_cols, dropna=False)
    ]

    for igroup, (_, group) in enumerate(groups, start=1):
        fig = plt.figure(figsize=(8.0, 7.0))
        gs = fig.add_gridspec(
            2, 1,
            height_ratios=[3.0, 1.15],
            hspace=0.05,
        )
        ax = fig.add_subplot(gs[0])
        pax = fig.add_subplot(gs[1], sharex=ax)

        ax.errorbar(
            group["phi_deg"],
            group["xs"],
            yerr=group["fit_sigma_default"],
            fmt="o",
            markersize=4.0,
            capsize=2,
            label="CLAS12 preliminary",
        )
        # The cross-section fit exists only at the measured CLAS12 points.
        # Do NOT connect these predictions with straight segments: that would
        # visually imply a continuous phi-dependent fit curve that was never
        # evaluated.  The continuous fitted objects are F1(t) and F2(t), shown
        # separately in the form-factor plots.
        ax.plot(
            group["phi_deg"],
            group["fit5_bh"],
            marker="o",
            linestyle="none",
            markersize=5.0,
            markerfacecolor="none",
            markeredgewidth=1.2,
            label="Fit 5 BH at CLAS12 points",
        )
        ax.plot(
            group["phi_deg"],
            group["fit8_bh"],
            marker="s",
            linestyle="none",
            markersize=4.8,
            markerfacecolor="none",
            markeredgewidth=1.2,
            label=r"Fit 8 BH at CLAS12 points ($F_2$ Kelly)",
        )
        ax.plot(
            group["phi_deg"],
            group["km15_bh"],
            marker="^",
            linestyle="none",
            markersize=4.6,
            markerfacecolor="none",
            markeredgewidth=1.0,
            label="KM15 nominal BH at CLAS12 points",
        )

        pax.axhline(0.0, linewidth=0.8)
        for level in (-3.0, -2.0, -1.0, 1.0, 2.0, 3.0):
            pax.axhline(level, linewidth=0.5, linestyle="--", alpha=0.5)
        #endfor

        pax.plot(
            group["phi_deg"],
            group["fit5_pull"],
            "o-",
            markersize=3.5,
            linewidth=0.9,
            label="Fit 5 pull",
        )
        pax.plot(
            group["phi_deg"],
            group["fit8_pull"],
            "s--",
            markersize=3.2,
            linewidth=0.9,
            label="Fit 8 pull",
        )

        all_pulls = np.concatenate([
            group["fit5_pull"].to_numpy(float),
            group["fit8_pull"].to_numpy(float),
        ])
        max_pull = max(
            4.0,
            float(np.nanmax(np.abs(all_pulls))) + 0.5,
        )
        pax.set_ylim(-max_pull, max_pull)

        ax.set_ylabel(r"$d^4\sigma$")
        pax.set_ylabel("pull")
        pax.set_xlabel(r"$\phi$ (deg)")
        ax.grid(alpha=0.2)
        pax.grid(alpha=0.2)

        ax.set_title(
            rf"$\langle x_B\rangle={group['xB'].mean():.3f}$, "
            rf"$\langle Q^2\rangle={group['Q2'].mean():.2f}$ GeV$^2$, "
            rf"$\langle|t|\rangle={group['t_abs'].mean():.3f}$ GeV$^2$"
        )
        ax.legend(fontsize=8)
        pax.legend(fontsize=8, ncol=2)

        plt.setp(ax.get_xticklabels(), visible=False)
        fig.tight_layout()
        fig.savefig(
            plotdir / f"cross_section_fit_{igroup:03d}.png",
            dpi=180,
        )
        plt.close(fig)
    #endfor

    columns = [
        "_row", "xB", "Q2", "t_abs", "phi_deg", "xs", "xs_stat",
        "ptp_sys_abs", "comb_sys_frac", "fit_sigma_default",
        "fit5_bh", "fit5_bh_bare", "fit5_pull", "fit5_frac_residual",
        "fit8_bh", "fit8_bh_bare", "fit8_pull", "fit8_frac_residual",
        "R_BH", "bh_delta",
    ]
    columns = [c for c in columns if c in work.columns]
    outliers = work[columns].copy()
    outliers["abs_fit5_pull"] = np.abs(work["fit5_pull"].to_numpy(float))
    outliers = outliers.sort_values("abs_fit5_pull", ascending=False)
    outliers.to_csv(plotdir / "fit5_points_sorted_by_abs_pull.csv", index=False)

    variables = [
        ("phi_deg", r"$\phi$ (deg)"),
        ("t_abs", r"$|t|$ (GeV$^2$)"),
        ("xB", r"$x_B$"),
        ("Q2", r"$Q^2$ (GeV$^2$)"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for ax, (column, xlabel) in zip(axes.flat, variables):
        ax.scatter(work[column], work["fit5_pull"], s=14)
        ax.axhline(0.0, linewidth=0.8)
        ax.axhline(3.0, linewidth=0.6, linestyle="--")
        ax.axhline(-3.0, linewidth=0.6, linestyle="--")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Fit 5 pull")
        ax.grid(alpha=0.2)
    #endfor
    fig.tight_layout()
    fig.savefig(plotdir / "fit5_pulls_vs_kinematics.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    ax.hist(work["fit5_pull"], bins=35, histtype="step", linewidth=1.4)
    ax.axvline(0.0, linewidth=0.8)
    ax.set_xlabel("Fit 5 pull")
    ax.set_ylabel("CLAS12 points")
    ax.set_title(
        rf"Fit 5 pulls: RMS={np.sqrt(np.mean(work['fit5_pull']**2)):.2f}"
    )
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(plotdir / "fit5_pull_distribution.png", dpi=180)
    plt.close(fig)
#enddef



def save_pass1_cross_section_diagnostics(data: pd.DataFrame,
                                         fit5: FitResult,
                                         fit8: FitResult,
                                         outdir: Path) -> None:
    """Per-bin pass-1 data-vs-fit PNGs plus pull summaries."""
    plotdir = outdir / "03_cross_section_fits"
    plotdir.mkdir(parents=True, exist_ok=True)

    work = data.copy()
    work["fit5_bh"] = evaluate_fit_cross_section(fit5, work, True)
    work["fit8_bh"] = evaluate_fit_cross_section(fit8, work, True)
    work["fit5_pull"] = (
        (work["fit5_bh"] - work["xs"]) / work["fit_sigma_default"]
    )
    work["fit8_pull"] = (
        (work["fit8_bh"] - work["xs"]) / work["fit_sigma_default"]
    )

    group_cols = [
        c for c in
        ["xBmin", "xBmax", "Q2min", "Q2max", "t_abs_min", "t_abs_max"]
        if c in work.columns
    ]
    if len(group_cols) < 6:
        work["_gx"] = work["xB"].round(3)
        work["_gq"] = work["Q2"].round(2)
        work["_gt"] = work["t_abs"].round(3)
        group_cols = ["_gx", "_gq", "_gt"]
    #endif

    groups = [
        (key, group.sort_values("phi_deg"))
        for key, group in work.groupby(group_cols, dropna=False)
    ]

    for igroup, (_, group) in enumerate(groups, start=1):
        fig = plt.figure(figsize=(8.0, 7.0))
        gs = fig.add_gridspec(
            2, 1,
            height_ratios=[3.0, 1.15],
            hspace=0.05,
        )
        ax = fig.add_subplot(gs[0])
        pax = fig.add_subplot(gs[1], sharex=ax)

        ax.errorbar(
            group["phi_deg"],
            group["xs"],
            yerr=group["fit_sigma_default"],
            fmt="o",
            markersize=4.0,
            capsize=2,
            label="CLAS12 Lee 2026",
        )
        ax.plot(
            group["phi_deg"],
            group["fit5_bh"],
            marker="o",
            linestyle="none",
            markersize=5.0,
            markerfacecolor="none",
            label="Fit 5 BH at pass-1 points",
        )
        ax.plot(
            group["phi_deg"],
            group["fit8_bh"],
            marker="s",
            linestyle="none",
            markersize=4.8,
            markerfacecolor="none",
            label=r"Fit 8 BH at pass-1 points ($F_2$ Kelly)",
        )
        ax.plot(
            group["phi_deg"],
            group["km15_bh"],
            marker="^",
            linestyle="none",
            markersize=4.6,
            markerfacecolor="none",
            label="KM15 nominal BH",
        )

        pax.axhline(0.0, linewidth=0.8)
        for level in (-3.0, -2.0, -1.0, 1.0, 2.0, 3.0):
            pax.axhline(level, linewidth=0.5, linestyle="--", alpha=0.5)
        #endfor

        pax.plot(
            group["phi_deg"],
            group["fit5_pull"],
            "o-",
            markersize=3.5,
            linewidth=0.9,
            label="Fit 5 pull",
        )
        pax.plot(
            group["phi_deg"],
            group["fit8_pull"],
            "s--",
            markersize=3.2,
            linewidth=0.9,
            label="Fit 8 pull",
        )

        all_pulls = np.concatenate([
            group["fit5_pull"].to_numpy(float),
            group["fit8_pull"].to_numpy(float),
        ])
        max_pull = max(
            4.0,
            float(np.nanmax(np.abs(all_pulls))) + 0.5,
        )
        pax.set_ylim(-max_pull, max_pull)

        ax.set_ylabel(r"$d^4\sigma$")
        pax.set_ylabel("pull")
        pax.set_xlabel(r"$\phi$ (deg)")
        ax.grid(alpha=0.2)
        pax.grid(alpha=0.2)
        ax.set_title(
            rf"Pass 1: $\langle x_B\rangle={group['xB'].mean():.3f}$, "
            rf"$\langle Q^2\rangle={group['Q2'].mean():.2f}$ GeV$^2$, "
            rf"$\langle|t|\rangle={group['t_abs'].mean():.3f}$ GeV$^2$"
        )
        ax.legend(fontsize=8)
        pax.legend(fontsize=8, ncol=2)
        plt.setp(ax.get_xticklabels(), visible=False)

        fig.tight_layout()
        fig.savefig(
            plotdir / f"pass1_cross_section_fit_{igroup:03d}.png",
            dpi=180,
        )
        plt.close(fig)
    #endfor

    columns = [
        "_row", "xB", "Q2", "t_abs", "phi_deg", "xs", "xs_stat",
        "ptp_sys_abs", "fit_sigma_default",
        "fit5_bh", "fit5_pull",
        "fit8_bh", "fit8_pull",
        "R_BH", "bh_delta",
    ]
    columns = [c for c in columns if c in work.columns]
    outliers = work[columns].copy()
    outliers["abs_fit5_pull"] = np.abs(work["fit5_pull"].to_numpy(float))
    outliers.sort_values(
        "abs_fit5_pull",
        ascending=False,
    ).to_csv(
        plotdir / "pass1_fit5_points_sorted_by_abs_pull.csv",
        index=False,
    )

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for ax, (column, xlabel) in zip(
        axes.flat,
        [
            ("phi_deg", r"$\phi$ (deg)"),
            ("t_abs", r"$|t|$ (GeV$^2$)"),
            ("xB", r"$x_B$"),
            ("Q2", r"$Q^2$ (GeV$^2$)"),
        ],
    ):
        ax.scatter(work[column], work["fit5_pull"], s=14)
        ax.axhline(0.0, linewidth=0.8)
        ax.axhline(3.0, linewidth=0.6, linestyle="--")
        ax.axhline(-3.0, linewidth=0.6, linestyle="--")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Fit 5 pull")
        ax.grid(alpha=0.2)
    #endfor
    fig.tight_layout()
    fig.savefig(plotdir / "pass1_fit5_pulls_vs_kinematics.png", dpi=180)
    plt.close(fig)
#enddef


def save_fit1_to_fit5_plots(results: List[FitResult],
                            set5_data: pd.DataFrame,
                            outdir: Path) -> None:
    """Combined paper Fig. 4/5 analogue: F1 and F2 on one 1x2 canvas."""
    fits = [
        next(r for r in results if r.name == f"Fit {i}")
        for i in range(1, 6)
    ]
    fit5 = fits[-1]
    q = np.linspace(0.0, 1.0, 400)

    # Remove obsolete separate-panel outputs from older script versions.
    for obsolete in [
        outdir / "04_F1_fits_1_to_5.png",
        outdir / "05_F2_fits_1_to_5.png",
    ]:
        if obsolete.exists():
            obsolete.unlink()
        #endif
    #endfor

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.0))

    for ax, which, ylabel, index in [
        (axes[0], "F1", r"$F_1(t)$", 0),
        (axes[1], "F2", r"$F_2(t)$", 1),
    ]:
        for fr in fits:
            values = evaluate_fit_form_factors(fr, q)[index]
            ax.plot(q, values, linewidth=1.2, label=fr.name)
        #endfor

        central, sigma = form_factor_band(fit5, q, which)
        ax.fill_between(
            q,
            central - sigma,
            central + sigma,
            alpha=0.20,
            label="Fit 5 68% Hessian band",
        )

        draw_t_reference_lines(ax, set5_data)
        ax.set_xlabel(r"$|t|$ (GeV$^2$)")
        ax.set_ylabel(ylabel)
        ax.set_xlim(0.0, 1.0)
        ax.grid(alpha=0.2)
    #endfor

    axes[0].legend(fontsize=8)
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "04_F1_F2_fits_1_to_5.png", dpi=180)
    plt.close(fig)
#enddef


def save_fit5_sachs_plot(fit5: FitResult,
                         set5_data: pd.DataFrame,
                         outdir: Path) -> None:
    """Paper Fig. 6 analogue using CLAS12 Fit 5."""
    q = np.linspace(0.0, 1.0, 400)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))

    for ax, which, ylabel in zip(
        axes,
        ["GE", "GM"],
        [r"$G_E(t)$", r"$G_M(t)$"],
    ):
        central, sigma = form_factor_band(fit5, q, which)
        ax.plot(q, central, linewidth=1.5, label="Fit 5")
        ax.fill_between(
            q,
            central - sigma,
            central + sigma,
            alpha=0.20,
            label="68% Hessian band",
        )
        draw_t_reference_lines(ax, set5_data)
        ax.set_xlabel(r"$|t|$ (GeV$^2$)")
        ax.set_ylabel(ylabel)
        ax.set_xlim(0.0, 1.0)
        ax.grid(alpha=0.2)
    #endfor

    axes[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "06_GE_GM_fit5.png")
    plt.close(fig)
#enddef



def save_low_q2_ratio_plots(
        results: List[FitResult],
        set5_data: pd.DataFrame,
        outdir: Path) -> None:
    """
    Low-Q2 electric/magnetic shape diagnostics.

    mu_p GE/GM = 1 at Q2=0 and its initial slope is proportional to
    r_M^2-r_E^2.  GE/GD and GM/(mu_p GD) expose deviations from common
    dipole behavior.
    """
    qmax = max(0.60, float(set5_data["t_abs"].max()) * 1.10)
    q = np.linspace(0.0, qmax, 500)
    gd = (1.0 + q / 0.71) ** -2

    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.9))
    configs = [
        ("GE_over_GD", r"$G_E/G_D$"),
        ("GM_over_muGD", r"$G_M/(\mu_pG_D)$"),
        ("muGE_over_GM", r"$\mu_pG_E/G_M$"),
    ]

    for fit_name in ["Fit 5", "Fit 8"]:
        fr = next(r for r in results if r.name == fit_name)
        for ax, (quantity, ylabel) in zip(axes, configs):
            central, sigma = derived_form_factor_band(fr, q, quantity)
            line, = ax.plot(q, central, linewidth=1.5, label=fit_name)
            ax.fill_between(
                q, central - sigma, central + sigma,
                alpha=0.16, color=line.get_color()
            )
            ax.set_ylabel(ylabel)
        #endfor
    #endfor

    for label, (ge, gm) in elastic_reference_curves(q).items():
        axes[0].plot(q, ge / gd, linewidth=1.0, linestyle="--", label=label)
        axes[1].plot(q, gm / (MU_P * gd), linewidth=1.0, linestyle="--")
        axes[2].plot(q, MU_P * ge / gm, linewidth=1.0, linestyle="--")
    #endfor

    for ax in axes:
        draw_t_reference_lines(ax, set5_data)
        ax.axhline(1.0, linewidth=0.8, linestyle=":")
        ax.set_xlabel(r"$Q^2=|t|$ (GeV$^2$)")
        ax.set_xlim(0.0, qmax)
        ax.grid(alpha=0.2)
    #endfor

    axes[0].legend(fontsize=7)
    fig.suptitle(
        r"Low-$Q^2$ electric/magnetic shape: "
        r"$\mu_pG_E/G_M=1+\frac{r_M^2-r_E^2}{6(\hbar c)^2}Q^2+\cdots$",
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.savefig(outdir / "09_lowQ2_GE_GM_ratios.png", dpi=180)
    plt.close(fig)
#enddef


def save_elastic_reference_comparison(
        results: List[FitResult],
        set5_data: pd.DataFrame,
        outdir: Path,
        fit_name: str = "Fit 5") -> None:
    """Compare BH-extracted GE/GM directly with elastic fits and A1 data."""
    fr = next(r for r in results if r.name == fit_name)
    qmax = max(0.60, float(set5_data["t_abs"].max()) * 1.10)
    q = np.linspace(0.0, qmax, 500)
    a1 = bernauer_rosenbluth_data()
    a1 = a1.loc[a1["Q2"] <= qmax].copy()

    fig, axes = plt.subplots(
        2, 2, figsize=(12.0, 8.5),
        gridspec_kw={"height_ratios": [2.2, 1.0]},
        sharex="col",
    )

    for j, (which, col, ylabel) in enumerate([
        ("GE", "GE", r"$G_E$"),
        ("GM", "GM", r"$G_M$"),
    ]):
        central, sigma = form_factor_band(fr, q, which)
        line, = axes[0, j].plot(q, central, linewidth=1.6, label=f"BH {fit_name}")
        axes[0, j].fill_between(
            q, central - sigma, central + sigma,
            alpha=0.18, color=line.get_color(),
            label="68% Hessian band",
        )

        for label, (ge_ref, gm_ref) in elastic_reference_curves(q).items():
            vals = ge_ref if which == "GE" else gm_ref
            axes[0, j].plot(q, vals, linewidth=1.0, linestyle="--", label=label)
        #endfor

        axes[0, j].errorbar(
            a1["Q2"], a1[col],
            yerr=a1[f"{col}_err"],
            fmt="o", markersize=3.0, fillstyle="none",
            capsize=1.5, linewidth=0.8,
            label="A1/Bernauer Rosenbluth",
        )

        # Diagnostic pull: A1 point minus BH fit, divided by the quadrature
        # of the tabulated A1 error and our propagated Hessian uncertainty.
        qpts = a1["Q2"].to_numpy(float)
        bh_c, bh_s = form_factor_band(fr, qpts, which)
        denom = np.sqrt(a1[f"{col}_err"].to_numpy(float)**2 + bh_s**2)
        pull = (a1[col].to_numpy(float) - bh_c) / denom
        axes[1, j].axhline(0.0, linewidth=0.8)
        axes[1, j].plot(qpts, pull, "o", markersize=3.0, fillstyle="none")
        axes[1, j].set_ylabel("diagnostic pull")
        axes[1, j].set_xlabel(r"$Q^2=|t|$ (GeV$^2$)")
        axes[1, j].grid(alpha=0.2)

        draw_t_reference_lines(axes[0, j], set5_data)
        axes[0, j].set_ylabel(ylabel)
        axes[0, j].grid(alpha=0.2)
        axes[0, j].set_xlim(0.0, qmax)
    #endfor

    axes[0, 0].legend(fontsize=7)
    fig.suptitle(
        f"{fit_name}: BH extraction vs direct elastic proton form factors",
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    safe = fit_name.lower().replace(" ", "_")
    fig.savefig(outdir / f"10_elastic_reference_{safe}.png", dpi=180)
    plt.close(fig)

    # Save the diagnostic numbers too.
    rows = []
    for which, col in [("GE", "GE"), ("GM", "GM")]:
        qpts = a1["Q2"].to_numpy(float)
        bh_c, bh_s = form_factor_band(fr, qpts, which)
        denom = np.sqrt(a1[f"{col}_err"].to_numpy(float)**2 + bh_s**2)
        for i in range(len(a1)):
            rows.append({
                "fit": fit_name,
                "form_factor": col,
                "Q2": qpts[i],
                "A1_value": float(a1.iloc[i][col]),
                "A1_error": float(a1.iloc[i][f"{col}_err"]),
                "BH_value": float(bh_c[i]),
                "BH_hessian_error": float(bh_s[i]),
                "diagnostic_pull": float(
                    (a1.iloc[i][col] - bh_c[i]) / denom[i]
                ),
            })
        #endfor
    #endfor
    pd.DataFrame(rows).to_csv(
        outdir / f"10_elastic_reference_{safe}_pulls.csv",
        index=False,
    )
#enddef




def sachs_first_coefficient_from_radius(radius_fm: float, family: str) -> float:
    """First normalized-shape coefficient corresponding to an RMS radius."""
    slope = -(float(radius_fm) / HBARC)**2 / 6.0
    if family.startswith("D"):
        # Regular dipole: G/G(0) = (1 + a Q2)^(-2), hence slope = -2 a.
        return float(-0.5 * slope)
    #endif
    if family.startswith("P"):
        return float(slope)
    #endif
    # IP and CF both have dG/dQ2|0 = -a1.
    return float(-slope)
#enddef


def apply_sachs_radius_limits(
        minimizer,
        names_e: Sequence[str],
        names_m: Sequence[str],
        family: str) -> None:
    """Apply broad radius-equivalent limits to the first E/M slope terms."""
    c_lo = sachs_first_coefficient_from_radius(SACHS_MIN_RADIUS_FM, family)
    c_hi = sachs_first_coefficient_from_radius(SACHS_MAX_RADIUS_FM, family)
    lo, hi = min(c_lo, c_hi), max(c_lo, c_hi)
    for name in [names_e[0], names_m[0]]:
        minimizer.limits[name] = (lo, hi)
    #endfor
#enddef


def sachs_family_value(q: np.ndarray, coeffs: np.ndarray, family: str) -> np.ndarray:
    """Normalized G(Q2)/G(0) candidate used by the radius-bias study."""
    q = np.asarray(q, dtype=float)
    c = np.asarray(coeffs, dtype=float)

    if family.startswith("D"):
        # One-parameter regular dipole.  D1 is deliberately kept distinct from
        # IP1: the latter is monopole-like, whereas the dipole has a squared
        # denominator and therefore a fixed relation between slope and curvature.
        a = float(c[0])
        return 1.0 / (1.0 + a * q)**2
    #endif

    if family.startswith("P"):
        out = np.ones_like(q)
        for i, a in enumerate(c, start=1):
            out += a * q**i
        #endfor
        return out
    #endif

    if family.startswith("IP"):
        den = np.ones_like(q)
        for i, a in enumerate(c, start=1):
            den += a * q**i
        #endfor
        return 1.0 / den
    #endif

    if family.startswith("CF"):
        # 1 / (1 + a1 q / (1 + a2 q / (...))).
        tail = np.ones_like(q)
        for a in c[::-1]:
            tail = 1.0 + a * q / tail
        #endfor
        return 1.0 / tail
    #endif

    raise ValueError(f"unknown Sachs family {family}")
#enddef


def sachs_family_value_precomputed(q: np.ndarray, coeffs: np.ndarray, family: str,
                                    q_powers: Optional[np.ndarray] = None) -> np.ndarray:
    """Fast repeated P/IP evaluation using invariant precomputed q powers."""
    c = np.asarray(coeffs, dtype=float)
    if q_powers is not None and (family.startswith("P") or family.startswith("IP")):
        poly = 1.0 + np.dot(c, q_powers[:len(c)])
        return 1.0 / poly if family.startswith("IP") else poly
    #endif
    return sachs_family_value(q, c, family)
#enddef


def sachs_family_radius(coeffs: np.ndarray, family: str) -> float:
    """Radius from the Q2=0 slope of a normalized Sachs family."""
    def shape(q):
        return sachs_family_value(np.asarray(q, dtype=float), coeffs, family)
    #enddef
    return radius_from_shape(shape, 1.0)
#enddef



def decode_sachs_family_pair(family: str) -> Tuple[str, str]:
    """Decode a same-family name or an explicit E:<family>|M:<family> pair."""
    family = str(family)
    match = re.fullmatch(r"E:([A-Z]+[0-9]+)\|M:([A-Z]+[0-9]+)", family)
    if match:
        return match.group(1), match.group(2)
    #endif
    return family, family
#enddef


def encode_sachs_family_pair(family_e: str, family_m: str) -> str:
    return f"E:{family_e}|M:{family_m}"
#enddef


def fit_cross_sections_with_sachs_family(
        specs: Sequence[Dict[str, object]],
        family: str,
        bh_cut: float = 0.05,
        add_moradi_bh_systematic: bool = False,
        override_y: Optional[Dict[str, np.ndarray]] = None,
        statistical_only: bool = True,
        return_parameters: bool = False):
    """Fit independently selectable normalized Sachs GE and GM families through BH."""
    family_e, family_m = decode_sachs_family_pair(family)
    ne = int(re.findall(r"\d+", family_e)[0])
    nm = int(re.findall(r"\d+", family_m)[0])
    names_e = [f"e{i}" for i in range(1, ne + 1)]
    names_m = [f"m{i}" for i in range(1, nm + 1)]
    shape_names = names_e + names_m

    p0 = np.zeros(ne + nm, dtype=float)
    p0[0] = sachs_first_coefficient_from_radius(SACHS_INITIAL_RADIUS_FM, family_e)
    p0[ne] = sachs_first_coefficient_from_radius(SACHS_INITIAL_RADIUS_FM, family_m)

    nuisance_names = []
    nuisance_is_free = {}
    nuisance_fracs = {}
    for spec in specs:
        key = str(spec["key"])
        free_norm = bool(spec.get("unconstrained_norm", False))
        frac = 1.0 if free_norm else float(spec.get("norm_frac", 0.0))
        if frac > 0.0:
            nname = "beta_" + re.sub(r"[^A-Za-z0-9]+", "_", key)
            nuisance_names.append(nname)
            nuisance_is_free[nname] = free_norm
            nuisance_fracs[nname] = frac
        #endif
    #endfor
    fit_names = shape_names + nuisance_names
    fit_p0 = np.concatenate([p0, np.zeros(len(nuisance_names), dtype=float)])
    nuisance_index = {
        name: len(shape_names) + i for i, name in enumerate(nuisance_names)
    }

    prepared = []
    max_order = max(ne, nm)
    for spec in specs:
        d = spec["data"]
        y = override_y[str(spec["key"])] if override_y is not None else d["xs"].to_numpy(float)
        err = (
            dataset_statistical_errors(d, str(spec["kind"]))
            if statistical_only
            else dataset_point_errors(d, str(spec["kind"]), bh_cut, add_moradi_bh_systematic)
        )
        q = d["t_abs"].to_numpy(float)
        tau = q / (4.0 * MP2)
        key = str(spec["key"])
        free_norm = bool(spec.get("unconstrained_norm", False))
        frac = 1.0 if free_norm else float(spec.get("norm_frac", 0.0))
        prepared.append({
            "key": key, "norm_frac": frac, "unconstrained_norm": free_norm,
            "nuisance_name": ("beta_" + re.sub(r"[^A-Za-z0-9]+", "_", key) if frac > 0.0 else None),
            "q": q, "q_powers": np.vstack([q**i for i in range(1, max_order + 1)]),
            "tau": tau, "inv_one_plus_tau": 1.0 / (1.0 + tau),
            "y": np.asarray(y, dtype=float), "e": np.asarray(err, dtype=float),
            "A": d["bh_A"].to_numpy(float), "B": d["bh_B"].to_numpy(float),
            "C": d["bh_C"].to_numpy(float),
        })
    #endfor

    def chi2(*values):
        p = np.asarray(values, dtype=float)
        ce = p[:ne]
        cm = p[ne:ne + nm]
        total = 0.0
        for item in prepared:
            q = item["q"]
            ge = sachs_family_value_precomputed(q, ce, family_e, item["q_powers"])
            gm = MU_P * sachs_family_value_precomputed(q, cm, family_m, item["q_powers"])
            tau = item["tau"]
            inv = item["inv_one_plus_tau"]
            f1 = (ge + tau * gm) * inv
            f2 = (gm - ge) * inv
            pred = bh_from_f1f2(item["A"], item["B"], item["C"], f1, f2)
            if item["norm_frac"] > 0.0:
                beta = p[nuisance_index[item["nuisance_name"]]]
                pred *= 1.0 + item["norm_frac"] * beta
            #endif
            pull = (pred - item["y"]) / item["e"]
            total += float(np.dot(pull, pull))
        #endfor
        for nname in nuisance_names:
            if not nuisance_is_free.get(nname, False):
                total += float(p[nuisance_index[nname]]**2)
            #endif
        #endfor
        return total
    #enddef

    m = Minuit(chi2, *fit_p0, name=tuple(fit_names))
    m.errordef = Minuit.LEAST_SQUARES
    # Apply the radius-equivalent first-slope bounds separately because E and M
    # may now use different functional families.
    ce_lo = sachs_first_coefficient_from_radius(SACHS_MIN_RADIUS_FM, family_e)
    ce_hi = sachs_first_coefficient_from_radius(SACHS_MAX_RADIUS_FM, family_e)
    cm_lo = sachs_first_coefficient_from_radius(SACHS_MIN_RADIUS_FM, family_m)
    cm_hi = sachs_first_coefficient_from_radius(SACHS_MAX_RADIUS_FM, family_m)
    m.limits[names_e[0]] = (min(ce_lo, ce_hi), max(ce_lo, ce_hi))
    m.limits[names_m[0]] = (min(cm_lo, cm_hi), max(cm_lo, cm_hi))
    for nname in nuisance_names:
        m.limits[nname] = (
            (-0.50, 0.50) if nuisance_is_free.get(nname, False)
            else (-10.0, 10.0)
        )
    #endfor
    m.migrad()

    full = np.array([float(m.values[n]) for n in fit_names])
    pars = full[:ne + nm]
    rE = sachs_family_radius(pars[:ne], family_e)
    rM = sachs_family_radius(pars[ne:ne + nm], family_m)
    valid = bool(m.valid and np.isfinite(rE) and np.isfinite(rM))
    if return_parameters:
        return rE, rM, valid, pars
    #endif
    return rE, rM, valid
#enddef


# Per-process state for the radius bias/variance replica workers.  These are
# initialized once per truth model so the experimental DataFrames and BH
# coefficients are not serialized with every individual replica task.
_RADIUS_BIAS_WORKER_SPECS = None
_RADIUS_BIAS_WORKER_CENTRAL = None
_RADIUS_BIAS_WORKER_SIGMA = None
_RADIUS_BIAS_WORKER_TRUTH_GE = None
_RADIUS_BIAS_WORKER_TRUTH_GM = None


def _init_radius_bias_worker(
        specs: Sequence[Dict[str, object]],
        central_by_key: Dict[str, np.ndarray],
        sigma_by_key: Dict[str, np.ndarray],
        truth_ge_by_key: Dict[str, np.ndarray],
        truth_gm_by_key: Dict[str, np.ndarray]) -> None:
    global _RADIUS_BIAS_WORKER_SPECS
    global _RADIUS_BIAS_WORKER_CENTRAL
    global _RADIUS_BIAS_WORKER_SIGMA
    global _RADIUS_BIAS_WORKER_TRUTH_GE
    global _RADIUS_BIAS_WORKER_TRUTH_GM

    _RADIUS_BIAS_WORKER_SPECS = specs
    _RADIUS_BIAS_WORKER_CENTRAL = central_by_key
    _RADIUS_BIAS_WORKER_SIGMA = sigma_by_key
    _RADIUS_BIAS_WORKER_TRUTH_GE = truth_ge_by_key
    _RADIUS_BIAS_WORKER_TRUTH_GM = truth_gm_by_key
#enddef


def _radius_bias_replica_worker(
        task: Tuple[str, int]):
    """Generate/fit one replica and return radius plus full-range shape diagnostics."""
    family, seed = task
    rng = np.random.default_rng(int(seed))

    specs = _RADIUS_BIAS_WORKER_SPECS
    central_by_key = _RADIUS_BIAS_WORKER_CENTRAL
    sigma_by_key = _RADIUS_BIAS_WORKER_SIGMA
    truth_ge_by_key = _RADIUS_BIAS_WORKER_TRUTH_GE
    truth_gm_by_key = _RADIUS_BIAS_WORKER_TRUTH_GM
    if (specs is None or central_by_key is None or sigma_by_key is None
            or truth_ge_by_key is None or truth_gm_by_key is None):
        raise RuntimeError("radius-bias worker was not initialized")
    #endif

    replica_y = {}
    for spec in specs:
        key = str(spec["key"])
        central = np.asarray(central_by_key[key], dtype=float)
        free_norm = bool(spec.get("unconstrained_norm", False))
        frac = float(spec.get("norm_frac", 0.0))
        beta_true = (
            rng.normal(0.0, 1.0)
            if (frac > 0.0 and not free_norm) else 0.0
        )
        # The quoted experimental normalization is a cross-section scale.
        # Apply it to sigma_BH, not separately to GE or GM.  Because the BH
        # cross section is homogeneous quadratic in F1/F2, a common form-factor
        # amplitude c0 would correspond to sigma -> c0^2 sigma; using the
        # measured cross-section nuisance directly avoids an artificial
        # GE/GM normalization ambiguity.
        shifted_central = central * (1.0 + frac * beta_true)
        replica_y[key] = rng.normal(shifted_central, sigma_by_key[key])
    #endfor

    try:
        re_val, rm_val, valid, pars = fit_cross_sections_with_sachs_family(
            specs, family=family, bh_cut=0.05,
            add_moradi_bh_systematic=False, override_y=replica_y,
            statistical_only=True, return_parameters=True,
        )
    except Exception:
        return family, np.nan, np.nan, False, np.nan, np.nan, np.nan
    #endtry

    valid = bool(valid and np.isfinite(re_val) and np.isfinite(rm_val))
    ge_frac = []
    gm_frac = []
    if valid:
        family_e, family_m = decode_sachs_family_pair(family)
        ne = int(re.findall(r"\d+", family_e)[0])
        nm = int(re.findall(r"\d+", family_m)[0])
        ce = np.asarray(pars[:ne], dtype=float)
        cm = np.asarray(pars[ne:ne + nm], dtype=float)
        for spec in specs:
            key = str(spec["key"])
            q = spec["data"]["t_abs"].to_numpy(float)
            ge_fit = sachs_family_value(q, ce, family_e)
            gm_fit = sachs_family_value(q, cm, family_m)
            ge_true = np.asarray(truth_ge_by_key[key], dtype=float)
            gm_true = np.asarray(truth_gm_by_key[key], dtype=float)
            good_e = np.isfinite(ge_fit) & np.isfinite(ge_true) & (np.abs(ge_true) > 1.0e-12)
            good_m = np.isfinite(gm_fit) & np.isfinite(gm_true) & (np.abs(gm_true) > 1.0e-12)
            ge_frac.extend(((ge_fit[good_e] - ge_true[good_e]) / ge_true[good_e]).tolist())
            gm_frac.extend(((gm_fit[good_m] - gm_true[good_m]) / gm_true[good_m]).tolist())
        #endfor
    #endif
    ge_rms = float(np.sqrt(np.mean(np.asarray(ge_frac)**2))) if ge_frac else np.nan
    gm_rms = float(np.sqrt(np.mean(np.asarray(gm_frac)**2))) if gm_frac else np.nan
    combined = float(np.sqrt(0.5 * (ge_rms**2 + gm_rms**2))) if np.isfinite(ge_rms) and np.isfinite(gm_rms) else np.nan
    return family, float(re_val), float(rm_val), valid, ge_rms, gm_rms, combined
#enddef



def radius_to_normalized_slope(radius_fm: float) -> float:
    """Return d[G/G(0)]/dQ2 at Q2=0 for a requested RMS radius."""
    return -(float(radius_fm) / HBARC)**2 / 6.0
#enddef


def sachs_family_coefficients_with_radius(
        template_coeffs: np.ndarray,
        family: str,
        radius_fm: float) -> np.ndarray:
    """
    Copy a shape template but force its Q2=0 slope to the requested radius.

    For D1, slope=-2*c1. For Pn, slope=c1. For IPn and CFn, slope=-c1.
    Higher-order coefficients retain the template curvature.
    """
    out = np.asarray(template_coeffs, dtype=float).copy()
    slope = radius_to_normalized_slope(radius_fm)
    if family.startswith("D"):
        out[0] = -0.5 * slope
    elif family.startswith("P"):
        out[0] = slope
    elif family.startswith("IP") or family.startswith("CF"):
        out[0] = -slope
    else:
        raise ValueError(f"unknown Sachs family {family}")
    #endif
    return out
#enddef


def fit_sachs_family_shape_template(
        family: str,
        target_q: np.ndarray,
        target_shape: np.ndarray) -> np.ndarray:
    """Fit a normalized D/P/IP/CF family to a smooth reference shape."""
    q = np.asarray(target_q, dtype=float)
    y = np.asarray(target_shape, dtype=float)
    order = int(re.findall(r"\d+", family)[0])
    names = tuple(f"c{i}" for i in range(1, order + 1))

    p0 = np.zeros(order, dtype=float)
    if family.startswith("D"):
        p0[0] = 1.0 / 0.71
    else:
        p0[0] = -2.0 / 0.71 if family.startswith("P") else 2.0 / 0.71
    #endif

    def objective(*values):
        pred = sachs_family_value(q, np.asarray(values, dtype=float), family)
        scale = np.maximum(np.abs(y), 0.10)
        pull = (pred - y) / scale
        return float(np.dot(pull, pull))
    #enddef

    m = Minuit(objective, *p0, name=names)
    m.errordef = Minuit.LEAST_SQUARES
    m.migrad()
    if not m.valid:
        raise RuntimeError(
            f"Could not construct synthetic truth template for {family}"
        )
    #endif
    return np.asarray([float(m.values[n]) for n in names], dtype=float)
#enddef


def synthetic_truth_is_physical(
        qmax: float,
        family: str,
        coeffs_e: np.ndarray,
        coeffs_m: np.ndarray) -> bool:
    """
    Loose physical sanity screen for synthetic normalized Sachs truth shapes.
    """
    q = np.linspace(0.0, max(float(qmax), 1.0e-4), 500)
    ge = sachs_family_value(q, coeffs_e, family)
    gm = sachs_family_value(q, coeffs_m, family)
    if not np.all(np.isfinite(ge)) or not np.all(np.isfinite(gm)):
        return False
    #endif
    if np.any(ge <= 0.0) or np.any(gm <= 0.0):
        return False
    #endif
    if np.nanmax(ge) > 1.15 or np.nanmax(gm) > 1.15:
        return False
    #endif
    return True
#enddef


def parse_radius_bias_grid(value: str) -> List[float]:
    vals = [float(x.strip()) for x in str(value).split(",") if x.strip()]
    if len(vals) < 2:
        raise ValueError(
            "--radius-bias-radius-grid must contain at least two radii"
        )
    #endif
    if any((not np.isfinite(v)) or v <= 0.0 for v in vals):
        raise ValueError("radius-bias radius grid contains an invalid radius")
    #endif
    return vals
#enddef


def make_synthetic_sachs_truth_scenarios(
        families: Sequence[str],
        radius_values: Sequence[float],
        qmax: float,
        curvature_stress: str = "sectoral") -> List[Dict[str, object]]:
    """
    Build synthetic closure truths with independently controlled slope and curvature.

    Baseline grid:
      For every generating family, retain the Cartesian rE x rM grid used
      previously.  Its higher-order coefficients are obtained by fitting the
      normalized Kelly GE and GM shapes, after which c1 is replaced so the
      requested radius is exact.

    Curvature stress tests:
      At the central radius-grid value only, add four sector-isolated
      Hayward--Griffioen-style curvature challenges:
        * exponential-like GE, Kelly-like GM
        * dipole-like GE,      Kelly-like GM
        * Kelly-like GE,       exponential-like GM
        * Kelly-like GE,       dipole-like GM

      Thus the true radius is held fixed while the finite-Q2 curvature is
      changed independently in the electric or magnetic sector.  Restricting
      these extra shapes to the central radius avoids a full
      radius x curvature Cartesian explosion: for the default three-point
      radius grid this adds only four truths per generating family (13 rather
      than 9), a 44% increase rather than a factor of 5 or more.

    All scenarios from a given generating family retain one common truth_group,
    so adding curvature stress realizations does not give that family extra
    conceptual weight in the equal-truth-group closure objective.
    """
    q_template = np.linspace(0.0, max(float(qmax), 0.05), 400)
    k_ge, k_gm = kelly_sachs(q_template)
    k_ge = np.asarray(k_ge, dtype=float) / float(k_ge[0])
    k_gm = np.asarray(k_gm, dtype=float) / float(k_gm[0])

    # Fit each generating representation once to the baseline Kelly shapes.
    # These fits are cached and reused for all radius points.
    templates = {}
    for family in families:
        templates[family] = {
            "E": fit_sachs_family_shape_template(family, q_template, k_ge),
            "M": fit_sachs_family_shape_template(family, q_template, k_gm),
        }
    #endfor

    scenarios = []
    skipped = 0

    def append_scenario(family, rE, rM, ce, cm, curvature_mode):
        nonlocal skipped
        if not synthetic_truth_is_physical(qmax, family, ce, cm):
            skipped += 1
            return
        #endif

        def truth_fn(q, fam=family, e=ce.copy(), m=cm.copy()):
            q = np.asarray(q, dtype=float)
            ge = sachs_family_value(q, e, fam)
            gm = MU_P * sachs_family_value(q, m, fam)
            return ge, gm
        #enddef

        suffix = "" if curvature_mode == "kelly_template" else f"_{curvature_mode}"
        scenarios.append({
            "truth_model": (
                f"synthetic_{family}_rE{rE:.3f}_rM{rM:.3f}{suffix}"
            ),
            "truth_group": family,
            "truth_family": family,
            "truth_rE_fm": float(rE),
            "truth_rM_fm": float(rM),
            "truth_curvature_mode": curvature_mode,
            "truth_fn": truth_fn,
            "synthetic": True,
        })
    #enddef

    # Existing radius grid with Kelly-derived higher-order curvature.
    for family in families:
        for rE in radius_values:
            for rM in radius_values:
                ce = sachs_family_coefficients_with_radius(
                    templates[family]["E"], family, rE
                )
                cm = sachs_family_coefficients_with_radius(
                    templates[family]["M"], family, rM
                )
                append_scenario(
                    family, rE, rM, ce, cm, "kelly_template"
                )
            #endfor
        #endfor
    #endfor

    # Add independent curvature challenges without forming a large Cartesian
    # product.  The median grid radius is used so slope is fixed at a central,
    # physically relevant value while curvature alone is changed.
    curvature_added = 0
    if str(curvature_stress).lower() == "sectoral":
        radius_sorted = sorted(float(x) for x in radius_values)
        r0 = float(radius_sorted[len(radius_sorted) // 2])
        slope0 = radius_to_normalized_slope(r0)

        # Complete normalized reference shapes with the same exact slope/radius.
        # These are the two generators used in the Hayward--Griffioen spirit:
        # exponential and dipole differ in Q4,Q6,... while sharing G(0) and G'(0).
        exp_shape = np.exp(slope0 * q_template)
        dip_scale = -(2.0 / slope0) if slope0 < 0.0 else np.inf
        dip_shape = (1.0 + q_template / dip_scale) ** -2

        for family in families:
            try:
                exp_template = fit_sachs_family_shape_template(
                    family, q_template, exp_shape
                )
                dip_template = fit_sachs_family_shape_template(
                    family, q_template, dip_shape
                )
            except Exception as exc:
                print(
                    f"[radius-bias] curvature template warning: {family}: "
                    f"{type(exc).__name__}: {exc}; skipping sectoral curvature "
                    "stress truths for this generating family"
                )
                continue
            #endtry

            # Re-impose the radius analytically after the template fit.  This
            # guarantees that any closure difference is driven by higher-order
            # shape rather than a small slope mismatch from the template fit.
            e_k = sachs_family_coefficients_with_radius(
                templates[family]["E"], family, r0
            )
            m_k = sachs_family_coefficients_with_radius(
                templates[family]["M"], family, r0
            )
            e_exp = sachs_family_coefficients_with_radius(
                exp_template, family, r0
            )
            e_dip = sachs_family_coefficients_with_radius(
                dip_template, family, r0
            )
            m_exp = sachs_family_coefficients_with_radius(
                exp_template, family, r0
            )
            m_dip = sachs_family_coefficients_with_radius(
                dip_template, family, r0
            )

            before = len(scenarios)
            append_scenario(
                family, r0, r0, e_exp, m_k, "E_exponential_M_kelly"
            )
            append_scenario(
                family, r0, r0, e_dip, m_k, "E_dipole_M_kelly"
            )
            append_scenario(
                family, r0, r0, e_k, m_exp, "E_kelly_M_exponential"
            )
            append_scenario(
                family, r0, r0, e_k, m_dip, "E_kelly_M_dipole"
            )
            curvature_added += len(scenarios) - before
        #endfor
    elif str(curvature_stress).lower() != "none":
        raise ValueError(
            f"unknown radius-bias curvature stress mode: {curvature_stress}"
        )
    #endif

    print(
        f"[radius-bias] synthetic truth ensemble: {len(scenarios)} accepted "
        f"scenario(s), {skipped} rejected by physical-shape screen; "
        f"radii={list(radius_values)}; curvature_stress={curvature_stress}; "
        f"sectoral_curvature_truths_added={curvature_added}"
    )
    return scenarios
#enddef



def robust_upper_limit(
        values,
        quantile: float = 0.98,
        pad_fraction: float = 0.12,
        minimum: float | None = None) -> float:
    """
    Robust upper plotting limit that ignores a very small number of extreme
    failures while preserving the scale occupied by the bulk of the study.

    This changes only visualization; no values are clipped in CSV output or
    in any ranking/objective calculation.
    """
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if len(arr) == 0:
        return 1.0 if minimum is None else float(minimum)
    #endif

    q = float(np.nanquantile(arr, quantile))
    med = float(np.nanmedian(arr))
    upper = q + pad_fraction * max(abs(q - med), abs(q), 1.0e-6)
    if minimum is not None:
        upper = max(upper, float(minimum))
    #endif
    return upper
#enddef


def robust_symmetric_limit(
        values,
        quantile: float = 0.98,
        pad_fraction: float = 0.12,
        minimum: float = 0.02) -> float:
    """Robust symmetric +/- plotting range around zero."""
    arr = np.asarray(values, dtype=float)
    arr = np.abs(arr[np.isfinite(arr)])
    if len(arr) == 0:
        return float(minimum)
    #endif
    q = float(np.nanquantile(arr, quantile))
    return max(float(minimum), q * (1.0 + pad_fraction))
#enddef


def robust_radius_limits(
        values,
        truth_values=None,
        lower_quantile: float = 0.01,
        upper_quantile: float = 0.99,
        pad_fraction: float = 0.08) -> tuple[float, float]:
    """
    Robust y-range for extracted-radius plots.

    Extreme failed-fit radii are ignored for display only. Truth radii are
    explicitly included so the physically relevant reference range is never
    cropped.
    """
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if truth_values is not None:
        truth = np.asarray(truth_values, dtype=float)
        truth = truth[np.isfinite(truth)]
        arr = np.concatenate([arr, truth])
    #endif
    if len(arr) == 0:
        return 0.5, 1.1
    #endif

    lo = float(np.nanquantile(arr, lower_quantile))
    hi = float(np.nanquantile(arr, upper_quantile))
    span = max(hi - lo, 0.05)
    return lo - pad_fraction * span, hi + pad_fraction * span
#enddef



EMPIRICAL_VALIDATION_TRUTH_GROUP = "empirical_elastic"


def production_closure_truth_rows(table: pd.DataFrame) -> pd.DataFrame:
    """Return only controlled synthetic truths used for production closure."""
    if "truth_group" not in table.columns:
        return table.copy()
    #endif
    return table.loc[
        table["truth_group"].astype(str)
        != EMPIRICAL_VALIDATION_TRUTH_GROUP
    ].copy()
#enddef


def empirical_validation_truth_rows(table: pd.DataFrame) -> pd.DataFrame:
    """Return only Kelly/AMT/Bernauer validation rows."""
    if "truth_group" not in table.columns:
        return table.iloc[0:0].copy()
    #endif
    return table.loc[
        table["truth_group"].astype(str)
        == EMPIRICAL_VALIDATION_TRUTH_GROUP
    ].copy()
#enddef


def save_empirical_closure_validation_summary(
        table: pd.DataFrame,
        outdir: Path) -> pd.DataFrame:
    """Save Kelly/AMT/Bernauer recovery as validation, never as production votes."""
    empirical = empirical_validation_truth_rows(table)
    outdir.mkdir(parents=True, exist_ok=True)

    if len(empirical) == 0:
        empty = pd.DataFrame()
        empty.to_csv(outdir / "empirical_truth_validation_summary.csv", index=False)
        return empty
    #endif

    detail_cols = [
        c for c in [
            "truth_model", "truth_family", "family", "quantity",
            "truth_radius_fm", "fit_mean_radius_fm", "bias_fm",
            "stat_RMS_fm", "sqrt_stat2_plus_bias2_fm",
            "valid_replicas", "requested_replicas",
            "attempted_replicas", "attempt_efficiency", "target_reached",
        ] if c in empirical.columns
    ]
    empirical[detail_cols].copy().to_csv(
        outdir / "empirical_truth_validation_detail.csv", index=False
    )

    rows = []
    for family, fam in empirical.groupby("family", sort=False):
        row = {"family": str(family)}
        for quantity in ["rE", "rM"]:
            qpart = fam.loc[fam["quantity"].astype(str) == quantity]
            b = qpart["bias_fm"].to_numpy(float)
            b = b[np.isfinite(b)]
            rmse = qpart["sqrt_stat2_plus_bias2_fm"].to_numpy(float)
            rmse = rmse[np.isfinite(rmse)]
            row[f"{quantity}_empirical_validation_RMS_bias_fm"] = (
                float(np.sqrt(np.mean(b**2))) if len(b) else np.nan
            )
            row[f"{quantity}_empirical_validation_RMS_RMSE_fm"] = (
                float(np.sqrt(np.mean(rmse**2))) if len(rmse) else np.nan
            )
            row[f"{quantity}_N_empirical_truths"] = int(len(qpart))
        #endfor
        rows.append(row)
    #endfor

    summary = pd.DataFrame(rows)
    if len(summary):
        summary["combined_empirical_validation_RMSE_fm"] = np.sqrt(
            0.5 * (
                summary["rE_empirical_validation_RMS_RMSE_fm"]**2
                + summary["rM_empirical_validation_RMS_RMSE_fm"]**2
            )
        )
        summary = summary.sort_values(
            "combined_empirical_validation_RMSE_fm"
        ).reset_index(drop=True)
    #endif
    summary.to_csv(
        outdir / "empirical_truth_validation_summary.csv", index=False
    )

    print(
        "[radius-bias] Kelly/AMT/Bernauer retained as validation only; "
        "zero production weight"
    )
    return summary
#enddef


def save_extended_radius_bias_matrices(
        table: pd.DataFrame,
        families: Sequence[str],
        outdir: Path,
        minimum_global_valid_fraction: float,
        minimum_scenario_valid_fraction: float) -> None:
    """
    Make production closure matrices/rankings from controlled synthetic truths.

    Kelly/AMT/Bernauer are validation-only and cannot favor or veto a family.
    """
    table = production_closure_truth_rows(table)
    if len(table) == 0:
        raise RuntimeError("Production closure has zero synthetic truth rows.")
    #endif
    groups = list(dict.fromkeys(table["truth_group"].astype(str).tolist()))
    matrix_rows = []

    for quantity in ["rE", "rM"]:
        objective = np.full((len(groups), len(families)), np.nan)
        abs_bias = np.full_like(objective, np.nan)
        valid_frac = np.full_like(objective, np.nan)

        for ig, group in enumerate(groups):
            for jf, fit_family in enumerate(families):
                part = table.loc[
                    (table["truth_group"] == group)
                    & (table["family"] == fit_family)
                    & (table["quantity"] == quantity)
                ]
                obj = part["sqrt_stat2_plus_bias2_fm"].to_numpy(float)
                bias = part["bias_fm"].to_numpy(float)
                nvalid = part["valid_replicas"].to_numpy(float)
                nreq = part["requested_replicas"].to_numpy(float)

                finite_obj = obj[np.isfinite(obj)]
                finite_bias = bias[np.isfinite(bias)]
                if len(finite_obj):
                    objective[ig, jf] = float(
                        np.sqrt(np.mean(finite_obj**2))
                    )
                #endif
                if len(finite_bias):
                    abs_bias[ig, jf] = float(
                        np.sqrt(np.mean(finite_bias**2))
                    )
                #endif
                if len(nvalid):
                    valid_frac[ig, jf] = float(
                        np.sum(nvalid) / np.sum(nreq)
                    )
                #endif

                matrix_rows.append({
                    "quantity": quantity,
                    "truth_group": group,
                    "fit_family": fit_family,
                    "RMS_objective_fm": objective[ig, jf],
                    "RMS_abs_bias_fm": abs_bias[ig, jf],
                    "aggregate_valid_fraction": valid_frac[ig, jf],
                })
            #endfor
        #endfor

        for metric_name, values, filename, label in [
            (
                "objective", objective,
                f"07_{quantity}_truth_family_x_fit_family_objective.png",
                r"RMS $\sqrt{\sigma_{\rm stat}^2+b^2}$ (fm)",
            ),
            (
                "bias", abs_bias,
                f"08_{quantity}_truth_family_x_fit_family_bias.png",
                "RMS |bias| (fm)",
            ),
            (
                "validity", valid_frac,
                f"09_{quantity}_truth_family_x_fit_family_validity.png",
                "valid replica fraction",
            ),
        ]:
            fig, ax = plt.subplots(
                figsize=(1.0 + 0.85 * len(families), 2.0 + 0.55 * len(groups))
            )
            finite_values = values[np.isfinite(values)]
            if metric_name == "validity":
                vmin, vmax = 0.0, 1.0
                extend = "neither"
            else:
                vmin = 0.0
                vmax = robust_upper_limit(
                    finite_values,
                    quantile=0.98,
                    pad_fraction=0.08,
                )
                extend = "max"
            #endif
            image = ax.imshow(
                values,
                aspect="auto",
                vmin=vmin,
                vmax=vmax,
            )
            ax.set_xticks(np.arange(len(families)))
            ax.set_xticklabels(families, rotation=45, ha="right")
            ax.set_yticks(np.arange(len(groups)))
            ax.set_yticklabels(groups)
            ax.set_xlabel("fit family")
            ax.set_ylabel("generating truth family/model")
            ax.set_title(
                f"{quantity}: generating family vs fitting family ({metric_name})"
            )
            cbar = fig.colorbar(image, ax=ax, pad=0.02, extend=extend)
            cbar.set_label(label)
            fig.tight_layout()
            fig.savefig(outdir / filename, dpi=180)
            plt.close(fig)
        #endfor
    #endfor

    pd.DataFrame(matrix_rows).to_csv(
        outdir / "radius_bias_truth_family_fit_family_matrix.csv",
        index=False,
    )

    rank_rows = []
    for family in families:
        fam_all = table.loc[table["family"] == family]
        scenario_valid = (
            fam_all["valid_replicas"].to_numpy(float)
            / fam_all["requested_replicas"].to_numpy(float)
        )
        min_valid = (
            float(np.nanmin(scenario_valid))
            if len(scenario_valid)
            else np.nan
        )
        global_valid = float(
            np.nansum(fam_all["valid_replicas"])
            / np.nansum(fam_all["requested_replicas"])
        )

        row = {
            "family": family,
            "minimum_scenario_valid_fraction": min_valid,
            "global_valid_fraction": global_valid,
            "minimum_global_valid_fraction_required": float(
                minimum_global_valid_fraction
            ),
            "minimum_scenario_valid_fraction_required": float(
                minimum_scenario_valid_fraction
            ),
            "eligible": bool(
                np.isfinite(global_valid)
                and np.isfinite(min_valid)
                and global_valid >= float(minimum_global_valid_fraction)
                and min_valid >= float(minimum_scenario_valid_fraction)
            ),
            "convergence_requirement_applied": True,
        }

        for quantity in ["rE", "rM"]:
            part = fam_all.loc[fam_all["quantity"] == quantity]
            vals = part["sqrt_stat2_plus_bias2_fm"].to_numpy(float)
            vals = vals[np.isfinite(vals)]
            biases = part["bias_fm"].to_numpy(float)
            biases = biases[np.isfinite(biases)]
            row[f"{quantity}_RMS_objective_all_truths_fm"] = (
                equal_truth_group_rms(part, "sqrt_stat2_plus_bias2_fm")
            )
            row[f"{quantity}_max_objective_all_truths_fm"] = (
                float(np.max(vals)) if len(vals) else np.nan
            )
            row[f"{quantity}_RMS_bias_all_truths_fm"] = (
                equal_truth_group_rms(part, "bias_fm")
            )
            row[f"{quantity}_scenario_weighted_RMS_objective_fm"] = (
                float(np.sqrt(np.mean(vals**2))) if len(vals) else np.nan
            )
            row[f"{quantity}_scenario_weighted_RMS_bias_fm"] = (
                float(np.sqrt(np.mean(biases**2))) if len(biases) else np.nan
            )
        #endfor

        e = row["rE_RMS_objective_all_truths_fm"]
        m = row["rM_RMS_objective_all_truths_fm"]
        row["combined_RMS_objective_fm"] = (
            float(math.sqrt(0.5 * (e**2 + m**2)))
            if (
                bool(row["eligible"])
                and np.isfinite(e)
                and np.isfinite(m)
            )
            else np.nan
        )
        rank_rows.append(row)
    #endfor

    ranking = pd.DataFrame(rank_rows).sort_values(
        ["combined_RMS_objective_fm"],
        ascending=[True],
    )
    ranking.to_csv(
        outdir / "radius_bias_extended_eligibility_ranking.csv",
        index=False,
    )

    fig, ax = plt.subplots(figsize=(10.5, 5.3))
    order = ranking["family"].tolist()
    xmap = {f: i for i, f in enumerate(order)}
    finite = ranking.loc[np.isfinite(ranking["combined_RMS_objective_fm"])]
    ax.scatter(
        [xmap[f] for f in finite["family"]],
        finite["combined_RMS_objective_fm"],
        marker="o",
        label="eligible families; ranked by closure objective",
    )

    ax.set_xticks(np.arange(len(order)))
    ax.set_xticklabels(order, rotation=45)
    ax.set_ylabel("combined RMS bias-variance objective (fm)")
    ax.set_title("Extended closure-study ranking after convergence eligibility")
    yvals = ranking["combined_RMS_objective_fm"].to_numpy(float)
    ax.set_ylim(
        0.0,
        robust_upper_limit(
            yvals,
            quantile=0.90,
            pad_fraction=0.12,
            minimum=0.05,
        ),
    )
    ax.grid(alpha=0.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "10_extended_eligibility_ranking.png", dpi=180)
    plt.close(fig)
#enddef



def _radius_bias_family_batch_worker(task: Tuple[str, Sequence[int]]):
    """Run all replica seeds for one family in one multiprocessing dispatch."""
    family, seeds = task
    return family, [_radius_bias_replica_worker((family, int(seed))) for seed in seeds]
#enddef




def equal_truth_group_rms(
        table: pd.DataFrame,
        value_column: str) -> float:
    """
    Hierarchical RMS with equal weight for each truth group.

    First compute mean(value^2) within each truth_group, then average those
    group-level mean-squares with equal group weight. Production callers
    remove empirical validation rows first, so equal weighting applies only
    among controlled synthetic truth groups.
    """
    group_ms = []
    for _, group in table.groupby("truth_group", sort=False):
        vals = group[value_column].to_numpy(float)
        vals = vals[np.isfinite(vals)]
        if len(vals):
            group_ms.append(float(np.mean(vals**2)))
        #endif
    #endfor
    if not group_ms:
        return np.nan
    #endif
    return float(np.sqrt(np.mean(group_ms)))
#enddef




def equal_truth_group_mean(
        table: pd.DataFrame,
        value_column: str) -> float:
    """
    Hierarchical SIGNED mean with equal total weight per conceptual truth group.

    This is complementary to equal_truth_group_rms().  It preserves the
    direction of a systematic bias.  A positive value means the fitted radius
    tends to lie above the generated truth and therefore a first-order
    bias-corrected real-data radius would subtract this value.
    """
    group_means = []
    for _, group in table.groupby("truth_group", sort=False):
        vals = group[value_column].to_numpy(float)
        vals = vals[np.isfinite(vals)]
        if len(vals):
            group_means.append(float(np.mean(vals)))
        #endif
    #endfor
    if not group_means:
        return np.nan
    #endif
    return float(np.mean(group_means))
#enddef


def equal_truth_group_positive_fraction(
        table: pd.DataFrame,
        value_column: str) -> float:
    """Equal-truth-group fraction of scenarios whose signed value is positive."""
    group_fracs = []
    for _, group in table.groupby("truth_group", sort=False):
        vals = group[value_column].to_numpy(float)
        vals = vals[np.isfinite(vals)]
        if len(vals):
            group_fracs.append(float(np.mean(vals > 0.0)))
        #endif
    #endfor
    if not group_fracs:
        return np.nan
    #endif
    return float(np.mean(group_fracs))
#enddef


def save_mixed_family_closure_ranking(
        table: pd.DataFrame,
        candidate_pairs: Sequence[str],
        outdir: Path,
        minimum_global_valid_fraction: float,
        minimum_scenario_valid_fraction: float) -> None:
    """Rank E/M families using controlled synthetic truths only."""
    table = production_closure_truth_rows(table)
    if len(table) == 0:
        raise RuntimeError("Mixed-family ranking has zero synthetic truth rows.")
    #endif
    rows = []
    for pair in candidate_pairs:
        part = table.loc[table["family"] == pair]
        attempted = part["attempted_replicas"].to_numpy(float)
        valid = part["valid_replicas"].to_numpy(float)
        valid_frac = np.divide(
            valid, attempted, out=np.zeros_like(valid), where=attempted > 0.0
        )
        global_valid = float(
            part["valid_replicas"].sum() / part["attempted_replicas"].sum()
        ) if float(part["attempted_replicas"].sum()) > 0.0 else np.nan
        min_valid = float(np.nanmin(valid_frac)) if len(valid_frac) else np.nan
        all_targets_reached = bool(part["target_reached"].astype(bool).all())
        eligible = bool(
            all_targets_reached
            and np.isfinite(global_valid) and np.isfinite(min_valid)
            and global_valid >= minimum_global_valid_fraction
            and min_valid >= minimum_scenario_valid_fraction
        )
        row = {
            "family": pair,
            "family_E": decode_sachs_family_pair(pair)[0],
            "family_M": decode_sachs_family_pair(pair)[1],
            "global_valid_fraction": global_valid,
            "minimum_scenario_valid_fraction": min_valid,
            "all_targets_reached": all_targets_reached,
            "eligible": eligible,
        }
        for quantity in ["rE", "rM"]:
            qpart = part.loc[part["quantity"] == quantity]

            # Production definition: each conceptual truth_group gets one
            # equal vote, regardless of whether it contains 1 or 9 scenarios.
            row[f"{quantity}_RMS_bias_fm"] = equal_truth_group_rms(
                qpart, "bias_fm"
            )
            row[f"{quantity}_RMS_replica_std_fm"] = equal_truth_group_rms(
                qpart, "stat_RMS_fm"
            )
            row[f"{quantity}_RMS_RMSE_fm"] = equal_truth_group_rms(
                qpart, "sqrt_stat2_plus_bias2_fm"
            )
            row[f"{quantity}_mean_signed_bias_fm"] = equal_truth_group_mean(
                qpart, "bias_fm"
            )
            row[f"{quantity}_positive_bias_fraction"] = (
                equal_truth_group_positive_fraction(qpart, "bias_fm")
            )
            signed = float(row[f"{quantity}_mean_signed_bias_fm"])
            rms_bias = float(row[f"{quantity}_RMS_bias_fm"])
            row[f"{quantity}_residual_bias_RMS_after_mean_correction_fm"] = (
                float(np.sqrt(max(0.0, rms_bias**2 - signed**2)))
                if np.isfinite(signed) and np.isfinite(rms_bias) else np.nan
            )

            # Preserve the previous scenario-weighted quantities for auditing.
            bias = qpart["bias_fm"].to_numpy(float)
            stat = qpart["stat_RMS_fm"].to_numpy(float)
            obj = qpart["sqrt_stat2_plus_bias2_fm"].to_numpy(float)
            row[f"{quantity}_scenario_weighted_RMS_bias_fm"] = (
                float(np.sqrt(np.nanmean(bias**2))) if len(bias) else np.nan
            )
            row[f"{quantity}_scenario_weighted_RMS_replica_std_fm"] = (
                float(np.sqrt(np.nanmean(stat**2))) if len(stat) else np.nan
            )
            row[f"{quantity}_scenario_weighted_RMS_RMSE_fm"] = (
                float(np.sqrt(np.nanmean(obj**2))) if len(obj) else np.nan
            )
        #endfor
        row["combined_RMS_objective_fm"] = (
            float(np.sqrt(0.5 * (
                row["rE_RMS_RMSE_fm"]**2 + row["rM_RMS_RMSE_fm"]**2
            ))) if eligible else np.nan
        )
        rows.append(row)
    #endfor
    ranking = pd.DataFrame(rows).sort_values(
        ["eligible", "combined_RMS_objective_fm"],
        ascending=[False, True], na_position="last",
    )
    ranking.to_csv(outdir / "radius_bias_mixed_family_ranking.csv", index=False)

    eligible = ranking.loc[ranking["eligible"]].head(25).copy()
    if len(eligible):
        x = np.arange(len(eligible))
        fig, ax = plt.subplots(figsize=(max(10.0, 0.48 * len(eligible)), 5.5))
        ax.plot(x, eligible["rE_RMS_bias_fm"], marker="o", label=r"$r_E$ RMS bias")
        ax.plot(x, eligible["rE_RMS_replica_std_fm"], marker="s", label=r"$r_E$ replica RMS")
        ax.plot(x, eligible["rM_RMS_bias_fm"], marker="o", linestyle="--", label=r"$r_M$ RMS bias")
        ax.plot(x, eligible["rM_RMS_replica_std_fm"], marker="s", linestyle="--", label=r"$r_M$ replica RMS")
        ax.set_xticks(x)
        ax.set_xticklabels(eligible["family"], rotation=60, ha="right")
        ax.set_ylabel("fm")
        ax.set_title("Top mixed E/M families: closure bias and replica variance")
        ax.grid(alpha=0.2)
        ax.legend()
        fig.tight_layout()
        fig.savefig(outdir / "11_mixed_family_bias_vs_variance.png", dpi=180)
        plt.close(fig)

        fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.3))
        for ax, quantity, ylabel in [
            (axes[0], "rE", r"Signed $r_E$ bias [fm]"),
            (axes[1], "rM", r"Signed $r_M$ bias [fm]"),
        ]:
            vals = eligible[f"{quantity}_mean_signed_bias_fm"].to_numpy(float)
            rms = eligible[f"{quantity}_RMS_bias_fm"].to_numpy(float)
            ax.axhline(0.0, linewidth=0.9, linestyle="--")
            ax.errorbar(
                x, vals, yerr=rms, marker="o", linestyle="none", capsize=3,
                label="equal-truth-group signed mean ± RMS magnitude",
            )
            ax.set_xticks(x)
            ax.set_xticklabels(eligible["family"], rotation=60, ha="right")
            ax.set_ylabel(ylabel)
            ax.grid(axis="y", alpha=0.2)
        #endfor
        fig.suptitle(
            "Closure bias direction for the leading mixed Sachs families",
            y=0.995,
        )
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        fig.savefig(outdir / "12_mixed_family_signed_bias.png", dpi=220)
        plt.close(fig)
    #endif
#enddef




def save_empirical_truth_similarity_diagnostics(
        qmax: float,
        outdir: Path) -> pd.DataFrame:
    """
    Quantify how similar the three empirical elastic truth parameterizations are.

    This makes the empirical validation ensemble auditable. The three
    curves are external known-answer checks only and receive zero production
    closure/model-selection weight.
    """
    qhi = min(1.0, max(0.05, float(qmax)))
    q = np.linspace(0.0, qhi, 800)
    models = {
        "Kelly": kelly_sachs,
        "AMT2007": amt2007_sachs,
        "Bernauer": bernauer_polyxdipole_sachs,
    }
    values = {}
    for name, fn in models.items():
        ge, gm = fn(q)
        values[name] = (
            np.asarray(ge, dtype=float),
            np.asarray(gm, dtype=float) / MU_P,
        )
    #endfor

    rows = []
    names = list(models)
    for i, left in enumerate(names):
        for right in names[i + 1:]:
            for quantity, iq in [("GE", 0), ("GM_over_mu", 1)]:
                a = values[left][iq]
                b = values[right][iq]
                scale = np.maximum(1.0e-12, 0.5 * (np.abs(a) + np.abs(b)))
                frac = (a - b) / scale
                rows.append({
                    "model_a": left,
                    "model_b": right,
                    "quantity": quantity,
                    "Q2_max_GeV2": qhi,
                    "fractional_RMS_difference": float(
                        np.sqrt(np.mean(frac**2))
                    ),
                    "maximum_absolute_fractional_difference": float(
                        np.max(np.abs(frac))
                    ),
                })
            #endfor
        #endfor
    #endfor
    table = pd.DataFrame(rows)
    table.to_csv(outdir / "empirical_truth_pairwise_similarity.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.9))
    kge, kgm = values["Kelly"]
    for name in names:
        ge, gm = values[name]
        axes[0].plot(q, ge / kge, label=name)
        axes[1].plot(q, gm / kgm, label=name)
    #endfor
    for ax, ylabel in zip(
            axes,
            [r"$G_E/G_E^{\rm Kelly}$",
             r"$(G_M/\mu_p)/(G_M/\mu_p)^{\rm Kelly}$"]):
        ax.axhline(1.0, linewidth=0.8, linestyle="--")
        ax.set_xlim(0.0, 1.0)
        ax.set_xlabel(r"$Q^2$ [GeV$^2$]")
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.2)
    #endfor
    axes[0].legend(frameon=False)
    fig.suptitle(
        "Empirical elastic truth similarity (one conceptual truth group)",
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(outdir / "00_empirical_truth_similarity.png", dpi=240)
    plt.close(fig)
    return table
#enddef


def pilot_prune_radius_bias_families(
        specs: Sequence[Dict[str, object]],
        scenarios: Sequence[Dict[str, object]],
        families: Sequence[str],
        args,
        outdir: Path) -> List[str]:
    """Cheap first-stage closure screen before the full high-statistics study.

    The pilot uses a fixed number of ATTEMPTED replicas (no retry-to-success
    loop) on one representative empirical truth (Kelly) plus a few
    central-radius, high-curvature synthetic truths.  It is deliberately conservative: numerically unreliable
    pairs are rejected, while all pairs close to the leader and a protected
    minimum top cohort survive.  This prevents a 10-replica fluctuation from
    prematurely eliminating near-degenerate forms such as IP2 and CF2.
    """
    npilot = max(0, int(getattr(args, "radius_bias_pilot_replicas", 10)))
    if npilot <= 0 or int(args.radius_bias_replicas) <= npilot or len(families) <= 1:
        return list(families)
    #endif

    # Pilot pruning is based only on controlled synthetic truths.
    # Empirical Kelly/AMT/Bernauer cases are validation-only.
    synthetic = [
        sc for sc in scenarios
        if bool(sc.get("synthetic", False))
        and str(sc.get("truth_family", "")) in {"P4", "IP4", "CF4"}
        and abs(float(sc.get("truth_rE_fm", np.nan)) - 0.85) < 1.0e-9
        and abs(float(sc.get("truth_rM_fm", np.nan)) - 0.85) < 1.0e-9
    ]
    pilot_scenarios = synthetic
    if not pilot_scenarios:
        pilot_scenarios = [
            sc for sc in scenarios if bool(sc.get("synthetic", False))
        ][:min(3, len(scenarios))]
    #endif

    nworkers = getattr(args, "radius_bias_workers", None)
    nworkers = int(nworkers if nworkers is not None else args.workers)
    nworkers = max(1, nworkers)
    records = []
    print(
        f"[radius-bias pilot] {len(families)} family pairs x "
        f"{len(pilot_scenarios)} pilot truths x {npilot} attempted replicas; "
        f"workers={nworkers}"
    )

    for truth_index, scenario in enumerate(pilot_scenarios):
        truth_fn = scenario["truth_fn"]
        if bool(scenario.get("synthetic", False)):
            truth_re = float(scenario["truth_rE_fm"])
            truth_rm = float(scenario["truth_rM_fm"])
        else:
            truth_re = radius_from_shape(
                lambda qq: truth_fn(np.asarray(qq, dtype=float))[0], 1.0
            )
            truth_rm = radius_from_shape(
                lambda qq: truth_fn(np.asarray(qq, dtype=float))[1], MU_P
            )
        #endif

        central_by_key = {}
        sigma_by_key = {}
        truth_ge_by_key = {}
        truth_gm_by_key = {}
        for spec in specs:
            d = spec["data"]
            q = d["t_abs"].to_numpy(float)
            ge, gm = truth_fn(q)
            truth_ge_by_key[str(spec["key"])] = np.asarray(ge, dtype=float)
            truth_gm_by_key[str(spec["key"])] = np.asarray(gm, dtype=float) / MU_P
            tau = q / (4.0 * MP2)
            f1 = (ge + tau * gm) / (1.0 + tau)
            f2 = (gm - ge) / (1.0 + tau)
            central_by_key[str(spec["key"])] = bh_from_f1f2(
                d["bh_A"].to_numpy(float), d["bh_B"].to_numpy(float),
                d["bh_C"].to_numpy(float), f1, f2,
            )
            sigma_by_key[str(spec["key"])] = dataset_statistical_errors(
                d, str(spec["kind"])
            )
        #endfor

        seed_root = np.random.SeedSequence(
            [int(args.radius_bias_seed), 900000 + int(truth_index)]
        )
        family_roots = seed_root.spawn(len(families))
        tasks = []
        for ifam, family in enumerate(families):
            seeds = [
                int(ss.generate_state(1, dtype=np.uint64)[0])
                for ss in family_roots[ifam].spawn(npilot)
            ]
            tasks.append((family, seeds))
        #endfor

        with ProcessPoolExecutor(
                max_workers=nworkers,
                initializer=_init_radius_bias_worker,
                initargs=(specs, central_by_key, sigma_by_key,
                          truth_ge_by_key, truth_gm_by_key)) as pool:
            for family, family_results in pool.map(
                    _radius_bias_family_batch_worker, tasks, chunksize=1):
                good = [r for r in family_results if bool(r[3])]
                re_vals = np.asarray([r[1] for r in good], dtype=float)
                rm_vals = np.asarray([r[2] for r in good], dtype=float)
                def one_obj(vals, truth):
                    if len(vals) == 0:
                        return np.nan, np.nan, np.nan
                    #endif
                    mean = float(np.mean(vals))
                    rms = float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0
                    bias = mean - float(truth)
                    return bias, rms, float(np.hypot(bias, rms))
                #enddef
                be, se, oe = one_obj(re_vals, truth_re)
                bm, sm, om = one_obj(rm_vals, truth_rm)
                combo = float(np.sqrt(0.5 * (oe**2 + om**2))) if np.isfinite(oe) and np.isfinite(om) else np.nan
                records.append({
                    "truth_model": str(scenario["truth_model"]),
                    "family": family,
                    "attempted": npilot,
                    "valid": len(good),
                    "valid_fraction": len(good) / npilot,
                    "rE_bias_fm": be, "rE_RMS_fm": se,
                    "rM_bias_fm": bm, "rM_RMS_fm": sm,
                    "combined_objective_fm": combo,
                })
            #endfor
        #endwith
    #endfor

    detail = pd.DataFrame(records)
    detail.to_csv(outdir / "radius_bias_pilot_detail.csv", index=False)
    summary_rows = []
    for family in families:
        part = detail.loc[detail["family"] == family]
        vf = float(part["valid"].sum() / part["attempted"].sum()) if len(part) else 0.0
        objs = part["combined_objective_fm"].to_numpy(float)
        finite = objs[np.isfinite(objs)]
        score = float(np.sqrt(np.mean(finite**2))) if len(finite) else np.inf
        worst = float(np.max(finite)) if len(finite) else np.inf
        summary_rows.append({
            "family": family, "pilot_valid_fraction": vf,
            "pilot_RMS_objective_fm": score,
            "pilot_max_objective_fm": worst,
        })
    #endfor
    summary = pd.DataFrame(summary_rows).sort_values(
        ["pilot_RMS_objective_fm", "pilot_max_objective_fm"]
    ).reset_index(drop=True)

    pilot_min_valid = float(
        getattr(args, "radius_bias_pilot_min_valid_fraction", 0.50)
    )
    reliable = summary.loc[
        np.isfinite(summary["pilot_RMS_objective_fm"])
        & (summary["pilot_valid_fraction"] >= pilot_min_valid)
    ].copy()
    if len(reliable) == 0:
        reliable = summary.loc[np.isfinite(summary["pilot_RMS_objective_fm"])].copy()
    #endif
    if len(reliable) == 0:
        print("[radius-bias pilot] no finite pilot scores; disabling pruning")
        summary["survives_pilot"] = True
        summary.to_csv(outdir / "radius_bias_pilot_ranking.csv", index=False)
        return list(families)
    #endif

    best = float(reliable.iloc[0]["pilot_RMS_objective_fm"])
    reltol = max(0.0, float(getattr(args, "radius_bias_pilot_relative_tolerance", 0.35)))
    abstol = max(0.0, float(getattr(args, "radius_bias_pilot_absolute_tolerance_fm", 0.020)))
    cutoff = max(best * (1.0 + reltol), best + abstol)
    min_keep = max(1, int(getattr(args, "radius_bias_pilot_min_keep", 10)))
    max_keep = max(min_keep, int(getattr(args, "radius_bias_pilot_max_keep", 16)))

    survivor_set = set(reliable.loc[
        reliable["pilot_RMS_objective_fm"] <= cutoff, "family"
    ].astype(str))
    survivor_set.update(reliable.head(min(min_keep, len(reliable)))["family"].astype(str))
    ordered_survivors = [f for f in reliable["family"].astype(str) if f in survivor_set]
    if len(ordered_survivors) > max_keep:
        ordered_survivors = ordered_survivors[:max_keep]
    #endif
    # Explicitly protect the common low-order IP2/CF2 neighborhood whenever
    # those pairs are numerically reliable in the pilot.
    protected = {
        encode_sachs_family_pair(e, m)
        for e in ("IP2", "CF2") for m in ("IP2", "CF2")
    }
    reliable_names = set(reliable["family"].astype(str))
    for family in families:
        if family in protected and family in reliable_names and family not in ordered_survivors:
            ordered_survivors.append(family)
        #endif
    #endfor

    summary["survives_pilot"] = summary["family"].astype(str).isin(ordered_survivors)
    summary["pilot_cutoff_fm"] = cutoff
    summary.to_csv(outdir / "radius_bias_pilot_ranking.csv", index=False)
    print(
        f"[radius-bias pilot] retained {len(ordered_survivors)}/{len(families)} "
        f"family pairs for the full study; best={best:.5f} fm, cutoff={cutoff:.5f} fm"
    )
    print("[radius-bias pilot] survivors: " + ", ".join(ordered_survivors))
    return ordered_survivors
#enddef


def run_radius_bias_variance_study(
        bundles: Sequence[Dict[str, object]],
        args,
        outdir: Path) -> None:
    """
    Hayward/Griffioen-style extrapolation optimization adapted to BH data.

    The trial parameters are normalized Sachs GE(Q2) and GM(Q2)/mu_p shapes,
    but they are fitted through the BH cross section, not to direct GE/GM
    pseudodata. At each measured point the trial GE/GM are converted to F1/F2
    and inserted into the exact BH quadratic.

    Production fit candidates:
      P2/P3, ID1/P2/P3/IP1/IP2/IP3/CF2/CF3, independently paired for GE and GM.
      A 10-replica pilot prunes clearly noncompetitive pairs before the full
      requested-statistics closure stage. Higher-order P4/IP4/CF4 forms remain
      in the truth ensemble as curvature stress tests.

    Truth ensemble:
      empirical Kelly/AMT/Bernauer truths plus optional synthetic D/P/IP/CF
      closure truths spanning a Cartesian rE x rM grid.  The extended ensemble
      also varies exponential/dipole-like higher-order curvature independently
      in GE and GM at fixed central radius, without a costly full Cartesian
      radius x curvature expansion.  Kelly, AMT2007, and
      Bernauer are treated as three realizations of ONE empirical-elastic
      truth_group rather than three independent conceptual votes. Production
      closure aggregation is hierarchical: every truth_group receives equal
      total weight, independent of the number of scenarios in that group.
    """
    outdir.mkdir(parents=True, exist_ok=True)

    # Exploratory fit-candidate set. P1 is omitted because a linear polynomial
    # has not provided an adequate description. Fourth-order candidates remain
    # truth generators rather than production fits. The 64 E x M pairs below
    # first pass through a cheap pilot closure screen; only competitive,
    # numerically reliable survivors enter the expensive full-replica study.
    candidate_families = [
        "D1",
        "P2", "P3",
        "IP1", "IP2", "IP3",
        "CF2", "CF3",
    ]
    families = [
        encode_sachs_family_pair(fe, fm)
        for fe in candidate_families for fm in candidate_families
    ]

    # Keep the broader generating-family set for closure stress tests.  Truth
    # generation is cheap compared with fitting every candidate to every
    # replica, and retaining D1 and P4/IP4/CF4 here tests whether the reduced
    # candidate set remains robust against more complicated underlying shapes.
    truth_families = ["D1"]
    truth_families += [f"P{i}" for i in range(1, 5)]
    truth_families += [f"IP{i}" for i in range(1, 5)]
    truth_families += [f"CF{i}" for i in range(2, 5)]

    specs_all = [bundle_to_measurement_spec(b) for b in bundles]
    empty_specs = [spec for spec in specs_all if len(spec["data"]) == 0]
    specs = [spec for spec in specs_all if len(spec["data"]) > 0]

    for spec in empty_specs:
        print(
            f"[radius-bias] {spec['label']}: 0 selected points; "
            "omitting this measurement from the closure study"
        )
    #endfor

    if not specs:
        raise RuntimeError(
            "Radius-bias study has zero selected points across all measurements"
        )
    #endif

    qmax = max(
        float(np.nanmax(spec["data"]["t_abs"].to_numpy(float)))
        for spec in specs
    )

    scenarios = [
        {
            "truth_model": "Kelly",
            "truth_group": EMPIRICAL_VALIDATION_TRUTH_GROUP,
            "truth_family": "Kelly",
            "truth_fn": kelly_sachs,
            "synthetic": False,
        },
        {
            "truth_model": "AMT2007",
            "truth_group": EMPIRICAL_VALIDATION_TRUTH_GROUP,
            "truth_family": "AMT2007",
            "truth_fn": amt2007_sachs,
            "synthetic": False,
        },
        {
            "truth_model": "Bernauer_order8_polyxdipole",
            "truth_group": EMPIRICAL_VALIDATION_TRUTH_GROUP,
            "truth_family": "Bernauer",
            "truth_fn": bernauer_polyxdipole_sachs,
            "synthetic": False,
        },
    ]

    if args.radius_bias_extended_truths:
        radius_values = parse_radius_bias_grid(args.radius_bias_radius_grid)
        scenarios += make_synthetic_sachs_truth_scenarios(
            families=truth_families,
            radius_values=radius_values,
            qmax=qmax,
            curvature_stress=args.radius_bias_curvature_stress,
        )
    #endif

    save_empirical_truth_similarity_diagnostics(qmax, outdir)

    families_full = list(families)
    families = pilot_prune_radius_bias_families(
        specs, scenarios, families_full, args, outdir
    )

    n_empirical_validation = sum(
        str(sc.get("truth_group", ""))
        == EMPIRICAL_VALIDATION_TRUTH_GROUP
        for sc in scenarios
    )
    n_production_truths = len(scenarios) - int(n_empirical_validation)
    print(
        f"[radius-bias] total truth scenarios={len(scenarios)} "
        f"({n_production_truths} controlled synthetic production truths + "
        f"{n_empirical_validation} empirical validation truths); "
        f"fit families={len(families)}, replicas/scenario/family="
        f"{args.radius_bias_replicas}"
    )
    print(
        "[radius-bias] Kelly/AMT/Bernauer are validation only: "
        "zero weight in ranking, production bias, or P(r) model weights"
    )

    rows = []
    replica_rows = []

    for truth_index, scenario in enumerate(scenarios):
        truth_name = str(scenario["truth_model"])
        truth_group = str(scenario["truth_group"])
        truth_fn = scenario["truth_fn"]

        if bool(scenario.get("synthetic", False)):
            truth_radius = {
                "E": float(scenario["truth_rE_fm"]),
                "M": float(scenario["truth_rM_fm"]),
            }
        else:
            truth_radius = {}
            for quantity in ["E", "M"]:
                norm = 1.0 if quantity == "E" else MU_P
                truth_radius[quantity] = radius_from_shape(
                    lambda qq, qn=quantity: (
                        truth_fn(np.asarray(qq, dtype=float))[0]
                        if qn == "E"
                        else truth_fn(np.asarray(qq, dtype=float))[1]
                    ),
                    norm,
                )
            #endfor
        #endif

        central_by_key = {}
        sigma_by_key = {}
        truth_ge_by_key = {}
        truth_gm_by_key = {}
        for spec in specs:
            d = spec["data"]
            q = d["t_abs"].to_numpy(float)
            ge, gm = truth_fn(q)
            truth_ge_by_key[str(spec["key"])] = np.asarray(ge, dtype=float)
            truth_gm_by_key[str(spec["key"])] = np.asarray(gm, dtype=float) / MU_P
            tau = q / (4.0 * MP2)
            f1 = (ge + tau * gm) / (1.0 + tau)
            f2 = (gm - ge) / (1.0 + tau)
            central_by_key[str(spec["key"])] = bh_from_f1f2(
                d["bh_A"].to_numpy(float),
                d["bh_B"].to_numpy(float),
                d["bh_C"].to_numpy(float),
                f1, f2,
            )
            sigma_by_key[str(spec["key"])] = dataset_statistical_errors(
                d, str(spec["kind"])
            )
        #endfor

        nworkers = (
            args.radius_bias_workers
            if args.radius_bias_workers is not None
            else args.workers
        )
        nworkers = max(1, int(nworkers))

        # Treat --radius-bias-replicas as the requested number of SUCCESSFUL
        # replica fits, rather than the number merely attempted.  Failed fits
        # are replaced by fresh statistically independent replicas until the
        # target is reached or a finite attempt budget is exhausted.  The
        # attempt efficiency remains an explicit numerical-robustness metric.
        target_valid = max(1, int(args.radius_bias_replicas))
        min_eff = max(0.50, min(0.99, float(
            getattr(args, "radius_bias_min_scenario_valid_fraction", 0.80)
        )))
        # Do not spend up to 5x the target rescuing a pathological family. A
        # family that cannot achieve the configured acceptable convergence
        # fraction within a modest safety margin should fail the closure gate.
        max_attempts = max(
            target_valid + 10,
            int(math.ceil(1.10 * target_valid / min_eff)),
        )
        seed_root = np.random.SeedSequence(
            [int(args.radius_bias_seed), int(truth_index)]
        )
        family_seed_roots = seed_root.spawn(len(families))
        seed_queues = {
            family: list(family_seed_roots[ifam].spawn(max_attempts))
            for ifam, family in enumerate(families)
        }

        results_by_family = {
            family: {"rE": [], "rM": [], "shape_GE": [], "shape_GM": [],
                     "shape_combined": [], "nvalid": 0, "nattempted": 0,
                     "target_reached": False}
            for family in families
        }

        print(
            f"[radius-bias] truth={truth_name}: target={target_valid} valid "
            f"replicas/family, max_attempts={max_attempts}, "
            f"{len(families)} family batches with {nworkers} worker(s)..."
        )
        t0 = time.time()

        # Dispatch one retry round at a time.  Each family receives only as
        # many new replicas as it still needs, capped by its remaining attempt
        # budget.  This preserves family-level batching while avoiding wasted
        # fits after a family has already accumulated the requested successes.
        with ProcessPoolExecutor(
                max_workers=nworkers,
                initializer=_init_radius_bias_worker,
                initargs=(specs, central_by_key, sigma_by_key,
                          truth_ge_by_key, truth_gm_by_key)) as pool:
            while True:
                tasks = []
                task_attempt_starts = {}
                for family in families:
                    state = results_by_family[family]
                    need = target_valid - int(state["nvalid"])
                    remaining = max_attempts - int(state["nattempted"])
                    if need <= 0 or remaining <= 0:
                        continue
                    #endif
                    nsubmit = min(need, remaining)
                    i0 = int(state["nattempted"])
                    seeds = [
                        int(ss.generate_state(1, dtype=np.uint64)[0])
                        for ss in seed_queues[family][i0:i0 + nsubmit]
                    ]
                    tasks.append((family, seeds))
                    task_attempt_starts[family] = i0
                #endfor

                if not tasks:
                    break
                #endif

                for family, family_results in pool.map(
                        _radius_bias_family_batch_worker, tasks, chunksize=1):
                    state = results_by_family[family]
                    attempt_start = task_attempt_starts[family]
                    for ilocal, result in enumerate(family_results):
                        _, re_val, rm_val, valid, ge_shape, gm_shape, combined_shape = result
                        attempt_index = attempt_start + ilocal
                        state["nattempted"] += 1
                        accepted = bool(valid and int(state["nvalid"]) < target_valid)
                        replica_rows.append({
                            "truth_model": truth_name, "truth_group": truth_group,
                            "truth_rE_fm": truth_radius["E"], "truth_rM_fm": truth_radius["M"],
                            "family": family, "attempt": attempt_index,
                            "accepted_replica": int(state["nvalid"]) if accepted else np.nan,
                            "rE_fit_fm": re_val, "rM_fit_fm": rm_val, "valid": valid,
                            "GE_shape_frac_rms": ge_shape, "GM_shape_frac_rms": gm_shape,
                            "combined_shape_frac_rms": combined_shape,
                        })
                        if accepted:
                            state["rE"].append(re_val)
                            state["rM"].append(rm_val)
                            state["shape_GE"].append(ge_shape)
                            state["shape_GM"].append(gm_shape)
                            state["shape_combined"].append(combined_shape)
                            state["nvalid"] += 1
                        #endif
                    #endfor
                #endfor
            #endwhile
        #endwith

        for family in families:
            results_by_family[family]["target_reached"] = bool(
                int(results_by_family[family]["nvalid"]) >= target_valid
            )
        #endfor

        print(
            f"[radius-bias] truth={truth_name}: finished in "
            f"{time.time() - t0:.1f} s"
        )

        for family in families:
            re_vals = results_by_family[family]["rE"]
            rm_vals = results_by_family[family]["rM"]
            nvalid = int(results_by_family[family]["nvalid"])
            shape_ge = np.asarray(results_by_family[family]["shape_GE"], dtype=float)
            shape_gm = np.asarray(results_by_family[family]["shape_GM"], dtype=float)
            shape_combined = np.asarray(results_by_family[family]["shape_combined"], dtype=float)
            shape_summary = {
                "GE_shape_fractional_RMS_mean": float(np.mean(shape_ge)) if len(shape_ge) else np.nan,
                "GM_shape_fractional_RMS_mean": float(np.mean(shape_gm)) if len(shape_gm) else np.nan,
                "combined_shape_fractional_RMS_mean": float(np.mean(shape_combined)) if len(shape_combined) else np.nan,
                "combined_shape_fractional_RMS_median": float(np.median(shape_combined)) if len(shape_combined) else np.nan,
            }

            for quantity, vals, rtrue in [
                ("rE", re_vals, truth_radius["E"]),
                ("rM", rm_vals, truth_radius["M"]),
            ]:
                arr = np.asarray(vals, dtype=float)
                if len(arr) == 0:
                    mean = stat = bias = total = np.nan
                else:
                    mean = float(np.mean(arr))
                    stat = float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0
                    bias = float(mean - rtrue)
                    total = float(math.sqrt(stat**2 + bias**2))
                #endif

                rows.append({
                    "truth_model": truth_name,
                    "truth_group": truth_group,
                    "truth_family": str(
                        scenario.get("truth_family", truth_group)
                    ),
                    "synthetic_truth": bool(
                        scenario.get("synthetic", False)
                    ),
                    "family": family,
                    "quantity": quantity,
                    "truth_radius_fm": rtrue,
                    "truth_rE_fm": truth_radius["E"],
                    "truth_rM_fm": truth_radius["M"],
                    "mean_extracted_radius_fm": mean,
                    "stat_RMS_fm": stat,
                    "bias_fm": bias,
                    "sqrt_stat2_plus_bias2_fm": total,
                    "valid_replicas": nvalid,
                    "requested_replicas": target_valid,
                    "attempted_replicas": int(results_by_family[family]["nattempted"]),
                    "attempt_efficiency": (
                        nvalid / int(results_by_family[family]["nattempted"])
                        if int(results_by_family[family]["nattempted"]) > 0 else np.nan
                    ),
                    "target_reached": bool(results_by_family[family]["target_reached"]),
                    "valid_fraction": (nvalid / target_valid),
                    "max_attempts": max_attempts,
                    "workers": nworkers,
                    **shape_summary,
                })
            #endfor

            print(
                f"[radius-bias] truth={truth_name:42s} family={family:15s} "
                f"valid={nvalid}/{target_valid} "
                f"attempted={int(results_by_family[family]['nattempted'])} "
                f"eff={nvalid / max(1, int(results_by_family[family]['nattempted'])):.3f}"
            )
        #endfor
    #endfor

    table = pd.DataFrame(rows)
    table.to_csv(outdir / "radius_bias_variance_study.csv", index=False)

    save_empirical_closure_validation_summary(table, outdir)

    replica_table = pd.DataFrame(replica_rows)
    replica_table.to_csv(
        outdir / "radius_bias_replica_results.csv",
        index=False,
    )

    agg_rows = []
    for family in families:
        for quantity in ["rE", "rM"]:
            part = production_closure_truth_rows(table.loc[
                (table["family"] == family)
                & (table["quantity"] == quantity)
            ])
            vals = part["sqrt_stat2_plus_bias2_fm"].to_numpy(float)
            vals = vals[np.isfinite(vals)]
            agg_rows.append({
                "family": family,
                "quantity": quantity,
                "RMS_objective_across_truths_fm": (
                    float(np.sqrt(np.mean(vals**2))) if len(vals) else np.nan
                ),
                "max_objective_across_truths_fm": (
                    float(np.max(vals)) if len(vals) else np.nan
                ),
            })
        #endfor
    #endfor

    agg = pd.DataFrame(agg_rows)
    agg.to_csv(outdir / "radius_bias_variance_ranking.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.8))
    for ax, quantity in zip(axes, ["rE", "rM"]):
        part = agg.loc[agg["quantity"] == quantity].copy()
        x = np.arange(len(part))
        ax.plot(
            x,
            part["RMS_objective_across_truths_fm"],
            marker="o",
        )
        ax.set_xticks(x)
        ax.set_xticklabels(part["family"], rotation=45)
        ax.set_ylabel(r"$\sqrt{\sigma_{\rm stat}^2+b^2}$ (fm)")
        ax.set_title(quantity)
        ax.set_ylim(
            0.0,
            robust_upper_limit(
                part["RMS_objective_across_truths_fm"].to_numpy(float),
                quantile=0.90,
                pad_fraction=0.15,
                minimum=0.05,
            ),
        )
        ax.grid(alpha=0.2)
    #endfor
    fig.suptitle(
        "Radius extrapolation bias-variance optimization across truth models",
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(outdir / "01_radius_bias_variance_ranking.png", dpi=180)
    plt.close(fig)

    empirical_table = table.loc[~table["synthetic_truth"]].copy()
    for quantity in ["rE", "rM"]:
        fig, ax = plt.subplots(figsize=(10.5, 5.2))
        x = np.arange(len(families))
        for truth_name in empirical_table["truth_model"].unique():
            part = empirical_table.loc[
                (empirical_table["truth_model"] == truth_name)
                & (empirical_table["quantity"] == quantity)
            ].set_index("family").reindex(families)
            ax.plot(
                x,
                part["sqrt_stat2_plus_bias2_fm"],
                marker="o",
                label=truth_name,
            )
        #endfor
        ax.set_xticks(x)
        ax.set_xticklabels(families, rotation=45)
        ax.set_ylabel(r"$\sqrt{\sigma_{\rm stat}^2+b^2}$ (fm)")
        ax.set_title(f"{quantity}: empirical-truth objective")
        plotted_vals = empirical_table.loc[
            empirical_table["quantity"] == quantity,
            "sqrt_stat2_plus_bias2_fm",
        ].to_numpy(float)
        ax.set_ylim(
            0.0,
            robust_upper_limit(
                plotted_vals,
                quantile=0.95,
                pad_fraction=0.12,
                minimum=0.05,
            ),
        )
        ax.grid(alpha=0.2)
        ax.legend()
        fig.tight_layout()
        fig.savefig(
            outdir / f"02_{quantity}_objective_empirical_truths.png",
            dpi=180,
        )
        plt.close(fig)
    #endfor

    for quantity in ["rE", "rM"]:
        fig, ax = plt.subplots(figsize=(10.5, 5.2))
        x = np.arange(len(families))
        part_q = table.loc[table["quantity"] == quantity]
        stat_rms = []
        abs_bias_rms = []
        for family in families:
            fam = part_q.loc[part_q["family"] == family]
            stat_vals = fam["stat_RMS_fm"].to_numpy(float)
            bias_vals = np.abs(fam["bias_fm"].to_numpy(float))
            stat_vals = stat_vals[np.isfinite(stat_vals)]
            bias_vals = bias_vals[np.isfinite(bias_vals)]
            stat_rms.append(
                float(np.sqrt(np.mean(stat_vals**2)))
                if len(stat_vals) else np.nan
            )
            abs_bias_rms.append(
                float(np.sqrt(np.mean(bias_vals**2)))
                if len(bias_vals) else np.nan
            )
        #endfor
        ax.plot(x, stat_rms, marker="o", label="statistical RMS")
        ax.plot(x, abs_bias_rms, marker="s", label="RMS |bias|")
        ax.set_xticks(x)
        ax.set_xticklabels(families, rotation=45)
        ax.set_ylabel("fm")
        ax.set_title(
            f"{quantity}: statistical variance versus extrapolation bias"
        )
        ax.set_ylim(
            0.0,
            robust_upper_limit(
                np.asarray(stat_rms + abs_bias_rms, dtype=float),
                quantile=0.90,
                pad_fraction=0.15,
                minimum=0.05,
            ),
        )
        ax.grid(alpha=0.2)
        ax.legend()
        fig.tight_layout()
        fig.savefig(outdir / f"03_{quantity}_stat_vs_bias.png", dpi=180)
        plt.close(fig)
    #endfor

    fig, ax = plt.subplots(figsize=(10.5, 5.2))
    x = np.arange(len(families))
    min_frac = []
    global_frac = []
    for family in families:
        part = table.loc[
            (table["family"] == family)
            & (table["quantity"] == "rE")
        ]
        attempted = part["attempted_replicas"].to_numpy(float)
        valid = part["valid_replicas"].to_numpy(float)
        fractions = np.divide(
            valid, attempted, out=np.zeros_like(valid), where=attempted > 0.0
        )
        min_frac.append(float(np.nanmin(fractions)))
        global_frac.append(
            float(np.nansum(valid) / np.nansum(attempted))
            if np.nansum(attempted) > 0.0 else np.nan
        )
    #endfor

    ax.plot(x, global_frac, marker="o", label="global attempt efficiency")
    ax.plot(x, min_frac, marker="s", label="minimum truth-scenario efficiency")
    # Attempt efficiency is the numerical-robustness gate; closure bias and
    # replica variance determine the ranking among families that pass it.
    ax.set_xticks(x)
    ax.set_xticklabels(families, rotation=45)
    ax.set_ylim(0.0, 1.05)
    ax.set_ylabel("successful fits / attempted fits")
    ax.set_title("Replica-fit attempt efficiency by candidate family")
    ax.grid(alpha=0.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "04_replica_attempt_efficiency.png", dpi=180)
    plt.close(fig)

    for quantity in ["rE", "rM"]:
        fig, ax = plt.subplots(figsize=(11.0, 5.6))
        x = np.arange(len(families), dtype=float)
        empirical_names = empirical_table["truth_model"].unique().tolist()
        offsets = np.linspace(-0.18, 0.18, len(empirical_names))
        for offset, truth_name in zip(offsets, empirical_names):
            part = empirical_table.loc[
                (empirical_table["truth_model"] == truth_name)
                & (empirical_table["quantity"] == quantity)
            ].set_index("family").reindex(families)
            ax.errorbar(
                x + offset,
                part["mean_extracted_radius_fm"],
                yerr=part["stat_RMS_fm"],
                fmt="o",
                capsize=2,
                label=truth_name,
            )
        #endfor
        ax.set_xticks(x)
        ax.set_xticklabels(families, rotation=45)
        ax.set_ylabel(f"{quantity} extracted radius (fm)")
        ax.set_title(
            f"{quantity}: empirical truths, mean extracted radius ± replica RMS"
        )
        qpart = empirical_table.loc[
            empirical_table["quantity"] == quantity
        ]
        lo, hi = robust_radius_limits(
            qpart["mean_extracted_radius_fm"].to_numpy(float),
            truth_values=qpart["truth_radius_fm"].to_numpy(float),
            lower_quantile=0.02,
            upper_quantile=0.98,
            pad_fraction=0.10,
        )
        ax.set_ylim(lo, hi)
        ax.grid(alpha=0.2)
        ax.legend()
        fig.tight_layout()
        fig.savefig(
            outdir / f"05_{quantity}_mean_extracted_radius_empirical.png",
            dpi=180,
        )
        plt.close(fig)
    #endfor

    # Diagnostic only: full-range GE/GM shape recovery over the actual measured
    # |t| support.  This does NOT enter the baseline family ranking yet.
    shape_rows = []
    for family in families:
        fam = table.loc[(table["family"] == family) & (table["quantity"] == "rE")].copy()
        for col, label in [("GE_shape_fractional_RMS_mean", "GE"), ("GM_shape_fractional_RMS_mean", "GM"), ("combined_shape_fractional_RMS_mean", "combined")]:
            vals = fam[col].to_numpy(float) if col in fam.columns else np.asarray([])
            vals = vals[np.isfinite(vals)]
            shape_rows.append({
                "family": family, "quantity": label,
                "RMS_across_truth_scenarios": float(np.sqrt(np.mean(vals**2))) if len(vals) else np.nan,
                "mean_across_truth_scenarios": float(np.mean(vals)) if len(vals) else np.nan,
                "N_truth_scenarios_with_metric": int(len(vals)),
            })
        #endfor
    #endfor
    shape_table = pd.DataFrame(shape_rows)
    shape_table.to_csv(outdir / "11_full_range_shape_recovery_diagnostic.csv", index=False)
    fig, ax = plt.subplots(figsize=(10.5, 5.3))
    x = np.arange(len(families))
    for quantity, marker in [("GE", "o"), ("GM", "s"), ("combined", "^")]:
        part = shape_table.loc[shape_table["quantity"] == quantity].set_index("family").reindex(families)
        ax.plot(x, 100.0 * part["RMS_across_truth_scenarios"], marker=marker, label=quantity)
    #endfor
    ax.set_xticks(x)
    ax.set_xticklabels(families, rotation=45)
    ax.set_ylabel("RMS fractional shape error over measured |t| (%)")
    ax.set_title("Diagnostic full-range form-factor closure (not used for ranking)")
    ax.grid(alpha=0.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "11_full_range_shape_recovery_diagnostic.png", dpi=180)
    plt.close(fig)

    if args.radius_bias_extended_truths:
        save_mixed_family_closure_ranking(
            table=table,
            candidate_pairs=families,
            outdir=outdir,
            minimum_global_valid_fraction=args.radius_bias_min_global_valid_fraction,
            minimum_scenario_valid_fraction=args.radius_bias_min_scenario_valid_fraction,
        )
    #endif

    print(f"[radius-bias] tables and plots -> {outdir}")
#enddef


def save_bh_local_f1_f2_sensitivity(
        outdir: Path,
        reference_data: pd.DataFrame,
        reference_fit: FitResult) -> None:
    """
    Local logarithmic BH sensitivities at every measured Set-5 point.

      S1 = d ln sigma_BH / d ln F1
      S2 = d ln sigma_BH / d ln F2

    A value S1=1.8 means a 1% fractional change in F1 changes that point's
    BH prediction by about 1.8%.
    """
    q = reference_data["t_abs"].to_numpy(float)
    A = reference_data["bh_A"].to_numpy(float)
    B = reference_data["bh_B"].to_numpy(float)
    C = reference_data["bh_C"].to_numpy(float)
    f1, f2, _, _ = evaluate_fit_form_factors(reference_fit, q)
    sigma = bh_from_f1f2(A, B, C, f1, f2)

    s1 = f1 * (2.0 * A * f1 + B * f2) / sigma
    s2 = f2 * (B * f1 + 2.0 * C * f2) / sigma

    diag = reference_data[["xB", "Q2", "t_abs", "phi_deg"]].copy()
    diag["S_F1"] = s1
    diag["S_F2"] = s2
    diag["abs_S_F1"] = np.abs(s1)
    diag["abs_S_F2"] = np.abs(s2)
    diag["F2_to_F1_sensitivity_ratio"] = np.abs(s2) / np.maximum(np.abs(s1), 1e-12)
    diag.to_csv(outdir / "11_BH_local_F1_F2_sensitivity.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))
    sc1 = axes[0].scatter(
        reference_data["t_abs"], s1,
        c=reference_data["phi_deg"], s=16,
    )
    sc2 = axes[1].scatter(
        reference_data["t_abs"], s2,
        c=reference_data["phi_deg"], s=16,
    )
    axes[0].set_ylabel(r"$\partial\ln\sigma_{\rm BH}/\partial\ln F_1$")
    axes[1].set_ylabel(r"$\partial\ln\sigma_{\rm BH}/\partial\ln F_2$")
    for ax in axes:
        ax.set_xlabel(r"$|t|$ (GeV$^2$)")
        ax.axhline(0.0, linewidth=0.8)
        ax.grid(alpha=0.2)
    #endfor
    fig.subplots_adjust(
        top=0.90, bottom=0.12, left=0.08, right=0.84, wspace=0.25
    )
    cbar = fig.colorbar(
        sc2, ax=axes.ravel().tolist(), pad=0.06, fraction=0.030
    )
    cbar.set_label(r"$\phi$ (deg)")
    fig.suptitle("Local BH sensitivity to Dirac and Pauli form factors", y=0.995)
    fig.savefig(outdir / "11_BH_local_F1_F2_sensitivity.png", dpi=180)
    plt.close(fig)
#enddef


def save_bh_f1_f2_sensitivity(outdir: Path,
                              ebeam: float,
                              reference_data: pd.DataFrame) -> None:
    """
    Paper Fig. 7 analogue.

    Use a representative CLAS12 Set-5 kinematic point rather than importing
    elastic-scattering data. The point chosen is the Set-5 event closest to
    the medians in xB, Q2, and phi.
    """
    try:
        import gepard as g
        from gepard.fits import th_KM15
    except Exception:
        return
    #endtry

    med_x = float(reference_data["xB"].median())
    med_q2 = float(reference_data["Q2"].median())
    med_phi = float(reference_data["phi_deg"].median())

    scale_x = max(float(reference_data["xB"].std()), 1.0e-6)
    scale_q2 = max(float(reference_data["Q2"].std()), 1.0e-6)
    scale_phi = max(float(reference_data["phi_deg"].std()), 1.0e-6)

    metric = (
        ((reference_data["xB"] - med_x) / scale_x)**2
        + ((reference_data["Q2"] - med_q2) / scale_q2)**2
        + ((reference_data["phi_deg"] - med_phi) / scale_phi)**2
    )
    row = reference_data.loc[metric.idxmin()]

    qmax = max(0.45, float(reference_data["t_abs"].max()))
    q = np.linspace(0.01, qmax, 220)
    full = np.zeros_like(q)
    f1zero = np.zeros_like(q)
    f2zero = np.zeros_like(q)

    for i, qi in enumerate(q):
        th = th_KM15
        pt = make_gepard_point(
            g,
            float(row["xB"]),
            float(row["Q2"]),
            float(qi),
            float(row["phi_deg"]),
            ebeam,
        )
        pref = float(th.PreFacSigma(pt))
        model = th.m
        old_f1, old_f2 = model.F1, model.F2
        f1k, f2k = kelly_f1_f2(np.array([qi]))

        try:
            model.F1 = lambda p, value=float(f1k[0]): value
            model.F2 = lambda p, value=float(f2k[0]): value
            full[i] = pref * float(th.TBH2unp(pt))

            model.F1 = lambda p: 0.0
            f1zero[i] = pref * float(th.TBH2unp(pt))

            model.F1 = lambda p, value=float(f1k[0]): value
            model.F2 = lambda p: 0.0
            f2zero[i] = pref * float(th.TBH2unp(pt))
        finally:
            model.F1, model.F2 = old_f1, old_f2
        #endtry
    #endfor

    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    ax.plot(q, full, label="Full BH")
    ax.plot(q, f1zero, "--", label=r"$F_1=0$")
    ax.plot(q, f2zero, ":", label=r"$F_2=0$")
    ax.set_xlabel(r"$|t|$ (GeV$^2$)")
    ax.set_ylabel(r"$d^4\sigma$")
    ax.set_title(
        rf"BH sensitivity: $x_B={row['xB']:.3f}$, "
        rf"$Q^2={row['Q2']:.2f}$ GeV$^2$, "
        rf"$\phi={row['phi_deg']:.1f}^\circ$"
    )
    ax.legend()
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(outdir / "07_BH_F1_F2_sensitivity.png")
    plt.close(fig)
#enddef


def save_fit5_to_fit8_sachs_plot(results: List[FitResult],
                                 set5_data: pd.DataFrame,
                                 outdir: Path) -> None:
    """Paper Fig. 8 analogue: Fits 5-8 in GE and GM."""
    names = ["Fit 5", "Fit 6", "Fit 7", "Fit 8"]
    fits = [next(r for r in results if r.name == name) for name in names]
    q = np.linspace(0.0, 1.0, 400)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))

    for fr in fits:
        _, _, ge, gm = evaluate_fit_form_factors(fr, q)
        axes[0].plot(q, ge, linewidth=1.3, label=fr.name)
        axes[1].plot(q, gm, linewidth=1.3, label=fr.name)
    #endfor

    # Show only the Fit-8 uncertainty band to keep the comparison legible.
    for ax, which in zip(axes, ["GE", "GM"]):
        central, sigma = form_factor_band(fits[-1], q, which)
        ax.fill_between(
            q,
            central - sigma,
            central + sigma,
            alpha=0.18,
            label="Fit 8 68% Hessian band",
        )
        draw_t_reference_lines(ax, set5_data)
        ax.set_xlabel(r"$|t|$ (GeV$^2$)")
        ax.set_xlim(0.0, 1.0)
        ax.grid(alpha=0.2)
    #endfor

    axes[0].set_ylabel(r"$G_E(t)$")
    axes[1].set_ylabel(r"$G_M(t)$")
    axes[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "08_GE_GM_fits_5_to_8.png")
    plt.close(fig)
#enddef


def save_radii_plot(results: List[FitResult],
                    outdir: Path) -> None:
    """Paper Table-V-style visual summary for Fits 5-8."""
    rows = [
        next(r for r in results if r.name == f"Fit {i}")
        for i in range(5, 9)
    ]
    labels = [r.name for r in rows]
    y = np.arange(len(rows))

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.8), sharey=True)

    axes[0].errorbar(
        [r.rE_fm for r in rows],
        y,
        xerr=[r.rE_err_fm for r in rows],
        fmt="o",
        capsize=2,
    )
    axes[1].errorbar(
        [r.rM_fm for r in rows],
        y,
        xerr=[r.rM_err_fm for r in rows],
        fmt="o",
        capsize=2,
    )

    axes[0].set_yticks(y, labels)
    axes[0].invert_yaxis()
    axes[0].set_xlabel(r"$r_E$ (fm)")
    axes[1].set_xlabel(r"$r_M$ (fm)")
    axes[0].grid(alpha=0.2)
    axes[1].grid(alpha=0.2)

    fig.tight_layout()
    fig.savefig(outdir / "09_radii_fits_5_to_8.png")
    plt.close(fig)
#enddef


def save_summary_tables(results: List[FitResult],
                        selection_records: List[Dict[str, float]],
                        outdir: Path) -> None:
    """Write compact paper-like Tables I-V to CSV files."""
    pd.DataFrame(selection_records).to_csv(
        outdir / "table_I_bh_datasets.csv",
        index=False,
    )

    fits15 = [
        next(r for r in results if r.name == f"Fit {i}")
        for i in range(1, 6)
    ]
    pd.DataFrame([
        {
            "Fit": fr.name,
            "BH criterion": f"<= {100*fr.bh_cut:.0f}%",
            "Npts": fr.npts,
            "chi2": fr.chi2,
            "dof": fr.ndof,
            "chi2/dof": fr.chi2_ndof,
        }
        for fr in fits15
    ]).to_csv(outdir / "table_II_fits_1_to_5_chi2.csv", index=False)

    pd.DataFrame([
        {
            "Fit": fr.name,
            "aE_GeV2": fr.params[0],
            "aE_err_GeV2": math.sqrt(fr.covariance[0, 0]),
            "aM_GeV2": fr.params[1],
            "aM_err_GeV2": math.sqrt(fr.covariance[1, 1]),
        }
        for fr in fits15
    ]).to_csv(outdir / "table_III_dipole_parameters.csv", index=False)

    fit6 = next(r for r in results if r.name == "Fit 6")
    fit7 = next(r for r in results if r.name == "Fit 7")
    fit8 = next(r for r in results if r.name == "Fit 8")

    rows = []
    for fr in [fit6, fit7, fit8]:
        row = {
            "Fit": fr.name,
            "chi2/dof": fr.chi2_ndof,
        }
        for i, name in enumerate(fr.param_names):
            row[name] = fr.params[i]
            row[name + "_err"] = math.sqrt(
                max(fr.covariance[i, i], 0.0)
            )
        #endfor
        rows.append(row)
    #endfor
    pd.DataFrame(rows).to_csv(
        outdir / "table_IV_ppole_and_fit8_parameters.csv",
        index=False,
    )

    pd.DataFrame([
        {
            "Fit": fr.name,
            "rE_fm": fr.rE_fm,
            "rE_err_fm": fr.rE_err_fm,
            "rM_fm": fr.rM_fm,
            "rM_err_fm": fr.rM_err_fm,
        }
        for fr in [fits15[-1], fit6, fit7, fit8]
    ]).to_csv(outdir / "table_V_radii.csv", index=False)
#enddef





def get_named_fit(bundle: Dict[str, object], fit_name: str) -> FitResult:
    return next(r for r in bundle["results"] if r.name == fit_name)
#enddef


def save_cross_dataset_form_factor_comparison(
        bundles: Sequence[Dict[str, object]],
        outdir: Path,
        fit_name: str) -> None:
    """
    Overlay CLAS6, CLAS12 Lee 2026, and CLAS12 pass-2 fitted form factors
    with propagated 68% Hessian bands. No KM15/Gepard calls are made here.
    """
    observed_max = max(
        float(bundle["set5"]["t_abs"].max()) for bundle in bundles
    )
    q = np.linspace(0.0, max(1.0, observed_max), 500)

    fig, axes = plt.subplots(2, 2, figsize=(12.0, 9.0))
    panels = [
        ("F1", r"$F_1(t)$"),
        ("F2", r"$F_2(t)$"),
        ("GE", r"$G_E(t)$"),
        ("GM", r"$G_M(t)$"),
    ]

    for ax, (which, ylabel) in zip(axes.flat, panels):
        for bundle in bundles:
            fr = get_named_fit(bundle, fit_name)
            central, sigma = form_factor_band(fr, q, which)
            line, = ax.plot(q, central, linewidth=1.5, label=bundle["label"])
            ax.fill_between(
                q, central - sigma, central + sigma,
                alpha=0.18, color=line.get_color()
            )

            tmin = float(bundle["set5"]["t_abs"].min())
            tmax = float(bundle["set5"]["t_abs"].max())
            ax.axvline(tmin, linewidth=0.55, alpha=0.55,
                       color=line.get_color())
            ax.axvline(tmax, linewidth=0.55, alpha=0.55,
                       color=line.get_color())
        #endfor

        ax.axvline(A1_BERNAUER_Q2_MIN, linewidth=0.7, linestyle="--")
        ax.axvline(PRAD_Q2_MIN, linewidth=0.7, linestyle=":")
        ax.set_xlabel(r"$|t|$ (GeV$^2$)")
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.2)
    #endfor

    axes[0, 0].legend(fontsize=9)
    fig.suptitle(f"{fit_name}: BH-extracted proton form factors", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    safe = fit_name.lower().replace(" ", "_")
    fig.savefig(outdir / f"comparison_{safe}_form_factors.png", dpi=180)
    plt.close(fig)
#enddef


def save_cross_dataset_radius_comparison(
        bundles: Sequence[Dict[str, object]],
        outdir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 5.2))
    labels = [b["label"] for b in bundles]
    x = np.arange(len(labels), dtype=float)
    offsets = {"Fit 5": -0.10, "Fit 8": 0.10}

    for fit_name in ["Fit 5", "Fit 8"]:
        fits = [get_named_fit(b, fit_name) for b in bundles]
        axes[0].errorbar(
            x + offsets[fit_name],
            [f.rE_fm for f in fits],
            yerr=[f.rE_err_fm for f in fits],
            fmt="o", capsize=2, label=fit_name
        )
        axes[1].errorbar(
            x + offsets[fit_name],
            [f.rM_fm for f in fits],
            yerr=[f.rM_err_fm for f in fits],
            fmt="o", capsize=2, label=fit_name
        )
    #endfor

    for ax, ylabel in zip(axes, [r"$r_E$ (fm)", r"$r_M$ (fm)"]):
        ax.set_xticks(x, labels, rotation=15)
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.2)
    #endfor
    axes[0].legend()
    fig.tight_layout()
    fig.savefig(outdir / "comparison_radii_fit5_fit8.png", dpi=180)
    plt.close(fig)
#enddef


def save_cross_dataset_summary(
        bundles: Sequence[Dict[str, object]],
        outdir: Path) -> None:
    rows = []
    for bundle in bundles:
        for fit_name in ["Fit 5", "Fit 8"]:
            fr = get_named_fit(bundle, fit_name)
            rows.append({
                "dataset": bundle["label"],
                "fit": fit_name,
                "N_set5": len(bundle["set5"]),
                "chi2_ndof": fr.chi2_ndof,
                "rE_fm": fr.rE_fm,
                "rE_err_fm": fr.rE_err_fm,
                "rM_fm": fr.rM_fm,
                "rM_err_fm": fr.rM_err_fm,
                "tmin": float(bundle["set5"]["t_abs"].min()),
                "tmax": float(bundle["set5"]["t_abs"].max()),
            })
        #endfor
    #endfor
    pd.DataFrame(rows).to_csv(
        outdir / "comparison_fit5_fit8_summary.csv", index=False
    )
#enddef




def save_bh_cut_plateau_study_multi(
        bundles: Sequence[Dict[str, object]],
        outdir: Path,
        tag: str = "all_three") -> None:
    """
    Stability scan from 1% through 10% in |1-R_BH|.

    Every threshold is fit twice:
      * published_errors: published experimental point errors only;
      * moradi_errors: same errors plus an extra uncorrelated threshold*xsec
        uncertainty, following Moradi et al.

    This separation is essential because the latter mechanically inflates
    errors as the threshold is loosened and can drive chi2/ndf downward.
    """
    cuts = np.arange(0.01, 0.101, 0.01)
    rows = []

    for cut in cuts:
        selected_specs = []
        counts = {}
        for bundle in bundles:
            all_data = bundle["all_data"]
            selected = all_data.loc[all_data["bh_delta"] <= cut].copy()
            counts[str(bundle["key"])] = len(selected)
            selected_specs.append(
                bundle_to_measurement_spec(bundle, selected)
            )
        #endfor

        if any(len(spec["data"]) < 5 for spec in selected_specs):
            continue
        #endif

        for error_mode, add_moradi in [
            ("published_errors", False),
            ("moradi_errors", True),
        ]:
            for kind, fit_label in [
                ("dipole", "Fit 5 form"),
                ("fit8_f2_kelly", "Fit 8 form"),
            ]:
                fr = fit_multi_measurements(
                    selected_specs,
                    kind=kind,
                    fit_name=fit_label,
                    bh_cut=float(cut),
                    add_moradi_bh_systematic=add_moradi,
                )

                row = {
                    "BH_cut_fraction": float(cut),
                    "BH_cut_percent": 100.0 * float(cut),
                    "error_mode": error_mode,
                    "fit_form": fit_label,
                    "N_total": sum(counts.values()),
                    "chi2_ndof": fr.chi2_ndof,
                    "rE_fm": fr.rE_fm,
                    "rE_err_fm": fr.rE_err_fm,
                    "rM_fm": fr.rM_fm,
                    "rM_err_fm": fr.rM_err_fm,
                }
                for key, n in counts.items():
                    row["N_" + key] = n
                #endfor
                for mk, mv in fr.meta.items():
                    if mk.endswith("_scale_factor"):
                        row[mk] = mv
                    #endif
                #endfor
                rows.append(row)
            #endfor
        #endfor
    #endfor

    table = pd.DataFrame(rows)
    csv_name = f"BH_cut_plateau_{tag}.csv"
    table.to_csv(outdir / csv_name, index=False)

    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.2), sharex=True)
    line_styles = {
        "published_errors": "-",
        "moradi_errors": "--",
    }

    for error_mode in ["published_errors", "moradi_errors"]:
        for fit_label in ["Fit 5 form", "Fit 8 form"]:
            part = table.loc[
                (table["error_mode"] == error_mode)
                & (table["fit_form"] == fit_label)
            ].sort_values("BH_cut_percent")

            label = (
                f"{fit_label}, "
                + ("published errors" if error_mode == "published_errors"
                   else "Moradi + threshold sys.")
            )
            x = part["BH_cut_percent"].to_numpy(float)
            ls = line_styles[error_mode]

            axes[0, 0].errorbar(
                x, part["rE_fm"], yerr=part["rE_err_fm"],
                marker="o", markersize=3.5, capsize=2,
                linestyle=ls, label=label,
            )
            axes[0, 1].errorbar(
                x, part["rM_fm"], yerr=part["rM_err_fm"],
                marker="o", markersize=3.5, capsize=2,
                linestyle=ls, label=label,
            )
            axes[1, 0].plot(
                x, part["chi2_ndof"],
                marker="o", markersize=3.5, linestyle=ls, label=label,
            )
            axes[1, 1].plot(
                x, part["N_total"],
                marker="o", markersize=3.5, linestyle=ls, label=label,
            )
        #endfor
    #endfor

    axes[0, 0].set_ylabel(r"$r_E$ (fm)")
    axes[0, 1].set_ylabel(r"$r_M$ (fm)")
    axes[1, 0].set_ylabel(r"$\chi^2/{\rm ndf}$")
    axes[1, 1].set_ylabel("selected points")
    axes[1, 0].set_xlabel(r"$|1-R_{\rm BH}|$ threshold (%)")
    axes[1, 1].set_xlabel(r"$|1-R_{\rm BH}|$ threshold (%)")

    for ax in axes.ravel():
        ax.grid(alpha=0.2)
    #endfor
    axes[0, 0].legend(fontsize=7)
    fig.suptitle(
        "BH-purity plateau: published uncertainties vs Moradi prescription",
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(outdir / f"BH_cut_plateau_{tag}.png", dpi=180)
    plt.close(fig)
#enddef



def run_clas6_pass1_combined_from_bundles(
        clas6_bundle: Dict[str, object],
        pass1_bundle: Dict[str, object],
        args) -> Dict[str, object]:
    """Fit common FF parameters to the 5% CLAS6 and pass-1 samples."""
    outdir = (
        Path(args.outdir).expanduser().resolve()
        / "clas6_pass1_combined"
    )
    outdir.mkdir(parents=True, exist_ok=True)

    clas6 = clas6_bundle["set5"]
    pass1 = pass1_bundle["set5"]

    results = []
    for kind, fit_name in [
        ("dipole", "Fit 5"),
        ("fit6_same_a", "Fit 6"),
        ("fit7_same_p", "Fit 7"),
        ("fit8_f2_kelly", "Fit 8"),
    ]:
        fr = fit_combined_clas6_pass1(
            clas6=clas6,
            pass1=pass1,
            kind=kind,
            fit_name=fit_name,
            bh_cut=0.05,
        )
        results.append(fr)
        print(
            f"[COMBINED] {fit_name}: N={fr.npts} "
            f"chi2/dof={fr.chi2_ndof:.3f} "
            f"rE={fr.rE_fm:.5f}+/-{fr.rE_err_fm:.5f} fm "
            f"rM={fr.rM_fm:.5f}+/-{fr.rM_err_fm:.5f} fm "
            f"beta_pass1={fr.meta['beta_pass1']:+.3f}"
        )
    #endfor

    pd.DataFrame(
        [fitresult_to_record(fr) for fr in results]
    ).to_csv(outdir / "fit_results.csv", index=False)

    # The combined sample's t-range markers use the union.
    reference = pd.concat([clas6, pass1], ignore_index=True)
    save_fit5_to_fit8_sachs_plot(results, reference, outdir)
    save_combination_f1_f2_fits_5_to_8(results, reference, outdir)
    save_radii_plot(results, outdir)
    save_low_q2_ratio_plots(results, reference, outdir)
    save_elastic_reference_comparison(results, reference, outdir, "Fit 5")
    save_elastic_reference_comparison(results, reference, outdir, "Fit 8")
    save_bh_local_f1_f2_sensitivity(
        outdir, reference, next(r for r in results if r.name == "Fit 5")
    )

    clas6_all = clas6_bundle.get("all_data", clas6)
    pass1_all = pass1_bundle.get("all_data", pass1)
    save_bh_cut_plateau_study_multi([clas6_bundle, pass1_bundle], outdir, tag='jo2015_plus_pass1')

    return {
        "label": "CLAS6 + CLAS12 Lee 2026",
        "results": results,
        "set5": reference,
        "outdir": outdir,
    }
#enddef


def save_published_dataset_comparisons(
        bundles: Sequence[Dict[str, object]],
        outdir: Path) -> None:
    """CLAS6, CLAS12 Lee 2026, and their simultaneous-fit comparison."""
    outdir.mkdir(parents=True, exist_ok=True)
    for fit_name in ["Fit 5", "Fit 8"]:
        save_cross_dataset_form_factor_comparison(
            bundles, outdir, fit_name
        )
    #endfor
    save_cross_dataset_radius_comparison(bundles, outdir)
    save_cross_dataset_summary(bundles, outdir)
#enddef


def bundle_to_measurement_spec(
        bundle: Dict[str, object],
        data: Optional[pd.DataFrame] = None) -> Dict[str, object]:
    """Convert a result bundle into the generic multi-measurement fit spec."""
    selected = bundle["set5"] if data is None else data
    spec = measurement_spec(
        key=str(bundle["key"]),
        label=str(bundle["label"]),
        kind=str(bundle["kind"]),
        data=selected,
        norm_frac=float(bundle.get("norm_frac", 0.0)),
    )
    # A dataset may deliberately carry an unconstrained overall cross-section
    # scale.  This is distinct from setting an arbitrarily large Gaussian prior:
    # the nuisance is fitted directly as a fractional scale and receives no
    # beta^2 penalty.  The all-four status baseline uses this for Georges.
    spec["unconstrained_norm"] = bool(bundle.get("unconstrained_norm", False))
    return spec
#enddef


def audit_measurement_bundle_for_combination(
        bundle: Dict[str, object]) -> None:
    """Fail early with a complete schema/error audit before combined fits."""
    label = str(bundle["label"])
    kind = str(bundle["kind"])
    d = bundle["set5"]

    required = {"t_abs", "xs", "bh_A", "bh_B", "bh_C"}
    missing = sorted(required.difference(d.columns))
    if missing:
        raise KeyError(
            f"{label} is missing required 5% combined-fit columns: "
            + ", ".join(missing)
        )
    #endif

    stat = dataset_statistical_errors(d, kind)
    total = dataset_point_errors(d, kind, 0.05, True)
    for name, arr in [("statistical", stat), ("combined-fit", total)]:
        arr = np.asarray(arr, dtype=float)
        if (
            len(arr) != len(d)
            or not np.all(np.isfinite(arr))
            or np.any(arr <= 0.0)
        ):
            raise ValueError(
                f"{label} has invalid {name} errors: "
                f"Ndata={len(d)}, Nerr={len(arr)}, "
                f"Nnonfinite={int(np.sum(~np.isfinite(arr)))}, "
                f"Nnonpositive={int(np.sum(arr <= 0.0))}."
            )
        #endif
    #endfor

    print(
        f"[COMBINATION preflight] {label}: N={len(d)}, "
        f"kind={kind}, error schema OK"
    )
#enddef



def save_combination_f1_f2_fits_5_to_8(
        results: Sequence[FitResult],
        reference: pd.DataFrame,
        outdir: Path) -> None:
    """Explicit F1/F2 curves for Fits 5--8 for every combined fit."""
    chosen = [fr for fr in results if fr.name in {"Fit 5", "Fit 6", "Fit 7", "Fit 8"}]
    if not chosen:
        return
    #endif

    qmax = float(np.nanmax(reference["t_abs"]))
    q = np.linspace(0.0, max(qmax, 1e-3), 350)

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.0))
    for fr in chosen:
        f1, f2 = paper_model_f1f2(fr.model_kind, q, fr.params)
        axes[0].plot(q, f1, label=fr.name)
        axes[1].plot(q, f2, label=fr.name)
    #endfor

    kf1, kf2 = kelly_f1_f2(q)
    axes[0].plot(q, kf1, linestyle="--", label="Kelly")
    axes[1].plot(q, kf2, linestyle="--", label="Kelly")

    axes[0].set_ylabel(r"$F_1$")
    axes[1].set_ylabel(r"$F_2$")
    tlo = float(np.nanmin(reference["t_abs"]))
    thi = float(np.nanmax(reference["t_abs"]))
    for ax in axes:
        ax.set_xlabel(r"$|t|$ (GeV$^2$)")
        ax.grid(alpha=0.2)
        ax.axvline(tlo, linewidth=0.8, linestyle=":")
        ax.axvline(thi, linewidth=0.8, linestyle=":")
    #endfor

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=5, bbox_to_anchor=(0.5, 0.965))
    fig.suptitle("Combined BH fit: Dirac and Pauli form factors", y=0.995)
    fig.subplots_adjust(top=0.84, bottom=0.14, left=0.08, right=0.98, wspace=0.24)
    fig.savefig(outdir / "07_F1_F2_fits_5_to_8.png", dpi=180)
    plt.close(fig)
#enddef


def save_fit5_residual_chi2_diagnostics(
        data: pd.DataFrame,
        fit: FitResult,
        dataset_kind: str,
        dataset_key: str,
        outdir: Path,
        bh_cut: float = 0.05) -> None:
    """Residual/pull and per-point chi2 decomposition versus kinematics."""
    outdir.mkdir(parents=True, exist_ok=True)

    q = data["t_abs"].to_numpy(float)
    f1, f2 = paper_model_f1f2(fit.model_kind, q, fit.params)
    pred = bh_from_f1f2(
        data["bh_A"].to_numpy(float),
        data["bh_B"].to_numpy(float),
        data["bh_C"].to_numpy(float),
        f1, f2,
    )

    beta_key = "beta_" + re.sub(r"[^A-Za-z0-9]+", "_", dataset_key)
    pred *= float(fit.meta.get(beta_key + "_scale_factor", 1.0))

    err = dataset_point_errors(data, dataset_kind, bh_cut, True)
    y = data["xs"].to_numpy(float)
    residual = y - pred
    pull = residual / err
    chi2_i = pull**2

    diag = data.copy()
    diag["fit5_prediction"] = pred
    diag["fit5_residual"] = residual
    diag["fit5_pull"] = pull
    diag["fit5_chi2_contribution"] = chi2_i
    diag.to_csv(outdir / "fit5_residual_chi2_per_point.csv", index=False)

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.0))
    variables = [
        ("phi_deg", r"$\phi$ (deg)"),
        ("t_abs", r"$|t|$ (GeV$^2$)"),
        ("xB", r"$x_B$"),
        ("Q2", r"$Q^2$ (GeV$^2$)"),
    ]

    binned_rows = []
    for ax, (col, xlabel) in zip(axes.flat, variables):
        x = diag[col].to_numpy(float)
        ax.scatter(x, pull, s=16, facecolors="none", alpha=0.45)
        ax.axhline(0.0, linewidth=0.8)
        ax.axhline(+2.0, linewidth=0.7, linestyle=":")
        ax.axhline(-2.0, linewidth=0.7, linestyle=":")
        finite = np.isfinite(x) & np.isfinite(pull)

        if np.sum(finite) >= 10:
            edges = np.linspace(np.nanmin(x[finite]), np.nanmax(x[finite]), 11)
            centers = 0.5 * (edges[:-1] + edges[1:])
            means = []
            for ibin, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
                m = finite & (x >= lo) & (
                    (x < hi) if ibin < len(edges) - 2 else (x <= hi)
                )
                means.append(float(np.nanmean(pull[m])) if np.any(m) else np.nan)
                if np.any(m):
                    binned_rows.append({
                        "variable": col,
                        "bin": ibin,
                        "low": lo,
                        "high": hi,
                        "N": int(np.sum(m)),
                        "chi2_sum": float(np.sum(chi2_i[m])),
                        "chi2_per_point": float(np.mean(chi2_i[m])),
                        "mean_pull": float(np.mean(pull[m])),
                        "rms_pull": float(np.sqrt(np.mean(pull[m]**2))),
                    })
                #endif
            #endfor
            ax.plot(centers, means, marker="o", linewidth=1.0, label="binned mean pull")
        #endif
        ax.set_xlabel(xlabel)
        ax.set_ylabel("pull")
        ax.grid(alpha=0.2)
    #endfor

    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.965))
    #endif
    fig.suptitle(
        f"{dataset_key}: Fit-5 residual structure; point chi2={np.sum(chi2_i):.1f}",
        y=0.995,
    )
    fig.subplots_adjust(top=0.92, bottom=0.08, left=0.08, right=0.98, hspace=0.28, wspace=0.24)
    fig.savefig(outdir / "fit5_pulls_vs_kinematics.png", dpi=180)
    plt.close(fig)

    pd.DataFrame(binned_rows).to_csv(
        outdir / "fit5_chi2_binned_vs_kinematics.csv", index=False
    )
#enddef


def run_measurement_combination(
        bundles: Sequence[Dict[str, object]],
        args,
        tag: str) -> Dict[str, object]:
    """
    Fit one pair/triple measurement combination with common FF parameters.

    Fits 5--8 use the common 5% BH-selected sample for each measurement.
    Each experiment retains its own pointwise errors and correlated
    normalization nuisance.
    """
    outdir = (
        Path(args.outdir).expanduser().resolve()
        / "measurement_combinations"
        / tag
    )
    outdir.mkdir(parents=True, exist_ok=True)

    specs = [bundle_to_measurement_spec(bundle) for bundle in bundles]
    labels = [str(bundle["label"]) for bundle in bundles]

    results = []
    for kind, fit_name in [
        ("dipole", "Fit 5"),
        ("fit6_same_a", "Fit 6"),
        ("fit7_same_p", "Fit 7"),
        ("fit8_f2_kelly", "Fit 8"),
    ]:
        fr = fit_multi_measurements(
            datasets=specs,
            kind=kind,
            fit_name=fit_name,
            bh_cut=0.05,
            add_moradi_bh_systematic=True,
        )
        results.append(fr)

        nuisance_bits = []
        for spec in specs:
            key = re.sub(r"[^A-Za-z0-9]+", "_", str(spec["key"]))
            beta_key = f"beta_{key}"
            scale_key = beta_key + "_scale_factor"
            if beta_key in fr.meta:
                nuisance_bits.append(f"{beta_key}={fr.meta[beta_key]:+.3f}")
            #endif
            if scale_key in fr.meta:
                nuisance_bits.append(f"{scale_key}={fr.meta[scale_key]:.3f}")
            #endif
        #endfor

        suffix = ("; " + ", ".join(nuisance_bits)) if nuisance_bits else ""
        print(
            f"[COMBINATION] {' + '.join(labels)} | {fit_name}: "
            f"N={fr.npts} chi2/dof={fr.chi2_ndof:.3f} "
            f"rE={fr.rE_fm:.5f}+/-{fr.rE_err_fm:.5f} fm "
            f"rM={fr.rM_fm:.5f}+/-{fr.rM_err_fm:.5f} fm"
            f"{suffix}"
        )
    #endfor

    pd.DataFrame(
        [fitresult_to_record(fr) for fr in results]
    ).to_csv(outdir / "fit_results.csv", index=False)

    reference = pd.concat(
        [bundle["set5"] for bundle in bundles],
        ignore_index=True,
    )

    save_fit5_to_fit8_sachs_plot(results, reference, outdir)
    save_radii_plot(results, outdir)
    save_low_q2_ratio_plots(results, reference, outdir)
    save_elastic_reference_comparison(
        results, reference, outdir, "Fit 5"
    )
    save_elastic_reference_comparison(
        results, reference, outdir, "Fit 8"
    )

    return {
        "label": " + ".join(labels),
        "key": tag,
        "kind": "combined",
        "norm_frac": 0.0,
        "results": results,
        "set5": reference,
        "outdir": outdir,
    }
#enddef


def save_separate_vs_combined_radius_summary(
        bundles: Sequence[Dict[str, object]],
        outdir: Path) -> None:
    """Compare Fit-5 and Fit-8 radii for all separate/combined cases."""
    outdir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(
        2, 2, figsize=(13.5, 8.5),
        sharex=False,
    )

    for irow, fit_name in enumerate(["Fit 5", "Fit 8"]):
        labels = []
        re_vals = []
        re_errs = []
        rm_vals = []
        rm_errs = []

        for bundle in bundles:
            matches = [
                fr for fr in bundle["results"]
                if fr.name == fit_name
            ]
            if not matches:
                continue
            #endif

            fr = matches[0]
            labels.append(str(bundle["label"]))
            re_vals.append(fr.rE_fm)
            re_errs.append(fr.rE_err_fm)
            rm_vals.append(fr.rM_fm)
            rm_errs.append(fr.rM_err_fm)
        #endfor

        y = np.arange(len(labels))
        axes[irow, 0].errorbar(
            re_vals, y, xerr=re_errs,
            fmt="o", fillstyle="none", capsize=2,
        )
        axes[irow, 1].errorbar(
            rm_vals, y, xerr=rm_errs,
            fmt="o", fillstyle="none", capsize=2,
        )

        for icol in range(2):
            axes[irow, icol].set_yticks(y)
            axes[irow, icol].set_yticklabels(labels, fontsize=8)
            axes[irow, icol].grid(alpha=0.2)
        #endfor

        axes[irow, 0].set_ylabel(fit_name)
    #endfor

    axes[0, 0].set_title(r"$r_E$")
    axes[0, 1].set_title(r"$r_M$")
    axes[1, 0].set_xlabel("radius (fm)")
    axes[1, 1].set_xlabel("radius (fm)")
    fig.suptitle(
        "BH radius extraction: separate measurements and combinations",
        y=0.995,
    )
    fig.subplots_adjust(
        top=0.93,
        bottom=0.08,
        left=0.28,
        right=0.98,
        hspace=0.28,
        wspace=0.28,
    )
    fig.savefig(
        outdir / "separate_vs_combined_radii.png",
        dpi=180,
    )
    plt.close(fig)
#enddef


def save_separate_vs_combined_table(
        bundles: Sequence[Dict[str, object]],
        outdir: Path) -> None:
    """CSV summary of Fit-5/Fit-8 separate and combined extractions."""
    rows = []
    for bundle in bundles:
        for fr in bundle["results"]:
            if fr.name not in {"Fit 5", "Fit 8"}:
                continue
            #endif
            row = {
                "dataset_or_combination": bundle["label"],
                "fit": fr.name,
                "N": fr.npts,
                "chi2_ndof": fr.chi2_ndof,
                "rE_fm": fr.rE_fm,
                "rE_err_fm": fr.rE_err_fm,
                "rM_fm": fr.rM_fm,
                "rM_err_fm": fr.rM_err_fm,
            }
            for key, value in fr.meta.items():
                if (
                    key.startswith("beta_")
                    or key.endswith("_scale_factor")
                ):
                    row[key] = value
                #endif
            #endfor
            rows.append(row)
        #endfor
    #endfor

    pd.DataFrame(rows).to_csv(
        outdir / "separate_vs_combined_summary.csv",
        index=False,
    )
#enddef



def save_published_workflow_manifest(
        bundles: Sequence[Dict[str, object]],
        combinations: Sequence[Dict[str, object]],
        outdir: Path) -> None:
    """Write a simple machine- and human-readable completion manifest."""
    rows = []
    for bundle in bundles:
        rows.append({
            "category": "single",
            "key": str(bundle["key"]),
            "label": str(bundle["label"]),
            "output_directory": str(Path(bundle["outdir"]).resolve()),
            "set5_points": int(len(bundle["set5"])),
        })
    #endfor

    for bundle in combinations:
        rows.append({
            "category": "combination",
            "key": str(bundle["key"]),
            "label": str(bundle["label"]),
            "output_directory": str(Path(bundle["outdir"]).resolve()),
            "set5_points": int(len(bundle["set5"])),
        })
    #endfor

    manifest = pd.DataFrame(rows)
    manifest.to_csv(outdir / "published_workflow_manifest.csv", index=False)

    with open(outdir / "published_workflow_manifest.txt", "w") as fout:
        fout.write("Published BH form-factor/radius workflow products\n")
        fout.write("=" * 72 + "\n")
        for row in rows:
            fout.write(
                f"{row['category']:11s} | {row['label']} | "
                f"N5={row['set5_points']} | {row['output_directory']}\n"
            )
        #endfor
    #endwith
#enddef




BH_MODEL_PHI_SOURCE = {
    "jo2015": (
        "Gepard dataset 98 (BMK-stored)",
        "180-minus",
    ),
    "halla_defurne2015": (
        "Gepard datasets 107/108/112/113 (BMK-stored)",
        "180-minus",
    ),
    "halla_defurne2017": (
        "Gepard datasets 129/130/131/132/133/134 (BMK-stored)",
        "180-minus",
    ),
    "saylor2018": (
        "published supplemental table (direct)",
        "identity",
    ),
    "halla_georges2022": (
        "published E12-06-114 spreadsheet (Trento)",
        "identity",
    ),
    "pass1": (
        "CLAS12 Lee table",
        "identity",
    ),
    "pass2": (
        "CLAS12 Hayward table",
        "identity",
    ),
}


def export_bh_model_selection_kinematics(
        bundles: Sequence[Dict[str, object]],
        root_outdir: Path) -> Path:
    """
    Export every exact experimental point once for external VGG/GK evaluation.

    The external model stage must not rerun KM15.  This table therefore carries
    the KM15 EP/BH decomposition and the exact kinematics already used by the
    production analysis.  Later model-selection results can be merged back by
    the deterministic point_id.
    """
    rows = []
    for bundle in bundles:
        data = bundle.get("all_data", bundle["set5"]).copy()
        key = str(bundle["key"])

        if key not in BH_MODEL_PHI_SOURCE:
            raise RuntimeError(
                f"{bundle['label']}: no explicit phi-source convention has "
                "been configured for external PARTONS evaluation"
            )
        #endif

        required = ["xB", "Q2", "t_abs", "phi_deg", "ebeam"]
        missing = [col for col in required if col not in data.columns]
        if missing:
            raise KeyError(
                f"{bundle['label']}: cannot export BH-model kinematics; "
                f"missing {missing}"
            )
        #endif

        for local_row, (_, row) in enumerate(data.iterrows()):
            source_repr, expected_mapping = BH_MODEL_PHI_SOURCE.get(
                key, ("UNCONFIGURED", "UNCONFIGURED")
            )
            record = {
                "point_id": f"{key}:{local_row}",
                "dataset": key,
                "source_row": local_row,
                "xB": float(row["xB"]),
                "Q2": float(row["Q2"]),
                "t_abs": float(row["t_abs"]),
                "phi_deg": float(row["phi_deg"]),
                "ebeam": float(row["ebeam"]),
                "phi_source_representation": source_repr,
                "phi_partons_expected_mapping": expected_mapping,
            }
            for col in [
                "km15_ep", "km15_bh", "km15_dvcs", "km15_int",
                "rbh", "delta_bh", "bh_A", "bh_B", "bh_C",
            ]:
                if col in data.columns:
                    record[col] = row[col]
                #endif
            #endfor
            rows.append(record)
        #endfor
    #endfor

    outdir = root_outdir / "bh_model_selection"
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / "bh_model_kinematics.csv"
    table = pd.DataFrame(rows)
    if table["point_id"].duplicated().any():
        raise RuntimeError("Duplicate point_id while exporting BH-model kinematics")
    #endif
    table.to_csv(path, index=False)
    print(
        f"[BH-model selection] exported {len(table)} exact points from "
        f"{len(bundles)} measurement(s) -> {path}"
    )
    return path
#enddef




FINAL_MODEL_SELECTION_DEFAULT = (
    "/work/clas12/thayward/CLAS12_exclusive/dvcs/model_predictions/"
    "bh_model_selection_all_points.csv"
)
FINAL_MODEL_NAMES = ("km15", "vgg99", "gk16")
FINAL_NOMINAL_MODEL = "km15"
FINAL_PHYSICS_DATASETS = ("jo2015", "pass1")

# Publication-facing chronological ordering.  Unknown/future datasets are
# placed after the known sequence while preserving their input order.
DATASET_CHRONOLOGY = {
    "jo2015": 0,
    "halla_defurne2015": 1,
    "halla_defurne2017": 2,
    "saylor2018": 3,
    "halla_georges2022": 4,
    "pass1": 5,
    "pass2": 6,
}


def sort_bundles_chronologically(
        bundles: Sequence[Dict[str, object]]) -> List[Dict[str, object]]:
    """Return measurement bundles in the standard publication chronology."""
    indexed = list(enumerate(bundles))
    indexed.sort(
        key=lambda item: (
            DATASET_CHRONOLOGY.get(str(item[1].get("key", "")), 10_000),
            item[0],
        )
    )
    return [bundle for _, bundle in indexed]
#enddef


def load_external_bh_model_selection(path: Path) -> pd.DataFrame:
    """Load and validate the completed KM15/VGG99/GK16 model-selection table."""
    if not path.exists():
        raise FileNotFoundError(
            "Completed external BH-model selection table not found:\n"
            f"  {path}\n"
            "Run evaluate_bh_model_selection.py --run-partons first."
        )
    #endif

    table = pd.read_csv(path)
    required = {
        "point_id", "dataset", "source_row",
        "delta_bh_km15", "delta_bh_vgg99", "delta_bh_gk16",
    }
    missing = sorted(required.difference(table.columns))
    if missing:
        raise KeyError(
            "External BH-model selection table is missing required columns: "
            + ", ".join(missing)
        )
    #endif

    if table["point_id"].duplicated().any():
        raise RuntimeError("Duplicate point_id values in BH-model selection table")
    #endif

    for col in ["source_row", "delta_bh_km15", "delta_bh_vgg99", "delta_bh_gk16"]:
        vals = pd.to_numeric(table[col], errors="coerce")
        if not np.all(np.isfinite(vals.to_numpy(float))):
            raise RuntimeError(f"Nonfinite values in external selection column {col}")
        #endif
        table[col] = vals
    #endfor

    table["source_row"] = table["source_row"].astype(int)
    print(
        f"[final model analysis] loaded {len(table)} externally evaluated points "
        f"from {path}"
    )
    return table
#enddef


def select_bundle_from_external_model(
        bundle: Dict[str, object],
        selection_table: pd.DataFrame,
        model: str,
        threshold: float) -> pd.DataFrame:
    """
    Reconstruct one model-specific selected sample from deterministic source_row.

    The external PARTONS stage never carries the measured cross section/error
    columns.  Those remain authoritative in bundle["all_data"]; source_row maps
    the external purity decision back onto that exact table.
    """
    if model not in FINAL_MODEL_NAMES:
        raise ValueError(f"Unknown BH-purity model {model}")
    #endif

    key = str(bundle["key"])
    all_data = bundle["all_data"].reset_index(drop=True)

    # Saylor is evaluated directly with the same Gepard/KM15 decomposition
    # already constructed by run_saylor_validation().  It is not present in
    # the external five-dataset PARTONS selection table.  A bundle may also
    # carry a fixed lower-|t| requirement, used for the nominal surviving
    # Saylor sample motivated independently by arXiv:2607.04481.
    if bool(bundle.get("direct_km15_selection", False)):
        if str(model) != "km15":
            raise ValueError(
                f"{bundle['label']}: direct selection is KM15-only"
            )
        #endif
        delta_col = str(bundle.get("direct_km15_delta_column", "bh_delta"))
        if delta_col not in all_data.columns and delta_col == "bh_delta":
            delta_col = "delta_bh"
        #endif
        if delta_col not in all_data.columns:
            raise KeyError(
                f"{bundle['label']}: no direct KM15 delta column "
                f"'{delta_col}'"
            )
        #endif
        delta = pd.to_numeric(
            all_data[delta_col], errors="coerce"
        ).to_numpy(float)
        mask = np.isfinite(delta) & (delta <= float(threshold))
        tmin = bundle.get("fixed_t_abs_min_GeV2", None)
        if tmin is not None:
            mask &= (
                pd.to_numeric(
                    all_data["t_abs"], errors="coerce"
                ).to_numpy(float) >= float(tmin) - 1.0e-12
            )
        #endif
        selected = all_data.loc[mask].copy()
        selected["external_delta_bh"] = delta[mask]
        selected["external_bh_model"] = "km15_direct"
        selected["external_bh_threshold"] = float(threshold)
        if tmin is not None:
            selected["fixed_t_abs_min_GeV2"] = float(tmin)
        #endif
        return selected
    #endif

    ext = selection_table.loc[selection_table["dataset"].astype(str) == key].copy()

    if len(ext) != len(all_data):
        raise RuntimeError(
            f"{bundle['label']}: external selection has {len(ext)} points but "
            f"the current bundle has {len(all_data)}. Regenerate "
            "bh_model_kinematics.csv and the PARTONS cache before fitting."
        )
    #endif

    ext = ext.sort_values("source_row")
    expected = np.arange(len(all_data), dtype=int)
    got = ext["source_row"].to_numpy(int)
    if not np.array_equal(got, expected):
        raise RuntimeError(
            f"{bundle['label']}: source_row mapping is not exactly 0..N-1"
        )
    #endif

    # Strong kinematic audit when the external table carries these columns.
    for col in ["xB", "Q2", "t_abs", "phi_deg", "ebeam"]:
        if col in ext.columns and col in all_data.columns:
            a = all_data[col].to_numpy(float)
            b = ext[col].to_numpy(float)
            scale = np.maximum(1.0, np.maximum(np.abs(a), np.abs(b)))
            if np.any(np.abs(a - b) > 1.0e-10 * scale):
                raise RuntimeError(
                    f"{bundle['label']}: external/current kinematic mismatch in {col}"
                )
            #endif
        #endif
    #endfor

    mask = ext[f"delta_bh_{model}"].to_numpy(float) <= float(threshold)
    selected = all_data.loc[mask].copy()
    selected["external_delta_bh"] = ext.loc[mask, f"delta_bh_{model}"].to_numpy(float)
    selected["external_bh_model"] = model
    selected["external_bh_threshold"] = float(threshold)
    return selected
#enddef


def fit_sachs_family_multi_measurements(
        datasets: Sequence[Dict[str, object]],
        family: str,
        bh_cut: float = 0.05,
        add_moradi_bh_systematic: bool = True,
        bh_systematic_fraction: Optional[float] = None,
        unconstrained_norm_keys: Optional[Sequence[str]] = None
        ) -> Dict[str, object]:
    """
    Production fit of one Sachs family to multiple BH-selected measurements.

    This is the production counterpart of fit_cross_sections_with_sachs_family:
    GE and GM/mu_p are parameterized directly, converted point-by-point to
    F1/F2, and evaluated with each point's exact BH quadratic.  Unlike the
    replica fitter, this version retains each experiment's correlated
    normalization nuisance and returns Hessian uncertainties.
    """
    family_e, family_m = decode_sachs_family_pair(family)
    ne = int(re.findall(r"\d+", family_e)[0])
    nm = int(re.findall(r"\d+", family_m)[0])
    names_e = [f"e{i}" for i in range(1, ne + 1)]
    names_m = [f"m{i}" for i in range(1, nm + 1)]
    shape_names = names_e + names_m

    p0 = np.zeros(ne + nm, dtype=float)
    p0[0] = sachs_first_coefficient_from_radius(SACHS_INITIAL_RADIUS_FM, family_e)
    p0[ne] = sachs_first_coefficient_from_radius(SACHS_INITIAL_RADIUS_FM, family_m)

    unconstrained_norm_keys = set(
        str(k) for k in (unconstrained_norm_keys or [])
    )
    unconstrained_norm_keys.update(
        str(spec["key"]) for spec in datasets
        if bool(spec.get("unconstrained_norm", False))
    )
    nuisance_names = []
    nuisance_fracs = {}
    nuisance_is_free = {}
    for spec in datasets:
        key = str(spec["key"])
        free_norm = key in unconstrained_norm_keys
        frac = 1.0 if free_norm else float(spec.get("norm_frac", 0.0))
        if frac > 0.0:
            nname = "beta_" + re.sub(r"[^A-Za-z0-9]+", "_", key)
            nuisance_names.append(nname)
            nuisance_fracs[nname] = frac
            nuisance_is_free[nname] = bool(free_norm)
        #endif
    #endfor

    fit_names = shape_names + nuisance_names
    fit_p0 = np.concatenate([p0, np.zeros(len(nuisance_names), dtype=float)])
    nuisance_index = {
        name: len(shape_names) + i for i, name in enumerate(nuisance_names)
    }

    prepared = []
    for spec in datasets:
        d = spec["data"]
        if len(d) == 0:
            print(
                f"[final fit] {spec['label']}: 0 selected points; "
                "omitting this measurement from this fit"
            )
            continue
        #endif
        error_bh_fraction = (
            float(bh_cut) if bh_systematic_fraction is None
            else float(bh_systematic_fraction)
        )
        err = dataset_point_errors(
            d, str(spec["kind"]), error_bh_fraction, add_moradi_bh_systematic
        )
        if (
            len(err) != len(d)
            or not np.all(np.isfinite(err))
            or np.any(err <= 0.0)
        ):
            raise RuntimeError(f"{spec['label']}: invalid production point errors")
        #endif
        q = d["t_abs"].to_numpy(float)
        tau = q / (4.0 * MP2)
        key = str(spec["key"])
        frac = 1.0 if key in unconstrained_norm_keys else float(spec.get("norm_frac", 0.0))
        prepared.append({
            "key": key, "norm_frac": frac,
            "nuisance_name": ("beta_" + re.sub(r"[^A-Za-z0-9]+", "_", key) if frac > 0.0 else None),
            "q": q, "q_powers": np.vstack([q**i for i in range(1, max(ne, nm) + 1)]),
            "tau": tau, "inv_one_plus_tau": 1.0 / (1.0 + tau),
            "y": d["xs"].to_numpy(float), "e": np.asarray(err, dtype=float),
            "A": d["bh_A"].to_numpy(float), "B": d["bh_B"].to_numpy(float),
            "C": d["bh_C"].to_numpy(float), "N": len(d),
        })
    #endfor

    if not prepared:
        raise RuntimeError(
            "Production Sachs fit has zero selected points across all measurements"
        )
    #endif

    def chi2_minuit(*values):
        p = np.asarray(values, dtype=float)
        ce = p[:ne]
        cm = p[ne:ne + nm]
        total = 0.0

        for item in prepared:
            q = item["q"]
            ge = sachs_family_value_precomputed(q, ce, family_e, item["q_powers"])
            gm = MU_P * sachs_family_value_precomputed(q, cm, family_m, item["q_powers"])
            tau = item["tau"]; inv = item["inv_one_plus_tau"]
            f1 = (ge + tau * gm) * inv
            f2 = (gm - ge) * inv
            pred = bh_from_f1f2(item["A"], item["B"], item["C"], f1, f2)

            frac = item["norm_frac"]
            if frac > 0.0:
                beta = p[nuisance_index[item["nuisance_name"]]]
                pred = pred * (1.0 + frac * beta)
            #endif

            pull = (pred - item["y"]) / item["e"]
            total += float(np.dot(pull, pull))
        #endfor

        for nname in nuisance_names:
            if not nuisance_is_free.get(nname, False):
                total += float(p[nuisance_index[nname]] ** 2)
            #endif
        #endfor
        return total
    #enddef

    m = Minuit(chi2_minuit, *fit_p0, name=tuple(fit_names))
    m.errordef = Minuit.LEAST_SQUARES

    # Broad radius-equivalent slope guards suppress optimizer excursions to
    # zero-radius or multi-fm solutions without choosing the physical answer.
    ce_lo = sachs_first_coefficient_from_radius(SACHS_MIN_RADIUS_FM, family_e)
    ce_hi = sachs_first_coefficient_from_radius(SACHS_MAX_RADIUS_FM, family_e)
    cm_lo = sachs_first_coefficient_from_radius(SACHS_MIN_RADIUS_FM, family_m)
    cm_hi = sachs_first_coefficient_from_radius(SACHS_MAX_RADIUS_FM, family_m)
    m.limits[names_e[0]] = (min(ce_lo, ce_hi), max(ce_lo, ce_hi))
    m.limits[names_m[0]] = (min(cm_lo, cm_hi), max(cm_lo, cm_hi))
    for nname in nuisance_names:
        if nuisance_is_free.get(nname, False):
            # beta is directly the fractional scale shift when frac=1.
            # Keep the diagnostic physically positive and broad.
            m.limits[nname] = (-0.50, 0.50)
        else:
            m.limits[nname] = (-10.0, 10.0)
        #endif
    #endfor

    m.migrad()
    m.hesse()

    full = np.array([float(m.values[n]) for n in fit_names], dtype=float)
    shape = full[:ne + nm]
    cov = np.full((ne + nm, ne + nm), np.nan)
    if m.covariance is not None:
        for i, ni in enumerate(shape_names):
            for j, nj in enumerate(shape_names):
                cov[i, j] = float(m.covariance[ni, nj])
            #endfor
        #endfor
    #endif

    def radius_e(p):
        return sachs_family_radius(np.asarray(p[:ne]), family_e)
    #enddef

    def radius_m(p):
        return sachs_family_radius(np.asarray(p[ne:ne + nm]), family_m)
    #enddef

    rE, rE_err = propagate_scalar(radius_e, shape, cov)
    rM, rM_err = propagate_scalar(radius_m, shape, cov)

    ndata = int(sum(item["N"] for item in prepared))
    nfree = len(fit_names)
    ndof = max(1, ndata - nfree)
    result = {
        "family": family,
        "N": ndata,
        "n_parameters": nfree,
        "chi2": float(m.fval),
        "ndof": int(ndof),
        "chi2_ndof": float(m.fval / ndof),
        "valid": bool(m.valid and np.isfinite(rE) and np.isfinite(rM)),
        "accurate": bool(m.accurate),
        "rE_fm": float(rE),
        "rE_fit_err_fm": float(rE_err),
        "rM_fm": float(rM),
        "rM_fit_err_fm": float(rM_err),
    }

    for name in shape_names:
        result[name] = float(m.values[name])
        result[name + "_err"] = float(m.errors[name])
    #endfor
    # Preserve the complete shape covariance for final GE/GM uncertainty
    # bands.  JSON keeps the result CSV self-contained without creating a
    # separate covariance file for every fit.
    result["shape_covariance_json"] = json.dumps(cov.tolist())

    for nname in nuisance_names:
        result[nname] = float(m.values[nname])
        result[nname + "_err"] = float(m.errors[nname])
        result[nname + "_norm_fraction"] = float(nuisance_fracs[nname])
        result[nname + "_is_unconstrained"] = bool(
            nuisance_is_free.get(nname, False)
        )
        result[nname + "_scale_factor"] = float(
            1.0 + nuisance_fracs[nname] * float(m.values[nname])
        )
    #endfor
    return result
#enddef


def summarize_bias_for_family(
        study_csv: Path,
        family: str) -> Dict[str, float]:
    """
    Bias estimates for one fitted family using equal truth-group weighting.

    The production RMS first averages bias^2 within each truth_group and then
    gives every truth_group equal weight.  Raw scenario-weighted RMS values are
    retained in parallel for audit/backward comparison.
    """
    table = pd.read_csv(study_csv)
    table = production_closure_truth_rows(table)
    fam = table.loc[table["family"].astype(str) == family].copy()
    if len(fam) == 0:
        raise RuntimeError(f"No {family} rows found in {study_csv}")
    #endif

    out = {}
    for quantity in ["rE", "rM"]:
        qpart = fam.loc[fam["quantity"].astype(str) == quantity].copy()
        vals = qpart["bias_fm"].to_numpy(float)
        vals = vals[np.isfinite(vals)]
        if len(vals) == 0:
            raise RuntimeError(f"No finite {quantity} bias values for {family}")
        #endif

        # Signed mean is also made hierarchical so dense synthetic grids do not
        # dominate its interpretation.
        group_means = []
        for _, group in qpart.groupby("truth_group", sort=False):
            gv = group["bias_fm"].to_numpy(float)
            gv = gv[np.isfinite(gv)]
            if len(gv):
                group_means.append(float(np.mean(gv)))
            #endif
        #endfor

        out[f"{quantity}_mean_signed_bias_fm"] = (
            float(np.mean(group_means)) if group_means else np.nan
        )
        out[f"{quantity}_RMS_bias_fm"] = equal_truth_group_rms(
            qpart, "bias_fm"
        )
        out[f"{quantity}_scenario_weighted_RMS_bias_fm"] = float(
            np.sqrt(np.mean(vals**2))
        )
        out[f"{quantity}_max_abs_bias_fm"] = float(np.max(np.abs(vals)))
        out[f"{quantity}_truth_group_count"] = int(
            qpart["truth_group"].astype(str).nunique()
        )
    #endfor
    return out
#enddef


def aggregate_model_bias_rankings(
        bias_dirs: Dict[str, Path],
        outdir: Path) -> Tuple[str, pd.DataFrame]:
    """
    Select the production Sachs family from the KM15 closure study only.

    KM15 is the nominal/production BH-purity prescription. VGG99 and GK16 are
    external diagnostics only and therefore do not enter closure ranking,
    family vetoes, or the extrapolation-bias systematic.
    """
    if FINAL_NOMINAL_MODEL not in bias_dirs:
        raise RuntimeError("KM15 closure directory is required for production ranking")
    #endif

    path = (
        bias_dirs[FINAL_NOMINAL_MODEL]
        / "radius_bias_extended_eligibility_ranking.csv"
    )
    if not path.exists():
        raise FileNotFoundError(f"Missing KM15 extended bias ranking: {path}")
    #endif

    table = pd.read_csv(path).copy()
    required = [
        "family",
        "combined_RMS_objective_fm",
        "minimum_scenario_valid_fraction",
        "global_valid_fraction",
    ]
    missing = [c for c in required if c not in table.columns]
    if missing:
        raise RuntimeError(
            f"KM15 closure ranking is missing required columns: {missing}"
        )
    #endif

    ranking = table.copy()
    ranking = ranking.rename(columns={
        "combined_RMS_objective_fm": "km15_combined_RMS_objective_fm",
        "minimum_scenario_valid_fraction":
            "km15_minimum_scenario_valid_fraction",
        "global_valid_fraction": "km15_global_valid_fraction",
    })
    ranking["production_eligible_km15"] = np.isfinite(
        ranking["km15_combined_RMS_objective_fm"].to_numpy(float)
    )
    ranking = ranking.sort_values(
        ["production_eligible_km15", "km15_combined_RMS_objective_fm"],
        ascending=[False, True],
    ).reset_index(drop=True)

    eligible = ranking.loc[ranking["production_eligible_km15"]]
    if len(eligible) == 0:
        raise RuntimeError("No Sachs family has a finite KM15 closure objective.")
    #endif

    rank1 = str(eligible.iloc[0]["family"])
    outdir.mkdir(parents=True, exist_ok=True)
    ranking.to_csv(outdir / "km15_production_family_ranking.csv", index=False)

    with open(outdir / "bias_rank1_sachs_family.txt", "w") as fout:
        fout.write(f"bias_rank1_family={rank1}\n")
        fout.write(
            "criterion=finite KM15 closure objective; VGG99/GK16 are external "
            "diagnostics only; minimum KM15 combined RMS bias-variance objective\n"
        )
        fout.write(
            "note=This is the KM15 closure rank-1 family. Final production "
            "selection may fall through to the next-ranked family if the "
            "nominal KM15 5% production fit is invalid. Threshold-scan failures "
            "do not veto it.\n"
        )
    #endwith

    print(
        f"[final model analysis] KM15 closure rank-1 Sachs family is {rank1}; "
        "KM15 production-fit validation will determine the final family"
    )
    return rank1, ranking
#enddef



MODEL_BH_COLUMN = {
    "km15": "km15_bh",
    "vgg99": "partons_bh_vgg99",
    "gk16": "partons_bh_gv08",
}
MODEL_EP_COLUMN = {
    "km15": "km15_ep",
    "vgg99": "partons_ep_vgg99",
    "gk16": "partons_ep_gk16",
}
MODEL_DISPLAY = {
    "km15": "KM15 / Gepard",
    "vgg99": "VGG99 / PARTONS",
    "gk16": "GK16 / PARTONS",
}


def _external_table_with_measurements(
        selection: pd.DataFrame,
        physics_bundles: Sequence[Dict[str, object]]) -> pd.DataFrame:
    """
    Attach measured cross sections to the externally evaluated model table.

    This table is diagnostic only.  Measured cross sections are NEVER used to
    decide whether a model is accepted, to rank the purity models, or to alter
    the BH-purity selections.  That separation is deliberate to avoid circular
    logic in the extraction.
    """
    blocks = []
    for bundle in physics_bundles:
        key = str(bundle["key"])
        label = str(bundle["label"])
        data = bundle["all_data"].reset_index(drop=True)
        ext = selection.loc[
            selection["dataset"].astype(str) == key
        ].copy().sort_values("source_row").reset_index(drop=True)

        if len(ext) != len(data):
            raise RuntimeError(
                f"{label}: external diagnostic table has {len(ext)} points, "
                f"but measurement bundle has {len(data)}"
            )
        #endif

        if not np.array_equal(
                ext["source_row"].to_numpy(int),
                np.arange(len(data), dtype=int)):
            raise RuntimeError(
                f"{label}: source_row mismatch while building model diagnostics"
            )
        #endif

        ext["dataset_label"] = label
        if "xs" in data.columns:
            ext["measured_xs"] = data["xs"].to_numpy(float)
        #endif

        # Preserve common measured kinematics even if an older evaluator table
        # did not carry all of them.
        for col in ["xB", "Q2", "t_abs", "phi_deg", "ebeam"]:
            if col not in ext.columns and col in data.columns:
                ext[col] = data[col].to_numpy(float)
            #endif
        #endfor

        blocks.append(ext)
    #endfor

    if not blocks:
        raise RuntimeError("No physics measurements available for model diagnostics")
    #endif
    return pd.concat(blocks, ignore_index=True)
#enddef


def _model_pass_mask(df: pd.DataFrame, model: str, threshold: float) -> np.ndarray:
    return (
        pd.to_numeric(df[f"delta_bh_{model}"], errors="coerce").to_numpy(float)
        <= float(threshold)
    )
#enddef


def select_bundle_from_external_intersection(
        bundle: Dict[str, object],
        selection_table: pd.DataFrame,
        threshold: float,
        models: Sequence[str] = FINAL_MODEL_NAMES) -> pd.DataFrame:
    """Select points that satisfy the same BH-purity threshold in every model."""
    key = str(bundle["key"])
    all_data = bundle["all_data"].reset_index(drop=True)
    ext = selection_table.loc[
        selection_table["dataset"].astype(str) == key
    ].copy().sort_values("source_row").reset_index(drop=True)

    if len(ext) != len(all_data):
        raise RuntimeError(
            f"{bundle['label']}: cannot construct model intersection because "
            "external/current point counts differ"
        )
    #endif

    mask = np.ones(len(ext), dtype=bool)
    for model in models:
        mask &= _model_pass_mask(ext, model, threshold)
    #endfor

    selected = all_data.loc[mask].copy()
    selected["external_bh_model"] = "intersection_" + "_".join(models)
    selected["external_bh_threshold"] = float(threshold)
    return selected
#enddef


def write_model_selection_diagnostics(
        selection: pd.DataFrame,
        physics_bundles: Sequence[Dict[str, object]],
        final_dir: Path) -> None:
    """
    Write diagnostics for KM15/VGG99/GK16 BH-purity selections.

    Nothing produced here changes the selected samples or favors a model based
    on agreement with the measured cross section.  In particular, data/BH
    ratios are diagnostic only.  The physics selection remains Moradi's
        |1 - sigma_BH/sigma_EP| <= threshold,
    so cancellation between DVCS and interference is allowed when their NET
    contribution makes the observed unpolarized cross section BH-like.
    """
    ddir = final_dir / "model_selection_diagnostics"
    ddir.mkdir(parents=True, exist_ok=True)

    df = _external_table_with_measurements(selection, physics_bundles)

    # Signed net non-BH contribution:
    #   sigma_EP - sigma_BH = sigma_DVCS + sigma_INT.
    # This is the quantity whose cancellation is physically relevant here.
    for model in FINAL_MODEL_NAMES:
        bh_col = MODEL_BH_COLUMN[model]
        ep_col = MODEL_EP_COLUMN[model]
        if bh_col in df.columns and ep_col in df.columns:
            bh = pd.to_numeric(df[bh_col], errors="coerce").to_numpy(float)
            ep = pd.to_numeric(df[ep_col], errors="coerce").to_numpy(float)
            good = np.isfinite(bh) & np.isfinite(ep) & (bh != 0.0) & (ep != 0.0)

            net_over_bh = np.full(len(df), np.nan, dtype=float)
            net_over_ep = np.full(len(df), np.nan, dtype=float)
            signed_delta = np.full(len(df), np.nan, dtype=float)
            net_over_bh[good] = (ep[good] - bh[good]) / bh[good]
            net_over_ep[good] = (ep[good] - bh[good]) / ep[good]
            signed_delta[good] = 1.0 - bh[good] / ep[good]

            df[f"net_nonbh_over_bh_{model}"] = net_over_bh
            df[f"net_nonbh_over_ep_{model}"] = net_over_ep
            df[f"signed_delta_bh_{model}"] = signed_delta
        #endif
    #endfor

    # Cross-process BH comparisons.  These diagnose process/elastic-FF
    # implementation differences independently of the GPD-dependent full EP.
    if {
        "km15_bh", "partons_bh_gv08", "partons_bh_vgg99"
    }.issubset(df.columns):
        km = pd.to_numeric(df["km15_bh"], errors="coerce").to_numpy(float)
        gv = pd.to_numeric(df["partons_bh_gv08"], errors="coerce").to_numpy(float)
        vg = pd.to_numeric(df["partons_bh_vgg99"], errors="coerce").to_numpy(float)
        good = np.isfinite(km) & np.isfinite(gv) & np.isfinite(vg) & (km != 0.0) & (gv != 0.0)

        r_gv_km = np.full(len(df), np.nan)
        r_vg_km = np.full(len(df), np.nan)
        r_vg_gv = np.full(len(df), np.nan)
        r_gv_km[good] = gv[good] / km[good]
        r_vg_km[good] = vg[good] / km[good]
        r_vg_gv[good] = vg[good] / gv[good]
        df["bh_gv08_over_gepard"] = r_gv_km
        df["bh_vgg99_over_gepard"] = r_vg_km
        df["bh_vgg99_over_gv08"] = r_vg_gv
    #endif

    # Diagnostic data/BH ratios.  These are deliberately never fed back into
    # selection, model ranking, weighting, or the final systematic.
    if "measured_xs" in df.columns:
        y = pd.to_numeric(df["measured_xs"], errors="coerce").to_numpy(float)
        for model in FINAL_MODEL_NAMES:
            bh_col = MODEL_BH_COLUMN[model]
            if bh_col in df.columns:
                bh = pd.to_numeric(df[bh_col], errors="coerce").to_numpy(float)
                ratio = np.full(len(df), np.nan)
                good = np.isfinite(y) & np.isfinite(bh) & (bh != 0.0)
                ratio[good] = y[good] / bh[good]
                df[f"data_over_reference_bh_{model}"] = ratio
            #endif
        #endfor
    #endif

    # Point-level membership table at the nominal 5% cut.
    membership_cols = [
        c for c in [
            "point_id", "dataset", "dataset_label", "source_row",
            "xB", "Q2", "t_abs", "phi_deg", "ebeam",
            "delta_bh_km15", "delta_bh_vgg99", "delta_bh_gk16",
            "signed_delta_bh_km15", "signed_delta_bh_vgg99",
            "signed_delta_bh_gk16",
            "net_nonbh_over_bh_km15", "net_nonbh_over_bh_vgg99",
            "net_nonbh_over_bh_gk16",
            "km15_bh", "partons_bh_vgg99", "partons_bh_gv08",
            "partons_ep_vgg99", "partons_ep_gk16",
            "bh_gv08_over_gepard", "bh_vgg99_over_gepard",
            "bh_vgg99_over_gv08",
            "measured_xs",
            "data_over_reference_bh_km15",
            "data_over_reference_bh_vgg99",
            "data_over_reference_bh_gk16",
        ] if c in df.columns
    ]
    membership = df[membership_cols].copy()
    for model in FINAL_MODEL_NAMES:
        membership[f"pass_05pct_{model}"] = _model_pass_mask(df, model, 0.05)
    #endfor
    membership["pass_05pct_all_models"] = np.logical_and.reduce([
        membership[f"pass_05pct_{model}"].to_numpy(bool)
        for model in FINAL_MODEL_NAMES
    ])
    membership.to_csv(ddir / "model_selection_membership_05pct.csv", index=False)

    # Per-dataset model counts/fractions at nominal threshold.
    count_rows = []
    status_rows = []
    groups = list(df.groupby(["dataset", "dataset_label"], sort=False))
    groups.append((("ALL", "All final measurements"), df))
    for (dataset, label), d in groups:
        row = {
            "dataset": dataset,
            "dataset_label": label,
            "N_available": int(len(d)),
        }
        for model in FINAL_MODEL_NAMES:
            n = int(_model_pass_mask(d, model, 0.05).sum())
            row[f"N_{model}"] = n
            row[f"fraction_{model}"] = float(n / len(d)) if len(d) else np.nan
        #endfor
        row["N_all_model_intersection"] = int(np.logical_and.reduce([
            _model_pass_mask(d, model, 0.05)
            for model in FINAL_MODEL_NAMES
        ]).sum())
        count_rows.append(row)

        if dataset != "ALL":
            zeros = [
                model for model in FINAL_MODEL_NAMES
                if row[f"N_{model}"] == 0
            ]
            status_rows.append({
                "dataset": dataset,
                "dataset_label": label,
                "nominal_threshold": 0.05,
                "status": "WARNING" if zeros else "OK",
                "zero_selected_models": ",".join(zeros),
                "message": (
                    "One or more purity prescriptions select zero points at 5%; "
                    "inspect purity distributions before interpreting model spread."
                    if zeros else
                    "All purity prescriptions retain at least one point at 5%."
                ),
            })
        #endif
    #endfor
    count_table = pd.DataFrame(count_rows)
    count_table.to_csv(ddir / "model_selection_counts_05pct.csv", index=False)
    pd.DataFrame(status_rows).to_csv(
        ddir / "model_selection_validation_status.csv", index=False
    )

    # Pairwise and three-way overlap versus threshold.  This directly tests
    # whether models are making small boundary migrations or selecting
    # qualitatively different regions.
    pair_rows = []
    triple_rows = []
    thresholds = np.arange(0.01, 0.101, 0.01)
    model_pairs = list(itertools.combinations(FINAL_MODEL_NAMES, 2))
    for (dataset, label), d in groups:
        for threshold in thresholds:
            masks = {
                model: _model_pass_mask(d, model, float(threshold))
                for model in FINAL_MODEL_NAMES
            }
            for model_a, model_b in model_pairs:
                ma = masks[model_a]
                mb = masks[model_b]
                na = int(ma.sum())
                nb = int(mb.sum())
                inter = int((ma & mb).sum())
                union = int((ma | mb).sum())
                pair_rows.append({
                    "dataset": dataset,
                    "dataset_label": label,
                    "threshold": float(threshold),
                    "threshold_percent": 100.0 * float(threshold),
                    "model_a": model_a,
                    "model_b": model_b,
                    "N_a": na,
                    "N_b": nb,
                    "N_intersection": inter,
                    "N_union": union,
                    "jaccard": float(inter / union) if union else np.nan,
                    "fraction_a_shared": float(inter / na) if na else np.nan,
                    "fraction_b_shared": float(inter / nb) if nb else np.nan,
                })
            #endfor

            mall = np.logical_and.reduce([
                masks[model] for model in FINAL_MODEL_NAMES
            ])
            many = np.logical_or.reduce([
                masks[model] for model in FINAL_MODEL_NAMES
            ])
            triple_rows.append({
                "dataset": dataset,
                "dataset_label": label,
                "threshold": float(threshold),
                "threshold_percent": 100.0 * float(threshold),
                "N_km15": int(masks["km15"].sum()),
                "N_vgg99": int(masks["vgg99"].sum()),
                "N_gk16": int(masks["gk16"].sum()),
                "N_all_model_intersection": int(mall.sum()),
                "N_any_model_union": int(many.sum()),
                "three_way_jaccard": (
                    float(mall.sum() / many.sum()) if many.sum() else np.nan
                ),
            })
        #endfor
    #endfor
    pd.DataFrame(pair_rows).to_csv(
        ddir / "model_pairwise_overlap_threshold_scan.csv", index=False
    )
    pd.DataFrame(triple_rows).to_csv(
        ddir / "model_three_way_overlap_threshold_scan.csv", index=False
    )

    # Purity-distribution summaries.
    purity_rows = []
    for (dataset, label), d in groups:
        for model in FINAL_MODEL_NAMES:
            vals = pd.to_numeric(
                d[f"delta_bh_{model}"], errors="coerce"
            ).to_numpy(float)
            vals = vals[np.isfinite(vals)]
            if len(vals) == 0:
                continue
            #endif
            q = np.quantile(vals, [0.0, 0.10, 0.25, 0.50, 0.75, 0.90, 1.0])
            purity_rows.append({
                "dataset": dataset,
                "dataset_label": label,
                "model": model,
                "N": int(len(vals)),
                "delta_min": float(q[0]),
                "delta_q10": float(q[1]),
                "delta_q25": float(q[2]),
                "delta_median": float(q[3]),
                "delta_q75": float(q[4]),
                "delta_q90": float(q[5]),
                "delta_max": float(q[6]),
                "N_le_05pct": int((vals <= 0.05).sum()),
            })
        #endfor
    #endfor
    pd.DataFrame(purity_rows).to_csv(
        ddir / "model_purity_distribution_summary.csv", index=False
    )

    # Summaries of PARTONS-vs-Gepard BH process differences and signed net
    # non-BH terms.  These are central diagnostics for understanding whether
    # the large sample migration is GPD/CFF-driven or process/BH-driven.
    diagnostic_rows = []
    ratio_columns = [
        "bh_gv08_over_gepard",
        "bh_vgg99_over_gepard",
        "bh_vgg99_over_gv08",
    ]
    signed_columns = [
        f"net_nonbh_over_bh_{model}" for model in FINAL_MODEL_NAMES
    ]
    data_ratio_columns = [
        f"data_over_reference_bh_{model}" for model in FINAL_MODEL_NAMES
    ]
    for (dataset, label), d in groups:
        for col in ratio_columns + signed_columns + data_ratio_columns:
            if col not in d.columns:
                continue
            #endif
            vals = pd.to_numeric(d[col], errors="coerce").to_numpy(float)
            vals = vals[np.isfinite(vals)]
            if len(vals) == 0:
                continue
            #endif
            diagnostic_rows.append({
                "dataset": dataset,
                "dataset_label": label,
                "quantity": col,
                "N": int(len(vals)),
                "mean": float(np.mean(vals)),
                "median": float(np.median(vals)),
                "std": float(np.std(vals)),
                "q05": float(np.quantile(vals, 0.05)),
                "q25": float(np.quantile(vals, 0.25)),
                "q75": float(np.quantile(vals, 0.75)),
                "q95": float(np.quantile(vals, 0.95)),
                "min": float(np.min(vals)),
                "max": float(np.max(vals)),
            })
        #endfor
    #endfor
    pd.DataFrame(diagnostic_rows).to_csv(
        ddir / "partons_gepard_process_diagnostic_summary.csv", index=False
    )

    # Selected-only data/reference-BH summaries.  Again, diagnostic only:
    # observed agreement is not allowed to decide which purity model is used.
    selected_data_rows = []
    if "measured_xs" in df.columns:
        for (dataset, label), d in groups:
            for model in FINAL_MODEL_NAMES:
                col = f"data_over_reference_bh_{model}"
                if col not in d.columns:
                    continue
                #endif
                mask = _model_pass_mask(d, model, 0.05)
                vals = pd.to_numeric(
                    d.loc[mask, col], errors="coerce"
                ).to_numpy(float)
                vals = vals[np.isfinite(vals)]
                if len(vals) == 0:
                    continue
                #endif
                selected_data_rows.append({
                    "dataset": dataset,
                    "dataset_label": label,
                    "model": model,
                    "threshold": 0.05,
                    "N": int(len(vals)),
                    "mean_data_over_reference_bh": float(np.mean(vals)),
                    "median_data_over_reference_bh": float(np.median(vals)),
                    "std_data_over_reference_bh": float(np.std(vals)),
                    "q16": float(np.quantile(vals, 0.16)),
                    "q84": float(np.quantile(vals, 0.84)),
                })
            #endfor
        #endfor
    #endif
    pd.DataFrame(selected_data_rows).to_csv(
        ddir / "selected_data_over_reference_bh_05pct.csv", index=False
    )

    # --------------------------------------------------------------
    # Plots: one directory per experiment.  These are intentionally
    # selection diagnostics, not fit-quality/model-ranking plots.
    # --------------------------------------------------------------
    for (dataset, label), d in df.groupby(
            ["dataset", "dataset_label"], sort=False):
        pdir = ddir / str(dataset)
        pdir.mkdir(parents=True, exist_ok=True)

        # Pairwise delta scatter.
        fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8))
        pairs = [
            ("km15", "vgg99"),
            ("km15", "gk16"),
            ("vgg99", "gk16"),
        ]
        finite_all = np.concatenate([
            pd.to_numeric(
                d[f"delta_bh_{m}"], errors="coerce"
            ).to_numpy(float)
            for m in FINAL_MODEL_NAMES
        ])
        finite_all = finite_all[np.isfinite(finite_all)]
        upper = (
            max(0.10, float(np.quantile(finite_all, 0.98)))
            if len(finite_all) else 0.10
        )
        for ax, (ma, mb) in zip(axes, pairs):
            x = 100.0 * d[f"delta_bh_{ma}"].to_numpy(float)
            y = 100.0 * d[f"delta_bh_{mb}"].to_numpy(float)
            lim = 100.0 * upper
            ax.scatter(x, y, s=10, alpha=0.55)
            ax.plot([0.0, lim], [0.0, lim], linestyle="--", linewidth=1.0)
            ax.axvline(5.0, linestyle=":", linewidth=1.0)
            ax.axhline(5.0, linestyle=":", linewidth=1.0)
            ax.set_xlim(0.0, lim)
            ax.set_ylim(0.0, lim)
            ax.set_xlabel(f"{MODEL_DISPLAY[ma]} deviation (%)")
            ax.set_ylabel(f"{MODEL_DISPLAY[mb]} deviation (%)")
            ax.grid(alpha=0.2)
        #endfor
        fig.suptitle(f"{label}: pairwise BH-purity comparison")
        fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.94])
        fig.savefig(pdir / "01_pairwise_delta_bh.png", dpi=180)
        plt.close(fig)

        # Purity distributions.
        fig, ax = plt.subplots(figsize=(8.4, 5.6))
        bins = np.linspace(0.0, 100.0 * upper, 70)
        for model in FINAL_MODEL_NAMES:
            ax.hist(
                100.0 * d[f"delta_bh_{model}"].to_numpy(float),
                bins=bins,
                histtype="step",
                linewidth=1.6,
                label=MODEL_DISPLAY[model],
            )
        #endfor
        ax.axvline(5.0, linestyle=":", linewidth=1.2)
        ax.set_xlabel(r"$|1-\sigma_{\rm BH}/\sigma_{\rm EP}|$ (%)")
        ax.set_ylabel("points")
        ax.set_title(f"{label}: BH-purity distributions")
        ax.legend()
        ax.grid(alpha=0.2)
        fig.tight_layout()
        fig.savefig(pdir / "02_delta_bh_distributions.png", dpi=180)
        plt.close(fig)

        # Signed NET non-BH term.  Cancellation is allowed and physically
        # relevant: zero means sigma_DVCS + sigma_INT = 0.
        available_signed = [
            model for model in FINAL_MODEL_NAMES
            if f"net_nonbh_over_bh_{model}" in d.columns
        ]
        if available_signed:
            fig, ax = plt.subplots(figsize=(8.4, 5.6))
            vals_all = np.concatenate([
                d[f"net_nonbh_over_bh_{model}"].to_numpy(float)
                for model in available_signed
            ])
            vals_all = vals_all[np.isfinite(vals_all)]
            if len(vals_all):
                lo = float(np.quantile(vals_all, 0.02))
                hi = float(np.quantile(vals_all, 0.98))
                if not np.isfinite(lo) or not np.isfinite(hi) or lo == hi:
                    lo, hi = -0.2, 0.2
                #endif
                bins_signed = np.linspace(lo, hi, 80)
                for model in available_signed:
                    ax.hist(
                        d[f"net_nonbh_over_bh_{model}"].to_numpy(float),
                        bins=bins_signed,
                        histtype="step",
                        linewidth=1.6,
                        label=MODEL_DISPLAY[model],
                    )
                #endfor
                ax.axvline(0.0, linestyle=":", linewidth=1.2)
                ax.set_xlim(lo, hi)
            #endif
            ax.set_xlabel(
                r"$(\sigma_{\rm EP}-\sigma_{\rm BH})/\sigma_{\rm BH}$"
                r" $=(\sigma_{\rm DVCS}+\sigma_{\rm INT})/\sigma_{\rm BH}$"
            )
            ax.set_ylabel("points")
            ax.set_title(f"{label}: signed net non-BH contribution")
            ax.legend()
            ax.grid(alpha=0.2)
            fig.tight_layout()
            fig.savefig(pdir / "03_signed_net_nonbh_over_bh.png", dpi=180)
            plt.close(fig)
        #endif

        # PARTONS/Gepard BH process ratios versus kinematics.
        if "bh_gv08_over_gepard" in d.columns:
            kine = [
                ("xB", r"$x_B$"),
                ("Q2", r"$Q^2$ (GeV$^2$)"),
                ("t_abs", r"$|t|$ (GeV$^2$)"),
                ("phi_deg", r"$\phi$ (deg)"),
            ]
            kine = [(c, lab) for c, lab in kine if c in d.columns]
            if kine:
                fig, axes = plt.subplots(
                    1, len(kine), figsize=(4.4 * len(kine), 4.5),
                    squeeze=False
                )
                axes = axes[0]
                for ax, (col, xlabel) in zip(axes, kine):
                    ax.scatter(
                        d[col], d["bh_gv08_over_gepard"],
                        s=8, alpha=0.45, label="GV08 BH / Gepard BH"
                    )
                    if "bh_vgg99_over_gepard" in d.columns:
                        ax.scatter(
                            d[col], d["bh_vgg99_over_gepard"],
                            s=8, alpha=0.45, label="VGG99 BH / Gepard BH"
                        )
                    #endif
                    ax.axhline(1.0, linestyle=":", linewidth=1.0)
                    ax.set_xlabel(xlabel)
                    ax.set_ylabel("BH cross-section ratio")
                    ax.grid(alpha=0.2)
                #endfor
                handles, labels = axes[0].get_legend_handles_labels()
                fig.suptitle(f"{label}: process-level BH implementation comparison")
                if handles:
                    fig.legend(
                        handles, labels, loc="upper center",
                        bbox_to_anchor=(0.5, 0.94), ncol=2
                    )
                #endif
                fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.86])
                fig.savefig(
                    pdir / "04_partons_over_gepard_bh_vs_kinematics.png",
                    dpi=180,
                )
                plt.close(fig)
            #endif
        #endif

        # For each model, show where the 5%-accepted points lie in the signed
        # net non-BH variable.  This is especially useful for Hall A.
        fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.7))
        plotted = False
        for ax, model in zip(axes, FINAL_MODEL_NAMES):
            col = f"net_nonbh_over_bh_{model}"
            if col not in d.columns:
                ax.set_axis_off()
                continue
            #endif
            mask = _model_pass_mask(d, model, 0.05)
            allv = d[col].to_numpy(float)
            sel = allv[mask]
            finite = allv[np.isfinite(allv)]
            if len(finite) == 0:
                ax.set_axis_off()
                continue
            #endif
            lo = float(np.quantile(finite, 0.02))
            hi = float(np.quantile(finite, 0.98))
            if lo == hi:
                lo, hi = -0.2, 0.2
            #endif
            bins_local = np.linspace(lo, hi, 65)
            ax.hist(
                finite, bins=bins_local, histtype="step",
                linewidth=1.2, label="all points"
            )
            if len(sel):
                ax.hist(
                    sel[np.isfinite(sel)], bins=bins_local, histtype="step",
                    linewidth=1.8, label="5% selected"
                )
            #endif
            ax.axvline(0.0, linestyle=":", linewidth=1.0)
            ax.set_xlabel(r"$(\sigma_{\rm EP}-\sigma_{\rm BH})/\sigma_{\rm BH}$")
            ax.set_ylabel("points")
            ax.set_title(MODEL_DISPLAY[model])
            ax.legend()
            ax.grid(alpha=0.2)
            plotted = True
        #endfor
        if plotted:
            fig.suptitle(f"{label}: net non-BH term for nominal 5% selections")
            fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.93])
            fig.savefig(pdir / "05_selected_net_nonbh_05pct.png", dpi=180)
        #endif
        plt.close(fig)
    #endfor

    # Human-readable statement of the non-circular diagnostic policy.
    with open(ddir / "README.txt", "w") as fout:
        fout.write(
            "MODEL-SELECTION DIAGNOSTICS\n"
            "===========================\n\n"
            "The extraction selection remains Moradi's observable-level criterion:\n"
            "  |1 - sigma_BH/sigma_EP| <= threshold.\n\n"
            "Therefore a cancellation sigma_DVCS + sigma_INT ~= 0 is allowed and\n"
            "is physically relevant: the total unpolarized observable is then\n"
            "BH-like even if the two non-BH terms are individually nonzero.\n\n"
            "Measured data/reference-BH ratios in this directory are diagnostics\n"
            "ONLY. They do not choose, reject, rank, or weight KM15/VGG99/GK16.\n"
            "This avoids circularly preferring the model that best describes the\n"
            "same measured cross sections used in the form-factor extraction.\n\n"
            "PARTONS/Gepard BH ratios diagnose process/elastic-form-factor\n"
            "implementation differences separately from the GPD-dependent full\n"
            "EP prediction. Large differences should be understood before the\n"
            "cross-model radius envelope is interpreted as a final systematic.\n"
        )
    #endwith

    print("\n[model diagnostics] nominal 5% selected counts")
    print(count_table.to_string(index=False))
    warning_table = pd.DataFrame(status_rows)
    warnings = warning_table.loc[warning_table["status"] != "OK"]
    if len(warnings):
        print("\n[model diagnostics] WARNING: zero-coverage model/dataset combinations")
        print(warnings.to_string(index=False))
    #endif
    print(f"[model diagnostics] outputs -> {ddir}")
#enddef



def validate_external_partons_bh_consistency(
        selection: pd.DataFrame,
        physics_bundles: Sequence[Dict[str, object]],
        final_dir: Path) -> None:
    """
    Fail-fast validation of the PARTONS BH kinematic/convention mapping.

    The BH subprocess is independent of the GPD/CFF model.  Therefore the
    PARTONS BH calculations and the Gepard BH calculation need not be exactly
    identical, but they must describe the same BH angular structure point by
    point.  Gross discrepancies are an implementation/convention failure, not
    a legitimate GPD-model uncertainty.

    This validation deliberately uses BH-vs-BH comparisons only.  It does not
    use the measured cross section, so it cannot circularly prefer the model
    that happens to describe the data better.

    The bounds below are intentionally very loose compared with the few-percent
    Jo/Hall-A differences observed in validated comparisons:
      * >= 95% of finite points must have 0.5 < PARTONS_BH/Gepard_BH < 2.0
      * the median ratio must lie in 0.8--1.2
      * >= 99% of points must have finite, positive BH values
    These are sanity checks, not physics cuts.
    """
    ddir = final_dir / "model_selection_diagnostics"
    ddir.mkdir(parents=True, exist_ok=True)

    df = _external_table_with_measurements(selection, physics_bundles)
    required = {"km15_bh", "partons_bh_gv08", "partons_bh_vgg99"}
    missing = sorted(required.difference(df.columns))
    if missing:
        raise RuntimeError(
            "Cannot validate PARTONS/Gepard BH consistency; missing columns: "
            + ", ".join(missing)
        )
    #endif

    rows = []
    failures = []
    comparisons = [
        ("gv08", "partons_bh_gv08"),
        ("vgg99", "partons_bh_vgg99"),
    ]

    for (dataset, label), d in df.groupby(
            ["dataset", "dataset_label"], sort=False):
        gepard = pd.to_numeric(
            d["km15_bh"], errors="coerce"
        ).to_numpy(float)

        for process_name, partons_col in comparisons:
            partons = pd.to_numeric(
                d[partons_col], errors="coerce"
            ).to_numpy(float)

            valid = (
                np.isfinite(gepard)
                & np.isfinite(partons)
                & (gepard > 0.0)
                & (partons > 0.0)
            )
            valid_fraction = float(valid.mean()) if len(d) else 0.0

            ratio = np.full(len(d), np.nan, dtype=float)
            ratio[valid] = partons[valid] / gepard[valid]
            finite_ratio = ratio[np.isfinite(ratio)]

            if len(finite_ratio):
                median = float(np.median(finite_ratio))
                q01 = float(np.quantile(finite_ratio, 0.01))
                q05 = float(np.quantile(finite_ratio, 0.05))
                q95 = float(np.quantile(finite_ratio, 0.95))
                q99 = float(np.quantile(finite_ratio, 0.99))
                rmin = float(np.min(finite_ratio))
                rmax = float(np.max(finite_ratio))
                broad_fraction = float(np.mean(
                    (finite_ratio > 0.5) & (finite_ratio < 2.0)
                ))
            else:
                median = np.nan
                q01 = np.nan
                q05 = np.nan
                q95 = np.nan
                q99 = np.nan
                rmin = np.nan
                rmax = np.nan
                broad_fraction = 0.0
            #endif

            reasons = []
            if valid_fraction < 0.99:
                reasons.append(
                    f"finite-positive fraction={valid_fraction:.4f}<0.99"
                )
            #endif
            if (
                not np.isfinite(median)
                or median < 0.8
                or median > 1.2
            ):
                reasons.append(
                    f"median PARTONS/Gepard BH={median:.6g} outside [0.8,1.2]"
                )
            #endif
            if broad_fraction < 0.95:
                reasons.append(
                    f"fraction in 0.5--2.0={broad_fraction:.4f}<0.95"
                )
            #endif

            status = "FAIL" if reasons else "PASS"
            row = {
                "dataset": dataset,
                "dataset_label": label,
                "partons_bh_process": process_name,
                "N_total": int(len(d)),
                "N_valid_positive": int(valid.sum()),
                "valid_positive_fraction": valid_fraction,
                "ratio_median": median,
                "ratio_q01": q01,
                "ratio_q05": q05,
                "ratio_q95": q95,
                "ratio_q99": q99,
                "ratio_min": rmin,
                "ratio_max": rmax,
                "fraction_ratio_between_0p5_and_2": broad_fraction,
                "status": status,
                "failure_reason": "; ".join(reasons),
            }
            rows.append(row)

            if reasons:
                failures.append(
                    f"{label} / {process_name}: " + "; ".join(reasons)
                )
            #endif
        #endfor
    #endfor

    table = pd.DataFrame(rows)
    table.to_csv(
        ddir / "partons_gepard_bh_preflight_validation.csv", index=False
    )

    print("\n[PARTONS BH preflight] process-level BH consistency")
    print(table.to_string(index=False))

    if failures:
        message = (
            "\nPARTONS/Gepard BH preflight FAILED.\n"
            "The alternative-model radius analysis has been stopped BEFORE "
            "closure studies or model-systematic fits because the BH subprocess "
            "does not have a consistent point-by-point kinematic/convention "
            "mapping.\n\n"
            "This is not evidence for a large GPD uncertainty: BH itself is "
            "already inconsistent before the GPD-dependent EP prediction is "
            "considered.\n\n"
            "For the current CLAS12 failure, validate the CLAS12-specific phi "
            "mapping in the PARTONS evaluator with a small BH-only scan before "
            "regenerating its PARTONS predictions. Do not automatically apply "
            "the Jo/Hall-A phi transformation to CLAS12 unless that BH-only "
            "comparison demonstrates it is correct.\n\nFailures:\n  - "
            + "\n  - ".join(failures)
            + "\n\nDiagnostic table: "
            + str(ddir / "partons_gepard_bh_preflight_validation.csv")
        )
        raise RuntimeError(message)
    #endif

    print("[PARTONS BH preflight] PASS: all datasets/processes are convention-consistent")
#enddef



def production_sachs_coefficients(
        fit_result: Dict[str, object],
        which: str) -> Tuple[np.ndarray, np.ndarray]:
    """Extract one production Sachs coefficient vector and its covariance."""
    family = str(fit_result["family"])
    order = int(re.findall(r"\d+", family)[0])
    prefix = "e" if which.upper() == "GE" else "m"
    coeffs = np.asarray(
        [float(fit_result[f"{prefix}{i}"]) for i in range(1, order + 1)],
        dtype=float,
    )

    raw_cov = fit_result.get("shape_covariance_json")
    if raw_cov is None:
        cov = np.diag([
            float(fit_result.get(f"{prefix}{i}_err", np.nan)) ** 2
            for i in range(1, order + 1)
        ])
        return coeffs, cov
    #endif

    full_cov = np.asarray(json.loads(str(raw_cov)), dtype=float)
    if which.upper() == "GE":
        cov = full_cov[:order, :order]
    else:
        cov = full_cov[order:2 * order, order:2 * order]
    #endif
    return coeffs, cov
#enddef


def production_sachs_band(
        fit_result: Dict[str, object],
        q: np.ndarray,
        which: str) -> Tuple[np.ndarray, np.ndarray]:
    """Central GE/GM curve and full-Hessian one-sigma shape band."""
    q = np.asarray(q, dtype=float)
    family = str(fit_result["family"])
    coeffs, cov = production_sachs_coefficients(fit_result, which)
    norm = 1.0 if which.upper() == "GE" else MU_P

    central = norm * sachs_family_value(q, coeffs, family)
    if (
        cov.shape != (len(coeffs), len(coeffs))
        or not np.all(np.isfinite(cov))
    ):
        return central, np.full_like(central, np.nan)
    #endif

    gradients = np.empty((len(q), len(coeffs)), dtype=float)
    for j in range(len(coeffs)):
        step = 1.0e-5 * max(1.0, abs(float(coeffs[j])))
        plus = coeffs.copy()
        minus = coeffs.copy()
        plus[j] += step
        minus[j] -= step
        gradients[:, j] = norm * (
            sachs_family_value(q, plus, family)
            - sachs_family_value(q, minus, family)
        ) / (2.0 * step)
    #endfor

    variance = np.einsum("ij,jk,ik->i", gradients, cov, gradients)
    sigma = np.sqrt(np.maximum(variance, 0.0))
    return central, sigma
#enddef


def save_preliminary_radius_summary_figure(
        fit_table: pd.DataFrame,
        final_summary: pd.DataFrame,
        figures_dir: Path) -> None:
    """Compact preliminary rE/rM comparison across the three purity models."""
    figures_dir.mkdir(parents=True, exist_ok=True)
    models = list(FINAL_MODEL_NAMES)
    labels = {"km15": "KM15", "vgg99": "VGG99", "gk16": "GK16"}

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.8))
    for ax, quantity, error_col, ylabel in [
        (axes[0], "rE_fm", "rE_fit_err_fm", r"$r_E$ (fm)"),
        (axes[1], "rM_fm", "rM_fit_err_fm", r"$r_M$ (fm)"),
    ]:
        for i, model in enumerate(models):
            row = fit_table.loc[fit_table["model"] == model].iloc[0]
            ax.errorbar(
                i, float(row[quantity]), yerr=float(row[error_col]),
                fmt="o", capsize=3.0, markersize=6.0,
            )
        #endfor
        ax.set_xticks(range(len(models)), [labels[m] for m in models])
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", alpha=0.25)
    #endfor

    summary = final_summary.iloc[0]
    axes[0].axhline(float(summary["rE_fm"]), linewidth=1.0, linestyle="--")
    axes[1].axhline(float(summary["rM_fm"]), linewidth=1.0, linestyle="--")
    fig.suptitle("Jo 2015 + Lee 2026: 5% BH-purity model comparison", y=0.98)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(figures_dir / "01_preliminary_radii.png", dpi=180)
    plt.close(fig)
#enddef


def save_final_threshold_stability_figure(
        threshold_table: pd.DataFrame,
        figures_dir: Path) -> None:
    """Chosen-family radius stability versus BH-purity threshold."""
    figures_dir.mkdir(parents=True, exist_ok=True)
    d = threshold_table.loc[
        (threshold_table["error_mode"].astype(str) == "published_errors")
        & threshold_table["valid"].astype(bool)
    ].copy()

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.8), sharex=True)
    for model in FINAL_MODEL_NAMES:
        m = d.loc[d["model"].astype(str) == model].sort_values("threshold")
        x = 100.0 * m["threshold"].to_numpy(float)
        axes[0].errorbar(
            x, m["rE_fm"], yerr=m["rE_fit_err_fm"],
            fmt="o-", markersize=3.5, linewidth=1.0, label=model.upper(),
        )
        axes[1].errorbar(
            x, m["rM_fm"], yerr=m["rM_fit_err_fm"],
            fmt="o-", markersize=3.5, linewidth=1.0, label=model.upper(),
        )
    #endfor

    for ax, ylabel in zip(axes, [r"$r_E$ (fm)", r"$r_M$ (fm)"]):
        ax.axvspan(3.0, 7.0, alpha=0.08)
        ax.axvline(5.0, linewidth=0.9, linestyle="--")
        ax.set_xlabel(r"BH-purity threshold $\delta_{\rm BH}$ (%)")
        ax.set_ylabel(ylabel)
        ax.set_xlim(1.0, 10.0)
        ax.grid(alpha=0.22)
    #endfor
    axes[0].legend(fontsize=8)
    fig.suptitle(
        "Jo 2015 + Lee 2026: chosen-family threshold stability "
        "(published experimental errors)",
        y=0.98,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(figures_dir / "02_threshold_stability.png", dpi=180)
    plt.close(fig)
#enddef


def save_final_form_factor_comparison(
        fit_rows: Sequence[Dict[str, object]],
        figures_dir: Path) -> None:
    """
    Compare the chosen Jo+Lee Sachs parameterization with elastic references.

    The top row shows GE and GM.  The bottom row shows GE/GD and
    GM/(mu_p GD).  A1/Bernauer Rosenbluth points are displayed directly.
    The KM15 nominal curve receives its full Hessian band; VGG99/GK16 are
    shown as central alternative-purity curves.
    """
    figures_dir.mkdir(parents=True, exist_ok=True)
    qmax = 0.60
    q = np.linspace(0.0, qmax, 600)
    gd = (1.0 + q / 0.71) ** -2
    a1 = bernauer_rosenbluth_data()
    a1 = a1.loc[a1["Q2"] <= qmax].copy()

    by_model = {str(r["model"]): r for r in fit_rows}
    fig, axes = plt.subplots(
        2, 2, figsize=(12.0, 8.3), sharex="col",
        gridspec_kw={"height_ratios": [1.0, 1.0]},
    )

    for model in FINAL_MODEL_NAMES:
        row = by_model[model]
        ge, ge_sigma = production_sachs_band(row, q, "GE")
        gm, gm_sigma = production_sachs_band(row, q, "GM")
        label = model.upper()

        line_ge, = axes[0, 0].plot(q, ge, linewidth=1.5, label=label)
        line_gm, = axes[0, 1].plot(q, gm, linewidth=1.5, label=label)
        axes[1, 0].plot(q, ge / gd, linewidth=1.5)
        axes[1, 1].plot(q, gm / (MU_P * gd), linewidth=1.5)

        if model == FINAL_NOMINAL_MODEL:
            axes[0, 0].fill_between(
                q, ge - ge_sigma, ge + ge_sigma,
                alpha=0.16, color=line_ge.get_color(),
                label="KM15 68% Hessian band",
            )
            axes[0, 1].fill_between(
                q, gm - gm_sigma, gm + gm_sigma,
                alpha=0.16, color=line_gm.get_color(),
            )
            axes[1, 0].fill_between(
                q, (ge - ge_sigma) / gd, (ge + ge_sigma) / gd,
                alpha=0.16, color=line_ge.get_color(),
            )
            axes[1, 1].fill_between(
                q, (gm - gm_sigma) / (MU_P * gd),
                (gm + gm_sigma) / (MU_P * gd),
                alpha=0.16, color=line_gm.get_color(),
            )
        #endif
    #endfor

    for label, (ge_ref, gm_ref) in elastic_reference_curves(q).items():
        axes[0, 0].plot(q, ge_ref, linewidth=1.0, linestyle="--", label=label)
        axes[0, 1].plot(q, gm_ref, linewidth=1.0, linestyle="--")
        axes[1, 0].plot(q, ge_ref / gd, linewidth=1.0, linestyle="--")
        axes[1, 1].plot(q, gm_ref / (MU_P * gd), linewidth=1.0, linestyle="--")
    #endfor

    axes[0, 0].errorbar(
        a1["Q2"], a1["GE"], yerr=a1["GE_err"],
        fmt="o", fillstyle="none", markersize=3.2, capsize=1.5,
        linewidth=0.8, label="A1/Bernauer Rosenbluth",
    )
    axes[0, 1].errorbar(
        a1["Q2"], a1["GM"], yerr=a1["GM_err"],
        fmt="o", fillstyle="none", markersize=3.2, capsize=1.5,
        linewidth=0.8,
    )
    axes[1, 0].errorbar(
        a1["Q2"], a1["GE"] / ((1.0 + a1["Q2"] / 0.71) ** -2),
        yerr=a1["GE_err"] / ((1.0 + a1["Q2"] / 0.71) ** -2),
        fmt="o", fillstyle="none", markersize=3.2, capsize=1.5,
        linewidth=0.8,
    )
    axes[1, 1].errorbar(
        a1["Q2"], a1["GM"] / (
            MU_P * ((1.0 + a1["Q2"] / 0.71) ** -2)
        ),
        yerr=a1["GM_err"] / (
            MU_P * ((1.0 + a1["Q2"] / 0.71) ** -2)
        ),
        fmt="o", fillstyle="none", markersize=3.2, capsize=1.5,
        linewidth=0.8,
    )

    axes[0, 0].set_ylabel(r"$G_E$")
    axes[0, 1].set_ylabel(r"$G_M$")
    axes[1, 0].set_ylabel(r"$G_E/G_D$")
    axes[1, 1].set_ylabel(r"$G_M/(\mu_p G_D)$")
    for ax in axes.ravel():
        ax.set_xlim(0.0, qmax)
        ax.grid(alpha=0.22)
    #endfor
    for ax in axes[1, :]:
        ax.set_xlabel(r"$Q^2=|t|$ (GeV$^2$)")
        ax.axhline(1.0, linewidth=0.8, linestyle=":")
    #endfor

    axes[0, 0].legend(fontsize=7, ncol=2)
    family = str(fit_rows[0]["family"])
    fig.suptitle(
        f"Jo 2015 + Lee 2026: chosen Sachs family {family} vs elastic data/fits",
        y=0.985,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.955))
    fig.savefig(figures_dir / "03_form_factors_vs_elastic.png", dpi=180)
    plt.close(fig)

    # Ratio most directly sensitive to the difference between electric and
    # magnetic slopes.
    fig, ax = plt.subplots(figsize=(7.0, 5.1))
    for model in FINAL_MODEL_NAMES:
        row = by_model[model]
        ge, _ = production_sachs_band(row, q, "GE")
        gm, _ = production_sachs_band(row, q, "GM")
        ax.plot(q, MU_P * ge / gm, linewidth=1.5, label=model.upper())
    #endfor
    for label, (ge_ref, gm_ref) in elastic_reference_curves(q).items():
        ax.plot(q, MU_P * ge_ref / gm_ref, linewidth=1.0, linestyle="--",
                label=label)
    #endfor
    a1_ratio = MU_P * a1["GE"] / a1["GM"]
    a1_ratio_err = np.abs(a1_ratio) * np.sqrt(
        (a1["GE_err"] / a1["GE"]) ** 2
        + (a1["GM_err"] / a1["GM"]) ** 2
    )
    ax.errorbar(
        a1["Q2"], a1_ratio, yerr=a1_ratio_err,
        fmt="o", fillstyle="none", markersize=3.2, capsize=1.5,
        linewidth=0.8, label="A1/Bernauer Rosenbluth",
    )
    ax.axhline(1.0, linewidth=0.8, linestyle=":")
    ax.set_xlim(0.0, qmax)
    ax.set_xlabel(r"$Q^2=|t|$ (GeV$^2$)")
    ax.set_ylabel(r"$\mu_p G_E/G_M$")
    ax.grid(alpha=0.22)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(figures_dir / "04_muGE_over_GM_vs_elastic.png", dpi=180)
    plt.close(fig)
#enddef


def save_final_family_ranking_figure(
        ranking: pd.DataFrame,
        figures_dir: Path) -> None:
    """KM15-only production closure-family ranking."""
    figures_dir.mkdir(parents=True, exist_ok=True)
    objective_col = "km15_combined_RMS_objective_fm"
    d = ranking.copy().sort_values(objective_col)
    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    y = np.arange(len(d))
    ax.barh(y, d[objective_col].to_numpy(float))
    ax.set_yticks(y, d["family"].astype(str).tolist())
    ax.invert_yaxis()
    ax.set_xlabel("KM15 combined RMS bias-variance objective (fm)")
    ax.set_title("Jo 2015 + Lee 2026 KM15 closure-family ranking")
    ax.grid(axis="x", alpha=0.22)
    for i, (_, row) in enumerate(d.iterrows()):
        if not bool(row["production_eligible_km15"]):
            x = float(row[objective_col])
            if not np.isfinite(x):
                x = 0.0
            #endif
            ax.text(x, i, "  ineligible", va="center", fontsize=8)
        #endif
    #endfor
    fig.tight_layout()
    fig.savefig(figures_dir / "05_closure_family_ranking.png", dpi=180)
    plt.close(fig)
#enddef


def write_final_analysis_readme(
        summary: pd.DataFrame,
        fit_table: pd.DataFrame,
        figures_dir: Path,
        summary_dir: Path,
        diagnostics_dir: Path,
        output_path: Path) -> None:
    """Human-readable map of the deliberately consolidated final products."""
    s = summary.iloc[0]
    lines = [
        "PRELIMINARY CLAS-ONLY BH FORM-FACTOR / RADIUS ANALYSIS",
        "=" * 62,
        "",
        "Production ensemble: CLAS6 Jo 2015 + CLAS12 Lee 2026.",
        "Excluded by default: CLAS6 Saylor 2018 and Hall A Defurne 2015.",
        "Nominal BH-purity prescription: KM15 at 5%.",
        (
            f"Chosen Sachs family: {s['chosen_family']} "
            f"(closure rank {int(s['chosen_family_closure_rank'])}; "
            f"bias-rank-1 family {s['bias_rank1_family']}; "
            f"fallback used: {bool(s['family_fallback_used'])})."
        ),
        (
            "Family selection: eligible families are ordered by the closure "
            "bias+variance objective. The script tries them in that order and "
            "accepts the first family giving valid nominal 5% production fits "
            "for all three purity prescriptions. Threshold-scan failures do "
            "not disqualify a family; invalid scan points are recorded and "
            "skipped when building the 3--7% selection envelope. Every rejected "
            "nominal-family attempt is recorded in "
            "summary/family_selection_attempts.csv."
        ),
        (
            "Closure-family eligibility: minimum 33% valid replicas in every "
            "truth scenario for every purity-model closure study."
        ),
        (
            "Candidate families: P1--P4, IP1--IP4, CF2--CF4. Failed replicas "
            "are retained as diagnostics rather than automatically excluding "
            "a family unless validity falls below the 33% threshold."
        ),
        "",
        "Preliminary radii",
        "-----------------",
        (
            f"rE = {float(s['rE_fm']):.5f} +/- "
            f"{float(s['rE_fit_err_fm']):.5f} (fit) +/- "
            f"{float(s['rE_method_sys_fm']):.5f} (method) fm"
        ),
        (
            f"rM = {float(s['rM_fm']):.5f} +/- "
            f"{float(s['rM_fit_err_fm']):.5f} (fit) +/- "
            f"{float(s['rM_method_sys_fm']):.5f} (method) fm"
        ),
        "",
        (
            "Method systematic = quadrature of the KM15 3--7% BH-selection "
            "systematic and the KM15 closure/extrapolation bias systematic."
        ),
        (
            "BH-selection systematic = maximum absolute shift from the nominal "
            "KM15 5% published-errors result using KM15 only over the 3--7% "
            "BH-purity window. VGG99/GK16 are external diagnostics only."
        ),
        "",
        "Primary files",
        "-------------",
        f"Summary table: {summary_dir / 'final_results.csv'}",
        f"Model x threshold table: {summary_dir / 'model_threshold_summary.csv'}",
        f"Closure ranking: {summary_dir / 'closure_summary.csv'}",
        f"Figures: {figures_dir}",
        f"Detailed diagnostics: {diagnostics_dir}",
        "",
        "The detailed diagnostics are retained for auditability but should not "
        "be the first place to look for the result.",
    ]
    output_path.write_text("\n".join(lines) + "\n")
#enddef



def save_nominal_model_selection_kinematic_histograms(
        physics_bundles: Sequence[Dict[str, object]],
        selected: Dict[str, Dict[float, Sequence[Dict[str, object]]]],
        outdir: Path) -> None:
    """Plot nominal-5% selected counts and selection fractions vs kinematics."""
    outdir.mkdir(parents=True, exist_ok=True)
    bundle_by_key = {str(b["key"]): b for b in physics_bundles}
    row_defs = [
        ("jo2015", "CLAS6 Jo 2015"),
        ("pass1", "CLAS12 Lee 2026"),
        ("combined", "Jo + Lee"),
    ]
    variables = [
        ("Q2", r"$Q^2$ (GeV$^2$)"),
        ("xB", r"$x_B$"),
        ("t_abs", r"$|t|$ (GeV$^2$)"),
    ]

    all_combined = pd.concat(
        [bundle_by_key[k]["all_data"] for k in ["jo2015", "pass1"]],
        ignore_index=True,
    )
    bins = {}
    for var, _ in variables:
        vals = all_combined[var].to_numpy(float)
        vals = vals[np.isfinite(vals)]
        bins[var] = np.linspace(float(np.min(vals)), float(np.max(vals)), 16)
    #endfor

    spec5 = {}
    for model in FINAL_MODEL_NAMES:
        spec5[model] = {str(sp["key"]): sp["data"] for sp in selected[model][0.05]}
    #endfor

    for mode in ["counts", "fraction"]:
        fig, axes = plt.subplots(3, 3, figsize=(14.0, 10.0), sharex="col")
        for irow, (row_key, row_label) in enumerate(row_defs):
            if row_key == "combined":
                denominator = all_combined
            else:
                denominator = bundle_by_key[row_key]["all_data"]
            #endif
            for icol, (var, xlabel) in enumerate(variables):
                ax = axes[irow, icol]
                edges = bins[var]
                centers = 0.5 * (edges[:-1] + edges[1:])
                denom_counts, _ = np.histogram(denominator[var].to_numpy(float), bins=edges)
                for model in FINAL_MODEL_NAMES:
                    if row_key == "combined":
                        data = pd.concat([spec5[model]["jo2015"], spec5[model]["pass1"]], ignore_index=True)
                    else:
                        data = spec5[model][row_key]
                    #endif
                    selected_counts, _ = np.histogram(data[var].to_numpy(float), bins=edges)
                    if mode == "counts":
                        y = selected_counts.astype(float)
                    else:
                        y = np.divide(selected_counts, denom_counts, out=np.full_like(selected_counts, np.nan, dtype=float), where=denom_counts > 0)
                    #endif
                    ax.step(centers, y, where="mid", linewidth=1.5, label=model.upper())
                #endfor
                if irow == 0:
                    ax.set_title(xlabel)
                #endif
                if irow == 2:
                    ax.set_xlabel(xlabel)
                #endif
                if icol == 0:
                    ax.set_ylabel(row_label + ("\nSelected points / bin" if mode == "counts" else "\nSelected / all"))
                #endif
                if mode == "fraction":
                    ax.set_ylim(0.0, 1.05)
                #endif
                ax.grid(alpha=0.2)
        #endfor
        handles, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 0.965))
        title = "Nominal 5% BH-purity selection: accepted points vs kinematics" if mode == "counts" else "Nominal 5% BH-purity selection fraction vs kinematics"
        fig.suptitle(title, y=0.995)
        fig.tight_layout(rect=(0, 0, 1, 0.925))
        filename = "05_selected_counts_vs_kinematics_05pct.png" if mode == "counts" else "06_selection_fraction_vs_kinematics_05pct.png"
        fig.savefig(outdir / filename, dpi=180)
        plt.close(fig)
    #endfor
#enddef



def build_data_driven_bh_selection_diagnostics(
        physics_bundles: Sequence[Dict[str, object]],
        nominal_specs: Sequence[Dict[str, object]],
        chosen_family: str,
        diagnostics_dir: Path) -> Dict[str, object]:
    """Compare nominal KM15 selection with data-driven pure-BH compatibility.

    Kelly GE/GM are used only as a fixed external reference to construct the BH
    prediction.  The measured cross section is never refit before selection.
    Two matched-statistics selectors are formed independently in each dataset:
    smallest absolute fractional BH deviation and smallest absolute standardized
    BH pull.  The latter uses statistical + point-to-point errors but excludes
    correlated normalization and the Moradi BH-purity term from the selector.
    """
    nominal_by_key = {str(spec["key"]): spec for spec in nominal_specs}
    point_rows = []
    alt_specs = {"kelly_fractional_matchedN": [], "kelly_pull_matchedN": []}
    count_rows = []

    for bundle in physics_bundles:
        key = str(bundle["key"])
        d = bundle.get("all_data", bundle.get("data", bundle.get("set5"))).copy()
        if d is None or len(d) == 0:
            continue
        #endif
        q = d["t_abs"].to_numpy(float)
        ge, gm = kelly_sachs(q)
        tau = q / (4.0 * MP2)
        f1 = (ge + tau * gm) / (1.0 + tau)
        f2 = (gm - ge) / (1.0 + tau)
        bh = bh_from_f1f2(
            d["bh_A"].to_numpy(float), d["bh_B"].to_numpy(float),
            d["bh_C"].to_numpy(float), f1, f2,
        )
        xs = d["xs"].to_numpy(float)
        point_err = dataset_point_errors(d, str(bundle["kind"]), 0.05, False)
        frac = np.divide(xs - bh, bh, out=np.full_like(xs, np.nan), where=np.abs(bh) > 1.0e-15)
        pull = np.divide(xs - bh, point_err, out=np.full_like(xs, np.nan), where=point_err > 0.0)
        d["kelly_bh_xs"] = bh
        d["kelly_bh_fractional_deviation"] = frac
        d["kelly_bh_pull"] = pull

        n_nom = int(len(nominal_by_key[key]["data"]))
        finite_frac = np.flatnonzero(np.isfinite(frac))
        finite_pull = np.flatnonzero(np.isfinite(pull))
        idx_frac = finite_frac[np.argsort(np.abs(frac[finite_frac]))[:min(n_nom, len(finite_frac))]]
        idx_pull = finite_pull[np.argsort(np.abs(pull[finite_pull]))[:min(n_nom, len(finite_pull))]]
        sel_frac = d.iloc[np.sort(idx_frac)].copy()
        sel_pull = d.iloc[np.sort(idx_pull)].copy()
        alt_specs["kelly_fractional_matchedN"].append(bundle_to_measurement_spec(bundle, sel_frac))
        alt_specs["kelly_pull_matchedN"].append(bundle_to_measurement_spec(bundle, sel_pull))

        nominal_index = set(nominal_by_key[key]["data"].index.tolist())
        frac_index = set(sel_frac.index.tolist())
        pull_index = set(sel_pull.index.tolist())
        def overlap(a, b):
            union = a | b
            return len(a & b), (len(a & b) / len(union) if union else np.nan)
        #enddef
        nf, jf = overlap(nominal_index, frac_index)
        npull, jpull = overlap(nominal_index, pull_index)
        count_rows.extend([
            {"dataset": key, "selector": "KM15_5pct", "N": n_nom, "N_overlap_KM15": n_nom, "Jaccard_with_KM15": 1.0},
            {"dataset": key, "selector": "kelly_fractional_matchedN", "N": len(sel_frac), "N_overlap_KM15": nf, "Jaccard_with_KM15": jf, "effective_abs_cut": float(np.nanmax(np.abs(sel_frac["kelly_bh_fractional_deviation"].to_numpy(float)))) if len(sel_frac) else np.nan},
            {"dataset": key, "selector": "kelly_pull_matchedN", "N": len(sel_pull), "N_overlap_KM15": npull, "Jaccard_with_KM15": jpull, "effective_abs_cut": float(np.nanmax(np.abs(sel_pull["kelly_bh_pull"].to_numpy(float)))) if len(sel_pull) else np.nan},
        ])
        for irow, row in d.iterrows():
            point_rows.append({
                "dataset": key, "row_index": irow,
                "t_abs": float(row["t_abs"]), "xs": float(row["xs"]),
                "kelly_bh_xs": float(row["kelly_bh_xs"]),
                "kelly_bh_fractional_deviation": float(row["kelly_bh_fractional_deviation"]),
                "kelly_bh_pull": float(row["kelly_bh_pull"]),
                "selected_KM15_5pct": irow in nominal_index,
                "selected_kelly_fractional_matchedN": irow in frac_index,
                "selected_kelly_pull_matchedN": irow in pull_index,
            })
        #endfor
    #endfor

    point_table = pd.DataFrame(point_rows)
    count_table = pd.DataFrame(count_rows)
    point_table.to_csv(diagnostics_dir / "data_driven_bh_deviance_points.csv", index=False)
    count_table.to_csv(diagnostics_dir / "data_driven_bh_deviance_overlap.csv", index=False)

    if len(point_table):
        fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.0))
        for dataset, part in point_table.groupby("dataset"):
            axes[0].hist(
                part["kelly_bh_fractional_deviation"].to_numpy(float),
                bins=60, histtype="step", density=True, label=str(dataset),
            )
            axes[1].hist(
                part["kelly_bh_pull"].to_numpy(float),
                bins=60, histtype="step", density=True, label=str(dataset),
            )
        #endfor
        axes[0].set_xlabel(r"$(\sigma_{\rm data}-\sigma_{\rm BH}^{\rm Kelly})/\sigma_{\rm BH}^{\rm Kelly}$")
        axes[1].set_xlabel(r"$(\sigma_{\rm data}-\sigma_{\rm BH}^{\rm Kelly})/\delta\sigma_{\rm ptp}$")
        axes[0].set_ylabel("Normalized entries")
        axes[0].set_title("Fractional deviation from fixed Kelly BH")
        axes[1].set_title("Standardized deviation from fixed Kelly BH")
        axes[0].legend()
        axes[1].legend()
        fig.suptitle("Data-driven pure-BH compatibility diagnostics")
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        fig.savefig(diagnostics_dir / "data_driven_bh_deviance_distributions.png", dpi=180)
        plt.close(fig)
    #endif

    fit_rows = []
    for selector, specs in alt_specs.items():
        result = fit_sachs_family_multi_measurements(
            specs, family=chosen_family, bh_cut=0.05,
            add_moradi_bh_systematic=True,
        )
        result["selector"] = selector
        fit_rows.append(result)
        print(
            f"[BH-deviance diagnostic] {selector}: N={result['N']} "
            f"chi2/ndf={result['chi2_ndof']:.4f} "
            f"rE={result['rE_fm']:.5f} fm rM={result['rM_fm']:.5f} fm"
        )
    #endfor
    fit_table = pd.DataFrame(fit_rows)
    fit_table.to_csv(diagnostics_dir / "data_driven_bh_deviance_fits.csv", index=False)
    return {"fit_table": fit_table, "alternative_specs": alt_specs}
#enddef


def run_final_model_selected_analysis(
        bundles: Sequence[Dict[str, object]],
        args,
        root_outdir: Path) -> None:
    """
    End-to-end final analysis after the external PARTONS model-selection stage.

    Preliminary CLAS-only production analysis using CLAS6 Jo 2015 + CLAS12
    Lee 2026.  Saylor 2018 and Hall A Defurne 2015 are intentionally excluded
    from the default/preliminary ensemble and remain opt-in diagnostics.

    The three BH-purity models define alternative selected samples.  They are
    not averaged point-by-point.  A common Sachs family is selected by
    model-specific closure studies and then fit independently to each
    model-selected Jo+Lee sample.  KM15 remains the nominal central
    prescription.  For the preliminary uncertainty, model and threshold
    dependence are treated jointly over the 3--7% BH-purity window so that
    strongly correlated sample-migration effects are not double counted.
    """
    selection_path = Path(args.bh_model_selection_results).expanduser().resolve()
    selection = load_external_bh_model_selection(selection_path)

    physics_bundles = [
        b for b in bundles if str(b["key"]) in FINAL_PHYSICS_DATASETS
    ]
    found = {str(b["key"]) for b in physics_bundles}
    missing = [key for key in FINAL_PHYSICS_DATASETS if key not in found]
    if missing:
        raise RuntimeError(
            "Preliminary final model-selected analysis requires "
            "CLAS6 Jo 2015 + CLAS12 Lee 2026. Missing: " + ", ".join(missing)
        )
    #endif

    final_dir = root_outdir / "final_analysis"
    summary_dir = final_dir / "summary"
    figures_dir = final_dir / "figures"
    diagnostics_dir = final_dir / "diagnostics"
    closure_dir = diagnostics_dir / "closure"
    for directory in [
        final_dir, summary_dir, figures_dir, diagnostics_dir, closure_dir
    ]:
        directory.mkdir(parents=True, exist_ok=True)
    #endfor

    # Build all model/threshold selections once and save the counts.
    selected = {}
    count_rows = []
    thresholds = np.arange(0.01, 0.101, 0.01)
    for model in FINAL_MODEL_NAMES:
        selected[model] = {}
        for threshold in thresholds:
            specs = []
            row = {
                "model": model,
                "threshold": float(threshold),
                "threshold_percent": 100.0 * float(threshold),
            }
            for bundle in physics_bundles:
                d = select_bundle_from_external_model(
                    bundle, selection, model, float(threshold)
                )
                specs.append(bundle_to_measurement_spec(bundle, d))
                row["N_" + str(bundle["key"])] = int(len(d))
            #endfor
            row["N_total"] = int(sum(len(spec["data"]) for spec in specs))
            selected[model][round(float(threshold), 2)] = specs
            count_rows.append(row)
        #endfor
    #endfor
    count_table = pd.DataFrame(count_rows)
    count_table.to_csv(
        diagnostics_dir / "model_threshold_selected_counts.csv", index=False
    )
    save_nominal_model_selection_kinematic_histograms(
        physics_bundles, selected, diagnostics_dir
    )

    # Comprehensive purity/process diagnostics are written BEFORE any bias
    # study or production fit.  They never feed measured-data agreement back
    # into model selection, avoiding circular model preference.
    write_model_selection_diagnostics(selection, physics_bundles, diagnostics_dir)

    # Fail loudly before doing expensive closure studies or quoting a model
    # envelope if the PARTONS BH subprocess does not reproduce the same
    # point-by-point BH angular structure as Gepard.  This catches phi/convention
    # mistakes without using the measured cross sections.
    validate_external_partons_bh_consistency(
        selection, physics_bundles, diagnostics_dir
    )

    # Run the expensive production closure/bias study for KM15 only.
    # VGG99/GK16 remain external selection diagnostics and do not consume
    # replica statistics or influence the production family/bias systematic.
    model = FINAL_NOMINAL_MODEL
    model_dir = closure_dir / model
    bias_dirs = {model: model_dir}
    print(f"\n[final model analysis] starting production radius-bias study for {model}")

    model_bundles = []
    specs5 = selected[model][0.05]
    spec_by_key = {str(spec["key"]): spec for spec in specs5}
    for bundle in physics_bundles:
        bcopy = dict(bundle)
        bcopy["set5"] = spec_by_key[str(bundle["key"])]["data"].copy()
        model_bundles.append(bcopy)
    #endfor

    if args.run_radius_bias_study:
        run_radius_bias_variance_study(model_bundles, args, model_dir)
    elif not (model_dir / "radius_bias_extended_eligibility_ranking.csv").exists():
        raise RuntimeError(
            f"KM15 bias ranking missing: {model_dir}. Run with "
            "--run-radius-bias-study --radius-bias-extended-truths."
        )
    #endif

    bias_rank1_family, family_ranking = aggregate_model_bias_rankings(
        bias_dirs, summary_dir
    )

    # ------------------------------------------------------------------
    # Ordered production-validation fallback.
    #
    # The closure ranking supplies an ordered list of rankable families.
    # Replica convergence fractions are diagnostic only and impose no veto.
    # Starting from the best bias+variance score, require ONLY that the family
    # produce a valid nominal 5% KM15 Moradi-error production fit. VGG99/GK16
    # remain external diagnostics and cannot veto a family. If the nominal fit
    # fails, retain the failure in the audit table
    # and automatically try the next-best eligible family.
    #
    # Threshold-scan failures do NOT disqualify a family. They are retained as
    # diagnostics and skipped when constructing the local selection envelope.
    # ------------------------------------------------------------------
    eligible_ranking = family_ranking.loc[
        family_ranking["production_eligible_km15"].astype(bool)
    ].copy().reset_index(drop=True)
    if len(eligible_ranking) == 0:
        raise RuntimeError(
            "No closure-rankable Sachs family is available for production "
            "validation."
        )
    #endif

    family_attempt_rows = []
    chosen_family = None
    chosen_family_rank = None
    chosen_5pct_fit_rows = None

    for candidate_index, rank_row in eligible_ranking.iterrows():
        candidate_family = str(rank_row["family"])
        candidate_rank = int(candidate_index) + 1
        candidate_score = float(
            rank_row["km15_combined_RMS_objective_fm"]
        )
        print(
            f"\n[family fallback] trying rank {candidate_rank}: "
            f"{candidate_family} "
            f"(worst objective={candidate_score:.6f} fm)"
        )

        attempt = {
            "closure_rank": candidate_rank,
            "family": candidate_family,
            "km15_combined_RMS_objective_fm": candidate_score,
            "passed_05pct_production": False,
            "selected_for_final_analysis": False,
            "failure_details": "",
        }

        candidate_fit_rows = []
        failure_details = []

        # Required production validation is KM15 only. VGG99/GK16 are
        # external diagnostics and cannot veto the production family.
        model = FINAL_NOMINAL_MODEL
        try:
            result = fit_sachs_family_multi_measurements(
                selected[model][0.05],
                family=candidate_family,
                bh_cut=0.05,
                add_moradi_bh_systematic=True,
            )
            result["model"] = model
            result["threshold"] = 0.05
            result["family_validation_candidate"] = candidate_family
            candidate_fit_rows.append(result)
            if not bool(result.get("valid", False)):
                failure_details.append(f"invalid 5% Moradi fit: {model}")
            #endif
        except Exception as exc:
            failure_details.append(
                f"exception in 5% Moradi fit {model}: "
                f"{type(exc).__name__}: {exc}"
            )
        #endtry

        passed_05pct = (
            len(candidate_fit_rows) == 1
            and bool(candidate_fit_rows[0].get("valid", False))
        )
        attempt["passed_05pct_production"] = bool(passed_05pct)

        if passed_05pct:
            chosen_family = candidate_family
            chosen_family_rank = candidate_rank
            chosen_5pct_fit_rows = candidate_fit_rows
            attempt["selected_for_final_analysis"] = True
            print(
                f"[family fallback] ACCEPTED {candidate_family} at closure "
                f"rank {candidate_rank}: nominal 5% production fits valid"
            )
        else:
            print(
                f"[family fallback] REJECTED {candidate_family} at closure "
                f"rank {candidate_rank}: "
                + "; ".join(failure_details)
            )
        #endif

        attempt["failure_details"] = "; ".join(failure_details)
        family_attempt_rows.append(attempt)

        # Rewrite on every attempt so a job interrupted during fallback still
        # leaves a complete record of everything tested so far.
        pd.DataFrame(family_attempt_rows).to_csv(
            summary_dir / "family_selection_attempts.csv", index=False
        )

        if chosen_family is not None:
            break
        #endif
    #endfor

    if chosen_family is None:
        raise RuntimeError(
            "Every closure-ranked Sachs family failed at least one required "
            "production-validation fit. See "
            f"{summary_dir / 'family_selection_attempts.csv'} for the ordered "
            "fallback audit. No preliminary result is quoted."
        )
    #endif

    fallback_used = bool(chosen_family_rank > 1)
    with open(summary_dir / "chosen_sachs_family.txt", "w") as fout:
        fout.write(f"chosen_family={chosen_family}\n")
        fout.write(f"closure_rank={chosen_family_rank}\n")
        fout.write(f"bias_rank1_family={bias_rank1_family}\n")
        fout.write(f"fallback_used={fallback_used}\n")
        fout.write(
            "criterion=first KM15-closure-ranked family, ordered by increasing "
            "KM15 combined RMS bias-variance objective, that yields a valid "
            "nominal 5% KM15 Moradi-error production fit. VGG99/GK16 are "
            "external diagnostics only. Threshold-scan failures do not "
            "disqualify the family.\n"
        )
    #endwith

    # Preserve the extrapolation-bias estimates for the family that actually
    # survived production validation, not automatically the closure rank-1
    # family.
    bias = summarize_bias_for_family(
        bias_dirs[FINAL_NOMINAL_MODEL] / "radius_bias_variance_study.csv",
        chosen_family,
    )
    bias_table = pd.DataFrame([{
        "model": FINAL_NOMINAL_MODEL,
        "family": chosen_family,
        **bias,
    }])
    bias_table.to_csv(
        summary_dir / "chosen_family_bias_estimates.csv", index=False
    )

    km15_bias_row = bias_table.loc[
        bias_table["model"].astype(str) == FINAL_NOMINAL_MODEL
    ]
    if len(km15_bias_row) != 1:
        raise RuntimeError("Could not identify unique KM15 bias row")
    #endif
    conservative_bias = {
        "rE_bias_systematic_fm": float(
            km15_bias_row.iloc[0]["rE_RMS_bias_fm"]
        ),
        "rM_bias_systematic_fm": float(
            km15_bias_row.iloc[0]["rM_RMS_bias_fm"]
        ),
    }

    # Production fit is the already validated KM15 5% fit. Add VGG99/GK16
    # 5% fits only as external diagnostics for the comparison outputs.
    fit_rows = list(chosen_5pct_fit_rows)
    for model in FINAL_MODEL_NAMES:
        if model == FINAL_NOMINAL_MODEL:
            continue
        #endif
        result = fit_sachs_family_multi_measurements(
            selected[model][0.05],
            family=chosen_family,
            bh_cut=0.05,
            add_moradi_bh_systematic=True,
        )
        result["model"] = model
        result["threshold"] = 0.05
        result["external_diagnostic_only"] = True
        fit_rows.append(result)
    #endfor
    for result in fit_rows:
        model = str(result["model"])
        print(
            f"[final fit] {model:5s} {chosen_family}: "
            f"N={result['N']} chi2/ndf={result['chi2_ndof']:.4f} "
            f"rE={result['rE_fm']:.5f} +/- "
            f"{result['rE_fit_err_fm']:.5f} fm, "
            f"rM={result['rM_fm']:.5f} +/- "
            f"{result['rM_fit_err_fm']:.5f} fm"
        )
    #endfor
    fit_table = pd.DataFrame(fit_rows)
    fit_table.to_csv(
        summary_dir / "model_selected_05pct_fits.csv", index=False
    )

    # Independent data-driven selection robustness test. Kelly is used only as
    # a fixed BH reference; matched-N selectors avoid conflating sample size
    # with selection criterion. These shifts are diagnostics and are not folded
    # into the quoted method systematic automatically.
    bh_deviance = build_data_driven_bh_selection_diagnostics(
        physics_bundles, selected[FINAL_NOMINAL_MODEL][0.05],
        chosen_family, diagnostics_dir
    )

    # Common-support diagnostic.  This is NOT a fourth purity prescription and
    # does not replace the model-specific fits.  It asks what the chosen Sachs
    # family gives on exactly the points classified BH-like by all three models,
    # helping distinguish sample-composition effects from the nominal result.
    intersection_specs = []
    intersection_count_row = {"threshold": 0.05}
    for bundle in physics_bundles:
        d_intersection = select_bundle_from_external_intersection(
            bundle, selection, threshold=0.05
        )
        intersection_specs.append(
            bundle_to_measurement_spec(bundle, d_intersection)
        )
        intersection_count_row["N_" + str(bundle["key"])] = int(
            len(d_intersection)
        )
    #endfor
    intersection_count_row["N_total"] = int(sum(
        len(spec["data"]) for spec in intersection_specs
    ))
    pd.DataFrame([intersection_count_row]).to_csv(
        diagnostics_dir / "all_model_intersection_05pct_counts.csv", index=False
    )
    if intersection_count_row["N_total"] > 0:
        intersection_fit = fit_sachs_family_multi_measurements(
            intersection_specs,
            family=chosen_family,
            bh_cut=0.05,
            add_moradi_bh_systematic=True,
        )
        intersection_fit["sample"] = "all_model_intersection"
        intersection_fit["threshold"] = 0.05
        pd.DataFrame([intersection_fit]).to_csv(
            diagnostics_dir / "all_model_intersection_05pct_fit.csv", index=False
        )
        print(
            f"[final fit diagnostic] all-model 5% intersection {chosen_family}: "
            f"N={intersection_fit['N']} "
            f"chi2/ndf={intersection_fit['chi2_ndof']:.4f} "
            f"rE={intersection_fit['rE_fm']:.5f} fm, "
            f"rM={intersection_fit['rM_fm']:.5f} fm"
        )
    else:
        print(
            "[final fit diagnostic] all-model 5% intersection is empty; "
            "no common-support fit written"
        )
    #endif

    # Full 1--10% threshold scan for every purity model in two error modes:
    #   published_errors : no extra threshold*xsec uncertainty;
    #   moradi_errors    : includes the Moradi c_BH*xsec term.
    #
    # The threshold systematic is derived from published_errors so changing the
    # BH cut changes the selected sample but does not simultaneously and
    # mechanically change every retained point's weight.
    threshold_rows = []
    for error_mode, add_moradi in [
        ("published_errors", False),
        ("moradi_errors", True),
    ]:
        for model in FINAL_MODEL_NAMES:
            for threshold in thresholds:
                result = fit_sachs_family_multi_measurements(
                    selected[model][round(float(threshold), 2)],
                    family=chosen_family,
                    bh_cut=float(threshold),
                    add_moradi_bh_systematic=add_moradi,
                )
                result["model"] = model
                result["threshold"] = float(threshold)
                result["error_mode"] = error_mode
                threshold_rows.append(result)
            #endfor
        #endfor
    #endfor
    threshold_table = pd.DataFrame(threshold_rows)
    threshold_table.to_csv(
        summary_dir / "chosen_family_threshold_scan_01_to_10pct.csv",
        index=False,
    )

    # Production prescription:
    #   * KM15 5% is the nominal central fit.
    #   * The BH-selection systematic is the maximum absolute KM15-only
    #     excursion from that 5% baseline over thresholds 3--7%.
    #   * VGG99/GK16 are external diagnostics only and never enter the quoted
    #     production uncertainty.
    #   * Keep the complete 1--10% scans for all models as diagnostics.
    nominal = fit_table.loc[fit_table["model"] == FINAL_NOMINAL_MODEL].iloc[0]

    published_scan = threshold_table.loc[
        threshold_table["error_mode"].astype(str) == "published_errors"
    ].copy()
    nominal_scan_05 = published_scan.loc[
        (published_scan["model"].astype(str) == FINAL_NOMINAL_MODEL)
        & np.isclose(published_scan["threshold"].to_numpy(float), 0.05)
    ]
    if len(nominal_scan_05) != 1:
        raise RuntimeError(
            "Could not identify unique KM15 5% published-errors baseline"
        )
    #endif
    selection_baseline_e = float(nominal_scan_05.iloc[0]["rE_fm"])
    selection_baseline_m = float(nominal_scan_05.iloc[0]["rM_fm"])

    local_scan = published_scan.loc[
        (published_scan["model"].astype(str) == FINAL_NOMINAL_MODEL)
        & (published_scan["threshold"].to_numpy(float) >= 0.03)
        & (published_scan["threshold"].to_numpy(float) <= 0.07)
    ].copy()

    invalid_local = local_scan.loc[
        ~local_scan["valid"].astype(bool)
    ].copy()
    valid_local = local_scan.loc[
        local_scan["valid"].astype(bool)
    ].copy()

    if len(invalid_local):
        bad = ", ".join(
            f"{row.model}@{100.0 * float(row.threshold):.0f}%"
            for row in invalid_local.itertuples()
        )
        print(
            "[final model analysis] WARNING: skipping invalid KM15 3--7% "
            f"threshold-scan fits when building the production envelope: {bad}"
        )
    #endif

    if len(valid_local) == 0:
        raise RuntimeError(
            "No valid KM15 3--7% published-error threshold fits remain for "
            "the production BH-selection envelope."
        )
    #endif

    km15_selection_sys_e = float(np.max(np.abs(
        valid_local["rE_fm"].to_numpy(float) - selection_baseline_e
    )))
    km15_selection_sys_m = float(np.max(np.abs(
        valid_local["rM_fm"].to_numpy(float) - selection_baseline_m
    )))

    # Retain fixed-cut and nominal-model threshold components only as
    # diagnostics so it remains possible to understand the origin of the joint
    # envelope without double counting them in the preliminary quoted error.
    alternatives = fit_table.loc[fit_table["model"] != FINAL_NOMINAL_MODEL]
    fixed5_model_sys_e = float(np.max(np.abs(
        alternatives["rE_fm"].to_numpy(float) - nominal["rE_fm"]
    )))
    fixed5_model_sys_m = float(np.max(np.abs(
        alternatives["rM_fm"].to_numpy(float) - nominal["rM_fm"]
    )))

    km15_scan = published_scan.loc[
        published_scan["model"].astype(str) == FINAL_NOMINAL_MODEL
    ].copy()
    km15_1to10_sys_e = float(np.max(np.abs(
        km15_scan["rE_fm"].to_numpy(float) - selection_baseline_e
    )))
    km15_1to10_sys_m = float(np.max(np.abs(
        km15_scan["rM_fm"].to_numpy(float) - selection_baseline_m
    )))

    bias_sys_e = conservative_bias["rE_bias_systematic_fm"]
    bias_sys_m = conservative_bias["rM_bias_systematic_fm"]
    method_sys_e = float(np.hypot(km15_selection_sys_e, bias_sys_e))
    method_sys_m = float(np.hypot(km15_selection_sys_m, bias_sys_m))

    devfits = bh_deviance["fit_table"]
    valid_devfits = devfits.loc[devfits["valid"].astype(bool)] if len(devfits) else devfits
    data_selection_sys_e = float(np.max(np.abs(
        valid_devfits["rE_fm"].to_numpy(float) - float(nominal["rE_fm"])
    ))) if len(valid_devfits) else np.nan
    data_selection_sys_m = float(np.max(np.abs(
        valid_devfits["rM_fm"].to_numpy(float) - float(nominal["rM_fm"])
    ))) if len(valid_devfits) else np.nan

    # Compact radius-level comparison of the nominal model selector against
    # the two matched-statistics, data-driven pure-BH selectors.
    radius_compare_rows = [{
        "selector": "KM15_5pct_nominal",
        "rE_fm": float(nominal["rE_fm"]),
        "rE_err_fm": float(nominal["rE_fit_err_fm"]),
        "rM_fm": float(nominal["rM_fm"]),
        "rM_err_fm": float(nominal["rM_fit_err_fm"]),
    }]
    for row in valid_devfits.itertuples():
        radius_compare_rows.append({
            "selector": str(row.selector),
            "rE_fm": float(row.rE_fm),
            "rE_err_fm": float(row.rE_fit_err_fm),
            "rM_fm": float(row.rM_fm),
            "rM_err_fm": float(row.rM_fit_err_fm),
        })
    #endfor
    radius_compare = pd.DataFrame(radius_compare_rows)
    radius_compare.to_csv(
        diagnostics_dir / "data_driven_bh_deviance_radius_comparison.csv",
        index=False,
    )
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.8))
    xx = np.arange(len(radius_compare))
    axes[0].errorbar(xx, radius_compare["rE_fm"], yerr=radius_compare["rE_err_fm"], fmt="o")
    axes[1].errorbar(xx, radius_compare["rM_fm"], yerr=radius_compare["rM_err_fm"], fmt="o")
    labels = [str(x).replace("kelly_", "").replace("_matchedN", "") for x in radius_compare["selector"]]
    for ax in axes:
        ax.set_xticks(xx)
        ax.set_xticklabels(labels, rotation=20, ha="right")
        ax.grid(alpha=0.2)
    #endfor
    axes[0].set_ylabel(r"$r_E$ [fm]")
    axes[1].set_ylabel(r"$r_M$ [fm]")
    axes[0].set_title("Electric radius")
    axes[1].set_title("Magnetic radius")
    fig.suptitle("KM15 selection vs data-driven BH-compatible selections")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(diagnostics_dir / "data_driven_bh_deviance_radius_comparison.png", dpi=180)
    plt.close(fig)

    summary = pd.DataFrame([{
        "analysis_ensemble": "CLAS6 Jo 2015 + CLAS12 Lee 2026",
        "excluded_default_datasets": (
            "CLAS6 Saylor 2018; Hall A Defurne 2015"
        ),
        "status": "preliminary",
        "nominal_model": FINAL_NOMINAL_MODEL,
        "chosen_family": chosen_family,
        "chosen_family_closure_rank": int(chosen_family_rank),
        "bias_rank1_family": bias_rank1_family,
        "family_fallback_used": bool(fallback_used),
        "nominal_threshold": 0.05,
        "selection_systematic_window_min": 0.03,
        "selection_systematic_window_max": 0.07,
        "selection_scan_invalid_rows_skipped": int(len(invalid_local)),
        "selection_scan_valid_rows_used": int(len(valid_local)),
        "closure_min_global_valid_fraction": float(
            args.radius_bias_min_global_valid_fraction
        ),
        "closure_min_scenario_valid_fraction": float(
            args.radius_bias_min_scenario_valid_fraction
        ),
        "rE_fm": float(nominal["rE_fm"]),
        "rE_fit_err_fm": float(nominal["rE_fit_err_fm"]),
        "rE_km15_3to7_bh_selection_sys_fm": km15_selection_sys_e,
        "rE_extrapolation_bias_sys_fm": bias_sys_e,
        "rE_method_sys_fm": method_sys_e,
        "rE_fixed5_model_spread_diagnostic_fm": fixed5_model_sys_e,
        "rE_km15_1to10_threshold_envelope_diagnostic_fm": km15_1to10_sys_e,
        "rE_data_driven_bh_selection_diagnostic_fm": data_selection_sys_e,
        "rM_fm": float(nominal["rM_fm"]),
        "rM_fit_err_fm": float(nominal["rM_fit_err_fm"]),
        "rM_km15_3to7_bh_selection_sys_fm": km15_selection_sys_m,
        "rM_extrapolation_bias_sys_fm": bias_sys_m,
        "rM_method_sys_fm": method_sys_m,
        "rM_fixed5_model_spread_diagnostic_fm": fixed5_model_sys_m,
        "rM_km15_1to10_threshold_envelope_diagnostic_fm": km15_1to10_sys_m,
        "rM_data_driven_bh_selection_diagnostic_fm": data_selection_sys_m,
        "threshold_systematic_error_mode": "published_errors",
        "note": (
            "Preliminary CLAS-only result. KM15 5% is nominal. "
            "The quoted method systematic combines in quadrature the KM15-only "
            "3--7% maximum BH-selection excursion and the KM15 closure/"
            "extrapolation RMS-bias systematic. VGG99/GK16 fixed-cut spreads "
            "and their threshold scans are external diagnostics only and are "
            "not included in the production uncertainty. A matched-statistics "
            "Kelly pure-BH fractional-deviation/pull selection is also retained "
            "as an independent data-selection diagnostic and is not yet folded "
            "into the quoted method uncertainty. The KM15 1--10% threshold "
            "envelope is also retained as a diagnostic. Fit "
            "uncertainty includes "
            "the configured point-error prescription and correlated "
            "normalization nuisances."
        ),
    }])

    # Deliberately small primary summary set.
    summary.to_csv(summary_dir / "final_results.csv", index=False)
    summary.to_csv(summary_dir / "final_radius_prescription.csv", index=False)
    family_ranking.to_csv(summary_dir / "closure_summary.csv", index=False)

    # One complete model x threshold table: selected counts and fitted radii in
    # one place, instead of requiring the user to cross-reference CSVs.
    threshold_summary = threshold_table.merge(
        count_table,
        on=["model", "threshold"],
        how="left",
        validate="many_to_one",
    )
    threshold_summary.to_csv(
        summary_dir / "model_threshold_summary.csv", index=False
    )

    save_preliminary_radius_summary_figure(
        fit_table, summary, figures_dir
    )
    save_final_threshold_stability_figure(
        threshold_table, figures_dir
    )
    save_final_form_factor_comparison(
        fit_rows, figures_dir
    )
    save_final_family_ranking_figure(
        family_ranking, figures_dir
    )
    write_final_analysis_readme(
        summary=summary,
        fit_table=fit_table,
        figures_dir=figures_dir,
        summary_dir=summary_dir,
        diagnostics_dir=diagnostics_dir,
        output_path=final_dir / "README.txt",
    )

    print("\n[final model analysis] PRELIMINARY CLAS-ONLY PRESCRIPTION")
    print(
        f"  ensemble            : CLAS6 Jo 2015 + CLAS12 Lee 2026\n"
        f"  common Sachs family : {chosen_family} "
        f"(closure rank {chosen_family_rank})\n"
        f"  bias-rank-1 family  : {bias_rank1_family}\n"
        f"  fallback used       : {fallback_used}\n"
        f"  nominal purity model: {FINAL_NOMINAL_MODEL}\n"
        f"  nominal BH cut      : 5%\n"
        f"  selection sys window: KM15 only, 3--7% maximum excursion "
        f"(invalid scan fits skipped)\n"
        f"  rE = {float(nominal['rE_fm']):.5f} +/- "
        f"{float(nominal['rE_fit_err_fm']):.5f} (fit) +/- "
        f"{method_sys_e:.5f} (method) fm\n"
        f"       KM15 selection={km15_selection_sys_e:.5f}, "
        f"closure bias={bias_sys_e:.5f}\n"
        f"  rM = {float(nominal['rM_fm']):.5f} +/- "
        f"{float(nominal['rM_fit_err_fm']):.5f} (fit) +/- "
        f"{method_sys_m:.5f} (method) fm\n"
        f"       KM15 selection={km15_selection_sys_m:.5f}, "
        f"closure bias={bias_sys_m:.5f}\n"
        f"  summary -> {summary_dir / 'final_results.csv'}\n"
        f"  figures -> {figures_dir}\n"
        f"  diagnostics -> {diagnostics_dir}"
    )

#enddef



def _fit5_prediction_and_pull_for_bundle(
        bundle: Dict[str, object],
        bh_cut: float = 0.05) -> pd.DataFrame:
    """Return Fit-5 prediction, residual, pull, and chi2 for one selected bundle."""
    data = bundle["set5"].copy().reset_index(drop=True)
    fit5 = next(r for r in bundle["results"] if r.name == "Fit 5")
    q = data["t_abs"].to_numpy(float)
    f1, f2 = paper_model_f1f2(fit5.model_kind, q, fit5.params)
    pred = bh_from_f1f2(
        data["bh_A"].to_numpy(float),
        data["bh_B"].to_numpy(float),
        data["bh_C"].to_numpy(float),
        f1, f2,
    )
    beta_key = "beta_" + re.sub(
        r"[^A-Za-z0-9]+", "_", str(bundle["key"])
    )
    pred *= float(fit5.meta.get(beta_key + "_scale_factor", 1.0))
    err = dataset_point_errors(
        data, str(bundle["kind"]), float(bh_cut), True
    )
    y = data["xs"].to_numpy(float)
    out = data.copy()
    out["fit5_prediction"] = pred
    out["fit5_error"] = err
    out["fit5_residual"] = y - pred
    out["fit5_pull"] = (y - pred) / err
    out["fit5_chi2_contribution"] = out["fit5_pull"].to_numpy(float) ** 2
    return out
#enddef


def save_jo_saylor_direct_comparison(
        jo_bundle: Dict[str, object],
        saylor_bundle: Dict[str, object],
        outdir: Path) -> None:
    """
    Directly compare why Jo 2015 and Saylor 2018 behave differently in pure-BH fits.

    The comparison uses each experiment's own nominal 5% KM15-selected sample and
    its own Fit-5 prediction.  It therefore does not require identical kinematics.
    Binned mean pulls and chi2 per point expose whether the Saylor discrepancy is
    localized in Q2, xB, |t|, or phi rather than merely reflecting a global offset.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    prepared = {
        "CLAS6 Jo 2015": _fit5_prediction_and_pull_for_bundle(jo_bundle),
        "CLAS6 Saylor 2018": _fit5_prediction_and_pull_for_bundle(saylor_bundle),
    }
    variables = [
        ("Q2", r"$Q^2$ (GeV$^2$)"),
        ("xB", r"$x_B$"),
        ("t_abs", r"$|t|$ (GeV$^2$)"),
        ("phi_deg", r"$\phi$ (deg)"),
    ]

    rows = []
    fig, axes = plt.subplots(2, 4, figsize=(17.0, 8.5), sharex="col")
    for icol, (var, xlabel) in enumerate(variables):
        allvals = np.concatenate([
            d[var].to_numpy(float) for d in prepared.values()
        ])
        finite_all = allvals[np.isfinite(allvals)]
        if len(finite_all) < 2:
            continue
        #endif
        edges = np.linspace(
            float(np.nanmin(finite_all)),
            float(np.nanmax(finite_all)),
            11,
        )
        centers = 0.5 * (edges[:-1] + edges[1:])

        for label, data in prepared.items():
            x = data[var].to_numpy(float)
            pull = data["fit5_pull"].to_numpy(float)
            chi2_i = data["fit5_chi2_contribution"].to_numpy(float)
            mean_pull = np.full(len(centers), np.nan)
            chi2_per_point = np.full(len(centers), np.nan)

            for ibin, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
                mask = (
                    np.isfinite(x)
                    & np.isfinite(pull)
                    & (x >= lo)
                    & ((x < hi) if ibin < len(edges) - 2 else (x <= hi))
                )
                if not np.any(mask):
                    continue
                #endif
                mean_pull[ibin] = float(np.mean(pull[mask]))
                chi2_per_point[ibin] = float(np.mean(chi2_i[mask]))
                rows.append({
                    "dataset": label,
                    "variable": var,
                    "bin": int(ibin),
                    "low": float(lo),
                    "high": float(hi),
                    "N": int(np.sum(mask)),
                    "mean_pull": float(mean_pull[ibin]),
                    "rms_pull": float(np.sqrt(np.mean(pull[mask] ** 2))),
                    "chi2_sum": float(np.sum(chi2_i[mask])),
                    "chi2_per_point": float(chi2_per_point[ibin]),
                })
            #endfor

            axes[0, icol].plot(
                centers, mean_pull, marker="o", linewidth=1.2, label=label
            )
            axes[1, icol].plot(
                centers, chi2_per_point, marker="o", linewidth=1.2, label=label
            )
        #endfor

        axes[0, icol].axhline(0.0, linewidth=0.8)
        axes[0, icol].set_title(xlabel)
        axes[1, icol].set_xlabel(xlabel)
        axes[0, icol].grid(alpha=0.2)
        axes[1, icol].grid(alpha=0.2)
    #endfor

    axes[0, 0].set_ylabel("Binned mean Fit-5 pull")
    axes[1, 0].set_ylabel(r"Binned $\chi^2$/point")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles, labels, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 0.965)
    )
    fig.suptitle(
        "CLAS6 Jo 2015 vs Saylor 2018: where the nominal pure-BH fit disagrees",
        y=0.995,
    )
    fig.subplots_adjust(
        top=0.89, bottom=0.10, left=0.07, right=0.985,
        hspace=0.22, wspace=0.24,
    )
    fig.savefig(outdir / "01_jo2015_vs_saylor2018_fit5_discrepancy.png", dpi=180)
    plt.close(fig)

    pd.DataFrame(rows).to_csv(
        outdir / "jo2015_vs_saylor2018_binned_fit5_diagnostics.csv",
        index=False,
    )

    summary_rows = []
    for label, data in prepared.items():
        pull = data["fit5_pull"].to_numpy(float)
        finite = np.isfinite(pull)
        summary_rows.append({
            "dataset": label,
            "N": int(np.sum(finite)),
            "chi2_sum": float(np.sum(pull[finite] ** 2)),
            "chi2_per_point": float(np.mean(pull[finite] ** 2)),
            "mean_pull": float(np.mean(pull[finite])),
            "rms_pull": float(np.sqrt(np.mean(pull[finite] ** 2))),
            "fraction_abs_pull_gt_2": float(np.mean(np.abs(pull[finite]) > 2.0)),
            "fraction_abs_pull_gt_3": float(np.mean(np.abs(pull[finite]) > 3.0)),
        })
    #endfor
    pd.DataFrame(summary_rows).to_csv(
        outdir / "jo2015_vs_saylor2018_fit5_summary.csv",
        index=False,
    )
#enddef


def _profile_model_normalization(
        measured: np.ndarray,
        predicted: np.ndarray,
        error: np.ndarray,
        norm_frac: float) -> Tuple[float, float, float]:
    """Profile one Gaussian correlated normalization nuisance analytically."""
    y = np.asarray(measured, dtype=float)
    m = np.asarray(predicted, dtype=float)
    e = np.asarray(error, dtype=float)
    finite = (
        np.isfinite(y) & np.isfinite(m) & np.isfinite(e)
        & (e > 0.0) & (m > 0.0)
    )
    y = y[finite]
    m = m[finite]
    e = e[finite]
    if len(y) == 0:
        return np.nan, np.nan, np.nan
    #endif

    frac = float(norm_frac)
    if frac > 0.0:
        w = 1.0 / (e ** 2)
        denom = 1.0 + frac * frac * float(np.sum(m * m * w))
        beta = (
            frac * float(np.sum(m * (y - m) * w)) / denom
        )
    else:
        beta = 0.0
    #endif
    scale = 1.0 + frac * beta
    pull = (scale * m - y) / e
    chi2 = float(np.dot(pull, pull) + beta * beta)
    return beta, scale, chi2
#enddef



def make_pass2_diagnostic_bundle(args) -> Dict[str, object]:
    """
    Build a diagnostic-only CLAS12 pass-2 bundle for all-point model scoring.

    This does not enter the production Jo+Lee form-factor fit.  It exists only
    so the published-point model-agreement comparison can include the current
    pass-2 cross sections alongside the historical measurements.
    """
    csv_path = Path(args.csv).expanduser().resolve()
    df = load_clas12_csv(csv_path)
    outdir = (
        Path(args.outdir).expanduser().resolve()
        / "clas12_pass2_model_diagnostic"
    )
    outdir.mkdir(parents=True, exist_ok=True)

    # Pass-2 uses the same bin boundaries as pass-1, but its stored average
    # kinematics can differ.  Evaluate KM15 at the actual pass-2 averages and
    # export these points for independent VGG99/GK16 PARTONS evaluation.
    df = evaluate_km15_dataframe(
        df,
        float(args.ebeam),
        max(1, int(args.workers)),
        outdir / "km15_bh_decomposition_pass2_diagnostic.csv",
        bool(args.force_km15),
    )
    finite = (
        np.isfinite(df["km15_ep"])
        & np.isfinite(df["km15_bh"])
        & np.isfinite(df["R_BH"])
        & (df["km15_ep"] > 0.0)
        & (df["km15_bh"] > 0.0)
    )
    df = df.loc[finite].copy().reset_index(drop=True)

    # The pass-2 CSV stores the measured average kinematics but not the
    # run-wide beam energy.  Add it explicitly so these points can be exported
    # to the external KM15/VGG99/GK16 evaluator just like the other datasets.
    if "ebeam" not in df.columns:
        df["ebeam"] = float(args.ebeam)
    #endif

    scale = pd.to_numeric(
        df.get("scale_sys_frac", pd.Series(dtype=float)),
        errors="coerce",
    ).to_numpy(float)
    scale = scale[np.isfinite(scale) & (scale > 0.0)]
    norm_frac = float(np.median(scale)) if len(scale) else 0.0

    print(
        f"[PASS2 diagnostic] Loaded {len(df)} current pass-2 points; "
        f"representative correlated scale={100.0*norm_frac:.2f}%."
    )
    print(
        "[PASS2 diagnostic] KM15 evaluated at pass-2 average kinematics; "
        "VGG99/GK16 will use the exported pass-2 points."
    )
    return {
        "label": "CLAS12 Hayward (unpublished)",
        "key": "pass2",
        "kind": "pass2",
        "norm_frac": norm_frac,
        "results": [],
        "set5": pd.DataFrame(),
        "all_data": df.copy(),
        "outdir": outdir,
    }
#enddef



def match_pass2_to_pass1_model_rows(
        pass2_data: pd.DataFrame,
        selection: pd.DataFrame,
        atol: float = 5.0e-9) -> pd.DataFrame:
    """
    Reuse pass-1 model predictions for pass-2 after an explicit kinematic match.

    The pass-1 and pass-2 analyses use the same CLAS12 binning.  Matching by
    kinematics rather than source_row also handles the case where one analysis
    contains one fewer valid cross-section point.
    """
    pass1 = selection.loc[
        selection["dataset"].astype(str) == "pass1"
    ].copy().reset_index(drop=True)
    if len(pass1) == 0:
        raise RuntimeError(
            "Cannot reuse pass-1 predictions for pass-2: the external model "
            "selection table contains no pass1 rows."
        )
    #endif

    kin_cols = ["xB", "Q2", "t_abs", "phi_deg", "ebeam"]
    missing_ext = [c for c in kin_cols if c not in pass1.columns]
    if missing_ext:
        raise KeyError(
            "Pass-2/pass-1 model reuse requires the external pass-1 table to "
            f"contain {kin_cols}; missing external={missing_ext}"
        )
    #endif

    # The pass-2 analysis CSV does not need to carry Ebeam because the beam
    # energy is a run-wide analysis setting.  Infer it from the already
    # evaluated pass-1 rows, which use the identical CLAS12 kinematic grid.
    pass2_data = pass2_data.copy()
    if "ebeam" not in pass2_data.columns:
        pass1_ebeams = np.unique(
            np.round(pass1["ebeam"].to_numpy(float), 10)
        )
        pass1_ebeams = pass1_ebeams[np.isfinite(pass1_ebeams)]
        if len(pass1_ebeams) != 1:
            raise RuntimeError(
                "Cannot infer a unique CLAS12 beam energy for pass-2 reuse; "
                f"pass-1 external table contains {pass1_ebeams.tolist()}."
            )
        #endif
        pass2_data["ebeam"] = float(pass1_ebeams[0])
        print(
            "[PASS2 model reuse] pass-2 CSV has no ebeam column; "
            f"using pass-1 Ebeam={float(pass1_ebeams[0]):.6f} GeV."
        )
    #endif

    missing_data = [c for c in kin_cols if c not in pass2_data.columns]
    if missing_data:
        raise KeyError(
            "Pass-2/pass-1 model reuse requires kinematic columns "
            f"{kin_cols}; missing pass2={missing_data}"
        )
    #endif

    # Build a deterministic nearest-exact lookup.  Rounding is only a hash
    # accelerator; every accepted match is subsequently checked numerically.
    def _key(row) -> tuple:
        return tuple(round(float(row[c]), 8) for c in kin_cols)
    #enddef

    lookup = {}
    for _, row in pass1.iterrows():
        key = _key(row)
        lookup.setdefault(key, []).append(row)
    #endfor

    matched_rows = []
    maxdiff = {c: 0.0 for c in kin_cols}
    unmatched = []
    ambiguous = []
    for irow, (_, drow) in enumerate(pass2_data.iterrows()):
        key = _key(drow)
        candidates = lookup.get(key, [])
        if len(candidates) != 1:
            if len(candidates) == 0:
                unmatched.append(irow)
            else:
                ambiguous.append(irow)
            #endif
            continue
        #endif

        erow = candidates[0].copy()
        diffs = {
            c: abs(float(erow[c]) - float(drow[c]))
            for c in kin_cols
        }
        if any(diff > atol for diff in diffs.values()):
            unmatched.append(irow)
            continue
        #endif
        for c, diff in diffs.items():
            maxdiff[c] = max(maxdiff[c], diff)
        #endfor

        erow["dataset"] = "pass2"
        erow["point_id"] = f"pass2:{irow}"
        erow["source_row"] = irow
        matched_rows.append(erow)
    #endfor

    if unmatched or ambiguous or len(matched_rows) != len(pass2_data):
        raise RuntimeError(
            "Pass-2 model reuse failed: "
            f"matched {len(matched_rows)}/{len(pass2_data)}, "
            f"unmatched={len(unmatched)}, ambiguous={len(ambiguous)}. "
            "The pass-1/pass-2 kinematic grids are not identical under the "
            "current matching tolerance."
        )
    #endif

    print(
        f"[PASS2 model reuse] {len(matched_rows)}/{len(pass2_data)} pass-2 "
        "points matched to pass-1 kinematics."
    )
    print(
        "[PASS2 model reuse] max |differences|: "
        + ", ".join(f"{c}={maxdiff[c]:.3e}" for c in kin_cols)
    )
    print(
        "[PASS2 model reuse] Reusing KM15/VGG99/GK16 predictions; "
        "no model reevaluation required."
    )
    return pd.DataFrame(matched_rows).reset_index(drop=True)
#enddef


def save_all_point_model_agreement_diagnostics(
        selection: pd.DataFrame,
        bundles: Sequence[Dict[str, object]],
        outdir: Path) -> None:
    """
    Compare every published point with the full EP prediction of all three models.

    This is a diagnostic model-validation score only.  It deliberately includes
    the central-phi/DVCS-sensitive region and is NEVER fed back into BH-purity
    selection, family ranking, or the quoted model systematic.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    selection_datasets = set(selection["dataset"].astype(str))
    usable = sort_bundles_chronologically([
        b for b in bundles if str(b["key"]) in selection_datasets
    ])
    missing = [
        str(b["key"]) for b in bundles
        if str(b["key"]) not in selection_datasets
    ]
    if missing:
        print(
            "[model agreement diagnostic] external model table does not contain "
            + ", ".join(missing)
            + "; those datasets are omitted. Use the exported all-enabled "
              "diagnostic kinematics file to evaluate them with KM15/VGG99/GK16."
        )
    #endif
    if not usable:
        print("[model agreement diagnostic] no enabled dataset has external model predictions")
        return
    #endif

    score_rows = []
    binned_rows = []
    variables = [
        ("Q2", r"$Q^2$ (GeV$^2$)"),
        ("xB", r"$x_B$"),
        ("t_abs", r"$|t|$ (GeV$^2$)"),
        ("phi_deg", r"$\phi$ (deg)"),
    ]

    for bundle in usable:
        key = str(bundle["key"])
        label = str(bundle["label"])
        kind = str(bundle["kind"])
        data = bundle["all_data"].reset_index(drop=True)
        ext = selection.loc[
            selection["dataset"].astype(str) == key
        ].copy().sort_values("source_row").reset_index(drop=True)
        if len(ext) != len(data):
            print(
                f"[model agreement diagnostic] {label}: external/current point "
                f"count mismatch ({len(ext)} vs {len(data)}); skipping"
            )
            continue
        #endif
        if not np.array_equal(
                ext["source_row"].to_numpy(int),
                np.arange(len(data), dtype=int)):
            print(
                f"[model agreement diagnostic] {label}: source_row mismatch; skipping"
            )
            continue
        #endif

        y = data["xs"].to_numpy(float)
        stat_err = dataset_statistical_errors(data, kind)
        err = dataset_point_errors(data, kind, 0.0, False)
        norm_frac = float(bundle.get("norm_frac", 0.0))

        precision_finite = (
            np.isfinite(y) & (y > 0.0)
            & np.isfinite(stat_err) & (stat_err > 0.0)
            & np.isfinite(err) & (err > 0.0)
        )
        stat_frac = np.full(len(data), np.nan)
        total_frac = np.full(len(data), np.nan)
        stat_frac[precision_finite] = stat_err[precision_finite] / y[precision_finite]
        total_frac[precision_finite] = err[precision_finite] / y[precision_finite]

        def _pct_summary(arr: np.ndarray) -> Tuple[float, float, float]:
            vals = np.asarray(arr, dtype=float)
            vals = vals[np.isfinite(vals)]
            if len(vals) == 0:
                return np.nan, np.nan, np.nan
            #endif
            q16, q50, q84 = np.percentile(vals, [16.0, 50.0, 84.0])
            return float(q16), float(q50), float(q84)
        #enddef

        stat16, stat50, stat84 = _pct_summary(stat_frac)
        total16, total50, total84 = _pct_summary(total_frac)

        fig, axes = plt.subplots(2, 2, figsize=(12.8, 9.2))
        for model in FINAL_MODEL_NAMES:
            ep_col = MODEL_EP_COLUMN[model]
            if ep_col not in ext.columns:
                print(
                    f"[model agreement diagnostic] {label}: missing {ep_col}; "
                    f"cannot score {model}"
                )
                continue
            #endif
            pred = pd.to_numeric(ext[ep_col], errors="coerce").to_numpy(float)
            finite = (
                np.isfinite(y) & np.isfinite(err) & (err > 0.0)
                & np.isfinite(pred) & (pred > 0.0)
            )
            if not np.any(finite):
                continue
            #endif

            raw_pull = (pred[finite] - y[finite]) / err[finite]
            raw_chi2 = float(np.dot(raw_pull, raw_pull))
            beta, scale, prof_chi2 = _profile_model_normalization(
                y[finite], pred[finite], err[finite], norm_frac
            )
            prof_pull_full = np.full(len(data), np.nan)
            prof_pull_full[finite] = (
                scale * pred[finite] - y[finite]
            ) / err[finite]
            ratio = np.full(len(data), np.nan)
            ratio[finite] = y[finite] / pred[finite]
            frac_residual = np.full(len(data), np.nan)
            frac_residual[finite] = (
                scale * pred[finite] - y[finite]
            ) / y[finite]
            abs_frac_residual = np.abs(frac_residual)

            # Central phi is intentionally retained in the all-point score.
            phi = np.mod(data["phi_deg"].to_numpy(float), 360.0)
            central = finite & (np.abs(((phi - 180.0 + 180.0) % 360.0) - 180.0) <= 60.0)
            central_chi2 = (
                float(np.sum(prof_pull_full[central] ** 2))
                if np.any(central) else np.nan
            )

            score_rows.append({
                "dataset": key,
                "dataset_label": label,
                "model": model,
                "model_display": MODEL_DISPLAY[model],
                "N_all": int(np.sum(finite)),
                "raw_chi2": raw_chi2,
                "raw_chi2_per_point": float(raw_chi2 / np.sum(finite)),
                "profiled_beta": float(beta),
                "profiled_scale": float(scale),
                "profiled_chi2": float(prof_chi2),
                "profiled_chi2_per_point": float(prof_chi2 / np.sum(finite)),
                "rms_profiled_pull": float(np.sqrt(np.mean(prof_pull_full[finite] ** 2))),
                "median_data_over_model": float(np.nanmedian(ratio[finite])),
                "median_abs_fractional_model_residual": float(
                    np.nanmedian(abs_frac_residual[finite])
                ),
                "rms_fractional_model_residual": float(
                    np.sqrt(np.nanmean(frac_residual[finite] ** 2))
                ),
                "median_stat_fraction": stat50,
                "stat_fraction_p16": stat16,
                "stat_fraction_p84": stat84,
                "median_pointwise_total_fraction": total50,
                "pointwise_total_fraction_p16": total16,
                "pointwise_total_fraction_p84": total84,
                "N_central_phi_pm60": int(np.sum(central)),
                "central_phi_chi2_per_point": (
                    float(central_chi2 / np.sum(central))
                    if np.any(central) else np.nan
                ),
                "normalization_fraction": norm_frac,
            })

            for ax, (var, xlabel) in zip(axes.flat, variables):
                x = data[var].to_numpy(float)
                f = finite & np.isfinite(x)
                if np.sum(f) < 4:
                    continue
                #endif
                edges = np.linspace(
                    float(np.nanmin(x[f])),
                    float(np.nanmax(x[f])),
                    13,
                )
                centers = 0.5 * (edges[:-1] + edges[1:])
                means = np.full(len(centers), np.nan)
                rms = np.full(len(centers), np.nan)
                for ibin, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
                    mask = (
                        f & (x >= lo)
                        & ((x < hi) if ibin < len(edges) - 2 else (x <= hi))
                    )
                    if not np.any(mask):
                        continue
                    #endif
                    means[ibin] = float(np.mean(prof_pull_full[mask]))
                    rms[ibin] = float(np.sqrt(np.mean(prof_pull_full[mask] ** 2)))
                    binned_rows.append({
                        "dataset": key,
                        "dataset_label": label,
                        "model": model,
                        "variable": var,
                        "bin": int(ibin),
                        "low": float(lo),
                        "high": float(hi),
                        "N": int(np.sum(mask)),
                        "mean_profiled_pull": float(means[ibin]),
                        "rms_profiled_pull": float(rms[ibin]),
                    })
                #endfor
                ax.plot(
                    centers, means, marker="o", linewidth=1.1,
                    label=MODEL_DISPLAY[model],
                )
                ax.axhline(0.0, linewidth=0.8)
                ax.set_xlabel(xlabel)
                ax.set_ylabel("Binned mean pull")
                ax.grid(alpha=0.2)
            #endfor
        #endfor

        handles, labels = axes[0, 0].get_legend_handles_labels()
        if handles:
            fig.legend(
                handles, labels, loc="upper center", ncol=3,
                bbox_to_anchor=(0.5, 0.965),
            )
        #endif
        fig.suptitle(
            f"{label}: all published points vs full electroproduction models",
            y=0.995,
        )
        fig.subplots_adjust(
            top=0.90, bottom=0.08, left=0.08, right=0.98,
            hspace=0.30, wspace=0.24,
        )
        safe_key = re.sub(r"[^A-Za-z0-9]+", "_", key)
        fig.savefig(
            outdir / f"all_points_model_pulls_{safe_key}.png",
            dpi=180,
        )
        plt.close(fig)
    #endfor

    scores = pd.DataFrame(score_rows)
    scores.to_csv(outdir / "all_point_model_agreement_scores.csv", index=False)

    # Compact note-ready model comparison.  Keep both profiled chi2/point and
    # the fitted scale factor so a reader can distinguish shape disagreement
    # from a simple overall normalization shift.
    if len(scores):
        note_rows = []
        for dataset_label in list(dict.fromkeys(scores["dataset_label"].astype(str))):
            row = {"dataset": dataset_label}
            ds = scores.loc[scores["dataset_label"].astype(str) == dataset_label]
            if len(ds):
                row["N"] = int(ds.iloc[0]["N_all"])
                row["median_pointwise_uncertainty_pct"] = (
                    100.0 * float(ds.iloc[0]["median_pointwise_total_fraction"])
                )
            #endif
            for model in FINAL_MODEL_NAMES:
                mr = ds.loc[ds["model"].astype(str) == model]
                prefix = re.sub(r"[^A-Za-z0-9]+", "_", MODEL_DISPLAY[model]).strip("_")
                if len(mr):
                    row[prefix + "_chi2_per_point"] = float(mr.iloc[0]["profiled_chi2_per_point"])
                    row[prefix + "_scale"] = float(mr.iloc[0]["profiled_scale"])
                    row[prefix + "_median_abs_residual_pct"] = (
                        100.0 * float(mr.iloc[0]["median_abs_fractional_model_residual"])
                    )
                #endif
            #endfor
            note_rows.append(row)
        #endfor
        pd.DataFrame(note_rows).to_csv(
            outdir / "all_point_model_agreement_note_table.csv", index=False
        )
    #endif

    pd.DataFrame(binned_rows).to_csv(
        outdir / "all_point_model_agreement_binned_pulls.csv",
        index=False,
    )

    if len(scores):
        # Enforce the same chronology used everywhere else in publication plots.
        label_to_key = {
            str(b["label"]): str(b["key"]) for b in usable
        }
        datasets = sorted(
            list(dict.fromkeys(scores["dataset_label"].astype(str))),
            key=lambda label: DATASET_CHRONOLOGY.get(
                label_to_key.get(str(label), ""), 10_000
            ),
        )
        x = np.arange(len(datasets), dtype=float)
        width = 0.24

        # Dynamic range is large for rejected models, so chi2/N is shown on a
        # logarithmic axis.  The lower panel shows the actual fractional
        # cross-section disagreement, which separates model mismatch from the
        # effect of experimental precision.
        fig, axes = plt.subplots(
            2, 1, figsize=(12.8, 9.4), sharex=True,
            gridspec_kw={"height_ratios": [1.0, 1.0]},
        )
        ax_chi2, ax_frac = axes

        for imodel, model in enumerate(FINAL_MODEL_NAMES):
            chi2_vals = []
            frac_vals = []
            for label in datasets:
                row = scores.loc[
                    (scores["dataset_label"].astype(str) == label)
                    & (scores["model"].astype(str) == model)
                ]
                if len(row):
                    chi2_vals.append(
                        float(row.iloc[0]["profiled_chi2_per_point"])
                    )
                    frac_vals.append(
                        100.0 * float(
                            row.iloc[0]["median_abs_fractional_model_residual"]
                        )
                    )
                else:
                    chi2_vals.append(np.nan)
                    frac_vals.append(np.nan)
                #endif
            #endfor

            ax_chi2.bar(
                x + (imodel - 1) * width,
                chi2_vals,
                width=width,
                label=MODEL_DISPLAY[model],
            )
            ax_frac.bar(
                x + (imodel - 1) * width,
                frac_vals,
                width=width,
                label=MODEL_DISPLAY[model],
            )
        #endfor

        # Experimental pointwise precision is a useful scale for interpreting
        # the fractional model residuals.
        precision_vals = []
        for label in datasets:
            row = scores.loc[
                scores["dataset_label"].astype(str) == label
            ]
            precision_vals.append(
                100.0 * float(row.iloc[0]["median_pointwise_total_fraction"])
                if len(row) else np.nan
            )
        #endfor
        ax_frac.plot(
            x,
            precision_vals,
            marker="D",
            linestyle="none",
            markersize=6,
            label="Median pointwise experimental uncertainty",
        )

        ax_chi2.set_yscale("log")
        ax_chi2.set_ylabel(r"Profiled $\chi^2$/point")
        ax_chi2.grid(axis="y", alpha=0.2, which="both")
        ax_chi2.legend(
            loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.20)
        )
        ax_chi2.set_title(
            "Full electroproduction model agreement with published data"
        )

        ax_frac.set_ylabel("Median absolute fractional model residual (%)")
        ax_frac.grid(axis="y", alpha=0.2)
        ax_frac.legend(
            loc="upper center", ncol=4, bbox_to_anchor=(0.5, 1.14)
        )
        ax_frac.set_xticks(x)
        ax_frac.set_xticklabels(datasets, rotation=20, ha="right")

        fig.subplots_adjust(
            top=0.88, bottom=0.16, left=0.09, right=0.98, hspace=0.28
        )
        fig.savefig(outdir / "07_all_point_model_agreement_chi2.png", dpi=300)
        plt.close(fig)

        # Compact dataset-level table combining experimental precision with
        # model agreement.  Precision columns are properties of the data and
        # are therefore repeated only once per dataset in the rendered table.
        summary_rows = []
        for label in datasets:
            drows = scores.loc[
                scores["dataset_label"].astype(str) == label
            ]
            if len(drows) == 0:
                continue
            #endif
            first = drows.iloc[0]
            row = {
                "dataset": label,
                "N": int(first["N_all"]),
                "median_stat_pct": 100.0 * float(first["median_stat_fraction"]),
                "stat_p16_pct": 100.0 * float(first["stat_fraction_p16"]),
                "stat_p84_pct": 100.0 * float(first["stat_fraction_p84"]),
                "median_pointwise_total_pct": (
                    100.0 * float(first["median_pointwise_total_fraction"])
                ),
                "pointwise_total_p16_pct": (
                    100.0 * float(first["pointwise_total_fraction_p16"])
                ),
                "pointwise_total_p84_pct": (
                    100.0 * float(first["pointwise_total_fraction_p84"])
                ),
            }
            for model in FINAL_MODEL_NAMES:
                mr = drows.loc[drows["model"].astype(str) == model]
                if len(mr):
                    rr = mr.iloc[0]
                    tag = MODEL_DISPLAY[model]
                    row[f"{tag}_chi2_per_point"] = float(
                        rr["profiled_chi2_per_point"]
                    )
                    row[f"{tag}_median_abs_frac_residual_pct"] = (
                        100.0 * float(rr["median_abs_fractional_model_residual"])
                    )
                    row[f"{tag}_rms_frac_residual_pct"] = (
                        100.0 * float(rr["rms_fractional_model_residual"])
                    )
                    row[f"{tag}_central_phi_chi2_per_point"] = float(
                        rr["central_phi_chi2_per_point"]
                    )
                #endif
            #endfor
            summary_rows.append(row)
        #endfor

        summary = pd.DataFrame(summary_rows)
        if len(summary):
            summary["_chronology"] = summary["dataset"].map(
                lambda label: DATASET_CHRONOLOGY.get(
                    label_to_key.get(str(label), ""), 10_000
                )
            )
            summary = (
                summary.sort_values("_chronology")
                .drop(columns=["_chronology"])
                .reset_index(drop=True)
            )
        #endif
        summary.to_csv(
            outdir / "08_all_point_model_agreement_with_precision.csv",
            index=False,
        )

        if len(summary):
            display_cols = [
                "Dataset", "N", "Stat. precision", "Pointwise total",
                "KM15 chi2/N", "VGG99 chi2/N", "GK16 chi2/N",
                "KM15 |frac.|", "VGG99 |frac.|", "GK16 |frac.|",
            ]
            table_rows = []
            for _, r in summary.iterrows():
                table_rows.append([
                    str(r["dataset"]),
                    f'{int(r["N"])}',
                    (
                        f'{r["median_stat_pct"]:.1f}% '
                        f'[{r["stat_p16_pct"]:.1f},{r["stat_p84_pct"]:.1f}]'
                    ),
                    (
                        f'{r["median_pointwise_total_pct"]:.1f}% '
                        f'[{r["pointwise_total_p16_pct"]:.1f},'
                        f'{r["pointwise_total_p84_pct"]:.1f}]'
                    ),
                    f'{r.get("KM15_chi2_per_point", np.nan):.2f}',
                    f'{r.get("VGG99_chi2_per_point", np.nan):.2f}',
                    f'{r.get("GK16_chi2_per_point", np.nan):.2f}',
                    f'{r.get("KM15_median_abs_frac_residual_pct", np.nan):.1f}%',
                    f'{r.get("VGG99_median_abs_frac_residual_pct", np.nan):.1f}%',
                    f'{r.get("GK16_median_abs_frac_residual_pct", np.nan):.1f}%',
                ])
            #endfor

            fig_h = max(4.2, 0.55 * len(table_rows) + 2.3)
            fig, ax = plt.subplots(figsize=(18.0, fig_h))
            ax.axis("off")
            tbl = ax.table(
                cellText=table_rows,
                colLabels=display_cols,
                loc="center",
                cellLoc="center",
            )
            tbl.auto_set_font_size(False)
            tbl.set_fontsize(9.5)
            tbl.scale(1.0, 1.55)
            for (irow, icol), cell in tbl.get_celld().items():
                if irow == 0:
                    cell.set_text_props(weight="bold")
                #endif
                if icol == 0:
                    cell.set_text_props(ha="left")
                #endif
            #endfor
            ax.set_title(
                "Experimental precision and full-electroproduction model agreement\n"
                "median [16th,84th] relative experimental errors; "
                "median absolute fractional model residual",
                pad=18.0,
            )
            fig.tight_layout()
            fig.savefig(
                outdir / "08_all_point_model_agreement_with_precision.png",
                dpi=190,
                bbox_inches="tight",
            )
            plt.close(fig)

            # A second diagnostic isolates model-data fractional disagreement
            # from experimental precision.  Unlike chi2, this quantity is not
            # amplified simply because an experiment has smaller error bars.
            x = np.arange(len(summary), dtype=float)
            width = 0.24
            fig, ax = plt.subplots(figsize=(12.0, 6.8))
            for imodel, model in enumerate(FINAL_MODEL_NAMES):
                tag = MODEL_DISPLAY[model]
                vals = summary[
                    f"{tag}_median_abs_frac_residual_pct"
                ].to_numpy(float)
                ax.bar(
                    x + (imodel - 1) * width,
                    vals,
                    width=width,
                    label=tag,
                )
            #endfor
            ax.set_xticks(x)
            ax.set_xticklabels(
                summary["dataset"].astype(str), rotation=20, ha="right"
            )
            ax.set_ylabel("Median absolute fractional model residual (%)")
            ax.grid(axis="y", alpha=0.2)
            ax.legend(loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.12))
            ax.set_title(
                "Model-data discrepancy independent of experimental error size"
            )
            fig.tight_layout()
            fig.savefig(
                outdir / "09_all_point_model_fractional_residual.png",
                dpi=180,
            )
            plt.close(fig)
        #endif
    #endif
#enddef




def _km15_selected_specs_for_bundles(
        bundles: Sequence[Dict[str, object]],
        selection: pd.DataFrame,
        threshold: float) -> Tuple[List[Dict[str, object]], Dict[str, int]]:
    """Build measurement specs for a source-defined KM15 BH-purity selection."""
    specs = []
    counts = {}
    for bundle in bundles:
        selected = select_bundle_from_external_model(
            bundle, selection, "km15", float(threshold)
        )
        specs.append(bundle_to_measurement_spec(bundle, selected))
        counts[str(bundle["key"])] = int(len(selected))
    #endfor
    return specs, counts
#enddef


def run_km15_ensemble_rotation(
        bundles: Sequence[Dict[str, object]],
        args,
        root_outdir: Path,
        family: str = "IP2") -> None:
    """
    Compare the requested KM15-only production ensembles.

    The three primary configurations are:
      1) CLAS12 Lee 2026 alone;
      2) CLAS6 Jo 2015 + CLAS12 Lee 2026;
      3) Jo + Hall A Defurne + Hall A Georges + Lee.

    Each measurement retains its own published uncertainty treatment and its
    own correlated normalization nuisance where an authoritative decomposition
    is available.  Georges currently uses its published total systematic
    pointwise and no additional normalization nuisance.
    """
    selection = load_external_bh_model_selection(
        Path(args.bh_model_selection_results).expanduser().resolve()
    )
    by_key = {str(b["key"]): b for b in bundles}
    required = ["jo2015", "halla_defurne2015", "halla_georges2022", "pass1"]
    missing = [k for k in required if k not in by_key]
    if missing:
        raise RuntimeError(
            "KM15 ensemble rotation requires Jo, Defurne, Georges, and Lee. "
            "Missing: " + ", ".join(missing)
        )
    #endif

    outdir = root_outdir / "final_analysis" / "km15_ensemble_rotation"
    outdir.mkdir(parents=True, exist_ok=True)

    configurations = [
        ("lee2026_only", [by_key["pass1"]], "CLAS12 Lee 2026"),
        (
            "jo2015_plus_lee2026",
            [by_key["jo2015"], by_key["pass1"]],
            "CLAS6 Jo 2015 + CLAS12 Lee 2026",
        ),
        (
            "all_four",
            [
                by_key["jo2015"],
                by_key["halla_defurne2015"],
                by_key["halla_georges2022"],
                by_key["pass1"],
            ],
            "Jo 2015 + Defurne 2015 + Georges 2022 + Lee 2026",
        ),
    ]

    rows = []
    thresholds = np.arange(0.01, 0.101, 0.01)
    for tag, cfg_bundles, label in configurations:
        cfg_dir = outdir / tag
        cfg_dir.mkdir(parents=True, exist_ok=True)
        threshold_rows = []

        for threshold in thresholds:
            specs, counts = _km15_selected_specs_for_bundles(
                cfg_bundles, selection, float(threshold)
            )
            fit = fit_sachs_family_multi_measurements(
                specs,
                family=family,
                bh_cut=float(threshold),
                add_moradi_bh_systematic=True,
            )
            row = {
                "configuration": tag,
                "configuration_label": label,
                "family": family,
                "threshold": float(threshold),
                "threshold_percent": 100.0 * float(threshold),
                **counts,
                **fit,
            }
            threshold_rows.append(row)
            if abs(float(threshold) - 0.05) < 1.0e-12:
                rows.append(row)
            #endif
        #endfor

        table = pd.DataFrame(threshold_rows)
        table.to_csv(cfg_dir / "km15_threshold_scan.csv", index=False)

        # Distinguish finite-|t| form-factor stability from the amplified
        # derivative-at-zero radius sensitivity.
        save_threshold_form_factor_band_study(
            bundles=cfg_bundles,
            selection=selection,
            family=family,
            threshold_table=table,
            label=label,
            outdir=cfg_dir / "form_factor_band_stability",
        )

        fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.7))
        axes[0].errorbar(
            table["threshold_percent"], table["rE_fm"],
            yerr=table["rE_fit_err_fm"], marker="o", linestyle="-",
        )
        axes[1].errorbar(
            table["threshold_percent"], table["rM_fm"],
            yerr=table["rM_fit_err_fm"], marker="o", linestyle="-",
        )
        axes[2].plot(
            table["threshold_percent"], table["chi2_ndof"],
            marker="o", linestyle="-",
        )
        axes[0].set_ylabel(r"$r_E$ (fm)")
        axes[1].set_ylabel(r"$r_M$ (fm)")
        axes[2].set_ylabel(r"$\chi^2/\mathrm{dof}$")
        for ax in axes:
            ax.set_xlabel(r"$|1-R_{\rm BH}^{\rm KM15}|$ threshold (%)")
            ax.set_xscale("log")
            ax.axvline(5.0, linewidth=0.8, linestyle=":")
            ax.grid(alpha=0.2)
        #endfor
        fig.suptitle(f"{label}: KM15-only threshold stability ({family})", y=0.995)
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        fig.savefig(cfg_dir / "01_km15_threshold_stability.png", dpi=260)
        plt.close(fig)
    #endfor

    nominal = pd.DataFrame(rows)
    nominal.to_csv(outdir / "km15_5pct_ensemble_summary.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.0))
    x = np.arange(len(nominal))
    axes[0].errorbar(
        x, nominal["rE_fm"], yerr=nominal["rE_fit_err_fm"],
        marker="o", linestyle="none", capsize=3,
    )
    axes[1].errorbar(
        x, nominal["rM_fm"], yerr=nominal["rM_fit_err_fm"],
        marker="o", linestyle="none", capsize=3,
    )
    labels = nominal["configuration_label"].astype(str).tolist()
    for ax, ylabel in zip(axes, [r"$r_E$ (fm)", r"$r_M$ (fm)"]):
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=18, ha="right")
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", alpha=0.2)
    #endfor
    fig.suptitle(f"KM15 5% BH-purity ensemble comparison ({family})", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(outdir / "00_km15_5pct_ensemble_radii.png", dpi=180)
    plt.close(fig)

    print("\n[KM15 ensemble rotation] nominal 5% results")
    print(
        nominal[[
            "configuration_label", "N", "chi2_ndof",
            "rE_fm", "rE_fit_err_fm", "rM_fm", "rM_fit_err_fm",
        ]].to_string(index=False)
    )
    print(f"[KM15 ensemble rotation] outputs -> {outdir}")
#enddef



def run_global_saylor_tmin_scan(
        production_bundles: Sequence[Dict[str, object]],
        saylor_bundle: Dict[str, object],
        selection: pd.DataFrame,
        family: str,
        root_outdir: Path) -> pd.DataFrame:
    """
    Global six-dataset analogue of the Goharipour/Moradi/Azizi CLAS-2018
    low-|t| removal study.

    All five production datasets remain fixed at the nominal KM15 5% BH-purity
    selection.  Saylor 2018 is also first restricted to its KM15 5% sample,
    after which progressively larger LOW-|t| portions of Saylor alone are
    removed.  For each threshold t_cut, Saylor points satisfy |t| >= t_cut.

    The Sachs family is held fixed to the family selected by the all-five
    bias/variance closure machinery.  This is deliberate: otherwise changes in
    chi2 and radii would conflate Saylor-data removal with changing functional
    form.

    Saylor's published cross sections already include the experiment's elastic
    normalization correction.  The paper quotes a 4.0% uncertainty on that
    correction, which is treated here as one correlated Gaussian normalization
    nuisance; it is removed in quadrature from the tabulated total systematic
    uncertainty by load_saylor_supplement().
    """
    outdir = (
        root_outdir / "final_analysis" / "diagnostics"
        / "saylor_global_tmin_scan"
    )
    outdir.mkdir(parents=True, exist_ok=True)

    # Nominal 5% samples for the five production datasets.
    base_specs = []
    base_counts = {}
    for bundle in production_bundles:
        selected = select_bundle_from_external_model(
            bundle, selection, "km15", 0.05
        )
        base_counts[str(bundle["key"])] = int(len(selected))
        base_specs.append(measurement_spec(
            str(bundle["key"]),
            str(bundle["label"]),
            str(bundle["kind"]),
            selected,
            float(bundle.get("norm_frac", 0.0)),
            unconstrained_norm=bool(bundle.get("unconstrained_norm", False)),
        ))
    #endfor

    # Saylor is diagnostic-only and is intentionally absent from the standard
    # five-dataset external PARTONS selection table.  Its KM15 decomposition,
    # however, is already evaluated directly by run_saylor_validation() and
    # stored in saylor_bundle["all_data"].  The present scan is KM15-only, so
    # selecting Saylor from that authoritative in-memory decomposition avoids
    # any unnecessary PARTONS regeneration and keeps the production external
    # table five-dataset by construction.
    saylor_all = saylor_bundle["all_data"].copy().reset_index(drop=True)
    if "bh_delta" in saylor_all.columns:
        saylor5 = saylor_all.loc[
            pd.to_numeric(saylor_all["bh_delta"], errors="coerce") <= 0.05
        ].copy()
    elif "delta_bh" in saylor_all.columns:
        saylor5 = saylor_all.loc[
            pd.to_numeric(saylor_all["delta_bh"], errors="coerce") <= 0.05
        ].copy()
    else:
        raise RuntimeError(
            "Saylor bundle has no KM15 BH-purity column (expected bh_delta "
            "or delta_bh). Rebuild the Saylor KM15 decomposition cache."
        )
    #endif
    saylor5 = saylor5.reset_index(drop=True)
    if len(saylor5) == 0:
        raise RuntimeError("Saylor has zero points in the direct KM15 5% sample.")
    #endif
    print(
        f"[global Saylor |t|min scan] direct KM15 selection: "
        f"{len(saylor5)}/{len(saylor_all)} Saylor points at 5%; "
        "no external PARTONS table required for Saylor"
    )

    # Scan actual Saylor |t| bin boundaries, matching the logic of Ref.
    # arXiv:2607.04481: points below the threshold are removed.
    tvals = np.sort(np.unique(np.round(
        saylor5["t_abs"].to_numpy(float), 6
    )))
    # Baseline t_cut=0 keeps the complete 5%-selected Saylor sample.
    cuts = np.concatenate(([0.0], tvals))
    rows = []

    for tcut in cuts:
        kept = saylor5.loc[
            saylor5["t_abs"].to_numpy(float) >= float(tcut) - 1.0e-10
        ].copy()
        if len(kept) < 5:
            continue
        #endif

        specs = list(base_specs)
        specs.append(measurement_spec(
            str(saylor_bundle["key"]),
            str(saylor_bundle["label"]),
            str(saylor_bundle["kind"]),
            kept,
            SAYLOR_GLOBAL_SCALE_FRAC,
        ))

        fit = fit_sachs_family_multi_measurements(
            specs,
            family=family,
            bh_cut=0.05,
            add_moradi_bh_systematic=True,
        )

        beta_key = "beta_" + re.sub(
            r"[^A-Za-z0-9]+", "_", str(saylor_bundle["key"])
        )
        rows.append({
            "saylor_tmin_GeV2": float(tcut),
            "saylor_N": int(len(kept)),
            "saylor_N_removed": int(len(saylor5) - len(kept)),
            "global_N": int(fit["N"]),
            "family": str(family),
            "valid": bool(fit["valid"]),
            "chi2": float(fit["chi2"]),
            "ndof": int(fit["ndof"]),
            "chi2_ndof": float(fit["chi2_ndof"]),
            "rE_fm": float(fit["rE_fm"]),
            "rE_fit_err_fm": float(fit["rE_fit_err_fm"]),
            "rM_fm": float(fit["rM_fm"]),
            "rM_fit_err_fm": float(fit["rM_fit_err_fm"]),
            "saylor_beta_norm": float(fit.get(beta_key, np.nan)),
            "saylor_scale_factor": float(
                fit.get(beta_key + "_scale_factor", np.nan)
            ),
            "saylor_norm_prior_fraction": SAYLOR_GLOBAL_SCALE_FRAC,
        })
    #endfor

    # Add the Saylor-excluded all-five endpoint explicitly.
    fit_exc = fit_sachs_family_multi_measurements(
        base_specs,
        family=family,
        bh_cut=0.05,
        add_moradi_bh_systematic=True,
    )
    rows.append({
        "saylor_tmin_GeV2": np.nan,
        "saylor_N": 0,
        "saylor_N_removed": int(len(saylor5)),
        "global_N": int(fit_exc["N"]),
        "family": str(family),
        "valid": bool(fit_exc["valid"]),
        "chi2": float(fit_exc["chi2"]),
        "ndof": int(fit_exc["ndof"]),
        "chi2_ndof": float(fit_exc["chi2_ndof"]),
        "rE_fm": float(fit_exc["rE_fm"]),
        "rE_fit_err_fm": float(fit_exc["rE_fit_err_fm"]),
        "rM_fm": float(fit_exc["rM_fm"]),
        "rM_fit_err_fm": float(fit_exc["rM_fit_err_fm"]),
        "saylor_beta_norm": np.nan,
        "saylor_scale_factor": np.nan,
        "saylor_norm_prior_fraction": SAYLOR_GLOBAL_SCALE_FRAC,
    })

    table = pd.DataFrame(rows)
    table.to_csv(outdir / "global_six_dataset_saylor_tmin_scan.csv", index=False)

    def scan_saylor_lower_tail(
            column: str,
            output_stem: str,
            display_label: str) -> pd.DataFrame:
        """Remove progressively larger lower tails of one Saylor variable."""
        vals = np.sort(np.unique(np.round(
            pd.to_numeric(saylor5[column], errors="coerce").to_numpy(float),
            6,
        )))
        vals = vals[np.isfinite(vals)]
        if len(vals) == 0:
            return pd.DataFrame()
        #endif

        scan_rows = []
        baseline = float(vals[0]) - 1.0e-9
        for cut in np.concatenate(([baseline], vals)):
            kept = saylor5.loc[
                pd.to_numeric(
                    saylor5[column], errors="coerce"
                ).to_numpy(float) >= float(cut) - 1.0e-10
            ].copy()
            if len(kept) < 5:
                continue
            #endif

            specs = list(base_specs)
            specs.append(measurement_spec(
                str(saylor_bundle["key"]),
                str(saylor_bundle["label"]),
                str(saylor_bundle["kind"]),
                kept,
                SAYLOR_GLOBAL_SCALE_FRAC,
            ))
            fit = fit_sachs_family_multi_measurements(
                specs,
                family=family,
                bh_cut=0.05,
                add_moradi_bh_systematic=True,
            )
            beta_key = "beta_" + re.sub(
                r"[^A-Za-z0-9]+", "_", str(saylor_bundle["key"])
            )
            scan_rows.append({
                "cut_value": float(cut),
                "saylor_N": int(len(kept)),
                "saylor_N_removed": int(len(saylor5) - len(kept)),
                "global_N": int(fit["N"]),
                "family": str(family),
                "valid": bool(fit["valid"]),
                "chi2": float(fit["chi2"]),
                "ndof": int(fit["ndof"]),
                "chi2_ndof": float(fit["chi2_ndof"]),
                "rE_fm": float(fit["rE_fm"]),
                "rE_fit_err_fm": float(fit["rE_fit_err_fm"]),
                "rM_fm": float(fit["rM_fm"]),
                "rM_fit_err_fm": float(fit["rM_fit_err_fm"]),
                "saylor_beta_norm": float(fit.get(beta_key, np.nan)),
                "saylor_scale_factor": float(
                    fit.get(beta_key + "_scale_factor", np.nan)
                ),
            })
        #endfor

        scan_rows.append({
            "cut_value": np.nan,
            "saylor_N": 0,
            "saylor_N_removed": int(len(saylor5)),
            "global_N": int(fit_exc["N"]),
            "family": str(family),
            "valid": bool(fit_exc["valid"]),
            "chi2": float(fit_exc["chi2"]),
            "ndof": int(fit_exc["ndof"]),
            "chi2_ndof": float(fit_exc["chi2_ndof"]),
            "rE_fm": float(fit_exc["rE_fm"]),
            "rE_fit_err_fm": float(fit_exc["rE_fit_err_fm"]),
            "rM_fm": float(fit_exc["rM_fm"]),
            "rM_fit_err_fm": float(fit_exc["rM_fit_err_fm"]),
            "saylor_beta_norm": np.nan,
            "saylor_scale_factor": np.nan,
        })
        scan_table = pd.DataFrame(scan_rows)
        scan_table.to_csv(
            outdir / f"global_six_dataset_saylor_{output_stem}_scan.csv",
            index=False,
        )

        finite = scan_table.loc[
            scan_table["cut_value"].notna()
        ].copy()
        excluded_row = scan_table.loc[
            scan_table["cut_value"].isna()
        ].iloc[0]

        # Use approximately a dozen equal-retained-fraction points for the
        # displayed figure. Full bin-boundary resolution stays in the CSV.
        if len(finite) > 12:
            target_n = np.linspace(
                float(finite["saylor_N"].max()),
                float(finite["saylor_N"].min()),
                12,
            )
            take = []
            narr = finite["saylor_N"].to_numpy(float)
            for nt in target_n:
                take.append(int(np.argmin(np.abs(narr - nt))))
            #endfor
            shown = finite.iloc[sorted(set(take))].copy()
        else:
            shown = finite.copy()
        #endif

        fig, axes = plt.subplots(3, 1, figsize=(8.4, 9.8), sharex=True)
        axes[0].plot(
            shown["cut_value"], shown["chi2_ndof"], "o-", lw=1.4, ms=4
        )
        axes[0].axhline(
            float(excluded_row["chi2_ndof"]), ls="--", lw=1.0
        )
        axes[0].set_ylabel(r"Global $\chi^2/\mathrm{d.o.f.}$")

        axes[1].errorbar(
            shown["cut_value"], shown["rE_fm"],
            yerr=shown["rE_fit_err_fm"],
            fmt="o-", lw=1.2, ms=4, capsize=2,
        )
        axes[1].axhline(float(excluded_row["rE_fm"]), ls="--", lw=1.0)
        axes[1].set_ylabel(r"$r_E$ [fm]")

        axes[2].errorbar(
            shown["cut_value"], shown["rM_fm"],
            yerr=shown["rM_fit_err_fm"],
            fmt="o-", lw=1.2, ms=4, capsize=2,
        )
        axes[2].axhline(float(excluded_row["rM_fm"]), ls="--", lw=1.0)
        axes[2].set_ylabel(r"$r_M$ [fm]")
        axes[2].set_xlabel(display_label)

        for ax in axes:
            ax.grid(alpha=0.2)
        #endfor
        fig.suptitle(
            f"Global Saylor lower-tail removal scan: {display_label} ({family})",
            y=0.992,
        )
        fig.tight_layout(rect=(0, 0, 1, 0.965))
        fig.savefig(
            outdir / f"03_global_chi2_and_radii_vs_saylor_{output_stem}.png",
            dpi=270,
        )
        plt.close(fig)
        return scan_table
    #enddef

    finite_cut = table["saylor_tmin_GeV2"].notna()
    scan = table.loc[finite_cut].copy()
    excluded = table.loc[~finite_cut].iloc[0]

    # Main publication/note diagnostic: preserve every threshold in the CSV,
    # but display a compact representative subset.  Explicitly retain the two
    # literature-motivated Saylor cuts used in arXiv:2607.04481 when a nearby
    # actual bin boundary exists.
    if len(scan) > 14:
        target_n = np.linspace(
            float(scan["saylor_N"].max()),
            float(scan["saylor_N"].min()),
            12,
        )
        idx = []
        narr = scan["saylor_N"].to_numpy(float)
        for nt in target_n:
            idx.append(int(np.argmin(np.abs(narr - nt))))
        #endfor
        tarr = scan["saylor_tmin_GeV2"].to_numpy(float)
        for target_t in [0.265, 0.343]:
            idx.append(int(np.argmin(np.abs(tarr - target_t))))
        #endfor
        shown_scan = scan.iloc[sorted(set(idx))].copy()
    else:
        shown_scan = scan.copy()
    #endif

    fig, axes = plt.subplots(3, 1, figsize=(8.2, 10.0), sharex=True)

    axes[0].plot(
        shown_scan["saylor_tmin_GeV2"], shown_scan["chi2_ndof"],
        "o-", lw=1.5, ms=4
    )
    axes[0].axhline(
        float(excluded["chi2_ndof"]), ls="--", lw=1.2,
        label=f"Saylor excluded: {float(excluded['chi2_ndof']):.3f}"
    )
    axes[0].set_ylabel(r"Global $\chi^2/\mathrm{d.o.f.}$")
    axes[0].legend(frameon=False)

    axes[1].errorbar(
        shown_scan["saylor_tmin_GeV2"], shown_scan["rE_fm"],
        yerr=shown_scan["rE_fit_err_fm"], fmt="o-", lw=1.3, ms=4, capsize=2
    )
    axes[1].axhline(float(excluded["rE_fm"]), ls="--", lw=1.2)
    axes[1].set_ylabel(r"$r_E$ [fm]")

    axes[2].errorbar(
        shown_scan["saylor_tmin_GeV2"], shown_scan["rM_fm"],
        yerr=shown_scan["rM_fit_err_fm"], fmt="o-", lw=1.3, ms=4, capsize=2
    )
    axes[2].axhline(float(excluded["rM_fm"]), ls="--", lw=1.2)
    axes[2].set_ylabel(r"$r_M$ [fm]")
    axes[2].set_xlabel(
        r"Minimum retained Saylor $|t|$ [GeV$^2$]"
    )

    fig.suptitle(
        f"Global six-dataset Saylor low-$|t|$ removal scan ({family})"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.savefig(outdir / "01_global_chi2_and_radii_vs_saylor_tmin.png", dpi=260)
    plt.close(fig)

    # A second diagnostic shows exactly how much Saylor is retained and whether
    # its correlated normalization nuisance attempts to absorb the discrepancy.
    fig, ax1 = plt.subplots(figsize=(8.2, 5.4))
    ax1.plot(
        scan["saylor_tmin_GeV2"], scan["saylor_N"], "o-", lw=1.5, ms=4
    )
    ax1.set_xlabel(r"Minimum retained Saylor $|t|$ [GeV$^2$]")
    ax1.set_ylabel("Retained Saylor points")
    ax2 = ax1.twinx()
    ax2.plot(
        scan["saylor_tmin_GeV2"], scan["saylor_scale_factor"],
        "s--", lw=1.3, ms=4
    )
    ax2.axhline(1.0, ls=":", lw=1.0)
    ax2.set_ylabel("Fitted Saylor normalization scale")
    fig.suptitle("Saylor retention and correlated normalization")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(outdir / "02_saylor_retention_and_normalization.png", dpi=260)
    plt.close(fig)

    q2_scan = scan_saylor_lower_tail(
        "Q2", "q2min", r"Minimum retained Saylor $Q^2$ [GeV$^2$]"
    )
    xb_scan = scan_saylor_lower_tail(
        "xB", "xbmin", r"Minimum retained Saylor $x_B$"
    )

    print("\n[global Saylor |t|min scan]")
    print(
        table[[
            "saylor_tmin_GeV2", "saylor_N", "global_N", "chi2_ndof",
            "rE_fm", "rE_fit_err_fm", "rM_fm", "rM_fit_err_fm",
            "saylor_scale_factor",
        ]].to_string(index=False)
    )
    print(f"[global Saylor |t|min scan] outputs -> {outdir}")
    return table
#enddef



def run_saylor_recovery_study(
        saylor_bundle: Dict[str, object],
        args,
        root_outdir: Path,
        family: str = "IP2") -> None:
    """
    Diagnose whether a restricted Saylor kinematic region is internally usable.

    This is deliberately a diagnostic, NOT an automatic prescription for adding
    Saylor to the final radius ensemble.  It scans Q2 in two complementary ways:
      * fixed equal-population Q2 bins over all published points, testing the
        full KM15 electroproduction prediction;
      * sliding Q2 windows, testing both full-EP agreement and the quality of
        the KM15-5%-selected BH form-factor fit.

    A region should only be promoted to a production selection if its boundary
    can be justified independently of obtaining a favorable radius/chi2.
    """
    outdir = root_outdir / "final_analysis" / "diagnostics" / "saylor_recovery"
    outdir.mkdir(parents=True, exist_ok=True)

    d = saylor_bundle["all_data"].copy()
    finite = (
        np.isfinite(d["Q2"])
        & np.isfinite(d["xs"])
        & np.isfinite(d["km15_ep"])
        & (d["xs"] > 0.0)
    )
    d = d.loc[finite].copy().reset_index(drop=True)

    stat = dataset_statistical_errors(d, "saylor2018")
    ptp = d["ptp_sys"].to_numpy(float)
    err = np.sqrt(stat**2 + ptp**2)
    d["_ep_err"] = err

    def score_region(region: pd.DataFrame) -> Dict[str, float]:
        y = region["xs"].to_numpy(float)
        m = region["km15_ep"].to_numpy(float)
        e = region["_ep_err"].to_numpy(float)
        beta, scale, chi2 = _profile_model_normalization(
            y, m, e, SAYLOR_GLOBAL_SCALE_FRAC
        )
        pull = (scale * m - y) / e
        return {
            "N_all": int(len(region)),
            "Q2_mean": float(np.mean(region["Q2"])),
            "xB_mean": float(np.mean(region["xB"])),
            "t_mean": float(np.mean(region["t_abs"])),
            "ep_beta_norm": float(beta),
            "ep_scale": float(scale),
            "ep_chi2_per_point": float(chi2 / max(1, len(region))),
            "ep_mean_pull": float(np.mean(pull)),
            "ep_rms_pull": float(np.sqrt(np.mean(pull**2))),
            "ep_median_abs_frac_residual": float(
                np.median(np.abs(scale * m - y) / y)
            ),
        }
    #enddef

    # Equal-population Q2 bins: robust overview without hand-picked boundaries.
    q = d["Q2"].to_numpy(float)
    edges = np.unique(np.quantile(q, np.linspace(0.0, 1.0, 9)))
    qbin_rows = []
    for ibin, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        mask = (d["Q2"] >= lo) & (
            (d["Q2"] < hi) if ibin < len(edges) - 2 else (d["Q2"] <= hi)
        )
        region = d.loc[mask].copy()
        if len(region) == 0:
            continue
        #endif
        row = {"bin": ibin, "Q2_low": float(lo), "Q2_high": float(hi)}
        row.update(score_region(region))
        qbin_rows.append(row)
    #endfor
    qbins = pd.DataFrame(qbin_rows)
    qbins.to_csv(outdir / "01_saylor_km15_ep_equal_population_q2_bins.csv", index=False)

    # Sliding windows.  Width is intentionally broad enough to retain useful
    # statistics and stepped finely enough to reveal whether a central plateau
    # exists rather than selecting one favorable bin.
    qmin = float(np.nanmin(q))
    qmax = float(np.nanmax(q))
    width = 1.0
    step = 0.25
    starts = np.arange(qmin, max(qmin, qmax - width) + 0.5 * step, step)
    sliding_rows = []
    for lo in starts:
        hi = min(float(lo + width), qmax + 1.0e-12)
        region = d.loc[(d["Q2"] >= lo) & (d["Q2"] < hi)].copy()
        if len(region) < 40:
            continue
        #endif
        row = {"Q2_low": float(lo), "Q2_high": float(hi)}
        row.update(score_region(region))

        bh = region.loc[region["bh_delta"] <= 0.05].copy()
        row["N_bh5"] = int(len(bh))
        if len(bh) >= 20:
            try:
                spec = measurement_spec(
                    key="saylor2018",
                    label="CLAS6 Saylor 2018",
                    kind="saylor2018",
                    data=bh,
                    norm_frac=SAYLOR_GLOBAL_SCALE_FRAC,
                )
                fit = fit_sachs_family_multi_measurements(
                    [spec], family=family, bh_cut=0.05,
                    add_moradi_bh_systematic=True,
                )
                row["bh5_valid"] = bool(fit["valid"])
                row["bh5_chi2_ndof"] = float(fit["chi2_ndof"])
                row["bh5_rE_fm"] = float(fit["rE_fm"])
                row["bh5_rE_fit_err_fm"] = float(fit["rE_fit_err_fm"])
                row["bh5_rM_fm"] = float(fit["rM_fm"])
                row["bh5_rM_fit_err_fm"] = float(fit["rM_fit_err_fm"])
            except Exception as exc:
                row["bh5_valid"] = False
                row["bh5_failure"] = str(exc)
            #endtry
        else:
            row["bh5_valid"] = False
            row["bh5_failure"] = "fewer than 20 KM15-5%-selected points"
        #endif
        sliding_rows.append(row)
    #endfor
    sliding = pd.DataFrame(sliding_rows)
    sliding.to_csv(outdir / "02_saylor_q2_sliding_window_scan.csv", index=False)

    # Overview plots.
    if len(qbins):
        centers = 0.5 * (qbins["Q2_low"] + qbins["Q2_high"])
        fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
        axes[0].plot(centers, qbins["ep_chi2_per_point"], marker="o")
        axes[1].plot(
            centers, 100.0 * qbins["ep_median_abs_frac_residual"], marker="o"
        )
        axes[0].set_ylabel(r"KM15 full-EP profiled $\chi^2/N$")
        axes[1].set_ylabel("Median absolute fractional residual (%)")
        for ax in axes:
            ax.set_xlabel(r"$Q^2$ (GeV$^2$)")
            ax.grid(alpha=0.2)
        #endfor
        fig.suptitle("Saylor 2018: KM15 agreement versus Q2", y=0.995)
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        fig.savefig(outdir / "01_saylor_km15_ep_vs_q2.png", dpi=180)
        plt.close(fig)
    #endif

    if len(sliding):
        centers = 0.5 * (sliding["Q2_low"] + sliding["Q2_high"])
        fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
        axes[0].plot(centers, sliding["ep_chi2_per_point"], marker="o")
        valid = sliding.get("bh5_valid", pd.Series(False, index=sliding.index)).astype(bool)
        axes[1].plot(
            centers[valid], sliding.loc[valid, "bh5_chi2_ndof"], marker="o"
        )
        axes[0].set_ylabel(r"Full-EP KM15 profiled $\chi^2/N$")
        axes[1].set_ylabel(r"KM15 5% BH-fit $\chi^2/\mathrm{dof}$")
        for ax in axes:
            ax.set_xlabel(r"center of 1-GeV$^2$ $Q^2$ window")
            ax.grid(alpha=0.2)
        #endfor
        fig.suptitle(
            "Saylor 2018 recoverability scan: full electroproduction and BH extraction",
            y=0.995,
        )
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        fig.savefig(outdir / "02_saylor_q2_recovery_scan.png", dpi=180)
        plt.close(fig)
    #endif

    print(f"[Saylor recovery] outputs -> {outdir}")
#enddef



def save_nominal_family_failure_diagnostic(
        audit: pd.DataFrame,
        configuration: str,
        configuration_label: str,
        outdir: Path) -> None:
    """
    Visualize the case where closure-ranked GE/GM Sachs families fail on the
    actual nominal 5% BH-selected dataset.

    This is a numerical/identifiability diagnostic, not a physics fit result.
    It deliberately shows both the closure ranking and the outcome of the
    corresponding nominal-data minimization attempts.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    part = audit.loc[
        audit["configuration"].astype(str) == str(configuration)
    ].copy()
    if len(part) == 0:
        return
    #endif

    part = part.sort_values("closure_rank").reset_index(drop=True)
    part.to_csv(
        outdir / "01_nominal_family_failure_audit.csv",
        index=False,
    )

    ranks = part["closure_rank"].to_numpy(float)
    objective = part["closure_objective_fm"].to_numpy(float)
    global_eff = part["global_valid_fraction"].to_numpy(float)
    worst_eff = part["minimum_scenario_valid_fraction"].to_numpy(float)
    nominal_valid = part["nominal_fit_valid"].astype(bool).to_numpy()

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))

    axes[0].plot(
        ranks, objective,
        marker="o", linestyle="none",
    )
    if np.any(nominal_valid):
        axes[0].plot(
            ranks[nominal_valid], objective[nominal_valid],
            marker="o", linestyle="none",
            markersize=9,
            fillstyle="none",
            label="Nominal fit valid",
        )
    #endif
    if np.any(~nominal_valid):
        axes[0].plot(
            ranks[~nominal_valid], objective[~nominal_valid],
            marker="x", linestyle="none",
            markersize=8,
            label="Nominal fit invalid",
        )
    #endif
    axes[0].set_xlabel("Closure rank")
    axes[0].set_ylabel("Closure combined RMSE objective [fm]")
    axes[0].set_title("Closure ranking vs nominal-fit outcome")
    axes[0].grid(alpha=0.2)
    axes[0].legend()

    axes[1].plot(
        ranks, global_eff,
        marker="o", linestyle="none",
        label="Global convergence efficiency",
    )
    axes[1].plot(
        ranks, worst_eff,
        marker="s", linestyle="none",
        fillstyle="none",
        label="Worst-truth convergence efficiency",
    )
    axes[1].set_xlabel("Closure rank")
    axes[1].set_ylabel("Convergence fraction")
    axes[1].set_ylim(-0.03, 1.03)
    axes[1].set_title("Closure numerical robustness")
    axes[1].grid(alpha=0.2)
    axes[1].legend()

    fig.suptitle(
        f"{configuration_label}: no valid free-$G_E/G_M$ nominal fit",
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(
        outdir / "01_nominal_family_failure_summary.png",
        dpi=180,
    )
    plt.close(fig)

    # A second compact plot shows any finite radius values returned by invalid
    # fits.  They are intentionally labeled invalid and are not interpreted as
    # measurements; this helps diagnose whether failure comes from a runaway
    # electric, magnetic, or both sectors.
    rE = pd.to_numeric(part["nominal_rE_fm"], errors="coerce").to_numpy(float)
    rM = pd.to_numeric(part["nominal_rM_fm"], errors="coerce").to_numpy(float)
    finite = np.isfinite(rE) | np.isfinite(rM)
    if np.any(finite):
        fig, ax = plt.subplots(figsize=(8.2, 5.0))
        ax.plot(
            ranks[np.isfinite(rE)],
            rE[np.isfinite(rE)],
            marker="o", linestyle="none",
            label=r"Returned $r_E$ (invalid fit)",
        )
        ax.plot(
            ranks[np.isfinite(rM)],
            rM[np.isfinite(rM)],
            marker="s", linestyle="none",
            fillstyle="none",
            label=r"Returned $r_M$ (invalid fit)",
        )
        ax.axhline(SACHS_MIN_RADIUS_FM, linewidth=0.8, linestyle=":")
        ax.axhline(SACHS_MAX_RADIUS_FM, linewidth=0.8, linestyle=":")
        ax.set_xlabel("Closure rank")
        ax.set_ylabel("Returned radius [fm]")
        ax.set_title(
            f"{configuration_label}: invalid nominal-fit radius returns"
        )
        ax.grid(alpha=0.2)
        ax.legend()
        fig.tight_layout()
        fig.savefig(
            outdir / "02_invalid_nominal_returned_radii.png",
            dpi=180,
        )
        plt.close(fig)
    #endif
#enddef


def save_lee_fit8_kelly_threshold_scan(
        lee_bundle: Dict[str, object],
        selection: pd.DataFrame,
        summary_dir: Path,
        outdir: Path) -> pd.DataFrame:
    """
    Moradi Fit-8-style Lee-only fallback/diagnostic:
      F1 is fitted with the paper's two-parameter power form;
      F2 is fixed to Kelly.

    This is run even when the simultaneous free-GE/free-GM Sachs-family
    extraction fails.  It therefore directly tests whether Lee still contains
    useful electric information once the poorly constrained magnetic sector is
    supplied externally.

    The threshold scan holds the additional Moradi BH-method uncertainty fixed
    at 5% by rescaling the bh_cut passed to the fitter only for the uncertainty
    term.  Selection itself is performed independently at each threshold.
    """
    outdir.mkdir(parents=True, exist_ok=True)

    scan_thresholds = [
        0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.10,
        0.15, 0.20, 0.30, 0.50, 0.75, 1.00,
    ]
    key = str(lee_bundle["key"])
    # External model-selection schema uses dataset and delta_bh_km15.
    sel = selection.loc[
        (selection["dataset"].astype(str) == key)
        & np.isfinite(selection["delta_bh_km15"].to_numpy(float))
    ]
    if len(sel):
        scan_thresholds.append(float(sel["delta_bh_km15"].max()) * 1.001)
    #endif
    scan_thresholds = sorted(set(float(x) for x in scan_thresholds))

    rows = []
    for threshold in scan_thresholds:
        selected = select_bundle_from_external_model(
            lee_bundle, selection, "km15", float(threshold)
        )
        if len(selected) < 4:
            continue
        #endif

        # fit_paper_model currently uses bh_cut for the added Moradi uncertainty.
        # To keep that method uncertainty fixed at 5% during the loose-threshold
        # contamination test, make a temporary copy whose xs errors are handled
        # by a local wrapper below.
        fit = fit_paper_model(
            selected,
            kind="fit8_f2_kelly",
            fit_name=f"Lee2026_Fit8_F2_Kelly_{100.0*threshold:.3f}pct",
            bh_cut=0.05,
            include_clas12_ptp_sys=True,
            include_bh_sys=True,
            include_combination_norm_sys=False,
            include_global_norm_sys=True,
            global_scale_frac=PASS1_GLOBAL_SCALE_FRAC,
        )

        row = fitresult_to_record(fit)
        row.update({
            "selection_threshold": float(threshold),
            "selection_threshold_percent": 100.0 * float(threshold),
            "N_selected": int(len(selected)),
            "F2_treatment": "Kelly fixed",
            "fit_kind": "fit8_f2_kelly",
        })
        rows.append(row)
    #endfor

    table = pd.DataFrame(rows)
    table.to_csv(
        outdir / "01_lee2026_fit8_kelly_threshold_scan.csv",
        index=False,
    )

    nominal = table.loc[
        np.isclose(table["selection_threshold"].to_numpy(float), 0.05)
    ].copy()
    nominal.to_csv(
        summary_dir / "lee2026_fit8_f2_kelly.csv",
        index=False,
    )

    if len(table):
        fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.8))

        axes[0].errorbar(
            table["selection_threshold_percent"],
            table["rE_fm"],
            yerr=table["rE_err_fm"],
            marker="o",
            linestyle="none",
            capsize=2,
        )
        axes[0].axvline(5.0, linewidth=0.8, linestyle=":")
        axes[0].set_xscale("log")
        axes[0].set_xlabel(
            r"$|1-R_{\rm BH}^{\rm KM15}|$ selection threshold (%)"
        )
        axes[0].set_ylabel(r"$r_E$ [fm]")
        axes[0].set_title(r"Lee-only Fit 8: $F_2$ fixed to Kelly")
        axes[0].grid(alpha=0.2)

        axes[1].plot(
            table["selection_threshold_percent"],
            table["chi2_ndof"],
            marker="o",
            linestyle="none",
        )
        axes[1].axvline(5.0, linewidth=0.8, linestyle=":")
        axes[1].set_xscale("log")
        axes[1].set_xlabel(
            r"$|1-R_{\rm BH}^{\rm KM15}|$ selection threshold (%)"
        )
        axes[1].set_ylabel(r"$\chi^2/\mathrm{dof}$")
        axes[1].set_title("Fit quality")
        axes[1].grid(alpha=0.2)

        fig.suptitle(
            "CLAS12 Lee 2026: electric-radius extraction with external magnetic input",
            y=0.995,
        )
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        fig.savefig(
            outdir / "01_lee2026_fit8_kelly_threshold_scan.png",
            dpi=180,
        )
        plt.close(fig)
    #endif

    return table
#enddef



def save_selected_experiment_kinematic_coverage(
        bundles: Sequence[Dict[str, object]],
        selection: pd.DataFrame,
        outdir: Path,
        threshold: float = 0.05) -> pd.DataFrame:
    """Compare the KM15-selected kinematic support of the production datasets."""
    outdir.mkdir(parents=True, exist_ok=True)
    selected_by_label = {}
    summary_rows = []
    for bundle in bundles:
        d = select_bundle_from_external_model(
            bundle, selection, "km15", float(threshold)
        )
        if len(d) == 0:
            continue
        #endif
        label = str(bundle["label"])
        selected_by_label[label] = d
        stat = dataset_statistical_errors(d, str(bundle["kind"]))
        point = dataset_point_errors(
            d, str(bundle["kind"]), float(threshold), False
        )
        xs = np.abs(d["xs"].to_numpy(float))
        rel_stat = np.divide(
            stat, xs, out=np.full_like(stat, np.nan), where=xs > 0.0
        )
        rel_point = np.divide(
            point, xs, out=np.full_like(point, np.nan), where=xs > 0.0
        )
        summary_rows.append({
            "dataset": label,
            "N_selected": int(len(d)),
            "Q2_min": float(np.nanmin(d["Q2"])),
            "Q2_median": float(np.nanmedian(d["Q2"])),
            "Q2_max": float(np.nanmax(d["Q2"])),
            "xB_min": float(np.nanmin(d["xB"])),
            "xB_median": float(np.nanmedian(d["xB"])),
            "xB_max": float(np.nanmax(d["xB"])),
            "t_abs_min": float(np.nanmin(d["t_abs"])),
            "t_abs_median": float(np.nanmedian(d["t_abs"])),
            "t_abs_max": float(np.nanmax(d["t_abs"])),
            "phi_deg_min": float(np.nanmin(d["phi_deg"])),
            "phi_deg_median": float(np.nanmedian(d["phi_deg"])),
            "phi_deg_max": float(np.nanmax(d["phi_deg"])),
            "median_relative_stat_error": float(np.nanmedian(rel_stat)),
            "median_relative_point_error": float(np.nanmedian(rel_point)),
        })
    #endfor

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(outdir / "selected_kinematic_coverage_summary.csv", index=False)

    if selected_by_label:
        fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.0))
        variables = [
            ("Q2", r"$Q^2$ [GeV$^2$]"),
            ("xB", r"$x_B$"),
            ("t_abs", r"$|t|$ [GeV$^2$]"),
            ("phi_deg", r"$\phi$ [deg]"),
        ]
        for ax, (column, xlabel) in zip(axes.flat, variables):
            all_values = np.concatenate([
                d[column].to_numpy(float) for d in selected_by_label.values()
            ])
            finite = all_values[np.isfinite(all_values)]
            if len(finite) == 0:
                continue
            #endif
            edges = np.linspace(float(np.min(finite)), float(np.max(finite)), 31)
            for label, d in selected_by_label.items():
                values = d[column].to_numpy(float)
                values = values[np.isfinite(values)]
                ax.hist(
                    values, bins=edges, histtype="step", density=True,
                    linewidth=1.4, label=label,
                )
            #endfor
            ax.set_xlabel(xlabel)
            ax.set_ylabel("Normalized entries")
            ax.grid(alpha=0.2)
        #endfor
        handles, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(
            handles, labels, loc="upper center", ncol=2,
            bbox_to_anchor=(0.5, 0.955),
        )
        fig.suptitle(
            f"KM15 {100.0*threshold:.0f}% selected points: experimental kinematic coverage",
            y=0.995,
        )
        fig.tight_layout(rect=(0, 0, 1, 0.90))
        fig.savefig(outdir / "selected_kinematic_coverage_4variables.png", dpi=180)
        plt.close(fig)

        fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2))
        for label, d in selected_by_label.items():
            axes[0].scatter(
                d["xB"], d["Q2"], s=10, alpha=0.55, label=label,
            )
            axes[1].scatter(
                d["t_abs"], d["phi_deg"], s=10, alpha=0.55, label=label,
            )
        #endfor
        axes[0].set_xlabel(r"$x_B$")
        axes[0].set_ylabel(r"$Q^2$ [GeV$^2$]")
        axes[1].set_xlabel(r"$|t|$ [GeV$^2$]")
        axes[1].set_ylabel(r"$\phi$ [deg]")
        for ax in axes:
            ax.grid(alpha=0.2)
        #endfor
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(
            handles, labels, loc="upper center", ncol=2,
            bbox_to_anchor=(0.5, 0.955),
        )
        fig.suptitle(
            f"KM15 {100.0*threshold:.0f}% selected points: correlated kinematic support",
            y=0.995,
        )
        fig.tight_layout(rect=(0, 0, 1, 0.90))
        fig.savefig(outdir / "selected_kinematic_coverage_correlations.png", dpi=180)
        plt.close(fig)
    #endif
    return summary
#enddef


def _copy_specs_with_norm_fraction(
        specs: Sequence[Dict[str, object]],
        norm_fraction_by_key: Dict[str, float]) -> List[Dict[str, object]]:
    out = []
    for spec in specs:
        item = dict(spec)
        key = str(item["key"])
        if key in norm_fraction_by_key:
            item["norm_frac"] = float(norm_fraction_by_key[key])
        #endif
        out.append(item)
    #endfor
    return out
#enddef


def run_normalization_tension_diagnostics(
        bundles: Sequence[Dict[str, object]],
        selection: pd.DataFrame,
        family: str,
        outdir: Path) -> pd.DataFrame:
    """Compare fixed, published-constrained, and freely floating normalizations."""
    outdir.mkdir(parents=True, exist_ok=True)
    specs, counts = _km15_selected_specs_for_bundles(
        bundles, selection, 0.05
    )
    keys = [str(spec["key"]) for spec in specs]

    def clear_free_flags(input_specs):
        out = [dict(spec) for spec in input_specs]
        for spec in out:
            spec["unconstrained_norm"] = False
        #endfor
        return out
    #enddef

    modes = []
    constrained_specs = clear_free_flags(specs)
    modes.append((
        "published_constrained_georges_fixed",
        constrained_specs,
        [],
        "Published Gaussian priors where known; Georges fixed only as a comparison diagnostic.",
    ))
    modes.append((
        "all_fixed",
        _copy_specs_with_norm_fraction(
            clear_free_flags(specs), {k: 0.0 for k in keys}
        ),
        [],
        "All absolute normalizations fixed; comparison diagnostic.",
    ))
    if "halla_georges2022" in keys:
        georges_specs = _copy_specs_with_norm_fraction(
            clear_free_flags(specs), {"halla_georges2022": 1.0}
        )
        for spec in georges_specs:
            if str(spec["key"]) == "halla_georges2022":
                spec["unconstrained_norm"] = True
            #endif
        #endfor
        modes.append((
            "georges_free_baseline",
            georges_specs,
            ["halla_georges2022"],
            "Current baseline: published constraints retained for other datasets; Georges normalization unconstrained.",
        ))
    #endif
    all_free_specs = _copy_specs_with_norm_fraction(
        clear_free_flags(specs), {k: 1.0 for k in keys}
    )
    for spec in all_free_specs:
        spec["unconstrained_norm"] = True
    #endfor
    modes.append((
        "all_free_diagnostic",
        all_free_specs,
        keys,
        "Every dataset normalization unconstrained; diagnostic only.",
    ))

    rows = []
    for mode, mode_specs, free_keys, note in modes:
        fit = fit_sachs_family_multi_measurements(
            mode_specs,
            family=family,
            bh_cut=0.05,
            add_moradi_bh_systematic=True,
            bh_systematic_fraction=0.05,
            unconstrained_norm_keys=free_keys,
        )
        row = {
            "normalization_mode": mode,
            "family": family,
            "note": note,
            **counts,
            **fit,
        }
        rows.append(row)
    #endfor
    table = pd.DataFrame(rows)
    table.to_csv(outdir / "normalization_treatment_comparison.csv", index=False)

    # A compact scale-factor table is easier to read than the wide fit table.
    scale_rows = []
    for _, row in table.iterrows():
        for key in keys:
            bname = "beta_" + re.sub(r"[^A-Za-z0-9]+", "_", key)
            if bname in table.columns and np.isfinite(row.get(bname, np.nan)):
                scale_rows.append({
                    "normalization_mode": row["normalization_mode"],
                    "dataset": key,
                    "beta": float(row[bname]),
                    "beta_err": float(row.get(bname + "_err", np.nan)),
                    "prior_fraction": float(row.get(
                        bname + "_norm_fraction", np.nan
                    )),
                    "scale_factor": float(row.get(
                        bname + "_scale_factor", np.nan
                    )),
                    "unconstrained": bool(row.get(
                        bname + "_is_unconstrained", False
                    )),
                })
            else:
                scale_rows.append({
                    "normalization_mode": row["normalization_mode"],
                    "dataset": key,
                    "beta": np.nan,
                    "beta_err": np.nan,
                    "prior_fraction": 0.0,
                    "scale_factor": 1.0,
                    "unconstrained": False,
                })
            #endif
        #endfor
    #endfor
    pd.DataFrame(scale_rows).to_csv(
        outdir / "normalization_scale_factors.csv", index=False
    )

    print("\n[normalization tension] KM15 5% all-four comparison")
    display_cols = [
        "normalization_mode", "N", "chi2_ndof",
        "rE_fm", "rE_fit_err_fm", "rM_fm", "rM_fit_err_fm",
    ]
    print(table[display_cols].to_string(index=False))
    return table
#enddef


def run_fixed_family_ensemble_decomposition(
        jo_bundle: Dict[str, object],
        defurne_bundle: Dict[str, object],
        georges_bundle: Dict[str, object],
        lee_bundle: Dict[str, object],
        selection: pd.DataFrame,
        family: str,
        outdir: Path) -> pd.DataFrame:
    """Diagnose which Hall-A dataset drives the all-four radius rotation."""
    outdir.mkdir(parents=True, exist_ok=True)
    configurations = [
        ("jo_plus_lee", [jo_bundle, lee_bundle], "Jo + Lee"),
        (
            "jo_plus_lee_plus_defurne",
            [jo_bundle, lee_bundle, defurne_bundle],
            "Jo + Lee + Defurne",
        ),
        (
            "jo_plus_lee_plus_georges",
            [jo_bundle, lee_bundle, georges_bundle],
            "Jo + Lee + Georges",
        ),
        (
            "all_four",
            [jo_bundle, defurne_bundle, georges_bundle, lee_bundle],
            "Jo + Defurne + Georges + Lee",
        ),
        (
            "hall_a_only",
            [defurne_bundle, georges_bundle],
            "Defurne + Georges",
        ),
    ]
    rows = []
    for tag, bundles, label in configurations:
        specs, counts = _km15_selected_specs_for_bundles(
            bundles, selection, 0.05
        )
        fit = fit_sachs_family_multi_measurements(
            specs, family=family, bh_cut=0.05,
            add_moradi_bh_systematic=True,
            bh_systematic_fraction=0.05,
        )
        rows.append({
            "configuration": tag,
            "configuration_label": label,
            "family": family,
            **counts,
            **fit,
        })
    #endfor
    table = pd.DataFrame(rows)
    table.to_csv(outdir / "fixed_family_ensemble_decomposition.csv", index=False)

    print(f"\n[Hall-A decomposition] fixed family {family}")
    cols = [
        "configuration_label", "N", "chi2_ndof",
        "rE_fm", "rE_fit_err_fm", "rM_fm", "rM_fit_err_fm",
    ]
    print(table[cols].to_string(index=False))

    # Publication/meeting-facing visualization of the fixed-family ensemble
    # decomposition.  The same Sachs family is held fixed in every row, so the
    # movement isolates dataset leverage rather than a simultaneous family
    # change.
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.2))
    xx = np.arange(len(table))
    axes[0].errorbar(
        xx, table["rE_fm"], yerr=table["rE_fit_err_fm"],
        marker="o", linestyle="none", capsize=3,
    )
    axes[1].errorbar(
        xx, table["rM_fm"], yerr=table["rM_fit_err_fm"],
        marker="o", linestyle="none", capsize=3,
    )
    labels = table["configuration_label"].astype(str).tolist()
    for ax, ylabel in zip(axes, [r"$r_E$ [fm]", r"$r_M$ [fm]"]):
        ax.set_xticks(xx)
        ax.set_xticklabels(labels, rotation=20, ha="right")
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", alpha=0.2)
    #endfor
    fig.suptitle(
        f"Hall-A ensemble decomposition at KM15 5% (fixed {family})",
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(outdir / "01_fixed_family_ensemble_decomposition.png", dpi=180)
    plt.close(fig)

    # Normalization-scale companion plot for configurations in which a
    # published correlated nuisance exists. Georges remains at unity in the
    # production prescription because no authoritative prior decomposition is
    # supplied.
    scale_rows = []
    for _, row in table.iterrows():
        for key, short in [
            ("jo2015", "Jo"),
            ("halla_defurne2015", "Defurne"),
            ("halla_georges2022", "Georges"),
            ("pass1", "Lee"),
        ]:
            bname = "beta_" + re.sub(r"[^A-Za-z0-9]+", "_", key)
            scale = row.get(bname + "_scale_factor", np.nan)
            if not np.isfinite(scale):
                scale = 1.0
            #endif
            scale_rows.append({
                "configuration_label": row["configuration_label"],
                "dataset": short,
                "scale_factor": float(scale),
            })
        #endfor
    #endfor
    pd.DataFrame(scale_rows).to_csv(
        outdir / "fixed_family_ensemble_normalization_scales.csv", index=False
    )
    return table
#enddef


def run_data_driven_bh_deviance_for_ensemble(
        bundles: Sequence[Dict[str, object]],
        selection: pd.DataFrame,
        family: str,
        label: str,
        outdir: Path) -> pd.DataFrame:
    """
    Run the observed-data/Kelly-BH matched-statistics selector beside KM15 5%.

    WARNING: these selectors rank points using the same measured cross section
    that is subsequently fit.  They therefore condition on the dependent
    variable and can produce artificially small chi2 values.  They are retained
    only as overlap/selection-geometry diagnostics and MUST NOT be interpreted
    as an independent BH-selection systematic.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    nominal_specs, _ = _km15_selected_specs_for_bundles(
        bundles, selection, 0.05
    )
    nominal = fit_sachs_family_multi_measurements(
        nominal_specs, family=family, bh_cut=0.05,
        add_moradi_bh_systematic=True,
        bh_systematic_fraction=0.05,
    )
    diag = build_data_driven_bh_selection_diagnostics(
        physics_bundles=bundles,
        nominal_specs=nominal_specs,
        chosen_family=family,
        diagnostics_dir=outdir,
    )
    rows = [{
        "selector": "KM15_5pct",
        "ensemble_label": label,
        **nominal,
    }]
    for _, r in diag["fit_table"].iterrows():
        d = r.to_dict()
        d["ensemble_label"] = label
        rows.append(d)
    #endfor
    table = pd.DataFrame(rows)
    base_e = float(nominal["rE_fm"])
    base_m = float(nominal["rM_fm"])
    table["delta_rE_from_KM15_fm"] = table["rE_fm"].to_numpy(float) - base_e
    table["delta_rM_from_KM15_fm"] = table["rM_fm"].to_numpy(float) - base_m
    table["selection_uses_observed_xs"] = (
        table["selector"].astype(str) != "KM15_5pct"
    )
    table["production_systematic_eligible"] = (
        table["selector"].astype(str) == "KM15_5pct"
    )
    table["methodological_note"] = np.where(
        table["selection_uses_observed_xs"].to_numpy(bool),
        (
            "Diagnostic only: selector conditions on measured cross section, "
            "so fitted chi2 and radius shifts are not an independent "
            "BH-selection systematic."
        ),
        "Nominal production selector.",
    )
    table.to_csv(outdir / "data_driven_bh_radius_comparison.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))
    x = np.arange(len(table))
    axes[0].errorbar(
        x, table["rE_fm"], yerr=table["rE_fit_err_fm"],
        marker="o", linestyle="none", capsize=3,
    )
    axes[1].errorbar(
        x, table["rM_fm"], yerr=table["rM_fit_err_fm"],
        marker="o", linestyle="none", capsize=3,
    )
    labels = table["selector"].astype(str).tolist()
    for ax, ylabel in zip(axes, [r"$r_E$ [fm]", r"$r_M$ [fm]"]):
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=18, ha="right")
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", alpha=0.2)
    #endfor
    fig.suptitle(
        f"{label}: KM15 versus data-driven Kelly-BH selection ({family})",
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(outdir / "data_driven_bh_radius_comparison.png", dpi=180)
    plt.close(fig)
    return table
#enddef



def run_moradi_fit5_fit8_ensemble_benchmarks(
        jo_bundle: Dict[str, object],
        defurne_bundle: Dict[str, object],
        georges_bundle: Dict[str, object],
        lee_bundle: Dict[str, object],
        selection: pd.DataFrame,
        outdir: Path) -> pd.DataFrame:
    """
    Run the original Moradi Fit-5 and Fit-8 F1/F2 parameterizations on the same
    KM15-selected production ensembles used by the Sachs-family analysis.

    These are benchmark/interpretation fits, not substitutes for the
    closure-ranked production family.  Fit 5 has independent dipole F1/F2
    scales.  Fit 8 fits the F1 power form while fixing F2 to Kelly.  The latter
    therefore reports an electric-radius constraint with the magnetic sector
    externally supplied.

    Every experiment keeps its own point-error prescription and independent
    correlated normalization nuisance.  The added Moradi BH-method uncertainty
    is held at 5% while the selection threshold is scanned, matching the
    threshold-stability convention used elsewhere in this script.
    """
    outdir.mkdir(parents=True, exist_ok=True)

    configurations = [
        ("jo_only", [jo_bundle], "Jo 2015"),
        ("defurne_only", [defurne_bundle], "Defurne 2015"),
        ("georges_only", [georges_bundle], "Georges 2022"),
        ("lee_only", [lee_bundle], "Lee 2026"),
        ("jo_plus_lee", [jo_bundle, lee_bundle], "Jo + Lee"),
        (
            "jo_plus_lee_plus_defurne",
            [jo_bundle, lee_bundle, defurne_bundle],
            "Jo + Lee + Defurne",
        ),
        (
            "jo_plus_lee_plus_georges",
            [jo_bundle, lee_bundle, georges_bundle],
            "Jo + Lee + Georges",
        ),
        (
            "hall_a_only",
            [defurne_bundle, georges_bundle],
            "Defurne + Georges",
        ),
        (
            "all_four",
            [jo_bundle, defurne_bundle, georges_bundle, lee_bundle],
            "Jo + Defurne + Georges + Lee",
        ),
    ]

    model_defs = [
        ("Fit 5", "dipole"),
        ("Fit 8", "fit8_f2_kelly"),
    ]
    scan_thresholds = [
        0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.075,
        0.10, 0.15, 0.20,
    ]

    all_rows = []
    nominal_rows = []

    for tag, bundles, label in configurations:
        cfg_dir = outdir / tag
        cfg_dir.mkdir(parents=True, exist_ok=True)
        rows = []

        for threshold in scan_thresholds:
            specs, counts = _km15_selected_specs_for_bundles(
                bundles, selection, float(threshold)
            )
            if sum(counts.values()) < 4:
                continue
            #endif

            for fit_label, kind in model_defs:
                try:
                    fit = fit_multi_measurements(
                        specs,
                        kind=kind,
                        fit_name=fit_label,
                        bh_cut=0.05,
                        add_moradi_bh_systematic=True,
                    )
                    row = {
                        "configuration": tag,
                        "configuration_label": label,
                        "fit_label": fit_label,
                        "model_kind": kind,
                        "selection_threshold": float(threshold),
                        "selection_threshold_percent": 100.0 * float(threshold),
                        "bh_method_uncertainty_fraction": 0.05,
                        "fit8_F2_fixed_to_Kelly": bool(
                            kind == "fit8_f2_kelly"
                        ),
                        **counts,
                        "N": int(fit.npts),
                        "chi2": float(fit.chi2),
                        "ndof": int(fit.ndof),
                        "chi2_ndof": float(fit.chi2_ndof),
                        "success": bool(fit.success),
                        "rE_fm": float(fit.rE_fm),
                        "rE_fit_err_fm": float(fit.rE_err_fm),
                        "rM_fm": float(fit.rM_fm),
                        "rM_fit_err_fm": float(fit.rM_err_fm),
                        **dict(fit.meta),
                    }
                except Exception as exc:
                    row = {
                        "configuration": tag,
                        "configuration_label": label,
                        "fit_label": fit_label,
                        "model_kind": kind,
                        "selection_threshold": float(threshold),
                        "selection_threshold_percent": 100.0 * float(threshold),
                        "bh_method_uncertainty_fraction": 0.05,
                        "fit8_F2_fixed_to_Kelly": bool(
                            kind == "fit8_f2_kelly"
                        ),
                        **counts,
                        "N": int(sum(counts.values())),
                        "success": False,
                        "failure": str(exc),
                        "chi2_ndof": np.nan,
                        "rE_fm": np.nan,
                        "rE_fit_err_fm": np.nan,
                        "rM_fm": np.nan,
                        "rM_fit_err_fm": np.nan,
                    }
                #endtry
                rows.append(row)
                all_rows.append(row)
                if abs(float(threshold) - 0.05) < 1.0e-12:
                    nominal_rows.append(row)
                #endif
            #endfor
        #endfor

        table = pd.DataFrame(rows)
        table.to_csv(cfg_dir / "moradi_fit5_fit8_threshold_scan.csv", index=False)

        if len(table):
            fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.8))
            for fit_label, marker in [("Fit 5", "o"), ("Fit 8", "s")]:
                part = table.loc[
                    (table["fit_label"] == fit_label)
                    & table["success"].astype(bool)
                ]
                if len(part) == 0:
                    continue
                #endif
                axes[0].errorbar(
                    part["selection_threshold_percent"],
                    part["rE_fm"],
                    yerr=part["rE_fit_err_fm"],
                    marker=marker, linestyle="-", capsize=2,
                    label=fit_label,
                )
                axes[1].errorbar(
                    part["selection_threshold_percent"],
                    part["rM_fm"],
                    yerr=part["rM_fit_err_fm"],
                    marker=marker, linestyle="-", capsize=2,
                    label=(
                        fit_label if fit_label == "Fit 5"
                        else "Fit 8 (F2 fixed Kelly)"
                    ),
                )
                axes[2].plot(
                    part["selection_threshold_percent"],
                    part["chi2_ndof"],
                    marker=marker, linestyle="-",
                    label=fit_label,
                )
            #endfor
            axes[0].set_ylabel(r"$r_E$ [fm]")
            axes[1].set_ylabel(r"$r_M$ [fm]")
            axes[2].set_ylabel(r"$\chi^2/\mathrm{dof}$")
            for ax in axes:
                ax.set_xlabel(r"KM15 BH-purity threshold (%)")
                ax.axvline(5.0, linewidth=0.8, linestyle=":")
                ax.grid(alpha=0.2)
            #endfor
            handles, labels = axes[0].get_legend_handles_labels()
            if handles:
                fig.legend(
                    handles, labels, loc="upper center", ncol=2,
                    bbox_to_anchor=(0.5, 0.955),
                )
            #endif
            fig.suptitle(
                f"{label}: Moradi Fit 5 / Fit 8 benchmarks",
                y=0.995,
            )
            fig.tight_layout(rect=(0, 0, 1, 0.90))
            fig.savefig(cfg_dir / "01_moradi_fit5_fit8_threshold_scan.png", dpi=180)
            plt.close(fig)
        #endif
    #endfor

    all_table = pd.DataFrame(all_rows)
    all_table.to_csv(outdir / "moradi_fit5_fit8_all_thresholds.csv", index=False)

    nominal = pd.DataFrame(nominal_rows)
    nominal.to_csv(outdir / "moradi_fit5_fit8_5pct_summary.csv", index=False)

    if len(nominal):
        labels_order = [x[2] for x in configurations]
        fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.4))
        x = np.arange(len(labels_order), dtype=float)
        offsets = {"Fit 5": -0.12, "Fit 8": +0.12}
        for fit_label in ["Fit 5", "Fit 8"]:
            part = nominal.loc[
                (nominal["fit_label"] == fit_label)
                & nominal["success"].astype(bool)
            ].copy()
            part = part.set_index("configuration_label").reindex(labels_order)
            ok = np.isfinite(part["rE_fm"].to_numpy(float))
            axes[0].errorbar(
                x[ok] + offsets[fit_label],
                part.loc[np.asarray(labels_order)[ok], "rE_fm"],
                yerr=part.loc[np.asarray(labels_order)[ok], "rE_fit_err_fm"],
                marker="o" if fit_label == "Fit 5" else "s",
                linestyle="none", capsize=3,
                label=fit_label,
            )
            okm = np.isfinite(part["rM_fm"].to_numpy(float))
            axes[1].errorbar(
                x[okm] + offsets[fit_label],
                part.loc[np.asarray(labels_order)[okm], "rM_fm"],
                yerr=part.loc[np.asarray(labels_order)[okm], "rM_fit_err_fm"],
                marker="o" if fit_label == "Fit 5" else "s",
                linestyle="none", capsize=3,
                label=(
                    fit_label if fit_label == "Fit 5"
                    else "Fit 8 (F2 fixed Kelly)"
                ),
            )
        #endfor
        axes[0].set_ylabel(r"$r_E$ [fm]")
        axes[1].set_ylabel(r"$r_M$ [fm]")
        for ax in axes:
            ax.set_xticks(x)
            ax.set_xticklabels(labels_order, rotation=23, ha="right")
            ax.grid(axis="y", alpha=0.2)
        #endfor
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(
            handles, labels, loc="upper center", ncol=2,
            bbox_to_anchor=(0.5, 0.955),
        )
        fig.suptitle("KM15 5% Moradi Fit 5 / Fit 8 ensemble benchmarks", y=0.995)
        fig.tight_layout(rect=(0, 0, 1, 0.90))
        fig.savefig(outdir / "00_moradi_fit5_fit8_5pct_ensemble_summary.png", dpi=180)
        plt.close(fig)

        print("\n[Moradi ensemble benchmarks] KM15 5%")
        show = [
            "configuration_label", "fit_label", "N", "chi2_ndof",
            "rE_fm", "rE_fit_err_fm", "rM_fm", "rM_fit_err_fm",
        ]
        print(nominal[show].to_string(index=False))
    #endif
    return nominal
#enddef


def run_competitive_family_real_data_diagnostics(
        ensemble_defs: Sequence[Tuple[str, Sequence[Dict[str, object]], str]],
        chosen: Dict[str, Optional[str]],
        closure_root: Path,
        selection: pd.DataFrame,
        outdir: Path,
        max_families: int = 8) -> pd.DataFrame:
    """
    Fit the real KM15-5% data with the leading closure-qualified Sachs families.

    This measures the discrete real-data spread among families that closure
    regards as competitive.  It is intentionally diagnostic: the spread is
    saved and printed but is NOT folded into the production method uncertainty.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    all_rows = []
    spread_rows = []

    for tag, bundles, label in ensemble_defs:
        chosen_family = chosen.get(tag)
        ranking_candidates = [
            closure_root / "all_six_saylor_t343"
            / "radius_bias_mixed_family_ranking.csv",
            closure_root / "all_five"
            / "radius_bias_mixed_family_ranking.csv",
        ]
        ranking_path = next(
            (p for p in ranking_candidates if p.exists()),
            ranking_candidates[0],
        )
        if chosen_family is None or not ranking_path.exists():
            continue
        #endif
        ranking = pd.read_csv(ranking_path)
        usable = ranking.loc[
            ranking["eligible"].astype(bool)
            & np.isfinite(
                ranking["combined_RMS_objective_fm"].to_numpy(float)
            )
        ].copy()
        usable = usable.sort_values(
            "combined_RMS_objective_fm", ascending=True
        ).head(int(max_families)).reset_index(drop=True)
        specs, counts = _km15_selected_specs_for_bundles(
            bundles, selection, 0.05
        )

        rows = []
        for irank, rr in usable.iterrows():
            family = str(rr["family"])
            try:
                fit = fit_sachs_family_multi_measurements(
                    specs,
                    family=family,
                    bh_cut=0.05,
                    add_moradi_bh_systematic=True,
                    bh_systematic_fraction=0.05,
                )
                nobs = int(fit["N"])
                kpar = int(fit["n_parameters"])
                chi2 = float(fit["chi2"])
                aic = chi2 + 2.0 * kpar
                aicc = (
                    aic
                    + 2.0 * kpar * (kpar + 1.0)
                    / max(1.0, nobs - kpar - 1.0)
                )
                bic = chi2 + kpar * math.log(max(1, nobs))
                row = {
                    "configuration": tag,
                    "configuration_label": label,
                    "closure_rank": int(irank) + 1,
                    "family": family,
                    "chosen_family": chosen_family,
                    "is_chosen_family": bool(family == chosen_family),
                    "closure_objective_fm": float(
                        rr["combined_RMS_objective_fm"]
                    ),
                    "AIC": float(aic),
                    "AICc": float(aicc),
                    "BIC": float(bic),
                    **counts,
                    **fit,
                }
            except Exception as exc:
                row = {
                    "configuration": tag,
                    "configuration_label": label,
                    "closure_rank": int(irank) + 1,
                    "family": family,
                    "chosen_family": chosen_family,
                    "is_chosen_family": bool(family == chosen_family),
                    "closure_objective_fm": float(
                        rr["combined_RMS_objective_fm"]
                    ),
                    **counts,
                    "valid": False,
                    "failure": str(exc),
                    "rE_fm": np.nan,
                    "rE_fit_err_fm": np.nan,
                    "rM_fm": np.nan,
                    "rM_fit_err_fm": np.nan,
                    "chi2_ndof": np.nan,
                }
            #endtry
            rows.append(row)
            all_rows.append(row)
        #endfor

        table = pd.DataFrame(rows)
        valid = table.loc[
            table.get("valid", pd.Series(False, index=table.index)).astype(bool)
            & np.isfinite(table["rE_fm"].to_numpy(float))
            & np.isfinite(table["rM_fm"].to_numpy(float))
        ].copy()
        chosen_row = valid.loc[valid["family"] == chosen_family]
        if len(chosen_row) == 1:
            base_e = float(chosen_row.iloc[0]["rE_fm"])
            base_m = float(chosen_row.iloc[0]["rM_fm"])
            valid["delta_rE_from_chosen_fm"] = (
                valid["rE_fm"].to_numpy(float) - base_e
            )
            valid["delta_rM_from_chosen_fm"] = (
                valid["rM_fm"].to_numpy(float) - base_m
            )
            spread_rows.append({
                "configuration": tag,
                "configuration_label": label,
                "chosen_family": chosen_family,
                "N_competitive_valid": int(len(valid)),
                "max_abs_delta_rE_fm": float(np.max(np.abs(
                    valid["delta_rE_from_chosen_fm"].to_numpy(float)
                ))),
                "max_abs_delta_rM_fm": float(np.max(np.abs(
                    valid["delta_rM_from_chosen_fm"].to_numpy(float)
                ))),
                "rE_min_fm": float(valid["rE_fm"].min()),
                "rE_max_fm": float(valid["rE_fm"].max()),
                "rM_min_fm": float(valid["rM_fm"].min()),
                "rM_max_fm": float(valid["rM_fm"].max()),
                "production_systematic_eligible": False,
                "note": (
                    "Diagnostic spread among top closure-qualified families; "
                    "not folded into quoted method uncertainty."
                ),
            })
        #endif

        cfg_dir = outdir / tag
        cfg_dir.mkdir(parents=True, exist_ok=True)
        table.to_csv(cfg_dir / "competitive_family_real_data_fits.csv", index=False)

        if len(valid):
            fit_quality = valid[[
                "family", "closure_rank", "closure_objective_fm",
                "N", "n_parameters", "chi2", "ndof", "chi2_ndof",
                "AIC", "AICc", "BIC",
                "rE_fm", "rE_fit_err_fm", "rM_fm", "rM_fit_err_fm",
            ]].copy()
            for metric in ["AIC", "AICc", "BIC"]:
                fit_quality["delta_" + metric] = (
                    fit_quality[metric].to_numpy(float)
                    - float(fit_quality[metric].min())
                )
            #endfor
            fit_quality.to_csv(
                cfg_dir / "competitive_family_fit_quality.csv",
                index=False,
            )

            fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.0))
            xxq = np.arange(len(fit_quality))
            axes[0].plot(
                xxq, fit_quality["chi2_ndof"],
                marker="o", linestyle="none",
            )
            axes[0].set_ylabel(r"$\chi^2/\mathrm{d.o.f.}$")
            axes[1].plot(
                xxq, fit_quality["delta_AICc"],
                marker="o", linestyle="none", label=r"$\Delta$AICc",
            )
            axes[1].plot(
                xxq, fit_quality["delta_BIC"],
                marker="s", linestyle="none", label=r"$\Delta$BIC",
            )
            axes[1].axhline(0.0, linewidth=0.8, linestyle="--")
            axes[1].set_ylabel("Information-criterion difference")
            axes[1].legend(frameon=False)
            qlabs = [
                f"{f}\nclosure rank {r}"
                for f, r in zip(
                    fit_quality["family"], fit_quality["closure_rank"]
                )
            ]
            for ax in axes:
                ax.set_xticks(xxq)
                ax.set_xticklabels(qlabs, rotation=18, ha="right")
                ax.grid(axis="y", alpha=0.2)
            #endfor
            fig.suptitle(
                f"{label}: closure ranking vs real-data fit quality",
                y=0.995,
            )
            fig.tight_layout(rect=(0, 0, 1, 0.94))
            fig.savefig(
                cfg_dir / "03_competitive_family_fit_quality.png",
                dpi=280,
            )
            plt.close(fig)
            fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.0))
            xx = np.arange(len(valid))
            axes[0].errorbar(
                xx, valid["rE_fm"], yerr=valid["rE_fit_err_fm"],
                marker="o", linestyle="none", capsize=3,
            )
            axes[1].errorbar(
                xx, valid["rM_fm"], yerr=valid["rM_fit_err_fm"],
                marker="o", linestyle="none", capsize=3,
            )
            labs = [
                f"{f}\n(rank {r})"
                for f, r in zip(valid["family"], valid["closure_rank"])
            ]
            for ax, ylabel in zip(axes, [r"$r_E$ [fm]", r"$r_M$ [fm]"]):
                ax.set_xticks(xx)
                ax.set_xticklabels(labs, rotation=18, ha="right")
                ax.set_ylabel(ylabel)
                ax.grid(axis="y", alpha=0.2)
            #endfor
            fig.suptitle(
                f"{label}: real-data spread among top closure-qualified families",
                y=0.995,
            )
            fig.tight_layout(rect=(0, 0, 1, 0.94))
            fig.savefig(cfg_dir / "01_competitive_family_real_data_radii.png", dpi=260)
            plt.close(fig)

            # Radius spread can be much larger than finite-|t| disagreement.
            # Quantify that distinction with full Hessian GE/GM bands.
            save_competitive_family_form_factor_bands(
                valid=valid,
                chosen_family=chosen_family,
                label=label,
                selected_specs=specs,
                outdir=cfg_dir,
            )
        #endif
    #endfor

    all_table = pd.DataFrame(all_rows)
    all_table.to_csv(outdir / "competitive_family_real_data_fits_all.csv", index=False)
    spread = pd.DataFrame(spread_rows)
    spread.to_csv(outdir / "competitive_family_spread_summary.csv", index=False)
    if len(spread):
        print("\n[competitive-family real-data diagnostic]")
        print(spread[[
            "configuration_label", "chosen_family",
            "max_abs_delta_rE_fm", "max_abs_delta_rM_fm",
        ]].to_string(index=False))
    #endif
    return spread
#enddef




def _closure_bootstrap_win_probabilities(
        closure_table: pd.DataFrame,
        families: Sequence[str],
        n_bootstrap: int = 20000,
        seed: int = 20260903) -> pd.DataFrame:
    """
    Bootstrap conceptual truth groups and within-group truth realizations.

    For every bootstrap draw, determine which candidate family minimizes the
    same combined E/M closure RMSE objective.  The resulting win fractions are
    a direct measure of how stable the discrete family ranking is to the finite
    truth ensemble, rather than interpreting tiny objective differences as
    deterministic.
    """
    fams = [str(f) for f in families]
    part = closure_table.loc[
        closure_table["family"].astype(str).isin(fams)
    ].copy()
    part = production_closure_truth_rows(part)
    groups = list(dict.fromkeys(part["truth_group"].astype(str).tolist()))
    if not groups:
        raise RuntimeError(
            "Truth-menu resampling has zero controlled synthetic groups."
        )
    #endif
    rng = np.random.default_rng(int(seed))

    cache = {}
    for fam in fams:
        cache[fam] = {}
        fp = part.loc[part["family"].astype(str) == fam]
        for group in groups:
            gp = fp.loc[fp["truth_group"].astype(str) == group]
            truths = list(dict.fromkeys(gp["truth_model"].astype(str).tolist()))
            vals = []
            for truth in truths:
                tp = gp.loc[gp["truth_model"].astype(str) == truth]
                erow = tp.loc[tp["quantity"].astype(str) == "rE"]
                mrow = tp.loc[tp["quantity"].astype(str) == "rM"]
                if len(erow) == 1 and len(mrow) == 1:
                    vals.append((
                        float(erow.iloc[0]["sqrt_stat2_plus_bias2_fm"]),
                        float(mrow.iloc[0]["sqrt_stat2_plus_bias2_fm"]),
                    ))
                #endif
            #endfor
            cache[fam][group] = vals
        #endfor
    #endfor

    wins = {fam: 0 for fam in fams}
    score_sum = {fam: 0.0 for fam in fams}
    score_sq = {fam: 0.0 for fam in fams}
    valid_draws = 0

    for _ in range(max(1, int(n_bootstrap))):
        drawn_groups = rng.choice(groups, size=len(groups), replace=True)
        scores = {}
        for fam in fams:
            e2 = []
            m2 = []
            okay = True
            for group in drawn_groups:
                vals = cache[fam].get(str(group), [])
                if not vals:
                    okay = False
                    break
                #endif
                e, m = vals[int(rng.integers(0, len(vals)))]
                if not np.isfinite(e) or not np.isfinite(m):
                    okay = False
                    break
                #endif
                e2.append(e * e)
                m2.append(m * m)
            #endfor
            if okay and e2 and m2:
                scores[fam] = float(np.sqrt(
                    0.5 * (float(np.mean(e2)) + float(np.mean(m2)))
                ))
            #endif
        #endfor

        if not scores:
            continue
        #endif
        valid_draws += 1
        winner = min(scores, key=scores.get)
        wins[winner] += 1
        for fam, score in scores.items():
            score_sum[fam] += score
            score_sq[fam] += score * score
        #endfor
    #endfor

    rows = []
    for fam in fams:
        mean = (
            score_sum[fam] / valid_draws if valid_draws > 0 else np.nan
        )
        var = (
            score_sq[fam] / valid_draws - mean**2
            if valid_draws > 0 and np.isfinite(mean) else np.nan
        )
        rows.append({
            "family": fam,
            "closure_bootstrap_wins": int(wins[fam]),
            "closure_bootstrap_valid_draws": int(valid_draws),
            "closure_bootstrap_win_fraction": (
                wins[fam] / valid_draws if valid_draws > 0 else np.nan
            ),
            "closure_bootstrap_mean_objective_fm": mean,
            "closure_bootstrap_std_objective_fm": (
                float(np.sqrt(max(0.0, var)))
                if np.isfinite(var) else np.nan
            ),
        })
    #endfor
    return pd.DataFrame(rows).sort_values(
        "closure_bootstrap_win_fraction", ascending=False
    ).reset_index(drop=True)
#enddef


def _mixture_density_and_interval(
        centers: np.ndarray,
        sigmas: np.ndarray,
        weights: np.ndarray,
        lo: float,
        hi: float,
        ngrid: int = 5000) -> Tuple[np.ndarray, np.ndarray, Dict[str, float]]:
    """Normalized Gaussian-mixture density plus central 68.3% interval."""
    x = np.linspace(float(lo), float(hi), int(ngrid))
    p = np.zeros_like(x)
    for mu, sig, w in zip(centers, sigmas, weights):
        sig = max(float(sig), 1.0e-4)
        p += float(w) * np.exp(-0.5 * ((x - float(mu)) / sig)**2) / (
            np.sqrt(2.0 * np.pi) * sig
        )
    #endfor
    norm = float(np.trapezoid(p, x))
    if norm > 0.0:
        p /= norm
    #endif
    dx = x[1] - x[0]
    cdf = np.cumsum(p) * dx
    cdf = np.clip(cdf, 0.0, 1.0)
    qlo = float(np.interp(0.1585, cdf, x))
    qmed = float(np.interp(0.5000, cdf, x))
    qhi = float(np.interp(0.8415, cdf, x))
    mean = float(np.trapezoid(x * p, x))
    mode_index = int(np.argmax(p))
    mode = float(x[mode_index])

    target = 0.683
    best_lo = 0
    best_hi = len(x) - 1
    best_width = np.inf
    j = 0
    for i in range(len(x)):
        if j < i:
            j = i
        #endif
        cdf_before = cdf[i - 1] if i > 0 else 0.0
        while j < len(x) - 1 and (cdf[j] - cdf_before) < target:
            j += 1
        #endwhile
        mass = cdf[j] - cdf_before
        if mass >= target and i <= mode_index <= j:
            width = x[j] - x[i]
            if width < best_width:
                best_width = width
                best_lo = i
                best_hi = j
            #endif
        #endif
    #endfor
    hdi_lo = float(x[best_lo])
    hdi_hi = float(x[best_hi])

    return x, p, {
        "mean_fm": mean,
        "median_fm": qmed,
        "mode_fm": mode,
        "central68_low_fm": qlo,
        "central68_high_fm": qhi,
        "central68_minus_fm": qmed - qlo,
        "central68_plus_fm": qhi - qmed,
        "mode_HDI68_low_fm": hdi_lo,
        "mode_HDI68_high_fm": hdi_hi,
        "mode_HDI68_minus_fm": mode - hdi_lo,
        "mode_HDI68_plus_fm": hdi_hi - mode,
    }
#enddef


def save_hayward_griffioen_model_averaging_diagnostics(
        closure_table_path: Path,
        closure_ranking_path: Path,
        competitive_fit_path: Path,
        selected_specs: Sequence[Dict[str, object]],
        outdir: Path,
        max_families: int = 8,
        n_bootstrap: int = 20000) -> pd.DataFrame:
    """
    Modernized Hayward-Griffioen-style P(r) construction.

    1. Closure determines the admissible leading model set.
    2. Truth-group bootstrap win fractions quantify ranking stability.
    3. Real-data AICc enters only AFTER closure as an empirical-adequacy weight.
    4. Each model Gaussian is centered at r_fit - <signed closure bias>.
    5. Its width combines Hessian fit uncertainty with the residual closure
       bias RMS after removal of the signed mean bias.
    6. Both closure-only and closure×AICc mixtures are saved, so the impact of
       incorporating real-data goodness-of-fit is explicit rather than hidden
       in a single arbitrary scalar objective.
    7. A companion 2x2 diagnostic plots every candidate family's weighted
       contribution underneath the summed P(r), separately for rE/rM and for
       the closure-only versus closure×AICc weighting prescriptions.

    This follows the spirit of Eq. (33) of Hayward & Griffioen (2018), while
    making the closure and empirical-fit contributions separately auditable.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    closure = pd.read_csv(closure_table_path)
    ranking = pd.read_csv(closure_ranking_path)
    fits = pd.read_csv(competitive_fit_path)

    usable = ranking.loc[
        ranking["eligible"].astype(bool)
        & np.isfinite(ranking["combined_RMS_objective_fm"].to_numpy(float))
    ].sort_values("combined_RMS_objective_fm").head(int(max_families)).copy()

    fams = usable["family"].astype(str).tolist()
    fits = fits.loc[
        fits["family"].astype(str).isin(fams)
        & fits["valid"].astype(bool)
    ].copy()
    fams = [f for f in fams if f in set(fits["family"].astype(str))]
    usable = usable.loc[usable["family"].astype(str).isin(fams)].copy()
    if not fams:
        return pd.DataFrame()
    #endif

    boot = _closure_bootstrap_win_probabilities(
        closure, fams, n_bootstrap=n_bootstrap
    )
    model = fits.merge(
        usable[[
            "family", "combined_RMS_objective_fm",
            "rE_RMS_bias_fm", "rM_RMS_bias_fm",
            "rE_mean_signed_bias_fm", "rM_mean_signed_bias_fm",
            "rE_residual_bias_RMS_after_mean_correction_fm",
            "rM_residual_bias_RMS_after_mean_correction_fm",
        ]],
        on="family", how="inner",
    ).merge(boot, on="family", how="left")

    nobs = model["N"].to_numpy(float)
    kpar = model["n_parameters"].to_numpy(float)
    chi2 = model["chi2"].to_numpy(float)
    model["AICc"] = (
        chi2 + 2.0 * kpar
        + 2.0 * kpar * (kpar + 1.0)
        / np.maximum(1.0, nobs - kpar - 1.0)
    )
    delta_aicc = model["AICc"].to_numpy(float) - float(model["AICc"].min())
    model["akaike_relative_support"] = np.exp(-0.5 * delta_aicc)

    closure_w = model["closure_bootstrap_win_fraction"].to_numpy(float)
    # A model may never be the strict winner while remaining essentially tied.
    # Add a tiny floor only to keep diagnostic mixture components visible.
    closure_w = np.maximum(closure_w, 1.0e-6)
    closure_w /= np.sum(closure_w)
    model["weight_closure_only"] = closure_w

    combined_w = closure_w * model["akaike_relative_support"].to_numpy(float)
    if np.sum(combined_w) <= 0.0:
        combined_w = closure_w.copy()
    #endif
    combined_w /= np.sum(combined_w)
    model["weight_closure_times_AICc"] = combined_w

    for quantity in ["rE", "rM"]:
        signed = model[f"{quantity}_mean_signed_bias_fm"].to_numpy(float)
        model[f"{quantity}_bias_corrected_center_fm"] = (
            model[f"{quantity}_fm"].to_numpy(float) - signed
        )
        residual = model[
            f"{quantity}_residual_bias_RMS_after_mean_correction_fm"
        ].to_numpy(float)
        model[f"{quantity}_mixture_sigma_fm"] = np.hypot(
            model[f"{quantity}_fit_err_fm"].to_numpy(float),
            np.nan_to_num(residual, nan=0.0),
        )
    #endfor

    model.to_csv(outdir / "radius_model_components_and_weights.csv", index=False)
    boot.to_csv(outdir / "closure_bootstrap_family_win_probabilities.csv", index=False)

    summary_rows = []
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.1))
    for ax, quantity, label in [
        (axes[0], "rE", r"$r_E$ (fm)"),
        (axes[1], "rM", r"$r_M$ (fm)"),
    ]:
        centers = model[f"{quantity}_bias_corrected_center_fm"].to_numpy(float)
        sigmas = model[f"{quantity}_mixture_sigma_fm"].to_numpy(float)
        lo = max(0.20, float(np.nanmin(centers - 5.0 * sigmas)) - 0.02)
        hi = min(1.50, float(np.nanmax(centers + 5.0 * sigmas)) + 0.02)

        for weight_col, tag, linestyle in [
            ("weight_closure_only", "closure bootstrap only", "--"),
            ("weight_closure_times_AICc", "closure × AICc", "-"),
        ]:
            x, p, stats = _mixture_density_and_interval(
                centers, sigmas, model[weight_col].to_numpy(float), lo, hi
            )
            ax.plot(x, p / np.max(p), linestyle=linestyle, label=tag)
            for key, val in stats.items():
                summary_rows.append({
                    "quantity": quantity,
                    "mixture": tag,
                    "statistic": key,
                    "value_fm": val,
                })
            #endfor
            pd.DataFrame({
                "radius_fm": x,
                "density": p,
                "mixture": tag,
                "quantity": quantity,
            }).to_csv(
                outdir / f"P_{quantity}_{weight_col}.csv",
                index=False,
            )
        #endfor
        ax.set_xlabel(label)
        ax.set_ylabel(r"$P(r)$ (normalized to peak)")
        ax.set_xlim(0.70, 1.00)
        ax.grid(alpha=0.2)
        ax.legend(frameon=False)
    #endfor
    fig.suptitle(
        "Model-averaged radius distributions: closure stability + fit adequacy",
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(outdir / "01_model_averaged_radius_probabilities.png", dpi=280)
    plt.close(fig)

    # ------------------------------------------------------------------
    # Component-resolved P(r) diagnostic.
    #
    # The main 01 plot intentionally shows only the two final mixtures.  This
    # companion figure exposes exactly how those smooth distributions are
    # assembled from the individual candidate-family contributions.
    #
    # Each thin curve below is WEIGHTED before plotting:
    #
    #     contribution_i(r) = w_i * Gaussian_i(r)
    #
    # so its visible area reflects the actual weight carried by that family.
    # The thick curve is the sum of those weighted components.  All curves in
    # each panel are divided by the peak of the summed distribution, preserving
    # the relative importance of the components while keeping a common vertical
    # scale.
    # ------------------------------------------------------------------
    fig_components, axes_components = plt.subplots(
        2, 2, figsize=(14.2, 10.0), sharey=False
    )
    component_rows = []

    weight_specs = [
        (
            "weight_closure_only",
            "Known-answer tests only",
            "closure bootstrap only",
        ),
        (
            "weight_closure_times_AICc",
            "Known-answer tests + real-data AICc",
            "closure × AICc",
        ),
    ]

    for irow, (weight_col, row_title, mixture_tag) in enumerate(weight_specs):
        for icol, (quantity, xlabel) in enumerate([
            ("rE", r"$r_E$ (fm)"),
            ("rM", r"$r_M$ (fm)"),
        ]):
            ax = axes_components[irow, icol]
            centers = model[
                f"{quantity}_bias_corrected_center_fm"
            ].to_numpy(float)
            sigmas = model[
                f"{quantity}_mixture_sigma_fm"
            ].to_numpy(float)
            weights = model[weight_col].to_numpy(float)

            lo = max(
                0.20,
                float(np.nanmin(centers - 5.0 * sigmas)) - 0.02,
            )
            hi = min(
                1.50,
                float(np.nanmax(centers + 5.0 * sigmas)) + 0.02,
            )
            x = np.linspace(lo, hi, 5000)

            weighted_components = []
            for (_, row), mu, sig, weight in zip(
                    model.iterrows(), centers, sigmas, weights):
                sig = max(float(sig), 1.0e-4)
                weight = float(weight)
                component = (
                    weight
                    * np.exp(-0.5 * ((x - float(mu)) / sig)**2)
                    / (np.sqrt(2.0 * np.pi) * sig)
                )
                weighted_components.append(component)
            #endfor

            weighted_components = np.asarray(weighted_components)
            total = np.sum(weighted_components, axis=0)
            peak = float(np.nanmax(total))
            if not np.isfinite(peak) or peak <= 0.0:
                peak = 1.0
            #endif

            # Plot individual weighted family contributions in decreasing
            # weight order so the legend is immediately interpretable.
            order = np.argsort(weights)[::-1]
            for idx in order:
                family = str(model.iloc[int(idx)]["family"])
                weight = float(weights[int(idx)])
                contribution = weighted_components[int(idx)] / peak
                ax.plot(
                    x,
                    contribution,
                    linewidth=1.15,
                    alpha=0.72,
                    label=f"{family}  ({100.0 * weight:.1f}%)",
                )

                for xx, yy in zip(x, contribution):
                    component_rows.append({
                        "quantity": quantity,
                        "mixture": mixture_tag,
                        "family": family,
                        "family_weight": weight,
                        "radius_fm": float(xx),
                        "weighted_component_normalized_to_total_peak": float(yy),
                    })
                #endfor
            #endfor

            total_norm = total / peak
            ax.plot(
                x,
                total_norm,
                linewidth=2.8,
                label="summed $P(r)$",
            )
            for xx, yy in zip(x, total_norm):
                component_rows.append({
                    "quantity": quantity,
                    "mixture": mixture_tag,
                    "family": "__SUM__",
                    "family_weight": 1.0,
                    "radius_fm": float(xx),
                    "weighted_component_normalized_to_total_peak": float(yy),
                })
            #endfor

            ax.set_xlabel(xlabel)
            ax.set_ylabel(r"Contribution to $P(r)$ (sum peak = 1)")
            ax.set_xlim(0.70, 1.00)
            ax.set_title(
                f"{row_title}: "
                + (r"$r_E$" if quantity == "rE" else r"$r_M$")
            )
            ax.grid(alpha=0.20)
            ax.legend(
                fontsize=8.0,
                frameon=False,
                ncol=1,
                loc="upper right",
            )
        #endfor
    #endfor

    fig_components.suptitle(
        "Radius-probability decomposition by candidate fit family",
        y=0.995,
    )
    fig_components.tight_layout(rect=(0, 0, 1, 0.965))
    fig_components.savefig(
        outdir / "03_radius_probability_family_components.png",
        dpi=300,
    )
    plt.close(fig_components)

    pd.DataFrame(component_rows).to_csv(
        outdir / "radius_probability_family_component_curves.csv",
        index=False,
    )

    # Compact table containing exactly the numbers needed to interpret the
    # component plot without inspecting the full competitive-fit table.
    component_summary = model[[
        "family",
        "combined_RMS_objective_fm",
        "closure_bootstrap_win_fraction",
        "AICc",
        "akaike_relative_support",
        "weight_closure_only",
        "weight_closure_times_AICc",
        "rE_fm",
        "rE_fit_err_fm",
        "rE_mean_signed_bias_fm",
        "rE_bias_corrected_center_fm",
        "rE_mixture_sigma_fm",
        "rM_fm",
        "rM_fit_err_fm",
        "rM_mean_signed_bias_fm",
        "rM_bias_corrected_center_fm",
        "rM_mixture_sigma_fm",
    ]].copy()
    component_summary = component_summary.sort_values(
        "weight_closure_only", ascending=False
    )
    component_summary.to_csv(
        outdir / "radius_probability_family_component_summary.csv",
        index=False,
    )

    # Primary interpolation result: model-average GE and GM over 0 <= Q2 <= 1.
    q = np.linspace(0.0, 1.0, 900)
    qdata = np.concatenate([
        spec["data"]["t_abs"].to_numpy(float)
        for spec in selected_specs if len(spec["data"])
    ])
    qlo = float(np.nanmin(qdata))
    qhi = float(np.nanmax(qdata))
    interpolation_rows = []

    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.1))
    for ax, which, ylabel in [
        (axes[0], "GE", r"$G_E^p$"),
        (axes[1], "GM", r"$G_M^p$"),
    ]:
        curves = []
        variances = []
        weights = model["weight_closure_times_AICc"].to_numpy(float)
        for _, row in model.iterrows():
            fam = str(row["family"])
            result = row.to_dict()
            central, err = _sachs_band_from_result(result, fam, q, which)
            curves.append(np.asarray(central, dtype=float))
            variances.append(np.asarray(err, dtype=float)**2)
        #endfor
        curves = np.asarray(curves)
        variances = np.asarray(variances)
        mean_curve = np.sum(weights[:, None] * curves, axis=0)
        within = np.sum(weights[:, None] * variances, axis=0)
        between = np.sum(
            weights[:, None] * (curves - mean_curve[None, :])**2,
            axis=0,
        )
        total = np.sqrt(np.maximum(0.0, within + between))
        ax.plot(q, mean_curve, linewidth=2.0, label="model-averaged BH extraction")
        ax.fill_between(
            q, mean_curve - total, mean_curve + total,
            alpha=0.20, linewidth=0.0,
            label="within + between-model 68% band",
        )
        ax.axvspan(qlo, qhi, alpha=0.07, label="selected-data support")
        ax.set_xlim(0.0, 1.0)
        ax.set_xlabel(r"$Q^2=|t|$ [GeV$^2$]")
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.2)
        ax.legend(fontsize=8.5, frameon=False)
        for iq, qq in enumerate(q):
            interpolation_rows.append({
                "form_factor": which,
                "Q2_GeV2": float(qq),
                "model_average": float(mean_curve[iq]),
                "within_model_sigma": float(np.sqrt(max(0.0, within[iq]))),
                "between_model_sigma": float(np.sqrt(max(0.0, between[iq]))),
                "total_model_averaged_sigma": float(total[iq]),
                "inside_selected_data_support": bool(qlo <= qq <= qhi),
            })
        #endfor
    #endfor
    fig.suptitle(
        "Primary finite-$Q^2$ result: closure/AICc model-averaged form factors",
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(outdir / "02_model_averaged_GE_GM_interpolation.png", dpi=280)
    plt.close(fig)
    pd.DataFrame(interpolation_rows).to_csv(
        outdir / "model_averaged_GE_GM_interpolation.csv", index=False
    )
    radius_summary = pd.DataFrame(summary_rows)
    radius_summary.to_csv(
        outdir / "radius_probability_summary.csv", index=False
    )

    preferred_rows = []
    for quantity in ["rE", "rM"]:
        s = radius_summary.loc[
            (radius_summary["quantity"] == quantity)
            & (radius_summary["mixture"] == "closure bootstrap only")
        ].set_index("statistic")["value_fm"]
        preferred_rows.append({
            "quantity": quantity,
            "preferred_prescription": "known_answer_tests_only",
            "central_estimator": "mode",
            "radius_fm": float(s.get("mode_fm", np.nan)),
            "minus_68p3_fm": float(s.get("mode_HDI68_minus_fm", np.nan)),
            "plus_68p3_fm": float(s.get("mode_HDI68_plus_fm", np.nan)),
            "interval_low_fm": float(s.get("mode_HDI68_low_fm", np.nan)),
            "interval_high_fm": float(s.get("mode_HDI68_high_fm", np.nan)),
            "interval_definition": (
                "shortest contiguous 68.3% probability interval containing "
                "the mode; candidate-family centers are signed-bias corrected"
            ),
        })
    #endfor
    pd.DataFrame(preferred_rows).to_csv(
        outdir / "preferred_known_answer_model_averaged_radii.csv",
        index=False,
    )

    print("\n[Hayward-Griffioen-style model averaging]")
    print(model[[
        "family", "combined_RMS_objective_fm",
        "closure_bootstrap_win_fraction", "AICc",
        "weight_closure_only", "weight_closure_times_AICc",
        "rE_fm", "rE_mean_signed_bias_fm",
        "rE_bias_corrected_center_fm",
        "rM_fm", "rM_mean_signed_bias_fm",
        "rM_bias_corrected_center_fm",
    ]].sort_values(
        "weight_closure_times_AICc", ascending=False
    ).to_string(index=False))
    return model
#enddef


def fit_sachs_ge_family_fixed_kelly_f2_multi_measurements(
        datasets: Sequence[Dict[str, object]],
        ge_family: str,
        bh_cut: float = 0.05,
        add_moradi_bh_systematic: bool = True
        ) -> Dict[str, object]:
    """
    Global BH fit with a flexible GE family while the Pauli form factor F2 is
    fixed point-by-point to Kelly 2004.

    This is the direct multi-dataset analogue of the Moradi Fit-8 idea, but it
    keeps the same GE family used in the closure machinery.  It isolates how
    much of the radius instability comes from simultaneous F1/F2 (GE/GM)
    freedom rather than from GE extrapolation itself.
    """
    ge_family = str(ge_family)
    ne = int(re.findall(r"\d+", ge_family)[0])
    names_e = [f"e{i}" for i in range(1, ne + 1)]
    p0 = np.zeros(ne, dtype=float)
    p0[0] = sachs_first_coefficient_from_radius(
        SACHS_INITIAL_RADIUS_FM, ge_family
    )

    nuisance_names = []
    nuisance_fracs = {}
    nuisance_is_free = {}
    for spec in datasets:
        key = str(spec["key"])
        free_norm = bool(spec.get("unconstrained_norm", False))
        frac = 1.0 if free_norm else float(spec.get("norm_frac", 0.0))
        if frac > 0.0:
            nname = "beta_" + re.sub(r"[^A-Za-z0-9]+", "_", key)
            nuisance_names.append(nname)
            nuisance_fracs[nname] = frac
            nuisance_is_free[nname] = free_norm
        #endif
    #endfor

    fit_names = names_e + nuisance_names
    fit_p0 = np.concatenate([p0, np.zeros(len(nuisance_names))])
    nuisance_index = {
        n: ne + i for i, n in enumerate(nuisance_names)
    }

    prepared = []
    for spec in datasets:
        d = spec["data"]
        if len(d) == 0:
            continue
        #endif
        err = dataset_point_errors(
            d, str(spec["kind"]), bh_cut, add_moradi_bh_systematic
        )
        q = d["t_abs"].to_numpy(float)
        tau = q / (4.0 * MP2)
        _, f2_kelly = kelly_f1_f2(q)
        key = str(spec["key"])
        free_norm = bool(spec.get("unconstrained_norm", False))
        frac = 1.0 if free_norm else float(spec.get("norm_frac", 0.0))
        prepared.append({
            "key": key,
            "q": q,
            "q_powers": np.vstack([q**i for i in range(1, ne + 1)]),
            "tau": tau,
            "f2_kelly": np.asarray(f2_kelly, dtype=float),
            "y": d["xs"].to_numpy(float),
            "e": np.asarray(err, dtype=float),
            "A": d["bh_A"].to_numpy(float),
            "B": d["bh_B"].to_numpy(float),
            "C": d["bh_C"].to_numpy(float),
            "N": int(len(d)),
            "norm_frac": frac,
            "nuisance_name": (
                "beta_" + re.sub(r"[^A-Za-z0-9]+", "_", key)
                if frac > 0.0 else None
            ),
        })
    #endfor
    if not prepared:
        raise RuntimeError("Kelly-F2 global fit has zero selected points.")
    #endif

    def chi2_minuit(*values):
        p = np.asarray(values, dtype=float)
        ce = p[:ne]
        total = 0.0
        for item in prepared:
            ge = sachs_family_value_precomputed(
                item["q"], ce, ge_family, item["q_powers"]
            )
            f2 = item["f2_kelly"]
            f1 = ge + item["tau"] * f2
            pred = bh_from_f1f2(
                item["A"], item["B"], item["C"], f1, f2
            )
            if item["norm_frac"] > 0.0:
                beta = p[nuisance_index[item["nuisance_name"]]]
                pred *= 1.0 + item["norm_frac"] * beta
            #endif
            pull = (pred - item["y"]) / item["e"]
            total += float(np.dot(pull, pull))
        #endfor
        for nname in nuisance_names:
            if not nuisance_is_free.get(nname, False):
                total += float(p[nuisance_index[nname]]**2)
            #endif
        #endfor
        return total
    #enddef

    mn = Minuit(chi2_minuit, *fit_p0, name=tuple(fit_names))
    mn.errordef = Minuit.LEAST_SQUARES
    ce_lo = sachs_first_coefficient_from_radius(
        SACHS_MIN_RADIUS_FM, ge_family
    )
    ce_hi = sachs_first_coefficient_from_radius(
        SACHS_MAX_RADIUS_FM, ge_family
    )
    mn.limits[names_e[0]] = (min(ce_lo, ce_hi), max(ce_lo, ce_hi))
    for nname in nuisance_names:
        mn.limits[nname] = (
            (-0.50, 0.50)
            if nuisance_is_free.get(nname, False)
            else (-10.0, 10.0)
        )
    #endfor
    mn.migrad()
    mn.hesse()

    shape = np.array([float(mn.values[n]) for n in names_e])
    cov = np.full((ne, ne), np.nan)
    if mn.covariance is not None:
        for i, ni in enumerate(names_e):
            for j, nj in enumerate(names_e):
                cov[i, j] = float(mn.covariance[ni, nj])
            #endfor
        #endfor
    #endif

    def radius_e(pars):
        return sachs_family_radius(np.asarray(pars), ge_family)
    #enddef

    def radius_m_derived(pars):
        pars = np.asarray(pars, dtype=float)
        def gm_shape(q):
            q = np.asarray(q, dtype=float)
            ge = sachs_family_value(q, pars, ge_family)
            _, f2 = kelly_f1_f2(q)
            tau = q / (4.0 * MP2)
            return ge + (1.0 + tau) * f2
        #enddef
        return radius_from_shape(gm_shape, MU_P)
    #enddef

    rE, rEerr = propagate_scalar(radius_e, shape, cov)
    rM, rMerr = propagate_scalar(radius_m_derived, shape, cov)

    ndata = int(sum(item["N"] for item in prepared))
    nfree = int(len(fit_names))
    ndof = max(1, ndata - nfree)
    result = {
        "model": f"GE:{ge_family}|F2:Kelly2004",
        "GE_family": ge_family,
        "F2_treatment": "Kelly 2004 fixed",
        "N": ndata,
        "n_parameters": nfree,
        "chi2": float(mn.fval),
        "ndof": int(ndof),
        "chi2_ndof": float(mn.fval / ndof),
        "valid": bool(
            mn.valid and np.isfinite(rE) and np.isfinite(rM)
        ),
        "rE_fm": float(rE),
        "rE_fit_err_fm": float(rEerr),
        "rM_fm": float(rM),
        "rM_fit_err_fm": float(rMerr),
        "shape_covariance_json": json.dumps(cov.tolist()),
    }
    for name in names_e:
        result[name] = float(mn.values[name])
        result[name + "_err"] = float(mn.errors[name])
    #endfor
    for nname in nuisance_names:
        result[nname] = float(mn.values[nname])
        result[nname + "_scale_factor"] = float(
            1.0 + nuisance_fracs[nname] * float(mn.values[nname])
        )
    #endfor
    return result
#enddef

def fit_sachs_ge_family_fixed_yahl18_f2_multi_measurements(
        datasets: Sequence[Dict[str, object]],
        ge_family: str,
        bh_cut: float = 0.05,
        add_moradi_bh_systematic: bool = True
        ) -> Dict[str, object]:
    """
    Global BH fit with a flexible GE family while the Pauli form factor F2 is
    fixed point-by-point to YAHL18.

    This is the direct multi-dataset analogue of the Moradi Fit-8 idea, but it
    keeps the same GE family used in the closure machinery.  It isolates how
    much of the radius instability comes from simultaneous F1/F2 (GE/GM)
    freedom rather than from GE extrapolation itself.
    """
    ge_family = str(ge_family)
    ne = int(re.findall(r"\d+", ge_family)[0])
    names_e = [f"e{i}" for i in range(1, ne + 1)]
    p0 = np.zeros(ne, dtype=float)
    p0[0] = sachs_first_coefficient_from_radius(
        SACHS_INITIAL_RADIUS_FM, ge_family
    )

    nuisance_names = []
    nuisance_fracs = {}
    nuisance_is_free = {}
    for spec in datasets:
        key = str(spec["key"])
        free_norm = bool(spec.get("unconstrained_norm", False))
        frac = 1.0 if free_norm else float(spec.get("norm_frac", 0.0))
        if frac > 0.0:
            nname = "beta_" + re.sub(r"[^A-Za-z0-9]+", "_", key)
            nuisance_names.append(nname)
            nuisance_fracs[nname] = frac
            nuisance_is_free[nname] = free_norm
        #endif
    #endfor

    fit_names = names_e + nuisance_names
    fit_p0 = np.concatenate([p0, np.zeros(len(nuisance_names))])
    nuisance_index = {
        n: ne + i for i, n in enumerate(nuisance_names)
    }

    prepared = []
    for spec in datasets:
        d = spec["data"]
        if len(d) == 0:
            continue
        #endif
        err = dataset_point_errors(
            d, str(spec["kind"]), bh_cut, add_moradi_bh_systematic
        )
        q = d["t_abs"].to_numpy(float)
        tau = q / (4.0 * MP2)
        _, f2_yahl18 = yahl18_f1_f2(q)
        key = str(spec["key"])
        free_norm = bool(spec.get("unconstrained_norm", False))
        frac = 1.0 if free_norm else float(spec.get("norm_frac", 0.0))
        prepared.append({
            "key": key,
            "q": q,
            "q_powers": np.vstack([q**i for i in range(1, ne + 1)]),
            "tau": tau,
            "f2_yahl18": np.asarray(f2_yahl18, dtype=float),
            "y": d["xs"].to_numpy(float),
            "e": np.asarray(err, dtype=float),
            "A": d["bh_A"].to_numpy(float),
            "B": d["bh_B"].to_numpy(float),
            "C": d["bh_C"].to_numpy(float),
            "N": int(len(d)),
            "norm_frac": frac,
            "nuisance_name": (
                "beta_" + re.sub(r"[^A-Za-z0-9]+", "_", key)
                if frac > 0.0 else None
            ),
        })
    #endfor
    if not prepared:
        raise RuntimeError("YAHL18-F2 global fit has zero selected points.")
    #endif

    def chi2_minuit(*values):
        p = np.asarray(values, dtype=float)
        ce = p[:ne]
        total = 0.0
        for item in prepared:
            ge = sachs_family_value_precomputed(
                item["q"], ce, ge_family, item["q_powers"]
            )
            f2 = item["f2_yahl18"]
            f1 = ge + item["tau"] * f2
            pred = bh_from_f1f2(
                item["A"], item["B"], item["C"], f1, f2
            )
            if item["norm_frac"] > 0.0:
                beta = p[nuisance_index[item["nuisance_name"]]]
                pred *= 1.0 + item["norm_frac"] * beta
            #endif
            pull = (pred - item["y"]) / item["e"]
            total += float(np.dot(pull, pull))
        #endfor
        for nname in nuisance_names:
            if not nuisance_is_free.get(nname, False):
                total += float(p[nuisance_index[nname]]**2)
            #endif
        #endfor
        return total
    #enddef

    mn = Minuit(chi2_minuit, *fit_p0, name=tuple(fit_names))
    mn.errordef = Minuit.LEAST_SQUARES
    ce_lo = sachs_first_coefficient_from_radius(
        SACHS_MIN_RADIUS_FM, ge_family
    )
    ce_hi = sachs_first_coefficient_from_radius(
        SACHS_MAX_RADIUS_FM, ge_family
    )
    mn.limits[names_e[0]] = (min(ce_lo, ce_hi), max(ce_lo, ce_hi))
    for nname in nuisance_names:
        mn.limits[nname] = (
            (-0.50, 0.50)
            if nuisance_is_free.get(nname, False)
            else (-10.0, 10.0)
        )
    #endfor
    mn.migrad()
    mn.hesse()

    shape = np.array([float(mn.values[n]) for n in names_e])
    cov = np.full((ne, ne), np.nan)
    if mn.covariance is not None:
        for i, ni in enumerate(names_e):
            for j, nj in enumerate(names_e):
                cov[i, j] = float(mn.covariance[ni, nj])
            #endfor
        #endfor
    #endif

    def radius_e(pars):
        return sachs_family_radius(np.asarray(pars), ge_family)
    #enddef

    def radius_m_derived(pars):
        pars = np.asarray(pars, dtype=float)
        def gm_shape(q):
            q = np.asarray(q, dtype=float)
            ge = sachs_family_value(q, pars, ge_family)
            _, f2 = yahl18_f1_f2(q)
            tau = q / (4.0 * MP2)
            return ge + (1.0 + tau) * f2
        #enddef
        return radius_from_shape(gm_shape, MU_P)
    #enddef

    rE, rEerr = propagate_scalar(radius_e, shape, cov)
    rM, rMerr = propagate_scalar(radius_m_derived, shape, cov)

    ndata = int(sum(item["N"] for item in prepared))
    nfree = int(len(fit_names))
    ndof = max(1, ndata - nfree)
    result = {
        "model": f"GE:{ge_family}|F2:YAHL18",
        "GE_family": ge_family,
        "F2_treatment": "YAHL18 fixed",
        "N": ndata,
        "n_parameters": nfree,
        "chi2": float(mn.fval),
        "ndof": int(ndof),
        "chi2_ndof": float(mn.fval / ndof),
        "valid": bool(
            mn.valid and np.isfinite(rE) and np.isfinite(rM)
        ),
        "rE_fm": float(rE),
        "rE_fit_err_fm": float(rEerr),
        "rM_fm": float(rM),
        "rM_fit_err_fm": float(rMerr),
        "shape_covariance_json": json.dumps(cov.tolist()),
    }
    for name in names_e:
        result[name] = float(mn.values[name])
        result[name + "_err"] = float(mn.errors[name])
    #endfor
    for nname in nuisance_names:
        result[nname] = float(mn.values[nname])
        result[nname + "_scale_factor"] = float(
            1.0 + nuisance_fracs[nname] * float(mn.values[nname])
        )
    #endfor
    return result
#enddef


def run_all_dataset_subset_and_external_f2_diagnostic(
        production_bundles: Sequence[Dict[str, object]],
        saylor_bundle: Dict[str, object],
        selection: pd.DataFrame,
        family: str,
        outdir: Path) -> pd.DataFrame:
    """
    Fit every non-empty subset of the six experiments at KM15 5%, using the nominal restricted Saylor sample.

    For every subset, run:
      (1) the closure-selected free-GE/free-GM family;
      (2) the same GE family with F2 fixed to Kelly 2004;
      (3) the same GE family with F2 fixed to YAHL18.

    This is diagnostic bookkeeping intended to expose exactly which dataset
    combinations determine the electric/magnetic solution.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    print(
        "[external-F2 preflight] starting restricted-Saylor subset diagnostic "
        f"with production family={family}"
    )

    selected_bundles = []
    for bundle in production_bundles:
        selected = select_bundle_from_external_model(
            bundle, selection, "km15", 0.05
        )
        selected_bundles.append((
            str(bundle["key"]),
            str(bundle["label"]),
            measurement_spec(
                str(bundle["key"]),
                str(bundle["label"]),
                str(bundle["kind"]),
                selected,
                float(bundle.get("norm_frac", 0.0)),
                unconstrained_norm=bool(
                    bundle.get("unconstrained_norm", False)
                ),
            ),
        ))
    #endfor

    saylor_all = saylor_bundle["all_data"].copy()
    delta_col = "bh_delta" if "bh_delta" in saylor_all.columns else "delta_bh"
    t_abs = np.abs(pd.to_numeric(
        saylor_all["t"], errors="coerce"
    ).to_numpy(float))
    saylor5 = saylor_all.loc[
        (
            pd.to_numeric(
                saylor_all[delta_col], errors="coerce"
            ).to_numpy(float) <= 0.05
        )
        & (t_abs >= 0.343)
    ].copy()

    if len(saylor5) == 0:
        raise RuntimeError(
            "Nominal restricted-Saylor fixed-F2 diagnostic selected zero points."
        )
    #endif
    if np.nanmin(np.abs(pd.to_numeric(
            saylor5["t"], errors="coerce"
        ).to_numpy(float))) < 0.343 - 1.0e-12:
        raise RuntimeError(
            "Restricted-Saylor preflight failed: a selected point has |t| < 0.343 GeV^2."
        )
    #endif
    print(
        "[external-F2 preflight] restricted Saylor: "
        f"N={len(saylor5)}, min|t|="
        f"{np.nanmin(np.abs(pd.to_numeric(saylor5['t'], errors='coerce').to_numpy(float))):.6f} GeV^2"
    )
    selected_bundles.append((
        str(saylor_bundle["key"]),
        str(saylor_bundle["label"]),
        measurement_spec(
            str(saylor_bundle["key"]),
            str(saylor_bundle["label"]),
            str(saylor_bundle["kind"]),
            saylor5,
            SAYLOR_GLOBAL_SCALE_FRAC,
        ),
    ))

    ge_family, _ = decode_sachs_family_pair(family)
    rows = []
    nset = len(selected_bundles)
    for mask in range(1, 1 << nset):
        picked = [
            selected_bundles[i]
            for i in range(nset) if mask & (1 << i)
        ]
        specs = [x[2] for x in picked]
        keys = [x[0] for x in picked]
        labels = [x[1] for x in picked]
        configuration = " + ".join(labels)

        for mode in [
            "free_GE_GM",
            "kelly_F2_fixed",
            "yahl18_F2_fixed",
        ]:
            try:
                if mode == "free_GE_GM":
                    fit = fit_sachs_family_multi_measurements(
                        specs,
                        family=family,
                        bh_cut=0.05,
                        add_moradi_bh_systematic=True,
                    )
                elif mode == "kelly_F2_fixed":
                    fit = fit_sachs_ge_family_fixed_kelly_f2_multi_measurements(
                        specs,
                        ge_family=ge_family,
                        bh_cut=0.05,
                        add_moradi_bh_systematic=True,
                    )
                else:
                    fit = fit_sachs_ge_family_fixed_yahl18_f2_multi_measurements(
                        specs,
                        ge_family=ge_family,
                        bh_cut=0.05,
                        add_moradi_bh_systematic=True,
                    )
                #endif
                rows.append({
                    "subset_mask": int(mask),
                    "dataset_count": int(len(specs)),
                    "dataset_keys": "|".join(keys),
                    "configuration_label": configuration,
                    "fit_mode": mode,
                    "production_family": family,
                    **fit,
                })
            except Exception as exc:
                rows.append({
                    "subset_mask": int(mask),
                    "dataset_count": int(len(specs)),
                    "dataset_keys": "|".join(keys),
                    "configuration_label": configuration,
                    "fit_mode": mode,
                    "production_family": family,
                    "valid": False,
                    "failure": str(exc),
                })
            #endtry

            # Persist progress after every fit.  This diagnostic is cheap
            # compared with the closure study, and the checkpoint guarantees
            # that a late failure cannot leave an empty directory.
            partial = pd.DataFrame(rows)
            partial.to_csv(
                outdir / "externalF2_partial_results.csv",
                index=False,
            )
        #endfor
    #endfor

    table = pd.DataFrame(rows)
    table.to_csv(
        outdir / "externalF2_partial_results.csv",
        index=False,
    )
    if "valid" in table.columns:
        failures = table.loc[
            ~table["valid"].fillna(False).astype(bool)
        ].copy()
        failures.to_csv(
            outdir / "externalF2_failures.csv",
            index=False,
        )
    else:
        pd.DataFrame().to_csv(
            outdir / "externalF2_failures.csv",
            index=False,
        )
    #endif

    # Fail fast if the nominal all-six external-F2 comparison did not actually
    # execute all three intended fit modes. This prevents a long production run
    # from silently falling back to the old Kelly-only bookkeeping.
    all_six_keys = "|".join(x[0] for x in selected_bundles)
    all_six_rows = table.loc[
        table["dataset_keys"].astype(str) == all_six_keys
    ].copy()
    expected_modes = {
        "free_GE_GM",
        "kelly_F2_fixed",
        "yahl18_F2_fixed",
    }
    present_modes = set(all_six_rows["fit_mode"].astype(str).tolist())
    if present_modes != expected_modes:
        raise RuntimeError(
            "Nominal all-six external-F2 diagnostic did not produce all three "
            f"fit modes. Present={sorted(present_modes)}, "
            f"expected={sorted(expected_modes)}"
        )
    #endif

    if "valid" in all_six_rows.columns:
        bad = all_six_rows.loc[
            ~all_six_rows["valid"].fillna(False).astype(bool)
        ]
        if len(bad):
            cols = [c for c in ["fit_mode", "failure"] if c in bad.columns]
            raise RuntimeError(
                "Nominal all-six external-F2 diagnostic contains failed fits:\n"
                + bad[cols].to_string(index=False)
            )
        #endif
    #endif

    raw_all_six_n = int(sum(
        len(spec["data"]) for _, _, spec in selected_bundles
    ))
    fitted_n = set()
    if "N" in all_six_rows.columns:
        fitted_n = set(
            pd.to_numeric(all_six_rows["N"], errors="coerce")
            .dropna().astype(int).tolist()
        )
        if len(fitted_n) != 1:
            raise RuntimeError(
                "Nominal all-six external-F2 modes did not use the same "
                f"number of fitted points: {sorted(fitted_n)}"
            )
        #endif
    #endif
    common_fit_n = (
        next(iter(fitted_n)) if fitted_n else raw_all_six_n
    )

    print(
        "[external-F2 preflight] nominal all-six comparison PASS: "
        f"Nfit={common_fit_n}, Nraw={raw_all_six_n}, "
        f"modes={sorted(present_modes)}"
    )

    table.to_csv(
        outdir / "all_six_dataset_subsets_restricted_saylor_free_vs_externalF2.csv",
        index=False,
    )

    # Canonical all-five/all-six comparison for immediate inspection.
    canonical = table.loc[
        table["dataset_keys"].astype(str).isin([
            "|".join(x[0] for x in selected_bundles[:5]),
            all_six_keys,
        ])
    ].copy()
    canonical.to_csv(
        outdir / "canonical_allfive_allsix_restricted_saylor_free_vs_externalF2.csv",
        index=False,
    )
    if len(canonical):
        print("\\n[all-subset free-vs-external-F2 diagnostic: canonical]")
        show = [
            "configuration_label", "fit_mode", "N", "chi2_ndof",
            "rE_fm", "rE_fit_err_fm", "rM_fm", "rM_fit_err_fm",
        ]
        print(canonical[[c for c in show if c in canonical.columns]].to_string(index=False))
    #endif
    return table
#enddef



def run_model_only_kelly_bh_selector_diagnostic(
        bundles: Sequence[Dict[str, object]],
        selection: pd.DataFrame,
        family: str,
        label: str,
        outdir: Path) -> pd.DataFrame:
    """
    Data-independent elastic-input BH-selection diagnostics.

    KM15 in Gepard inherits KellyEFF.  Therefore the historical
    Kelly-BH/KM15-EP matched-N selector is *not* an independent model test:
    apart from the EP-vs-BH denominator convention, it reuses the same elastic
    F1/F2 input as the nominal KM15 BH numerator.  It is retained as an
    algebraic/cache-consistency check.

    Genuine elastic-input variations are added with AMT 2007 and the
    A1/Bernauer polynomial-times-dipole reference.  For each reference,

      score = |sigma_EP^KM15 - sigma_BH^reference| / |sigma_BH^reference|,

    and the number of retained points is matched separately for each
    experiment to its nominal KM15-5% count.  No measured cross section enters
    selection, so these avoid the dependent-variable conditioning of the older
    observed-data residual/pull selectors.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    nominal_specs, nominal_counts = _km15_selected_specs_for_bundles(
        bundles, selection, 0.05
    )
    nominal_fit = fit_sachs_family_multi_measurements(
        nominal_specs,
        family=family,
        bh_cut=0.05,
        add_moradi_bh_systematic=True,
        bh_systematic_fraction=0.05,
    )

    reference_defs = [
        ("Kelly2004", "Kelly 2004", False),
        ("AMT2007", "AMT 2007", True),
        ("YAHL18", "YAHL18", True),
        (
            "BernauerA1",
            "A1/Bernauer order-8 poly×dipole",
            True,
        ),
    ]

    comparison_rows = [{
        "selector": "KM15_5pct",
        "elastic_reference": "KM15/Gepard KellyEFF",
        "selection_uses_observed_xs": False,
        "independent_elastic_input_from_KM15": False,
        "production_systematic_eligible": True,
        "note": "Nominal production selector.",
        **nominal_fit,
    }]
    overlap_rows = []
    point_rows = []

    for ref_tag, ref_label, independent in reference_defs:
        alt_specs = []
        for bundle in bundles:
            key = str(bundle["key"])
            data = bundle["all_data"].copy().reset_index(drop=True)
            data["source_row"] = np.arange(len(data), dtype=int)
            required = ["t_abs", "bh_A", "bh_B", "bh_C", "km15_ep"]
            missing = [c for c in required if c not in data.columns]
            if missing:
                raise KeyError(
                    f"{bundle['label']} missing model-only selector columns: "
                    f"{missing}"
                )
            #endif

            q = data["t_abs"].to_numpy(float)
            refs = elastic_reference_curves(q)
            if ref_label not in refs:
                raise KeyError(
                    f"Elastic reference {ref_label!r} is unavailable."
                )
            #endif
            ge_ref, gm_ref = refs[ref_label]
            f1_ref, f2_ref = sachs_to_f1f2(q, ge_ref, gm_ref)
            bh_ref = bh_from_f1f2(
                data["bh_A"].to_numpy(float),
                data["bh_B"].to_numpy(float),
                data["bh_C"].to_numpy(float),
                f1_ref, f2_ref,
            )
            km15ep = data["km15_ep"].to_numpy(float)
            score = (
                np.abs(km15ep - bh_ref)
                / np.maximum(np.abs(bh_ref), 1.0e-30)
            )
            score_col = f"model_only_{ref_tag}_vs_km15ep_score"
            data[score_col] = score
            data[f"{ref_tag}_bh_reference"] = bh_ref

            n_target = int(nominal_counts[key])
            finite = np.isfinite(score)
            ranked = data.loc[finite].sort_values(
                score_col, ascending=True
            )
            selected_alt = ranked.head(n_target).copy()
            if len(selected_alt) != n_target:
                raise RuntimeError(
                    f"{bundle['label']}: only {len(selected_alt)} finite "
                    f"{ref_label} model-only scores for matched N={n_target}"
                )
            #endif

            nominal_selected = select_bundle_from_external_model(
                bundle, selection, model="km15", threshold=0.05
            ).reset_index(drop=True)
            # Nominal point IDs use deterministic source-row positions.
            nominal_ids = set(
                selection.loc[
                    (selection["dataset"].astype(str) == key)
                    & (
                        selection["delta_bh_km15"].to_numpy(float) <= 0.05
                    ),
                    "source_row",
                ].astype(int).tolist()
            )
            alt_ids = set(selected_alt["source_row"].astype(int).tolist())
            inter = nominal_ids & alt_ids
            union = nominal_ids | alt_ids
            overlap_rows.append({
                "elastic_reference": ref_label,
                "independent_elastic_input_from_KM15": bool(independent),
                "dataset": key,
                "dataset_label": str(bundle["label"]),
                "N_nominal": int(len(nominal_ids)),
                "N_alternate": int(len(alt_ids)),
                "N_overlap": int(len(inter)),
                "nominal_recovered_fraction": (
                    float(len(inter) / len(nominal_ids))
                    if nominal_ids else np.nan
                ),
                "jaccard": (
                    float(len(inter) / len(union)) if union else np.nan
                ),
            })

            alt_specs.append({
                "key": key,
                "label": str(bundle["label"]),
                "kind": str(bundle["kind"]),
                "norm_frac": float(bundle.get("norm_frac", 0.0)),
                "data": selected_alt,
            })

            keep = selected_alt[[
                "source_row", "Q2", "xB", "t_abs", "phi_deg",
                "xs", score_col, f"{ref_tag}_bh_reference", "km15_ep",
            ]].copy()
            keep.insert(0, "dataset", key)
            keep.insert(1, "dataset_label", str(bundle["label"]))
            keep.insert(2, "elastic_reference", ref_label)
            keep["selected_by_nominal_km15_5pct"] = (
                keep["source_row"].astype(int).isin(nominal_ids)
            )
            point_rows.append(keep)
        #endfor

        alt_fit = fit_sachs_family_multi_measurements(
            alt_specs,
            family=family,
            bh_cut=0.05,
            add_moradi_bh_systematic=True,
            bh_systematic_fraction=0.05,
        )
        comparison_rows.append({
            "selector": f"{ref_tag}BH_vs_KM15EP_model_only_matchedN",
            "elastic_reference": ref_label,
            "selection_uses_observed_xs": False,
            "independent_elastic_input_from_KM15": bool(independent),
            "production_systematic_eligible": False,
            "note": (
                "Kelly is an algebraic/cache consistency check because "
                "Gepard KM15 inherits KellyEFF."
                if not independent else
                "Data-independent elastic-input variation; diagnostic pending "
                "systematic prescription."
            ),
            **alt_fit,
        })
    #endfor

    comparison = pd.DataFrame(comparison_rows)
    base_e = float(nominal_fit["rE_fm"])
    base_m = float(nominal_fit["rM_fm"])
    comparison["delta_rE_from_KM15_fm"] = (
        comparison["rE_fm"].to_numpy(float) - base_e
    )
    comparison["delta_rM_from_KM15_fm"] = (
        comparison["rM_fm"].to_numpy(float) - base_m
    )
    comparison.to_csv(
        outdir / "model_only_bh_radius_comparison.csv", index=False
    )
    pd.DataFrame(overlap_rows).to_csv(
        outdir / "model_only_bh_selection_overlap.csv", index=False
    )
    if point_rows:
        pd.concat(point_rows, ignore_index=True).to_csv(
            outdir / "model_only_bh_selection_points.csv", index=False
        )
    #endif

    fig, axes = plt.subplots(1, 2, figsize=(12.2, 5.1))
    xx = np.arange(len(comparison))
    axes[0].errorbar(
        xx, comparison["rE_fm"], yerr=comparison["rE_fit_err_fm"],
        marker="o", linestyle="none", capsize=3,
    )
    axes[1].errorbar(
        xx, comparison["rM_fm"], yerr=comparison["rM_fit_err_fm"],
        marker="o", linestyle="none", capsize=3,
    )
    labels = [
        "KM15 5%",
        "Kelly\n(same KM15 EFF)",
        "AMT 2007",
        "YAHL18",
        "A1/Bernauer",
    ]
    for ax, ylabel in zip(axes, [r"$r_E$ [fm]", r"$r_M$ [fm]"]):
        ax.set_xticks(xx)
        ax.set_xticklabels(labels)
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", alpha=0.2)
    #endfor
    fig.suptitle(
        f"{label}: data-independent elastic-input BH-selection diagnostics\n"
        "matched N per experiment; Kelly is not independent of KM15",
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(outdir / "model_only_bh_radius_comparison.png", dpi=260)
    plt.close(fig)
    return comparison
#enddef



def _shape_parameters_and_covariance_from_result(
        fit: Dict[str, object],
        family: str) -> Tuple[np.ndarray, np.ndarray, int, int, str, str]:
    """Reconstruct Sachs shape parameters/covariance from a saved fit dict."""
    family_e, family_m = decode_sachs_family_pair(family)
    ne = int(re.findall(r"\d+", family_e)[0])
    nm = int(re.findall(r"\d+", family_m)[0])
    names = (
        [f"e{i}" for i in range(1, ne + 1)]
        + [f"m{i}" for i in range(1, nm + 1)]
    )
    params = np.asarray([float(fit[n]) for n in names], dtype=float)
    raw_cov = fit.get("shape_covariance_json", "")
    try:
        cov = np.asarray(json.loads(str(raw_cov)), dtype=float)
    except Exception:
        cov = np.full((ne + nm, ne + nm), np.nan, dtype=float)
    #endtry
    if cov.shape != (ne + nm, ne + nm):
        cov = np.full((ne + nm, ne + nm), np.nan, dtype=float)
    #endif
    return params, cov, ne, nm, family_e, family_m
#enddef


def _sachs_band_from_result(
        fit: Dict[str, object],
        family: str,
        q: np.ndarray,
        which: str) -> Tuple[np.ndarray, np.ndarray]:
    """Evaluate GE or GM and its full Hessian covariance band."""
    q = np.asarray(q, dtype=float)
    params, cov, ne, nm, family_e, family_m = (
        _shape_parameters_and_covariance_from_result(fit, family)
    )

    def evaluate(p):
        ce = np.asarray(p[:ne], dtype=float)
        cm = np.asarray(p[ne:ne + nm], dtype=float)
        ge = sachs_family_value(q, ce, family_e)
        gm = MU_P * sachs_family_value(q, cm, family_m)
        return ge if which == "GE" else gm
    #enddef

    central = np.asarray(evaluate(params), dtype=float)
    sigma = np.full_like(central, np.nan)
    if not np.all(np.isfinite(cov)):
        return central, sigma
    #endif

    for iq in range(len(q)):
        def scalar(p):
            return float(evaluate(np.asarray(p, dtype=float))[iq])
        #enddef
        grad = numerical_gradient(scalar, params)
        var = float(grad @ cov @ grad)
        if np.isfinite(var):
            sigma[iq] = math.sqrt(max(var, 0.0))
        #endif
    #endfor
    return central, sigma
#enddef


def _band_overlap_metrics(
        ca: np.ndarray,
        sa: np.ndarray,
        cb: np.ndarray,
        sb: np.ndarray) -> Dict[str, float]:
    """Pointwise 68% band-overlap and standardized-separation diagnostics."""
    ca = np.asarray(ca, dtype=float)
    sa = np.asarray(sa, dtype=float)
    cb = np.asarray(cb, dtype=float)
    sb = np.asarray(sb, dtype=float)
    finite = (
        np.isfinite(ca) & np.isfinite(sa)
        & np.isfinite(cb) & np.isfinite(sb)
        & (sa >= 0.0) & (sb >= 0.0)
    )
    if not np.any(finite):
        return {
            "band_overlap_fraction": np.nan,
            "central_a_inside_band_b_fraction": np.nan,
            "central_b_inside_band_a_fraction": np.nan,
            "median_standardized_curve_separation": np.nan,
            "max_standardized_curve_separation": np.nan,
            "median_absolute_fractional_curve_difference": np.nan,
            "max_absolute_fractional_curve_difference": np.nan,
        }
    #endif
    ca = ca[finite]; sa = sa[finite]; cb = cb[finite]; sb = sb[finite]
    overlap = np.maximum(ca - sa, cb - sb) <= np.minimum(ca + sa, cb + sb)
    a_in_b = (ca >= cb - sb) & (ca <= cb + sb)
    b_in_a = (cb >= ca - sa) & (cb <= ca + sa)
    denom = np.sqrt(sa**2 + sb**2)
    d = np.abs(ca - cb) / np.maximum(denom, 1.0e-30)
    frac = np.abs(ca - cb) / np.maximum(0.5 * (np.abs(ca) + np.abs(cb)), 1.0e-30)
    return {
        "band_overlap_fraction": float(np.mean(overlap)),
        "central_a_inside_band_b_fraction": float(np.mean(a_in_b)),
        "central_b_inside_band_a_fraction": float(np.mean(b_in_a)),
        "median_standardized_curve_separation": float(np.median(d)),
        "max_standardized_curve_separation": float(np.max(d)),
        "median_absolute_fractional_curve_difference": float(np.median(frac)),
        "max_absolute_fractional_curve_difference": float(np.max(frac)),
    }
#enddef


def save_threshold_form_factor_band_study(
        bundles: Sequence[Dict[str, object]],
        selection: pd.DataFrame,
        family: str,
        threshold_table: pd.DataFrame,
        label: str,
        outdir: Path) -> pd.DataFrame:
    """
    Compare finite-|t| GE/GM shapes across the 1--7% KM15 selections.

    Radius changes can strongly amplify small differences in the extrapolated
    Q2->0 slope.  This diagnostic therefore compares the fitted form-factor
    curves themselves, with full Hessian covariance bands, over both the
    extrapolation region and the common data-supported |t| interval.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    wanted = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07]
    rows_by_thr = {}
    supports = {}
    for thr in wanted:
        part = threshold_table.loc[
            np.isclose(threshold_table["threshold"].to_numpy(float), thr)
            & threshold_table["valid"].astype(bool)
        ]
        if len(part) != 1:
            continue
        #endif
        rows_by_thr[thr] = dict(part.iloc[0])
        specs, _ = _km15_selected_specs_for_bundles(bundles, selection, thr)
        vals = []
        for spec in specs:
            if len(spec["data"]):
                vals.extend(spec["data"]["t_abs"].to_numpy(float).tolist())
            #endif
        #endfor
        vals = np.asarray(vals, dtype=float)
        vals = vals[np.isfinite(vals)]
        if len(vals):
            supports[thr] = (float(np.min(vals)), float(np.max(vals)))
        #endif
    #endfor

    if 0.05 not in rows_by_thr or not supports:
        return pd.DataFrame()
    #endif

    qmax = max(v[1] for v in supports.values())
    q = np.linspace(0.0, max(0.05, 1.03 * qmax), 500)
    common_min = max(v[0] for v in supports.values())
    common_max = min(v[1] for v in supports.values())
    if common_max <= common_min:
        common_min = supports[0.05][0]
        common_max = supports[0.05][1]
    #endif
    common_mask = (q >= common_min) & (q <= common_max)
    extrap_mask = (q >= 0.0) & (q < common_min)

    # High-information figure: central curves plus 68% Hessian bands.
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.3))
    for thr in sorted(rows_by_thr):
        fit = rows_by_thr[thr]
        for ax, which in zip(axes, ["GE", "GM"]):
            c, s = _sachs_band_from_result(fit, family, q, which)
            if which == "GM":
                c = c / MU_P
                s = s / MU_P
            #endif
            lw = 2.2 if abs(thr - 0.05) < 1.0e-12 else 1.15
            alpha = 0.16 if abs(thr - 0.05) < 1.0e-12 else 0.055
            line, = ax.plot(
                q, c, linewidth=lw,
                label=f"{100.0*thr:.0f}% selection",
            )
            ax.fill_between(
                q, c - s, c + s,
                alpha=alpha, color=line.get_color(), linewidth=0.0,
            )
        #endfor
    #endfor
    for ax in axes:
        ax.axvspan(
            common_min, common_max, alpha=0.08,
            label="common selected-data |t| support" if ax is axes[0] else None,
        )
        ax.axvline(common_min, linewidth=0.8, linestyle=":")
        ax.set_xlabel(r"$|t|$ (GeV$^2$)")
        ax.grid(alpha=0.20)
    #endfor
    axes[0].set_ylabel(r"$G_E$")
    axes[1].set_ylabel(r"$G_M/\mu_p$")
    handles, labs = axes[0].get_legend_handles_labels()
    fig.legend(
        handles, labs, loc="upper center", ncol=4,
        bbox_to_anchor=(0.5, 0.965), frameon=True,
    )
    fig.suptitle(
        f"{label}: finite-$|t|$ stability across KM15 BH-purity selections\n"
        f"closure-selected {family}; shaded curves are 68% Hessian bands",
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    fig.savefig(outdir / "01_GE_GM_hessian_bands_thresholds_1to7.png", dpi=260)
    plt.close(fig)

    base = rows_by_thr[0.05]
    metric_rows = []
    curve_rows = []
    for thr in sorted(rows_by_thr):
        fit = rows_by_thr[thr]
        for which in ["GE", "GM"]:
            c, s = _sachs_band_from_result(fit, family, q, which)
            c5, s5 = _sachs_band_from_result(base, family, q, which)
            scale = MU_P if which == "GM" else 1.0
            c = c / scale; s = s / scale
            c5 = c5 / scale; s5 = s5 / scale

            for region, mask in [
                ("common_data_support", common_mask),
                ("extrapolation_to_zero", extrap_mask),
                ("zero_to_common_max", q <= common_max),
            ]:
                if not np.any(mask):
                    continue
                #endif
                metrics = _band_overlap_metrics(
                    c[mask], s[mask], c5[mask], s5[mask]
                )
                metric_rows.append({
                    "configuration_label": label,
                    "family": family,
                    "form_factor": which,
                    "threshold": thr,
                    "threshold_percent": 100.0 * thr,
                    "reference_threshold": 0.05,
                    "region": region,
                    "region_qmin": float(np.min(q[mask])),
                    "region_qmax": float(np.max(q[mask])),
                    **metrics,
                })
            #endfor

            for iq in range(len(q)):
                curve_rows.append({
                    "form_factor": which,
                    "threshold": thr,
                    "threshold_percent": 100.0 * thr,
                    "q": float(q[iq]),
                    "central": float(c[iq]),
                    "hessian_sigma": float(s[iq]),
                    "common_data_support": bool(common_mask[iq]),
                    "extrapolation_to_zero": bool(extrap_mask[iq]),
                })
            #endfor
        #endfor
    #endfor

    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(outdir / "GE_GM_threshold_band_overlap_metrics.csv", index=False)
    pd.DataFrame(curve_rows).to_csv(
        outdir / "GE_GM_threshold_hessian_band_curves.csv", index=False
    )

    # Compact summary of how radius excursions compare with finite-t curve motion.
    finite_part = metrics.loc[
        (metrics["region"] == "common_data_support")
        & (metrics["threshold"] != 0.05)
    ].copy()
    if len(finite_part):
        fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
        for ax, which in zip(axes, ["GE", "GM"]):
            p = finite_part.loc[finite_part["form_factor"] == which]
            ax.plot(
                p["threshold_percent"],
                100.0 * p["band_overlap_fraction"],
                marker="o", label="68% band overlap",
            )
            ax.plot(
                p["threshold_percent"],
                100.0 * p["central_a_inside_band_b_fraction"],
                marker="s", label="alternate central inside 5% band",
            )
            ax.set_ylim(0.0, 105.0)
            ax.set_xlabel("KM15 BH-purity threshold (%)")
            ax.set_ylabel("common-support overlap (%)")
            ax.set_title(which)
            ax.grid(alpha=0.2)
        #endfor
        handles, labs = axes[0].get_legend_handles_labels()
        fig.legend(handles, labs, loc="upper center", ncol=2,
                   bbox_to_anchor=(0.5, 0.96))
        fig.suptitle(
            f"{label}: are threshold-dependent radius shifts visible at finite $|t|$?",
            y=0.995,
        )
        fig.tight_layout(rect=(0, 0, 1, 0.88))
        fig.savefig(outdir / "02_GE_GM_common_support_band_overlap.png", dpi=260)
        plt.close(fig)
    #endif
    return metrics
#enddef


def save_competitive_family_form_factor_bands(
        valid: pd.DataFrame,
        chosen_family: str,
        label: str,
        selected_specs: Sequence[Dict[str, object]],
        outdir: Path) -> pd.DataFrame:
    """Compare closure-competitive real-data GE/GM curves and their ratios."""
    outdir.mkdir(parents=True, exist_ok=True)
    qvals = []
    for spec in selected_specs:
        if len(spec["data"]):
            qvals.extend(spec["data"]["t_abs"].to_numpy(float).tolist())
        #endif
    #endfor
    qvals = np.asarray(qvals, dtype=float)
    qvals = qvals[np.isfinite(qvals)]
    if len(qvals) == 0:
        return pd.DataFrame()
    #endif
    qmin = float(np.min(qvals)); qmax = float(np.max(qvals))
    # Extend the displayed form-factor curves beyond the selected-data support
    # to make the extrapolation behavior of closure-competitive families visible.
    # The current proton samples end below |t|=1.2 GeV^2, so use a common
    # 0--1.2 GeV^2 presentation range for both GE and GM panels.
    q_plot_max = 1.20
    q = np.linspace(0.0, q_plot_max, 600)
    support = (q >= qmin) & (q <= qmax)
    extrap = q < qmin

    chosen_rows = valid.loc[valid["family"].astype(str) == chosen_family]
    if len(chosen_rows) != 1:
        return pd.DataFrame()
    #endif
    chosen = dict(chosen_rows.iloc[0])
    chosen_curves = {}
    for which in ["GE", "GM"]:
        c0, s0 = _sachs_band_from_result(chosen, chosen_family, q, which)
        scale = MU_P if which == "GM" else 1.0
        chosen_curves[which] = (c0 / scale, s0 / scale)
    #endfor

    # Top row: absolute normalized Sachs functions with 68% Hessian bands.
    # Bottom row: central-function ratio to the closure-selected family.  To
    # avoid overloading the ratio panels, show propagated Hessian ratio bands
    # only for the two non-chosen families whose closure objective is closest
    # to that of the chosen family.  The inter-family covariance is not
    # available from separate Hessian fits, so those ratio bands use the
    # conservative independence approximation.
    chosen_objective = float(chosen.get("closure_objective_fm", np.nan))
    ratio_band_families = []
    if np.isfinite(chosen_objective) and "closure_objective_fm" in valid.columns:
        nearest = valid.loc[
            valid["family"].astype(str) != chosen_family
        ].copy()
        nearest["_objective_distance"] = np.abs(
            nearest["closure_objective_fm"].to_numpy(float) - chosen_objective
        )
        nearest = nearest.loc[
            np.isfinite(nearest["_objective_distance"].to_numpy(float))
        ].sort_values("_objective_distance", ascending=True)
        ratio_band_families = nearest["family"].astype(str).head(2).tolist()
    #endif

    fig, axes = plt.subplots(
        2, 2, figsize=(13.4, 8.3), sharex="col",
        gridspec_kw={"height_ratios": [1.55, 0.85]},
    )
    metric_rows = []
    ratio_band_rows = []
    ratio_extrema = {"GE": [], "GM": []}
    handles = []
    labs = []
    for _, rr in valid.iterrows():
        family = str(rr["family"])
        fit = dict(rr)
        for j, which in enumerate(["GE", "GM"]):
            c, serr = _sachs_band_from_result(fit, family, q, which)
            c0, s0 = chosen_curves[which]
            scale = MU_P if which == "GM" else 1.0
            c = c / scale; serr = serr / scale
            lw = 2.4 if family == chosen_family else 1.15
            alpha = 0.16 if family == chosen_family else 0.045
            line, = axes[0, j].plot(
                q, c, linewidth=lw,
                label=(f"{family} (chosen)" if family == chosen_family else family),
            )
            axes[0, j].fill_between(
                q, c - serr, c + serr, alpha=alpha,
                color=line.get_color(), linewidth=0.0,
            )
            denom = np.maximum(np.abs(c0), 1.0e-30)
            ratio = c / denom
            axes[1, j].plot(q, ratio, linewidth=lw, color=line.get_color())

            # For the two closure-nearest alternative families, propagate both
            # Hessian bands into the ratio.  With no cross-family covariance
            # matrix available, use
            #   Var(c/c0) = (serr/c0)^2 + (c*s0/c0^2)^2.
            # Do not draw a band for chosen/chosen: that ratio is identically 1
            # for the same realization and therefore has exactly zero width.
            if family in ratio_band_families:
                ratio_sigma = np.sqrt(
                    (serr / denom) ** 2
                    + ((c * s0) / (denom ** 2)) ** 2
                )
                good_ratio_band = (
                    np.isfinite(ratio)
                    & np.isfinite(ratio_sigma)
                    & (ratio_sigma >= 0.0)
                )
                axes[1, j].fill_between(
                    q,
                    np.where(good_ratio_band, ratio - ratio_sigma, np.nan),
                    np.where(good_ratio_band, ratio + ratio_sigma, np.nan),
                    alpha=0.12,
                    color=line.get_color(),
                    linewidth=0.0,
                )
                for iq in np.flatnonzero(good_ratio_band):
                    ratio_band_rows.append({
                        "configuration_label": label,
                        "family": family,
                        "chosen_family": chosen_family,
                        "form_factor": which,
                        "t_abs_GeV2": float(q[iq]),
                        "ratio": float(ratio[iq]),
                        "ratio_hessian_sigma_independent": float(ratio_sigma[iq]),
                        "closure_objective_fm": float(
                            fit.get("closure_objective_fm", np.nan)
                        ),
                        "chosen_closure_objective_fm": float(chosen_objective),
                        "inter_family_covariance_included": False,
                    })
                #endfor
            #endif
            ratio_extrema[which].extend(ratio[np.isfinite(ratio)].tolist())
            if j == 0:
                handles.append(line); labs.append(line.get_label())
            #endif
            for region, mask in [
                ("selected_data_support", support),
                ("extrapolation_to_zero", extrap),
            ]:
                if np.any(mask):
                    metric_rows.append({
                        "configuration_label": label,
                        "family": family,
                        "chosen_family": chosen_family,
                        "form_factor": which,
                        "region": region,
                        **_band_overlap_metrics(
                            c[mask], serr[mask], c0[mask], s0[mask]
                        ),
                    })
                #endif
            #endfor
        #endfor
    #endfor

    axes[0, 0].set_ylabel(r"$G_E$")
    axes[0, 1].set_ylabel(r"$G_M/\mu_p$")
    axes[1, 0].set_ylabel(r"$G_E/G_E^{\rm chosen}$")
    axes[1, 1].set_ylabel(r"$G_M/G_M^{\rm chosen}$")
    for j, which in enumerate(["GE", "GM"]):
        for ax in axes[:, j]:
            ax.axvspan(qmin, qmax, alpha=0.07)
            ax.axvline(qmin, linewidth=0.8, linestyle=":")
            ax.grid(alpha=0.2)
        #endfor
        axes[1, j].axhline(1.0, linewidth=0.8, linestyle="--")
        axes[1, j].set_xlabel(r"$|t|$ (GeV$^2$)")
        # Use identical ratio scales so GE and GM family differences can be
        # compared directly by eye.
        axes[1, j].set_ylim(0.70, 1.30)
        axes[0, j].set_xlim(0.0, q_plot_max)
        axes[1, j].set_xlim(0.0, q_plot_max)
    #endfor

    # Reserve a dedicated header region.  The legend spans two rows and is kept
    # well separated from both the title and the explanatory line.
    fig.suptitle(
        f"{label}: closure-competitive Sachs functions",
        y=0.992, fontsize=15,
    )
    fig.legend(
        handles, labs, loc="upper center", ncol=4,
        bbox_to_anchor=(0.5, 0.952), frameon=True, fontsize=9,
        borderaxespad=0.0,
    )
    ratio_band_text = (
        ", ".join(ratio_band_families) if len(ratio_band_families)
        else "none"
    )
    fig.text(
        0.5, 0.865,
        (
            "Top bands: 68% Hessian uncertainty; bottom ratio bands: "
            f"{ratio_band_text} only (inter-family covariance neglected); "
            "shaded region: selected-data |t| support"
        ),
        ha="center", va="center", fontsize=9.2,
    )
    fig.subplots_adjust(
        top=0.80, bottom=0.09, left=0.08, right=0.985,
        hspace=0.08, wspace=0.16,
    )
    fig.savefig(outdir / "02_competitive_family_GE_GM_hessian_bands.png", dpi=300)
    plt.close(fig)

    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(
        outdir / "competitive_family_GE_GM_band_overlap_metrics.csv",
        index=False,
    )
    pd.DataFrame(ratio_band_rows).to_csv(
        outdir / "competitive_family_GE_GM_ratio_hessian_bands.csv",
        index=False,
    )
    if len(ratio_band_families):
        print(
            f"[competitive-family ratio bands] {label}: "
            f"{', '.join(ratio_band_families)} "
            "(two closest closure objectives to chosen; "
            "inter-family covariance neglected)"
        )
    #endif
    return metrics
#enddef

def save_f1_f2_bh_sensitivity_diagnostics(
        bundles: Sequence[Dict[str, object]],
        selection: pd.DataFrame,
        fit: Dict[str, object],
        family: str,
        outdir: Path) -> pd.DataFrame:
    """
    Exact quadratic BH decomposition and logarithmic F1/F2 sensitivities.

    sigma_BH = A F1^2 + B F1 F2 + C F2^2
    S_F1 = d ln sigma / d ln F1
    S_F2 = d ln sigma / d ln F2

    The exact term fractions avoid the ambiguity of calling F1=0/F2=0 curves
    "additive contributions", while those turn-off ratios are also tabulated
    for intuitive comparison with the original Moradi-style visualization.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    params, _, ne, nm, family_e, family_m = (
        _shape_parameters_and_covariance_from_result(fit, family)
    )
    ce = params[:ne]; cm = params[ne:ne + nm]

    point_rows = []
    for bundle in bundles:
        specs, _ = _km15_selected_specs_for_bundles(
            [bundle], selection, 0.05
        )
        if not specs or len(specs[0]["data"]) == 0:
            continue
        #endif
        d = specs[0]["data"].copy()
        q = d["t_abs"].to_numpy(float)
        ge = sachs_family_value(q, ce, family_e)
        gm = MU_P * sachs_family_value(q, cm, family_m)
        tau = q / (4.0 * MP2)
        f1 = (ge + tau * gm) / (1.0 + tau)
        f2 = (gm - ge) / (1.0 + tau)
        A = d["bh_A"].to_numpy(float)
        B = d["bh_B"].to_numpy(float)
        C = d["bh_C"].to_numpy(float)
        t11 = A * f1**2
        t12 = B * f1 * f2
        t22 = C * f2**2
        full = t11 + t12 + t22
        sf1 = (2.0 * t11 + t12) / np.maximum(np.abs(full), 1.0e-30)
        sf2 = (t12 + 2.0 * t22) / np.maximum(np.abs(full), 1.0e-30)
        f1_only = t11
        f2_only = t22
        for irow in range(len(d)):
            point_rows.append({
                "dataset": str(bundle["key"]),
                "dataset_label": str(bundle["label"]),
                "t_abs": float(q[irow]),
                "F1": float(f1[irow]),
                "F2": float(f2[irow]),
                "sigma_full_bh": float(full[irow]),
                "fraction_F1sq": float(t11[irow] / full[irow]),
                "fraction_F1F2": float(t12[irow] / full[irow]),
                "fraction_F2sq": float(t22[irow] / full[irow]),
                "S_F1": float(sf1[irow]),
                "S_F2": float(sf2[irow]),
                "abs_SF2_over_abs_SF1": float(
                    abs(sf2[irow]) / max(abs(sf1[irow]), 1.0e-30)
                ),
                "sigma_F2_zero_over_full": float(f1_only[irow] / full[irow]),
                "sigma_F1_zero_over_full": float(f2_only[irow] / full[irow]),
            })
        #endfor
    #endfor

    pts = pd.DataFrame(point_rows)
    pts.to_csv(outdir / "f1_f2_bh_sensitivity_points.csv", index=False)
    if len(pts) == 0:
        return pts
    #endif

    summary_rows = []
    for key, d in list(pts.groupby("dataset", sort=False)) + [("ALL", pts)]:
        label = "All four combined" if key == "ALL" else str(d.iloc[0]["dataset_label"])
        summary_rows.append({
            "dataset": key,
            "dataset_label": label,
            "N": int(len(d)),
            "median_S_F1": float(np.median(d["S_F1"])),
            "median_S_F2": float(np.median(d["S_F2"])),
            "median_abs_S_F1": float(np.median(np.abs(d["S_F1"]))),
            "median_abs_S_F2": float(np.median(np.abs(d["S_F2"]))),
            "median_abs_SF2_over_abs_SF1": float(
                np.median(d["abs_SF2_over_abs_SF1"])
            ),
            "median_fraction_F1sq": float(np.median(d["fraction_F1sq"])),
            "median_fraction_F1F2": float(np.median(d["fraction_F1F2"])),
            "median_fraction_F2sq": float(np.median(d["fraction_F2sq"])),
        })
    #endfor
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(outdir / "f1_f2_bh_sensitivity_summary.csv", index=False)

    # Five panels: one per experiment plus the combined ensemble.
    panel_defs = [
        (str(b["key"]), str(b["label"])) for b in bundles
    ] + [("ALL", "All four combined")]
    fig, axes = plt.subplots(2, 3, figsize=(14.5, 8.8), sharex=False)
    for ax, (key, label) in zip(axes.flat, panel_defs):
        d = pts if key == "ALL" else pts.loc[pts["dataset"] == key]
        ax.scatter(
            d["t_abs"], d["fraction_F1sq"],
            s=13, alpha=0.50, label=r"$A F_1^2/\sigma_{\rm BH}$",
        )
        ax.scatter(
            d["t_abs"], d["fraction_F1F2"],
            s=13, alpha=0.50, label=r"$B F_1F_2/\sigma_{\rm BH}$",
        )
        ax.scatter(
            d["t_abs"], d["fraction_F2sq"],
            s=13, alpha=0.50, label=r"$C F_2^2/\sigma_{\rm BH}$",
        )
        ax.axhline(0.0, linewidth=0.7)
        ax.set_title(label)
        ax.set_xlabel(r"$|t|$ (GeV$^2$)")
        ax.set_ylabel("exact BH term / total BH")
        ax.grid(alpha=0.18)
    #endfor
    if len(panel_defs) < len(axes.flat):
        axes.flat[-1].axis("off")
    #endif
    handles, labs = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labs, loc="upper center", ncol=3,
               bbox_to_anchor=(0.5, 0.965))
    fig.suptitle(
        f"BH electromagnetic decomposition at KM15 5% selection ({family})",
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    fig.savefig(outdir / "01_F1_F2_exact_BH_term_fractions.png", dpi=260)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.2, 6.8))
    for key, label in panel_defs[:-1]:
        d = pts.loc[pts["dataset"] == key]
        ax.scatter(
            d["S_F1"], d["S_F2"], s=18, alpha=0.55, label=label
        )
    #endfor
    ax.axhline(0.0, linewidth=0.7)
    ax.axvline(0.0, linewidth=0.7)
    ax.set_xlabel(r"$S_{F_1}=\partial\ln\sigma_{\rm BH}/\partial\ln F_1$")
    ax.set_ylabel(r"$S_{F_2}=\partial\ln\sigma_{\rm BH}/\partial\ln F_2$")
    ax.set_title(
        "Complementary electromagnetic sensitivity of the four proton datasets"
    )
    ax.grid(alpha=0.2)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "02_F1_F2_log_sensitivity_plane.png", dpi=260)
    plt.close(fig)

    # Compact note/slide-ready medians.
    fig, ax = plt.subplots(figsize=(9.2, 5.4))
    ss = summary.loc[summary["dataset"] != "ALL"].copy()
    x = np.arange(len(ss))
    width = 0.36
    ax.bar(x - width/2, ss["median_abs_S_F1"], width, label=r"median $|S_{F_1}|$")
    ax.bar(x + width/2, ss["median_abs_S_F2"], width, label=r"median $|S_{F_2}|$")
    ax.set_xticks(x)
    ax.set_xticklabels(ss["dataset_label"], rotation=18, ha="right")
    ax.set_ylabel("median absolute logarithmic sensitivity")
    ax.set_title("Why adding Hall A improves simultaneous $F_1/F_2$ leverage")
    ax.grid(axis="y", alpha=0.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "03_F1_F2_sensitivity_medians.png", dpi=300)
    plt.close(fig)

    # Intuitive Moradi-style turn-off demonstration.
    #
    # For every selected point we have the exact quadratic BH decomposition
    #
    #   sigma_BH = A F1^2 + B F1 F2 + C F2^2.
    #
    # Setting F1=0 leaves C F2^2; setting F2=0 leaves A F1^2.  These are not
    # additive "fractions" because the full cross section also contains the
    # interference term B F1 F2.  The upper panel shows the actual cross-section
    # scale of the three cases across the measured |t| support.  The lower panel
    # shows the ratio of the F2-only to F1-only cases, which is a compact visual
    # measure of relative electromagnetic leverage.
    #
    # Because xB, Q2 and phi are not fixed as |t| changes across a dataset,
    # individual points are shown explicitly.  A running median in |t| bins is
    # overlaid only as a visual guide; it is not a model curve.
    def _turnoff_plot_for_points(d: pd.DataFrame, label: str, stem: str) -> None:
        if len(d) == 0:
            return
        #endif

        work = d.sort_values("t_abs").copy()
        t = work["t_abs"].to_numpy(float)
        full = work["sigma_full_bh"].to_numpy(float)
        f1zero = work["sigma_F1_zero_over_full"].to_numpy(float) * full
        f2zero = work["sigma_F2_zero_over_full"].to_numpy(float) * full

        with np.errstate(divide="ignore", invalid="ignore"):
            leverage_ratio = np.divide(
                f1zero,
                f2zero,
                out=np.full_like(f1zero, np.nan, dtype=float),
                where=np.abs(f2zero) > 1.0e-30,
            )
        #endwith

        fig, axes = plt.subplots(
            2, 1,
            figsize=(9.2, 7.8),
            sharex=True,
            gridspec_kw={"height_ratios": [2.1, 1.0]},
        )
        ax, rax = axes

        # Measured-point values.
        ax.scatter(t, full, s=15, alpha=0.30, label="Full BH")
        ax.scatter(t, f1zero, s=15, alpha=0.30, marker="s",
                   label=r"$F_1=0$ ($F_2$ only)")
        ax.scatter(t, f2zero, s=15, alpha=0.30, marker="^",
                   label=r"$F_2=0$ ($F_1$ only)")
        rax.scatter(t, leverage_ratio, s=15, alpha=0.35)

        # Running medians in equal-population |t| bins.  This keeps the display
        # readable without pretending the data constitute a fixed-kinematics
        # continuous trajectory in t.
        nbin = min(14, max(4, int(round(math.sqrt(len(work))))))
        try:
            groups = pd.qcut(
                work["t_abs"],
                q=nbin,
                duplicates="drop",
            )
            tmp = pd.DataFrame({
                "t": t,
                "full": full,
                "f1zero": f1zero,
                "f2zero": f2zero,
                "ratio": leverage_ratio,
                "_grp": groups,
            })
            med = tmp.groupby("_grp", observed=True).median(numeric_only=True)
            med = med.sort_values("t")
            ax.plot(med["t"], med["full"], linewidth=2.0)
            ax.plot(med["t"], med["f1zero"], linewidth=2.0, linestyle="--")
            ax.plot(med["t"], med["f2zero"], linewidth=2.0, linestyle=":")
            rax.plot(med["t"], med["ratio"], linewidth=2.0)
        except Exception:
            pass
        #endtry

        ax.set_ylabel(r"BH cross section")
        ax.set_title(
            f"{label}: electromagnetic leverage across the selected $|t|$ range"
        )
        ax.grid(alpha=0.18)
        ax.legend(fontsize=9, loc="best")

        rax.axhline(1.0, linewidth=0.8, linestyle="--")
        rax.set_xlabel(r"$|t|$ (GeV$^2$)")
        rax.set_ylabel(
            r"$\sigma_{\rm BH}(F_1=0)\,/\,\sigma_{\rm BH}(F_2=0)$"
        )
        rax.grid(alpha=0.18)

        fig.subplots_adjust(
            top=0.93, bottom=0.10, left=0.10, right=0.985,
            hspace=0.08,
        )
        fig.savefig(
            outdir / f"04_F1_F2_turnoff_{stem}.png",
            dpi=300,
        )
        plt.close(fig)
    #enddef

    for key, label in panel_defs[:-1]:
        _turnoff_plot_for_points(
            pts.loc[pts["dataset"] == key].copy(),
            label,
            str(key),
        )
    #endfor
    _turnoff_plot_for_points(
        pts.copy(),
        "All four combined",
        "all_four_combined",
    )

    # Keep one compact all-dataset summary of the *relative* turn-off leverage.
    # This is useful for a slide where five separate absolute-cross-section
    # panels would be too dense.
    fig, ax = plt.subplots(figsize=(9.4, 6.0))
    for key, label in panel_defs[:-1]:
        d = pts.loc[pts["dataset"] == key].copy()
        denom = d["sigma_F2_zero_over_full"].to_numpy(float)
        numer = d["sigma_F1_zero_over_full"].to_numpy(float)
        ratio = np.divide(
            numer, denom,
            out=np.full(len(d), np.nan, dtype=float),
            where=np.abs(denom) > 1.0e-30,
        )
        ax.scatter(
            d["t_abs"], ratio,
            s=16, alpha=0.45, label=label,
        )
    #endfor
    ax.axhline(1.0, linewidth=0.8, linestyle="--")
    ax.set_xlabel(r"$|t|$ (GeV$^2$)")
    ax.set_ylabel(
        r"$\sigma_{\rm BH}(F_1=0)\,/\,\sigma_{\rm BH}(F_2=0)$"
    )
    ax.set_title(
        r"Relative $F_2$-only versus $F_1$-only BH leverage at KM15 5%"
    )
    ax.grid(alpha=0.18)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(
        outdir / "05_F1_F2_turnoff_leverage_ratio_all_datasets.png",
        dpi=300,
    )
    plt.close(fig)

    return summary
#enddef



def save_preferred_sachs_vs_elastic_data(
        fit: Dict[str, object],
        family: str,
        selected_specs: Sequence[Dict[str, object]],
        outdir: Path) -> None:
    """
    Compare the preferred BH-extracted Sachs functions with A1 and PRad data.

    Layout:
      top row    : absolute GE and GM; A1 on both, PRad on GE;
      middle row : A1/BH ratios for GE and GM;
      bottom row : PRad GE/BH ratio spanning the full figure width.

    PRad measured GE only, so no artificial PRad GM panel is constructed.
    The bottom PRad panel uses the actual PRad low-Q2 support rather than the
    full BH-extraction x range, making the percent-level shape comparison
    visible.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    qdata = np.concatenate([
        spec["data"]["t_abs"].to_numpy(float)
        for spec in selected_specs if len(spec["data"])
    ])
    qmax = 1.0
    q = np.linspace(0.0, qmax, 900)

    ge, ge_err = _sachs_band_from_result(fit, family, q, "GE")
    gm, gm_err = _sachs_band_from_result(fit, family, q, "GM")
    yahl_ge, yahl_gm = yahl18_sachs(q)

    a1 = bernauer_rosenbluth_data()
    a1 = a1.loc[a1["Q2"] <= qmax].copy()
    prad = load_prad_normalized_ge()
    prad = prad.loc[prad["Q2_GeV2"] <= qmax].copy()

    fig = plt.figure(figsize=(13.2, 11.2))
    gs = fig.add_gridspec(
        3, 2,
        height_ratios=[2.15, 1.0, 1.15],
        hspace=0.10,
        wspace=0.18,
    )
    ax_ge = fig.add_subplot(gs[0, 0])
    ax_gm = fig.add_subplot(gs[0, 1])
    ax_a1_ge = fig.add_subplot(gs[1, 0], sharex=ax_ge)
    ax_a1_gm = fig.add_subplot(gs[1, 1], sharex=ax_gm)
    ax_prad = fig.add_subplot(gs[2, :])

    qlo = float(np.nanmin(qdata))
    qhi = float(np.nanmax(qdata))

    # --------------------------
    # Absolute GE
    # --------------------------
    line_ge, = ax_ge.plot(
        q, ge, linewidth=2.0, label="BH extraction"
    )
    ax_ge.fill_between(
        q, ge - ge_err, ge + ge_err,
        alpha=0.18, linewidth=0.0,
        label="68% Hessian band",
    )
    ax_ge.plot(
        q, yahl_ge, linewidth=1.5, linestyle="--", label="YAHL18"
    )
    ax_ge.errorbar(
        a1["Q2"], a1["GE"], yerr=a1["GE_err"],
        fmt="o", fillstyle="none", markersize=4.0, capsize=2,
        linewidth=0.9, label="A1/Bernauer Rosenbluth",
    )
    for ebeam, group in prad.groupby("beam_energy_MeV", sort=True):
        ax_ge.errorbar(
            group["Q2_GeV2"], group["GE"], yerr=group["GE_total"],
            fmt="o", fillstyle="none", markersize=4.2, capsize=2,
            linewidth=0.9,
            label=f"PRad {float(ebeam)/1000.0:.3g} GeV",
        )
    #endfor
    ax_ge.axvspan(qlo, qhi, alpha=0.06)
    ax_ge.set_xlim(0.0, qmax)
    ax_ge.set_ylabel(r"$G_E^p$")
    ax_ge.grid(alpha=0.2)
    ax_ge.legend(fontsize=8.5, loc="best")

    # --------------------------
    # Absolute GM
    # --------------------------
    ax_gm.plot(q, gm, linewidth=2.0, label="BH extraction")
    ax_gm.fill_between(
        q, gm - gm_err, gm + gm_err,
        alpha=0.18, linewidth=0.0,
        label="68% Hessian band",
    )
    ax_gm.plot(
        q, yahl_gm, linewidth=1.5, linestyle="--", label="YAHL18"
    )
    ax_gm.errorbar(
        a1["Q2"], a1["GM"], yerr=a1["GM_err"],
        fmt="o", fillstyle="none", markersize=4.0, capsize=2,
        linewidth=0.9, label="A1/Bernauer Rosenbluth",
    )
    ax_gm.axvspan(qlo, qhi, alpha=0.06)
    ax_gm.set_xlim(0.0, qmax)
    ax_gm.set_ylabel(r"$G_M^p$")
    ax_gm.grid(alpha=0.2)

    # --------------------------
    # A1/BH ratios
    # --------------------------
    for which, col, errcol, rax in [
        ("GE", "GE", "GE_err", ax_a1_ge),
        ("GM", "GM", "GM_err", ax_a1_gm),
    ]:
        qpts = a1["Q2"].to_numpy(float)
        bh_pts, _ = _sachs_band_from_result(
            fit, family, qpts, which
        )
        ratio = a1[col].to_numpy(float) / bh_pts
        ratio_err = a1[errcol].to_numpy(float) / bh_pts
        rax.errorbar(
            qpts, ratio, yerr=ratio_err,
            fmt="o", fillstyle="none", markersize=4.0,
            capsize=2, linewidth=0.9, label="A1/BH",
        )
        yahl_curve = yahl_ge if which == "GE" else yahl_gm
        bh_curve = ge if which == "GE" else gm
        rax.plot(
            q, yahl_curve / bh_curve,
            linewidth=1.5, linestyle="--", label="YAHL18/BH",
        )

        central = ge if which == "GE" else gm
        sigma = ge_err if which == "GE" else gm_err
        rel_bh = np.divide(
            sigma, central,
            out=np.full_like(sigma, np.nan, dtype=float),
            where=np.abs(central) > 1.0e-12,
        )
        rax.fill_between(
            q, 1.0 - rel_bh, 1.0 + rel_bh,
            alpha=0.18, linewidth=0.0,
        )
        rax.axhline(1.0, linewidth=0.9, linestyle="--")
        rax.axvspan(qlo, qhi, alpha=0.06)
        rax.set_xlim(0.0, qmax)
        rax.set_ylim(0.90, 1.10)
        rax.grid(alpha=0.2)
        rax.set_xlabel(r"$Q^2=|t|$ (GeV$^2$)")
    #endfor
    ax_a1_ge.set_ylabel(r"$G_E^{\rm elastic}/G_E^{\rm BH}$")
    ax_a1_ge.legend(fontsize=8.0, loc="best")
    ax_a1_gm.set_ylabel(r"$G_M^{\rm elastic}/G_M^{\rm BH}$")
    ax_a1_gm.legend(fontsize=8.0, loc="best")

    # --------------------------
    # PRad/BH electric ratio, full-width low-Q2 panel
    # --------------------------
    q_prad_max = 1.0
    q_prad_min = 1.0e-4
    q_prad_band = np.geomspace(q_prad_min, q_prad_max, 500)
    bh_prad_band, bh_prad_band_err = _sachs_band_from_result(
        fit, family, q_prad_band, "GE"
    )
    rel_prad_band = np.divide(
        bh_prad_band_err,
        bh_prad_band,
        out=np.full_like(bh_prad_band_err, np.nan, dtype=float),
        where=np.abs(bh_prad_band) > 1.0e-12,
    )
    ax_prad.fill_between(
        q_prad_band,
        1.0 - rel_prad_band,
        1.0 + rel_prad_band,
        alpha=0.18,
        linewidth=0.0,
        label="BH 68% Hessian band",
    )

    prad_ratio_rows = []
    for ebeam, group in prad.groupby("beam_energy_MeV", sort=True):
        qpts = group["Q2_GeV2"].to_numpy(float)
        bh_pts, bh_pts_err = _sachs_band_from_result(
            fit, family, qpts, "GE"
        )
        ratio = group["GE"].to_numpy(float) / bh_pts
        ratio_stat = group["GE_stat"].to_numpy(float) / bh_pts
        ratio_syst = group["GE_syst"].to_numpy(float) / bh_pts
        ratio_total = np.hypot(ratio_stat, ratio_syst)

        ax_prad.errorbar(
            qpts, ratio, yerr=ratio_total,
            fmt="o", fillstyle="none", markersize=4.3,
            capsize=2, linewidth=0.9,
            label=f"PRad {float(ebeam)/1000.0:.3g} GeV / BH",
        )

        for i in range(len(group)):
            prad_ratio_rows.append({
                "beam_energy_MeV": int(ebeam),
                "Q2_GeV2": float(qpts[i]),
                "PRad_GE": float(group.iloc[i]["GE"]),
                "PRad_GE_stat": float(group.iloc[i]["GE_stat"]),
                "PRad_GE_syst": float(group.iloc[i]["GE_syst"]),
                "PRad_GE_total": float(group.iloc[i]["GE_total"]),
                "BH_GE": float(bh_pts[i]),
                "BH_GE_hessian_error": float(bh_pts_err[i]),
                "PRad_over_BH": float(ratio[i]),
                "PRad_over_BH_stat": float(ratio_stat[i]),
                "PRad_over_BH_syst": float(ratio_syst[i]),
                "PRad_over_BH_total": float(ratio_total[i]),
            })
        #endfor
    #endfor

    ax_prad.axhline(1.0, linewidth=0.9, linestyle="--")
    ax_prad.set_xscale("log")
    ax_prad.set_xlim(q_prad_min, q_prad_max)
    ax_prad.set_ylim(0.985, 1.040)
    ax_prad.set_ylabel(r"$G_E^{\rm PRad}/G_E^{\rm BH}$")
    ax_prad.set_xlabel(r"$Q^2$ (GeV$^2$)")
    ax_prad.grid(alpha=0.2, which="both")
    ax_prad.legend(fontsize=8.5, ncol=3, loc="best")

    family_label = str(family)
    fig.suptitle(
        f"Preferred BH extraction vs YAHL18, A1 and PRad elastic form factors "
        f"({family_label})",
        y=0.992,
    )
    fig.subplots_adjust(
        top=0.94,
        bottom=0.07,
        left=0.075,
        right=0.985,
    )
    fig.savefig(
        outdir / "01_preferred_GE_GM_vs_YAHL18_A1_PRad.png",
        dpi=300,
    )
    plt.close(fig)

    # Save A1 ratios used in the middle panels.
    a1_rows = []
    for which, col, errcol in [
        ("GE", "GE", "GE_err"),
        ("GM", "GM", "GM_err"),
    ]:
        qpts = a1["Q2"].to_numpy(float)
        bh_pts, bh_pts_err = _sachs_band_from_result(
            fit, family, qpts, which
        )
        for i in range(len(a1)):
            a1_rows.append({
                "form_factor": col,
                "Q2": float(qpts[i]),
                "A1_value": float(a1.iloc[i][col]),
                "A1_error": float(a1.iloc[i][errcol]),
                "BH_value": float(bh_pts[i]),
                "BH_hessian_error": float(bh_pts_err[i]),
                "A1_over_BH": float(a1.iloc[i][col] / bh_pts[i]),
                "A1_over_BH_error": float(
                    a1.iloc[i][errcol] / bh_pts[i]
                ),
            })
        #endfor
    #endfor
    pd.DataFrame(a1_rows).to_csv(
        outdir / "preferred_GE_GM_vs_A1_Bernauer_ratios.csv",
        index=False,
    )
    pd.DataFrame(prad_ratio_rows).to_csv(
        outdir / "preferred_GE_vs_PRad_ratios.csv",
        index=False,
    )

    yahl_rows = []
    for which, yahl_curve, bh_curve, bh_err in [
        ("GE", yahl_ge, ge, ge_err),
        ("GM", yahl_gm, gm, gm_err),
    ]:
        for i in range(len(q)):
            yahl_rows.append({
                "form_factor": which,
                "Q2": float(q[i]),
                "YAHL18_value": float(yahl_curve[i]),
                "BH_value": float(bh_curve[i]),
                "BH_hessian_error": float(bh_err[i]),
                "YAHL18_over_BH": float(yahl_curve[i] / bh_curve[i]),
            })
        #endfor
    #endfor
    pd.DataFrame(yahl_rows).to_csv(
        outdir / "preferred_GE_GM_vs_YAHL18_ratios.csv",
        index=False,
    )
#enddef



def run_georges_normalization_prior_scan(
        bundles: Sequence[Dict[str, object]],
        georges_key: str,
        selection: pd.DataFrame,
        family: str,
        outdir: Path) -> pd.DataFrame:
    """
    Scan Georges normalization assumptions around the production baseline.

    The production/status baseline leaves the common Georges scale
    unconstrained.  Fixed and 3--5% Gaussian-prior cases are retained only as
    diagnostics showing how the radii evolve as external normalization
    information is imposed.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    base_specs, counts = _km15_selected_specs_for_bundles(
        bundles, selection, 0.05
    )
    scan = [
        ("fixed", 0.0, False),
        ("3pct", 0.030, False),
        ("3p5pct", 0.035, False),
        ("4pct", 0.040, False),
        ("5pct", 0.050, False),
        ("free", 1.0, True),
    ]
    rows = []
    for tag, frac, free in scan:
        specs = [dict(s) for s in base_specs]
        for spec in specs:
            if str(spec["key"]) == str(georges_key):
                spec["norm_frac"] = float(frac)
                spec["unconstrained_norm"] = bool(free)
            #endif
        #endfor
        fit = fit_sachs_family_multi_measurements(
            specs,
            family=family,
            bh_cut=0.05,
            add_moradi_bh_systematic=True,
            bh_systematic_fraction=0.05,
            unconstrained_norm_keys=[georges_key] if free else None,
        )
        nname = "beta_" + re.sub(r"[^A-Za-z0-9]+", "_", str(georges_key))
        rows.append({
            "mode": tag,
            "georges_prior_fraction": np.nan if free else float(frac),
            "georges_unconstrained": bool(free),
            **counts,
            **fit,
            "georges_scale_factor": float(
                fit.get(nname + "_scale_factor", 1.0)
            ),
            "georges_beta": float(fit.get(nname, 0.0)),
        })
    #endfor
    table = pd.DataFrame(rows)
    table.to_csv(outdir / "georges_normalization_prior_scan.csv", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(14.2, 4.8))
    labels = table["mode"].astype(str).tolist()
    x = np.arange(len(table))
    axes[0].errorbar(
        x, table["rE_fm"], yerr=table["rE_fit_err_fm"],
        marker="o", linestyle="-", capsize=3,
    )
    axes[1].errorbar(
        x, table["rM_fm"], yerr=table["rM_fit_err_fm"],
        marker="o", linestyle="-", capsize=3,
    )
    axes[2].plot(x, table["georges_scale_factor"], marker="o")
    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=20, ha="right")
        ax.grid(alpha=0.2)
    #endfor
    axes[0].set_ylabel(r"$r_E$ [fm]")
    axes[1].set_ylabel(r"$r_M$ [fm]")
    axes[2].set_ylabel(r"$N_{\rm Georges}$")
    axes[2].axhline(1.0, linewidth=0.8, linestyle=":")
    fig.suptitle(
        f"All-four KM15 5% fit: Georges normalization-prior scan ({family})",
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(outdir / "01_georges_normalization_prior_scan.png", dpi=260)
    plt.close(fig)
    return table
#enddef


def run_normalization_shape_degeneracy_matrix(
        configurations: Sequence[
            Tuple[str, Sequence[Dict[str, object]], str]
        ],
        selection: pd.DataFrame,
        family: str,
        outdir: Path) -> pd.DataFrame:
    """
    Controlled normalization x form-factor-model matrix.

    Each ensemble is fit with the same Sachs family plus Moradi Fit 5 and Fit 8
    under fixed, published-constrained, and fully free experiment scales.
    This explicitly diagnoses when a surprising standalone radius is a
    normalization/G_E/G_M degeneracy rather than a robust slope measurement.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    for tag, bundles, label in configurations:
        specs0, counts = _km15_selected_specs_for_bundles(
            bundles, selection, 0.05
        )
        keys = [str(b["key"]) for b in bundles]
        for norm_mode in ["fixed", "published_constrained", "free"]:
            specs = [dict(s) for s in specs0]
            free_keys = None
            if norm_mode == "fixed":
                for spec in specs:
                    spec["norm_frac"] = 0.0
                #endfor
            elif norm_mode == "free":
                free_keys = keys
            #endif

            # Closure-family Sachs fit.
            try:
                fit = fit_sachs_family_multi_measurements(
                    specs, family=family, bh_cut=0.05,
                    add_moradi_bh_systematic=True,
                    bh_systematic_fraction=0.05,
                    unconstrained_norm_keys=free_keys,
                )
                rows.append({
                    "configuration": tag,
                    "configuration_label": label,
                    "normalization_mode": norm_mode,
                    "fit_model": f"Sachs {family}",
                    **counts, **fit,
                })
            except Exception as exc:
                rows.append({
                    "configuration": tag,
                    "configuration_label": label,
                    "normalization_mode": norm_mode,
                    "fit_model": f"Sachs {family}",
                    **counts, "valid": False, "failure": str(exc),
                })
            #endtry

            # Moradi Fit 5 and Fit 8 benchmarks under the same norm treatment.
            for fit_label, kind in [("Moradi Fit 5", "dipole"),
                                    ("Moradi Fit 8", "fit8_f2_kelly")]:
                try:
                    fr = fit_multi_measurements(
                        specs, kind=kind, fit_name=fit_label,
                        bh_cut=0.05, add_moradi_bh_systematic=True,
                        unconstrained_norm_keys=free_keys,
                    )
                    rows.append({
                        "configuration": tag,
                        "configuration_label": label,
                        "normalization_mode": norm_mode,
                        "fit_model": fit_label,
                        **counts,
                        "valid": bool(fr.success),
                        "N": int(fr.npts),
                        "chi2": float(fr.chi2),
                        "ndof": int(fr.ndof),
                        "chi2_ndof": float(fr.chi2_ndof),
                        "rE_fm": float(fr.rE_fm),
                        "rE_fit_err_fm": float(fr.rE_err_fm),
                        "rM_fm": float(fr.rM_fm),
                        "rM_fit_err_fm": float(fr.rM_err_fm),
                        **dict(fr.meta),
                    })
                except Exception as exc:
                    rows.append({
                        "configuration": tag,
                        "configuration_label": label,
                        "normalization_mode": norm_mode,
                        "fit_model": fit_label,
                        **counts, "valid": False, "failure": str(exc),
                    })
                #endtry
            #endfor
        #endfor
    #endfor
    table = pd.DataFrame(rows)
    table.to_csv(outdir / "normalization_shape_degeneracy_matrix.csv", index=False)

    # One concise figure per configuration.
    for tag, _, label in configurations:
        p = table.loc[
            (table["configuration"] == tag)
            & table.get("valid", pd.Series(False, index=table.index)).fillna(False).astype(bool)
        ].copy()
        if len(p) == 0:
            continue
        #endif
        modes = ["fixed", "published_constrained", "free"]
        models = [f"Sachs {family}", "Moradi Fit 5", "Moradi Fit 8"]
        fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.0))
        x = np.arange(len(modes))
        offsets = np.linspace(-0.20, 0.20, len(models))
        for off, model in zip(offsets, models):
            pm = p.loc[p["fit_model"] == model].set_index("normalization_mode")
            vals = pm.reindex(modes)
            ok = np.isfinite(pd.to_numeric(vals["rE_fm"], errors="coerce"))
            axes[0].errorbar(
                x[ok] + off, vals.loc[np.asarray(modes)[ok], "rE_fm"],
                yerr=vals.loc[np.asarray(modes)[ok], "rE_fit_err_fm"],
                marker="o", linestyle="none", capsize=3, label=model,
            )
            okm = np.isfinite(pd.to_numeric(vals["rM_fm"], errors="coerce"))
            axes[1].errorbar(
                x[okm] + off, vals.loc[np.asarray(modes)[okm], "rM_fm"],
                yerr=vals.loc[np.asarray(modes)[okm], "rM_fit_err_fm"],
                marker="o", linestyle="none", capsize=3, label=model,
            )
        #endfor
        for ax, ylabel in zip(axes, [r"$r_E$ [fm]", r"$r_M$ [fm]"]):
            ax.set_xticks(x)
            ax.set_xticklabels(["all fixed", "published priors", "all free"])
            ax.set_ylabel(ylabel)
            ax.grid(axis="y", alpha=0.2)
        #endfor
        handles, labs = axes[0].get_legend_handles_labels()
        fig.legend(handles, labs, loc="upper center", ncol=3,
                   bbox_to_anchor=(0.5, 0.96))
        fig.suptitle(
            f"{label}: normalization--shape degeneracy at KM15 5%",
            y=0.995,
        )
        fig.tight_layout(rect=(0, 0, 1, 0.88))
        fig.savefig(
            outdir / f"{tag}_normalization_shape_degeneracy.png", dpi=260
        )
        plt.close(fig)
    #endfor
    return table
#enddef



def run_paired_nested_threshold_correlation_study(
        bundles: Sequence[Dict[str, object]],
        selection: pd.DataFrame,
        family: str,
        outdir: Path,
        nreplicas: int = 200,
        seed: int = 20260904) -> pd.DataFrame:
    """
    Quantify statistical covariance among the nested 3--7% BH-purity samples.

    A single fluctuated realization of every underlying cross-section point is
    generated per replica.  The SAME point fluctuation is reused at every
    threshold in which that point appears.  This preserves the strong
    correlations implied by 3% subset 4% subset ... subset 7%.

    The central threshold excursion remains a methodological stability
    variation.  This study answers the separate question: how statistically
    significant is each excursion once the nesting covariance is respected?
    """
    outdir.mkdir(parents=True, exist_ok=True)
    thresholds = np.asarray([0.03, 0.04, 0.05, 0.06, 0.07], dtype=float)
    rng = np.random.default_rng(int(seed))

    selected_by_t = {}
    experimental_sigma = {}
    master_z_size = {}
    for bundle in bundles:
        key = str(bundle["key"])
        master_z_size[key] = len(bundle["all_data"])
        experimental_sigma[key] = {}
        for t in thresholds:
            d = select_bundle_from_external_model(
                bundle, selection, "km15", float(t)
            )
            selected_by_t[(key, float(t))] = d
            experimental_sigma[key][float(t)] = dataset_point_errors(
                d, str(bundle["kind"]), 0.0, False
            )
        #endfor
    #endfor

    # Observed central fits.
    central = {}
    for t in thresholds:
        specs, _ = _km15_selected_specs_for_bundles(
            bundles, selection, float(t)
        )
        central[float(t)] = fit_sachs_family_multi_measurements(
            specs, family=family, bh_cut=float(t),
            add_moradi_bh_systematic=True,
            bh_systematic_fraction=0.05,
        )
    #endfor

    rows = []
    for irep in range(max(1, int(nreplicas))):
        z_by_key = {
            key: rng.normal(size=n)
            for key, n in master_z_size.items()
        }
        for t in thresholds:
            specs = []
            for bundle in bundles:
                key = str(bundle["key"])
                d0 = selected_by_t[(key, float(t))]
                d = d0.copy()
                # selected frames preserve source-row indices from all_data.
                idx = d.index.to_numpy(int)
                sigma = experimental_sigma[key][float(t)]
                d["xs"] = (
                    d0["xs"].to_numpy(float)
                    + sigma * z_by_key[key][idx]
                )
                specs.append(bundle_to_measurement_spec(bundle, d))
            #endfor
            try:
                fit = fit_sachs_family_multi_measurements(
                    specs, family=family, bh_cut=float(t),
                    add_moradi_bh_systematic=True,
                    bh_systematic_fraction=0.05,
                )
                rows.append({
                    "replica": int(irep),
                    "threshold": float(t),
                    "threshold_percent": 100.0 * float(t),
                    "valid": bool(fit.get("valid", False)),
                    "rE_fm": float(fit.get("rE_fm", np.nan)),
                    "rM_fm": float(fit.get("rM_fm", np.nan)),
                })
            except Exception:
                rows.append({
                    "replica": int(irep),
                    "threshold": float(t),
                    "threshold_percent": 100.0 * float(t),
                    "valid": False,
                    "rE_fm": np.nan,
                    "rM_fm": np.nan,
                })
            #endtry
        #endfor
    #endfor

    rep = pd.DataFrame(rows)
    rep.to_csv(outdir / "paired_nested_threshold_replicas.csv", index=False)

    summary_rows = []
    for quantity in ["rE", "rM"]:
        piv = rep.loc[rep["valid"].astype(bool)].pivot(
            index="replica", columns="threshold", values=f"{quantity}_fm"
        ).reindex(columns=thresholds)
        corr = piv.corr()
        cov = piv.cov()
        corr.to_csv(outdir / f"{quantity}_threshold_correlation_matrix.csv")
        cov.to_csv(outdir / f"{quantity}_threshold_covariance_matrix_fm2.csv")

        nominal = 0.05
        for t in thresholds:
            paired = piv[[nominal, float(t)]].dropna()
            delta_rep = (
                paired[float(t)].to_numpy(float)
                - paired[nominal].to_numpy(float)
            )
            observed_delta = (
                float(central[float(t)][f"{quantity}_fm"])
                - float(central[nominal][f"{quantity}_fm"])
            )
            sigma_delta = (
                float(np.std(delta_rep, ddof=1))
                if len(delta_rep) > 1 else np.nan
            )
            summary_rows.append({
                "quantity": quantity,
                "threshold": float(t),
                "threshold_percent": 100.0 * float(t),
                "N_paired_replicas": int(len(paired)),
                "observed_delta_from_5pct_fm": observed_delta,
                "paired_stat_sigma_of_delta_fm": sigma_delta,
                "quadrature_excess_shift_fm": (
                    math.copysign(
                        math.sqrt(max(0.0, observed_delta**2 - sigma_delta**2)),
                        observed_delta,
                    )
                    if np.isfinite(sigma_delta) else np.nan
                ),
                "abs_quadrature_excess_shift_fm": (
                    math.sqrt(max(0.0, observed_delta**2 - sigma_delta**2))
                    if np.isfinite(sigma_delta) else np.nan
                ),
                "observed_delta_over_paired_stat_sigma": (
                    observed_delta / sigma_delta
                    if np.isfinite(sigma_delta) and sigma_delta > 0 else np.nan
                ),
                "correlation_with_5pct": (
                    float(corr.loc[float(t), nominal])
                    if float(t) in corr.index and nominal in corr.columns
                    else np.nan
                ),
            })
        #endfor
    #endfor

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(
        outdir / "paired_nested_threshold_difference_significance.csv",
        index=False,
    )

    residual_rows = []
    for quantity in ["rE", "rM"]:
        q = summary.loc[
            (summary["quantity"] == quantity)
            & summary["threshold"].between(0.03, 0.07)
            & ~np.isclose(summary["threshold"].to_numpy(float), 0.05)
        ].copy()
        residual_rows.append({
            "quantity": quantity,
            "conservative_raw_envelope_fm": (
                float(np.nanmax(np.abs(
                    q["observed_delta_from_5pct_fm"].to_numpy(float)
                ))) if len(q) else np.nan
            ),
            "covariance_aware_quadrature_excess_envelope_fm": (
                float(np.nanmax(
                    q["abs_quadrature_excess_shift_fm"].to_numpy(float)
                )) if len(q) else np.nan
            ),
            "prescription_note": (
                "Diagnostic only: at each nested threshold use "
                "sqrt(max(observed_shift^2 - paired_stat_sigma^2, 0)), then "
                "take the 3--7% envelope. The raw envelope remains production "
                "default until the prescription is finalized."
            ),
        })
    #endfor
    pd.DataFrame(residual_rows).to_csv(
        outdir / "threshold_systematic_raw_vs_covariance_aware.csv",
        index=False,
    )

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.9))
    for ax, quantity, ylabel in [
        (axes[0], "rE", r"$\Delta r_E$ from 5% (fm)"),
        (axes[1], "rM", r"$\Delta r_M$ from 5% (fm)"),
    ]:
        s = summary.loc[summary["quantity"] == quantity]
        ax.errorbar(
            s["threshold_percent"],
            s["observed_delta_from_5pct_fm"],
            yerr=s["paired_stat_sigma_of_delta_fm"],
            marker="o", linestyle="-", capsize=3,
        )
        ax.axhline(0.0, linewidth=0.8, linestyle="--")
        ax.set_xlabel("BH-purity threshold (%)")
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.2)
    #endfor
    fig.suptitle(
        "Nested BH-threshold shifts with paired statistical covariance",
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(outdir / "01_paired_threshold_shift_covariance.png", dpi=260)
    plt.close(fig)
    return summary
#enddef


def run_unified_km15_final_analysis(
        jo_bundle: Dict[str, object],
        defurne_bundle: Dict[str, object],
        defurne2017_bundle: Dict[str, object],
        georges_bundle: Dict[str, object],
        lee_bundle: Dict[str, object],
        saylor_bundle: Dict[str, object],
        args,
        root_outdir: Path) -> Dict[str, str]:
    yahl_check = validate_yahl18_implementation()
    print(
        "[YAHL18 validation] "
        f"GE(0)={yahl_check['GE0']:.8f}, "
        f"GM(0)/mu_p={yahl_check['GM0_over_mu']:.8f}, "
        f"rE={yahl_check['rE_fm']:.6f} fm, "
        f"rM={yahl_check['rM_fm']:.6f} fm"
    )
    """
    Unified KM15 workflow with an INDEPENDENT functional-form determination for
    every requested measurement ensemble.

    The expensive Hayward/Griffioen closure study is performed ONLY for the
    complete nominal all-six production ensemble.  The resulting production family is
    then held fixed for the reduced-dataset ablation fits.  This makes changes
    in radii/Hessian information interpretable as DATASET effects rather than
    a mixture of dataset changes and independently re-selected fit families.
    """
    final_dir = root_outdir / "final_analysis"
    diagnostics_dir = final_dir / "diagnostics"
    summary_dir = final_dir / "summary"
    closure_root = diagnostics_dir / "closure" / "km15_by_ensemble"
    for directory in [final_dir, diagnostics_dir, summary_dir, closure_root]:
        directory.mkdir(parents=True, exist_ok=True)
    #endfor

    selection = load_external_bh_model_selection(
        Path(args.bh_model_selection_results).expanduser().resolve()
    )

    # Current-status baseline: Georges carries an unconstrained overall
    # cross-section normalization.  A 3% fitted shift is physically modest,
    # while fixing this dataset to exactly unit normalization is not defensible
    # without a published decomposition of its common scale uncertainty.
    georges_baseline = dict(georges_bundle)
    georges_baseline["unconstrained_norm"] = True
    saylor_nominal = dict(saylor_bundle)
    saylor_nominal["direct_km15_selection"] = True
    saylor_nominal["direct_km15_delta_column"] = "bh_delta"
    saylor_nominal["fixed_t_abs_min_GeV2"] = 0.343
    saylor_nominal["label"] = (
        str(saylor_bundle["label"]) + r" ($|t|\geq0.343$ GeV$^2$)"
    )

    pre12_bundles = sort_bundles_chronologically(
        [jo_bundle, defurne_bundle, defurne2017_bundle]
    )
    pre12_with_saylor = sort_bundles_chronologically(
        [jo_bundle, defurne_bundle, defurne2017_bundle, saylor_nominal]
    )
    pre_lee_bundles = sort_bundles_chronologically(
        [jo_bundle, defurne_bundle, defurne2017_bundle,
         saylor_nominal, georges_baseline]
    )
    five_bundles = sort_bundles_chronologically(
        [jo_bundle, defurne_bundle, defurne2017_bundle,
         georges_baseline, lee_bundle]
    )
    production_bundles = sort_bundles_chronologically(
        [jo_bundle, defurne_bundle, defurne2017_bundle, saylor_nominal,
         georges_baseline, lee_bundle]
    )
    save_selected_experiment_kinematic_coverage(
        production_bundles,
        selection,
        diagnostics_dir / "selected_experiment_kinematic_coverage",
        threshold=0.05,
    )
    ensemble_defs = [
        ("lee2026_only", [lee_bundle], "CLAS12 Lee 2026"),
        ("jo2015_only", [jo_bundle], "CLAS6 Jo 2015"),
        (
            "jo2015_plus_defurne2015",
            [jo_bundle, defurne_bundle],
            "Jo 2015 + Defurne 2015",
        ),
        (
            "complete_6gev_pre_saylor",
            pre12_bundles,
            "Jo 2015 + Defurne 2015 + Defurne 2017",
        ),
        (
            "complete_6gev_with_saylor_t343",
            pre12_with_saylor,
            r"Jo + Defurne15 + Defurne17 + Saylor ($|t|\geq0.343$)",
        ),
        (
            "pre_lee_with_georges",
            pre_lee_bundles,
            r"6-GeV ensemble + Saylor ($|t|\geq0.343$) + Georges 2022",
        ),
        (
            "jo2015_plus_lee2026",
            [jo_bundle, lee_bundle],
            "CLAS6 Jo 2015 + CLAS12 Lee 2026",
        ),
        (
            "all_five",
            five_bundles,
            "Jo 2015 + Defurne 2015 + Defurne 2017 + Georges 2022 + Lee 2026",
        ),
        (
            "all_six_saylor_t343",
            production_bundles,
            r"All six: Saylor restricted to $|t|\geq0.343$ GeV$^2$",
        ),
    ]

    print("\n" + "=" * 78)
    print("[UNIFIED KM15 FINAL ANALYSIS]")
    print("[BH-purity prescription] KM15 only")
    print("[functional forms] closure-ranked ONCE on nominal all-six (Saylor |t|>=0.343); fixed for ablations")
    print(
        f"[optimizer] initial rE=rM={SACHS_INITIAL_RADIUS_FM:.2f} fm; "
        f"broad slope-equivalent radius guards "
        f"{SACHS_MIN_RADIUS_FM:.2f}--{SACHS_MAX_RADIUS_FM:.2f} fm"
    )
    print("[normalization] Georges free overall scale in all-five baseline; published nuisances retained for Jo/Defurne15/Defurne17/Lee")
    print("[benchmarks] Moradi Fit 5 and Fit 8 run across production ensembles")
    print("[family robustness] top closure-qualified families refit to real 5% data")
    print("[alternate BH selection] Kelly consistency plus AMT/Bernauer elastic-input diagnostics")
    print("[Saylor] nominal surviving sample uses KM15 5% and |t|>=0.343 GeV^2; full-sample scans retained as diagnostics")
    print("=" * 78)

    chosen = {}
    closure_bias_by_tag = {}
    family_audit_rows = []

    # ------------------------------------------------------------------
    # EXPENSIVE STEP: perform closure exactly once, on the complete
    # production ensemble.  Reduced ensembles below are ablation studies
    # and deliberately inherit this same production family.
    # ------------------------------------------------------------------
    production_tag = "all_six_saylor_t343"
    production_label = next(
        label for tag, _, label in ensemble_defs if tag == production_tag
    )
    production_cfg_bundles = next(
        bundles for tag, bundles, _ in ensemble_defs if tag == production_tag
    )
    production_closure_dir = closure_root / production_tag
    production_closure_dir.mkdir(parents=True, exist_ok=True)

    production_closure_bundles = []
    for bundle in production_cfg_bundles:
        selected5 = select_bundle_from_external_model(
            bundle, selection, "km15", 0.05
        )
        bcopy = dict(bundle)
        bcopy["set5"] = selected5.copy()
        production_closure_bundles.append(bcopy)
    #endfor

    production_ranking_path = (
        production_closure_dir / "radius_bias_mixed_family_ranking.csv"
    )

    if args.run_radius_bias_study:
        if production_ranking_path.exists():
            print(
                f"\n[production closure] {production_label}: existing closure "
                "ranking found; reusing it instead of recomputing replicas"
            )
        else:
            print(
                f"\n[production closure] {production_label}: running the ONE "
                "full bias/variance functional-form determination"
            )
            run_radius_bias_variance_study(
                production_closure_bundles, args, production_closure_dir
            )
        #endif
    #endif

    if not production_ranking_path.exists():
        raise RuntimeError(
            f"Missing production closure ranking: {production_ranking_path}. "
            "Run with --run-radius-bias-study "
            "--radius-bias-extended-truths."
        )
    #endif

    production_ranking = pd.read_csv(production_ranking_path)
    production_usable = production_ranking.loc[
        production_ranking["eligible"].astype(bool)
        & np.isfinite(
            production_ranking["combined_RMS_objective_fm"].to_numpy(float)
        )
    ].copy()
    production_usable = production_usable.sort_values(
        "combined_RMS_objective_fm", ascending=True
    ).reset_index(drop=True)
    if len(production_usable) == 0:
        raise RuntimeError(
            "No convergence-eligible Sachs family in the all-five production "
            "closure ranking."
        )
    #endif

    production_family = None
    for irank, rank_row in production_usable.iterrows():
        family = str(rank_row["family"])
        specs, counts = _km15_selected_specs_for_bundles(
            production_cfg_bundles, selection, 0.05
        )
        try:
            fit = fit_sachs_family_multi_measurements(
                specs,
                family=family,
                bh_cut=0.05,
                add_moradi_bh_systematic=True,
            )
            valid = bool(fit.get("valid", False))
        except Exception as exc:
            valid = False
            fit = {"failure": str(exc)}
        #endtry

        family_audit_rows.append({
            "configuration": production_tag,
            "configuration_label": production_label,
            "closure_rank": int(irank) + 1,
            "family": family,
            "closure_objective_fm": float(
                rank_row["combined_RMS_objective_fm"]
            ),
            "minimum_scenario_valid_fraction": float(
                rank_row["minimum_scenario_valid_fraction"]
            ),
            "global_valid_fraction": float(
                rank_row["global_valid_fraction"]
            ),
            "nominal_fit_valid": valid,
            "nominal_N": int(sum(counts.values())),
            "nominal_chi2_ndof": float(
                fit.get("chi2_ndof", np.nan)
            ),
            "nominal_rE_fm": float(fit.get("rE_fm", np.nan)),
            "nominal_rM_fm": float(fit.get("rM_fm", np.nan)),
            "nominal_failure": str(fit.get("failure", "")),
            "family_source": "all_five_production_closure",
        })
        if valid:
            production_family = family
            break
        #endif
    #endfor

    if production_family is None:
        raise RuntimeError(
            "No closure-ranked all-five family produced a valid nominal "
            "simultaneous free-GE/free-GM fit."
        )
    #endif

    production_bias_summary = summarize_bias_for_family(
        production_closure_dir / "radius_bias_variance_study.csv",
        production_family,
    )
    print(
        f"[production closure] selected family = {production_family}; "
        f"closure bias sys rE="
        f"{production_bias_summary['rE_RMS_bias_fm']:.5f} fm, "
        f"rM={production_bias_summary['rM_RMS_bias_fm']:.5f} fm"
    )

    # ------------------------------------------------------------------
    # CHEAP STEP: hold the production family fixed and refit each dataset
    # combination.  These are controlled ablations, not independent
    # functional-form determinations.
    # ------------------------------------------------------------------
    for tag, cfg_bundles, label in ensemble_defs:
        chosen[tag] = production_family

        # Use the all-five closure bias only for the actual production result.
        # For reduced ensembles it would not be an independently established
        # method uncertainty, so leave it undefined rather than mislabel it.
        if tag == production_tag:
            closure_bias_by_tag[tag] = dict(production_bias_summary)
        else:
            closure_bias_by_tag[tag] = {
                "rE_RMS_bias_fm": np.nan,
                "rM_RMS_bias_fm": np.nan,
            }
        #endif

        specs, counts = _km15_selected_specs_for_bundles(
            cfg_bundles, selection, 0.05
        )
        try:
            fit = fit_sachs_family_multi_measurements(
                specs,
                family=production_family,
                bh_cut=0.05,
                add_moradi_bh_systematic=True,
            )
            valid = bool(fit.get("valid", False))
        except Exception as exc:
            valid = False
            fit = {"failure": str(exc)}
        #endtry

        family_audit_rows.append({
            "configuration": tag,
            "configuration_label": label,
            "closure_rank": 1,
            "family": production_family,
            "closure_objective_fm": float(
                production_usable.iloc[0]["combined_RMS_objective_fm"]
            ),
            "minimum_scenario_valid_fraction": float(
                production_usable.iloc[0]["minimum_scenario_valid_fraction"]
            ),
            "global_valid_fraction": float(
                production_usable.iloc[0]["global_valid_fraction"]
            ),
            "nominal_fit_valid": valid,
            "nominal_N": int(sum(counts.values())),
            "nominal_chi2_ndof": float(fit.get("chi2_ndof", np.nan)),
            "nominal_rE_fm": float(fit.get("rE_fm", np.nan)),
            "nominal_rM_fm": float(fit.get("rM_fm", np.nan)),
            "nominal_failure": str(fit.get("failure", "")),
            "family_source": (
                "all_five_production_closure"
                if tag == production_tag
                else "fixed_from_all_five_for_ablation"
            ),
        })

        print(
            f"[fixed-family ablation] {label}: family={production_family}, "
            f"N={int(sum(counts.values()))}, valid={valid}, "
            f"chi2/dof={float(fit.get('chi2_ndof', np.nan)):.3f}, "
            f"rE={float(fit.get('rE_fm', np.nan)):.5f}, "
            f"rM={float(fit.get('rM_fm', np.nan)):.5f}"
        )
    #endfor

    pd.DataFrame(family_audit_rows).to_csv(
        summary_dir / "km15_family_selection_by_ensemble.csv",
        index=False,
    )
    with open(
            summary_dir / "chosen_sachs_family_by_ensemble.txt", "w"
    ) as fout:
        for tag, family in chosen.items():
            value = family if family is not None else "UNRESOLVED"
            fout.write(f"{tag}={value}\n")
        #endfor
    #endwith

    # Run each threshold scan with the common all-five production family.
    rotation_root = final_dir / "km15_ensemble_rotation"
    rotation_root.mkdir(parents=True, exist_ok=True)
    nominal_rows = []
    for tag, cfg_bundles, label in ensemble_defs:
        family = chosen[tag]
        cfg_dir = rotation_root / tag
        cfg_dir.mkdir(parents=True, exist_ok=True)

        if family is None:
            pd.DataFrame([{
                "configuration": tag,
                "configuration_label": label,
                "status": "UNRESOLVED_FREE_GE_GM",
                "reason": (
                    "No closure-ranked mixed Sachs family produced a valid "
                    "nominal 5% simultaneous free-GE/free-GM fit."
                ),
            }]).to_csv(
                cfg_dir / "free_ge_gm_status.csv",
                index=False,
            )
            print(
                f"[KM15 ensemble rotation] {label}: skipping simultaneous "
                "free-GE/free-GM threshold scan (unresolved nominal fit)."
            )
            continue
        #endif

        threshold_rows = []

        # Dense near the nominal BH-pure region, sparse into the DVCS-contaminated
        # region.  The final point is chosen above the largest finite KM15 delta
        # so it includes essentially every point for which model selection exists.
        scan_thresholds = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07,
                           0.075, 0.10, 0.15, 0.20, 0.30, 0.50, 0.75, 1.00]
        finite_deltas = []
        for bundle in cfg_bundles:
            key = str(bundle["key"])
            sel = selection.loc[
                (selection["dataset"] == key)
                & np.isfinite(selection["delta_bh_km15"].to_numpy(float))
            ]
            finite_deltas.extend(sel["delta_bh_km15"].to_numpy(float).tolist())
        #endfor
        if finite_deltas:
            scan_thresholds.append(max(finite_deltas) * 1.001)
        #endif
        scan_thresholds = sorted(set(float(x) for x in scan_thresholds))

        for threshold in scan_thresholds:
            specs, counts = _km15_selected_specs_for_bundles(
                cfg_bundles, selection, float(threshold)
            )
            fit = fit_sachs_family_multi_measurements(
                specs,
                family=family,
                bh_cut=float(threshold),
                add_moradi_bh_systematic=True,
                bh_systematic_fraction=0.05,
            )
            bias_summary = closure_bias_by_tag.get(tag, {})
            row = {
                "configuration": tag,
                "configuration_label": label,
                "selected_family": family,
                "threshold": float(threshold),
                "threshold_percent": 100.0 * float(threshold),
                "rE_closure_bias_sys_fm": float(
                    bias_summary.get("rE_RMS_bias_fm", np.nan)
                ),
                "rM_closure_bias_sys_fm": float(
                    bias_summary.get("rM_RMS_bias_fm", np.nan)
                ),
                **counts,
                **fit,
            }
            threshold_rows.append(row)
            if abs(float(threshold) - 0.05) < 1.0e-12:
                nominal_rows.append(row)
            #endif
        #endfor

        table = pd.DataFrame(threshold_rows)

        # Quote the two immediately useful systematic scales alongside every
        # ensemble result: closure/extrapolation bias and the KM15-only maximum
        # threshold excursion over the explicit 3--7% scan.
        nominal5 = table.loc[
            np.isclose(table["threshold"].to_numpy(float), 0.05)
            & table["valid"].astype(bool)
        ]
        local = table.loc[
            (table["threshold"].to_numpy(float) >= 0.03)
            & (table["threshold"].to_numpy(float) <= 0.07)
            & table["valid"].astype(bool)
        ]
        if len(nominal5) == 1 and len(local):
            base_e = float(nominal5.iloc[0]["rE_fm"])
            base_m = float(nominal5.iloc[0]["rM_fm"])
            sel_sys_e = float(np.max(np.abs(
                local["rE_fm"].to_numpy(float) - base_e
            )))
            sel_sys_m = float(np.max(np.abs(
                local["rM_fm"].to_numpy(float) - base_m
            )))
        else:
            sel_sys_e = np.nan
            sel_sys_m = np.nan
        #endif
        bias_summary = closure_bias_by_tag.get(tag, {})
        bias_e = float(bias_summary.get("rE_RMS_bias_fm", np.nan))
        bias_m = float(bias_summary.get("rM_RMS_bias_fm", np.nan))
        table["rE_km15_3to7_selection_sys_fm"] = sel_sys_e
        table["rM_km15_3to7_selection_sys_fm"] = sel_sys_m
        table["rE_method_sys_fm"] = (
            float(np.hypot(sel_sys_e, bias_e))
            if np.isfinite(sel_sys_e) and np.isfinite(bias_e) else np.nan
        )
        table["rM_method_sys_fm"] = (
            float(np.hypot(sel_sys_m, bias_m))
            if np.isfinite(sel_sys_m) and np.isfinite(bias_m) else np.nan
        )
        for nrow in nominal_rows:
            if str(nrow.get("configuration")) == tag:
                nrow["rE_km15_3to7_selection_sys_fm"] = sel_sys_e
                nrow["rM_km15_3to7_selection_sys_fm"] = sel_sys_m
                nrow["rE_method_sys_fm"] = table["rE_method_sys_fm"].iloc[0]
                nrow["rM_method_sys_fm"] = table["rM_method_sys_fm"].iloc[0]
            #endif
        #endfor

        table.to_csv(cfg_dir / "km15_threshold_scan.csv", index=False)

        # Compare the fitted finite-|t| Sachs form factors themselves across
        # the 1--7% KM15 selections.  This is deliberately attached to the
        # ACTIVE independently closure-selected ensemble workflow (rather than
        # the legacy fixed-family rotation helper) so the analysis-note
        # Hessian-band products are always produced alongside the radius scan.
        band_dir = cfg_dir / "form_factor_band_stability"
        band_metrics = save_threshold_form_factor_band_study(
            bundles=cfg_bundles,
            selection=selection,
            family=family,
            threshold_table=table,
            label=label,
            outdir=band_dir,
        )
        if len(band_metrics):
            print(
                f"[KM15 form-factor stability] {label}: "
                f"{len(band_metrics)} overlap rows -> {band_dir}"
            )
        else:
            print(
                f"[KM15 form-factor stability] WARNING: {label}: "
                "no threshold-band products were generated; inspect the "
                "1--7% valid-fit rows and selected-data support."
            )
        #endif

        fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.7))
        axes[0].errorbar(
            table["threshold_percent"], table["rE_fm"],
            yerr=table["rE_fit_err_fm"], marker="o", linestyle="-",
        )
        axes[1].errorbar(
            table["threshold_percent"], table["rM_fm"],
            yerr=table["rM_fit_err_fm"], marker="o", linestyle="-",
        )
        axes[2].plot(
            table["threshold_percent"], table["chi2_ndof"],
            marker="o", linestyle="-",
        )
        axes[0].set_ylabel(r"$r_E$ (fm)")
        axes[1].set_ylabel(r"$r_M$ (fm)")
        axes[2].set_ylabel(r"$\chi^2/\mathrm{dof}$")
        for ax in axes:
            ax.set_xlabel(r"$|1-R_{\rm BH}^{\rm KM15}|$ threshold (%)")
            ax.axvline(5.0, linewidth=0.8, linestyle=":")
            ax.grid(alpha=0.2)
        #endfor
        fig.suptitle(
            f"{label}: KM15 threshold stability; closure-selected {family}",
            y=0.995,
        )
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        fig.savefig(cfg_dir / "01_km15_threshold_stability.png", dpi=260)
        plt.close(fig)
    #endfor

    nominal = pd.DataFrame(nominal_rows)
    nominal.to_csv(
        rotation_root / "km15_5pct_ensemble_summary.csv", index=False
    )

    if len(nominal):
        incremental_order = [
            "jo2015_plus_defurne2015",
            "complete_6gev_pre_saylor",
            "complete_6gev_with_saylor_t343",
            "pre_lee_with_georges",
            "all_six_saylor_t343",
        ]
        inc = nominal.loc[
            nominal["configuration"].astype(str).isin(incremental_order)
        ].copy()
        inc["_order"] = inc["configuration"].map(
            {k: i for i, k in enumerate(incremental_order)}
        )
        inc = inc.sort_values("_order").drop(columns=["_order"])
        if len(inc):
            inc["family_policy"] = (
                "fixed to nominal all-six closure-selected production family"
            )
            inc["delta_rE_from_previous_fm"] = inc["rE_fm"].diff()
            inc["delta_rM_from_previous_fm"] = inc["rM_fm"].diff()
            inc["rE_fit_error_ratio_to_previous"] = (
                inc["rE_fit_err_fm"] / inc["rE_fit_err_fm"].shift(1)
            )
            inc["rM_fit_error_ratio_to_previous"] = (
                inc["rM_fit_err_fm"] / inc["rM_fit_err_fm"].shift(1)
            )
            inc.to_csv(
                summary_dir / "incremental_dataset_radius_impact.csv",
                index=False,
            )
            print("\n[incremental dataset impact at KM15 5%]")
            show = [
                "configuration_label", "selected_family", "N",
                "rE_fm", "rE_fit_err_fm",
                "rM_fm", "rM_fit_err_fm",
                "delta_rE_from_previous_fm",
                "delta_rM_from_previous_fm",
            ]
            show = [c for c in show if c in inc.columns]
            print(inc[show].to_string(index=False))
        #endif
    #endif

    # The all-five row is now the actual free-Georges production/status
    # baseline.  Save a narrow artifact that carries central values, fit
    # uncertainties, closure bias, 3--7% selection terms, and their current
    # quadrature method combination under one self-consistent prescription.
    if len(nominal):
        preferred = nominal.loc[
            nominal["configuration"].astype(str) == "all_six_saylor_t343"
        ].copy()
        if len(preferred) == 1:
            status_dir = diagnostics_dir / "preferred_current_status"
            status_dir.mkdir(parents=True, exist_ok=True)
            preferred["georges_normalization_treatment"] = "unconstrained"
            preferred["status_role"] = "preferred_nominal_all_six_saylor_t343"
            preferred.to_csv(
                status_dir / "preferred_all_six_saylor_t343_georges_free_with_method_systematics.csv",
                index=False,
            )
        #endif
    #endif

    if len(nominal):
        fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.2))
        x = np.arange(len(nominal))
        for ax, radius, fiterr, biaserr, ylabel in [
            (
                axes[0], "rE_fm", "rE_fit_err_fm",
                "rE_closure_bias_sys_fm", r"$r_E$ [fm]",
            ),
            (
                axes[1], "rM_fm", "rM_fit_err_fm",
                "rM_closure_bias_sys_fm", r"$r_M$ [fm]",
            ),
        ]:
            outer = np.hypot(
                nominal[fiterr].to_numpy(float),
                nominal[biaserr].to_numpy(float),
            )
            ax.errorbar(
                x, nominal[radius], yerr=outer,
                marker="o", linestyle="none", capsize=4,
                label="fit + closure-bias (quadrature)",
            )
            ax.errorbar(
                x, nominal[radius], yerr=nominal[fiterr],
                marker="o", linestyle="none", capsize=2,
                label="fit uncertainty",
            )
            ax.set_xticks(x)
            ax.set_xticklabels(
                nominal["configuration_label"].astype(str),
                rotation=18, ha="right",
            )
            ax.set_ylabel(ylabel)
            ax.grid(axis="y", alpha=0.2)
        #endfor
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(
            handles, labels, loc="upper center", ncol=2,
            bbox_to_anchor=(0.5, 0.955),
        )
        fig.suptitle(
            "KM15 5% ensemble radii: common nominal all-six family",
            y=0.995,
        )
        fig.tight_layout(rect=(0, 0, 1, 0.90))
        fig.savefig(
            rotation_root / "00_km15_5pct_ensemble_radii_with_bias.png",
            dpi=180,
        )
        plt.close(fig)
    #endif

    print("\n[KM15 ensemble rotation] common nominal all-six production family")
    cols = [
        "configuration_label", "selected_family", "N", "chi2_ndof",
        "rE_fm", "rE_fit_err_fm", "rE_closure_bias_sys_fm",
        "rE_km15_3to7_selection_sys_fm", "rE_method_sys_fm",
        "rM_fm", "rM_fit_err_fm", "rM_closure_bias_sys_fm",
        "rM_km15_3to7_selection_sys_fm", "rM_method_sys_fm",
    ]
    if len(nominal):
        present = [c for c in cols if c in nominal.columns]
        print(nominal[present].to_string(index=False))
    else:
        print("[KM15 ensemble rotation] no resolved nominal free-GE/free-GM fits.")
    #endif

    # Compare the leading closure-qualified families directly on the real
    # nominal data. This is a discrete-family robustness diagnostic only and is
    # deliberately not folded into the production method uncertainty.
    run_competitive_family_real_data_diagnostics(
        ensemble_defs=ensemble_defs,
        chosen=chosen,
        closure_root=closure_root,
        selection=selection,
        outdir=diagnostics_dir / "competitive_family_real_data",
        max_families=8,
    )

    if chosen.get("all_six_saylor_t343") is not None:
        production_specs, _ = _km15_selected_specs_for_bundles(
            production_bundles, selection, 0.05
        )
        try:
            save_hayward_griffioen_model_averaging_diagnostics(
                closure_table_path=(
                    closure_root / "all_six_saylor_t343"
                    / "radius_bias_variance_study.csv"
                ),
                closure_ranking_path=(
                    closure_root / "all_six_saylor_t343"
                    / "radius_bias_mixed_family_ranking.csv"
                ),
                competitive_fit_path=(
                    diagnostics_dir / "competitive_family_real_data"
                    / "all_six_saylor_t343"
                    / "competitive_family_real_data_fits.csv"
                ),
                selected_specs=production_specs,
                outdir=diagnostics_dir / "radius_model_averaging",
                max_families=8,
                n_bootstrap=20000,
            )
        except Exception as exc:
            print(
                "[radius model averaging] WARNING: "
                f"{type(exc).__name__}: {exc}"
            )
        #endtry

        try:
            run_paired_nested_threshold_correlation_study(
                bundles=production_bundles,
                selection=selection,
                family=chosen["all_six_saylor_t343"],
                outdir=diagnostics_dir / "paired_threshold_covariance",
                nreplicas=int(
                    getattr(args, "threshold_correlation_replicas", 200)
                ),
                seed=int(args.radius_bias_seed) + 707,
            )
        except Exception as exc:
            print(
                "[paired threshold covariance] WARNING: "
                f"{type(exc).__name__}: {exc}"
            )
        #endtry
    #endif

    # Preserve the original Moradi Fit-5 and Fit-8 parameterizations as a
    # common benchmark across every relevant single/combined ensemble.
    run_moradi_fit5_fit8_ensemble_benchmarks(
        jo_bundle=jo_bundle,
        defurne_bundle=defurne_bundle,
        georges_bundle=georges_baseline,
        lee_bundle=lee_bundle,
        selection=selection,
        outdir=diagnostics_dir / "moradi_fit5_fit8_ensembles",
    )

    # Diagnose the Hall-A-induced radius rotation with one common family so
    # changes are attributable to the dataset ensemble rather than to a
    # simultaneous change of extrapolation parameterization.
    if chosen.get("all_six_saylor_t343") is not None:
        run_fixed_family_ensemble_decomposition(
            jo_bundle=jo_bundle,
            defurne_bundle=defurne_bundle,
            georges_bundle=georges_baseline,
            lee_bundle=lee_bundle,
            selection=selection,
            family=chosen["all_six_saylor_t343"],
            outdir=diagnostics_dir / "hall_a_ensemble_decomposition",
        )
        norm_tension = run_normalization_tension_diagnostics(
            bundles=production_bundles,
            selection=selection,
            family=chosen["all_six_saylor_t343"],
            outdir=diagnostics_dir / "normalization_tension",
        )

        # Current-status headline: use the completely free Georges scale.
        # This does NOT alter the closure ranking or erase the fixed/constrained
        # normalization diagnostics; it simply records the user's preferred
        # provisional all-five result in a dedicated, unambiguous artifact.
        # The free scale is especially useful for the status note because the
        # fit itself determines whether Georges wants a pathological rescaling.
        try:
            if isinstance(norm_tension, pd.DataFrame):
                nt = norm_tension
            elif isinstance(norm_tension, dict) and "fit_table" in norm_tension:
                nt = norm_tension["fit_table"]
            else:
                nt = pd.DataFrame()
            #endif
            free_row = nt.loc[
                nt.get("normalization_mode", pd.Series(dtype=str)).astype(str)
                == "georges_free_baseline"
            ].copy()
            if len(free_row) == 1:
                status_dir = diagnostics_dir / "preferred_current_status"
                status_dir.mkdir(parents=True, exist_ok=True)
                free_row = free_row.copy()
                free_row["status_role"] = "preferred_provisional_all_five"
                free_row["georges_normalization_treatment"] = "unconstrained"
                free_row["note"] = (
                    "Preferred current-status result: Georges normalization "
                    "floats freely. Fixed and finite-prior Georges treatments "
                    "remain available only as normalization diagnostics."
                )
                free_row.to_csv(
                    status_dir / "preferred_all_four_georges_free.csv",
                    index=False,
                )
                print(
                    "[preferred current status] all-five with Georges free -> "
                    f"{status_dir / 'preferred_all_four_georges_free.csv'}"
                )
            #endif
        except Exception as exc:
            print(
                "[preferred current status] WARNING: could not write "
                f"free-Georges status artifact: {type(exc).__name__}: {exc}"
            )
        #endtry

        # Exact BH electromagnetic decomposition and derivative sensitivities.
        all5_specs, _ = _km15_selected_specs_for_bundles(
            production_bundles, selection, 0.05
        )
        all5_fit = fit_sachs_family_multi_measurements(
            all5_specs,
            family=chosen["all_six_saylor_t343"],
            bh_cut=0.05,
            add_moradi_bh_systematic=True,
            bh_systematic_fraction=0.05,
        )
        save_f1_f2_bh_sensitivity_diagnostics(
            bundles=production_bundles,
            selection=selection,
            fit=all5_fit,
            family=chosen["all_six_saylor_t343"],
            outdir=diagnostics_dir / "f1_f2_sensitivity",
        )
        all5_specs, _ = _km15_selected_specs_for_bundles(
            production_bundles, selection, 0.05
        )
        try:
            save_preferred_sachs_vs_elastic_data(
                fit=all5_fit, family=chosen["all_six_saylor_t343"],
                selected_specs=all5_specs,
                outdir=diagnostics_dir / "preferred_form_factors_vs_elastic",
            )
        except Exception as exc:
            print(
                "[preferred elastic comparison] WARNING: diagnostic failed; "
                "continuing without aborting the completed production fit: "
                f"{type(exc).__name__}: {exc}"
            )
        #endtry

        # E12-06-114 / Georges normalization study: fixed, plausible few-percent
        # Gaussian priors, and an unconstrained diagnostic endpoint.
        run_georges_normalization_prior_scan(
            bundles=production_bundles,
            georges_key=str(georges_bundle["key"]),
            selection=selection,
            family=chosen["all_six_saylor_t343"],
            outdir=diagnostics_dir / "georges_normalization_prior_scan",
        )

        # Explicitly expose the normalization--GE/GM--shape degeneracy that can
        # produce very small standalone radii.  Use one common Sachs family so
        # changes are not confounded by simultaneous family rotation.
        run_normalization_shape_degeneracy_matrix(
            configurations=[
                ("jo_only", [jo_bundle], "Jo 2015"),
                ("jo_plus_lee", [jo_bundle, lee_bundle], "Jo + Lee"),
                (
                    "hall_a_only",
                    [defurne_bundle, defurne2017_bundle, georges_bundle],
                    "Defurne 2015 + Defurne 2017 + Georges",
                ),
                (
                    "all_five",
                    production_bundles,
                    "All six; Saylor |t|>=0.343",
                ),
            ],
            selection=selection,
            family=chosen["all_six_saylor_t343"],
            outdir=diagnostics_dir / "normalization_shape_degeneracy",
        )

        run_model_only_kelly_bh_selector_diagnostic(
            bundles=production_bundles,
            selection=selection,
            family=chosen["all_six_saylor_t343"],
            label="All six; Saylor |t|>=0.343",
            outdir=diagnostics_dir / "bh_model_only_alternate" / "all_five",
        )
        run_data_driven_bh_deviance_for_ensemble(
            bundles=production_bundles,
            selection=selection,
            family=chosen["all_six_saylor_t343"],
            label="All six; Saylor |t|>=0.343",
            outdir=diagnostics_dir / "bh_deviance" / "all_five",
        )
    #endif

    if chosen.get("jo2015_plus_lee2026") is not None:
        run_model_only_kelly_bh_selector_diagnostic(
            bundles=[jo_bundle, lee_bundle],
            selection=selection,
            family=chosen["jo2015_plus_lee2026"],
            label="Jo + Lee",
            outdir=diagnostics_dir / "bh_model_only_alternate" / "jo_plus_lee",
        )
        run_data_driven_bh_deviance_for_ensemble(
            bundles=[jo_bundle, lee_bundle],
            selection=selection,
            family=chosen["jo2015_plus_lee2026"],
            label="Jo + Lee",
            outdir=diagnostics_dir / "bh_deviance" / "jo_plus_lee",
        )
    #endif

    # Moradi Fit 8 / Kelly-F2 Lee-only diagnostic.  This is intentionally
    # independent of whether the simultaneous free-GE/free-GM Lee extraction
    # succeeded.  In particular, when Lee-only is unresolved above, this
    # becomes the scientifically relevant test of whether the data still
    # constrain the electric sector once F2 is supplied externally.
    lee_fit8_dir = diagnostics_dir / "lee2026_fit8_kelly"
    lee_fit8_table = save_lee_fit8_kelly_threshold_scan(
        lee_bundle=lee_bundle,
        selection=selection,
        summary_dir=summary_dir,
        outdir=lee_fit8_dir,
    )
    nominal_fit8 = lee_fit8_table.loc[
        np.isclose(
            lee_fit8_table["selection_threshold"].to_numpy(float),
            0.05,
        )
    ]
    if len(nominal_fit8):
        row = nominal_fit8.iloc[0]
        print(
            f"[Lee Fit8 / Kelly F2] N={int(row['N_selected'])}, "
            f"chi2/ndof={float(row['chi2_ndof']):.6f}, "
            f"rE={float(row['rE_fm']):.6f}+/-"
            f"{float(row['rE_err_fm']):.6f} fm"
        )
    else:
        print(
            "[Lee Fit8 / Kelly F2] no nominal 5% fit was available; "
            "inspect diagnostics/lee2026_fit8_kelly."
        )
    #endif

    # Global six-dataset Saylor low-|t| removal study, following the logic
    # of arXiv:2607.04481 but adding Georges 2022 and Lee 2026 and using the
    # closure-selected production Sachs family.
    if chosen.get("all_six_saylor_t343") is not None:
        try:
            run_global_saylor_tmin_scan(
                five_bundles,
                saylor_bundle,
                selection,
                chosen["all_six_saylor_t343"],
                root_outdir,
            )
        except Exception as exc:
            print(
                "[global Saylor |t|min scan] WARNING: diagnostic failed; "
                "continuing without aborting the completed production fit: "
                f"{type(exc).__name__}: {exc}"
            )
        #endtry
    #endif

    if chosen.get("all_six_saylor_t343") is not None:
        try:
            run_all_dataset_subset_and_external_f2_diagnostic(
                production_bundles=five_bundles,
                saylor_bundle=saylor_bundle,
                selection=selection,
                family=chosen["all_six_saylor_t343"],
                outdir=diagnostics_dir / "all_dataset_subsets_free_vs_externalF2",
            )
        except Exception as exc:
            print(
                "[all-subset free-vs-external-F2] WARNING: diagnostic failed; "
                f"continuing: {type(exc).__name__}: {exc}"
            )
            print(
                "[all-subset free-vs-external-F2] Inspect "
                "diagnostics/all_dataset_subsets_free_vs_externalF2/"
                "externalF2_partial_results.csv and "
                "externalF2_failures.csv for the exact failed mode/subset."
            )
        #endtry
    else:
        print(
            "[all-subset free-vs-external-F2] skipped because the nominal "
            "restricted-Saylor production family is unresolved."
        )
    #endif

    # Saylor: keep both model-based and direct Jo-vs-Saylor diagnostics.
    # For the BH-fit portion use the all-five selected family only as a
    # diagnostic reference; it does not affect any production selection.
    if chosen.get("all_six_saylor_t343") is not None:
        run_saylor_recovery_study(
            saylor_bundle,
            args,
            root_outdir,
            family=chosen["all_six_saylor_t343"],
        )
    else:
        print(
            "[Saylor recovery] skipped family-referenced BH-fit diagnostic "
            "because all-five free-GE/free-GM extraction was unresolved."
        )
    #endif

    return chosen
#enddef




def save_jo_saylor_matched_kinematics(
        jo_bundle: Dict[str, object],
        saylor_bundle: Dict[str, object],
        outdir: Path) -> None:
    """
    Direct Jo-vs-Saylor comparison in overlapping multidimensional kinematics.

    For every Saylor point, find the nearest Jo point in standardized
    (Q2, xB, |t|, phi) space.  Two ratios are retained:

      raw ratio
          sigma_Saylor / sigma_Jo

      locally transported ratio
          [sigma_Saylor * KM15(Jo kinematics) / KM15(Saylor kinematics)]
          / sigma_Jo

    The second quantity uses KM15 only as a local interpolation/transport
    between nearby measured kinematic points.  It does NOT use KM15 to decide
    which dataset is correct.  Point-to-point experimental uncertainties from
    both measurements are propagated to the transported ratio and to the
    normalized Jo-Saylor difference.

    The presentation-facing plots emphasize the best-overlap 50% of matches,
    while CSV tables retain the full sample and fixed 25/50/75/100% overlap
    quantiles.  This makes the central tendency near unity and the substantial
    point-to-point spread visually explicit without hiding matching quality.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    jo = jo_bundle["all_data"].copy().reset_index(drop=True)
    sa = saylor_bundle["all_data"].copy().reset_index(drop=True)

    cols = ["Q2", "xB", "t_abs", "phi_deg"]
    for frame in [jo, sa]:
        for col in cols + ["xs", "km15_ep"]:
            frame[col] = pd.to_numeric(frame[col], errors="coerce")
        #endfor
    #endfor

    # Canonical point-to-point errors; correlated overall normalization
    # uncertainties are intentionally not folded into this point-pair diagnostic.
    jo["_point_err"] = dataset_point_errors(
        jo, str(jo_bundle["kind"]), 0.0, False
    )
    sa["_point_err"] = dataset_point_errors(
        sa, str(saylor_bundle["kind"]), 0.0, False
    )

    jo = jo.dropna(subset=cols + ["xs", "km15_ep", "_point_err"]).copy()
    sa = sa.dropna(subset=cols + ["xs", "km15_ep", "_point_err"]).copy()
    jo = jo.loc[
        (jo["xs"] > 0.0) & (jo["km15_ep"] > 0.0) & (jo["_point_err"] > 0.0)
    ].reset_index(drop=True)
    sa = sa.loc[
        (sa["xs"] > 0.0) & (sa["km15_ep"] > 0.0) & (sa["_point_err"] > 0.0)
    ].reset_index(drop=True)

    # Circular phi distance is represented by sin/cos coordinates.
    def features(frame):
        phi = np.deg2rad(frame["phi_deg"].to_numpy(float))
        return np.column_stack([
            frame["Q2"].to_numpy(float),
            frame["xB"].to_numpy(float),
            frame["t_abs"].to_numpy(float),
            np.sin(phi),
            np.cos(phi),
        ])
    #enddef

    def circular_phi_delta_deg(a, b):
        return np.abs((float(a) - float(b) + 180.0) % 360.0 - 180.0)
    #enddef

    combined = np.vstack([features(jo), features(sa)])
    scale = np.nanstd(combined, axis=0)
    scale[~np.isfinite(scale) | (scale <= 0.0)] = 1.0
    fj = features(jo) / scale
    fs = features(sa) / scale

    rows = []
    # Chunked brute-force matching avoids an optional scipy spatial dependency.
    for i0 in range(0, len(sa), 250):
        block = fs[i0:i0 + 250]
        d2 = np.sum(
            (block[:, None, :] - fj[None, :, :])**2, axis=2
        )
        jj = np.argmin(d2, axis=1)
        dd = np.sqrt(d2[np.arange(len(jj)), jj])

        for local, (j, dist) in enumerate(zip(jj, dd)):
            si = i0 + local
            srow = sa.iloc[si]
            jrow = jo.iloc[int(j)]

            s_xs = float(srow["xs"])
            j_xs = float(jrow["xs"])
            s_err = float(srow["_point_err"])
            j_err = float(jrow["_point_err"])
            s_model = float(srow["km15_ep"])
            j_model = float(jrow["km15_ep"])

            raw_ratio = s_xs / j_xs if j_xs != 0.0 else np.nan

            transport_factor = (
                j_model / s_model
                if np.isfinite(j_model) and np.isfinite(s_model) and s_model > 0.0
                else np.nan
            )
            s_trans = s_xs * transport_factor
            s_trans_err = s_err * abs(transport_factor)

            transported_ratio = (
                s_trans / j_xs if np.isfinite(s_trans) and j_xs != 0.0 else np.nan
            )
            transported_ratio_err = (
                abs(transported_ratio)
                * np.sqrt((s_err / s_xs)**2 + (j_err / j_xs)**2)
                if (
                    np.isfinite(transported_ratio)
                    and s_xs > 0.0 and j_xs > 0.0
                    and s_err > 0.0 and j_err > 0.0
                )
                else np.nan
            )
            denom = np.sqrt(s_trans_err**2 + j_err**2)
            normalized_difference = (
                (s_trans - j_xs) / denom
                if np.isfinite(denom) and denom > 0.0
                else np.nan
            )

            rs = s_xs / s_model
            rj = j_xs / j_model

            rows.append({
                "saylor_index": int(si),
                "jo_index": int(j),
                "match_distance": float(dist),
                "saylor_Q2": float(srow["Q2"]),
                "jo_Q2": float(jrow["Q2"]),
                "delta_Q2_abs": abs(float(srow["Q2"]) - float(jrow["Q2"])),
                "saylor_xB": float(srow["xB"]),
                "jo_xB": float(jrow["xB"]),
                "delta_xB_abs": abs(float(srow["xB"]) - float(jrow["xB"])),
                "saylor_t_abs": float(srow["t_abs"]),
                "jo_t_abs": float(jrow["t_abs"]),
                "delta_t_abs": abs(float(srow["t_abs"]) - float(jrow["t_abs"])),
                "saylor_phi_deg": float(srow["phi_deg"]),
                "jo_phi_deg": float(jrow["phi_deg"]),
                "delta_phi_deg_abs": circular_phi_delta_deg(
                    srow["phi_deg"], jrow["phi_deg"]
                ),
                "saylor_xs": s_xs,
                "jo_xs": j_xs,
                "saylor_point_err": s_err,
                "jo_point_err": j_err,
                "saylor_km15_ep": s_model,
                "jo_km15_ep": j_model,
                "saylor_data_over_km15": rs,
                "jo_data_over_km15": rj,
                "raw_ratio_saylor_over_jo": raw_ratio,
                "km15_transport_factor_saylor_to_jo": transport_factor,
                "saylor_xs_transported_to_jo": s_trans,
                "saylor_transported_point_err": s_trans_err,
                "ratio_saylor_over_jo": transported_ratio,
                "ratio_saylor_over_jo_err": transported_ratio_err,
                "normalized_difference_transported": normalized_difference,
            })
        #endfor
    #endfor

    matches = pd.DataFrame(rows)
    matches.to_csv(
        outdir / "02_jo_saylor_matched_kinematics.csv", index=False
    )

    if len(matches) == 0:
        print("[Jo/Saylor matched diagnostic] no valid matches")
        return
    #endif

    # Fixed overlap quantiles.  These are defined solely by matching distance,
    # not by the observed Jo/Saylor cross-section agreement.
    quantiles = [0.25, 0.50, 0.75, 1.00]
    summary = []
    for q in quantiles:
        cut = float(matches["match_distance"].quantile(q))
        part = matches.loc[matches["match_distance"] <= cut].copy()
        rr = part["ratio_saylor_over_jo"].to_numpy(float)
        raw = part["raw_ratio_saylor_over_jo"].to_numpy(float)
        z = part["normalized_difference_transported"].to_numpy(float)
        rr = rr[np.isfinite(rr)]
        raw = raw[np.isfinite(raw)]
        z = z[np.isfinite(z)]

        if len(rr):
            p16, p50, p84 = np.percentile(rr, [16.0, 50.0, 84.0])
        else:
            p16 = p50 = p84 = np.nan
        #endif

        summary.append({
            "match_distance_quantile": q,
            "max_match_distance": cut,
            "N": len(part),
            "median_raw_saylor_over_jo": (
                float(np.nanmedian(raw)) if len(raw) else np.nan
            ),
            "median_transported_saylor_over_jo": float(p50),
            "transported_ratio_p16": float(p16),
            "transported_ratio_p84": float(p84),
            "transported_ratio_half_68_width": (
                float(0.5 * (p84 - p16))
                if np.isfinite(p16) and np.isfinite(p84) else np.nan
            ),
            "rms_transported_ratio_minus_one": (
                float(np.sqrt(np.nanmean((rr - 1.0)**2))) if len(rr) else np.nan
            ),
            "fraction_transported_within_10pct": (
                float(np.mean(np.abs(rr - 1.0) <= 0.10)) if len(rr) else np.nan
            ),
            "fraction_transported_within_20pct": (
                float(np.mean(np.abs(rr - 1.0) <= 0.20)) if len(rr) else np.nan
            ),
            "mean_normalized_difference": (
                float(np.mean(z)) if len(z) else np.nan
            ),
            "rms_normalized_difference": (
                float(np.sqrt(np.mean(z**2))) if len(z) else np.nan
            ),
            "fraction_abs_normalized_difference_gt_2": (
                float(np.mean(np.abs(z) > 2.0)) if len(z) else np.nan
            ),
            "fraction_abs_normalized_difference_gt_3": (
                float(np.mean(np.abs(z) > 3.0)) if len(z) else np.nan
            ),
            "median_abs_delta_Q2": float(np.nanmedian(part["delta_Q2_abs"])),
            "median_abs_delta_xB": float(np.nanmedian(part["delta_xB_abs"])),
            "median_abs_delta_t": float(np.nanmedian(part["delta_t_abs"])),
            "median_abs_delta_phi_deg": float(
                np.nanmedian(part["delta_phi_deg_abs"])
            ),
        })
    #endfor
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(
        outdir / "03_jo_saylor_matched_overlap_summary.csv", index=False
    )

    # Presentation-facing subset: best-overlap half, selected only by geometry.
    q50_cut = float(matches["match_distance"].quantile(0.50))
    tight = matches.loc[matches["match_distance"] <= q50_cut].copy()
    ratio = tight["ratio_saylor_over_jo"].to_numpy(float)
    finite_ratio = ratio[np.isfinite(ratio)]
    if len(finite_ratio):
        rp16, rmed, rp84 = np.percentile(finite_ratio, [16.0, 50.0, 84.0])
    else:
        rp16 = rmed = rp84 = np.nan
    #endif

    def add_equal_population_summary(ax, x, y, nbins=8):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        good = np.isfinite(x) & np.isfinite(y)
        x = x[good]
        y = y[good]
        if len(x) < max(8, nbins):
            return
        #endif
        order = np.argsort(x)
        groups = np.array_split(order, nbins)
        xc, ym, ylo, yhi = [], [], [], []
        for g in groups:
            if len(g) == 0:
                continue
            #endif
            xx = x[g]
            yy = y[g]
            q16, q50, q84 = np.percentile(yy, [16.0, 50.0, 84.0])
            xc.append(float(np.median(xx)))
            ym.append(float(q50))
            ylo.append(float(q50 - q16))
            yhi.append(float(q84 - q50))
        #endfor
        ax.errorbar(
            xc, ym, yerr=np.vstack([ylo, yhi]),
            fmt="o-", linewidth=1.5, markersize=4.5, capsize=2.5,
            label="equal-population median and 16--84%",
        )
    #enddef

    # ------------------------------------------------------------------
    # Figure 02: the main direct-comparison figure.
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.4))
    ax = axes[0, 0]
    ax.scatter(
        tight["jo_t_abs"], tight["ratio_saylor_over_jo"],
        s=11, alpha=0.22, label="matched points",
    )
    add_equal_population_summary(
        ax, tight["jo_t_abs"], tight["ratio_saylor_over_jo"]
    )
    ax.axhline(1.0, linewidth=1.0, linestyle="--")
    ax.set_ylim(0.0, 2.5)
    ratio_t = tight["ratio_saylor_over_jo"].to_numpy(float)
    n_ratio_t_out = int(np.sum(np.isfinite(ratio_t) & ((ratio_t < 0.0) | (ratio_t > 2.5))))
    if n_ratio_t_out > 0:
        ax.text(
            0.98, 0.96, f"{n_ratio_t_out} points outside displayed y range",
            transform=ax.transAxes, ha="right", va="top", fontsize=8,
        )
    #endif
    ax.set_xlabel(r"Jo matched $|t|$ (GeV$^2$)")
    ax.set_ylabel(r"transported $\sigma_{\rm Saylor}/\sigma_{\rm Jo}$")
    ax.set_title(r"Ratio versus $|t|$")
    ax.grid(alpha=0.2)
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    ax.scatter(
        tight["jo_phi_deg"], tight["ratio_saylor_over_jo"],
        s=11, alpha=0.22,
    )
    add_equal_population_summary(
        ax, tight["jo_phi_deg"], tight["ratio_saylor_over_jo"]
    )
    ax.axhline(1.0, linewidth=1.0, linestyle="--")
    ax.set_ylim(0.0, 2.5)
    ratio_phi = tight["ratio_saylor_over_jo"].to_numpy(float)
    n_ratio_phi_out = int(np.sum(np.isfinite(ratio_phi) & ((ratio_phi < 0.0) | (ratio_phi > 2.5))))
    if n_ratio_phi_out > 0:
        ax.text(
            0.98, 0.96, f"{n_ratio_phi_out} points outside displayed y range",
            transform=ax.transAxes, ha="right", va="top", fontsize=8,
        )
    #endif
    ax.set_xlabel(r"Jo matched $\phi$ (deg)")
    ax.set_ylabel(r"transported $\sigma_{\rm Saylor}/\sigma_{\rm Jo}$")
    ax.set_title(r"Ratio versus $\phi$")
    ax.grid(alpha=0.2)

    ax = axes[1, 0]
    z = tight["normalized_difference_transported"].to_numpy(float)
    ax.scatter(tight["jo_t_abs"], z, s=11, alpha=0.24)
    ax.axhline(0.0, linewidth=1.0)
    ax.axhline(+2.0, linewidth=0.9, linestyle="--")
    ax.axhline(-2.0, linewidth=0.9, linestyle="--")
    ax.axhline(+3.0, linewidth=0.8, linestyle=":")
    ax.axhline(-3.0, linewidth=0.8, linestyle=":")
    ax.set_ylim(-5.0, 5.0)
    n_pull_out = int(np.sum(np.isfinite(z) & ((z < -5.0) | (z > 5.0))))
    if n_pull_out > 0:
        ax.text(
            0.98, 0.96, f"{n_pull_out} points outside displayed y range",
            transform=ax.transAxes, ha="right", va="top", fontsize=8,
        )
    #endif
    ax.set_xlabel(r"Jo matched $|t|$ (GeV$^2$)")
    ax.set_ylabel(
        r"$(\sigma_{\rm S}^{\rm tr}-\sigma_{\rm J})/"
        r"\sqrt{\delta_{\rm S,tr}^2+\delta_{\rm J}^2}$"
    )
    ax.set_title("Normalized point-pair difference")
    ax.grid(alpha=0.2)

    ax = axes[1, 1]
    goodr = finite_ratio
    if len(goodr):
        lo = max(0.0, float(np.percentile(goodr, 1.0)))
        hi = float(np.percentile(goodr, 99.0))
        if not np.isfinite(hi) or hi <= lo:
            hi = float(np.nanmax(goodr))
        #endif
        shown = goodr[(goodr >= lo) & (goodr <= hi)]
        ax.hist(shown, bins=35, histtype="step", linewidth=1.5)
        ax.axvline(1.0, linewidth=1.0, linestyle="--", label="unity")
        ax.axvline(rmed, linewidth=1.3, label=f"median = {rmed:.3f}")
        ax.axvspan(rp16, rp84, alpha=0.12, label="16--84% interval")
        ax.set_xlim(lo, hi)
    #endif
    ax.set_xlabel(r"transported $\sigma_{\rm Saylor}/\sigma_{\rm Jo}$")
    ax.set_ylabel("Matched pairs")
    ax.set_title("Ratio distribution (central 98% shown)")
    ax.grid(alpha=0.2)
    ax.legend(fontsize=8)

    fig.suptitle(
        "Jo 2015 vs Saylor 2018: direct best-overlap matched-kinematics comparison",
        y=0.992,
    )
    subtitle = (
        f"best-overlap 50% selected only by kinematic distance: N={len(tight)}; "
        f"median ratio={rmed:.3f}, 16--84%=[{rp16:.3f}, {rp84:.3f}]"
    )
    fig.text(0.5, 0.953, subtitle, ha="center", va="top", fontsize=10)
    fig.text(
        0.5, 0.928,
        "Saylor cross sections are locally transported to the matched Jo "
        "kinematics using only the KM15 point-to-point kinematic dependence.",
        ha="center", va="top", fontsize=9,
    )
    fig.tight_layout(rect=(0.02, 0.02, 0.99, 0.90))
    fig.savefig(
        outdir / "02_jo_saylor_matched_kinematics.png", dpi=220
    )
    plt.close(fig)

    # ------------------------------------------------------------------
    # Figure 03: make explicit what the local transport changes.
    # ------------------------------------------------------------------
    raw_all = tight["raw_ratio_saylor_over_jo"].to_numpy(float)
    tr_all = tight["ratio_saylor_over_jo"].to_numpy(float)
    finite_both = np.isfinite(raw_all) & np.isfinite(tr_all)
    raw_all = raw_all[finite_both]
    tr_all = tr_all[finite_both]
    dist_all = tight["match_distance"].to_numpy(float)[finite_both]

    combined_ratio = np.concatenate([raw_all, tr_all]) if len(raw_all) else np.array([])
    if len(combined_ratio):
        ylo = max(0.0, float(np.percentile(combined_ratio, 1.0)))
        yhi = float(np.percentile(combined_ratio, 99.0))
    else:
        ylo, yhi = 0.0, 2.0
    #endif

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.1), sharey=True)
    axes[0].scatter(dist_all, raw_all, s=11, alpha=0.25)
    axes[0].axhline(1.0, linewidth=1.0, linestyle="--")
    axes[0].set_xlabel("standardized nearest-Jo distance")
    axes[0].set_ylabel(r"$\sigma_{\rm Saylor}/\sigma_{\rm Jo}$")
    axes[0].set_title(
        f"Raw cross-section ratio\nmedian={np.nanmedian(raw_all):.3f}"
        if len(raw_all) else "Raw cross-section ratio"
    )
    axes[0].grid(alpha=0.2)

    axes[1].scatter(dist_all, tr_all, s=11, alpha=0.25)
    axes[1].axhline(1.0, linewidth=1.0, linestyle="--")
    axes[1].set_xlabel("standardized nearest-Jo distance")
    axes[1].set_title(
        f"After local KM15 transport\nmedian={np.nanmedian(tr_all):.3f}"
        if len(tr_all) else "After local KM15 transport"
    )
    axes[1].grid(alpha=0.2)
    # Fixed display range chosen to resolve the bulk ratio distribution.
    # The numerical analysis and CSV retain every point.
    axes[0].set_ylim(0.0, 2.5)
    n_raw_out = int(np.sum(np.isfinite(raw_all) & ((raw_all < 0.0) | (raw_all > 2.5))))
    n_tr_out = int(np.sum(np.isfinite(tr_all) & ((tr_all < 0.0) | (tr_all > 2.5))))
    if n_raw_out > 0:
        axes[0].text(
            0.98, 0.96, f"{n_raw_out} outside y range",
            transform=axes[0].transAxes, ha="right", va="top", fontsize=8,
        )
    #endif
    if n_tr_out > 0:
        axes[1].text(
            0.98, 0.96, f"{n_tr_out} outside y range",
            transform=axes[1].transAxes, ha="right", va="top", fontsize=8,
        )
    #endif
    fig.suptitle(
        "Jo 2015 vs Saylor 2018: raw and locally kinematic-transported ratios",
        y=0.985,
    )
    fig.tight_layout(rect=(0.02, 0.02, 0.99, 0.92))
    fig.savefig(
        outdir / "03_jo_saylor_raw_vs_transport_ratio.png", dpi=220
    )
    plt.close(fig)

    # ------------------------------------------------------------------
    # Figure 04: matching quality itself, so the reader can judge overlap.
    # ------------------------------------------------------------------
    quality_vars = [
        ("delta_Q2_abs", r"$|\Delta Q^2|$ (GeV$^2$)"),
        ("delta_xB_abs", r"$|\Delta x_B|$"),
        ("delta_t_abs", r"$|\Delta |t||$ (GeV$^2$)"),
        ("delta_phi_deg_abs", r"$|\Delta\phi|$ (deg)"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.4))
    for ax, (col, xlabel) in zip(axes.flat, quality_vars):
        vals = tight[col].to_numpy(float)
        vals = vals[np.isfinite(vals)]
        ax.hist(vals, bins=30, histtype="step", linewidth=1.5)
        if len(vals):
            med = float(np.median(vals))
            p84 = float(np.percentile(vals, 84.0))
            ax.axvline(med, linewidth=1.1, linestyle="--")
            ax.text(
                0.97, 0.94,
                f"median={med:.4g}\n84%<{p84:.4g}",
                transform=ax.transAxes, ha="right", va="top", fontsize=9,
            )
        #endif
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Matched pairs")
        ax.grid(alpha=0.2)
    #endfor
    fig.suptitle(
        "Jo 2015 vs Saylor 2018: kinematic quality of the best-overlap 50%",
        y=0.985,
    )
    fig.tight_layout(rect=(0.02, 0.02, 0.99, 0.94))
    fig.savefig(
        outdir / "04_jo_saylor_matching_quality.png", dpi=220
    )
    plt.close(fig)

    # Compact presentation/note table for the best-overlap half.
    zgood = z[np.isfinite(z)]
    note_summary = pd.DataFrame([{
        "subset": "best_overlap_50pct",
        "N": int(len(tight)),
        "max_standardized_match_distance": q50_cut,
        "median_transported_ratio": float(rmed),
        "ratio_p16": float(rp16),
        "ratio_p84": float(rp84),
        "rms_ratio_minus_one": float(
            np.sqrt(np.nanmean((finite_ratio - 1.0)**2))
        ) if len(finite_ratio) else np.nan,
        "fraction_ratio_within_10pct": float(
            np.mean(np.abs(finite_ratio - 1.0) <= 0.10)
        ) if len(finite_ratio) else np.nan,
        "fraction_ratio_within_20pct": float(
            np.mean(np.abs(finite_ratio - 1.0) <= 0.20)
        ) if len(finite_ratio) else np.nan,
        "mean_normalized_difference": float(np.mean(zgood)) if len(zgood) else np.nan,
        "rms_normalized_difference": float(
            np.sqrt(np.mean(zgood**2))
        ) if len(zgood) else np.nan,
        "fraction_abs_normalized_difference_gt_2": float(
            np.mean(np.abs(zgood) > 2.0)
        ) if len(zgood) else np.nan,
        "fraction_abs_normalized_difference_gt_3": float(
            np.mean(np.abs(zgood) > 3.0)
        ) if len(zgood) else np.nan,
    }])
    note_summary.to_csv(
        outdir / "04_jo_saylor_best_overlap_note_summary.csv", index=False
    )

    print(
        "[Jo/Saylor matched diagnostic] "
        f"best-overlap 50%: N={len(tight)}, median transported ratio="
        f"{rmed:.4f}, 16--84%=[{rp16:.4f},{rp84:.4f}]"
    )
#enddef


def save_five_dataset_bh_selected_consistency(
        bundles: Sequence[Dict[str, object]],
        outdir: Path,
        threshold: float = 0.05) -> None:
    """
    Five-dataset context table/figure for the nominal KM15 BH-like region.

    This diagnostic answers a different question from the full-EP model audit:
    once the same |1-BH/EP| <= threshold condition used by the radius analysis
    is imposed, how broad are the residuals of each published dataset relative
    to KM15 full electroproduction?

    Overall normalization is profiled using each dataset's published correlated
    normalization prior when one is available.  The table reports both raw and
    profiled scores.  This remains diagnostic and is never fed back into the
    production selection or functional-form ranking.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    pull_store = {}

    for bundle in sort_bundles_chronologically(list(bundles)):
        data = bundle["all_data"].copy().reset_index(drop=True)
        kind = str(bundle["kind"])
        key = str(bundle["key"])
        label = str(bundle["label"])

        needed = ["xs", "km15_bh", "km15_ep"]
        if any(c not in data.columns for c in needed):
            print(
                f"[five-dataset BH consistency] {label}: missing KM15 columns; skipped"
            )
            continue
        #endif

        y = pd.to_numeric(data["xs"], errors="coerce").to_numpy(float)
        bh = pd.to_numeric(data["km15_bh"], errors="coerce").to_numpy(float)
        ep = pd.to_numeric(data["km15_ep"], errors="coerce").to_numpy(float)
        err = dataset_point_errors(data, kind, 0.0, False)

        finite = (
            np.isfinite(y) & (y > 0.0)
            & np.isfinite(bh) & (bh > 0.0)
            & np.isfinite(ep) & (ep > 0.0)
            & np.isfinite(err) & (err > 0.0)
        )
        purity = np.full(len(data), np.nan)
        purity[finite] = np.abs(1.0 - bh[finite] / ep[finite])
        selected = finite & (purity <= float(threshold))
        if not np.any(selected):
            continue
        #endif

        yy = y[selected]
        mm = ep[selected]
        ee = err[selected]
        raw_pull = (mm - yy) / ee
        raw_chi2 = float(np.sum(raw_pull**2))

        # Dataset normalization priors are diagnostic here.  A bundle with
        # norm_frac <= 0 remains fixed to unity.
        norm_frac = float(bundle.get("norm_frac", 0.0))
        beta, scale, prof_chi2 = _profile_model_normalization(
            yy, mm, ee, norm_frac
        )
        prof_pull = (scale * mm - yy) / ee
        ratio = yy / mm

        # Unbinned Gaussian maximum-likelihood estimates.  For a normal
        # distribution these are simply the sample mean and the MLE width
        # (ddof=0).  Fit the FULL finite pull sample; do not let the plotting
        # range or histogram binning bias mu or sigma.
        prof_pull_finite = prof_pull[np.isfinite(prof_pull)]
        gaussian_mu = (
            float(np.mean(prof_pull_finite))
            if len(prof_pull_finite) else np.nan
        )
        gaussian_sigma = (
            float(np.std(prof_pull_finite, ddof=0))
            if len(prof_pull_finite) else np.nan
        )

        pull_store[label] = prof_pull.copy()
        rows.append({
            "dataset": key,
            "dataset_label": label,
            "N_5pct": int(np.sum(selected)),
            "normalization_prior_fraction": norm_frac,
            "profiled_beta": float(beta),
            "profiled_scale": float(scale),
            "raw_chi2_per_point": float(raw_chi2 / np.sum(selected)),
            "profiled_chi2_per_point": float(prof_chi2 / np.sum(selected)),
            "mean_profiled_pull": float(np.mean(prof_pull)),
            "rms_profiled_pull": float(np.sqrt(np.mean(prof_pull**2))),
            "gaussian_fit_mu": gaussian_mu,
            "gaussian_fit_sigma": gaussian_sigma,
            "median_data_over_km15_ep": float(np.median(ratio)),
            "fraction_abs_profiled_pull_gt_2": float(
                np.mean(np.abs(prof_pull) > 2.0)
            ),
            "fraction_abs_profiled_pull_gt_3": float(
                np.mean(np.abs(prof_pull) > 3.0)
            ),
        })
    #endfor

    table = pd.DataFrame(rows)
    table.to_csv(
        outdir / "05_five_dataset_km15_5pct_consistency_summary.csv",
        index=False,
    )
    if table.empty:
        return
    #endif

    # Human-readable LaTeX fragment for direct inclusion in the note.
    tex_cols = [
        "dataset_label", "N_5pct", "profiled_chi2_per_point",
        "rms_profiled_pull", "fraction_abs_profiled_pull_gt_2",
        "fraction_abs_profiled_pull_gt_3",
    ]
    latex_table = table[tex_cols].copy()
    latex_table.columns = [
        "Dataset", "$N_{5\\%}$", "$\\chi^2/N$",
        "pull RMS", "$f(|p|>2)$", "$f(|p|>3)$",
    ]
    latex_table["$\\chi^2/N$"] = latex_table["$\\chi^2/N$"].map(
        lambda x: f"{x:.2f}"
    )
    latex_table["pull RMS"] = latex_table["pull RMS"].map(
        lambda x: f"{x:.2f}"
    )
    latex_table["$f(|p|>2)$"] = latex_table["$f(|p|>2)$"].map(
        lambda x: f"{100.0*x:.1f}\\%"
    )
    latex_table["$f(|p|>3)$"] = latex_table["$f(|p|>3)$"].map(
        lambda x: f"{100.0*x:.1f}\\%"
    )
    with open(
        outdir / "05_five_dataset_km15_5pct_consistency_table.tex",
        "w",
        encoding="utf-8",
    ) as fout:
        fout.write(
            latex_table.to_latex(
                index=False,
                escape=False,
                caption=(
                    "Consistency of the five proton cross-section products "
                    "with KM15 full electroproduction within the common "
                    "$5\\%$ BH-like selection. Overall normalization is "
                    "profiled using each dataset's published correlated "
                    "normalization prior where available."
                ),
                label="tab:five_dataset_km15_bh_consistency",
                float_format="%.3f",
            )
        )
    #endwith

    # One compact context figure: profiled pull distributions with identical
    # binning, so Saylor can be judged against all other datasets immediately.
    finite_pulls = np.concatenate([
        p[np.isfinite(p)] for p in pull_store.values() if np.any(np.isfinite(p))
    ])
    if len(finite_pulls):
        lim = max(4.0, float(np.percentile(np.abs(finite_pulls), 99.0)))
        lim = min(lim, 8.0)
    else:
        lim = 5.0
    #endif
    bins = np.linspace(-lim, lim, 51)

    fig, ax = plt.subplots(figsize=(10.5, 6.4))
    xx = np.linspace(-lim, lim, 600)

    for label, pull in pull_store.items():
        pp_full = pull[np.isfinite(pull)]
        pp_plot = pp_full[(pp_full >= -lim) & (pp_full <= lim)]

        # Draw the histogram first and retain its automatically assigned color
        # so the corresponding Gaussian fit uses exactly the same color.
        hist_out = ax.hist(
            pp_plot, bins=bins, histtype="step", linewidth=1.4,
            density=True,
            label="_nolegend_",
        )
        hist_color = hist_out[2][0].get_edgecolor()

        if len(pp_full):
            mu_fit = float(np.mean(pp_full))
            sigma_fit = float(np.std(pp_full, ddof=0))
        else:
            mu_fit = np.nan
            sigma_fit = np.nan
        #endif

        # The fit is unbinned and uses all finite pulls, including any tails
        # outside the displayed x range.  Only the visualization is clipped.
        if np.isfinite(mu_fit) and np.isfinite(sigma_fit) and sigma_fit > 0.0:
            yy_fit = (
                np.exp(-0.5 * ((xx - mu_fit) / sigma_fit)**2)
                / (sigma_fit * np.sqrt(2.0 * np.pi))
            )
            ax.plot(
                xx, yy_fit,
                color=hist_color,
                linewidth=1.15,
                alpha=0.85,
                label=(
                    f"{label}: "
                    rf"$\mu={mu_fit:+.2f},\ \sigma={sigma_fit:.2f}$"
                ),
            )
        else:
            ax.plot(
                [], [], color=hist_color, linewidth=1.4,
                label=f"{label}: Gaussian fit unavailable",
            )
        #endif
    #endfor

    ax.plot(
        xx, np.exp(-0.5 * xx**2) / np.sqrt(2.0 * np.pi),
        linestyle="--", linewidth=1.2, label=r"unit Gaussian: $\mu=0,\ \sigma=1$",
    )
    ax.axvline(0.0, linewidth=0.9)
    ax.set_xlabel("profiled KM15 residual pull")
    ax.set_ylabel("normalized density")
    ax.set_title(
        "Five proton datasets in the common KM15 5% BH-like region"
    )
    ax.grid(alpha=0.2)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(
        outdir / "05_five_dataset_km15_5pct_pull_distributions.png", dpi=220
    )
    plt.close(fig)

    # Ranked summary plot makes the relative anomaly visually immediate.
    plot_table = table.sort_values("rms_profiled_pull", ascending=True).copy()
    ypos = np.arange(len(plot_table))
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.4))
    axes[0].barh(ypos, plot_table["rms_profiled_pull"])
    axes[0].axvline(1.0, linewidth=1.0, linestyle="--")
    axes[0].set_yticks(ypos)
    axes[0].set_yticklabels(plot_table["dataset_label"])
    axes[0].set_xlabel("profiled pull RMS")
    axes[0].grid(axis="x", alpha=0.2)

    axes[1].barh(
        ypos, 100.0 * plot_table["fraction_abs_profiled_pull_gt_2"]
    )
    axes[1].set_yticks(ypos)
    axes[1].set_yticklabels([])
    axes[1].set_xlabel(r"fraction with $|p|>2$ (\%)")
    axes[1].grid(axis="x", alpha=0.2)

    fig.suptitle(
        "Common 5% BH-like selection: cross-dataset KM15 residual comparison",
        y=0.985,
    )
    fig.tight_layout(rect=(0.02, 0.02, 0.99, 0.93))
    fig.savefig(
        outdir / "06_five_dataset_km15_5pct_residual_summary.png", dpi=220
    )
    plt.close(fig)

    print(
        "[five-dataset BH consistency] outputs -> "
        f"{outdir}"
    )
#enddef


def run_published_default(args) -> int:
    """
    Publication-facing workflow.

    A normal --run-final-model-analysis invocation now assembles all material
    needed for the current study automatically:
      * Jo 2015, Defurne 2015, Georges 2022, Lee 2026 production candidates;
      * Saylor 2018 diagnostic sample;
      * five-dataset all-point model-agreement diagnostics;
      * KM15-only all-four closure/radius-bias family selection;
      * Lee-only, Jo+Lee, and all-four KM15 radius/threshold results;
      * Saylor recoverability scans.

    The replica count and worker count remain controlled only by the existing
    --radius-bias-replicas and --radius-bias-workers options.
    """
    print("\n" + "=" * 78)
    print("[DEFAULT MODE] Unified published BH form-factor/radius study")
    print("[production candidates] Jo 2015; Defurne 2015; Defurne 2017; Georges 2022; Lee 2026")
    print("[BH-purity prescription] KM15")
    print("[Saylor 2018] diagnostic/recoverability study")
    print("=" * 78)

    root_outdir = Path(args.outdir).expanduser().resolve()

    jo = run_clas6_validation(args, return_results=True)
    lee = run_pass1_validation(args, return_results=True)
    defurne = run_halla_defurne_validation(args, return_results=True)
    defurne2017 = run_halla_defurne2017_validation(args, return_results=True)
    georges = make_georges_diagnostic_bundle(args)

    # Saylor is always loaded for the unified diagnostic study.  No extra
    # --include-saylor switch is required in the current final workflow.
    saylor_path = resolve_script_relative_path(args.saylor_file)
    if args.download_saylor:
        saylor_path = maybe_download_saylor_supplement(
            saylor_path, force=args.force_saylor_download
        )
        args.saylor_file = str(saylor_path)
    #endif
    if not saylor_path.exists():
        raise FileNotFoundError(
            "Unified final analysis requires the Saylor supplemental table "
            "for its diagnostic/recoverability study.\n"
            f"Resolved path: {saylor_path}"
        )
    #endif
    args.saylor_file = str(saylor_path)
    saylor = run_saylor_validation(args, return_results=True)

    production_bundles = sort_bundles_chronologically(
        [jo, defurne, defurne2017, georges, lee]
    )
    diagnostic_bundles = sort_bundles_chronologically(
        [jo, defurne, defurne2017, saylor, georges, lee]
    )

    # Export the exact five-dataset table expected by the external evaluator.
    diagnostic_input_root = (
        root_outdir / "final_analysis" / "diagnostics"
        / "all_point_model_agreement_input"
    )
    export_bh_model_selection_kinematics(
        diagnostic_bundles, diagnostic_input_root
    )

    # Also export the four production candidates to the standard selection
    # input location.  The external result table may remain a superset.
    export_bh_model_selection_kinematics(production_bundles, root_outdir)

    jo_saylor_dir = (
        root_outdir / "final_analysis" / "diagnostics"
        / "jo2015_vs_saylor2018"
    )
    save_jo_saylor_direct_comparison(jo, saylor, jo_saylor_dir)
    save_jo_saylor_matched_kinematics(jo, saylor, jo_saylor_dir)

    if args.run_final_model_analysis:
        selection = load_external_bh_model_selection(
            Path(args.bh_model_selection_results).expanduser().resolve()
        )

        # Full five-dataset electroproduction model-agreement audit remains
        # available even though KM15 alone defines the BH-purity selection.
        save_all_point_model_agreement_diagnostics(
            selection,
            diagnostic_bundles,
            root_outdir / "final_analysis" / "diagnostics"
            / "all_point_model_agreement",
        )

        # Dedicated context for the Saylor decision: compare all five proton
        # products under the SAME nominal KM15 5% BH-like criterion.
        save_five_dataset_bh_selected_consistency(
            diagnostic_bundles,
            root_outdir / "final_analysis" / "diagnostics"
            / "jo2015_vs_saylor2018",
            threshold=0.05,
        )

        run_unified_km15_final_analysis(
            jo,
            defurne,
            defurne2017,
            georges,
            lee,
            saylor,
            args,
            root_outdir,
        )
    #endif

    print(f"\nDone. Unified results are in {root_outdir}")
    return 0
#enddef




def run_all_three(args) -> int:
    print("\\n" + "=" * 78)
    print("[DEFAULT MODE] Running pass 2, pass 1, and CLAS6")
    print("=" * 78)

    pass2 = run_pass2_analysis(args, return_results=True)
    pass1 = run_pass1_validation(args, return_results=True)
    clas6 = run_clas6_validation(args, return_results=True)

    comparison_dir = Path(args.outdir).expanduser().resolve() / "comparisons"
    comparison_dir.mkdir(parents=True, exist_ok=True)

    bundles = [clas6, pass1, pass2]
    save_cross_dataset_form_factor_comparison(
        bundles, comparison_dir, "Fit 5"
    )
    save_cross_dataset_form_factor_comparison(
        bundles, comparison_dir, "Fit 8"
    )
    save_cross_dataset_radius_comparison(bundles, comparison_dir)
    save_cross_dataset_summary(bundles, comparison_dir)

    print(f"\\n[comparison] Results written to {comparison_dir}")
    return 0
#enddef


def run_pass1_validation(args, return_results: bool = False):
    """
    Run the same BH-dominance/form-factor diagnostic on the legacy CLAS12
    pass-1 cross sections.

    This mode bypasses the pass-2 CSV and uses the authoritative CLAS Physics Database E214M1 release by default:
      * released statistical uncertainty
      * released total systematic uncertainty, decomposed into pointwise and 31% correlated normalization pieces
      * paper's 1--5% BH-selection uncertainty
      * one correlated 31% overall normalization nuisance
      * no pass-2 combination-systematic nuisance
    """
    csv_path = Path(args.pass1_csv).expanduser().resolve()
    outdir = (
        Path(args.outdir).expanduser().resolve()
        / "clas12_lee2026"
    )
    outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 78)
    print("[CLAS12 LEE 2026 VALIDATION MODE]")
    print(f"[input]  {csv_path}")
    print(f"[output] {outdir}")
    print("=" * 78)

    df = load_clas12_pass1_csv(csv_path)

    cache = outdir / "km15_bh_decomposition_pass1.csv"
    df = evaluate_km15_dataframe(
        df,
        10.604,
        max(1, args.workers),
        cache,
        args.force_km15,
    )

    finite = (
        np.isfinite(df["km15_ep"])
        & np.isfinite(df["km15_bh"])
        & np.isfinite(df["R_BH"])
        & (df["km15_ep"] > 0.0)
        & (df["km15_bh"] > 0.0)
    )
    df = df.loc[finite].copy().reset_index(drop=True)

    max_ep_relerr = float(np.nanmax(df["km15_ep_decomp_relerr"]))
    max_bh_relerr = float(np.nanmax(df["bh_quad_relerr"]))
    print(
        f"[check] max KM15 EP decomposition rel. error = {max_ep_relerr:.3e}"
    )
    print(
        f"[check] max BH quadratic reconstruction rel. error = {max_bh_relerr:.3e}"
    )

    if max_ep_relerr > 1.0e-7 or max_bh_relerr > 1.0e-7:
        raise RuntimeError(
            "Pass-1 validation failed internal KM15 decomposition checks."
        )
    #endif

    bh_cuts = [0.01, 0.02, 0.03, 0.04, 0.05]
    selected_sets = {}
    selection_rows = []

    print("\n[PASS1 BH selections]")
    for cut in bh_cuts:
        selected = df.loc[df["bh_delta"] <= cut].copy()
        selected_sets[cut] = selected
        print(
            f"  |1-R_BH| <= {100*cut:.0f}% : {len(selected)} points"
        )
        selection_rows.append({
            "BH_cut": cut,
            "Npts": len(selected),
            "xB_min": selected["xB"].min() if len(selected) else np.nan,
            "xB_max": selected["xB"].max() if len(selected) else np.nan,
            "Q2_min": selected["Q2"].min() if len(selected) else np.nan,
            "Q2_max": selected["Q2"].max() if len(selected) else np.nan,
            "t_min": selected["t_abs"].min() if len(selected) else np.nan,
            "t_max": selected["t_abs"].max() if len(selected) else np.nan,
        })
    #endfor

    pd.DataFrame(selection_rows).to_csv(
        outdir / "bh_selection_summary.csv",
        index=False,
    )

    fit_results = []

    for fit_index, cut in enumerate(bh_cuts, start=1):
        data = selected_sets[cut]
        fr = fit_paper_model(
            data=data,
            kind="dipole",
            fit_name=f"Fit {fit_index}",
            bh_cut=cut,
            include_clas12_ptp_sys=True,
            include_bh_sys=True,
            include_combination_norm_sys=False,
            include_global_norm_sys=True,
            global_scale_frac=PASS1_GLOBAL_SCALE_FRAC,
        )
        fit_results.append(fr)

        print(
            f"[PASS1 fit] Fit {fit_index}: "
            f"N={fr.npts:4d} "
            f"chi2/dof={fr.chi2_ndof:.3f} "
            f"aE={fr.params[0]:.4f} "
            f"aM={fr.params[1]:.4f} "
            f"beta_global={fr.meta.get('beta_global', np.nan):+.3f}"
        )
    #endfor

    set5 = selected_sets[0.05].copy()
    set5["fit_sigma_default"] = fit_sigma_errors(
        set5,
        0.05,
        include_clas12_ptp_sys=True,
        include_bh_sys=True,
    )
    set5.to_csv(outdir / "set5_selected_points.csv", index=False)

    for kind, fit_name in [
        ("fit6_same_a", "Fit 6"),
        ("fit7_same_p", "Fit 7"),
        ("fit8_f2_kelly", "Fit 8"),
    ]:
        fr = fit_paper_model(
            data=set5,
            kind=kind,
            fit_name=fit_name,
            bh_cut=0.05,
            include_clas12_ptp_sys=True,
            include_bh_sys=True,
            include_combination_norm_sys=False,
            include_global_norm_sys=True,
            global_scale_frac=PASS1_GLOBAL_SCALE_FRAC,
        )
        fit_results.append(fr)
        print(
            f"[PASS1 fit] {fit_name}: "
            f"chi2/dof={fr.chi2_ndof:.3f} "
            f"beta_global={fr.meta.get('beta_global', np.nan):+.3f}"
        )
    #endfor

    pd.DataFrame(
        [fitresult_to_record(fr) for fr in fit_results]
    ).to_csv(outdir / "fit_results.csv", index=False)

    fit5 = next(r for r in fit_results if r.name == "Fit 5")
    fit8 = next(r for r in fit_results if r.name == "Fit 8")

    save_fit1_to_fit5_plots(fit_results, set5, outdir)
    save_fit5_to_fit8_sachs_plot(fit_results, set5, outdir)
    save_radii_plot(fit_results, outdir)
    save_low_q2_ratio_plots(fit_results, set5, outdir)
    save_elastic_reference_comparison(fit_results, set5, outdir, "Fit 5")
    save_elastic_reference_comparison(fit_results, set5, outdir, "Fit 8")
    save_bh_local_f1_f2_sensitivity(outdir, set5, fit5)
    save_pass1_cross_section_diagnostics(set5, fit5, fit8, outdir)

    print("\n[PASS1 radii]")
    for fit_name in ["Fit 5", "Fit 6", "Fit 7", "Fit 8"]:
        fr = next(r for r in fit_results if r.name == fit_name)
        print(
            f"  {fit_name}: "
            f"rE={fr.rE_fm:.5f}+/-{fr.rE_err_fm:.5f} fm, "
            f"rM={fr.rM_fm:.5f}+/-{fr.rM_err_fm:.5f} fm, "
            f"beta_global={fr.meta.get('beta_global', np.nan):+.3f}, "
            f"chi2/dof={fr.chi2_ndof:.3f}"
        )
    #endfor

    print(f"\nDone. Pass-1 validation results are in {outdir}")
    if return_results:
        return {
            "label": "CLAS12 Lee 2026",
            "key": "pass1",
            "kind": "pass1",
            "norm_frac": PASS1_GLOBAL_SCALE_FRAC,
            "results": fit_results,
            "set5": set5,
            "all_data": df.copy(),
            "outdir": outdir,
        }
    #endif
    return 0
#enddef



def run_saylor_validation(args, return_results: bool = False):
    """Standalone Saylor 2018 BH/form-factor analysis."""
    outdir = (
        Path(args.outdir).expanduser().resolve()
        / "clas6_saylor2018"
    )
    outdir.mkdir(parents=True, exist_ok=True)

    saylor_path = resolve_script_relative_path(args.saylor_file)
    if args.download_saylor:
        saylor_path = maybe_download_saylor_supplement(
            saylor_path, force=args.force_saylor_download
        )
    #endif

    print("=" * 78)
    print("[CLAS6 SAYLOR 2018 MODE]")
    print(f"[input]  {saylor_path.resolve()}")
    print(f"[output] {outdir}")
    print("=" * 78)

    df = load_saylor_supplement(str(saylor_path))
    cache = outdir / "km15_bh_decomposition_saylor2018.csv"
    df = evaluate_km15_saylor_dataframe(
        df,
        workers=args.workers,
        cache_path=cache,
        force=args.force_km15,
    )

    finite = (
        np.isfinite(df["R_BH"])
        & np.isfinite(df["bh_A"])
        & np.isfinite(df["bh_B"])
        & np.isfinite(df["bh_C"])
    )
    df = df.loc[finite].reset_index(drop=True)

    fit_results = []
    sets = {}
    for i, cut in enumerate(args.bh_cuts, start=1):
        selected = df.loc[df["bh_delta"] <= cut].copy()
        sets[float(cut)] = selected
        if len(selected) < 5:
            continue
        #endif
        fr = fit_multi_measurements(
            [measurement_spec(
                "saylor2018", "CLAS6 Saylor 2018", "saylor2018",
                selected, SAYLOR_GLOBAL_SCALE_FRAC
            )],
            kind="dipole",
            fit_name=f"Fit {i}",
            bh_cut=float(cut),
            add_moradi_bh_systematic=True,
        )
        fit_results.append(fr)
        print(
            f"[SAYLOR fit] Fit {i}: N={fr.npts:4d} "
            f"chi2/dof={fr.chi2_ndof:.3f}"
        )
    #endfor

    set5 = sets.get(0.05)
    if set5 is None or len(set5) < 5:
        raise RuntimeError("Saylor 5% BH-selected sample has too few points.")
    #endif

    for kind, name in [
        ("fit6_same_a", "Fit 6"),
        ("fit7_same_p", "Fit 7"),
        ("fit8_f2_kelly", "Fit 8"),
    ]:
        fit_results.append(
            fit_multi_measurements(
                [measurement_spec(
                    "saylor2018", "CLAS6 Saylor 2018", "saylor2018",
                    set5, SAYLOR_GLOBAL_SCALE_FRAC
                )],
                kind=kind,
                fit_name=name,
                bh_cut=0.05,
                add_moradi_bh_systematic=True,
            )
        )
    #endfor

    pd.DataFrame(
        [fitresult_to_record(fr) for fr in fit_results]
    ).to_csv(outdir / "fit_results.csv", index=False)

    save_fit1_to_fit5_plots(fit_results, set5, outdir)
    save_fit5_to_fit8_sachs_plot(fit_results, set5, outdir)
    save_radii_plot(fit_results, outdir)
    save_low_q2_ratio_plots(fit_results, set5, outdir)
    save_elastic_reference_comparison(fit_results, set5, outdir, "Fit 5")
    save_elastic_reference_comparison(fit_results, set5, outdir, "Fit 8")
    fit5 = next(r for r in fit_results if r.name == "Fit 5")
    save_bh_local_f1_f2_sensitivity(outdir, set5, fit5)
    save_fit5_residual_chi2_diagnostics(
        set5, fit5, "saylor2018", "saylor2018",
        outdir / "12_fit5_residual_diagnostics", bh_cut=0.05,
    )

    bundle = {
        "label": "CLAS6 Saylor 2018",
        "key": "saylor2018",
        "kind": "saylor2018",
        "norm_frac": SAYLOR_GLOBAL_SCALE_FRAC,
        "results": fit_results,
        "set5": set5,
        "all_data": df.copy(),
        "outdir": outdir,
    }
    if return_results:
        return bundle
    #endif
    return 0
#enddef




def run_halla_defurne_validation(args, return_results: bool = False):
    """Standalone Hall A Defurne 2015 BH/form-factor analysis."""
    outdir = Path(args.outdir).expanduser().resolve() / "halla_defurne2015"
    outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 78)
    print("[HALL A DEFURNE 2015 MODE]")
    print(f"[Gepard datasets] {HALLA_DEFURNE_GEPARD_DATASET_IDS}")
    print(f"[output] {outdir}")
    print("=" * 78)

    df = load_halla_defurne_gepard_datasets()
    df = evaluate_km15_halla_dataframe(
        df,
        workers=args.workers,
        cache_path=outdir / "km15_bh_decomposition_halla_defurne2015.csv",
        force=args.force_km15,
    )

    finite = (
        np.isfinite(df["R_BH"])
        & np.isfinite(df["bh_A"])
        & np.isfinite(df["bh_B"])
        & np.isfinite(df["bh_C"])
    )
    df = df.loc[finite].reset_index(drop=True)

    fit_results = []
    sets = {}
    print("\n[HALLA BH selections]")
    for i, cut in enumerate(args.bh_cuts, start=1):
        selected = df.loc[df["bh_delta"] <= cut].copy()
        sets[float(cut)] = selected
        print(f"  |1-R_BH| <= {100*cut:.0f}% : {len(selected)} points")
        if len(selected) < 5:
            continue
        #endif
        fr = fit_multi_measurements(
            [measurement_spec(
                "halla_defurne2015", "Hall A Defurne 2015",
                "halla_defurne2015", selected,
                HALLA_DEFURNE_GLOBAL_SCALE_FRAC,
            )],
            kind="dipole",
            fit_name=f"Fit {i}",
            bh_cut=float(cut),
            add_moradi_bh_systematic=True,
        )
        fit_results.append(fr)
        print(f"[HALLA fit] Fit {i}: N={fr.npts:4d} chi2/dof={fr.chi2_ndof:.3f}")
    #endfor

    set5 = sets.get(0.05)
    if set5 is None or len(set5) < 5:
        raise RuntimeError("Hall A Defurne 5% BH-selected sample has too few points.")
    #endif

    for kind, name in [
        ("fit6_same_a", "Fit 6"),
        ("fit7_same_p", "Fit 7"),
        ("fit8_f2_kelly", "Fit 8"),
    ]:
        fit_results.append(
            fit_multi_measurements(
                [measurement_spec(
                    "halla_defurne2015", "Hall A Defurne 2015",
                    "halla_defurne2015", set5,
                    HALLA_DEFURNE_GLOBAL_SCALE_FRAC,
                )],
                kind=kind,
                fit_name=name,
                bh_cut=0.05,
                add_moradi_bh_systematic=True,
            )
        )
    #endfor

    pd.DataFrame([fitresult_to_record(fr) for fr in fit_results]).to_csv(
        outdir / "fit_results.csv", index=False
    )
    set5.to_csv(outdir / "set5_selected_points.csv", index=False)

    save_fit1_to_fit5_plots(fit_results, set5, outdir)
    save_fit5_to_fit8_sachs_plot(fit_results, set5, outdir)
    save_combination_f1_f2_fits_5_to_8(fit_results, set5, outdir)
    save_radii_plot(fit_results, outdir)
    save_low_q2_ratio_plots(fit_results, set5, outdir)
    save_elastic_reference_comparison(fit_results, set5, outdir, "Fit 5")
    save_elastic_reference_comparison(fit_results, set5, outdir, "Fit 8")
    fit5 = next(r for r in fit_results if r.name == "Fit 5")
    save_bh_local_f1_f2_sensitivity(outdir, set5, fit5)
    save_fit5_residual_chi2_diagnostics(
        set5, fit5, "halla_defurne2015", "halla_defurne2015",
        outdir / "12_fit5_residual_diagnostics", bh_cut=0.05,
    )

    bundle = {
        "label": "Hall A Defurne 2015",
        "key": "halla_defurne2015",
        "kind": "halla_defurne2015",
        "norm_frac": HALLA_DEFURNE_GLOBAL_SCALE_FRAC,
        "results": fit_results,
        "set5": set5,
        "all_data": df.copy(),
        "outdir": outdir,
    }
    if return_results:
        return bundle
    #endif
    return 0
#enddef




def run_halla_defurne2017_validation(args, return_results: bool = False):
    """Standalone Hall A E07-007 / Defurne 2017 BH/form-factor analysis."""
    outdir = Path(args.outdir).expanduser().resolve() / "halla_defurne2017"
    outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 78)
    print("[HALL A DEFURNE 2017 / E07-007 MODE]")
    print(f"[Gepard datasets] {HALLA_DEFURNE2017_GEPARD_DATASET_IDS}")
    print(
        "[systematics] 3.2% point-to-point; "
        f"{100*HALLA_DEFURNE2017_GLOBAL_SCALE_FRAC:.3f}% correlated scale"
    )
    print(f"[output] {outdir}")
    print("=" * 78)

    df = load_halla_defurne2017_gepard_datasets()
    df = evaluate_km15_halla_dataframe(
        df,
        workers=args.workers,
        cache_path=outdir / "km15_bh_decomposition_halla_defurne2017.csv",
        force=args.force_km15,
    )

    finite = (
        np.isfinite(df["R_BH"])
        & np.isfinite(df["bh_A"])
        & np.isfinite(df["bh_B"])
        & np.isfinite(df["bh_C"])
    )
    df = df.loc[finite].reset_index(drop=True)

    fit_results = []
    sets = {}
    print("\n[HALLA 2017 BH selections]")
    for i, cut in enumerate(args.bh_cuts, start=1):
        selected = df.loc[df["bh_delta"] <= cut].copy()
        sets[float(cut)] = selected
        print(f"  |1-R_BH| <= {100*cut:.0f}% : {len(selected)} points")
        if len(selected) < 5:
            continue
        #endif
        fr = fit_multi_measurements(
            [measurement_spec(
                "halla_defurne2017", "Hall A Defurne 2017",
                "halla_defurne2017", selected,
                HALLA_DEFURNE2017_GLOBAL_SCALE_FRAC,
            )],
            kind="dipole",
            fit_name=f"Fit {i}",
            bh_cut=float(cut),
            add_moradi_bh_systematic=True,
        )
        fit_results.append(fr)
        print(
            f"[HALLA 2017 fit] Fit {i}: N={fr.npts:4d} "
            f"chi2/dof={fr.chi2_ndof:.3f}"
        )
    #endfor

    set5 = sets.get(0.05)
    if set5 is None or len(set5) < 5:
        raise RuntimeError(
            "Hall A Defurne 2017 5% BH-selected sample has too few points."
        )
    #endif

    for kind, name in [
        ("fit6_same_a", "Fit 6"),
        ("fit7_same_p", "Fit 7"),
        ("fit8_f2_kelly", "Fit 8"),
    ]:
        fit_results.append(
            fit_multi_measurements(
                [measurement_spec(
                    "halla_defurne2017", "Hall A Defurne 2017",
                    "halla_defurne2017", set5,
                    HALLA_DEFURNE2017_GLOBAL_SCALE_FRAC,
                )],
                kind=kind,
                fit_name=name,
                bh_cut=0.05,
                add_moradi_bh_systematic=True,
            )
        )
    #endfor

    pd.DataFrame([fitresult_to_record(fr) for fr in fit_results]).to_csv(
        outdir / "fit_results.csv", index=False
    )
    set5.to_csv(outdir / "set5_selected_points.csv", index=False)

    save_fit1_to_fit5_plots(fit_results, set5, outdir)
    save_fit5_to_fit8_sachs_plot(fit_results, set5, outdir)
    save_combination_f1_f2_fits_5_to_8(fit_results, set5, outdir)
    save_radii_plot(fit_results, outdir)
    save_low_q2_ratio_plots(fit_results, set5, outdir)
    fit5 = next(r for r in fit_results if r.name == "Fit 5")
    save_bh_local_f1_f2_sensitivity(outdir, set5, fit5)
    save_fit5_residual_chi2_diagnostics(
        set5, fit5, "halla_defurne2017", "halla_defurne2017",
        outdir / "12_fit5_residual_diagnostics", bh_cut=0.05,
    )

    bundle = {
        "label": "Hall A Defurne 2017",
        "key": "halla_defurne2017",
        "kind": "halla_defurne2017",
        "norm_frac": HALLA_DEFURNE2017_GLOBAL_SCALE_FRAC,
        "results": fit_results,
        "set5": set5,
        "all_data": df.copy(),
        "outdir": outdir,
    }
    if return_results:
        return bundle
    #endif
    return 0
#enddef


def run_clas6_validation(args, return_results: bool = False):
    """
    Reproduce the Moradi et al. CLAS6 analysis using Gepard's bundled
    CLAS:2015uuo XUU data.

    This is an override/validation mode:
      * no CLAS12 CSV is read
      * no CLAS12 point-to-point systematic is added
      * no CLAS12 correlated normalization nuisances are used
      * the published Gepard point.err is used as the experimental error
      * only the paper's additional 1--5% BH-selection uncertainty is added

    Expected paper benchmarks:
      selection counts = 98, 212, 292, 404, 473
      chi2/dof Fits 1--5 = 0.82, 0.75, 0.73, 0.68, 0.69
    """
    outdir = (
        Path(args.outdir).expanduser().resolve()
        / "clas6_jo2015"
    )
    outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 78)
    print("[CLAS6 VALIDATION MODE]")
    print(f"[output] {outdir}")
    print(
        f"[source] Gepard dataset {args.clas6_dataset_id} "
        "(CLAS 2015 XUU, Jo et al.)"
    )
    print("=" * 78)

    df = load_clas6_gepard_dataset(args.clas6_dataset_id)
    if len(df) != 2640 and args.clas6_dataset_id == 98:
        warnings.warn(
            f"Expected 2640 points in Gepard dataset 98, found {len(df)}."
        )
    #endif

    cache = outdir / "km15_bh_decomposition_clas6.csv"
    df = evaluate_km15_clas6_dataframe(
        df,
        max(1, args.workers),
        cache,
        args.force_km15,
    )

    finite = (
        np.isfinite(df["km15_ep"])
        & np.isfinite(df["km15_bh"])
        & np.isfinite(df["R_BH"])
        & (df["km15_ep"] > 0.0)
        & (df["km15_bh"] > 0.0)
    )
    df = df.loc[finite].copy().reset_index(drop=True)

    max_ep_relerr = float(np.nanmax(df["km15_ep_decomp_relerr"]))
    max_bh_relerr = float(np.nanmax(df["bh_quad_relerr"]))
    print(
        f"[check] max KM15 EP decomposition rel. error = {max_ep_relerr:.3e}"
    )
    print(
        f"[check] max BH quadratic reconstruction rel. error = {max_bh_relerr:.3e}"
    )

    if max_ep_relerr > 1.0e-7 or max_bh_relerr > 1.0e-7:
        raise RuntimeError(
            "CLAS6 validation failed the internal Gepard decomposition checks."
        )
    #endif

    bh_cuts = [0.01, 0.02, 0.03, 0.04, 0.05]
    selected_sets = {}
    selection_rows = []

    print("\n[CLAS6 selection benchmark]")
    for cut in bh_cuts:
        selected = df.loc[df["bh_delta"] <= cut].copy()
        selected_sets[cut] = selected

        expected = PAPER_SELECTION_COUNTS[cut]
        delta = len(selected) - expected
        status = "PASS" if delta == 0 else "MISMATCH"

        print(
            f"  {100*cut:.0f}%: got {len(selected):4d}, "
            f"paper {expected:4d}, delta={delta:+d}  [{status}]"
        )

        selection_rows.append({
            "BH_cut": cut,
            "our_Npts": len(selected),
            "paper_Npts": expected,
            "delta_Npts": delta,
            "match": delta == 0,
        })
    #endfor

    pd.DataFrame(selection_rows).to_csv(
        outdir / "selection_count_validation.csv",
        index=False,
    )

    # If these counts do not match, continuing the fit is still useful for
    # diagnosis, but loudly flag that we have not reproduced the paper's
    # selection/conventions yet.
    counts_match = all(
        len(selected_sets[c]) == PAPER_SELECTION_COUNTS[c]
        for c in bh_cuts
    )
    if not counts_match:
        warnings.warn(
            "CLAS6 BH-selected counts do not exactly reproduce Table I. "
            "Inspect phi/frame/observable conventions before interpreting "
            "the chi2 validation."
        )
    #endif

    fit_results = []

    # Fits 1--5 exactly as the paper: published total experimental error +
    # threshold-matched uncorrelated BH-selection systematic.
    for fit_index, cut in enumerate(bh_cuts, start=1):
        fr = fit_paper_model(
            data=selected_sets[cut],
            kind="dipole",
            fit_name=f"Fit {fit_index}",
            bh_cut=cut,
            include_clas12_ptp_sys=False,
            include_bh_sys=True,
            include_combination_norm_sys=False,
            include_global_norm_sys=False,
        )
        fit_results.append(fr)

        paper_chi = PAPER_CHI2_NDOF[fit_index]
        paper_aE, paper_daE, paper_aM, paper_daM = PAPER_DIPOLE_PARAMS[fit_index]

        print(
            f"[CLAS6 fit] Fit {fit_index}: "
            f"chi2/dof={fr.chi2_ndof:.3f} "
            f"(paper {paper_chi:.2f}); "
            f"aE={fr.params[0]:.4f} (paper {paper_aE:.3f}); "
            f"aM={fr.params[1]:.4f} (paper {paper_aM:.3f})"
        )
    #endfor

    set5 = selected_sets[0.05].copy()
    set5["fit_sigma_default"] = fit_sigma_errors(
        set5,
        0.05,
        include_clas12_ptp_sys=False,
        include_bh_sys=True,
    )
    set5.to_csv(outdir / "set5_selected_points.csv", index=False)

    for kind, fit_name, fit_index in [
        ("fit6_same_a", "Fit 6", 6),
        ("fit7_same_p", "Fit 7", 7),
        ("fit8_f2_kelly", "Fit 8", 8),
    ]:
        fr = fit_paper_model(
            data=set5,
            kind=kind,
            fit_name=fit_name,
            bh_cut=0.05,
            include_clas12_ptp_sys=False,
            include_bh_sys=True,
            include_combination_norm_sys=False,
            include_global_norm_sys=False,
        )
        fit_results.append(fr)
        print(
            f"[CLAS6 fit] {fit_name}: "
            f"chi2/dof={fr.chi2_ndof:.3f} "
            f"(paper {PAPER_CHI2_NDOF[fit_index]:.2f})"
        )
    #endfor

    # Direct comparison table to the published benchmarks.
    rows = []
    for fit_index in range(1, 9):
        fr = next(r for r in fit_results if r.name == f"Fit {fit_index}")
        row = {
            "Fit": fit_index,
            "our_chi2_ndof": fr.chi2_ndof,
            "paper_chi2_ndof": PAPER_CHI2_NDOF[fit_index],
            "delta_chi2_ndof": (
                fr.chi2_ndof - PAPER_CHI2_NDOF[fit_index]
            ),
            "rE_fm": fr.rE_fm,
            "rE_err_fm": fr.rE_err_fm,
            "rM_fm": fr.rM_fm,
            "rM_err_fm": fr.rM_err_fm,
        }
        if fit_index <= 5:
            row["our_aE"] = fr.params[0]
            row["our_aM"] = fr.params[1]
            row["paper_aE"] = PAPER_DIPOLE_PARAMS[fit_index][0]
            row["paper_aM"] = PAPER_DIPOLE_PARAMS[fit_index][2]
        #endif
        rows.append(row)
    #endfor
    comparison = pd.DataFrame(rows)
    comparison.to_csv(outdir / "paper_benchmark_comparison.csv", index=False)

    pd.DataFrame(
        [fitresult_to_record(fr) for fr in fit_results]
    ).to_csv(outdir / "fit_results.csv", index=False)

    # Paper-style consolidated plots.
    save_fit1_to_fit5_plots(fit_results, set5, outdir)
    fit5 = next(r for r in fit_results if r.name == "Fit 5")
    save_fit5_to_fit8_sachs_plot(fit_results, set5, outdir)
    save_radii_plot(fit_results, outdir)
    save_low_q2_ratio_plots(fit_results, set5, outdir)
    save_elastic_reference_comparison(fit_results, set5, outdir, "Fit 5")
    save_elastic_reference_comparison(fit_results, set5, outdir, "Fit 8")
    save_bh_local_f1_f2_sensitivity(outdir, set5, fit5)

    # Validation plot: our chi2/dof against the paper values.
    fig, ax = plt.subplots(figsize=(7.4, 5.0))
    fit_numbers = np.arange(1, 9)
    ours = [
        next(r for r in fit_results if r.name == f"Fit {i}").chi2_ndof
        for i in fit_numbers
    ]
    paper = [PAPER_CHI2_NDOF[i] for i in fit_numbers]
    ax.plot(fit_numbers, ours, "o-", label="This script / Gepard data")
    ax.plot(fit_numbers, paper, "s--", label="Moradi et al.")
    ax.set_xlabel("Fit")
    ax.set_ylabel(r"$\chi^2/\mathrm{dof}$")
    ax.set_xticks(fit_numbers)
    ax.grid(alpha=0.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "00_chi2_validation.png", dpi=180)
    plt.close(fig)

    print("\n[CLAS6 validation summary]")
    print(comparison[[
        "Fit", "our_chi2_ndof", "paper_chi2_ndof", "delta_chi2_ndof"
    ]].to_string(index=False))
    print(f"\nDone. CLAS6 validation results are in {outdir}")
    if return_results:
        return {
            "label": "CLAS6 Jo 2015",
            "key": "jo2015",
            "kind": "jo2015",
            "norm_frac": JO_GLOBAL_SCALE_FRAC,
            "results": fit_results,
            "set5": set5,
            "all_data": df.copy(),
            "outdir": outdir,
        }
    #endif
    return 0
#enddef


# -----------------------------------------------------------------------------
# CLI / main
# -----------------------------------------------------------------------------

def load_ndvcs_preliminary_unpolarized(path: Path) -> pd.DataFrame:
    """
    Load the 72 preliminary CLAS12 neutron-DVCS unpolarized points from
    Tables A.1-A.9 of the analysis note.

    The import file preserves the note's pb/GeV^4 values and also carries
    convenience nb/GeV^4 columns.  Internally this script uses nb/GeV^4, as in
    the existing Gepard-based proton workflow.
    """
    names = [
        "kin_bin", "phi_deg", "Q2_GeV2", "xB", "t_abs_GeV2",
        "xs_pb_GeV4", "stat_pb_GeV4", "sys_pb_GeV4",
        "xs_nb_GeV4", "stat_nb_GeV4", "sys_nb_GeV4",
    ]
    data = pd.read_csv(
        path,
        sep=r"\s+",
        comment="#",
        names=names,
        engine="python",
    )
    if len(data) != 72:
        raise RuntimeError(
            f"Expected 72 unpolarized neutron-DVCS points in {path}, "
            f"found {len(data)}."
        )
    #endif

    out = pd.DataFrame({
        "kin_bin": data["kin_bin"].astype(int),
        "phi": np.deg2rad(data["phi_deg"].to_numpy(float)),
        "phi_deg": data["phi_deg"].to_numpy(float),
        "Q2": data["Q2_GeV2"].to_numpy(float),
        "xB": data["xB"].to_numpy(float),
        "t_abs": data["t_abs_GeV2"].to_numpy(float),
        "xs": data["xs_nb_GeV4"].to_numpy(float),
        "stat": data["stat_nb_GeV4"].to_numpy(float),
        "sys": data["sys_nb_GeV4"].to_numpy(float),
    })
    out["t"] = -out["t_abs"]
    return out
#enddef



NEUTRON_PHI_CONVENTIONS = (
    "identity",
    "180-minus",
    "180-shift",
)


def transform_neutron_phi_deg(phi_deg: float, convention: str) -> float:
    """
    Transform the preliminary neutron-table phi value into the convention
    passed to Gepard/PARTONS.

    identity:
        phi_model = phi_table

    180-minus:
        phi_model = (180 - phi_table) mod 360

    180-shift:
        phi_model = (180 + phi_table) mod 360

    For an unpolarized cross section with the expected phi -> 360-phi symmetry,
    the latter two can be numerically degenerate.  They are nevertheless kept
    separate in the audit so the convention choice is explicit and documented.
    """
    phi = float(phi_deg) % 360.0
    if convention == "identity":
        return phi
    elif convention == "180-minus":
        return (180.0 - phi) % 360.0
    elif convention == "180-shift":
        return (180.0 + phi) % 360.0
    #endif
    raise ValueError(
        f"Unknown neutron phi convention '{convention}'. "
        f"Allowed: {NEUTRON_PHI_CONVENTIONS}"
    )
#enddef


def make_gepard_neutron_point(g, xB: float, Q2: float, t_abs: float,
                              phi_deg: float, ebeam: float,
                              phi_convention: str = "identity"):
    """
    Construct an unpolarized e n -> e n gamma point after applying the chosen
    neutron phi-convention transformation.
    """
    phi_model_deg = transform_neutron_phi_deg(phi_deg, phi_convention)
    pt = g.DataPoint(
        xB=float(xB),
        t=-abs(float(t_abs)),
        Q2=float(Q2),
        phi=math.radians(float(phi_model_deg)),
        observable="XS",
        frame="trento",
        process="en2engamma",
        exptype="fixed target",
        in1energy=float(ebeam),
        in1charge=-1,
        in1polarization=0,
        in2particle="n",
    )
    pt.prepare()
    return pt
#enddef


def evaluate_km15_neutron_point(
        args: Tuple[int, float, float, float, float, float, str]) -> Dict[str, float]:
    """
    Evaluate KM15/Gepard for a neutron point and build the exact BH quadratic.

    Gepard's KellyEFF provides neutron F1/F2 when in2particle == "n".  KM15 is
    used here only as an exploratory BH-purity classifier; the neutron CFF
    content of KM15 was not fitted to these preliminary CLAS12 neutron data.
    """
    idx, xB, Q2, t_abs, phi_deg, ebeam, phi_convention = args
    try:
        import gepard as g
        from gepard.fits import th_KM15
    except Exception as exc:
        raise RuntimeError(
            "Could not import gepard/KM15. Activate the same Python "
            "environment used for the proton analysis."
        ) from exc
    #endtry

    th = th_KM15
    pt = make_gepard_neutron_point(
        g, xB, Q2, t_abs, phi_deg, ebeam,
        phi_convention=phi_convention,
    )

    pref = float(th.PreFacSigma(pt))
    sigma_bh = pref * float(th.TBH2unp(pt))
    sigma_int = pref * float(th.TINTunp(pt))
    sigma_dvcs = pref * float(th.TDVCS2unp(pt))
    sigma_ep = sigma_bh + sigma_int + sigma_dvcs
    sigma_predict = float(th.predict(pt))
    ep_relerr = abs(sigma_ep - sigma_predict) / max(abs(sigma_predict), 1e-30)

    model = th.m
    try:
        old_f1, old_f2 = _set_eff_basis(model, 1.0, 0.0)
        A = pref * float(th.TBH2unp(pt))
        model.F1 = lambda p: 0.0
        model.F2 = lambda p: 1.0
        C = pref * float(th.TBH2unp(pt))
        model.F1 = lambda p: 1.0
        model.F2 = lambda p: 1.0
        one_one = pref * float(th.TBH2unp(pt))
        B = one_one - A - C
    finally:
        _restore_eff(model, old_f1, old_f2)
    #endtry

    f1_nom = float(model.F1(pt))
    f2_nom = float(model.F2(pt))
    sigma_bh_reco = A * f1_nom**2 + B * f1_nom * f2_nom + C * f2_nom**2
    quad_relerr = abs(sigma_bh_reco - sigma_bh) / max(abs(sigma_bh), 1e-30)

    rbh = sigma_bh / sigma_ep if sigma_ep != 0.0 else np.nan
    return {
        "_row": int(idx),
        "km15_ep": sigma_ep,
        "km15_bh": sigma_bh,
        "km15_dvcs": sigma_dvcs,
        "km15_int": sigma_int,
        "R_BH": rbh,
        "bh_delta": abs(1.0 - rbh) if np.isfinite(rbh) else np.nan,
        "bh_A": A,
        "bh_B": B,
        "bh_C": C,
        "bh_quad_relerr": quad_relerr,
        "km15_F1": f1_nom,
        "km15_F2": f2_nom,
        "km15_ep_predict": sigma_predict,
        "km15_ep_decomp_relerr": ep_relerr,
        "phi_model_deg": float(np.degrees(float(pt.phi)) % 360.0),
        "phi_convention": str(phi_convention),
    }
#enddef


def evaluate_km15_neutron_dataframe(
        df: pd.DataFrame,
        ebeam: float,
        workers: int,
        cache_path: Path,
        phi_convention: str = "identity",
        force: bool = False) -> pd.DataFrame:
    """Evaluate/cache KM15 neutron EP and exact BH coefficients for all points."""
    cache_cols = [
        "_row", "km15_ep", "km15_bh", "km15_dvcs", "km15_int",
        "R_BH", "bh_delta", "bh_A", "bh_B", "bh_C",
        "bh_quad_relerr", "km15_F1", "km15_F2",
        "km15_ep_predict", "km15_ep_decomp_relerr",
        "phi_model_deg", "phi_convention",
    ]
    if cache_path.exists() and not force:
        cached = pd.read_csv(cache_path)
        cache_matches_convention = (
            "phi_convention" in cached
            and len(cached) > 0
            and str(cached["phi_convention"].iloc[0]) == str(phi_convention)
        )
        if (
            len(cached) == len(df)
            and all(c in cached for c in cache_cols)
            and cache_matches_convention
        ):
            out = df.reset_index(drop=True).copy()
            for col in cache_cols:
                if col != "_row":
                    out[col] = cached[col].to_numpy()
                #endif
            #endfor
            return out
        #endif
    #endif

    tasks = [
        (
            i,
            float(row.xB),
            float(row.Q2),
            float(row.t_abs),
            float(row.phi_deg),
            float(ebeam),
            str(phi_convention),
        )
        for i, row in df.reset_index(drop=True).iterrows()
    ]
    if int(workers) <= 1:
        rows = [evaluate_km15_neutron_point(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=int(workers)) as ex:
            rows = list(ex.map(evaluate_km15_neutron_point, tasks, chunksize=8))
        #endwith
    #endif

    evaluated = pd.DataFrame(rows).sort_values("_row")
    out = df.reset_index(drop=True).copy()
    for col in cache_cols:
        if col != "_row":
            out[col] = evaluated[col].to_numpy()
        #endif
    #endfor

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    evaluated.to_csv(cache_path, index=False)
    return out
#enddef




def summarize_neutron_phi_convention_agreement(
        evaluated: pd.DataFrame,
        convention: str) -> Dict[str, float]:
    """
    Quantify absolute data/model agreement for one neutron phi convention.

    This is intentionally based on the measured cross sections, not only on
    R_BH = sigma_BH/sigma_EP.  The quoted systematic uncertainty is treated as
    point-to-point here because the preliminary note does not provide a
    covariance decomposition.
    """
    y = evaluated["xs"].to_numpy(float)
    err = np.sqrt(
        evaluated["stat"].to_numpy(float)**2
        + evaluated["sys"].to_numpy(float)**2
    )
    bh = evaluated["km15_bh"].to_numpy(float)
    ep = evaluated["km15_ep"].to_numpy(float)

    finite_bh = np.isfinite(y) & np.isfinite(err) & np.isfinite(bh) & (err > 0)
    finite_ep = np.isfinite(y) & np.isfinite(err) & np.isfinite(ep) & (err > 0)

    chi2_bh = float(np.sum(((y[finite_bh] - bh[finite_bh]) / err[finite_bh])**2))
    chi2_ep = float(np.sum(((y[finite_ep] - ep[finite_ep]) / err[finite_ep])**2))

    ratio_bh = y / bh
    ratio_ep = y / ep

    return {
        "phi_convention": convention,
        "N": int(len(evaluated)),
        "chi2_per_point_bh": chi2_bh / max(int(np.sum(finite_bh)), 1),
        "chi2_per_point_ep": chi2_ep / max(int(np.sum(finite_ep)), 1),
        "median_data_over_bh": float(np.nanmedian(ratio_bh)),
        "median_data_over_ep": float(np.nanmedian(ratio_ep)),
        "median_abs_frac_resid_bh": float(
            np.nanmedian(np.abs((y - bh) / np.where(y != 0.0, y, np.nan)))
        ),
        "median_abs_frac_resid_ep": float(
            np.nanmedian(np.abs((y - ep) / np.where(y != 0.0, y, np.nan)))
        ),
        "min_bh_delta": float(np.nanmin(evaluated["bh_delta"].to_numpy(float))),
        "median_bh_delta": float(np.nanmedian(evaluated["bh_delta"].to_numpy(float))),
    }
#enddef


def save_neutron_phi_convention_audit(
        convention_tables: Dict[str, pd.DataFrame],
        outdir: Path) -> str:
    """
    Compare identity, 180-minus, and 180-shift neutron phi mappings.

    The preferred convention is selected by the full-electroproduction
    chi2/point, with BH-only chi2/point used as a diagnostic.  This convention
    is then propagated into the exact BH quadratic, R_BH selection, threshold
    scan, and neutron-radius fits.
    """
    diagdir = outdir / "phi_convention_audit"
    diagdir.mkdir(parents=True, exist_ok=True)

    summary_rows = [
        summarize_neutron_phi_convention_agreement(table, convention)
        for convention, table in convention_tables.items()
    ]
    summary = pd.DataFrame(summary_rows).sort_values(
        ["chi2_per_point_ep", "chi2_per_point_bh"]
    ).reset_index(drop=True)
    summary.to_csv(diagdir / "01_neutron_phi_convention_summary.csv", index=False)

    preferred = str(summary.iloc[0]["phi_convention"])
    print("[neutron-phi] convention audit:")
    for _, row in summary.iterrows():
        print(
            f"[neutron-phi] {row['phi_convention']:10s} "
            f"EP chi2/N={row['chi2_per_point_ep']:.3f} "
            f"BH chi2/N={row['chi2_per_point_bh']:.3f} "
            f"median data/EP={row['median_data_over_ep']:.3f} "
            f"median data/BH={row['median_data_over_bh']:.3f} "
            f"min |1-RBH|={100.0*row['min_bh_delta']:.2f}%"
        )
    #endfor
    print(f"[neutron-phi] preferred convention = {preferred}")

    # 3x3 visual audit: data are identical in all three conventions, while the
    # model markers move according to the phi transformation.  Each panel
    # contains all three full-EP predictions and BH predictions at measured
    # phi_table locations so the angular reflection/shift is obvious.
    reference = next(iter(convention_tables.values()))
    fig, axes = plt.subplots(3, 3, figsize=(15.5, 11.8), sharex=True)
    axes = np.asarray(axes).ravel()

    for ax, kin_bin in zip(axes, sorted(reference["kin_bin"].unique())):
        data_part = reference.loc[reference["kin_bin"] == kin_bin].sort_values("phi_deg")
        phi = data_part["phi_deg"].to_numpy(float)
        yerr = np.sqrt(
            data_part["stat"].to_numpy(float)**2
            + data_part["sys"].to_numpy(float)**2
        )
        ax.errorbar(
            phi,
            data_part["xs"].to_numpy(float),
            yerr=yerr,
            marker="o",
            linestyle="none",
            capsize=2,
            label="CLAS12 preliminary",
        )

        marker_map = {
            "identity": "s",
            "180-minus": "^",
            "180-shift": "D",
        }
        for convention in NEUTRON_PHI_CONVENTIONS:
            p = convention_tables[convention]
            p = p.loc[p["kin_bin"] == kin_bin].sort_values("phi_deg")
            ax.plot(
                p["phi_deg"],
                p["km15_ep"],
                marker=marker_map[convention],
                linestyle="none",
                fillstyle="none",
                label=f"EP {convention}",
            )
        #endfor

        ax.set_title(
            rf"Bin {int(kin_bin)}: "
            rf"$Q^2={data_part['Q2'].mean():.2f}$, "
            rf"$x_B={data_part['xB'].mean():.3f}$, "
            rf"$|t|={data_part['t_abs'].mean():.2f}$"
        )
        ax.grid(alpha=0.2)
    #endfor

    for ax in axes[6:]:
        ax.set_xlabel(r"$\phi_{\rm table}$ [deg]")
    #endfor
    for ax in axes[::3]:
        ax.set_ylabel(r"$d^4\sigma$ [nb/GeV$^4$]")
    #endfor

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles, labels,
        loc="upper center",
        ncol=4,
        bbox_to_anchor=(0.5, 0.965),
    )
    fig.suptitle(
        "Neutron phi-convention audit: full KM15 electroproduction",
        y=0.995,
        fontsize=16,
    )
    fig.subplots_adjust(top=0.90, hspace=0.28, wspace=0.22)
    fig.savefig(diagdir / "02_neutron_phi_convention_ep_3x3.png", dpi=180)
    plt.close(fig)

    # Same audit for BH-only predictions.
    fig, axes = plt.subplots(3, 3, figsize=(15.5, 11.8), sharex=True)
    axes = np.asarray(axes).ravel()
    for ax, kin_bin in zip(axes, sorted(reference["kin_bin"].unique())):
        data_part = reference.loc[reference["kin_bin"] == kin_bin].sort_values("phi_deg")
        phi = data_part["phi_deg"].to_numpy(float)
        yerr = np.sqrt(
            data_part["stat"].to_numpy(float)**2
            + data_part["sys"].to_numpy(float)**2
        )
        ax.errorbar(
            phi,
            data_part["xs"].to_numpy(float),
            yerr=yerr,
            marker="o",
            linestyle="none",
            capsize=2,
            label="CLAS12 preliminary",
        )
        marker_map = {
            "identity": "s",
            "180-minus": "^",
            "180-shift": "D",
        }
        for convention in NEUTRON_PHI_CONVENTIONS:
            p = convention_tables[convention]
            p = p.loc[p["kin_bin"] == kin_bin].sort_values("phi_deg")
            ax.plot(
                p["phi_deg"],
                p["km15_bh"],
                marker=marker_map[convention],
                linestyle="none",
                fillstyle="none",
                label=f"BH {convention}",
            )
        #endfor
        ax.set_title(f"Bin {int(kin_bin)}")
        ax.grid(alpha=0.2)
    #endfor

    for ax in axes[6:]:
        ax.set_xlabel(r"$\phi_{\rm table}$ [deg]")
    #endfor
    for ax in axes[::3]:
        ax.set_ylabel(r"$d^4\sigma_{\rm BH}$ [nb/GeV$^4$]")
    #endfor

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles, labels,
        loc="upper center",
        ncol=4,
        bbox_to_anchor=(0.5, 0.965),
    )
    fig.suptitle(
        "Neutron phi-convention audit: KM15 Bethe-Heitler",
        y=0.995,
        fontsize=16,
    )
    fig.subplots_adjust(top=0.90, hspace=0.28, wspace=0.22)
    fig.savefig(diagdir / "03_neutron_phi_convention_bh_3x3.png", dpi=180)
    plt.close(fig)

    # Ratio audit for compact comparison.
    fig, axes = plt.subplots(3, 3, figsize=(15.5, 11.8), sharex=True)
    axes = np.asarray(axes).ravel()
    for ax, kin_bin in zip(axes, sorted(reference["kin_bin"].unique())):
        for convention in NEUTRON_PHI_CONVENTIONS:
            p = convention_tables[convention]
            p = p.loc[p["kin_bin"] == kin_bin].sort_values("phi_deg")
            ratio = p["xs"].to_numpy(float) / p["km15_bh"].to_numpy(float)
            ax.plot(
                p["phi_deg"],
                ratio,
                marker={"identity":"s","180-minus":"^","180-shift":"D"}[convention],
                linestyle="none",
                fillstyle="none",
                label=convention,
            )
        #endfor
        ax.axhline(1.0, linewidth=0.9, linestyle="--")
        ax.set_title(f"Bin {int(kin_bin)}")
        ax.grid(alpha=0.2)
    #endfor
    for ax in axes[6:]:
        ax.set_xlabel(r"$\phi_{\rm table}$ [deg]")
    #endfor
    for ax in axes[::3]:
        ax.set_ylabel("Data / BH")
    #endfor
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles, labels,
        loc="upper center",
        ncol=3,
        bbox_to_anchor=(0.5, 0.965),
    )
    fig.suptitle(
        "Neutron phi-convention audit: data/BH",
        y=0.995,
        fontsize=16,
    )
    fig.subplots_adjust(top=0.90, hspace=0.28, wspace=0.22)
    fig.savefig(diagdir / "04_neutron_phi_convention_data_over_bh_3x3.png", dpi=180)
    plt.close(fig)

    # Save one all-point table with each convention side-by-side.
    base_cols = ["kin_bin", "phi_deg", "Q2", "xB", "t_abs", "xs", "stat", "sys"]
    merged = reference[base_cols].reset_index(drop=True).copy()
    for convention in NEUTRON_PHI_CONVENTIONS:
        p = convention_tables[convention].reset_index(drop=True)
        suffix = convention.replace("-", "_")
        merged[f"phi_model_deg_{suffix}"] = p["phi_model_deg"]
        merged[f"km15_bh_{suffix}"] = p["km15_bh"]
        merged[f"km15_ep_{suffix}"] = p["km15_ep"]
        merged[f"R_BH_{suffix}"] = p["R_BH"]
        merged[f"bh_delta_{suffix}"] = p["bh_delta"]
        merged[f"data_over_bh_{suffix}"] = p["xs"] / p["km15_bh"]
        merged[f"data_over_ep_{suffix}"] = p["xs"] / p["km15_ep"]
    #endfor
    merged.to_csv(diagdir / "05_neutron_phi_convention_all_points.csv", index=False)

    # Explicitly document whether the two 180-degree mappings are degenerate
    # for this unpolarized observable.
    minus = convention_tables["180-minus"].sort_values(
        ["kin_bin", "phi_deg"]
    ).reset_index(drop=True)
    shift = convention_tables["180-shift"].sort_values(
        ["kin_bin", "phi_deg"]
    ).reset_index(drop=True)
    bh_rel = np.nanmax(
        np.abs(minus["km15_bh"].to_numpy(float) - shift["km15_bh"].to_numpy(float))
        / np.maximum(np.abs(minus["km15_bh"].to_numpy(float)), 1e-30)
    )
    ep_rel = np.nanmax(
        np.abs(minus["km15_ep"].to_numpy(float) - shift["km15_ep"].to_numpy(float))
        / np.maximum(np.abs(minus["km15_ep"].to_numpy(float)), 1e-30)
    )
    with (diagdir / "06_neutron_phi_convention_note.txt").open("w") as handle:
        handle.write(f"Preferred convention: {preferred}\n")
        handle.write(
            "Maximum relative 180-minus vs 180-shift difference:\n"
        )
        handle.write(f"  BH: {bh_rel:.6e}\n")
        handle.write(f"  EP: {ep_rel:.6e}\n")
        if bh_rel < 1e-10 and ep_rel < 1e-10:
            handle.write(
                "The two 180-degree mappings are numerically degenerate for "
                "these unpolarized points; XUU alone cannot distinguish them.\n"
            )
        #endif
    #endwith

    print(
        "[neutron-phi] 180-minus vs 180-shift max relative difference: "
        f"BH={bh_rel:.3e}, EP={ep_rel:.3e}"
    )
    print(f"[neutron-phi] diagnostics -> {diagdir}")
    return preferred
#enddef


def save_neutron_model_comparison_diagnostics(
        evaluated: pd.DataFrame,
        outdir: Path) -> None:
    """
    Compare the preliminary CLAS12 neutron cross sections directly with the
    unfitted KM15 full-electroproduction and BH-only predictions.

    All model markers are evaluated only at the measured kinematic points.
    The lower panels show data/model ratios so an absolute-normalization
    mismatch cannot be confused with R_BH = sigma_BH/sigma_EP being near one.
    """
    diagdir = outdir / "model_comparison"
    diagdir.mkdir(parents=True, exist_ok=True)

    work = evaluated.copy()
    work["data_over_bh"] = work["xs"] / work["km15_bh"]
    work["data_over_ep"] = work["xs"] / work["km15_ep"]
    work["data_total_err"] = np.sqrt(work["stat"]**2 + work["sys"]**2)
    work["data_over_bh_err"] = work["data_total_err"] / np.abs(work["km15_bh"])
    work["data_over_ep_err"] = work["data_total_err"] / np.abs(work["km15_ep"])
    work.to_csv(diagdir / "01_neutron_data_vs_km15_all_points.csv", index=False)

    # Per-bin edge audit: use the lowest- and highest-phi measured point in
    # each of the nine (Q2,xB,t) bins.
    edge_rows = []
    for kin_bin, part in work.groupby("kin_bin", sort=True):
        p = part.sort_values("phi_deg").reset_index(drop=True)
        for edge_name, row in [("low_phi", p.iloc[0]), ("high_phi", p.iloc[-1])]:
            edge_rows.append({
                "kin_bin": int(kin_bin),
                "edge": edge_name,
                "phi_deg": float(row["phi_deg"]),
                "Q2_GeV2": float(row["Q2"]),
                "xB": float(row["xB"]),
                "t_abs_GeV2": float(row["t_abs"]),
                "data_nb_GeV4": float(row["xs"]),
                "km15_bh_nb_GeV4": float(row["km15_bh"]),
                "km15_ep_nb_GeV4": float(row["km15_ep"]),
                "data_over_bh": float(row["data_over_bh"]),
                "data_over_ep": float(row["data_over_ep"]),
                "R_BH": float(row["R_BH"]),
                "bh_delta": float(row["bh_delta"]),
            })
        #endfor
    #endfor
    edge = pd.DataFrame(edge_rows)
    edge.to_csv(diagdir / "02_neutron_edge_point_audit.csv", index=False)

    print("[neutron-model] direct edge-point data/model audit:")
    for _, row in edge.iterrows():
        print(
            f"[neutron-model] bin={int(row['kin_bin'])} "
            f"{row['edge']:8s} phi={row['phi_deg']:6.1f} deg "
            f"data/BH={row['data_over_bh']:.3f} "
            f"data/EP={row['data_over_ep']:.3f} "
            f"R_BH={row['R_BH']:.3f}"
        )
    #endfor

    # Individual two-panel plots, one for each kinematic bin.
    for kin_bin, part in work.groupby("kin_bin", sort=True):
        p = part.sort_values("phi_deg")
        phi = p["phi_deg"].to_numpy(float)
        y = p["xs"].to_numpy(float)
        yerr = p["data_total_err"].to_numpy(float)

        fig, (ax_top, ax_ratio) = plt.subplots(
            2, 1, figsize=(7.4, 6.6), sharex=True,
            gridspec_kw={"height_ratios": [2.2, 1.0], "hspace": 0.05},
        )

        ax_top.errorbar(
            phi, y, yerr=yerr,
            marker="o", linestyle="none", capsize=2,
            label="CLAS12 preliminary",
        )
        ax_top.plot(
            phi, p["km15_bh"].to_numpy(float),
            marker="s", linestyle="none", fillstyle="none",
            label="KM15 BH only",
        )
        ax_top.plot(
            phi, p["km15_ep"].to_numpy(float),
            marker="^", linestyle="none", fillstyle="none",
            label="KM15 full EP",
        )
        ax_top.set_ylabel(r"$d^4\sigma$ [nb/GeV$^4$]")
        ax_top.set_title(
            rf"Neutron bin {int(kin_bin)}: "
            rf"$Q^2={p['Q2'].mean():.3f}$ GeV$^2$, "
            rf"$x_B={p['xB'].mean():.3f}$, "
            rf"$|t|={p['t_abs'].mean():.3f}$ GeV$^2$"
        )
        ax_top.grid(alpha=0.2)
        ax_top.legend()

        ax_ratio.errorbar(
            phi,
            p["data_over_bh"].to_numpy(float),
            yerr=p["data_over_bh_err"].to_numpy(float),
            marker="s", linestyle="none", capsize=2,
            label="Data / BH",
        )
        ax_ratio.errorbar(
            phi,
            p["data_over_ep"].to_numpy(float),
            yerr=p["data_over_ep_err"].to_numpy(float),
            marker="^", linestyle="none", capsize=2,
            label="Data / full EP",
        )
        ax_ratio.axhline(1.0, linewidth=0.9, linestyle="--")
        ax_ratio.set_xlabel(r"$\phi$ [deg]")
        ax_ratio.set_ylabel("Data / model")
        ax_ratio.grid(alpha=0.2)
        ax_ratio.legend(ncol=2, fontsize=9)

        fig.subplots_adjust(top=0.90)
        fig.savefig(
            diagdir / f"03_bin_{int(kin_bin):02d}_data_vs_km15.png",
            dpi=180,
        )
        plt.close(fig)
    #endfor

    # 3x3 cross-section summary matching the nine analysis-note kinematic bins.
    fig, axes = plt.subplots(3, 3, figsize=(15.0, 11.5), sharex=True)
    axes = np.asarray(axes).ravel()
    for ax, (kin_bin, part) in zip(
            axes, work.groupby("kin_bin", sort=True)):
        p = part.sort_values("phi_deg")
        phi = p["phi_deg"].to_numpy(float)
        ax.errorbar(
            phi, p["xs"].to_numpy(float),
            yerr=p["data_total_err"].to_numpy(float),
            marker="o", linestyle="none", capsize=2,
            label="CLAS12 preliminary",
        )
        ax.plot(
            phi, p["km15_bh"].to_numpy(float),
            marker="s", linestyle="none", fillstyle="none",
            label="KM15 BH",
        )
        ax.plot(
            phi, p["km15_ep"].to_numpy(float),
            marker="^", linestyle="none", fillstyle="none",
            label="KM15 full EP",
        )
        ax.set_title(
            rf"Bin {int(kin_bin)}: "
            rf"$Q^2={p['Q2'].mean():.2f}$, "
            rf"$x_B={p['xB'].mean():.3f}$, "
            rf"$|t|={p['t_abs'].mean():.2f}$"
        )
        ax.grid(alpha=0.2)
    #endfor
    for ax in axes[6:]:
        ax.set_xlabel(r"$\phi$ [deg]")
    #endfor
    for ax in axes[::3]:
        ax.set_ylabel(r"$d^4\sigma$ [nb/GeV$^4$]")
    #endfor

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles, labels, loc="upper center", ncol=3,
        bbox_to_anchor=(0.5, 0.965),
    )
    fig.suptitle(
        "Preliminary CLAS12 neutron cross sections vs unfitted KM15",
        y=0.995,
        fontsize=16,
    )
    fig.subplots_adjust(top=0.90, hspace=0.28, wspace=0.22)
    fig.savefig(diagdir / "04_neutron_data_vs_km15_3x3.png", dpi=180)
    plt.close(fig)

    # Separate 3x3 ratio summary.  This is deliberately not folded into the
    # cross-section canvas so the absolute scale and the ratio scale are both
    # readable.
    fig, axes = plt.subplots(3, 3, figsize=(15.0, 11.5), sharex=True)
    axes = np.asarray(axes).ravel()
    for ax, (kin_bin, part) in zip(
            axes, work.groupby("kin_bin", sort=True)):
        p = part.sort_values("phi_deg")
        phi = p["phi_deg"].to_numpy(float)
        ax.errorbar(
            phi,
            p["data_over_bh"].to_numpy(float),
            yerr=p["data_over_bh_err"].to_numpy(float),
            marker="s", linestyle="none", capsize=2,
            label="Data / BH",
        )
        ax.errorbar(
            phi,
            p["data_over_ep"].to_numpy(float),
            yerr=p["data_over_ep_err"].to_numpy(float),
            marker="^", linestyle="none", capsize=2,
            label="Data / full EP",
        )
        ax.axhline(1.0, linewidth=0.9, linestyle="--")
        ax.set_title(f"Bin {int(kin_bin)}")
        ax.grid(alpha=0.2)
    #endfor
    for ax in axes[6:]:
        ax.set_xlabel(r"$\phi$ [deg]")
    #endfor
    for ax in axes[::3]:
        ax.set_ylabel("Data / model")
    #endfor

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles, labels, loc="upper center", ncol=2,
        bbox_to_anchor=(0.5, 0.965),
    )
    fig.suptitle(
        "Preliminary CLAS12 neutron data/KM15 ratios",
        y=0.995,
        fontsize=16,
    )
    fig.subplots_adjust(top=0.90, hspace=0.28, wspace=0.22)
    fig.savefig(diagdir / "05_neutron_data_over_km15_3x3.png", dpi=180)
    plt.close(fig)

    # Compact numerical summary of the absolute normalization discrepancy.
    summary_rows = []
    for kin_bin, part in work.groupby("kin_bin", sort=True):
        summary_rows.append({
            "kin_bin": int(kin_bin),
            "median_data_over_bh": float(np.nanmedian(part["data_over_bh"])),
            "median_data_over_ep": float(np.nanmedian(part["data_over_ep"])),
            "min_bh_delta": float(np.nanmin(part["bh_delta"])),
            "max_bh_delta": float(np.nanmax(part["bh_delta"])),
        })
    #endfor
    pd.DataFrame(summary_rows).to_csv(
        diagdir / "06_neutron_model_ratio_summary_by_bin.csv",
        index=False,
    )

    print(f"[neutron-model] diagnostics -> {diagdir}")
#enddef


def neutron_atac_ge(q: np.ndarray, A: float, B: float, C: float) -> np.ndarray:
    """
    Atac et al. (Nature Communications 12, 1759 (2021)) neutron GE form:

      GE_n(Q2) = (1 + Q2/A)^(-2) * B*tau/(1 + C*tau).

    B controls the Q2->0 slope and therefore <r_E,n^2>.
    """
    q = np.asarray(q, dtype=float)
    tau = q / (4.0 * MP2)
    return (1.0 + q / float(A))**(-2) * (float(B) * tau) / (
        1.0 + float(C) * tau
    )
#enddef


def neutron_charge_radius_sq_from_atac_B(B: float) -> float:
    """Return <r_E,n^2> in fm^2 from the Atac-form low-Q2 slope."""
    slope_gev2 = float(B) / (4.0 * MP2)
    return -6.0 * slope_gev2 * HBARC**2
#enddef


def fit_neutron_atac_bh(
        data: pd.DataFrame,
        magnetic_mode: str = "kelly",
        bh_systematic_fraction: float = 0.05) -> Dict[str, object]:
    """
    Fit neutron BH-dominated cross sections with the Atac GE_n form.

    magnetic_mode="kelly":
      neutron GM is fixed to Gepard/Kelly.  This is analogous in spirit to the
      proton Moradi Fit 8: use external magnetic information to isolate the
      electric slope.

    magnetic_mode="ip1":
      GM_n/mu_n = 1/(1+m1*Q2) is fitted simultaneously, yielding an exploratory
      neutron magnetic radius as well.
    """
    if magnetic_mode not in {"kelly", "ip1"}:
        raise ValueError(f"unknown neutron magnetic_mode={magnetic_mode}")
    #endif

    q = data["t_abs"].to_numpy(float)
    tau = q / (4.0 * MP2)
    inv = 1.0 / (1.0 + tau)
    y = data["xs"].to_numpy(float)
    stat = data["stat"].to_numpy(float)
    sys = data["sys"].to_numpy(float)
    err = np.sqrt(stat**2 + sys**2 + (float(bh_systematic_fraction) * y)**2)

    bh_A = data["bh_A"].to_numpy(float)
    bh_B = data["bh_B"].to_numpy(float)
    bh_C = data["bh_C"].to_numpy(float)

    # Nominal Gepard neutron EFF is Kelly.  Recover Sachs GM from F1/F2.
    f1_kelly = data["km15_F1"].to_numpy(float)
    f2_kelly = data["km15_F2"].to_numpy(float)
    gm_kelly = f1_kelly + f2_kelly

    names = ["A_GE", "B_GE", "C_GE"]
    p0 = [0.505, 1.655, 0.909]
    if magnetic_mode == "ip1":
        names.append("m1_GM")
        p0.append(
            sachs_first_coefficient_from_radius(
                SACHS_INITIAL_RADIUS_FM, "IP1"
            )
        )
    #endif

    def chi2(*values):
        A_GE = float(values[0])
        B_GE = float(values[1])
        C_GE = float(values[2])
        ge = neutron_atac_ge(q, A_GE, B_GE, C_GE)

        if magnetic_mode == "kelly":
            gm = gm_kelly
        else:
            gm = MU_N * sachs_family_value(
                q, np.asarray([values[3]], dtype=float), "IP1"
            )
        #endif

        f1 = (ge + tau * gm) * inv
        f2 = (gm - ge) * inv
        pred = bh_from_f1f2(bh_A, bh_B, bh_C, f1, f2)
        pull = (pred - y) / err
        return float(np.dot(pull, pull))
    #enddef

    m = Minuit(chi2, *p0, name=tuple(names))
    m.errordef = Minuit.LEAST_SQUARES

    # Broad numerical guards.  The Nature fit is near A~0.5, B~1.7, C~0.9;
    # these limits are intentionally much wider and are not physics priors.
    m.limits["A_GE"] = (0.05, 10.0)
    m.limits["B_GE"] = (-10.0, 10.0)
    m.limits["C_GE"] = (-20.0, 30.0)
    if magnetic_mode == "ip1":
        lo = sachs_first_coefficient_from_radius(SACHS_MIN_RADIUS_FM, "IP1")
        hi = sachs_first_coefficient_from_radius(SACHS_MAX_RADIUS_FM, "IP1")
        m.limits["m1_GM"] = (min(lo, hi), max(lo, hi))
    #endif

    m.migrad()
    m.hesse()

    pars = np.asarray([float(m.values[n]) for n in names], dtype=float)
    cov = np.full((len(names), len(names)), np.nan)
    if m.covariance is not None:
        for i, ni in enumerate(names):
            for j, nj in enumerate(names):
                cov[i, j] = float(m.covariance[ni, nj])
            #endfor
        #endfor
    #endif

    rn2 = neutron_charge_radius_sq_from_atac_B(pars[1])
    rn2_err = np.nan
    if np.isfinite(cov[1, 1]) and cov[1, 1] >= 0.0:
        scale = 6.0 * HBARC**2 / (4.0 * MP2)
        rn2_err = scale * math.sqrt(cov[1, 1])
    #endif

    rM = np.nan
    rM_err = np.nan
    if magnetic_mode == "ip1":
        def rM_fun(p):
            return sachs_family_radius(np.asarray([p[3]], dtype=float), "IP1")
        #enddef
        rM, rM_err = propagate_scalar(rM_fun, pars, cov)
    #endif

    ndof = max(1, len(data) - len(names))
    return {
        "magnetic_mode": magnetic_mode,
        "N": int(len(data)),
        "n_parameters": int(len(names)),
        "chi2": float(m.fval),
        "ndof": int(ndof),
        "chi2_ndof": float(m.fval / ndof),
        "valid": bool(m.valid and np.isfinite(rn2)),
        "accurate": bool(m.accurate),
        "rn2_fm2": float(rn2),
        "rn2_fit_err_fm2": float(rn2_err),
        "rM_fm": float(rM),
        "rM_fit_err_fm": float(rM_err),
        "A_GE_GeV2": float(pars[0]),
        "B_GE": float(pars[1]),
        "C_GE": float(pars[2]),
        "m1_GM_GeVminus2": (
            float(pars[3]) if magnetic_mode == "ip1" else np.nan
        ),
    }
#enddef



def save_neutron_bh_sensitivity_study(
        evaluated: pd.DataFrame,
        outdir: Path) -> None:
    """
    Diagnose what the measured neutron BH cross sections actually constrain.

    The exact precomputed quadratic
        sigma_BH = A F1^2 + B F1 F2 + C F2^2
    is differentiated numerically with respect to Sachs GE_n and GM_n at each
    measured point, using the nominal Kelly neutron elastic form factors as
    the expansion point.

    Two complementary sensitivities are saved:
      1) logarithmic/fractional:
           S_GE = d ln(sigma) / d ln(|GE|)
           S_GM = d ln(sigma) / d ln(|GM|)
         These answer how a fractional FF change modifies the cross section.
      2) absolute:
           D_GE = d sigma / d GE
           D_GM = d sigma / d GM
         These remain meaningful when GE_n is small and the logarithmic GE
         sensitivity is correspondingly suppressed.

    The local two-parameter Fisher matrix in each kinematic bin is also saved.
    Its eigenvalues/condition number quantify the GE/GM degeneracy without
    attempting a Q2->0 radius extrapolation.
    """
    diagdir = outdir / "bh_sensitivity"
    diagdir.mkdir(parents=True, exist_ok=True)

    work = evaluated.copy().reset_index(drop=True)
    q = work["t_abs"].to_numpy(float)
    tau = q / (4.0 * MP2)
    inv = 1.0 / (1.0 + tau)

    f1 = work["km15_F1"].to_numpy(float)
    f2 = work["km15_F2"].to_numpy(float)
    ge0 = f1 - tau * f2
    gm0 = f1 + f2

    A = work["bh_A"].to_numpy(float)
    B = work["bh_B"].to_numpy(float)
    C = work["bh_C"].to_numpy(float)

    def sigma_from_sachs(ge, gm):
        ff1 = (ge + tau * gm) * inv
        ff2 = (gm - ge) * inv
        return bh_from_f1f2(A, B, C, ff1, ff2)
    #enddef

    sigma0 = sigma_from_sachs(ge0, gm0)

    # Central finite differences.  The step is relative to the natural neutron
    # FF scale, with a floor so GE_n~0 does not lead to a vanishing step.
    h_ge = 1.0e-5 * np.maximum(np.abs(ge0), 0.05)
    h_gm = 1.0e-5 * np.maximum(np.abs(gm0), 0.20)

    dsdge = (
        sigma_from_sachs(ge0 + h_ge, gm0)
        - sigma_from_sachs(ge0 - h_ge, gm0)
    ) / (2.0 * h_ge)
    dsdgm = (
        sigma_from_sachs(ge0, gm0 + h_gm)
        - sigma_from_sachs(ge0, gm0 - h_gm)
    ) / (2.0 * h_gm)

    with np.errstate(divide="ignore", invalid="ignore"):
        s_ge = ge0 * dsdge / sigma0
        s_gm = gm0 * dsdgm / sigma0
        ratio_abs = np.abs(dsdge) / np.maximum(np.abs(dsdgm), 1.0e-30)
        ratio_frac = np.abs(s_ge) / np.maximum(np.abs(s_gm), 1.0e-30)
    #endwith

    err = np.sqrt(
        work["stat"].to_numpy(float)**2
        + work["sys"].to_numpy(float)**2
        + (0.05 * work["xs"].to_numpy(float))**2
    )
    info_ge = (dsdge / err)**2
    info_gm = (dsdgm / err)**2
    info_cross = dsdge * dsdgm / err**2

    work["GE_n_Kelly"] = ge0
    work["GM_n_Kelly"] = gm0
    work["bh_sigma_Kelly"] = sigma0
    work["dSigma_dGE"] = dsdge
    work["dSigma_dGM"] = dsdgm
    work["dlnSigma_dlnAbsGE"] = s_ge
    work["dlnSigma_dlnAbsGM"] = s_gm
    work["abs_derivative_ratio_GE_over_GM"] = ratio_abs
    work["fractional_sensitivity_ratio_GE_over_GM"] = ratio_frac
    work["point_error_for_sensitivity"] = err
    work["fisher_GE_GE"] = info_ge
    work["fisher_GE_GM"] = info_cross
    work["fisher_GM_GM"] = info_gm
    work.to_csv(diagdir / "01_neutron_bh_point_sensitivities.csv", index=False)

    # Moradi-style neutron F1/F2 turn-off demonstration over the full measured
    # neutron sample.  Here the nominal Kelly neutron F1/F2 values are used only
    # as the electromagnetic expansion point for visualizing BH leverage.
    t11_n = A * f1**2
    t12_n = B * f1 * f2
    t22_n = C * f2**2
    full_n = t11_n + t12_n + t22_n

    with np.errstate(divide="ignore", invalid="ignore"):
        f2only_over_f1only_n = np.divide(
            t22_n,
            t11_n,
            out=np.full_like(t22_n, np.nan, dtype=float),
            where=np.abs(t11_n) > 1.0e-30,
        )
    #endwith

    fig, axes = plt.subplots(
        2, 1,
        figsize=(9.2, 7.8),
        sharex=True,
        gridspec_kw={"height_ratios": [2.1, 1.0]},
    )
    ax, rax = axes

    ax.scatter(q, full_n, s=22, alpha=0.45, label="Full BH")
    ax.scatter(q, t22_n, s=22, alpha=0.45, marker="s",
               label=r"$F_1^n=0$ ($F_2^n$ only)")
    ax.scatter(q, t11_n, s=22, alpha=0.45, marker="^",
               label=r"$F_2^n=0$ ($F_1^n$ only)")

    rax.scatter(q, f2only_over_f1only_n, s=22, alpha=0.50)

    # The neutron table consists of repeated phi measurements in a small number
    # of (Q2,xB,t) bins.  Plot the per-bin medians as connected guides.
    nplot = pd.DataFrame({
        "kin_bin": work["kin_bin"].to_numpy(int),
        "t": q,
        "full": full_n,
        "f2only": t22_n,
        "f1only": t11_n,
        "ratio": f2only_over_f1only_n,
    })
    nmed = (
        nplot.groupby("kin_bin", sort=True)
        .median(numeric_only=True)
        .sort_values("t")
    )
    ax.plot(nmed["t"], nmed["full"], linewidth=2.0)
    ax.plot(nmed["t"], nmed["f2only"], linewidth=2.0, linestyle="--")
    ax.plot(nmed["t"], nmed["f1only"], linewidth=2.0, linestyle=":")
    rax.plot(nmed["t"], nmed["ratio"], linewidth=2.0)

    ax.set_ylabel(r"BH cross section")
    ax.set_title(
        r"Neutron electromagnetic leverage across the measured $|t|$ range"
    )
    ax.grid(alpha=0.18)
    ax.legend(fontsize=9, loc="best")

    rax.axhline(1.0, linewidth=0.8, linestyle="--")
    rax.set_xlabel(r"$|t|$ (GeV$^2$)")
    rax.set_ylabel(
        r"$\sigma_{\rm BH}(F_1^n=0)\,/\,\sigma_{\rm BH}(F_2^n=0)$"
    )
    rax.grid(alpha=0.18)

    fig.subplots_adjust(
        top=0.93, bottom=0.10, left=0.10, right=0.985,
        hspace=0.08,
    )
    fig.savefig(
        diagdir / "03_neutron_F1_F2_turnoff_vs_t.png",
        dpi=300,
    )
    plt.close(fig)

    pd.DataFrame({
        "kin_bin": work["kin_bin"].to_numpy(int),
        "t_abs": q,
        "sigma_full_bh": full_n,
        "sigma_F1_zero_F2_only": t22_n,
        "sigma_F2_zero_F1_only": t11_n,
        "sigma_F2only_over_F1only": f2only_over_f1only_n,
    }).to_csv(
        diagdir / "03_neutron_F1_F2_turnoff_vs_t.csv",
        index=False,
    )

    # Per-bin Fisher diagnostics.  This deliberately treats GE and GM as two
    # local amplitudes within one (Q2,xB,t) bin, avoiding any radius model.
    fisher_rows = []
    for kin_bin, part in work.groupby("kin_bin", sort=True):
        I = np.array([
            [np.nansum(part["fisher_GE_GE"]),
             np.nansum(part["fisher_GE_GM"])],
            [np.nansum(part["fisher_GE_GM"]),
             np.nansum(part["fisher_GM_GM"])],
        ], dtype=float)
        eig = np.linalg.eigvalsh(I)
        eig = np.sort(eig)
        condition = (
            float(eig[-1] / eig[0])
            if eig[0] > 0.0 and np.isfinite(eig[0])
            else np.inf
        )
        covariance = np.full((2, 2), np.nan)
        corr = np.nan
        sigma_ge = np.nan
        sigma_gm = np.nan
        try:
            covariance = np.linalg.inv(I)
            sigma_ge = math.sqrt(max(0.0, covariance[0, 0]))
            sigma_gm = math.sqrt(max(0.0, covariance[1, 1]))
            denom = sigma_ge * sigma_gm
            if denom > 0.0:
                corr = float(covariance[0, 1] / denom)
            #endif
        except np.linalg.LinAlgError:
            pass
        #endtry

        fisher_rows.append({
            "kin_bin": int(kin_bin),
            "N": int(len(part)),
            "Q2_mean_GeV2": float(np.nanmean(part["Q2"])),
            "xB_mean": float(np.nanmean(part["xB"])),
            "t_abs_mean_GeV2": float(np.nanmean(part["t_abs"])),
            "GE_n_Kelly_mean": float(np.nanmean(part["GE_n_Kelly"])),
            "GM_n_Kelly_mean": float(np.nanmean(part["GM_n_Kelly"])),
            "median_abs_dlnSigma_dlnAbsGE": float(
                np.nanmedian(np.abs(part["dlnSigma_dlnAbsGE"]))
            ),
            "median_abs_dlnSigma_dlnAbsGM": float(
                np.nanmedian(np.abs(part["dlnSigma_dlnAbsGM"]))
            ),
            "median_fractional_sensitivity_ratio_GE_over_GM": float(
                np.nanmedian(part["fractional_sensitivity_ratio_GE_over_GM"])
            ),
            "fisher_eigenvalue_small": float(eig[0]),
            "fisher_eigenvalue_large": float(eig[-1]),
            "fisher_condition_number": condition,
            "local_sigma_GE": sigma_ge,
            "local_sigma_GM": sigma_gm,
            "local_GE_GM_correlation": corr,
        })
    #endfor

    fisher = pd.DataFrame(fisher_rows)
    fisher.to_csv(diagdir / "02_neutron_bh_fisher_by_kinematic_bin.csv", index=False)

    # 3x3 fractional-sensitivity canvas.  This is the most direct visualization
    # of whether GE_n or GM_n controls each measured angular point.
    fig, axes = plt.subplots(3, 3, figsize=(15.0, 11.5), sharex=True)
    axes = np.asarray(axes).ravel()
    for ax, (kin_bin, part) in zip(
            axes, work.groupby("kin_bin", sort=True)):
        p = part.sort_values("phi_deg")
        phi = p["phi_deg"].to_numpy(float)
        ax.plot(
            phi,
            np.abs(p["dlnSigma_dlnAbsGE"].to_numpy(float)),
            marker="o", linestyle="none",
            label=r"$|\partial\ln\sigma/\partial\ln|G_E^n||$",
        )
        ax.plot(
            phi,
            np.abs(p["dlnSigma_dlnAbsGM"].to_numpy(float)),
            marker="s", linestyle="none", fillstyle="none",
            label=r"$|\partial\ln\sigma/\partial\ln|G_M^n||$",
        )
        ax.set_title(
            rf"Bin {int(kin_bin)}: "
            rf"$|t|={p['t_abs'].mean():.3f}$ GeV$^2$"
        )
        ax.grid(alpha=0.2)
    #endfor
    for ax in axes[6:]:
        ax.set_xlabel(r"$\phi$ [deg]")
    #endfor
    for ax in axes[::3]:
        ax.set_ylabel("Absolute fractional sensitivity")
    #endfor
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles, labels, loc="upper center", ncol=2,
        bbox_to_anchor=(0.5, 0.965),
    )
    fig.suptitle(
        r"Neutron BH local sensitivity to $G_E^n$ and $G_M^n$",
        y=0.995, fontsize=16,
    )
    fig.subplots_adjust(top=0.90, hspace=0.28, wspace=0.22)
    fig.savefig(
        diagdir / "03_neutron_bh_fractional_sensitivity_3x3.png",
        dpi=180,
    )
    plt.close(fig)

    # Absolute derivative canvas is essential because GE_n is small: a small
    # logarithmic GE sensitivity can arise simply from the GE prefactor.
    fig, axes = plt.subplots(3, 3, figsize=(15.0, 11.5), sharex=True)
    axes = np.asarray(axes).ravel()
    for ax, (kin_bin, part) in zip(
            axes, work.groupby("kin_bin", sort=True)):
        p = part.sort_values("phi_deg")
        phi = p["phi_deg"].to_numpy(float)
        ax.plot(
            phi, np.abs(p["dSigma_dGE"].to_numpy(float)),
            marker="o", linestyle="none",
            label=r"$|\partial\sigma/\partial G_E^n|$",
        )
        ax.plot(
            phi, np.abs(p["dSigma_dGM"].to_numpy(float)),
            marker="s", linestyle="none", fillstyle="none",
            label=r"$|\partial\sigma/\partial G_M^n|$",
        )
        ax.set_title(f"Bin {int(kin_bin)}")
        ax.grid(alpha=0.2)
    #endfor
    for ax in axes[6:]:
        ax.set_xlabel(r"$\phi$ [deg]")
    #endfor
    for ax in axes[::3]:
        ax.set_ylabel(r"Absolute derivative [cross section / FF]")
    #endfor
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles, labels, loc="upper center", ncol=2,
        bbox_to_anchor=(0.5, 0.965),
    )
    fig.suptitle(
        r"Neutron BH absolute $G_E^n/G_M^n$ leverage",
        y=0.995, fontsize=16,
    )
    fig.subplots_adjust(top=0.90, hspace=0.28, wspace=0.22)
    fig.savefig(
        diagdir / "04_neutron_bh_absolute_sensitivity_3x3.png",
        dpi=180,
    )
    plt.close(fig)

    # Fisher condition number and GE-GM correlation summarize how nearly
    # singular a two-amplitude extraction is in each kinematic bin.
    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    ax.plot(
        fisher["kin_bin"],
        fisher["fisher_condition_number"],
        marker="o", linestyle="none",
    )
    ax.set_yscale("log")
    ax.set_xlabel("Kinematic bin")
    ax.set_ylabel("Fisher condition number")
    ax.set_title(r"Local neutron $G_E^n/G_M^n$ identifiability")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(
        diagdir / "05_neutron_bh_fisher_condition_by_bin.png",
        dpi=180,
    )
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    ax.plot(
        fisher["kin_bin"],
        fisher["local_GE_GM_correlation"],
        marker="o", linestyle="none",
    )
    ax.axhline(0.0, linewidth=0.8, linestyle=":")
    ax.set_ylim(-1.05, 1.05)
    ax.set_xlabel("Kinematic bin")
    ax.set_ylabel(r"Local $G_E^n$-$G_M^n$ correlation")
    ax.set_title("Two-amplitude Fisher correlation")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(
        diagdir / "06_neutron_bh_GE_GM_correlation_by_bin.png",
        dpi=180,
    )
    plt.close(fig)

    # Compact threshold summary: how the sensitivity composition changes under
    # the same KM15 BH-purity cuts used by the exploratory radius scan.
    thresholds = [
        0.05, 0.075, 0.10, 0.15, 0.20, 0.30, 0.50, 0.75, 1.00
    ]
    finite_delta = work.loc[
        np.isfinite(work["bh_delta"]), "bh_delta"
    ].to_numpy(float)
    if len(finite_delta):
        thresholds.append(float(np.max(finite_delta) * 1.001))
    #endif
    threshold_rows = []
    for threshold in sorted(set(thresholds)):
        p = work.loc[
            np.isfinite(work["bh_delta"])
            & (work["bh_delta"] <= float(threshold))
        ]
        if not len(p):
            continue
        #endif
        threshold_rows.append({
            "bh_threshold": float(threshold),
            "bh_threshold_percent": 100.0 * float(threshold),
            "N": int(len(p)),
            "median_abs_fractional_GE_sensitivity": float(
                np.nanmedian(np.abs(p["dlnSigma_dlnAbsGE"]))
            ),
            "median_abs_fractional_GM_sensitivity": float(
                np.nanmedian(np.abs(p["dlnSigma_dlnAbsGM"]))
            ),
            "median_fractional_GE_over_GM": float(
                np.nanmedian(p["fractional_sensitivity_ratio_GE_over_GM"])
            ),
            "sum_fisher_GE_GE": float(np.nansum(p["fisher_GE_GE"])),
            "sum_fisher_GM_GM": float(np.nansum(p["fisher_GM_GM"])),
        })
    #endfor
    pd.DataFrame(threshold_rows).to_csv(
        diagdir / "07_neutron_bh_sensitivity_by_threshold.csv",
        index=False,
    )

    print(f"[neutron-sensitivity] diagnostics -> {diagdir}")
#enddef



def fit_neutron_magnetic_only_bh(
        data: pd.DataFrame,
        ge_mode: str = "atac2021",
        bh_systematic_fraction: float = 0.05) -> Dict[str, object]:
    """
    Exploratory one-unknown neutron BH fit.

    GE_n is fixed externally and only GM_n is fitted:
        GM_n(Q2) / mu_n = 1 / (1 + m1 Q2).

    ge_mode="atac2021" fixes GE_n to the central Atac et al. 2021 form
        A=0.505 GeV^2, B=1.655, C=0.909.
    ge_mode="zero" sets GE_n identically to zero as a limiting diagnostic.

    This is deliberately a cross-section-level fit rather than a point-by-point
    inversion, so the measured uncertainties remain in their native space.
    """
    if ge_mode not in {"atac2021", "zero"}:
        raise ValueError(f"unknown neutron magnetic-only ge_mode={ge_mode}")
    #endif

    q = data["t_abs"].to_numpy(float)
    tau = q / (4.0 * MP2)
    inv = 1.0 / (1.0 + tau)
    y = data["xs"].to_numpy(float)
    stat = data["stat"].to_numpy(float)
    sys = data["sys"].to_numpy(float)
    err = np.sqrt(
        stat**2 + sys**2 + (float(bh_systematic_fraction) * y)**2
    )

    bh_A = data["bh_A"].to_numpy(float)
    bh_B = data["bh_B"].to_numpy(float)
    bh_C = data["bh_C"].to_numpy(float)

    if ge_mode == "atac2021":
        ge_fixed = neutron_atac_ge(q, 0.505, 1.655, 0.909)
    else:
        ge_fixed = np.zeros_like(q)
    #endif

    m1_initial = sachs_first_coefficient_from_radius(
        SACHS_INITIAL_RADIUS_FM, "IP1"
    )

    def chi2(m1_GM):
        gm = MU_N * sachs_family_value(
            q, np.asarray([float(m1_GM)], dtype=float), "IP1"
        )
        f1 = (ge_fixed + tau * gm) * inv
        f2 = (gm - ge_fixed) * inv
        pred = bh_from_f1f2(bh_A, bh_B, bh_C, f1, f2)
        pull = (pred - y) / err
        return float(np.dot(pull, pull))
    #enddef

    m = Minuit(chi2, m1_GM=m1_initial)
    m.errordef = Minuit.LEAST_SQUARES
    lo = sachs_first_coefficient_from_radius(SACHS_MIN_RADIUS_FM, "IP1")
    hi = sachs_first_coefficient_from_radius(SACHS_MAX_RADIUS_FM, "IP1")
    m.limits["m1_GM"] = (min(lo, hi), max(lo, hi))
    m.migrad()
    m.hesse()

    m1 = float(m.values["m1_GM"])
    m1_err = float(m.errors["m1_GM"]) if np.isfinite(m.errors["m1_GM"]) else np.nan

    pars = np.asarray([m1], dtype=float)
    cov = np.asarray([[m1_err**2]], dtype=float)
    def rM_fun(p):
        return sachs_family_radius(np.asarray([p[0]], dtype=float), "IP1")
    #enddef
    rM, rM_err = propagate_scalar(rM_fun, pars, cov)

    ndof = max(1, len(data) - 1)
    return {
        "ge_mode": ge_mode,
        "gm_family": "IP1",
        "N": int(len(data)),
        "n_parameters": 1,
        "chi2": float(m.fval),
        "ndof": int(ndof),
        "chi2_ndof": float(m.fval / ndof),
        "valid": bool(m.valid and np.isfinite(rM)),
        "accurate": bool(m.accurate),
        "m1_GM_GeVminus2": m1,
        "m1_GM_err_GeVminus2": m1_err,
        "rM_fm": float(rM),
        "rM_fit_err_fm": float(rM_err),
    }
#enddef



def fit_neutron_magnetic_family_bh(
        data: pd.DataFrame,
        family: str,
        override_y: Optional[np.ndarray] = None,
        statistical_only: bool = True,
        bh_systematic_fraction: float = 0.05,
        return_parameters: bool = False,
        return_uncertainty: bool = False):
    """
    Fit one normalized neutron magnetic Sachs family through the exact BH
    cross section while fixing GE_n to the central Atac et al. 2021 form.

        GM_n(Q2) = mu_n * g_M(Q2),  g_M(0)=1.

    This is the one-form-factor analogue of the proton Sachs-family fitter.
    """
    nshape = int(re.findall(r"\d+", family)[0])
    names = [f"m{i}" for i in range(1, nshape + 1)]
    p0 = np.zeros(nshape, dtype=float)
    p0[0] = sachs_first_coefficient_from_radius(
        SACHS_INITIAL_RADIUS_FM, family
    )

    q = data["t_abs"].to_numpy(float)
    tau = q / (4.0 * MP2)
    inv = 1.0 / (1.0 + tau)
    ge = neutron_atac_ge(q, 0.505, 1.655, 0.909)
    y = (
        np.asarray(override_y, dtype=float)
        if override_y is not None
        else data["xs"].to_numpy(float)
    )
    if statistical_only:
        err = data["stat"].to_numpy(float)
    else:
        err = np.sqrt(
            data["stat"].to_numpy(float)**2
            + data["sys"].to_numpy(float)**2
            + (float(bh_systematic_fraction)
               * data["xs"].to_numpy(float))**2
        )
    #endif
    err = np.maximum(err, 1.0e-15)

    A = data["bh_A"].to_numpy(float)
    B = data["bh_B"].to_numpy(float)
    C = data["bh_C"].to_numpy(float)
    q_powers = np.vstack([q**i for i in range(1, nshape + 1)])

    def chi2(*values):
        coeffs = np.asarray(values, dtype=float)
        gm = MU_N * sachs_family_value_precomputed(
            q, coeffs, family, q_powers
        )
        f1 = (ge + tau * gm) * inv
        f2 = (gm - ge) * inv
        pred = bh_from_f1f2(A, B, C, f1, f2)
        pull = (pred - y) / err
        return float(np.dot(pull, pull))
    #enddef

    m = Minuit(chi2, *p0, name=tuple(names))
    m.errordef = Minuit.LEAST_SQUARES
    c_lo = sachs_first_coefficient_from_radius(
        SACHS_MIN_RADIUS_FM, family
    )
    c_hi = sachs_first_coefficient_from_radius(
        SACHS_MAX_RADIUS_FM, family
    )
    m.limits[names[0]] = (min(c_lo, c_hi), max(c_lo, c_hi))
    m.migrad()
    try:
        m.hesse()
    except Exception:
        pass
    #endtry

    pars = np.asarray([float(m.values[n]) for n in names], dtype=float)
    radius = sachs_family_radius(pars, family)
    radius_err = np.nan
    if return_uncertainty and m.valid:
        try:
            cov = np.asarray(m.covariance, dtype=float)
            grad = np.zeros(nshape, dtype=float)
            for ip in range(nshape):
                step = max(1.0e-6, 1.0e-5 * max(1.0, abs(pars[ip])))
                pp = pars.copy(); pm = pars.copy()
                pp[ip] += step; pm[ip] -= step
                grad[ip] = (
                    sachs_family_radius(pp, family)
                    - sachs_family_radius(pm, family)
                ) / (2.0 * step)
            #endfor
            var_r = float(grad @ cov @ grad)
            if np.isfinite(var_r) and var_r >= 0.0:
                radius_err = math.sqrt(var_r)
            #endif
        except Exception:
            radius_err = np.nan
        #endtry
    #endif
    valid = bool(m.valid and np.isfinite(radius))

    # Reject numerically converged but obviously pathological magnetic shapes
    # over the measured neutron |t| interval.
    shape = sachs_family_value(q, pars, family)
    if (
        np.any(~np.isfinite(shape))
        or np.any(shape <= 0.0)
        or np.any(shape > 1.15)
    ):
        valid = False
    #endif

    if return_parameters and return_uncertainty:
        return radius, radius_err, valid, pars, float(m.fval)
    #endif
    if return_parameters:
        return radius, valid, pars, float(m.fval)
    #endif
    if return_uncertainty:
        return radius, radius_err, valid
    #endif
    return radius, valid
#enddef



def _neutron_closure_family_batch_worker(task):
    """Fit one neutron GM family to a batch of pseudodata seeds."""
    evaluated, family, central, sigma, seeds, target = task
    successes = []
    attempts = []
    for attempted, seed in enumerate(seeds):
        rng = np.random.default_rng(int(seed))
        pseudo = np.asarray(central, dtype=float) + rng.normal(
            0.0, np.asarray(sigma, dtype=float)
        )
        try:
            radius, valid, pars, chi2 = fit_neutron_magnetic_family_bh(
                evaluated, family, override_y=pseudo,
                statistical_only=True, return_parameters=True,
            )
        except Exception:
            radius, valid, chi2 = np.nan, False, np.nan
        #endtry
        accepted = bool(valid and len(successes) < int(target))
        attempts.append((attempted, float(radius), float(chi2), bool(valid), accepted))
        if accepted:
            successes.append(float(radius))
        #endif
        if len(successes) >= int(target):
            break
        #endif
    #endfor
    return family, successes, attempts
#enddef


def run_neutron_magnetic_function_closure(
        evaluated: pd.DataFrame,
        args,
        outdir: Path) -> Optional[str]:
    """
    Select the neutron GM_n functional family by the same closure principle
    used for the proton analysis: generate pseudodata at the exact measured
    kinematics, fit competing P/IP/CF families through the exact BH quadratic,
    and rank by sqrt(variance^2 + bias^2) across a broad truth ensemble.

    GE_n is fixed to Atac et al. 2021 in both generation and fitting, matching
    the intended one-unknown production extraction.

    --radius-bias-replicas is interpreted as SUCCESSFUL fits per
    truth/family. Failed fits are retried with fresh seeds up to a finite cap.
    """
    closuredir = outdir / "magnetic_only" / "function_closure"
    closuredir.mkdir(parents=True, exist_ok=True)

    # The neutron closure is expensive and unchanged by plotting-only updates.
    # Reuse a completed ranking when available; delete this file explicitly if
    # a fresh closure calculation is desired.
    cached_ranking_path = closuredir / "03_neutron_GM_closure_ranking.csv"
    if cached_ranking_path.exists():
        cached = pd.read_csv(cached_ranking_path)
        picked = cached.loc[cached.get("selected", False).astype(bool)]
        if len(picked) == 1:
            chosen = str(picked.iloc[0]["family"])
            print(
                f"[neutron-GM-closure] reusing existing completed ranking; "
                f"selected family={chosen}"
            )
            return chosen
        #endif
    #endif

    candidate_families = [
        "D1",
        "P2", "P3",
        "IP1", "IP2", "IP3",
        "CF2", "CF3",
    ]
    truth_families = ["D1"]
    truth_families += [f"P{i}" for i in range(1, 5)]
    truth_families += [f"IP{i}" for i in range(1, 5)]
    truth_families += [f"CF{i}" for i in range(2, 5)]

    q = evaluated["t_abs"].to_numpy(float)
    qmax = float(np.nanmax(q))
    tau = q / (4.0 * MP2)
    inv = 1.0 / (1.0 + tau)
    ge_fixed = neutron_atac_ge(q, 0.505, 1.655, 0.909)
    A = evaluated["bh_A"].to_numpy(float)
    B = evaluated["bh_B"].to_numpy(float)
    C = evaluated["bh_C"].to_numpy(float)
    sigma = np.maximum(evaluated["stat"].to_numpy(float), 1.0e-15)

    # Empirical Kelly neutron GM truth.
    gm_kelly = (
        evaluated["km15_F1"].to_numpy(float)
        + evaluated["km15_F2"].to_numpy(float)
    )
    gm_kelly_norm = gm_kelly / MU_N

    def kelly_norm_shape(qq):
        qq = np.asarray(qq, dtype=float)
        # Interpolate only for radius evaluation/template construction below;
        # exact measured-point Kelly values are retained for pseudodata.
        q_aug = np.concatenate([[0.0], q])
        g_aug = np.concatenate([[1.0], gm_kelly_norm])
        order = np.argsort(q_aug)
        return np.interp(qq, q_aug[order], g_aug[order])
    #enddef

    # Estimate the Kelly truth radius from a very-low-Q2 finite difference of
    # the script's native Kelly Sachs implementation when available.  This
    # avoids inferring the Q2=0 slope from the neutron data range.
    h = 1.0e-5
    try:
        _, gm0_pair = kelly_sachs(np.asarray([0.0, h], dtype=float))
        kelly_rM = radius_from_shape(
            lambda qq: kelly_sachs(np.asarray(qq, dtype=float))[1],
            MU_P,
        )
        # kelly_sachs is proton-specific; do not use its radius for neutron.
        # The branch above is intentionally discarded after confirming that
        # this helper is not the neutron Kelly implementation.
        _ = gm0_pair, kelly_rM
        kelly_rM = np.nan
    except Exception:
        kelly_rM = np.nan
    #endtry

    # Fit a high-order smooth template directly to the measured neutron Kelly
    # GM values plus the exact normalization point. Its slope supplies a
    # reproducible empirical-Kelly reference radius for closure bookkeeping.
    q_template_fit = np.concatenate([[0.0], q])
    g_template_fit = np.concatenate([[1.0], gm_kelly_norm])
    coeff_kelly = np.polyfit(q_template_fit, g_template_fit, deg=4)
    deriv0 = float(np.polyval(np.polyder(coeff_kelly), 0.0))
    if deriv0 < 0.0:
        kelly_rM = HBARC * math.sqrt(-6.0 * deriv0)
    #endif

    scenarios = [{
        "truth_model": "Kelly_neutron_measured_point_template",
        "truth_group": "Kelly",
        "truth_family": "empirical",
        "truth_radius_fm": float(kelly_rM),
        "truth_shape_measured": gm_kelly_norm.copy(),
        "synthetic": False,
    }]

    # Broad synthetic truth ensemble.  As in the proton study, higher-order
    # P4/IP4/CF4 shapes remain in truth generation even though they are not
    # production candidates.
    radius_values = parse_radius_bias_grid(args.neutron_radius_bias_radius_grid)
    q_template = np.linspace(0.0, max(qmax, 0.05), 400)

    # Use a smooth IP3 approximation to the measured Kelly neutron GM as the
    # curvature template from which all generating families are constructed.
    template_base = fit_sachs_family_shape_template(
        "IP3",
        q_template_fit,
        g_template_fit,
    )
    smooth_kelly = sachs_family_value(q_template, template_base, "IP3")

    templates = {}
    for family in truth_families:
        templates[family] = fit_sachs_family_shape_template(
            family, q_template, smooth_kelly
        )
    #endfor

    skipped = 0
    for family in truth_families:
        for radius in radius_values:
            coeffs = sachs_family_coefficients_with_radius(
                templates[family], family, radius
            )
            shape_grid = sachs_family_value(q_template, coeffs, family)
            if (
                np.any(~np.isfinite(shape_grid))
                or np.any(shape_grid <= 0.0)
                or np.any(shape_grid > 1.15)
            ):
                skipped += 1
                continue
            #endif
            scenarios.append({
                "truth_model": f"synthetic_{family}_rM{radius:.3f}",
                "truth_group": family,
                "truth_family": family,
                "truth_radius_fm": float(radius),
                "truth_shape_measured": sachs_family_value(
                    q, coeffs, family
                ),
                "synthetic": True,
            })
        #endfor
    #endfor

    print(
        f"[neutron-GM-closure] truths={len(scenarios)} "
        f"({skipped} synthetic rejected), "
        f"families={len(candidate_families)}, "
        f"successful replicas/truth/family={args.radius_bias_replicas}"
    )

    target = max(1, int(args.radius_bias_replicas))
    min_eff = max(0.50, min(0.99, float(
        getattr(args, "radius_bias_min_scenario_valid_fraction", 0.80)
    )))
    max_attempts = max(target + 10, int(math.ceil(1.10 * target / min_eff)))
    rows = []
    replica_rows = []

    nworkers = getattr(args, "radius_bias_workers", None)
    nworkers = int(nworkers if nworkers is not None else args.workers)
    nworkers = max(1, min(nworkers, len(candidate_families)))

    with ProcessPoolExecutor(max_workers=nworkers) as neutron_pool:
        for itruth, scenario in enumerate(scenarios):
            gm_truth = MU_N * np.asarray(
                scenario["truth_shape_measured"], dtype=float
            )
            f1_truth = (ge_fixed + tau * gm_truth) * inv
            f2_truth = (gm_truth - ge_fixed) * inv
            central = bh_from_f1f2(A, B, C, f1_truth, f2_truth)
            rtrue = float(scenario["truth_radius_fm"])

            tasks = []
            for ifam, family in enumerate(candidate_families):
                seed_seq = np.random.SeedSequence([
                    int(args.radius_bias_seed), 918273, int(itruth), int(ifam)
                ])
                seed_states = seed_seq.spawn(max_attempts)
                seeds = [
                    int(ss.generate_state(1, dtype=np.uint64)[0])
                    for ss in seed_states
                ]
                tasks.append((evaluated, family, central, sigma, seeds, target))
            #endfor

            for family, successes, attempts in neutron_pool.map(
                    _neutron_closure_family_batch_worker, tasks, chunksize=1):
                for attempt_index, radius, chi2, valid, accepted in attempts:
                    replica_rows.append({
                        "truth_model": scenario["truth_model"],
                        "truth_group": scenario["truth_group"],
                        "truth_family": scenario["truth_family"],
                        "synthetic_truth": bool(scenario["synthetic"]),
                        "truth_rM_fm": rtrue,
                        "family": family,
                        "attempt": attempt_index,
                        "accepted_replica": (
                            sum(1 for a in attempts[:attempt_index + 1] if a[4]) - 1
                            if accepted else np.nan
                        ),
                        "rM_fit_fm": radius,
                        "chi2": chi2,
                        "valid": bool(valid),
                    })
                #endfor

                arr = np.asarray(successes, dtype=float)
                if len(arr):
                    mean = float(np.mean(arr))
                    stat = float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0
                    bias = float(mean - rtrue)
                    objective = float(math.sqrt(stat**2 + bias**2))
                else:
                    mean = stat = bias = objective = np.nan
                #endif
                attempted_count = len(attempts)
                rows.append({
                    "truth_model": scenario["truth_model"],
                    "truth_group": scenario["truth_group"],
                    "truth_family": scenario["truth_family"],
                    "synthetic_truth": bool(scenario["synthetic"]),
                    "truth_rM_fm": rtrue,
                    "family": family,
                    "mean_extracted_rM_fm": mean,
                    "stat_RMS_fm": stat,
                    "bias_fm": bias,
                    "sqrt_stat2_plus_bias2_fm": objective,
                    "valid_replicas": int(len(arr)),
                    "requested_replicas": target,
                    "attempted_replicas": attempted_count,
                    "attempt_efficiency": (
                        len(arr) / attempted_count if attempted_count else np.nan
                    ),
                    "target_reached": bool(len(arr) >= target),
                })
                print(
                    f"[neutron-GM-closure] "
                    f"truth={str(scenario['truth_model']):28s} "
                    f"family={family:3s} valid={len(arr)}/{target} "
                    f"attempted={attempted_count}"
                )
            #endfor
        #endfor
    #endwith

    table = pd.DataFrame(rows)
    replicas = pd.DataFrame(replica_rows)
    table.to_csv(
        closuredir / "01_neutron_GM_closure_by_truth.csv", index=False
    )
    replicas.to_csv(
        closuredir / "02_neutron_GM_closure_replicas.csv", index=False
    )

    ranking_rows = []
    for family in candidate_families:
        part = table.loc[table["family"] == family].copy()
        obj = part["sqrt_stat2_plus_bias2_fm"].to_numpy(float)
        finite = np.isfinite(obj)
        target_fraction = float(np.mean(part["target_reached"].astype(float)))
        ranking_rows.append({
            "family": family,
            "RMS_objective_across_truths_fm": (
                float(np.sqrt(np.mean(obj[finite]**2)))
                if np.any(finite) else np.nan
            ),
            "max_objective_across_truths_fm": (
                float(np.max(obj[finite])) if np.any(finite) else np.nan
            ),
            "mean_abs_bias_fm": float(
                np.nanmean(np.abs(part["bias_fm"].to_numpy(float)))
            ),
            "mean_stat_RMS_fm": float(
                np.nanmean(part["stat_RMS_fm"].to_numpy(float))
            ),
            "truth_scenarios": int(len(part)),
            "target_reached_fraction": target_fraction,
            "mean_attempt_efficiency": float(
                np.nanmean(part["attempt_efficiency"].to_numpy(float))
            ),
        })
    #endfor

    ranking = pd.DataFrame(ranking_rows)
    eligible = ranking.loc[
        ranking["target_reached_fraction"] >= 0.999999
    ].copy()
    if not len(eligible):
        print(
            "[neutron-GM-closure] WARNING: no family reached the requested "
            "successful-replica target for every truth scenario"
        )
    #endif

    # Second gate: a closure-eligible family must also give a physically valid
    # fit to the actual neutron data.  Audit every candidate on all evaluated
    # points so a good synthetic-closure score cannot hide a pathological
    # nominal-data extrapolation.
    nominal_rows = []
    for family in candidate_families:
        npars = int(re.findall(r"\d+", family)[0])
        radius, valid, pars, chi2 = fit_neutron_magnetic_family_bh(
            evaluated,
            family,
            statistical_only=False,
            bh_systematic_fraction=0.05,
            return_parameters=True,
        )
        ndof = max(1, len(evaluated) - npars)
        rank_row = ranking.loc[ranking["family"] == family].iloc[0]
        closure_eligible = bool(
            float(rank_row["target_reached_fraction"]) >= 0.999999
        )
        row = {
            "family": family,
            "closure_eligible": closure_eligible,
            "closure_RMS_objective_fm": float(
                rank_row["RMS_objective_across_truths_fm"]
            ),
            "closure_max_objective_fm": float(
                rank_row["max_objective_across_truths_fm"]
            ),
            "nominal_valid": bool(valid),
            "nominal_rM_fm": float(radius) if valid else np.nan,
            "nominal_chi2": float(chi2),
            "nominal_ndof": int(ndof),
            "nominal_chi2_ndof": float(chi2 / ndof),
        }
        for i, value in enumerate(pars, start=1):
            row[f"m{i}"] = float(value)
        #endfor
        nominal_rows.append(row)
        print(
            f"[neutron-GM-nominal-audit] family={family:3s} "
            f"closure={'PASS' if closure_eligible else 'FAIL'} "
            f"nominal={'PASS' if valid else 'FAIL'} "
            f"chi2/ndof={chi2/ndof:.3f} "
            + (
                f"rM={radius:.4f} fm"
                if valid else
                "rM=UNRESOLVED"
            )
        )
    #endfor

    nominal_audit = pd.DataFrame(nominal_rows)
    nominal_audit.to_csv(
        closuredir / "03b_neutron_GM_nominal_family_audit.csv",
        index=False,
    )

    # Select in closure-rank order, but only among families passing both gates.
    ranked_eligible = nominal_audit.loc[
        nominal_audit["closure_eligible"]
        & nominal_audit["nominal_valid"]
    ].sort_values(
        ["closure_RMS_objective_fm", "closure_max_objective_fm"]
    )
    if not len(ranked_eligible):
        chosen = None
        print(
            "[neutron-GM-closure] WARNING: no family passes both closure "
            "completion and nominal-data physical-validity gates"
        )
    else:
        chosen = str(ranked_eligible.iloc[0]["family"])
    #endif

    ranking["nominal_valid"] = ranking["family"].map(
        nominal_audit.set_index("family")["nominal_valid"]
    )
    ranking["passes_both_gates"] = (
        (ranking["target_reached_fraction"] >= 0.999999)
        & ranking["nominal_valid"].fillna(False)
    )
    ranking["selected"] = ranking["family"].eq(chosen)
    ranking.to_csv(
        closuredir / "03_neutron_GM_closure_ranking.csv", index=False
    )

    fig, ax = plt.subplots(figsize=(9.0, 5.4))
    ordered = ranking.sort_values(
        "RMS_objective_across_truths_fm"
    ).reset_index(drop=True)
    ax.plot(
        np.arange(len(ordered)),
        ordered["RMS_objective_across_truths_fm"],
        marker="o", linestyle="none",
    )
    ax.set_xticks(np.arange(len(ordered)))
    ax.set_xticklabels(ordered["family"])
    ax.set_ylabel(
        r"RMS across truths of $\sqrt{\sigma_{\rm stat}^2+b^2}$ [fm]"
    )
    ax.set_xlabel(r"$G_M^n/\mu_n$ fit family")
    ax.set_title(
        r"Neutron $G_M^n$ functional-form closure ranking"
    )
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(
        closuredir / "04_neutron_GM_closure_ranking.png", dpi=180
    )
    plt.close(fig)

    # Bias and variance are shown separately so the closure objective is not
    # mistaken for a systematic uncertainty.
    fig, ax = plt.subplots(figsize=(9.0, 5.4))
    x = np.arange(len(ordered))
    ax.plot(
        x, ordered["mean_abs_bias_fm"],
        marker="o", linestyle="none", label="mean |bias|",
    )
    ax.plot(
        x, ordered["mean_stat_RMS_fm"],
        marker="s", linestyle="none", fillstyle="none",
        label="mean replica RMS",
    )
    ax.set_xticks(x)
    ax.set_xticklabels(ordered["family"])
    ax.set_ylabel("Radius scale [fm]")
    ax.set_xlabel(r"$G_M^n/\mu_n$ fit family")
    ax.set_title(r"Neutron $r_{M,n}$ closure: bias and variance")
    ax.grid(alpha=0.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(
        closuredir / "05_neutron_GM_closure_bias_variance.png", dpi=180
    )
    plt.close(fig)

    pd.DataFrame([{
        "selected_family": chosen if chosen is not None else "UNRESOLVED",
        "selection_rule": (
            "best closure-ranked family passing both successful-replica "
            "completion and nominal-data physical-validity gates"
        ),
        "fixed_GE_model": "Atac2021 central",
        "candidate_families": ",".join(candidate_families),
        "truth_families": ",".join(truth_families),
        "radius_grid_fm": str(args.neutron_radius_bias_radius_grid),
        "successful_replicas_per_truth_family": target,
    }]).to_csv(
        closuredir / "06_neutron_GM_selected_family.csv", index=False
    )

    print(
        f"[neutron-GM-closure] selected family = "
        f"{chosen if chosen is not None else 'UNRESOLVED'}"
    )
    return chosen
#enddef


def save_neutron_magnetic_selected_family_study(
        evaluated: pd.DataFrame,
        family: Optional[str],
        outdir: Path) -> None:
    """Run the neutron magnetic-radius threshold scan with the closure-selected family."""
    diagdir = outdir / "magnetic_only"
    if family is None:
        return
    #endif

    finite_delta = evaluated.loc[
        np.isfinite(evaluated["bh_delta"]), "bh_delta"
    ].to_numpy(float)
    thresholds = [
        0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.10,
        0.15, 0.20, 0.30, 0.50, 0.75, 1.00,
    ]
    if len(finite_delta):
        thresholds.append(float(np.max(finite_delta) * 1.001))
    #endif
    thresholds = sorted(set(thresholds))

    rows = []
    npars = int(re.findall(r"\d+", family)[0])
    for threshold in thresholds:
        selected = evaluated.loc[
            np.isfinite(evaluated["bh_delta"])
            & (evaluated["bh_delta"] <= float(threshold))
        ].copy()
        if len(selected) <= npars:
            continue
        #endif
        radius, radius_err, valid, pars, chi2 = fit_neutron_magnetic_family_bh(
            selected,
            family,
            statistical_only=False,
            bh_systematic_fraction=0.05,
            return_parameters=True,
            return_uncertainty=True,
        )
        ndof = max(1, len(selected) - npars)
        row = {
            "family": family,
            "bh_threshold": float(threshold),
            "bh_threshold_percent": 100.0 * float(threshold),
            "N": int(len(selected)),
            "n_parameters": npars,
            # Never expose a radius from a fit that failed the physical-shape
            # validity gate. Keep the raw fit coefficients/chi2 below for
            # diagnosis, but mark the radius unresolved.
            "rM_fm": float(radius) if valid else np.nan,
            "rM_fit_err_fm": float(radius_err) if valid else np.nan,
            "chi2": float(chi2),
            "ndof": int(ndof),
            "chi2_ndof": float(chi2 / ndof),
            "valid": bool(valid),
        }
        for i, value in enumerate(pars, start=1):
            row[f"m{i}"] = float(value)
        #endfor
        rows.append(row)
    #endfor

    scan = pd.DataFrame(rows)
    scan.to_csv(
        diagdir / "05_neutron_selected_family_threshold_scan.csv",
        index=False,
    )
    if not len(scan):
        return
    #endif

    valid_scan = scan.loc[
        scan["valid"] & np.isfinite(scan["rM_fm"])
    ].copy()
    if not len(valid_scan):
        print(
            f"[neutron-GM-selected] family={family}: no BH-threshold fit "
            "passed the physical-validity gate; radius plot suppressed"
        )
        return
    #endif

    fig, ax = plt.subplots(figsize=(8.4, 5.6))
    ax.errorbar(
        valid_scan["bh_threshold_percent"],
        valid_scan["rM_fm"],
        yerr=valid_scan["rM_fit_err_fm"],
        marker="o", linestyle="none", capsize=3, linewidth=1.0,
    )
    ax.axhline(
        0.864, linewidth=1.0, linestyle="--",
        label=r"PDG $r_{M,n}=0.864$ fm",
    )
    ax.axvline(5.0, linewidth=0.8, linestyle=":")
    ax.set_xscale("log")
    ax.set_xlabel(r"$|1-R_{\rm BH}^{\rm KM15}|$ threshold (%)")
    ax.set_ylabel(r"$r_{M,n}$ [fm]")
    ax.set_ylim(0.70, 1.00)
    ax.set_title(
        rf"Neutron magnetic radius with closure-selected {family}"
    )
    ax.grid(alpha=0.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(
        diagdir / "05_neutron_selected_family_radius_threshold.png",
        dpi=300,
    )
    plt.close(fig)
#enddef



def save_neutron_phi_shape_bh_compatibility_study(
        evaluated: pd.DataFrame,
        outdir: Path) -> None:
    """
    Test BH compatibility at the level of each complete eight-point phi
    distribution rather than by selecting individual cross-section points.

    In each (Q2,xB,t) bin, GE_n is fixed to the Atac et al. central form and a
    single local magnetic amplitude gM_scale is fitted:
        GM_n(Q2) = gM_scale * GM_n^Kelly(Q2).

    The fit therefore lets the magnetic normalization appropriate to that
    finite-|t| bin float while testing whether the measured phi dependence is
    compatible with the exact BH angular structure.  With eight phi points and
    one fitted amplitude, the shape test has seven degrees of freedom.

    This is a diagnostic of non-BH azimuthal structure.  It does not use KM15
    EP/BH purity to select points and it does not extrapolate to Q2=0.
    """
    diagdir = outdir / "phi_shape_bh_compatibility"
    diagdir.mkdir(parents=True, exist_ok=True)

    rows = []
    predictions = []

    for kin_bin, part in evaluated.groupby("kin_bin", sort=True):
        p = part.sort_values("phi_deg").copy()
        q = p["t_abs"].to_numpy(float)
        tau = q / (4.0 * MP2)
        inv = 1.0 / (1.0 + tau)
        ge = neutron_atac_ge(q, 0.505, 1.655, 0.909)
        # Use the same pointwise neutron Kelly/Gepard Sachs GM already
        # carried by the evaluated neutron table: GM = F1 + F2.
        gm_ref = (
            p["km15_F1"].to_numpy(float)
            + p["km15_F2"].to_numpy(float)
        )
        y = p["xs"].to_numpy(float)
        err = np.sqrt(
            p["stat"].to_numpy(float)**2 + p["sys"].to_numpy(float)**2
        )
        err = np.maximum(err, 1.0e-15)
        A = p["bh_A"].to_numpy(float)
        B = p["bh_B"].to_numpy(float)
        C = p["bh_C"].to_numpy(float)

        def chi2(gm_scale):
            gm = float(gm_scale) * gm_ref
            f1 = (ge + tau * gm) * inv
            f2 = (gm - ge) * inv
            pred = bh_from_f1f2(A, B, C, f1, f2)
            pull = (pred - y) / err
            return float(np.dot(pull, pull))
        #enddef

        m = Minuit(chi2, gm_scale=1.0)
        m.errordef = Minuit.LEAST_SQUARES
        m.limits["gm_scale"] = (0.20, 1.80)
        m.migrad()

        gm_scale = float(m.values["gm_scale"])
        gm_scale_err = float(m.errors["gm_scale"])
        chi2_min = float(m.fval)
        ndof = max(1, len(p) - 1)
        chi2_ndof = chi2_min / ndof

        # Survival probability for integer ndof without introducing scipy:
        # Q(nu/2, chi2/2).  Here nu=7 for the complete 8-point bins, but use
        # mpmath only if it is already available would add a dependency.
        # Store chi2/ndof as the primary compatibility statistic; the familiar
        # 95% critical value is supplied below for the expected 7 dof.
        chi2_95 = 14.067 if ndof == 7 else np.nan
        compatible_95 = bool(
            np.isfinite(chi2_95) and chi2_min <= chi2_95
        )

        gm = gm_scale * gm_ref
        f1 = (ge + tau * gm) * inv
        f2 = (gm - ge) * inv
        pred = bh_from_f1f2(A, B, C, f1, f2)
        pulls = (y - pred) / err

        rows.append({
            "kin_bin": int(kin_bin),
            "N_phi": int(len(p)),
            "Q2_mean_GeV2": float(np.nanmean(p["Q2"])),
            "xB_mean": float(np.nanmean(p["xB"])),
            "t_abs_mean_GeV2": float(np.nanmean(p["t_abs"])),
            "gm_scale_vs_Kelly": gm_scale,
            "gm_scale_err": gm_scale_err,
            "chi2": chi2_min,
            "ndof": int(ndof),
            "chi2_ndof": chi2_ndof,
            "chi2_95_critical": chi2_95,
            "BH_shape_compatible_95pct": compatible_95,
            "median_abs_pull": float(np.nanmedian(np.abs(pulls))),
            "max_abs_pull": float(np.nanmax(np.abs(pulls))),
            "median_km15_bh_delta": float(np.nanmedian(p["bh_delta"])),
            "max_km15_bh_delta": float(np.nanmax(p["bh_delta"])),
        })

        for j, (_, point) in enumerate(p.iterrows()):
            predictions.append({
                "kin_bin": int(kin_bin),
                "phi_deg": float(point["phi_deg"]),
                "Q2_GeV2": float(point["Q2"]),
                "xB": float(point["xB"]),
                "t_abs_GeV2": float(point["t_abs"]),
                "xs_data": float(point["xs"]),
                "xs_err_total": float(err[j]),
                "xs_BH_local_GM_fit": float(pred[j]),
                "pull_BH_local_GM_fit": float(pulls[j]),
                "gm_scale_vs_Kelly": gm_scale,
                "km15_bh_delta": float(point["bh_delta"]),
            })
        #endfor
    #endfor

    summary = pd.DataFrame(rows).sort_values("kin_bin")
    points = pd.DataFrame(predictions).sort_values(["kin_bin", "phi_deg"])
    summary.to_csv(
        diagdir / "01_neutron_phi_shape_BH_compatibility_by_bin.csv",
        index=False,
    )
    points.to_csv(
        diagdir / "02_neutron_phi_shape_BH_fit_points.csv",
        index=False,
    )

    # Main 3x3 observable-space diagnostic: data versus the best local-GM BH
    # shape in each independent kinematic bin.
    fig, axes = plt.subplots(
        3, 3, figsize=(15.0, 11.5), sharex=True
    )
    axes = np.asarray(axes).ravel()
    for ax, (_, row) in zip(axes, summary.iterrows()):
        kin_bin = int(row["kin_bin"])
        pp = points.loc[points["kin_bin"] == kin_bin].sort_values("phi_deg")
        ax.errorbar(
            pp["phi_deg"], pp["xs_data"], yerr=pp["xs_err_total"],
            fmt="o", linestyle="none", label="Data",
        )
        ax.plot(
            pp["phi_deg"], pp["xs_BH_local_GM_fit"],
            marker="s", fillstyle="none", linestyle="none",
            label=r"BH, local $G_M^n$ fit",
        )
        verdict = "compatible" if row["BH_shape_compatible_95pct"] else "rejected"
        ax.set_title(
            rf"Bin {kin_bin}: $\chi^2/\nu={row['chi2_ndof']:.2f}$, "
            rf"$\nu={int(row['ndof'])}$ ({verdict})"
        )
        ax.grid(alpha=0.2)
    #endfor
    for ax in axes[6:]:
        ax.set_xlabel(r"$\phi$ [deg]")
    #endfor
    for ax in axes[::3]:
        ax.set_ylabel(r"$d^4\sigma$")
    #endfor
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles, labels, loc="upper center", ncol=2,
        bbox_to_anchor=(0.5, 0.965),
    )
    fig.suptitle(
        r"Neutron BH $\phi$-shape compatibility with local $G_M^n$ floated",
        y=0.995, fontsize=16,
    )
    fig.subplots_adjust(top=0.90, hspace=0.30, wspace=0.24)
    fig.savefig(
        diagdir / "03_neutron_phi_shape_BH_compatibility_3x3.png",
        dpi=180,
    )
    plt.close(fig)

    # Per-bin pulls make coherent non-BH angular residuals much easier to see.
    fig, axes = plt.subplots(
        3, 3, figsize=(15.0, 11.5), sharex=True, sharey=True
    )
    axes = np.asarray(axes).ravel()
    for ax, (_, row) in zip(axes, summary.iterrows()):
        kin_bin = int(row["kin_bin"])
        pp = points.loc[points["kin_bin"] == kin_bin].sort_values("phi_deg")
        ax.axhline(0.0, linewidth=0.8)
        ax.axhline(2.0, linewidth=0.6, linestyle=":")
        ax.axhline(-2.0, linewidth=0.6, linestyle=":")
        ax.plot(
            pp["phi_deg"], pp["pull_BH_local_GM_fit"],
            marker="o", linestyle="none",
        )
        ax.set_title(
            rf"Bin {kin_bin}: max $|z|={row['max_abs_pull']:.2f}$"
        )
        ax.grid(alpha=0.2)
    #endfor
    for ax in axes[6:]:
        ax.set_xlabel(r"$\phi$ [deg]")
    #endfor
    for ax in axes[::3]:
        ax.set_ylabel("BH-fit pull")
    #endfor
    fig.suptitle(
        r"Residual angular structure after local-$G_M^n$ BH fit",
        y=0.995, fontsize=16,
    )
    fig.subplots_adjust(top=0.92, hspace=0.28, wspace=0.18)
    fig.savefig(
        diagdir / "04_neutron_phi_shape_BH_pull_3x3.png",
        dpi=180,
    )
    plt.close(fig)

    # Compare whole-bin BH-shape compatibility with KM15's point-level purity
    # expectation, without using either quantity to select the other.
    fig, ax = plt.subplots(figsize=(7.2, 5.5))
    ax.scatter(
        100.0 * summary["median_km15_bh_delta"],
        summary["chi2_ndof"],
        s=42,
    )
    for _, row in summary.iterrows():
        ax.annotate(
            str(int(row["kin_bin"])),
            (100.0 * row["median_km15_bh_delta"], row["chi2_ndof"]),
            xytext=(4, 4), textcoords="offset points", fontsize=8,
        )
    #endfor
    ax.axhline(14.067 / 7.0, linewidth=0.8, linestyle="--")
    ax.set_xlabel(
        r"Median KM15 $|1-\sigma_{\rm BH}/\sigma_{\rm EP}|$ in bin (%)"
    )
    ax.set_ylabel(r"Local-$G_M^n$ BH shape $\chi^2/\nu$")
    ax.set_title("Whole-bin BH shape test vs KM15 purity expectation")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(
        diagdir / "05_neutron_phi_shape_chi2_vs_km15_purity.png",
        dpi=180,
    )
    plt.close(fig)

    n_ok = int(np.sum(summary["BH_shape_compatible_95pct"]))
    print(
        f"[neutron-phi-shape-BH] {n_ok}/{len(summary)} kinematic bins "
        "compatible with a local-GM BH phi shape at the 95% level"
    )
    print(f"[neutron-phi-shape-BH] outputs -> {diagdir}")
#enddef


def save_neutron_empirical_bh_compatibility_study(
        evaluated: pd.DataFrame,
        family: Optional[str],
        outdir: Path) -> None:
    """
    Diagnostic study of empirical compatibility with the reference BH
    prediction.

    This is deliberately NOT used to define the production sample because the
    reference BH calculation contains an assumed neutron form factor.  It asks
    whether the measured cross section itself looks BH-like and compares that
    empirical statement with the KM15 model-purity classifier.

    Two observed-compatibility variables are used:
        delta_obs = |1 - sigma_BH_ref / sigma_data|
        pull_BH   = (sigma_data - sigma_BH_ref) / sqrt(stat^2 + sys^2)

    The closure-selected GM family is then refitted after progressively tighter
    cuts on each observed variable.  Those radius scans are diagnostics of the
    KM15 purity criterion, not alternate production extractions.
    """
    diagdir = outdir / "empirical_bh_compatibility"
    diagdir.mkdir(parents=True, exist_ok=True)

    work = evaluated.copy()
    y = work["xs"].to_numpy(float)
    bh = work["km15_bh"].to_numpy(float)
    exp_err = np.sqrt(
        work["stat"].to_numpy(float)**2
        + work["sys"].to_numpy(float)**2
    )
    safe_y = np.maximum(np.abs(y), 1.0e-15)
    safe_err = np.maximum(exp_err, 1.0e-15)

    work["bh_signed_fractional_residual"] = (y - bh) / safe_y
    work["bh_abs_fractional_residual"] = np.abs(
        work["bh_signed_fractional_residual"].to_numpy(float)
    )
    work["bh_pull_exp"] = (y - bh) / safe_err
    work["bh_abs_pull_exp"] = np.abs(work["bh_pull_exp"].to_numpy(float))
    work["data_total_err"] = exp_err
    work["km15_bh_purity_percent"] = 100.0 * work["bh_delta"]
    work["observed_bh_difference_percent"] = (
        100.0 * work["bh_abs_fractional_residual"]
    )
    work.to_csv(
        diagdir / "01_neutron_empirical_bh_compatibility_points.csv",
        index=False,
    )

    # Direct comparison between the model-predicted non-BH importance and what
    # the data actually do relative to the reference BH calculation.
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    ax.scatter(
        100.0 * work["bh_delta"],
        100.0 * work["bh_abs_fractional_residual"],
        s=30,
    )
    lim = max(
        5.0,
        float(np.nanmax(100.0 * work["bh_delta"])) * 1.08,
        float(np.nanpercentile(
            100.0 * work["bh_abs_fractional_residual"], 97
        )) * 1.08,
    )
    ax.plot([0.0, lim], [0.0, lim], linewidth=0.9, linestyle="--")
    ax.set_xlim(0.0, lim)
    ax.set_ylim(0.0, lim)
    ax.set_xlabel(
        r"KM15 predicted $|1-\sigma_{\rm BH}/\sigma_{\rm EP}|$ (%)"
    )
    ax.set_ylabel(
        r"Observed $|1-\sigma_{\rm BH}^{\rm ref}/\sigma_{\rm data}|$ (%)"
    )
    ax.set_title("Neutron: predicted non-BH importance vs observed BH difference")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(
        diagdir / "02_neutron_observed_vs_km15_bh_difference.png", dpi=180
    )
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    ax.scatter(
        100.0 * work["bh_delta"],
        work["bh_pull_exp"],
        s=30,
    )
    for level in [-3.0, -2.0, -1.0, 1.0, 2.0, 3.0]:
        ax.axhline(level, linewidth=0.7, linestyle=":")
    #endfor
    ax.axhline(0.0, linewidth=0.9)
    ax.set_xlabel(
        r"KM15 predicted $|1-\sigma_{\rm BH}/\sigma_{\rm EP}|$ (%)"
    )
    ax.set_ylabel(
        r"BH pull $(\sigma_{\rm data}-\sigma_{\rm BH}^{\rm ref})/"
        r"\sqrt{\delta_{\rm stat}^2+\delta_{\rm sys}^2}$"
    )
    ax.set_title("Neutron: empirical BH pull vs KM15 purity classifier")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(
        diagdir / "03_neutron_bh_pull_vs_km15_purity.png", dpi=180
    )
    plt.close(fig)

    # Summarize how often the two criteria agree/disagree.
    model_thresholds = [0.05, 0.10, 0.20, 0.30, 0.50, 0.75]
    observed_thresholds = [0.05, 0.10, 0.15, 0.20, 0.30, 0.50]
    agreement_rows = []
    for model_cut in model_thresholds:
        for obs_cut in observed_thresholds:
            model_bh = work["bh_delta"] <= model_cut
            observed_bh = work["bh_abs_fractional_residual"] <= obs_cut
            agreement_rows.append({
                "km15_delta_cut": model_cut,
                "km15_delta_cut_percent": 100.0 * model_cut,
                "observed_delta_cut": obs_cut,
                "observed_delta_cut_percent": 100.0 * obs_cut,
                "N_both_BH_like": int(np.sum(model_bh & observed_bh)),
                "N_KM15_only_BH_like": int(np.sum(model_bh & ~observed_bh)),
                "N_observed_only_BH_like": int(np.sum(~model_bh & observed_bh)),
                "N_neither_BH_like": int(np.sum(~model_bh & ~observed_bh)),
            })
        #endfor
    #endfor
    pd.DataFrame(agreement_rows).to_csv(
        diagdir / "04_neutron_km15_vs_observed_bh_agreement.csv",
        index=False,
    )

    if family is None:
        print(
            "[neutron-empirical-BH] no closure-selected GM family; "
            "compatibility plots saved but radius diagnostics skipped"
        )
        return
    #endif

    npars = int(re.findall(r"\d+", family)[0])
    scan_rows = []

    # Fractional central-value compatibility scan. This is intentionally
    # diagnostic because the reference BH contains an assumed neutron FF.
    for cut in observed_thresholds:
        selected = work.loc[
            np.isfinite(work["bh_abs_fractional_residual"])
            & (work["bh_abs_fractional_residual"] <= cut)
        ].copy()
        if len(selected) <= npars:
            continue
        #endif
        radius, valid, pars, chi2 = fit_neutron_magnetic_family_bh(
            selected,
            family,
            statistical_only=False,
            bh_systematic_fraction=0.05,
            return_parameters=True,
        )
        ndof = max(1, len(selected) - npars)
        row = {
            "selection_type": "observed_fractional_BH_difference",
            "selection_cut": float(cut),
            "selection_cut_display": 100.0 * float(cut),
            "selection_units": "percent",
            "family": family,
            "N": int(len(selected)),
            "rM_fm": float(radius) if valid else np.nan,
            "chi2": float(chi2),
            "ndof": int(ndof),
            "chi2_ndof": float(chi2 / ndof),
            "valid": bool(valid),
            "median_km15_bh_delta": float(np.nanmedian(selected["bh_delta"])),
            "max_km15_bh_delta": float(np.nanmax(selected["bh_delta"])),
        }
        for i, value in enumerate(pars, start=1):
            row[f"m{i}"] = float(value)
        #endfor
        scan_rows.append(row)
    #endfor

    # Statistical-compatibility scan. Unlike the fractional cut, this accounts
    # for the fact that a large central-value difference may be insignificant
    # for a low-precision point.
    pull_thresholds = [0.5, 1.0, 1.5, 2.0, 3.0]
    for cut in pull_thresholds:
        selected = work.loc[
            np.isfinite(work["bh_abs_pull_exp"])
            & (work["bh_abs_pull_exp"] <= cut)
        ].copy()
        if len(selected) <= npars:
            continue
        #endif
        radius, valid, pars, chi2 = fit_neutron_magnetic_family_bh(
            selected,
            family,
            statistical_only=False,
            bh_systematic_fraction=0.05,
            return_parameters=True,
        )
        ndof = max(1, len(selected) - npars)
        row = {
            "selection_type": "observed_BH_pull",
            "selection_cut": float(cut),
            "selection_cut_display": float(cut),
            "selection_units": "sigma",
            "family": family,
            "N": int(len(selected)),
            "rM_fm": float(radius) if valid else np.nan,
            "chi2": float(chi2),
            "ndof": int(ndof),
            "chi2_ndof": float(chi2 / ndof),
            "valid": bool(valid),
            "median_km15_bh_delta": float(np.nanmedian(selected["bh_delta"])),
            "max_km15_bh_delta": float(np.nanmax(selected["bh_delta"])),
        }
        for i, value in enumerate(pars, start=1):
            row[f"m{i}"] = float(value)
        #endfor
        scan_rows.append(row)
    #endfor

    scan = pd.DataFrame(scan_rows)
    scan.to_csv(
        diagdir / "05_neutron_empirical_bh_radius_diagnostic.csv",
        index=False,
    )

    # Radius vs observed fractional BH compatibility.
    frac = scan.loc[
        (scan["selection_type"] == "observed_fractional_BH_difference")
        & scan["valid"]
        & np.isfinite(scan["rM_fm"])
    ].copy()
    if len(frac):
        fig, ax = plt.subplots(figsize=(7.2, 5.4))
        ax.plot(
            frac["selection_cut_display"],
            frac["rM_fm"],
            marker="o", linestyle="none",
        )
        ax.axhline(
            0.864, linewidth=1.0, linestyle="--",
            label=r"PDG $r_{M,n}=0.864$ fm",
        )
        for _, row in frac.iterrows():
            ax.annotate(
                f"N={int(row['N'])}",
                (row["selection_cut_display"], row["rM_fm"]),
                xytext=(4, 4), textcoords="offset points", fontsize=8,
            )
        #endfor
        ax.set_xlabel(
            r"Observed $|1-\sigma_{\rm BH}^{\rm ref}/\sigma_{\rm data}|$ "
            r"cut (%)"
        )
        ax.set_ylabel(r"$r_{M,n}$ [fm]")
        ax.set_title(
            rf"Diagnostic BH-compatible selection, {family} $G_M^n$"
        )
        ax.grid(alpha=0.2)
        ax.legend()
        fig.tight_layout()
        fig.savefig(
            diagdir / "06_neutron_radius_vs_observed_bh_difference.png",
            dpi=180,
        )
        plt.close(fig)
    #endif

    pullscan = scan.loc[
        (scan["selection_type"] == "observed_BH_pull")
        & scan["valid"]
        & np.isfinite(scan["rM_fm"])
    ].copy()
    if len(pullscan):
        fig, ax = plt.subplots(figsize=(7.2, 5.4))
        ax.plot(
            pullscan["selection_cut_display"],
            pullscan["rM_fm"],
            marker="o", linestyle="none",
        )
        ax.axhline(
            0.864, linewidth=1.0, linestyle="--",
            label=r"PDG $r_{M,n}=0.864$ fm",
        )
        for _, row in pullscan.iterrows():
            ax.annotate(
                f"N={int(row['N'])}",
                (row["selection_cut_display"], row["rM_fm"]),
                xytext=(4, 4), textcoords="offset points", fontsize=8,
            )
        #endfor
        ax.set_xlabel(r"Observed BH compatibility cut, $|z_{\rm BH}|$")
        ax.set_ylabel(r"$r_{M,n}$ [fm]")
        ax.set_title(
            rf"Diagnostic BH-pull selection, {family} $G_M^n$"
        )
        ax.grid(alpha=0.2)
        ax.legend()
        fig.tight_layout()
        fig.savefig(
            diagdir / "07_neutron_radius_vs_bh_pull.png", dpi=180
        )
        plt.close(fig)
    #endif

    print(
        "[neutron-empirical-BH] diagnostic only: observed BH compatibility "
        "uses a reference-BH FF assumption and is not a production selection"
    )
    print(
        f"[neutron-empirical-BH] outputs -> {diagdir}"
    )
#enddef


def save_neutron_magnetic_only_study(
        evaluated: pd.DataFrame,
        outdir: Path) -> None:
    """
    Explore a neutron magnetic-radius extraction with GE_n fixed externally.

    The primary result fixes GE_n to Atac et al. 2021.  GE_n=0 is retained only
    as a limiting diagnostic of how much the small electric contribution moves
    the extracted magnetic radius.

    The BH-purity threshold is scanned while the added Moradi-style BH-method
    uncertainty remains fixed at 5%.
    """
    diagdir = outdir / "magnetic_only"
    diagdir.mkdir(parents=True, exist_ok=True)

    finite_delta = evaluated.loc[
        np.isfinite(evaluated["bh_delta"]), "bh_delta"
    ].to_numpy(float)
    thresholds = [
        0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.10,
        0.15, 0.20, 0.30, 0.50, 0.75, 1.00,
    ]
    if len(finite_delta):
        thresholds.append(float(np.max(finite_delta) * 1.001))
    #endif
    thresholds = sorted(set(thresholds))

    rows = []
    for threshold in thresholds:
        selected = evaluated.loc[
            np.isfinite(evaluated["bh_delta"])
            & (evaluated["bh_delta"] <= float(threshold))
        ].copy()
        if len(selected) < 3:
            continue
        #endif

        for ge_mode in ["atac2021", "zero"]:
            fit = fit_neutron_magnetic_only_bh(
                selected,
                ge_mode=ge_mode,
                bh_systematic_fraction=0.05,
            )
            fit["bh_threshold"] = float(threshold)
            fit["bh_threshold_percent"] = 100.0 * float(threshold)
            fit["max_selected_delta"] = float(
                np.nanmax(selected["bh_delta"].to_numpy(float))
            )
            rows.append(fit)
            print(
                f"[neutron-GM-only] threshold={100.0*threshold:7.2f}% "
                f"GE={ge_mode:8s} N={fit['N']:2d} "
                f"chi2/ndof={fit['chi2_ndof']:.3f} "
                f"rM={fit['rM_fm']:.4f}+/-{fit['rM_fit_err_fm']:.4f} fm"
            )
        #endfor
    #endfor

    scan = pd.DataFrame(rows)
    scan.to_csv(
        diagdir / "01_neutron_magnetic_only_threshold_scan.csv",
        index=False,
    )

    if not len(scan):
        print("[neutron-GM-only] no thresholds contained enough points")
        return
    #endif

    # Direct comparison of the primary fixed-Atac result with the GE=0 limit.
    fig, ax = plt.subplots(figsize=(8.0, 5.4))
    for ge_mode, label in [
        ("atac2021", r"$G_E^n$ fixed: Atac et al. 2021"),
        ("zero", r"$G_E^n=0$ diagnostic"),
    ]:
        p = scan.loc[scan["ge_mode"] == ge_mode].sort_values(
            "bh_threshold_percent"
        )
        ax.errorbar(
            p["bh_threshold_percent"],
            p["rM_fm"],
            yerr=p["rM_fit_err_fm"],
            marker="o" if ge_mode == "atac2021" else "s",
            linestyle="none",
            fillstyle="full" if ge_mode == "atac2021" else "none",
            capsize=2,
            label=label,
        )
    #endfor
    ax.axvline(5.0, linewidth=0.8, linestyle=":")
    ax.set_xscale("log")
    ax.set_xlabel(r"$|1-R_{\rm BH}^{\rm KM15}|$ threshold (%)")
    ax.set_ylabel(r"$r_{M,n}$ [fm]")
    ax.set_title(
        r"Exploratory neutron magnetic-radius extraction: fixed $G_E^n$"
    )
    ax.grid(alpha=0.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(
        diagdir / "02_neutron_magnetic_radius_threshold_scan.png",
        dpi=180,
    )
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.0, 5.4))
    for ge_mode, label in [
        ("atac2021", r"$G_E^n$ fixed: Atac et al. 2021"),
        ("zero", r"$G_E^n=0$ diagnostic"),
    ]:
        p = scan.loc[scan["ge_mode"] == ge_mode].sort_values(
            "bh_threshold_percent"
        )
        ax.plot(
            p["bh_threshold_percent"],
            p["chi2_ndof"],
            marker="o" if ge_mode == "atac2021" else "s",
            linestyle="none",
            fillstyle="full" if ge_mode == "atac2021" else "none",
            label=label,
        )
    #endfor
    ax.axvline(5.0, linewidth=0.8, linestyle=":")
    ax.set_xscale("log")
    ax.set_xlabel(r"$|1-R_{\rm BH}^{\rm KM15}|$ threshold (%)")
    ax.set_ylabel(r"$\chi^2/\mathrm{ndof}$")
    ax.set_title("Neutron magnetic-only BH fit quality")
    ax.grid(alpha=0.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(
        diagdir / "03_neutron_magnetic_only_chi2_threshold_scan.png",
        dpi=180,
    )
    plt.close(fig)

    # At each of the nine fixed-|t| bins, extract a single GM_n value directly
    # from all phi points with GE_n fixed to Atac.  This is not used to obtain
    # the radius; it visualizes the actual finite-|t| magnetic information.
    point_rows = []
    for kin_bin, part in evaluated.groupby("kin_bin", sort=True):
        qbar = float(np.nanmean(part["t_abs"]))
        ge = float(neutron_atac_ge(
            np.asarray([qbar]), 0.505, 1.655, 0.909
        )[0])

        q = part["t_abs"].to_numpy(float)
        tau = q / (4.0 * MP2)
        inv = 1.0 / (1.0 + tau)
        y = part["xs"].to_numpy(float)
        err = np.sqrt(
            part["stat"].to_numpy(float)**2
            + part["sys"].to_numpy(float)**2
            + (0.05 * y)**2
        )
        A = part["bh_A"].to_numpy(float)
        B = part["bh_B"].to_numpy(float)
        C = part["bh_C"].to_numpy(float)
        ge_arr = neutron_atac_ge(q, 0.505, 1.655, 0.909)

        # One GM amplitude at this narrow-|t| bin; initialize from Kelly.
        gm_kelly = (
            part["km15_F1"].to_numpy(float)
            + part["km15_F2"].to_numpy(float)
        )
        gm0 = float(np.nanmean(gm_kelly))

        def chi2_bin(GM_n):
            gm_arr = np.full_like(q, float(GM_n))
            f1 = (ge_arr + tau * gm_arr) * inv
            f2 = (gm_arr - ge_arr) * inv
            pred = bh_from_f1f2(A, B, C, f1, f2)
            pull = (pred - y) / err
            return float(np.dot(pull, pull))
        #enddef

        m = Minuit(chi2_bin, GM_n=gm0)
        m.errordef = Minuit.LEAST_SQUARES
        # Neutron GM is known to be negative over this low-|t| region.  This
        # guard prevents the mathematically equivalent positive square-root
        # branch from being selected by the BH-only cross section.
        m.limits["GM_n"] = (-3.0, -0.01)
        m.migrad()
        m.hesse()

        gm_fit = float(m.values["GM_n"])
        gm_err = float(m.errors["GM_n"])
        point_rows.append({
            "kin_bin": int(kin_bin),
            "N_phi": int(len(part)),
            "Q2_mean_GeV2": float(np.nanmean(part["Q2"])),
            "xB_mean": float(np.nanmean(part["xB"])),
            "t_abs_mean_GeV2": qbar,
            "GE_n_fixed_Atac2021": ge,
            "GM_n_fit": gm_fit,
            "GM_n_fit_err": gm_err,
            "GM_n_over_mu_n": gm_fit / MU_N,
            "GM_n_Kelly_mean": gm0,
            "chi2": float(m.fval),
            "ndof": max(1, int(len(part) - 1)),
            "chi2_ndof": float(m.fval / max(1, len(part) - 1)),
            "valid": bool(m.valid),
        })
    #endfor

    finite_gm = pd.DataFrame(point_rows)
    finite_gm.to_csv(
        diagdir / "04_neutron_finite_t_GM_extraction.csv",
        index=False,
    )

    fig, ax = plt.subplots(figsize=(8.0, 5.4))
    ax.errorbar(
        finite_gm["t_abs_mean_GeV2"],
        finite_gm["GM_n_over_mu_n"],
        yerr=finite_gm["GM_n_fit_err"] / abs(MU_N),
        marker="o", linestyle="none", capsize=2,
        label=r"BH extraction, $G_E^n$ fixed to Atac",
    )
    kelly_norm = finite_gm["GM_n_Kelly_mean"] / MU_N
    ax.plot(
        finite_gm["t_abs_mean_GeV2"],
        kelly_norm,
        marker="s", linestyle="none", fillstyle="none",
        label="Kelly reference at measured bins",
    )
    ax.set_xlabel(r"$|t|$ [GeV$^2$]")
    ax.set_ylabel(r"$G_M^n/\mu_n$")
    ax.set_title(
        r"Finite-$|t|$ neutron magnetic form factor from preliminary BH data"
    )
    ax.grid(alpha=0.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(
        diagdir / "04_neutron_finite_t_GM_extraction.png",
        dpi=180,
    )
    plt.close(fig)

    print(f"[neutron-GM-only] diagnostics -> {diagdir}")
#enddef


def run_neutron_analysis(args) -> int:
    """
    Exploratory neutron-radius extension.

    This intentionally does not reuse the proton closure machinery verbatim:
    GE_p(0)=1 whereas GE_n(0)=0, so the proton normalized-Sachs radius families
    are mathematically inappropriate for the neutron electric radius.  The
    neutron branch instead fits <r_E,n^2> through the low-Q2 slope of the Atac
    GE_n functional form and optionally fits GM_n with a normalized IP1 shape.
    """
    input_path = Path(args.neutron_file)
    outdir = Path(args.outdir) / "neutron_clas12_preliminary"
    outdir.mkdir(parents=True, exist_ok=True)

    data = load_ndvcs_preliminary_unpolarized(input_path)

    # Evaluate every plausible neutron phi mapping before doing any BH-purity
    # selection or radius fit.  The audit is based on absolute data/model
    # agreement, not merely on R_BH within the model.
    convention_tables: Dict[str, pd.DataFrame] = {}
    for convention in NEUTRON_PHI_CONVENTIONS:
        cache_name = (
            "km15_neutron_bh_decomposition_"
            + convention.replace("-", "_")
            + ".csv"
        )
        convention_tables[convention] = evaluate_km15_neutron_dataframe(
            data,
            ebeam=float(args.neutron_ebeam),
            workers=int(args.workers),
            cache_path=outdir / cache_name,
            phi_convention=convention,
            force=bool(args.force_km15),
        )
    #endfor

    if args.neutron_phi_convention == "auto":
        preferred_phi_convention = save_neutron_phi_convention_audit(
            convention_tables, outdir
        )
    else:
        # Still write the audit even when the user overrides the convention, so
        # the consequence of the override remains visible and reproducible.
        auto_preferred = save_neutron_phi_convention_audit(
            convention_tables, outdir
        )
        preferred_phi_convention = str(args.neutron_phi_convention)
        print(
            "[neutron-phi] user override: "
            f"{preferred_phi_convention} "
            f"(automatic audit preferred {auto_preferred})"
        )
    #endif

    evaluated = convention_tables[preferred_phi_convention].copy()

    # Preserve the traditional cache/output filenames for downstream tools,
    # but they now contain the explicitly selected phi convention.
    evaluated.to_csv(outdir / "km15_neutron_bh_decomposition.csv", index=False)

    # Save the full evaluated table in both the script's nb units and the
    # original-note pb units for easy cross-checking.
    evaluated["xs_pb_GeV4"] = 1000.0 * evaluated["xs"]
    evaluated["stat_pb_GeV4"] = 1000.0 * evaluated["stat"]
    evaluated["sys_pb_GeV4"] = 1000.0 * evaluated["sys"]
    evaluated.to_csv(outdir / "neutron_evaluated_points.csv", index=False)

    # Compare the unfitted KM15 BH/full-EP predictions directly with the data
    # using the selected convention before attempting any neutron-radius fit.
    save_neutron_model_comparison_diagnostics(evaluated, outdir)

    # Before extrapolating to a neutron charge or magnetic radius, quantify the
    # actual finite-|t| GE_n/GM_n leverage of these measured BH cross sections.
    save_neutron_bh_sensitivity_study(evaluated, outdir)

    # The sensitivity/Fisher study shows that these cross sections are
    # overwhelmingly magnetic and nearly singular in simultaneous GE_n/GM_n.
    # Therefore also run the complementary one-unknown extraction in which
    # GE_n is supplied externally and only GM_n is fitted.
    save_neutron_magnetic_only_study(evaluated, outdir)

    # Select the GM_n extrapolation family by pseudodata closure on this exact
    # neutron phase space, then repeat the threshold scan with that family.
    neutron_gm_family = run_neutron_magnetic_function_closure(
        evaluated, args, outdir
    )
    save_neutron_magnetic_selected_family_study(
        evaluated, neutron_gm_family, outdir
    )

    # Empirically test whether the data themselves are compatible with the
    # reference BH prediction, and compare that with KM15's model-defined
    # BH-purity ranking. These selections are diagnostics only because the
    # reference BH contains an assumed neutron form factor.
    save_neutron_empirical_bh_compatibility_study(
        evaluated, neutron_gm_family, outdir
    )

    # Whole-bin angular-shape test: float one finite-|t| GM amplitude in each
    # kinematic bin and test the complete eight-point phi dependence against BH.
    save_neutron_phi_shape_bh_compatibility_study(
        evaluated, outdir
    )

    print(
        "[neutron] preliminary CLAS12 nDVCS: "
        f"N={len(evaluated)}, Ebeam={args.neutron_ebeam:.3f} GeV, "
        f"phi={preferred_phi_convention}"
    )
    print(
        "[neutron] KM15/Gepard decomposition checks: "
        f"max EP relerr={np.nanmax(evaluated['km15_ep_decomp_relerr']):.3e}, "
        f"max BH quadratic relerr={np.nanmax(evaluated['bh_quad_relerr']):.3e}"
    )

    finite_delta = evaluated.loc[
        np.isfinite(evaluated["bh_delta"]), "bh_delta"
    ].to_numpy(float)
    thresholds = [
        0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.10,
        0.15, 0.20, 0.30, 0.50, 0.75, 1.00,
    ]
    if len(finite_delta):
        thresholds.append(float(np.max(finite_delta) * 1.001))
    #endif
    thresholds = sorted(set(thresholds))

    scan_rows = []
    for threshold in thresholds:
        selected = evaluated.loc[
            np.isfinite(evaluated["bh_delta"])
            & (evaluated["bh_delta"] <= float(threshold))
        ].copy()
        if len(selected) < 8:
            continue
        #endif

        for mode in ["kelly", "ip1"]:
            fit = fit_neutron_atac_bh(
                selected,
                magnetic_mode=mode,
                # Keep the Moradi-style additional BH-method uncertainty fixed
                # at 5% while the selection threshold is scanned.  Otherwise a
                # 100% threshold would artificially give every point a 100%
                # uncertainty and conceal DVCS-contamination effects.
                bh_systematic_fraction=0.05,
            )
            fit["bh_threshold"] = float(threshold)
            fit["bh_threshold_percent"] = 100.0 * float(threshold)
            fit["max_selected_delta"] = float(
                np.nanmax(selected["bh_delta"].to_numpy(float))
            )
            scan_rows.append(fit)
            print(
                f"[neutron] threshold={100.0*threshold:7.2f}% "
                f"mode={mode:5s} N={fit['N']:2d} "
                f"chi2/ndof={fit['chi2_ndof']:.3f} "
                f"<rn^2>={fit['rn2_fm2']:+.4f}+/-"
                f"{fit['rn2_fit_err_fm2']:.4f} fm^2 "
                + (
                    f"rM={fit['rM_fm']:.4f}+/-{fit['rM_fit_err_fm']:.4f} fm"
                    if mode == "ip1" else
                    "GM=Kelly"
                )
            )
        #endfor
    #endfor

    scan = pd.DataFrame(scan_rows)
    scan.to_csv(outdir / "neutron_radius_threshold_scan.csv", index=False)

    # Nominal 5% summary.
    nominal = scan.loc[np.isclose(scan["bh_threshold"], 0.05)].copy()
    nominal.to_csv(outdir / "neutron_radius_nominal_5pct.csv", index=False)

    # Reference from Atac et al. Nature Communications 12, 1759 (2021).
    reference = pd.DataFrame([{
        "source": "Atac et al., Nature Communications 12, 1759 (2021)",
        "rn2_fm2": -0.110,
        "rn2_err_fm2": 0.008,
        "A_GeV2": 0.505,
        "A_err_GeV2": 0.079,
        "B": 1.655,
        "B_err": 0.126,
        "C": 0.909,
        "C_err": 0.583,
    }])
    reference.to_csv(outdir / "neutron_charge_radius_reference_atac2021.csv", index=False)

    if len(scan):
        for mode, label in [
            ("kelly", r"$G_M^n$ fixed to Kelly"),
            ("ip1", r"$G_M^n$ fitted (IP1)"),
        ]:
            part = scan.loc[scan["magnetic_mode"] == mode].sort_values(
                "bh_threshold_percent"
            )
            if not len(part):
                continue
            #endif
            fig, ax = plt.subplots(figsize=(7.5, 5.0))
            ax.errorbar(
                part["bh_threshold_percent"],
                part["rn2_fm2"],
                yerr=part["rn2_fit_err_fm2"],
                marker="o",
                linestyle="none",
                capsize=2,
                label=label,
            )
            ax.axhline(-0.110, linewidth=1.0, linestyle="--",
                       label=r"Atac et al. $-0.110$ fm$^2$")
            ax.axvline(5.0, linewidth=0.8, linestyle=":")
            ax.set_xscale("log")
            ax.set_xlabel(r"$|1-R_{\rm BH}^{\rm KM15}|$ threshold (%)")
            ax.set_ylabel(r"$\langle r_{E,n}^2\rangle$ [fm$^2$]")
            ax.set_title("Preliminary CLAS12 neutron charge-radius threshold scan")
            ax.grid(alpha=0.2)
            ax.legend()
            fig.tight_layout()
            fig.savefig(
                outdir / f"01_neutron_rn2_threshold_{mode}.png", dpi=180
            )
            plt.close(fig)
        #endfor

        ip1 = scan.loc[scan["magnetic_mode"] == "ip1"].sort_values(
            "bh_threshold_percent"
        )
        if len(ip1):
            fig, ax = plt.subplots(figsize=(7.5, 5.0))
            ax.errorbar(
                ip1["bh_threshold_percent"],
                ip1["rM_fm"],
                yerr=ip1["rM_fit_err_fm"],
                marker="o",
                linestyle="none",
                capsize=2,
            )
            ax.axvline(5.0, linewidth=0.8, linestyle=":")
            ax.set_xscale("log")
            ax.set_xlabel(r"$|1-R_{\rm BH}^{\rm KM15}|$ threshold (%)")
            ax.set_ylabel(r"$r_{M,n}$ [fm]")
            ax.set_title("Exploratory CLAS12 neutron magnetic-radius threshold scan")
            ax.grid(alpha=0.2)
            fig.tight_layout()
            fig.savefig(outdir / "02_neutron_rM_threshold_ip1.png", dpi=180)
            plt.close(fig)
        #endif
    #endif

    print(
        "[neutron] NOTE: this is exploratory.  The preliminary analysis note "
        "states that radiative corrections were not yet included, and its "
        "bin-integral correction uses a BH calculation."
    )
    print(f"[neutron] outputs -> {outdir}")
    return 0
#enddef


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    mode = p.add_mutually_exclusive_group()
    mode.add_argument(
        "--only-pass2",
        action="store_true",
        help="Run only the CLAS12 pass-2 analysis",
    )
    mode.add_argument(
        "--only-pass1",
        action="store_true",
        help="Run only the CLAS12 Lee 2026 diagnostic",
    )
    mode.add_argument(
        "--only-clas6",
        action="store_true",
        help="Run only the CLAS6/Gepard validation",
    )
    mode.add_argument(
        "--only-saylor",
        action="store_true",
        help="Run only the CLAS6 Hirlinger Saylor et al. 2018 measurement",
    )
    mode.add_argument(
        "--only-halla-defurne",
        action="store_true",
        help="Run only Hall A Defurne et al. 2015",
    )
    mode.add_argument(
        "--only-clas6-pass1-combined",
        action="store_true",
        help=(
            "Run CLAS6 and CLAS12 Lee 2026 inputs and then fit them "
            "simultaneously with common form-factor parameters"
        ),
    )
    mode.add_argument(
        "--all-including-pass2",
        action="store_true",
        help="Run CLAS6, pass 1, combined CLAS6+pass1, and pass 2",
    )
    p.add_argument(
        "--pass1-csv",
        default=DEFAULT_PASS1_CSV,
        help="CLAS12 Lee 2026 input table; default is authoritative CLAS Physics Database clasdb_E214M1.txt",
    )
    p.add_argument(
        "--saylor-file",
        default=DEFAULT_SAYLOR_FILE,
        help="CLAS6 Saylor 2018 saylor_CLAS6.txt or normalized CSV",
    )
    p.add_argument(
        "--georges-file",
        default=DEFAULT_GEORGES_FILE,
        help="Hall A Georges 2022 E12-06-114 supplemental file (default .xlsx is read without optional Excel packages)",
    )
    p.add_argument(
        "--download-saylor",
        action="store_true",
        help="Attempt to download the public Saylor supplemental TXT",
    )
    p.add_argument(
        "--force-saylor-download",
        action="store_true",
        help="Redownload the Saylor supplemental even if the local file exists",
    )
    p.add_argument(
        "--include-saylor",
        action="store_true",
        help=(
            "Opt in to CLAS6 Saylor 2018 diagnostics. It is excluded from the "
            "preliminary Jo+Lee production ensemble by default."
        ),
    )
    p.add_argument(
        "--include-halla-defurne",
        action="store_true",
        help=(
            "Opt in to Hall A Defurne 2015 diagnostics. It is excluded from "
            "the preliminary Jo+Lee production ensemble by default."
        ),
    )
    p.add_argument(
        "--include-georges",
        action="store_true",
        help=(
            "Opt in to Hall A Georges 2022 E12-06-114 diagnostics. Georges is "
            "diagnostic-only and is never added to the Jo+Lee production fit."
        ),
    )
    p.add_argument(
        "--skip-saylor",
        action="store_true",
        help=(
            "Deprecated compatibility flag; Saylor is already excluded by "
            "default. Use --include-saylor to enable it."
        ),
    )
    p.add_argument(
        "--skip-halla-defurne",
        action="store_true",
        help=(
            "Deprecated compatibility flag; Hall A is already excluded by "
            "default. Use --include-halla-defurne to enable it."
        ),
    )
    p.add_argument(
        "--run-final-model-analysis",
        action="store_true",
        help=(
            "Consume the completed KM15/VGG99/GK16 external purity table, "
            "run model-specific 5%% selections, closure-based Sachs-family "
            "selection, production fits, and 1--10%% threshold systematics"
        ),
    )
    p.add_argument(
        "--bh-model-selection-results",
        default=FINAL_MODEL_SELECTION_DEFAULT,
        help=(
            "Completed evaluate_bh_model_selection.py all-point table; default: "
            + FINAL_MODEL_SELECTION_DEFAULT
        ),
    )
    p.add_argument(
        "--run-radius-bias-study",
        action="store_true",
        help=(
            "Run the expensive GE/GM polynomial, inverse-polynomial, and "
            "continued-fraction radius bias/variance pseudodata study"
        ),
    )
    p.add_argument(
        "--radius-bias-replicas",
        type=int,
        default=300,
        help=("Target number of SUCCESSFUL Gaussian pseudodata replica fits per "
              "truth model and fit family; failed fits are retried up to a "
              "finite automatic attempt cap"),
    )
    p.add_argument(
        "--threshold-correlation-replicas",
        type=int,
        default=200,
        help=(
            "Paired replicas used to estimate covariance among nested 3--7%% "
            "BH-purity threshold fits (default: 200). The same fluctuation of "
            "a shared event is reused at every threshold."
        ),
    )
    p.add_argument(
        "--radius-bias-workers",
        type=int,
        default=None,
        help=(
            "Process workers for the radius bias/variance replicas; "
            "defaults to --workers"
        ),
    )
    p.add_argument(
        "--radius-bias-seed",
        type=int,
        default=20260830,
        help="Random seed for the radius bias/variance pseudodata study",
    )
    p.add_argument(
        "--radius-bias-extended-truths",
        action="store_true",
        help=(
            "Add synthetic P/IP/CF generating families across an rE x rM "
            "grid to the empirical Kelly/AMT/Bernauer truth ensemble"
        ),
    )
    p.add_argument(
        "--radius-bias-radius-grid",
        default="0.81,0.85,0.89",
        help=(
            "Comma-separated proton synthetic truth radii in fm; the extended study "
            "uses the Cartesian rE x rM grid for every generating family"
        ),
    )
    p.add_argument(
        "--radius-bias-curvature-stress",
        choices=("sectoral", "none"),
        default="sectoral",
        help=(
            "Independent higher-order curvature stress truths for the proton "
            "closure ensemble. 'sectoral' (default) adds exponential- and "
            "dipole-like curvature in GE or GM separately at the central "
            "radius-grid point, while preserving the requested radius. "
            "'none' reproduces the older Kelly-template-only synthetic grid."
        ),
    )
    p.add_argument(
        "--radius-bias-pilot-replicas",
        type=int,
        default=10,
        help=(
            "Fixed-attempt replicas per pilot truth/family pair used to prune "
            "clearly noncompetitive Sachs-family pairs before the expensive "
            "full closure pass (default: 10). Set to 0 to disable pruning."
        ),
    )
    p.add_argument(
        "--radius-bias-pilot-min-keep",
        type=int,
        default=10,
        help=(
            "Minimum number of numerically reliable family pairs protected "
            "through the full replica study after pilot pruning (default: 10)."
        ),
    )
    p.add_argument(
        "--radius-bias-pilot-max-keep",
        type=int,
        default=16,
        help=(
            "Maximum family-pair cohort normally carried into the full replica "
            "study after pilot pruning (default: 16)."
        ),
    )
    p.add_argument(
        "--radius-bias-pilot-relative-tolerance",
        type=float,
        default=0.35,
        help=(
            "Retain pilot families whose objective is within this fractional "
            "amount of the pilot leader, subject to the protected cohort "
            "rules (default: 0.35)."
        ),
    )
    p.add_argument(
        "--radius-bias-pilot-min-valid-fraction",
        type=float,
        default=0.50,
        help=(
            "Minimum valid-fit fraction required for a family pair to be treated "
            "as numerically reliable in the cheap pilot stage (default: 0.50). "
            "This does not weaken the stricter full-closure convergence gates."
        ),
    )
    p.add_argument(
        "--radius-bias-pilot-absolute-tolerance-fm",
        type=float,
        default=0.020,
        help=(
            "Absolute objective tolerance in fm used together with the relative "
            "pilot tolerance (default: 0.020 fm)."
        ),
    )
    p.add_argument(
        "--neutron-radius-bias-radius-grid",
        default="0.82,0.86,0.90",
        help=(
            "Comma-separated neutron magnetic-radius truth grid in fm. The "
            "default spans a representative range around the accepted neutron "
            "magnetic radius without inheriting the proton grid."
        ),
    )
    p.add_argument(
        "--radius-bias-min-global-valid-fraction",
        type=float,
        default=0.95,
        help=(
            "Minimum aggregate valid-replica fraction required for a Sachs "
            "family to enter the closure ranking (default: 0.95)."
        ),
    )
    p.add_argument(
        "--radius-bias-min-scenario-valid-fraction",
        type=float,
        default=0.80,
        help=(
            "Minimum valid-replica fraction in every individual truth scenario "
            "required for a Sachs family to enter the closure ranking "
            "(default: 0.80)."
        ),
    )
    p.add_argument(
        "--clas6-dataset-id",
        type=int,
        default=CLAS6_GEPARD_DATASET_ID,
        help="Gepard dataset ID for CLAS6 validation (paper data: 98)",
    )
    p.add_argument(
        "--neutron",
        action="store_true",
        help=(
            "Run the exploratory CLAS12 neutron-DVCS radius workflow instead "
            "of the proton analysis.  This uses the preliminary unpolarized "
            "Tables A.1-A.9 and neutron-specific GE_n/GM_n parameterizations."
        ),
    )
    p.add_argument(
        "--neutron-file",
        default=DEFAULT_NDVCS_UNPOLARIZED,
        help=(
            "Machine-readable preliminary CLAS12 neutron-DVCS unpolarized "
            "cross-section table; default: " + DEFAULT_NDVCS_UNPOLARIZED
        ),
    )
    p.add_argument(
        "--neutron-ebeam",
        type=float,
        default=DEFAULT_NDVCS_EBEAM,
        help=(
            "Effective beam energy in GeV for the combined preliminary "
            "neutron table (default: 10.45 GeV, as quoted approximately in "
            "the analysis note)."
        ),
    )
    p.add_argument(
        "--neutron-phi-convention",
        choices=("auto",) + NEUTRON_PHI_CONVENTIONS,
        default="auto",
        help=(
            "Neutron table-to-Gepard phi mapping.  Default 'auto' evaluates "
            "identity, 180-minus, and 180-shift and selects the convention "
            "with the best full-KM15 absolute chi2/point.  An explicit choice "
            "overrides the automatic selection but still writes the audit."
        ),
    )
    p.add_argument(
        "--csv",
        default=DEFAULT_CSV,
        help="CLAS12 DVCS cross-section CSV",
    )
    p.add_argument(
        "--outdir",
        default=DEFAULT_OUTDIR,
        help="Output directory",
    )
    p.add_argument(
        "--ebeam",
        type=float,
        default=10.604,
        help="Beam energy in GeV",
    )
    p.add_argument(
        "--workers",
        type=int,
        default=min(6, os.cpu_count() or 1),
        help="KM15 worker processes",
    )
    p.add_argument(
        "--bh-cuts",
        type=float,
        nargs="+",
        default=[0.01, 0.02, 0.03, 0.04, 0.05],
        help="Nested |1-R_BH| thresholds; nominal paper values are 1-5%%",
    )
    p.add_argument(
        "--force-km15",
        action="store_true",
        help="Ignore the cached KM15 decomposition",
    )
    p.add_argument(
        "--no-clas12-point-to-point-systematic",
        action="store_true",
        help=(
            "Disable the preliminary CLAS12 absolute point-to-point "
            "systematic. It is included by default."
        ),
    )
    p.add_argument(
        "--no-combination-normalization-systematic",
        action="store_true",
        help=(
            "Disable the correlated normalization nuisance from the CSV "
            "combination systematic. It is included by default."
        ),
    )
    p.add_argument(
        "--no-global-normalization-systematic",
        action="store_true",
        help=(
            "Disable the correlated 4.76%% target-thickness/absolute-charge "
            "normalization nuisance. It is included by default."
        ),
    )
    p.add_argument(
        "--no-bh-selection-systematic",
        action="store_true",
        help="Disable the paper's additional 1-5%% uncorrelated uncertainty",
    )
    p.add_argument(
        "--tmin",
        type=float,
        default=None,
        help="Optional minimum |t|; normally leave unset",
    )
    p.add_argument(
        "--tmax",
        type=float,
        default=None,
        help="Optional maximum |t|; normally leave unset",
    )


    return p
#enddef


def apply_optional_t_range(data: pd.DataFrame,
                           tmin: Optional[float],
                           tmax: Optional[float]) -> pd.DataFrame:
    result = data.copy()
    if tmin is not None:
        result = result.loc[result["t_abs"] >= tmin]
    #endif
    if tmax is not None:
        result = result.loc[result["t_abs"] <= tmax]
    #endif
    return result.copy()
#enddef


def run_pass2_analysis(args, return_results: bool = False):
    csv_path = Path(args.csv).expanduser().resolve()
    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"[input]  {csv_path}")
    print(f"[output] {outdir}")

    df = load_clas12_csv(csv_path)
    print(
        f"[data]   {len(df)} valid 10.6-GeV unpolarized "
        "cross-section points"
    )

    # Reuse the validated cache from the earlier analysis if present.
    cache = outdir / "km15_bh_decomposition.csv"
    old_cache_candidates = [
        Path("output/emff_from_bh_paper_method/km15_bh_decomposition.csv").resolve(),
        Path("output/emff_from_bh/km15_bh_decomposition.csv").resolve(),
    ]
    if (not cache.exists()) and (not args.force_km15):
        for old_cache in old_cache_candidates:
            if old_cache.exists():
                cache.parent.mkdir(parents=True, exist_ok=True)
                import shutil
                shutil.copy2(old_cache, cache)
                print(f"[KM15] Copied existing validated cache: {old_cache}")
                break
            #endif
        #endfor
    #endif

    df = evaluate_km15_dataframe(
        df,
        args.ebeam,
        max(1, args.workers),
        cache,
        args.force_km15,
    )

    finite_model = (
        np.isfinite(df["km15_ep"])
        & np.isfinite(df["km15_bh"])
        & np.isfinite(df["R_BH"])
        & (df["km15_ep"] > 0)
        & (df["km15_bh"] > 0)
    )
    df = df.loc[finite_model].copy().reset_index(drop=True)
    print(f"[KM15]   {len(df)} points have finite positive model predictions")

    max_ep_relerr = float(np.nanmax(df["km15_ep_decomp_relerr"]))
    max_bh_relerr = float(np.nanmax(df["bh_quad_relerr"]))
    print(
        "[check]  max KM15 EP decomposition relative error "
        f"= {max_ep_relerr:.3e}"
    )
    print(
        "[check]  max BH quadratic reconstruction relative error "
        f"= {max_bh_relerr:.3e}"
    )

    if max_ep_relerr > 1.0e-7:
        raise RuntimeError(
            "KM15 BH+INT+DVCS decomposition does not reproduce "
            "th_KM15.predict."
        )
    #endif

    if max_bh_relerr > 1.0e-7:
        raise RuntimeError(
            "Precomputed BH quadratic coefficients do not reproduce "
            "the nominal BH cross section."
        )
    #endif

    bh_cuts = sorted(float(cut) for cut in args.bh_cuts)
    if bh_cuts != [0.01, 0.02, 0.03, 0.04, 0.05]:
        warnings.warn(
            "The paper defines Fits 1-5 with 1%%,2%%,3%%,4%%,5%% cuts. "
            "Nonstandard --bh-cuts change that correspondence."
        )
    #endif

    selection_records = []
    selected_sets: Dict[float, pd.DataFrame] = {}

    for cut in bh_cuts:
        selected = df.loc[df["bh_delta"] <= cut].copy()
        selected = apply_optional_t_range(selected, args.tmin, args.tmax)
        selected_sets[cut] = selected

        flag = df["bh_delta"] <= cut
        df[f"BH_{int(round(100*cut))}pct"] = flag

        selection_records.append({
            "Set": int(round(100 * cut)),
            "BH_cut": cut,
            "criterion": f"|1-R_BH| <= {100*cut:.0f}%",
            "Npts": len(selected),
            "xB_min": float(selected["xB"].min()) if len(selected) else np.nan,
            "xB_max": float(selected["xB"].max()) if len(selected) else np.nan,
            "Q2_min": float(selected["Q2"].min()) if len(selected) else np.nan,
            "Q2_max": float(selected["Q2"].max()) if len(selected) else np.nan,
            "t_min": float(selected["t_abs"].min()) if len(selected) else np.nan,
            "t_max": float(selected["t_abs"].max()) if len(selected) else np.nan,
            "phi_min": float(selected["phi_deg"].min()) if len(selected) else np.nan,
            "phi_max": float(selected["phi_deg"].max()) if len(selected) else np.nan,
        })

        print(
            f"[select] |1-R_BH| <= {100*cut:.0f}% : "
            f"{len(selected)} points"
        )
    #endfor

    df.to_csv(
        outdir / "clas12_with_km15_bh_decomposition.csv",
        index=False,
    )
    save_bh_selection_plots(df, bh_cuts, outdir)

    include_ptp = not args.no_clas12_point_to_point_systematic
    include_bh_sys = not args.no_bh_selection_systematic
    include_comb_norm = not args.no_combination_normalization_systematic
    include_global_norm = not args.no_global_normalization_systematic

    print(
        "[errors] point uncertainty = "
        "sqrt(stat^2 + point-to-point^2 + BH-selection^2)"
        if include_ptp
        else "[errors] point uncertainty = sqrt(stat^2 + BH-selection^2)"
    )
    print(
        "[errors] CLAS12 point-to-point systematic: "
        + ("INCLUDED" if include_ptp else "NOT included")
    )
    print(
        "[norm]   row-dependent combination systematic nuisance: "
        + ("INCLUDED" if include_comb_norm else "NOT included")
    )
    print(
        f"[norm]   global {100*GLOBAL_SCALE_FRAC:.2f}% normalization nuisance: "
        + ("INCLUDED" if include_global_norm else "NOT included")
    )

    fit_results: List[FitResult] = []

    # Fits 1-5: dipole F1 and F2 on the five nested BH samples.
    for fit_index, cut in enumerate(bh_cuts, start=1):
        data = selected_sets[cut]
        if len(data) < 10:
            raise RuntimeError(
                f"Fit {fit_index} has only {len(data)} selected points."
            )
        #endif

        fr = fit_paper_model(
            data=data,
            kind="dipole",
            fit_name=f"Fit {fit_index}",
            bh_cut=cut,
            include_clas12_ptp_sys=include_ptp,
            include_bh_sys=include_bh_sys,
            include_combination_norm_sys=include_comb_norm,
            include_global_norm_sys=include_global_norm,
        )
        fit_results.append(fr)

        aE = fr.params[0]
        aM = fr.params[1]
        daE = math.sqrt(max(fr.covariance[0, 0], 0.0))
        daM = math.sqrt(max(fr.covariance[1, 1], 0.0))

        print(
            f"[fit] Fit {fit_index}: N={fr.npts:4d} "
            f"chi2/dof={fr.chi2_ndof:.3f} "
            f"aE={aE:.4f}+/-{daE:.4f} "
            f"aM={aM:.4f}+/-{daM:.4f} "
            f"beta_global={fr.meta.get('beta_global', np.nan):+.3f} "
            f"beta_comb={fr.meta.get('beta_comb', np.nan):+.3f}"
        )
    #endfor

    if 0.05 not in selected_sets:
        raise RuntimeError(
            "The exact paper-method Fit 6-8 stage requires the 5%% Set 5."
        )
    #endif

    set5 = selected_sets[0.05].copy()
    set5["fit_sigma_default"] = fit_sigma_errors(
        set5,
        0.05,
        include_ptp,
        include_bh_sys,
    )
    set5.to_csv(outdir / "set5_selected_points.csv", index=False)

    # Fits 6-8 exactly as in the paper.
    for kind, fit_name in [
        ("fit6_same_a", "Fit 6"),
        ("fit7_same_p", "Fit 7"),
        ("fit8_f2_kelly", "Fit 8"),
    ]:
        fr = fit_paper_model(
            data=set5,
            kind=kind,
            fit_name=fit_name,
            bh_cut=0.05,
            include_clas12_ptp_sys=include_ptp,
            include_bh_sys=include_bh_sys,
            include_combination_norm_sys=include_comb_norm,
            include_global_norm_sys=include_global_norm,
        )
        fit_results.append(fr)

        parameter_text = " ".join(
            f"{name}={value:.4f}+/-"
            f"{math.sqrt(max(fr.covariance[i, i], 0.0)):.4f}"
            for i, (name, value) in enumerate(
                zip(fr.param_names, fr.params)
            )
        )
        print(
            f"[fit] {fit_name}: N={fr.npts:4d} "
            f"chi2/dof={fr.chi2_ndof:.3f} {parameter_text} "
            f"beta_global={fr.meta.get('beta_global', np.nan):+.3f} "
            f"beta_comb={fr.meta.get('beta_comb', np.nan):+.3f}"
        )
    #endfor

    fit5 = next(r for r in fit_results if r.name == "Fit 5")
    fit8 = next(r for r in fit_results if r.name == "Fit 8")

    # Save machine-readable outputs.
    pd.DataFrame(
        [fitresult_to_record(fr) for fr in fit_results]
    ).to_csv(outdir / "fit_results.csv", index=False)

    arrays = {}
    metadata = []
    for i, fr in enumerate(fit_results):
        key = f"fit_{i + 1}"
        arrays[key + "_params"] = fr.params
        arrays[key + "_cov"] = fr.covariance
        metadata.append({
            "key": key,
            "name": fr.name,
            "param_names": fr.param_names,
            "bh_cut": fr.bh_cut,
            "model_kind": fr.model_kind,
            "meta": fr.meta,
        })
    #endfor

    np.savez_compressed(outdir / "fit_covariances.npz", **arrays)
    with open(outdir / "fit_metadata.json", "w") as handle:
        json.dump(metadata, handle, indent=2)
    #endwith

    save_summary_tables(fit_results, selection_records, outdir)

    # Paper-style and observable-space plots.
    save_selected_cross_section_pages(
        set5,
        fit5,
        fit8,
        outdir,
    )
    save_fit1_to_fit5_plots(fit_results, set5, outdir)
    save_fit5_to_fit8_sachs_plot(fit_results, set5, outdir)
    save_radii_plot(fit_results, outdir)
    save_low_q2_ratio_plots(fit_results, set5, outdir)
    save_elastic_reference_comparison(fit_results, set5, outdir, "Fit 5")
    save_elastic_reference_comparison(fit_results, set5, outdir, "Fit 8")
    save_bh_local_f1_f2_sensitivity(outdir, set5, fit5)

    print("\n[paper-method radii]")
    for fit_name in ["Fit 5", "Fit 6", "Fit 7", "Fit 8"]:
        fr = next(r for r in fit_results if r.name == fit_name)
        print(
            f"  {fit_name}: "
            f"rE={fr.rE_fm:.5f}+/-{fr.rE_err_fm:.5f} fm, "
            f"rM={fr.rM_fm:.5f}+/-{fr.rM_err_fm:.5f} fm, "
            f"beta_global={fr.meta.get('beta_global', np.nan):+.3f}, "
            f"beta_comb={fr.meta.get('beta_comb', np.nan):+.3f}, "
            f"chi2/dof={fr.chi2_ndof:.3f}"
        )
    #endfor

    print(f"\nDone. Results are in {outdir}")

#enddef
    if return_results:
        return {"label": "CLAS12 pass 2", "results": fit_results,
                "set5": set5, "outdir": outdir}
    #endif
    return 0
#enddef


def main() -> int:
    args = build_parser().parse_args()
    print_run_plan(args)

    if Minuit is None:
        raise RuntimeError(
            "iminuit is required. Install it in the active Python environment."
        )
    #endif

    if args.neutron:
        return run_neutron_analysis(args)
    #endif

    if args.only_pass2:
        return run_pass2_analysis(args)
    #endif
    if args.only_pass1:
        return run_pass1_validation(args)
    #endif
    if args.only_clas6:
        return run_clas6_validation(args)
    #endif
    if args.only_saylor:
        return run_saylor_validation(args)
    #endif
    if args.only_halla_defurne:
        return run_halla_defurne_validation(args)
    #endif
    if args.only_clas6_pass1_combined:
        print(
            "[NOTE] --only-clas6-pass1-combined is a restricted legacy/debug "
            "mode. It intentionally does NOT run Saylor or the other "
            "measurement combinations."
        )
        jo = run_clas6_validation(args, return_results=True)
        pass1 = run_pass1_validation(args, return_results=True)
        run_measurement_combination(
            [jo, pass1], args, "jo2015_plus_pass1"
        )
        return 0
    #endif
    if args.all_including_pass2:
        status = run_published_default(args)
        if status != 0:
            return status
        #endif
        return run_pass2_analysis(args)
    #endif

    return run_published_default(args)
#enddef


if __name__ == "__main__":
    raise SystemExit(main())
#endif
