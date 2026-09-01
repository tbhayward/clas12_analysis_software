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

Thus the 4.76% target-thickness/absolute-charge uncertainty moves every
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
  * CLAS12 Lee 2026 (all_bin_v3.csv)

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

SAYLOR_EBEAM = 5.88
SAYLOR_GLOBAL_SCALE_FRAC = 0.05
DEFAULT_SAYLOR_FILE = "import/saylor_CLAS6.txt"
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
DEFAULT_PASS1_CSV = "../imports/all_bin_v3.csv"
PASS1_GLOBAL_SCALE_FRAC = 0.31

# Georges 2022 supplemental tables provide statistical and TOTAL systematic
# uncertainties, while the publication states that the latter contain both
# point-to-point and correlated components.  The public table does not provide
# a decomposition that can be reconstructed unambiguously here.  Therefore the
# current combined-radius implementation treats the quoted total systematic as
# pointwise and does NOT add a separate Georges normalization nuisance.  This is
# conservative for point errors and, critically, avoids double counting a
# correlated contribution.  Keep this explicit until an authoritative
# point-to-point/correlated decomposition is supplied.

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
    Load the legacy/pass-1 CLAS12 DVCS cross-section table all_bin_v3.csv.

    Relevant columns in the attached file:
      xBavg, Q2avg, t_abs_avg, phiavg
      cross sections, ep->epg, exp
      cross sections, ep->epg, exp, stat. unc.
      cross sections, ep->epg, exp, syst. unc. (up/down)
      valid bin

    The pass-1 systematic column is treated as an absolute pointwise
    systematic uncertainty. The large 31% overall scale uncertainty is kept
    separate and represented by one correlated normalization nuisance.
    """
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
            "Pass-1 CSV is missing required columns: "
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

    # Use the mean absolute up/down uncertainty as the symmetric pointwise
    # systematic entering chi2. In the attached CSV up and down are in fact
    # equal for the inspected rows, but this remains robust if some bins differ.
    outdf["ptp_sys_abs"] = 0.5 * (
        np.abs(outdf["pass1_sys_up"].to_numpy(float))
        + np.abs(outdf["pass1_sys_down"].to_numpy(float))
    )

    # Canonical aliases used only by generic diagnostics.  The authoritative
    # pass-1 columns remain xs_stat and ptp_sys_abs.
    outdf["stat"] = outdf["xs_stat"]
    outdf["ptp_sys"] = outdf["ptp_sys_abs"]

    outdf["comb_sys_frac"] = 0.0
    outdf["scale_sys_frac"] = PASS1_GLOBAL_SCALE_FRAC
    outdf["ebeam"] = 10.604

    # Carry bin-boundary columns through for grouped phi plots when available.
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
        f"[PASS1] Loaded {len(outdf)} valid pass-1 CLAS12 "
        f"cross-section points from {csv_path}"
    )

    rel_stat = outdf["xs_stat"] / outdf["xs"]
    rel_sys = outdf["ptp_sys_abs"] / outdf["xs"]
    print(
        f"[PASS1] median stat/xs={100*np.nanmedian(rel_stat):.1f}%, "
        f"median pointwise syst/xs={100*np.nanmedian(rel_sys):.1f}%, "
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


def evaluate_km15_point(args: Tuple[int, float, float, float, float, float]) -> Dict[str, float]:
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

    tasks = [
        (int(i), float(r.xB), float(r.Q2), float(r.t_abs), float(r.phi_deg), float(ebeam))
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


def elastic_reference_curves(q: np.ndarray) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    """Named direct-elastic reference curves used consistently in all plots."""
    return {
        "Kelly 2004": kelly_sachs(q),
        "AMT 2007": amt2007_sachs(q),
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
    elif dataset_kind == "halla_defurne2015":
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
        # Gepard pt.err is the published total point uncertainty used in the
        # Moradi reproduction.  The Jo loader canonicalizes that quantity as
        # ``clas6_err_total`` and also copies it to ``xs_stat`` for the
        # standalone fit machinery.  Do not require a raw ``err`` column here:
        # it is not present in the dataframe returned by
        # load_clas6_gepard_dataset().
        if "clas6_err_total" in data.columns:
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
    elif dataset_kind == "halla_defurne2015":
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
        add_moradi_bh_systematic: bool = True) -> FitResult:
    """
    Simultaneous common-form-factor fit to any subset of Jo 2015,
    Saylor 2018, and CLAS12 Lee 2026.

    Each measurement has its own cross-section data and uncertainty model.
    Correlated normalization nuisances are independent between experiments:
      Jo 2015: none, matching the Moradi/Gepard treatment;
      Saylor 2018: 5% global-normalization nuisance;
      Hall A Defurne 2015: 2.8% global-normalization nuisance;
      CLAS12 Lee 2026: 31% conservative global normalization nuisance.
    """
    names, p0 = paper_model_setup(kind)

    nuisance_names = []
    nuisance_fracs = {}
    for spec in datasets:
        frac = float(spec.get("norm_frac", 0.0))
        if frac > 0.0:
            nname = "beta_" + re.sub(r"[^A-Za-z0-9]+", "_", str(spec["key"]))
            nuisance_names.append(nname)
            nuisance_fracs[nname] = frac
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

        prepared.append({
            "key": spec["key"],
            "kind": spec["kind"],
            "norm_frac": float(spec.get("norm_frac", 0.0)),
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
            beta = p[nuisance_index[nname]]
            total += float(beta**2)
        #endfor
        return total
    #enddef

    m = Minuit(chi2_minuit, *fit_p0, name=tuple(fit_names))
    m.errordef = Minuit.LEAST_SQUARES
    for name in names:
        m.limits[name] = (1.0e-6, None)
    #endfor
    for nname in nuisance_names:
        m.limits[nname] = (-10.0, 10.0)
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
        norm_frac: float = 0.0) -> Dict[str, object]:
    return {
        "key": key,
        "label": label,
        "kind": kind,
        "data": data,
        "norm_frac": float(norm_frac),
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
    ax.legend(fontsize=8, ncol=2)
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
    for spec in specs:
        frac = float(spec.get("norm_frac", 0.0))
        if frac > 0.0:
            nuisance_names.append(
                "beta_" + re.sub(r"[^A-Za-z0-9]+", "_", str(spec["key"]))
            )
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
        frac = float(spec.get("norm_frac", 0.0))
        prepared.append({
            "key": key, "norm_frac": frac,
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
            total += float(p[nuisance_index[nname]]**2)
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
        m.limits[nname] = (-10.0, 10.0)
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
        frac = float(spec.get("norm_frac", 0.0))
        beta_true = rng.normal(0.0, 1.0) if frac > 0.0 else 0.0
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

    For Pn, slope=c1. For IPn and CFn, slope=-c1. Higher-order coefficients
    retain the template curvature.
    """
    out = np.asarray(template_coeffs, dtype=float).copy()
    slope = radius_to_normalized_slope(radius_fm)
    if family.startswith("P"):
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
    """Fit a normalized P/IP/CF family to a smooth reference shape."""
    q = np.asarray(target_q, dtype=float)
    y = np.asarray(target_shape, dtype=float)
    order = int(re.findall(r"\d+", family)[0])
    names = tuple(f"c{i}" for i in range(1, order + 1))

    p0 = np.zeros(order, dtype=float)
    p0[0] = -2.0 / 0.71 if family.startswith("P") else 2.0 / 0.71

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
        qmax: float) -> List[Dict[str, object]]:
    """
    Build cross-family closure truths spanning a Cartesian rE x rM grid.

    Each family first approximates Kelly over the measured Q2 range to obtain
    a smooth higher-order curvature template. Then c1 is replaced so the
    requested radius is exact. This varies slope and functional family
    independently while retaining smooth proton-like curvature.
    """
    q_template = np.linspace(0.0, max(float(qmax), 0.05), 400)
    k_ge, k_gm = kelly_sachs(q_template)
    k_ge = k_ge / k_ge[0]
    k_gm = k_gm / k_gm[0]

    templates = {}
    for family in families:
        templates[family] = {
            "E": fit_sachs_family_shape_template(family, q_template, k_ge),
            "M": fit_sachs_family_shape_template(family, q_template, k_gm),
        }
    #endfor

    scenarios = []
    skipped = 0
    for family in families:
        for rE in radius_values:
            for rM in radius_values:
                ce = sachs_family_coefficients_with_radius(
                    templates[family]["E"], family, rE
                )
                cm = sachs_family_coefficients_with_radius(
                    templates[family]["M"], family, rM
                )
                if not synthetic_truth_is_physical(qmax, family, ce, cm):
                    skipped += 1
                    continue
                #endif

                def truth_fn(q, fam=family, e=ce.copy(), m=cm.copy()):
                    q = np.asarray(q, dtype=float)
                    ge = sachs_family_value(q, e, fam)
                    gm = MU_P * sachs_family_value(q, m, fam)
                    return ge, gm
                #enddef

                scenarios.append({
                    "truth_model": (
                        f"synthetic_{family}_rE{rE:.3f}_rM{rM:.3f}"
                    ),
                    "truth_group": family,
                    "truth_family": family,
                    "truth_rE_fm": float(rE),
                    "truth_rM_fm": float(rM),
                    "truth_fn": truth_fn,
                    "synthetic": True,
                })
            #endfor
        #endfor
    #endfor

    print(
        f"[radius-bias] synthetic truth ensemble: {len(scenarios)} accepted "
        f"scenario(s), {skipped} rejected by physical-shape screen; "
        f"radii={list(radius_values)}"
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


def save_extended_radius_bias_matrices(
        table: pd.DataFrame,
        families: Sequence[str],
        outdir: Path,
        minimum_global_valid_fraction: float,
        minimum_scenario_valid_fraction: float) -> None:
    """
    Make cross-family closure matrices and an eligibility-aware ranking.

    Eligibility requires both strong aggregate convergence and no catastrophic
    truth-scenario failure.  Convergence is a gate only; eligible families are
    still ranked purely by the bias-variance objective.
    """
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
                float(np.sqrt(np.mean(vals**2))) if len(vals) else np.nan
            )
            row[f"{quantity}_max_objective_all_truths_fm"] = (
                float(np.max(vals)) if len(vals) else np.nan
            )
            row[f"{quantity}_RMS_bias_all_truths_fm"] = (
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



def save_mixed_family_closure_ranking(
        table: pd.DataFrame,
        candidate_pairs: Sequence[str],
        outdir: Path,
        minimum_global_valid_fraction: float,
        minimum_scenario_valid_fraction: float) -> None:
    """Rank E/M family pairs and retain bias and replica variance separately."""
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
            bias = qpart["bias_fm"].to_numpy(float)
            stat = qpart["stat_RMS_fm"].to_numpy(float)
            obj = qpart["sqrt_stat2_plus_bias2_fm"].to_numpy(float)
            row[f"{quantity}_RMS_bias_fm"] = float(np.sqrt(np.nanmean(bias**2)))
            row[f"{quantity}_RMS_replica_std_fm"] = float(np.sqrt(np.nanmean(stat**2)))
            row[f"{quantity}_RMS_RMSE_fm"] = float(np.sqrt(np.nanmean(obj**2)))
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
    #endif
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

    Fitting families:
      P1..P4, IP1..IP4, CF2..CF4.

    Truth ensemble:
      empirical Kelly, AMT2007, Bernauer/A1, plus optional synthetic P/IP/CF
      closure truths spanning a Cartesian rE x rM grid.
    """
    outdir.mkdir(parents=True, exist_ok=True)

    # Exploratory fit-candidate set.  P1 is omitted because a linear
    # polynomial has not provided an adequate description in the studies so
    # far.  Fourth-order P4/IP4/CF4 candidates are temporarily omitted to keep
    # the independent GE x GM closure scan computationally manageable.  They
    # remain in the synthetic truth ensemble below, so the retained candidates
    # are still tested against higher-curvature underlying form factors.
    candidate_families = [
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
    # replica, and retaining P4/IP4/CF4 here tests whether the reduced
    # candidate set remains robust against more complicated underlying shapes.
    truth_families = [f"P{i}" for i in range(1, 5)]
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
            "truth_group": "Kelly",
            "truth_family": "empirical",
            "truth_fn": kelly_sachs,
            "synthetic": False,
        },
        {
            "truth_model": "AMT2007",
            "truth_group": "AMT2007",
            "truth_family": "empirical",
            "truth_fn": amt2007_sachs,
            "synthetic": False,
        },
        {
            "truth_model": "Bernauer_order8_polyxdipole",
            "truth_group": "Bernauer",
            "truth_family": "empirical",
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
        )
    #endif

    print(
        f"[radius-bias] total truth scenarios={len(scenarios)}, "
        f"fit families={len(families)}, replicas/scenario/family="
        f"{args.radius_bias_replicas}"
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
        max_attempts = max(5 * target_valid, target_valid + 50)
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

    replica_table = pd.DataFrame(replica_rows)
    replica_table.to_csv(
        outdir / "radius_bias_replica_results.csv",
        index=False,
    )

    agg_rows = []
    for family in families:
        for quantity in ["rE", "rM"]:
            part = table.loc[
                (table["family"] == family)
                & (table["quantity"] == quantity)
            ]
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
    return measurement_spec(
        key=str(bundle["key"]),
        label=str(bundle["label"]),
        kind=str(bundle["kind"]),
        data=selected,
        norm_frac=float(bundle.get("norm_frac", 0.0)),
    )
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
    "saylor2018": 2,
    "halla_georges2022": 3,
    "pass1": 4,
    "pass2": 5,
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
        bh_systematic_fraction: Optional[float] = None) -> Dict[str, object]:
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

    nuisance_names = []
    nuisance_fracs = {}
    for spec in datasets:
        frac = float(spec.get("norm_frac", 0.0))
        if frac > 0.0:
            nname = "beta_" + re.sub(r"[^A-Za-z0-9]+", "_", str(spec["key"]))
            nuisance_names.append(nname)
            nuisance_fracs[nname] = frac
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
        key = str(spec["key"]); frac = float(spec.get("norm_frac", 0.0))
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
            total += float(p[nuisance_index[nname]] ** 2)
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
        m.limits[nname] = (-10.0, 10.0)
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
    #endfor
    return result
#enddef


def summarize_bias_for_family(
        study_csv: Path,
        family: str) -> Dict[str, float]:
    """Preserve signed, RMS, and worst-case bias estimates for one family."""
    table = pd.read_csv(study_csv)
    fam = table.loc[table["family"].astype(str) == family].copy()
    if len(fam) == 0:
        raise RuntimeError(f"No {family} rows found in {study_csv}")
    #endif

    out = {}
    for quantity in ["rE", "rM"]:
        vals = fam.loc[
            fam["quantity"].astype(str) == quantity, "bias_fm"
        ].to_numpy(float)
        vals = vals[np.isfinite(vals)]
        if len(vals) == 0:
            raise RuntimeError(f"No finite {quantity} bias values for {family}")
        #endif
        out[f"{quantity}_mean_signed_bias_fm"] = float(np.mean(vals))
        out[f"{quantity}_RMS_bias_fm"] = float(np.sqrt(np.mean(vals**2)))
        out[f"{quantity}_max_abs_bias_fm"] = float(np.max(np.abs(vals)))
    #endfor
    return out
#enddef


def aggregate_model_bias_rankings(
        bias_dirs: Dict[str, Path],
        outdir: Path) -> Tuple[str, pd.DataFrame]:
    """
    Select one common Sachs family across all three purity-model samples.

    Replica convergence fractions are diagnostic only and do not veto a family.
    Rank every family with a finite closure objective in all purity models by
    the worst (largest) combined RMS bias-variance objective across
    KM15/VGG99/GK16.  This avoids
    choosing a different extrapolation form for each purity model and avoids
    allowing one especially favorable selection to dominate the choice.
    """
    rows = []
    by_model = {}
    for model, bdir in bias_dirs.items():
        path = bdir / "radius_bias_extended_eligibility_ranking.csv"
        if not path.exists():
            raise FileNotFoundError(f"Missing extended bias ranking: {path}")
        #endif
        table = pd.read_csv(path)
        by_model[model] = table.set_index("family")
    #endfor

    families = sorted(
        set.intersection(
            *[set(t.index.astype(str)) for t in by_model.values()]
        )
    )
    for family in families:
        row = {"family": family}
        objectives = []
        rankable_all = True
        for model in FINAL_MODEL_NAMES:
            t = by_model[model]
            obj = float(t.loc[family, "combined_RMS_objective_fm"])
            min_valid = float(t.loc[family, "minimum_scenario_valid_fraction"])
            global_valid = float(t.loc[family, "global_valid_fraction"])
            row[f"{model}_minimum_scenario_valid_fraction"] = min_valid
            row[f"{model}_global_valid_fraction"] = global_valid
            row[f"{model}_combined_RMS_objective_fm"] = obj
            row[f"{model}_eligible"] = bool(np.isfinite(obj))
            rankable_all = rankable_all and np.isfinite(obj)
            objectives.append(obj)
        #endfor
        row["eligible_all_models"] = bool(rankable_all)
        finite_obj = np.asarray(objectives, dtype=float)
        finite_obj = finite_obj[np.isfinite(finite_obj)]
        row["mean_combined_RMS_objective_fm"] = float(np.mean(finite_obj)) if len(finite_obj) else np.nan
        row["worst_combined_RMS_objective_fm"] = float(np.max(finite_obj)) if len(finite_obj) else np.nan
        rows.append(row)
    #endfor

    ranking = pd.DataFrame(rows).sort_values(
        ["eligible_all_models", "worst_combined_RMS_objective_fm",
         "mean_combined_RMS_objective_fm"],
        ascending=[False, True, True],
    )
    eligible = ranking.loc[ranking["eligible_all_models"]]
    if len(eligible) == 0:
        raise RuntimeError(
            "No Sachs family has a finite closure objective for all three "
            "model-selected samples."
        )
    #endif

    rank1 = str(eligible.iloc[0]["family"])
    outdir.mkdir(parents=True, exist_ok=True)
    ranking.to_csv(outdir / "model_aggregated_family_ranking.csv", index=False)

    with open(outdir / "bias_rank1_sachs_family.txt", "w") as fout:
        fout.write(f"bias_rank1_family={rank1}\n")
        fout.write(
            "criterion=finite closure objective in all KM15/VGG99/GK16 studies; "
            "convergence fractions are diagnostic only; minimum worst combined "
            "RMS bias-variance objective\n"
        )
        fout.write(
            "note=This is the closure rank-1 family. Final production selection "
            "may fall through to the next-ranked family if the nominal 5% "
            "production fits are invalid. Threshold-scan failures do not veto it.\n"
        )
    #endwith

    print(
        f"[final model analysis] closure rank-1 Sachs family is {rank1}; "
        "production-fit validation will determine the final family"
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
    """Single consolidated closure-family ranking for the final-analysis folder."""
    figures_dir.mkdir(parents=True, exist_ok=True)
    d = ranking.copy().sort_values("worst_combined_RMS_objective_fm")
    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    y = np.arange(len(d))
    ax.barh(y, d["worst_combined_RMS_objective_fm"].to_numpy(float))
    ax.set_yticks(y, d["family"].astype(str).tolist())
    ax.invert_yaxis()
    ax.set_xlabel("Worst-model combined RMS bias-variance objective (fm)")
    ax.set_title("Jo 2015 + Lee 2026 closure-family ranking")
    ax.grid(axis="x", alpha=0.22)
    for i, (_, row) in enumerate(d.iterrows()):
        if not bool(row["eligible_all_models"]):
            ax.text(
                float(row["worst_combined_RMS_objective_fm"]),
                i, "  ineligible", va="center", fontsize=8,
            )
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
            "Method systematic = quadrature of the joint BH-selection "
            "systematic and closure/extrapolation bias systematic."
        ),
        (
            "Joint BH-selection systematic = maximum shift from the nominal "
            "KM15 5% published-errors result over KM15/VGG99/GK16 and the "
            "3--7% BH-purity window."
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

    # Run a closure/bias study for each 5% model-selected sample.  This is
    # intentionally three studies: the extrapolation family must be robust to
    # the modestly different kinematic support induced by purity-model choice.
    bias_dirs = {}
    for model in FINAL_MODEL_NAMES:
        model_dir = closure_dir / model
        bias_dirs[model] = model_dir
        print(f"\n[final model analysis] starting radius-bias study for {model}")

        # Rebuild lightweight bundle copies whose set5 is the external 5% sample
        # so the existing closure machinery can be reused unchanged.
        model_bundles = []
        specs5 = selected[model][0.05]
        spec_by_key = {str(s["key"]): s for s in specs5}
        for bundle in physics_bundles:
            bcopy = dict(bundle)
            bcopy["set5"] = spec_by_key[str(bundle["key"])]["data"].copy()
            model_bundles.append(bcopy)
        #endfor

        if args.run_radius_bias_study:
            run_radius_bias_variance_study(model_bundles, args, model_dir)
        elif not (model_dir / "radius_bias_extended_eligibility_ranking.csv").exists():
            raise RuntimeError(
                f"Bias ranking missing for {model}: {model_dir}. "
                "Run with --run-radius-bias-study "
                "--radius-bias-extended-truths."
            )
        #endif
    #endfor

    bias_rank1_family, family_ranking = aggregate_model_bias_rankings(
        bias_dirs, summary_dir
    )

    # ------------------------------------------------------------------
    # Ordered production-validation fallback.
    #
    # The closure ranking supplies an ordered list of rankable families.
    # Replica convergence fractions are diagnostic only and impose no veto.
    # Starting from the best bias+variance score, require ONLY that the family
    # produce valid nominal 5% Moradi-error production fits for KM15, VGG99,
    # and GK16. If the nominal fit fails, retain the failure in the audit table
    # and automatically try the next-best eligible family.
    #
    # Threshold-scan failures do NOT disqualify a family. They are retained as
    # diagnostics and skipped when constructing the local selection envelope.
    # ------------------------------------------------------------------
    eligible_ranking = family_ranking.loc[
        family_ranking["eligible_all_models"].astype(bool)
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
            rank_row["worst_combined_RMS_objective_fm"]
        )
        print(
            f"\n[family fallback] trying rank {candidate_rank}: "
            f"{candidate_family} "
            f"(worst objective={candidate_score:.6f} fm)"
        )

        attempt = {
            "closure_rank": candidate_rank,
            "family": candidate_family,
            "worst_combined_RMS_objective_fm": candidate_score,
            "mean_combined_RMS_objective_fm": float(
                rank_row["mean_combined_RMS_objective_fm"]
            ),
            "passed_05pct_production": False,
            "selected_for_final_analysis": False,
            "failure_details": "",
        }

        candidate_fit_rows = []
        failure_details = []

        # Required 5% production fits with Moradi's additional 5% pointwise
        # systematic.
        for model in FINAL_MODEL_NAMES:
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
                    failure_details.append(
                        f"invalid 5% Moradi fit: {model}"
                    )
                #endif
            except Exception as exc:
                failure_details.append(
                    f"exception in 5% Moradi fit {model}: "
                    f"{type(exc).__name__}: {exc}"
                )
            #endtry
        #endfor

        passed_05pct = (
            len(candidate_fit_rows) == len(FINAL_MODEL_NAMES)
            and all(bool(row.get("valid", False)) for row in candidate_fit_rows)
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
            "criterion=first closure-ranked family with finite closure objectives, ordered by increasing "
            "worst combined RMS bias-variance objective, that yields valid "
            "nominal 5% Moradi-error production fits for KM15/VGG99/GK16. "
            "Threshold-scan failures do not disqualify the family.\n"
        )
    #endwith

    # Preserve the extrapolation-bias estimates for the family that actually
    # survived production validation, not automatically the closure rank-1
    # family.
    bias_rows = []
    for model in FINAL_MODEL_NAMES:
        bias = summarize_bias_for_family(
            bias_dirs[model] / "radius_bias_variance_study.csv",
            chosen_family,
        )
        bias_rows.append({"model": model, "family": chosen_family, **bias})
    #endfor
    bias_table = pd.DataFrame(bias_rows)
    bias_table.to_csv(
        summary_dir / "chosen_family_bias_estimates.csv", index=False
    )

    conservative_bias = {
        "rE_bias_systematic_fm": float(
            bias_table["rE_RMS_bias_fm"].max()
        ),
        "rM_bias_systematic_fm": float(
            bias_table["rM_RMS_bias_fm"].max()
        ),
    }

    # Reuse the already validated 5% fits for the accepted family.
    fit_rows = chosen_5pct_fit_rows
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

    # Preliminary prescription:
    #   * KM15 5% remains the nominal central fit.
    #   * Model and threshold dependence are NOT treated as independent
    #     systematics.  Both alter which measured points enter the fit, and the
    #     smoke test showed strong correlation, especially for rM.
    #   * Define one joint BH-selection systematic as the largest shift from
    #     the nominal KM15 5% published-errors baseline over all three purity
    #     models and thresholds 3--7%.
    #   * Keep the complete 1--10% scan as a diagnostic.
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
        (published_scan["threshold"].to_numpy(float) >= 0.03)
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
            "[final model analysis] WARNING: skipping invalid 3--7% "
            f"threshold-scan fits when building the joint envelope: {bad}"
        )
    #endif

    if len(valid_local) == 0:
        raise RuntimeError(
            "No valid 3--7% published-error threshold fits remain for the "
            "joint BH-selection envelope."
        )
    #endif

    joint_selection_sys_e = float(np.max(np.abs(
        valid_local["rE_fm"].to_numpy(float) - selection_baseline_e
    )))
    joint_selection_sys_m = float(np.max(np.abs(
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
    method_sys_e = float(np.hypot(joint_selection_sys_e, bias_sys_e))
    method_sys_m = float(np.hypot(joint_selection_sys_m, bias_sys_m))

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
        "closure_min_valid_fraction": float(
            args.radius_bias_min_valid_fraction
        ),
        "rE_fm": float(nominal["rE_fm"]),
        "rE_fit_err_fm": float(nominal["rE_fit_err_fm"]),
        "rE_joint_bh_selection_sys_fm": joint_selection_sys_e,
        "rE_extrapolation_bias_sys_fm": bias_sys_e,
        "rE_method_sys_fm": method_sys_e,
        "rE_fixed5_model_spread_diagnostic_fm": fixed5_model_sys_e,
        "rE_km15_1to10_threshold_envelope_diagnostic_fm": km15_1to10_sys_e,
        "rM_fm": float(nominal["rM_fm"]),
        "rM_fit_err_fm": float(nominal["rM_fit_err_fm"]),
        "rM_joint_bh_selection_sys_fm": joint_selection_sys_m,
        "rM_extrapolation_bias_sys_fm": bias_sys_m,
        "rM_method_sys_fm": method_sys_m,
        "rM_fixed5_model_spread_diagnostic_fm": fixed5_model_sys_m,
        "rM_km15_1to10_threshold_envelope_diagnostic_fm": km15_1to10_sys_m,
        "threshold_systematic_error_mode": "published_errors",
        "note": (
            "Preliminary CLAS-only result. KM15 5% is nominal. "
            "The quoted method systematic combines in quadrature a JOINT "
            "BH-selection envelope (KM15/VGG99/GK16 over 3--7%) and the "
            "closure/extrapolation RMS-bias systematic. Fixed-5% model spread "
            "and the KM15 1--10% threshold envelope are retained as diagnostics "
            "but are not separately added, avoiding double counting of "
            "correlated sample-selection migration. Fit uncertainty includes "
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
        f"  selection sys window: 3--7%, all three purity models "
        f"(invalid scan fits skipped)\n"
        f"  rE = {float(nominal['rE_fm']):.5f} +/- "
        f"{float(nominal['rE_fit_err_fm']):.5f} (fit) +/- "
        f"{method_sys_e:.5f} (method) fm\n"
        f"       joint selection={joint_selection_sys_e:.5f}, "
        f"closure bias={bias_sys_e:.5f}\n"
        f"  rM = {float(nominal['rM_fm']):.5f} +/- "
        f"{float(nominal['rM_fit_err_fm']):.5f} (fit) +/- "
        f"{method_sys_m:.5f} (method) fm\n"
        f"       joint selection={joint_selection_sys_m:.5f}, "
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
            "Full electroproduction model agreement (diagnostic only)"
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
        fig.savefig(outdir / "07_all_point_model_agreement_chi2.png", dpi=180)
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
        fig.savefig(cfg_dir / "01_km15_threshold_stability.png", dpi=180)
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


def run_unified_km15_final_analysis(
        jo_bundle: Dict[str, object],
        defurne_bundle: Dict[str, object],
        georges_bundle: Dict[str, object],
        lee_bundle: Dict[str, object],
        saylor_bundle: Dict[str, object],
        args,
        root_outdir: Path) -> Dict[str, str]:
    """
    Unified KM15 workflow with an INDEPENDENT functional-form determination for
    every requested measurement ensemble.

    This follows the Hayward/Griffioen bias-variance logic: the preferred
    extrapolation form depends on the actual kinematic range and precision of
    the data being fit, so Lee-only, Jo+Lee, and all-four do not inherit one
    another's preferred family.
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

    production_bundles = sort_bundles_chronologically(
        [jo_bundle, defurne_bundle, georges_bundle, lee_bundle]
    )
    ensemble_defs = [
        ("lee2026_only", [lee_bundle], "CLAS12 Lee 2026"),
        ("jo2015_only", [jo_bundle], "CLAS6 Jo 2015"),
        (
            "jo2015_plus_lee2026",
            [jo_bundle, lee_bundle],
            "CLAS6 Jo 2015 + CLAS12 Lee 2026",
        ),
        (
            "all_four",
            production_bundles,
            "Jo 2015 + Defurne 2015 + Georges 2022 + Lee 2026",
        ),
    ]

    print("\n" + "=" * 78)
    print("[UNIFIED KM15 FINAL ANALYSIS]")
    print("[BH-purity prescription] KM15 only")
    print("[functional forms] independently closure-ranked for EACH ensemble")
    print(
        f"[optimizer] initial rE=rM={SACHS_INITIAL_RADIUS_FM:.2f} fm; "
        f"broad slope-equivalent radius guards "
        f"{SACHS_MIN_RADIUS_FM:.2f}--{SACHS_MAX_RADIUS_FM:.2f} fm"
    )
    print("[normalization] correlated scale nuisances included in data and closure")
    print("[Saylor] KM15 diagnostic plus direct Jo comparison; not production")
    print("=" * 78)

    chosen = {}
    family_audit_rows = []

    for tag, cfg_bundles, label in ensemble_defs:
        closure_dir = closure_root / tag
        closure_dir.mkdir(parents=True, exist_ok=True)

        closure_bundles = []
        for bundle in cfg_bundles:
            selected5 = select_bundle_from_external_model(
                bundle, selection, "km15", 0.05
            )
            bcopy = dict(bundle)
            bcopy["set5"] = selected5.copy()
            closure_bundles.append(bcopy)
        #endfor

        ranking_path = (
            closure_dir / "radius_bias_mixed_family_ranking.csv"
        )

        if args.run_radius_bias_study:
            if ranking_path.exists():
                print(
                    f"\n[ensemble closure] {label}: existing closure ranking "
                    "found; reusing it instead of recomputing replicas"
                )
            else:
                print(
                    f"\n[ensemble closure] {label}: running independent "
                    "functional-form determination"
                )
                run_radius_bias_variance_study(
                    closure_bundles, args, closure_dir
                )
            #endif
        #endif
        if not ranking_path.exists():
            raise RuntimeError(
                f"Missing closure ranking for {label}: {ranking_path}. "
                "Run with --run-radius-bias-study "
                "--radius-bias-extended-truths."
            )
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
        ).reset_index(drop=True)
        if len(usable) == 0:
            raise RuntimeError(
                f"No convergence-eligible Sachs family for {label}. "
                "Inspect the ensemble closure before production fitting."
            )
        #endif

        selected_family = None
        for irank, rank_row in usable.iterrows():
            family = str(rank_row["family"])
            specs, counts = _km15_selected_specs_for_bundles(
                cfg_bundles, selection, 0.05
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
                "configuration": tag,
                "configuration_label": label,
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
            })
            # Persist the audit after every attempted production family so a
            # later fatal exception never hides which closure-ranked fits were
            # tried or why they failed.
            pd.DataFrame(family_audit_rows).to_csv(
                summary_dir / "km15_family_selection_by_ensemble.csv",
                index=False,
            )
            if valid:
                selected_family = family
                break
            #endif
        #endfor

        if selected_family is None:
            chosen[tag] = None
            print(
                f"[ensemble closure] {label}: NO valid simultaneous "
                "free-GE/free-GM nominal fit from the closure-ranked families."
            )
            print(
                "[ensemble closure] This ensemble will be recorded as "
                "unresolved and skipped in the free-GE/free-GM threshold "
                "rotation; external-magnetic-input diagnostics will still run "
                "where defined."
            )
            audit_now = pd.DataFrame(family_audit_rows)
            failure_dir = (
                diagnostics_dir / "nominal_family_failures" / tag
            )
            save_nominal_family_failure_diagnostic(
                audit_now,
                configuration=tag,
                configuration_label=label,
                outdir=failure_dir,
            )
        else:
            chosen[tag] = selected_family
            print(
                f"[ensemble closure] {label}: selected family = "
                f"{selected_family}"
            )
        #endif
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

    # Run each threshold scan with that ensemble's own selected family.
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
        scan_thresholds = [0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.10,
                           0.15, 0.20, 0.30, 0.50, 0.75, 1.00]
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
            row = {
                "configuration": tag,
                "configuration_label": label,
                "selected_family": family,
                "threshold": float(threshold),
                "threshold_percent": 100.0 * float(threshold),
                **counts,
                **fit,
            }
            threshold_rows.append(row)
            if abs(float(threshold) - 0.05) < 1.0e-12:
                nominal_rows.append(row)
            #endif
        #endfor

        table = pd.DataFrame(threshold_rows)
        table.to_csv(cfg_dir / "km15_threshold_scan.csv", index=False)

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
        fig.savefig(cfg_dir / "01_km15_threshold_stability.png", dpi=180)
        plt.close(fig)
    #endfor

    nominal = pd.DataFrame(nominal_rows)
    nominal.to_csv(
        rotation_root / "km15_5pct_ensemble_summary.csv", index=False
    )

    print("\n[KM15 ensemble rotation] independently selected families")
    cols = [
        "configuration_label", "selected_family", "N", "chi2_ndof",
        "rE_fm", "rE_fit_err_fm", "rM_fm", "rM_fit_err_fm",
    ]
    if len(nominal):
        present = [c for c in cols if c in nominal.columns]
        print(nominal[present].to_string(index=False))
    else:
        print("[KM15 ensemble rotation] no resolved nominal free-GE/free-GM fits.")
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

    # Saylor: keep both model-based and direct Jo-vs-Saylor diagnostics.
    # For the BH-fit portion use the all-four selected family only as a
    # diagnostic reference; it does not affect any production selection.
    if chosen.get("all_four") is not None:
        run_saylor_recovery_study(
            saylor_bundle,
            args,
            root_outdir,
            family=chosen["all_four"],
        )
    else:
        print(
            "[Saylor recovery] skipped family-referenced BH-fit diagnostic "
            "because all-four free-GE/free-GM extraction was unresolved."
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
    (Q2, xB, |t|, phi) space.  Compare DATA/KM15 rather than raw cross sections,
    so the known kinematic dependence is divided out before the two
    measurements are compared.  This is diagnostic only.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    jo = jo_bundle["all_data"].copy()
    sa = saylor_bundle["all_data"].copy()

    cols = ["Q2", "xB", "t_abs", "phi_deg"]
    for frame in [jo, sa]:
        for col in cols + ["xs", "km15_ep"]:
            frame[col] = pd.to_numeric(frame[col], errors="coerce")
        #endfor
    #endfor
    jo = jo.dropna(subset=cols + ["xs", "km15_ep"]).copy()
    sa = sa.dropna(subset=cols + ["xs", "km15_ep"]).copy()

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
            rs = float(srow["xs"] / srow["km15_ep"])
            rj = float(jrow["xs"] / jrow["km15_ep"])
            rows.append({
                "saylor_index": int(si),
                "jo_index": int(j),
                "match_distance": float(dist),
                "saylor_Q2": float(srow["Q2"]),
                "jo_Q2": float(jrow["Q2"]),
                "saylor_xB": float(srow["xB"]),
                "jo_xB": float(jrow["xB"]),
                "saylor_t_abs": float(srow["t_abs"]),
                "jo_t_abs": float(jrow["t_abs"]),
                "saylor_phi_deg": float(srow["phi_deg"]),
                "jo_phi_deg": float(jrow["phi_deg"]),
                "saylor_data_over_km15": rs,
                "jo_data_over_km15": rj,
                "ratio_saylor_over_jo": rs / rj if rj != 0.0 else np.nan,
            })
        #endfor
    #endfor

    matches = pd.DataFrame(rows)
    matches.to_csv(
        outdir / "02_jo_saylor_matched_kinematics.csv", index=False
    )

    # Show increasingly strict overlap quality without choosing a posteriori
    # kinematic cuts.
    quantiles = [0.25, 0.50, 0.75, 1.00]
    summary = []
    for q in quantiles:
        cut = float(matches["match_distance"].quantile(q))
        part = matches.loc[matches["match_distance"] <= cut]
        summary.append({
            "match_distance_quantile": q,
            "max_match_distance": cut,
            "N": len(part),
            "median_saylor_over_jo": float(
                np.nanmedian(part["ratio_saylor_over_jo"])
            ),
            "rms_saylor_over_jo_minus_one": float(
                np.sqrt(np.nanmean(
                    (part["ratio_saylor_over_jo"] - 1.0)**2
                ))
            ),
        })
    #endfor
    pd.DataFrame(summary).to_csv(
        outdir / "03_jo_saylor_matched_overlap_summary.csv", index=False
    )

    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    ax.scatter(
        matches["match_distance"],
        matches["ratio_saylor_over_jo"],
        s=10, alpha=0.35,
    )
    ax.axhline(1.0, linewidth=0.9, linestyle="--")
    ax.set_xlabel(
        r"nearest-Jo distance in standardized $(Q^2,x_B,|t|,\phi)$ space"
    )
    ax.set_ylabel(
        r"$(\sigma_{\rm Saylor}/\sigma_{\rm KM15})/"
        r"(\sigma_{\rm Jo}/\sigma_{\rm KM15})$"
    )
    ax.set_title("Direct Jo 2015 vs Saylor 2018 matched-kinematics diagnostic")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(
        outdir / "02_jo_saylor_matched_kinematics.png", dpi=180
    )
    plt.close(fig)
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
    print("[production candidates] Jo 2015; Defurne 2015; Georges 2022; Lee 2026")
    print("[BH-purity prescription] KM15")
    print("[Saylor 2018] diagnostic/recoverability study")
    print("=" * 78)

    root_outdir = Path(args.outdir).expanduser().resolve()

    jo = run_clas6_validation(args, return_results=True)
    lee = run_pass1_validation(args, return_results=True)
    defurne = run_halla_defurne_validation(args, return_results=True)
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
        [jo, defurne, georges, lee]
    )
    diagnostic_bundles = sort_bundles_chronologically(
        [jo, defurne, saylor, georges, lee]
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

        run_unified_km15_final_analysis(
            jo,
            defurne,
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

    This mode bypasses the pass-2 CSV and uses:
      * pass-1 stat uncertainty
      * pass-1 per-point systematic uncertainty
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
    print("[CLAS12 PASS-1 VALIDATION MODE]")
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
            "norm_frac": 0.0,
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
        return_parameters: bool = False):
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

    pars = np.asarray([float(m.values[n]) for n in names], dtype=float)
    radius = sachs_family_radius(pars, family)
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

    if return_parameters:
        return radius, valid, pars, float(m.fval)
    #endif
    return radius, valid
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

    candidate_families = [
        "P2", "P3",
        "IP1", "IP2", "IP3",
        "CF2", "CF3",
    ]
    truth_families = [f"P{i}" for i in range(1, 5)]
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
    radius_values = parse_radius_bias_grid(args.radius_bias_radius_grid)
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
    max_attempts = max(5 * target, target + 50)
    rows = []
    replica_rows = []

    for itruth, scenario in enumerate(scenarios):
        gm_truth = MU_N * np.asarray(
            scenario["truth_shape_measured"], dtype=float
        )
        f1_truth = (ge_fixed + tau * gm_truth) * inv
        f2_truth = (gm_truth - ge_fixed) * inv
        central = bh_from_f1f2(A, B, C, f1_truth, f2_truth)
        rtrue = float(scenario["truth_radius_fm"])

        for ifam, family in enumerate(candidate_families):
            successes = []
            attempted = 0
            seed_seq = np.random.SeedSequence([
                int(args.radius_bias_seed), 918273, int(itruth), int(ifam)
            ])
            seeds = seed_seq.spawn(max_attempts)

            while len(successes) < target and attempted < max_attempts:
                rng = np.random.default_rng(
                    int(seeds[attempted].generate_state(
                        1, dtype=np.uint64
                    )[0])
                )
                pseudo = central + rng.normal(0.0, sigma)
                radius, valid, pars, chi2 = fit_neutron_magnetic_family_bh(
                    evaluated,
                    family,
                    override_y=pseudo,
                    statistical_only=True,
                    return_parameters=True,
                )
                accepted = bool(valid and len(successes) < target)
                replica_rows.append({
                    "truth_model": scenario["truth_model"],
                    "truth_group": scenario["truth_group"],
                    "truth_family": scenario["truth_family"],
                    "synthetic_truth": bool(scenario["synthetic"]),
                    "truth_rM_fm": rtrue,
                    "family": family,
                    "attempt": attempted,
                    "accepted_replica": (
                        len(successes) if accepted else np.nan
                    ),
                    "rM_fit_fm": radius,
                    "chi2": chi2,
                    "valid": bool(valid),
                })
                attempted += 1
                if accepted:
                    successes.append(radius)
                #endif
            #endwhile

            arr = np.asarray(successes, dtype=float)
            if len(arr):
                mean = float(np.mean(arr))
                stat = (
                    float(np.std(arr, ddof=1))
                    if len(arr) > 1 else 0.0
                )
                bias = float(mean - rtrue)
                objective = float(math.sqrt(stat**2 + bias**2))
            else:
                mean = stat = bias = objective = np.nan
            #endif

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
                "attempted_replicas": attempted,
                "attempt_efficiency": (
                    len(arr) / attempted if attempted else np.nan
                ),
                "target_reached": bool(len(arr) >= target),
            })
            print(
                f"[neutron-GM-closure] "
                f"truth={str(scenario['truth_model']):28s} "
                f"family={family:3s} "
                f"valid={len(arr)}/{target} attempted={attempted}"
            )
        #endfor
    #endfor

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
        "radius_grid_fm": str(args.radius_bias_radius_grid),
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
        radius, valid, pars, chi2 = fit_neutron_magnetic_family_bh(
            selected,
            family,
            statistical_only=False,
            bh_systematic_fraction=0.05,
            return_parameters=True,
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

    fig, ax = plt.subplots(figsize=(8.0, 5.4))
    ax.plot(
        valid_scan["bh_threshold_percent"],
        valid_scan["rM_fm"],
        marker="o", linestyle="none",
    )
    ax.axhline(
        0.864, linewidth=1.0, linestyle="--",
        label=r"PDG $r_{M,n}=0.864$ fm",
    )
    ax.axvline(5.0, linewidth=0.8, linestyle=":")
    ax.set_xscale("log")
    ax.set_xlabel(r"$|1-R_{\rm BH}^{\rm KM15}|$ threshold (%)")
    ax.set_ylabel(r"$r_{M,n}$ [fm]")
    ax.set_title(
        rf"Neutron magnetic radius with closure-selected {family}"
    )
    ax.grid(alpha=0.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(
        diagdir / "05_neutron_selected_family_radius_threshold.png",
        dpi=180,
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
        help="Legacy CLAS12 Lee 2026 cross-section CSV",
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
            "run model-specific 5% selections, closure-based Sachs-family "
            "selection, production fits, and 1--10% threshold systematics"
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
        default="0.75,0.80,0.85,0.90,0.92",
        help=(
            "Comma-separated synthetic truth radii in fm; the extended study "
            "uses the Cartesian rE x rM grid for every generating family"
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
        help="Nested |1-R_BH| thresholds; nominal paper values are 1-5%",
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
            "Disable the correlated 4.76% target-thickness/absolute-charge "
            "normalization nuisance. It is included by default."
        ),
    )
    p.add_argument(
        "--no-bh-selection-systematic",
        action="store_true",
        help="Disable the paper's additional 1-5% uncorrelated uncertainty",
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
            "The paper defines Fits 1-5 with 1%,2%,3%,4%,5% cuts. "
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
            "The exact paper-method Fit 6-8 stage requires the 5% Set 5."
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
