#!/usr/bin/env python3
#
# dvcs_current_dependence.py
#
# DATA:
#   - Read each DVCS ROOT file
#   - Count events per run using the runnum branch
#   - Read accumulated charge per run from global.csv
#   - Resolve beam current per run using the same hard-coded mapping logic
#   - Aggregate total counts and total charge by (period, current)
#   - Compute counts per nC with Poisson uncertainties
#   - Fit y = m*x + b with weighted least squares
#
# MC:
#   - Scan the DVCSGEN MC directory for gen_*.root and rec_*.root files
#   - Pair generated and reconstructed ROOT files by name
#   - Identify period and current from the filename
#   - Treat "nobkg" as 0 nA
#   - Compute efficiency = N_rec / N_gen
#   - Compute uncertainty = sqrt(N_rec) / N_gen
#   - Fit y = m*x + b with weighted least squares
#
# IMPORTANT OPERATIONAL POINT:
#   The production acceptance correction does NOT use all MC currents.
#   Instead, it uses only one MC current per period:
#
#     Sp18 Inb  -> 50 nA
#     Sp18 Out  -> 45 nA
#     Fa18 Inb  -> 50 nA
#     Fa18 Out  -> 50 nA
#     Sp19 Inb  -> 50 nA
#
#   Therefore, the quantity relevant for efficiency.json is NOT
#     epsilon_data_rel(I) / epsilon_mc_rel(I)
#   but rather
#     epsilon_data_rel(I) / epsilon_mc_rel(I_ref)
#   where I_ref is the single MC current actually used in the production
#   acceptance for that period.
#
#   This script keeps the diagnostic CSVs and the residual-correction table,
#   but the plotting has been streamlined into four integrated figures saved to:
#
#     output/dvcs_current_dependence/integrated/
#
#   Integrated plots produced:
#     (1) six-panel DATA counts/nC vs current
#     (2) six-panel DATA percent of fitted 0 nA vs current
#     (3) six-panel absolute DATA and MC vs current together
#         - data shown as counts/nC
#         - MC shown as rec/gen
#         - this remains a shape-only comparison because the absolute units differ
#     (4) six-panel DATA and MC percent of fitted 0 nA vs current together,
#         with a lower ratio pad in each panel
#
#   Requested plotting conventions:
#     - all x-axes fixed to 0-80 nA
#     - percent panels fixed to 40-120 (%)
#     - MC shown with dashed lines and open-circle markers
#     - DATA shown with solid lines and filled-circle markers
#
# Output figures:
#   output/dvcs_current_dependence/integrated/dvcs_data_counts_per_nc_integrated.png
#   output/dvcs_current_dependence/integrated/dvcs_data_percent_of_zero_integrated.png
#   output/dvcs_current_dependence/integrated/dvcs_absolute_data_mc_integrated.png
#   output/dvcs_current_dependence/integrated/dvcs_percent_data_mc_with_ratio_integrated.png
#
# Diagnostic CSVs:
#   output/dvcs_current_dependence/dvcs_current_dependence_run_table.csv
#   output/dvcs_current_dependence/dvcs_current_dependence_current_table.csv
#   output/dvcs_current_dependence/dvcs_current_dependence_mc_table.csv
#   output/dvcs_current_dependence/dvcs_current_dependence_residual_table.csv
#
# Temporary runtime optional override:
#   --skip_temp_heavy_mc
#     currently disabled in code below unless you uncomment skip_pairs
#

import os
import math
import csv
import argparse
from collections import defaultdict

import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec


# -----------------------------------------------------------------------------
# Input configuration
# -----------------------------------------------------------------------------

CSV_FILE = "/u/home/thayward/clas12_analysis_software/analysis_scripts/dvcs_cross_section/imports/integrated_luminosity/global.csv"

DATA_PERIOD_FILES = [
    ("Sp18 Inb", "rga_sp18_inb", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp18_inb_epgamma.root"),
    ("Sp18 Out", "rga_sp18_out", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp18_out_epgamma.root"),
    ("Fa18 Inb", "rga_fa18_inb", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_inb_epgamma.root"),
    ("Fa18 Out", "rga_fa18_out", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_out_epgamma.root"),
    ("Sp19 Inb", "rga_sp19_inb", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp19_inb_epgamma.root"),
]

PERIOD_ORDER = ["Sp18 Inb", "Sp18 Out", "Fa18 Inb", "Fa18 Out", "Sp19 Inb"]

PERIOD_DISPLAY_FROM_INTERNAL = {
    "rga_sp18_inb": "Sp18 Inb",
    "rga_sp18_out": "Sp18 Out",
    "rga_fa18_inb": "Fa18 Inb",
    "rga_fa18_out": "Fa18 Out",
    "rga_sp19_inb": "Sp19 Inb",
}

PERIOD_INTERNAL_FROM_DISPLAY = {
    "Sp18 Inb": "rga_sp18_inb",
    "Sp18 Out": "rga_sp18_out",
    "Fa18 Inb": "rga_fa18_inb",
    "Fa18 Out": "rga_fa18_out",
    "Sp19 Inb": "rga_sp19_inb",
}

# The single MC current actually used in the production acceptance correction.
MC_REFERENCE_CURRENT = {
    "rga_sp18_inb": 50,
    "rga_sp18_out": 45,
    "rga_fa18_inb": 50,
    "rga_fa18_out": 50,
    "rga_sp19_inb": 50,
}

MC_DIR = "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen"
MC_TREE_NAME = "PhysicsEvents"

OUTPUT_DIR = "output/dvcs_current_dependence"
INTEGRATED_DIR = os.path.join(OUTPUT_DIR, "integrated")

XMIN = 0.0
XMAX = 80.0

PERCENT_YMIN = 40.0
PERCENT_YMAX = 120.0


# -----------------------------------------------------------------------------
# Run -> current maps
# -----------------------------------------------------------------------------

FA18_INB_CURRENT = {
    5418: 5, 5419: 5,

    5335: 40, 5339: 40, 5341: 40,
    5340: 40, 5342: 40, 5343: 40, 5344: 40,

    5032: 45, 5036: 45, 5038: 45, 5039: 45, 5040: 45, 5041: 45,
    5043: 45, 5045: 45, 5052: 45, 5053: 45, 5116: 45, 5117: 45,
    5119: 45, 5120: 45, 5124: 45, 5125: 45, 5126: 45, 5127: 45,
    5139: 45, 5153: 45, 5158: 45, 5162: 45, 5163: 45, 5164: 45,
    5181: 45, 5191: 45, 5193: 45, 5195: 45, 5196: 45, 5197: 45,
    5198: 45, 5199: 45, 5200: 45, 5201: 45, 5202: 45, 5203: 45,
    5204: 45, 5205: 45, 5206: 45, 5208: 45, 5211: 45, 5212: 45,
    5215: 45, 5216: 45, 5219: 45, 5220: 45, 5221: 45, 5222: 45,
    5223: 45, 5230: 45, 5231: 45, 5232: 45, 5233: 45, 5234: 45,
    5235: 45, 5237: 45, 5238: 45, 5247: 45, 5248: 45, 5249: 45,
    5252: 45, 5253: 45, 5257: 45, 5258: 45, 5259: 45, 5261: 45,
    5262: 45, 5303: 45, 5304: 45, 5305: 45, 5306: 45, 5307: 45,
    5310: 45, 5311: 45, 5315: 45, 5317: 45, 5318: 45, 5319: 45,
    5320: 45, 5323: 45, 5324: 45, 5333: 45, 5334: 45, 5346: 45,
    5347: 45, 5349: 45, 5351: 45, 5354: 45, 5355: 45, 5367: 45,

    5046: 45, 5047: 45, 5051: 45,
    5128: 45, 5129: 45, 5130: 45,
    5159: 45, 5160: 45,
    5165: 45, 5166: 45, 5167: 45, 5168: 45, 5169: 45,
    5180: 45, 5182: 45, 5183: 45,
    5190: 45,
    5239: 45,
    5336: 45,

    5356: 50, 5357: 50, 5358: 50, 5359: 50, 5360: 50, 5361: 50,
    5362: 50, 5366: 50,

    5368: 55, 5369: 55, 5372: 55, 5373: 55, 5374: 55, 5375: 55,
    5376: 55, 5377: 55, 5378: 55, 5379: 55, 5380: 55, 5381: 55,
    5382: 55, 5383: 55, 5386: 55, 5390: 55, 5391: 55, 5392: 55,
    5393: 55, 5398: 55, 5400: 55, 5401: 55, 5403: 55, 5404: 55,
    5406: 55, 5407: 55,
}

FA18_OUT_CURRENT = {
    5443: 5,

    5444: 20,

    5423: 40, 5424: 40, 5425: 40, 5426: 40, 5428: 40, 5429: 40,
    5430: 40, 5432: 40, 5434: 40, 5435: 40, 5436: 40, 5437: 40,
    5438: 40, 5440: 40, 5441: 40, 5442: 40, 5445: 40, 5447: 40,
    5449: 40, 5450: 40, 5451: 40, 5452: 40, 5453: 40, 5454: 40,
    5455: 40, 5460: 40, 5464: 40, 5465: 40, 5466: 40, 5467: 40,
    5468: 40, 5469: 40, 5470: 40, 5471: 40, 5472: 40, 5473: 40,
    5474: 40, 5475: 40, 5476: 40, 5478: 40, 5479: 40, 5480: 40,
    5481: 40, 5482: 40, 5483: 40, 5485: 40, 5486: 40, 5487: 40,
    5497: 40, 5498: 40, 5499: 40, 5500: 40, 5504: 40,

    5448: 40, 5495: 40, 5496: 40,

    5507: 50, 5516: 50, 5517: 50, 5518: 50, 5519: 50, 5520: 50,
    5521: 50, 5522: 50, 5523: 50, 5524: 50, 5525: 50, 5526: 50,
    5527: 50, 5528: 50, 5530: 50, 5532: 50, 5533: 50, 5534: 50,
    5535: 50, 5536: 50, 5537: 50, 5538: 50, 5540: 50, 5541: 50,
    5543: 50, 5544: 50, 5545: 50, 5546: 50, 5547: 50, 5548: 50,
    5549: 50, 5550: 50, 5551: 50, 5552: 50, 5555: 50, 5556: 50,
    5557: 50, 5558: 50, 5559: 50, 5562: 50, 5569: 50, 5570: 50,
    5571: 50, 5572: 50, 5573: 50, 5574: 50, 5577: 50, 5578: 50,
    5591: 50, 5592: 50, 5594: 50, 5597: 50, 5598: 50, 5600: 50,
    5601: 50, 5602: 50, 5603: 50, 5604: 50, 5606: 50, 5607: 50,
    5611: 50, 5612: 50, 5613: 50, 5614: 50, 5616: 50, 5618: 50,
    5619: 50, 5624: 50, 5625: 50, 5626: 50, 5627: 50, 5628: 50,
    5629: 50, 5630: 50, 5631: 50, 5632: 50, 5633: 50, 5635: 50,
    5637: 50, 5638: 50, 5639: 50, 5641: 50, 5643: 50, 5644: 50,
    5645: 50, 5646: 50, 5647: 50, 5648: 50, 5649: 50, 5650: 50,
    5651: 50, 5652: 50, 5654: 50, 5655: 50, 5656: 50, 5662: 50,
    5663: 50, 5664: 50, 5665: 50, 5666: 50,
    5615: 50,

    5505: 50, 5567: 50, 5617: 50, 5621: 50, 5623: 50,

    5610: 5,
}


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

def resolve_current(period_label_internal, runnum):
    """
    Resolve beam current (nA) for a given period and run number.

    Returns:
      (True, current_nA) on success
      (False, None) otherwise
    """
    label = period_label_internal.lower()

    if label == "rga_fa18_inb":
        if runnum in FA18_INB_CURRENT:
            return True, FA18_INB_CURRENT[runnum]
        #endif
        return False, None
    #endif

    if label == "rga_fa18_out":
        if runnum in FA18_OUT_CURRENT:
            return True, FA18_OUT_CURRENT[runnum]
        #endif
        return False, None
    #endif

    if label == "rga_sp18_out":
        if 3211 <= runnum <= 3293:
            return True, 30
        #endif
        if 3867 <= runnum <= 3987:
            return True, 45
        #endif
        return False, None
    #endif

    if label == "rga_sp18_inb":
        if runnum == 3418:
            return True, 70
        #endif
        if runnum == 3421 or runnum == 3422:
            return True, 35
        #endif
        if runnum == 3429:
            return True, 50
        #endif

        if 3306 <= runnum <= 3411:
            return True, 35
        #endif

        if 3431 <= runnum <= 4325:
            return True, 50
        #endif

        return False, None
    #endif

    if label == "rga_sp19_inb":
        if runnum == 6616:
            return True, 5
        #endif
        return True, 50
    #endif

    return False, None
#enddef


def f5(val):
    """
    Format a float with 5 digits after the decimal, no scientific notation.
    """
    try:
        return f"{float(val):.5f}"
    except Exception:
        return "nan"
    #endif
#enddef


def weighted_linear_fit(x, y, sy):
    """
    Weighted least squares fit for y = m*x + b.

    Inputs:
      x, y, sy: numpy arrays

    Returns dict:
      m, b, sm, sb, cov_mb, chi2, ndof
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    sy = np.asarray(sy, dtype=float)

    if x.size < 2:
        return {
            "m": float("nan"),
            "b": float("nan"),
            "sm": float("nan"),
            "sb": float("nan"),
            "cov_mb": float("nan"),
            "chi2": float("nan"),
            "ndof": 0,
        }
    #endif

    if np.any(sy <= 0.0):
        raise RuntimeError("Encountered non-positive uncertainty in weighted_linear_fit.")
    #endif

    w = 1.0 / (sy * sy)

    S = np.sum(w)
    Sx = np.sum(w * x)
    Sy = np.sum(w * y)
    Sxx = np.sum(w * x * x)
    Sxy = np.sum(w * x * y)

    D = S * Sxx - Sx * Sx
    if D == 0.0:
        return {
            "m": float("nan"),
            "b": float("nan"),
            "sm": float("nan"),
            "sb": float("nan"),
            "cov_mb": float("nan"),
            "chi2": float("nan"),
            "ndof": 0,
        }
    #endif

    m = (S * Sxy - Sx * Sy) / D
    b = (Sxx * Sy - Sx * Sxy) / D

    var_m = S / D
    var_b = Sxx / D
    cov_mb = -Sx / D

    sm = math.sqrt(var_m) if var_m >= 0.0 else float("nan")
    sb = math.sqrt(var_b) if var_b >= 0.0 else float("nan")

    yfit = m * x + b
    chi2 = np.sum(((y - yfit) / sy) ** 2)
    ndof = int(x.size - 2)

    return {
        "m": m,
        "b": b,
        "sm": sm,
        "sb": sb,
        "cov_mb": cov_mb,
        "chi2": chi2,
        "ndof": ndof,
    }
#enddef


def read_charge_map(csv_file):
    """
    Read global.csv and return:
      run_charge_map[runnum] = charge_nC
    """
    if not os.path.exists(csv_file):
        raise RuntimeError(f"Charge CSV not found: {csv_file}")
    #endif

    run_info_df = pd.read_csv(
        csv_file,
        comment="#",
        header=None,
        names=["runnum", "charge_nC", "col2", "col3", "col4", "col5"],
    )

    if "runnum" not in run_info_df.columns or "charge_nC" not in run_info_df.columns:
        raise RuntimeError("Failed to read required columns from charge CSV.")
    #endif

    run_charge_map = {}

    for _, row in run_info_df.iterrows():
        runnum = int(row["runnum"])
        charge = float(row["charge_nC"])
        run_charge_map[runnum] = charge
    #endfor

    return run_charge_map
#enddef


def read_run_counts_from_root(root_path):
    """
    Read PhysicsEvents/runnum from a ROOT file and return:
      run_counts[runnum] = number of entries in that run
    """
    if not os.path.exists(root_path):
        raise RuntimeError(f"ROOT file not found: {root_path}")
    #endif

    root_file = uproot.open(root_path)

    if "PhysicsEvents" not in root_file:
        raise RuntimeError(f"'PhysicsEvents' tree not found in {root_path}")
    #endif

    tree = root_file["PhysicsEvents"]

    if "runnum" not in tree.keys():
        raise RuntimeError(f"'runnum' branch not found in {root_path}:PhysicsEvents")
    #endif

    runnum_data = tree["runnum"].array(library="np")

    unique_runs, event_counts = np.unique(runnum_data, return_counts=True)

    run_counts = {}
    for run_num, count in zip(unique_runs, event_counts):
        run_counts[int(run_num)] = int(count)
    #endfor

    return run_counts
#enddef


def count_tree_entries(root_path, tree_name):
    """
    Count total entries in a ROOT tree.
    """
    if not os.path.exists(root_path):
        raise RuntimeError(f"ROOT file not found: {root_path}")
    #endif

    root_file = uproot.open(root_path)
    if tree_name not in root_file:
        raise RuntimeError(f"'{tree_name}' tree not found in {root_path}")
    #endif

    tree = root_file[tree_name]
    return int(tree.num_entries)
#enddef


def parse_mc_filename(basename):
    """
    Parse names like:
      gen_dvcsgen_rga_fa18_inb_45nA_10604MeV.root
      rec_dvcsgen_rga_sp18_out_nobkg_10594MeV.root

    Returns:
      kind, period_internal, current_nA
    """
    if not basename.endswith(".root"):
        raise RuntimeError(f"Not a ROOT file: {basename}")
    #endif

    stem = basename[:-5]
    tokens = stem.split("_")

    if len(tokens) != 7:
        raise RuntimeError(f"Unexpected MC filename format: {basename}")
    #endif

    kind = tokens[0]
    if kind != "gen" and kind != "rec":
        raise RuntimeError(f"Unexpected MC file prefix in {basename}")
    #endif

    if tokens[1] != "dvcsgen":
        raise RuntimeError(f"Unexpected MC generator tag in {basename}")
    #endif

    period_internal = "_".join(tokens[2:5])
    current_token = tokens[5]

    if current_token == "nobkg":
        current_nA = 0
    else:
        if not current_token.endswith("nA"):
            raise RuntimeError(f"Could not parse current token in {basename}")
        #endif
        current_nA = int(current_token[:-2])
    #endif

    return kind, period_internal, current_nA
#enddef


def get_mc_reference_current(period_name_or_internal):
    """
    Return the single MC current actually used in the production acceptance
    for this period.
    """
    if period_name_or_internal in MC_REFERENCE_CURRENT:
        return MC_REFERENCE_CURRENT[period_name_or_internal]
    #endif

    if period_name_or_internal in PERIOD_INTERNAL_FROM_DISPLAY:
        internal = PERIOD_INTERNAL_FROM_DISPLAY[period_name_or_internal]
        return MC_REFERENCE_CURRENT[internal]
    #endif

    raise RuntimeError(f"No MC reference current configured for period: {period_name_or_internal}")
#enddef


def add_reference_current_text(ax, period_name):
    """
    Annotate a subplot with the production MC reference current used in the
    acceptance correction.
    """
    ref_current = get_mc_reference_current(period_name)
    text = f"MC ref in acceptance: {ref_current} nA"

    ax.text(
        0.03,
        0.03,
        text,
        transform=ax.transAxes,
        fontsize=8.5,
        verticalalignment="bottom",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.75, edgecolor="0.7"),
    )
#enddef


def build_period_aggregation(period_display_name, period_internal_name, root_path, run_charge_map):
    """
    For one data period:
      - read run counts from ROOT
      - resolve charge and current for each run
      - build run-level rows
      - build current-level aggregation

    Returns:
      run_rows, current_rows
    """
    run_counts = read_run_counts_from_root(root_path)

    missing_charge_runs = []
    nonpositive_charge_runs = []
    unknown_current_runs = []

    run_rows = []
    per_current_counts = defaultdict(int)
    per_current_charge = defaultdict(float)
    per_current_runs = defaultdict(list)

    for runnum in sorted(run_counts.keys()):
        count = run_counts[runnum]

        if runnum not in run_charge_map:
            missing_charge_runs.append(runnum)
            continue
        #endif

        charge = float(run_charge_map[runnum])
        if charge <= 0.0:
            nonpositive_charge_runs.append((runnum, charge))
            continue
        #endif

        ok_cur, current_nA = resolve_current(period_internal_name, runnum)
        if not ok_cur:
            unknown_current_runs.append(runnum)
            continue
        #endif

        rate = float(count) / charge
        rate_err = math.sqrt(float(count)) / charge if count > 0 else 0.0

        run_rows.append({
            "period": period_display_name,
            "period_internal": period_internal_name,
            "runnum": runnum,
            "current_nA": int(current_nA),
            "counts": int(count),
            "charge_nC": float(charge),
            "counts_per_nC": float(rate),
            "counts_per_nC_err": float(rate_err),
        })

        per_current_counts[current_nA] += int(count)
        per_current_charge[current_nA] += float(charge)
        per_current_runs[current_nA].append(runnum)
    #endfor

    if len(missing_charge_runs) > 0 or len(nonpositive_charge_runs) > 0 or len(unknown_current_runs) > 0:
        msg = []
        msg.append(f"Fatal bookkeeping problem while processing {period_display_name}:")
        if len(missing_charge_runs) > 0:
            msg.append(f"  Runs missing from charge CSV ({len(missing_charge_runs)}): {missing_charge_runs}")
        #endif
        if len(nonpositive_charge_runs) > 0:
            msg.append(f"  Runs with non-positive charge ({len(nonpositive_charge_runs)}): {nonpositive_charge_runs}")
        #endif
        if len(unknown_current_runs) > 0:
            msg.append(f"  Runs with no current mapping ({len(unknown_current_runs)}): {unknown_current_runs}")
        #endif
        raise RuntimeError("\n".join(msg))
    #endif

    current_rows = []
    for current_nA in sorted(per_current_counts.keys()):
        total_counts = int(per_current_counts[current_nA])
        total_charge = float(per_current_charge[current_nA])

        if total_charge <= 0.0:
            raise RuntimeError(
                f"Total charge <= 0 for {period_display_name}, current {current_nA} nA"
            )
        #endif

        counts_per_nC = float(total_counts) / total_charge
        counts_per_nC_err = math.sqrt(float(total_counts)) / total_charge if total_counts > 0 else 0.0

        current_rows.append({
            "period": period_display_name,
            "period_internal": period_internal_name,
            "current_nA": int(current_nA),
            "n_runs": len(per_current_runs[current_nA]),
            "counts": int(total_counts),
            "charge_nC": float(total_charge),
            "counts_per_nC": float(counts_per_nC),
            "counts_per_nC_err": float(counts_per_nC_err),
        })
    #endfor

    return run_rows, current_rows
#enddef


def build_mc_aggregation(mc_dir, skip_temp_heavy_mc=False):
    """
    Scan the MC directory, pair gen and rec files, and compute MC efficiencies.

    Returns:
      mc_rows
    """
    if not os.path.isdir(mc_dir):
        raise RuntimeError(f"MC directory not found: {mc_dir}")
    #endif

    skip_pairs = set()
    if skip_temp_heavy_mc:
        # skip_pairs.add(("rga_sp18_out", 45))
        # skip_pairs.add(("rga_fa18_inb", 50))
        # skip_pairs.add(("rga_fa18_out", 50))
        pass
    #endif

    grouped = {}

    for basename in sorted(os.listdir(mc_dir)):
        if not basename.endswith(".root"):
            continue
        #endif

        kind, period_internal, current_nA = parse_mc_filename(basename)

        if period_internal not in PERIOD_DISPLAY_FROM_INTERNAL:
            raise RuntimeError(f"Unknown MC period in filename: {basename}")
        #endif

        if (period_internal, current_nA) in skip_pairs:
            continue
        #endif

        key = (period_internal, current_nA)
        if key not in grouped:
            grouped[key] = {"gen": None, "rec": None}
        #endif

        full_path = os.path.join(mc_dir, basename)

        if grouped[key][kind] is not None:
            raise RuntimeError(f"Duplicate MC file for {key}: {basename}")
        #endif

        grouped[key][kind] = full_path
    #endfor

    mc_rows = []

    for key in sorted(grouped.keys(), key=lambda item: (PERIOD_ORDER.index(PERIOD_DISPLAY_FROM_INTERNAL[item[0]]), item[1])):
        period_internal, current_nA = key
        period_display_name = PERIOD_DISPLAY_FROM_INTERNAL[period_internal]

        gen_path = grouped[key]["gen"]
        rec_path = grouped[key]["rec"]

        if gen_path is None or rec_path is None:
            raise RuntimeError(
                f"Missing gen/rec MC pair for period={period_internal}, current={current_nA} nA"
            )
        #endif

        n_gen = count_tree_entries(gen_path, MC_TREE_NAME)
        n_rec = count_tree_entries(rec_path, MC_TREE_NAME)

        if n_gen <= 0:
            raise RuntimeError(
                f"Generated MC count <= 0 for period={period_internal}, current={current_nA} nA"
            )
        #endif

        efficiency = float(n_rec) / float(n_gen)
        efficiency_err = math.sqrt(float(n_rec)) / float(n_gen) if n_rec > 0 else 0.0

        mc_rows.append({
            "period": period_display_name,
            "period_internal": period_internal,
            "current_nA": int(current_nA),
            "n_gen": int(n_gen),
            "n_rec": int(n_rec),
            "efficiency": float(efficiency),
            "efficiency_err": float(efficiency_err),
            "gen_file": gen_path,
            "rec_file": rec_path,
        })
    #endfor

    return mc_rows
#enddef


def period_points_from_current_rows(current_rows, period_name):
    """
    Extract x, y, sy arrays for one data period from the aggregated current rows.
    """
    rows = [r for r in current_rows if r["period"] == period_name]
    rows = sorted(rows, key=lambda r: r["current_nA"])

    if len(rows) == 0:
        raise RuntimeError(f"No current rows found for period {period_name}")
    #endif

    x = np.asarray([float(r["current_nA"]) for r in rows], dtype=float)
    y = np.asarray([float(r["counts_per_nC"]) for r in rows], dtype=float)
    sy = np.asarray([float(r["counts_per_nC_err"]) for r in rows], dtype=float)

    return x, y, sy, rows
#enddef


def period_points_from_mc_rows(mc_rows, period_name):
    """
    Extract x, y, sy arrays for one MC period from the MC table.
    """
    rows = [r for r in mc_rows if r["period"] == period_name]
    rows = sorted(rows, key=lambda r: r["current_nA"])

    if len(rows) == 0:
        raise RuntimeError(f"No MC rows found for period {period_name}")
    #endif

    x = np.asarray([float(r["current_nA"]) for r in rows], dtype=float)
    y = np.asarray([float(r["efficiency"]) for r in rows], dtype=float)
    sy = np.asarray([float(r["efficiency_err"]) for r in rows], dtype=float)

    return x, y, sy, rows
#enddef


def write_run_table_csv(path, run_rows):
    """
    Write a diagnostic run-level CSV.
    """
    fieldnames = [
        "period",
        "period_internal",
        "runnum",
        "current_nA",
        "counts",
        "charge_nC",
        "counts_per_nC",
        "counts_per_nC_err",
    ]

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in run_rows:
            writer.writerow(row)
        #endfor
    #endwith
#enddef


def write_current_table_csv(path, current_rows):
    """
    Write a diagnostic current-level CSV for data.
    """
    fieldnames = [
        "period",
        "period_internal",
        "current_nA",
        "n_runs",
        "counts",
        "charge_nC",
        "counts_per_nC",
        "counts_per_nC_err",
    ]

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in current_rows:
            writer.writerow(row)
        #endfor
    #endwith
#enddef


def write_mc_table_csv(path, mc_rows):
    """
    Write a diagnostic MC table.
    """
    fieldnames = [
        "period",
        "period_internal",
        "current_nA",
        "n_gen",
        "n_rec",
        "efficiency",
        "efficiency_err",
        "gen_file",
        "rec_file",
    ]

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in mc_rows:
            writer.writerow(row)
        #endfor
    #endwith
#enddef


def write_residual_table_csv(path, residual_rows):
    """
    Write the operational residual current-dependence correction table.

    Stored quantity:
      epsilon_resid(I) = epsilon_data_rel(I) / epsilon_mc_rel(I_ref)

    cross_sections.cpp should then apply:
      1 / epsilon_resid(I)
    """
    fieldnames = [
        "period",
        "period_internal",
        "data_current_nA",
        "mc_reference_current_nA",
        "data_rel_percent",
        "data_rel_fraction",
        "mc_rel_at_reference_percent",
        "mc_rel_at_reference_fraction",
        "residual_percent",
        "residual_fraction",
        "applied_scale_percent",
        "applied_scale_fraction",
    ]

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in residual_rows:
            writer.writerow(row)
        #endfor
    #endwith
#enddef


def compute_percent_curve(xfit, fit_result):
    """
    Convert y = m*x + b into percent-of-intercept:
      100 * (m*x + b) / b

    The uncertainty band is the same slope-only band used elsewhere:
      lower: 100 * (((m - sm) * x + b) / b)
      upper: 100 * (((m + sm) * x + b) / b)

    So we are varying only the slope by +/- 1 sigma while keeping the
    intercept fixed.
    """
    m = fit_result["m"]
    b = fit_result["b"]
    sm = fit_result["sm"]

    pct_fit = 100.0 * ((m * xfit + b) / b)
    pct_fit_lo = 100.0 * (((m - sm) * xfit + b) / b)
    pct_fit_hi = 100.0 * (((m + sm) * xfit + b) / b)

    return pct_fit, pct_fit_lo, pct_fit_hi
#enddef


def compute_percent_slope_and_error(fit_result):
    """
    For p = 100*m/b, propagate uncertainty using m, b, cov(m,b).
    """
    m = fit_result["m"]
    b = fit_result["b"]
    sm = fit_result["sm"]
    sb = fit_result["sb"]
    cov_mb = fit_result["cov_mb"]

    slope_pct = 100.0 * (m / b)

    dp_dm = 100.0 / b
    dp_db = -100.0 * m / (b * b)
    var_p = dp_dm * dp_dm * (sm * sm) + dp_db * dp_db * (sb * sb) + 2.0 * dp_dm * dp_db * cov_mb
    slope_pct_err = math.sqrt(var_p) if var_p >= 0.0 else float("nan")

    return slope_pct, slope_pct_err
#enddef


def compute_relative_value_at_current(current_nA, fit_result):
    """
    Compute relative survival at a given current from the fitted line:
      rel(I) = (m*I + b) / b
    """
    m = fit_result["m"]
    b = fit_result["b"]

    return (m * float(current_nA) + b) / b
#enddef


def build_residual_correction_rows(all_current_rows, data_fit_results, mc_fit_results):
    """
    Build the operational residual correction table for the case where the
    production acceptance uses only one MC current per period.

    For each period and each DATA current I:
      epsilon_data_rel(I) = R_data(I) / R_data(0)
      epsilon_mc_rel(I_ref) = epsilon_mc(I_ref) / epsilon_mc(0)

      epsilon_resid(I) = epsilon_data_rel(I) / epsilon_mc_rel(I_ref)

    The scale actually applied to the cross section is:
      1 / epsilon_resid(I)
    """
    residual_rows = []

    for period in PERIOD_ORDER:
        period_internal = PERIOD_INTERNAL_FROM_DISPLAY[period]
        ref_current = get_mc_reference_current(period_internal)

        data_fit = data_fit_results[period]
        mc_fit = mc_fit_results[period]

        data_rows = [r for r in all_current_rows if r["period"] == period]
        data_rows = sorted(data_rows, key=lambda r: r["current_nA"])

        mc_rel_at_ref = compute_relative_value_at_current(ref_current, mc_fit)

        for row in data_rows:
            current_nA = int(row["current_nA"])
            data_rel = compute_relative_value_at_current(current_nA, data_fit)
            residual = data_rel / mc_rel_at_ref
            applied_scale = 1.0 / residual

            residual_rows.append({
                "period": period,
                "period_internal": period_internal,
                "data_current_nA": current_nA,
                "mc_reference_current_nA": int(ref_current),
                "data_rel_percent": 100.0 * data_rel,
                "data_rel_fraction": data_rel,
                "mc_rel_at_reference_percent": 100.0 * mc_rel_at_ref,
                "mc_rel_at_reference_fraction": mc_rel_at_ref,
                "residual_percent": 100.0 * residual,
                "residual_fraction": residual,
                "applied_scale_percent": 100.0 * applied_scale,
                "applied_scale_fraction": applied_scale,
            })
        #endfor
    #endfor

    return residual_rows
#enddef


def make_absolute_band(fit_result, xfit):
    """
    Absolute y = m*x + b band with slope varied by +/- sm.
    """
    yfit = fit_result["m"] * xfit + fit_result["b"]
    yfit_lo = (fit_result["m"] - fit_result["sm"]) * xfit + fit_result["b"]
    yfit_hi = (fit_result["m"] + fit_result["sm"]) * xfit + fit_result["b"]
    return yfit, yfit_lo, yfit_hi
#enddef


def prepare_period_cache(all_current_rows, mc_rows, data_fit_results, mc_fit_results):
    """
    Build a compact cache used by the plotting code, so the plotting section
    stays streamlined and does not repeatedly recompute the same arrays.
    """
    cache = {}

    for period in PERIOD_ORDER:
        xd, yd, syd, rowsd = period_points_from_current_rows(all_current_rows, period)
        xm, ym, sym, rowsm = period_points_from_mc_rows(mc_rows, period)

        frd = data_fit_results[period]
        frm = mc_fit_results[period]

        cache[period] = {
            "xd": xd,
            "yd": yd,
            "syd": syd,
            "rowsd": rowsd,
            "xm": xm,
            "ym": ym,
            "sym": sym,
            "rowsm": rowsm,
            "frd": frd,
            "frm": frm,
        }
    #endfor

    return cache
#enddef


def set_common_style(ax):
    """
    Apply common style to an axis.
    """
    ax.grid(True, alpha=0.3)
#enddef


def draw_data_absolute_panel(ax, period, cache, color, xfit, show_legend=True):
    """
    Draw one panel for DATA counts/nC vs current.
    """
    frd = cache[period]["frd"]
    xd = cache[period]["xd"]
    yd = cache[period]["yd"]
    syd = cache[period]["syd"]

    data_fit, data_fit_lo, data_fit_hi = make_absolute_band(frd, xfit)

    ax.fill_between(xfit, data_fit_lo, data_fit_hi, color=color, alpha=0.20, linewidth=0)
    ax.errorbar(
        xd,
        yd,
        yerr=syd,
        fmt="o",
        capsize=3,
        color=color,
        label=f"m={f5(frd['m'])} +/- {f5(frd['sm'])}\n b={f5(frd['b'])} +/- {f5(frd['sb'])}",
    )
    ax.plot(xfit, data_fit, color=color)

    ax.set_title(period)
    ax.set_xlim(XMIN, XMAX)
    ax.set_xlabel("Beam current (nA)")
    ax.set_ylabel("DVCS counts per charge (1/nC)")
    set_common_style(ax)
    if show_legend:
        ax.legend(frameon=True, fontsize=9)
    #endif
    add_reference_current_text(ax, period)
#enddef


def draw_data_percent_panel(ax, period, cache, color, xfit, show_legend=True):
    """
    Draw one panel for DATA percent of fitted 0 nA vs current.
    """
    frd = cache[period]["frd"]
    xd = cache[period]["xd"]
    yd = cache[period]["yd"]
    syd = cache[period]["syd"]

    data_pct = 100.0 * (yd / frd["b"])
    data_pct_err = 100.0 * np.sqrt((syd / frd["b"]) ** 2 + ((yd * frd["sb"]) / (frd["b"] * frd["b"])) ** 2)
    data_pct_fit, data_pct_fit_lo, data_pct_fit_hi = compute_percent_curve(xfit, frd)
    slope_pct, slope_pct_err = compute_percent_slope_and_error(frd)

    ax.fill_between(xfit, data_pct_fit_lo, data_pct_fit_hi, color=color, alpha=0.20, linewidth=0)
    ax.errorbar(
        xd,
        data_pct,
        yerr=data_pct_err,
        fmt="o",
        capsize=3,
        color=color,
        label=f"slope={f5(slope_pct)} +/- {f5(slope_pct_err)} (%/nA)",
    )
    ax.plot(xfit, data_pct_fit, color=color)

    ax.set_title(period)
    ax.set_xlim(XMIN, XMAX)
    ax.set_ylim(PERCENT_YMIN, PERCENT_YMAX)
    ax.set_xlabel("Beam current (nA)")
    ax.set_ylabel("Efficiency relative to fitted 0 nA (%)")
    set_common_style(ax)
    if show_legend:
        ax.legend(frameon=True, fontsize=9)
    #endif
    add_reference_current_text(ax, period)
#enddef


def draw_absolute_overlay_panel(ax, period, cache, color, xfit, show_legend=True):
    """
    Draw one panel for absolute DATA and MC together.
    DATA: counts/nC
    MC:   rec/gen
    """
    frd = cache[period]["frd"]
    frm = cache[period]["frm"]

    xd = cache[period]["xd"]
    yd = cache[period]["yd"]
    syd = cache[period]["syd"]

    xm = cache[period]["xm"]
    ym = cache[period]["ym"]
    sym = cache[period]["sym"]

    data_fit, data_fit_lo, data_fit_hi = make_absolute_band(frd, xfit)
    mc_fit, mc_fit_lo, mc_fit_hi = make_absolute_band(frm, xfit)

    ax.fill_between(xfit, data_fit_lo, data_fit_hi, color=color, alpha=0.14, linewidth=0)
    ax.errorbar(
        xd,
        yd,
        yerr=syd,
        fmt="o",
        capsize=3,
        color=color,
        label="Data (counts/nC)",
    )
    ax.plot(xfit, data_fit, color=color, linestyle="-")

    ax.fill_between(xfit, mc_fit_lo, mc_fit_hi, color=color, alpha=0.08, linewidth=0)
    ax.errorbar(
        xm,
        ym,
        yerr=sym,
        fmt="o",
        capsize=3,
        linestyle="none",
        markerfacecolor="none",
        markeredgecolor=color,
        ecolor=color,
        color=color,
        label="MC (rec/gen)",
    )
    ax.plot(xfit, mc_fit, color=color, linestyle="--")

    ax.set_title(period)
    ax.set_xlim(XMIN, XMAX)
    ax.set_xlabel("Beam current (nA)")
    ax.set_ylabel("Absolute value")
    set_common_style(ax)
    if show_legend:
        ax.legend(frameon=True, fontsize=9)
    #endif
    add_reference_current_text(ax, period)
#enddef


def evaluate_percent_fit_and_band(xvals, fit_result):
    """
    Evaluate percent curve and slope-only percent band at specific x values.
    """
    xvals = np.asarray(xvals, dtype=float)

    pct = 100.0 * ((fit_result["m"] * xvals + fit_result["b"]) / fit_result["b"])
    pct_lo = 100.0 * ((((fit_result["m"] - fit_result["sm"]) * xvals) + fit_result["b"]) / fit_result["b"])
    pct_hi = 100.0 * ((((fit_result["m"] + fit_result["sm"]) * xvals) + fit_result["b"]) / fit_result["b"])

    return pct, pct_lo, pct_hi
#enddef


def draw_percent_overlay_with_ratio(top_ax, bot_ax, period, cache, color, xfit, show_legend=True):
    """
    Draw one panel with:
      - top: DATA and MC percent of fitted 0 nA
      - bottom: ratio = DATA / MC_fit_evaluated_at_DATA_currents

    This follows the style the user requested:
      - larger top panel
      - smaller ratio panel underneath
      - MC dashed and open circles
      - DATA solid and filled circles

    Ratio points are evaluated at the DATA current values:
      ratio_i = DATA_percent(I_i) / MC_percent_fit(I_i)

    Ratio band is built from the fit bands:
      ratio_low  = DATA_fit_low / MC_fit_high
      ratio_high = DATA_fit_high / MC_fit_low
    """
    frd = cache[period]["frd"]
    frm = cache[period]["frm"]

    xd = cache[period]["xd"]
    yd = cache[period]["yd"]
    syd = cache[period]["syd"]

    xm = cache[period]["xm"]
    ym = cache[period]["ym"]
    sym = cache[period]["sym"]

    data_pct = 100.0 * (yd / frd["b"])
    data_pct_err = 100.0 * np.sqrt((syd / frd["b"]) ** 2 + ((yd * frd["sb"]) / (frd["b"] * frd["b"])) ** 2)
    data_pct_fit, data_pct_fit_lo, data_pct_fit_hi = compute_percent_curve(xfit, frd)

    mc_pct = 100.0 * (ym / frm["b"])
    mc_pct_err = 100.0 * np.sqrt((sym / frm["b"]) ** 2 + ((ym * frm["sb"]) / (frm["b"] * frm["b"])) ** 2)
    mc_pct_fit, mc_pct_fit_lo, mc_pct_fit_hi = compute_percent_curve(xfit, frm)

    # -------------------------
    # Top panel
    # -------------------------

    top_ax.fill_between(xfit, data_pct_fit_lo, data_pct_fit_hi, color=color, alpha=0.12, linewidth=0)
    top_ax.errorbar(
        xd,
        data_pct,
        yerr=data_pct_err,
        fmt="o",
        capsize=3,
        color=color,
        label="Data",
    )
    top_ax.plot(xfit, data_pct_fit, color=color, linestyle="-")

    top_ax.fill_between(xfit, mc_pct_fit_lo, mc_pct_fit_hi, color=color, alpha=0.08, linewidth=0)
    top_ax.errorbar(
        xm,
        mc_pct,
        yerr=mc_pct_err,
        fmt="o",
        capsize=3,
        linestyle="none",
        markerfacecolor="none",
        markeredgecolor=color,
        ecolor=color,
        color=color,
        label="MC",
    )
    top_ax.plot(xfit, mc_pct_fit, color=color, linestyle="--")

    top_ax.set_title(period)
    top_ax.set_xlim(XMIN, XMAX)
    top_ax.set_ylim(PERCENT_YMIN, PERCENT_YMAX)
    top_ax.set_ylabel("Efficiency relative to fitted 0 nA (%)")
    set_common_style(top_ax)
    top_ax.tick_params(labelbottom=False)
    if show_legend:
        top_ax.legend(frameon=True, fontsize=9, loc="upper right")
    #endif
    add_reference_current_text(top_ax, period)

    # -------------------------
    # Ratio panel
    # -------------------------

    mc_pct_at_xd, mc_pct_at_xd_lo, mc_pct_at_xd_hi = evaluate_percent_fit_and_band(xd, frm)
    ratio_pts = data_pct / mc_pct_at_xd

    # Approximate point uncertainty:
    #   dominant visible uncertainty is DATA point uncertainty combined with the
    #   MC fit-band half-width at the same x.
    mc_pct_sigma_at_xd = 0.5 * np.abs(mc_pct_at_xd_hi - mc_pct_at_xd_lo)
    ratio_err = ratio_pts * np.sqrt((data_pct_err / data_pct) ** 2 + (mc_pct_sigma_at_xd / mc_pct_at_xd) ** 2)

    ratio_fit = data_pct_fit / mc_pct_fit
    ratio_fit_lo = data_pct_fit_lo / mc_pct_fit_hi
    ratio_fit_hi = data_pct_fit_hi / mc_pct_fit_lo

    bot_ax.fill_between(xfit, ratio_fit_lo, ratio_fit_hi, color=color, alpha=0.12, linewidth=0)
    bot_ax.errorbar(
        xd,
        ratio_pts,
        yerr=ratio_err,
        fmt="o",
        capsize=3,
        color=color,
    )
    bot_ax.plot(xfit, ratio_fit, color=color, linestyle="-")
    bot_ax.axhline(1.0, color="gray", linestyle="--", linewidth=1.0)

    bot_ax.set_xlim(XMIN, XMAX)
    bot_ax.set_xlabel("Beam current (nA)")
    bot_ax.set_ylabel("Data/MC")
    set_common_style(bot_ax)

    finite_ratio = ratio_pts[np.isfinite(ratio_pts)]
    finite_ratio_err = ratio_err[np.isfinite(ratio_err)]

    if finite_ratio.size > 0:
        rmin = np.min(finite_ratio - finite_ratio_err)
        rmax = np.max(finite_ratio + finite_ratio_err)
        pad = max(0.01, 0.15 * (rmax - rmin))
        if rmax - rmin < 0.03:
            pad = max(pad, 0.01)
        #endif
        bot_ax.set_ylim(rmin - pad, rmax + pad)
    else:
        bot_ax.set_ylim(0.8, 1.2)
    #endif
#enddef


def plot_six_panel_simple(period_order, cache, period_color, drawer, outfile, figsize=(18, 10)):
    """
    Draw a simple 2x3 figure with one axis per slot.
    """
    fig, axs = plt.subplots(2, 3, figsize=figsize, constrained_layout=True)
    axs = axs.flatten()

    for i, period in enumerate(period_order):
        drawer(axs[i], period, cache, period_color[period], np.linspace(XMIN, XMAX, 300), True)
    #endfor

    # Overlay panel
    axc = axs[5]
    axc.clear()
    axc.set_title("All periods (overlay)")
    axc.set_xlim(XMIN, XMAX)
    set_common_style(axc)

    xfit = np.linspace(XMIN, XMAX, 300)

    if drawer == draw_data_absolute_panel:
        axc.set_xlabel("Beam current (nA)")
        axc.set_ylabel("DVCS counts per charge (1/nC)")
        for period in period_order:
            c = period_color[period]
            frd = cache[period]["frd"]
            xd = cache[period]["xd"]
            yd = cache[period]["yd"]
            syd = cache[period]["syd"]
            data_fit, data_fit_lo, data_fit_hi = make_absolute_band(frd, xfit)

            axc.fill_between(xfit, data_fit_lo, data_fit_hi, color=c, alpha=0.10, linewidth=0)
            axc.errorbar(xd, yd, yerr=syd, fmt="o", capsize=3, color=c, label=f"{period}")
            axc.plot(xfit, data_fit, color=c)
        #endfor
    elif drawer == draw_data_percent_panel:
        axc.set_xlabel("Beam current (nA)")
        axc.set_ylabel("Efficiency relative to fitted 0 nA (%)")
        axc.set_ylim(PERCENT_YMIN, PERCENT_YMAX)
        for period in period_order:
            c = period_color[period]
            frd = cache[period]["frd"]
            xd = cache[period]["xd"]
            yd = cache[period]["yd"]
            syd = cache[period]["syd"]

            data_pct = 100.0 * (yd / frd["b"])
            data_pct_err = 100.0 * np.sqrt((syd / frd["b"]) ** 2 + ((yd * frd["sb"]) / (frd["b"] * frd["b"])) ** 2)
            data_pct_fit, data_pct_fit_lo, data_pct_fit_hi = compute_percent_curve(xfit, frd)

            axc.fill_between(xfit, data_pct_fit_lo, data_pct_fit_hi, color=c, alpha=0.10, linewidth=0)
            axc.errorbar(xd, data_pct, yerr=data_pct_err, fmt="o", capsize=3, color=c, label=f"{period}")
            axc.plot(xfit, data_pct_fit, color=c)
        #endfor
    elif drawer == draw_absolute_overlay_panel:
        axc.set_xlabel("Beam current (nA)")
        axc.set_ylabel("Absolute value")
        for period in period_order:
            c = period_color[period]

            frd = cache[period]["frd"]
            frm = cache[period]["frm"]

            xd = cache[period]["xd"]
            yd = cache[period]["yd"]
            syd = cache[period]["syd"]

            xm = cache[period]["xm"]
            ym = cache[period]["ym"]
            sym = cache[period]["sym"]

            data_fit, data_fit_lo, data_fit_hi = make_absolute_band(frd, xfit)
            mc_fit, mc_fit_lo, mc_fit_hi = make_absolute_band(frm, xfit)

            axc.fill_between(xfit, data_fit_lo, data_fit_hi, color=c, alpha=0.08, linewidth=0)
            axc.errorbar(xd, yd, yerr=syd, fmt="o", capsize=3, color=c, label=f"{period} data")
            axc.plot(xfit, data_fit, color=c, linestyle="-")

            axc.fill_between(xfit, mc_fit_lo, mc_fit_hi, color=c, alpha=0.05, linewidth=0)
            axc.errorbar(
                xm,
                ym,
                yerr=sym,
                fmt="o",
                capsize=3,
                linestyle="none",
                markerfacecolor="none",
                markeredgecolor=c,
                ecolor=c,
                color=c,
                label=f"{period} MC",
            )
            axc.plot(xfit, mc_fit, color=c, linestyle="--")
        #endfor
    else:
        raise RuntimeError("Unknown drawer passed to plot_six_panel_simple.")
    #endif

    axc.legend(frameon=True, fontsize=8)
    fig.savefig(outfile, dpi=200)
    plt.close(fig)
#enddef


def plot_percent_overlay_with_ratio(period_order, cache, period_color, outfile):
    """
    Draw the requested 2x3 integrated figure where each slot contains:
      - top panel: DATA and MC percent of fitted 0 nA
      - bottom panel: ratio Data/MC
    """
    fig = plt.figure(figsize=(18, 11), constrained_layout=False)
    outer = GridSpec(2, 3, figure=fig, wspace=0.22, hspace=0.28)

    xfit = np.linspace(XMIN, XMAX, 300)

    # Five period panels
    for i, period in enumerate(period_order):
        row = i // 3
        col = i % 3

        inner = GridSpecFromSubplotSpec(
            2,
            1,
            subplot_spec=outer[row, col],
            height_ratios=[3.2, 1.0],
            hspace=0.0,
        )

        top_ax = fig.add_subplot(inner[0, 0])
        bot_ax = fig.add_subplot(inner[1, 0], sharex=top_ax)

        draw_percent_overlay_with_ratio(
            top_ax,
            bot_ax,
            period,
            cache,
            period_color[period],
            xfit,
            show_legend=True,
        )
    #endfor

    # Sixth overlay panel
    inner = GridSpecFromSubplotSpec(
        2,
        1,
        subplot_spec=outer[1, 2],
        height_ratios=[3.2, 1.0],
        hspace=0.0,
    )

    top_ax = fig.add_subplot(inner[0, 0])
    bot_ax = fig.add_subplot(inner[1, 0], sharex=top_ax)

    top_ax.set_title("All periods (overlay)")
    top_ax.set_xlim(XMIN, XMAX)
    top_ax.set_ylim(PERCENT_YMIN, PERCENT_YMAX)
    top_ax.set_ylabel("Efficiency relative to fitted 0 nA (%)")
    set_common_style(top_ax)
    top_ax.tick_params(labelbottom=False)

    bot_ax.set_xlim(XMIN, XMAX)
    bot_ax.set_xlabel("Beam current (nA)")
    bot_ax.set_ylabel("Data/MC")
    set_common_style(bot_ax)
    bot_ax.axhline(1.0, color="gray", linestyle="--", linewidth=1.0)

    all_ratio_y = []
    all_ratio_err = []

    for period in period_order:
        c = period_color[period]
        frd = cache[period]["frd"]
        frm = cache[period]["frm"]

        xd = cache[period]["xd"]
        yd = cache[period]["yd"]
        syd = cache[period]["syd"]

        xm = cache[period]["xm"]
        ym = cache[period]["ym"]
        sym = cache[period]["sym"]

        data_pct = 100.0 * (yd / frd["b"])
        data_pct_err = 100.0 * np.sqrt((syd / frd["b"]) ** 2 + ((yd * frd["sb"]) / (frd["b"] * frd["b"])) ** 2)
        data_pct_fit, data_pct_fit_lo, data_pct_fit_hi = compute_percent_curve(xfit, frd)

        mc_pct = 100.0 * (ym / frm["b"])
        mc_pct_err = 100.0 * np.sqrt((sym / frm["b"]) ** 2 + ((ym * frm["sb"]) / (frm["b"] * frm["b"])) ** 2)
        mc_pct_fit, mc_pct_fit_lo, mc_pct_fit_hi = compute_percent_curve(xfit, frm)

        top_ax.fill_between(xfit, data_pct_fit_lo, data_pct_fit_hi, color=c, alpha=0.08, linewidth=0)
        top_ax.errorbar(xd, data_pct, yerr=data_pct_err, fmt="o", capsize=3, color=c, label=f"{period} data")
        top_ax.plot(xfit, data_pct_fit, color=c, linestyle="-")

        top_ax.fill_between(xfit, mc_pct_fit_lo, mc_pct_fit_hi, color=c, alpha=0.05, linewidth=0)
        top_ax.errorbar(
            xm,
            mc_pct,
            yerr=mc_pct_err,
            fmt="o",
            capsize=3,
            linestyle="none",
            markerfacecolor="none",
            markeredgecolor=c,
            ecolor=c,
            color=c,
            label=f"{period} MC",
        )
        top_ax.plot(xfit, mc_pct_fit, color=c, linestyle="--")

        mc_pct_at_xd, mc_pct_at_xd_lo, mc_pct_at_xd_hi = evaluate_percent_fit_and_band(xd, frm)
        ratio_pts = data_pct / mc_pct_at_xd
        mc_pct_sigma_at_xd = 0.5 * np.abs(mc_pct_at_xd_hi - mc_pct_at_xd_lo)
        ratio_err = ratio_pts * np.sqrt((data_pct_err / data_pct) ** 2 + (mc_pct_sigma_at_xd / mc_pct_at_xd) ** 2)

        ratio_fit = data_pct_fit / mc_pct_fit
        ratio_fit_lo = data_pct_fit_lo / mc_pct_fit_hi
        ratio_fit_hi = data_pct_fit_hi / mc_pct_fit_lo

        bot_ax.fill_between(xfit, ratio_fit_lo, ratio_fit_hi, color=c, alpha=0.08, linewidth=0)
        bot_ax.errorbar(xd, ratio_pts, yerr=ratio_err, fmt="o", capsize=3, color=c)
        bot_ax.plot(xfit, ratio_fit, color=c, linestyle="-")

        all_ratio_y.extend(list(ratio_pts))
        all_ratio_err.extend(list(ratio_err))
    #endfor

    if len(all_ratio_y) > 0:
        all_ratio_y = np.asarray(all_ratio_y, dtype=float)
        all_ratio_err = np.asarray(all_ratio_err, dtype=float)

        rmin = np.min(all_ratio_y - all_ratio_err)
        rmax = np.max(all_ratio_y + all_ratio_err)
        pad = max(0.01, 0.15 * (rmax - rmin))
        if rmax - rmin < 0.03:
            pad = max(pad, 0.01)
        #endif
        bot_ax.set_ylim(rmin - pad, rmax + pad)
    else:
        bot_ax.set_ylim(0.8, 1.2)
    #endif

    top_ax.legend(frameon=True, fontsize=7.5, ncol=1)
    fig.tight_layout()
    fig.savefig(outfile, dpi=200, bbox_inches="tight")
    plt.close(fig)
#enddef


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--skip_temp_heavy_mc",
        action="store_true",
        help="Temporarily skip selected MC points if skip_pairs is uncommented in build_mc_aggregation.",
    )
    args = parser.parse_args()

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(INTEGRATED_DIR, exist_ok=True)

    # -------------------------------------------------------------------------
    # Read charge map once.
    # -------------------------------------------------------------------------
    print("")
    print("Reading charge map...")
    run_charge_map = read_charge_map(CSV_FILE)
    print(f"Loaded charge entries for {len(run_charge_map)} runs from:")
    print(f"  {CSV_FILE}")

    # -------------------------------------------------------------------------
    # Build data run-level and current-level tables.
    # -------------------------------------------------------------------------
    all_run_rows = []
    all_current_rows = []

    print("")
    print("Processing DATA ROOT files and aggregating by current...")
    for period_display_name, period_internal_name, root_path in DATA_PERIOD_FILES:
        print("")
        print("=" * 90)
        print(f"DATA period: {period_display_name}")
        print(f"Internal label: {period_internal_name}")
        print(f"ROOT file: {root_path}")
        print("=" * 90)

        run_rows, current_rows = build_period_aggregation(
            period_display_name=period_display_name,
            period_internal_name=period_internal_name,
            root_path=root_path,
            run_charge_map=run_charge_map,
        )

        all_run_rows.extend(run_rows)
        all_current_rows.extend(current_rows)

        print(f"Run-level rows: {len(run_rows)}")
        print(f"Current groups: {len(current_rows)}")

        for row in current_rows:
            print(
                f"  {period_display_name:8s}  "
                f"{int(row['current_nA']):3d} nA  "
                f"runs={int(row['n_runs']):3d}  "
                f"counts={int(row['counts']):8d}  "
                f"charge={float(row['charge_nC']):12.5f} nC  "
                f"counts/nC={float(row['counts_per_nC']):.8f} +/- {float(row['counts_per_nC_err']):.8f}"
            )
        #endfor
    #endfor

    # -------------------------------------------------------------------------
    # Build MC table.
    # -------------------------------------------------------------------------
    print("")
    print("Processing MC ROOT files and computing efficiencies...")
    if args.skip_temp_heavy_mc:
        print("Temporary MC skip override flag is ON.")
    #endif

    mc_rows = build_mc_aggregation(MC_DIR, skip_temp_heavy_mc=args.skip_temp_heavy_mc)

    for row in mc_rows:
        print(
            f"  {row['period']:8s}  "
            f"{int(row['current_nA']):3d} nA  "
            f"N_gen={int(row['n_gen']):8d}  "
            f"N_rec={int(row['n_rec']):8d}  "
            f"eff={float(row['efficiency']):.8f} +/- {float(row['efficiency_err']):.8f}"
        )
    #endfor

    # -------------------------------------------------------------------------
    # Write diagnostic CSVs.
    # -------------------------------------------------------------------------
    run_table_csv = os.path.join(OUTPUT_DIR, "dvcs_current_dependence_run_table.csv")
    current_table_csv = os.path.join(OUTPUT_DIR, "dvcs_current_dependence_current_table.csv")
    mc_table_csv = os.path.join(OUTPUT_DIR, "dvcs_current_dependence_mc_table.csv")

    write_run_table_csv(run_table_csv, all_run_rows)
    write_current_table_csv(current_table_csv, all_current_rows)
    write_mc_table_csv(mc_table_csv, mc_rows)

    print("")
    print(f"[saved] {run_table_csv}")
    print(f"[saved] {current_table_csv}")
    print(f"[saved] {mc_table_csv}")

    # -------------------------------------------------------------------------
    # Choose a consistent color per period.
    # -------------------------------------------------------------------------
    default_colors = plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
    if len(default_colors) < len(PERIOD_ORDER):
        raise RuntimeError("Matplotlib default color cycle is shorter than number of periods.")
    #endif

    period_color = {}
    for i, period in enumerate(PERIOD_ORDER):
        period_color[period] = default_colors[i]
    #endfor

    # -------------------------------------------------------------------------
    # Fit DATA.
    # -------------------------------------------------------------------------
    data_fit_results = {}

    print("")
    print("=== DATA fits: y = m*x + b for counts_per_nC vs current ===")

    for period in PERIOD_ORDER:
        x, y, sy, rows = period_points_from_current_rows(all_current_rows, period)
        fr = weighted_linear_fit(x, y, sy)
        data_fit_results[period] = fr

        print(
            f"{period}: m = {fr['m']:.10f} +/- {fr['sm']:.10f}, "
            f"b = {fr['b']:.10f} +/- {fr['sb']:.10f}, "
            f"chi2/ndof = {fr['chi2']:.2f}/{fr['ndof']}"
        )
    #endfor

    # -------------------------------------------------------------------------
    # Fit MC.
    # -------------------------------------------------------------------------
    mc_fit_results = {}

    print("")
    print("=== MC fits: y = m*x + b for efficiency vs current ===")

    for period in PERIOD_ORDER:
        x, y, sy, rows = period_points_from_mc_rows(mc_rows, period)
        fr = weighted_linear_fit(x, y, sy)
        mc_fit_results[period] = fr

        print(
            f"{period}: m = {fr['m']:.10f} +/- {fr['sm']:.10f}, "
            f"b = {fr['b']:.10f} +/- {fr['sb']:.10f}, "
            f"chi2/ndof = {fr['chi2']:.2f}/{fr['ndof']}"
        )
    #endfor

    # -------------------------------------------------------------------------
    # Build and write the operational residual correction table.
    # -------------------------------------------------------------------------
    residual_rows = build_residual_correction_rows(
        all_current_rows=all_current_rows,
        data_fit_results=data_fit_results,
        mc_fit_results=mc_fit_results,
    )

    residual_table_csv = os.path.join(OUTPUT_DIR, "dvcs_current_dependence_residual_table.csv")
    write_residual_table_csv(residual_table_csv, residual_rows)

    print("")
    print("[saved] " + residual_table_csv)
    print("")
    print("=== Operational residual current-dependence correction ===")
    print("Stored quantity for efficiency.json / cross_sections.cpp usage:")
    print("  epsilon_resid(I) = epsilon_data_rel(I) / epsilon_mc_rel(I_ref)")
    print("  applied scale    = 1 / epsilon_resid(I)")
    print("where I_ref is the single MC current actually used in production acceptance.")
    print("")

    for row in residual_rows:
        print(
            f"  {row['period']:8s}  "
            f"data I={int(row['data_current_nA']):3d} nA  "
            f"MC ref={int(row['mc_reference_current_nA']):3d} nA  "
            f"data_rel={float(row['data_rel_fraction']):.6f}  "
            f"mc_rel_ref={float(row['mc_rel_at_reference_fraction']):.6f}  "
            f"epsilon_resid={float(row['residual_fraction']):.6f}  "
            f"scale={float(row['applied_scale_fraction']):.6f}"
        )
    #endfor

    # -------------------------------------------------------------------------
    # Build compact plotting cache.
    # -------------------------------------------------------------------------
    cache = prepare_period_cache(
        all_current_rows=all_current_rows,
        mc_rows=mc_rows,
        data_fit_results=data_fit_results,
        mc_fit_results=mc_fit_results,
    )

    # -------------------------------------------------------------------------
    # Streamlined integrated plots requested by the user.
    # -------------------------------------------------------------------------

    out1 = os.path.join(INTEGRATED_DIR, "dvcs_data_counts_per_nc_integrated.png")
    plot_six_panel_simple(
        period_order=PERIOD_ORDER,
        cache=cache,
        period_color=period_color,
        drawer=draw_data_absolute_panel,
        outfile=out1,
    )
    print("")
    print(f"[saved] {out1}")

    out2 = os.path.join(INTEGRATED_DIR, "dvcs_data_percent_of_zero_integrated.png")
    plot_six_panel_simple(
        period_order=PERIOD_ORDER,
        cache=cache,
        period_color=period_color,
        drawer=draw_data_percent_panel,
        outfile=out2,
    )
    print("")
    print(f"[saved] {out2}")

    out3 = os.path.join(INTEGRATED_DIR, "dvcs_absolute_data_mc_integrated.png")
    plot_six_panel_simple(
        period_order=PERIOD_ORDER,
        cache=cache,
        period_color=period_color,
        drawer=draw_absolute_overlay_panel,
        outfile=out3,
    )
    print("")
    print(f"[saved] {out3}")

    out4 = os.path.join(INTEGRATED_DIR, "dvcs_percent_data_mc_with_ratio_integrated.png")
    plot_percent_overlay_with_ratio(
        period_order=PERIOD_ORDER,
        cache=cache,
        period_color=period_color,
        outfile=out4,
    )
    print("")
    print(f"[saved] {out4}")
    print("")
#enddef


if __name__ == "__main__":
    main()
#endif