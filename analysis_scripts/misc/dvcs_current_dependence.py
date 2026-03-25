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
# Integrated plots produced in:
#   output/dvcs_current_dependence/integrated
#
#   (1) six-panel data counts per accumulated charge vs current
#   (2) six-panel data drop from assumed 100 percent efficiency at 0 current
#   (3) six-panel data counts per accumulated charge + MC rec/gen together
#   (4) six-panel data + MC drop from assumed 100 percent efficiency at 0 current,
#       with a ratio pad below each period subplot and in the overlay panel
#
# Diagnostic CSVs:
#   output/dvcs_current_dependence/dvcs_current_dependence_run_table.csv
#   output/dvcs_current_dependence/dvcs_current_dependence_current_table.csv
#   output/dvcs_current_dependence/dvcs_current_dependence_mc_table.csv
#   output/dvcs_current_dependence/dvcs_current_dependence_residual_table.csv
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
INTEGRATED_OUTPUT_DIR = os.path.join(OUTPUT_DIR, "integrated")


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
    try:
        return f"{float(val):.5f}"
    except Exception:
        return "nan"
    #endif
#enddef


def weighted_linear_fit(x, y, sy):

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

    ref_current = get_mc_reference_current(period_name)
    text = f"MC ref in acceptance: {ref_current} nA"

    ax.text(
        0.03,
        0.03,
        text,
        transform=ax.transAxes,
        fontsize=9,
        verticalalignment="bottom",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.75, edgecolor="0.7"),
    )
#enddef


def build_period_aggregation(period_display_name, period_internal_name, root_path, run_charge_map):

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

    if not os.path.isdir(mc_dir):
        raise RuntimeError(f"MC directory not found: {mc_dir}")
    #endif

    skip_pairs = set()
    if skip_temp_heavy_mc:
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

    m = fit_result["m"]
    b = fit_result["b"]
    sm = fit_result["sm"]

    pct_fit = 100.0 * ((m * xfit + b) / b)
    pct_fit_lo = 100.0 * (((m - sm) * xfit + b) / b)
    pct_fit_hi = 100.0 * (((m + sm) * xfit + b) / b)

    return pct_fit, pct_fit_lo, pct_fit_hi
#enddef


def compute_percent_slope_and_error(fit_result):

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

    m = fit_result["m"]
    b = fit_result["b"]

    return (m * float(current_nA) + b) / b
#enddef


def build_residual_correction_rows(all_current_rows, data_fit_results, mc_fit_results):

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


def style_percent_axis(ax, ylabel):
    ax.set_xlim(0.0, 80.0)
    ax.set_ylim(40.0, 120.0)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)
#enddef


def style_absolute_axis(ax, ylabel):
    ax.set_xlim(0.0, 80.0)
    ax.set_xlabel("Beam current (nA)")
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)
#enddef


def style_ratio_axis(ax):
    ax.set_xlim(0.0, 80.0)
    ax.set_ylim(0.80, 1.05)
    ax.set_xlabel("Beam current (nA)")
    ax.set_ylabel("Data/MC")
    ax.grid(True, alpha=0.3)
    ax.axhline(1.0, color="0.5", linestyle="--", linewidth=1.0)

    ticks = [0.80, 0.85, 0.90, 0.95, 1.00, 1.05]
    ax.set_yticks(ticks)
#enddef


def create_simple_2x3_figure():
    fig = plt.figure(figsize=(18, 10), constrained_layout=True)
    gs = GridSpec(2, 3, figure=fig)
    axes = [fig.add_subplot(gs[i // 3, i % 3]) for i in range(6)]
    return fig, axes
#enddef


def create_doublepad_2x3_figure():
    fig = plt.figure(figsize=(18, 11), constrained_layout=True)
    outer = GridSpec(2, 3, figure=fig)

    top_axes = []
    bottom_axes = []

    for i in range(6):
        r = i // 3
        c = i % 3
        inner = GridSpecFromSubplotSpec(
            2, 1,
            subplot_spec=outer[r, c],
            height_ratios=[3.2, 1.2],
            hspace=0.0,
        )
        ax_top = fig.add_subplot(inner[0])
        ax_bot = fig.add_subplot(inner[1], sharex=ax_top)
        plt.setp(ax_top.get_xticklabels(), visible=False)
        top_axes.append(ax_top)
        bottom_axes.append(ax_bot)
    #endfor

    return fig, top_axes, bottom_axes
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
    os.makedirs(INTEGRATED_OUTPUT_DIR, exist_ok=True)

    print("")
    print("Reading charge map...")
    run_charge_map = read_charge_map(CSV_FILE)
    print(f"Loaded charge entries for {len(run_charge_map)} runs from:")
    print(f"  {CSV_FILE}")

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

    default_colors = plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
    if len(default_colors) < len(PERIOD_ORDER):
        raise RuntimeError("Matplotlib default color cycle is shorter than number of periods.")
    #endif

    period_color = {}
    for i, period in enumerate(PERIOD_ORDER):
        period_color[period] = default_colors[i]
    #endfor

    band_alpha = 0.20
    xfit = np.linspace(0.0, 80.0, 300)

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

    residual_rows = build_residual_correction_rows(
        all_current_rows=all_current_rows,
        data_fit_results=data_fit_results,
        mc_fit_results=mc_fit_results,
    )

    residual_table_csv = os.path.join(OUTPUT_DIR, "dvcs_current_dependence_residual_table.csv")
    write_residual_table_csv(residual_table_csv, residual_rows)

    print("")
    print("[saved] " + residual_table_csv)

    # -------------------------------------------------------------------------
    # Plot A: data counts per accumulated charge vs current
    # -------------------------------------------------------------------------
    fig_a, axes_a = create_simple_2x3_figure()

    for i, period in enumerate(PERIOD_ORDER):
        ax = axes_a[i]
        c = period_color[period]

        x, y, sy, rows = period_points_from_current_rows(all_current_rows, period)
        fr = data_fit_results[period]

        yfit = fr["m"] * xfit + fr["b"]
        yfit_lo = (fr["m"] - fr["sm"]) * xfit + fr["b"]
        yfit_hi = (fr["m"] + fr["sm"]) * xfit + fr["b"]

        ax.fill_between(xfit, yfit_lo, yfit_hi, color=c, alpha=band_alpha, linewidth=0)
        ax.errorbar(x, y, yerr=sy, fmt="o", capsize=3, color=c, label="Data")
        ax.plot(xfit, yfit, color=c)

        ax.set_title(period)
        style_absolute_axis(ax, "Counts / accumulated charge (1/nC)")
        ax.legend(frameon=True)
        add_reference_current_text(ax, period)
    #endfor

    ax = axes_a[5]
    ax.set_title("All periods (overlay)")
    style_absolute_axis(ax, "Counts / accumulated charge (1/nC)")

    for period in PERIOD_ORDER:
        c = period_color[period]
        x, y, sy, rows = period_points_from_current_rows(all_current_rows, period)
        fr = data_fit_results[period]

        yfit = fr["m"] * xfit + fr["b"]
        yfit_lo = (fr["m"] - fr["sm"]) * xfit + fr["b"]
        yfit_hi = (fr["m"] + fr["sm"]) * xfit + fr["b"]

        ax.fill_between(xfit, yfit_lo, yfit_hi, color=c, alpha=0.08, linewidth=0)
        ax.errorbar(x, y, yerr=sy, fmt="o", capsize=3, color=c, label=period)
        ax.plot(xfit, yfit, color=c)
    #endfor

    ax.legend(frameon=True, fontsize=9)

    out_a = os.path.join(INTEGRATED_OUTPUT_DIR, "dvcs_counts_per_charge_data_integrated.png")
    fig_a.savefig(out_a, dpi=200)
    print("")
    print(f"[saved] {out_a}")

    # -------------------------------------------------------------------------
    # Plot B: data drop from assumed 100 percent efficiency at 0 current
    # -------------------------------------------------------------------------
    fig_b, axes_b = create_simple_2x3_figure()

    for i, period in enumerate(PERIOD_ORDER):
        ax = axes_b[i]
        c = period_color[period]

        x, y, sy, rows = period_points_from_current_rows(all_current_rows, period)
        fr = data_fit_results[period]

        pct = 100.0 * (y / fr["b"])
        pct_err = 100.0 * np.sqrt((sy / fr["b"]) ** 2 + ((y * fr["sb"]) / (fr["b"] * fr["b"])) ** 2)
        pct_fit, pct_fit_lo, pct_fit_hi = compute_percent_curve(xfit, fr)

        ax.fill_between(xfit, pct_fit_lo, pct_fit_hi, color=c, alpha=band_alpha, linewidth=0)
        ax.errorbar(x, pct, yerr=pct_err, fmt="o", capsize=3, color=c, label="Data")
        ax.plot(xfit, pct_fit, color=c)

        ax.set_title(period)
        style_percent_axis(ax, "Efficiency relative to fitted 0 nA (%)")
        ax.set_xlabel("Beam current (nA)")
        ax.legend(frameon=True)
        add_reference_current_text(ax, period)
    #endfor

    ax = axes_b[5]
    ax.set_title("All periods (overlay)")
    style_percent_axis(ax, "Efficiency relative to fitted 0 nA (%)")
    ax.set_xlabel("Beam current (nA)")

    for period in PERIOD_ORDER:
        c = period_color[period]
        x, y, sy, rows = period_points_from_current_rows(all_current_rows, period)
        fr = data_fit_results[period]

        pct = 100.0 * (y / fr["b"])
        pct_err = 100.0 * np.sqrt((sy / fr["b"]) ** 2 + ((y * fr["sb"]) / (fr["b"] * fr["b"])) ** 2)
        pct_fit, pct_fit_lo, pct_fit_hi = compute_percent_curve(xfit, fr)

        ax.fill_between(xfit, pct_fit_lo, pct_fit_hi, color=c, alpha=0.08, linewidth=0)
        ax.errorbar(x, pct, yerr=pct_err, fmt="o", capsize=3, color=c, label=period)
        ax.plot(xfit, pct_fit, color=c)
    #endfor

    ax.legend(frameon=True, fontsize=9)

    out_b = os.path.join(INTEGRATED_OUTPUT_DIR, "dvcs_percent_of_zero_data_integrated.png")
    fig_b.savefig(out_b, dpi=200)
    print(f"[saved] {out_b}")

    # -------------------------------------------------------------------------
    # Plot C: data counts per accumulated charge and MC rec/gen together
    # -------------------------------------------------------------------------
    fig_c, axes_c = create_simple_2x3_figure()

    for i, period in enumerate(PERIOD_ORDER):
        ax = axes_c[i]
        c = period_color[period]

        xd, yd, syd, rowsd = period_points_from_current_rows(all_current_rows, period)
        frd = data_fit_results[period]
        data_fit = frd["m"] * xfit + frd["b"]
        data_fit_lo = (frd["m"] - frd["sm"]) * xfit + frd["b"]
        data_fit_hi = (frd["m"] + frd["sm"]) * xfit + frd["b"]

        xm, ym, sym, rowsm = period_points_from_mc_rows(mc_rows, period)
        frm = mc_fit_results[period]
        mc_fit = frm["m"] * xfit + frm["b"]
        mc_fit_lo = (frm["m"] - frm["sm"]) * xfit + frm["b"]
        mc_fit_hi = (frm["m"] + frm["sm"]) * xfit + frm["b"]

        ax.fill_between(xfit, data_fit_lo, data_fit_hi, color=c, alpha=0.12, linewidth=0)
        ax.errorbar(xd, yd, yerr=syd, fmt="o", capsize=3, color=c, label="Data")
        ax.plot(xfit, data_fit, color=c)

        ax.fill_between(xfit, mc_fit_lo, mc_fit_hi, color=c, alpha=0.08, linewidth=0)
        ax.errorbar(
            xm, ym, yerr=sym,
            fmt="o",
            capsize=3,
            linestyle="none",
            markerfacecolor="none",
            markeredgecolor=c,
            ecolor=c,
            color=c,
            label="MC",
        )
        ax.plot(xfit, mc_fit, color=c, linestyle="--")

        ax.set_title(period)
        style_absolute_axis(ax, "Absolute value")
        ax.legend(frameon=True)
        add_reference_current_text(ax, period)
    #endfor

    ax = axes_c[5]
    ax.set_title("All periods (overlay)")
    style_absolute_axis(ax, "Absolute value")

    for period in PERIOD_ORDER:
        c = period_color[period]

        xd, yd, syd, rowsd = period_points_from_current_rows(all_current_rows, period)
        frd = data_fit_results[period]
        data_fit = frd["m"] * xfit + frd["b"]
        data_fit_lo = (frd["m"] - frd["sm"]) * xfit + frd["b"]
        data_fit_hi = (frd["m"] + frd["sm"]) * xfit + frd["b"]

        xm, ym, sym, rowsm = period_points_from_mc_rows(mc_rows, period)
        frm = mc_fit_results[period]
        mc_fit = frm["m"] * xfit + frm["b"]
        mc_fit_lo = (frm["m"] - frm["sm"]) * xfit + frm["b"]
        mc_fit_hi = (frm["m"] + frm["sm"]) * xfit + frm["b"]

        ax.fill_between(xfit, data_fit_lo, data_fit_hi, color=c, alpha=0.08, linewidth=0)
        ax.errorbar(xd, yd, yerr=syd, fmt="o", capsize=3, color=c, label=f"{period} data")
        ax.plot(xfit, data_fit, color=c)

        ax.fill_between(xfit, mc_fit_lo, mc_fit_hi, color=c, alpha=0.05, linewidth=0)
        ax.errorbar(
            xm, ym, yerr=sym,
            fmt="o",
            capsize=3,
            linestyle="none",
            markerfacecolor="none",
            markeredgecolor=c,
            ecolor=c,
            color=c,
            label=f"{period} MC",
        )
        ax.plot(xfit, mc_fit, color=c, linestyle="--")
    #endfor

    ax.legend(frameon=True, fontsize=8)

    out_c = os.path.join(INTEGRATED_OUTPUT_DIR, "dvcs_absolute_data_vs_mc_integrated.png")
    fig_c.savefig(out_c, dpi=200)
    print(f"[saved] {out_c}")

    # -------------------------------------------------------------------------
    # Plot D: data and MC drop from assumed 100 percent efficiency at 0 current
    #         with ratio pad below each subplot
    # -------------------------------------------------------------------------
    fig_d, top_axes_d, bottom_axes_d = create_doublepad_2x3_figure()

    for i, period in enumerate(PERIOD_ORDER):
        ax_top = top_axes_d[i]
        ax_bot = bottom_axes_d[i]
        c = period_color[period]

        xd, yd, syd, rowsd = period_points_from_current_rows(all_current_rows, period)
        frd = data_fit_results[period]
        data_pct = 100.0 * (yd / frd["b"])
        data_pct_err = 100.0 * np.sqrt((syd / frd["b"]) ** 2 + ((yd * frd["sb"]) / (frd["b"] * frd["b"])) ** 2)
        data_pct_fit, data_pct_fit_lo, data_pct_fit_hi = compute_percent_curve(xfit, frd)

        xm, ym, sym, rowsm = period_points_from_mc_rows(mc_rows, period)
        frm = mc_fit_results[period]
        mc_pct = 100.0 * (ym / frm["b"])
        mc_pct_err = 100.0 * np.sqrt((sym / frm["b"]) ** 2 + ((ym * frm["sb"]) / (frm["b"] * frm["b"])) ** 2)
        mc_pct_fit, mc_pct_fit_lo, mc_pct_fit_hi = compute_percent_curve(xfit, frm)

        ax_top.fill_between(xfit, data_pct_fit_lo, data_pct_fit_hi, color=c, alpha=0.12, linewidth=0)
        ax_top.errorbar(xd, data_pct, yerr=data_pct_err, fmt="o", capsize=3, color=c, label="Data")
        ax_top.plot(xfit, data_pct_fit, color=c)

        ax_top.fill_between(xfit, mc_pct_fit_lo, mc_pct_fit_hi, color=c, alpha=0.08, linewidth=0)
        ax_top.errorbar(
            xm, mc_pct, yerr=mc_pct_err,
            fmt="o",
            capsize=3,
            linestyle="none",
            markerfacecolor="none",
            markeredgecolor=c,
            ecolor=c,
            color=c,
            label="MC",
        )
        ax_top.plot(xfit, mc_pct_fit, color=c, linestyle="--")

        ax_top.set_title(period)
        style_percent_axis(ax_top, "Efficiency relative to fitted 0 nA (%)")
        ax_top.set_xlabel("")
        ax_top.legend(frameon=True)
        add_reference_current_text(ax_top, period)

        ratio_x = np.intersect1d(xd, xm)
        ratio_y = []
        ratio_sy = []

        for xr in ratio_x:
            idd = np.where(xd == xr)[0][0]
            idm = np.where(xm == xr)[0][0]

            rv = data_pct[idd] / mc_pct[idm]

            rel_err_sq = 0.0
            if data_pct[idd] != 0.0:
                rel_err_sq += (data_pct_err[idd] / data_pct[idd]) ** 2
            #endif
            if mc_pct[idm] != 0.0:
                rel_err_sq += (mc_pct_err[idm] / mc_pct[idm]) ** 2
            #endif

            re = rv * math.sqrt(rel_err_sq)

            ratio_y.append(rv)
            ratio_sy.append(re)
        #endfor

        ratio_x = np.asarray(ratio_x, dtype=float)
        ratio_y = np.asarray(ratio_y, dtype=float)
        ratio_sy = np.asarray(ratio_sy, dtype=float)

        data_fit_curve = 100.0 * ((frd["m"] * xfit + frd["b"]) / frd["b"])
        mc_fit_curve = 100.0 * ((frm["m"] * xfit + frm["b"]) / frm["b"])
        ratio_fit_curve = data_fit_curve / mc_fit_curve

        data_fit_lo_curve = 100.0 * (((frd["m"] - frd["sm"]) * xfit + frd["b"]) / frd["b"])
        data_fit_hi_curve = 100.0 * (((frd["m"] + frd["sm"]) * xfit + frd["b"]) / frd["b"])
        mc_fit_lo_curve = 100.0 * (((frm["m"] - frm["sm"]) * xfit + frm["b"]) / frm["b"])
        mc_fit_hi_curve = 100.0 * (((frm["m"] + frm["sm"]) * xfit + frm["b"]) / frm["b"])

        ratio_fit_lo_curve = data_fit_lo_curve / mc_fit_hi_curve
        ratio_fit_hi_curve = data_fit_hi_curve / mc_fit_lo_curve

        ax_bot.fill_between(xfit, ratio_fit_lo_curve, ratio_fit_hi_curve, color=c, alpha=0.12, linewidth=0)
        ax_bot.errorbar(ratio_x, ratio_y, yerr=ratio_sy, fmt="o", capsize=3, color=c)
        ax_bot.plot(xfit, ratio_fit_curve, color=c)

        style_ratio_axis(ax_bot)
    #endfor

    ax_top = top_axes_d[5]
    ax_bot = bottom_axes_d[5]
    ax_top.set_title("All periods (overlay)")
    style_percent_axis(ax_top, "Efficiency relative to fitted 0 nA (%)")
    ax_top.set_xlabel("")
    style_ratio_axis(ax_bot)

    for period in PERIOD_ORDER:
        c = period_color[period]

        xd, yd, syd, rowsd = period_points_from_current_rows(all_current_rows, period)
        frd = data_fit_results[period]
        data_pct = 100.0 * (yd / frd["b"])
        data_pct_err = 100.0 * np.sqrt((syd / frd["b"]) ** 2 + ((yd * frd["sb"]) / (frd["b"] * frd["b"])) ** 2)
        data_pct_fit, data_pct_fit_lo, data_pct_fit_hi = compute_percent_curve(xfit, frd)

        xm, ym, sym, rowsm = period_points_from_mc_rows(mc_rows, period)
        frm = mc_fit_results[period]
        mc_pct = 100.0 * (ym / frm["b"])
        mc_pct_err = 100.0 * np.sqrt((sym / frm["b"]) ** 2 + ((ym * frm["sb"]) / (frm["b"] * frm["b"])) ** 2)
        mc_pct_fit, mc_pct_fit_lo, mc_pct_fit_hi = compute_percent_curve(xfit, frm)

        ax_top.fill_between(xfit, data_pct_fit_lo, data_pct_fit_hi, color=c, alpha=0.08, linewidth=0)
        ax_top.errorbar(xd, data_pct, yerr=data_pct_err, fmt="o", capsize=3, color=c, label=f"{period} data")
        ax_top.plot(xfit, data_pct_fit, color=c)

        ax_top.fill_between(xfit, mc_pct_fit_lo, mc_pct_fit_hi, color=c, alpha=0.05, linewidth=0)
        ax_top.errorbar(
            xm, mc_pct, yerr=mc_pct_err,
            fmt="o",
            capsize=3,
            linestyle="none",
            markerfacecolor="none",
            markeredgecolor=c,
            ecolor=c,
            color=c,
            label=f"{period} MC",
        )
        ax_top.plot(xfit, mc_pct_fit, color=c, linestyle="--")

        ratio_x = np.intersect1d(xd, xm)
        ratio_y = []
        ratio_sy = []

        for xr in ratio_x:
            idd = np.where(xd == xr)[0][0]
            idm = np.where(xm == xr)[0][0]

            rv = data_pct[idd] / mc_pct[idm]

            rel_err_sq = 0.0
            if data_pct[idd] != 0.0:
                rel_err_sq += (data_pct_err[idd] / data_pct[idd]) ** 2
            #endif
            if mc_pct[idm] != 0.0:
                rel_err_sq += (mc_pct_err[idm] / mc_pct[idm]) ** 2
            #endif

            re = rv * math.sqrt(rel_err_sq)

            ratio_y.append(rv)
            ratio_sy.append(re)
        #endfor

        ratio_x = np.asarray(ratio_x, dtype=float)
        ratio_y = np.asarray(ratio_y, dtype=float)
        ratio_sy = np.asarray(ratio_sy, dtype=float)

        data_fit_curve = 100.0 * ((frd["m"] * xfit + frd["b"]) / frd["b"])
        mc_fit_curve = 100.0 * ((frm["m"] * xfit + frm["b"]) / frm["b"])
        ratio_fit_curve = data_fit_curve / mc_fit_curve

        ax_bot.errorbar(ratio_x, ratio_y, yerr=ratio_sy, fmt="o", capsize=3, color=c)
        ax_bot.plot(xfit, ratio_fit_curve, color=c)
    #endfor

    ax_top.legend(frameon=True, fontsize=8)

    out_d = os.path.join(INTEGRATED_OUTPUT_DIR, "dvcs_percent_of_zero_data_vs_mc_with_ratio_integrated.png")
    fig_d.savefig(out_d, dpi=200)
    print(f"[saved] {out_d}")

    print("")
#enddef


if __name__ == "__main__":
    main()
#endif