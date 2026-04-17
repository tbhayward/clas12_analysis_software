#!/usr/bin/env python3
#
# dvcs_current_dependence.py
#
# Purpose
# -------
# Study current dependence in reconstruction using:
#
#   DATA:
#     reconstructed counts / accumulated charge
#
#   MC:
#     reconstructed / generated
#
# Operational logic for the production cross section correction
# -------------------------------------------------------------
# The production acceptance correction uses only one MC current per period:
#
#   Sp18 Inb -> 50 nA
#   Sp18 Out -> 45 nA
#   Fa18 Inb -> 50 nA
#   Fa18 Out -> 50 nA
#   Sp19 Inb -> 50 nA
#
# Therefore the period-dependent divisor relevant for cross_sections.cpp is:
#
#   divisor(I) = epsilon_data_rel(I) / epsilon_mc_rel(I_ref)
#
# where
#
#   epsilon_data_rel(I)   = R_data(I) / R_data(0)
#   epsilon_mc_rel(I_ref) = epsilon_mc(I_ref) / epsilon_mc(0)
#
# and cross_sections.cpp should divide by divisor(I), equivalently multiply by:
#
#   applied_scale(I) = 1.0 / divisor(I)
#
# In this script we also build a single representative period-level number:
#
#   weighted_data_rel
#     = event-weighted average of epsilon_data_rel(I)
#       using reconstructed DATA counts as weights across the actual current mix
#
#   mc_ref_rel
#     = epsilon_mc_rel(I_ref)
#
#   divisor_weighted
#     = weighted_data_rel / mc_ref_rel
#
# This is the representative value to divide by if you want one integrated
# period-level divisor compatible with cross_sections.cpp.
#
# Detector angle dependence
# -------------------------
# We also repeat the study in detector-level angular bins using:
#
#   e_theta   (electron) : 8 to 35 deg
#   p1_theta  (proton)   : 8 to 70 deg
#   p2_theta  (photon)   : 0 to 35 deg
#
# These branches are stored in radians in the ROOT trees and are converted
# to degrees before binning.
#
# For DATA in an angular bin:
#   counts are restricted to that angular bin
#   charge remains the full accumulated charge for the runs/current settings
#
# For MC in an angular bin:
#   just use generated and reconstructed counts in that same angular bin
#
# Fine histogram plots
# --------------------
# For each detector angle variable we also make a 1x2 overlay figure:
#
#   left  : fine-binned counts / accumulated charge for DATA
#   right : the same DATA histograms normalized to their own integral
#
# We also make an MC reconstructed version:
#
#   left  : fine-binned reconstructed MC counts
#   right : the same reconstructed MC histograms normalized to their own integral
#
# All run periods are drawn on the same axes in each panel.
#
# Additional event cut
# --------------------
# Any event is skipped if the cone angle between the scattered electron and
# the photon is less than 7 degrees.
#
# Additional photon histogram version
# ----------------------------------
# For the photon theta 1x2 distribution plot only, an additional version is
# produced requiring p2_phi > 60 degrees.
#
# Additional photon sector panel
# ------------------------------
# For photon theta we also produce a 3x2 panel of normalized distributions,
# one pad per CLAS12 sector, where sectors are defined from p2_phi:
#
#   sector 1 : [330,360) U [0,30)
#   sector 2 : [30,90)
#   sector 3 : [90,150)
#   sector 4 : [150,210)
#   sector 5 : [210,270)
#   sector 6 : [270,330)
#
# For the sector panel, nothing is drawn below 5 degrees.

import os
import math
import csv
from collections import defaultdict
import concurrent.futures as cf

import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec


# -----------------------------------------------------------------------------
# Input configuration
# -----------------------------------------------------------------------------

CSV_FILE = "/u/home/thayward/clas12_analysis_software/analysis_scripts/dvcs_cross_section/imports/integrated_luminosity/global.csv"

PERIOD_ORDER = ["Sp18 Inb", "Sp18 Out", "Fa18 Inb", "Fa18 Out", "Sp19 Inb"]
MAX_WORKERS = 5
ITERATE_STEP_SIZE = "200 MB"
MC_TREE_NAME = "PhysicsEvents"
MIN_E_GAMMA_CONE_ANGLE_DEG = 7.0
PHOTON_HIST_PHI_MIN_DEG = 60.0

COS_MIN_E_GAMMA_CONE = math.cos(math.radians(MIN_E_GAMMA_CONE_ANGLE_DEG))

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

DATA_PERIOD_FILES = [
    ("Sp18 Inb", "rga_sp18_inb", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp18_inb_epgamma.root"),
    ("Sp18 Out", "rga_sp18_out", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp18_out_epgamma.root"),
    ("Fa18 Inb", "rga_fa18_inb", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_inb_epgamma.root"),
    ("Fa18 Out", "rga_fa18_out", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_out_epgamma.root"),
    ("Sp19 Inb", "rga_sp19_inb", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp19_inb_epgamma.root"),
]

MC_DIR = "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen"

ANGLE_DEPENDENCE_CONFIG = [
    {
        "key": "e_theta",
        "display_name": "electron theta",
        "branch": "e_theta",
        "min_deg": 8.0,
        "max_deg": 35.0,
        "n_bins": 7,
        "plot_n_bins": 120,
    },
    {
        "key": "p1_theta",
        "display_name": "proton theta",
        "branch": "p1_theta",
        "min_deg": 8.0,
        "max_deg": 70.0,
        "n_bins": 7,
        "plot_n_bins": 140,
    },
    {
        "key": "p2_theta",
        "display_name": "photon theta",
        "branch": "p2_theta",
        "min_deg": 0.0,
        "max_deg": 35.0,
        "n_bins": 7,
        "plot_n_bins": 120,
    },
]

CLAS12_SECTOR_ORDER = [1, 2, 3, 4, 5, 6]


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
# Uniform angle-bin definitions
# -----------------------------------------------------------------------------

def build_angle_bins(min_deg, max_deg, n_bins):

    edges = np.linspace(min_deg, max_deg, n_bins + 1)
    bins = []

    for i in range(n_bins):
        bins.append((float(edges[i]), float(edges[i + 1])))
    #endfor

    return bins
#enddef


def get_summary_bins(var_key):

    if var_key == "e_theta":
        base = build_angle_bins(8.0, 35.0, 7)
        return [
            base[0],
            base[1],
            base[2],
            base[3],
            base[4],
            (base[5][0], base[6][1]),
        ]
    #endif

    if var_key == "p1_theta":
        base = build_angle_bins(8.0, 70.0, 7)
        return [
            (base[0][0], base[1][1]),
            base[2],
            base[3],
            base[4],
            base[5],
            base[6],
        ]
    #endif

    if var_key == "p2_theta":
        return [
            (0.0, 5.0),
            (5.0, 10.0),
            (10.0, 15.0),
            (15.0, 20.0),
            (20.0, 25.0),
            (25.0, 30.0),
            (30.0, 35.0),
        ]
    #endif

    raise RuntimeError(f"Unknown angle variable for summary bins: {var_key}")
#enddef


ANGLE_BINS = {}
PLOT_ANGLE_EDGES = {}

for var_cfg in ANGLE_DEPENDENCE_CONFIG:
    ANGLE_BINS[var_cfg["key"]] = get_summary_bins(var_cfg["key"])
    PLOT_ANGLE_EDGES[var_cfg["key"]] = np.linspace(
        float(var_cfg["min_deg"]),
        float(var_cfg["max_deg"]),
        int(var_cfg["plot_n_bins"]) + 1,
    )
#endfor


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

def get_output_paths():

    output_dir = "output/dvcs_current_dependence"
    integrated_output_dir = os.path.join(output_dir, "integrated")
    angle_output_dir = os.path.join(output_dir, "angle_dependence")

    return output_dir, integrated_output_dir, angle_output_dir
#enddef


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


def weighted_linear_fit(x, y, sy):

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    sy = np.asarray(sy, dtype=float)

    if x.size == 0:
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

    if x.size == 1:
        sb = float(sy[0]) if np.isfinite(sy[0]) and sy[0] > 0.0 else 0.0
        return {
            "m": 0.0,
            "b": float(y[0]),
            "sm": 0.0,
            "sb": sb,
            "cov_mb": 0.0,
            "chi2": 0.0,
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
    beam_energy_token = tokens[6]

    if not beam_energy_token.endswith("MeV"):
        raise RuntimeError(f"Could not parse beam energy token in {basename}")
    #endif

    if current_token == "nobkg":
        current_nA = 0
    else:
        if not current_token.endswith("nA"):
            raise RuntimeError(f"Could not parse current token in {basename}")
        #endif
        current_nA = int(current_token[:-2])
    #endif

    return kind, period_internal, current_nA, beam_energy_token
#enddef


def is_candidate_mc_filename(basename):
    if not basename.endswith(".root"):
        return False
    #endif
    if basename.startswith("gen_") or basename.startswith("rec_"):
        return True
    #endif
    return False
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


def angle_bin_to_dirname(variable_key, low, high):
    return f"{variable_key}_{low:.2f}_{high:.2f}".replace(".", "p")
#enddef


def angle_bin_to_label(low, high, is_last=False):
    if is_last:
        return f"[{low:.2f}, {high:.2f}]"
    #endif
    return f"[{low:.2f}, {high:.2f})"
#enddef


def angle_bin_masks(theta_deg_array, bins):

    masks = []
    n_bins = len(bins)

    for i, (low, high) in enumerate(bins):
        if i == n_bins - 1:
            mask = (theta_deg_array >= low) & (theta_deg_array <= high)
        else:
            mask = (theta_deg_array >= low) & (theta_deg_array < high)
        #endif
        masks.append(mask)
    #endfor

    return masks
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


def wrap_phi_deg(phi_deg):

    phi_wrapped = np.mod(phi_deg, 360.0)
    phi_wrapped = np.where(phi_wrapped < 0.0, phi_wrapped + 360.0, phi_wrapped)
    return phi_wrapped
#enddef


def get_clas12_sector_from_phi_deg(phi_deg_array):

    phi = wrap_phi_deg(phi_deg_array)
    sector = np.zeros(phi.shape, dtype=np.int32)

    sector[((phi >= 330.0) & (phi < 360.0)) | ((phi >= 0.0) & (phi < 30.0))] = 1
    sector[(phi >= 30.0) & (phi < 90.0)] = 2
    sector[(phi >= 90.0) & (phi < 150.0)] = 3
    sector[(phi >= 150.0) & (phi < 210.0)] = 4
    sector[(phi >= 210.0) & (phi < 270.0)] = 5
    sector[(phi >= 270.0) & (phi < 330.0)] = 6

    return sector
#enddef


def compute_e_gamma_dot(arrays):

    required = ["e_theta", "e_phi", "p2_theta", "p2_phi"]
    missing = [name for name in required if name not in arrays]
    if len(missing) > 0:
        raise RuntimeError(f"Missing branches required for cone-angle cut: {missing}")
    #endif

    e_theta = arrays["e_theta"]
    e_phi = arrays["e_phi"]
    g_theta = arrays["p2_theta"]
    g_phi = arrays["p2_phi"]

    dot = (
        np.sin(e_theta) * np.sin(g_theta) * np.cos(e_phi - g_phi) +
        np.cos(e_theta) * np.cos(g_theta)
    )

    dot = np.clip(dot, -1.0, 1.0)
    return dot
#enddef


def apply_event_mask(arrays):

    dot = compute_e_gamma_dot(arrays)
    mask = dot <= COS_MIN_E_GAMMA_CONE
    return mask
#enddef


def get_angle_axis_label(var_key):

    if var_key == "e_theta":
        return r"$\theta_{e}$ (deg)"
    #endif
    if var_key == "p1_theta":
        return r"$\theta_{p}$ (deg)"
    #endif
    if var_key == "p2_theta":
        return r"$\theta_{\gamma}$ (deg)"
    #endif

    return f"{var_key} (deg)"
#enddef


def classify_run(period_display_name, period_internal_name, runnum, run_charge_map):

    if runnum not in run_charge_map:
        return {
            "status": "missing_charge",
            "period": period_display_name,
            "period_internal": period_internal_name,
            "runnum": int(runnum),
        }
    #endif

    charge = float(run_charge_map[runnum])

    ok_cur, current_nA = resolve_current(period_internal_name, runnum)

    if charge <= 0.0:
        current_out = int(current_nA) if ok_cur else None
        return {
            "status": "nonpositive_charge",
            "period": period_display_name,
            "period_internal": period_internal_name,
            "runnum": int(runnum),
            "current_nA": current_out,
            "charge_nC": float(charge),
            "reason": "nonpositive_charge",
        }
    #endif

    if not ok_cur:
        return {
            "status": "unknown_current",
            "period": period_display_name,
            "period_internal": period_internal_name,
            "runnum": int(runnum),
        }
    #endif

    return {
        "status": "valid",
        "period": period_display_name,
        "period_internal": period_internal_name,
        "runnum": int(runnum),
        "current_nA": int(current_nA),
        "charge_nC": float(charge),
    }
#enddef


def count_mc_fine_angle_histograms(root_path, tree_name):

    if not os.path.exists(root_path):
        raise RuntimeError(f"ROOT file not found: {root_path}")
    #endif

    root_file = uproot.open(root_path)
    if tree_name not in root_file:
        raise RuntimeError(f"'{tree_name}' tree not found in {root_path}")
    #endif

    tree = root_file[tree_name]

    needed = {"e_theta", "e_phi", "p1_theta", "p2_theta", "p2_phi"}
    missing = [name for name in needed if name not in tree.keys()]
    if len(missing) > 0:
        raise RuntimeError(f"Missing required branches in {root_path}:{tree_name}: {missing}")
    #endif

    hist_counts = {}
    for var_cfg in ANGLE_DEPENDENCE_CONFIG:
        var_key = var_cfg["key"]
        hist_counts[var_key] = np.zeros(len(PLOT_ANGLE_EDGES[var_key]) - 1, dtype=np.int64)
    #endfor

    photon_hist_counts_phi_gt_60 = np.zeros(len(PLOT_ANGLE_EDGES["p2_theta"]) - 1, dtype=np.int64)
    photon_sector_hist_counts = {
        sector: np.zeros(len(PLOT_ANGLE_EDGES["p2_theta"]) - 1, dtype=np.int64)
        for sector in CLAS12_SECTOR_ORDER
    }

    iterate_branches = ["e_theta", "e_phi", "p1_theta", "p2_theta", "p2_phi"]

    for arrays in tree.iterate(iterate_branches, library="np", step_size=ITERATE_STEP_SIZE):
        event_mask = apply_event_mask(arrays)

        if not np.any(event_mask):
            continue
        #endif

        theta_deg_map = {
            "e_theta": np.degrees(arrays["e_theta"][event_mask]),
            "p1_theta": np.degrees(arrays["p1_theta"][event_mask]),
            "p2_theta": np.degrees(arrays["p2_theta"][event_mask]),
        }
        p2_phi_deg = np.degrees(arrays["p2_phi"][event_mask])

        for var_cfg in ANGLE_DEPENDENCE_CONFIG:
            var_key = var_cfg["key"]
            counts_chunk, _ = np.histogram(theta_deg_map[var_key], bins=PLOT_ANGLE_EDGES[var_key])
            hist_counts[var_key] += counts_chunk.astype(np.int64)
        #endfor

        phi_gt_60_mask = p2_phi_deg > PHOTON_HIST_PHI_MIN_DEG
        if np.any(phi_gt_60_mask):
            counts_chunk_phi_gt_60, _ = np.histogram(
                theta_deg_map["p2_theta"][phi_gt_60_mask],
                bins=PLOT_ANGLE_EDGES["p2_theta"],
            )
            photon_hist_counts_phi_gt_60 += counts_chunk_phi_gt_60.astype(np.int64)
        #endif

        photon_sector = get_clas12_sector_from_phi_deg(p2_phi_deg)
        photon_theta_deg = theta_deg_map["p2_theta"]

        for sector in CLAS12_SECTOR_ORDER:
            sector_mask = photon_sector == sector
            if not np.any(sector_mask):
                continue
            #endif

            counts_chunk_sector, _ = np.histogram(
                photon_theta_deg[sector_mask],
                bins=PLOT_ANGLE_EDGES["p2_theta"],
            )
            photon_sector_hist_counts[sector] += counts_chunk_sector.astype(np.int64)
        #endfor
    #endfor

    out = {}
    for var_cfg in ANGLE_DEPENDENCE_CONFIG:
        var_key = var_cfg["key"]
        out[var_key] = {
            "edges": PLOT_ANGLE_EDGES[var_key].copy(),
            "counts": hist_counts[var_key].copy(),
        }
    #endfor

    out["p2_theta"]["counts_phi_gt_60"] = photon_hist_counts_phi_gt_60.copy()
    out["p2_theta"]["sector_counts"] = {
        sector: photon_sector_hist_counts[sector].copy()
        for sector in CLAS12_SECTOR_ORDER
    }

    return out
#enddef


def build_data_period_aggregations(period_display_name, period_internal_name, root_path, run_charge_map):

    if not os.path.exists(root_path):
        raise RuntimeError(f"ROOT file not found: {root_path}")
    #endif

    root_file = uproot.open(root_path)
    if "PhysicsEvents" not in root_file:
        raise RuntimeError(f"'PhysicsEvents' tree not found in {root_path}")
    #endif

    tree = root_file["PhysicsEvents"]

    needed = {"runnum", "e_theta", "e_phi", "p1_theta", "p2_theta", "p2_phi"}
    missing = [name for name in needed if name not in tree.keys()]
    if len(missing) > 0:
        raise RuntimeError(f"Missing required branches in {root_path}: {missing}")
    #endif

    total_run_counts = defaultdict(int)
    angle_run_counts = {}
    angle_theta_sum = {}
    angle_theta_n = {}
    fine_hist_counts = {}

    for var_cfg in ANGLE_DEPENDENCE_CONFIG:
        var_key = var_cfg["key"]
        bins = ANGLE_BINS[var_key]
        angle_run_counts[var_key] = [defaultdict(int) for _ in bins]
        angle_theta_sum[var_key] = [0.0 for _ in bins]
        angle_theta_n[var_key] = [0 for _ in bins]
        fine_hist_counts[var_key] = np.zeros(len(PLOT_ANGLE_EDGES[var_key]) - 1, dtype=np.int64)
    #endfor

    photon_hist_counts_phi_gt_60 = np.zeros(len(PLOT_ANGLE_EDGES["p2_theta"]) - 1, dtype=np.int64)
    photon_sector_hist_counts = {
        sector: np.zeros(len(PLOT_ANGLE_EDGES["p2_theta"]) - 1, dtype=np.int64)
        for sector in CLAS12_SECTOR_ORDER
    }

    run_meta = {}
    current_charge_totals = defaultdict(float)
    skipped_nonpositive_charge_rows = []

    missing_charge_runs = set()
    unknown_current_runs = set()
    seen_runs = set()

    iterate_branches = ["runnum", "e_theta", "e_phi", "p1_theta", "p2_theta", "p2_phi"]

    for arrays in tree.iterate(iterate_branches, library="np", step_size=ITERATE_STEP_SIZE):
        base_mask = apply_event_mask(arrays)

        if not np.any(base_mask):
            continue
        #endif

        runnum_chunk = arrays["runnum"][base_mask]
        p2_phi_deg_chunk = np.degrees(arrays["p2_phi"][base_mask])

        unique_runs, inverse, counts = np.unique(runnum_chunk, return_inverse=True, return_counts=True)
        run_is_valid = np.zeros(len(unique_runs), dtype=bool)

        for i_run, runnum in enumerate(unique_runs):
            runnum_int = int(runnum)

            if runnum_int not in seen_runs:
                seen_runs.add(runnum_int)

                info = classify_run(
                    period_display_name=period_display_name,
                    period_internal_name=period_internal_name,
                    runnum=runnum_int,
                    run_charge_map=run_charge_map,
                )

                if info["status"] == "valid":
                    run_meta[runnum_int] = {
                        "period": info["period"],
                        "period_internal": info["period_internal"],
                        "runnum": info["runnum"],
                        "current_nA": info["current_nA"],
                        "charge_nC": info["charge_nC"],
                    }
                    current_charge_totals[int(info["current_nA"])] += float(info["charge_nC"])
                elif info["status"] == "nonpositive_charge":
                    skipped_nonpositive_charge_rows.append({
                        "period": info["period"],
                        "period_internal": info["period_internal"],
                        "runnum": info["runnum"],
                        "current_nA": info["current_nA"],
                        "charge_nC": info["charge_nC"],
                        "reason": info["reason"],
                    })
                elif info["status"] == "missing_charge":
                    missing_charge_runs.add(runnum_int)
                elif info["status"] == "unknown_current":
                    unknown_current_runs.add(runnum_int)
                #endif
            #endif

            if runnum_int in run_meta:
                run_is_valid[i_run] = True
                total_run_counts[runnum_int] += int(counts[i_run])
            #endif
        #endfor

        valid_event_mask = run_is_valid[inverse]

        if not np.any(valid_event_mask):
            continue
        #endif

        runnum_valid = runnum_chunk[valid_event_mask]
        p2_phi_deg_valid = p2_phi_deg_chunk[valid_event_mask]

        theta_deg_map = {
            "e_theta": np.degrees(arrays["e_theta"][base_mask][valid_event_mask]),
            "p1_theta": np.degrees(arrays["p1_theta"][base_mask][valid_event_mask]),
            "p2_theta": np.degrees(arrays["p2_theta"][base_mask][valid_event_mask]),
        }

        for var_cfg in ANGLE_DEPENDENCE_CONFIG:
            var_key = var_cfg["key"]
            bins = ANGLE_BINS[var_key]
            theta_deg = theta_deg_map[var_key]

            masks = angle_bin_masks(theta_deg, bins)
            counts_store = angle_run_counts[var_key]
            sum_store = angle_theta_sum[var_key]
            n_store = angle_theta_n[var_key]

            for ibin, mask in enumerate(masks):
                if not np.any(mask):
                    continue
                #endif

                selected_runs = runnum_valid[mask]
                unique_sel, counts_sel = np.unique(selected_runs, return_counts=True)

                for r, c in zip(unique_sel, counts_sel):
                    counts_store[ibin][int(r)] += int(c)
                #endfor

                sum_store[ibin] += float(np.sum(theta_deg[mask]))
                n_store[ibin] += int(np.count_nonzero(mask))
            #endfor

            counts_chunk_hist, _ = np.histogram(theta_deg, bins=PLOT_ANGLE_EDGES[var_key])
            fine_hist_counts[var_key] += counts_chunk_hist.astype(np.int64)
        #endfor

        phi_gt_60_mask = p2_phi_deg_valid > PHOTON_HIST_PHI_MIN_DEG
        if np.any(phi_gt_60_mask):
            counts_chunk_phi_gt_60, _ = np.histogram(
                theta_deg_map["p2_theta"][phi_gt_60_mask],
                bins=PLOT_ANGLE_EDGES["p2_theta"],
            )
            photon_hist_counts_phi_gt_60 += counts_chunk_phi_gt_60.astype(np.int64)
        #endif

        photon_sector = get_clas12_sector_from_phi_deg(p2_phi_deg_valid)
        photon_theta_deg = theta_deg_map["p2_theta"]

        for sector in CLAS12_SECTOR_ORDER:
            sector_mask = photon_sector == sector
            if not np.any(sector_mask):
                continue
            #endif

            counts_chunk_sector, _ = np.histogram(
                photon_theta_deg[sector_mask],
                bins=PLOT_ANGLE_EDGES["p2_theta"],
            )
            photon_sector_hist_counts[sector] += counts_chunk_sector.astype(np.int64)
        #endfor
    #endfor

    if len(missing_charge_runs) > 0 or len(unknown_current_runs) > 0:
        msg = []
        msg.append(f"Fatal bookkeeping problem while processing {period_display_name}:")
        if len(missing_charge_runs) > 0:
            msg.append(
                f"  Runs missing from charge CSV ({len(missing_charge_runs)}): "
                f"{sorted(missing_charge_runs)}"
            )
        #endif
        if len(unknown_current_runs) > 0:
            msg.append(
                f"  Runs with no current mapping ({len(unknown_current_runs)}): "
                f"{sorted(unknown_current_runs)}"
            )
        #endif
        raise RuntimeError("\n".join(msg))
    #endif

    integrated_run_rows = []
    for runnum in sorted(run_meta.keys()):
        count = int(total_run_counts.get(runnum, 0))
        charge = float(run_meta[runnum]["charge_nC"])
        current_nA = int(run_meta[runnum]["current_nA"])

        rate = float(count) / charge
        rate_err = math.sqrt(float(count)) / charge if count > 0 else 0.0

        integrated_run_rows.append({
            "period": run_meta[runnum]["period"],
            "period_internal": run_meta[runnum]["period_internal"],
            "runnum": int(runnum),
            "current_nA": int(current_nA),
            "counts": int(count),
            "charge_nC": float(charge),
            "counts_per_nC": float(rate),
            "counts_per_nC_err": float(rate_err),
        })
    #endfor

    integrated_current_rows = []
    per_current_counts = defaultdict(int)
    per_current_n_runs = defaultdict(int)

    for runnum in sorted(run_meta.keys()):
        current_nA = int(run_meta[runnum]["current_nA"])
        count = int(total_run_counts.get(runnum, 0))
        per_current_counts[current_nA] += int(count)
        per_current_n_runs[current_nA] += 1
    #endfor

    for current_nA in sorted(current_charge_totals.keys()):
        total_charge = float(current_charge_totals[current_nA])
        total_counts = int(per_current_counts.get(current_nA, 0))
        n_runs = int(per_current_n_runs.get(current_nA, 0))

        if total_charge <= 0.0:
            raise RuntimeError(f"Total charge <= 0 for {period_display_name}, current {current_nA} nA")
        #endif

        if total_counts <= 0:
            continue
        #endif

        counts_per_nC = float(total_counts) / total_charge
        counts_per_nC_err = math.sqrt(float(total_counts)) / total_charge

        integrated_current_rows.append({
            "period": period_display_name,
            "period_internal": period_internal_name,
            "current_nA": int(current_nA),
            "n_runs": int(n_runs),
            "counts": int(total_counts),
            "charge_nC": float(total_charge),
            "counts_per_nC": float(counts_per_nC),
            "counts_per_nC_err": float(counts_per_nC_err),
        })
    #endfor

    angle_current_rows = {}
    for var_cfg in ANGLE_DEPENDENCE_CONFIG:
        var_key = var_cfg["key"]
        bins = ANGLE_BINS[var_key]
        angle_current_rows[var_key] = {
            "bins": bins,
            "rows": [],
            "x_means": [],
        }

        for ibin in range(len(bins)):
            s = angle_theta_sum[var_key][ibin]
            n = angle_theta_n[var_key][ibin]

            if n > 0:
                angle_current_rows[var_key]["x_means"].append(float(s) / float(n))
            else:
                angle_current_rows[var_key]["x_means"].append(float("nan"))
            #endif

            per_current_counts_bin = defaultdict(int)
            per_current_n_runs_bin = defaultdict(int)

            for runnum in sorted(run_meta.keys()):
                current_nA = int(run_meta[runnum]["current_nA"])
                count = int(angle_run_counts[var_key][ibin].get(runnum, 0))
                per_current_counts_bin[current_nA] += int(count)
                per_current_n_runs_bin[current_nA] += 1
            #endfor

            rows_bin = []
            for current_nA in sorted(current_charge_totals.keys()):
                total_charge = float(current_charge_totals[current_nA])
                total_counts = int(per_current_counts_bin.get(current_nA, 0))
                n_runs = int(per_current_n_runs_bin.get(current_nA, 0))

                if total_charge <= 0.0:
                    raise RuntimeError(f"Total charge <= 0 for {period_display_name}, current {current_nA} nA")
                #endif

                if total_counts <= 0:
                    continue
                #endif

                counts_per_nC = float(total_counts) / total_charge
                counts_per_nC_err = math.sqrt(float(total_counts)) / total_charge

                rows_bin.append({
                    "period": period_display_name,
                    "period_internal": period_internal_name,
                    "current_nA": int(current_nA),
                    "n_runs": int(n_runs),
                    "counts": int(total_counts),
                    "charge_nC": float(total_charge),
                    "counts_per_nC": float(counts_per_nC),
                    "counts_per_nC_err": float(counts_per_nC_err),
                })
            #endfor

            angle_current_rows[var_key]["rows"].append(rows_bin)
        #endfor
    #endfor

    fine_angle_histograms = {}
    for var_cfg in ANGLE_DEPENDENCE_CONFIG:
        var_key = var_cfg["key"]
        fine_angle_histograms[var_key] = {
            "edges": PLOT_ANGLE_EDGES[var_key].copy(),
            "counts": fine_hist_counts[var_key].copy(),
        }
    #endfor

    fine_angle_histograms["p2_theta"]["counts_phi_gt_60"] = photon_hist_counts_phi_gt_60.copy()
    fine_angle_histograms["p2_theta"]["sector_counts"] = {
        sector: photon_sector_hist_counts[sector].copy()
        for sector in CLAS12_SECTOR_ORDER
    }

    total_valid_charge_nC = float(sum(current_charge_totals.values()))

    return (
        integrated_run_rows,
        integrated_current_rows,
        angle_current_rows,
        skipped_nonpositive_charge_rows,
        fine_angle_histograms,
        total_valid_charge_nC,
    )
#enddef


def process_data_period_worker(args):

    period_display_name, period_internal_name, root_path, run_charge_map = args

    (
        run_rows,
        current_rows,
        angle_current_rows,
        skipped_nonpositive_charge_rows,
        fine_angle_histograms,
        total_valid_charge_nC,
    ) = build_data_period_aggregations(
        period_display_name=period_display_name,
        period_internal_name=period_internal_name,
        root_path=root_path,
        run_charge_map=run_charge_map,
    )

    return {
        "period": period_display_name,
        "period_internal": period_internal_name,
        "run_rows": run_rows,
        "current_rows": current_rows,
        "angle_current_rows": angle_current_rows,
        "skipped_nonpositive_charge_rows": skipped_nonpositive_charge_rows,
        "fine_angle_histograms": fine_angle_histograms,
        "total_valid_charge_nC": total_valid_charge_nC,
        "root_path": root_path,
    }
#enddef


def count_mc_total_and_angle_entries(root_path, tree_name):

    if not os.path.exists(root_path):
        raise RuntimeError(f"ROOT file not found: {root_path}")
    #endif

    root_file = uproot.open(root_path)
    if tree_name not in root_file:
        raise RuntimeError(f"'{tree_name}' tree not found in {root_path}")
    #endif

    tree = root_file[tree_name]

    needed = {"e_theta", "e_phi", "p1_theta", "p2_theta", "p2_phi"}
    missing = [name for name in needed if name not in tree.keys()]
    if len(missing) > 0:
        raise RuntimeError(f"Missing required branches in {root_path}:{tree_name}: {missing}")
    #endif

    iterate_branches = ["e_theta", "e_phi", "p1_theta", "p2_theta", "p2_phi"]

    total_count = 0
    angle_counts = {}

    for var_cfg in ANGLE_DEPENDENCE_CONFIG:
        var_key = var_cfg["key"]
        bins = ANGLE_BINS[var_key]
        angle_counts[var_key] = [0 for _ in bins]
    #endfor

    for arrays in tree.iterate(iterate_branches, library="np", step_size=ITERATE_STEP_SIZE):
        event_mask = apply_event_mask(arrays)

        if not np.any(event_mask):
            continue
        #endif

        total_count += int(np.count_nonzero(event_mask))

        theta_deg_map = {
            "e_theta": np.degrees(arrays["e_theta"][event_mask]),
            "p1_theta": np.degrees(arrays["p1_theta"][event_mask]),
            "p2_theta": np.degrees(arrays["p2_theta"][event_mask]),
        }

        for var_cfg in ANGLE_DEPENDENCE_CONFIG:
            var_key = var_cfg["key"]
            theta_deg = theta_deg_map[var_key]
            bins = ANGLE_BINS[var_key]
            masks = angle_bin_masks(theta_deg, bins)
            store = angle_counts[var_key]

            for ibin, mask in enumerate(masks):
                store[ibin] += int(np.count_nonzero(mask))
            #endfor
        #endfor
    #endfor

    return total_count, angle_counts
#enddef


def process_mc_pair_worker(args):

    period_internal, current_nA, gen_path, rec_path = args
    period_display_name = PERIOD_DISPLAY_FROM_INTERNAL[period_internal]

    n_gen, angle_counts_gen = count_mc_total_and_angle_entries(
        root_path=gen_path,
        tree_name=MC_TREE_NAME,
    )
    n_rec, angle_counts_rec = count_mc_total_and_angle_entries(
        root_path=rec_path,
        tree_name=MC_TREE_NAME,
    )

    rec_fine_histograms = count_mc_fine_angle_histograms(
        root_path=rec_path,
        tree_name=MC_TREE_NAME,
    )

    integrated_row = None
    if n_gen > 0:
        eff = float(n_rec) / float(n_gen)
        eff_err = math.sqrt(float(n_rec)) / float(n_gen) if n_rec > 0 else 0.0

        integrated_row = {
            "period": period_display_name,
            "period_internal": period_internal,
            "current_nA": int(current_nA),
            "n_gen": int(n_gen),
            "n_rec": int(n_rec),
            "efficiency": float(eff),
            "efficiency_err": float(eff_err),
            "gen_file": gen_path,
            "rec_file": rec_path,
        }
    #endif

    angle_rows_by_variable = {}
    for var_cfg in ANGLE_DEPENDENCE_CONFIG:
        var_key = var_cfg["key"]
        bins = ANGLE_BINS[var_key]
        rows = []

        for ibin in range(len(bins)):
            n_gen_bin = int(angle_counts_gen[var_key][ibin])
            n_rec_bin = int(angle_counts_rec[var_key][ibin])

            if n_gen_bin <= 0:
                rows.append(None)
                continue
            #endif
            if n_rec_bin <= 0:
                rows.append(None)
                continue
            #endif

            eff_bin = float(n_rec_bin) / float(n_gen_bin)
            eff_err_bin = math.sqrt(float(n_rec_bin)) / float(n_gen_bin)

            rows.append({
                "period": period_display_name,
                "period_internal": period_internal,
                "current_nA": int(current_nA),
                "n_gen": int(n_gen_bin),
                "n_rec": int(n_rec_bin),
                "efficiency": float(eff_bin),
                "efficiency_err": float(eff_err_bin),
                "gen_file": gen_path,
                "rec_file": rec_path,
            })
        #endfor

        angle_rows_by_variable[var_key] = rows
    #endfor

    return {
        "period": period_display_name,
        "period_internal": period_internal,
        "current_nA": int(current_nA),
        "integrated_row": integrated_row,
        "angle_rows_by_variable": angle_rows_by_variable,
        "rec_fine_histograms": rec_fine_histograms,
    }
#enddef


def build_mc_aggregation(mc_dir):

    if not os.path.isdir(mc_dir):
        raise RuntimeError(f"MC directory not found: {mc_dir}")
    #endif

    grouped = {}
    skipped_non_mc_files = []

    for basename in sorted(os.listdir(mc_dir)):
        full_path = os.path.join(mc_dir, basename)

        if os.path.isdir(full_path):
            continue
        #endif

        if not basename.endswith(".root"):
            continue
        #endif

        if not is_candidate_mc_filename(basename):
            skipped_non_mc_files.append(basename)
            continue
        #endif

        kind, period_internal, current_nA, beam_energy_token = parse_mc_filename(basename)

        if period_internal not in PERIOD_DISPLAY_FROM_INTERNAL:
            raise RuntimeError(f"Unknown MC period in filename: {basename}")
        #endif

        key = (period_internal, current_nA)
        if key not in grouped:
            grouped[key] = {"gen": None, "rec": None}
        #endif

        if grouped[key][kind] is not None:
            raise RuntimeError(f"Duplicate MC file for {key}: {basename}")
        #endif

        grouped[key][kind] = full_path
    #endfor

    if len(skipped_non_mc_files) > 0:
        print("")
        print("Skipping non-MC ROOT files found in MC directory:")
        for name in skipped_non_mc_files:
            print(f"  {name}")
        #endfor
    #endif

    tasks = []
    for key in sorted(grouped.keys(), key=lambda item: (PERIOD_ORDER.index(PERIOD_DISPLAY_FROM_INTERNAL[item[0]]), item[1])):
        period_internal, current_nA = key
        gen_path = grouped[key]["gen"]
        rec_path = grouped[key]["rec"]

        if gen_path is None or rec_path is None:
            raise RuntimeError(f"Missing gen/rec MC pair for period={period_internal}, current={current_nA} nA")
        #endif

        tasks.append((period_internal, current_nA, gen_path, rec_path))
    #endfor

    integrated_mc_rows = []
    angle_mc_rows_by_variable = {}
    rec_fine_angle_histograms_by_period = {}

    for var_cfg in ANGLE_DEPENDENCE_CONFIG:
        var_key = var_cfg["key"]
        bins = ANGLE_BINS[var_key]
        angle_mc_rows_by_variable[var_key] = {
            "bins": bins,
            "rows": {period_name: [[] for _ in bins] for period_name in PERIOD_ORDER},
        }
    #endfor

    for period_name in PERIOD_ORDER:
        rec_fine_angle_histograms_by_period[period_name] = {}
        for var_cfg in ANGLE_DEPENDENCE_CONFIG:
            var_key = var_cfg["key"]
            rec_fine_angle_histograms_by_period[period_name][var_key] = {
                "edges": PLOT_ANGLE_EDGES[var_key].copy(),
                "counts": np.zeros(len(PLOT_ANGLE_EDGES[var_key]) - 1, dtype=np.int64),
            }
        #endfor

        rec_fine_angle_histograms_by_period[period_name]["p2_theta"]["counts_phi_gt_60"] = np.zeros(
            len(PLOT_ANGLE_EDGES["p2_theta"]) - 1,
            dtype=np.int64,
        )
        rec_fine_angle_histograms_by_period[period_name]["p2_theta"]["sector_counts"] = {
            sector: np.zeros(len(PLOT_ANGLE_EDGES["p2_theta"]) - 1, dtype=np.int64)
            for sector in CLAS12_SECTOR_ORDER
        }
    #endfor

    n_workers = min(MAX_WORKERS, max(1, len(tasks)))

    with cf.ProcessPoolExecutor(max_workers=n_workers) as ex:
        futures = [ex.submit(process_mc_pair_worker, task) for task in tasks]

        for fut in cf.as_completed(futures):
            result = fut.result()

            if result["integrated_row"] is not None:
                integrated_mc_rows.append(result["integrated_row"])
            #endif

            period_name = result["period"]

            for var_cfg in ANGLE_DEPENDENCE_CONFIG:
                var_key = var_cfg["key"]
                rows = result["angle_rows_by_variable"][var_key]
                store = angle_mc_rows_by_variable[var_key]["rows"][period_name]

                for ibin, row in enumerate(rows):
                    if row is not None:
                        store[ibin].append(row)
                    #endif
                #endfor
            #endfor

            rec_hist = result["rec_fine_histograms"]
            for var_cfg in ANGLE_DEPENDENCE_CONFIG:
                var_key = var_cfg["key"]
                rec_fine_angle_histograms_by_period[period_name][var_key]["counts"] += np.asarray(
                    rec_hist[var_key]["counts"],
                    dtype=np.int64,
                )
            #endfor

            rec_fine_angle_histograms_by_period[period_name]["p2_theta"]["counts_phi_gt_60"] += np.asarray(
                rec_hist["p2_theta"]["counts_phi_gt_60"],
                dtype=np.int64,
            )

            for sector in CLAS12_SECTOR_ORDER:
                rec_fine_angle_histograms_by_period[period_name]["p2_theta"]["sector_counts"][sector] += np.asarray(
                    rec_hist["p2_theta"]["sector_counts"][sector],
                    dtype=np.int64,
                )
            #endfor
        #endfor
    #endwith

    integrated_mc_rows = sorted(
        integrated_mc_rows,
        key=lambda row: (PERIOD_ORDER.index(row["period"]), row["current_nA"])
    )

    for var_cfg in ANGLE_DEPENDENCE_CONFIG:
        var_key = var_cfg["key"]
        for period_name in PERIOD_ORDER:
            bins_rows = angle_mc_rows_by_variable[var_key]["rows"][period_name]
            for ibin in range(len(bins_rows)):
                bins_rows[ibin] = sorted(bins_rows[ibin], key=lambda row: row["current_nA"])
            #endfor
        #endfor
    #endfor

    return integrated_mc_rows, angle_mc_rows_by_variable, rec_fine_angle_histograms_by_period
#enddef


def period_points_from_current_rows(current_rows, period_name):

    rows = [r for r in current_rows if r["period"] == period_name]
    rows = sorted(rows, key=lambda r: r["current_nA"])

    if len(rows) == 0:
        return np.asarray([]), np.asarray([]), np.asarray([]), []
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
        return np.asarray([]), np.asarray([]), np.asarray([]), []
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


def write_skipped_run_table_csv(path, skipped_rows):

    fieldnames = [
        "period",
        "period_internal",
        "runnum",
        "current_nA",
        "charge_nC",
        "reason",
    ]

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in skipped_rows:
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


def write_period_summary_csv(path, summary_rows):

    fieldnames = [
        "period",
        "period_internal",
        "weighted_data_rel_percent",
        "weighted_data_rel_fraction",
        "weighted_data_rel_stat_err_percent",
        "weighted_data_rel_stat_err_fraction",
        "data_rel_at_ref_percent",
        "data_rel_at_ref_fraction",
        "data_rel_at_ref_stat_err_percent",
        "data_rel_at_ref_stat_err_fraction",
        "mc_rel_at_ref_percent",
        "mc_rel_at_ref_fraction",
        "mc_rel_at_ref_stat_err_percent",
        "mc_rel_at_ref_stat_err_fraction",
        "data_over_mc_at_ref_percent",
        "data_over_mc_at_ref_fraction",
        "data_over_mc_at_ref_stat_err_percent",
        "data_over_mc_at_ref_stat_err_fraction",
        "divisor_to_divide_by_percent",
        "divisor_to_divide_by_fraction",
        "divisor_to_divide_by_stat_err_percent",
        "divisor_to_divide_by_stat_err_fraction",
        "applied_scale_percent",
        "applied_scale_fraction",
        "applied_scale_stat_err_percent",
        "applied_scale_stat_err_fraction",
    ]

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in summary_rows:
            writer.writerow(row)
        #endfor
    #endwith
#enddef


def write_angle_summary_csv(path, summary_rows):

    fieldnames = [
        "period",
        "period_internal",
        "angle_variable",
        "angle_bin_index",
        "angle_min_deg",
        "angle_max_deg",
        "angle_x_mean_deg",
        "weighted_data_rel_percent",
        "weighted_data_rel_fraction",
        "weighted_data_rel_stat_err_percent",
        "weighted_data_rel_stat_err_fraction",
        "data_rel_at_ref_percent",
        "data_rel_at_ref_fraction",
        "data_rel_at_ref_stat_err_percent",
        "data_rel_at_ref_stat_err_fraction",
        "mc_rel_at_ref_percent",
        "mc_rel_at_ref_fraction",
        "mc_rel_at_ref_stat_err_percent",
        "mc_rel_at_ref_stat_err_fraction",
        "data_over_mc_at_ref_percent",
        "data_over_mc_at_ref_fraction",
        "data_over_mc_at_ref_stat_err_percent",
        "data_over_mc_at_ref_stat_err_fraction",
        "divisor_to_divide_by_percent",
        "divisor_to_divide_by_fraction",
        "divisor_to_divide_by_stat_err_percent",
        "divisor_to_divide_by_stat_err_fraction",
        "applied_scale_percent",
        "applied_scale_fraction",
        "applied_scale_stat_err_percent",
        "applied_scale_stat_err_fraction",
    ]

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in summary_rows:
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


def compute_relative_value_at_current(current_nA, fit_result):

    m = fit_result["m"]
    b = fit_result["b"]

    return (m * float(current_nA) + b) / b
#enddef


def compute_relative_value_error_at_current(current_nA, fit_result):

    m = float(fit_result["m"])
    b = float(fit_result["b"])
    sm = float(fit_result["sm"])
    sb = float(fit_result["sb"])
    cov_mb = float(fit_result["cov_mb"])

    if not np.isfinite(m) or not np.isfinite(b) or b == 0.0:
        return float("nan")
    #endif

    if not np.isfinite(sm) or not np.isfinite(sb) or not np.isfinite(cov_mb):
        return float("nan")
    #endif

    var_m = sm * sm
    var_b = sb * sb
    I = float(current_nA)

    dr_dm = I / b
    dr_db = -(m * I) / (b * b)

    var_r = dr_dm * dr_dm * var_m + dr_db * dr_db * var_b + 2.0 * dr_dm * dr_db * cov_mb

    if var_r < 0.0 and abs(var_r) < 1.0e-15:
        var_r = 0.0
    #endif
    if var_r < 0.0:
        return float("nan")
    #endif

    return math.sqrt(var_r)
#enddef


def divide_with_error(num, num_err, den, den_err):

    if not np.isfinite(num) or not np.isfinite(den) or den == 0.0:
        return float("nan"), float("nan")
    #endif

    val = num / den

    rel2 = 0.0
    if np.isfinite(num_err) and num != 0.0:
        rel2 += (num_err / num) ** 2
    #endif
    if np.isfinite(den_err) and den != 0.0:
        rel2 += (den_err / den) ** 2
    #endif

    err = abs(val) * math.sqrt(rel2)
    return val, err
#enddef


def reciprocal_with_error(val, val_err):

    if not np.isfinite(val) or val == 0.0:
        return float("nan"), float("nan")
    #endif

    out = 1.0 / val
    if not np.isfinite(val_err):
        return out, float("nan")
    #endif

    out_err = abs(val_err / (val * val))
    return out, out_err
#enddef


def compute_weighted_data_rel_and_error(rows, fit_result):

    m = float(fit_result["m"])
    b = float(fit_result["b"])
    sm = float(fit_result["sm"])
    sb = float(fit_result["sb"])
    cov_mb = float(fit_result["cov_mb"])

    if not np.isfinite(m) or not np.isfinite(b) or b == 0.0:
        return float("nan"), float("nan")
    #endif
    if not np.isfinite(sm) or not np.isfinite(sb) or not np.isfinite(cov_mb):
        return float("nan"), float("nan")
    #endif

    total_counts = float(sum(int(r["counts"]) for r in rows))
    if total_counts <= 0.0:
        return float("nan"), float("nan")
    #endif

    Ibar = 0.0
    for row in rows:
        w = float(row["counts"]) / total_counts
        Ibar += w * float(row["current_nA"])
    #endfor

    value = 1.0 + (m * Ibar) / b

    var_m = sm * sm
    var_b = sb * sb

    dr_dm = Ibar / b
    dr_db = -(m * Ibar) / (b * b)

    var_r = dr_dm * dr_dm * var_m + dr_db * dr_db * var_b + 2.0 * dr_dm * dr_db * cov_mb

    if var_r < 0.0 and abs(var_r) < 1.0e-15:
        var_r = 0.0
    #endif
    if var_r < 0.0:
        return value, float("nan")
    #endif

    return value, math.sqrt(var_r)
#enddef


def build_residual_correction_rows(all_current_rows, data_fit_results, mc_fit_results):

    residual_rows = []

    for period in PERIOD_ORDER:
        period_internal = PERIOD_INTERNAL_FROM_DISPLAY[period]

        if period not in data_fit_results or period not in mc_fit_results:
            continue
        #endif

        data_fit = data_fit_results[period]
        mc_fit = mc_fit_results[period]

        if not np.isfinite(data_fit["b"]) or not np.isfinite(mc_fit["b"]):
            continue
        #endif

        ref_current = get_mc_reference_current(period_internal)
        mc_rel_at_ref = compute_relative_value_at_current(ref_current, mc_fit)

        data_rows = [r for r in all_current_rows if r["period"] == period]
        data_rows = sorted(data_rows, key=lambda r: r["current_nA"])

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


def build_period_summary_rows(current_rows, data_fit_results, mc_fit_results):

    summary_rows = []

    for period in PERIOD_ORDER:
        period_internal = PERIOD_INTERNAL_FROM_DISPLAY[period]

        if period not in data_fit_results or period not in mc_fit_results:
            continue
        #endif

        data_fit = data_fit_results[period]
        mc_fit = mc_fit_results[period]

        if not np.isfinite(data_fit["b"]) or not np.isfinite(mc_fit["b"]):
            continue
        #endif

        rows = [r for r in current_rows if r["period"] == period]
        rows = sorted(rows, key=lambda r: r["current_nA"])

        if len(rows) == 0:
            continue
        #endif

        weighted_data_rel, weighted_data_rel_err = compute_weighted_data_rel_and_error(rows, data_fit)

        ref_current = get_mc_reference_current(period_internal)

        data_rel_at_ref = compute_relative_value_at_current(ref_current, data_fit)
        data_rel_at_ref_err = compute_relative_value_error_at_current(ref_current, data_fit)

        mc_rel_at_ref = compute_relative_value_at_current(ref_current, mc_fit)
        mc_rel_at_ref_err = compute_relative_value_error_at_current(ref_current, mc_fit)

        data_over_mc_at_ref, data_over_mc_at_ref_err = divide_with_error(
            data_rel_at_ref, data_rel_at_ref_err,
            mc_rel_at_ref, mc_rel_at_ref_err,
        )

        divisor_to_divide_by, divisor_to_divide_by_err = divide_with_error(
            weighted_data_rel, weighted_data_rel_err,
            mc_rel_at_ref, mc_rel_at_ref_err,
        )

        applied_scale, applied_scale_err = reciprocal_with_error(
            divisor_to_divide_by,
            divisor_to_divide_by_err,
        )

        summary_rows.append({
            "period": period,
            "period_internal": period_internal,

            "weighted_data_rel_percent": 100.0 * weighted_data_rel,
            "weighted_data_rel_fraction": weighted_data_rel,
            "weighted_data_rel_stat_err_percent": 100.0 * weighted_data_rel_err,
            "weighted_data_rel_stat_err_fraction": weighted_data_rel_err,

            "data_rel_at_ref_percent": 100.0 * data_rel_at_ref,
            "data_rel_at_ref_fraction": data_rel_at_ref,
            "data_rel_at_ref_stat_err_percent": 100.0 * data_rel_at_ref_err,
            "data_rel_at_ref_stat_err_fraction": data_rel_at_ref_err,

            "mc_rel_at_ref_percent": 100.0 * mc_rel_at_ref,
            "mc_rel_at_ref_fraction": mc_rel_at_ref,
            "mc_rel_at_ref_stat_err_percent": 100.0 * mc_rel_at_ref_err,
            "mc_rel_at_ref_stat_err_fraction": mc_rel_at_ref_err,

            "data_over_mc_at_ref_percent": 100.0 * data_over_mc_at_ref,
            "data_over_mc_at_ref_fraction": data_over_mc_at_ref,
            "data_over_mc_at_ref_stat_err_percent": 100.0 * data_over_mc_at_ref_err,
            "data_over_mc_at_ref_stat_err_fraction": data_over_mc_at_ref_err,

            "divisor_to_divide_by_percent": 100.0 * divisor_to_divide_by,
            "divisor_to_divide_by_fraction": divisor_to_divide_by,
            "divisor_to_divide_by_stat_err_percent": 100.0 * divisor_to_divide_by_err,
            "divisor_to_divide_by_stat_err_fraction": divisor_to_divide_by_err,

            "applied_scale_percent": 100.0 * applied_scale,
            "applied_scale_fraction": applied_scale,
            "applied_scale_stat_err_percent": 100.0 * applied_scale_err,
            "applied_scale_stat_err_fraction": applied_scale_err,
        })
    #endfor

    return summary_rows
#enddef


def build_angle_summary_rows(angle_rows_by_period, mc_rows_by_variable, var_cfg):

    summary_rows = []
    var_key = var_cfg["key"]

    for period_name in PERIOD_ORDER:
        bins = angle_rows_by_period[period_name]["bins"]
        current_rows_list = angle_rows_by_period[period_name]["rows"]
        mc_rows_list = mc_rows_by_variable["rows"][period_name]
        x_means = angle_rows_by_period[period_name]["x_means"]

        for ibin, (low, high) in enumerate(bins):
            current_rows = current_rows_list[ibin]
            mc_rows = mc_rows_list[ibin]

            data_fit_results = {}
            mc_fit_results = {}

            for loop_period in PERIOD_ORDER:
                if loop_period == period_name:
                    xd, yd, syd, rowsd = period_points_from_current_rows(current_rows, loop_period)
                    xm, ym, sym, rowsm = period_points_from_mc_rows(mc_rows, loop_period)
                else:
                    xd = np.asarray([])
                    yd = np.asarray([])
                    syd = np.asarray([])
                    rowsd = []
                    xm = np.asarray([])
                    ym = np.asarray([])
                    sym = np.asarray([])
                    rowsm = []
                #endif

                data_fit_results[loop_period] = weighted_linear_fit(xd, yd, syd)
                mc_fit_results[loop_period] = weighted_linear_fit(xm, ym, sym)
            #endfor

            period_rows = build_period_summary_rows(current_rows, data_fit_results, mc_fit_results)

            for row in period_rows:
                row_out = dict(row)
                row_out["angle_variable"] = var_key
                row_out["angle_bin_index"] = int(ibin)
                row_out["angle_min_deg"] = float(low)
                row_out["angle_max_deg"] = float(high)
                row_out["angle_x_mean_deg"] = float(x_means[ibin])
                summary_rows.append(row_out)
            #endfor
        #endfor
    #endfor

    return summary_rows
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
    ax.set_ylim(0.70, 1.10)
    ax.set_xlabel("Beam current (nA)")
    ax.set_ylabel("Data/MC")
    ax.grid(True, alpha=0.3)
    ax.axhline(1.0, color="0.5", linestyle="--", linewidth=1.0)
    ax.set_yticks([0.70, 0.80, 0.90, 1.00, 1.10])
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


def plot_four_panel_set(output_dir, tag, data_current_rows, mc_rows, data_fit_results, mc_fit_results, period_color, title_suffix=""):

    os.makedirs(output_dir, exist_ok=True)
    xfit = np.linspace(0.0, 80.0, 300)
    band_alpha = 0.20

    fig_a, axes_a = create_simple_2x3_figure()

    for i, period in enumerate(PERIOD_ORDER):
        ax = axes_a[i]
        c = period_color[period]

        xd, yd, syd, rowsd = period_points_from_current_rows(data_current_rows, period)

        if len(rowsd) > 0:
            ax.errorbar(xd, yd, yerr=syd, fmt="o", capsize=3, color=c, label="Data")
        #endif

        frd = data_fit_results.get(period, None)
        if frd is not None and np.isfinite(frd["m"]) and np.isfinite(frd["b"]):
            data_fit = frd["m"] * xfit + frd["b"]
            data_fit_lo = (frd["m"] - frd["sm"]) * xfit + frd["b"]
            data_fit_hi = (frd["m"] + frd["sm"]) * xfit + frd["b"]
            ax.fill_between(xfit, data_fit_lo, data_fit_hi, color=c, alpha=band_alpha, linewidth=0)
            ax.plot(xfit, data_fit, color=c)
        #endif

        ax.set_title(period + title_suffix)
        style_absolute_axis(ax, "Counts / accumulated charge (1/nC)")
        if len(rowsd) > 0:
            ax.legend(frameon=True)
        #endif
        add_reference_current_text(ax, period)
    #endfor

    ax = axes_a[5]
    ax.set_title("All periods" + title_suffix)
    style_absolute_axis(ax, "Counts / accumulated charge (1/nC)")

    for period in PERIOD_ORDER:
        c = period_color[period]
        xd, yd, syd, rowsd = period_points_from_current_rows(data_current_rows, period)

        if len(rowsd) > 0:
            ax.errorbar(xd, yd, yerr=syd, fmt="o", capsize=3, color=c, label=period)
        #endif

        frd = data_fit_results.get(period, None)
        if frd is not None and np.isfinite(frd["m"]) and np.isfinite(frd["b"]):
            data_fit = frd["m"] * xfit + frd["b"]
            data_fit_lo = (frd["m"] - frd["sm"]) * xfit + frd["b"]
            data_fit_hi = (frd["m"] + frd["sm"]) * xfit + frd["b"]
            ax.fill_between(xfit, data_fit_lo, data_fit_hi, color=c, alpha=0.08, linewidth=0)
            ax.plot(xfit, data_fit, color=c)
        #endif
    #endfor

    ax.legend(frameon=True, fontsize=9)

    out_a = os.path.join(output_dir, f"{tag}_counts_per_charge_data.png")
    fig_a.savefig(out_a, dpi=200)
    plt.close(fig_a)

    fig_b, axes_b = create_simple_2x3_figure()

    for i, period in enumerate(PERIOD_ORDER):
        ax = axes_b[i]
        c = period_color[period]

        xd, yd, syd, rowsd = period_points_from_current_rows(data_current_rows, period)
        frd = data_fit_results.get(period, None)

        if frd is not None and np.isfinite(frd["b"]):
            if len(rowsd) > 0:
                data_pct = 100.0 * (yd / frd["b"])
                data_pct_err = 100.0 * np.sqrt((syd / frd["b"]) ** 2 + ((yd * frd["sb"]) / (frd["b"] * frd["b"])) ** 2)
                ax.errorbar(xd, data_pct, yerr=data_pct_err, fmt="o", capsize=3, color=c, label="Data")
            #endif

            pct_fit, pct_fit_lo, pct_fit_hi = compute_percent_curve(xfit, frd)
            ax.fill_between(xfit, pct_fit_lo, pct_fit_hi, color=c, alpha=band_alpha, linewidth=0)
            ax.plot(xfit, pct_fit, color=c)
        #endif

        ax.set_title(period + title_suffix)
        style_percent_axis(ax, "Efficiency relative to fitted 0 nA (%)")
        ax.set_xlabel("Beam current (nA)")
        if len(rowsd) > 0:
            ax.legend(frameon=True)
        #endif
        add_reference_current_text(ax, period)
    #endfor

    ax = axes_b[5]
    ax.set_title("All periods" + title_suffix)
    style_percent_axis(ax, "Efficiency relative to fitted 0 nA (%)")
    ax.set_xlabel("Beam current (nA)")

    for period in PERIOD_ORDER:
        c = period_color[period]
        xd, yd, syd, rowsd = period_points_from_current_rows(data_current_rows, period)
        frd = data_fit_results.get(period, None)

        if frd is not None and np.isfinite(frd["b"]):
            if len(rowsd) > 0:
                data_pct = 100.0 * (yd / frd["b"])
                data_pct_err = 100.0 * np.sqrt((syd / frd["b"]) ** 2 + ((yd * frd["sb"]) / (frd["b"] * frd["b"])) ** 2)
                ax.errorbar(xd, data_pct, yerr=data_pct_err, fmt="o", capsize=3, color=c, label=period)
            #endif

            pct_fit, pct_fit_lo, pct_fit_hi = compute_percent_curve(xfit, frd)
            ax.fill_between(xfit, pct_fit_lo, pct_fit_hi, color=c, alpha=0.08, linewidth=0)
            ax.plot(xfit, pct_fit, color=c)
        #endif
    #endfor

    ax.legend(frameon=True, fontsize=9)

    out_b = os.path.join(output_dir, f"{tag}_percent_of_zero_data.png")
    fig_b.savefig(out_b, dpi=200)
    plt.close(fig_b)

    fig_c, axes_c = create_simple_2x3_figure()

    for i, period in enumerate(PERIOD_ORDER):
        ax = axes_c[i]
        c = period_color[period]

        xd, yd, syd, rowsd = period_points_from_current_rows(data_current_rows, period)
        xm, ym, sym, rowsm = period_points_from_mc_rows(mc_rows, period)

        if len(rowsd) > 0:
            ax.errorbar(xd, yd, yerr=syd, fmt="o", capsize=3, color=c, label="Data")
        #endif
        if len(rowsm) > 0:
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
        #endif

        frd = data_fit_results.get(period, None)
        if frd is not None and np.isfinite(frd["m"]) and np.isfinite(frd["b"]):
            data_fit = frd["m"] * xfit + frd["b"]
            data_fit_lo = (frd["m"] - frd["sm"]) * xfit + frd["b"]
            data_fit_hi = (frd["m"] + frd["sm"]) * xfit + frd["b"]
            ax.fill_between(xfit, data_fit_lo, data_fit_hi, color=c, alpha=0.12, linewidth=0)
            ax.plot(xfit, data_fit, color=c)
        #endif

        frm = mc_fit_results.get(period, None)
        if frm is not None and np.isfinite(frm["m"]) and np.isfinite(frm["b"]):
            mc_fit = frm["m"] * xfit + frm["b"]
            mc_fit_lo = (frm["m"] - frm["sm"]) * xfit + frm["b"]
            mc_fit_hi = (frm["m"] + frm["sm"]) * xfit + frm["b"]
            ax.fill_between(xfit, mc_fit_lo, mc_fit_hi, color=c, alpha=0.08, linewidth=0)
            ax.plot(xfit, mc_fit, color=c, linestyle="--")
        #endif

        ax.set_title(period + title_suffix)
        style_absolute_axis(ax, "Absolute value")
        if len(rowsd) > 0 or len(rowsm) > 0:
            ax.legend(frameon=True)
        #endif
        add_reference_current_text(ax, period)
    #endfor

    ax = axes_c[5]
    ax.set_title("All periods" + title_suffix)
    style_absolute_axis(ax, "Absolute value")

    for period in PERIOD_ORDER:
        c = period_color[period]
        xd, yd, syd, rowsd = period_points_from_current_rows(data_current_rows, period)
        xm, ym, sym, rowsm = period_points_from_mc_rows(mc_rows, period)

        if len(rowsd) > 0:
            ax.errorbar(xd, yd, yerr=syd, fmt="o", capsize=3, color=c, label=f"{period} data")
        #endif
        if len(rowsm) > 0:
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
        #endif

        frd = data_fit_results.get(period, None)
        if frd is not None and np.isfinite(frd["m"]) and np.isfinite(frd["b"]):
            data_fit = frd["m"] * xfit + frd["b"]
            data_fit_lo = (frd["m"] - frd["sm"]) * xfit + frd["b"]
            data_fit_hi = (frd["m"] + frd["sm"]) * xfit + frd["b"]
            ax.fill_between(xfit, data_fit_lo, data_fit_hi, color=c, alpha=0.08, linewidth=0)
            ax.plot(xfit, data_fit, color=c)
        #endif

        frm = mc_fit_results.get(period, None)
        if frm is not None and np.isfinite(frm["m"]) and np.isfinite(frm["b"]):
            mc_fit = frm["m"] * xfit + frm["b"]
            mc_fit_lo = (frm["m"] - frm["sm"]) * xfit + frm["b"]
            mc_fit_hi = (frm["m"] + frm["sm"]) * xfit + frm["b"]
            ax.fill_between(xfit, mc_fit_lo, mc_fit_hi, color=c, alpha=0.05, linewidth=0)
            ax.plot(xfit, mc_fit, color=c, linestyle="--")
        #endif
    #endfor

    ax.legend(frameon=True, fontsize=8)

    out_c = os.path.join(output_dir, f"{tag}_absolute_data_vs_mc.png")
    fig_c.savefig(out_c, dpi=200)
    plt.close(fig_c)

    fig_d, top_axes_d, bottom_axes_d = create_doublepad_2x3_figure()

    for i, period in enumerate(PERIOD_ORDER):
        ax_top = top_axes_d[i]
        ax_bot = bottom_axes_d[i]
        c = period_color[period]

        xd, yd, syd, rowsd = period_points_from_current_rows(data_current_rows, period)
        xm, ym, sym, rowsm = period_points_from_mc_rows(mc_rows, period)

        frd = data_fit_results.get(period, None)
        frm = mc_fit_results.get(period, None)

        if frd is not None and np.isfinite(frd["b"]):
            if len(rowsd) > 0:
                data_pct = 100.0 * (yd / frd["b"])
                data_pct_err = 100.0 * np.sqrt((syd / frd["b"]) ** 2 + ((yd * frd["sb"]) / (frd["b"] * frd["b"])) ** 2)
                ax_top.errorbar(xd, data_pct, yerr=data_pct_err, fmt="o", capsize=3, color=c, label="Data")
            #endif

            data_pct_fit, data_pct_fit_lo, data_pct_fit_hi = compute_percent_curve(xfit, frd)
            ax_top.fill_between(xfit, data_pct_fit_lo, data_pct_fit_hi, color=c, alpha=0.12, linewidth=0)
            ax_top.plot(xfit, data_pct_fit, color=c)
        #endif

        if frm is not None and np.isfinite(frm["b"]):
            if len(rowsm) > 0:
                mc_pct = 100.0 * (ym / frm["b"])
                mc_pct_err = 100.0 * np.sqrt((sym / frm["b"]) ** 2 + ((ym * frm["sb"]) / (frm["b"] * frm["b"])) ** 2)
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
            #endif

            mc_pct_fit, mc_pct_fit_lo, mc_pct_fit_hi = compute_percent_curve(xfit, frm)
            ax_top.fill_between(xfit, mc_pct_fit_lo, mc_pct_fit_hi, color=c, alpha=0.08, linewidth=0)
            ax_top.plot(xfit, mc_pct_fit, color=c, linestyle="--")
        #endif

        ax_top.set_title(period + title_suffix)
        style_percent_axis(ax_top, "Efficiency relative to fitted 0 nA (%)")
        ax_top.set_xlabel("")
        if len(rowsd) > 0 or len(rowsm) > 0:
            ax_top.legend(frameon=True)
        #endif
        add_reference_current_text(ax_top, period)

        if frd is not None and frm is not None and np.isfinite(frd["b"]) and np.isfinite(frm["b"]) and len(rowsd) > 0 and len(rowsm) > 0:
            data_pct = 100.0 * (yd / frd["b"])
            data_pct_err = 100.0 * np.sqrt((syd / frd["b"]) ** 2 + ((yd * frd["sb"]) / (frd["b"] * frd["b"])) ** 2)

            mc_pct = 100.0 * (ym / frm["b"])
            mc_pct_err = 100.0 * np.sqrt((sym / frm["b"]) ** 2 + ((ym * frm["sb"]) / (frm["b"] * frm["b"])) ** 2)

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
        #endif

        style_ratio_axis(ax_bot)
    #endfor

    ax_top = top_axes_d[5]
    ax_bot = bottom_axes_d[5]
    ax_top.set_title("All periods" + title_suffix)
    style_percent_axis(ax_top, "Efficiency relative to fitted 0 nA (%)")
    ax_top.set_xlabel("")
    style_ratio_axis(ax_bot)

    for period in PERIOD_ORDER:
        c = period_color[period]

        xd, yd, syd, rowsd = period_points_from_current_rows(data_current_rows, period)
        xm, ym, sym, rowsm = period_points_from_mc_rows(mc_rows, period)

        frd = data_fit_results.get(period, None)
        frm = mc_fit_results.get(period, None)

        if frd is not None and np.isfinite(frd["b"]):
            if len(rowsd) > 0:
                data_pct = 100.0 * (yd / frd["b"])
                data_pct_err = 100.0 * np.sqrt((syd / frd["b"]) ** 2 + ((yd * frd["sb"]) / (frd["b"] * frd["b"])) ** 2)
                ax_top.errorbar(xd, data_pct, yerr=data_pct_err, fmt="o", capsize=3, color=c, label=f"{period} data")
            #endif

            data_pct_fit, data_pct_fit_lo, data_pct_fit_hi = compute_percent_curve(xfit, frd)
            ax_top.fill_between(xfit, data_pct_fit_lo, data_pct_fit_hi, color=c, alpha=0.08, linewidth=0)
            ax_top.plot(xfit, data_pct_fit, color=c)
        #endif

        if frm is not None and np.isfinite(frm["b"]):
            if len(rowsm) > 0:
                mc_pct = 100.0 * (ym / frm["b"])
                mc_pct_err = 100.0 * np.sqrt((sym / frm["b"]) ** 2 + ((ym * frm["sb"]) / (frm["b"] * frm["b"])) ** 2)
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
            #endif

            mc_pct_fit, mc_pct_fit_lo, mc_pct_fit_hi = compute_percent_curve(xfit, frm)
            ax_top.fill_between(xfit, mc_pct_fit_lo, mc_pct_fit_hi, color=c, alpha=0.05, linewidth=0)
            ax_top.plot(xfit, mc_pct_fit, color=c, linestyle="--")
        #endif

        if frd is not None and frm is not None and np.isfinite(frd["b"]) and np.isfinite(frm["b"]) and len(rowsd) > 0 and len(rowsm) > 0:
            data_pct = 100.0 * (yd / frd["b"])
            data_pct_err = 100.0 * np.sqrt((syd / frd["b"]) ** 2 + ((yd * frd["sb"]) / (frd["b"] * frd["b"])) ** 2)

            mc_pct = 100.0 * (ym / frm["b"])
            mc_pct_err = 100.0 * np.sqrt((sym / frm["b"]) ** 2 + ((ym * frm["sb"]) / (frm["b"] * frm["b"])) ** 2)

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
        #endif
    #endfor

    ax_top.legend(frameon=True, fontsize=8)

    out_d = os.path.join(output_dir, f"{tag}_percent_of_zero_data_vs_mc_with_ratio.png")
    fig_d.savefig(out_d, dpi=200)
    plt.close(fig_d)

    return out_a, out_b, out_c, out_d
#enddef


def make_fit_result_map_for_data(current_rows):

    fit_results = {}

    for period in PERIOD_ORDER:
        x, y, sy, rows = period_points_from_current_rows(current_rows, period)
        fit_results[period] = weighted_linear_fit(x, y, sy)
    #endfor

    return fit_results
#enddef


def make_fit_result_map_for_mc(mc_rows):

    fit_results = {}

    for period in PERIOD_ORDER:
        x, y, sy, rows = period_points_from_mc_rows(mc_rows, period)
        fit_results[period] = weighted_linear_fit(x, y, sy)
    #endfor

    return fit_results
#enddef


def make_angle_hist_overlay_plot(angle_histograms_by_period, output_dir, period_color, var_cfg, title_suffix=""):

    os.makedirs(output_dir, exist_ok=True)

    var_key = var_cfg["key"]
    axis_label = get_angle_axis_label(var_key)
    saved_paths = []

    fig, axes = plt.subplots(1, 2, figsize=(16, 6), constrained_layout=True)
    ax_left = axes[0]
    ax_right = axes[1]

    for period in PERIOD_ORDER:
        if period not in angle_histograms_by_period:
            continue
        #endif

        entry = angle_histograms_by_period[period]
        edges = np.asarray(entry["edges"], dtype=float)
        counts = np.asarray(entry["counts"], dtype=float)
        charge_nC = float(entry["charge_nC"])

        if charge_nC <= 0.0:
            continue
        #endif

        counts_per_nC = counts / charge_nC

        ax_left.stairs(
            counts_per_nC,
            edges,
            label=period,
            color=period_color[period],
            linewidth=1.5,
        )

        bin_widths = np.diff(edges)
        integral = float(np.sum(counts_per_nC * bin_widths))

        if integral > 0.0:
            normalized = counts_per_nC / integral
            ax_right.stairs(
                normalized,
                edges,
                label=period,
                color=period_color[period],
                linewidth=1.5,
            )
        #endif
    #endfor

    left_title = f"{var_cfg['display_name']} histograms: counts / accumulated charge"
    right_title = f"{var_cfg['display_name']} histograms: normalized to integral"

    if title_suffix != "":
        left_title += title_suffix
        right_title += title_suffix
    #endif

    ax_left.set_title(left_title)
    ax_left.set_xlabel(axis_label)
    ax_left.set_ylabel("Counts / accumulated charge (1/nC)")
    ax_left.set_xlim(float(var_cfg["min_deg"]), float(var_cfg["max_deg"]))
    ax_left.grid(True, alpha=0.3)
    ax_left.legend(frameon=True)

    ax_right.set_title(right_title)
    ax_right.set_xlabel(axis_label)
    ax_right.set_ylabel("Normalized distribution")
    ax_right.set_xlim(float(var_cfg["min_deg"]), float(var_cfg["max_deg"]))
    ax_right.grid(True, alpha=0.3)
    ax_right.legend(frameon=True)

    out_path = os.path.join(output_dir, f"{var_key}_dependence_histograms_1x2.png")
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    saved_paths.append(out_path)

    if var_key == "p2_theta":
        fig2, axes2 = plt.subplots(1, 2, figsize=(16, 6), constrained_layout=True)
        ax2_left = axes2[0]
        ax2_right = axes2[1]

        for period in PERIOD_ORDER:
            if period not in angle_histograms_by_period:
                continue
            #endif

            entry = angle_histograms_by_period[period]
            edges = np.asarray(entry["edges"], dtype=float)
            counts_phi_gt_60 = np.asarray(entry.get("counts_phi_gt_60", np.zeros(len(edges) - 1)), dtype=float)
            charge_nC = float(entry["charge_nC"])

            if charge_nC <= 0.0:
                continue
            #endif

            counts_per_nC_phi_gt_60 = counts_phi_gt_60 / charge_nC

            ax2_left.stairs(
                counts_per_nC_phi_gt_60,
                edges,
                label=period,
                color=period_color[period],
                linewidth=1.5,
            )

            bin_widths = np.diff(edges)
            integral_phi_gt_60 = float(np.sum(counts_per_nC_phi_gt_60 * bin_widths))

            if integral_phi_gt_60 > 0.0:
                normalized_phi_gt_60 = counts_per_nC_phi_gt_60 / integral_phi_gt_60
                ax2_right.stairs(
                    normalized_phi_gt_60,
                    edges,
                    label=period,
                    color=period_color[period],
                    linewidth=1.5,
                )
            #endif
        #endfor

        left_title2 = (
            f"{var_cfg['display_name']} histograms: counts / accumulated charge"
            f" (p2_phi > {PHOTON_HIST_PHI_MIN_DEG:.0f} deg)"
        )
        right_title2 = (
            f"{var_cfg['display_name']} histograms: normalized to integral"
            f" (p2_phi > {PHOTON_HIST_PHI_MIN_DEG:.0f} deg)"
        )

        if title_suffix != "":
            left_title2 += title_suffix
            right_title2 += title_suffix
        #endif

        ax2_left.set_title(left_title2)
        ax2_left.set_xlabel(axis_label)
        ax2_left.set_ylabel("Counts / accumulated charge (1/nC)")
        ax2_left.set_xlim(float(var_cfg["min_deg"]), float(var_cfg["max_deg"]))
        ax2_left.grid(True, alpha=0.3)
        ax2_left.legend(frameon=True)

        ax2_right.set_title(right_title2)
        ax2_right.set_xlabel(axis_label)
        ax2_right.set_ylabel("Normalized distribution")
        ax2_right.set_xlim(float(var_cfg["min_deg"]), float(var_cfg["max_deg"]))
        ax2_right.grid(True, alpha=0.3)
        ax2_right.legend(frameon=True)

        out_path2 = os.path.join(
            output_dir,
            f"{var_key}_dependence_histograms_1x2_p2phi_gt_{int(PHOTON_HIST_PHI_MIN_DEG)}deg.png"
        )
        fig2.savefig(out_path2, dpi=200)
        plt.close(fig2)
        saved_paths.append(out_path2)
    #endif

    return saved_paths
#enddef


def make_mc_rec_angle_hist_overlay_plot(angle_histograms_by_period, output_dir, period_color, var_cfg, title_suffix=""):

    os.makedirs(output_dir, exist_ok=True)

    var_key = var_cfg["key"]
    axis_label = get_angle_axis_label(var_key)
    saved_paths = []

    fig, axes = plt.subplots(1, 2, figsize=(16, 6), constrained_layout=True)
    ax_left = axes[0]
    ax_right = axes[1]

    for period in PERIOD_ORDER:
        if period not in angle_histograms_by_period:
            continue
        #endif

        entry = angle_histograms_by_period[period]
        edges = np.asarray(entry["edges"], dtype=float)
        counts = np.asarray(entry["counts"], dtype=float)

        ax_left.stairs(
            counts,
            edges,
            label=period,
            color=period_color[period],
            linewidth=1.5,
        )

        bin_widths = np.diff(edges)
        integral = float(np.sum(counts * bin_widths))

        if integral > 0.0:
            normalized = counts / integral
            ax_right.stairs(
                normalized,
                edges,
                label=period,
                color=period_color[period],
                linewidth=1.5,
            )
        #endif
    #endfor

    left_title = f"{var_cfg['display_name']} histograms: reconstructed MC counts"
    right_title = f"{var_cfg['display_name']} histograms: reconstructed MC normalized to integral"

    if title_suffix != "":
        left_title += title_suffix
        right_title += title_suffix
    #endif

    ax_left.set_title(left_title)
    ax_left.set_xlabel(axis_label)
    ax_left.set_ylabel("Reconstructed MC counts")
    ax_left.set_xlim(float(var_cfg["min_deg"]), float(var_cfg["max_deg"]))
    ax_left.grid(True, alpha=0.3)
    ax_left.legend(frameon=True)

    ax_right.set_title(right_title)
    ax_right.set_xlabel(axis_label)
    ax_right.set_ylabel("Normalized distribution")
    ax_right.set_xlim(float(var_cfg["min_deg"]), float(var_cfg["max_deg"]))
    ax_right.grid(True, alpha=0.3)
    ax_right.legend(frameon=True)

    out_path = os.path.join(output_dir, f"{var_key}_dependence_histograms_1x2_rec_mc.png")
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    saved_paths.append(out_path)

    if var_key == "p2_theta":
        fig2, axes2 = plt.subplots(1, 2, figsize=(16, 6), constrained_layout=True)
        ax2_left = axes2[0]
        ax2_right = axes2[1]

        for period in PERIOD_ORDER:
            if period not in angle_histograms_by_period:
                continue
            #endif

            entry = angle_histograms_by_period[period]
            edges = np.asarray(entry["edges"], dtype=float)
            counts_phi_gt_60 = np.asarray(entry.get("counts_phi_gt_60", np.zeros(len(edges) - 1)), dtype=float)

            ax2_left.stairs(
                counts_phi_gt_60,
                edges,
                label=period,
                color=period_color[period],
                linewidth=1.5,
            )

            bin_widths = np.diff(edges)
            integral_phi_gt_60 = float(np.sum(counts_phi_gt_60 * bin_widths))

            if integral_phi_gt_60 > 0.0:
                normalized_phi_gt_60 = counts_phi_gt_60 / integral_phi_gt_60
                ax2_right.stairs(
                    normalized_phi_gt_60,
                    edges,
                    label=period,
                    color=period_color[period],
                    linewidth=1.5,
                )
            #endif
        #endfor

        left_title2 = (
            f"{var_cfg['display_name']} histograms: reconstructed MC counts"
            f" (p2_phi > {PHOTON_HIST_PHI_MIN_DEG:.0f} deg)"
        )
        right_title2 = (
            f"{var_cfg['display_name']} histograms: reconstructed MC normalized to integral"
            f" (p2_phi > {PHOTON_HIST_PHI_MIN_DEG:.0f} deg)"
        )

        if title_suffix != "":
            left_title2 += title_suffix
            right_title2 += title_suffix
        #endif

        ax2_left.set_title(left_title2)
        ax2_left.set_xlabel(axis_label)
        ax2_left.set_ylabel("Reconstructed MC counts")
        ax2_left.set_xlim(float(var_cfg["min_deg"]), float(var_cfg["max_deg"]))
        ax2_left.grid(True, alpha=0.3)
        ax2_left.legend(frameon=True)

        ax2_right.set_title(right_title2)
        ax2_right.set_xlabel(axis_label)
        ax2_right.set_ylabel("Normalized distribution")
        ax2_right.set_xlim(float(var_cfg["min_deg"]), float(var_cfg["max_deg"]))
        ax2_right.grid(True, alpha=0.3)
        ax2_right.legend(frameon=True)

        out_path2 = os.path.join(
            output_dir,
            f"{var_key}_dependence_histograms_1x2_rec_mc_p2phi_gt_{int(PHOTON_HIST_PHI_MIN_DEG)}deg.png"
        )
        fig2.savefig(out_path2, dpi=200)
        plt.close(fig2)
        saved_paths.append(out_path2)
    #endif

    return saved_paths
#enddef


def make_photon_sector_normalized_panel(angle_histograms_by_period, output_dir, period_color, title_suffix=""):

    os.makedirs(output_dir, exist_ok=True)

    fig = plt.figure(figsize=(18, 10), constrained_layout=True)
    gs = GridSpec(2, 3, figure=fig)
    axes = [fig.add_subplot(gs[i // 3, i % 3]) for i in range(6)]

    for i, sector in enumerate(CLAS12_SECTOR_ORDER):
        ax = axes[i]

        for period in PERIOD_ORDER:
            if period not in angle_histograms_by_period:
                continue
            #endif

            entry = angle_histograms_by_period[period]
            edges = np.asarray(entry["edges"], dtype=float)
            charge_nC = float(entry["charge_nC"])

            if charge_nC <= 0.0:
                continue
            #endif

            sector_counts_map = entry.get("sector_counts", {})
            if sector not in sector_counts_map:
                continue
            #endif

            counts_sector = np.asarray(sector_counts_map[sector], dtype=float)
            counts_per_nC_sector = counts_sector / charge_nC

            bin_widths = np.diff(edges)
            integral_sector = float(np.sum(counts_per_nC_sector * bin_widths))

            if integral_sector <= 0.0:
                continue
            #endif

            normalized_sector = counts_per_nC_sector / integral_sector

            normalized_sector_plot = normalized_sector.copy()
            low_edge_per_bin = edges[:-1]
            normalized_sector_plot[low_edge_per_bin < 5.0] = 0.0

            ax.stairs(
                normalized_sector_plot,
                edges,
                label=period,
                color=period_color[period],
                linewidth=1.5,
            )
        #endfor

        title = f"Photon theta, sector {sector}"
        if title_suffix != "":
            title += title_suffix
        #endif

        ax.set_title(title)
        ax.set_xlabel(get_angle_axis_label("p2_theta"))
        ax.set_ylabel("Normalized distribution")
        ax.set_xlim(5.0, 35.0)
        ax.grid(True, alpha=0.3)
        ax.legend(frameon=True, fontsize=9)
    #endfor

    out_path = os.path.join(output_dir, "p2_theta_dependence_histograms_sector_panel_normalized.png")
    fig.savefig(out_path, dpi=200)
    plt.close(fig)

    return out_path
#enddef


def make_mc_rec_photon_sector_normalized_panel(angle_histograms_by_period, output_dir, period_color, title_suffix=""):

    os.makedirs(output_dir, exist_ok=True)

    fig = plt.figure(figsize=(18, 10), constrained_layout=True)
    gs = GridSpec(2, 3, figure=fig)
    axes = [fig.add_subplot(gs[i // 3, i % 3]) for i in range(6)]

    for i, sector in enumerate(CLAS12_SECTOR_ORDER):
        ax = axes[i]

        for period in PERIOD_ORDER:
            if period not in angle_histograms_by_period:
                continue
            #endif

            entry = angle_histograms_by_period[period]
            edges = np.asarray(entry["edges"], dtype=float)

            sector_counts_map = entry.get("sector_counts", {})
            if sector not in sector_counts_map:
                continue
            #endif

            counts_sector = np.asarray(sector_counts_map[sector], dtype=float)

            bin_widths = np.diff(edges)
            integral_sector = float(np.sum(counts_sector * bin_widths))

            if integral_sector <= 0.0:
                continue
            #endif

            normalized_sector = counts_sector / integral_sector

            normalized_sector_plot = normalized_sector.copy()
            low_edge_per_bin = edges[:-1]
            normalized_sector_plot[low_edge_per_bin < 5.0] = 0.0

            ax.stairs(
                normalized_sector_plot,
                edges,
                label=period,
                color=period_color[period],
                linewidth=1.5,
            )
        #endfor

        title = f"Photon theta, sector {sector} reconstructed MC"
        if title_suffix != "":
            title += title_suffix
        #endif

        ax.set_title(title)
        ax.set_xlabel(get_angle_axis_label("p2_theta"))
        ax.set_ylabel("Normalized distribution")
        ax.set_xlim(5.0, 35.0)
        ax.grid(True, alpha=0.3)
        ax.legend(frameon=True, fontsize=9)
    #endfor

    out_path = os.path.join(output_dir, "p2_theta_dependence_histograms_sector_panel_normalized_rec_mc.png")
    fig.savefig(out_path, dpi=200)
    plt.close(fig)

    return out_path
#enddef


def make_angle_summary_plot(angle_summary_rows, output_dir, period_color, var_cfg):

    fig = plt.figure(figsize=(18, 10), constrained_layout=True)
    gs = GridSpec(2, 3, figure=fig)
    axes = [fig.add_subplot(gs[i // 3, i % 3]) for i in range(6)]

    all_xmeans = []
    for row in angle_summary_rows:
        if row["angle_variable"] != var_cfg["key"]:
            continue
        #endif

        xmean = float(row["angle_x_mean_deg"])
        val = float(row["divisor_to_divide_by_fraction"])

        if not np.isfinite(xmean):
            continue
        #endif
        if not np.isfinite(val):
            continue
        #endif

        all_xmeans.append(xmean)
    #endfor

    if len(all_xmeans) > 0:
        xmin_global = min(all_xmeans)
        xmax_global = max(all_xmeans)

        xmin_global = 5.0 * math.floor(xmin_global / 5.0)
        xmax_global = 5.0 * math.ceil(xmax_global / 5.0)

        if xmin_global == xmax_global:
            xmin_global -= 5.0
            xmax_global += 5.0
        #endif
    else:
        xmin_global = 0.0
        xmax_global = 5.0
    #endif

    for i, period in enumerate(PERIOD_ORDER):
        ax = axes[i]
        c = period_color[period]

        rows = [r for r in angle_summary_rows if r["period"] == period and r["angle_variable"] == var_cfg["key"]]
        rows = sorted(rows, key=lambda r: r["angle_x_mean_deg"])

        x = []
        y = []
        yerr = []

        for row in rows:
            xmean = float(row["angle_x_mean_deg"])
            val = float(row["divisor_to_divide_by_fraction"])
            err = float(row["divisor_to_divide_by_stat_err_fraction"])

            if not np.isfinite(val):
                continue
            #endif
            if not np.isfinite(xmean):
                continue
            #endif

            x.append(xmean)
            y.append(val)
            yerr.append(err if np.isfinite(err) else 0.0)
        #endfor

        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        yerr = np.asarray(yerr, dtype=float)

        if len(x) > 0:
            ax.errorbar(x, y, yerr=yerr, fmt="o", capsize=3, color=c)
        #endif

        ax.set_title(period)
        ax.set_xlabel(get_angle_axis_label(var_cfg["key"]))
        ax.set_ylabel("Divisor for cross_sections.cpp")
        ax.set_xlim(xmin_global, xmax_global)
        ax.set_ylim(0.4, 1.3)
        ax.grid(True, alpha=0.3)
    #endfor

    axes[5].axis("off")

    out_path = os.path.join(output_dir, f"{var_cfg['key']}_dependence_normalization_divisor.png")
    fig.savefig(out_path, dpi=200)
    plt.close(fig)

    return out_path
#enddef


def print_summary_table(title, rows):

    print("")
    print(title)
    header = (
        f"{'Period':12s}  "
        f"{'weighted_data_rel':>18s}  "
        f"{'mc_rel@ref':>12s}  "
        f"{'data/mc@ref':>12s}  "
        f"{'divide_by':>12s}  "
        f"{'applied_scale':>13s}"
    )
    print(header)
    print("-" * len(header))

    for row in rows:
        print(
            f"{row['period']:12s}  "
            f"{float(row['weighted_data_rel_fraction']):18.6f}  "
            f"{float(row['mc_rel_at_ref_fraction']):12.6f}  "
            f"{float(row['data_over_mc_at_ref_fraction']):12.6f}  "
            f"{float(row['divisor_to_divide_by_fraction']):12.6f}  "
            f"{float(row['applied_scale_fraction']):13.6f}"
        )
    #endfor
#enddef


def run_selection_analysis(run_charge_map, period_color):

    OUTPUT_DIR, INTEGRATED_OUTPUT_DIR, ANGLE_OUTPUT_DIR = get_output_paths()

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(INTEGRATED_OUTPUT_DIR, exist_ok=True)
    os.makedirs(ANGLE_OUTPUT_DIR, exist_ok=True)

    print("")
    print("------------------------------------------------------------")
    print("Selection: overall")
    print(f"Output dir: {OUTPUT_DIR}")
    print("------------------------------------------------------------")

    all_run_rows = []
    all_current_rows = []
    all_skipped_run_rows = []

    angle_current_rows_by_variable = {}
    fine_angle_histograms_by_variable = {}

    for var_cfg in ANGLE_DEPENDENCE_CONFIG:
        angle_current_rows_by_variable[var_cfg["key"]] = {}
        fine_angle_histograms_by_variable[var_cfg["key"]] = {}
    #endfor

    print("")
    print("Processing DATA ROOT files and aggregating by current...")

    data_tasks = []
    for period_display_name, period_internal_name, root_path in DATA_PERIOD_FILES:
        data_tasks.append((period_display_name, period_internal_name, root_path, run_charge_map))
    #endfor

    n_workers = min(MAX_WORKERS, max(1, len(data_tasks)))

    with cf.ProcessPoolExecutor(max_workers=n_workers) as ex:
        futures = [ex.submit(process_data_period_worker, task) for task in data_tasks]

        results_by_period = {}
        for fut in cf.as_completed(futures):
            result = fut.result()
            results_by_period[result["period"]] = result
        #endfor
    #endwith

    for period_display_name, period_internal_name, root_path in DATA_PERIOD_FILES:
        result = results_by_period[period_display_name]

        print("")
        print("=" * 90)
        print(f"DATA period: {period_display_name}")
        print(f"Internal label: {period_internal_name}")
        print(f"ROOT file: {root_path}")
        print("=" * 90)

        run_rows = result["run_rows"]
        current_rows = result["current_rows"]
        angle_current_rows_this_period = result["angle_current_rows"]
        skipped_nonpositive_charge_rows = result["skipped_nonpositive_charge_rows"]
        fine_angle_histograms_this_period = result["fine_angle_histograms"]
        total_valid_charge_nC = float(result["total_valid_charge_nC"])

        all_run_rows.extend(run_rows)
        all_current_rows.extend(current_rows)
        all_skipped_run_rows.extend(skipped_nonpositive_charge_rows)

        for var_cfg in ANGLE_DEPENDENCE_CONFIG:
            var_key = var_cfg["key"]
            hist_entry = {
                "edges": np.asarray(fine_angle_histograms_this_period[var_key]["edges"], dtype=float),
                "counts": np.asarray(fine_angle_histograms_this_period[var_key]["counts"], dtype=float),
                "charge_nC": float(total_valid_charge_nC),
            }

            if var_key == "p2_theta":
                if "counts_phi_gt_60" in fine_angle_histograms_this_period[var_key]:
                    hist_entry["counts_phi_gt_60"] = np.asarray(
                        fine_angle_histograms_this_period[var_key]["counts_phi_gt_60"],
                        dtype=float,
                    )
                #endif

                if "sector_counts" in fine_angle_histograms_this_period[var_key]:
                    hist_entry["sector_counts"] = {
                        sector: np.asarray(
                            fine_angle_histograms_this_period[var_key]["sector_counts"][sector],
                            dtype=float,
                        )
                        for sector in CLAS12_SECTOR_ORDER
                    }
                #endif
            #endif

            angle_current_rows_by_variable[var_key][period_display_name] = angle_current_rows_this_period[var_key]
            fine_angle_histograms_by_variable[var_key][period_display_name] = hist_entry
        #endfor

        print(f"Run-level rows kept: {len(run_rows)}")
        print(f"Integrated current groups: {len(current_rows)}")
        print(f"Skipped non-positive-charge runs: {len(skipped_nonpositive_charge_rows)}")
        print(f"Total valid charge used for fine histograms: {total_valid_charge_nC:.5f} nC")

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

    integrated_mc_rows, angle_mc_rows_by_variable, rec_fine_angle_histograms_by_period = build_mc_aggregation(
        mc_dir=MC_DIR,
    )

    for row in integrated_mc_rows:
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
    skipped_run_table_csv = os.path.join(OUTPUT_DIR, "dvcs_current_dependence_skipped_runs.csv")

    write_run_table_csv(run_table_csv, all_run_rows)
    write_current_table_csv(current_table_csv, all_current_rows)
    write_mc_table_csv(mc_table_csv, integrated_mc_rows)
    write_skipped_run_table_csv(skipped_run_table_csv, all_skipped_run_rows)

    print("")
    print(f"[saved] {run_table_csv}")
    print(f"[saved] {current_table_csv}")
    print(f"[saved] {mc_table_csv}")
    print(f"[saved] {skipped_run_table_csv}")

    integrated_data_fit_results = make_fit_result_map_for_data(all_current_rows)
    integrated_mc_fit_results = make_fit_result_map_for_mc(integrated_mc_rows)

    print("")
    print("=== DATA fits: y = m*x + b for counts_per_nC vs current ===")
    for period in PERIOD_ORDER:
        fr = integrated_data_fit_results[period]
        print(
            f"{period}: m = {fr['m']:.10f} +/- {fr['sm']:.10f}, "
            f"b = {fr['b']:.10f} +/- {fr['sb']:.10f}, "
            f"chi2/ndof = {fr['chi2']:.2f}/{fr['ndof']}"
        )
    #endfor

    print("")
    print("=== MC fits: y = m*x + b for efficiency vs current ===")
    for period in PERIOD_ORDER:
        fr = integrated_mc_fit_results[period]
        print(
            f"{period}: m = {fr['m']:.10f} +/- {fr['sm']:.10f}, "
            f"b = {fr['b']:.10f} +/- {fr['sb']:.10f}, "
            f"chi2/ndof = {fr['chi2']:.2f}/{fr['ndof']}"
        )
    #endfor

    residual_rows = build_residual_correction_rows(
        all_current_rows=all_current_rows,
        data_fit_results=integrated_data_fit_results,
        mc_fit_results=integrated_mc_fit_results,
    )

    residual_table_csv = os.path.join(OUTPUT_DIR, "dvcs_current_dependence_residual_table.csv")
    write_residual_table_csv(residual_table_csv, residual_rows)
    print("")
    print(f"[saved] {residual_table_csv}")

    period_summary_rows = build_period_summary_rows(
        current_rows=all_current_rows,
        data_fit_results=integrated_data_fit_results,
        mc_fit_results=integrated_mc_fit_results,
    )

    period_summary_csv = os.path.join(OUTPUT_DIR, "dvcs_current_dependence_period_summary.csv")
    write_period_summary_csv(period_summary_csv, period_summary_rows)
    print(f"[saved] {period_summary_csv}")

    out_a, out_b, out_c, out_d = plot_four_panel_set(
        output_dir=INTEGRATED_OUTPUT_DIR,
        tag="integrated",
        data_current_rows=all_current_rows,
        mc_rows=integrated_mc_rows,
        data_fit_results=integrated_data_fit_results,
        mc_fit_results=integrated_mc_fit_results,
        period_color=period_color,
        title_suffix="",
    )

    print("")
    print(f"[saved] {out_a}")
    print(f"[saved] {out_b}")
    print(f"[saved] {out_c}")
    print(f"[saved] {out_d}")

    for var_cfg in ANGLE_DEPENDENCE_CONFIG:
        var_key = var_cfg["key"]
        var_output_dir = os.path.join(ANGLE_OUTPUT_DIR, var_key)
        os.makedirs(var_output_dir, exist_ok=True)

        angle_summary_rows = build_angle_summary_rows(
            angle_rows_by_period=angle_current_rows_by_variable[var_key],
            mc_rows_by_variable=angle_mc_rows_by_variable[var_key],
            var_cfg=var_cfg,
        )

        angle_summary_csv = os.path.join(var_output_dir, f"dvcs_current_dependence_{var_key}_summary.csv")
        write_angle_summary_csv(angle_summary_csv, angle_summary_rows)
        print(f"[saved] {angle_summary_csv}")

        for period_name in PERIOD_ORDER:
            bins = angle_current_rows_by_variable[var_key][period_name]["bins"]

            for ibin, (low, high) in enumerate(bins):
                subdir = os.path.join(var_output_dir, period_name.replace(" ", "_"), angle_bin_to_dirname(var_key, low, high))
                os.makedirs(subdir, exist_ok=True)

                current_rows_bin = angle_current_rows_by_variable[var_key][period_name]["rows"][ibin]
                mc_rows_bin = angle_mc_rows_by_variable[var_key]["rows"][period_name][ibin]

                data_fit_bin = make_fit_result_map_for_data(current_rows_bin)
                mc_fit_bin = make_fit_result_map_for_mc(mc_rows_bin)

                label = angle_bin_to_label(low, high, is_last=(ibin == len(bins) - 1))
                angle_title_suffix = f"  {var_cfg['display_name']} {label}"

                out_a_bin, out_b_bin, out_c_bin, out_d_bin = plot_four_panel_set(
                    output_dir=subdir,
                    tag=angle_bin_to_dirname(var_key, low, high),
                    data_current_rows=current_rows_bin,
                    mc_rows=mc_rows_bin,
                    data_fit_results=data_fit_bin,
                    mc_fit_results=mc_fit_bin,
                    period_color=period_color,
                    title_suffix=angle_title_suffix,
                )

                print(f"[saved] {out_a_bin}")
                print(f"[saved] {out_b_bin}")
                print(f"[saved] {out_c_bin}")
                print(f"[saved] {out_d_bin}")
            #endfor
        #endfor

        angle_hist_plots = make_angle_hist_overlay_plot(
            angle_histograms_by_period=fine_angle_histograms_by_variable[var_key],
            output_dir=var_output_dir,
            period_color=period_color,
            var_cfg=var_cfg,
            title_suffix="",
        )
        for path in angle_hist_plots:
            print(f"[saved] {path}")
        #endfor

        mc_rec_angle_hist_plots = make_mc_rec_angle_hist_overlay_plot(
            angle_histograms_by_period=rec_fine_angle_histograms_by_period,
            output_dir=var_output_dir,
            period_color=period_color,
            var_cfg=var_cfg,
            title_suffix="",
        )
        for path in mc_rec_angle_hist_plots:
            print(f"[saved] {path}")
        #endfor

        if var_key == "p2_theta":
            photon_sector_panel_path = make_photon_sector_normalized_panel(
                angle_histograms_by_period=fine_angle_histograms_by_variable[var_key],
                output_dir=var_output_dir,
                period_color=period_color,
                title_suffix="",
            )
            print(f"[saved] {photon_sector_panel_path}")

            mc_rec_photon_sector_panel_path = make_mc_rec_photon_sector_normalized_panel(
                angle_histograms_by_period=rec_fine_angle_histograms_by_period,
                output_dir=var_output_dir,
                period_color=period_color,
                title_suffix="",
            )
            print(f"[saved] {mc_rec_photon_sector_panel_path}")
        #endif

        angle_summary_plot = make_angle_summary_plot(
            angle_summary_rows=angle_summary_rows,
            output_dir=var_output_dir,
            period_color=period_color,
            var_cfg=var_cfg,
        )
        print(f"[saved] {angle_summary_plot}")
    #endfor

    if len(all_skipped_run_rows) > 0:
        print("")
        print("Skipped runs with non-positive accumulated charge:")
        print(f"  total skipped = {len(all_skipped_run_rows)}")
        print(f"  details saved to: {skipped_run_table_csv}")
    #endif

    return period_summary_rows
#enddef


def main():

    default_colors = plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
    if len(default_colors) < len(PERIOD_ORDER):
        raise RuntimeError("Matplotlib default color cycle is shorter than number of periods.")
    #endif

    period_color = {}
    for i, period in enumerate(PERIOD_ORDER):
        period_color[period] = default_colors[i]
    #endfor

    print("")
    print("============================================================")
    print("Running channel: epgamma")
    print(f"MC directory:    {MC_DIR}")
    print(f"Max workers:     {MAX_WORKERS}")
    print(f"e-gamma cone cut: {MIN_E_GAMMA_CONE_ANGLE_DEG:.1f} deg")
    print(f"extra photon histogram phi cut: p2_phi > {PHOTON_HIST_PHI_MIN_DEG:.1f} deg")
    print("============================================================")

    print("")
    print("Reading charge map...")
    run_charge_map = read_charge_map(CSV_FILE)
    print(f"Loaded charge entries for {len(run_charge_map)} runs from:")
    print(f"  {CSV_FILE}")

    summary_rows = run_selection_analysis(
        run_charge_map=run_charge_map,
        period_color=period_color,
    )

    title = "Representative integrated period-level normalization summary: overall channel epgamma"
    print_summary_table(title, summary_rows)

    print("")
#enddef


if __name__ == "__main__":
    main()
#endif