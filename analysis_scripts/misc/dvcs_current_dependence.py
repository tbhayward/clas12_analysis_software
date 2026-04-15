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
# Default channel
# ---------------
# If no channel is specified, this runs for:
#
#   epgamma
#
# Additional runtime channel option
# ---------------------------------
# You can also switch channels at runtime, for example:
#
#   python dvcs_current_dependence.py --channel eX
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

import os
import math
import csv
import argparse
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

MC_TREE_NAME = "PhysicsEvents"

ANGLE_DEPENDENCE_CONFIG = [
    {
        "key": "e_theta",
        "display_name": "electron theta",
        "branch": "e_theta",
        "min_deg": 8.0,
        "max_deg": 35.0,
        "n_bins": 7,
    },
    {
        "key": "p1_theta",
        "display_name": "proton theta",
        "branch": "p1_theta",
        "min_deg": 8.0,
        "max_deg": 70.0,
        "n_bins": 7,
    },
    {
        "key": "p2_theta",
        "display_name": "photon theta",
        "branch": "p2_theta",
        "min_deg": 0.0,
        "max_deg": 35.0,
        "n_bins": 7,
    },
]

CHANNEL_CONFIG = {
    "epgamma": {
        "display_name": "epgamma",
        "mc_dir": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen",
        "mc_channel_tag": None,
        "supports_topology": True,
        "data_period_files": [
            ("Sp18 Inb", "rga_sp18_inb", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp18_inb_epgamma.root"),
            ("Sp18 Out", "rga_sp18_out", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp18_out_epgamma.root"),
            ("Fa18 Inb", "rga_fa18_inb", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_inb_epgamma.root"),
            ("Fa18 Out", "rga_fa18_out", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_out_epgamma.root"),
            ("Sp19 Inb", "rga_sp19_inb", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp19_inb_epgamma.root"),
        ],
    },
    "eX": {
        "display_name": "eX",
        "mc_dir": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/eX",
        "mc_channel_tag": "eX",
        "supports_topology": False,
        "data_period_files": [
            ("Sp18 Inb", "rga_sp18_inb", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/eX/rga_sp18_inb_eX.root"),
            ("Sp18 Out", "rga_sp18_out", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/eX/rga_sp18_out_eX.root"),
            ("Fa18 Inb", "rga_fa18_inb", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/eX/rga_fa18_inb_eX.root"),
            ("Fa18 Out", "rga_fa18_out", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/eX/rga_fa18_out_eX.root"),
            ("Sp19 Inb", "rga_sp19_inb", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/eX/rga_sp19_inb_eX.root"),
        ],
    },
}

TOPOLOGY_CONFIG = [
    {
        "key": "overall",
        "display_name": "overall",
        "dir_name": None,
        "cuts": None,
    },
    {
        "key": "CD_FT",
        "display_name": "CD,FT",
        "dir_name": "CD_FT",
        "cuts": {
            "p1_theta_min_deg": 30.0,
            "p1_theta_max_deg": 70.0,
            "p2_theta_min_deg": 0.0,
            "p2_theta_max_deg": 6.0,
        },
    },
    {
        "key": "CD_FD",
        "display_name": "CD,FD",
        "dir_name": "CD_FD",
        "cuts": {
            "p1_theta_min_deg": 30.0,
            "p1_theta_max_deg": 70.0,
            "p2_theta_min_deg": 6.0,
            "p2_theta_max_deg": 30.0,
        },
    },
    {
        "key": "FD_FD",
        "display_name": "FD,FD",
        "dir_name": "FD_FD",
        "cuts": {
            "p1_theta_min_deg": 6.0,
            "p1_theta_max_deg": 30.0,
            "p2_theta_min_deg": 6.0,
            "p2_theta_max_deg": 30.0,
        },
    },
]

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
# Bin cache
# -----------------------------------------------------------------------------

SUMMARY_BINS_CACHE = {}
for var_cfg in ANGLE_DEPENDENCE_CONFIG:
    SUMMARY_BINS_CACHE[var_cfg["key"]] = {}
    for period_name in PERIOD_ORDER:
        SUMMARY_BINS_CACHE[var_cfg["key"]][period_name] = get_summary_bins_for_period(var_cfg, period_name)
    #endfor
#endfor

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

def get_channel_config(channel):

    if channel not in CHANNEL_CONFIG:
        raise RuntimeError(
            f"Unknown channel '{channel}'. Available channels: {sorted(CHANNEL_CONFIG.keys())}"
        )
    #endif

    return CHANNEL_CONFIG[channel]
#enddef


def get_summary_bins(var_key, period_name):
    return SUMMARY_BINS_CACHE[var_key][period_name]
#enddef


def read_data_run_counts_and_angle_counts(root_path, topology_cuts=None):

    if not os.path.exists(root_path):
        raise RuntimeError(f"ROOT file not found: {root_path}")
    #endif

    root_file = uproot.open(root_path)
    if "PhysicsEvents" not in root_file:
        raise RuntimeError(f"'PhysicsEvents' tree not found in {root_path}")
    #endif

    tree = root_file["PhysicsEvents"]

    needed = {"runnum", "e_theta", "p1_theta", "p2_theta"}
    missing = [name for name in needed if name not in tree.keys()]
    if len(missing) > 0:
        raise RuntimeError(f"Missing required branches in {root_path}: {missing}")
    #endif

    total_run_counts = defaultdict(int)
    angle_run_counts = {}
    angle_theta_sum = {}
    angle_theta_n = {}

    for var_cfg in ANGLE_DEPENDENCE_CONFIG:
        var_key = var_cfg["key"]
        angle_run_counts[var_key] = {}
        angle_theta_sum[var_key] = {}
        angle_theta_n[var_key] = {}
        for period_name in PERIOD_ORDER:
            bins = get_summary_bins(var_key, period_name)
            angle_run_counts[var_key][period_name] = [defaultdict(int) for _ in bins]
            angle_theta_sum[var_key][period_name] = [0.0 for _ in bins]
            angle_theta_n[var_key][period_name] = [0 for _ in bins]
        #endfor
    #endfor

    iterate_branches = ["runnum", "e_theta", "p1_theta", "p2_theta"]

    for arrays in tree.iterate(iterate_branches, library="np", step_size=ITERATE_STEP_SIZE):
        event_mask = apply_topology_mask(arrays, topology_cuts)

        if not np.any(event_mask):
            continue
        #endif

        runnum = arrays["runnum"][event_mask]

        unique_runs, counts = np.unique(runnum, return_counts=True)
        for r, c in zip(unique_runs, counts):
            total_run_counts[int(r)] += int(c)
        #endfor

        theta_deg_map = {
            "e_theta": np.degrees(arrays["e_theta"][event_mask]),
            "p1_theta": np.degrees(arrays["p1_theta"][event_mask]),
            "p2_theta": np.degrees(arrays["p2_theta"][event_mask]),
        }

        for var_cfg in ANGLE_DEPENDENCE_CONFIG:
            var_key = var_cfg["key"]
            theta_deg = theta_deg_map[var_key]

            for period_name in PERIOD_ORDER:
                bins = get_summary_bins(var_key, period_name)
                masks = angle_bin_masks(theta_deg, bins)

                counts_store = angle_run_counts[var_key][period_name]
                sum_store = angle_theta_sum[var_key][period_name]
                n_store = angle_theta_n[var_key][period_name]

                for ibin, mask in enumerate(masks):
                    if not np.any(mask):
                        continue
                    #endif

                    selected_runs = runnum[mask]
                    unique_sel, counts_sel = np.unique(selected_runs, return_counts=True)

                    for r, c in zip(unique_sel, counts_sel):
                        counts_store[ibin][int(r)] += int(c)
                    #endfor

                    sum_store[ibin] += float(np.sum(theta_deg[mask]))
                    n_store[ibin] += int(np.count_nonzero(mask))
                #endfor
            #endfor
        #endfor
    #endfor

    angle_run_counts_out = {}
    angle_x_means_out = {}

    for var_cfg in ANGLE_DEPENDENCE_CONFIG:
        var_key = var_cfg["key"]
        angle_run_counts_out[var_key] = {}
        angle_x_means_out[var_key] = {}

        for period_name in PERIOD_ORDER:
            angle_run_counts_out[var_key][period_name] = [dict(d) for d in angle_run_counts[var_key][period_name]]
            means = []
            for s, n in zip(angle_theta_sum[var_key][period_name], angle_theta_n[var_key][period_name]):
                if n > 0:
                    means.append(float(s) / float(n))
                else:
                    means.append(float("nan"))
                #endif
            #endfor
            angle_x_means_out[var_key][period_name] = means
        #endfor
    #endfor

    return dict(total_run_counts), angle_run_counts_out, angle_x_means_out
#enddef


def count_mc_total_and_angle_entries(root_path, tree_name, topology_cuts=None):

    if not os.path.exists(root_path):
        raise RuntimeError(f"ROOT file not found: {root_path}")
    #endif

    root_file = uproot.open(root_path)
    if tree_name not in root_file:
        raise RuntimeError(f"'{tree_name}' tree not found in {root_path}")
    #endif

    tree = root_file[tree_name]

    needed = {"e_theta", "p1_theta", "p2_theta"}
    missing = [name for name in needed if name not in tree.keys()]
    if len(missing) > 0:
        raise RuntimeError(f"Missing required branches in {root_path}:{tree_name}: {missing}")
    #endif

    iterate_branches = ["e_theta", "p1_theta", "p2_theta"]

    total_count = 0
    angle_counts = {}

    for var_cfg in ANGLE_DEPENDENCE_CONFIG:
        var_key = var_cfg["key"]
        angle_counts[var_key] = {}
        for period_name in PERIOD_ORDER:
            bins = get_summary_bins(var_key, period_name)
            angle_counts[var_key][period_name] = [0 for _ in bins]
        #endfor
    #endfor

    for arrays in tree.iterate(iterate_branches, library="np", step_size=ITERATE_STEP_SIZE):
        event_mask = apply_topology_mask(arrays, topology_cuts)

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

            for period_name in PERIOD_ORDER:
                bins = get_summary_bins(var_key, period_name)
                masks = angle_bin_masks(theta_deg, bins)
                store = angle_counts[var_key][period_name]

                for ibin, mask in enumerate(masks):
                    store[ibin] += int(np.count_nonzero(mask))
                #endfor
            #endfor
        #endfor
    #endfor

    return total_count, angle_counts
#enddef


def build_run_metadata(period_display_name, period_internal_name, run_list, run_charge_map):

    missing_charge_runs = []
    unknown_current_runs = []

    run_meta = {}
    current_charge_totals = defaultdict(float)
    skipped_nonpositive_charge_rows = []

    for runnum in sorted(run_list):
        if runnum not in run_charge_map:
            missing_charge_runs.append(runnum)
            continue
        #endif

        charge = float(run_charge_map[runnum])
        if charge <= 0.0:
            ok_cur, current_nA = resolve_current(period_internal_name, runnum)
            current_out = int(current_nA) if ok_cur else None

            skipped_nonpositive_charge_rows.append({
                "period": period_display_name,
                "period_internal": period_internal_name,
                "runnum": int(runnum),
                "current_nA": current_out,
                "charge_nC": float(charge),
                "reason": "nonpositive_charge",
            })
            continue
        #endif

        ok_cur, current_nA = resolve_current(period_internal_name, runnum)
        if not ok_cur:
            unknown_current_runs.append(runnum)
            continue
        #endif

        run_meta[runnum] = {
            "period": period_display_name,
            "period_internal": period_internal_name,
            "runnum": int(runnum),
            "current_nA": int(current_nA),
            "charge_nC": float(charge),
        }

        current_charge_totals[int(current_nA)] += float(charge)
    #endfor

    if len(missing_charge_runs) > 0 or len(unknown_current_runs) > 0:
        msg = []
        msg.append(f"Fatal bookkeeping problem while processing {period_display_name}:")
        if len(missing_charge_runs) > 0:
            msg.append(f"  Runs missing from charge CSV ({len(missing_charge_runs)}): {missing_charge_runs}")
        #endif
        if len(unknown_current_runs) > 0:
            msg.append(f"  Runs with no current mapping ({len(unknown_current_runs)}): {unknown_current_runs}")
        #endif
        raise RuntimeError("\n".join(msg))
    #endif

    return run_meta, dict(current_charge_totals), skipped_nonpositive_charge_rows
#enddef


def build_run_rows_from_counts(run_meta, run_counts):

    run_rows = []

    for runnum in sorted(run_meta.keys()):
        count = int(run_counts.get(runnum, 0))
        charge = float(run_meta[runnum]["charge_nC"])
        current_nA = int(run_meta[runnum]["current_nA"])

        rate = float(count) / charge
        rate_err = math.sqrt(float(count)) / charge if count > 0 else 0.0

        run_rows.append({
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

    return run_rows
#enddef


def build_current_rows_from_counts(period_display_name, period_internal_name, run_meta, current_charge_totals, run_counts):

    per_current_counts = defaultdict(int)
    per_current_n_runs = defaultdict(int)

    for runnum in sorted(run_meta.keys()):
        current_nA = int(run_meta[runnum]["current_nA"])
        count = int(run_counts.get(runnum, 0))
        per_current_counts[current_nA] += int(count)
        per_current_n_runs[current_nA] += 1
    #endfor

    current_rows = []

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

        current_rows.append({
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

    return current_rows
#enddef


def build_data_period_aggregations(period_display_name, period_internal_name, root_path, run_charge_map, topology_cuts=None):

    total_run_counts, angle_run_counts, angle_x_means = read_data_run_counts_and_angle_counts(
        root_path=root_path,
        topology_cuts=topology_cuts,
    )

    run_meta, current_charge_totals, skipped_nonpositive_charge_rows = build_run_metadata(
        period_display_name=period_display_name,
        period_internal_name=period_internal_name,
        run_list=sorted(total_run_counts.keys()),
        run_charge_map=run_charge_map,
    )

    integrated_run_rows = build_run_rows_from_counts(run_meta, total_run_counts)
    integrated_current_rows = build_current_rows_from_counts(
        period_display_name=period_display_name,
        period_internal_name=period_internal_name,
        run_meta=run_meta,
        current_charge_totals=current_charge_totals,
        run_counts=total_run_counts,
    )

    angle_current_rows = {}
    for var_cfg in ANGLE_DEPENDENCE_CONFIG:
        var_key = var_cfg["key"]
        bins = get_summary_bins(var_key, period_display_name)
        angle_current_rows[var_key] = {
            "bins": bins,
            "rows": [],
            "x_means": angle_x_means[var_key][period_display_name],
        }

        for ibin in range(len(bins)):
            rows_bin = build_current_rows_from_counts(
                period_display_name=period_display_name,
                period_internal_name=period_internal_name,
                run_meta=run_meta,
                current_charge_totals=current_charge_totals,
                run_counts=angle_run_counts[var_key][period_display_name][ibin],
            )
            angle_current_rows[var_key]["rows"].append(rows_bin)
        #endfor
    #endfor

    return integrated_run_rows, integrated_current_rows, angle_current_rows, skipped_nonpositive_charge_rows
#enddef


def process_data_period_worker(args):

    period_display_name, period_internal_name, root_path, run_charge_map, topology_cuts = args

    run_rows, current_rows, angle_current_rows, skipped_nonpositive_charge_rows = build_data_period_aggregations(
        period_display_name=period_display_name,
        period_internal_name=period_internal_name,
        root_path=root_path,
        run_charge_map=run_charge_map,
        topology_cuts=topology_cuts,
    )

    return {
        "period": period_display_name,
        "period_internal": period_internal_name,
        "run_rows": run_rows,
        "current_rows": current_rows,
        "angle_current_rows": angle_current_rows,
        "skipped_nonpositive_charge_rows": skipped_nonpositive_charge_rows,
        "root_path": root_path,
    }
#enddef


def process_mc_pair_worker(args):

    period_internal, current_nA, gen_path, rec_path, topology_cuts = args
    period_display_name = PERIOD_DISPLAY_FROM_INTERNAL[period_internal]

    n_gen, angle_counts_gen = count_mc_total_and_angle_entries(
        root_path=gen_path,
        tree_name=MC_TREE_NAME,
        topology_cuts=topology_cuts,
    )
    n_rec, angle_counts_rec = count_mc_total_and_angle_entries(
        root_path=rec_path,
        tree_name=MC_TREE_NAME,
        topology_cuts=topology_cuts,
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
        bins = get_summary_bins(var_key, period_display_name)
        rows = []

        for ibin in range(len(bins)):
            n_gen_bin = int(angle_counts_gen[var_key][period_display_name][ibin])
            n_rec_bin = int(angle_counts_rec[var_key][period_display_name][ibin])

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
    }
#enddef


def build_mc_aggregation(mc_dir, requested_channel_tag=None, topology_cuts=None, skip_temp_heavy_mc=False):

    if not os.path.isdir(mc_dir):
        raise RuntimeError(f"MC directory not found: {mc_dir}")
    #endif

    skip_pairs = set()
    if skip_temp_heavy_mc:
        pass
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

        kind, period_internal, current_nA, beam_energy_token, channel_tag = parse_mc_filename(basename)

        if requested_channel_tag is None:
            if channel_tag is not None:
                continue
            #endif
        else:
            if channel_tag != requested_channel_tag:
                continue
            #endif
        #endif

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

        tasks.append((period_internal, current_nA, gen_path, rec_path, topology_cuts))
    #endfor

    integrated_mc_rows = []
    angle_mc_rows_by_variable = {}

    for var_cfg in ANGLE_DEPENDENCE_CONFIG:
        var_key = var_cfg["key"]
        angle_mc_rows_by_variable[var_key] = {}
        for period_name in PERIOD_ORDER:
            bins = get_summary_bins(var_key, period_name)
            angle_mc_rows_by_variable[var_key][period_name] = {
                "bins": bins,
                "rows": [[] for _ in bins],
            }
        #endfor
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
                store = angle_mc_rows_by_variable[var_key][period_name]["rows"]

                for ibin, row in enumerate(rows):
                    if row is not None:
                        store[ibin].append(row)
                    #endif
                #endfor
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
            bins = angle_mc_rows_by_variable[var_key][period_name]["rows"]
            for ibin in range(len(bins)):
                bins[ibin] = sorted(bins[ibin], key=lambda row: row["current_nA"])
            #endfor
        #endfor
    #endfor

    return integrated_mc_rows, angle_mc_rows_by_variable
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


def build_angle_summary_rows(angle_rows_by_period, mc_rows_by_period, var_cfg):

    summary_rows = []
    var_key = var_cfg["key"]

    for period_name in PERIOD_ORDER:
        bins = angle_rows_by_period[period_name]["bins"]
        current_rows_list = angle_rows_by_period[period_name]["rows"]
        mc_rows_list = mc_rows_by_period[period_name]["rows"]
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
        ax.set_xlabel(f"{var_cfg['display_name']} (deg)")
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


def run_selection_analysis(
    channel,
    channel_cfg,
    topology_cfg,
    run_charge_map,
    period_color,
    skip_temp_heavy_mc=False,
):
    topology_key = topology_cfg["key"]
    topology_display_name = topology_cfg["display_name"]
    topology_dir_name = topology_cfg["dir_name"]
    topology_cuts = topology_cfg["cuts"]

    DATA_PERIOD_FILES = channel_cfg["data_period_files"]
    MC_DIR = channel_cfg["mc_dir"]
    MC_CHANNEL_TAG = channel_cfg["mc_channel_tag"]

    OUTPUT_DIR, INTEGRATED_OUTPUT_DIR, ANGLE_OUTPUT_DIR = get_output_paths(channel, topology_dir_name)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(INTEGRATED_OUTPUT_DIR, exist_ok=True)
    os.makedirs(ANGLE_OUTPUT_DIR, exist_ok=True)

    print("")
    print("------------------------------------------------------------")
    print(f"Selection: {topology_display_name}")
    print(f"Output dir: {OUTPUT_DIR}")
    print("------------------------------------------------------------")

    all_run_rows = []
    all_current_rows = []
    all_skipped_run_rows = []

    angle_current_rows_by_variable = {}
    for var_cfg in ANGLE_DEPENDENCE_CONFIG:
        angle_current_rows_by_variable[var_cfg["key"]] = {}
    #endfor

    print("")
    print("Processing DATA ROOT files and aggregating by current...")

    data_tasks = []
    for period_display_name, period_internal_name, root_path in DATA_PERIOD_FILES:
        data_tasks.append((period_display_name, period_internal_name, root_path, run_charge_map, topology_cuts))
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

        all_run_rows.extend(run_rows)
        all_current_rows.extend(current_rows)
        all_skipped_run_rows.extend(skipped_nonpositive_charge_rows)

        for var_cfg in ANGLE_DEPENDENCE_CONFIG:
            var_key = var_cfg["key"]
            angle_current_rows_by_variable[var_key][period_display_name] = angle_current_rows_this_period[var_key]
        #endfor

        print(f"Run-level rows kept: {len(run_rows)}")
        print(f"Integrated current groups: {len(current_rows)}")
        print(f"Skipped non-positive-charge runs: {len(skipped_nonpositive_charge_rows)}")

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
    if skip_temp_heavy_mc:
        print("Temporary MC skip override flag is ON.")
    #endif

    integrated_mc_rows, angle_mc_rows_by_variable = build_mc_aggregation(
        mc_dir=MC_DIR,
        requested_channel_tag=MC_CHANNEL_TAG,
        topology_cuts=topology_cuts,
        skip_temp_heavy_mc=skip_temp_heavy_mc,
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

    title_suffix = ""
    if topology_key != "overall":
        title_suffix = f"  topology {topology_display_name}"
    #endif

    out_a, out_b, out_c, out_d = plot_four_panel_set(
        output_dir=INTEGRATED_OUTPUT_DIR,
        tag="integrated",
        data_current_rows=all_current_rows,
        mc_rows=integrated_mc_rows,
        data_fit_results=integrated_data_fit_results,
        mc_fit_results=integrated_mc_fit_results,
        period_color=period_color,
        title_suffix=title_suffix,
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
            mc_rows_by_period=angle_mc_rows_by_variable[var_key],
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
                mc_rows_bin = angle_mc_rows_by_variable[var_key][period_name]["rows"][ibin]

                data_fit_bin = make_fit_result_map_for_data(current_rows_bin)
                mc_fit_bin = make_fit_result_map_for_mc(mc_rows_bin)

                label = angle_bin_to_label(low, high, is_last=(ibin == len(bins) - 1))
                angle_title_suffix = title_suffix + f"  {var_cfg['display_name']} {label}"

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
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--channel",
        default="epgamma",
        help="Channel to analyze. Default: epgamma. Currently supported: epgamma, eX",
    )
    parser.add_argument(
        "--skip_temp_heavy_mc",
        action="store_true",
        help="Temporarily skip selected MC points if skip_pairs is uncommented in build_mc_aggregation.",
    )
    args = parser.parse_args()

    channel_cfg = get_channel_config(args.channel)

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
    print(f"Running channel: {args.channel}")
    print(f"MC directory:    {channel_cfg['mc_dir']}")
    print(f"Max workers:     {MAX_WORKERS}")
    print("============================================================")

    print("")
    print("Reading charge map...")
    run_charge_map = read_charge_map(CSV_FILE)
    print(f"Loaded charge entries for {len(run_charge_map)} runs from:")
    print(f"  {CSV_FILE}")

    topology_list = [TOPOLOGY_CONFIG[0]]
    if channel_cfg["supports_topology"]:
        topology_list.extend(TOPOLOGY_CONFIG[1:])
    #endif

    all_summary_outputs = []

    for topology_cfg in topology_list:
        summary_rows = run_selection_analysis(
            channel=args.channel,
            channel_cfg=channel_cfg,
            topology_cfg=topology_cfg,
            run_charge_map=run_charge_map,
            period_color=period_color,
            skip_temp_heavy_mc=args.skip_temp_heavy_mc,
        )
        all_summary_outputs.append((topology_cfg, summary_rows))
    #endfor

    for topology_cfg, summary_rows in all_summary_outputs:
        if topology_cfg["key"] == "overall":
            title = f"Representative integrated period-level normalization summary: overall channel {args.channel}"
        else:
            title = f"Representative integrated period-level normalization summary: {args.channel} topology {topology_cfg['display_name']}"
        #endif

        print_summary_table(title, summary_rows)
    #endfor

    print("")
#enddef


if __name__ == "__main__":
    main()
#endif