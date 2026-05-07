#!/usr/bin/env python3

import os
import sys
import csv
import math
import time
import argparse
import multiprocessing as mp
import ROOT

# -----------------------------------------------------------------------------
# User-editable constants
# -----------------------------------------------------------------------------

GLOBAL_CHARGE_CSV = "/u/home/thayward/clas12_analysis_software/analysis_scripts/dvcs_cross_section/imports/integrated_luminosity/global.csv"

DEFAULT_AAOGEN_SIGMA_INTEGRATED_MICROBARN = 6.04885572e-4

CHARGE_TO_MC_FACTOR = 1.0e-6

RGA_LUMINOSITY_FACTOR_PB_INV_PER_MC = 1316.875

OUTPUT_ROOT_DIR = "output/data_mc_normalization"
DEFAULT_OUTPUT_TAG = "default"

# p1 region definitions:
#   FD: 0 <= p1_theta < 40 deg
#   CD: 40 <= p1_theta < 70 deg
P1_FD_THETA_MIN_DEG = 0.0
P1_FD_THETA_MAX_DEG = 40.0

P1_CD_THETA_MIN_DEG = 40.0
P1_CD_THETA_MAX_DEG = 70.0

# p2 region definitions:
#   FT: 0 <= p2_theta < 5 deg
#   FD: 5 <= p2_theta < 40 deg
P2_FT_THETA_MIN_DEG = 0.0
P2_FT_THETA_MAX_DEG = 5.0

P2_FD_THETA_MIN_DEG = 5.0
P2_FD_THETA_MAX_DEG = 40.0

N_BINS_THETA_FIRST_REGION = 12
N_BINS_THETA_SECOND_REGION = 23

N_BINS_P_FIRST_REGION = 12
N_BINS_P_SECOND_REGION = 23

N_BINS_PHI_FIRST_REGION = 12
N_BINS_PHI_SECOND_REGION = 24

P1_P_MIN_GEV = 0.3
P1_P_MAX_GEV = 1.3

P2_P_MIN_GEV = 2.0
P2_P_MAX_GEV = 6.0

PHI_MIN_DEG = 0.0
PHI_MAX_DEG = 360.0

DEFAULT_STATUS_EVERY = 250000
DEFAULT_MAX_WORKERS = 5

LOG_Y_MIN = 0.5
LINEAR_Y_MIN = 0.0
Y_PADDING_LINEAR = 1.30
Y_PADDING_LOG = 30.0

# Ratio plots are linear only.
RATIO_LINEAR_Y_MIN = 0.2
RATIO_Y_MAX = 1.5

TREE_NAME = "PhysicsEvents"

# -----------------------------------------------------------------------------
# Plot variable configuration
# -----------------------------------------------------------------------------

PLOT_VARIABLES = [
    {
        "particle": "p1",
        "key": "p1_theta",
        "branch": "p1_theta",
        "theta_branch": "p1_theta",
        "phi_branch": "p1_phi",
        "p_branch": "p1_p",
        "title": "p_{1} #theta",
        "x_title": "p_{1} #theta (deg)",
        "first_region_name": "FD",
        "second_region_name": "CD",
        "first_region_theta_min_deg": P1_FD_THETA_MIN_DEG,
        "first_region_theta_max_deg": P1_FD_THETA_MAX_DEG,
        "second_region_theta_min_deg": P1_CD_THETA_MIN_DEG,
        "second_region_theta_max_deg": P1_CD_THETA_MAX_DEG,
        "first_region_n_bins": N_BINS_THETA_FIRST_REGION,
        "second_region_n_bins": N_BINS_THETA_SECOND_REGION,
        "first_region_min": P1_FD_THETA_MIN_DEG,
        "first_region_max": P1_FD_THETA_MAX_DEG,
        "second_region_min": P1_CD_THETA_MIN_DEG,
        "second_region_max": P1_CD_THETA_MAX_DEG,
        "convert": "rad_to_deg",
        "layout": "sector_2x4",
        "n_panels": 7,
    },
    {
        "particle": "p1",
        "key": "p1_p",
        "branch": "p1_p",
        "theta_branch": "p1_theta",
        "phi_branch": "p1_phi",
        "p_branch": "p1_p",
        "title": "p_{1} momentum",
        "x_title": "p_{1} momentum (GeV)",
        "first_region_name": "FD",
        "second_region_name": "CD",
        "first_region_theta_min_deg": P1_FD_THETA_MIN_DEG,
        "first_region_theta_max_deg": P1_FD_THETA_MAX_DEG,
        "second_region_theta_min_deg": P1_CD_THETA_MIN_DEG,
        "second_region_theta_max_deg": P1_CD_THETA_MAX_DEG,
        "first_region_n_bins": N_BINS_P_FIRST_REGION,
        "second_region_n_bins": N_BINS_P_SECOND_REGION,
        "first_region_min": P1_P_MIN_GEV,
        "first_region_max": P1_P_MAX_GEV,
        "second_region_min": P1_P_MIN_GEV,
        "second_region_max": P1_P_MAX_GEV,
        "convert": "identity",
        "layout": "sector_2x4",
        "n_panels": 7,
    },
    {
        "particle": "p1",
        "key": "p1_phi",
        "branch": "p1_phi",
        "theta_branch": "p1_theta",
        "phi_branch": "p1_phi",
        "p_branch": "p1_p",
        "title": "p_{1} #phi",
        "x_title": "p_{1} #phi (deg)",
        "first_region_name": "FD",
        "second_region_name": "CD",
        "first_region_theta_min_deg": P1_FD_THETA_MIN_DEG,
        "first_region_theta_max_deg": P1_FD_THETA_MAX_DEG,
        "second_region_theta_min_deg": P1_CD_THETA_MIN_DEG,
        "second_region_theta_max_deg": P1_CD_THETA_MAX_DEG,
        "first_region_n_bins": N_BINS_PHI_FIRST_REGION,
        "second_region_n_bins": N_BINS_PHI_SECOND_REGION,
        "first_region_min": PHI_MIN_DEG,
        "first_region_max": PHI_MAX_DEG,
        "second_region_min": PHI_MIN_DEG,
        "second_region_max": PHI_MAX_DEG,
        "convert": "rad_to_deg_wrapped",
        "layout": "phi_1x3",
        "n_panels": 2,
    },
    {
        "particle": "p2",
        "key": "p2_theta",
        "branch": "p2_theta",
        "theta_branch": "p2_theta",
        "phi_branch": "p2_phi",
        "p_branch": "p2_p",
        "title": "p_{2} #theta",
        "x_title": "p_{2} #theta (deg)",
        "first_region_name": "FD",
        "second_region_name": "FT",
        "first_region_theta_min_deg": P2_FD_THETA_MIN_DEG,
        "first_region_theta_max_deg": P2_FD_THETA_MAX_DEG,
        "second_region_theta_min_deg": P2_FT_THETA_MIN_DEG,
        "second_region_theta_max_deg": P2_FT_THETA_MAX_DEG,
        "first_region_n_bins": N_BINS_THETA_FIRST_REGION,
        "second_region_n_bins": N_BINS_THETA_SECOND_REGION,
        "first_region_min": P2_FD_THETA_MIN_DEG,
        "first_region_max": P2_FD_THETA_MAX_DEG,
        "second_region_min": P2_FT_THETA_MIN_DEG,
        "second_region_max": P2_FT_THETA_MAX_DEG,
        "convert": "rad_to_deg",
        "layout": "sector_2x4",
        "n_panels": 7,
    },
    {
        "particle": "p2",
        "key": "p2_p",
        "branch": "p2_p",
        "theta_branch": "p2_theta",
        "phi_branch": "p2_phi",
        "p_branch": "p2_p",
        "title": "p_{2} momentum",
        "x_title": "p_{2} momentum (GeV)",
        "first_region_name": "FD",
        "second_region_name": "FT",
        "first_region_theta_min_deg": P2_FD_THETA_MIN_DEG,
        "first_region_theta_max_deg": P2_FD_THETA_MAX_DEG,
        "second_region_theta_min_deg": P2_FT_THETA_MIN_DEG,
        "second_region_theta_max_deg": P2_FT_THETA_MAX_DEG,
        "first_region_n_bins": N_BINS_P_FIRST_REGION,
        "second_region_n_bins": N_BINS_P_SECOND_REGION,
        "first_region_min": P2_P_MIN_GEV,
        "first_region_max": P2_P_MAX_GEV,
        "second_region_min": P2_P_MIN_GEV,
        "second_region_max": P2_P_MAX_GEV,
        "convert": "identity",
        "layout": "sector_2x4",
        "n_panels": 7,
    },
    {
        "particle": "p2",
        "key": "p2_phi",
        "branch": "p2_phi",
        "theta_branch": "p2_theta",
        "phi_branch": "p2_phi",
        "p_branch": "p2_p",
        "title": "p_{2} #phi",
        "x_title": "p_{2} #phi (deg)",
        "first_region_name": "FD",
        "second_region_name": "FT",
        "first_region_theta_min_deg": P2_FD_THETA_MIN_DEG,
        "first_region_theta_max_deg": P2_FD_THETA_MAX_DEG,
        "second_region_theta_min_deg": P2_FT_THETA_MIN_DEG,
        "second_region_theta_max_deg": P2_FT_THETA_MAX_DEG,
        "first_region_n_bins": N_BINS_PHI_FIRST_REGION,
        "second_region_n_bins": N_BINS_PHI_SECOND_REGION,
        "first_region_min": PHI_MIN_DEG,
        "first_region_max": PHI_MAX_DEG,
        "second_region_min": PHI_MIN_DEG,
        "second_region_max": PHI_MAX_DEG,
        "convert": "rad_to_deg_wrapped",
        "layout": "phi_1x3",
        "n_panels": 2,
    },
]

# -----------------------------------------------------------------------------
# Basic helpers
# -----------------------------------------------------------------------------

def fatal(message):
    print("ERROR: {}".format(message), file=sys.stderr, flush=True)
    sys.exit(1)


def status(message):
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print("[{}] {}".format(timestamp, message), flush=True)


def format_percent(done, total):
    if total <= 0:
        return "0.0%"
    #endif

    return "{:.1f}%".format(100.0 * float(done) / float(total))


def sanitize_output_tag(output_tag):
    sanitized = output_tag.strip()

    if sanitized == "":
        fatal("output_tag cannot be an empty string")
    #endif

    if "/" in sanitized or "\\" in sanitized:
        fatal("output_tag must be a single directory name, not a path: {}".format(output_tag))
    #endif

    return sanitized


def ensure_output_directory(output_path):
    output_dir = os.path.dirname(output_path)

    if output_dir != "":
        os.makedirs(output_dir, exist_ok=True)
    #endif


def build_output_paths(output_tag, variable_config):
    tag = sanitize_output_tag(output_tag)
    particle = variable_config["particle"]
    variable_key = variable_config["key"]

    output_dir = os.path.join(OUTPUT_ROOT_DIR, tag, particle)
    base = os.path.join(output_dir, "output_{}".format(variable_key))

    return {
        "output_dir": output_dir,
        "comparison_log": base + "_log.pdf",
        "comparison_linear": base + "_linear.pdf",
        "ratio_linear": base + "_ratio_linear.pdf",
    }


def open_root_file(path, label):
    status("Opening {} ROOT file: {}".format(label, path))

    if not os.path.isfile(path):
        fatal("{} ROOT file does not exist: {}".format(label, path))
    #endif

    root_file = ROOT.TFile.Open(path, "READ")

    if root_file is None or root_file.IsZombie():
        fatal("could not open {} ROOT file: {}".format(label, path))
    #endif

    status("Opened {} ROOT file.".format(label))

    return root_file


def get_tree(root_file, tree_name, label):
    status("Loading tree '{}' from {} ROOT file.".format(tree_name, label))

    tree = root_file.Get(tree_name)

    if tree is None:
        fatal("{} ROOT file does not contain tree '{}': {}".format(label, tree_name, root_file.GetName()))
    #endif

    status("Loaded tree '{}' from {} ROOT file with {} entries.".format(tree_name, label, tree.GetEntries()))

    return tree


def has_branch(tree, branch_name):
    return tree.GetBranch(branch_name) is not None


def require_branches(tree, branch_names, label):
    status("Checking required branches for {} tree.".format(label))

    missing = []

    for branch_name in branch_names:
        if not has_branch(tree, branch_name):
            missing.append(branch_name)
        #endif
    #endfor

    if len(missing) > 0:
        fatal("{} tree is missing required branches: {}".format(label, ", ".join(missing)))
    #endif

    status("All required branches are present for {} tree.".format(label))


def get_entry_count(root_path, tree_name, label):
    root_file = open_root_file(root_path, label)
    tree = get_tree(root_file, tree_name, label)
    n_entries = tree.GetEntries()
    root_file.Close()

    return n_entries


def check_input_tree(root_path, tree_name, branch_names, label):
    root_file = open_root_file(root_path, label)
    tree = get_tree(root_file, tree_name, label)
    require_branches(tree, branch_names, label)
    n_entries = tree.GetEntries()
    root_file.Close()

    return n_entries


def wrap_phi_0_360(phi_deg):
    phi = math.fmod(phi_deg, 360.0)

    if phi < 0.0:
        phi += 360.0
    #endif

    return phi


def convert_value(raw_value, conversion_mode):
    if conversion_mode == "identity":
        return raw_value
    #endif

    if conversion_mode == "rad_to_deg":
        return math.degrees(raw_value)
    #endif

    if conversion_mode == "rad_to_deg_wrapped":
        return wrap_phi_0_360(math.degrees(raw_value))
    #endif

    fatal("unknown conversion mode: {}".format(conversion_mode))


def get_fd_sector_from_phi_deg(phi_deg):
    phi = wrap_phi_0_360(phi_deg)

    if phi >= 330.0 or phi < 30.0:
        return 1
    #endif

    if phi >= 30.0 and phi < 90.0:
        return 2
    #endif

    if phi >= 90.0 and phi < 150.0:
        return 3
    #endif

    if phi >= 150.0 and phi < 210.0:
        return 4
    #endif

    if phi >= 210.0 and phi < 270.0:
        return 5
    #endif

    if phi >= 270.0 and phi < 330.0:
        return 6
    #endif

    return 0


def get_region_from_theta(variable_config, theta_rad):
    theta_deg = math.degrees(theta_rad)

    first_min = variable_config["first_region_theta_min_deg"]
    first_max = variable_config["first_region_theta_max_deg"]
    second_min = variable_config["second_region_theta_min_deg"]
    second_max = variable_config["second_region_theta_max_deg"]

    if theta_deg >= first_min and theta_deg < first_max:
        return "first"
    #endif

    if theta_deg >= second_min and theta_deg < second_max:
        return "second"
    #endif

    return "none"


def get_sector_panel_index_from_theta_and_phi(variable_config, theta_rad, phi_rad):
    region = get_region_from_theta(variable_config, theta_rad)

    if region == "first":
        phi_deg = math.degrees(phi_rad)
        sector = get_fd_sector_from_phi_deg(phi_deg)

        if sector >= 1 and sector <= 6:
            return sector - 1
        #endif

        return -1
    #endif

    if region == "second":
        return 6
    #endif

    return -1


def get_variable_panel_index(variable_config, theta_rad, phi_rad):
    if variable_config["layout"] == "sector_2x4":
        return get_sector_panel_index_from_theta_and_phi(variable_config, theta_rad, phi_rad)
    #endif

    if variable_config["layout"] == "phi_1x3":
        region = get_region_from_theta(variable_config, theta_rad)

        if region == "first":
            return 0
        #endif

        if region == "second":
            return 1
        #endif

        return -1
    #endif

    fatal("unknown layout: {}".format(variable_config["layout"]))


def get_plot_range_for_panel(variable_config, i_panel):
    if variable_config["layout"] == "sector_2x4":
        if i_panel >= 0 and i_panel <= 5:
            return variable_config["first_region_min"], variable_config["first_region_max"]
        #endif

        if i_panel == 6:
            return variable_config["second_region_min"], variable_config["second_region_max"]
        #endif
    #endif

    if variable_config["layout"] == "phi_1x3":
        if i_panel == 0:
            return variable_config["first_region_min"], variable_config["first_region_max"]
        #endif

        if i_panel == 1:
            return variable_config["second_region_min"], variable_config["second_region_max"]
        #endif
    #endif

    fatal("invalid panel index {} for variable {}".format(i_panel, variable_config["key"]))


def get_n_bins_for_panel(variable_config, i_panel):
    if variable_config["layout"] == "sector_2x4":
        if i_panel >= 0 and i_panel <= 5:
            return variable_config["first_region_n_bins"]
        #endif

        if i_panel == 6:
            return variable_config["second_region_n_bins"]
        #endif
    #endif

    if variable_config["layout"] == "phi_1x3":
        if i_panel == 0:
            return variable_config["first_region_n_bins"]
        #endif

        if i_panel == 1:
            return variable_config["second_region_n_bins"]
        #endif
    #endif

    fatal("invalid panel index {} for variable {}".format(i_panel, variable_config["key"]))


def get_panel_labels(variable_config):
    first_name = variable_config["first_region_name"]
    second_name = variable_config["second_region_name"]

    if variable_config["layout"] == "sector_2x4":
        return [
            "{} sector 1".format(first_name),
            "{} sector 2".format(first_name),
            "{} sector 3".format(first_name),
            "{} sector 4".format(first_name),
            "{} sector 5".format(first_name),
            "{} sector 6".format(first_name),
            second_name,
        ]
    #endif

    if variable_config["layout"] == "phi_1x3":
        return [
            "{} combined".format(first_name),
            second_name,
        ]
    #endif

    fatal("unknown layout: {}".format(variable_config["layout"]))


def get_value_bin_for_panel(variable_config, i_panel, value):
    value_min, value_max = get_plot_range_for_panel(variable_config, i_panel)
    n_bins = get_n_bins_for_panel(variable_config, i_panel)

    if value < value_min or value >= value_max:
        return -1
    #endif

    bin_width = (value_max - value_min) / float(n_bins)
    i_bin = int((value - value_min) / bin_width)

    if i_bin < 0 or i_bin >= n_bins:
        return -1
    #endif

    return i_bin


def make_empty_counts_for_variable(variable_config):
    counts = []

    for i_panel in range(variable_config["n_panels"]):
        n_bins = get_n_bins_for_panel(variable_config, i_panel)
        counts.append([0.0 for _ in range(n_bins)])
    #endfor

    return counts


def make_empty_sumw2_for_variable(variable_config):
    sumw2 = []

    for i_panel in range(variable_config["n_panels"]):
        n_bins = get_n_bins_for_panel(variable_config, i_panel)
        sumw2.append([0.0 for _ in range(n_bins)])
    #endfor

    return sumw2


def add_counts_for_variable(total_counts, chunk_counts, variable_config):
    for i_panel in range(variable_config["n_panels"]):
        n_bins = get_n_bins_for_panel(variable_config, i_panel)

        for i_bin in range(n_bins):
            total_counts[i_panel][i_bin] += chunk_counts[i_panel][i_bin]
        #endfor
    #endfor


def make_ranges(n_entries, n_workers):
    ranges = []

    if n_entries <= 0:
        return ranges
    #endif

    n_workers = max(1, min(n_workers, n_entries))
    chunk_size = int(math.ceil(float(n_entries) / float(n_workers)))

    start = 0

    while start < n_entries:
        end = min(start + chunk_size, n_entries)
        ranges.append((start, end))
        start = end
    #endwhile

    return ranges


def get_variable_keys():
    keys = []

    for variable_config in PLOT_VARIABLES:
        keys.append(variable_config["key"])
    #endfor

    return keys


def get_required_plot_branches():
    required = set()

    for variable_config in PLOT_VARIABLES:
        required.add(variable_config["branch"])
        required.add(variable_config["theta_branch"])
        required.add(variable_config["phi_branch"])
        required.add(variable_config["p_branch"])
    #endfor

    return sorted(required)


# -----------------------------------------------------------------------------
# Charge helpers
# -----------------------------------------------------------------------------

def load_charge_map(csv_path):
    status("Loading charge map from CSV: {}".format(csv_path))

    if not os.path.isfile(csv_path):
        fatal("charge CSV does not exist: {}".format(csv_path))
    #endif

    charge_by_run = {}
    n_rows_read = 0
    n_rows_used = 0

    with open(csv_path, "r") as f:
        reader = csv.reader(f)

        for row in reader:
            n_rows_read += 1

            if len(row) < 2:
                continue
            #endif

            try:
                runnum = int(row[0])
                total_charge_raw = float(row[1])
            except Exception:
                continue
            #endtry

            charge_by_run[runnum] = total_charge_raw
            n_rows_used += 1
        #endfor
    #endwith

    if len(charge_by_run) == 0:
        fatal("no valid run-charge rows were read from {}".format(csv_path))
    #endif

    status("Loaded charge map: {} valid runs from {} CSV rows.".format(n_rows_used, n_rows_read))

    return charge_by_run


def sum_charge_for_runs(unique_runs, charge_by_run):
    status("Summing accumulated charge for {} unique data runs.".format(len(unique_runs)))

    missing_runs = []
    total_charge_raw = 0.0

    for runnum in sorted(unique_runs):
        if runnum not in charge_by_run:
            missing_runs.append(runnum)
            continue
        #endif

        total_charge_raw += charge_by_run[runnum]
    #endfor

    if len(missing_runs) > 0:
        fatal("these run numbers from the data file were not found in the charge CSV: {}".format(
            ", ".join(str(runnum) for runnum in missing_runs)
        ))
    #endif

    status("Finished charge sum. Raw accumulated charge sum: {:.12g}".format(total_charge_raw))

    return total_charge_raw


# -----------------------------------------------------------------------------
# Worker functions
# -----------------------------------------------------------------------------

def configure_tree_for_data(tree):
    tree.SetBranchStatus("*", 0)
    tree.SetBranchStatus("runnum", 1)

    for branch_name in get_required_plot_branches():
        tree.SetBranchStatus(branch_name, 1)
    #endfor


def configure_tree_for_mc(tree):
    tree.SetBranchStatus("*", 0)

    for branch_name in get_required_plot_branches():
        tree.SetBranchStatus(branch_name, 1)
    #endfor


def make_empty_result_maps():
    counts_by_variable = {}
    sumw2_by_variable = {}

    for variable_config in PLOT_VARIABLES:
        key = variable_config["key"]
        counts_by_variable[key] = make_empty_counts_for_variable(variable_config)
        sumw2_by_variable[key] = make_empty_sumw2_for_variable(variable_config)
    #endfor

    return counts_by_variable, sumw2_by_variable


def fill_variable_counts(counts_by_variable, sumw2_by_variable, variable_config, theta_rad, phi_rad, raw_value, event_weight):
    key = variable_config["key"]
    i_panel = get_variable_panel_index(variable_config, theta_rad, phi_rad)

    if i_panel < 0:
        return False, "panel"
    #endif

    value = convert_value(raw_value, variable_config["convert"])
    i_bin = get_value_bin_for_panel(variable_config, i_panel, value)

    if i_bin < 0:
        return False, "range"
    #endif

    counts_by_variable[key][i_panel][i_bin] += event_weight
    sumw2_by_variable[key][i_panel][i_bin] += event_weight * event_weight

    return True, "filled"


def initialize_variable_counters():
    variable_fills = {}
    variable_panel_skips = {}
    variable_range_skips = {}

    for key in get_variable_keys():
        variable_fills[key] = 0
        variable_panel_skips[key] = 0
        variable_range_skips[key] = 0
    #endfor

    return variable_fills, variable_panel_skips, variable_range_skips


def initialize_region_counters():
    region_counts = {
        "p1_first": 0,
        "p1_second": 0,
        "p1_none": 0,
        "p2_first": 0,
        "p2_second": 0,
        "p2_none": 0,
    }

    return region_counts


def update_region_counters(region_counts, particle, region):
    key = "{}_{}".format(particle, region)

    if key not in region_counts:
        fatal("internal error: unknown region counter key {}".format(key))
    #endif

    region_counts[key] += 1


def data_worker(args):
    root_path = args["root_path"]
    tree_name = args["tree_name"]
    start_entry = args["start_entry"]
    end_entry = args["end_entry"]
    worker_id = args["worker_id"]
    status_every = args["status_every"]

    root_file = ROOT.TFile.Open(root_path, "READ")

    if root_file is None or root_file.IsZombie():
        raise RuntimeError("worker {} could not open data ROOT file: {}".format(worker_id, root_path))
    #endif

    tree = root_file.Get(tree_name)

    if tree is None:
        root_file.Close()
        raise RuntimeError("worker {} could not load tree {}".format(worker_id, tree_name))
    #endif

    configure_tree_for_data(tree)

    counts_by_variable, sumw2_by_variable = make_empty_result_maps()
    unique_runs = set()

    region_counts = initialize_region_counters()
    variable_fills, variable_panel_skips, variable_range_skips = initialize_variable_counters()

    n_total = end_entry - start_entry

    print("[{}] data worker {} starting entries {} to {}.".format(
        time.strftime("%Y-%m-%d %H:%M:%S"),
        worker_id,
        start_entry,
        end_entry
    ), flush=True)

    for local_index, i_entry in enumerate(range(start_entry, end_entry), start=1):
        tree.GetEntry(i_entry)

        runnum = int(getattr(tree, "runnum"))
        unique_runs.add(runnum)

        raw_values = {
            "p1_theta": float(getattr(tree, "p1_theta")),
            "p1_phi": float(getattr(tree, "p1_phi")),
            "p1_p": float(getattr(tree, "p1_p")),
            "p2_theta": float(getattr(tree, "p2_theta")),
            "p2_phi": float(getattr(tree, "p2_phi")),
            "p2_p": float(getattr(tree, "p2_p")),
        }

        for variable_config in PLOT_VARIABLES:
            key = variable_config["key"]
            particle = variable_config["particle"]
            theta_rad = raw_values[variable_config["theta_branch"]]
            phi_rad = raw_values[variable_config["phi_branch"]]
            raw_value = raw_values[variable_config["branch"]]
            event_weight = 1.0

            region = get_region_from_theta(variable_config, theta_rad)
            update_region_counters(region_counts, particle, region)

            filled, reason = fill_variable_counts(
                counts_by_variable,
                sumw2_by_variable,
                variable_config,
                theta_rad,
                phi_rad,
                raw_value,
                event_weight
            )

            if filled:
                variable_fills[key] += 1
            elif reason == "panel":
                variable_panel_skips[key] += 1
            elif reason == "range":
                variable_range_skips[key] += 1
            #endif
        #endfor

        if local_index % status_every == 0:
            print("[{}] data worker {} progress: {}/{} chunk entries ({}) | p1 first/second/none = {}/{}/{} | p2 first/second/none = {}/{}/{} | unique runs: {}".format(
                time.strftime("%Y-%m-%d %H:%M:%S"),
                worker_id,
                local_index,
                n_total,
                format_percent(local_index, n_total),
                region_counts["p1_first"],
                region_counts["p1_second"],
                region_counts["p1_none"],
                region_counts["p2_first"],
                region_counts["p2_second"],
                region_counts["p2_none"],
                len(unique_runs)
            ), flush=True)
        #endif
    #endfor

    root_file.Close()

    print("[{}] data worker {} finished: p1 first/second/none = {}/{}/{}, p2 first/second/none = {}/{}/{}, unique runs = {}.".format(
        time.strftime("%Y-%m-%d %H:%M:%S"),
        worker_id,
        region_counts["p1_first"],
        region_counts["p1_second"],
        region_counts["p1_none"],
        region_counts["p2_first"],
        region_counts["p2_second"],
        region_counts["p2_none"],
        len(unique_runs)
    ), flush=True)

    return {
        "counts_by_variable": counts_by_variable,
        "sumw2_by_variable": sumw2_by_variable,
        "unique_runs": sorted(unique_runs),
        "region_counts": region_counts,
        "variable_fills": variable_fills,
        "variable_panel_skips": variable_panel_skips,
        "variable_range_skips": variable_range_skips,
    }


def mc_worker(args):
    root_path = args["root_path"]
    tree_name = args["tree_name"]
    start_entry = args["start_entry"]
    end_entry = args["end_entry"]
    worker_id = args["worker_id"]
    status_every = args["status_every"]
    normalization_factor = args["normalization_factor"]

    root_file = ROOT.TFile.Open(root_path, "READ")

    if root_file is None or root_file.IsZombie():
        raise RuntimeError("worker {} could not open MC ROOT file: {}".format(worker_id, root_path))
    #endif

    tree = root_file.Get(tree_name)

    if tree is None:
        root_file.Close()
        raise RuntimeError("worker {} could not load tree {}".format(worker_id, tree_name))
    #endif

    configure_tree_for_mc(tree)

    counts_by_variable, sumw2_by_variable = make_empty_result_maps()

    region_counts = initialize_region_counters()
    variable_fills, variable_panel_skips, variable_range_skips = initialize_variable_counters()

    n_total = end_entry - start_entry
    filled_event_count = 0
    scaled_weight_sum = 0.0

    print("[{}] MC worker {} starting entries {} to {}.".format(
        time.strftime("%Y-%m-%d %H:%M:%S"),
        worker_id,
        start_entry,
        end_entry
    ), flush=True)

    for local_index, i_entry in enumerate(range(start_entry, end_entry), start=1):
        tree.GetEntry(i_entry)

        event_weight = normalization_factor

        raw_values = {
            "p1_theta": float(getattr(tree, "p1_theta")),
            "p1_phi": float(getattr(tree, "p1_phi")),
            "p1_p": float(getattr(tree, "p1_p")),
            "p2_theta": float(getattr(tree, "p2_theta")),
            "p2_phi": float(getattr(tree, "p2_phi")),
            "p2_p": float(getattr(tree, "p2_p")),
        }

        any_filled = False

        for variable_config in PLOT_VARIABLES:
            key = variable_config["key"]
            particle = variable_config["particle"]
            theta_rad = raw_values[variable_config["theta_branch"]]
            phi_rad = raw_values[variable_config["phi_branch"]]
            raw_value = raw_values[variable_config["branch"]]

            region = get_region_from_theta(variable_config, theta_rad)
            update_region_counters(region_counts, particle, region)

            filled, reason = fill_variable_counts(
                counts_by_variable,
                sumw2_by_variable,
                variable_config,
                theta_rad,
                phi_rad,
                raw_value,
                event_weight
            )

            if filled:
                variable_fills[key] += 1
                any_filled = True
            elif reason == "panel":
                variable_panel_skips[key] += 1
            elif reason == "range":
                variable_range_skips[key] += 1
            #endif
        #endfor

        if any_filled:
            filled_event_count += 1
            scaled_weight_sum += event_weight
        #endif

        if local_index % status_every == 0:
            print("[{}] MC worker {} progress: {}/{} chunk entries ({}) | p1 first/second/none = {}/{}/{} | p2 first/second/none = {}/{}/{} | any-filled entries = {} | constant event weight = {:.12g}".format(
                time.strftime("%Y-%m-%d %H:%M:%S"),
                worker_id,
                local_index,
                n_total,
                format_percent(local_index, n_total),
                region_counts["p1_first"],
                region_counts["p1_second"],
                region_counts["p1_none"],
                region_counts["p2_first"],
                region_counts["p2_second"],
                region_counts["p2_none"],
                filled_event_count,
                event_weight
            ), flush=True)
        #endif
    #endfor

    root_file.Close()

    print("[{}] MC worker {} finished: p1 first/second/none = {}/{}/{}, p2 first/second/none = {}/{}/{}, any-filled entries = {}, scaled weight sum = {:.12g}, constant event weight = {:.12g}.".format(
        time.strftime("%Y-%m-%d %H:%M:%S"),
        worker_id,
        region_counts["p1_first"],
        region_counts["p1_second"],
        region_counts["p1_none"],
        region_counts["p2_first"],
        region_counts["p2_second"],
        region_counts["p2_none"],
        filled_event_count,
        scaled_weight_sum,
        normalization_factor
    ), flush=True)

    return {
        "counts_by_variable": counts_by_variable,
        "sumw2_by_variable": sumw2_by_variable,
        "region_counts": region_counts,
        "filled_event_count": filled_event_count,
        "scaled_weight_sum": scaled_weight_sum,
        "constant_event_weight": normalization_factor,
        "variable_fills": variable_fills,
        "variable_panel_skips": variable_panel_skips,
        "variable_range_skips": variable_range_skips,
    }


def merge_region_counts(total_region_counts, chunk_region_counts):
    for key in total_region_counts:
        total_region_counts[key] += chunk_region_counts[key]
    #endfor


def run_data_parallel(root_path, tree_name, n_entries, max_workers, status_every):
    n_workers = min(max_workers, n_entries)

    if n_workers < 1:
        fatal("data tree has zero entries")
    #endif

    ranges = make_ranges(n_entries, n_workers)

    tasks = []

    for i_worker, entry_range in enumerate(ranges):
        start_entry, end_entry = entry_range

        tasks.append({
            "root_path": root_path,
            "tree_name": tree_name,
            "start_entry": start_entry,
            "end_entry": end_entry,
            "worker_id": i_worker + 1,
            "status_every": status_every,
        })
    #endfor

    status("Starting parallel data pass with {} workers.".format(len(tasks)))

    total_counts_by_variable, total_sumw2_by_variable = make_empty_result_maps()
    unique_runs = set()
    total_region_counts = initialize_region_counters()
    variable_fills, variable_panel_skips, variable_range_skips = initialize_variable_counters()

    with mp.Pool(processes=len(tasks)) as pool:
        results = pool.map(data_worker, tasks)
    #endwith

    for result in results:
        for variable_config in PLOT_VARIABLES:
            key = variable_config["key"]

            add_counts_for_variable(
                total_counts_by_variable[key],
                result["counts_by_variable"][key],
                variable_config
            )

            add_counts_for_variable(
                total_sumw2_by_variable[key],
                result["sumw2_by_variable"][key],
                variable_config
            )

            variable_fills[key] += result["variable_fills"][key]
            variable_panel_skips[key] += result["variable_panel_skips"][key]
            variable_range_skips[key] += result["variable_range_skips"][key]
        #endfor

        unique_runs.update(result["unique_runs"])
        merge_region_counts(total_region_counts, result["region_counts"])
    #endfor

    status("Finished parallel data pass.")
    status("Data p1 region totals: FD = {}, CD = {}, none = {}.".format(
        total_region_counts["p1_first"],
        total_region_counts["p1_second"],
        total_region_counts["p1_none"]
    ))
    status("Data p2 region totals: FD = {}, FT = {}, none = {}.".format(
        total_region_counts["p2_first"],
        total_region_counts["p2_second"],
        total_region_counts["p2_none"]
    ))
    status("Data unique runs = {}.".format(len(unique_runs)))

    for key in get_variable_keys():
        status("Data variable {}: filled = {}, skipped panel = {}, skipped range = {}.".format(
            key,
            variable_fills[key],
            variable_panel_skips[key],
            variable_range_skips[key]
        ))
    #endfor

    return {
        "counts_by_variable": total_counts_by_variable,
        "sumw2_by_variable": total_sumw2_by_variable,
        "unique_runs": unique_runs,
        "region_counts": total_region_counts,
        "variable_fills": variable_fills,
        "variable_panel_skips": variable_panel_skips,
        "variable_range_skips": variable_range_skips,
    }


def run_mc_parallel(root_path, tree_name, n_entries, max_workers, status_every, normalization_factor):
    n_workers = min(max_workers, n_entries)

    if n_workers < 1:
        fatal("reconstructed MC tree has zero entries")
    #endif

    ranges = make_ranges(n_entries, n_workers)

    tasks = []

    for i_worker, entry_range in enumerate(ranges):
        start_entry, end_entry = entry_range

        tasks.append({
            "root_path": root_path,
            "tree_name": tree_name,
            "start_entry": start_entry,
            "end_entry": end_entry,
            "worker_id": i_worker + 1,
            "status_every": status_every,
            "normalization_factor": normalization_factor,
        })
    #endfor

    status("Starting parallel reconstructed MC pass with {} workers.".format(len(tasks)))

    total_counts_by_variable, total_sumw2_by_variable = make_empty_result_maps()
    total_region_counts = initialize_region_counters()
    variable_fills, variable_panel_skips, variable_range_skips = initialize_variable_counters()

    filled_event_count = 0
    scaled_weight_sum = 0.0

    with mp.Pool(processes=len(tasks)) as pool:
        results = pool.map(mc_worker, tasks)
    #endwith

    for result in results:
        for variable_config in PLOT_VARIABLES:
            key = variable_config["key"]

            add_counts_for_variable(
                total_counts_by_variable[key],
                result["counts_by_variable"][key],
                variable_config
            )

            add_counts_for_variable(
                total_sumw2_by_variable[key],
                result["sumw2_by_variable"][key],
                variable_config
            )

            variable_fills[key] += result["variable_fills"][key]
            variable_panel_skips[key] += result["variable_panel_skips"][key]
            variable_range_skips[key] += result["variable_range_skips"][key]
        #endfor

        merge_region_counts(total_region_counts, result["region_counts"])

        filled_event_count += result["filled_event_count"]
        scaled_weight_sum += result["scaled_weight_sum"]
    #endfor

    status("Finished parallel reconstructed MC pass.")
    status("Reco MC p1 region totals: FD = {}, CD = {}, none = {}.".format(
        total_region_counts["p1_first"],
        total_region_counts["p1_second"],
        total_region_counts["p1_none"]
    ))
    status("Reco MC p2 region totals: FD = {}, FT = {}, none = {}.".format(
        total_region_counts["p2_first"],
        total_region_counts["p2_second"],
        total_region_counts["p2_none"]
    ))
    status("Reco MC normalization summary over any-filled entries: any-filled entries = {}, scaled weight sum = {:.12g}, constant event weight = {:.12g}.".format(
        filled_event_count,
        scaled_weight_sum,
        normalization_factor
    ))

    for key in get_variable_keys():
        status("Reco MC variable {}: filled = {}, skipped panel = {}, skipped range = {}.".format(
            key,
            variable_fills[key],
            variable_panel_skips[key],
            variable_range_skips[key]
        ))
    #endfor

    return {
        "counts_by_variable": total_counts_by_variable,
        "sumw2_by_variable": total_sumw2_by_variable,
        "region_counts": total_region_counts,
        "filled_event_count": filled_event_count,
        "scaled_weight_sum": scaled_weight_sum,
        "constant_event_weight": normalization_factor,
        "variable_fills": variable_fills,
        "variable_panel_skips": variable_panel_skips,
        "variable_range_skips": variable_range_skips,
    }


# -----------------------------------------------------------------------------
# Histogram and plotting helpers
# -----------------------------------------------------------------------------

def arrays_to_histograms(prefix, variable_config, counts, sumw2):
    histograms = []
    panel_labels = get_panel_labels(variable_config)

    for i_panel in range(variable_config["n_panels"]):
        value_min, value_max = get_plot_range_for_panel(variable_config, i_panel)
        n_bins = get_n_bins_for_panel(variable_config, i_panel)

        hist_name = "{}_{}_panel_{}".format(prefix, variable_config["key"], i_panel + 1)
        hist_title = "{};{};Counts".format(panel_labels[i_panel], variable_config["x_title"])

        hist = ROOT.TH1D(
            hist_name,
            hist_title,
            n_bins,
            value_min,
            value_max
        )

        hist.Sumw2()

        for i_bin in range(n_bins):
            root_bin = i_bin + 1
            content = counts[i_panel][i_bin]
            error = math.sqrt(sumw2[i_panel][i_bin])

            hist.SetBinContent(root_bin, content)
            hist.SetBinError(root_bin, error)
        #endfor

        histograms.append(hist)
    #endfor

    return histograms


def make_ratio_graphs(variable_config, data_histograms, mc_histograms):
    ratio_graphs = []

    for i_panel in range(variable_config["n_panels"]):
        value_min, value_max = get_plot_range_for_panel(variable_config, i_panel)
        n_bins = get_n_bins_for_panel(variable_config, i_panel)
        bin_width = (value_max - value_min) / float(n_bins)

        graph = ROOT.TGraphErrors()
        graph.SetName("ratio_{}_panel_{}".format(variable_config["key"], i_panel + 1))
        graph.SetTitle("ratio_{}_panel_{};{};data / MC".format(
            variable_config["key"],
            i_panel + 1,
            variable_config["x_title"]
        ))

        h_data = data_histograms[i_panel]
        h_mc = mc_histograms[i_panel]

        point_index = 0

        for i_bin in range(1, n_bins + 1):
            data_content = h_data.GetBinContent(i_bin)
            data_error = h_data.GetBinError(i_bin)
            mc_content = h_mc.GetBinContent(i_bin)
            mc_error = h_mc.GetBinError(i_bin)

            if data_content > 0.0 and mc_content > 0.0:
                ratio = data_content / mc_content
                rel_data = data_error / data_content
                rel_mc = mc_error / mc_content
                ratio_error = ratio * math.sqrt(rel_data * rel_data + rel_mc * rel_mc)

                x_center = value_min + (float(i_bin) - 0.5) * bin_width

                graph.SetPoint(point_index, x_center, ratio)
                graph.SetPointError(point_index, 0.0, ratio_error)
                point_index += 1
            #endif
        #endfor

        ratio_graphs.append(graph)
    #endfor

    return ratio_graphs


def style_histograms(data_histograms, mc_histograms, ratio_graphs):
    status("Styling histograms and ratio graphs.")

    for hist in data_histograms:
        hist.SetLineColor(ROOT.kBlue)
        hist.SetMarkerColor(ROOT.kBlue)
        hist.SetMarkerStyle(20)
        hist.SetMarkerSize(0.0)
        hist.SetLineWidth(2)
        hist.SetStats(False)
    #endfor

    for hist in mc_histograms:
        hist.SetLineColor(ROOT.kRed)
        hist.SetMarkerColor(ROOT.kRed)
        hist.SetMarkerStyle(24)
        hist.SetMarkerSize(0.0)
        hist.SetLineWidth(2)
        hist.SetStats(False)
    #endfor

    for graph in ratio_graphs:
        graph.SetLineColor(ROOT.kBlack)
        graph.SetMarkerColor(ROOT.kBlack)
        graph.SetMarkerStyle(20)
        graph.SetMarkerSize(0.95)
        graph.SetLineWidth(1)
    #endfor


def get_hist_positive_min_and_max(hist):
    min_positive = None
    max_value = 0.0

    for i_bin in range(1, hist.GetNbinsX() + 1):
        content = hist.GetBinContent(i_bin)
        error = hist.GetBinError(i_bin)
        upper = content + error
        lower = content - error

        if upper > max_value:
            max_value = upper
        #endif

        if content > 0.0:
            candidate = content

            if lower > 0.0:
                candidate = min(candidate, lower)
            #endif

            if min_positive is None or candidate < min_positive:
                min_positive = candidate
            #endif
        #endif
    #endfor

    return min_positive, max_value


def get_common_y_range_for_hist_list(histograms, log_y):
    global_min_positive = None
    global_max = 0.0

    for hist in histograms:
        min_positive, max_value = get_hist_positive_min_and_max(hist)

        if max_value > global_max:
            global_max = max_value
        #endif

        if min_positive is not None:
            if global_min_positive is None or min_positive < global_min_positive:
                global_min_positive = min_positive
            #endif
        #endif
    #endfor

    if global_max <= 0.0:
        global_max = 1.0
    #endif

    if log_y:
        y_min = LOG_Y_MIN

        if global_min_positive is not None:
            y_min = max(LOG_Y_MIN, 0.5 * global_min_positive)
        #endif

        y_max = Y_PADDING_LOG * global_max
    else:
        y_min = LINEAR_Y_MIN
        y_max = Y_PADDING_LINEAR * global_max
    #endif

    if y_max <= y_min:
        y_max = y_min + 1.0
    #endif

    return y_min, y_max


def get_common_y_ranges_for_comparison(variable_config, data_histograms, mc_histograms, log_y):
    if variable_config["layout"] == "sector_2x4":
        first_region_histograms = []
        second_region_histograms = []

        for i_panel in range(6):
            first_region_histograms.append(data_histograms[i_panel])
            first_region_histograms.append(mc_histograms[i_panel])
        #endfor

        second_region_histograms.append(data_histograms[6])
        second_region_histograms.append(mc_histograms[6])

        first_y_min, first_y_max = get_common_y_range_for_hist_list(first_region_histograms, log_y)
        second_y_min, second_y_max = get_common_y_range_for_hist_list(second_region_histograms, log_y)

        return {
            "first": (first_y_min, first_y_max),
            "second": (second_y_min, second_y_max),
        }
    #endif

    if variable_config["layout"] == "phi_1x3":
        first_region_histograms = [data_histograms[0], mc_histograms[0]]
        second_region_histograms = [data_histograms[1], mc_histograms[1]]

        first_y_min, first_y_max = get_common_y_range_for_hist_list(first_region_histograms, log_y)
        second_y_min, second_y_max = get_common_y_range_for_hist_list(second_region_histograms, log_y)

        return {
            "first": (first_y_min, first_y_max),
            "second": (second_y_min, second_y_max),
        }
    #endif

    fatal("unknown layout: {}".format(variable_config["layout"]))


def get_comparison_y_range_for_panel(variable_config, y_ranges, i_panel):
    if variable_config["layout"] == "sector_2x4":
        if i_panel >= 0 and i_panel <= 5:
            return y_ranges["first"]
        #endif

        if i_panel == 6:
            return y_ranges["second"]
        #endif
    #endif

    if variable_config["layout"] == "phi_1x3":
        if i_panel == 0:
            return y_ranges["first"]
        #endif

        if i_panel == 1:
            return y_ranges["second"]
        #endif
    #endif

    fatal("invalid panel {} for variable {}".format(i_panel, variable_config["key"]))


def draw_normalization_pad(output_tag, variable_config, total_charge_mc, integrated_luminosity_pb_inv, n_gen, normalization_factor, ratio_mode):
    pad = ROOT.gPad
    pad.Clear()
    pad.SetFillColor(ROOT.kWhite)

    latex = ROOT.TLatex()
    latex.SetNDC(True)

    if variable_config["layout"] == "phi_1x3":
        latex.SetTextSize(0.040)
        x0 = 0.10
    else:
        latex.SetTextSize(0.034)
        x0 = 0.10
    #endif

    if ratio_mode:
        title = "Ratio normalization"
    else:
        title = "MC normalization"
    #endif

    latex.DrawLatex(x0, 0.90, output_tag)
    latex.DrawLatex(x0, 0.80, variable_config["particle"])
    latex.DrawLatex(x0, 0.70, variable_config["title"])
    latex.DrawLatex(x0, 0.60, title)
    latex.DrawLatex(x0, 0.49, "Q = %.6g mC" % total_charge_mc)
    latex.DrawLatex(x0, 0.39, "L_{int} = %.6g pb^{-1}" % integrated_luminosity_pb_inv)
    latex.DrawLatex(x0, 0.29, "N_{gen} = %.6g" % n_gen)
    latex.DrawLatex(x0, 0.19, "event weight = %.6g" % normalization_factor)

    if ratio_mode:
        latex.DrawLatex(x0, 0.09, "ratio: data / MC")
    else:
        latex.DrawLatex(x0, 0.09, "MC fill: const AAOgen norm")
    #endif


def make_comparison_legend(h_data, h_mc):
    legend = ROOT.TLegend(0.64, 0.76, 0.92, 0.89)
    legend.SetBorderSize(1)
    legend.SetFillStyle(1001)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetTextSize(0.032)
    legend.AddEntry(h_data, "data", "l")
    legend.AddEntry(h_mc, "AAOgen MC", "l")

    return legend


def draw_sector_comparison_canvas(output_tag, variable_config, data_histograms, mc_histograms, output_pdf, normalization_factor, total_charge_mc, integrated_luminosity_pb_inv, n_gen, log_y):
    if log_y:
        status("Drawing {} sector data-vs-MC output canvas with log y scale.".format(variable_config["key"]))
    else:
        status("Drawing {} sector data-vs-MC output canvas with linear y scale.".format(variable_config["key"]))
    #endif

    ROOT.gStyle.SetOptStat(0)

    if log_y:
        canvas_name = "canvas_{}_comparison_log".format(variable_config["key"])
        canvas_title = "{} {} data vs MC log y".format(output_tag, variable_config["key"])
    else:
        canvas_name = "canvas_{}_comparison_linear".format(variable_config["key"])
        canvas_title = "{} {} data vs MC linear y".format(output_tag, variable_config["key"])
    #endif

    canvas = ROOT.TCanvas(canvas_name, canvas_title, 1600, 900)
    canvas.Divide(4, 2)

    panel_labels = get_panel_labels(variable_config)

    canvas_pad_for_panel = [
        1,
        2,
        3,
        5,
        6,
        7,
        8,
    ]

    y_ranges = get_common_y_ranges_for_comparison(variable_config, data_histograms, mc_histograms, log_y)

    status("{} comparison canvas {} common y-range: [{:.12g}, {:.12g}]".format(
        variable_config["key"],
        variable_config["first_region_name"],
        y_ranges["first"][0],
        y_ranges["first"][1]
    ))
    status("{} comparison canvas {} y-range: [{:.12g}, {:.12g}]".format(
        variable_config["key"],
        variable_config["second_region_name"],
        y_ranges["second"][0],
        y_ranges["second"][1]
    ))

    for i_panel in range(7):
        canvas.cd(canvas_pad_for_panel[i_panel])

        pad = ROOT.gPad
        pad.SetLeftMargin(0.16)
        pad.SetRightMargin(0.06)
        pad.SetTopMargin(0.12)
        pad.SetBottomMargin(0.14)
        pad.SetLogy(log_y)

        h_data = data_histograms[i_panel]
        h_mc = mc_histograms[i_panel]

        y_min, y_max = get_comparison_y_range_for_panel(variable_config, y_ranges, i_panel)

        h_data.SetMaximum(y_max)
        h_data.SetMinimum(y_min)

        h_data.SetTitle("{}  {}".format(output_tag, panel_labels[i_panel]))
        h_data.GetXaxis().SetTitle(variable_config["x_title"])
        h_data.GetYaxis().SetTitle("Counts")
        h_data.GetXaxis().CenterTitle(True)
        h_data.GetYaxis().CenterTitle(True)
        h_data.GetXaxis().SetTitleSize(0.050)
        h_data.GetYaxis().SetTitleSize(0.050)
        h_data.GetXaxis().SetLabelSize(0.045)
        h_data.GetYaxis().SetLabelSize(0.045)
        h_data.GetYaxis().SetTitleOffset(1.45)

        h_data.Draw("HIST")
        h_mc.Draw("HIST SAME")

        legend = make_comparison_legend(h_data, h_mc)
        legend.Draw()
    #endfor

    canvas.cd(4)
    draw_normalization_pad(output_tag, variable_config, total_charge_mc, integrated_luminosity_pb_inv, n_gen, normalization_factor, False)

    ensure_output_directory(output_pdf)

    status("Saving output PDF: {}".format(output_pdf))
    canvas.SaveAs(output_pdf)
    status("Saved output PDF.")


def draw_phi_comparison_canvas(output_tag, variable_config, data_histograms, mc_histograms, output_pdf, normalization_factor, total_charge_mc, integrated_luminosity_pb_inv, n_gen, log_y):
    if log_y:
        status("Drawing {} combined first-region phi data-vs-MC output canvas with log y scale.".format(variable_config["key"]))
    else:
        status("Drawing {} combined first-region phi data-vs-MC output canvas with linear y scale.".format(variable_config["key"]))
    #endif

    ROOT.gStyle.SetOptStat(0)

    if log_y:
        canvas_name = "canvas_{}_comparison_log".format(variable_config["key"])
        canvas_title = "{} {} data vs MC log y".format(output_tag, variable_config["key"])
    else:
        canvas_name = "canvas_{}_comparison_linear".format(variable_config["key"])
        canvas_title = "{} {} data vs MC linear y".format(output_tag, variable_config["key"])
    #endif

    canvas = ROOT.TCanvas(canvas_name, canvas_title, 1600, 500)
    canvas.Divide(3, 1)

    panel_labels = get_panel_labels(variable_config)
    y_ranges = get_common_y_ranges_for_comparison(variable_config, data_histograms, mc_histograms, log_y)

    status("{} comparison canvas {} combined y-range: [{:.12g}, {:.12g}]".format(
        variable_config["key"],
        variable_config["first_region_name"],
        y_ranges["first"][0],
        y_ranges["first"][1]
    ))
    status("{} comparison canvas {} y-range: [{:.12g}, {:.12g}]".format(
        variable_config["key"],
        variable_config["second_region_name"],
        y_ranges["second"][0],
        y_ranges["second"][1]
    ))

    for i_panel in range(2):
        canvas.cd(i_panel + 1)

        pad = ROOT.gPad
        pad.SetLeftMargin(0.16)
        pad.SetRightMargin(0.06)
        pad.SetTopMargin(0.12)
        pad.SetBottomMargin(0.14)
        pad.SetLogy(log_y)

        h_data = data_histograms[i_panel]
        h_mc = mc_histograms[i_panel]

        y_min, y_max = get_comparison_y_range_for_panel(variable_config, y_ranges, i_panel)

        h_data.SetMaximum(y_max)
        h_data.SetMinimum(y_min)

        h_data.SetTitle("{}  {}".format(output_tag, panel_labels[i_panel]))
        h_data.GetXaxis().SetTitle(variable_config["x_title"])
        h_data.GetYaxis().SetTitle("Counts")
        h_data.GetXaxis().CenterTitle(True)
        h_data.GetYaxis().CenterTitle(True)
        h_data.GetXaxis().SetTitleSize(0.050)
        h_data.GetYaxis().SetTitleSize(0.050)
        h_data.GetXaxis().SetLabelSize(0.045)
        h_data.GetYaxis().SetLabelSize(0.045)
        h_data.GetYaxis().SetTitleOffset(1.45)

        h_data.Draw("HIST")
        h_mc.Draw("HIST SAME")

        legend = make_comparison_legend(h_data, h_mc)
        legend.Draw()
    #endfor

    canvas.cd(3)
    draw_normalization_pad(output_tag, variable_config, total_charge_mc, integrated_luminosity_pb_inv, n_gen, normalization_factor, False)

    ensure_output_directory(output_pdf)

    status("Saving output PDF: {}".format(output_pdf))
    canvas.SaveAs(output_pdf)
    status("Saved output PDF.")


def draw_comparison_canvas(output_tag, variable_config, data_histograms, mc_histograms, output_pdf, normalization_factor, total_charge_mc, integrated_luminosity_pb_inv, n_gen, log_y):
    if variable_config["layout"] == "sector_2x4":
        draw_sector_comparison_canvas(
            output_tag,
            variable_config,
            data_histograms,
            mc_histograms,
            output_pdf,
            normalization_factor,
            total_charge_mc,
            integrated_luminosity_pb_inv,
            n_gen,
            log_y
        )

        return
    #endif

    if variable_config["layout"] == "phi_1x3":
        draw_phi_comparison_canvas(
            output_tag,
            variable_config,
            data_histograms,
            mc_histograms,
            output_pdf,
            normalization_factor,
            total_charge_mc,
            integrated_luminosity_pb_inv,
            n_gen,
            log_y
        )

        return
    #endif

    fatal("unknown layout: {}".format(variable_config["layout"]))


def configure_ratio_frame(frame):
    frame.SetStats(False)
    frame.SetMinimum(RATIO_LINEAR_Y_MIN)
    frame.SetMaximum(RATIO_Y_MAX)
    frame.GetXaxis().CenterTitle(True)
    frame.GetYaxis().CenterTitle(True)
    frame.GetXaxis().SetTitleSize(0.050)
    frame.GetYaxis().SetTitleSize(0.050)
    frame.GetXaxis().SetLabelSize(0.045)
    frame.GetYaxis().SetLabelSize(0.045)
    frame.GetYaxis().SetTitleOffset(1.45)


def draw_ratio_title(title_text):
    latex = ROOT.TLatex()
    latex.SetNDC(True)
    latex.SetTextSize(0.055)
    latex.SetTextAlign(22)
    latex.DrawLatex(0.50, 0.94, title_text)


def draw_sector_ratio_canvas(output_tag, variable_config, ratio_graphs, output_pdf, normalization_factor, total_charge_mc, integrated_luminosity_pb_inv, n_gen):
    status("Drawing {} sector data/MC ratio output canvas with linear y scale.".format(variable_config["key"]))

    ROOT.gStyle.SetOptStat(0)

    canvas_name = "canvas_{}_ratio_linear".format(variable_config["key"])
    canvas_title = "{} {} data over MC linear y".format(output_tag, variable_config["key"])

    canvas = ROOT.TCanvas(canvas_name, canvas_title, 1600, 900)
    canvas.Divide(4, 2)

    panel_labels = get_panel_labels(variable_config)

    canvas_pad_for_panel = [
        1,
        2,
        3,
        5,
        6,
        7,
        8,
    ]

    frame_histograms = []
    unity_lines = []

    for i_panel in range(7):
        canvas.cd(canvas_pad_for_panel[i_panel])

        pad = ROOT.gPad
        pad.SetLeftMargin(0.16)
        pad.SetRightMargin(0.06)
        pad.SetTopMargin(0.16)
        pad.SetBottomMargin(0.14)
        pad.SetLogy(False)

        value_min, value_max = get_plot_range_for_panel(variable_config, i_panel)
        n_bins = get_n_bins_for_panel(variable_config, i_panel)

        frame = ROOT.TH1D(
            "frame_{}_panel_{}".format(variable_config["key"], i_panel + 1),
            ";{};data / MC".format(variable_config["x_title"]),
            n_bins,
            value_min,
            value_max
        )

        configure_ratio_frame(frame)
        frame.Draw("AXIS")

        draw_ratio_title("{}  {}".format(output_tag, panel_labels[i_panel]))

        frame_histograms.append(frame)

        unity_line = ROOT.TLine(value_min, 1.0, value_max, 1.0)
        unity_line.SetLineColor(ROOT.kRed + 1)
        unity_line.SetLineStyle(2)
        unity_line.SetLineWidth(2)
        unity_line.Draw("SAME")
        unity_lines.append(unity_line)

        graph = ratio_graphs[i_panel]

        if graph.GetN() > 0:
            graph.SetMarkerSize(0.95)
            graph.SetLineWidth(1)
            graph.Draw("PZ SAME")
        #endif
    #endfor

    canvas.cd(4)
    draw_normalization_pad(output_tag, variable_config, total_charge_mc, integrated_luminosity_pb_inv, n_gen, normalization_factor, True)

    ensure_output_directory(output_pdf)

    status("Saving output PDF: {}".format(output_pdf))
    canvas.SaveAs(output_pdf)
    status("Saved output PDF.")


def draw_phi_ratio_canvas(output_tag, variable_config, ratio_graphs, output_pdf, normalization_factor, total_charge_mc, integrated_luminosity_pb_inv, n_gen):
    status("Drawing {} combined first-region phi data/MC ratio output canvas with linear y scale.".format(variable_config["key"]))

    ROOT.gStyle.SetOptStat(0)

    canvas_name = "canvas_{}_ratio_linear".format(variable_config["key"])
    canvas_title = "{} {} data over MC linear y".format(output_tag, variable_config["key"])

    canvas = ROOT.TCanvas(canvas_name, canvas_title, 1600, 500)
    canvas.Divide(3, 1)

    panel_labels = get_panel_labels(variable_config)
    frame_histograms = []
    unity_lines = []

    for i_panel in range(2):
        canvas.cd(i_panel + 1)

        pad = ROOT.gPad
        pad.SetLeftMargin(0.16)
        pad.SetRightMargin(0.06)
        pad.SetTopMargin(0.16)
        pad.SetBottomMargin(0.14)
        pad.SetLogy(False)

        value_min, value_max = get_plot_range_for_panel(variable_config, i_panel)
        n_bins = get_n_bins_for_panel(variable_config, i_panel)

        frame = ROOT.TH1D(
            "frame_{}_panel_{}".format(variable_config["key"], i_panel + 1),
            ";{};data / MC".format(variable_config["x_title"]),
            n_bins,
            value_min,
            value_max
        )

        configure_ratio_frame(frame)
        frame.Draw("AXIS")

        draw_ratio_title("{}  {}".format(output_tag, panel_labels[i_panel]))

        frame_histograms.append(frame)

        unity_line = ROOT.TLine(value_min, 1.0, value_max, 1.0)
        unity_line.SetLineColor(ROOT.kRed + 1)
        unity_line.SetLineStyle(2)
        unity_line.SetLineWidth(2)
        unity_line.Draw("SAME")
        unity_lines.append(unity_line)

        graph = ratio_graphs[i_panel]

        if graph.GetN() > 0:
            graph.SetMarkerSize(0.95)
            graph.SetLineWidth(1)
            graph.Draw("PZ SAME")
        #endif
    #endfor

    canvas.cd(3)
    draw_normalization_pad(output_tag, variable_config, total_charge_mc, integrated_luminosity_pb_inv, n_gen, normalization_factor, True)

    ensure_output_directory(output_pdf)

    status("Saving output PDF: {}".format(output_pdf))
    canvas.SaveAs(output_pdf)
    status("Saved output PDF.")


def draw_ratio_canvas(output_tag, variable_config, ratio_graphs, output_pdf, normalization_factor, total_charge_mc, integrated_luminosity_pb_inv, n_gen):
    status("{} ratio canvas y-range: [{:.12g}, {:.12g}]".format(
        variable_config["key"],
        RATIO_LINEAR_Y_MIN,
        RATIO_Y_MAX
    ))

    if variable_config["layout"] == "sector_2x4":
        draw_sector_ratio_canvas(
            output_tag,
            variable_config,
            ratio_graphs,
            output_pdf,
            normalization_factor,
            total_charge_mc,
            integrated_luminosity_pb_inv,
            n_gen
        )

        return
    #endif

    if variable_config["layout"] == "phi_1x3":
        draw_phi_ratio_canvas(
            output_tag,
            variable_config,
            ratio_graphs,
            output_pdf,
            normalization_factor,
            total_charge_mc,
            integrated_luminosity_pb_inv,
            n_gen
        )

        return
    #endif

    fatal("unknown layout: {}".format(variable_config["layout"]))


def build_and_save_plots_for_variable(output_tag, variable_config, data_result, mc_result, normalization_factor, total_charge_mc, integrated_luminosity_pb_inv, n_gen):
    output_paths = build_output_paths(output_tag, variable_config)

    status("Building ROOT histograms for variable {}.".format(variable_config["key"]))

    data_histograms = arrays_to_histograms(
        "data",
        variable_config,
        data_result["counts_by_variable"][variable_config["key"]],
        data_result["sumw2_by_variable"][variable_config["key"]]
    )

    mc_histograms = arrays_to_histograms(
        "mc",
        variable_config,
        mc_result["counts_by_variable"][variable_config["key"]],
        mc_result["sumw2_by_variable"][variable_config["key"]]
    )

    ratio_graphs = make_ratio_graphs(variable_config, data_histograms, mc_histograms)

    style_histograms(data_histograms, mc_histograms, ratio_graphs)

    draw_comparison_canvas(
        output_tag,
        variable_config,
        data_histograms,
        mc_histograms,
        output_paths["comparison_log"],
        normalization_factor,
        total_charge_mc,
        integrated_luminosity_pb_inv,
        n_gen,
        True
    )

    draw_comparison_canvas(
        output_tag,
        variable_config,
        data_histograms,
        mc_histograms,
        output_paths["comparison_linear"],
        normalization_factor,
        total_charge_mc,
        integrated_luminosity_pb_inv,
        n_gen,
        False
    )

    draw_ratio_canvas(
        output_tag,
        variable_config,
        ratio_graphs,
        output_paths["ratio_linear"],
        normalization_factor,
        total_charge_mc,
        integrated_luminosity_pb_inv,
        n_gen
    )

    return output_paths


# -----------------------------------------------------------------------------
# Normalization diagnostics
# -----------------------------------------------------------------------------

def sum_counts_for_variable(result, variable_key):
    total = 0.0

    for panel_counts in result["counts_by_variable"][variable_key]:
        for value in panel_counts:
            total += value
        #endfor
    #endfor

    return total


def print_implied_cross_section_diagnostics(data_result, mc_result, integrated_luminosity_pb_inv, n_gen, sigma_used_microbarn, constant_event_weight):
    print("")
    print("AAOgen normalization diagnostic")
    print("-------------------------------")
    print("This diagnostic asks what integrated cross section would be required for the")
    print("constant-weight reconstructed MC prediction to match the data normalization.")
    print("It is a consistency check, not a replacement cross-section measurement.")
    print("Formula: sigma_implied = (N_data / L_int) * (N_gen / N_rec) / 1e6, in microbarn")
    print("")

    for variable_key in get_variable_keys():
        n_data = sum_counts_for_variable(data_result, variable_key)
        n_mc = sum_counts_for_variable(mc_result, variable_key)

        if constant_event_weight > 0.0:
            n_rec_equivalent = n_mc / constant_event_weight
        else:
            n_rec_equivalent = 0.0
        #endif

        if n_data > 0.0 and n_mc > 0.0 and n_rec_equivalent > 0.0 and integrated_luminosity_pb_inv > 0.0:
            data_over_mc = n_data / n_mc
            sigma_implied = (n_data / integrated_luminosity_pb_inv) * (float(n_gen) / n_rec_equivalent) / MICROBARN_TO_PB
            sigma_ratio = sigma_implied / sigma_used_microbarn

            print("{}:".format(variable_key))
            print("  N_data(hist total) = {:.12g}".format(n_data))
            print("  N_MC(pred total)   = {:.12g}".format(n_mc))
            print("  N_rec equivalent   = {:.12g}".format(n_rec_equivalent))
            print("  data / MC          = {:.12g}".format(data_over_mc))
            print("  sigma_used         = {:.12g} microbarn".format(sigma_used_microbarn))
            print("  sigma_implied      = {:.12g} microbarn".format(sigma_implied))
            print("  sigma_implied / sigma_used = {:.12g}".format(sigma_ratio))
        else:
            print("{}:".format(variable_key))
            print("  skipped diagnostic because one of N_data, N_MC, N_rec, or L_int is non-positive.")
        #endif
    #endfor


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compare data and reconstructed AAOgen MC for p1 and p2 variables using a constant integrated-cross-section normalization."
    )

    parser.add_argument("data_root", help="Input data ROOT file containing PhysicsEvents")
    parser.add_argument("reco_mc_root", help="Input reconstructed AAOgen MC ROOT file containing PhysicsEvents")
    parser.add_argument("gen_mc_root", help="Input generated MC ROOT file containing PhysicsEvents")
    parser.add_argument("output_tag", nargs="?", default=DEFAULT_OUTPUT_TAG, help="Output tag / run-period string, e.g. Fa18_Inb. Plots are saved under output/data_mc_normalization/<output_tag>/p1 and /p2.")
    parser.add_argument("--charge-csv", default=GLOBAL_CHARGE_CSV, help="CSV containing run number and accumulated charge")
    parser.add_argument("--charge-to-mc-factor", type=float, default=CHARGE_TO_MC_FACTOR, help="Conversion factor from CSV charge units to mC")
    parser.add_argument("--status-every", type=int, default=DEFAULT_STATUS_EVERY, help="Print loop progress every N entries per worker")
    parser.add_argument("--max-workers", type=int, default=DEFAULT_MAX_WORKERS, help="Maximum number of worker processes. Hard capped at 5.")
    parser.add_argument("--n-gen-override", type=float, default=None, help="Override N_GEN with the true number of thrown/generated MC events.")
    parser.add_argument("--sigma-integrated-microbarn", type=float, default=DEFAULT_AAOGEN_SIGMA_INTEGRATED_MICROBARN, help="AAOgen integrated cross section in microbarns. Default is a temporary 10.6 GeV test value.")

    args = parser.parse_args()

    if args.status_every <= 0:
        fatal("--status-every must be positive")
    #endif

    if args.max_workers <= 0:
        fatal("--max-workers must be positive")
    #endif

    if args.n_gen_override is not None and args.n_gen_override <= 0.0:
        fatal("--n-gen-override must be positive")
    #endif

    if args.sigma_integrated_microbarn <= 0.0:
        fatal("--sigma-integrated-microbarn must be positive")
    #endif

    max_workers = min(args.max_workers, 5)

    output_tag = sanitize_output_tag(args.output_tag)

    start_time = time.time()

    required_data_branches = ["runnum"] + get_required_plot_branches()
    required_mc_branches = get_required_plot_branches()

    status("Starting data/MC normalization script.")
    status("Output tag: {}".format(output_tag))
    status("Data file: {}".format(args.data_root))
    status("Reco MC file: {}".format(args.reco_mc_root))
    status("Generated MC file: {}".format(args.gen_mc_root))
    status("Output base directory: {}".format(os.path.join(OUTPUT_ROOT_DIR, output_tag)))
    status("p1 output directory: {}".format(os.path.join(OUTPUT_ROOT_DIR, output_tag, "p1")))
    status("p2 output directory: {}".format(os.path.join(OUTPUT_ROOT_DIR, output_tag, "p2")))
    status("Maximum worker processes: {}".format(max_workers))
    status("Using constant AAOgen integrated-cross-section MC event weight; MC::Event.weight is not used.")
    status("p1 detector definition: FD = 0 <= p1_theta < 40 deg; CD = 40 <= p1_theta < 70 deg.")
    status("p2 detector definition: FT = 0 <= p2_theta < 5 deg; FD = 5 <= p2_theta < 40 deg.")
    status("Comparison y scales: common across first-region sectors where applicable; second region uses its own scale.")
    status("Ratio y scale: linear only, fixed from 0.2 to 1.5.")
    status("Ratio canvases include explicit subplot titles for sectors and CD/FT panels.")
    status("Data and MC comparison plots are drawn as lines.")
    status("Ratio points use vertical statistical error bars only, with horizontal errors set to zero.")
    status("Variables to plot: {}".format(", ".join(get_variable_keys())))
    status("Special phi behavior: 1x3 canvas with first region combined, second region, and normalization pad.")
    status("First-region bins: theta = {}, p = {}, phi = {}.".format(
        N_BINS_THETA_FIRST_REGION,
        N_BINS_P_FIRST_REGION,
        N_BINS_PHI_FIRST_REGION
    ))
    status("Second-region bins: theta = {}, p = {}, phi = {}.".format(
        N_BINS_THETA_SECOND_REGION,
        N_BINS_P_SECOND_REGION,
        N_BINS_PHI_SECOND_REGION
    ))

    ROOT.gROOT.SetBatch(True)

    n_data_entries = check_input_tree(
        args.data_root,
        TREE_NAME,
        required_data_branches,
        "data"
    )

    n_reco_mc_entries = check_input_tree(
        args.reco_mc_root,
        TREE_NAME,
        required_mc_branches,
        "reconstructed MC"
    )

    n_gen_from_file = get_entry_count(args.gen_mc_root, TREE_NAME, "generated MC")

    if args.n_gen_override is not None:
        n_gen = args.n_gen_override
        status("Using user-provided N_GEN override: {:.12g}".format(n_gen))
        status("Generated MC ROOT file entry count was: {}".format(n_gen_from_file))
    else:
        n_gen = float(n_gen_from_file)
    #endif

    status("Input entry counts:")
    status("  data entries = {}".format(n_data_entries))
    status("  reco MC entries = {}".format(n_reco_mc_entries))
    status("  generated MC entries from file = {}".format(n_gen_from_file))
    status("  generated MC normalization N_GEN = {:.12g}".format(n_gen))

    if n_gen <= 0.0:
        fatal("generated MC normalization N_GEN is zero")
    #endif

    if n_gen < float(n_reco_mc_entries):
        fatal("N_GEN = {:.12g} is smaller than the reconstructed MC entry count = {}. This usually means the generated MC file is not a valid generated-event tree, or you need to use --n-gen-override.".format(
            n_gen,
            n_reco_mc_entries
        ))
    #endif

    charge_by_run = load_charge_map(args.charge_csv)

    data_result = run_data_parallel(
        args.data_root,
        TREE_NAME,
        n_data_entries,
        max_workers,
        args.status_every
    )

    unique_runs = data_result["unique_runs"]
    total_charge_raw = sum_charge_for_runs(unique_runs, charge_by_run)
    total_charge_mc = total_charge_raw * args.charge_to_mc_factor

    integrated_luminosity_pb_inv = RGA_LUMINOSITY_FACTOR_PB_INV_PER_MC * total_charge_mc

    sigma_integrated_microbarn = args.sigma_integrated_microbarn
    sigma_integrated_pb = sigma_integrated_microbarn * MICROBARN_TO_PB
    expected_generated_yield = integrated_luminosity_pb_inv * sigma_integrated_pb
    normalization_factor = expected_generated_yield / float(n_gen)

    status("Computed AAOgen integrated-cross-section MC normalization:")
    status("  Raw charge sum = {:.12g}".format(total_charge_raw))
    status("  Q = {:.12g} mC".format(total_charge_mc))
    status("  L_int = {:.12g} pb^-1".format(integrated_luminosity_pb_inv))
    status("  sigma_int = {:.12g} microbarn = {:.12g} pb".format(sigma_integrated_microbarn, sigma_integrated_pb))
    status("  N_GEN = {:.12g}".format(n_gen))
    status("  expected generated yield = L_int * sigma_int = {:.12g}".format(expected_generated_yield))
    status("  constant reconstructed-MC event weight = expected generated yield / N_GEN = {:.12g}".format(normalization_factor))
    status("  MC::Event.weight is intentionally not used for AAOgen normalization.")

    mc_result = run_mc_parallel(
        args.reco_mc_root,
        TREE_NAME,
        n_reco_mc_entries,
        max_workers,
        args.status_every,
        normalization_factor
    )

    saved_paths_by_variable = {}

    for variable_config in PLOT_VARIABLES:
        output_paths = build_and_save_plots_for_variable(
            output_tag,
            variable_config,
            data_result,
            mc_result,
            normalization_factor,
            total_charge_mc,
            integrated_luminosity_pb_inv,
            n_gen
        )

        saved_paths_by_variable[variable_config["key"]] = output_paths
    #endfor

    print_implied_cross_section_diagnostics(
        data_result,
        mc_result,
        integrated_luminosity_pb_inv,
        n_gen,
        sigma_integrated_microbarn,
        normalization_factor
    )

    elapsed_time = time.time() - start_time

    print("")
    print("Normalization summary")
    print("---------------------")
    print("Output tag: {}".format(output_tag))
    print("Data ROOT file: {}".format(args.data_root))
    print("Reco MC ROOT file: {}".format(args.reco_mc_root))
    print("Gen MC ROOT file: {}".format(args.gen_mc_root))
    print("Charge CSV: {}".format(args.charge_csv))
    print("Maximum workers used: {}".format(max_workers))
    print("Unique data runs found: {}".format(len(unique_runs)))
    print("p1 FD definition: 0 <= p1_theta < 40 deg")
    print("p1 CD definition: 40 <= p1_theta < 70 deg")
    print("p2 FT definition: 0 <= p2_theta < 5 deg")
    print("p2 FD definition: 5 <= p2_theta < 40 deg")
    print("p1_theta FD plotting range: {:.1f} to {:.1f} deg".format(P1_FD_THETA_MIN_DEG, P1_FD_THETA_MAX_DEG))
    print("p1_theta CD plotting range: {:.1f} to {:.1f} deg".format(P1_CD_THETA_MIN_DEG, P1_CD_THETA_MAX_DEG))
    print("p1_p plotting range: {:.1f} to {:.1f} GeV".format(P1_P_MIN_GEV, P1_P_MAX_GEV))
    print("p1_phi plotting range: {:.1f} to {:.1f} deg".format(PHI_MIN_DEG, PHI_MAX_DEG))
    print("p2_theta FT plotting range: {:.1f} to {:.1f} deg".format(P2_FT_THETA_MIN_DEG, P2_FT_THETA_MAX_DEG))
    print("p2_theta FD plotting range: {:.1f} to {:.1f} deg".format(P2_FD_THETA_MIN_DEG, P2_FD_THETA_MAX_DEG))
    print("p2_p plotting range: {:.1f} to {:.1f} GeV".format(P2_P_MIN_GEV, P2_P_MAX_GEV))
    print("p2_phi plotting range: {:.1f} to {:.1f} deg".format(PHI_MIN_DEG, PHI_MAX_DEG))
    print("phi canvas layout: 1x3, first region combined, second region, normalization pad")
    print("Ratio canvases: subplot titles are explicitly drawn above each ratio pad")
    print("First-region theta bins: {}".format(N_BINS_THETA_FIRST_REGION))
    print("Second-region theta bins: {}".format(N_BINS_THETA_SECOND_REGION))
    print("First-region p bins: {}".format(N_BINS_P_FIRST_REGION))
    print("Second-region p bins: {}".format(N_BINS_P_SECOND_REGION))
    print("First-region phi bins: {}".format(N_BINS_PHI_FIRST_REGION))
    print("Second-region phi bins: {}".format(N_BINS_PHI_SECOND_REGION))
    print("Comparison y scales: common across first-region sectors where applicable; second region uses its own scale")
    print("Ratio y range linear: {:.12g} to {:.12g}".format(RATIO_LINEAR_Y_MIN, RATIO_Y_MAX))
    print("Total accumulated charge from CSV raw units: {:.12g}".format(total_charge_raw))
    print("Charge conversion factor to mC: {:.12g}".format(args.charge_to_mc_factor))
    print("Total accumulated charge Q: {:.12g} mC".format(total_charge_mc))
    print("Integrated luminosity: {:.12g} pb^-1".format(integrated_luminosity_pb_inv))
    print("Generated MC entries from file: {}".format(n_gen_from_file))
    print("N_GEN used for normalization: {:.12g}".format(n_gen))
    print("AAOgen integrated cross section used: {:.12g} microbarn".format(sigma_integrated_microbarn))
    print("AAOgen integrated cross section used: {:.12g} pb".format(sigma_integrated_pb))
    print("Expected generated yield L_int * sigma_int: {:.12g}".format(expected_generated_yield))
    print("Constant reconstructed-MC event weight: {:.12g}".format(normalization_factor))
    print("Data p1 FD entries: {}".format(data_result["region_counts"]["p1_first"]))
    print("Data p1 CD entries: {}".format(data_result["region_counts"]["p1_second"]))
    print("Data p1 outside entries: {}".format(data_result["region_counts"]["p1_none"]))
    print("Data p2 FD entries: {}".format(data_result["region_counts"]["p2_first"]))
    print("Data p2 FT entries: {}".format(data_result["region_counts"]["p2_second"]))
    print("Data p2 outside entries: {}".format(data_result["region_counts"]["p2_none"]))
    print("Reco MC p1 FD entries: {}".format(mc_result["region_counts"]["p1_first"]))
    print("Reco MC p1 CD entries: {}".format(mc_result["region_counts"]["p1_second"]))
    print("Reco MC p1 outside entries: {}".format(mc_result["region_counts"]["p1_none"]))
    print("Reco MC p2 FD entries: {}".format(mc_result["region_counts"]["p2_first"]))
    print("Reco MC p2 FT entries: {}".format(mc_result["region_counts"]["p2_second"]))
    print("Reco MC p2 outside entries: {}".format(mc_result["region_counts"]["p2_none"]))
    print("Reco MC any-filled entries: {}".format(mc_result["filled_event_count"]))
    print("Reco MC scaled weight sum over any-filled entries: {:.12g}".format(mc_result["scaled_weight_sum"]))
    print("Reco MC constant event weight: {:.12g}".format(mc_result["constant_event_weight"]))
    print("Output base directory: {}".format(os.path.join(OUTPUT_ROOT_DIR, output_tag)))

    for variable_key in get_variable_keys():
        print("{} comparison log PDF: {}".format(variable_key, saved_paths_by_variable[variable_key]["comparison_log"]))
        print("{} comparison linear PDF: {}".format(variable_key, saved_paths_by_variable[variable_key]["comparison_linear"]))
        print("{} ratio linear PDF: {}".format(variable_key, saved_paths_by_variable[variable_key]["ratio_linear"]))
    #endfor

    print("Elapsed time: {:.2f} seconds".format(elapsed_time))

    status("Finished data/MC normalization script.")


if __name__ == "__main__":
    main()
#endif