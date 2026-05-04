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

CHARGE_TO_MC_FACTOR = 1.0e-6

RGA_LUMINOSITY_FACTOR_PB_INV_PER_MC = 1316.875

OUTPUT_ROOT_DIR = "output/data_mc_normalization"
DEFAULT_OUTPUT_TAG = "default"

FD_THETA_MIN_DEG = 0.0
FD_THETA_MAX_DEG = 40.0

CD_THETA_MIN_DEG = 40.0
CD_THETA_MAX_DEG = 70.0

# Previous value was 70. This is approximately one third of that.
N_BINS_THETA = 23
N_BINS_P = 23
N_BINS_PHI = 24

P1_P_MIN_GEV = 0.3
P1_P_MAX_GEV = 1.3

P1_PHI_MIN_DEG = 0.0
P1_PHI_MAX_DEG = 360.0

DEFAULT_STATUS_EVERY = 250000
DEFAULT_MAX_WORKERS = 5

LOG_Y_MIN = 0.5
LINEAR_Y_MIN = 0.0
Y_PADDING_LINEAR = 1.30
Y_PADDING_LOG = 30.0

RATIO_LINEAR_Y_MIN = 0.0
RATIO_Y_MAX = 2.0

# -----------------------------------------------------------------------------
# Plot variable configuration
# -----------------------------------------------------------------------------

PLOT_VARIABLES = [
    {
        "key": "p1_theta",
        "branch": "p1_theta",
        "title": "p_{1} #theta",
        "x_title": "p_{1} #theta (deg)",
        "n_bins": N_BINS_THETA,
        "fd_min": FD_THETA_MIN_DEG,
        "fd_max": FD_THETA_MAX_DEG,
        "cd_min": CD_THETA_MIN_DEG,
        "cd_max": CD_THETA_MAX_DEG,
        "convert": "rad_to_deg",
    },
    {
        "key": "p1_p",
        "branch": "p1_p",
        "title": "p_{1} momentum",
        "x_title": "p_{1} momentum (GeV)",
        "n_bins": N_BINS_P,
        "fd_min": P1_P_MIN_GEV,
        "fd_max": P1_P_MAX_GEV,
        "cd_min": P1_P_MIN_GEV,
        "cd_max": P1_P_MAX_GEV,
        "convert": "identity",
    },
    {
        "key": "p1_phi",
        "branch": "p1_phi",
        "title": "p_{1} #phi",
        "x_title": "p_{1} #phi (deg)",
        "n_bins": N_BINS_PHI,
        "fd_min": P1_PHI_MIN_DEG,
        "fd_max": P1_PHI_MAX_DEG,
        "cd_min": P1_PHI_MIN_DEG,
        "cd_max": P1_PHI_MAX_DEG,
        "convert": "rad_to_deg_wrapped",
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


def build_output_paths(output_tag, variable_key):
    tag = sanitize_output_tag(output_tag)
    output_dir = os.path.join(OUTPUT_ROOT_DIR, tag)
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


def get_panel_index_from_theta_and_phi(p1_theta_rad, p1_phi_rad):
    p1_theta_deg = math.degrees(p1_theta_rad)

    if p1_theta_deg < FD_THETA_MIN_DEG:
        return -1
    #endif

    if p1_theta_deg < FD_THETA_MAX_DEG:
        p1_phi_deg = math.degrees(p1_phi_rad)
        sector = get_fd_sector_from_phi_deg(p1_phi_deg)

        if sector >= 1 and sector <= 6:
            return sector - 1
        #endif

        return -1
    #endif

    if p1_theta_deg >= CD_THETA_MIN_DEG and p1_theta_deg < CD_THETA_MAX_DEG:
        return 6
    #endif

    return -1


def get_plot_range_for_panel(variable_config, i_panel):
    if i_panel >= 0 and i_panel <= 5:
        return variable_config["fd_min"], variable_config["fd_max"]
    #endif

    if i_panel == 6:
        return variable_config["cd_min"], variable_config["cd_max"]
    #endif

    fatal("invalid panel index in get_plot_range_for_panel: {}".format(i_panel))


def get_value_bin_for_panel(variable_config, i_panel, value):
    value_min, value_max = get_plot_range_for_panel(variable_config, i_panel)
    n_bins = variable_config["n_bins"]

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
    n_bins = variable_config["n_bins"]
    return [[0.0 for _ in range(n_bins)] for _ in range(7)]


def make_empty_sumw2_for_variable(variable_config):
    n_bins = variable_config["n_bins"]
    return [[0.0 for _ in range(n_bins)] for _ in range(7)]


def add_counts_for_variable(total_counts, chunk_counts, variable_config):
    n_bins = variable_config["n_bins"]

    for i_panel in range(7):
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
    tree.SetBranchStatus("p1_phi", 1)
    tree.SetBranchStatus("p1_theta", 1)
    tree.SetBranchStatus("p1_p", 1)


def configure_tree_for_mc(tree):
    tree.SetBranchStatus("*", 0)
    tree.SetBranchStatus("p1_phi", 1)
    tree.SetBranchStatus("p1_theta", 1)
    tree.SetBranchStatus("p1_p", 1)
    tree.SetBranchStatus("weight", 1)


def make_empty_result_maps():
    counts_by_variable = {}
    sumw2_by_variable = {}

    for variable_config in PLOT_VARIABLES:
        key = variable_config["key"]
        counts_by_variable[key] = make_empty_counts_for_variable(variable_config)
        sumw2_by_variable[key] = make_empty_sumw2_for_variable(variable_config)
    #endfor

    return counts_by_variable, sumw2_by_variable


def fill_variable_counts(counts_by_variable, sumw2_by_variable, variable_config, i_panel, raw_value, event_weight):
    key = variable_config["key"]
    value = convert_value(raw_value, variable_config["convert"])
    i_bin = get_value_bin_for_panel(variable_config, i_panel, value)

    if i_bin < 0:
        return False
    #endif

    counts_by_variable[key][i_panel][i_bin] += event_weight
    sumw2_by_variable[key][i_panel][i_bin] += event_weight * event_weight

    return True


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

    n_filled_panel = 0
    n_fd = 0
    n_cd = 0
    n_skipped_panel = 0
    n_total = end_entry - start_entry

    variable_fills = {}
    variable_skips = {}

    for key in get_variable_keys():
        variable_fills[key] = 0
        variable_skips[key] = 0
    #endfor

    print("[{}] data worker {} starting entries {} to {}.".format(
        time.strftime("%Y-%m-%d %H:%M:%S"),
        worker_id,
        start_entry,
        end_entry
    ), flush=True)

    for local_index, i_entry in enumerate(range(start_entry, end_entry), start=1):
        tree.GetEntry(i_entry)

        runnum = int(getattr(tree, "runnum"))
        p1_phi_rad = float(getattr(tree, "p1_phi"))
        p1_theta_rad = float(getattr(tree, "p1_theta"))
        p1_p = float(getattr(tree, "p1_p"))

        unique_runs.add(runnum)

        i_panel = get_panel_index_from_theta_and_phi(p1_theta_rad, p1_phi_rad)

        if i_panel < 0:
            n_skipped_panel += 1
        else:
            n_filled_panel += 1

            if i_panel >= 0 and i_panel <= 5:
                n_fd += 1
            elif i_panel == 6:
                n_cd += 1
            #endif

            raw_values = {
                "p1_theta": p1_theta_rad,
                "p1_p": p1_p,
                "p1_phi": p1_phi_rad,
            }

            for variable_config in PLOT_VARIABLES:
                key = variable_config["key"]
                event_weight = 1.0
                filled = fill_variable_counts(
                    counts_by_variable,
                    sumw2_by_variable,
                    variable_config,
                    i_panel,
                    raw_values[key],
                    event_weight
                )

                if filled:
                    variable_fills[key] += 1
                else:
                    variable_skips[key] += 1
                #endif
            #endfor
        #endif

        if local_index % status_every == 0:
            print("[{}] data worker {} progress: {}/{} chunk entries ({}) | panel-filled: {} | FD(theta<40): {} | CD(theta>=40): {} | unique runs: {}".format(
                time.strftime("%Y-%m-%d %H:%M:%S"),
                worker_id,
                local_index,
                n_total,
                format_percent(local_index, n_total),
                n_filled_panel,
                n_fd,
                n_cd,
                len(unique_runs)
            ), flush=True)
        #endif
    #endfor

    root_file.Close()

    print("[{}] data worker {} finished: panel-filled = {}, FD(theta<40) = {}, CD(theta>=40) = {}, skipped panel = {}, unique runs = {}.".format(
        time.strftime("%Y-%m-%d %H:%M:%S"),
        worker_id,
        n_filled_panel,
        n_fd,
        n_cd,
        n_skipped_panel,
        len(unique_runs)
    ), flush=True)

    return {
        "counts_by_variable": counts_by_variable,
        "sumw2_by_variable": sumw2_by_variable,
        "unique_runs": sorted(unique_runs),
        "n_filled_panel": n_filled_panel,
        "n_fd": n_fd,
        "n_cd": n_cd,
        "n_skipped_panel": n_skipped_panel,
        "variable_fills": variable_fills,
        "variable_skips": variable_skips,
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

    n_filled_panel = 0
    n_fd = 0
    n_cd = 0
    n_skipped_panel = 0
    n_total = end_entry - start_entry
    raw_weight_sum = 0.0
    scaled_weight_sum = 0.0
    max_raw_weight = None

    variable_fills = {}
    variable_skips = {}

    for key in get_variable_keys():
        variable_fills[key] = 0
        variable_skips[key] = 0
    #endfor

    print("[{}] MC worker {} starting entries {} to {}.".format(
        time.strftime("%Y-%m-%d %H:%M:%S"),
        worker_id,
        start_entry,
        end_entry
    ), flush=True)

    for local_index, i_entry in enumerate(range(start_entry, end_entry), start=1):
        tree.GetEntry(i_entry)

        p1_phi_rad = float(getattr(tree, "p1_phi"))
        p1_theta_rad = float(getattr(tree, "p1_theta"))
        p1_p = float(getattr(tree, "p1_p"))
        raw_weight = float(getattr(tree, "weight"))
        event_weight = normalization_factor * raw_weight

        i_panel = get_panel_index_from_theta_and_phi(p1_theta_rad, p1_phi_rad)

        if i_panel < 0:
            n_skipped_panel += 1
        else:
            n_filled_panel += 1
            raw_weight_sum += raw_weight
            scaled_weight_sum += event_weight

            if max_raw_weight is None or raw_weight > max_raw_weight:
                max_raw_weight = raw_weight
            #endif

            if i_panel >= 0 and i_panel <= 5:
                n_fd += 1
            elif i_panel == 6:
                n_cd += 1
            #endif

            raw_values = {
                "p1_theta": p1_theta_rad,
                "p1_p": p1_p,
                "p1_phi": p1_phi_rad,
            }

            for variable_config in PLOT_VARIABLES:
                key = variable_config["key"]
                filled = fill_variable_counts(
                    counts_by_variable,
                    sumw2_by_variable,
                    variable_config,
                    i_panel,
                    raw_values[key],
                    event_weight
                )

                if filled:
                    variable_fills[key] += 1
                else:
                    variable_skips[key] += 1
                #endif
            #endfor
        #endif

        if local_index % status_every == 0:
            print("[{}] MC worker {} progress: {}/{} chunk entries ({}) | panel-filled: {} | FD(theta<40): {} | CD(theta>=40): {} | raw weight sum in accepted panels: {:.12g}".format(
                time.strftime("%Y-%m-%d %H:%M:%S"),
                worker_id,
                local_index,
                n_total,
                format_percent(local_index, n_total),
                n_filled_panel,
                n_fd,
                n_cd,
                raw_weight_sum
            ), flush=True)
        #endif
    #endfor

    root_file.Close()

    if max_raw_weight is None:
        max_raw_weight = 0.0
    #endif

    print("[{}] MC worker {} finished: panel-filled = {}, FD(theta<40) = {}, CD(theta>=40) = {}, skipped panel = {}, raw weight sum = {:.12g}, scaled weight sum = {:.12g}, max raw weight = {:.12g}.".format(
        time.strftime("%Y-%m-%d %H:%M:%S"),
        worker_id,
        n_filled_panel,
        n_fd,
        n_cd,
        n_skipped_panel,
        raw_weight_sum,
        scaled_weight_sum,
        max_raw_weight
    ), flush=True)

    return {
        "counts_by_variable": counts_by_variable,
        "sumw2_by_variable": sumw2_by_variable,
        "n_filled_panel": n_filled_panel,
        "n_fd": n_fd,
        "n_cd": n_cd,
        "n_skipped_panel": n_skipped_panel,
        "raw_weight_sum": raw_weight_sum,
        "scaled_weight_sum": scaled_weight_sum,
        "max_raw_weight": max_raw_weight,
        "variable_fills": variable_fills,
        "variable_skips": variable_skips,
    }


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

    n_filled_panel = 0
    n_fd = 0
    n_cd = 0
    n_skipped_panel = 0

    variable_fills = {}
    variable_skips = {}

    for key in get_variable_keys():
        variable_fills[key] = 0
        variable_skips[key] = 0
    #endfor

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
            variable_skips[key] += result["variable_skips"][key]
        #endfor

        unique_runs.update(result["unique_runs"])
        n_filled_panel += result["n_filled_panel"]
        n_fd += result["n_fd"]
        n_cd += result["n_cd"]
        n_skipped_panel += result["n_skipped_panel"]
    #endfor

    status("Finished parallel data pass.")
    status("Data total: panel-filled = {}, FD(theta<40) = {}, CD(theta>=40) = {}, skipped panel = {}, unique runs = {}.".format(
        n_filled_panel,
        n_fd,
        n_cd,
        n_skipped_panel,
        len(unique_runs)
    ))

    for key in get_variable_keys():
        status("Data variable {}: filled bins entries = {}, skipped range = {}.".format(
            key,
            variable_fills[key],
            variable_skips[key]
        ))
    #endfor

    return {
        "counts_by_variable": total_counts_by_variable,
        "sumw2_by_variable": total_sumw2_by_variable,
        "unique_runs": unique_runs,
        "n_filled_panel": n_filled_panel,
        "n_fd": n_fd,
        "n_cd": n_cd,
        "n_skipped_panel": n_skipped_panel,
        "variable_fills": variable_fills,
        "variable_skips": variable_skips,
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

    n_filled_panel = 0
    n_fd = 0
    n_cd = 0
    n_skipped_panel = 0
    raw_weight_sum = 0.0
    scaled_weight_sum = 0.0
    max_raw_weight = None

    variable_fills = {}
    variable_skips = {}

    for key in get_variable_keys():
        variable_fills[key] = 0
        variable_skips[key] = 0
    #endfor

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
            variable_skips[key] += result["variable_skips"][key]
        #endfor

        n_filled_panel += result["n_filled_panel"]
        n_fd += result["n_fd"]
        n_cd += result["n_cd"]
        n_skipped_panel += result["n_skipped_panel"]
        raw_weight_sum += result["raw_weight_sum"]
        scaled_weight_sum += result["scaled_weight_sum"]

        if max_raw_weight is None or result["max_raw_weight"] > max_raw_weight:
            max_raw_weight = result["max_raw_weight"]
        #endif
    #endfor

    if max_raw_weight is None:
        max_raw_weight = 0.0
    #endif

    status("Finished parallel reconstructed MC pass.")
    status("Reco MC total: panel-filled = {}, FD(theta<40) = {}, CD(theta>=40) = {}, skipped panel = {}.".format(
        n_filled_panel,
        n_fd,
        n_cd,
        n_skipped_panel
    ))
    status("Reco MC weight summary over accepted panels: raw weight sum = {:.12g}, scaled weight sum = {:.12g}, max raw weight = {:.12g}.".format(
        raw_weight_sum,
        scaled_weight_sum,
        max_raw_weight
    ))

    for key in get_variable_keys():
        status("Reco MC variable {}: filled bins entries = {}, skipped range = {}.".format(
            key,
            variable_fills[key],
            variable_skips[key]
        ))
    #endfor

    return {
        "counts_by_variable": total_counts_by_variable,
        "sumw2_by_variable": total_sumw2_by_variable,
        "n_filled_panel": n_filled_panel,
        "n_fd": n_fd,
        "n_cd": n_cd,
        "n_skipped_panel": n_skipped_panel,
        "raw_weight_sum": raw_weight_sum,
        "scaled_weight_sum": scaled_weight_sum,
        "max_raw_weight": max_raw_weight,
        "variable_fills": variable_fills,
        "variable_skips": variable_skips,
    }


# -----------------------------------------------------------------------------
# Histogram and plotting helpers
# -----------------------------------------------------------------------------

def arrays_to_histograms(prefix, variable_config, counts, sumw2):
    histograms = []

    panel_names = [
        "FD sector 1",
        "FD sector 2",
        "FD sector 3",
        "FD sector 4",
        "FD sector 5",
        "FD sector 6",
        "CD",
    ]

    for i_panel in range(7):
        value_min, value_max = get_plot_range_for_panel(variable_config, i_panel)

        hist_name = "{}_{}_panel_{}".format(prefix, variable_config["key"], i_panel + 1)
        hist_title = "{};{};Counts".format(panel_names[i_panel], variable_config["x_title"])

        hist = ROOT.TH1D(
            hist_name,
            hist_title,
            variable_config["n_bins"],
            value_min,
            value_max
        )

        hist.Sumw2()

        for i_bin in range(variable_config["n_bins"]):
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

    for i_panel in range(7):
        value_min, value_max = get_plot_range_for_panel(variable_config, i_panel)
        n_bins = variable_config["n_bins"]
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
        graph.SetMarkerSize(0.65)
        graph.SetLineWidth(2)
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


def get_common_y_ranges_for_comparison(data_histograms, mc_histograms, log_y):
    fd_histograms = []
    cd_histograms = []

    for i_panel in range(6):
        fd_histograms.append(data_histograms[i_panel])
        fd_histograms.append(mc_histograms[i_panel])
    #endfor

    cd_histograms.append(data_histograms[6])
    cd_histograms.append(mc_histograms[6])

    fd_y_min, fd_y_max = get_common_y_range_for_hist_list(fd_histograms, log_y)
    cd_y_min, cd_y_max = get_common_y_range_for_hist_list(cd_histograms, log_y)

    return fd_y_min, fd_y_max, cd_y_min, cd_y_max


def draw_normalization_pad(output_tag, variable_config, total_charge_mc, integrated_luminosity_pb_inv, n_gen, normalization_factor, log_y, ratio_mode):
    pad4 = ROOT.gPad
    pad4.Clear()
    pad4.SetFillColor(ROOT.kWhite)

    latex = ROOT.TLatex()
    latex.SetNDC(True)
    latex.SetTextSize(0.034)

    if ratio_mode:
        title = "Ratio normalization"
    else:
        title = "MC normalization"
    #endif

    if log_y:
        scale_label = "log y"
    else:
        scale_label = "linear y"
    #endif

    latex.DrawLatex(0.10, 0.90, output_tag)
    latex.DrawLatex(0.10, 0.81, variable_config["title"])
    latex.DrawLatex(0.10, 0.72, title)
    latex.DrawLatex(0.10, 0.63, "scale: {}".format(scale_label))
    latex.DrawLatex(0.10, 0.54, "Q = %.6g mC" % total_charge_mc)
    latex.DrawLatex(0.10, 0.45, "L_{int} = %.6g pb^{-1}" % integrated_luminosity_pb_inv)
    latex.DrawLatex(0.10, 0.36, "N_{gen} = %.6g" % n_gen)
    latex.DrawLatex(0.10, 0.27, "scale = L_{int}/N_{gen}")
    latex.DrawLatex(0.10, 0.18, "scale = %.6g pb^{-1}" % normalization_factor)

    if ratio_mode:
        latex.DrawLatex(0.10, 0.09, "ratio: data / MC")
    else:
        latex.DrawLatex(0.10, 0.09, "MC fill: scale #times weight")
    #endif


def draw_comparison_canvas(output_tag, variable_config, data_histograms, mc_histograms, output_pdf, normalization_factor, total_charge_mc, integrated_luminosity_pb_inv, n_gen, log_y):
    if log_y:
        status("Drawing {} data-vs-MC output canvas with log y scale.".format(variable_config["key"]))
    else:
        status("Drawing {} data-vs-MC output canvas with linear y scale.".format(variable_config["key"]))
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

    panel_labels = [
        "FD sector 1",
        "FD sector 2",
        "FD sector 3",
        "FD sector 4",
        "FD sector 5",
        "FD sector 6",
        "Central detector",
    ]

    canvas_pad_for_panel = [
        1,
        2,
        3,
        5,
        6,
        7,
        8,
    ]

    fd_y_min, fd_y_max, cd_y_min, cd_y_max = get_common_y_ranges_for_comparison(data_histograms, mc_histograms, log_y)

    status("{} comparison canvas FD common y-range: [{:.12g}, {:.12g}]".format(variable_config["key"], fd_y_min, fd_y_max))
    status("{} comparison canvas CD y-range: [{:.12g}, {:.12g}]".format(variable_config["key"], cd_y_min, cd_y_max))

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

        if i_panel >= 0 and i_panel <= 5:
            h_data.SetMaximum(fd_y_max)
            h_data.SetMinimum(fd_y_min)
        else:
            h_data.SetMaximum(cd_y_max)
            h_data.SetMinimum(cd_y_min)
        #endif

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

        legend = ROOT.TLegend(0.64, 0.76, 0.92, 0.89)
        legend.SetBorderSize(1)
        legend.SetFillStyle(1001)
        legend.SetFillColor(ROOT.kWhite)
        legend.SetTextSize(0.032)
        legend.AddEntry(h_data, "data", "l")
        legend.AddEntry(h_mc, "MC weighted", "l")
        legend.Draw()
    #endfor

    canvas.cd(4)
    draw_normalization_pad(output_tag, variable_config, total_charge_mc, integrated_luminosity_pb_inv, n_gen, normalization_factor, log_y, False)

    ensure_output_directory(output_pdf)

    status("Saving output PDF: {}".format(output_pdf))
    canvas.SaveAs(output_pdf)
    status("Saved output PDF.")


def draw_ratio_canvas(output_tag, variable_config, ratio_graphs, output_pdf, normalization_factor, total_charge_mc, integrated_luminosity_pb_inv, n_gen):
    status("Drawing {} data/MC ratio output canvas with linear y scale.".format(variable_config["key"]))

    ROOT.gStyle.SetOptStat(0)

    canvas_name = "canvas_{}_ratio_linear".format(variable_config["key"])
    canvas_title = "{} {} data over MC linear y".format(output_tag, variable_config["key"])

    canvas = ROOT.TCanvas(canvas_name, canvas_title, 1600, 900)
    canvas.Divide(4, 2)

    panel_labels = [
        "FD sector 1",
        "FD sector 2",
        "FD sector 3",
        "FD sector 4",
        "FD sector 5",
        "FD sector 6",
        "Central detector",
    ]

    canvas_pad_for_panel = [
        1,
        2,
        3,
        5,
        6,
        7,
        8,
    ]

    status("{} ratio canvas y-range: [{:.12g}, {:.12g}]".format(
        variable_config["key"],
        RATIO_LINEAR_Y_MIN,
        RATIO_Y_MAX
    ))

    frame_histograms = []

    for i_panel in range(7):
        canvas.cd(canvas_pad_for_panel[i_panel])

        pad = ROOT.gPad
        pad.SetLeftMargin(0.16)
        pad.SetRightMargin(0.06)
        pad.SetTopMargin(0.12)
        pad.SetBottomMargin(0.14)
        pad.SetLogy(False)

        value_min, value_max = get_plot_range_for_panel(variable_config, i_panel)

        frame = ROOT.TH1D(
            "frame_{}_panel_{}".format(variable_config["key"], i_panel + 1),
            "{}  {};{};data / MC".format(output_tag, panel_labels[i_panel], variable_config["x_title"]),
            variable_config["n_bins"],
            value_min,
            value_max
        )

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
        frame.Draw("AXIS")

        frame_histograms.append(frame)

        graph = ratio_graphs[i_panel]

        if graph.GetN() > 0:
            graph.Draw("PZ SAME")
        #endif

        line = ROOT.TLine(value_min, 1.0, value_max, 1.0)
        line.SetLineColor(ROOT.kRed + 1)
        line.SetLineStyle(2)
        line.SetLineWidth(2)
        line.Draw("SAME")
    #endfor

    canvas.cd(4)
    draw_normalization_pad(output_tag, variable_config, total_charge_mc, integrated_luminosity_pb_inv, n_gen, normalization_factor, False, True)

    ensure_output_directory(output_pdf)

    status("Saving output PDF: {}".format(output_pdf))
    canvas.SaveAs(output_pdf)
    status("Saved output PDF.")


def build_and_save_plots_for_variable(output_tag, variable_config, data_result, mc_result, normalization_factor, total_charge_mc, integrated_luminosity_pb_inv, n_gen):
    output_paths = build_output_paths(output_tag, variable_config["key"])

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
# Main
# -----------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compare data and reconstructed MC for p1_theta, p1_p, and p1_phi using theta-defined FD/CD regions with event-by-event MC::Event.weight normalization."
    )

    parser.add_argument("data_root", help="Input data ROOT file containing PhysicsEvents")
    parser.add_argument("reco_mc_root", help="Input reconstructed MC ROOT file containing PhysicsEvents and weight")
    parser.add_argument("gen_mc_root", help="Input generated MC ROOT file containing PhysicsEvents")
    parser.add_argument("output_tag", nargs="?", default=DEFAULT_OUTPUT_TAG, help="Output tag / run-period string, e.g. Fa18_Inb. Plots are saved under output/data_mc_normalization/<output_tag>/.")
    parser.add_argument("--charge-csv", default=GLOBAL_CHARGE_CSV, help="CSV containing run number and accumulated charge")
    parser.add_argument("--charge-to-mc-factor", type=float, default=CHARGE_TO_MC_FACTOR, help="Conversion factor from CSV charge units to mC")
    parser.add_argument("--status-every", type=int, default=DEFAULT_STATUS_EVERY, help="Print loop progress every N entries per worker")
    parser.add_argument("--max-workers", type=int, default=DEFAULT_MAX_WORKERS, help="Maximum number of worker processes. Hard capped at 5.")
    parser.add_argument("--n-gen-override", type=float, default=None, help="Override N_GEN with the true number of thrown/generated MC events.")

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

    max_workers = min(args.max_workers, 5)

    output_tag = sanitize_output_tag(args.output_tag)

    start_time = time.time()

    status("Starting data/MC normalization script.")
    status("Output tag: {}".format(output_tag))
    status("Data file: {}".format(args.data_root))
    status("Reco MC file: {}".format(args.reco_mc_root))
    status("Generated MC file: {}".format(args.gen_mc_root))
    status("Output directory: {}".format(os.path.join(OUTPUT_ROOT_DIR, output_tag)))
    status("Maximum worker processes: {}".format(max_workers))
    status("Using event-by-event reconstructed MC branch: weight")
    status("Detector definition override: FD = p1_theta < 40 deg; CD = p1_theta >= 40 deg.")
    status("Comparison y scales: common across FD sectors; CD uses its own scale.")
    status("Ratio y scale: linear only, fixed from 0 to 2.")
    status("Data and MC comparison plots are drawn as lines.")
    status("Ratio points use vertical statistical error bars only, with horizontal errors set to zero.")
    status("Variables to plot: {}".format(", ".join(get_variable_keys())))

    ROOT.gROOT.SetBatch(True)

    n_data_entries = check_input_tree(
        args.data_root,
        "PhysicsEvents",
        ["runnum", "p1_phi", "p1_theta", "p1_p"],
        "data"
    )

    n_reco_mc_entries = check_input_tree(
        args.reco_mc_root,
        "PhysicsEvents",
        ["p1_phi", "p1_theta", "p1_p", "weight"],
        "reconstructed MC"
    )

    n_gen_from_file = get_entry_count(args.gen_mc_root, "PhysicsEvents", "generated MC")

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
        "PhysicsEvents",
        n_data_entries,
        max_workers,
        args.status_every
    )

    unique_runs = data_result["unique_runs"]
    total_charge_raw = sum_charge_for_runs(unique_runs, charge_by_run)
    total_charge_mc = total_charge_raw * args.charge_to_mc_factor

    integrated_luminosity_pb_inv = RGA_LUMINOSITY_FACTOR_PB_INV_PER_MC * total_charge_mc

    normalization_factor = integrated_luminosity_pb_inv / float(n_gen)

    status("Computed event-by-event MC normalization:")
    status("  Raw charge sum = {:.12g}".format(total_charge_raw))
    status("  Q = {:.12g} mC".format(total_charge_mc))
    status("  L_int = {:.12g} pb^-1".format(integrated_luminosity_pb_inv))
    status("  N_GEN = {:.12g}".format(n_gen))
    status("  normalization factor L_int / N_GEN = {:.12g} pb^-1".format(normalization_factor))
    status("  MC event fill weight = (L_int / N_GEN) * weight")

    mc_result = run_mc_parallel(
        args.reco_mc_root,
        "PhysicsEvents",
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
    print("FD definition: p1_theta < 40 deg")
    print("CD definition: p1_theta >= 40 deg")
    print("p1_theta FD plotting range: {:.1f} to {:.1f} deg".format(FD_THETA_MIN_DEG, FD_THETA_MAX_DEG))
    print("p1_theta CD plotting range: {:.1f} to {:.1f} deg".format(CD_THETA_MIN_DEG, CD_THETA_MAX_DEG))
    print("p1_p plotting range: {:.1f} to {:.1f} GeV".format(P1_P_MIN_GEV, P1_P_MAX_GEV))
    print("p1_phi plotting range: {:.1f} to {:.1f} deg".format(P1_PHI_MIN_DEG, P1_PHI_MAX_DEG))
    print("Number of p1_theta bins: {}".format(N_BINS_THETA))
    print("Number of p1_p bins: {}".format(N_BINS_P))
    print("Number of p1_phi bins: {}".format(N_BINS_PHI))
    print("Comparison y scales: common across FD sectors; CD uses its own scale")
    print("Ratio y range linear: {:.12g} to {:.12g}".format(RATIO_LINEAR_Y_MIN, RATIO_Y_MAX))
    print("Total accumulated charge from CSV raw units: {:.12g}".format(total_charge_raw))
    print("Charge conversion factor to mC: {:.12g}".format(args.charge_to_mc_factor))
    print("Total accumulated charge Q: {:.12g} mC".format(total_charge_mc))
    print("Integrated luminosity: {:.12g} pb^-1".format(integrated_luminosity_pb_inv))
    print("Generated MC entries from file: {}".format(n_gen_from_file))
    print("N_GEN used for normalization: {:.12g}".format(n_gen))
    print("Normalization factor L_int / N_GEN: {:.12g} pb^-1".format(normalization_factor))
    print("Data panel entries filled: {}".format(data_result["n_filled_panel"]))
    print("Data FD panel entries filled: {}".format(data_result["n_fd"]))
    print("Data CD panel entries filled: {}".format(data_result["n_cd"]))
    print("Reco MC panel entries filled: {}".format(mc_result["n_filled_panel"]))
    print("Reco MC FD panel entries filled: {}".format(mc_result["n_fd"]))
    print("Reco MC CD panel entries filled: {}".format(mc_result["n_cd"]))
    print("Reco MC raw weight sum over accepted panels: {:.12g}".format(mc_result["raw_weight_sum"]))
    print("Reco MC scaled weight sum over accepted panels: {:.12g}".format(mc_result["scaled_weight_sum"]))
    print("Reco MC maximum raw weight over accepted panels: {:.12g}".format(mc_result["max_raw_weight"]))
    print("Output directory: {}".format(os.path.join(OUTPUT_ROOT_DIR, output_tag)))

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