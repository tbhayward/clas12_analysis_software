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

# The global.csv charge values look like nC in the examples you showed.
# Since the normalization formula wants Q in mC, the default conversion is:
#
#   Q(mC) = Q(nC) / 1.0e6
#
# If your CSV is already in mC, change this to 1.0 or use:
#
#   --charge-to-mc-factor 1.0
#
CHARGE_TO_MC_FACTOR = 1.0e-6

# RGA integrated luminosity:
#
#   L_int = 1316.875 * Q(mC) pb^{-1}
#
RGA_LUMINOSITY_FACTOR_PB_INV_PER_MC = 1316.875

OUTPUT_ROOT_DIR = "output/data_mc_normalization"
DEFAULT_OUTPUT_TAG = "default"

# Panels 1-6 are FD sectors, defined only by p1_theta < 40 degrees.
# Panel 7 is CD, defined only by p1_theta >= 40 degrees.
FD_THETA_MIN_DEG = 0.0
FD_THETA_MAX_DEG = 40.0

CD_THETA_MIN_DEG = 40.0
CD_THETA_MAX_DEG = 70.0

N_BINS_THETA = 70

DEFAULT_STATUS_EVERY = 250000
DEFAULT_MAX_WORKERS = 5

LOG_Y_MIN = 0.5
LINEAR_Y_MIN = 0.0

# Ratio plots are requested to be limited from 0 to 2.
# A log-y plot cannot have y_min = 0, so the log-y ratio canvas uses this
# small positive lower bound while still using y_max = 2.
RATIO_LOG_Y_MIN = 1.0e-4
RATIO_LINEAR_Y_MIN = 0.0
RATIO_Y_MAX = 2.0

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


def build_output_paths(output_tag):
    tag = sanitize_output_tag(output_tag)
    output_dir = os.path.join(OUTPUT_ROOT_DIR, tag)
    base = os.path.join(output_dir, "output")

    return {
        "output_dir": output_dir,
        "comparison_log": base + "_log.pdf",
        "comparison_linear": base + "_linear.pdf",
        "ratio_log": base + "_ratio_log.pdf",
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


def get_theta_bin_for_panel(i_panel, p1_theta_rad):
    p1_theta_deg = math.degrees(p1_theta_rad)

    if i_panel >= 0 and i_panel <= 5:
        theta_min = FD_THETA_MIN_DEG
        theta_max = FD_THETA_MAX_DEG
    elif i_panel == 6:
        theta_min = CD_THETA_MIN_DEG
        theta_max = CD_THETA_MAX_DEG
    else:
        return -1
    #endif

    if p1_theta_deg < theta_min or p1_theta_deg >= theta_max:
        return -1
    #endif

    bin_width = (theta_max - theta_min) / float(N_BINS_THETA)
    i_bin = int((p1_theta_deg - theta_min) / bin_width)

    if i_bin < 0 or i_bin >= N_BINS_THETA:
        return -1
    #endif

    return i_bin


def get_theta_range_for_panel(i_panel):
    if i_panel >= 0 and i_panel <= 5:
        return FD_THETA_MIN_DEG, FD_THETA_MAX_DEG
    #endif

    if i_panel == 6:
        return CD_THETA_MIN_DEG, CD_THETA_MAX_DEG
    #endif

    fatal("invalid panel index in get_theta_range_for_panel: {}".format(i_panel))


def make_empty_counts():
    return [[0.0 for _ in range(N_BINS_THETA)] for _ in range(7)]


def make_empty_sumw2():
    return [[0.0 for _ in range(N_BINS_THETA)] for _ in range(7)]


def add_counts(total_counts, chunk_counts):
    for i_panel in range(7):
        for i_bin in range(N_BINS_THETA):
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


def configure_tree_for_mc(tree):
    tree.SetBranchStatus("*", 0)
    tree.SetBranchStatus("p1_phi", 1)
    tree.SetBranchStatus("p1_theta", 1)
    tree.SetBranchStatus("weight", 1)


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

    counts = make_empty_counts()
    sumw2 = make_empty_sumw2()
    unique_runs = set()

    n_filled = 0
    n_fd = 0
    n_cd = 0
    n_skipped_panel = 0
    n_skipped_theta = 0
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
        p1_phi_rad = float(getattr(tree, "p1_phi"))
        p1_theta_rad = float(getattr(tree, "p1_theta"))

        unique_runs.add(runnum)

        i_panel = get_panel_index_from_theta_and_phi(p1_theta_rad, p1_phi_rad)

        if i_panel < 0:
            n_skipped_panel += 1
        else:
            i_bin = get_theta_bin_for_panel(i_panel, p1_theta_rad)

            if i_bin < 0:
                n_skipped_theta += 1
            else:
                event_weight = 1.0
                counts[i_panel][i_bin] += event_weight
                sumw2[i_panel][i_bin] += event_weight * event_weight
                n_filled += 1

                if i_panel >= 0 and i_panel <= 5:
                    n_fd += 1
                elif i_panel == 6:
                    n_cd += 1
                #endif
            #endif
        #endif

        if local_index % status_every == 0:
            print("[{}] data worker {} progress: {}/{} chunk entries ({}) | filled: {} | FD(theta<40): {} | CD(theta>=40): {} | unique runs: {}".format(
                time.strftime("%Y-%m-%d %H:%M:%S"),
                worker_id,
                local_index,
                n_total,
                format_percent(local_index, n_total),
                n_filled,
                n_fd,
                n_cd,
                len(unique_runs)
            ), flush=True)
        #endif
    #endfor

    root_file.Close()

    print("[{}] data worker {} finished: filled = {}, FD(theta<40) = {}, CD(theta>=40) = {}, skipped panel = {}, skipped theta = {}, unique runs = {}.".format(
        time.strftime("%Y-%m-%d %H:%M:%S"),
        worker_id,
        n_filled,
        n_fd,
        n_cd,
        n_skipped_panel,
        n_skipped_theta,
        len(unique_runs)
    ), flush=True)

    return {
        "counts": counts,
        "sumw2": sumw2,
        "unique_runs": sorted(unique_runs),
        "n_filled": n_filled,
        "n_fd": n_fd,
        "n_cd": n_cd,
        "n_skipped_panel": n_skipped_panel,
        "n_skipped_theta": n_skipped_theta,
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

    counts = make_empty_counts()
    sumw2 = make_empty_sumw2()

    n_filled = 0
    n_fd = 0
    n_cd = 0
    n_skipped_panel = 0
    n_skipped_theta = 0
    n_total = end_entry - start_entry
    raw_weight_sum = 0.0
    scaled_weight_sum = 0.0
    max_raw_weight = None

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
        raw_weight = float(getattr(tree, "weight"))
        event_weight = normalization_factor * raw_weight

        i_panel = get_panel_index_from_theta_and_phi(p1_theta_rad, p1_phi_rad)

        if i_panel < 0:
            n_skipped_panel += 1
        else:
            i_bin = get_theta_bin_for_panel(i_panel, p1_theta_rad)

            if i_bin < 0:
                n_skipped_theta += 1
            else:
                counts[i_panel][i_bin] += event_weight
                sumw2[i_panel][i_bin] += event_weight * event_weight
                raw_weight_sum += raw_weight
                scaled_weight_sum += event_weight

                if max_raw_weight is None or raw_weight > max_raw_weight:
                    max_raw_weight = raw_weight
                #endif

                n_filled += 1

                if i_panel >= 0 and i_panel <= 5:
                    n_fd += 1
                elif i_panel == 6:
                    n_cd += 1
                #endif
            #endif
        #endif

        if local_index % status_every == 0:
            print("[{}] MC worker {} progress: {}/{} chunk entries ({}) | filled: {} | FD(theta<40): {} | CD(theta>=40): {} | raw weight sum in filled bins: {:.12g}".format(
                time.strftime("%Y-%m-%d %H:%M:%S"),
                worker_id,
                local_index,
                n_total,
                format_percent(local_index, n_total),
                n_filled,
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

    print("[{}] MC worker {} finished: filled = {}, FD(theta<40) = {}, CD(theta>=40) = {}, skipped panel = {}, skipped theta = {}, raw weight sum = {:.12g}, scaled weight sum = {:.12g}, max raw weight = {:.12g}.".format(
        time.strftime("%Y-%m-%d %H:%M:%S"),
        worker_id,
        n_filled,
        n_fd,
        n_cd,
        n_skipped_panel,
        n_skipped_theta,
        raw_weight_sum,
        scaled_weight_sum,
        max_raw_weight
    ), flush=True)

    return {
        "counts": counts,
        "sumw2": sumw2,
        "n_filled": n_filled,
        "n_fd": n_fd,
        "n_cd": n_cd,
        "n_skipped_panel": n_skipped_panel,
        "n_skipped_theta": n_skipped_theta,
        "raw_weight_sum": raw_weight_sum,
        "scaled_weight_sum": scaled_weight_sum,
        "max_raw_weight": max_raw_weight,
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

    total_counts = make_empty_counts()
    total_sumw2 = make_empty_sumw2()
    unique_runs = set()

    n_filled = 0
    n_fd = 0
    n_cd = 0
    n_skipped_panel = 0
    n_skipped_theta = 0

    with mp.Pool(processes=len(tasks)) as pool:
        results = pool.map(data_worker, tasks)
    #endwith

    for result in results:
        add_counts(total_counts, result["counts"])
        add_counts(total_sumw2, result["sumw2"])
        unique_runs.update(result["unique_runs"])
        n_filled += result["n_filled"]
        n_fd += result["n_fd"]
        n_cd += result["n_cd"]
        n_skipped_panel += result["n_skipped_panel"]
        n_skipped_theta += result["n_skipped_theta"]
    #endfor

    status("Finished parallel data pass.")
    status("Data total: filled = {}, FD(theta<40) = {}, CD(theta>=40) = {}, skipped panel = {}, skipped theta = {}, unique runs = {}.".format(
        n_filled,
        n_fd,
        n_cd,
        n_skipped_panel,
        n_skipped_theta,
        len(unique_runs)
    ))

    return {
        "counts": total_counts,
        "sumw2": total_sumw2,
        "unique_runs": unique_runs,
        "n_filled": n_filled,
        "n_fd": n_fd,
        "n_cd": n_cd,
        "n_skipped_panel": n_skipped_panel,
        "n_skipped_theta": n_skipped_theta,
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

    total_counts = make_empty_counts()
    total_sumw2 = make_empty_sumw2()

    n_filled = 0
    n_fd = 0
    n_cd = 0
    n_skipped_panel = 0
    n_skipped_theta = 0
    raw_weight_sum = 0.0
    scaled_weight_sum = 0.0
    max_raw_weight = None

    with mp.Pool(processes=len(tasks)) as pool:
        results = pool.map(mc_worker, tasks)
    #endwith

    for result in results:
        add_counts(total_counts, result["counts"])
        add_counts(total_sumw2, result["sumw2"])
        n_filled += result["n_filled"]
        n_fd += result["n_fd"]
        n_cd += result["n_cd"]
        n_skipped_panel += result["n_skipped_panel"]
        n_skipped_theta += result["n_skipped_theta"]
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
    status("Reco MC total: filled = {}, FD(theta<40) = {}, CD(theta>=40) = {}, skipped panel = {}, skipped theta = {}.".format(
        n_filled,
        n_fd,
        n_cd,
        n_skipped_panel,
        n_skipped_theta
    ))
    status("Reco MC weight summary over filled entries: raw weight sum = {:.12g}, scaled weight sum = {:.12g}, max raw weight = {:.12g}.".format(
        raw_weight_sum,
        scaled_weight_sum,
        max_raw_weight
    ))

    return {
        "counts": total_counts,
        "sumw2": total_sumw2,
        "n_filled": n_filled,
        "n_fd": n_fd,
        "n_cd": n_cd,
        "n_skipped_panel": n_skipped_panel,
        "n_skipped_theta": n_skipped_theta,
        "raw_weight_sum": raw_weight_sum,
        "scaled_weight_sum": scaled_weight_sum,
        "max_raw_weight": max_raw_weight,
    }


# -----------------------------------------------------------------------------
# Histogram and plotting helpers
# -----------------------------------------------------------------------------

def arrays_to_histograms(prefix, counts, sumw2):
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
        theta_min, theta_max = get_theta_range_for_panel(i_panel)

        hist_name = "{}_panel_{}".format(prefix, i_panel + 1)
        hist_title = "{};p_{{1}} #theta (deg);Counts".format(panel_names[i_panel])

        hist = ROOT.TH1D(
            hist_name,
            hist_title,
            N_BINS_THETA,
            theta_min,
            theta_max
        )

        hist.Sumw2()

        for i_bin in range(N_BINS_THETA):
            root_bin = i_bin + 1
            content = counts[i_panel][i_bin]
            error = math.sqrt(sumw2[i_panel][i_bin])

            hist.SetBinContent(root_bin, content)
            hist.SetBinError(root_bin, error)
        #endfor

        histograms.append(hist)
    #endfor

    return histograms


def make_ratio_histograms(data_histograms, mc_histograms):
    ratio_histograms = []

    for i_panel in range(7):
        theta_min, theta_max = get_theta_range_for_panel(i_panel)

        ratio_hist = ROOT.TH1D(
            "ratio_panel_{}".format(i_panel + 1),
            "ratio_panel_{};p_{{1}} #theta (deg);data / MC".format(i_panel + 1),
            N_BINS_THETA,
            theta_min,
            theta_max
        )

        ratio_hist.Sumw2()

        h_data = data_histograms[i_panel]
        h_mc = mc_histograms[i_panel]

        for i_bin in range(1, N_BINS_THETA + 1):
            data_content = h_data.GetBinContent(i_bin)
            data_error = h_data.GetBinError(i_bin)
            mc_content = h_mc.GetBinContent(i_bin)
            mc_error = h_mc.GetBinError(i_bin)

            if data_content > 0.0 and mc_content > 0.0:
                ratio = data_content / mc_content
                rel_data = data_error / data_content
                rel_mc = mc_error / mc_content
                ratio_error = ratio * math.sqrt(rel_data * rel_data + rel_mc * rel_mc)

                ratio_hist.SetBinContent(i_bin, ratio)
                ratio_hist.SetBinError(i_bin, ratio_error)
            else:
                ratio_hist.SetBinContent(i_bin, 0.0)
                ratio_hist.SetBinError(i_bin, 0.0)
            #endif
        #endfor

        ratio_histograms.append(ratio_hist)
    #endfor

    return ratio_histograms


def style_histograms(data_histograms, mc_histograms, ratio_histograms):
    status("Styling histograms.")

    for hist in data_histograms:
        hist.SetLineColor(ROOT.kBlue)
        hist.SetMarkerColor(ROOT.kBlue)
        hist.SetMarkerStyle(20)
        hist.SetMarkerSize(0.65)
        hist.SetLineWidth(2)
        hist.SetStats(False)
    #endfor

    for hist in mc_histograms:
        hist.SetLineColor(ROOT.kRed)
        hist.SetMarkerColor(ROOT.kRed)
        hist.SetMarkerStyle(24)
        hist.SetMarkerSize(0.65)
        hist.SetLineWidth(2)
        hist.SetStats(False)
    #endfor

    for hist in ratio_histograms:
        hist.SetLineColor(ROOT.kBlack)
        hist.SetMarkerColor(ROOT.kBlack)
        hist.SetMarkerStyle(20)
        hist.SetMarkerSize(0.65)
        hist.SetLineWidth(2)
        hist.SetStats(False)
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

        y_max = 10.0 * global_max
    else:
        y_min = LINEAR_Y_MIN
        y_max = 1.15 * global_max
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


def get_ratio_y_range(log_y):
    if log_y:
        return RATIO_LOG_Y_MIN, RATIO_Y_MAX
    #endif

    return RATIO_LINEAR_Y_MIN, RATIO_Y_MAX


def draw_normalization_pad(output_tag, total_charge_mc, integrated_luminosity_pb_inv, n_gen, normalization_factor, log_y, ratio_mode):
    pad4 = ROOT.gPad
    pad4.Clear()
    pad4.SetFillColor(ROOT.kWhite)

    latex = ROOT.TLatex()
    latex.SetNDC(True)
    latex.SetTextSize(0.035)

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

    latex.DrawLatex(0.10, 0.88, output_tag)
    latex.DrawLatex(0.10, 0.78, title)
    latex.DrawLatex(0.10, 0.68, "scale: {}".format(scale_label))
    latex.DrawLatex(0.10, 0.58, "Q = %.6g mC" % total_charge_mc)
    latex.DrawLatex(0.10, 0.48, "L_{int} = %.6g pb^{-1}" % integrated_luminosity_pb_inv)
    latex.DrawLatex(0.10, 0.38, "N_{gen} = %.6g" % n_gen)
    latex.DrawLatex(0.10, 0.28, "scale = L_{int}/N_{gen}")
    latex.DrawLatex(0.10, 0.18, "scale = %.6g pb^{-1}" % normalization_factor)

    if ratio_mode:
        latex.DrawLatex(0.10, 0.08, "ratio: data / MC")
    else:
        latex.DrawLatex(0.10, 0.08, "MC fill: scale #times weight")
    #endif


def draw_comparison_canvas(output_tag, data_histograms, mc_histograms, output_pdf, normalization_factor, total_charge_mc, integrated_luminosity_pb_inv, n_gen, log_y):
    if log_y:
        status("Drawing data-vs-MC output canvas with log y scale.")
    else:
        status("Drawing data-vs-MC output canvas with linear y scale.")
    #endif

    ROOT.gStyle.SetOptStat(0)

    if log_y:
        canvas_name = "canvas_comparison_log"
        canvas_title = "{} data vs MC p1_theta log y".format(output_tag)
    else:
        canvas_name = "canvas_comparison_linear"
        canvas_title = "{} data vs MC p1_theta linear y".format(output_tag)
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

    status("Comparison canvas FD common y-range: [{:.12g}, {:.12g}]".format(fd_y_min, fd_y_max))
    status("Comparison canvas CD y-range: [{:.12g}, {:.12g}]".format(cd_y_min, cd_y_max))

    for i_panel in range(7):
        canvas.cd(canvas_pad_for_panel[i_panel])

        pad = ROOT.gPad
        pad.SetLeftMargin(0.14)
        pad.SetRightMargin(0.05)
        pad.SetTopMargin(0.10)
        pad.SetBottomMargin(0.13)
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
        h_data.GetXaxis().SetTitle("p_{1} #theta (deg)")
        h_data.GetYaxis().SetTitle("Counts")
        h_data.GetXaxis().CenterTitle(True)
        h_data.GetYaxis().CenterTitle(True)
        h_data.GetXaxis().SetTitleSize(0.050)
        h_data.GetYaxis().SetTitleSize(0.050)
        h_data.GetXaxis().SetLabelSize(0.045)
        h_data.GetYaxis().SetLabelSize(0.045)

        h_data.Draw("E1")
        h_mc.Draw("HIST SAME")

        legend = ROOT.TLegend(0.66, 0.76, 0.92, 0.89)
        legend.SetBorderSize(1)
        legend.SetFillStyle(1001)
        legend.SetFillColor(ROOT.kWhite)
        legend.SetTextSize(0.032)
        legend.AddEntry(h_data, "data", "lep")
        legend.AddEntry(h_mc, "MC weighted", "l")
        legend.Draw()
    #endfor

    canvas.cd(4)
    draw_normalization_pad(output_tag, total_charge_mc, integrated_luminosity_pb_inv, n_gen, normalization_factor, log_y, False)

    ensure_output_directory(output_pdf)

    status("Saving output PDF: {}".format(output_pdf))
    canvas.SaveAs(output_pdf)
    status("Saved output PDF.")


def draw_ratio_canvas(output_tag, ratio_histograms, output_pdf, normalization_factor, total_charge_mc, integrated_luminosity_pb_inv, n_gen, log_y):
    if log_y:
        status("Drawing data/MC ratio output canvas with log y scale.")
    else:
        status("Drawing data/MC ratio output canvas with linear y scale.")
    #endif

    ROOT.gStyle.SetOptStat(0)

    if log_y:
        canvas_name = "canvas_ratio_log"
        canvas_title = "{} data over MC p1_theta log y".format(output_tag)
    else:
        canvas_name = "canvas_ratio_linear"
        canvas_title = "{} data over MC p1_theta linear y".format(output_tag)
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

    ratio_y_min, ratio_y_max = get_ratio_y_range(log_y)

    status("Ratio canvas y-range: [{:.12g}, {:.12g}]".format(ratio_y_min, ratio_y_max))

    for i_panel in range(7):
        canvas.cd(canvas_pad_for_panel[i_panel])

        pad = ROOT.gPad
        pad.SetLeftMargin(0.14)
        pad.SetRightMargin(0.05)
        pad.SetTopMargin(0.10)
        pad.SetBottomMargin(0.13)
        pad.SetLogy(log_y)

        h_ratio = ratio_histograms[i_panel]

        h_ratio.SetMaximum(ratio_y_max)
        h_ratio.SetMinimum(ratio_y_min)

        h_ratio.SetTitle("{}  {}".format(output_tag, panel_labels[i_panel]))
        h_ratio.GetXaxis().SetTitle("p_{1} #theta (deg)")
        h_ratio.GetYaxis().SetTitle("data / MC")
        h_ratio.GetXaxis().CenterTitle(True)
        h_ratio.GetYaxis().CenterTitle(True)
        h_ratio.GetXaxis().SetTitleSize(0.050)
        h_ratio.GetYaxis().SetTitleSize(0.050)
        h_ratio.GetXaxis().SetLabelSize(0.045)
        h_ratio.GetYaxis().SetLabelSize(0.045)

        h_ratio.Draw("E1")

        if (not log_y) and ratio_y_min < 1.0 and ratio_y_max > 1.0:
            line = ROOT.TLine(h_ratio.GetXaxis().GetXmin(), 1.0, h_ratio.GetXaxis().GetXmax(), 1.0)
            line.SetLineColor(ROOT.kRed + 1)
            line.SetLineStyle(2)
            line.SetLineWidth(2)
            line.Draw("SAME")
        #endif
    #endfor

    canvas.cd(4)
    draw_normalization_pad(output_tag, total_charge_mc, integrated_luminosity_pb_inv, n_gen, normalization_factor, log_y, True)

    ensure_output_directory(output_pdf)

    status("Saving output PDF: {}".format(output_pdf))
    canvas.SaveAs(output_pdf)
    status("Saved output PDF.")


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compare data and reconstructed MC p1_theta using theta-defined FD/CD regions with event-by-event MC::Event.weight normalization."
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
    output_paths = build_output_paths(output_tag)

    start_time = time.time()

    status("Starting data/MC normalization script.")
    status("Output tag: {}".format(output_tag))
    status("Data file: {}".format(args.data_root))
    status("Reco MC file: {}".format(args.reco_mc_root))
    status("Generated MC file: {}".format(args.gen_mc_root))
    status("Output directory: {}".format(output_paths["output_dir"]))
    status("Output comparison log PDF: {}".format(output_paths["comparison_log"]))
    status("Output comparison linear PDF: {}".format(output_paths["comparison_linear"]))
    status("Output ratio log PDF: {}".format(output_paths["ratio_log"]))
    status("Output ratio linear PDF: {}".format(output_paths["ratio_linear"]))
    status("Maximum worker processes: {}".format(max_workers))
    status("Using event-by-event reconstructed MC branch: weight")
    status("Detector definition override: FD = p1_theta < 40 deg; CD = p1_theta >= 40 deg.")
    status("FD histograms use theta range 0 to 40 deg; CD histogram uses theta range 40 to 70 deg.")
    status("Comparison y scales: common across FD sectors; CD uses its own scale.")
    status("Ratio y scales: fixed to 0 to 2 for linear, and positive lower bound to 2 for log.")

    ROOT.gROOT.SetBatch(True)

    n_data_entries = check_input_tree(
        args.data_root,
        "PhysicsEvents",
        ["runnum", "p1_phi", "p1_theta"],
        "data"
    )

    n_reco_mc_entries = check_input_tree(
        args.reco_mc_root,
        "PhysicsEvents",
        ["p1_phi", "p1_theta", "weight"],
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

    # Weighted-MC normalization:
    #
    #   per-event plot weight = L_int * MC::Event.weight / N_GEN
    #
    # The reconstructed ROOT tree now carries MC::Event.weight as the branch:
    #
    #   weight
    #
    # Therefore the worker fills MC histograms using:
    #
    #   event_weight = normalization_factor * weight
    #
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

    status("Building ROOT histograms from accumulated bin arrays.")

    data_histograms = arrays_to_histograms("data", data_result["counts"], data_result["sumw2"])
    mc_histograms = arrays_to_histograms("mc", mc_result["counts"], mc_result["sumw2"])
    ratio_histograms = make_ratio_histograms(data_histograms, mc_histograms)

    style_histograms(data_histograms, mc_histograms, ratio_histograms)

    draw_comparison_canvas(
        output_tag,
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
        ratio_histograms,
        output_paths["ratio_log"],
        normalization_factor,
        total_charge_mc,
        integrated_luminosity_pb_inv,
        n_gen,
        True
    )

    draw_ratio_canvas(
        output_tag,
        ratio_histograms,
        output_paths["ratio_linear"],
        normalization_factor,
        total_charge_mc,
        integrated_luminosity_pb_inv,
        n_gen,
        False
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
    print("FD definition: p1_theta < 40 deg")
    print("CD definition: p1_theta >= 40 deg")
    print("FD theta plotting range: {:.1f} to {:.1f} deg".format(FD_THETA_MIN_DEG, FD_THETA_MAX_DEG))
    print("CD theta plotting range: {:.1f} to {:.1f} deg".format(CD_THETA_MIN_DEG, CD_THETA_MAX_DEG))
    print("Comparison y scales: common across FD sectors; CD uses its own scale")
    print("Ratio y range linear: {:.12g} to {:.12g}".format(RATIO_LINEAR_Y_MIN, RATIO_Y_MAX))
    print("Ratio y range log: {:.12g} to {:.12g}".format(RATIO_LOG_Y_MIN, RATIO_Y_MAX))
    print("Total accumulated charge from CSV raw units: {:.12g}".format(total_charge_raw))
    print("Charge conversion factor to mC: {:.12g}".format(args.charge_to_mc_factor))
    print("Total accumulated charge Q: {:.12g} mC".format(total_charge_mc))
    print("Integrated luminosity: {:.12g} pb^-1".format(integrated_luminosity_pb_inv))
    print("Generated MC entries from file: {}".format(n_gen_from_file))
    print("N_GEN used for normalization: {:.12g}".format(n_gen))
    print("Normalization factor L_int / N_GEN: {:.12g} pb^-1".format(normalization_factor))
    print("Data entries filled: {}".format(data_result["n_filled"]))
    print("Data FD entries filled: {}".format(data_result["n_fd"]))
    print("Data CD entries filled: {}".format(data_result["n_cd"]))
    print("Reco MC entries filled: {}".format(mc_result["n_filled"]))
    print("Reco MC FD entries filled: {}".format(mc_result["n_fd"]))
    print("Reco MC CD entries filled: {}".format(mc_result["n_cd"]))
    print("Reco MC raw weight sum over filled entries: {:.12g}".format(mc_result["raw_weight_sum"]))
    print("Reco MC scaled weight sum over filled entries: {:.12g}".format(mc_result["scaled_weight_sum"]))
    print("Reco MC maximum raw weight over filled entries: {:.12g}".format(mc_result["max_raw_weight"]))
    print("Output directory: {}".format(output_paths["output_dir"]))
    print("Output comparison log PDF: {}".format(output_paths["comparison_log"]))
    print("Output comparison linear PDF: {}".format(output_paths["comparison_linear"]))
    print("Output ratio log PDF: {}".format(output_paths["ratio_log"]))
    print("Output ratio linear PDF: {}".format(output_paths["ratio_linear"]))
    print("Elapsed time: {:.2f} seconds".format(elapsed_time))

    status("Finished data/MC normalization script.")


if __name__ == "__main__":
    main()
#endif