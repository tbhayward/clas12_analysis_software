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

# Temporary approximation:
#
# The generated MC has event-by-event MC::Event.weight values, but the current
# reconstructed ROOT tree does not yet carry those weights. For now, use the
# average generated MC::Event.weight as an effective constant generator weight.
#
# You reported:
#
#   sum(MC::Event.weight) = 4 * 3.508383992386323E11
#   N_GEN                 = 127630261
#
# Therefore:
#
#   average generated weight = sum(MC::Event.weight) / N_GEN
#
DEFAULT_AVG_GEN_WEIGHT_PB = (4.0 * 3.508383992386323e11) / 127630261.0

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

OUTPUT_PDF = "output/data_mc_normalization/output.pdf"

# Panels 1-6 are FD sectors, now defined only by p1_theta < 40 degrees.
# Panel 7 is CD, now defined only by p1_theta >= 40 degrees.
FD_THETA_MIN_DEG = 0.0
FD_THETA_MAX_DEG = 40.0

CD_THETA_MIN_DEG = 40.0
CD_THETA_MAX_DEG = 70.0

N_BINS_THETA = 70

DEFAULT_STATUS_EVERY = 250000
DEFAULT_MAX_WORKERS = 5

LOG_Y_MIN = 0.5

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


def ensure_output_directory(output_path):
    output_dir = os.path.dirname(output_path)

    if output_dir != "":
        os.makedirs(output_dir, exist_ok=True)
    #endif


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
                counts[i_panel][i_bin] += 1.0
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

    n_filled = 0
    n_fd = 0
    n_cd = 0
    n_skipped_panel = 0
    n_skipped_theta = 0
    n_total = end_entry - start_entry

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

        i_panel = get_panel_index_from_theta_and_phi(p1_theta_rad, p1_phi_rad)

        if i_panel < 0:
            n_skipped_panel += 1
        else:
            i_bin = get_theta_bin_for_panel(i_panel, p1_theta_rad)

            if i_bin < 0:
                n_skipped_theta += 1
            else:
                counts[i_panel][i_bin] += 1.0
                n_filled += 1

                if i_panel >= 0 and i_panel <= 5:
                    n_fd += 1
                elif i_panel == 6:
                    n_cd += 1
                #endif
            #endif
        #endif

        if local_index % status_every == 0:
            print("[{}] MC worker {} progress: {}/{} chunk entries ({}) | filled: {} | FD(theta<40): {} | CD(theta>=40): {}".format(
                time.strftime("%Y-%m-%d %H:%M:%S"),
                worker_id,
                local_index,
                n_total,
                format_percent(local_index, n_total),
                n_filled,
                n_fd,
                n_cd
            ), flush=True)
        #endif
    #endfor

    root_file.Close()

    print("[{}] MC worker {} finished: filled = {}, FD(theta<40) = {}, CD(theta>=40) = {}, skipped panel = {}, skipped theta = {}.".format(
        time.strftime("%Y-%m-%d %H:%M:%S"),
        worker_id,
        n_filled,
        n_fd,
        n_cd,
        n_skipped_panel,
        n_skipped_theta
    ), flush=True)

    return {
        "counts": counts,
        "n_filled": n_filled,
        "n_fd": n_fd,
        "n_cd": n_cd,
        "n_skipped_panel": n_skipped_panel,
        "n_skipped_theta": n_skipped_theta,
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
        "unique_runs": unique_runs,
        "n_filled": n_filled,
        "n_fd": n_fd,
        "n_cd": n_cd,
        "n_skipped_panel": n_skipped_panel,
        "n_skipped_theta": n_skipped_theta,
    }


def run_mc_parallel(root_path, tree_name, n_entries, max_workers, status_every):
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
        })
    #endfor

    status("Starting parallel reconstructed MC pass with {} workers.".format(len(tasks)))

    total_counts = make_empty_counts()

    n_filled = 0
    n_fd = 0
    n_cd = 0
    n_skipped_panel = 0
    n_skipped_theta = 0

    with mp.Pool(processes=len(tasks)) as pool:
        results = pool.map(mc_worker, tasks)
    #endwith

    for result in results:
        add_counts(total_counts, result["counts"])
        n_filled += result["n_filled"]
        n_fd += result["n_fd"]
        n_cd += result["n_cd"]
        n_skipped_panel += result["n_skipped_panel"]
        n_skipped_theta += result["n_skipped_theta"]
    #endfor

    status("Finished parallel reconstructed MC pass.")
    status("Reco MC total: filled = {}, FD(theta<40) = {}, CD(theta>=40) = {}, skipped panel = {}, skipped theta = {}.".format(
        n_filled,
        n_fd,
        n_cd,
        n_skipped_panel,
        n_skipped_theta
    ))

    return {
        "counts": total_counts,
        "n_filled": n_filled,
        "n_fd": n_fd,
        "n_cd": n_cd,
        "n_skipped_panel": n_skipped_panel,
        "n_skipped_theta": n_skipped_theta,
    }


# -----------------------------------------------------------------------------
# Histogram and plotting helpers
# -----------------------------------------------------------------------------

def counts_to_histograms(prefix, counts):
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
            hist.SetBinContent(root_bin, content)

            if content >= 0.0:
                hist.SetBinError(root_bin, math.sqrt(content))
            else:
                hist.SetBinError(root_bin, 0.0)
            #endif
        #endfor

        histograms.append(hist)
    #endfor

    return histograms


def scaled_counts_to_histograms(prefix, counts, scale):
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
            raw_count = counts[i_panel][i_bin]
            content = scale * raw_count
            error = abs(scale) * math.sqrt(raw_count)

            hist.SetBinContent(root_bin, content)
            hist.SetBinError(root_bin, error)
        #endfor

        histograms.append(hist)
    #endfor

    return histograms


def style_histograms(data_histograms, mc_histograms):
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


def draw_canvas(data_histograms, mc_histograms, output_pdf, mc_event_weight, avg_gen_weight_pb, total_charge_mc, integrated_luminosity_pb_inv, n_gen):
    status("Drawing output canvas.")

    ROOT.gStyle.SetOptStat(0)

    canvas = ROOT.TCanvas("canvas", "data vs MC p1_theta", 1600, 900)
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

    for i_panel in range(7):
        canvas.cd(canvas_pad_for_panel[i_panel])

        pad = ROOT.gPad
        pad.SetLeftMargin(0.14)
        pad.SetRightMargin(0.05)
        pad.SetTopMargin(0.10)
        pad.SetBottomMargin(0.13)
        pad.SetLogy(True)

        h_data = data_histograms[i_panel]
        h_mc = mc_histograms[i_panel]

        max_y = max(h_data.GetMaximum(), h_mc.GetMaximum())

        if max_y <= 0.0:
            max_y = 1.0
        #endif

        h_data.SetMaximum(10.0 * max_y)
        h_data.SetMinimum(LOG_Y_MIN)

        h_data.SetTitle(panel_labels[i_panel])
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
        legend.AddEntry(h_mc, "MC scaled", "l")
        legend.Draw()
    #endfor

    canvas.cd(4)
    pad4 = ROOT.gPad
    pad4.Clear()
    pad4.SetFillColor(ROOT.kWhite)

    latex = ROOT.TLatex()
    latex.SetNDC(True)
    latex.SetTextSize(0.037)
    latex.DrawLatex(0.10, 0.84, "MC normalization")
    latex.DrawLatex(0.10, 0.74, "Q = %.6g mC" % total_charge_mc)
    latex.DrawLatex(0.10, 0.64, "L_{int} = %.6g pb^{-1}" % integrated_luminosity_pb_inv)
    latex.DrawLatex(0.10, 0.54, "#LT w_{gen} #GT = %.6g pb" % avg_gen_weight_pb)
    latex.DrawLatex(0.10, 0.44, "N_{gen} = %d" % n_gen)
    latex.DrawLatex(0.10, 0.34, "event weight = %.6g" % mc_event_weight)
    latex.DrawLatex(0.10, 0.24, "temporary average-weight approx.")
    latex.DrawLatex(0.10, 0.14, "FD: #theta < 40 deg, CD: #theta #geq 40 deg")

    ensure_output_directory(output_pdf)

    status("Saving output PDF: {}".format(output_pdf))
    canvas.SaveAs(output_pdf)
    status("Saved output PDF.")


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compare data and reconstructed MC p1_theta using theta-defined FD/CD regions with temporary average-weight MC normalization."
    )

    parser.add_argument("data_root", help="Input data ROOT file containing PhysicsEvents")
    parser.add_argument("reco_mc_root", help="Input reconstructed MC ROOT file containing PhysicsEvents")
    parser.add_argument("gen_mc_root", help="Input generated MC ROOT file containing PhysicsEvents")
    parser.add_argument("--output", default=OUTPUT_PDF, help="Output PDF path")
    parser.add_argument("--avg-gen-weight-pb", type=float, default=DEFAULT_AVG_GEN_WEIGHT_PB, help="Temporary average generated MC::Event.weight value in pb")
    parser.add_argument("--charge-csv", default=GLOBAL_CHARGE_CSV, help="CSV containing run number and accumulated charge")
    parser.add_argument("--charge-to-mc-factor", type=float, default=CHARGE_TO_MC_FACTOR, help="Conversion factor from CSV charge units to mC")
    parser.add_argument("--status-every", type=int, default=DEFAULT_STATUS_EVERY, help="Print loop progress every N entries per worker")
    parser.add_argument("--max-workers", type=int, default=DEFAULT_MAX_WORKERS, help="Maximum number of worker processes. Hard capped at 5.")

    args = parser.parse_args()

    if args.status_every <= 0:
        fatal("--status-every must be positive")
    #endif

    if args.max_workers <= 0:
        fatal("--max-workers must be positive")
    #endif

    max_workers = min(args.max_workers, 5)

    start_time = time.time()

    status("Starting data/MC normalization script.")
    status("Data file: {}".format(args.data_root))
    status("Reco MC file: {}".format(args.reco_mc_root))
    status("Generated MC file: {}".format(args.gen_mc_root))
    status("Output PDF: {}".format(args.output))
    status("Maximum worker processes: {}".format(max_workers))
    status("Temporary average generated weight: {:.12g} pb".format(args.avg_gen_weight_pb))
    status("Detector definition override: FD = p1_theta < 40 deg; CD = p1_theta >= 40 deg.")
    status("FD histograms use theta range 0 to 40 deg; CD histogram uses theta range 40 to 70 deg.")

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
        ["p1_phi", "p1_theta"],
        "reconstructed MC"
    )

    n_gen = get_entry_count(args.gen_mc_root, "PhysicsEvents", "generated MC")

    status("Input entry counts:")
    status("  data entries = {}".format(n_data_entries))
    status("  reco MC entries = {}".format(n_reco_mc_entries))
    status("  generated MC entries N_GEN = {}".format(n_gen))

    if n_gen <= 0:
        fatal("generated MC tree has zero entries")
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

    # Temporary approximation:
    #
    # Proper weighted MC should use:
    #
    #   per-event plot weight = L_int * mc_weight / N_GEN
    #
    # where mc_weight is MC::Event.weight carried into the reconstructed ROOT tree.
    #
    # For now, because the reconstructed tree does not carry mc_weight, use:
    #
    #   per-event plot weight = L_int * <mc_weight> / N_GEN
    #
    mc_event_weight = integrated_luminosity_pb_inv * args.avg_gen_weight_pb / float(n_gen)

    status("Computed temporary average-weight normalization:")
    status("  Raw charge sum = {:.12g}".format(total_charge_raw))
    status("  Q = {:.12g} mC".format(total_charge_mc))
    status("  L_int = {:.12g} pb^-1".format(integrated_luminosity_pb_inv))
    status("  <MC::Event.weight> = {:.12g} pb".format(args.avg_gen_weight_pb))
    status("  N_GEN = {}".format(n_gen))
    status("  MC event weight = {:.12g}".format(mc_event_weight))

    mc_result = run_mc_parallel(
        args.reco_mc_root,
        "PhysicsEvents",
        n_reco_mc_entries,
        max_workers,
        args.status_every
    )

    status("Building ROOT histograms from accumulated bin arrays.")

    data_histograms = counts_to_histograms("data", data_result["counts"])
    mc_histograms = scaled_counts_to_histograms("mc", mc_result["counts"], mc_event_weight)

    style_histograms(data_histograms, mc_histograms)

    draw_canvas(
        data_histograms,
        mc_histograms,
        args.output,
        mc_event_weight,
        args.avg_gen_weight_pb,
        total_charge_mc,
        integrated_luminosity_pb_inv,
        n_gen
    )

    elapsed_time = time.time() - start_time

    print("")
    print("Normalization summary")
    print("---------------------")
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
    print("Total accumulated charge from CSV raw units: {:.12g}".format(total_charge_raw))
    print("Charge conversion factor to mC: {:.12g}".format(args.charge_to_mc_factor))
    print("Total accumulated charge Q: {:.12g} mC".format(total_charge_mc))
    print("Integrated luminosity: {:.12g} pb^-1".format(integrated_luminosity_pb_inv))
    print("Temporary average MC::Event.weight: {:.12g} pb".format(args.avg_gen_weight_pb))
    print("N_GEN: {}".format(n_gen))
    print("MC event weight: {:.12g}".format(mc_event_weight))
    print("Data entries filled: {}".format(data_result["n_filled"]))
    print("Data FD entries filled: {}".format(data_result["n_fd"]))
    print("Data CD entries filled: {}".format(data_result["n_cd"]))
    print("Reco MC entries filled: {}".format(mc_result["n_filled"]))
    print("Reco MC FD entries filled: {}".format(mc_result["n_fd"]))
    print("Reco MC CD entries filled: {}".format(mc_result["n_cd"]))
    print("Output PDF: {}".format(args.output))
    print("Elapsed time: {:.2f} seconds".format(elapsed_time))

    status("Finished data/MC normalization script.")


if __name__ == "__main__":
    main()
#endif