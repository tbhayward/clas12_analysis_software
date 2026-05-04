#!/usr/bin/env python3

import os
import sys
import csv
import math
import argparse
import ROOT

# -----------------------------------------------------------------------------
# User-editable constants
# -----------------------------------------------------------------------------

GLOBAL_CHARGE_CSV = "/u/home/thayward/clas12_analysis_software/analysis_scripts/dvcs_cross_section/imports/integrated_luminosity/global.csv"

# Generator total cross section in pb.
DEFAULT_SIGMA_GEN_TOT_PB = 294734.8125

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

N_BINS_THETA = 70
THETA_MIN_DEG = 0.0
THETA_MAX_DEG = 70.0

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

def fatal(message):
    print("ERROR: {}".format(message), file=sys.stderr)
    sys.exit(1)


def ensure_output_directory(output_path):
    output_dir = os.path.dirname(output_path)

    if output_dir != "":
        os.makedirs(output_dir, exist_ok=True)
    #endif


def open_root_file(path, label):
    if not os.path.isfile(path):
        fatal("{} ROOT file does not exist: {}".format(label, path))
    #endif

    root_file = ROOT.TFile.Open(path, "READ")

    if root_file is None or root_file.IsZombie():
        fatal("could not open {} ROOT file: {}".format(label, path))
    #endif

    return root_file


def get_tree(root_file, tree_name, label):
    tree = root_file.Get(tree_name)

    if tree is None:
        fatal("{} ROOT file does not contain tree '{}': {}".format(label, tree_name, root_file.GetName()))
    #endif

    return tree


def has_branch(tree, branch_name):
    return tree.GetBranch(branch_name) is not None


def require_branches(tree, branch_names, label):
    missing = []

    for branch_name in branch_names:
        if not has_branch(tree, branch_name):
            missing.append(branch_name)
        #endif
    #endfor

    if len(missing) > 0:
        fatal("{} tree is missing required branches: {}".format(label, ", ".join(missing)))
    #endif


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


def get_panel_index(detector1, p1_phi_rad):
    if detector1 == 1:
        p1_phi_deg = math.degrees(p1_phi_rad)
        sector = get_fd_sector_from_phi_deg(p1_phi_deg)

        if sector >= 1 and sector <= 6:
            return sector - 1
        #endif

        return -1
    #endif

    if detector1 == 2:
        return 6
    #endif

    return -1


def load_charge_map(csv_path):
    if not os.path.isfile(csv_path):
        fatal("charge CSV does not exist: {}".format(csv_path))
    #endif

    charge_by_run = {}

    with open(csv_path, "r") as f:
        reader = csv.reader(f)

        for row in reader:
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
        #endfor
    #endwith

    if len(charge_by_run) == 0:
        fatal("no valid run-charge rows were read from {}".format(csv_path))
    #endif

    return charge_by_run


def collect_unique_data_runs(data_tree):
    unique_runs = set()

    n_entries = data_tree.GetEntries()

    for i_entry in range(n_entries):
        data_tree.GetEntry(i_entry)

        runnum = int(getattr(data_tree, "runnum"))
        unique_runs.add(runnum)
    #endfor

    return unique_runs


def sum_charge_for_runs(unique_runs, charge_by_run):
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

    return total_charge_raw


def make_histograms(prefix):
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
        hist_name = "{}_panel_{}".format(prefix, i_panel + 1)
        hist_title = "{};p_{{1}} #theta (deg);Counts".format(panel_names[i_panel])

        hist = ROOT.TH1D(
            hist_name,
            hist_title,
            N_BINS_THETA,
            THETA_MIN_DEG,
            THETA_MAX_DEG
        )

        hist.Sumw2()
        histograms.append(hist)
    #endfor

    return histograms


def fill_histograms(tree, histograms, event_weight):
    n_entries = tree.GetEntries()
    n_filled = 0

    for i_entry in range(n_entries):
        tree.GetEntry(i_entry)

        detector1 = int(getattr(tree, "detector1"))
        p1_phi_rad = float(getattr(tree, "p1_phi"))
        p1_theta_rad = float(getattr(tree, "p1_theta"))

        i_panel = get_panel_index(detector1, p1_phi_rad)

        if i_panel < 0 or i_panel >= len(histograms):
            continue
        #endif

        p1_theta_deg = math.degrees(p1_theta_rad)

        if p1_theta_deg < THETA_MIN_DEG or p1_theta_deg >= THETA_MAX_DEG:
            continue
        #endif

        histograms[i_panel].Fill(p1_theta_deg, event_weight)
        n_filled += 1
    #endfor

    return n_filled


def style_histograms(data_histograms, mc_histograms):
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


def draw_canvas(data_histograms, mc_histograms, output_pdf, mc_event_weight, sigma_gen_tot_pb):
    ROOT.gStyle.SetOptStat(0)

    canvas = ROOT.TCanvas("canvas", "data vs MC p1_theta", 1500, 1200)
    canvas.Divide(3, 3)

    panel_labels = [
        "FD sector 1",
        "FD sector 2",
        "FD sector 3",
        "FD sector 4",
        "FD sector 5",
        "FD sector 6",
        "Central detector",
    ]

    for i_panel in range(7):
        canvas.cd(i_panel + 1)

        pad = ROOT.gPad
        pad.SetLeftMargin(0.14)
        pad.SetRightMargin(0.05)
        pad.SetTopMargin(0.10)
        pad.SetBottomMargin(0.13)
        pad.SetGrid()

        h_data = data_histograms[i_panel]
        h_mc = mc_histograms[i_panel]

        max_y = max(h_data.GetMaximum(), h_mc.GetMaximum())

        if max_y <= 0.0:
            max_y = 1.0
        #endif

        h_data.SetMaximum(1.25 * max_y)
        h_data.SetMinimum(0.0)

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

        legend = ROOT.TLegend(0.58, 0.72, 0.92, 0.88)
        legend.SetBorderSize(1)
        legend.SetFillStyle(1001)
        legend.SetFillColor(ROOT.kWhite)
        legend.SetTextSize(0.040)
        legend.AddEntry(h_data, "data", "lep")
        legend.AddEntry(h_mc, "MC scaled", "l")
        legend.Draw()
    #endfor

    canvas.cd(8)
    pad8 = ROOT.gPad
    pad8.Clear()
    pad8.SetFillColor(ROOT.kWhite)

    latex = ROOT.TLatex()
    latex.SetNDC(True)
    latex.SetTextSize(0.045)
    latex.DrawLatex(0.12, 0.80, "MC normalization")
    latex.DrawLatex(0.12, 0.70, "#sigma_{gen}^{tot} = %.6g pb" % sigma_gen_tot_pb)
    latex.DrawLatex(0.12, 0.60, "event weight = %.6g" % mc_event_weight)

    canvas.cd(9)
    pad9 = ROOT.gPad
    pad9.Clear()
    pad9.SetFillColor(ROOT.kWhite)

    ensure_output_directory(output_pdf)
    canvas.SaveAs(output_pdf)


def main():
    parser = argparse.ArgumentParser(
        description="Compare data and reconstructed MC p1_theta by FD sector and CD, with MC normalized to expected yield."
    )

    parser.add_argument("data_root", help="Input data ROOT file containing PhysicsEvents")
    parser.add_argument("reco_mc_root", help="Input reconstructed MC ROOT file containing PhysicsEvents")
    parser.add_argument("gen_mc_root", help="Input generated MC ROOT file containing PhysicsEvents")
    parser.add_argument("--output", default=OUTPUT_PDF, help="Output PDF path")
    parser.add_argument("--sigma-gen-tot-pb", type=float, default=DEFAULT_SIGMA_GEN_TOT_PB, help="Generator total cross section in pb")
    parser.add_argument("--charge-csv", default=GLOBAL_CHARGE_CSV, help="CSV containing run number and accumulated charge")
    parser.add_argument("--charge-to-mc-factor", type=float, default=CHARGE_TO_MC_FACTOR, help="Conversion factor from CSV charge units to mC")

    args = parser.parse_args()

    data_file = open_root_file(args.data_root, "data")
    reco_mc_file = open_root_file(args.reco_mc_root, "reconstructed MC")
    gen_mc_file = open_root_file(args.gen_mc_root, "generated MC")

    data_tree = get_tree(data_file, "PhysicsEvents", "data")
    reco_mc_tree = get_tree(reco_mc_file, "PhysicsEvents", "reconstructed MC")
    gen_mc_tree = get_tree(gen_mc_file, "PhysicsEvents", "generated MC")

    require_branches(data_tree, ["runnum", "detector1", "p1_phi", "p1_theta"], "data")
    require_branches(reco_mc_tree, ["detector1", "p1_phi", "p1_theta"], "reconstructed MC")

    n_gen = gen_mc_tree.GetEntries()

    if n_gen <= 0:
        fatal("generated MC tree has zero entries")
    #endif

    charge_by_run = load_charge_map(args.charge_csv)
    unique_runs = collect_unique_data_runs(data_tree)
    total_charge_raw = sum_charge_for_runs(unique_runs, charge_by_run)
    total_charge_mc = total_charge_raw * args.charge_to_mc_factor

    integrated_luminosity_pb_inv = RGA_LUMINOSITY_FACTOR_PB_INV_PER_MC * total_charge_mc
    mc_event_weight = integrated_luminosity_pb_inv * args.sigma_gen_tot_pb / float(n_gen)

    data_histograms = make_histograms("data")
    mc_histograms = make_histograms("mc")

    n_data_filled = fill_histograms(data_tree, data_histograms, 1.0)
    n_mc_filled = fill_histograms(reco_mc_tree, mc_histograms, mc_event_weight)

    style_histograms(data_histograms, mc_histograms)
    draw_canvas(
        data_histograms,
        mc_histograms,
        args.output,
        mc_event_weight,
        args.sigma_gen_tot_pb
    )

    print("")
    print("Normalization summary")
    print("---------------------")
    print("Data ROOT file: {}".format(args.data_root))
    print("Reco MC ROOT file: {}".format(args.reco_mc_root))
    print("Gen MC ROOT file: {}".format(args.gen_mc_root))
    print("Charge CSV: {}".format(args.charge_csv))
    print("Unique data runs found: {}".format(len(unique_runs)))
    print("Total accumulated charge from CSV raw units: {:.12g}".format(total_charge_raw))
    print("Charge conversion factor to mC: {:.12g}".format(args.charge_to_mc_factor))
    print("Total accumulated charge Q: {:.12g} mC".format(total_charge_mc))
    print("Integrated luminosity: {:.12g} pb^-1".format(integrated_luminosity_pb_inv))
    print("sigma_GEN_TOT: {:.12g} pb".format(args.sigma_gen_tot_pb))
    print("N_GEN: {}".format(n_gen))
    print("MC event weight: {:.12g}".format(mc_event_weight))
    print("Data entries filled: {}".format(n_data_filled))
    print("Reco MC entries filled: {}".format(n_mc_filled))
    print("Output PDF: {}".format(args.output))

    data_file.Close()
    reco_mc_file.Close()
    gen_mc_file.Close()


if __name__ == "__main__":
    main()
#endif