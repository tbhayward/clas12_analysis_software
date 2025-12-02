#!/usr/bin/env python3

import os
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)


# ----------------------------------------------------------------------
# File configuration
# ----------------------------------------------------------------------

PERIODS = ["Sp18 Inb", "Sp18 Out", "Fa18 Inb", "Fa18 Out", "Sp19 Inb"]

DVCS_FILES = {
    "Sp18 Inb": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_sp18_inb_10594MeV.root",
    "Sp18 Out": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_sp18_out_10594MeV.root",
    "Fa18 Inb": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_fa18_inb_10604MeV.root",
    "Fa18 Out": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_fa18_out_10604MeV.root",
    "Sp19 Inb": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_sp19_inb_10200MeV.root",
}

PI0_FILES = {
    "Sp18 Inb": "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/gen_aaogen_norad_sp18_inb_10594MeV.root",
    "Sp18 Out": "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/gen_aaogen_norad_sp18_out_10594MeV.root",
    "Fa18 Inb": "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/gen_aaogen_norad_fa18_inb_10604MeV.root",
    "Fa18 Out": "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/gen_aaogen_norad_fa18_out_10604MeV.root",
    "Sp19 Inb": "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/gen_aaogen_norad_sp19_inb_10200MeV.root",
}

TREE_NAME = "PhysicsEvents"
MAX_EVENTS = 100000


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

def load_trees(file_map):
    files = {}
    trees = {}
    for period in PERIODS:
        path = file_map[period]
        f = ROOT.TFile.Open(path)
        if not f or f.IsZombie():
            raise RuntimeError("Failed to open file: {0}".format(path))
        #endif
        tree = f.Get(TREE_NAME)
        if not tree:
            raise RuntimeError("Tree {0} not found in file {1}".format(TREE_NAME, path))
        #endif
        files[period] = f
        trees[period] = tree
    #endfor
    return files, trees
#enddef


def make_2d_hist(tree, name, title,
                 x_branch, y_branch,
                 x_bins, x_min, x_max,
                 y_bins, y_min, y_max,
                 max_events=MAX_EVENTS,
                 negate_y=False):
    hist = ROOT.TH2D(name, title,
                     x_bins, x_min, x_max,
                     y_bins, y_min, y_max)
    hist.Sumw2()

    n_entries = int(min(max_events, tree.GetEntries()))

    for i in range(n_entries):
        tree.GetEntry(i)
        x_val = getattr(tree, x_branch)
        y_raw = getattr(tree, y_branch)
        y_val = -y_raw if negate_y else y_raw
        hist.Fill(x_val, y_val)
    #endfor

    return hist
#enddef


def make_p2p_hists(tree, name_prefix, max_events=MAX_EVENTS):
    # You can tune binning/range if needed
    h_all = ROOT.TH1D(
        "{0}_all".format(name_prefix),
        "p2_p; p2_p (GeV); Counts",
        100, 0.0, 10.0
    )
    h_cut = ROOT.TH1D(
        "{0}_cut".format(name_prefix),
        "p2_p with cuts; p2_p (GeV); Counts",
        100, 0.0, 10.0
    )
    h_all.Sumw2()
    h_cut.Sumw2()

    n_entries = int(min(max_events, tree.GetEntries()))

    for i in range(n_entries):
        tree.GetEntry(i)

        p2_p = getattr(tree, "p2_p")
        t1 = getattr(tree, "t1")
        open_angle_ep2 = getattr(tree, "open_angle_ep2")
        pTmiss = getattr(tree, "pTmiss")
        y = getattr(tree, "y")
        W = getattr(tree, "W")
        Q2 = getattr(tree, "Q2")

        h_all.Fill(p2_p)

        # Cuts:
        # -t1 < 1
        # open_angle_ep2 < 5
        # pTmiss < 0.2
        # y < 0.8
        # W > 2
        # Q2 > 1
        if (-t1 < 1.0 and
            open_angle_ep2 < 5.0 and
            pTmiss < 0.2 and
            y < 0.8 and
            W > 2.0 and
            Q2 > 1.0):
            h_cut.Fill(p2_p)
        #endif
    #endfor

    return h_all, h_cut
#enddef


def sanitize_name(label):
    return label.replace(" ", "").replace("/", "_")
#enddef


# ----------------------------------------------------------------------
# Plot 1: y vs W (2D) for DVCS (top) and e p -> e p pi0 (bottom)
# ----------------------------------------------------------------------

def plot_y_vs_W(dvcs_trees, pi0_trees, output_dir):
    canvas = ROOT.TCanvas("c_y_vs_W", "y vs W (generated MC)", 2500, 900)
    canvas.Divide(5, 2)

    dvcs_hists = []
    pi0_hists = []

    for idx, period in enumerate(PERIODS):
        # DVCS row (top)
        pad_index_top = idx + 1
        canvas.cd(pad_index_top)

        hist_name_dvcs = "h_y_vs_W_dvcs_{0}".format(sanitize_name(period))
        title_dvcs = "DVCS {0}; W (GeV); y".format(period)
        h_dvcs = make_2d_hist(
            dvcs_trees[period],
            hist_name_dvcs,
            title_dvcs,
            x_branch="W",
            y_branch="y",
            x_bins=100,
            x_min=0.0,
            x_max=6.0,
            y_bins=100,
            y_min=0.0,
            y_max=1.0,
            max_events=MAX_EVENTS,
            negate_y=False
        )
        dvcs_hists.append(h_dvcs)
        h_dvcs.Draw("COLZ")

        # pi0 row (bottom)
        pad_index_bottom = idx + 1 + 5
        canvas.cd(pad_index_bottom)

        hist_name_pi0 = "h_y_vs_W_pi0_{0}".format(sanitize_name(period))
        title_pi0 = "e p -> e p #pi_{0} {0}; W (GeV); y".format(period)
        h_pi0 = make_2d_hist(
            pi0_trees[period],
            hist_name_pi0,
            title_pi0,
            x_branch="W",
            y_branch="y",
            x_bins=100,
            x_min=0.0,
            x_max=6.0,
            y_bins=100,
            y_min=0.0,
            y_max=1.0,
            max_events=MAX_EVENTS,
            negate_y=False
        )
        pi0_hists.append(h_pi0)
        h_pi0.Draw("COLZ")
    #endfor

    out_path = os.path.join(output_dir, "generated_yields_y_vs_W.png")
    canvas.SaveAs(out_path)
#enddef


# ----------------------------------------------------------------------
# Plot 2: y vs Q2 (2D)
# ----------------------------------------------------------------------

def plot_y_vs_Q2(dvcs_trees, pi0_trees, output_dir):
    canvas = ROOT.TCanvas("c_y_vs_Q2", "y vs Q2 (generated MC)", 2500, 900)
    canvas.Divide(5, 2)

    dvcs_hists = []
    pi0_hists = []

    for idx, period in enumerate(PERIODS):
        # DVCS row (top)
        pad_index_top = idx + 1
        canvas.cd(pad_index_top)

        hist_name_dvcs = "h_y_vs_Q2_dvcs_{0}".format(sanitize_name(period))
        title_dvcs = "DVCS {0}; Q^{{2}} (GeV^{{2}}); y".format(period)
        h_dvcs = make_2d_hist(
            dvcs_trees[period],
            hist_name_dvcs,
            title_dvcs,
            x_branch="Q2",
            y_branch="y",
            x_bins=100,
            x_min=0.0,
            x_max=8.0,
            y_bins=100,
            y_min=0.0,
            y_max=1.0,
            max_events=MAX_EVENTS,
            negate_y=False
        )
        dvcs_hists.append(h_dvcs)
        h_dvcs.Draw("COLZ")

        # pi0 row (bottom)
        pad_index_bottom = idx + 1 + 5
        canvas.cd(pad_index_bottom)

        hist_name_pi0 = "h_y_vs_Q2_pi0_{0}".format(sanitize_name(period))
        title_pi0 = "e p -> e p #pi_{0} {0}; Q^{{2}} (GeV^{{2}}); y".format(period)
        h_pi0 = make_2d_hist(
            pi0_trees[period],
            hist_name_pi0,
            title_pi0,
            x_branch="Q2",
            y_branch="y",
            x_bins=100,
            x_min=0.0,
            x_max=8.0,
            y_bins=100,
            y_min=0.0,
            y_max=1.0,
            max_events=MAX_EVENTS,
            negate_y=False
        )
        pi0_hists.append(h_pi0)
        h_pi0.Draw("COLZ")
    #endfor

    out_path = os.path.join(output_dir, "generated_yields_y_vs_Q2.png")
    canvas.SaveAs(out_path)
#enddef


# ----------------------------------------------------------------------
# Plot 3: -t1 vs Q2 (2D), y-axis is -t1, x-axis is Q2
# ----------------------------------------------------------------------

def plot_minust_vs_Q2(dvcs_trees, pi0_trees, output_dir):
    canvas = ROOT.TCanvas("c_minust_vs_Q2", "-t1 vs Q2 (generated MC)", 2500, 900)
    canvas.Divide(5, 2)

    dvcs_hists = []
    pi0_hists = []

    for idx, period in enumerate(PERIODS):
        # DVCS row (top)
        pad_index_top = idx + 1
        canvas.cd(pad_index_top)

        hist_name_dvcs = "h_minust_vs_Q2_dvcs_{0}".format(sanitize_name(period))
        title_dvcs = "DVCS {0}; Q^{{2}} (GeV^{{2}}); -t1 (GeV^{{2}})".format(period)
        h_dvcs = make_2d_hist(
            dvcs_trees[period],
            hist_name_dvcs,
            title_dvcs,
            x_branch="Q2",
            y_branch="t1",
            x_bins=100,
            x_min=0.0,
            x_max=8.0,
            y_bins=100,
            y_min=-1.0,
            y_max=3.0,
            max_events=MAX_EVENTS,
            negate_y=True
        )
        dvcs_hists.append(h_dvcs)
        h_dvcs.Draw("COLZ")

        # pi0 row (bottom)
        pad_index_bottom = idx + 1 + 5
        canvas.cd(pad_index_bottom)

        hist_name_pi0 = "h_minust_vs_Q2_pi0_{0}".format(sanitize_name(period))
        title_pi0 = "e p -> e p #pi_{0} {0}; Q^{{2}} (GeV^{{2}}); -t1 (GeV^{{2}})".format(period)
        h_pi0 = make_2d_hist(
            pi0_trees[period],
            hist_name_pi0,
            title_pi0,
            x_branch="Q2",
            y_branch="t1",
            x_bins=100,
            x_min=0.0,
            x_max=8.0,
            y_bins=100,
            y_min=-1.0,
            y_max=3.0,
            max_events=MAX_EVENTS,
            negate_y=True
        )
        pi0_hists.append(h_pi0)
        h_pi0.Draw("COLZ")
    #endfor

    out_path = os.path.join(output_dir, "generated_yields_minust_vs_Q2.png")
    canvas.SaveAs(out_path)
#enddef


# ----------------------------------------------------------------------
# Plot 4: p2_p distributions before/after cuts, 2x5
# ----------------------------------------------------------------------

def plot_p2p_before_after(dvcs_trees, pi0_trees, output_dir):
    canvas = ROOT.TCanvas("c_p2p", "p2_p before/after cuts (generated MC)", 2500, 900)
    canvas.Divide(5, 2)

    dvcs_all_hists = []
    dvcs_cut_hists = []
    pi0_all_hists = []
    pi0_cut_hists = []

    for idx, period in enumerate(PERIODS):
        # DVCS row (top)
        pad_index_top = idx + 1
        canvas.cd(pad_index_top)

        name_prefix_dvcs = "h_p2p_dvcs_{0}".format(sanitize_name(period))
        h_all_dvcs, h_cut_dvcs = make_p2p_hists(dvcs_trees[period], name_prefix_dvcs)
        dvcs_all_hists.append(h_all_dvcs)
        dvcs_cut_hists.append(h_cut_dvcs)

        h_all_dvcs.SetLineColor(ROOT.kBlack)
        h_all_dvcs.SetLineWidth(2)
        h_cut_dvcs.SetLineColor(ROOT.kRed)
        h_cut_dvcs.SetLineWidth(2)

        title = "DVCS {0}; p2_p (GeV); Counts".format(period)
        h_all_dvcs.SetTitle(title)

        max_y = max(h_all_dvcs.GetMaximum(), h_cut_dvcs.GetMaximum())
        if max_y <= 0.0:
            max_y = 1.0
        #endif
        h_all_dvcs.SetMaximum(1.2 * max_y)

        h_all_dvcs.Draw("HIST")
        h_cut_dvcs.Draw("HIST SAME")

        leg_top = ROOT.TLegend(0.60, 0.70, 0.90, 0.88)
        leg_top.SetBorderSize(0)
        leg_top.SetFillStyle(0)
        leg_top.AddEntry(h_all_dvcs, "Before cuts", "l")
        leg_top.AddEntry(h_cut_dvcs, "After cuts", "l")
        leg_top.Draw()

        # pi0 row (bottom)
        pad_index_bottom = idx + 1 + 5
        canvas.cd(pad_index_bottom)

        name_prefix_pi0 = "h_p2p_pi0_{0}".format(sanitize_name(period))
        h_all_pi0, h_cut_pi0 = make_p2p_hists(pi0_trees[period], name_prefix_pi0)
        pi0_all_hists.append(h_all_pi0)
        pi0_cut_hists.append(h_cut_pi0)

        h_all_pi0.SetLineColor(ROOT.kBlack)
        h_all_pi0.SetLineWidth(2)
        h_cut_pi0.SetLineColor(ROOT.kRed)
        h_cut_pi0.SetLineWidth(2)

        title_pi0 = "e p -> e p #pi_{0} {0}; p2_p (GeV); Counts".format(period)
        h_all_pi0.SetTitle(title_pi0)

        max_y_pi0 = max(h_all_pi0.GetMaximum(), h_cut_pi0.GetMaximum())
        if max_y_pi0 <= 0.0:
            max_y_pi0 = 1.0
        #endif
        h_all_pi0.SetMaximum(1.2 * max_y_pi0)

        h_all_pi0.Draw("HIST")
        h_cut_pi0.Draw("HIST SAME")

        leg_bottom = ROOT.TLegend(0.60, 0.70, 0.90, 0.88)
        leg_bottom.SetBorderSize(0)
        leg_bottom.SetFillStyle(0)
        leg_bottom.AddEntry(h_all_pi0, "Before cuts", "l")
        leg_bottom.AddEntry(h_cut_pi0, "After cuts", "l")
        leg_bottom.Draw()
    #endfor

    out_path = os.path.join(output_dir, "generated_yields_p2p_before_after_cuts.png")
    canvas.SaveAs(out_path)
#enddef


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():
    output_dir = "output/generated_yields_check"
    os.makedirs(output_dir, exist_ok=True)

    # Load trees
    dvcs_files, dvcs_trees = load_trees(DVCS_FILES)
    pi0_files, pi0_trees = load_trees(PI0_FILES)

    # Make plots
    plot_y_vs_W(dvcs_trees, pi0_trees, output_dir)
    plot_y_vs_Q2(dvcs_trees, pi0_trees, output_dir)
    plot_minust_vs_Q2(dvcs_trees, pi0_trees, output_dir)
    plot_p2p_before_after(dvcs_trees, pi0_trees, output_dir)

    print("All plots written to {0}".format(output_dir))
#enddef


if __name__ == "__main__":
    main()
#endif