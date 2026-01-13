#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
fermi_motion_double_check.py

Make a 2x2 canvas comparing delta variables with/without Fermi motion smearing.

Layout:
  (row 0, col 0): clasdis,  x    - mc_x
  (row 0, col 1): aaogen,   x    - mc_x
  (row 1, col 0): clasdis,  p_p  - mc_p_p
  (row 1, col 1): aaogen,   p_p  - mc_p_p

Each pad overlays:
  - Fermi motion (red)
  - No motion    (black)

Saves: output/fermi_motion_double_check.png
"""

import os
import sys
import argparse

import ROOT


DEFAULT_CLASDIS_FERMI = "/volatile/clas12/thayward/fermi_motion_study/clasdis_fermi_motion.root"
DEFAULT_CLASDIS_NO    = "/volatile/clas12/thayward/fermi_motion_study/clasdis_no_motion.root"
DEFAULT_AAOGEN_FERMI  = "/volatile/clas12/thayward/fermi_motion_study/aaogen_fermi_motion.root"
DEFAULT_AAOGEN_NO     = "/volatile/clas12/thayward/fermi_motion_study/aaogen_no_motion.root"

TREE_NAME = "PhysicsEvents"

# Histogram configuration (explicit and deterministic)
NBINS_DX   = 240
DX_MIN    = -0.30
DX_MAX    =  0.30

NBINS_DP   = 240
DP_MIN    = -1.50
DP_MAX    =  1.50

# If True, scale each histogram to unit area (shape comparison).
# If False, keep raw counts.
NORMALIZE = True

# If True, set log-y in each pad.
DO_LOGY = False

OUTPUT_PNG = "output/fermi_motion_double_check.png"


def fatal(msg):
    raise RuntimeError(msg)


def require_file(path):
    if not os.path.isfile(path):
        fatal("Missing file: " + path)
    #endif


def open_tree(path, tree_name):
    f = ROOT.TFile.Open(path, "READ")
    if not f or f.IsZombie():
        fatal("Failed to open ROOT file: " + path)
    #endif

    t = f.Get(tree_name)
    if not t:
        fatal("Tree '" + tree_name + "' not found in file: " + path)
    #endif

    # Keep file handle alive by returning it too.
    return f, t


def require_branches(tree, branches, file_label):
    for b in branches:
        if not tree.GetListOfBranches().FindObject(b):
            fatal("Missing branch '" + b + "' in " + file_label)
        #endif
    #endfor


def make_hist(name, title, nbins, xmin, xmax):
    h = ROOT.TH1F(name, title, nbins, xmin, xmax)
    h.Sumw2()
    return h


def fill_hist(tree, expr, hist, max_entries):
    draw_opt = "goff"
    if max_entries >= 0:
        tree.SetEstimate(max_entries)
    #endif

    if max_entries >= 0:
        n = tree.Draw(expr + ">>" + hist.GetName(), "", draw_opt, max_entries, 0)
    else:
        n = tree.Draw(expr + ">>" + hist.GetName(), "", draw_opt)
    #endif

    if n <= 0:
        fatal("No entries filled for expression '" + expr + "' into hist '" + hist.GetName() + "'.")
    #endif

    if NORMALIZE:
        integral = hist.Integral(0, hist.GetNbinsX() + 1)
        if integral <= 0.0:
            fatal("Histogram '" + hist.GetName() + "' has non-positive integral; cannot normalize.")
        #endif
        hist.Scale(1.0 / integral)
    #endif


def style_hist(h, color, width):
    h.SetLineColor(color)
    h.SetLineWidth(width)
    h.SetMarkerColor(color)
    h.SetMarkerSize(0.0)


def draw_overlay_pad(pad, h_no, h_fermi, xlabel, ylabel, pad_title):
    pad.cd()
    pad.SetGrid(1, 1)
    if DO_LOGY:
        pad.SetLogy(1)
    else:
        pad.SetLogy(0)
    #endif

    h_no.GetXaxis().SetTitle(xlabel)
    h_no.GetYaxis().SetTitle(ylabel)

    max_no = h_no.GetMaximum()
    max_fm = h_fermi.GetMaximum()
    ymax = max(max_no, max_fm)
    h_no.SetMaximum(1.25 * ymax)

    h_no.Draw("hist")
    h_fermi.Draw("hist same")

    leg = ROOT.TLegend(0.58, 0.74, 0.90, 0.90)
    leg.SetBorderSize(1)
    leg.SetFillStyle(1001)
    leg.SetFillColor(ROOT.kWhite)
    leg.AddEntry(h_no, "no motion", "l")
    leg.AddEntry(h_fermi, "fermi motion", "l")
    leg.Draw()

    tex = ROOT.TLatex()
    tex.SetNDC(True)
    tex.SetTextSize(0.045)
    tex.DrawLatex(0.12, 0.92, pad_title)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clasdis-fermi", default=DEFAULT_CLASDIS_FERMI)
    ap.add_argument("--clasdis-no",    default=DEFAULT_CLASDIS_NO)
    ap.add_argument("--aaogen-fermi",  default=DEFAULT_AAOGEN_FERMI)
    ap.add_argument("--aaogen-no",     default=DEFAULT_AAOGEN_NO)
    ap.add_argument("--max-entries",   type=int, default=-1,
                    help="If >= 0, process only the first N entries of each tree.")
    args = ap.parse_args()

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)

    require_file(args.clasdis_fermi)
    require_file(args.clasdis_no)
    require_file(args.aaogen_fermi)
    require_file(args.aaogen_no)

    f_cd_fm, t_cd_fm = open_tree(args.clasdis_fermi, TREE_NAME)
    f_cd_no, t_cd_no = open_tree(args.clasdis_no, TREE_NAME)
    f_aa_fm, t_aa_fm = open_tree(args.aaogen_fermi, TREE_NAME)
    f_aa_no, t_aa_no = open_tree(args.aaogen_no, TREE_NAME)

    need_x  = ["x", "mc_x"]
    need_pp = ["p_p", "mc_p_p"]

    require_branches(t_cd_fm, need_x + need_pp, "clasdis_fermi_motion")
    require_branches(t_cd_no, need_x + need_pp, "clasdis_no_motion")
    require_branches(t_aa_fm, need_x + need_pp, "aaogen_fermi_motion")
    require_branches(t_aa_no, need_x + need_pp, "aaogen_no_motion")

    outdir = os.path.dirname(OUTPUT_PNG)
    if outdir != "":
        os.makedirs(outdir, exist_ok=True)
    #endif

    # Top row: delta x
    h_cd_dx_no = make_hist("h_cd_dx_no", "clasdis: x - mc_x (no)", NBINS_DX, DX_MIN, DX_MAX)
    h_cd_dx_fm = make_hist("h_cd_dx_fm", "clasdis: x - mc_x (fm)", NBINS_DX, DX_MIN, DX_MAX)

    h_aa_dx_no = make_hist("h_aa_dx_no", "aaogen: x - mc_x (no)", NBINS_DX, DX_MIN, DX_MAX)
    h_aa_dx_fm = make_hist("h_aa_dx_fm", "aaogen: x - mc_x (fm)", NBINS_DX, DX_MIN, DX_MAX)

    # Bottom row: delta p_p
    h_cd_dp_no = make_hist("h_cd_dp_no", "clasdis: p_p - mc_p_p (no)", NBINS_DP, DP_MIN, DP_MAX)
    h_cd_dp_fm = make_hist("h_cd_dp_fm", "clasdis: p_p - mc_p_p (fm)", NBINS_DP, DP_MIN, DP_MAX)

    h_aa_dp_no = make_hist("h_aa_dp_no", "aaogen: p_p - mc_p_p (no)", NBINS_DP, DP_MIN, DP_MAX)
    h_aa_dp_fm = make_hist("h_aa_dp_fm", "aaogen: p_p - mc_p_p (fm)", NBINS_DP, DP_MIN, DP_MAX)

    style_hist(h_cd_dx_no, ROOT.kBlack, 2)
    style_hist(h_cd_dx_fm, ROOT.kRed,   2)
    style_hist(h_aa_dx_no, ROOT.kBlack, 2)
    style_hist(h_aa_dx_fm, ROOT.kRed,   2)

    style_hist(h_cd_dp_no, ROOT.kBlack, 2)
    style_hist(h_cd_dp_fm, ROOT.kRed,   2)
    style_hist(h_aa_dp_no, ROOT.kBlack, 2)
    style_hist(h_aa_dp_fm, ROOT.kRed,   2)

    fill_hist(t_cd_no, "x - mc_x", h_cd_dx_no, args.max_entries)
    fill_hist(t_cd_fm, "x - mc_x", h_cd_dx_fm, args.max_entries)

    fill_hist(t_aa_no, "x - mc_x", h_aa_dx_no, args.max_entries)
    fill_hist(t_aa_fm, "x - mc_x", h_aa_dx_fm, args.max_entries)

    fill_hist(t_cd_no, "p_p - mc_p_p", h_cd_dp_no, args.max_entries)
    fill_hist(t_cd_fm, "p_p - mc_p_p", h_cd_dp_fm, args.max_entries)

    fill_hist(t_aa_no, "p_p - mc_p_p", h_aa_dp_no, args.max_entries)
    fill_hist(t_aa_fm, "p_p - mc_p_p", h_aa_dp_fm, args.max_entries)

    c = ROOT.TCanvas("c_fermi_motion_double_check", "Fermi motion double check", 1400, 1000)
    c.Divide(2, 2)

    ylab = "Counts (arb.)" if NORMALIZE else "Counts"

    draw_overlay_pad(c.cd(1), h_cd_dx_no, h_cd_dx_fm, "x - mc_x", ylab, "clasdis: x - mc_x")
    draw_overlay_pad(c.cd(2), h_aa_dx_no, h_aa_dx_fm, "x - mc_x", ylab, "aaogen:  x - mc_x")

    draw_overlay_pad(c.cd(3), h_cd_dp_no, h_cd_dp_fm, "p_p - mc_p_p (GeV)", ylab, "clasdis: p_p - mc_p_p")
    draw_overlay_pad(c.cd(4), h_aa_dp_no, h_aa_dp_fm, "p_p - mc_p_p (GeV)", ylab, "aaogen:  p_p - mc_p_p")

    c.SaveAs(OUTPUT_PNG)

    f_cd_fm.Close()
    f_cd_no.Close()
    f_aa_fm.Close()
    f_aa_no.Close()


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.stderr.write("FATAL: " + str(e) + "\n")
        sys.exit(1)
    #endif
#endif