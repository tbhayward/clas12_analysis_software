#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
mx2_twofile_hist.py

Usage:
  python mx2_twofile_hist.py <file1.root> <file2.root>

Reads TTree "PhysicsEvents" and branch "Mx2" from two ROOT files and overlays
their histograms on a single canvas.

Output:
  output/this_had_better_work_please_god_let_this_work_im_going_insane.png
"""

import os
import sys

import ROOT


def die(msg):
    print(f"FATAL: {msg}", file=sys.stderr)
    sys.exit(1)


def ensure_outdir(path):
    d = os.path.dirname(path)
    if d and (not os.path.isdir(d)):
        os.makedirs(d, exist_ok=True)
    #endif


def open_root_file(path):
    f = ROOT.TFile.Open(path, "READ")
    if (not f) or f.IsZombie():
        die(f"Could not open ROOT file: {path}")
    #endif
    return f


def get_tree(f, name):
    t = f.Get(name)
    if not t:
        die(f"Could not find TTree '{name}' in file: {f.GetName()}")
    #endif
    return t


def hist_from_tree(tree, hist_name, title, nbins, xmin, xmax):
    h = ROOT.TH1F(hist_name, title, nbins, xmin, xmax)
    h.Sumw2()
    h.SetDirectory(0)

    # Use Draw to avoid manual branch binding pitfalls.
    draw_expr = "Mx2"
    draw_cut = ""  # add cuts here if desired
    draw_opt = "goff"

    n_drawn = tree.Draw(f"{draw_expr}>>{hist_name}", draw_cut, draw_opt)
    if n_drawn < 0:
        die("TTree.Draw failed (returned < 0). Check that branch 'Mx2' exists and is readable.")
    #endif

    if h.GetEntries() <= 0:
        die("Histogram has zero entries. Check that 'Mx2' exists and that values fall in [0.4, 2.0].")
    #endif

    return h


def normalize_unit_area_in_range(h, xmin, xmax):
    bmin = h.FindBin(xmin + 1e-9)
    bmax = h.FindBin(xmax - 1e-9)
    integral = h.Integral(bmin, bmax)
    if integral > 0.0:
        h.Scale(1.0 / integral)
    else:
        die("Histogram integral over the plotted range is zero; cannot normalize.")
    #endif


def main():
    if len(sys.argv) != 3:
        die("Usage: python mx2_twofile_hist.py <file1.root> <file2.root>")
    #endif

    file1_path = sys.argv[1]
    file2_path = sys.argv[2]

    out_png = "output/this_had_better_work_please_god_let_this_work_im_going_insane.png"
    ensure_outdir(out_png)

    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.SetDefaultSumw2(True)

    f1 = open_root_file(file1_path)
    f2 = open_root_file(file2_path)

    t1 = get_tree(f1, "PhysicsEvents")
    t2 = get_tree(f2, "PhysicsEvents")

    nbins = 200
    xmin = 0.4
    xmax = 2.0

    h1 = hist_from_tree(t1, "h_mx2_1", "Mx2 comparison", nbins, xmin, xmax)
    h2 = hist_from_tree(t2, "h_mx2_2", "Mx2 comparison", nbins, xmin, xmax)

    normalize_unit_area_in_range(h1, xmin, xmax)
    normalize_unit_area_in_range(h2, xmin, xmax)

    # Style
    h1.SetLineWidth(2)
    h2.SetLineWidth(2)
    h1.SetLineColor(ROOT.kBlack)
    h2.SetLineColor(ROOT.kRed)

    h1.GetXaxis().SetTitle("M_{x}^{2} (GeV^{2})")
    h1.GetYaxis().SetTitle("Normalized counts")
    h1.GetXaxis().CenterTitle(True)
    h1.GetYaxis().CenterTitle(True)

    c = ROOT.TCanvas("c_mx2", "c_mx2", 900, 700)
    c.SetLeftMargin(0.14)
    c.SetRightMargin(0.05)
    c.SetBottomMargin(0.12)
    c.SetTopMargin(0.08)

    maxy = max(h1.GetMaximum(), h2.GetMaximum())
    h1.SetMaximum(1.15 * maxy)

    h1.Draw("HIST")
    h2.Draw("HIST SAME")

    leg = ROOT.TLegend(0.60, 0.75, 0.93, 0.92)
    leg.SetFillStyle(1001)
    leg.SetFillColor(ROOT.kWhite)
    leg.SetBorderSize(1)
    leg.AddEntry(h1, os.path.basename(file1_path), "l")
    leg.AddEntry(h2, os.path.basename(file2_path), "l")
    leg.Draw()

    c.SaveAs(out_png)

    # Clean up
    f1.Close()
    f2.Close()

    print(f"Wrote: {out_png}")


if __name__ == "__main__":
    main()
#endif