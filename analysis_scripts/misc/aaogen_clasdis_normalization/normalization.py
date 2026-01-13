#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
yields.py

Inputs:
  --data    <data.root>
  --aaogen  <aaogen.root>
  --clasdis <clasdis.root>

All files must contain a TTree named "PhysicsEvents" with branches:
  - x
  - tprime
  - Mx2

Binning:
  xB edges     = [0.10, 0.25, 0.35, 0.45, 0.60]  -> 4 rows
  -t' edges    = [0.05, 0.25, 0.45, 0.65, 0.85, 1.05, 1.25] -> 6 columns
We bin using: (-tprime) in those edges.

In each (xB, -t') bin, we plot the Mx2 distribution (0 to 4) for:
  data (black), aaogen (red), clasdis (blue),
each normalized to unit area.

Output:
  output/yields.png
"""

import os
import sys
import argparse

import ROOT


TREE_NAME = "PhysicsEvents"

XB_EDGES  = [0.10, 0.25, 0.35, 0.45, 0.60]  # 4 bins -> 4 rows
TNEG_EDGES = [0.05, 0.25, 0.45, 0.65, 0.85, 1.05, 1.25]  # 6 bins -> 6 cols

MX2_MIN   = 0.0
MX2_MAX   = 4.0
MX2_NBINS = 200

OUTPUT_PNG = "output/yields.png"


def fatal(msg):
    raise RuntimeError(msg)


def require_file(path):
    if path is None or path.strip() == "":
        fatal("Missing required input path.")
    #endif
    if not os.path.isfile(path):
        fatal("File not found: " + path)
    #endif


def open_tree(path, tree_name):
    f = ROOT.TFile.Open(path, "READ")
    if not f or f.IsZombie():
        fatal("Failed to open ROOT file: " + path)
    #endif

    t = f.Get(tree_name)
    if not t:
        fatal("Tree '" + tree_name + "' not found in: " + path)
    #endif

    return f, t


def require_branches(t, needed, label):
    blist = t.GetListOfBranches()
    if not blist:
        fatal("No branch list available for tree in: " + label)
    #endif
    missing = []
    for b in needed:
        if not blist.FindObject(b):
            missing.append(b)
        #endif
    #endfor
    if len(missing) > 0:
        fatal("Missing required branches in " + label + ": " + ", ".join(missing))
    #endif


def normalize_unit_area(h):
    integ = h.Integral(0, h.GetNbinsX() + 1)  # include under/overflow
    if integ > 0.0:
        h.Scale(1.0 / integ)
    #endif


def style_hist(h, color, width, linestyle):
    h.SetLineColor(color)
    h.SetLineWidth(width)
    h.SetLineStyle(linestyle)
    h.SetMarkerSize(0.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True, help="Path to data ROOT file")
    ap.add_argument("--aaogen", required=True, help="Path to aaogen ROOT file")
    ap.add_argument("--clasdis", required=True, help="Path to clasdis ROOT file")
    args = ap.parse_args()

    require_file(args.data)
    require_file(args.aaogen)
    require_file(args.clasdis)

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)

    f_data, t_data = open_tree(args.data, TREE_NAME)
    f_aao,  t_aao  = open_tree(args.aaogen, TREE_NAME)
    f_dis,  t_dis  = open_tree(args.clasdis, TREE_NAME)

    needed = ["x", "tprime", "Mx2"]
    require_branches(t_data, needed, "data")
    require_branches(t_aao,  needed, "aaogen")
    require_branches(t_dis,  needed, "clasdis")

    outdir = os.path.dirname(OUTPUT_PNG)
    if outdir != "":
        os.makedirs(outdir, exist_ok=True)
    #endif

    nrows = len(XB_EDGES) - 1
    ncols = len(TNEG_EDGES) - 1

    c = ROOT.TCanvas("c_yields", "Mx2 yields in xB and -tprime bins", 2400, 1400)
    c.Divide(ncols, nrows, 0.001, 0.001)

    # Colors: data black, aaogen red, clasdis blue
    col_data = ROOT.kBlack
    col_aao  = ROOT.kRed
    col_dis  = ROOT.kBlue

    # Loop over bins: row-major (ROOT pads are 1-indexed, left-to-right then top-to-bottom)
    pad_idx = 1
    for ixb in range(nrows):
        xb_lo = XB_EDGES[ixb]
        xb_hi = XB_EDGES[ixb + 1]

        for it in range(ncols):
            t_lo = TNEG_EDGES[it]
            t_hi = TNEG_EDGES[it + 1]

            pad = c.cd(pad_idx)
            pad_idx += 1

            pad.SetGrid(1, 1)
            pad.SetLeftMargin(0.13)
            pad.SetRightMargin(0.04)
            pad.SetBottomMargin(0.14)
            pad.SetTopMargin(0.08)

            # Cut uses (-tprime) because your tprime is (t - tmin), typically negative.
            cut = (
                f"(x >= {xb_lo} && x < {xb_hi})"
                f" && (-tprime >= {t_lo} && -tprime < {t_hi})"
            )

            # Unique histogram names per pad
            hname_data = f"h_mx2_data_{ixb}_{it}"
            hname_aao  = f"h_mx2_aao_{ixb}_{it}"
            hname_dis  = f"h_mx2_dis_{ixb}_{it}"

            h_data = ROOT.TH1F(hname_data, "", MX2_NBINS, MX2_MIN, MX2_MAX)
            h_aao  = ROOT.TH1F(hname_aao,  "", MX2_NBINS, MX2_MIN, MX2_MAX)
            h_dis  = ROOT.TH1F(hname_dis,  "", MX2_NBINS, MX2_MIN, MX2_MAX)

            h_data.Sumw2()
            h_aao.Sumw2()
            h_dis.Sumw2()

            # Fill (goff so nothing is drawn yet)
            n_data = t_data.Draw(f"Mx2>>{hname_data}", cut, "goff")
            n_aao  = t_aao.Draw (f"Mx2>>{hname_aao}",  cut, "goff")
            n_dis  = t_dis.Draw (f"Mx2>>{hname_dis}",  cut, "goff")

            # Normalize each to unit area (per your request)
            normalize_unit_area(h_data)
            normalize_unit_area(h_aao)
            normalize_unit_area(h_dis)

            # Styling (use different linestyles so overlap is still visible)
            style_hist(h_data, col_data, 2, 1)  # solid
            style_hist(h_aao,  col_aao,  2, 2)  # dashed
            style_hist(h_dis,  col_dis,  2, 3)  # dotted

            # Axes
            h_data.GetXaxis().SetTitle("Mx2 (GeV^2)")
            h_data.GetYaxis().SetTitle("Normalized yield")
            h_data.GetXaxis().SetTitleSize(0.06)
            h_data.GetYaxis().SetTitleSize(0.06)
            h_data.GetXaxis().SetLabelSize(0.05)
            h_data.GetYaxis().SetLabelSize(0.05)

            # Common y-max so all three are visible
            ymax = max(h_data.GetMaximum(), h_aao.GetMaximum(), h_dis.GetMaximum())
            if ymax <= 0.0:
                ymax = 1.0
            #endif
            h_data.SetMaximum(1.25 * ymax)

            # Draw
            h_data.Draw("hist")
            h_aao.Draw("hist same")
            h_dis.Draw("hist same")

            # Legend (in every pad, for clarity)
            leg = ROOT.TLegend(0.55, 0.68, 0.94, 0.92)
            leg.SetBorderSize(1)
            leg.SetFillStyle(1001)
            leg.SetFillColor(ROOT.kWhite)
            leg.SetTextSize(0.045)
            leg.AddEntry(h_data, f"data (N={int(n_data)})", "l")
            leg.AddEntry(h_aao,  f"aaogen (N={int(n_aao)})", "l")
            leg.AddEntry(h_dis,  f"clasdis (N={int(n_dis)})", "l")
            leg.Draw()

            # Bin label
            tex = ROOT.TLatex()
            tex.SetNDC(True)
            tex.SetTextSize(0.05)
            tex.DrawLatex(0.14, 0.93, f"xB [{xb_lo:.2f}, {xb_hi:.2f})   -t' [{t_lo:.2f}, {t_hi:.2f})")
        #endfor
    #endfor

    c.SaveAs(OUTPUT_PNG)

    f_data.Close()
    f_aao.Close()
    f_dis.Close()
#enddef


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.stderr.write("FATAL: " + str(e) + "\n")
        sys.exit(1)
    #endif
#endif