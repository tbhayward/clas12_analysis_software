#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
normalization.py

Make a 4x6 canvas of Mx2 distributions binned in (xB, -t') where -t' is computed
on-the-fly from electron and pion kinematics (no recoil neutron detected).

Inputs:
  --data    <data.root>
  --aaogen  <aaogen.root>
  --clasdis <clasdis.root>

Tree:
  PhysicsEvents

Required branches (must exist in each file):
  runnum, e_p, e_theta, e_phi, p_p, p_theta, p_phi, x, Q2, Mx2

Binning:
  xB edges   = [0.10, 0.25, 0.35, 0.45, 0.60]  -> 4 rows
  -t' edges  = [0.05, 0.25, 0.45, 0.65, 0.85, 1.05, 1.25] -> 6 cols

Plot:
  Mx2 in [0, 4], normalized to unit area for each dataset in each bin.

Output:
  output/yields.png
"""

import os
import sys
import math
import argparse

import ROOT


TREE_NAME = "PhysicsEvents"

XB_EDGES   = [0.10, 0.25, 0.35, 0.45, 0.60]
TNEG_EDGES = [0.05, 0.25, 0.45, 0.65, 0.85, 1.05, 1.25]

MX2_MIN   = 0.0
MX2_MAX   = 4.0
MX2_NBINS = 200

OUTPUT_PNG = "output/yields.png"

# Masses (GeV)
MASS_E  = 0.000511
MASS_PI = 0.139570
MASS_N  = 0.9382720813


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


def beam_energy(runnum):
    # Matches your C++ beamEnergy(int runnum)
    if runnum >= 6616 and runnum <= 6783:
        return 10.1998
    if runnum >= 16042 and runnum <= 17065:
        return 10.5473
    if runnum >= 17067 and runnum <= 17724:
        return 10.5563
    if runnum >= 17725 and runnum <= 17811:
        return 10.5593
    return 10.5563
#enddef


def compute_t_scalar(runnum, e_p, e_theta, e_phi, p_p, p_theta, p_phi):
    Eb = beam_energy(int(runnum))
    if Eb <= 0.0:
        return 1.0e9
    #endif

    E_e = math.sqrt(max(0.0, e_p * e_p + MASS_E * MASS_E))
    se  = math.sin(e_theta)
    ce  = math.cos(e_theta)
    ex  = e_p * se * math.cos(e_phi)
    ey  = e_p * se * math.sin(e_phi)
    ez  = e_p * ce

    E_pi = math.sqrt(max(0.0, p_p * p_p + MASS_PI * MASS_PI))
    sp   = math.sin(p_theta)
    cp   = math.cos(p_theta)
    px   = p_p * sp * math.cos(p_phi)
    py   = p_p * sp * math.sin(p_phi)
    pz   = p_p * cp

    # q = k - k' = (Eb - E_e, -ex, -ey, Eb - ez)
    # then (q - p_pi)^2 = ( (Eb - E_e) - E_pi )^2 - | (-e_vec + (0,0,Eb)) - p_vec |^2
    dE = (Eb - E_e) - E_pi
    dx = -ex - px
    dy = -ey - py
    dz = (Eb - ez) - pz

    return dE * dE - (dx * dx + dy * dy + dz * dz)
#enddef


def compute_tmin_exact(xB, Q2):
    xb_ok = (xB > 0.0 and xB < 1.0)
    if Q2 <= 0.0 or (not xb_ok):
        if xb_ok:
            denom = (1.0 - xB)
            if denom > 0.0:
                return - (MASS_N * xB) * (MASS_N * xB) / denom
            #endif
        #endif
        return 0.0
    #endif

    eps2 = 4.0 * MASS_N * MASS_N * xB * xB / Q2
    root = math.sqrt(1.0 + eps2)
    num  = Q2 * (2.0 * (1.0 - xB) * (1.0 - root) + eps2)
    den  = 4.0 * xB * (1.0 - xB) + eps2
    if den == 0.0:
        return 0.0
    #endif
    return - num / den
#enddef


def normalize_unit_area(h):
    integ = h.Integral(0, h.GetNbinsX() + 1)
    if integ > 0.0:
        h.Scale(1.0 / integ)
    #endif
#enddef


def style_hist(h, color, width, linestyle):
    h.SetLineColor(color)
    h.SetLineWidth(width)
    h.SetLineStyle(linestyle)
    h.SetMarkerSize(0.0)
#enddef


def fill_mx2_hists_in_bin(tree, xb_lo, xb_hi, tneg_lo, tneg_hi, hist):
    """
    Loop over the tree, compute -tprime, and fill Mx2 into `hist`
    for events in:
      xb_lo <= x < xb_hi
      tneg_lo <= (-tprime) < tneg_hi
    """
    # Bind branch readers
    runnum  = ROOT.Long64_t(0)
    e_p     = ROOT.Double(0.0)
    e_theta = ROOT.Double(0.0)
    e_phi   = ROOT.Double(0.0)
    p_p     = ROOT.Double(0.0)
    p_theta = ROOT.Double(0.0)
    p_phi   = ROOT.Double(0.0)
    xB      = ROOT.Double(0.0)
    Q2      = ROOT.Double(0.0)
    Mx2     = ROOT.Double(0.0)

    tree.SetBranchStatus("*", 0)
    for bn in ["runnum", "e_p", "e_theta", "e_phi", "p_p", "p_theta", "p_phi", "x", "Q2", "Mx2"]:
        tree.SetBranchStatus(bn, 1)
    #endfor

    tree.SetBranchAddress("runnum",  runnum)
    tree.SetBranchAddress("e_p",     e_p)
    tree.SetBranchAddress("e_theta", e_theta)
    tree.SetBranchAddress("e_phi",   e_phi)
    tree.SetBranchAddress("p_p",     p_p)
    tree.SetBranchAddress("p_theta", p_theta)
    tree.SetBranchAddress("p_phi",   p_phi)
    tree.SetBranchAddress("x",       xB)
    tree.SetBranchAddress("Q2",      Q2)
    tree.SetBranchAddress("Mx2",     Mx2)

    n = tree.GetEntries()
    n_filled = 0

    for i in range(int(n)):
        tree.GetEntry(i)

        xb_val = float(xB)
        if xb_val < xb_lo or xb_val >= xb_hi:
            continue
        #endif

        Q2_val = float(Q2)
        t_val = compute_t_scalar(int(runnum), float(e_p), float(e_theta), float(e_phi),
                                 float(p_p), float(p_theta), float(p_phi))
        tmin_val = compute_tmin_exact(xb_val, Q2_val)
        tprime = t_val - tmin_val
        tneg = -tprime

        if tneg < tneg_lo or tneg >= tneg_hi:
            continue
        #endif

        hist.Fill(float(Mx2))
        n_filled += 1
    #endfor

    tree.SetBranchStatus("*", 1)

    return n_filled
#enddef


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

    needed = ["runnum", "e_p", "e_theta", "e_phi", "p_p", "p_theta", "p_phi", "x", "Q2", "Mx2"]
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

    col_data = ROOT.kBlack
    col_aao  = ROOT.kRed
    col_dis  = ROOT.kBlue

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

            hname_data = f"h_mx2_data_{ixb}_{it}"
            hname_aao  = f"h_mx2_aao_{ixb}_{it}"
            hname_dis  = f"h_mx2_dis_{ixb}_{it}"

            h_data = ROOT.TH1F(hname_data, "", MX2_NBINS, MX2_MIN, MX2_MAX)
            h_aao  = ROOT.TH1F(hname_aao,  "", MX2_NBINS, MX2_MIN, MX2_MAX)
            h_dis  = ROOT.TH1F(hname_dis,  "", MX2_NBINS, MX2_MIN, MX2_MAX)

            h_data.Sumw2()
            h_aao.Sumw2()
            h_dis.Sumw2()

            n_data = fill_mx2_hists_in_bin(t_data, xb_lo, xb_hi, t_lo, t_hi, h_data)
            n_aao  = fill_mx2_hists_in_bin(t_aao,  xb_lo, xb_hi, t_lo, t_hi, h_aao)
            n_dis  = fill_mx2_hists_in_bin(t_dis,  xb_lo, xb_hi, t_lo, t_hi, h_dis)

            normalize_unit_area(h_data)
            normalize_unit_area(h_aao)
            normalize_unit_area(h_dis)

            style_hist(h_data, col_data, 2, 1)  # solid
            style_hist(h_aao,  col_aao,  2, 2)  # dashed
            style_hist(h_dis,  col_dis,  2, 3)  # dotted

            h_data.GetXaxis().SetTitle("Mx2 (GeV^2)")
            h_data.GetYaxis().SetTitle("Normalized yield")
            h_data.GetXaxis().SetTitleSize(0.06)
            h_data.GetYaxis().SetTitleSize(0.06)
            h_data.GetXaxis().SetLabelSize(0.05)
            h_data.GetYaxis().SetLabelSize(0.05)

            ymax = max(h_data.GetMaximum(), h_aao.GetMaximum(), h_dis.GetMaximum())
            if ymax <= 0.0:
                ymax = 1.0
            #endif
            h_data.SetMaximum(1.25 * ymax)

            h_data.Draw("hist")
            h_aao.Draw("hist same")
            h_dis.Draw("hist same")

            leg = ROOT.TLegend(0.55, 0.68, 0.94, 0.92)
            leg.SetBorderSize(1)
            leg.SetFillStyle(1001)
            leg.SetFillColor(ROOT.kWhite)
            leg.SetTextSize(0.045)
            leg.AddEntry(h_data, f"data (N={int(n_data)})", "l")
            leg.AddEntry(h_aao,  f"aaogen (N={int(n_aao)})", "l")
            leg.AddEntry(h_dis,  f"clasdis (N={int(n_dis)})", "l")
            leg.Draw()

            tex = ROOT.TLatex()
            tex.SetNDC(True)
            tex.SetTextSize(0.05)
            tex.DrawLatex(0.14, 0.93,
                          f"xB [{xb_lo:.2f}, {xb_hi:.2f})   -t' [{t_lo:.2f}, {t_hi:.2f})")
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