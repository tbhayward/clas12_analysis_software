#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
normalization.py

Phase 1:
  - Fill 4x6 Mx2 histograms binned in (xB, -tprime) for data, aaogen, clasdis
  - Compute t, tmin, tprime on-the-fly
  - Normalize each histogram to unit area per (xB, -tprime) bin (shape-only)

Phase 2:
  - Determine a per-bin mixture fraction w[r,c] in [0,1]:
        H_mix(r,c) = w[r,c] * H_aaogen(r,c) + (1-w[r,c]) * H_clasdis(r,c)

  - Choose w[r,c] by minimizing the unweighted sum of squared differences
    between the data shape and mixture shape within an Mx2 fit window:
        SSE_{r,c}(w) = sum_{i in window} ( D_i - (w*A_i + (1-w)*C_i) )^2

    Exact minimizer in each pad (over bins i in window):
        X_i = (A_i - C_i),  Y_i = (D_i - C_i)
        w = sum_i (X_i Y_i) / sum_i (X_i^2), clipped to [0,1].

Additional diagnostic output:
  - output/yields_data_only.png is now a SINGLE-PAD integrated plot:
        (sum over all xB and -tprime bins) for data vs aaogen vs clasdis.

Outputs:
  output/yields.png           : data vs aaogen vs clasdis (shape-only) per pad
  output/yields_mix.png       : data vs per-bin mixture per pad, legend shows w[r,c]
  output/yields_data_only.png : integrated (all bins summed) data vs aaogen vs clasdis
  output/weights.txt          : per-bin w[r,c], w_unclipped, SSE summary
"""

import os
import sys
import math
import argparse
from array import array

import ROOT

TREE_NAME = "PhysicsEvents"

XB_EDGES = [0.10, 0.25, 0.35, 0.45, 0.60]  # 4 rows
TNEG_EDGES = [0.05, 0.25, 0.45, 0.65, 0.85, 1.05, 1.25]  # 6 cols in -tprime

# Histogram window (identical to fit window)
MX2_MIN = 0.7
MX2_MAX = 1.1
MX2_NBINS = 25

# Fit window used for solving w and computing SSE (identical to hist range)
MX2_FIT_MIN = 0.7
MX2_FIT_MAX = 1.1

OUTPUT_YIELDS_PNG = "output/yields.png"
OUTPUT_MIX_PNG = "output/yields_mix.png"
OUTPUT_DATAONLY_PNG = "output/yields_data_only.png"  # now integrated 1-pad three-way
OUTPUT_WEIGHTS_TXT = "output/weights.txt"

# Masses (GeV)
MASS_E = 0.000511
MASS_PI = 0.139570
MASS_N = 0.9382720813


def fatal(msg):
    raise RuntimeError(msg)
#enddef


def require_file(path):
    if path is None or path.strip() == "":
        fatal("Missing required input path.")
    #endif
    if not os.path.isfile(path):
        fatal("File not found: " + path)
    #endif
#enddef


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
#enddef


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
#enddef


def beam_energy(runnum):
    if runnum >= 6616 and runnum <= 6783:
        return 10.1998
    #endif
    if runnum >= 16042 and runnum <= 17065:
        return 10.5473
    #endif
    if runnum >= 17067 and runnum <= 17724:
        return 10.5563
    #endif
    if runnum >= 17725 and runnum <= 17811:
        return 10.5593
    #endif
    return 10.5563
#enddef


def compute_t_scalar(runnum, e_p, e_theta, e_phi, p_p, p_theta, p_phi):
    Eb = beam_energy(int(runnum))
    if Eb <= 0.0:
        return 1.0e9
    #endif

    E_e = math.sqrt(max(0.0, e_p * e_p + MASS_E * MASS_E))
    se = math.sin(e_theta)
    ce = math.cos(e_theta)
    ex = e_p * se * math.cos(e_phi)
    ey = e_p * se * math.sin(e_phi)
    ez = e_p * ce

    E_pi = math.sqrt(max(0.0, p_p * p_p + MASS_PI * MASS_PI))
    sp = math.sin(p_theta)
    cp = math.cos(p_theta)
    px = p_p * sp * math.cos(p_phi)
    py = p_p * sp * math.sin(p_phi)
    pz = p_p * cp

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
    num = Q2 * (2.0 * (1.0 - xB) * (1.0 - root) + eps2)
    den = 4.0 * xB * (1.0 - xB) + eps2
    if den == 0.0:
        return 0.0
    #endif
    return - num / den
#enddef


def find_bin(val, edges):
    for i in range(len(edges) - 1):
        if val >= edges[i] and val < edges[i + 1]:
            return i
        #endif
    #endfor
    return -1
#enddef


def normalize_unit_area(h):
    integ = h.Integral(1, h.GetNbinsX())
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


def make_hist_grid(prefix, nrows, ncols):
    grid = []
    for r in range(nrows):
        row = []
        for c in range(ncols):
            name = f"{prefix}_r{r}_c{c}"
            h = ROOT.TH1F(name, "", MX2_NBINS, MX2_MIN, MX2_MAX)
            h.Sumw2()
            row.append(h)
        #endfor
        grid.append(row)
    #endfor
    return grid
#enddef


def make_count_grid(nrows, ncols):
    return [[0 for _ in range(ncols)] for _ in range(nrows)]
#enddef


def fill_all_bins_single_pass(tree, hgrid, cgrid, max_events):
    runnum = array("i", [0])
    e_p = array("d", [0.0])
    e_theta = array("d", [0.0])
    e_phi = array("d", [0.0])
    p_p = array("d", [0.0])
    p_theta = array("d", [0.0])
    p_phi = array("d", [0.0])
    xB = array("d", [0.0])
    Q2 = array("d", [0.0])
    Mx2 = array("d", [0.0])

    tree.SetBranchStatus("*", 0)
    for bn in ["runnum", "e_p", "e_theta", "e_phi", "p_p", "p_theta", "p_phi", "x", "Q2", "Mx2"]:
        tree.SetBranchStatus(bn, 1)
    #endfor

    tree.SetBranchAddress("runnum", runnum)
    tree.SetBranchAddress("e_p", e_p)
    tree.SetBranchAddress("e_theta", e_theta)
    tree.SetBranchAddress("e_phi", e_phi)
    tree.SetBranchAddress("p_p", p_p)
    tree.SetBranchAddress("p_theta", p_theta)
    tree.SetBranchAddress("p_phi", p_phi)
    tree.SetBranchAddress("x", xB)
    tree.SetBranchAddress("Q2", Q2)
    tree.SetBranchAddress("Mx2", Mx2)

    n_entries = int(tree.GetEntries())
    if max_events is None:
        n_to_process = n_entries
    else:
        n_to_process = min(n_entries, int(max_events))
    #endif

    for i in range(n_to_process):
        tree.GetEntry(i)

        xb_val = float(xB[0])
        rb = find_bin(xb_val, XB_EDGES)
        if rb < 0:
            continue
        #endif

        Q2_val = float(Q2[0])

        t_val = compute_t_scalar(int(runnum[0]),
                                 float(e_p[0]), float(e_theta[0]), float(e_phi[0]),
                                 float(p_p[0]), float(p_theta[0]), float(p_phi[0]))
        tmin_val = compute_tmin_exact(xb_val, Q2_val)
        tprime = t_val - tmin_val
        tneg = -tprime

        cb = find_bin(tneg, TNEG_EDGES)
        if cb < 0:
            continue
        #endif

        hgrid[rb][cb].Fill(float(Mx2[0]))
        cgrid[rb][cb] += 1
    #endfor

    tree.SetBranchStatus("*", 1)
#enddef


def ensure_outdir(path):
    d = os.path.dirname(path)
    if d != "":
        os.makedirs(d, exist_ok=True)
    #endif
#enddef


def compute_best_w_unweighted_for_pad(hd, ha, hc):
    """
    Minimizes SSE(w) = sum_{i in fit window} ( D_i - (w*A_i + (1-w)*C_i) )^2.

    Exact minimizer:
      X_i = (A_i - C_i)
      Y_i = (D_i - C_i)
      w_unclipped = sum(X_i*Y_i)/sum(X_i^2)   (over i in window)
      w = clip(w_unclipped, 0, 1)
    """

    if hd.Integral(1, hd.GetNbinsX()) <= 0.0:
        return 0.0, 0.0, 0.0
    #endif

    nb = hd.GetNbinsX()
    num = 0.0
    den = 0.0
    used_bins = 0

    for i in range(1, nb + 1):
        xcen = hd.GetXaxis().GetBinCenter(i)
        if xcen < MX2_FIT_MIN or xcen > MX2_FIT_MAX:
            continue
        #endif

        D = hd.GetBinContent(i)
        A = ha.GetBinContent(i)
        C = hc.GetBinContent(i)
        X = (A - C)
        Y = (D - C)
        num += X * Y
        den += X * X
        used_bins += 1
    #endfor

    if used_bins == 0 or den <= 0.0:
        w_unclipped = 0.0
    else:
        w_unclipped = num / den
    #endif

    w = w_unclipped
    if w < 0.0:
        w = 0.0
    #endif
    if w > 1.0:
        w = 1.0
    #endif

    sse = 0.0
    for i in range(1, nb + 1):
        xcen = hd.GetXaxis().GetBinCenter(i)
        if xcen < MX2_FIT_MIN or xcen > MX2_FIT_MAX:
            continue
        #endif
        D = hd.GetBinContent(i)
        M = w * ha.GetBinContent(i) + (1.0 - w) * hc.GetBinContent(i)
        diff = D - M
        sse += diff * diff
    #endfor

    return w, w_unclipped, sse
#enddef


def compute_w_grid_and_mix(h_data, h_aao, h_dis):
    nrows = len(h_data)
    ncols = len(h_data[0])

    w_grid = [[0.0 for _ in range(ncols)] for _ in range(nrows)]
    wun_grid = [[0.0 for _ in range(ncols)] for _ in range(nrows)]
    sse_grid = [[0.0 for _ in range(ncols)] for _ in range(nrows)]
    h_mix = []

    for r in range(nrows):
        row = []
        for c in range(ncols):
            hd = h_data[r][c]
            ha = h_aao[r][c]
            hc = h_dis[r][c]

            w, w_unclipped, sse = compute_best_w_unweighted_for_pad(hd, ha, hc)
            w_grid[r][c] = w
            wun_grid[r][c] = w_unclipped
            sse_grid[r][c] = sse

            hm = ha.Clone(f"h_mix_r{r}_c{c}")
            hm.Reset("ICESM")
            hm.Add(ha, w)
            hm.Add(hc, (1.0 - w))
            row.append(hm)
        #endfor
        h_mix.append(row)
    #endfor

    return w_grid, wun_grid, sse_grid, h_mix
#enddef


def pad_set_margins(pad):
    pad.SetGrid(1, 1)
    pad.SetLeftMargin(0.20)
    pad.SetRightMargin(0.04)
    pad.SetBottomMargin(0.14)
    pad.SetTopMargin(0.08)
#enddef


def set_axes_and_range(h_frame, ymax):
    h_frame.GetXaxis().SetTitle("Mx2 (GeV^2)")
    h_frame.GetYaxis().SetTitle("Normalized yield")
    h_frame.GetXaxis().SetTitleSize(0.06)
    h_frame.GetYaxis().SetTitleSize(0.06)
    h_frame.GetXaxis().SetLabelSize(0.05)
    h_frame.GetYaxis().SetLabelSize(0.05)

    h_frame.SetMinimum(0.0)
    h_frame.SetMaximum(ymax)
#enddef


def sum_counts_grid(cgrid):
    nrows = len(cgrid)
    ncols = len(cgrid[0])
    tot = 0
    for r in range(nrows):
        for c in range(ncols):
            tot += int(cgrid[r][c])
        #endfor
    #endfor
    return tot
#enddef


def make_integrated_hist_from_grid(hgrid, name):
    nrows = len(hgrid)
    ncols = len(hgrid[0])
    hsum = hgrid[0][0].Clone(name)
    hsum.Reset("ICESM")
    for r in range(nrows):
        for c in range(ncols):
            hsum.Add(hgrid[r][c])
        #endfor
    #endfor
    return hsum
#enddef


def draw_canvas_integrated_threeway(hd_int, ha_int, hc_int, Nd, Na, Nc, outpng):
    canv = ROOT.TCanvas("c_integrated", "Mx2 integrated: data vs aaogen vs clasdis", 1200, 900)
    canv.cd()

    pad = canv.cd(1)
    pad.SetGrid(1, 1)
    pad.SetLeftMargin(0.18)
    pad.SetRightMargin(0.05)
    pad.SetBottomMargin(0.14)
    pad.SetTopMargin(0.08)

    ymax = max(hd_int.GetMaximum(), ha_int.GetMaximum(), hc_int.GetMaximum())
    if ymax <= 0.0:
        ymax = 1.0
    #endif
    ymax = 1.2 * ymax

    set_axes_and_range(hd_int, ymax)

    hd_int.Draw("hist")
    ha_int.Draw("hist same")
    hc_int.Draw("hist same")

    leg = ROOT.TLegend(0.55, 0.72, 0.94, 0.92)
    leg.SetBorderSize(1)
    leg.SetFillStyle(1001)
    leg.SetFillColor(ROOT.kWhite)
    leg.SetTextSize(0.042)
    leg.AddEntry(hd_int, f"data (N={int(Nd)})", "l")
    leg.AddEntry(ha_int, f"aaogen (N={int(Na)})", "l")
    leg.AddEntry(hc_int, f"clasdis (N={int(Nc)})", "l")
    leg.Draw()

    tex = ROOT.TLatex()
    tex.SetNDC(True)
    tex.SetTextSize(0.045)
    tex.DrawLatex(0.16, 0.93, "Integrated over all xB and -tprime bins")

    canv.SaveAs(outpng)
#enddef


def draw_canvas_threeway(h_data, h_aao, h_dis, c_data, c_aao, c_dis, outpng):
    nrows = len(h_data)
    ncols = len(h_data[0])

    canv = ROOT.TCanvas("c_yields", "Mx2 in xB and -tprime bins", 2400, 1400)
    canv.Divide(ncols, nrows, 0.001, 0.001)

    pad_idx = 1
    for r in range(nrows):
        xb_lo = XB_EDGES[r]
        xb_hi = XB_EDGES[r + 1]

        for c in range(ncols):
            t_lo = TNEG_EDGES[c]
            t_hi = TNEG_EDGES[c + 1]

            pad = canv.cd(pad_idx)
            pad_idx += 1
            pad_set_margins(pad)

            hd = h_data[r][c]
            ha = h_aao[r][c]
            hc = h_dis[r][c]

            ymax = max(hd.GetMaximum(), ha.GetMaximum(), hc.GetMaximum())
            if ymax <= 0.0:
                ymax = 1.0
            #endif
            ymax = 1.2 * ymax

            set_axes_and_range(hd, ymax)

            hd.Draw("hist")
            ha.Draw("hist same")
            hc.Draw("hist same")

            leg = ROOT.TLegend(0.55, 0.68, 0.94, 0.92)
            leg.SetBorderSize(1)
            leg.SetFillStyle(1001)
            leg.SetFillColor(ROOT.kWhite)
            leg.SetTextSize(0.045)
            leg.AddEntry(hd, f"data (N={int(c_data[r][c])})", "l")
            leg.AddEntry(ha, f"aaogen (N={int(c_aao[r][c])})", "l")
            leg.AddEntry(hc, f"clasdis (N={int(c_dis[r][c])})", "l")
            leg.Draw()

            tex = ROOT.TLatex()
            tex.SetNDC(True)
            tex.SetTextSize(0.05)
            tex.DrawLatex(0.14, 0.93,
                          f"xB [{xb_lo:.2f}, {xb_hi:.2f})   -tprime [{t_lo:.2f}, {t_hi:.2f})")
        #endfor
    #endfor

    canv.SaveAs(outpng)
#enddef


def draw_canvas_mix(h_data, h_mix, c_data, w_grid, outpng):
    nrows = len(h_data)
    ncols = len(h_data[0])

    canv = ROOT.TCanvas("c_mix", "Mx2: data vs per-bin mixture", 2400, 1400)
    canv.Divide(ncols, nrows, 0.001, 0.001)

    pad_idx = 1
    for r in range(nrows):
        xb_lo = XB_EDGES[r]
        xb_hi = XB_EDGES[r + 1]

        for c in range(ncols):
            t_lo = TNEG_EDGES[c]
            t_hi = TNEG_EDGES[c + 1]

            pad = canv.cd(pad_idx)
            pad_idx += 1
            pad_set_margins(pad)

            hd = h_data[r][c]
            hm = h_mix[r][c]
            w = w_grid[r][c]

            ymax = max(hd.GetMaximum(), hm.GetMaximum())
            if ymax <= 0.0:
                ymax = 1.0
            #endif
            ymax = 1.2 * ymax

            set_axes_and_range(hd, ymax)

            hd.Draw("hist")
            hm.Draw("hist same")

            leg = ROOT.TLegend(0.50, 0.70, 0.94, 0.92)
            leg.SetBorderSize(1)
            leg.SetFillStyle(1001)
            leg.SetFillColor(ROOT.kWhite)
            leg.SetTextSize(0.040)
            leg.AddEntry(hd, f"data (N={int(c_data[r][c])})", "l")
            leg.AddEntry(hm, f"mix (aaogen frac w={w:.4f})", "l")
            leg.AddEntry("", f"fit window: [{MX2_FIT_MIN:.1f}, {MX2_FIT_MAX:.1f}]", "")
            leg.Draw()

            tex = ROOT.TLatex()
            tex.SetNDC(True)
            tex.SetTextSize(0.05)
            tex.DrawLatex(0.14, 0.93,
                          f"xB [{xb_lo:.2f}, {xb_hi:.2f})   -tprime [{t_lo:.2f}, {t_hi:.2f})")
        #endfor
    #endfor

    canv.SaveAs(outpng)
#enddef


def write_weights_report(w_grid, wun_grid, sse_grid, c_data, c_aao, c_dis, path):
    ensure_outdir(path)
    nrows = len(w_grid)
    ncols = len(w_grid[0])

    total_sse = 0.0
    for r in range(nrows):
        for c in range(ncols):
            total_sse += sse_grid[r][c]
        #endfor
    #endfor

    with open(path, "w") as f:
        f.write("Per-bin mixture weights (shape-only, unweighted SSE)\n")
        f.write("Definition: H_mix = w * H_aaogen + (1-w) * H_clasdis\n")
        f.write(f"Fit window for w and SSE: Mx2 in [{MX2_FIT_MIN:.3f}, {MX2_FIT_MAX:.3f}]\n")
        f.write("Note: w is clipped to [0,1]. Report includes w_unclipped for diagnostics.\n\n")
        f.write(f"Total SSE (sum over pads) = {total_sse:.6e}\n\n")

        for r in range(nrows):
            xb_lo = XB_EDGES[r]
            xb_hi = XB_EDGES[r + 1]
            f.write(f"Row {r}: xB [{xb_lo:.2f}, {xb_hi:.2f})\n")
            for c in range(ncols):
                t_lo = TNEG_EDGES[c]
                t_hi = TNEG_EDGES[c + 1]
                w = w_grid[r][c]
                wun = wun_grid[r][c]
                sse = sse_grid[r][c]
                Nd = int(c_data[r][c])
                Na = int(c_aao[r][c])
                Nc = int(c_dis[r][c])
                f.write(
                    f"  Col {c}: -tprime [{t_lo:.2f}, {t_hi:.2f})  "
                    f"w={w:.6f}  w_unclipped={wun:.6f}  SSE={sse:.6e}  "
                    f"N(data,aao,dis)=({Nd},{Na},{Nc})\n"
                )
            #endfor
            f.write("\n")
        #endfor
#enddef


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True, help="Path to data ROOT file")
    ap.add_argument("--aaogen", required=True, help="Path to aaogen ROOT file")
    ap.add_argument("--clasdis", required=True, help="Path to clasdis ROOT file")
    ap.add_argument("--max_events", type=int, default=-1,
                    help="Max events per file (-1 means all events)")
    args = ap.parse_args()

    require_file(args.data)
    require_file(args.aaogen)
    require_file(args.clasdis)

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)

    f_data, t_data = open_tree(args.data, TREE_NAME)
    f_aao, t_aao = open_tree(args.aaogen, TREE_NAME)
    f_dis, t_dis = open_tree(args.clasdis, TREE_NAME)

    needed = ["runnum", "e_p", "e_theta", "e_phi", "p_p", "p_theta", "p_phi", "x", "Q2", "Mx2"]
    require_branches(t_data, needed, "data")
    require_branches(t_aao, needed, "aaogen")
    require_branches(t_dis, needed, "clasdis")

    ensure_outdir(OUTPUT_YIELDS_PNG)
    ensure_outdir(OUTPUT_MIX_PNG)
    ensure_outdir(OUTPUT_DATAONLY_PNG)
    ensure_outdir(OUTPUT_WEIGHTS_TXT)

    nrows = len(XB_EDGES) - 1
    ncols = len(TNEG_EDGES) - 1

    h_data = make_hist_grid("h_data", nrows, ncols)
    h_aao = make_hist_grid("h_aaogen", nrows, ncols)
    h_dis = make_hist_grid("h_clasdis", nrows, ncols)

    c_data = make_count_grid(nrows, ncols)
    c_aao = make_count_grid(nrows, ncols)
    c_dis = make_count_grid(nrows, ncols)

    max_events = None if args.max_events is None or int(args.max_events) < 0 else int(args.max_events)

    fill_all_bins_single_pass(t_data, h_data, c_data, max_events)
    fill_all_bins_single_pass(t_aao, h_aao, c_aao, max_events)
    fill_all_bins_single_pass(t_dis, h_dis, c_dis, max_events)

    # Build integrated histograms from RAW (unnormalized) per-pad histograms.
    h_data_int = make_integrated_hist_from_grid(h_data, "h_data_integrated")
    h_aao_int = make_integrated_hist_from_grid(h_aao, "h_aaogen_integrated")
    h_dis_int = make_integrated_hist_from_grid(h_dis, "h_clasdis_integrated")

    Nd_tot = sum_counts_grid(c_data)
    Na_tot = sum_counts_grid(c_aao)
    Nc_tot = sum_counts_grid(c_dis)

    # Normalize integrated shapes (like per-pad behavior).
    normalize_unit_area(h_data_int)
    normalize_unit_area(h_aao_int)
    normalize_unit_area(h_dis_int)

    col_data = ROOT.kBlack
    col_aao = ROOT.kRed
    col_dis = ROOT.kBlue

    style_hist(h_data_int, col_data, 2, 1)
    style_hist(h_aao_int, col_aao, 2, 2)
    style_hist(h_dis_int, col_dis, 2, 3)

    # Now normalize the per-pad histograms to unit area (shape-only).
    for r in range(nrows):
        for c in range(ncols):
            normalize_unit_area(h_data[r][c])
            normalize_unit_area(h_aao[r][c])
            normalize_unit_area(h_dis[r][c])
        #endfor
    #endfor

    for r in range(nrows):
        for c in range(ncols):
            style_hist(h_data[r][c], col_data, 2, 1)
            style_hist(h_aao[r][c], col_aao, 2, 2)
            style_hist(h_dis[r][c], col_dis, 2, 3)
        #endfor
    #endfor

    draw_canvas_threeway(h_data, h_aao, h_dis, c_data, c_aao, c_dis, OUTPUT_YIELDS_PNG)

    # Repurposed "data_only" output: 1-pad integrated three-way plot.
    draw_canvas_integrated_threeway(h_data_int, h_aao_int, h_dis_int, Nd_tot, Na_tot, Nc_tot, OUTPUT_DATAONLY_PNG)

    w_grid, wun_grid, sse_grid, h_mix = compute_w_grid_and_mix(h_data, h_aao, h_dis)

    col_mix = ROOT.kGreen + 2
    for r in range(nrows):
        for c in range(ncols):
            style_hist(h_mix[r][c], col_mix, 3, 1)
        #endfor
    #endfor

    total_sse = 0.0
    for r in range(nrows):
        for c in range(ncols):
            total_sse += sse_grid[r][c]
        #endfor
    #endfor
    print("Per-bin mixture fit (shape-only, unweighted SSE):")
    print(f"  fit window: Mx2 in [{MX2_FIT_MIN:.3f}, {MX2_FIT_MAX:.3f}]")
    print(f"  total SSE  = {total_sse:.6e}")
    print(f"  wrote weights report: {OUTPUT_WEIGHTS_TXT}")

    write_weights_report(w_grid, wun_grid, sse_grid, c_data, c_aao, c_dis, OUTPUT_WEIGHTS_TXT)

    draw_canvas_mix(h_data, h_mix, c_data, w_grid, OUTPUT_MIX_PNG)

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