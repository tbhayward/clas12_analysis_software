#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
normalization.py

Phase 1:
  - Fill 4x6 Mx2 histograms binned in (xB, -tprime) for data, aaogen, clasdis
  - Compute t, tmin, tprime on-the-fly (from runnum-based beam energy)
  - Normalize each histogram to unit area per (xB, -tprime) bin (shape-only)

Phase 2:
  - Determine a per-bin mixture fraction w[r,c] in [0,1]:
        H_mix(r,c) = w[r,c] * H_aaogen(r,c) + (1-w[r,c]) * H_clasdis(r,c)

  - Choose w[r,c] by minimizing a WEIGHTED chi2 between the data shape and mixture shape
    within an Mx2 fit window, where weights come from the DATA statistical uncertainty
    per Mx2 bin (derived from the RAW data counts in that pad):

        chi2_{r,c}(w) = sum_{i in window} ( D_i - (w*A_i + (1-w)*C_i) )^2 / sigma_i^2

    With:
        D_i = normalized data bin content
        A_i = normalized aaogen bin content
        C_i = normalized clasdis bin content
        n_i = raw (unnormalized) data bin count
        N   = total raw data counts in the pad
        sigma_i = sqrt(n_i)/N   (with a floor for n_i=0)

    Exact minimizer in each pad (over bins i in window):
        X_i = (A_i - C_i),  Y_i = (D_i - C_i),  w_i = 1/sigma_i^2
        w_unclipped = sum_i (w_i X_i Y_i) / sum_i (w_i X_i^2), clipped to [0,1].

Additional diagnostic output:
  - output/yields_data_only.png is a SINGLE-PAD integrated plot:
        (sum over all xB and -tprime bins) for data vs aaogen vs clasdis.

Optional Phase 3 (event-level mixing into a new ROOT file):
  - If --mc_aaogen and --mc_clasdis are provided, create an output ROOT file containing
    a per-(xB,-tprime) event mixture of those MC files using the weights w[r,c].
  - Strategy: write ALL clasdis events (max stats), then add aaogen events until each
    bin achieves the target mixture fraction:
        w = N_aao / (N_aao + N_dis)  =>  N_aao_target = w/(1-w) * N_dis
  - These MC files do NOT have runnum. Assume fixed beam energy Eb=10.55 GeV.
  - IMPORTANT: Any pre-existing t/tmin branches in the MC inputs are IGNORED.
    We recompute:
        t, tmin, tprime from reconstructed kinematics and (x,Q2)
        mc_t, mc_tmin, mc_tprime from generated kinematics and (mc_x,mc_Q2)
  - Bin assignment for mixing uses recomputed reco (x, tprime).

Outputs:
  output/yields.png           : data vs aaogen vs clasdis (shape-only) per pad
  output/yields_mix.png       : data vs per-bin mixture per pad, legend shows w[r,c]
  output/yields_data_only.png : integrated (all bins summed) data vs aaogen vs clasdis
  output/weights.txt          : per-bin w[r,c], w_unclipped, chi2 summary
  --out mixed MC root         : optional mixed MC ROOT file (Phase 3)

Debug outputs:
  output/mix_debug_report.txt : detailed per-bin and global counters for the mixed ROOT file
  output/mix_debug_mx2.png    : integrated Mx2 shapes of events WRITTEN to the output, split by source

Optional diagnostics:
  --force: override a few w[r,c] values (after the fit) for debugging:
      Row 0, Col 2 = 0.05
      Row 0, Col 3 = 0.08
      Row 0, Col 4 = 0.03
      Row 0, Col 5 = 0.03
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

# Histogram window.
MX2_MIN = 0.3
MX2_MAX = 2.0

# One third of the previous 100 Mx2 bins.
MX2_NBINS = 33

# Peak-focused fit window used for solving w and computing chi2.
MX2_FIT_MIN = 0.81
MX2_FIT_MAX = 1.10

OUTPUT_YIELDS_PNG = "output/yields.png"
OUTPUT_MIX_PNG = "output/yields_mix.png"
OUTPUT_DATAONLY_PNG = "output/yields_data_only.png"
OUTPUT_WEIGHTS_TXT = "output/weights.txt"

DEFAULT_MIXED_MC_ROOT = "output/mixed_mc.root"

# Debug outputs for Phase 3
OUTPUT_MIX_DEBUG_TXT = "output/mix_debug_report.txt"
OUTPUT_MIX_DEBUG_MX2_PNG = "output/mix_debug_mx2.png"

# Fixed beam energy for MC mixing stage (no runnum in MC files)
MC_EB_FIXED = 10.55

# Masses (GeV)
MASS_E = 0.000511
MASS_PI = 0.139570
MASS_N = 0.9382720813


def fatal(msg):
    raise RuntimeError(msg)
#enddef


def require_file(path):
    if path is None or str(path).strip() == "":
        fatal("Missing required input path.")
    #endif

    if not os.path.isfile(path):
        fatal("File not found: " + str(path))
    #endif
#enddef


def open_tree(path, tree_name):
    f = ROOT.TFile.Open(path, "READ")
    if not f or f.IsZombie():
        fatal("Failed to open ROOT file: " + str(path))
    #endif

    t = f.Get(tree_name)
    if not t:
        fatal("Tree '" + str(tree_name) + "' not found in: " + str(path))
    #endif

    return f, t
#enddef


def require_branches(t, needed, label):
    blist = t.GetListOfBranches()
    if not blist:
        fatal("No branch list available for tree in: " + str(label))
    #endif

    missing = []
    for b in needed:
        if not blist.FindObject(b):
            missing.append(b)
        #endif
    #endfor

    if len(missing) > 0:
        fatal("Missing required branches in " + str(label) + ": " + ", ".join(missing))
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


def compute_t_scalar_from_Eb(Eb, e_p, e_theta, e_phi, p_p, p_theta, p_phi):
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


def compute_t_scalar(runnum, e_p, e_theta, e_phi, p_p, p_theta, p_phi):
    Eb = beam_energy(int(runnum))
    return compute_t_scalar_from_Eb(Eb, e_p, e_theta, e_phi, p_p, p_theta, p_phi)
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

        t_val = compute_t_scalar(
            int(runnum[0]),
            float(e_p[0]),
            float(e_theta[0]),
            float(e_phi[0]),
            float(p_p[0]),
            float(p_theta[0]),
            float(p_phi[0])
        )

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
    d = os.path.dirname(str(path))
    if d != "":
        os.makedirs(d, exist_ok=True)
    #endif
#enddef


def sigma_from_raw_counts(n_i, N_pad):
    if N_pad <= 0:
        return 1.0
    #endif

    if n_i <= 0.0:
        return 1.0 / float(N_pad)
    #endif

    return math.sqrt(n_i) / float(N_pad)
#enddef


def compute_best_w_weighted_for_pad(hd_norm, ha_norm, hc_norm, hd_raw):
    nb = hd_norm.GetNbinsX()

    N_pad = hd_raw.Integral(1, nb)
    if N_pad <= 0.0 or hd_norm.Integral(1, nb) <= 0.0:
        return 0.0, 0.0, 0.0
    #endif

    num = 0.0
    den = 0.0
    used_bins = 0

    for i in range(1, nb + 1):
        xcen = hd_norm.GetXaxis().GetBinCenter(i)

        if xcen < MX2_FIT_MIN or xcen > MX2_FIT_MAX:
            continue
        #endif

        D = hd_norm.GetBinContent(i)
        A = ha_norm.GetBinContent(i)
        C = hc_norm.GetBinContent(i)

        n_i = hd_raw.GetBinContent(i)
        sig = sigma_from_raw_counts(n_i, N_pad)

        if sig <= 0.0:
            continue
        #endif

        wgt = 1.0 / (sig * sig)

        X = (A - C)
        Y = (D - C)

        num += wgt * X * Y
        den += wgt * X * X
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

    chi2 = 0.0

    for i in range(1, nb + 1):
        xcen = hd_norm.GetXaxis().GetBinCenter(i)

        if xcen < MX2_FIT_MIN or xcen > MX2_FIT_MAX:
            continue
        #endif

        D = hd_norm.GetBinContent(i)
        M = w * ha_norm.GetBinContent(i) + (1.0 - w) * hc_norm.GetBinContent(i)

        n_i = hd_raw.GetBinContent(i)
        sig = sigma_from_raw_counts(n_i, N_pad)

        if sig <= 0.0:
            continue
        #endif

        diff = D - M
        chi2 += (diff * diff) / (sig * sig)
    #endfor

    return w, w_unclipped, chi2
#enddef


def compute_chi2_for_pad_given_w(hd_norm, ha_norm, hc_norm, hd_raw, w):
    nb = hd_norm.GetNbinsX()

    N_pad = hd_raw.Integral(1, nb)
    if N_pad <= 0.0 or hd_norm.Integral(1, nb) <= 0.0:
        return 0.0
    #endif

    chi2 = 0.0

    for i in range(1, nb + 1):
        xcen = hd_norm.GetXaxis().GetBinCenter(i)

        if xcen < MX2_FIT_MIN or xcen > MX2_FIT_MAX:
            continue
        #endif

        D = hd_norm.GetBinContent(i)
        M = w * ha_norm.GetBinContent(i) + (1.0 - w) * hc_norm.GetBinContent(i)

        n_i = hd_raw.GetBinContent(i)
        sig = sigma_from_raw_counts(n_i, N_pad)

        if sig <= 0.0:
            continue
        #endif

        diff = D - M
        chi2 += (diff * diff) / (sig * sig)
    #endfor

    return chi2
#enddef


def compute_w_grid_and_mix(h_data, h_aao, h_dis, h_data_raw):
    nrows = len(h_data)
    ncols = len(h_data[0])

    w_grid = [[0.0 for _ in range(ncols)] for _ in range(nrows)]
    wun_grid = [[0.0 for _ in range(ncols)] for _ in range(nrows)]
    chi2_grid = [[0.0 for _ in range(ncols)] for _ in range(nrows)]
    h_mix = []

    for r in range(nrows):
        row = []

        for c in range(ncols):
            hd = h_data[r][c]
            ha = h_aao[r][c]
            hc = h_dis[r][c]
            hdr = h_data_raw[r][c]

            w, w_unclipped, chi2 = compute_best_w_weighted_for_pad(hd, ha, hc, hdr)

            w_grid[r][c] = w
            wun_grid[r][c] = w_unclipped
            chi2_grid[r][c] = chi2

            hm = ha.Clone(f"h_mix_r{r}_c{c}")
            hm.Reset("ICESM")
            hm.Add(ha, w)
            hm.Add(hc, (1.0 - w))
            row.append(hm)
        #endfor

        h_mix.append(row)
    #endfor

    return w_grid, wun_grid, chi2_grid, h_mix
#enddef


def pad_set_margins(pad):
    pad.SetGrid(1, 1)
    pad.SetLeftMargin(0.20)
    pad.SetRightMargin(0.04)
    pad.SetBottomMargin(0.14)
    pad.SetTopMargin(0.08)
#enddef


def set_axes_and_range(h_frame, ymax):
    h_frame.GetXaxis().SetTitle("M_{x}^{2} (GeV^{2})")
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
    tex.DrawLatex(0.16, 0.93, "Integrated over all x_{B} and #minus t^{#prime} bins")

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
            tex.DrawLatex(
                0.14,
                0.93,
                f"x_{{B}} [{xb_lo:.2f}, {xb_hi:.2f})   #minus t^{{#prime}} [{t_lo:.2f}, {t_hi:.2f}) (GeV^{{2}})"
            )
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
            leg.AddEntry("", f"peak window: [{MX2_FIT_MIN:.2f}, {MX2_FIT_MAX:.2f}]", "")
            leg.Draw()

            tex = ROOT.TLatex()
            tex.SetNDC(True)
            tex.SetTextSize(0.05)
            tex.DrawLatex(
                0.14,
                0.93,
                f"x_{{B}} [{xb_lo:.2f}, {xb_hi:.2f})   #minus t^{{#prime}} [{t_lo:.2f}, {t_hi:.2f}) (GeV^{{2}})"
            )
        #endfor
    #endfor

    canv.SaveAs(outpng)
#enddef


def write_weights_report(w_grid, wun_grid, chi2_grid, c_data, c_aao, c_dis, path, forced_map=None):
    ensure_outdir(path)

    nrows = len(w_grid)
    ncols = len(w_grid[0])

    total_chi2 = 0.0

    for r in range(nrows):
        for c in range(ncols):
            total_chi2 += chi2_grid[r][c]
        #endfor
    #endfor

    with open(path, "w") as f:
        f.write("Per-bin mixture weights (shape-only, WEIGHTED chi2)\n")
        f.write("Definition: H_mix = w * H_aaogen + (1-w) * H_clasdis\n")
        f.write(f"Peak window for w and chi2: Mx2 in [{MX2_FIT_MIN:.3f}, {MX2_FIT_MAX:.3f}]\n")
        f.write(f"Histogram binning: {MX2_NBINS} bins from {MX2_MIN:.3f} to {MX2_MAX:.3f} GeV^2\n")
        f.write("Weights: sigma_i from RAW data counts in the pad: sigma_i = sqrt(n_i)/N (floor for n_i=0)\n")
        f.write("Note: w is clipped to [0,1]. Report includes w_unclipped for diagnostics.\n")

        if forced_map is not None and len(forced_map) > 0:
            f.write("NOTE: Some w values were FORCED by --force and chi2 was recomputed for those pads.\n")
        #endif

        f.write("\n")
        f.write(f"Total chi2 (sum over pads) = {total_chi2:.6e}\n\n")

        for r in range(nrows):
            xb_lo = XB_EDGES[r]
            xb_hi = XB_EDGES[r + 1]
            f.write(f"Row {r}: xB [{xb_lo:.2f}, {xb_hi:.2f})\n")

            for c in range(ncols):
                t_lo = TNEG_EDGES[c]
                t_hi = TNEG_EDGES[c + 1]
                w = w_grid[r][c]
                wun = wun_grid[r][c]
                chi2 = chi2_grid[r][c]
                Nd = int(c_data[r][c])
                Na = int(c_aao[r][c])
                Nc = int(c_dis[r][c])

                tag = ""

                if forced_map is not None:
                    if (r, c) in forced_map:
                        tag = "  FORCED"
                    #endif
                #endif

                f.write(
                    f"  Col {c}: -tprime [{t_lo:.2f}, {t_hi:.2f})  "
                    f"w={w:.6f}  w_unclipped={wun:.6f}  chi2={chi2:.6e}  "
                    f"N(data,aao,dis)=({Nd},{Na},{Nc}){tag}\n"
                )
            #endfor

            f.write("\n")
        #endfor
    #endwith
#enddef


def compute_bin_indices_from_reco(x_val, Q2_val, e_p, e_theta, e_phi, p_p, p_theta, p_phi):
    """
    For MC mixing stage:
      - Use fixed Eb=MC_EB_FIXED
      - Recompute t, tmin, tprime from reco kinematics and (x,Q2)
      - Return (row, col, t, tmin, tprime)
    """
    t_val = compute_t_scalar_from_Eb(MC_EB_FIXED, e_p, e_theta, e_phi, p_p, p_theta, p_phi)
    tmin_val = compute_tmin_exact(x_val, Q2_val)
    tprime_val = t_val - tmin_val
    tneg = -tprime_val

    r = find_bin(x_val, XB_EDGES)
    if r < 0:
        return -1, -1, t_val, tmin_val, tprime_val
    #endif

    c = find_bin(tneg, TNEG_EDGES)
    if c < 0:
        return -1, -1, t_val, tmin_val, tprime_val
    #endif

    return r, c, t_val, tmin_val, tprime_val
#enddef


def compute_mc_t_quantities_from_gen(mc_x_val, mc_Q2_val, mc_e_p, mc_e_theta, mc_e_phi, mc_p_p, mc_p_theta, mc_p_phi):
    """
    Recompute generated-level t, tmin, tprime using fixed Eb=MC_EB_FIXED.
    """
    mc_t_val = compute_t_scalar_from_Eb(
        MC_EB_FIXED,
        mc_e_p,
        mc_e_theta,
        mc_e_phi,
        mc_p_p,
        mc_p_theta,
        mc_p_phi
    )

    mc_tmin_val = compute_tmin_exact(mc_x_val, mc_Q2_val)
    mc_tprime_val = mc_t_val - mc_tmin_val

    return mc_t_val, mc_tmin_val, mc_tprime_val
#enddef


def mix_mc_to_output_root(mc_aao_path, mc_dis_path, out_root_path, w_grid, max_events):
    require_file(mc_aao_path)
    require_file(mc_dis_path)

    f_aao, t_aao = open_tree(mc_aao_path, TREE_NAME)
    f_dis, t_dis = open_tree(mc_dis_path, TREE_NAME)

    mc_needed = [
        "e_p", "e_theta", "e_phi", "vz_e",
        "p_p", "p_theta", "p_phi", "vz_p",
        "Q2", "W", "Mx", "Mx2", "x", "y", "z", "xF", "pT", "xi", "eta", "phi",
        "DepA", "DepB", "DepC", "DepV", "DepW",
        "mc_e_p", "mc_e_theta", "mc_e_phi", "mc_vz_e",
        "mc_p_p", "mc_p_theta", "mc_p_phi", "mc_vz_p",
        "mc_Q2", "mc_W", "mc_Mx", "mc_Mx2", "mc_x", "mc_y", "mc_z",
        "mc_xF", "mc_pT", "mc_xi", "mc_eta", "mc_phi",
        "mc_DepA", "mc_DepB", "mc_DepC", "mc_DepV", "mc_DepW",
        "matching_e_pid", "matching_p1_pid", "mc_p1_parent"
    ]

    require_branches(t_aao, mc_needed, "mc_aaogen")
    require_branches(t_dis, mc_needed, "mc_clasdis")

    ensure_outdir(out_root_path)

    fout = ROOT.TFile.Open(out_root_path, "RECREATE")
    if not fout or fout.IsZombie():
        fatal("Failed to create output ROOT file: " + str(out_root_path))
    #endif

    tout = ROOT.TTree(TREE_NAME, "Mixed MC: write ALL clasdis; top up aaogen in-grid using w_grid")

    e_p = array("d", [0.0])
    e_theta = array("d", [0.0])
    e_phi = array("d", [0.0])
    vz_e = array("d", [0.0])
    p_p = array("d", [0.0])
    p_theta = array("d", [0.0])
    p_phi = array("d", [0.0])
    vz_p = array("d", [0.0])
    Q2 = array("d", [0.0])
    W = array("d", [0.0])
    Mx = array("d", [0.0])
    Mx2 = array("d", [0.0])
    x = array("d", [0.0])
    y = array("d", [0.0])
    z = array("d", [0.0])
    xF = array("d", [0.0])
    pT = array("d", [0.0])
    xi = array("d", [0.0])
    eta = array("d", [0.0])
    phi = array("d", [0.0])
    DepA = array("d", [0.0])
    DepB = array("d", [0.0])
    DepC = array("d", [0.0])
    DepV = array("d", [0.0])
    DepW = array("d", [0.0])

    mc_e_p = array("d", [0.0])
    mc_e_theta = array("d", [0.0])
    mc_e_phi = array("d", [0.0])
    mc_vz_e = array("d", [0.0])
    mc_p_p = array("d", [0.0])
    mc_p_theta = array("d", [0.0])
    mc_p_phi = array("d", [0.0])
    mc_vz_p = array("d", [0.0])
    mc_Q2 = array("d", [0.0])
    mc_W = array("d", [0.0])
    mc_Mx = array("d", [0.0])
    mc_Mx2 = array("d", [0.0])
    mc_x = array("d", [0.0])
    mc_y = array("d", [0.0])
    mc_z = array("d", [0.0])
    mc_xF = array("d", [0.0])
    mc_pT = array("d", [0.0])
    mc_xi = array("d", [0.0])
    mc_eta = array("d", [0.0])
    mc_phi = array("d", [0.0])
    mc_DepA = array("d", [0.0])
    mc_DepB = array("d", [0.0])
    mc_DepC = array("d", [0.0])
    mc_DepV = array("d", [0.0])
    mc_DepW = array("d", [0.0])

    matching_e_pid = array("i", [0])
    matching_p1_pid = array("i", [0])
    mc_p1_parent = array("i", [0])

    t = array("d", [0.0])
    tmin = array("d", [0.0])
    tprime = array("d", [0.0])
    mc_t = array("d", [0.0])
    mc_tmin = array("d", [0.0])
    mc_tprime = array("d", [0.0])

    tout.Branch("e_p", e_p, "e_p/D")
    tout.Branch("e_theta", e_theta, "e_theta/D")
    tout.Branch("e_phi", e_phi, "e_phi/D")
    tout.Branch("vz_e", vz_e, "vz_e/D")
    tout.Branch("p_p", p_p, "p_p/D")
    tout.Branch("p_theta", p_theta, "p_theta/D")
    tout.Branch("p_phi", p_phi, "p_phi/D")
    tout.Branch("vz_p", vz_p, "vz_p/D")
    tout.Branch("Q2", Q2, "Q2/D")
    tout.Branch("W", W, "W/D")
    tout.Branch("Mx", Mx, "Mx/D")
    tout.Branch("Mx2", Mx2, "Mx2/D")
    tout.Branch("x", x, "x/D")
    tout.Branch("y", y, "y/D")
    tout.Branch("z", z, "z/D")
    tout.Branch("xF", xF, "xF/D")
    tout.Branch("pT", pT, "pT/D")
    tout.Branch("xi", xi, "xi/D")
    tout.Branch("eta", eta, "eta/D")
    tout.Branch("phi", phi, "phi/D")
    tout.Branch("DepA", DepA, "DepA/D")
    tout.Branch("DepB", DepB, "DepB/D")
    tout.Branch("DepC", DepC, "DepC/D")
    tout.Branch("DepV", DepV, "DepV/D")
    tout.Branch("DepW", DepW, "DepW/D")

    tout.Branch("mc_e_p", mc_e_p, "mc_e_p/D")
    tout.Branch("mc_e_theta", mc_e_theta, "mc_e_theta/D")
    tout.Branch("mc_e_phi", mc_e_phi, "mc_e_phi/D")
    tout.Branch("mc_vz_e", mc_vz_e, "mc_vz_e/D")
    tout.Branch("mc_p_p", mc_p_p, "mc_p_p/D")
    tout.Branch("mc_p_theta", mc_p_theta, "mc_p_theta/D")
    tout.Branch("mc_p_phi", mc_p_phi, "mc_p_phi/D")
    tout.Branch("mc_vz_p", mc_vz_p, "mc_vz_p/D")
    tout.Branch("mc_Q2", mc_Q2, "mc_Q2/D")
    tout.Branch("mc_W", mc_W, "mc_W/D")
    tout.Branch("mc_Mx", mc_Mx, "mc_Mx/D")
    tout.Branch("mc_Mx2", mc_Mx2, "mc_Mx2/D")
    tout.Branch("mc_x", mc_x, "mc_x/D")
    tout.Branch("mc_y", mc_y, "mc_y/D")
    tout.Branch("mc_z", mc_z, "mc_z/D")
    tout.Branch("mc_xF", mc_xF, "mc_xF/D")
    tout.Branch("mc_pT", mc_pT, "mc_pT/D")
    tout.Branch("mc_xi", mc_xi, "mc_xi/D")
    tout.Branch("mc_eta", mc_eta, "mc_eta/D")
    tout.Branch("mc_phi", mc_phi, "mc_phi/D")
    tout.Branch("mc_DepA", mc_DepA, "mc_DepA/D")
    tout.Branch("mc_DepB", mc_DepB, "mc_DepB/D")
    tout.Branch("mc_DepC", mc_DepC, "mc_DepC/D")
    tout.Branch("mc_DepV", mc_DepV, "mc_DepV/D")
    tout.Branch("mc_DepW", mc_DepW, "mc_DepW/D")

    tout.Branch("matching_e_pid", matching_e_pid, "matching_e_pid/I")
    tout.Branch("matching_p1_pid", matching_p1_pid, "matching_p1_pid/I")
    tout.Branch("mc_p1_parent", mc_p1_parent, "mc_p1_parent/I")

    tout.Branch("t", t, "t/D")
    tout.Branch("tmin", tmin, "tmin/D")
    tout.Branch("tprime", tprime, "tprime/D")
    tout.Branch("mc_t", mc_t, "mc_t/D")
    tout.Branch("mc_tmin", mc_tmin, "mc_tmin/D")
    tout.Branch("mc_tprime", mc_tprime, "mc_tprime/D")

    nrows = len(XB_EDGES) - 1
    ncols = len(TNEG_EDGES) - 1

    Ndis = [[0 for _ in range(ncols)] for _ in range(nrows)]

    clasdis_written_total = 0
    clasdis_written_in_grid = 0
    clasdis_written_out_grid = 0

    aaogen_written_total = 0
    aaogen_written_in_grid = 0

    aaogen_skipped_out_grid = 0
    aaogen_skipped_quota_met = 0

    def bind_mc_tree(tree):
        tree.SetBranchStatus("*", 0)

        for bn in mc_needed:
            tree.SetBranchStatus(bn, 1)
        #endfor

        tree.SetBranchAddress("e_p", e_p)
        tree.SetBranchAddress("e_theta", e_theta)
        tree.SetBranchAddress("e_phi", e_phi)
        tree.SetBranchAddress("vz_e", vz_e)
        tree.SetBranchAddress("p_p", p_p)
        tree.SetBranchAddress("p_theta", p_theta)
        tree.SetBranchAddress("p_phi", p_phi)
        tree.SetBranchAddress("vz_p", vz_p)
        tree.SetBranchAddress("Q2", Q2)
        tree.SetBranchAddress("W", W)
        tree.SetBranchAddress("Mx", Mx)
        tree.SetBranchAddress("Mx2", Mx2)
        tree.SetBranchAddress("x", x)
        tree.SetBranchAddress("y", y)
        tree.SetBranchAddress("z", z)
        tree.SetBranchAddress("xF", xF)
        tree.SetBranchAddress("pT", pT)
        tree.SetBranchAddress("xi", xi)
        tree.SetBranchAddress("eta", eta)
        tree.SetBranchAddress("phi", phi)
        tree.SetBranchAddress("DepA", DepA)
        tree.SetBranchAddress("DepB", DepB)
        tree.SetBranchAddress("DepC", DepC)
        tree.SetBranchAddress("DepV", DepV)
        tree.SetBranchAddress("DepW", DepW)

        tree.SetBranchAddress("mc_e_p", mc_e_p)
        tree.SetBranchAddress("mc_e_theta", mc_e_theta)
        tree.SetBranchAddress("mc_e_phi", mc_e_phi)
        tree.SetBranchAddress("mc_vz_e", mc_vz_e)
        tree.SetBranchAddress("mc_p_p", mc_p_p)
        tree.SetBranchAddress("mc_p_theta", mc_p_theta)
        tree.SetBranchAddress("mc_p_phi", mc_p_phi)
        tree.SetBranchAddress("mc_vz_p", mc_vz_p)
        tree.SetBranchAddress("mc_Q2", mc_Q2)
        tree.SetBranchAddress("mc_W", mc_W)
        tree.SetBranchAddress("mc_Mx", mc_Mx)
        tree.SetBranchAddress("mc_Mx2", mc_Mx2)
        tree.SetBranchAddress("mc_x", mc_x)
        tree.SetBranchAddress("mc_y", mc_y)
        tree.SetBranchAddress("mc_z", mc_z)
        tree.SetBranchAddress("mc_xF", mc_xF)
        tree.SetBranchAddress("mc_pT", mc_pT)
        tree.SetBranchAddress("mc_xi", mc_xi)
        tree.SetBranchAddress("mc_eta", mc_eta)
        tree.SetBranchAddress("mc_phi", mc_phi)
        tree.SetBranchAddress("mc_DepA", mc_DepA)
        tree.SetBranchAddress("mc_DepB", mc_DepB)
        tree.SetBranchAddress("mc_DepC", mc_DepC)
        tree.SetBranchAddress("mc_DepV", mc_DepV)
        tree.SetBranchAddress("mc_DepW", mc_DepW)

        tree.SetBranchAddress("matching_e_pid", matching_e_pid)
        tree.SetBranchAddress("matching_p1_pid", matching_p1_pid)
        tree.SetBranchAddress("mc_p1_parent", mc_p1_parent)
    #enddef

    bind_mc_tree(t_dis)

    n_entries_dis = int(t_dis.GetEntries())

    if max_events is None:
        n_to_process_dis = n_entries_dis
    else:
        n_to_process_dis = min(n_entries_dis, int(max_events))
    #endif

    for i in range(n_to_process_dis):
        t_dis.GetEntry(i)

        r, c, t_val, tmin_val, tp_val = compute_bin_indices_from_reco(
            float(x[0]),
            float(Q2[0]),
            float(e_p[0]),
            float(e_theta[0]),
            float(e_phi[0]),
            float(p_p[0]),
            float(p_theta[0]),
            float(p_phi[0])
        )

        mc_t_val, mc_tmin_val, mc_tp_val = compute_mc_t_quantities_from_gen(
            float(mc_x[0]),
            float(mc_Q2[0]),
            float(mc_e_p[0]),
            float(mc_e_theta[0]),
            float(mc_e_phi[0]),
            float(mc_p_p[0]),
            float(mc_p_theta[0]),
            float(mc_p_phi[0])
        )

        t[0] = float(t_val)
        tmin[0] = float(tmin_val)
        tprime[0] = float(tp_val)
        mc_t[0] = float(mc_t_val)
        mc_tmin[0] = float(mc_tmin_val)
        mc_tprime[0] = float(mc_tp_val)

        tout.Fill()
        clasdis_written_total += 1

        if r >= 0 and c >= 0:
            Ndis[r][c] += 1
            clasdis_written_in_grid += 1
        else:
            clasdis_written_out_grid += 1
        #endif
    #endfor

    Naao_target = [[0 for _ in range(ncols)] for _ in range(nrows)]

    for r in range(nrows):
        for c in range(ncols):
            w = float(w_grid[r][c])
            nd = int(Ndis[r][c])

            if nd <= 0:
                Naao_target[r][c] = 0
                continue
            #endif

            if w <= 0.0:
                Naao_target[r][c] = 0
                continue
            #endif

            if w >= 1.0:
                Naao_target[r][c] = 10**18
                continue
            #endif

            need = (w / (1.0 - w)) * float(nd)
            Naao_target[r][c] = int(math.ceil(need))
        #endfor
    #endfor

    Naao_written = [[0 for _ in range(ncols)] for _ in range(nrows)]

    bind_mc_tree(t_aao)

    n_entries_aao = int(t_aao.GetEntries())

    if max_events is None:
        n_to_process_aao = n_entries_aao
    else:
        n_to_process_aao = min(n_entries_aao, int(max_events))
    #endif

    for i in range(n_to_process_aao):
        t_aao.GetEntry(i)

        r, c, t_val, tmin_val, tp_val = compute_bin_indices_from_reco(
            float(x[0]),
            float(Q2[0]),
            float(e_p[0]),
            float(e_theta[0]),
            float(e_phi[0]),
            float(p_p[0]),
            float(p_theta[0]),
            float(p_phi[0])
        )

        if r < 0 or c < 0:
            aaogen_skipped_out_grid += 1
            continue
        #endif

        if Naao_written[r][c] >= Naao_target[r][c]:
            aaogen_skipped_quota_met += 1
            continue
        #endif

        mc_t_val, mc_tmin_val, mc_tp_val = compute_mc_t_quantities_from_gen(
            float(mc_x[0]),
            float(mc_Q2[0]),
            float(mc_e_p[0]),
            float(mc_e_theta[0]),
            float(mc_e_phi[0]),
            float(mc_p_p[0]),
            float(mc_p_theta[0]),
            float(mc_p_phi[0])
        )

        t[0] = float(t_val)
        tmin[0] = float(tmin_val)
        tprime[0] = float(tp_val)
        mc_t[0] = float(mc_t_val)
        mc_tmin[0] = float(mc_tmin_val)
        mc_tprime[0] = float(mc_tp_val)

        tout.Fill()

        Naao_written[r][c] += 1

        aaogen_written_total += 1
        aaogen_written_in_grid += 1
    #endfor

    fout.cd()
    tout.Write()
    fout.Close()
    f_aao.Close()
    f_dis.Close()

    print("Phase 3 debug:")
    print("  Wrote ALL clasdis events (in-grid + out-of-grid). No Mx2 cuts were applied.")
    print(f"  clasdis written total       = {clasdis_written_total}")
    print(f"  clasdis written in-grid     = {clasdis_written_in_grid}")
    print(f"  clasdis written out-of-grid = {clasdis_written_out_grid}")
    print(f"  aaogen written total        = {aaogen_written_total}")
    print(f"  aaogen written in-grid      = {aaogen_written_in_grid}")
    print(f"  aaogen skipped out-of-grid  = {aaogen_skipped_out_grid}")
    print(f"  aaogen skipped quota met    = {aaogen_skipped_quota_met}")

    denom = float(aaogen_written_total + clasdis_written_total)

    if denom > 0.0:
        print(f"  achieved global aaogen fraction (written) = {aaogen_written_total / denom:.6f}")
    else:
        print("  achieved global aaogen fraction (written) = 0.0")
    #endif
#enddef


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("--data", required=True, help="Path to data ROOT file")
    ap.add_argument("--aaogen", required=True, help="Path to aaogen ROOT file")
    ap.add_argument("--clasdis", required=True, help="Path to clasdis ROOT file")
    ap.add_argument("--max_events", type=int, default=-1, help="Max events per file (-1 means all events)")
    ap.add_argument("--mc_aaogen", default=None, help="Optional: Path to MC aaogen ROOT file for event-level mixing")
    ap.add_argument("--mc_clasdis", default=None, help="Optional: Path to MC clasdis ROOT file for event-level mixing")
    ap.add_argument("--out", default=DEFAULT_MIXED_MC_ROOT, help="Optional: Output ROOT file for mixed MC")
    ap.add_argument("--force", action="store_true", help="Optional diagnostic: force specific w[r,c] values after fit")

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
    ensure_outdir(OUTPUT_MIX_DEBUG_TXT)
    ensure_outdir(OUTPUT_MIX_DEBUG_MX2_PNG)

    nrows = len(XB_EDGES) - 1
    ncols = len(TNEG_EDGES) - 1

    h_data = make_hist_grid("h_data", nrows, ncols)
    h_aao = make_hist_grid("h_aaogen", nrows, ncols)
    h_dis = make_hist_grid("h_clasdis", nrows, ncols)

    c_data = make_count_grid(nrows, ncols)
    c_aao = make_count_grid(nrows, ncols)
    c_dis = make_count_grid(nrows, ncols)

    if args.max_events is None or int(args.max_events) < 0:
        max_events = None
    else:
        max_events = int(args.max_events)
    #endif

    fill_all_bins_single_pass(t_data, h_data, c_data, max_events)
    fill_all_bins_single_pass(t_aao, h_aao, c_aao, max_events)
    fill_all_bins_single_pass(t_dis, h_dis, c_dis, max_events)

    h_data_raw = []

    for r in range(nrows):
        row = []

        for c in range(ncols):
            row.append(h_data[r][c].Clone(f"h_data_raw_r{r}_c{c}"))
        #endfor

        h_data_raw.append(row)
    #endfor

    h_data_int = make_integrated_hist_from_grid(h_data, "h_data_integrated")
    h_aao_int = make_integrated_hist_from_grid(h_aao, "h_aaogen_integrated")
    h_dis_int = make_integrated_hist_from_grid(h_dis, "h_clasdis_integrated")

    Nd_tot = sum_counts_grid(c_data)
    Na_tot = sum_counts_grid(c_aao)
    Nc_tot = sum_counts_grid(c_dis)

    normalize_unit_area(h_data_int)
    normalize_unit_area(h_aao_int)
    normalize_unit_area(h_dis_int)

    col_data = ROOT.kBlack
    col_aao = ROOT.kRed
    col_dis = ROOT.kBlue

    style_hist(h_data_int, col_data, 2, 1)
    style_hist(h_aao_int, col_aao, 2, 2)
    style_hist(h_dis_int, col_dis, 2, 3)

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
    draw_canvas_integrated_threeway(h_data_int, h_aao_int, h_dis_int, Nd_tot, Na_tot, Nc_tot, OUTPUT_DATAONLY_PNG)

    w_grid, wun_grid, chi2_grid, h_mix = compute_w_grid_and_mix(h_data, h_aao, h_dis, h_data_raw)

    forced_map = {}

    if args.force:
        forced_map[(0, 2)] = 0.05
        forced_map[(0, 3)] = 0.08
        forced_map[(0, 4)] = 0.03
        forced_map[(0, 5)] = 0.03

        print("FORCE mode enabled: overriding selected w[r,c] after fit:")

        for (rr, cc), wv in forced_map.items():
            old = float(w_grid[rr][cc])
            w_grid[rr][cc] = float(wv)

            hm = h_mix[rr][cc]
            hm.Reset("ICESM")
            hm.Add(h_aao[rr][cc], float(w_grid[rr][cc]))
            hm.Add(h_dis[rr][cc], float(1.0 - w_grid[rr][cc]))

            chi2_grid[rr][cc] = compute_chi2_for_pad_given_w(
                h_data[rr][cc],
                h_aao[rr][cc],
                h_dis[rr][cc],
                h_data_raw[rr][cc],
                float(w_grid[rr][cc])
            )

            print(f"  Row {rr}, Col {cc}: w {old:.6f} -> {w_grid[rr][cc]:.6f}")
        #endfor
    #endif

    col_mix = ROOT.kGreen + 2

    for r in range(nrows):
        for c in range(ncols):
            style_hist(h_mix[r][c], col_mix, 3, 1)
        #endfor
    #endfor

    total_chi2 = 0.0

    for r in range(nrows):
        for c in range(ncols):
            total_chi2 += chi2_grid[r][c]
        #endfor
    #endfor

    print("Per-bin mixture fit (shape-only, WEIGHTED chi2):")
    print(f"  histogram binning: {MX2_NBINS} bins from {MX2_MIN:.3f} to {MX2_MAX:.3f} GeV^2")
    print(f"  peak window: Mx2 in [{MX2_FIT_MIN:.3f}, {MX2_FIT_MAX:.3f}] GeV^2")
    print(f"  total chi2 = {total_chi2:.6e}")
    print(f"  wrote weights report: {OUTPUT_WEIGHTS_TXT}")

    write_weights_report(w_grid, wun_grid, chi2_grid, c_data, c_aao, c_dis, OUTPUT_WEIGHTS_TXT, forced_map)
    draw_canvas_mix(h_data, h_mix, c_data, w_grid, OUTPUT_MIX_PNG)

    if args.mc_aaogen is not None or args.mc_clasdis is not None:
        if args.mc_aaogen is None or args.mc_clasdis is None:
            fatal("If using Phase 3, you must provide BOTH --mc_aaogen and --mc_clasdis.")
        #endif

        require_file(args.mc_aaogen)
        require_file(args.mc_clasdis)

        out_path = str(args.out)

        if out_path.strip() == "":
            fatal("--out may not be empty when Phase 3 is enabled.")
        #endif

        print("Phase 3: writing mixed MC ROOT file")
        print("  mc_aaogen  =", args.mc_aaogen)
        print("  mc_clasdis =", args.mc_clasdis)
        print("  out        =", out_path)
        print(f"  recompute t,tmin,tprime and mc_t,mc_tmin,mc_tprime with Eb={MC_EB_FIXED:.3f} (ignore any existing t/tmin branches)")
        print(f"  will write debug report: {OUTPUT_MIX_DEBUG_TXT}")
        print(f"  will write debug Mx2 png: {OUTPUT_MIX_DEBUG_MX2_PNG}")

        mix_mc_to_output_root(args.mc_aaogen, args.mc_clasdis, out_path, w_grid, max_events)

        print("  wrote mixed MC: ", out_path)
    #endif

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