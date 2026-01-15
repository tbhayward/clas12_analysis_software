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

Phase 3 (optional):
  - If --mc_aaogen and --mc_clasdis are provided:
      * Build a mixed MC ROOT file by combining events from the provided MC aaogen/clasdis trees
        on a bin-by-bin basis in (xB, -tprime), using the computed w[r,c] from Phase 2.
      * In each (r,c) bin, use ALL clasdis events, then add enough aaogen events to approach
        the target aaogen fraction w[r,c]. If aaogen runs out in a bin, use all available.
      * The new MC files do NOT contain runnum, so beam energy is assumed constant Eb=10.55 GeV.
      * Copy the required branches into the new output tree, and additionally write:
            t, tmin, tprime, mc_t, mc_tmin, mc_tprime

Outputs:
  output/yields.png           : data vs aaogen vs clasdis (shape-only) per pad
  output/yields_mix.png       : data vs per-bin mixture per pad, legend shows w[r,c]
  output/yields_data_only.png : integrated (all bins summed) data vs aaogen vs clasdis
  output/weights.txt          : per-bin w[r,c], w_unclipped, chi2 summary
  output/mixed_mc.root        : (optional) mixed MC built from --mc_aaogen/--mc_clasdis
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

# Histogram window (also fit window)
MX2_MIN = 0.3
MX2_MAX = 2.0
MX2_NBINS = 100

# Fit window used for solving w and computing chi2 (kept explicit even if identical)
MX2_FIT_MIN = 0.4
MX2_FIT_MAX = 2.0

OUTPUT_YIELDS_PNG = "output/yields.png"
OUTPUT_MIX_PNG = "output/yields_mix.png"
OUTPUT_DATAONLY_PNG = "output/yields_data_only.png"  # integrated 1-pad three-way
OUTPUT_WEIGHTS_TXT = "output/weights.txt"

OUTPUT_MIXED_MC_ROOT = "output/mixed_mc.root"

# For the new MC mixing files (no runnum):
MIXED_MC_EB = 10.55  # (GeV)

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


def compute_t_scalar_Eb(Eb, e_p, e_theta, e_phi, p_p, p_theta, p_phi):
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

    # q = k - k' = (Eb - E_e, -ex, -ey, Eb - ez)
    # then (q - p_pi)^2 = ( (Eb - E_e) - E_pi )^2 - | (-e_vec + (0,0,Eb)) - p_vec |^2
    dE = (Eb - E_e) - E_pi
    dx = -ex - px
    dy = -ey - py
    dz = (Eb - ez) - pz

    return dE * dE - (dx * dx + dy * dy + dz * dz)
#enddef


def compute_t_scalar(runnum, e_p, e_theta, e_phi, p_p, p_theta, p_phi):
    Eb = beam_energy(int(runnum))
    return compute_t_scalar_Eb(Eb, e_p, e_theta, e_phi, p_p, p_theta, p_phi)
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


def clone_hist_grid(hgrid, prefix):
    nrows = len(hgrid)
    ncols = len(hgrid[0])
    out = []
    for r in range(nrows):
        row = []
        for c in range(ncols):
            h = hgrid[r][c].Clone(f"{see_safe(prefix)}_r{r}_c{c}")
            row.append(h)
        #endfor
        out.append(row)
    #endfor
    return out
#enddef


def see_safe(s):
    # Avoid surprises with ROOT object naming
    return str(s).replace(" ", "_")
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


def sigma_from_raw_counts(n_i, N_pad):
    """
    For normalized data D_i = n_i / N_pad, Poisson stat gives:
      sigma(D_i) = sqrt(n_i) / N_pad.
    We must avoid sigma=0 when n_i=0; use a floor ~ 1/N_pad.
    """
    if N_pad <= 0:
        return 1.0
    #endif
    if n_i <= 0.0:
        return 1.0 / float(N_pad)
    #endif
    return math.sqrt(n_i) / float(N_pad)
#enddef


def compute_best_w_weighted_for_pad(hd_norm, ha_norm, hc_norm, hd_raw):
    """
    Weighted least squares minimization in the fit window:

      chi2(w) = sum_i (D_i - (w*A_i + (1-w)*C_i))^2 / sigma_i^2

    sigma_i is derived from RAW data counts:
      D_i = (n_i / N), sigma_i = sqrt(n_i)/N, with floor for n_i=0.

    Closed form:
      X_i = (A_i - C_i), Y_i = (D_i - C_i), w_i = 1/sigma_i^2
      w_unclipped = sum(w_i X_i Y_i) / sum(w_i X_i^2), then clip to [0,1].
    """

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


def write_weights_report(w_grid, wun_grid, chi2_grid, c_data, c_aao, c_dis, path):
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
        f.write(f"Fit window for w and chi2: Mx2 in [{MX2_FIT_MIN:.3f}, {MX2_FIT_MAX:.3f}]\n")
        f.write("Weights: sigma_i from RAW data counts in the pad: sigma_i = sqrt(n_i)/N (floor for n_i=0)\n")
        f.write("Note: w is clipped to [0,1]. Report includes w_unclipped for diagnostics.\n\n")
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
                f.write(
                    f"  Col {c}: -tprime [{t_lo:.2f}, {t_hi:.2f})  "
                    f"w={w:.6f}  w_unclipped={wun:.6f}  chi2={chi2:.6e}  "
                    f"N(data,aao,dis)=({Nd},{Na},{Nc})\n"
                )
            #endfor
            f.write("\n")
        #endfor
#enddef


def _make_out_branch_double(tout, name):
    arr = array("d", [0.0])
    tout.Branch(name, arr, name + "/D")
    return arr
#enddef


def _make_out_branch_int(tout, name):
    arr = array("i", [0])
    tout.Branch(name, arr, name + "/I")
    return arr
#enddef


def _compute_rc_bin_from_reco(Eb, x_rec, Q2_rec, e_p, e_theta, e_phi, p_p, p_theta, p_phi):
    rb = find_bin(x_rec, XB_EDGES)
    if rb < 0:
        return -1, -1, 0.0, 0.0, 0.0
    #endif

    t_val = compute_t_scalar_Eb(Eb, e_p, e_theta, e_phi, p_p, p_theta, p_phi)
    tmin_val = compute_tmin_exact(x_rec, Q2_rec)
    tprime_val = t_val - tmin_val
    tneg = -tprime_val

    cb = find_bin(tneg, TNEG_EDGES)
    if cb < 0:
        return -1, -1, t_val, tmin_val, tprime_val
    #endif

    return rb, cb, t_val, tmin_val, tprime_val
#enddef


def write_mixed_mc_root(mc_aaogen_path, mc_clasdis_path, w_grid, out_root_path, max_events):
    """
    Build a mixed MC ROOT file using the per-(xB,-tprime) weights w_grid:

      - Write ALL clasdis events (subject to binning cuts)
      - Then write enough aaogen events per bin to approximate target aaogen fraction w

    Target relation per bin:
      w = N_aao / (N_aao + N_dis)  =>  N_aao_target = w/(1-w) * N_dis

    If w=0 => no aaogen added in that bin.
    If w is very close to 1 => extremely large N_aao_target; in that case we write all aaogen
    available and warn if the target is unreachable.

    Also adds derived branches:
      t, tmin, tprime, mc_t, mc_tmin, mc_tprime

    Eb is assumed constant MIXED_MC_EB.
    """

    require_file(mc_aaogen_path)
    require_file(mc_clasdis_path)

    f_aao, t_aao = open_tree(mc_aaogen_path, TREE_NAME)
    f_dis, t_dis = open_tree(mc_clasdis_path, TREE_NAME)

    # Required branches for these special MC files (no runnum).
    # NOTE: fail fast if any are missing.
    needed_mc = [
        "e_p", "mc_e_p",
        "e_theta", "mc_e_theta",
        "e_phi", "mc_e_phi",
        "vz_e", "mc_vz_e",
        "p_p", "mc_p_p",
        "p_theta", "mc_p_theta",
        "p_phi", "mc_p_phi",
        "vz_p", "mc_vz_p",
        "Q2", "mc_Q2",
        "W", "mc_W",
        "Mx", "mc_Mx",
        "Mx2", "mc_Mx2",
        "x", "mc_x",
        "y", "mc_y",
        "z", "mc_z",
        "xF", "mc_xF",
        "pT", "mc_pT",
        "xi", "mc_xi",
        "eta", "mc_eta",
        "trento_phi", "mc_trento_phi",
        "Depolarization_A", "mc_Depolarization_A",
        "Depolarization_B", "mc_Depolarization_B",
        "Depolarization_C", "mc_Depolarization_C",
        "Depolarization_V", "mc_Depolarization_V",
        "Depolarization_W", "mc_Depolarization_W",
        "matching_e_pid",
        "matching_p1_pid",
        "mc_p1_parent"
    ]
    require_branches(t_aao, needed_mc, "mc_aaogen")
    require_branches(t_dis, needed_mc, "mc_clasdis")

    ensure_outdir(out_root_path)

    fout = ROOT.TFile.Open(out_root_path, "RECREATE")
    if not fout or fout.IsZombie():
        fatal("Failed to create output ROOT file: " + out_root_path)
    #endif

    tout = ROOT.TTree(TREE_NAME, "Mixed MC: clasdis + aaogen (per-bin)")

    # Output branches (copy inputs + add derived).
    out = {}

    # Doubles
    for bn in [
        "e_p", "mc_e_p",
        "e_theta", "mc_e_theta",
        "e_phi", "mc_e_phi",
        "vz_e", "mc_vz_e",
        "p_p", "mc_p_p",
        "p_theta", "mc_p_theta",
        "p_phi", "mc_p_phi",
        "vz_p", "mc_vz_p",
        "Q2", "mc_Q2",
        "W", "mc_W",
        "Mx", "mc_Mx",
        "Mx2", "mc_Mx2",
        "x", "mc_x",
        "y", "mc_y",
        "z", "mc_z",
        "xF", "mc_xF",
        "pT", "mc_pT",
        "xi", "mc_xi",
        "eta", "mc_eta",
        "trento_phi", "mc_trento_phi",
        "Depolarization_A", "mc_Depolarization_A",
        "Depolarization_B", "mc_Depolarization_B",
        "Depolarization_C", "mc_Depolarization_C",
        "Depolarization_V", "mc_Depolarization_V",
        "Depolarization_W", "mc_Depolarization_W"
    ]:
        out[bn] = _make_out_branch_double(tout, bn)
    #endfor

    # Ints
    for bn in ["matching_e_pid", "matching_p1_pid", "mc_p1_parent"]:
        out[bn] = _make_out_branch_int(tout, bn)
    #endfor

    # Derived doubles
    out["t"] = _make_out_branch_double(tout, "t")
    out["tmin"] = _make_out_branch_double(tout, "tmin")
    out["tprime"] = _make_out_branch_double(tout, "tprime")
    out["mc_t"] = _make_out_branch_double(tout, "mc_t")
    out["mc_tmin"] = _make_out_branch_double(tout, "mc_tmin")
    out["mc_tprime"] = _make_out_branch_double(tout, "mc_tprime")

    # Reader buffers for input trees
    def _bind_inputs(tree):
        bufs = {}
        tree.SetBranchStatus("*", 0)
        for bn in needed_mc:
            tree.SetBranchStatus(bn, 1)
        #endfor

        # doubles
        for bn in [
            "e_p", "mc_e_p",
            "e_theta", "mc_e_theta",
            "e_phi", "mc_e_phi",
            "vz_e", "mc_vz_e",
            "p_p", "mc_p_p",
            "p_theta", "mc_p_theta",
            "p_phi", "mc_p_phi",
            "vz_p", "mc_vz_p",
            "Q2", "mc_Q2",
            "W", "mc_W",
            "Mx", "mc_Mx",
            "Mx2", "mc_Mx2",
            "x", "mc_x",
            "y", "mc_y",
            "z", "mc_z",
            "xF", "mc_xF",
            "pT", "mc_pT",
            "xi", "mc_xi",
            "eta", "mc_eta",
            "trento_phi", "mc_trento_phi",
            "Depolarization_A", "mc_Depolarization_A",
            "Depolarization_B", "mc_Depolarization_B",
            "Depolarization_C", "mc_Depolarization_C",
            "Depolarization_V", "mc_Depolarization_V",
            "Depolarization_W", "mc_Depolarization_W"
        ]:
            bufs[bn] = array("d", [0.0])
            tree.SetBranchAddress(bn, bufs[bn])
        #endfor

        # ints
        for bn in ["matching_e_pid", "matching_p1_pid", "mc_p1_parent"]:
            bufs[bn] = array("i", [0])
            tree.SetBranchAddress(bn, bufs[bn])
        #endfor

        return bufs
    #enddef

    buf_dis = _bind_inputs(t_dis)
    buf_aao = _bind_inputs(t_aao)

    nrows = len(XB_EDGES) - 1
    ncols = len(TNEG_EDGES) - 1

    # First pass: write ALL clasdis events and count per bin.
    Ndis_written = [[0 for _ in range(ncols)] for _ in range(nrows)]

    n_dis_entries = int(t_dis.GetEntries())
    if max_events is None:
        n_dis_to_process = n_dis_entries
    else:
        n_dis_to_process = min(n_dis_entries, int(max_events))
    #endif

    for i in range(n_dis_to_process):
        t_dis.GetEntry(i)

        x_rec = float(buf_dis["x"][0])
        Q2_rec = float(buf_dis["Q2"][0])

        rb, cb, t_val, tmin_val, tprime_val = _compute_rc_bin_from_reco(
            MIXED_MC_EB,
            x_rec, Q2_rec,
            float(buf_dis["e_p"][0]), float(buf_dis["e_theta"][0]), float(buf_dis["e_phi"][0]),
            float(buf_dis["p_p"][0]), float(buf_dis["p_theta"][0]), float(buf_dis["p_phi"][0])
        )
        if rb < 0 or cb < 0:
            continue
        #endif

        # Fill output copied branches
        for bn in out.keys():
            if bn in ["t", "tmin", "tprime", "mc_t", "mc_tmin", "mc_tprime"]:
                continue
            #endif
            if bn in ["matching_e_pid", "matching_p1_pid", "mc_p1_parent"]:
                out[bn][0] = int(buf_dis[bn][0])
            else:
                out[bn][0] = float(buf_dis[bn][0])
            #endif
        #endfor

        # Derived (reco)
        out["t"][0] = float(t_val)
        out["tmin"][0] = float(tmin_val)
        out["tprime"][0] = float(tprime_val)

        # Derived (mc)
        mc_x = float(buf_dis["mc_x"][0])
        mc_Q2 = float(buf_dis["mc_Q2"][0])
        mc_t_val = compute_t_scalar_Eb(
            MIXED_MC_EB,
            float(buf_dis["mc_e_p"][0]), float(buf_dis["mc_e_theta"][0]), float(buf_dis["mc_e_phi"][0]),
            float(buf_dis["mc_p_p"][0]), float(buf_dis["mc_p_theta"][0]), float(buf_dis["mc_p_phi"][0])
        )
        mc_tmin_val = compute_tmin_exact(mc_x, mc_Q2)
        mc_tprime_val = mc_t_val - mc_tmin_val

        out["mc_t"][0] = float(mc_t_val)
        out["mc_tmin"][0] = float(mc_tmin_val)
        out["mc_tprime"][0] = float(mc_tprime_val)

        tout.Fill()
        Ndis_written[rb][cb] += 1
    #endfor

    # Decide how many aaogen events to add per bin
    Naao_need = [[0 for _ in range(ncols)] for _ in range(nrows)]
    for r in range(nrows):
        for c in range(ncols):
            Nc = int(Ndis_written[r][c])
            w = float(w_grid[r][c])
            if Nc <= 0:
                Naao_need[r][c] = 0
                continue
            #endif
            if w <= 0.0:
                Naao_need[r][c] = 0
                continue
            #endif
            if w >= 0.999999:
                # Practically "all aaogen"; cannot be achieved with finite Nc unless Nc->0.
                # We will just try to write all aaogen events in that bin.
                Naao_need[r][c] = 2**31 - 1
                continue
            #endif

            # N_aao_target = w/(1-w)*Nc
            target = (w / (1.0 - w)) * float(Nc)
            Naao_need[r][c] = int(math.ceil(target))
        #endfor
    #endfor

    # Second pass: write aaogen until Naao_need per bin is satisfied (or aaogen exhausts).
    Naao_written = [[0 for _ in range(ncols)] for _ in range(nrows)]

    n_aao_entries = int(t_aao.GetEntries())
    if max_events is None:
        n_aao_to_process = n_aao_entries
    else:
        n_aao_to_process = min(n_aao_entries, int(max_events))
    #endif

    remaining_total = 0
    for r in range(nrows):
        for c in range(ncols):
            if Naao_need[r][c] > 0 and Naao_need[r][c] < (2**31 - 1):
                remaining_total += Naao_need[r][c]
            elif Naao_need[r][c] >= (2**31 - 1):
                remaining_total += 1
            #endif
        #endfor
    #endfor

    for i in range(n_aao_to_process):
        if remaining_total <= 0:
            break
        #endif

        t_aao.GetEntry(i)

        x_rec = float(buf_aao["x"][0])
        Q2_rec = float(buf_aao["Q2"][0])

        rb, cb, t_val, tmin_val, tprime_val = _compute_rc_bin_from_reco(
            MIXED_MC_EB,
            x_rec, Q2_rec,
            float(buf_aao["e_p"][0]), float(buf_aao["e_theta"][0]), float(buf_aao["e_phi"][0]),
            float(buf_aao["p_p"][0]), float(buf_aao["p_theta"][0]), float(buf_aao["p_phi"][0])
        )
        if rb < 0 or cb < 0:
            continue
        #endif

        need = int(Naao_need[rb][cb])
        have = int(Naao_written[rb][cb])
        if have >= need:
            continue
        #endif

        # Fill output copied branches
        for bn in out.keys():
            if bn in ["t", "tmin", "tprime", "mc_t", "mc_tmin", "mc_tprime"]:
                continue
            #endif
            if bn in ["matching_e_pid", "matching_p1_pid", "mc_p1_parent"]:
                out[bn][0] = int(buf_aao[bn][0])
            else:
                out[bn][0] = float(buf_aao[bn][0])
            #endif
        #endfor

        # Derived (reco)
        out["t"][0] = float(t_val)
        out["tmin"][0] = float(tmin_val)
        out["tprime"][0] = float(tprime_val)

        # Derived (mc)
        mc_x = float(buf_aao["mc_x"][0])
        mc_Q2 = float(buf_aao["mc_Q2"][0])
        mc_t_val = compute_t_scalar_Eb(
            MIXED_MC_EB,
            float(buf_aao["mc_e_p"][0]), float(buf_aao["mc_e_theta"][0]), float(buf_aao["mc_e_phi"][0]),
            float(buf_aao["mc_p_p"][0]), float(buf_aao["mc_p_theta"][0]), float(buf_aao["mc_p_phi"][0])
        )
        mc_tmin_val = compute_tmin_exact(mc_x, mc_Q2)
        mc_tprime_val = mc_t_val - mc_tmin_val

        out["mc_t"][0] = float(mc_t_val)
        out["mc_tmin"][0] = float(mc_tmin_val)
        out["mc_tprime"][0] = float(mc_tprime_val)

        tout.Fill()
        Naao_written[rb][cb] += 1

        # update remaining_total conservatively
        if need < (2**31 - 1):
            remaining_total -= 1
        #endif
    #endfor

    # Diagnostics: achieved fractions per bin
    print("Mixed MC writing summary (per bin):")
    for r in range(nrows):
        for c in range(ncols):
            Nc = int(Ndis_written[r][c])
            Na = int(Naao_written[r][c])
            w_target = float(w_grid[r][c])

            denom = float(Na + Nc)
            if denom > 0.0:
                w_ach = float(Na) / denom
            else:
                w_ach = 0.0
            #endif

            if Nc > 0 or Na > 0:
                msg = f"  bin (r={r}, c={c}): w_target={w_target:.6f}  Ndis={Nc}  Naao={Na}  w_ach={w_ach:.6f}"
                if w_target > 0.0 and Na < Naao_need[r][c] and Naao_need[r][c] < (2**31 - 1):
                    msg += "  WARNING: aaogen insufficient for target"
                #endif
                if Naao_need[r][c] >= (2**31 - 1) and w_target >= 0.999999:
                    msg += "  NOTE: w_target~1; wrote as many aaogen as available"
                #endif
                print(msg)
            #endif
        #endfor
    #endfor

    fout.cd()
    tout.Write()
    fout.Close()

    t_aao.SetBranchStatus("*", 1)
    t_dis.SetBranchStatus("*", 1)

    f_aao.Close()
    f_dis.Close()

    print("Wrote mixed MC ROOT file: " + out_root_path)
#enddef


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True, help="Path to data ROOT file")
    ap.add_argument("--aaogen", required=True, help="Path to aaogen ROOT file")
    ap.add_argument("--clasdis", required=True, help="Path to clasdis ROOT file")
    ap.add_argument("--max_events", type=int, default=-1,
                    help="Max events per file (-1 means all events)")
    ap.add_argument("--mc_aaogen", default="", help="(optional) Path to special MC aaogen ROOT file to be mixed")
    ap.add_argument("--mc_clasdis", default="", help="(optional) Path to special MC clasdis ROOT file to be mixed")
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
    ensure_outdir(OUTPUT_MIXED_MC_ROOT)

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

    # Keep RAW data histograms for computing sigma_i in the weighted fit.
    h_data_raw = []
    for r in range(nrows):
        row = []
        for c in range(ncols):
            row.append(h_data[r][c].Clone(f"h_data_raw_r{r}_c{c}"))
        #endfor
        h_data_raw.append(row)
    #endfor

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

    # 1-pad integrated three-way plot.
    draw_canvas_integrated_threeway(h_data_int, h_aao_int, h_dis_int, Nd_tot, Na_tot, Nc_tot, OUTPUT_DATAONLY_PNG)

    # Weighted per-pad mixture fit (weights from RAW data stats).
    w_grid, wun_grid, chi2_grid, h_mix = compute_w_grid_and_mix(h_data, h_aao, h_dis, h_data_raw)

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
    print(f"  fit window: Mx2 in [{MX2_FIT_MIN:.3f}, {MX2_FIT_MAX:.3f}]")
    print(f"  total chi2 = {total_chi2:.6e}")
    print(f"  wrote weights report: {OUTPUT_WEIGHTS_TXT}")

    write_weights_report(w_grid, wun_grid, chi2_grid, c_data, c_aao, c_dis, OUTPUT_WEIGHTS_TXT)

    draw_canvas_mix(h_data, h_mix, c_data, w_grid, OUTPUT_MIX_PNG)

    # Phase 3 (optional): write mixed MC ROOT using provided special MC inputs
    mc_aao = str(args.mc_aaogen).strip()
    mc_dis = str(args.mc_clasdis).strip()
    if (mc_aao != "") or (mc_dis != ""):
        if mc_aao == "" or mc_dis == "":
            fatal("If using Phase 3, you must provide BOTH --mc_aaogen and --mc_clasdis.")
        #endif
        print("Phase 3: building mixed MC ROOT from:")
        print("  mc_aaogen  = " + mc_aao)
        print("  mc_clasdis = " + mc_dis)
        print("  out_root   = " + OUTPUT_MIXED_MC_ROOT)
        print("  Eb (fixed) = %.3f GeV" % MIXED_MC_EB)
        write_mixed_mc_root(mc_aao, mc_dis, w_grid, OUTPUT_MIXED_MC_ROOT, max_events)
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