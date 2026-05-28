#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
normalization.py

Purpose
-------
Build a mixed aaogen + clasdis MC sample for the enpi+X missing-mass study.

Inputs
------
  --data     data ROOT file
  --aaogen   aaogen ROOT file
  --clasdis  clasdis ROOT file
  --out      output mixed MC ROOT file

Main workflow
-------------
Phase 1:
  - Fill 4x6 Mx2 histograms binned in (xB, -tprime) for data, aaogen, clasdis.
  - Compute t, tmin, and tprime on the fly.
      data:           use runnum-based beam energy
      aaogen/clasdis: use fixed MC_EB_FIXED
  - Normalize each histogram to unit area per (xB, -tprime) bin.

Phase 2:
  - Determine a per-bin mixture fraction w[row][col]:
        H_mix = w * H_aaogen + (1 - w) * H_clasdis
  - Choose w by minimizing a weighted chi2 in the peak Mx2 window.
  - The data statistical uncertainty comes from the raw data counts in each Mx2 bin.

Phase 3:
  - Always write a mixed output ROOT file using the same aaogen and clasdis inputs.
  - Write all clasdis events.
  - Add aaogen events in each (xB, -tprime) bin until the target fraction is reached:
        w = N_aao / (N_aao + N_dis)
        N_aao_target = w / (1 - w) * N_dis
  - Recompute and write:
        t, tmin, tprime
        mc_t, mc_tmin, mc_tprime
  - Use fixed MC_EB_FIXED for aaogen/clasdis and generated-level quantities.
  - No runnum is required for aaogen or clasdis.

Outputs
-------
  output/yields.png
  output/yields_data_only.png
  output/yields_mix.png
  output/weights.txt
  output/mix_debug_report.txt
  output/mix_debug_mx2.png
  --out mixed MC ROOT file

Optional diagnostics
--------------------
  --force:
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

XB_EDGES = [0.10, 0.25, 0.35, 0.45, 0.60]
TNEG_EDGES = [0.05, 0.25, 0.45, 0.65, 0.85, 1.05, 1.25]

MX2_MIN = 0.3
MX2_MAX = 2.0
MX2_NBINS = 33

MX2_FIT_MIN = 0.81
MX2_FIT_MAX = 1.10

OUTPUT_YIELDS_PNG = "output/yields.png"
OUTPUT_MIX_PNG = "output/yields_mix.png"
OUTPUT_DATAONLY_PNG = "output/yields_data_only.png"
OUTPUT_WEIGHTS_TXT = "output/weights.txt"

DEFAULT_MIXED_MC_ROOT = "output/mixed_mc.root"

OUTPUT_MIX_DEBUG_TXT = "output/mix_debug_report.txt"
OUTPUT_MIX_DEBUG_MX2_PNG = "output/mix_debug_mx2.png"

MC_EB_FIXED = 10.55

MASS_E = 0.000511
MASS_PI = 0.139570
MASS_N = 0.9382720813


PHASE1_BRANCHES_BASE = [
    "e_p",
    "e_theta",
    "e_phi",
    "p_p",
    "p_theta",
    "p_phi",
    "x",
    "Q2",
    "Mx2"
]


PHASE3_DOUBLE_BRANCHES = [
    "e_p", "e_theta", "e_phi", "vz_e",
    "p_p", "p_theta", "p_phi", "vz_p",
    "Q2", "W", "Mx", "Mx2", "x", "y", "z", "xF", "pT", "xi", "eta", "phi",
    "DepA", "DepB", "DepC", "DepV", "DepW",
    "mc_e_p", "mc_e_theta", "mc_e_phi", "mc_vz_e",
    "mc_p_p", "mc_p_theta", "mc_p_phi", "mc_vz_p",
    "mc_Q2", "mc_W", "mc_Mx", "mc_Mx2", "mc_x", "mc_y", "mc_z",
    "mc_xF", "mc_pT", "mc_xi", "mc_eta", "mc_phi",
    "mc_DepA", "mc_DepB", "mc_DepC", "mc_DepV", "mc_DepW",
    "t", "tmin", "tprime", "mc_t", "mc_tmin", "mc_tprime"
]


PHASE3_INPUT_DOUBLE_BRANCHES = [
    "e_p", "e_theta", "e_phi", "vz_e",
    "p_p", "p_theta", "p_phi", "vz_p",
    "Q2", "W", "Mx", "Mx2", "x", "y", "z", "xF", "pT", "xi", "eta", "phi",
    "DepA", "DepB", "DepC", "DepV", "DepW",
    "mc_e_p", "mc_e_theta", "mc_e_phi", "mc_vz_e",
    "mc_p_p", "mc_p_theta", "mc_p_phi", "mc_vz_p",
    "mc_Q2", "mc_W", "mc_Mx", "mc_Mx2", "mc_x", "mc_y", "mc_z",
    "mc_xF", "mc_pT", "mc_xi", "mc_eta", "mc_phi",
    "mc_DepA", "mc_DepB", "mc_DepC", "mc_DepV", "mc_DepW"
]


PHASE3_INT_BRANCHES = [
    "matching_e_pid",
    "matching_p1_pid",
    "mc_p1_parent"
]


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


def ensure_outdir(path):
    d = os.path.dirname(str(path))
    if d != "":
        os.makedirs(d, exist_ok=True)
    #endif
#enddef


def open_tree(path, tree_name):
    fobj = ROOT.TFile.Open(path, "READ")
    if not fobj or fobj.IsZombie():
        fatal("Failed to open ROOT file: " + str(path))
    #endif

    tree = fobj.Get(tree_name)
    if not tree:
        fatal("Tree '" + str(tree_name) + "' not found in: " + str(path))
    #endif

    return fobj, tree
#enddef


def require_branches(tree, needed, label):
    blist = tree.GetListOfBranches()
    if not blist:
        fatal("No branch list available for tree in: " + str(label))
    #endif

    missing = []
    for bname in needed:
        if not blist.FindObject(bname):
            missing.append(bname)
        #endif
    #endfor

    if len(missing) > 0:
        fatal("Missing required branches in " + str(label) + ": " + ", ".join(missing))
    #endif
#enddef


def parse_max_events(val):
    if val is None:
        return None
    #endif

    ival = int(val)
    if ival < 0:
        return None
    #endif

    return ival
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
            denom = 1.0 - xB
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


def normalize_unit_area(hist):
    integral = hist.Integral(1, hist.GetNbinsX())
    if integral > 0.0:
        hist.Scale(1.0 / integral)
    #endif
#enddef


def style_hist(hist, color, width, linestyle):
    hist.SetLineColor(color)
    hist.SetLineWidth(width)
    hist.SetLineStyle(linestyle)
    hist.SetMarkerSize(0.0)
#enddef


def make_hist_grid(prefix, nrows, ncols):
    grid = []

    for r in range(nrows):
        row = []

        for c in range(ncols):
            name = f"{prefix}_r{r}_c{c}"
            hist = ROOT.TH1F(name, "", MX2_NBINS, MX2_MIN, MX2_MAX)
            hist.Sumw2()
            row.append(hist)
        #endfor

        grid.append(row)
    #endfor

    return grid
#enddef


def make_count_grid(nrows, ncols):
    return [[0 for _ in range(ncols)] for _ in range(nrows)]
#enddef


def make_phase1_buffers(use_runnum):
    buffers = {
        "e_p": array("d", [0.0]),
        "e_theta": array("d", [0.0]),
        "e_phi": array("d", [0.0]),
        "p_p": array("d", [0.0]),
        "p_theta": array("d", [0.0]),
        "p_phi": array("d", [0.0]),
        "x": array("d", [0.0]),
        "Q2": array("d", [0.0]),
        "Mx2": array("d", [0.0])
    }

    if use_runnum:
        buffers["runnum"] = array("i", [0])
    #endif

    return buffers
#enddef


def bind_phase1_tree(tree, buffers, use_runnum):
    tree.SetBranchStatus("*", 0)

    needed = list(PHASE1_BRANCHES_BASE)
    if use_runnum:
        needed.append("runnum")
    #endif

    for bname in needed:
        tree.SetBranchStatus(bname, 1)
    #endfor

    for bname in needed:
        tree.SetBranchAddress(bname, buffers[bname])
    #endfor
#enddef


def fill_all_bins_single_pass(tree, hgrid, cgrid, max_events, use_runnum, fixed_Eb):
    buffers = make_phase1_buffers(use_runnum)
    bind_phase1_tree(tree, buffers, use_runnum)

    n_entries = int(tree.GetEntries())
    if max_events is None:
        n_to_process = n_entries
    else:
        n_to_process = min(n_entries, int(max_events))
    #endif

    for ientry in range(n_to_process):
        tree.GetEntry(ientry)

        xb_val = float(buffers["x"][0])
        row = find_bin(xb_val, XB_EDGES)
        if row < 0:
            continue
        #endif

        Q2_val = float(buffers["Q2"][0])

        if use_runnum:
            t_val = compute_t_scalar(
                int(buffers["runnum"][0]),
                float(buffers["e_p"][0]),
                float(buffers["e_theta"][0]),
                float(buffers["e_phi"][0]),
                float(buffers["p_p"][0]),
                float(buffers["p_theta"][0]),
                float(buffers["p_phi"][0])
            )
        else:
            t_val = compute_t_scalar_from_Eb(
                float(fixed_Eb),
                float(buffers["e_p"][0]),
                float(buffers["e_theta"][0]),
                float(buffers["e_phi"][0]),
                float(buffers["p_p"][0]),
                float(buffers["p_theta"][0]),
                float(buffers["p_phi"][0])
            )
        #endif

        tmin_val = compute_tmin_exact(xb_val, Q2_val)
        tprime_val = t_val - tmin_val
        tneg_val = -tprime_val

        col = find_bin(tneg_val, TNEG_EDGES)
        if col < 0:
            continue
        #endif

        hgrid[row][col].Fill(float(buffers["Mx2"][0]))
        cgrid[row][col] += 1
    #endfor

    tree.SetBranchStatus("*", 1)
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
    nbins = hd_norm.GetNbinsX()

    N_pad = hd_raw.Integral(1, nbins)
    if N_pad <= 0.0 or hd_norm.Integral(1, nbins) <= 0.0:
        return 0.0, 0.0, 0.0
    #endif

    num = 0.0
    den = 0.0
    used_bins = 0

    for ibin in range(1, nbins + 1):
        xcen = hd_norm.GetXaxis().GetBinCenter(ibin)

        if xcen < MX2_FIT_MIN or xcen > MX2_FIT_MAX:
            continue
        #endif

        D = hd_norm.GetBinContent(ibin)
        A = ha_norm.GetBinContent(ibin)
        C = hc_norm.GetBinContent(ibin)

        n_i = hd_raw.GetBinContent(ibin)
        sig = sigma_from_raw_counts(n_i, N_pad)

        if sig <= 0.0:
            continue
        #endif

        weight = 1.0 / (sig * sig)
        X = A - C
        Y = D - C

        num += weight * X * Y
        den += weight * X * X
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

    chi2 = compute_chi2_for_pad_given_w(hd_norm, ha_norm, hc_norm, hd_raw, w)

    return w, w_unclipped, chi2
#enddef


def compute_chi2_for_pad_given_w(hd_norm, ha_norm, hc_norm, hd_raw, w):
    nbins = hd_norm.GetNbinsX()

    N_pad = hd_raw.Integral(1, nbins)
    if N_pad <= 0.0 or hd_norm.Integral(1, nbins) <= 0.0:
        return 0.0
    #endif

    chi2 = 0.0

    for ibin in range(1, nbins + 1):
        xcen = hd_norm.GetXaxis().GetBinCenter(ibin)

        if xcen < MX2_FIT_MIN or xcen > MX2_FIT_MAX:
            continue
        #endif

        D = hd_norm.GetBinContent(ibin)
        M = w * ha_norm.GetBinContent(ibin) + (1.0 - w) * hc_norm.GetBinContent(ibin)

        n_i = hd_raw.GetBinContent(ibin)
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
            hm.Add(hc, 1.0 - w)

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
    total = 0

    for r in range(len(cgrid)):
        for c in range(len(cgrid[r])):
            total += int(cgrid[r][c])
        #endfor
    #endfor

    return total
#enddef


def make_integrated_hist_from_grid(hgrid, name):
    hsum = hgrid[0][0].Clone(name)
    hsum.Reset("ICESM")

    for r in range(len(hgrid)):
        for c in range(len(hgrid[r])):
            hsum.Add(hgrid[r][c])
        #endfor
    #endfor

    return hsum
#enddef


def draw_canvas_integrated_threeway(hd_int, ha_int, hc_int, Nd, Na, Nc, outpng):
    canv = ROOT.TCanvas("c_integrated", "Mx2 integrated: data vs aaogen vs clasdis", 1200, 900)
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

    ymax *= 1.2

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

            ymax *= 1.2

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

            ymax *= 1.2

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

    with open(path, "w") as fout:
        fout.write("Per-bin mixture weights (shape-only, WEIGHTED chi2)\n")
        fout.write("Definition: H_mix = w * H_aaogen + (1-w) * H_clasdis\n")
        fout.write(f"Peak window for w and chi2: Mx2 in [{MX2_FIT_MIN:.3f}, {MX2_FIT_MAX:.3f}]\n")
        fout.write(f"Histogram binning: {MX2_NBINS} bins from {MX2_MIN:.3f} to {MX2_MAX:.3f} GeV^2\n")
        fout.write("Weights: sigma_i from RAW data counts in the pad: sigma_i = sqrt(n_i)/N (floor for n_i=0)\n")
        fout.write("Note: w is clipped to [0,1]. Report includes w_unclipped for diagnostics.\n")

        if forced_map is not None and len(forced_map) > 0:
            fout.write("NOTE: Some w values were FORCED by --force and chi2 was recomputed for those pads.\n")
        #endif

        fout.write("\n")
        fout.write(f"Total chi2 (sum over pads) = {total_chi2:.6e}\n\n")

        for r in range(nrows):
            xb_lo = XB_EDGES[r]
            xb_hi = XB_EDGES[r + 1]
            fout.write(f"Row {r}: xB [{xb_lo:.2f}, {xb_hi:.2f})\n")

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
                if forced_map is not None and (r, c) in forced_map:
                    tag = "  FORCED"
                #endif

                fout.write(
                    f"  Col {c}: -tprime [{t_lo:.2f}, {t_hi:.2f})  "
                    f"w={w:.6f}  w_unclipped={wun:.6f}  chi2={chi2:.6e}  "
                    f"N(data,aao,dis)=({Nd},{Na},{Nc}){tag}\n"
                )
            #endfor

            fout.write("\n")
        #endfor
    #endwith
#enddef


def compute_bin_indices_from_reco(x_val, Q2_val, e_p, e_theta, e_phi, p_p, p_theta, p_phi):
    t_val = compute_t_scalar_from_Eb(MC_EB_FIXED, e_p, e_theta, e_phi, p_p, p_theta, p_phi)
    tmin_val = compute_tmin_exact(x_val, Q2_val)
    tprime_val = t_val - tmin_val
    tneg_val = -tprime_val

    row = find_bin(x_val, XB_EDGES)
    if row < 0:
        return -1, -1, t_val, tmin_val, tprime_val
    #endif

    col = find_bin(tneg_val, TNEG_EDGES)
    if col < 0:
        return -1, -1, t_val, tmin_val, tprime_val
    #endif

    return row, col, t_val, tmin_val, tprime_val
#enddef


def compute_mc_t_quantities_from_gen(mc_x_val, mc_Q2_val, mc_e_p, mc_e_theta, mc_e_phi, mc_p_p, mc_p_theta, mc_p_phi):
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


def make_phase3_buffers():
    buffers = {}

    for bname in PHASE3_DOUBLE_BRANCHES:
        buffers[bname] = array("d", [0.0])
    #endfor

    for bname in PHASE3_INT_BRANCHES:
        buffers[bname] = array("i", [0])
    #endfor

    return buffers
#enddef


def bind_phase3_input_tree(tree, buffers):
    tree.SetBranchStatus("*", 0)

    for bname in PHASE3_INPUT_DOUBLE_BRANCHES:
        tree.SetBranchStatus(bname, 1)
    #endfor

    for bname in PHASE3_INT_BRANCHES:
        tree.SetBranchStatus(bname, 1)
    #endfor

    for bname in PHASE3_INPUT_DOUBLE_BRANCHES:
        tree.SetBranchAddress(bname, buffers[bname])
    #endfor

    for bname in PHASE3_INT_BRANCHES:
        tree.SetBranchAddress(bname, buffers[bname])
    #endfor
#enddef


def make_phase3_output_tree():
    buffers = make_phase3_buffers()
    tout = ROOT.TTree(TREE_NAME, "Mixed MC: write all clasdis; top up aaogen in-grid using w_grid")

    for bname in PHASE3_DOUBLE_BRANCHES:
        tout.Branch(bname, buffers[bname], bname + "/D")
    #endfor

    for bname in PHASE3_INT_BRANCHES:
        tout.Branch(bname, buffers[bname], bname + "/I")
    #endfor

    return tout, buffers
#enddef


def fill_recomputed_t_buffers(buffers):
    row, col, t_val, tmin_val, tprime_val = compute_bin_indices_from_reco(
        float(buffers["x"][0]),
        float(buffers["Q2"][0]),
        float(buffers["e_p"][0]),
        float(buffers["e_theta"][0]),
        float(buffers["e_phi"][0]),
        float(buffers["p_p"][0]),
        float(buffers["p_theta"][0]),
        float(buffers["p_phi"][0])
    )

    mc_t_val, mc_tmin_val, mc_tprime_val = compute_mc_t_quantities_from_gen(
        float(buffers["mc_x"][0]),
        float(buffers["mc_Q2"][0]),
        float(buffers["mc_e_p"][0]),
        float(buffers["mc_e_theta"][0]),
        float(buffers["mc_e_phi"][0]),
        float(buffers["mc_p_p"][0]),
        float(buffers["mc_p_theta"][0]),
        float(buffers["mc_p_phi"][0])
    )

    buffers["t"][0] = float(t_val)
    buffers["tmin"][0] = float(tmin_val)
    buffers["tprime"][0] = float(tprime_val)

    buffers["mc_t"][0] = float(mc_t_val)
    buffers["mc_tmin"][0] = float(mc_tmin_val)
    buffers["mc_tprime"][0] = float(mc_tprime_val)

    return row, col
#enddef


def write_mix_debug_report(path, w_grid, Ndis, Naao_target, Naao_written,
                           clasdis_written_total, clasdis_written_in_grid, clasdis_written_out_grid,
                           aaogen_written_total, aaogen_written_in_grid,
                           aaogen_skipped_out_grid, aaogen_skipped_quota_met):
    ensure_outdir(path)

    nrows = len(w_grid)
    ncols = len(w_grid[0])

    with open(path, "w") as fout:
        fout.write("Phase 3 debug report\n")
        fout.write("All clasdis events were written. Aaogen was topped up in-grid according to w_grid.\n")
        fout.write("No Mx2 cuts were applied during output-tree writing.\n")
        fout.write("\n")
        fout.write(f"clasdis written total       = {clasdis_written_total}\n")
        fout.write(f"clasdis written in-grid     = {clasdis_written_in_grid}\n")
        fout.write(f"clasdis written out-of-grid = {clasdis_written_out_grid}\n")
        fout.write(f"aaogen written total        = {aaogen_written_total}\n")
        fout.write(f"aaogen written in-grid      = {aaogen_written_in_grid}\n")
        fout.write(f"aaogen skipped out-of-grid  = {aaogen_skipped_out_grid}\n")
        fout.write(f"aaogen skipped quota met    = {aaogen_skipped_quota_met}\n")
        fout.write("\n")

        denom = float(clasdis_written_total + aaogen_written_total)
        if denom > 0.0:
            fout.write(f"global achieved aaogen fraction = {aaogen_written_total / denom:.6f}\n")
        else:
            fout.write("global achieved aaogen fraction = 0.000000\n")
        #endif

        fout.write("\n")

        for r in range(nrows):
            xb_lo = XB_EDGES[r]
            xb_hi = XB_EDGES[r + 1]
            fout.write(f"Row {r}: xB [{xb_lo:.2f}, {xb_hi:.2f})\n")

            for c in range(ncols):
                t_lo = TNEG_EDGES[c]
                t_hi = TNEG_EDGES[c + 1]

                nd = int(Ndis[r][c])
                nt = int(Naao_target[r][c])
                nw = int(Naao_written[r][c])

                achieved = 0.0
                if (nd + nw) > 0:
                    achieved = float(nw) / float(nd + nw)
                #endif

                fout.write(
                    f"  Col {c}: -tprime [{t_lo:.2f}, {t_hi:.2f})  "
                    f"w_target={float(w_grid[r][c]):.6f}  "
                    f"Ndis={nd}  "
                    f"Naao_target={nt}  "
                    f"Naao_written={nw}  "
                    f"w_achieved={achieved:.6f}\n"
                )
            #endfor

            fout.write("\n")
        #endfor
    #endwith
#enddef


def draw_mix_debug_mx2(h_written_dis, h_written_aao, h_written_mix, outpng):
    h_dis_norm = h_written_dis.Clone("h_written_dis_norm")
    h_aao_norm = h_written_aao.Clone("h_written_aao_norm")
    h_mix_norm = h_written_mix.Clone("h_written_mix_norm")

    normalize_unit_area(h_dis_norm)
    normalize_unit_area(h_aao_norm)
    normalize_unit_area(h_mix_norm)

    style_hist(h_dis_norm, ROOT.kBlue, 2, 3)
    style_hist(h_aao_norm, ROOT.kRed, 2, 2)
    style_hist(h_mix_norm, ROOT.kBlack, 3, 1)

    canv = ROOT.TCanvas("c_mix_debug_mx2", "Written MC debug", 1200, 900)
    pad = canv.cd(1)
    pad.SetGrid(1, 1)
    pad.SetLeftMargin(0.18)
    pad.SetRightMargin(0.05)
    pad.SetBottomMargin(0.14)
    pad.SetTopMargin(0.08)

    ymax = max(h_dis_norm.GetMaximum(), h_aao_norm.GetMaximum(), h_mix_norm.GetMaximum())
    if ymax <= 0.0:
        ymax = 1.0
    #endif

    ymax *= 1.2

    set_axes_and_range(h_mix_norm, ymax)

    h_mix_norm.Draw("hist")
    h_dis_norm.Draw("hist same")
    h_aao_norm.Draw("hist same")

    leg = ROOT.TLegend(0.55, 0.72, 0.94, 0.92)
    leg.SetBorderSize(1)
    leg.SetFillStyle(1001)
    leg.SetFillColor(ROOT.kWhite)
    leg.SetTextSize(0.042)
    leg.AddEntry(h_mix_norm, f"written total (N={int(h_written_mix.Integral())})", "l")
    leg.AddEntry(h_dis_norm, f"written clasdis (N={int(h_written_dis.Integral())})", "l")
    leg.AddEntry(h_aao_norm, f"written aaogen (N={int(h_written_aao.Integral())})", "l")
    leg.Draw()

    tex = ROOT.TLatex()
    tex.SetNDC(True)
    tex.SetTextSize(0.045)
    tex.DrawLatex(0.16, 0.93, "Integrated M_{x}^{2} shapes of events written to mixed MC")

    ensure_outdir(outpng)
    canv.SaveAs(outpng)
#enddef


def mix_mc_to_output_root(aaogen_path, clasdis_path, out_root_path, w_grid, max_events):
    require_file(aaogen_path)
    require_file(clasdis_path)

    f_aao, t_aao = open_tree(aaogen_path, TREE_NAME)
    f_dis, t_dis = open_tree(clasdis_path, TREE_NAME)

    phase3_needed = list(PHASE3_INPUT_DOUBLE_BRANCHES) + list(PHASE3_INT_BRANCHES)

    require_branches(t_aao, phase3_needed, "aaogen phase3")
    require_branches(t_dis, phase3_needed, "clasdis phase3")

    ensure_outdir(out_root_path)

    fout = ROOT.TFile.Open(out_root_path, "RECREATE")
    if not fout or fout.IsZombie():
        fatal("Failed to create output ROOT file: " + str(out_root_path))
    #endif

    tout, buffers = make_phase3_output_tree()

    bind_phase3_input_tree(t_dis, buffers)

    nrows = len(XB_EDGES) - 1
    ncols = len(TNEG_EDGES) - 1

    Ndis = [[0 for _ in range(ncols)] for _ in range(nrows)]
    Naao_target = [[0 for _ in range(ncols)] for _ in range(nrows)]
    Naao_written = [[0 for _ in range(ncols)] for _ in range(nrows)]

    clasdis_written_total = 0
    clasdis_written_in_grid = 0
    clasdis_written_out_grid = 0

    aaogen_written_total = 0
    aaogen_written_in_grid = 0
    aaogen_skipped_out_grid = 0
    aaogen_skipped_quota_met = 0

    h_written_dis = ROOT.TH1F("h_written_dis", "", MX2_NBINS, MX2_MIN, MX2_MAX)
    h_written_aao = ROOT.TH1F("h_written_aao", "", MX2_NBINS, MX2_MIN, MX2_MAX)
    h_written_mix = ROOT.TH1F("h_written_mix", "", MX2_NBINS, MX2_MIN, MX2_MAX)

    n_entries_dis = int(t_dis.GetEntries())
    if max_events is None:
        n_to_process_dis = n_entries_dis
    else:
        n_to_process_dis = min(n_entries_dis, int(max_events))
    #endif

    for ientry in range(n_to_process_dis):
        t_dis.GetEntry(ientry)

        row, col = fill_recomputed_t_buffers(buffers)

        tout.Fill()
        clasdis_written_total += 1

        mx2_val = float(buffers["Mx2"][0])
        h_written_dis.Fill(mx2_val)
        h_written_mix.Fill(mx2_val)

        if row >= 0 and col >= 0:
            Ndis[row][col] += 1
            clasdis_written_in_grid += 1
        else:
            clasdis_written_out_grid += 1
        #endif
    #endfor

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
                Naao_target[r][c] = 10 ** 18
                continue
            #endif

            need = (w / (1.0 - w)) * float(nd)
            Naao_target[r][c] = int(math.ceil(need))
        #endfor
    #endfor

    bind_phase3_input_tree(t_aao, buffers)

    n_entries_aao = int(t_aao.GetEntries())
    if max_events is None:
        n_to_process_aao = n_entries_aao
    else:
        n_to_process_aao = min(n_entries_aao, int(max_events))
    #endif

    for ientry in range(n_to_process_aao):
        t_aao.GetEntry(ientry)

        row, col = fill_recomputed_t_buffers(buffers)

        if row < 0 or col < 0:
            aaogen_skipped_out_grid += 1
            continue
        #endif

        if Naao_written[row][col] >= Naao_target[row][col]:
            aaogen_skipped_quota_met += 1
            continue
        #endif

        tout.Fill()
        Naao_written[row][col] += 1

        aaogen_written_total += 1
        aaogen_written_in_grid += 1

        mx2_val = float(buffers["Mx2"][0])
        h_written_aao.Fill(mx2_val)
        h_written_mix.Fill(mx2_val)
    #endfor

    fout.cd()
    tout.Write()
    h_written_dis.Write()
    h_written_aao.Write()
    h_written_mix.Write()
    fout.Close()

    f_aao.Close()
    f_dis.Close()

    write_mix_debug_report(
        OUTPUT_MIX_DEBUG_TXT,
        w_grid,
        Ndis,
        Naao_target,
        Naao_written,
        clasdis_written_total,
        clasdis_written_in_grid,
        clasdis_written_out_grid,
        aaogen_written_total,
        aaogen_written_in_grid,
        aaogen_skipped_out_grid,
        aaogen_skipped_quota_met
    )

    draw_mix_debug_mx2(h_written_dis, h_written_aao, h_written_mix, OUTPUT_MIX_DEBUG_MX2_PNG)

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

    print(f"  wrote debug report: {OUTPUT_MIX_DEBUG_TXT}")
    print(f"  wrote debug Mx2 png: {OUTPUT_MIX_DEBUG_MX2_PNG}")
#enddef


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("--data", required=True, help="Path to data ROOT file")
    ap.add_argument("--aaogen", required=True, help="Path to aaogen ROOT file")
    ap.add_argument("--clasdis", required=True, help="Path to clasdis ROOT file")
    ap.add_argument("--max_events", type=int, default=-1, help="Max events per file (-1 means all events)")
    ap.add_argument("--out", default=DEFAULT_MIXED_MC_ROOT, help="Output ROOT file for mixed MC")
    ap.add_argument("--force", action="store_true", help="Optional diagnostic: force specific w[row][col] values after fit")

    args = ap.parse_args()

    require_file(args.data)
    require_file(args.aaogen)
    require_file(args.clasdis)

    out_path = str(args.out)
    if out_path.strip() == "":
        fatal("--out may not be empty.")
    #endif

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)

    ensure_outdir(OUTPUT_YIELDS_PNG)
    ensure_outdir(OUTPUT_MIX_PNG)
    ensure_outdir(OUTPUT_DATAONLY_PNG)
    ensure_outdir(OUTPUT_WEIGHTS_TXT)
    ensure_outdir(OUTPUT_MIX_DEBUG_TXT)
    ensure_outdir(OUTPUT_MIX_DEBUG_MX2_PNG)
    ensure_outdir(out_path)

    f_data, t_data = open_tree(args.data, TREE_NAME)
    f_aao, t_aao = open_tree(args.aaogen, TREE_NAME)
    f_dis, t_dis = open_tree(args.clasdis, TREE_NAME)

    needed_data = ["runnum"] + list(PHASE1_BRANCHES_BASE)
    needed_mc_like = list(PHASE1_BRANCHES_BASE)

    require_branches(t_data, needed_data, "data")
    require_branches(t_aao, needed_mc_like, "aaogen")
    require_branches(t_dis, needed_mc_like, "clasdis")

    nrows = len(XB_EDGES) - 1
    ncols = len(TNEG_EDGES) - 1

    h_data = make_hist_grid("h_data", nrows, ncols)
    h_aao = make_hist_grid("h_aaogen", nrows, ncols)
    h_dis = make_hist_grid("h_clasdis", nrows, ncols)

    c_data = make_count_grid(nrows, ncols)
    c_aao = make_count_grid(nrows, ncols)
    c_dis = make_count_grid(nrows, ncols)

    max_events = parse_max_events(args.max_events)

    fill_all_bins_single_pass(t_data, h_data, c_data, max_events, True, MC_EB_FIXED)
    fill_all_bins_single_pass(t_aao, h_aao, c_aao, max_events, False, MC_EB_FIXED)
    fill_all_bins_single_pass(t_dis, h_dis, c_dis, max_events, False, MC_EB_FIXED)

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
    col_mix = ROOT.kGreen + 2

    style_hist(h_data_int, col_data, 2, 1)
    style_hist(h_aao_int, col_aao, 2, 2)
    style_hist(h_dis_int, col_dis, 2, 3)

    for r in range(nrows):
        for c in range(ncols):
            normalize_unit_area(h_data[r][c])
            normalize_unit_area(h_aao[r][c])
            normalize_unit_area(h_dis[r][c])

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

        print("FORCE mode enabled: overriding selected w[row][col] after fit:")

        for key, wv in forced_map.items():
            rr = key[0]
            cc = key[1]

            old = float(w_grid[rr][cc])
            w_grid[rr][cc] = float(wv)

            h_mix[rr][cc].Reset("ICESM")
            h_mix[rr][cc].Add(h_aao[rr][cc], float(w_grid[rr][cc]))
            h_mix[rr][cc].Add(h_dis[rr][cc], float(1.0 - w_grid[rr][cc]))

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
    print(f"  fixed MC-like beam energy for aaogen/clasdis histogram filling = {MC_EB_FIXED:.3f} GeV")
    print(f"  total chi2 = {total_chi2:.6e}")
    print(f"  wrote weights report: {OUTPUT_WEIGHTS_TXT}")

    write_weights_report(w_grid, wun_grid, chi2_grid, c_data, c_aao, c_dis, OUTPUT_WEIGHTS_TXT, forced_map)
    draw_canvas_mix(h_data, h_mix, c_data, w_grid, OUTPUT_MIX_PNG)

    print("Phase 3: writing mixed MC ROOT file")
    print("  aaogen  =", args.aaogen)
    print("  clasdis =", args.clasdis)
    print("  out     =", out_path)
    print(f"  recompute t,tmin,tprime and mc_t,mc_tmin,mc_tprime with Eb={MC_EB_FIXED:.3f} (ignore any existing t/tmin branches)")

    mix_mc_to_output_root(args.aaogen, args.clasdis, out_path, w_grid, max_events)

    print("  wrote mixed MC:", out_path)

    f_data.Close()
    f_aao.Close()
    f_dis.Close()
#enddef


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        sys.stderr.write("FATAL: " + str(exc) + "\n")
        sys.exit(1)
    #endif
#endif