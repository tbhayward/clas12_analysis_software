#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
normalization.py

Purpose
-------
Combine aaogen and clasdis for the epi+X missing-mass study in a way that is
much easier to tune manually and much more stable when aaogen has far fewer
events than clasdis.

Main ideas
----------
1) Fill per-pad Mx2 histograms in (xB, -tprime).
2) Use fine histograms for plotting and coarse histograms for fitting.
3) Determine an aaogen fraction w[r,c] per pad from a peak-focused objective,
   with an additional sideband penalty.
4) Allow manual pad-by-pad overrides and simple row suppression.
5) Dump all histogram contents to CSV so the mixture can be tuned manually.
6) Always build a mixed MC ROOT file using aaogen and clasdis as inputs.
7) Recompute and write t, tmin, tprime, mc_t, mc_tmin, mc_tprime.

Mixture definition
------------------
For each pad:
    H_mix = w * H_aaogen + (1 - w) * H_clasdis

Automatic fit strategy
----------------------
The automatic fit no longer tries to match the entire Mx2 range with a single
number.

Instead, it minimizes:
    objective(w) = chi2_peak(w) + lambda_sideband * penalty_sideband(w)

where:
  - chi2_peak compares the normalized DATA shape to the normalized mixture
    in the peak-sensitive region only.
  - penalty_sideband forces the mixture to remain consistent with the DATA
    sideband fraction.

This is solved with a simple scan in w from 0 to 1.

Manual weights
--------------
You can provide a manual weights file with lines of the form:
    row col weight
Example:
    0 0 0.20
    0 1 0.08
    0 2 0.05

Anything not listed keeps its automatic value.

Recommended first-pass usage
----------------------------
python normalization.py \
  --data /path/to/data.root \
  --aaogen /path/to/aaogen.root \
  --clasdis /path/to/clasdis.root \
  --max_events_clasdis 1000000 \
  --zero_aaogen_rows 1,2,3

Outputs
-------
output/yields.png
    Fine-bin per-pad DATA vs aaogen vs clasdis

output/yields_data_only.png
    Fine-bin integrated DATA vs aaogen vs clasdis

output/yields_mix.png
    Fine-bin per-pad DATA vs mixed template

output/yields_mix_integrated.png
    Fine-bin integrated DATA vs mixed template

output/weights.txt
    Per-pad weights and diagnostics

output/hist_shapes.csv
    Bin-by-bin raw and normalized histogram contents for all pads

output/mixed_mc.root
    Event-level mixed MC

output/mix_debug_report.txt
    Debug counters and per-pad written fractions

output/mix_debug_mx2.png
    Integrated Mx2 shapes of written clasdis, written aaogen, and total mixture
"""

import os
import sys
import math
import csv
import argparse
from array import array

import ROOT

TREE_NAME = "PhysicsEvents"

XB_EDGES = [0.10, 0.25, 0.35, 0.45, 0.60]
TNEG_EDGES = [0.05, 0.25, 0.45, 0.65, 0.85, 1.05, 1.25]

# Plot histogram configuration
MX2_PLOT_MIN = 0.3
MX2_PLOT_MAX = 2.0
MX2_PLOT_NBINS_DEFAULT = 100

# Fit histogram configuration
MX2_FIT_MIN_GLOBAL = 0.3
MX2_FIT_MAX_GLOBAL = 2.0
MX2_FIT_NBINS_DEFAULT = 40

# Default fit windows
PEAK_MIN_DEFAULT = 0.82
PEAK_MAX_DEFAULT = 1.02
SIDEBAND_MIN_DEFAULT = 1.15
SIDEBAND_MAX_DEFAULT = 2.00

# Objective scan resolution
W_SCAN_STEPS_DEFAULT = 1000

# Sideband penalty strength
LAMBDA_SIDEBAND_DEFAULT = 1.0

# Default output paths
OUTPUT_YIELDS_PNG = "output/yields.png"
OUTPUT_DATAONLY_PNG = "output/yields_data_only.png"
OUTPUT_MIX_PNG = "output/yields_mix.png"
OUTPUT_MIX_INTEGRATED_PNG = "output/yields_mix_integrated.png"
OUTPUT_WEIGHTS_TXT = "output/weights.txt"
OUTPUT_HIST_CSV = "output/hist_shapes.csv"

DEFAULT_MIXED_MC_ROOT = "output/mixed_mc.root"
OUTPUT_MIX_DEBUG_TXT = "output/mix_debug_report.txt"
OUTPUT_MIX_DEBUG_MX2_PNG = "output/mix_debug_mx2.png"

# Fixed beam energy for MC mixing stage
MC_EB_FIXED = 10.556

# Default event limits
DEFAULT_MAX_EVENTS_DATA = -1
DEFAULT_MAX_EVENTS_AAOGEN = -1
DEFAULT_MAX_EVENTS_CLASDIS = 1000000
DEFAULT_MAX_EVENTS_MC_AAOGEN = -1
DEFAULT_MAX_EVENTS_MC_CLASDIS = -1

# Masses (GeV)
MASS_E = 0.000511
MASS_PI = 0.139570
MASS_N = 0.9382720813


def fatal(msg):
    raise RuntimeError(msg)
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


def require_file(path):
    if path is None or str(path).strip() == "":
        fatal("Missing required input path.")
    #endif
    if not os.path.isfile(path):
        fatal("File not found: " + str(path))
    #endif
#enddef


def ensure_outdir(path):
    dname = os.path.dirname(str(path))
    if dname != "":
        os.makedirs(dname, exist_ok=True)
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


def make_hist_grid(prefix, nrows, ncols, nbins, xmin, xmax):
    grid = []
    for r in range(nrows):
        row = []
        for c in range(ncols):
            hname = f"{prefix}_r{r}_c{c}"
            hist = ROOT.TH1F(hname, "", nbins, xmin, xmax)
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


def clone_grid(src_grid, prefix):
    out = []
    for r in range(len(src_grid)):
        row = []
        for c in range(len(src_grid[r])):
            hist = src_grid[r][c].Clone(f"{prefix}_r{r}_c{c}")
            row.append(hist)
        #endfor
        out.append(row)
    #endfor
    return out
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


def sigma_from_raw_counts(n_i, N_pad):
    if N_pad <= 0.0:
        return 1.0
    #endif
    if n_i <= 0.0:
        return 1.0 / float(N_pad)
    #endif
    return math.sqrt(n_i) / float(N_pad)
#enddef


def raw_count_in_window(hist_raw, xmin, xmax):
    total = 0.0
    for ibin in range(1, hist_raw.GetNbinsX() + 1):
        xcen = hist_raw.GetXaxis().GetBinCenter(ibin)
        if xcen >= xmin and xcen <= xmax:
            total += hist_raw.GetBinContent(ibin)
        #endif
    #endfor
    return total
#enddef


def normalized_fraction_in_window(hist_norm, xmin, xmax):
    total = 0.0
    for ibin in range(1, hist_norm.GetNbinsX() + 1):
        xcen = hist_norm.GetXaxis().GetBinCenter(ibin)
        if xcen >= xmin and xcen <= xmax:
            total += hist_norm.GetBinContent(ibin)
        #endif
    #endfor
    return total
#enddef


def compute_peak_chi2_for_w(hd_norm, ha_norm, hc_norm, hd_raw, w, peak_min, peak_max):
    nbins = hd_norm.GetNbinsX()
    N_pad = hd_raw.Integral(1, nbins)

    if N_pad <= 0.0:
        return 0.0
    #endif

    chi2 = 0.0
    used_bins = 0

    for ibin in range(1, nbins + 1):
        xcen = hd_norm.GetXaxis().GetBinCenter(ibin)
        if xcen < peak_min or xcen > peak_max:
            continue
        #endif

        D = hd_norm.GetBinContent(ibin)
        A = ha_norm.GetBinContent(ibin)
        C = hc_norm.GetBinContent(ibin)
        M = w * A + (1.0 - w) * C

        n_i = hd_raw.GetBinContent(ibin)
        sig = sigma_from_raw_counts(n_i, N_pad)
        if sig <= 0.0:
            continue
        #endif

        diff = D - M
        chi2 += (diff * diff) / (sig * sig)
        used_bins += 1
    #endfor

    if used_bins <= 0:
        return 0.0
    #endif

    return chi2
#enddef


def compute_sideband_penalty_for_w(hd_norm, ha_norm, hc_norm, hd_raw,
                                   w, side_min, side_max, lambda_sideband):
    nbins = hd_raw.GetNbinsX()
    N_pad = hd_raw.Integral(1, nbins)

    if N_pad <= 0.0:
        return 0.0, 0.0, 0.0, 1.0
    #endif

    S_data = normalized_fraction_in_window(hd_norm, side_min, side_max)
    S_aao = normalized_fraction_in_window(ha_norm, side_min, side_max)
    S_dis = normalized_fraction_in_window(hc_norm, side_min, side_max)
    S_mix = w * S_aao + (1.0 - w) * S_dis

    n_side = raw_count_in_window(hd_raw, side_min, side_max)
    sigma_side = sigma_from_raw_counts(n_side, N_pad)

    if sigma_side <= 0.0:
        sigma_side = 1.0 / float(max(1.0, N_pad))
    #endif

    pull = (S_mix - S_data) / sigma_side
    penalty = lambda_sideband * pull * pull

    return penalty, S_data, S_mix, sigma_side
#enddef


def compute_objective_for_w(hd_norm, ha_norm, hc_norm, hd_raw,
                            w, peak_min, peak_max,
                            side_min, side_max, lambda_sideband):
    chi2_peak = compute_peak_chi2_for_w(hd_norm, ha_norm, hc_norm, hd_raw,
                                        w, peak_min, peak_max)

    penalty_side, S_data, S_mix, sigma_side = compute_sideband_penalty_for_w(
        hd_norm, ha_norm, hc_norm, hd_raw,
        w, side_min, side_max, lambda_sideband
    )

    objective = chi2_peak + penalty_side
    return objective, chi2_peak, penalty_side, S_data, S_mix, sigma_side
#enddef


def scan_best_w_for_pad(hd_norm, ha_norm, hc_norm, hd_raw,
                        peak_min, peak_max,
                        side_min, side_max,
                        lambda_sideband, nsteps):
    best_w = 0.0
    best_obj = None
    best_peak_chi2 = 0.0
    best_side_penalty = 0.0
    best_S_data = 0.0
    best_S_mix = 0.0
    best_sigma_side = 1.0

    hist_ok = (
        hd_norm.Integral(1, hd_norm.GetNbinsX()) > 0.0 and
        ha_norm.Integral(1, ha_norm.GetNbinsX()) > 0.0 and
        hc_norm.Integral(1, hc_norm.GetNbinsX()) > 0.0 and
        hd_raw.Integral(1, hd_raw.GetNbinsX()) > 0.0
    )

    if not hist_ok:
        return {
            "w": 0.0,
            "w_unclipped": 0.0,
            "objective": 0.0,
            "peak_chi2": 0.0,
            "sideband_penalty": 0.0,
            "sideband_data_fraction": 0.0,
            "sideband_mix_fraction": 0.0,
            "sideband_sigma": 1.0
        }
    #endif

    for istep in range(nsteps + 1):
        w = float(istep) / float(nsteps)

        obj, chi2_peak, penalty_side, S_data, S_mix, sigma_side = compute_objective_for_w(
            hd_norm, ha_norm, hc_norm, hd_raw,
            w, peak_min, peak_max,
            side_min, side_max, lambda_sideband
        )

        if best_obj is None or obj < best_obj:
            best_obj = obj
            best_w = w
            best_peak_chi2 = chi2_peak
            best_side_penalty = penalty_side
            best_S_data = S_data
            best_S_mix = S_mix
            best_sigma_side = sigma_side
        #endif
    #endfor

    return {
        "w": best_w,
        "w_unclipped": best_w,
        "objective": best_obj,
        "peak_chi2": best_peak_chi2,
        "sideband_penalty": best_side_penalty,
        "sideband_data_fraction": best_S_data,
        "sideband_mix_fraction": best_S_mix,
        "sideband_sigma": best_sigma_side
    }
#enddef


def compute_w_grid_and_mix(h_data_fit_norm, h_aao_fit_norm, h_dis_fit_norm, h_data_fit_raw,
                           h_aao_plot_norm, h_dis_plot_norm,
                           peak_min, peak_max,
                           side_min, side_max,
                           lambda_sideband, nsteps):
    nrows = len(h_data_fit_norm)
    ncols = len(h_data_fit_norm[0])

    w_grid = [[0.0 for _ in range(ncols)] for _ in range(nrows)]
    wun_grid = [[0.0 for _ in range(ncols)] for _ in range(nrows)]
    objective_grid = [[0.0 for _ in range(ncols)] for _ in range(nrows)]
    peak_chi2_grid = [[0.0 for _ in range(ncols)] for _ in range(nrows)]
    sideband_penalty_grid = [[0.0 for _ in range(ncols)] for _ in range(nrows)]
    sideband_data_fraction_grid = [[0.0 for _ in range(ncols)] for _ in range(nrows)]
    sideband_mix_fraction_grid = [[0.0 for _ in range(ncols)] for _ in range(nrows)]
    sideband_sigma_grid = [[1.0 for _ in range(ncols)] for _ in range(nrows)]

    h_mix_plot = []
    for r in range(nrows):
        row = []
        for c in range(ncols):
            result = scan_best_w_for_pad(
                h_data_fit_norm[r][c],
                h_aao_fit_norm[r][c],
                h_dis_fit_norm[r][c],
                h_data_fit_raw[r][c],
                peak_min,
                peak_max,
                side_min,
                side_max,
                lambda_sideband,
                nsteps
            )

            w_grid[r][c] = result["w"]
            wun_grid[r][c] = result["w_unclipped"]
            objective_grid[r][c] = result["objective"]
            peak_chi2_grid[r][c] = result["peak_chi2"]
            sideband_penalty_grid[r][c] = result["sideband_penalty"]
            sideband_data_fraction_grid[r][c] = result["sideband_data_fraction"]
            sideband_mix_fraction_grid[r][c] = result["sideband_mix_fraction"]
            sideband_sigma_grid[r][c] = result["sideband_sigma"]

            hmix = h_aao_plot_norm[r][c].Clone(f"h_mix_plot_r{r}_c{c}")
            hmix.Reset("ICESM")
            hmix.Add(h_aao_plot_norm[r][c], float(w_grid[r][c]))
            hmix.Add(h_dis_plot_norm[r][c], float(1.0 - w_grid[r][c]))
            row.append(hmix)
        #endfor
        h_mix_plot.append(row)
    #endfor

    return {
        "w_grid": w_grid,
        "wun_grid": wun_grid,
        "objective_grid": objective_grid,
        "peak_chi2_grid": peak_chi2_grid,
        "sideband_penalty_grid": sideband_penalty_grid,
        "sideband_data_fraction_grid": sideband_data_fraction_grid,
        "sideband_mix_fraction_grid": sideband_mix_fraction_grid,
        "sideband_sigma_grid": sideband_sigma_grid,
        "h_mix_plot": h_mix_plot
    }
#enddef


def compute_objective_grid_for_existing_weights(w_grid,
                                                h_data_fit_norm, h_aao_fit_norm, h_dis_fit_norm, h_data_fit_raw,
                                                peak_min, peak_max,
                                                side_min, side_max,
                                                lambda_sideband):
    nrows = len(w_grid)
    ncols = len(w_grid[0])

    objective_grid = [[0.0 for _ in range(ncols)] for _ in range(nrows)]
    peak_chi2_grid = [[0.0 for _ in range(ncols)] for _ in range(nrows)]
    sideband_penalty_grid = [[0.0 for _ in range(ncols)] for _ in range(nrows)]
    sideband_data_fraction_grid = [[0.0 for _ in range(ncols)] for _ in range(nrows)]
    sideband_mix_fraction_grid = [[0.0 for _ in range(ncols)] for _ in range(nrows)]
    sideband_sigma_grid = [[1.0 for _ in range(ncols)] for _ in range(nrows)]

    for r in range(nrows):
        for c in range(ncols):
            obj, chi2_peak, penalty_side, S_data, S_mix, sigma_side = compute_objective_for_w(
                h_data_fit_norm[r][c],
                h_aao_fit_norm[r][c],
                h_dis_fit_norm[r][c],
                h_data_fit_raw[r][c],
                float(w_grid[r][c]),
                peak_min,
                peak_max,
                side_min,
                side_max,
                lambda_sideband
            )

            objective_grid[r][c] = obj
            peak_chi2_grid[r][c] = chi2_peak
            sideband_penalty_grid[r][c] = penalty_side
            sideband_data_fraction_grid[r][c] = S_data
            sideband_mix_fraction_grid[r][c] = S_mix
            sideband_sigma_grid[r][c] = sigma_side
        #endfor
    #endfor

    return {
        "objective_grid": objective_grid,
        "peak_chi2_grid": peak_chi2_grid,
        "sideband_penalty_grid": sideband_penalty_grid,
        "sideband_data_fraction_grid": sideband_data_fraction_grid,
        "sideband_mix_fraction_grid": sideband_mix_fraction_grid,
        "sideband_sigma_grid": sideband_sigma_grid
    }
#enddef


def rebuild_mix_plot_grid(w_grid, h_aao_plot_norm, h_dis_plot_norm):
    nrows = len(w_grid)
    ncols = len(w_grid[0])
    out = []

    for r in range(nrows):
        row = []
        for c in range(ncols):
            hmix = h_aao_plot_norm[r][c].Clone(f"h_mix_plot_manual_r{r}_c{c}")
            hmix.Reset("ICESM")
            hmix.Add(h_aao_plot_norm[r][c], float(w_grid[r][c]))
            hmix.Add(h_dis_plot_norm[r][c], float(1.0 - w_grid[r][c]))
            row.append(hmix)
        #endfor
        out.append(row)
    #endfor

    return out
#enddef


def draw_canvas_integrated_threeway(hd_int, ha_int, hc_int, Nd, Na, Nc, outpng):
    canv = ROOT.TCanvas("c_integrated_threeway", "Integrated Mx2", 1200, 900)
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


def draw_canvas_integrated_mix(hd_int, hm_int, Nd, outpng):
    canv = ROOT.TCanvas("c_integrated_mix", "Integrated Mx2 data vs mix", 1200, 900)
    pad = canv.cd(1)
    pad.SetGrid(1, 1)
    pad.SetLeftMargin(0.18)
    pad.SetRightMargin(0.05)
    pad.SetBottomMargin(0.14)
    pad.SetTopMargin(0.08)

    ymax = max(hd_int.GetMaximum(), hm_int.GetMaximum())
    if ymax <= 0.0:
        ymax = 1.0
    #endif
    ymax *= 1.2

    set_axes_and_range(hd_int, ymax)
    hd_int.Draw("hist")
    hm_int.Draw("hist same")

    leg = ROOT.TLegend(0.55, 0.76, 0.94, 0.92)
    leg.SetBorderSize(1)
    leg.SetFillStyle(1001)
    leg.SetFillColor(ROOT.kWhite)
    leg.SetTextSize(0.042)
    leg.AddEntry(hd_int, f"data (N={int(Nd)})", "l")
    leg.AddEntry(hm_int, "mixed template", "l")
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

    canv = ROOT.TCanvas("c_yields", "Mx2 by pad", 2400, 1400)
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


def draw_canvas_mix(h_data, h_mix, c_data, w_grid, outpng,
                    peak_min, peak_max, side_min, side_max):
    nrows = len(h_data)
    ncols = len(h_data[0])

    canv = ROOT.TCanvas("c_mix", "Mx2 data vs mix", 2400, 1400)
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

            leg = ROOT.TLegend(0.43, 0.63, 0.94, 0.92)
            leg.SetBorderSize(1)
            leg.SetFillStyle(1001)
            leg.SetFillColor(ROOT.kWhite)
            leg.SetTextSize(0.036)
            leg.AddEntry(hd, f"data (N={int(c_data[r][c])})", "l")
            leg.AddEntry(hm, f"mix (aaogen frac w={w:.4f})", "l")
            leg.AddEntry("", f"peak: [{peak_min:.2f}, {peak_max:.2f}]", "")
            leg.AddEntry("", f"sideband: [{side_min:.2f}, {side_max:.2f}]", "")
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


def write_weights_report(path, w_grid, wun_grid,
                         objective_grid, peak_chi2_grid, sideband_penalty_grid,
                         sideband_data_fraction_grid, sideband_mix_fraction_grid, sideband_sigma_grid,
                         c_data, c_aao, c_dis,
                         peak_min, peak_max, side_min, side_max, lambda_sideband,
                         forced_map=None, manual_override_map=None):
    ensure_outdir(path)

    total_obj = 0.0
    total_peak = 0.0
    total_side = 0.0

    for r in range(len(w_grid)):
        for c in range(len(w_grid[r])):
            total_obj += objective_grid[r][c]
            total_peak += peak_chi2_grid[r][c]
            total_side += sideband_penalty_grid[r][c]
        #endfor
    #endfor

    with open(path, "w") as fout:
        fout.write("Per-bin aaogen mixture weights\n")
        fout.write("Definition: H_mix = w * H_aaogen + (1-w) * H_clasdis\n")
        fout.write(f"Peak window used in fit objective: [{peak_min:.3f}, {peak_max:.3f}]\n")
        fout.write(f"Sideband window used in fit objective: [{side_min:.3f}, {side_max:.3f}]\n")
        fout.write(f"Sideband penalty strength lambda_sideband = {lambda_sideband:.6f}\n")
        fout.write("Automatic fit method: scan in w from 0 to 1\n")
        fout.write("\n")
        fout.write(f"Total objective = {total_obj:.6e}\n")
        fout.write(f"Total peak chi2 = {total_peak:.6e}\n")
        fout.write(f"Total sideband penalty = {total_side:.6e}\n")
        fout.write("\n")

        if forced_map is not None and len(forced_map) > 0:
            fout.write("NOTE: Some pads were changed by the built-in force map.\n")
        #endif
        if manual_override_map is not None and len(manual_override_map) > 0:
            fout.write("NOTE: Some pads were changed by the manual weights file.\n")
        #endif
        fout.write("\n")

        for r in range(len(w_grid)):
            xb_lo = XB_EDGES[r]
            xb_hi = XB_EDGES[r + 1]
            fout.write(f"Row {r}: xB [{xb_lo:.2f}, {xb_hi:.2f})\n")

            for c in range(len(w_grid[r])):
                t_lo = TNEG_EDGES[c]
                t_hi = TNEG_EDGES[c + 1]

                tag_bits = []
                if forced_map is not None and (r, c) in forced_map:
                    tag_bits.append("FORCED")
                #endif
                if manual_override_map is not None and (r, c) in manual_override_map:
                    tag_bits.append("MANUAL")
                #endif

                tag = ""
                if len(tag_bits) > 0:
                    tag = "  " + ",".join(tag_bits)
                #endif

                fout.write(
                    f"  Col {c}: -tprime [{t_lo:.2f}, {t_hi:.2f})  "
                    f"w={w_grid[r][c]:.6f}  "
                    f"w_scan={wun_grid[r][c]:.6f}  "
                    f"objective={objective_grid[r][c]:.6e}  "
                    f"peak_chi2={peak_chi2_grid[r][c]:.6e}  "
                    f"side_penalty={sideband_penalty_grid[r][c]:.6e}  "
                    f"Sdata={sideband_data_fraction_grid[r][c]:.6f}  "
                    f"Smix={sideband_mix_fraction_grid[r][c]:.6f}  "
                    f"Ssigma={sideband_sigma_grid[r][c]:.6e}  "
                    f"N(data,aao,dis)=({int(c_data[r][c])},{int(c_aao[r][c])},{int(c_dis[r][c])})"
                    f"{tag}\n"
                )
            #endfor
            fout.write("\n")
        #endfor
#enddef


def dump_hist_shapes_csv(path,
                         h_data_raw, h_aao_raw, h_dis_raw,
                         h_data_norm, h_aao_norm, h_dis_norm,
                         h_mix_norm, w_grid):
    ensure_outdir(path)

    with open(path, "w", newline="") as fout:
        writer = csv.writer(fout)
        writer.writerow([
            "row",
            "col",
            "xB_lo",
            "xB_hi",
            "tneg_lo",
            "tneg_hi",
            "mx2_bin",
            "mx2_low",
            "mx2_high",
            "mx2_center",
            "w",
            "data_raw",
            "aaogen_raw",
            "clasdis_raw",
            "data_norm",
            "aaogen_norm",
            "clasdis_norm",
            "mix_norm"
        ])

        for r in range(len(h_data_raw)):
            xb_lo = XB_EDGES[r]
            xb_hi = XB_EDGES[r + 1]

            for c in range(len(h_data_raw[r])):
                t_lo = TNEG_EDGES[c]
                t_hi = TNEG_EDGES[c + 1]

                hd_raw = h_data_raw[r][c]
                ha_raw = h_aao_raw[r][c]
                hc_raw = h_dis_raw[r][c]

                hd_norm = h_data_norm[r][c]
                ha_norm = h_aao_norm[r][c]
                hc_norm = h_dis_norm[r][c]
                hm_norm = h_mix_norm[r][c]

                nbins = hd_raw.GetNbinsX()

                for ibin in range(1, nbins + 1):
                    xlow = hd_raw.GetXaxis().GetBinLowEdge(ibin)
                    xhigh = hd_raw.GetXaxis().GetBinUpEdge(ibin)
                    xcen = hd_raw.GetXaxis().GetBinCenter(ibin)

                    writer.writerow([
                        r,
                        c,
                        xb_lo,
                        xb_hi,
                        t_lo,
                        t_hi,
                        ibin,
                        xlow,
                        xhigh,
                        xcen,
                        w_grid[r][c],
                        hd_raw.GetBinContent(ibin),
                        ha_raw.GetBinContent(ibin),
                        hc_raw.GetBinContent(ibin),
                        hd_norm.GetBinContent(ibin),
                        ha_norm.GetBinContent(ibin),
                        hc_norm.GetBinContent(ibin),
                        hm_norm.GetBinContent(ibin)
                    ])
                #endfor
            #endfor
        #endfor
#enddef


def parse_int_list(text):
    if text is None:
        return []
    #endif
    s = str(text).strip()
    if s == "":
        return []
    #endif

    out = []
    parts = s.split(",")
    for part in parts:
        p = part.strip()
        if p == "":
            continue
        #endif
        out.append(int(p))
    #endfor
    return out
#enddef


def parse_manual_weights(path):
    require_file(path)

    out = {}
    with open(path, "r") as fin:
        for lineno, line in enumerate(fin, start=1):
            s = line.strip()

            if s == "":
                continue
            #endif
            if s.startswith("#"):
                continue
            #endif

            tokens = s.replace(",", " ").split()
            if len(tokens) != 3:
                fatal(f"Malformed manual weights file line {lineno}: expected 'row col weight'")
            #endif

            row = int(tokens[0])
            col = int(tokens[1])
            weight = float(tokens[2])

            if weight < 0.0 or weight > 1.0:
                fatal(f"Manual weight out of range on line {lineno}: {weight}")
            #endif

            out[(row, col)] = weight
        #endfor
    #endwith

    return out
#enddef


def apply_zero_aaogen_rows(w_grid, row_list):
    for row in row_list:
        if row < 0 or row >= len(w_grid):
            fatal(f"Requested zero_aaogen_rows contains invalid row index: {row}")
        #endif
        for c in range(len(w_grid[row])):
            w_grid[row][c] = 0.0
        #endfor
    #endfor
#enddef


def apply_manual_weights_map(w_grid, manual_map):
    for key, val in manual_map.items():
        r = key[0]
        c = key[1]

        if r < 0 or r >= len(w_grid):
            fatal(f"Manual weight row out of range: {r}")
        #endif
        if c < 0 or c >= len(w_grid[r]):
            fatal(f"Manual weight col out of range: {c}")
        #endif

        w_grid[r][c] = float(val)
    #endfor
#enddef


def apply_force_map(w_grid):
    forced_map = {
        (0, 2): 0.05,
        (0, 3): 0.08,
        (0, 4): 0.03,
        (0, 5): 0.03
    }

    for key, val in forced_map.items():
        r = key[0]
        c = key[1]
        if r < len(w_grid) and c < len(w_grid[r]):
            w_grid[r][c] = float(val)
        #endif
    #endfor

    return forced_map
#enddef


def build_branch_buffers_phase1():
    return {
        "runnum": array("i", [0]),
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
#enddef


def bind_phase1_tree(tree, buffers):
    needed = ["runnum", "e_p", "e_theta", "e_phi", "p_p", "p_theta", "p_phi", "x", "Q2", "Mx2"]

    tree.SetBranchStatus("*", 0)
    for bname in needed:
        tree.SetBranchStatus(bname, 1)
    #endfor

    tree.SetBranchAddress("runnum", buffers["runnum"])
    tree.SetBranchAddress("e_p", buffers["e_p"])
    tree.SetBranchAddress("e_theta", buffers["e_theta"])
    tree.SetBranchAddress("e_phi", buffers["e_phi"])
    tree.SetBranchAddress("p_p", buffers["p_p"])
    tree.SetBranchAddress("p_theta", buffers["p_theta"])
    tree.SetBranchAddress("p_phi", buffers["p_phi"])
    tree.SetBranchAddress("x", buffers["x"])
    tree.SetBranchAddress("Q2", buffers["Q2"])
    tree.SetBranchAddress("Mx2", buffers["Mx2"])
#enddef


def fill_all_bins_single_pass(tree,
                              h_plot_grid, h_fit_grid, cgrid,
                              max_events):
    buffers = build_branch_buffers_phase1()
    bind_phase1_tree(tree, buffers)

    n_entries = int(tree.GetEntries())
    if max_events is None:
        n_to_process = n_entries
    else:
        n_to_process = min(n_entries, int(max_events))
    #endif

    for ientry in range(n_to_process):
        tree.GetEntry(ientry)

        xb_val = float(buffers["x"][0])
        rb = find_bin(xb_val, XB_EDGES)
        if rb < 0:
            continue
        #endif

        Q2_val = float(buffers["Q2"][0])

        t_val = compute_t_scalar(
            int(buffers["runnum"][0]),
            float(buffers["e_p"][0]),
            float(buffers["e_theta"][0]),
            float(buffers["e_phi"][0]),
            float(buffers["p_p"][0]),
            float(buffers["p_theta"][0]),
            float(buffers["p_phi"][0])
        )

        tmin_val = compute_tmin_exact(xb_val, Q2_val)
        tprime_val = t_val - tmin_val
        tneg_val = -tprime_val

        cb = find_bin(tneg_val, TNEG_EDGES)
        if cb < 0:
            continue
        #endif

        mx2_val = float(buffers["Mx2"][0])

        h_plot_grid[rb][cb].Fill(mx2_val)
        h_fit_grid[rb][cb].Fill(mx2_val)
        cgrid[rb][cb] += 1
    #endfor

    tree.SetBranchStatus("*", 1)
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


def compute_mc_t_quantities_from_gen(mc_x_val, mc_Q2_val,
                                     mc_e_p, mc_e_theta, mc_e_phi,
                                     mc_p_p, mc_p_theta, mc_p_phi):
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


def choose_phase3_input_paths(args):
    mc_aao_path = args.aaogen
    mc_dis_path = args.clasdis

    if args.mc_aaogen is not None and str(args.mc_aaogen).strip() != "":
        mc_aao_path = args.mc_aaogen
    #endif
    if args.mc_clasdis is not None and str(args.mc_clasdis).strip() != "":
        mc_dis_path = args.mc_clasdis
    #endif

    require_file(mc_aao_path)
    require_file(mc_dis_path)

    return mc_aao_path, mc_dis_path
#enddef


def mix_mc_to_output_root(mc_aao_path, mc_dis_path, out_root_path, w_grid,
                          max_events_mc_aao, max_events_mc_dis):
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

    require_branches(t_aao, mc_needed, "phase3_aaogen")
    require_branches(t_dis, mc_needed, "phase3_clasdis")

    ensure_outdir(out_root_path)

    fout = ROOT.TFile.Open(out_root_path, "RECREATE")
    if not fout or fout.IsZombie():
        fatal("Failed to create output ROOT file: " + str(out_root_path))
    #endif

    tout = ROOT.TTree(TREE_NAME, "Mixed MC")

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
    Naao_target = [[0 for _ in range(ncols)] for _ in range(nrows)]
    Naao_written = [[0 for _ in range(ncols)] for _ in range(nrows)]

    clasdis_written_total = 0
    clasdis_written_in_grid = 0
    clasdis_written_out_grid = 0

    aaogen_written_total = 0
    aaogen_written_in_grid = 0
    aaogen_skipped_out_grid = 0
    aaogen_skipped_quota_met = 0

    h_written_dis = ROOT.TH1F("h_written_dis", "", MX2_PLOT_NBINS_DEFAULT, MX2_PLOT_MIN, MX2_PLOT_MAX)
    h_written_aao = ROOT.TH1F("h_written_aao", "", MX2_PLOT_NBINS_DEFAULT, MX2_PLOT_MIN, MX2_PLOT_MAX)
    h_written_mix = ROOT.TH1F("h_written_mix", "", MX2_PLOT_NBINS_DEFAULT, MX2_PLOT_MIN, MX2_PLOT_MAX)

    def bind_mc_tree(tree):
        tree.SetBranchStatus("*", 0)
        for bname in mc_needed:
            tree.SetBranchStatus(bname, 1)
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
    if max_events_mc_dis is None:
        n_to_process_dis = n_entries_dis
    else:
        n_to_process_dis = min(n_entries_dis, int(max_events_mc_dis))
    #endif

    for ientry in range(n_to_process_dis):
        t_dis.GetEntry(ientry)

        row, col, t_val, tmin_val, tp_val = compute_bin_indices_from_reco(
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
        h_written_dis.Fill(float(Mx2[0]))
        h_written_mix.Fill(float(Mx2[0]))

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

    bind_mc_tree(t_aao)

    n_entries_aao = int(t_aao.GetEntries())
    if max_events_mc_aao is None:
        n_to_process_aao = n_entries_aao
    else:
        n_to_process_aao = min(n_entries_aao, int(max_events_mc_aao))
    #endif

    for ientry in range(n_to_process_aao):
        t_aao.GetEntry(ientry)

        row, col, t_val, tmin_val, tp_val = compute_bin_indices_from_reco(
            float(x[0]),
            float(Q2[0]),
            float(e_p[0]),
            float(e_theta[0]),
            float(e_phi[0]),
            float(p_p[0]),
            float(p_theta[0]),
            float(p_phi[0])
        )

        if row < 0 or col < 0:
            aaogen_skipped_out_grid += 1
            continue
        #endif

        if Naao_written[row][col] >= Naao_target[row][col]:
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
        Naao_written[row][col] += 1
        aaogen_written_total += 1
        aaogen_written_in_grid += 1

        h_written_aao.Fill(float(Mx2[0]))
        h_written_mix.Fill(float(Mx2[0]))
    #endfor

    fout.cd()
    tout.Write()
    fout.Close()
    f_aao.Close()
    f_dis.Close()

    ensure_outdir(OUTPUT_MIX_DEBUG_TXT)
    with open(OUTPUT_MIX_DEBUG_TXT, "w") as fout_txt:
        fout_txt.write("Phase 3 debug report\n")
        fout_txt.write("All clasdis events were written; aaogen was topped up in-grid according to w_grid.\n")
        fout_txt.write("\n")
        fout_txt.write(f"clasdis written total       = {clasdis_written_total}\n")
        fout_txt.write(f"clasdis written in-grid     = {clasdis_written_in_grid}\n")
        fout_txt.write(f"clasdis written out-of-grid = {clasdis_written_out_grid}\n")
        fout_txt.write(f"aaogen written total        = {aaogen_written_total}\n")
        fout_txt.write(f"aaogen written in-grid      = {aaogen_written_in_grid}\n")
        fout_txt.write(f"aaogen skipped out-of-grid  = {aaogen_skipped_out_grid}\n")
        fout_txt.write(f"aaogen skipped quota met    = {aaogen_skipped_quota_met}\n")
        fout_txt.write("\n")

        denom = float(clasdis_written_total + aaogen_written_total)
        if denom > 0.0:
            fout_txt.write(f"global achieved aaogen fraction = {aaogen_written_total / denom:.6f}\n")
        else:
            fout_txt.write("global achieved aaogen fraction = 0.000000\n")
        #endif
        fout_txt.write("\n")

        for r in range(nrows):
            xb_lo = XB_EDGES[r]
            xb_hi = XB_EDGES[r + 1]
            fout_txt.write(f"Row {r}: xB [{xb_lo:.2f}, {xb_hi:.2f})\n")
            for c in range(ncols):
                t_lo = TNEG_EDGES[c]
                t_hi = TNEG_EDGES[c + 1]

                nd = int(Ndis[r][c])
                na_target = int(Naao_target[r][c])
                na_written = int(Naao_written[r][c])

                frac = 0.0
                if (nd + na_written) > 0:
                    frac = float(na_written) / float(nd + na_written)
                #endif

                fout_txt.write(
                    f"  Col {c}: -tprime [{t_lo:.2f}, {t_hi:.2f})  "
                    f"w_target={float(w_grid[r][c]):.6f}  "
                    f"Ndis={nd}  "
                    f"Naao_target={na_target}  "
                    f"Naao_written={na_written}  "
                    f"w_achieved={frac:.6f}\n"
                )
            #endfor
            fout_txt.write("\n")
        #endfor
    #endwith

    h_dis_norm = h_written_dis.Clone("h_written_dis_norm")
    h_aao_norm = h_written_aao.Clone("h_written_aao_norm")
    h_mix_norm = h_written_mix.Clone("h_written_mix_norm")

    normalize_unit_area(h_dis_norm)
    normalize_unit_area(h_aao_norm)
    normalize_unit_area(h_mix_norm)

    style_hist(h_dis_norm, ROOT.kBlue, 2, 1)
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

    canv.SaveAs(OUTPUT_MIX_DEBUG_MX2_PNG)
#enddef


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("--data", required=True, help="Path to data ROOT file")
    ap.add_argument("--aaogen", required=True, help="Path to aaogen ROOT file")
    ap.add_argument("--clasdis", required=True, help="Path to clasdis ROOT file")

    ap.add_argument("--plot_nbins", type=int, default=MX2_PLOT_NBINS_DEFAULT,
                    help="Number of bins for plotting histograms")
    ap.add_argument("--fit_nbins", type=int, default=MX2_FIT_NBINS_DEFAULT,
                    help="Number of bins for fitting histograms")

    ap.add_argument("--peak_min", type=float, default=PEAK_MIN_DEFAULT,
                    help="Lower edge of peak window for automatic fit")
    ap.add_argument("--peak_max", type=float, default=PEAK_MAX_DEFAULT,
                    help="Upper edge of peak window for automatic fit")
    ap.add_argument("--sideband_min", type=float, default=SIDEBAND_MIN_DEFAULT,
                    help="Lower edge of sideband window for automatic fit")
    ap.add_argument("--sideband_max", type=float, default=SIDEBAND_MAX_DEFAULT,
                    help="Upper edge of sideband window for automatic fit")
    ap.add_argument("--lambda_sideband", type=float, default=LAMBDA_SIDEBAND_DEFAULT,
                    help="Strength of sideband penalty in objective")
    ap.add_argument("--w_scan_steps", type=int, default=W_SCAN_STEPS_DEFAULT,
                    help="Number of scan intervals in w from 0 to 1")

    ap.add_argument("--max_events_data", type=int, default=DEFAULT_MAX_EVENTS_DATA,
                    help="Max data events for phase 1 (-1 means all)")
    ap.add_argument("--max_events_aaogen", type=int, default=DEFAULT_MAX_EVENTS_AAOGEN,
                    help="Max aaogen events for phase 1 (-1 means all)")
    ap.add_argument("--max_events_clasdis", type=int, default=DEFAULT_MAX_EVENTS_CLASDIS,
                    help="Max clasdis events for phase 1 (-1 means all)")

    ap.add_argument("--manual_weights", default=None,
                    help="Optional manual weights file with lines: row col weight")
    ap.add_argument("--zero_aaogen_rows", default="",
                    help="Comma-separated row indices to force w=0, for example 1,2,3")
    ap.add_argument("--force", action="store_true",
                    help="Apply built-in diagnostic force map after automatic fit")

    ap.add_argument("--weights_out", default=OUTPUT_WEIGHTS_TXT,
                    help="Weights report output path")
    ap.add_argument("--hist_csv_out", default=OUTPUT_HIST_CSV,
                    help="Histogram CSV output path")
    ap.add_argument("--yields_out", default=OUTPUT_YIELDS_PNG,
                    help="Per-pad three-way plot output path")
    ap.add_argument("--yields_data_only_out", default=OUTPUT_DATAONLY_PNG,
                    help="Integrated three-way plot output path")
    ap.add_argument("--yields_mix_out", default=OUTPUT_MIX_PNG,
                    help="Per-pad mixed plot output path")
    ap.add_argument("--yields_mix_integrated_out", default=OUTPUT_MIX_INTEGRATED_PNG,
                    help="Integrated mixed plot output path")

    ap.add_argument("--mc_aaogen", default=None,
                    help="Optional override aaogen ROOT file for output-tree writing")
    ap.add_argument("--mc_clasdis", default=None,
                    help="Optional override clasdis ROOT file for output-tree writing")
    ap.add_argument("--max_events_mc_aaogen", type=int, default=DEFAULT_MAX_EVENTS_MC_AAOGEN,
                    help="Max aaogen events for output-tree writing (-1 means all)")
    ap.add_argument("--max_events_mc_clasdis", type=int, default=DEFAULT_MAX_EVENTS_MC_CLASDIS,
                    help="Max clasdis events for output-tree writing (-1 means all)")
    ap.add_argument("--out", default=DEFAULT_MIXED_MC_ROOT,
                    help="Output ROOT file for mixed MC")

    args = ap.parse_args()

    require_file(args.data)
    require_file(args.aaogen)
    require_file(args.clasdis)

    if args.plot_nbins <= 0:
        fatal("--plot_nbins must be positive.")
    #endif
    if args.fit_nbins <= 0:
        fatal("--fit_nbins must be positive.")
    #endif
    if args.peak_min >= args.peak_max:
        fatal("Invalid peak window: peak_min must be < peak_max.")
    #endif
    if args.sideband_min >= args.sideband_max:
        fatal("Invalid sideband window: sideband_min must be < sideband_max.")
    #endif
    if args.w_scan_steps <= 0:
        fatal("--w_scan_steps must be positive.")
    #endif
    if args.lambda_sideband < 0.0:
        fatal("--lambda_sideband must be >= 0.")
    #endif

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)

    ensure_outdir(args.yields_out)
    ensure_outdir(args.yields_data_only_out)
    ensure_outdir(args.yields_mix_out)
    ensure_outdir(args.yields_mix_integrated_out)
    ensure_outdir(args.weights_out)
    ensure_outdir(args.hist_csv_out)
    ensure_outdir(args.out)
    ensure_outdir(OUTPUT_MIX_DEBUG_TXT)
    ensure_outdir(OUTPUT_MIX_DEBUG_MX2_PNG)

    f_data, t_data = open_tree(args.data, TREE_NAME)
    f_aao, t_aao = open_tree(args.aaogen, TREE_NAME)
    f_dis, t_dis = open_tree(args.clasdis, TREE_NAME)

    needed_phase1 = ["runnum", "e_p", "e_theta", "e_phi", "p_p", "p_theta", "p_phi", "x", "Q2", "Mx2"]
    require_branches(t_data, needed_phase1, "data")
    require_branches(t_aao, needed_phase1, "aaogen")
    require_branches(t_dis, needed_phase1, "clasdis")

    nrows = len(XB_EDGES) - 1
    ncols = len(TNEG_EDGES) - 1

    h_data_plot_raw = make_hist_grid("h_data_plot_raw", nrows, ncols,
                                     args.plot_nbins, MX2_PLOT_MIN, MX2_PLOT_MAX)
    h_aao_plot_raw = make_hist_grid("h_aao_plot_raw", nrows, ncols,
                                    args.plot_nbins, MX2_PLOT_MIN, MX2_PLOT_MAX)
    h_dis_plot_raw = make_hist_grid("h_dis_plot_raw", nrows, ncols,
                                    args.plot_nbins, MX2_PLOT_MIN, MX2_PLOT_MAX)

    h_data_fit_raw = make_hist_grid("h_data_fit_raw", nrows, ncols,
                                    args.fit_nbins, MX2_FIT_MIN_GLOBAL, MX2_FIT_MAX_GLOBAL)
    h_aao_fit_raw = make_hist_grid("h_aao_fit_raw", nrows, ncols,
                                   args.fit_nbins, MX2_FIT_MIN_GLOBAL, MX2_FIT_MAX_GLOBAL)
    h_dis_fit_raw = make_hist_grid("h_dis_fit_raw", nrows, ncols,
                                   args.fit_nbins, MX2_FIT_MIN_GLOBAL, MX2_FIT_MAX_GLOBAL)

    c_data = make_count_grid(nrows, ncols)
    c_aao = make_count_grid(nrows, ncols)
    c_dis = make_count_grid(nrows, ncols)

    max_events_data = parse_max_events(args.max_events_data)
    max_events_aao = parse_max_events(args.max_events_aaogen)
    max_events_dis = parse_max_events(args.max_events_clasdis)

    fill_all_bins_single_pass(t_data, h_data_plot_raw, h_data_fit_raw, c_data, max_events_data)
    fill_all_bins_single_pass(t_aao, h_aao_plot_raw, h_aao_fit_raw, c_aao, max_events_aao)
    fill_all_bins_single_pass(t_dis, h_dis_plot_raw, h_dis_fit_raw, c_dis, max_events_dis)

    h_data_plot_norm = clone_grid(h_data_plot_raw, "h_data_plot_norm")
    h_aao_plot_norm = clone_grid(h_aao_plot_raw, "h_aao_plot_norm")
    h_dis_plot_norm = clone_grid(h_dis_plot_raw, "h_dis_plot_norm")

    h_data_fit_norm = clone_grid(h_data_fit_raw, "h_data_fit_norm")
    h_aao_fit_norm = clone_grid(h_aao_fit_raw, "h_aao_fit_norm")
    h_dis_fit_norm = clone_grid(h_dis_fit_raw, "h_dis_fit_norm")

    for r in range(nrows):
        for c in range(ncols):
            normalize_unit_area(h_data_plot_norm[r][c])
            normalize_unit_area(h_aao_plot_norm[r][c])
            normalize_unit_area(h_dis_plot_norm[r][c])

            normalize_unit_area(h_data_fit_norm[r][c])
            normalize_unit_area(h_aao_fit_norm[r][c])
            normalize_unit_area(h_dis_fit_norm[r][c])
        #endfor
    #endfor

    Nd_tot = sum_counts_grid(c_data)
    Na_tot = sum_counts_grid(c_aao)
    Nc_tot = sum_counts_grid(c_dis)

    h_data_int = make_integrated_hist_from_grid(h_data_plot_raw, "h_data_int_raw")
    h_aao_int = make_integrated_hist_from_grid(h_aao_plot_raw, "h_aao_int_raw")
    h_dis_int = make_integrated_hist_from_grid(h_dis_plot_raw, "h_dis_int_raw")

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
            style_hist(h_data_plot_norm[r][c], col_data, 2, 1)
            style_hist(h_aao_plot_norm[r][c], col_aao, 2, 2)
            style_hist(h_dis_plot_norm[r][c], col_dis, 2, 3)
        #endfor
    #endfor

    draw_canvas_threeway(
        h_data_plot_norm,
        h_aao_plot_norm,
        h_dis_plot_norm,
        c_data,
        c_aao,
        c_dis,
        args.yields_out
    )

    draw_canvas_integrated_threeway(
        h_data_int,
        h_aao_int,
        h_dis_int,
        Nd_tot,
        Na_tot,
        Nc_tot,
        args.yields_data_only_out
    )

    fit_result = compute_w_grid_and_mix(
        h_data_fit_norm,
        h_aao_fit_norm,
        h_dis_fit_norm,
        h_data_fit_raw,
        h_aao_plot_norm,
        h_dis_plot_norm,
        args.peak_min,
        args.peak_max,
        args.sideband_min,
        args.sideband_max,
        args.lambda_sideband,
        args.w_scan_steps
    )

    w_grid = fit_result["w_grid"]
    wun_grid = fit_result["wun_grid"]

    manual_override_map = {}
    forced_map = {}

    zero_rows = parse_int_list(args.zero_aaogen_rows)
    if len(zero_rows) > 0:
        apply_zero_aaogen_rows(w_grid, zero_rows)
    #endif

    if args.manual_weights is not None:
        manual_override_map = parse_manual_weights(args.manual_weights)
        apply_manual_weights_map(w_grid, manual_override_map)
    #endif

    if args.force:
        forced_map = apply_force_map(w_grid)
    #endif

    diag_after_overrides = compute_objective_grid_for_existing_weights(
        w_grid,
        h_data_fit_norm,
        h_aao_fit_norm,
        h_dis_fit_norm,
        h_data_fit_raw,
        args.peak_min,
        args.peak_max,
        args.sideband_min,
        args.sideband_max,
        args.lambda_sideband
    )

    h_mix_plot = rebuild_mix_plot_grid(w_grid, h_aao_plot_norm, h_dis_plot_norm)

    for r in range(nrows):
        for c in range(ncols):
            style_hist(h_mix_plot[r][c], col_mix, 3, 1)
        #endfor
    #endfor

    total_obj = 0.0
    total_peak = 0.0
    total_side = 0.0
    for r in range(nrows):
        for c in range(ncols):
            total_obj += diag_after_overrides["objective_grid"][r][c]
            total_peak += diag_after_overrides["peak_chi2_grid"][r][c]
            total_side += diag_after_overrides["sideband_penalty_grid"][r][c]
        #endfor
    #endfor

    print("Automatic mixture determination complete")
    print(f"  peak window     = [{args.peak_min:.3f}, {args.peak_max:.3f}]")
    print(f"  sideband window = [{args.sideband_min:.3f}, {args.sideband_max:.3f}]")
    print(f"  lambda_sideband = {args.lambda_sideband:.6f}")
    print(f"  total objective = {total_obj:.6e}")
    print(f"  total peak chi2 = {total_peak:.6e}")
    print(f"  total sideband penalty = {total_side:.6e}")
    print(f"  data events processed    = {Nd_tot}")
    print(f"  aaogen events processed  = {Na_tot}")
    print(f"  clasdis events processed = {Nc_tot}")
    print(f"  wrote weights report: {args.weights_out}")
    print(f"  wrote histogram CSV:  {args.hist_csv_out}")

    write_weights_report(
        args.weights_out,
        w_grid,
        wun_grid,
        diag_after_overrides["objective_grid"],
        diag_after_overrides["peak_chi2_grid"],
        diag_after_overrides["sideband_penalty_grid"],
        diag_after_overrides["sideband_data_fraction_grid"],
        diag_after_overrides["sideband_mix_fraction_grid"],
        diag_after_overrides["sideband_sigma_grid"],
        c_data,
        c_aao,
        c_dis,
        args.peak_min,
        args.peak_max,
        args.sideband_min,
        args.sideband_max,
        args.lambda_sideband,
        forced_map,
        manual_override_map
    )

    draw_canvas_mix(
        h_data_plot_norm,
        h_mix_plot,
        c_data,
        w_grid,
        args.yields_mix_out,
        args.peak_min,
        args.peak_max,
        args.sideband_min,
        args.sideband_max
    )

    h_mix_int = make_integrated_hist_from_grid(h_mix_plot, "h_mix_int")
    normalize_unit_area(h_mix_int)
    style_hist(h_mix_int, col_mix, 3, 1)

    draw_canvas_integrated_mix(
        h_data_int,
        h_mix_int,
        Nd_tot,
        args.yields_mix_integrated_out
    )

    dump_hist_shapes_csv(
        args.hist_csv_out,
        h_data_fit_raw,
        h_aao_fit_raw,
        h_dis_fit_raw,
        h_data_fit_norm,
        h_aao_fit_norm,
        h_dis_fit_norm,
        rebuild_mix_plot_grid(w_grid, h_aao_fit_norm, h_dis_fit_norm),
        w_grid
    )

    mc_aao_path, mc_dis_path = choose_phase3_input_paths(args)

    max_events_mc_aao = parse_max_events(args.max_events_mc_aaogen)
    max_events_mc_dis = parse_max_events(args.max_events_mc_clasdis)

    print("Phase 3: writing mixed MC ROOT file")
    print(f"  phase3 aaogen input = {mc_aao_path}")
    print(f"  phase3 clasdis input = {mc_dis_path}")
    print(f"  out = {args.out}")
    print(f"  max_events_mc_aao = {args.max_events_mc_aaogen}")
    print(f"  max_events_mc_dis = {args.max_events_mc_clasdis}")
    print(f"  fixed Eb for MC mix = {MC_EB_FIXED:.3f}")
    print(f"  debug report = {OUTPUT_MIX_DEBUG_TXT}")
    print(f"  debug png = {OUTPUT_MIX_DEBUG_MX2_PNG}")

    mix_mc_to_output_root(
        mc_aao_path,
        mc_dis_path,
        args.out,
        w_grid,
        max_events_mc_aao,
        max_events_mc_dis
    )

    print(f"  wrote mixed MC: {args.out}")

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