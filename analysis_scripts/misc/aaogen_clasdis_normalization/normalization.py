#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import math
import argparse

import ROOT


# ---------------------------
# Physics helpers (from your C++)
# ---------------------------

MASS_E  = 0.000511
MASS_PI = 0.139570
MASS_N  = 0.9382720813


def beam_energy(runnum):
    if 6616 <= runnum <= 6783:
        return 10.1998
    if 16042 <= runnum <= 17065:
        return 10.5473
    if 17067 <= runnum <= 17724:
        return 10.5563
    if 17725 <= runnum <= 17811:
        return 10.5593
    return 10.5563


def compute_t_scalar(runnum, e_p, e_theta, e_phi, p_p, p_theta, p_phi):
    Eb = beam_energy(runnum)
    if Eb <= 0.0:
        return 1.0e9

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


def compute_tmin_exact(xB, Q2):
    xb_ok = (xB > 0.0 and xB < 1.0)
    if Q2 <= 0.0 or (not xb_ok):
        if xb_ok:
            denom = (1.0 - xB)
            if denom > 0.0:
                return - (MASS_N * xB) * (MASS_N * xB) / denom
        return 0.0

    eps2 = 4.0 * MASS_N * MASS_N * xB * xB / Q2
    root = math.sqrt(1.0 + eps2)
    num = Q2 * (2.0 * (1.0 - xB) * (1.0 - root) + eps2)
    den = 4.0 * xB * (1.0 - xB) + eps2
    if den == 0.0:
        return 0.0
    return - num / den


def compute_tprime(runnum, e_p, e_theta, e_phi, p_p, p_theta, p_phi, xB, Q2):
    t_val = compute_t_scalar(runnum, e_p, e_theta, e_phi, p_p, p_theta, p_phi)
    tmin_val = compute_tmin_exact(xB, Q2)
    return t_val - tmin_val


# ---------------------------
# Binning helpers
# ---------------------------

def find_bin(x, edges):
    # returns i such that edges[i] <= x < edges[i+1], else -1
    for i in range(len(edges) - 1):
        if edges[i] <= x < edges[i + 1]:
            return i
    return -1


# ---------------------------
# Histogram / fit helpers
# ---------------------------

def hist_integral_in_range(h, xmin, xmax):
    b1 = h.FindBin(xmin + 1.0e-12)
    b2 = h.FindBin(xmax - 1.0e-12)
    if b2 < b1:
        return 0.0
    return h.Integral(b1, b2)


def vector_from_hist_in_range(h, xmin, xmax):
    b1 = h.FindBin(xmin + 1.0e-12)
    b2 = h.FindBin(xmax - 1.0e-12)
    v = []
    for b in range(b1, b2 + 1):
        v.append(h.GetBinContent(b))
    #endfor
    return v


def normalize_vector(v):
    s = 0.0
    for x in v:
        s += x
    #endfor
    if s <= 0.0:
        return None
    out = []
    for x in v:
        out.append(x / s)
    #endfor
    return out


def dot(u, v):
    s = 0.0
    for i in range(len(u)):
        s += u[i] * v[i]
    #endfor
    return s


def sse(u, v):
    s = 0.0
    for i in range(len(u)):
        d = u[i] - v[i]
        s += d * d
    #endfor
    return s


def best_w_shape_only(D, A, C):
    # Minimizes || D - (w A + (1-w) C) ||^2 for w in [0,1]
    # Equivalent: minimize || (D-C) - w (A-C) ||^2
    AC = []
    DC = []
    for i in range(len(D)):
        AC.append(A[i] - C[i])
        DC.append(D[i] - C[i])
    #endfor

    denom = dot(AC, AC)
    if denom <= 0.0:
        # A and C identical in this window -> any w is equivalent; pick 0
        return 0.0
    w = dot(DC, AC) / denom
    if w < 0.0:
        w = 0.0
    if w > 1.0:
        w = 1.0
    return w


def clone_and_scale(h, name, factor):
    out = h.Clone(name)
    out.SetDirectory(0)
    out.Scale(factor)
    return out


def area_normalize_hist(h, xmin, xmax):
    integ = hist_integral_in_range(h, xmin, xmax)
    if integ > 0.0:
        h.Scale(1.0 / integ)
    #endif


# ---------------------------
# Main
# ---------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", required=True, help="data ROOT file (PhysicsEvents)")
    parser.add_argument("--aaogen", required=True, help="aaogen ROOT file (PhysicsEvents)")
    parser.add_argument("--clasdis", required=True, help="clasdis ROOT file (PhysicsEvents)")
    args = parser.parse_args()

    # Plot settings
    x_edges = [0.10, 0.25, 0.35, 0.45, 0.60]
    tp_edges = [0.05, 0.25, 0.45, 0.65, 0.85, 1.05, 1.25]
    nrows = len(x_edges) - 1
    ncols = len(tp_edges) - 1

    mx2_plot_min = 0.0
    mx2_plot_max = 2.0

    fit_min = 0.4
    fit_max = 1.5

    nbins_mx2 = 80

    out_png = "output/yields_mix.png"
    out_txt = "output/weights.txt"
    os.makedirs("output", exist_ok=True)

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetTitleFontSize(0.06)

    # Open files/trees
    def open_tree(path, label):
        f = ROOT.TFile.Open(path, "READ")
        if (not f) or f.IsZombie():
            raise RuntimeError("FATAL: could not open %s file: %s" % (label, path))
        t = f.Get("PhysicsEvents")
        if not t:
            raise RuntimeError("FATAL: PhysicsEvents not found in %s: %s" % (label, path))
        return f, t

    f_data, t_data = open_tree(args.data, "data")
    f_aao, t_aao = open_tree(args.aaogen, "aaogen")
    f_dis, t_dis = open_tree(args.clasdis, "clasdis")

    # Required branches for computing tprime and filling
    needed = ["runnum", "e_p", "e_theta", "e_phi", "p_p", "p_theta", "p_phi", "x", "Q2", "Mx2"]

    def check_branches(t, label):
        missing = []
        for b in needed:
            if not t.GetBranch(b):
                missing.append(b)
            #endif
        #endfor
        if len(missing) > 0:
            raise RuntimeError("FATAL: Missing required branches in %s: %s" % (label, ", ".join(missing)))
        #endif

    check_branches(t_data, "data")
    check_branches(t_aao, "aaogen")
    check_branches(t_dis, "clasdis")

    # Allocate histograms: [row][col] for each sample
    def make_hist_set(prefix):
        H = []
        for r in range(nrows):
            row = []
            for c in range(ncols):
                name = "%s_r%d_c%d" % (prefix, r, c)
                h = ROOT.TH1F(name, name, nbins_mx2, mx2_plot_min, mx2_plot_max)
                h.SetDirectory(0)
                row.append(h)
            #endfor
            H.append(row)
        #endfor
        return H

    H_data = make_hist_set("h_data")
    H_aao = make_hist_set("h_aao")
    H_dis = make_hist_set("h_dis")

    # One-pass fill
    def fill_from_tree(t, H, label):
        n = t.GetEntries()
        for i in range(n):
            t.GetEntry(i)

            runnum = int(getattr(t, "runnum"))
            e_p = float(getattr(t, "e_p"))
            e_theta = float(getattr(t, "e_theta"))
            e_phi = float(getattr(t, "e_phi"))
            p_p = float(getattr(t, "p_p"))
            p_theta = float(getattr(t, "p_theta"))
            p_phi = float(getattr(t, "p_phi"))
            xB = float(getattr(t, "x"))
            Q2 = float(getattr(t, "Q2"))
            mx2 = float(getattr(t, "Mx2"))

            r = find_bin(xB, x_edges)
            if r < 0:
                continue
            #endif

            tprime = compute_tprime(runnum, e_p, e_theta, e_phi, p_p, p_theta, p_phi, xB, Q2)
            tp = -tprime

            c = find_bin(tp, tp_edges)
            if c < 0:
                continue
            #endif

            H[r][c].Fill(mx2)

        #endfor

    fill_from_tree(t_data, H_data, "data")
    fill_from_tree(t_aao, H_aao, "aaogen")
    fill_from_tree(t_dis, H_dis, "clasdis")

    # Determine per-bin weights (shape-only in fit window), then draw mix vs data
    weights = [[0.0 for _ in range(ncols)] for __ in range(nrows)]
    perbin_sse = [[0.0 for _ in range(ncols)] for __ in range(nrows)]
    total_sse = 0.0

    # Prepare canvas
    can_w = 350 * ncols + 220
    can_h = 260 * nrows + 220
    c = ROOT.TCanvas("c", "c", can_w, can_h)
    c.Divide(ncols, nrows, 0.001, 0.001)

    # For a single legend placed on the last pad
    legend_made = False

    for r in range(nrows):
        for cc in range(ncols):
            hD_raw = H_data[r][cc]
            hA_raw = H_aao[r][cc]
            hC_raw = H_dis[r][cc]

            # Build normalized vectors in the fit window
            vD = normalize_vector(vector_from_hist_in_range(hD_raw, fit_min, fit_max))
            vA = normalize_vector(vector_from_hist_in_range(hA_raw, fit_min, fit_max))
            vC = normalize_vector(vector_from_hist_in_range(hC_raw, fit_min, fit_max))

            # If any are empty in the fit window, weight is not meaningful; default to 0
            if (vD is None) or (vA is None) or (vC is None):
                w = 0.0
                s = 0.0
            else:
                w = best_w_shape_only(vD, vA, vC)
                vMix = []
                for k in range(len(vD)):
                    vMix.append(w * vA[k] + (1.0 - w) * vC[k])
                #endfor
                s = sse(vD, vMix)
            #endif

            weights[r][cc] = w
            perbin_sse[r][cc] = s
            total_sse += s

            # Build mixture histogram in raw counts, then area-normalize for plotting
            hMix = clone_and_scale(hA_raw, "h_mix_r%d_c%d" % (r, cc), w)
            hTmpC = clone_and_scale(hC_raw, "h_tmpc_r%d_c%d" % (r, cc), (1.0 - w))
            hMix.Add(hTmpC)

            # Normalize data + mix in the plotting range so shapes are comparable
            hD = hD_raw.Clone("hD_plot_r%d_c%d" % (r, cc))
            hD.SetDirectory(0)
            area_normalize_hist(hD, mx2_plot_min, mx2_plot_max)
            area_normalize_hist(hMix, mx2_plot_min, mx2_plot_max)

            # Styling
            hD.SetLineColor(ROOT.kBlack)
            hD.SetLineWidth(2)
            hMix.SetLineColor(ROOT.kGreen + 2)
            hMix.SetLineWidth(2)

            # Determine per-pad ymax from the two plotted histograms (avoid 0 / collapse)
            maxD = hD.GetMaximum()
            maxM = hMix.GetMaximum()
            ymax = 1.2 * max(maxD, maxM)
            if ymax <= 0.0:
                ymax = 1.0
            #endif

            pad = c.cd(r * ncols + cc + 1)
            pad.SetLeftMargin(0.18)
            pad.SetRightMargin(0.05)
            pad.SetBottomMargin(0.16)
            pad.SetTopMargin(0.18)
            pad.SetGrid(1, 1)

            title = "xB [%.2f, %.2f)  -tprime [%.2f, %.2f)" % (
                x_edges[r], x_edges[r + 1], tp_edges[cc], tp_edges[cc + 1]
            )

            # Always draw a frame FIRST so axes never disappear
            frame = pad.DrawFrame(mx2_plot_min, 0.0, mx2_plot_max, ymax)
            frame.SetTitle(title)
            frame.GetXaxis().SetTitle("Mx2 (GeV2)")
            frame.GetYaxis().SetTitle("Normalized yield")
            frame.GetXaxis().SetTitleSize(0.06)
            frame.GetYaxis().SetTitleSize(0.06)
            frame.GetXaxis().SetLabelSize(0.05)
            frame.GetYaxis().SetLabelSize(0.05)
            frame.GetYaxis().SetTitleOffset(1.45)

            # Draw the histograms on top of the frame
            hD.Draw("hist same")
            hMix.Draw("hist same")

            # Put a legend on the very last pad only
            if (not legend_made) and (r == nrows - 1) and (cc == ncols - 1):
                leg = ROOT.TLegend(0.52, 0.62, 0.93, 0.92)
                leg.SetFillStyle(1001)
                leg.SetFillColor(ROOT.kWhite)
                leg.SetBorderSize(1)
                leg.SetTextSize(0.04)
                leg.AddEntry(hD, "data", "l")
                leg.AddEntry(hMix, "mix: w=%.4f" % (w), "l")
                leg.AddEntry(0, "fit window: [%.1f, %.1f]" % (fit_min, fit_max), "")
                leg.Draw()
                legend_made = True
            #endif

        #endfor
    #endfor

    c.SaveAs(out_png)

    # Write weights report
    with open(out_txt, "w") as f:
        f.write("Per-bin mixture weights (shape-only, analytic w, unweighted SSE)\n")
        f.write("Definition: H_mix = w * H_aaogen + (1-w) * H_clasdis\n")
        f.write("\n")
        f.write("plot range: Mx2 in [%.3f, %.3f]\n" % (mx2_plot_min, mx2_plot_max))
        f.write("fit window: Mx2 in [%.3f, %.3f]\n" % (fit_min, fit_max))
        f.write("total SSE  = %.6e\n" % total_sse)
        f.write("\n")

        for r in range(nrows):
            f.write("Row %d: xB [%.2f, %.2f)\n" % (r, x_edges[r], x_edges[r + 1]))
            for cc in range(ncols):
                Nd = int(H_data[r][cc].GetEntries())
                Na = int(H_aao[r][cc].GetEntries())
                Nc = int(H_dis[r][cc].GetEntries())
                f.write("  Col %d: -tprime [%.2f, %.2f)  w=%.6f  SSE=%.6e  N(data,aao,dis)=(%d,%d,%d)\n" % (
                    cc, tp_edges[cc], tp_edges[cc + 1],
                    weights[r][cc], perbin_sse[r][cc],
                    Nd, Na, Nc
                ))
            #endfor
            f.write("\n")
        #endfor

    print("plot range: Mx2 in [%.3f, %.3f]" % (mx2_plot_min, mx2_plot_max))
    print("fit window: Mx2 in [%.3f, %.3f]" % (fit_min, fit_max))
    print("total SSE  = %.6e" % total_sse)
    print("wrote plot: %s" % out_png)
    print("wrote weights report: %s" % out_txt)


if __name__ == "__main__":
    main()
#endif