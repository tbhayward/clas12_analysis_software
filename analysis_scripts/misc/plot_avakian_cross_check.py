#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
avakian_cross_check_multipage.py

Build a multi-page PDF comparing Hayward vs Avakian across several derived quantities.

Inputs (runtime):
  Three Hayward per-bin export text files with header:
    "# x  -tprime  Df  AULsin  eAULsin  <cosphi>  e<cosphi>  aulsin2phi  eaulsin2phi  <V/A>  <B/A>"
  (See examples you posted.)

Hard-coded:
  Avakian's three 5-point blocks.

Output:
  output/enpi+/avakian_cross_check_multipage.pdf

Pages:
  1) 1x3: Df vs x_B (y:[0.2,0.6], x:[0.05,0.65]), titles carry -t' range.
     Markers: Avakian = CLOSED circles; Hayward = OPEN circles.
  2) 1x3: (V/A) and (B/A) vs x_B (y:[0,2], x:[0.05,0.65])
     Markers: Hayward = CLOSED circles; Avakian = OPEN circles.
  3) 1x3: <x> and <-t'> vs x_B (y:[0,1], x:[0.05,0.65])
     Markers: Hayward = CLOSED circles; Avakian = OPEN circles.
  4) 2x3: A_{UL}^{sinφ} (top), A_{UL}^{sin2φ} (bottom) vs x_B (auto symmetric y-lims per row)
     Markers: Hayward = CLOSED circles; Avakian = OPEN circles.
     NOTE: These are the *asymmetry amplitudes* (not ratios).
  5) 2x3: F_{UL}^{sinφ}/F_{UU} (top), F_{UL}^{sin2φ}/F_{UU} (bottom) vs x_B
     Computed from per-bin amplitudes as:
       ratio_sinφ   = AULsin   / (V/A)
       ratio_sin2φ  = AULsin2  / (B/A)
     Fixed y-lims:
       sinφ:   [-0.35, 0.25]
       sin2φ:  [-0.60, 0.15]
     Markers: Hayward = CLOSED circles; Avakian = OPEN circles.

Conventions:
  - t-bin titles (using -t′ ranges): ["0.05 < -t' < 0.45", "0.45 < -t' < 0.85", "0.85 < -t' < 1.225"]
  - Colors for t-bins: ["tab:blue", "tab:orange", "tab:green"]
  - Colors for Page 2 (V/A, B/A): V/A="tab:purple", B/A="tab:red"
  - Colors for Page 3 (<x>, <-t′>): <x>="tab:blue", <-t′>="tab:olive"

"""

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages

# ----------------------------
# Config / constants
# ----------------------------
X_MIN, X_MAX = 0.05, 0.65

DF_YMIN,   DF_YMAX   = 0.20, 0.60
VBA_YMIN,  VBA_YMAX  = 0.00, 2.00
XT_YMIN,   XT_YMAX   = 0.00, 1.00

R_SIN_YMIN,  R_SIN_YMAX  = -0.35, 0.25
R_SIN2_YMIN, R_SIN2_YMAX = -0.60, 0.15

OUTPATH = "output/enpi+/avakian_cross_check_multipage.pdf"

tbin_titles = ["0.05 < -t' < 0.45", "0.45 < -t' < 0.85", "0.85 < -t' < 1.225"]
tbin_colors = ["tab:blue", "tab:orange", "tab:green"]

# Page 2 palette (variables inside each t-panel)
COL_VA = "tab:purple"
COL_BA = "tab:red"

# Page 3 palette (variables inside each t-panel)
COL_X   = "tab:blue"
COL_TPR = "tab:olive"

# ----------------------------
# Avakian (colleague) blocks (5 rows each, 11 cols)
# Columns (1-based):
#   1: xB
#   2: <-t'>
#   3: Df
#   4: AULsinphi
#   5: eAULsinphi
#   6: (unused here)
#   7: (unused here)
#   8: AULsin2phi
#   9: eAULsin2phi
#  10: <V/A>     (interpreted)
#  11: <B/A>     (interpreted)
# ----------------------------
avak_block1 = np.array([
    [0.1709, 0.1721, 0.3714, -0.1960E-01, 0.1160E-01, 0.5400, 0.7950E-01, -0.2010E-01, 0.1060E-01, 1.568, 0.7172],
    [0.2538, 0.1675, 0.3817,  0.7600E-02, 0.6300E-02, 0.4999, 0.5580E-01, -0.2830E-01, 0.5800E-02, 1.713, 0.8107],
    [0.3463, 0.1485, 0.3933,  0.5680E-01, 0.6200E-02, 0.4600, 0.5660E-01, -0.1270E-01, 0.5800E-02, 1.790, 0.8610],
    [0.4326, 0.1058, 0.4041,  0.6920E-01, 0.1090E-01, 0.4200, 0.8000E-01,  0.1190E-01, 0.1040E-01, 1.797, 0.8656],
    [0.5092, 0.0000, 0.4137,  0.8080E-01, 0.8400E-01, 0.3800, 0.5690E-01, -0.4470E-01, 0.8210E-01, 1.781, 0.8552],
], dtype=float)

avak_block2 = np.array([
    [0.1679, 0.5951, 0.3710, -0.1415,      0.2250E-01, 0.5355, 0.7160E-01, -0.1272,      0.2140E-01, 1.554, 0.7082],
    [0.2588, 0.5381, 0.3824, -0.4270E-01,  0.1090E-01, 0.4200, 0.6350E-01, -0.1054,      0.1100E-01, 1.719, 0.8152],
    [0.3525, 0.4503, 0.3941,  0.2740E-01,  0.7600E-02, 0.4004, 0.4970E-01, -0.6440E-01,  0.9900E-02, 1.790, 0.8617],
    [0.4460, 0.3051, 0.4057,  0.6470E-01,  0.8700E-02, 0.4200, 0.6790E-01, -0.4740E-01,  0.8600E-02, 1.780, 0.8546],
    [0.5331, 0.1758, 0.4167,  0.6900E-02,  0.1460E-01, 0.3800, 0.5640E-01, -0.4720E-01,  0.1400E-01, 1.730, 0.8217],
], dtype=float)

avak_block3 = np.array([
    [0.2103, 0.9800, 0.3763,  0.2430E-01,  0.1940E-01, -0.4600, 0.7960E-01, -0.4390E-01, 0.2110E-01, 1.640, 0.7637],
    [0.3055, 0.9099, 0.3882,  0.6780E-01,  0.1250E-01, -0.4200, 0.6100E-01, -0.1189,     0.1330E-01, 1.762, 0.8433],
    [0.3995, 0.7925, 0.4000,  0.9340E-01,  0.1080E-01, -0.3800, 0.5640E-01, -0.1179,     0.1130E-01, 1.796, 0.8651],
    [0.4946, 0.5996, 0.4118,  0.8190E-01,  0.1340E-01, -0.3400, 0.7740E-01, -0.7610E-01, 0.1340E-01, 1.752, 0.8362],
    [0.5881, 0.3125, 0.4234,  0.4350E-01,  0.2130E-01, -0.3000, 0.5690E-01, -0.2750E-01, 0.2120E-01, 1.674, 0.7851],
], dtype=float)

avak_blocks = [avak_block1, avak_block2, avak_block3]

# ----------------------------
# Utilities
# ----------------------------
def ensure_outdir(path: str):
    outdir = os.path.dirname(path)
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    #endif
#endfor

def load_hayward_export(path):
    """
    Read one Hayward per-bin export file (your new text format).
    Returns a dict of numpy arrays with keys:
      x, tprime, Df, AULsin, eAULsin, cosphi, e_cosphi, AULsin2, eAULsin2, VoverA, BoverA
    """
    xs, tps, dfs = [], [], []
    a1, ea1, cph, ecph = [], [], [], []
    a2, ea2, rva, rba = [], [], [], []

    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            # split by whitespace; expect 11 columns
            cols = s.split()
            if len(cols) < 11:
                continue
            #endfor
            x        = float(cols[0])
            tprime   = float(cols[1])
            Df       = float(cols[2])
            AULsin   = float(cols[3])
            eAULsin  = float(cols[4])
            cosphi   = float(cols[5])
            e_cosphi = float(cols[6])
            AULsin2  = float(cols[7])
            eAULsin2 = float(cols[8])
            VoverA   = float(cols[9])
            BoverA   = float(cols[10])

            xs.append(x); tps.append(tprime); dfs.append(Df)
            a1.append(AULsin); ea1.append(eAULsin)
            cph.append(cosphi); ecph.append(e_cosphi)
            a2.append(AULsin2); ea2.append(eAULsin2)
            rva.append(VoverA); rba.append(BoverA)
        #endfor
    #endwith

    return dict(
        x=np.array(xs, float),
        tprime=np.array(tps, float),
        Df=np.array(dfs, float),
        AULsin=np.array(a1, float),
        eAULsin=np.array(ea1, float),
        cosphi=np.array(cph, float),
        e_cosphi=np.array(ecph, float),
        AULsin2=np.array(a2, float),
        eAULsin2=np.array(ea2, float),
        VoverA=np.array(rva, float),
        BoverA=np.array(rba, float),
    )
#endfor

def split_avak_dicts():
    """
    Convert Avakian blocks into list of dicts compatible with Hayward dicts.
    Keys: x, tprime, Df, AULsin, eAULsin, AULsin2, eAULsin2, VoverA, BoverA
    (cosphi not available; omitted.)
    """
    out = []
    for blk in avak_blocks:
        d = dict(
            x=blk[:,0],
            tprime=blk[:,1],
            Df=blk[:,2],
            AULsin=blk[:,3],
            eAULsin=np.abs(blk[:,4]),
            AULsin2=blk[:,7],
            eAULsin2=np.abs(blk[:,8]),
            VoverA=blk[:,9],
            BoverA=blk[:,10],
        )
        out.append(d)
    #endfor
    return out
#endfor

def symmetric_ylim_from_series(series_list, pad=0.05, min_halfspan=0.05):
    """
    Given list of 1D arrays, compute symmetric y-lims around 0 with a small pad.
    """
    vmax = 0.0
    for s in series_list:
        if s is None: 
            continue
        s_abs = np.nanmax(np.abs(np.asarray(s, float)))
        vmax = max(vmax, s_abs)
    #endfor
    half = max(min_halfspan, (1.0 + pad) * vmax)
    return (-half, half)
#endfor

def add_panel_title(ax, idx):
    ax.set_title(tbin_titles[idx], fontsize=11)
#endfor

def legend_handles(dataset_closed_color="black", dataset_open_color="black",
                   closed_label="Hayward", open_label="Avakian"):
    h_closed = Line2D([0], [0], marker="o", linestyle="None",
                      mfc=dataset_closed_color, mec="black", ms=6, label=closed_label)
    h_open   = Line2D([0], [0], marker="o", linestyle="None",
                      mfc="none", mec=dataset_open_color, ms=6, label=open_label)
    return h_closed, h_open
#endfor

# ----------------------------
# Plotting pages
# ----------------------------
def page1_compare_df(pdf, hayward_dicts, avak_dicts):
    """Page 1: 1x3 Df vs x_B; y:[0.2,0.6]; x:[0.05,0.65].
       Avakian CLOSED, Hayward OPEN."""
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), sharex=True, sharey=True)
    for i in range(3):
        ax = axes[i]
        H, A = hayward_dicts[i], avak_dicts[i]
        # Avakian = closed
        ax.plot(A["x"], A["Df"], "o", ms=6, color=tbin_colors[i], mfc=tbin_colors[i], mec="black")
        # Hayward = open
        ax.plot(H["x"], H["Df"], "o", ms=6, color=tbin_colors[i], mfc="none", mec=tbin_colors[i])
        add_panel_title(ax, i)
        ax.set_xlim(X_MIN, X_MAX)
        ax.set_ylim(DF_YMIN, DF_YMAX)
        if i == 0:
            ax.set_ylabel("Dilution factor $D_{f}$")
        #endif
        ax.set_xlabel(r"$x_{B}$")
    #endfor

    # Legend
    h_closed = Line2D([0],[0], marker="o", linestyle="None", mfc="black", mec="black", ms=6, label="Avakian")
    h_open   = Line2D([0],[0], marker="o", linestyle="None", mfc="none",  mec="black", ms=6, label="Hayward")
    axes[0].legend(handles=[h_closed, h_open], loc="lower right", frameon=False, title="Dataset")

    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)
#endfor

def page2_compare_va_ba(pdf, hayward_dicts, avak_dicts):
    """Page 2: 1x3 compare <V/A> and <B/A>, y:[0,2]; x:[0.05,0.65].
       Hayward CLOSED, Avakian OPEN."""
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), sharex=True, sharey=True)
    for i in range(3):
        ax = axes[i]
        H, A = hayward_dicts[i], avak_dicts[i]

        # V/A (color COL_VA)
        ax.plot(H["x"], H["VoverA"], "o", ms=6, color=COL_VA, mfc=COL_VA, mec="black", label=r"$\langle V/A\rangle$ (H)")  # closed
        ax.plot(A["x"], A["VoverA"], "o", ms=6, color=COL_VA, mfc="none",  mec=COL_VA, label=r"$\langle V/A\rangle$ (A)")  # open

        # B/A (color COL_BA)
        ax.plot(H["x"], H["BoverA"], "o", ms=6, color=COL_BA, mfc=COL_BA, mec="black", label=r"$\langle B/A\rangle$ (H)")
        ax.plot(A["x"], A["BoverA"], "o", ms=6, color=COL_BA, mfc="none",  mec=COL_BA, label=r"$\langle B/A\rangle$ (A)")

        add_panel_title(ax, i)
        ax.set_xlim(X_MIN, X_MAX)
        ax.set_ylim(VBA_YMIN, VBA_YMAX)
        if i == 0:
            ax.set_ylabel(r"$\langle V/A\rangle,\,\langle B/A\rangle$")
        #endif
        ax.set_xlabel(r"$x_{B}$")
    #endfor

    # Build a simple legend: color encodes (V/A vs B/A); fill encodes (Hayward vs Avakian)
    h_va_h, h_va_a = (Line2D([0],[0], marker="o", linestyle="None", mfc=COL_VA, mec="black", ms=6, label=r"$\langle V/A\rangle$ (Hayward)"),
                      Line2D([0],[0], marker="o", linestyle="None", mfc="none",  mec=COL_VA, ms=6, label=r"$\langle V/A\rangle$ (Avakian)"))
    h_ba_h, h_ba_a = (Line2D([0],[0], marker="o", linestyle="None", mfc=COL_BA, mec="black", ms=6, label=r"$\langle B/A\rangle$ (Hayward)"),
                      Line2D([0],[0], marker="o", linestyle="None", mfc="none",  mec=COL_BA, ms=6, label=r"$\langle B/A\rangle$ (Avakian)"))
    axes[0].legend(handles=[h_va_h, h_va_a, h_ba_h, h_ba_a], loc="upper right", frameon=False, title="Vars/Datasets")

    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)
#endfor

def page3_compare_x_and_tprime(pdf, hayward_dicts, avak_dicts):
    """Page 3: 1x3 plot <x> and <-t'> vs x_B; y:[0,1]; x:[0.05,0.65].
       Hayward CLOSED, Avakian OPEN."""
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), sharex=True, sharey=True)
    for i in range(3):
        ax = axes[i]
        H, A = hayward_dicts[i], avak_dicts[i]

        # <x> (which is just x vs x) – included for direct comparison
        ax.plot(H["x"], H["x"], "o", ms=6, color=COL_X, mfc=COL_X, mec="black", label=r"$\langle x\rangle$ (H)")
        ax.plot(A["x"], A["x"], "o", ms=6, color=COL_X, mfc="none",  mec=COL_X, label=r"$\langle x\rangle$ (A)")

        # <-t'>
        ax.plot(H["x"], H["tprime"], "o", ms=6, color=COL_TPR, mfc=COL_TPR, mec="black", label=r"$\langle -t'\rangle$ (H)")
        ax.plot(A["x"], A["tprime"], "o", ms=6, color=COL_TPR, mfc="none",  mec=COL_TPR, label=r"$\langle -t'\rangle$ (A)")

        add_panel_title(ax, i)
        ax.set_xlim(X_MIN, X_MAX)
        ax.set_ylim(XT_YMIN, XT_YMAX)
        if i == 0:
            ax.set_ylabel(r"$\langle x\rangle,\ \langle -t'\rangle$")
        #endif
        ax.set_xlabel(r"$x_{B}$")
    #endfor

    h_x_h, h_x_a = (Line2D([0],[0], marker="o", linestyle="None", mfc=COL_X,   mec="black", ms=6, label=r"$\langle x\rangle$ (Hayward)"),
                    Line2D([0],[0], marker="o", linestyle="None", mfc="none",   mec=COL_X,   ms=6, label=r"$\langle x\rangle$ (Avakian)"))
    h_t_h, h_t_a = (Line2D([0],[0], marker="o", linestyle="None", mfc=COL_TPR, mec="black", ms=6, label=r"$\langle -t'\rangle$ (Hayward)"),
                    Line2D([0],[0], marker="o", linestyle="None", mfc="none",   mec=COL_TPR, ms=6, label=r"$\langle -t'\rangle$ (Avakian)"))
    axes[0].legend(handles=[h_x_h, h_x_a, h_t_h, h_t_a], loc="upper left", frameon=False, title="Vars/Datasets")

    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)
#endfor

def page4_compare_asymmetries(pdf, hayward_dicts, avak_dicts):
    """Page 4: 2x3; top A_UL^{sinφ}, bottom A_UL^{sin2φ} (asymmetry amplitudes)."""
    fig, axes = plt.subplots(2, 3, figsize=(12, 7.6), sharex=True)
    # Determine symmetric y-lims per row from both datasets
    row1_series = []
    row2_series = []
    for i in range(3):
        H, A = hayward_dicts[i], avak_dicts[i]
        row1_series.extend([H["AULsin"], A["AULsin"]])
        row2_series.extend([H["AULsin2"], A["AULsin2"]])
    #endfor
    y1min, y1max = symmetric_ylim_from_series(row1_series, pad=0.1, min_halfspan=0.05)
    y2min, y2max = symmetric_ylim_from_series(row2_series, pad=0.1, min_halfspan=0.05)

    for i in range(3):
        # Top row: AUL sinφ
        axT = axes[0, i]
        H, A = hayward_dicts[i], avak_dicts[i]
        axT.errorbar(H["x"], H["AULsin"],  yerr=H["eAULsin"],  fmt="o", ms=6, lw=1.2, capsize=3,
                     color=tbin_colors[i], mfc=tbin_colors[i], mec="black", label="Hayward")
        axT.errorbar(A["x"], A["AULsin"],  yerr=A["eAULsin"],  fmt="o", ms=6, lw=1.2, capsize=3,
                     color=tbin_colors[i], mfc="none", mec=tbin_colors[i], label="Avakian")
        add_panel_title(axT, i)
        axT.set_xlim(X_MIN, X_MAX)
        axT.set_ylim(y1min, y1max)
        if i == 0:
            axT.set_ylabel(r"$A_{UL}^{\sin\phi}$")
        #endif

        # Bottom row: AUL sin2φ
        axB = axes[1, i]
        axB.errorbar(H["x"], H["AULsin2"], yerr=H["eAULsin2"], fmt="o", ms=6, lw=1.2, capsize=3,
                     color=tbin_colors[i], mfc=tbin_colors[i], mec="black", label="Hayward")
        axB.errorbar(A["x"], A["AULsin2"], yerr=A["eAULsin2"], fmt="o", ms=6, lw=1.2, capsize=3,
                     color=tbin_colors[i], mfc="none", mec=tbin_colors[i], label="Avakian")
        axB.set_xlim(X_MIN, X_MAX)
        axB.set_ylim(y2min, y2max)
        if i == 0:
            axB.set_ylabel(r"$A_{UL}^{\sin2\phi}$")
        #endif
        axB.set_xlabel(r"$x_{B}$")
    #endfor

    # Single legend (top-left panel)
    h_h = Line2D([0],[0], marker="o", linestyle="None", mfc="black", mec="black", ms=6, label="Hayward")
    h_a = Line2D([0],[0], marker="o", linestyle="None", mfc="none",  mec="black", ms=6, label="Avakian")
    axes[0,0].legend(handles=[h_h, h_a], loc="lower right", frameon=False, title="Dataset")

    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)
#endfor

def page5_compare_ratios(pdf, hayward_dicts, avak_dicts):
    """Page 5: 2x3; top F_UL^{sinφ}/F_UU, bottom F_UL^{sin2φ}/F_UU (ratios).
       Using ratio = AUL / (V/A) and AUL2 / (B/A)."""
    fig, axes = plt.subplots(2, 3, figsize=(12, 7.6), sharex=True, sharey=False)

    for i in range(3):
        H, A = hayward_dicts[i], avak_dicts[i]

        # Ratios and errors via straightforward scaling
        H_r1  = H["AULsin"]  / H["VoverA"]
        H_er1 = H["eAULsin"] / H["VoverA"]
        H_r2  = H["AULsin2"]  / H["BoverA"]
        H_er2 = H["eAULsin2"] / H["BoverA"]

        A_r1  = A["AULsin"]  / A["VoverA"]
        A_er1 = A["eAULsin"] / A["VoverA"]
        A_r2  = A["AULsin2"]  / A["BoverA"]
        A_er2 = A["eAULsin2"] / A["BoverA"]

        # Top row: sinφ ratios
        axT = axes[0, i]
        axT.errorbar(H["x"], H_r1, yerr=H_er1, fmt="o", ms=6, lw=1.2, capsize=3,
                     color=tbin_colors[i], mfc=tbin_colors[i], mec="black", label="Hayward")
        axT.errorbar(A["x"], A_r1, yerr=A_er1, fmt="o", ms=6, lw=1.2, capsize=3,
                     color=tbin_colors[i], mfc="none",    mec=tbin_colors[i], label="Avakian")
        add_panel_title(axT, i)
        axT.set_xlim(X_MIN, X_MAX)
        axT.set_ylim(R_SIN_YMIN, R_SIN_YMAX)
        if i == 0:
            axT.set_ylabel(r"$F_{UL}^{\sin\phi}/F_{UU}$")
        #endif

        # Bottom row: sin2φ ratios
        axB = axes[1, i]
        axB.errorbar(H["x"], H_r2, yerr=H_er2, fmt="o", ms=6, lw=1.2, capsize=3,
                     color=tbin_colors[i], mfc=tbin_colors[i], mec="black", label="Hayward")
        axB.errorbar(A["x"], A_r2, yerr=A_er2, fmt="o", ms=6, lw=1.2, capsize=3,
                     color=tbin_colors[i], mfc="none",    mec=tbin_colors[i], label="Avakian")
        axB.set_xlim(X_MIN, X_MAX)
        axB.set_ylim(R_SIN2_YMIN, R_SIN2_YMAX)
        if i == 0:
            axB.set_ylabel(r"$F_{UL}^{\sin2\phi}/F_{UU}$")
        #endif
        axB.set_xlabel(r"$x_{B}$")
    #endfor

    # Single legend (top-left)
    h_h = Line2D([0],[0], marker="o", linestyle="None", mfc="black", mec="black", ms=6, label="Hayward")
    h_a = Line2D([0],[0], marker="o", linestyle="None", mfc="none",  mec="black", ms=6, label="Avakian")
    axes[0,0].legend(handles=[h_h, h_a], loc="lower right", frameon=False, title="Dataset")

    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)
#endfor

# ----------------------------
# Main
# ----------------------------
def main():
    parser = argparse.ArgumentParser(description="Build multi-page Avakian/Hayward cross-check PDF.")
    parser.add_argument("--hay1", required=True, help="Hayward per-bin export file (t-bin 1: 0.05 < -t' < 0.45)")
    parser.add_argument("--hay2", required=True, help="Hayward per-bin export file (t-bin 2: 0.45 < -t' < 0.85)")
    parser.add_argument("--hay3", required=True, help="Hayward per-bin export file (t-bin 3: 0.85 < -t' < 1.225)")
    parser.add_argument("--out",  default=OUTPATH, help=f"Output PDF path (default: {OUTPATH})")
    args = parser.parse_args()

    ensure_outdir(args.out)

    # Load Hayward data (3 files)
    H1 = load_hayward_export(args.hay1)
    H2 = load_hayward_export(args.hay2)
    H3 = load_hayward_export(args.hay3)
    hayward_dicts = [H1, H2, H3]

    # Prepare Avakian dicts
    avak_dicts = split_avak_dicts()

    # Build multipage PDF
    with PdfPages(args.out) as pdf:
        page1_compare_df(pdf, hayward_dicts, avak_dicts)
        page2_compare_va_ba(pdf, hayward_dicts, avak_dicts)
        page3_compare_x_and_tprime(pdf, hayward_dicts, avak_dicts)
        page4_compare_asymmetries(pdf, hayward_dicts, avak_dicts)
        page5_compare_ratios(pdf, hayward_dicts, avak_dicts)
    #endwith

    print("Wrote:", args.out)
#endfor

if __name__ == "__main__":
    main()
#endif