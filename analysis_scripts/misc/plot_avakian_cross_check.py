#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
harut_cross_check.py

Make a 1x2 figure:
  Left:  F_UL^{sin phi} / F_UU  vs x_B
  Right: F_UL^{sin 2phi} / F_UU vs x_B

Data:
  • Hayward (your three t-bins): filled circle markers with error bars
  • Avakian (three t-bins): open circle markers with error bars,
    converted from A_UL by dividing by (Df * D(y)) per your column map

Legends:
  • Left subplot: colors and their -t ranges
  • Right subplot: marker meaning (filled = Hayward, open = Avakian)

Output:
  output/enpi+/harut_cross_check.pdf
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ----------------------------
# Config
# ----------------------------
# Per your request: different y-limits per subplot
Y_L_MIN, Y_L_MAX = -0.35, 0.25   # sin(phi) panel
Y_R_MIN, Y_R_MAX = -0.60, 0.15   # sin(2phi) panel (updated to 0.15)
X_MIN, X_MAX = 0.05, 0.65
OUTPATH = "output/enpi+/harut_cross_check.pdf"

colors = ["tab:blue", "tab:orange", "tab:green"]  # three t-bin colors
tbin_labels = ["0.05 < -t < 0.45", "0.45 < -t < 0.85", "0.85 < -t < 1.225"]

# ----------------------------
# Hayward (your) results — UPDATED ARRAYS
# Format per point: (xB, value, err)
# ----------------------------
enpiHarut1GEchi2FitsAULsinphi = np.array([
    (0.169464440, -0.033040862, 0.019637808),
    (0.255432176,  0.024405972, 0.011847844),
    (0.348042430,  0.081420649, 0.009146418),
    (0.442051735,  0.099088760, 0.014005917),
    (0.537668007,  0.047496516, 0.017856319),
], dtype=float)

enpiHarut1GEchi2FitsAULsin2phi = np.array([
    (0.169464440, -0.147836872, 0.050605779),
    (0.255432176, -0.089364439, 0.019479172),
    (0.348042430, -0.102451672, 0.018080756),
    (0.442051735, -0.046612422, 0.026358805),
    (0.537668007, -0.049588114, 0.036631835),
], dtype=float)

enpiHarut2GEchi2FitsAULsinphi = np.array([
    (0.166711789, -0.185373350, 0.035154993),
    (0.257666829, -0.029528304, 0.017115088),
    (0.352791737,  0.043233876, 0.015079326),
    (0.444353621,  0.124480576, 0.021166973),
    (0.540514011,  0.064941958, 0.022453429),
], dtype=float)

enpiHarut2GEchi2FitsAULsin2phi = np.array([
    (0.166711789, -0.409780182, 0.084752034),
    (0.257666829, -0.327936196, 0.034823519),
    (0.352791737, -0.215493468, 0.027052988),
    (0.444353621, -0.161014987, 0.034770771),
    (0.540514011, -0.136942728, 0.054149666),
], dtype=float)

enpiHarut3GEchi2FitsAULsinphi = np.array([
    (0.167355144,  0.108010583, 0.044185777),
    (0.256121660,  0.104709608, 0.017636483),
    (0.349415100,  0.104212356, 0.024030534),
    (0.447242357,  0.077923070, 0.035373521),
    (0.540723607,  0.161757957, 0.030181965),
], dtype=float)

enpiHarut3GEchi2FitsAULsin2phi = np.array([
    (0.167355144, -0.258947941, 0.123175313),
    (0.256121660, -0.202247019, 0.031047992),
    (0.349415100, -0.331442135, 0.041434944),
    (0.447242357, -0.245808267, 0.051171297),
    (0.540723607, -0.186831477, 0.064951121),
], dtype=float)

hayward_sets = [
    (enpiHarut1GEchi2FitsAULsinphi,  enpiHarut1GEchi2FitsAULsin2phi),
    (enpiHarut2GEchi2FitsAULsinphi,  enpiHarut2GEchi2FitsAULsin2phi),
    (enpiHarut3GEchi2FitsAULsinphi,  enpiHarut3GEchi2FitsAULsin2phi),
]

# ----------------------------
# Avakian (colleague) results (three blocks, 5 rows each, 11 columns)
# Corrected column mapping (1-based -> 0-based):
#   1: xB                 -> 0
#   3: Df                 -> 2
#   4: AUL_sinphi         -> 3
#   5: err_sinphi         -> 4
#   8: AUL_sin2phi        -> 7
#   9: err_sin2phi        -> 8
#  10: D(y) for sinphi    -> 9
#  11: D(y) for sin2phi   -> 10
#
# Structure-function ratios:
#   F_UL^{sinphi}/F_UU  = AUL_sinphi  / (Df * Dy_sinphi)
#   F_UL^{sin2phi}/F_UU = AUL_sin2phi / (Df * Dy_sin2phi)
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

def avak_to_ratios(block):
    """
    Convert Avakian block to (xB, R_phi, R_phi_err, R_2phi, R_2phi_err, denom_phi, denom_2phi).
    Corrected index map (0-based):
      0:xB, 2:Df, 3:AUL_sinphi, 4:err_sinphi,
      7:AUL_sin2phi, 8:err_sin2phi,
      9:Dy_sinphi, 10:Dy_sin2phi
    """
    xB = block[:, 0]
    Df = block[:, 2]
    A_sinphi = block[:, 3]
    E_sinphi = np.abs(block[:, 4])
    A_sin2   = block[:, 7]
    E_sin2   = np.abs(block[:, 8])
    Dy_phi   = block[:, 9]
    Dy_2phi  = block[:, 10]

    denom_phi = Df * Dy_phi
    denom_2phi = Df * Dy_2phi

    safe_phi = np.where(denom_phi != 0.0, denom_phi, np.nan)
    safe_2phi = np.where(denom_2phi != 0.0, denom_2phi, np.nan)

    R_phi = A_sinphi / safe_phi
    R_phi_err = E_sinphi / safe_phi

    R_2phi = A_sin2 / safe_2phi
    R_2phi_err = E_sin2 / safe_2phi

    return (xB, R_phi, R_phi_err, R_2phi, R_2phi_err, denom_phi, denom_2phi)
#endfor

# Convert all three Avakian blocks
avak_blocks = [avak_block1, avak_block2, avak_block3]
avak_ratios = [avak_to_ratios(b) for b in avak_blocks]
#endfor

# ----------------------------
# Pretty printing helpers
# ----------------------------
def print_table(title, x, y, yerr, col_y_label):
    print("=" * 78)
    print(title)
    print("-" * 78)
    print("{:>3s}  {:>10s}  {:>16s}  {:>16s}".format("i", "xB", col_y_label, "err"))
    print("-" * 78)
    for i in range(len(x)):
        print("{:3d}  {:10.6f}  {:16.8f}  {:16.8f}".format(i, float(x[i]), float(y[i]), float(yerr[i])))
    #endfor
    print("")
#endfor

def print_denominator_ranges():
    for idx, (_x, _r1, _e1, _r2, _e2, dphi, d2phi) in enumerate(avak_ratios):
        tlabel = tbin_labels[idx]
        dphi_finite = dphi[np.isfinite(dphi)]
        d2phi_finite = d2phi[np.isfinite(d2phi)]
        print("Denominators (Df*Dy) ranges for {}:".format(tlabel))
        print("  sinphi : min={:.6f}, max={:.6f}".format(np.min(dphi_finite), np.max(dphi_finite)))
        print("  sin2phi: min={:.6f}, max={:.6f}".format(np.min(d2phi_finite), np.max(d2phi_finite)))
    #endfor
    print("")
#endif

def print_all_results():
    # Hayward tables (already structure-function ratios as plotted)
    for idx, (arr_phi, arr_2phi) in enumerate(hayward_sets):
        tlabel = tbin_labels[idx]
        print_table(
            f"Hayward F_UL^{{sinphi}}/F_UU   (t-bin: {tlabel})",
            arr_phi[:, 0], arr_phi[:, 1], arr_phi[:, 2],
            "FUL^{sinphi}/FUU"
        )
        print_table(
            f"Hayward F_UL^{{sin2phi}}/F_UU  (t-bin: {tlabel})",
            arr_2phi[:, 0], arr_2phi[:, 1], arr_2phi[:, 2],
            "FUL^{sin2phi}/FUU"
        )
    #endfor

    # Avakian tables (converted to structure-function ratios)
    for idx, (xB, R_phi, R_phi_err, R_2phi, R_2phi_err, _dphi, _d2phi) in enumerate(avak_ratios):
        tlabel = tbin_labels[idx]
        print_table(
            f"Avakian F_UL^{{sinphi}}/F_UU   (t-bin: {tlabel})",
            xB, R_phi, R_phi_err,
            "FUL^{sinphi}/FUU"
        )
        print_table(
            f"Avakian F_UL^{{sin2phi}}/F_UU  (t-bin: {tlabel})",
            xB, R_2phi, R_2phi_err,
            "FUL^{sin2phi}/FUU"
        )
    #endfor
#endfor

# ----------------------------
# Plot
# ----------------------------
def main():
    os.makedirs(os.path.dirname(OUTPATH), exist_ok=True)

    # Diagnostics first so you can verify sane denominators
    print_denominator_ranges()

    # Then print all numeric tables
    print_all_results()

    # No shared y so ticks/labels appear for BOTH subplots
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharex=True)
    axL, axR = axes[0], axes[1]

    # Plot Hayward sets (filled markers) with error bars; colors by t-bin
    for i, (arr_phi, arr_2phi) in enumerate(hayward_sets):
        x_phi, y_phi, e_phi = arr_phi[:, 0], arr_phi[:, 1], arr_phi[:, 2]
        x_2phi, y_2phi, e_2phi = arr_2phi[:, 0], arr_2phi[:, 1], arr_2phi[:, 2]

        # Left: sin phi
        axL.errorbar(
            x_phi, y_phi, yerr=e_phi, fmt="o", ms=6, lw=1.2, capsize=3,
            color=colors[i], mfc=colors[i], mec="black", label=tbin_labels[i]
        )

        # Right: sin 2phi
        axR.errorbar(
            x_2phi, y_2phi, yerr=e_2phi, fmt="o", ms=6, lw=1.2, capsize=3,
            color=colors[i], mfc=colors[i], mec="black"
        )
    #endfor

    # Overlay Avakian (all three t-bins) as OPEN circles in matching colors
    for i, (xB, R_phi, R_phi_err, R_2phi, R_2phi_err, _dphi, _d2phi) in enumerate(avak_ratios):
        axL.errorbar(
            xB, R_phi, yerr=R_phi_err, fmt="o", ms=6, lw=1.2, capsize=3,
            color=colors[i], mfc="none", mec=colors[i]
        )
        axR.errorbar(
            xB, R_2phi, yerr=R_2phi_err, fmt="o", ms=6, lw=1.2, capsize=3,
            color=colors[i], mfc="none", mec=colors[i]
        )
    #endfor

    # Axes labels and limits
    axL.set_xlabel(r"$x_{B}$")
    axR.set_xlabel(r"$x_{B}$")
    axL.set_ylabel(r"$F_{UL}^{\sin\phi}/F_{UU}$")
    axR.set_ylabel(r"$F_{UL}^{\sin2\phi}/F_{UU}$")
    axL.set_xlim(X_MIN, X_MAX)
    axR.set_xlim(X_MIN, X_MAX)
    axL.set_ylim(Y_L_MIN, Y_L_MAX)  # sinphi panel limits
    axR.set_ylim(Y_R_MIN, Y_R_MAX)  # sin2phi panel limits

    # Ensure y ticks & labels show for BOTH subplots
    axL.tick_params(axis="y", which="both", labelleft=True)
    axR.tick_params(axis="y", which="both", labelleft=True)

    # Legends:
    # 1) Left: t-bin color legend (already labeled for Hayward series)
    leg1 = axL.legend(loc="lower left", frameon=False, title="-t bins")

    # 2) Right: marker meaning legend
    hayward_handle = Line2D([0], [0], marker="o", linestyle="None",
                            mfc="black", mec="black", ms=6, label="Hayward")
    avakian_handle = Line2D([0], [0], marker="o", linestyle="None",
                            mfc="none", mec="black", ms=6, label="Avakian")
    leg2 = axR.legend(handles=[hayward_handle, avakian_handle],
                      loc="lower right", frameon=False, title="Markers")

    # Layout and save
    fig.tight_layout()
    fig.savefig(OUTPATH)
    print("Wrote:", OUTPATH)
#endfor

if __name__ == "__main__":
    main()
#endif