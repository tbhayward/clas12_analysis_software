#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
harut_cross_check.py

Make a 1x2 figure:
  Left:  F_LU^{sin phi} / F_UU  vs x_B
  Right: F_LU^{sin 2phi} / F_UU vs x_B

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
Y_MIN, Y_MAX = -0.60, 0.25
X_MIN, X_MAX = 0.05, 0.65  
OUTPATH = "output/enpi+/harut_cross_check.pdf"

colors = ["tab:blue", "tab:orange", "tab:green"]  # three t-bin colors
tbin_labels = ["0.05 < -t < 0.45", "0.45 < -t < 0.85", "0.85 < -t < 1.225"]

# ----------------------------
# Hayward (your) results
# Format per point: (xB, value, err)
# ----------------------------
enpiHarut1GEchi2FitsAULsinphi = np.array([
    (0.148219132, -0.018428140, 0.010281730),
    (0.238723759,  0.013108866, 0.004384448),
    (0.345392328,  0.028923008, 0.004619240),
    (0.443749834,  0.040931397, 0.006674216),
    (0.544170246, -0.006558865, 0.015482225),
], dtype=float)

enpiHarut1GEchi2FitsAULsin2phi = np.array([
    (0.148219132, -0.059938776, 0.032469481),
    (0.238723759, -0.024996669, 0.007703952),
    (0.345392328, -0.028058999, 0.007292368),
    (0.443749834, -0.025447428, 0.015817982),
    (0.544170246, -0.014103236, 0.072110831),
], dtype=float)

enpiHarut2GEchi2FitsAULsinphi = np.array([
    (0.145303000, -0.076281939, 0.016150381),
    (0.243772763, -0.027950448, 0.006584280),
    (0.355394572,  0.017372853, 0.004129027),
    (0.461519734,  0.030582373, 0.004596055),
    (0.557834047,  0.008527575, 0.009339692),
], dtype=float)

enpiHarut2GEchi2FitsAULsin2phi = np.array([
    (0.145303000,  0.102963164, 0.051422033),
    (0.243772763, -0.137232106, 0.015075803),
    (0.355394572, -0.062615979, 0.009788863),
    (0.461519734, -0.040056204, 0.009909883),
    (0.557834047, -0.018670646, 0.019436705),
], dtype=float)

enpiHarut3GEchi2FitsAULsinphi = np.array([
    (0.146069979,  0.039829520, 0.021220302),
    (0.241307119,  0.004658975, 0.009210669),
    (0.351589535,  0.018260599, 0.008218727),
    (0.466428246,  0.045420074, 0.011250287),
    (0.571300813,  0.024085029, 0.009393540),
], dtype=float)

enpiHarut3GEchi2FitsAULsin2phi = np.array([
    (0.146069979, -0.162472154, 0.062153993),
    (0.241307119, -0.074837665, 0.023728363),
    (0.351589535, -0.131269298, 0.013533272),
    (0.466428246, -0.088536541, 0.013662723),
    (0.571300813, -0.031148374, 0.020804943),
], dtype=float)

hayward_sets = [
    (enpiHarut1GEchi2FitsAULsinphi,  enpiHarut1GEchi2FitsAULsin2phi),
    (enpiHarut2GEchi2FitsAULsinphi,  enpiHarut2GEchi2FitsAULsin2phi),
    (enpiHarut3GEchi2FitsAULsinphi,  enpiHarut3GEchi2FitsAULsin2phi),
]

# ----------------------------
# Avakian (colleague) results (three blocks, 5 rows each, 11 columns)
# Column mapping (1-based indices you provided):
#   1: xB
#   3: Df
#   4: AULsinphi
#   5: err_sinphi
#   7: AULsin2phi
#   8: err_sin2phi
#   9: D(y) for sinphi
#  10: D(y) for sin2phi
#
# We compute structure-function ratios as:
#   F_LU^{sinphi}/F_UU  = AULsinphi  / (Df * D_y_sinphi)
#   F_LU^{sin2phi}/F_UU = AULsin2phi / (Df * D_y_sin2phi)
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
    Convert Avakian block to (xB, R_phi, R_phi_err, R_2phi, R_2phi_err).
    Index map (0-based):
      0:xB, 2:Df, 3:AUL_sinphi, 4:err_sinphi, 6:AUL_sin2phi, 7:err_sin2phi,
      8:Dy_sinphi, 9:Dy_sin2phi
    """
    xB = block[:, 0]
    Df = block[:, 2]
    A_sinphi = block[:, 3]
    E_sinphi = np.abs(block[:, 4])
    A_sin2   = block[:, 6]
    E_sin2   = np.abs(block[:, 7])
    Dy_phi   = block[:, 8]
    Dy_2phi  = block[:, 9]

    denom_phi = Df * Dy_phi
    denom_2phi = Df * Dy_2phi

    # Avoid division by zero
    safe_phi = np.where(denom_phi != 0.0, denom_phi, np.nan)
    safe_2phi = np.where(denom_2phi != 0.0, denom_2phi, np.nan)

    R_phi = A_sinphi / safe_phi
    R_phi_err = E_sinphi / safe_phi

    R_2phi = A_sin2 / safe_2phi
    R_2phi_err = E_sin2 / safe_2phi

    return (xB, R_phi, R_phi_err, R_2phi, R_2phi_err)
#endfor

# Convert all three Avakian blocks
avak_blocks = [avak_block1, avak_block2, avak_block3]
avak_ratios = [avak_to_ratios(b) for b in avak_blocks]  # list of tuples
#endfor

# ----------------------------
# Plot
# ----------------------------
def main():
    os.makedirs(os.path.dirname(OUTPATH), exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
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
    for i, (xB, R_phi, R_phi_err, R_2phi, R_2phi_err) in enumerate(avak_ratios):
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
    axL.set_ylim(Y_MIN, Y_MAX)
    axR.set_ylim(Y_MIN, Y_MAX)

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