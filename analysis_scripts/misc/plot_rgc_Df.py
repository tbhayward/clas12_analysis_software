#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Make two PDFs of dilution factor D_f vs -t for ep -> en pi+:
  1) Integrated xB (single axes) with Su22, Fa22, Sp23 plotted separately
  2) 2x2 canvas with Low/MidLow/MidHigh/High xB bins (each subplot shows all periods)

Saves:
  output/enpi+/dilution_factor_integrated.pdf
  output/enpi+/dilution_factor_binned.pdf
"""

import os
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.lines import Line2D

# ─────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────
def sort_and_mask(x, y, e):
    """Sort by x ascending and drop points with non-finite y or e<=0."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    e = np.asarray(e, float)
    n = min(len(x), len(y), len(e))
    x, y, e = x[:n], y[:n], e[:n]
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(e) & (e > 0)
    if not np.any(mask):
        return np.array([]), np.array([]), np.array([])
    x, y, e = x[mask], y[mask], e[mask]
    order = np.argsort(x)
    return x[order], y[order], e[order]

# ─────────────────────────────────────────────────────────────────────
# Styling
# ─────────────────────────────────────────────────────────────────────
COLORS = {"Su22": "tab:blue", "Fa22": "tab:orange", "Sp23": "tab:green"}
MARKER = "o"
CAPSIZE = 3
X_LABEL = r"$-t\ \mathrm{(GeV^{2})}$"
Y_LABEL = r"$D_{f}$"
YFIX = (0.2, 0.6)  # fixed y-limits everywhere

COMMON_CUTS = r"$Q^{2}>1,\ W>2,\ y<0.75,\ 0.81<M_{x}^{2}<1.00\ \mathrm{GeV}^{2}$"
XB_LABELS = {
    "integrated": r"$0.10 < x_{B} < 0.60$",
    "low":        r"$0.10 < x_{B} < 0.25$",
    "midlow":     r"$0.25 < x_{B} < 0.35$",
    "midhigh":    r"$0.35 < x_{B} < 0.45$",
    "high":       r"$0.45 < x_{B} < 0.60$",
}

# ─────────────────────────────────────────────────────────────────────
# Hard-coded -t points (negating the provided t<0 means these are positive)
# ─────────────────────────────────────────────────────────────────────
# Integrated (12)
t_integrated = np.array([
    -1.199254858, -1.098419750, -0.998794627, -0.898830973,
    -0.798462126, -0.698657933, -0.598504529, -0.498563527,
    -0.398570196, -0.298828009, -0.200246087, -0.109031497
], dtype=float)
xt_integrated = -t_integrated

# Low xB (6)
t_low = np.array([
    -1.144558703, -0.944756340, -0.742782191,
    -0.542254429, -0.338530753, -0.133457081
], dtype=float)
xt_low = -t_low

# Mid-Low xB (6)
t_midlow = np.array([
    -1.145725849, -0.944360461, -0.743307062,
    -0.542090771, -0.337765072, -0.175366590
], dtype=float)
xt_midlow = -t_midlow

# Mid-High xB (6)
t_midhigh = np.array([
    -1.144822976, -0.944546718, -0.743696111,
    -0.542723444, -0.344232228, -0.214512059
], dtype=float)
xt_midhigh = -t_midhigh

# High xB (6)
t_high = np.array([
    -1.145182109, -0.945414028, -0.744747106,
    -0.549817284, -0.388625205, -0.248053003
], dtype=float)
xt_high = -t_high

# ─────────────────────────────────────────────────────────────────────
# Hard-coded dilution factor series (Su22, Fa22, Sp23), per group
# ─────────────────────────────────────────────────────────────────────
# Integrated (12 each)
Su22_int_val = np.array([0.48661, 0.449281, 0.418568, 0.466025, 0.447367, 0.483743, 0.428328, 0.451445, 0.424836, 0.43053, 0.433139, 0.480003])
Su22_int_err = np.array([0.0278322, 0.0284053, 0.0293135, 0.0218104, 0.022329, 0.0180305, 0.0215714, 0.0165681, 0.0156522, 0.0147116, 0.0130666, 0.0151448])

Fa22_int_val = np.array([0.476551, 0.479322, 0.463499, 0.477185, 0.474655, 0.473627, 0.463388, 0.45963, 0.463811, 0.453108, 0.444489, 0.494976])
Fa22_int_err = np.array([0.0111445, 0.0101659, 0.0100682, 0.00869438, 0.00822032, 0.00768402, 0.00694078, 0.0066358, 0.00583311, 0.00555336, 0.00556698, 0.00665952])

Sp23_int_val = np.array([0.482783, 0.478903, 0.506463, 0.475348, 0.483202, 0.477652, 0.479375, 0.464754, 0.439925, 0.440424, 0.447364, 0.49759])
Sp23_int_err = np.array([0.0230691, 0.0207842, 0.0172687, 0.0180106, 0.0168914, 0.0156568, 0.0136001, 0.0129023, 0.0123617, 0.011096, 0.0103078, 0.0111493])

# Low xB (6 each)
Su22_low_val = np.array([0.454878, 0.418518, 0.47636, 0.308913, 0.373904, 0.460977])
Su22_low_err = np.array([0.0640986, 0.0644511, 0.0458086, 0.0750217, 0.0371739, 0.0155668])

Fa22_low_val = np.array([0.42825, 0.455208, 0.451466, 0.437499, 0.436306, 0.467968])
Fa22_low_err = np.array([0.0265942, 0.0203284, 0.0190961, 0.0157999, 0.0116434, 0.00690595])

Sp23_low_val = np.array([0.490316, 0.454087, 0.512283, 0.373836, 0.421837, 0.46064])
Sp23_low_err = np.array([0.03844, 0.0391033, 0.0340766, 0.0357946, 0.0238771, 0.0122773])

# Mid-Low xB (6 each)
Su22_midlow_val = np.array([0.463875, 0.481922, 0.453652, 0.44336, 0.440043, 0.431272])
Su22_midlow_err = np.array([0.0385929, 0.0307211, 0.030964, 0.0254525, 0.0165836, 0.0152497])

Fa22_midlow_val = np.array([0.448325, 0.419726, 0.457637, 0.452798, 0.44676, 0.458582])
Fa22_midlow_err = np.array([0.0167498, 0.0153794, 0.011201, 0.00941558, 0.0068275, 0.00607505])

Sp23_midlow_val = np.array([0.454674, 0.477387, 0.440084, 0.483706, 0.438356, 0.472769])
Sp23_midlow_err = np.array([0.03273, 0.0280927, 0.0236535, 0.0168341, 0.0137165, 0.0108316])

# Mid-High xB (6 each)
Su22_midhigh_val = np.array([0.452209, 0.393156, 0.485142, 0.457775, 0.440511, 0.499196])
Su22_midhigh_err = np.array([0.0390229, 0.0343088, 0.0208589, 0.0191145, 0.0156473, 0.0228349])

Fa22_midhigh_val = np.array([0.478929, 0.49171, 0.469737, 0.458246, 0.468491, 0.468527])
Fa22_midhigh_err = np.array([0.0123799, 0.0102398, 0.00940869, 0.00780097, 0.00600854, 0.0118922])

Sp23_midhigh_val = np.array([0.459963, 0.484987, 0.472083, 0.479568, 0.435612, 0.473629])
Sp23_midhigh_err = np.array([0.0290533, 0.0207754, 0.0203631, 0.0147656, 0.0127638, 0.0211814])

# High xB (6 each)
Su22_high_val = np.array([0.484279, 0.483005, 0.450656, 0.45627, 0.371169, 0.999997])
Su22_high_err = np.array([0.0315924, 0.0289463, 0.0275378, 0.0265606, 0.0467666, 0.0])  # σ=0 -> masked out

Fa22_high_val = np.array([0.516664, 0.490963, 0.50077, 0.485811, 0.486819, 0.247805])
Fa22_high_err = np.array([0.012285, 0.0118878, 0.0101202, 0.009192, 0.0137284, 0.746155])

Sp23_high_val = np.array([0.514416, 0.517251, 0.511247, 0.48416, 0.500329, -0.0114814])
Sp23_high_err = np.array([0.0263089, 0.0213761, 0.0197661, 0.0195764, 0.0258954, 0.21053])

# ─────────────────────────────────────────────────────────────────────
# Prepare series (each period has its own x after masking)
# ─────────────────────────────────────────────────────────────────────
# Integrated
xI_Su22, yI_Su22, eI_Su22 = sort_and_mask(xt_integrated, Su22_int_val, Su22_int_err)
xI_Fa22, yI_Fa22, eI_Fa22 = sort_and_mask(xt_integrated, Fa22_int_val, Fa22_int_err)
xI_Sp23, yI_Sp23, eI_Sp23 = sort_and_mask(xt_integrated, Sp23_int_val, Sp23_int_err)

# Low
xL_Su22, yL_Su22, eL_Su22 = sort_and_mask(xt_low, Su22_low_val, Su22_low_err)
xL_Fa22, yL_Fa22, eL_Fa22 = sort_and_mask(xt_low, Fa22_low_val, Fa22_low_err)
xL_Sp23, yL_Sp23, eL_Sp23 = sort_and_mask(xt_low, Sp23_low_val, Sp23_low_err)

# Mid-Low
xMLo_Su22, yMLo_Su22, eMLo_Su22 = sort_and_mask(xt_midlow, Su22_midlow_val, Su22_midlow_err)
xMLo_Fa22, yMLo_Fa22, eMLo_Fa22 = sort_and_mask(xt_midlow, Fa22_midlow_val, Fa22_midlow_err)
xMLo_Sp23, yMLo_Sp23, eMLo_Sp23 = sort_and_mask(xt_midlow, Sp23_midlow_val, Sp23_midlow_err)

# Mid-High
xMHi_Su22, yMHi_Su22, eMHi_Su22 = sort_and_mask(xt_midhigh, Su22_midhigh_val, Su22_midhigh_err)
xMHi_Fa22, yMHi_Fa22, eMHi_Fa22 = sort_and_mask(xt_midhigh, Fa22_midhigh_val, Fa22_midhigh_err)
xMHi_Sp23, yMHi_Sp23, eMHi_Sp23 = sort_and_mask(xt_midhigh, Sp23_midhigh_val, Sp23_midhigh_err)

# High
xH_Su22, yH_Su22, eH_Su22 = sort_and_mask(xt_high, Su22_high_val, Su22_high_err)
xH_Fa22, yH_Fa22, eH_Fa22 = sort_and_mask(xt_high, Fa22_high_val, Fa22_high_err)
xH_Sp23, yH_Sp23, eH_Sp23 = sort_and_mask(xt_high, Sp23_high_val, Sp23_high_err)

# ─────────────────────────────────────────────────────────────────────
# Plotting
# ─────────────────────────────────────────────────────────────────────
outdir = os.path.join("output", "enpi+")
os.makedirs(outdir, exist_ok=True)

# 1) Integrated: one axes, 3 periods (each with its own x)
plt.figure(figsize=(7.5, 6.0))
title_int = rf"$ep \rightarrow en\pi^{{+}}$ — {XB_LABELS['integrated']}, {COMMON_CUTS}"
plt.suptitle(title_int, y=0.95, fontsize=13)
ax = plt.gca()

if yI_Su22.size:
    ax.errorbar(xI_Su22, yI_Su22, yerr=eI_Su22, fmt=MARKER, color=COLORS["Su22"], ecolor=COLORS["Su22"], capsize=CAPSIZE, label="Su22")
if yI_Fa22.size:
    ax.errorbar(xI_Fa22, yI_Fa22, yerr=eI_Fa22, fmt=MARKER, color=COLORS["Fa22"], ecolor=COLORS["Fa22"], capsize=CAPSIZE, label="Fa22")
if yI_Sp23.size:
    ax.errorbar(xI_Sp23, yI_Sp23, yerr=eI_Sp23, fmt=MARKER, color=COLORS["Sp23"], ecolor=COLORS["Sp23"], capsize=CAPSIZE, label="Sp23")

ax.set_xlabel(X_LABEL)
ax.set_ylabel(Y_LABEL)
ax.grid(True, linestyle="--", alpha=0.6)
ax.set_xlim(0.0, 1.30)
ax.set_ylim(*YFIX)
ax.legend(title="Run Period")

out1 = os.path.join(outdir, "dilution_factor_integrated.pdf")
plt.tight_layout(rect=[0, 0, 1, 0.92])
plt.savefig(out1)
plt.close()
print(f"Saved: {out1}")

# 2) 2x2 canvas: four xB bins; each subplot: Su22, Fa22, Sp23 with its own legend
fig, axes = plt.subplots(2, 2, figsize=(12, 9))
fig.suptitle(rf"$ep \rightarrow en\pi^{{+}}$, {COMMON_CUTS}", y=0.985, fontsize=13)

def plot_panel(ax, xS,yS,eS, xF,yF,eF, xP,yP,eP, title, show_xlabel=True, show_ylabel=True):
    h = []
    if yS.size:
        h1 = ax.errorbar(xS, yS, yerr=eS, fmt=MARKER, color=COLORS["Su22"], ecolor=COLORS["Su22"], capsize=CAPSIZE, label="Su22")
        h.append(Line2D([0],[0], marker=MARKER, color=COLORS["Su22"], linestyle='', label="Su22"))
    if yF.size:
        h2 = ax.errorbar(xF, yF, yerr=eF, fmt=MARKER, color=COLORS["Fa22"], ecolor=COLORS["Fa22"], capsize=CAPSIZE, label="Fa22")
        h.append(Line2D([0],[0], marker=MARKER, color=COLORS["Fa22"], linestyle='', label="Fa22"))
    if yP.size:
        h3 = ax.errorbar(xP, yP, yerr=eP, fmt=MARKER, color=COLORS["Sp23"], ecolor=COLORS["Sp23"], capsize=CAPSIZE, label="Sp23")
        h.append(Line2D([0],[0], marker=MARKER, color=COLORS["Sp23"], linestyle='', label="Sp23"))
    ax.set_title(title, fontsize=12)
    ax.set_xlim(0.0, 1.30)
    ax.set_ylim(*YFIX)
    ax.grid(True, linestyle="--", alpha=0.6)
    ax.set_xlabel(X_LABEL if show_xlabel else "")
    ax.set_ylabel(Y_LABEL if show_ylabel else "")
    if h:
        ax.legend(handles=h, title="Run Period", loc="best", fontsize=10, title_fontsize=10, frameon=True)

# Top row: no x-labels; right column: no y-labels
plot_panel(axes[0,0],
           xL_Su22,yL_Su22,eL_Su22, xL_Fa22,yL_Fa22,eL_Fa22, xL_Sp23,yL_Sp23,eL_Sp23,
           XB_LABELS["low"], show_xlabel=False, show_ylabel=True)

plot_panel(axes[0,1],
           xMLo_Su22,yMLo_Su22,eMLo_Su22, xMLo_Fa22,yMLo_Fa22,eMLo_Fa22, xMLo_Sp23,yMLo_Sp23,eMLo_Sp23,
           XB_LABELS["midlow"], show_xlabel=False, show_ylabel=False)

plot_panel(axes[1,0],
           xMHi_Su22,yMHi_Su22,eMHi_Su22, xMHi_Fa22,yMHi_Fa22,eMHi_Fa22, xMHi_Sp23,yMHi_Sp23,eMHi_Sp23,
           XB_LABELS["midhigh"], show_xlabel=True, show_ylabel=True)

plot_panel(axes[1,1],
           xH_Su22,yH_Su22,eH_Su22, xH_Fa22,yH_Fa22,eH_Fa22, xH_Sp23,yH_Sp23,eH_Sp23,
           XB_LABELS["high"], show_xlabel=True, show_ylabel=False)

plt.tight_layout(rect=[0, 0.03, 1, 0.95], w_pad=2.0, h_pad=2.0)
out2 = os.path.join(outdir, "dilution_factor_binned.pdf")
plt.savefig(out2)
plt.close()
print(f"Saved: {out2}")