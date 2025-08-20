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

# Removed explicit Mx^2 window from the title text
COMMON_CUTS = r"$Q^{2}>1,\ W>2,\ y<0.75$"
XB_LABELS = {
    "integrated": r"$0.10 < x_{B} < 0.60$",
    "low":        r"$0.10 < x_{B} < 0.25$",
    "midlow":     r"$0.25 < x_{B} < 0.35$",
    "midhigh":    r"$0.35 < x_{B} < 0.45$",
    "high":       r"$0.45 < x_{B} < 0.60$",
}

# ─────────────────────────────────────────────────────────────────────
# -t points for each grouping (positive -t, derived from provided t<0)
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
# NEW dilution factor series (Su22, Fa22, Sp23), per group
# ─────────────────────────────────────────────────────────────────────
# ---- Su22 ----
Su22_int_val = np.array([0.436756, 0.389743, 0.375151, 0.416042, 0.38615, 0.432918,
                         0.386351, 0.399457, 0.363955, 0.38188, 0.403243, 0.40839])
Su22_int_err = np.array([0.0252046, 0.0257974, 0.0252794, 0.0193817, 0.0216881, 0.0165798,
                         0.0187477, 0.0150859, 0.0145927, 0.0129786, 0.0113738, 0.01624])

Su22_low_val = np.array([0.371153, 0.411246, 0.442005, 0.24115, 0.345906, 0.386037])
Su22_low_err = np.array([0.0546634, 0.048527, 0.0434388, 0.0646981, 0.0298382, 0.0160896])

Su22_midlow_val = np.array([0.401264, 0.406516, 0.413401, 0.395165, 0.365751, 0.40021])
Su22_midlow_err = np.array([0.0336909, 0.027903, 0.0293736, 0.0228527, 0.0158932, 0.0140583])

Su22_midhigh_val = np.array([0.423179, 0.36488, 0.413048, 0.407872, 0.394538, 0.465601])
Su22_midhigh_err = np.array([0.0365177, 0.0306887, 0.0198313, 0.0167711, 0.0142019, 0.0186615])

Su22_high_val = np.array([0.42464, 0.417366, 0.399337, 0.416074, 0.341447, 0.602045])
Su22_high_err = np.array([0.029742, 0.0264294, 0.0258747, 0.0241789, 0.0386353, 0.4432])

# ---- Fa22 ----
Fa22_int_val = np.array([0.430085, 0.438061, 0.421113, 0.419502, 0.424511, 0.431267,
                         0.413241, 0.415529, 0.402449, 0.400425, 0.401987, 0.444058])
Fa22_int_err = np.array([0.00998986, 0.00918585, 0.00872031, 0.00801108, 0.00721993, 0.00663911,
                         0.00631867, 0.00586152, 0.00524712, 0.00489041, 0.00492321, 0.00602917])

Fa22_low_val = np.array([0.387409, 0.418177, 0.420803, 0.396076, 0.393667, 0.411351])
Fa22_low_err = np.array([0.0211724, 0.0162096, 0.0168978, 0.0138455, 0.0101413, 0.00596732])

Fa22_midlow_val = np.array([0.403978, 0.367673, 0.427534, 0.419153, 0.383641, 0.423595])
Fa22_midlow_err = np.array([0.0142024, 0.0132877, 0.00987297, 0.00831253, 0.00594896, 0.00549946])

Fa22_midhigh_val = np.array([0.461051, 0.448859, 0.421039, 0.409896, 0.418673, 0.415829])
Fa22_midhigh_err = np.array([0.0116705, 0.00960309, 0.00781199, 0.00682293, 0.00537743, 0.0114879])

Fa22_high_val = np.array([0.453268, 0.432248, 0.439673, 0.425048, 0.407705, 0.0568268])
Fa22_high_err = np.array([0.011452, 0.010685, 0.00920494, 0.00877271, 0.0132452, 0.512669])

# ---- Sp23 ----
# (Provided values match Fa22 block in your message)
Sp23_int_val = np.array([0.430085, 0.438061, 0.421113, 0.419502, 0.424511, 0.431267,
                         0.413241, 0.415529, 0.402449, 0.400425, 0.401987, 0.444058])
Sp23_int_err = np.array([0.00998986, 0.00918585, 0.00872031, 0.00801108, 0.00721993, 0.00663911,
                         0.00631867, 0.00586152, 0.00524712, 0.00489041, 0.00492321, 0.00602917])

Sp23_low_val = np.array([0.387409, 0.418177, 0.420803, 0.396076, 0.393667, 0.411351])
Sp23_low_err = np.array([0.0211724, 0.0162096, 0.0168978, 0.0138455, 0.0101413, 0.00596732])

Sp23_midlow_val = np.array([0.403978, 0.367673, 0.427534, 0.419153, 0.383641, 0.423595])
Sp23_midlow_err = np.array([0.0142024, 0.0132877, 0.00987297, 0.00831253, 0.00594896, 0.00549946])

Sp23_midhigh_val = np.array([0.461051, 0.448859, 0.421039, 0.409896, 0.418673, 0.415829])
Sp23_midhigh_err = np.array([0.0116705, 0.00960309, 0.00781199, 0.00682293, 0.00537743, 0.0114879])

Sp23_high_val = np.array([0.453268, 0.432248, 0.439673, 0.425048, 0.407705, 0.0568268])
Sp23_high_err = np.array([0.011452, 0.010685, 0.00920494, 0.00877271, 0.0132452, 0.512669])

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
    ax.errorbar(xI_Su22, yI_Su22, yerr=eI_Su22, fmt=MARKER, color=COLORS["Su22"],
                ecolor=COLORS["Su22"], capsize=CAPSIZE, label="Su22")
if yI_Fa22.size:
    ax.errorbar(xI_Fa22, yI_Fa22, yerr=eI_Fa22, fmt=MARKER, color=COLORS["Fa22"],
                ecolor=COLORS["Fa22"], capsize=CAPSIZE, label="Fa22")
if yI_Sp23.size:
    ax.errorbar(xI_Sp23, yI_Sp23, yerr=eI_Sp23, fmt=MARKER, color=COLORS["Sp23"],
                ecolor=COLORS["Sp23"], capsize=CAPSIZE, label="Sp23")

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
        ax.errorbar(xS, yS, yerr=eS, fmt=MARKER, color=COLORS["Su22"], ecolor=COLORS["Su22"], capsize=CAPSIZE)
        h.append(Line2D([0],[0], marker=MARKER, color=COLORS["Su22"], linestyle='', label="Su22"))
    if yF.size:
        ax.errorbar(xF, yF, yerr=eF, fmt=MARKER, color=COLORS["Fa22"], ecolor=COLORS["Fa22"], capsize=CAPSIZE)
        h.append(Line2D([0],[0], marker=MARKER, color=COLORS["Fa22"], linestyle='', label="Fa22"))
    if yP.size:
        ax.errorbar(xP, yP, yerr=eP, fmt=MARKER, color=COLORS["Sp23"], ecolor=COLORS["Sp23"], capsize=CAPSIZE)
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