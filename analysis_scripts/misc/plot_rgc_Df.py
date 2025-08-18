#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot dilution factor D_f vs -t for:
  (1) Integrated xB (single canvas)
  (2) Four xB bins on a 2x2 canvas (Low, MidLow, MidHigh, High)

Data are hard-coded from user-provided values:
- -t axis comes from the first numbers of the specified text blocks (negated, since t<0).
- D_f values come in three sets per group (three run periods); we combine them by
  inverse-variance weighting at each -t point, skipping invalid entries (NaN or σ<=0).

Outputs (PDF):
  output/enpi+/dilution_factor_integrated.pdf
  output/enpi+/dilution_factor_binned.pdf
"""

import os
import math
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# Helper: inverse-variance combiner
# -----------------------------
def combine_weighted(values, errors):
    """
    values, errors: list of lists (periods) each of length N, or None
    Returns combined (y, yerr) arrays of length N with NaNs where no valid contribs.
    We treat entries with non-finite y or non-finite/<=0 σ as invalid.
    """
    # Ensure list-of-lists
    series = [np.asarray(v, dtype=float) for v in values]
    sigmas = [np.asarray(e, dtype=float) for e in errors]
    N = min(len(s) for s in series)

    y_out  = np.full(N, np.nan)
    e_out  = np.full(N, np.nan)

    for i in range(N):
        contrib = []
        for s, se in zip(series, sigmas):
            yi, ei = s[i], se[i]
            if np.isfinite(yi) and np.isfinite(ei) and ei > 0.0:
                contrib.append((yi, ei))
        if not contrib:
            continue
        wsum = 0.0
        ywsum = 0.0
        for yi, ei in contrib:
            w = 1.0 / (ei * ei)
            wsum  += w
            ywsum += yi * w
        y_out[i] = ywsum / wsum
        e_out[i] = 1.0 / math.sqrt(wsum)
    return y_out, e_out

def sort_by_x(x, y, e):
    """Return arrays sorted by ascending x, keeping only entries with finite y,e."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    e = np.asarray(e, float)
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(e)
    x, y, e = x[mask], y[mask], e[mask]
    if len(x) == 0:
        return x, y, e
    order = np.argsort(x)
    return x[order], y[order], e[order]

def set_sensible_ylim(all_y, all_e, pad=0.05):
    """Compute a y-limit that captures all points ± errors, clamped to [0, 1.05]."""
    ymin, ymax = +np.inf, -np.inf
    for y, e in zip(all_y, all_e):
        if y is None or e is None:
            continue
        if len(y) == 0:
            continue
        ymin = min(ymin, np.nanmin(y - e))
        ymax = max(ymax, np.nanmax(y + e))
    if not np.isfinite(ymin) or not np.isfinite(ymax):
        return (0.0, 1.0)
    span = ymax - ymin
    if span <= 0:
        span = 0.1
    lo = max(0.0, ymin - pad * span)
    hi = min(1.05, ymax + pad * span)
    if hi - lo < 0.1:
        hi = min(1.05, lo + 0.1)
    return (lo, hi)

# -----------------------------
# Hard-coded -t axes (take t from lists and negate to plot -t)
# -----------------------------
# Integrated (12 points): enpiGEchi2FitsALUoffset first entries
t_integrated = np.array([
    -1.199254858, -1.098419750, -0.998794627, -0.898830973,
    -0.798462126, -0.698657933, -0.598504529, -0.498563527,
    -0.398570196, -0.298828009, -0.200246087, -0.109031497
], dtype=float)
xt_integrated = -t_integrated  # plot vs -t (GeV^2)

# Low xB (6 points): enpiLowxBGEchi2FitsALUoffset
t_low = np.array([
    -1.144558703, -0.944756340, -0.742782191,
    -0.542254429, -0.338530753, -0.133457081
], dtype=float)
xt_low = -t_low

# Mid-Low xB (6 points): enpiMidLowxBGEchi2FitsAULoffset (name typo OK: we only use t)
t_midlow = np.array([
    -1.145725849, -0.944360461, -0.743307062,
    -0.542090771, -0.337765072, -0.175366590
], dtype=float)
xt_midlow = -t_midlow

# Mid-High xB (6 points): enpiMidHighxBGEchi2FitsALUoffset
t_midhigh = np.array([
    -1.144822976, -0.944546718, -0.743696111,
    -0.542723444, -0.344232228, -0.214512059
], dtype=float)
xt_midhigh = -t_midhigh

# High xB (6 points): enpiHighxBGEchi2FitsALUoffset
t_high = np.array([
    -1.145182109, -0.945414028, -0.744747106,
    -0.549817284, -0.388625205, -0.248053003
], dtype=float)
xt_high = -t_high

# -----------------------------
# Hard-coded D_f (value, uncertainty) series for each group
#   Each group has up to 3 series (three run periods); we will combine them.
# -----------------------------
# Integrated (12)
int_vals_A = np.array([0.48661, 0.449281, 0.418568, 0.466025, 0.447367, 0.483743, 0.428328, 0.451445, 0.424836, 0.43053, 0.433139, 0.480003])
int_errs_A = np.array([0.0278322, 0.0284053, 0.0293135, 0.0218104, 0.022329, 0.0180305, 0.0215714, 0.0165681, 0.0156522, 0.0147116, 0.0130666, 0.0151448])

int_vals_B = np.array([0.476551, 0.479322, 0.463499, 0.477185, 0.474655, 0.473627, 0.463388, 0.45963, 0.463811, 0.453108, 0.444489, 0.494976])
int_errs_B = np.array([0.0111445, 0.0101659, 0.0100682, 0.00869438, 0.00822032, 0.00768402, 0.00694078, 0.0066358, 0.00583311, 0.00555336, 0.00556698, 0.00665952])

int_vals_C = np.array([0.482783, 0.478903, 0.506463, 0.475348, 0.483202, 0.477652, 0.479375, 0.464754, 0.439925, 0.440424, 0.447364, 0.49759])
int_errs_C = np.array([0.0230691, 0.0207842, 0.0172687, 0.0180106, 0.0168914, 0.0156568, 0.0136001, 0.0129023, 0.0123617, 0.011096, 0.0103078, 0.0111493])

# Low xB (6) — three period sets A/B/C
low_vals_A = np.array([0.454878, 0.418518, 0.47636, 0.308913, 0.373904, 0.460977])
low_errs_A = np.array([0.0640986, 0.0644511, 0.0458086, 0.0750217, 0.0371739, 0.0155668])

low_vals_B = np.array([0.42825, 0.455208, 0.451466, 0.437499, 0.436306, 0.467968])
low_errs_B = np.array([0.0265942, 0.0203284, 0.0190961, 0.0157999, 0.0116434, 0.00690595])

low_vals_C = np.array([0.490316, 0.454087, 0.512283, 0.373836, 0.421837, 0.46064])
low_errs_C = np.array([0.03844, 0.0391033, 0.0340766, 0.0357946, 0.0238771, 0.0122773])

# Mid-Low xB (6)
midlow_vals_A = np.array([0.463875, 0.481922, 0.453652, 0.44336, 0.440043, 0.431272])
midlow_errs_A = np.array([0.0385929, 0.0307211, 0.030964, 0.0254525, 0.0165836, 0.0152497])

midlow_vals_B = np.array([0.448325, 0.419726, 0.457637, 0.452798, 0.44676, 0.458582])
midlow_errs_B = np.array([0.0167498, 0.0153794, 0.011201, 0.00941558, 0.0068275, 0.00607505])

midlow_vals_C = np.array([0.454674, 0.477387, 0.440084, 0.483706, 0.438356, 0.472769])
midlow_errs_C = np.array([0.03273, 0.0280927, 0.0236535, 0.0168341, 0.0137165, 0.0108316])

# Mid-High xB (6)
midhigh_vals_A = np.array([0.452209, 0.393156, 0.485142, 0.457775, 0.440511, 0.499196])
midhigh_errs_A = np.array([0.0390229, 0.0343088, 0.0208589, 0.0191145, 0.0156473, 0.0228349])

midhigh_vals_B = np.array([0.478929, 0.49171, 0.469737, 0.458246, 0.468491, 0.468527])
midhigh_errs_B = np.array([0.0123799, 0.0102398, 0.00940869, 0.00780097, 0.00600854, 0.0118922])

midhigh_vals_C = np.array([0.459963, 0.484987, 0.472083, 0.479568, 0.435612, 0.473629])
midhigh_errs_C = np.array([0.0290533, 0.0207754, 0.0203631, 0.0147656, 0.0127638, 0.0211814])

# High xB (6) — two period sets provided (3rd missing)
high_vals_A = np.array([0.484279, 0.483005, 0.450656, 0.45627, 0.371169, 0.999997])
high_errs_A = np.array([0.0315924, 0.0289463, 0.0275378, 0.0265606, 0.0467666, 0.0])  # σ=0 -> ignored

high_vals_B = np.array([0.516664, 0.490963, 0.50077, 0.485811, 0.486819, 0.247805])
high_errs_B = np.array([0.012285, 0.0118878, 0.0101202, 0.009192, 0.0137284, 0.746155])

# -----------------------------
# Combine per group
# -----------------------------
# Integrated:
int_y, int_e = combine_weighted(
    [int_vals_A, int_vals_B, int_vals_C],
    [int_errs_A, int_errs_B, int_errs_C]
)

# Low:
low_y, low_e = combine_weighted(
    [low_vals_A, low_vals_B, low_vals_C],
    [low_errs_A, low_errs_B, low_errs_C]
)

# Mid-Low:
midlow_y, midlow_e = combine_weighted(
    [midlow_vals_A, midlow_vals_B, midlow_vals_C],
    [midlow_errs_A, midlow_errs_B, midlow_errs_C]
)

# Mid-High:
midhigh_y, midhigh_e = combine_weighted(
    [midhigh_vals_A, midhigh_vals_B, midhigh_vals_C],
    [midhigh_errs_A, midhigh_errs_B, midhigh_errs_C]
)

# High:
high_y, high_e = combine_weighted(
    [high_vals_A, high_vals_B],
    [high_errs_A, high_errs_B]
)

# Sort by -t (ascending), clean NaNs
xt_int,  y_int,  e_int  = sort_by_x(xt_integrated, int_y, int_e)
xt_low,  y_low,  e_low  = sort_by_x(xt_low,      low_y, low_e)
xt_mlo,  y_mlo,  e_mlo  = sort_by_x(xt_midlow,   midlow_y, midlow_e)
xt_mhi,  y_mhi,  e_mhi  = sort_by_x(xt_midhigh,  midhigh_y, midhigh_e)
xt_hi,   y_hi,   e_hi   = sort_by_x(xt_high,     high_y, high_e)

# -----------------------------
# Plotting
# -----------------------------
os.makedirs(os.path.join("output", "enpi+"), exist_ok=True)

COMMON_CUTS = r"$Q^{2}>1,\ W>2,\ y<0.75,\ 0.81<M_{x}^{2}<1.00\ \mathrm{GeV}^{2}$"
XB_LABELS = {
    "integrated": r"$0.10 < x_{B} < 0.60$",
    "low":       r"$0.10 < x_{B} < 0.25$",
    "midlow":    r"$0.25 < x_{B} < 0.35$",
    "midhigh":   r"$0.35 < x_{B} < 0.45$",
    "high":      r"$0.45 < x_{B} < 0.60$",
}
X_LABEL = r"$-t\ \mathrm{(GeV^{2})}$"
Y_LABEL = r"$D_{f}$"

# --- Canvas 1: integrated only ---
plt.figure(figsize=(7.5, 6.0))
title_int = rf"$ep \rightarrow en\pi^{{+}}$ — Dilution factor $D_{{f}}$ vs $-t$\n{XB_LABELS['integrated']}; {COMMON_CUTS}"
plt.suptitle(title_int, y=0.96, fontsize=14)
ax = plt.gca()

ax.errorbar(xt_int, y_int, yerr=e_int, fmt='o', color='black', ecolor='black', capsize=3, label="Combined")

ax.set_xlabel(X_LABEL)
ax.set_ylabel(Y_LABEL)
ax.grid(True, linestyle="--", alpha=0.6)
ax.set_xlim(0.0, 1.30)

ylim = set_sensible_ylim([y_int], [e_int])
ax.set_ylim(*ylim)

out1 = os.path.join("output", "enpi+", "dilution_factor_integrated.pdf")
plt.tight_layout(rect=[0, 0, 1, 0.93])
plt.savefig(out1)
plt.close()
print(f"Saved: {out1}")

# --- Canvas 2: 2x2 with four xB bins ---
fig, axes = plt.subplots(2, 2, figsize=(12, 9))
fig.suptitle(rf"$ep \rightarrow en\pi^{{+}}$ — Dilution factor $D_{{f}}$ vs $-t$\n{COMMON_CUTS}", y=0.98, fontsize=16)

panels = [
    (axes[0,0], "low",    xt_low, y_low, e_low),
    (axes[0,1], "midlow", xt_mlo, y_mlo, e_mlo),
    (axes[1,0], "midhigh",xt_mhi, y_mhi, e_mhi),
    (axes[1,1], "high",   xt_hi,  y_hi,  e_hi),
]

# Determine a common y-limit that fits everything
all_y = [p[3] for p in panels]
all_e = [p[4] for p in panels]
common_ylim = set_sensible_ylim(all_y, all_e)

for ax, key, x, y, e in panels:
    ax.errorbar(x, y, yerr=e, fmt='o', color='black', ecolor='black', capsize=3)
    ax.set_title(XB_LABELS[key], fontsize=13)
    ax.set_xlabel(X_LABEL)
    ax.set_ylabel(Y_LABEL)
    ax.grid(True, linestyle="--", alpha=0.6)
    ax.set_xlim(0.0, 1.30)
    ax.set_ylim(*common_ylim)

plt.tight_layout(rect=[0, 0, 1, 0.95], w_pad=2.0, h_pad=2.0)
out2 = os.path.join("output", "enpi+", "dilution_factor_binned.pdf")
plt.savefig(out2)
plt.close()
print(f"Saved: {out2}")