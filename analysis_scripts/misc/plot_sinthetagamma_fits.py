#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2x2 panels of F_UL^{sin(n phi)}/F_UU vs sin(theta_gamma) with uncertainties.

Panels (x_B, -t):
  [0,0] Low x_B & Low -t      : 0.10<x_B<0.30, 0.05<-t<0.35
  [0,1] Low x_B & High -t     : 0.10<x_B<0.30, 0.85<-t<1.35
  [1,0] High x_B & Low -t     : 0.35<x_B<0.60, 0.05<-t<0.35
  [1,1] High x_B & High -t    : 0.35<x_B<0.60, 0.85<-t<1.35

Each panel overlays two harmonics (markers only, no connecting lines):
  • n = 1  (AULsinphi)   — circles, tab:blue
  • n = 2  (AULsin2phi)  — squares, tab:orange

Also: fit each dataset to a constant (weighted mean) and plot that as a thin
dashed horizontal line in the same color. Legends include χ²/ndf of that fit.

Y range fixed to [-0.2, 0.2].
Output:  output/enpi+/sinthetagamma_fits.pdf
"""

import os
import numpy as np
import matplotlib.pyplot as plt

# ─────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────
def as_xyz(triples, drop_x_le=0.0):
    """
    Convert a list of {x, y, e} triples to clean numpy arrays.
    Drops any row with non-finite values or with x <= drop_x_le (sentinels).
    Sorts by x ascending.
    """
    if not triples:
        return np.array([]), np.array([]), np.array([])
    arr = np.asarray(triples, dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 3:
        return np.array([]), np.array([]), np.array([])
    x, y, e = arr[:,0], arr[:,1], arr[:,2]
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(e) & (x > drop_x_le)
    if not np.any(m):
        return np.array([]), np.array([]), np.array([])
    x, y, e = x[m], y[m], e[m]
    order = np.argsort(x)
    return x[order], y[order], e[order]

def fit_constant(y, e):
    """
    Weighted-mean fit (constant model):
      mu   = sum(w*y)/sum(w),  w=1/e^2
      chi2 = sum((y-mu)^2 / e^2)
      ndf  = N - 1
    Returns (mu, chi2, ndf). If no valid points, returns (nan, 0, 0).
    """
    y = np.asarray(y, float)
    e = np.asarray(e, float)
    m = np.isfinite(y) & np.isfinite(e) & (e > 0)
    n = int(np.count_nonzero(m))
    if n == 0:
        return np.nan, 0.0, 0
    w  = 1.0 / (e[m]**2)
    mu = float(np.sum(w * y[m]) / np.sum(w))
    chi2 = float(np.sum(((y[m] - mu)**2) * w))
    ndf = max(0, n - 1)
    return mu, chi2, ndf

def plot_series_points(ax, x, y, e, label, marker, color):
    """Markers + error bars only (no connecting lines)."""
    if x.size == 0:
        return None
    return ax.errorbar(
        x, y, yerr=e,
        fmt=marker, color=color, ecolor=color,
        linestyle='None',   # no line between points
        ms=4, lw=1.0, capsize=3,
        label=label
    )

def panel(ax, x1,y1,e1, x2,y2,e2, title):
    # Constant fits
    mu1, chi2_1, ndf_1 = fit_constant(y1, e1) if y1.size else (np.nan, 0.0, 0)
    mu2, chi2_2, ndf_2 = fit_constant(y2, e2) if y2.size else (np.nan, 0.0, 0)

    # Labels with χ²/ndf
    lab1 = rf"$n=1$  ($\chi^{{2}}/\mathrm{{ndf}}={chi2_1:.1f}/{ndf_1}$)"
    lab2 = rf"$n=2$  ($\chi^{{2}}/\mathrm{{ndf}}={chi2_2:.1f}/{ndf_2}$)"

    # Plot points
    h1 = plot_series_points(ax, x1, y1, e1, lab1, marker="o", color="tab:blue")
    h2 = plot_series_points(ax, x2, y2, e2, lab2, marker="s", color="tab:orange")

    # Axes cosmetics
    ax.set_title(title, fontsize=12)
    ax.set_xlabel(r"$\sin\theta_{\gamma}$")
    ax.set_ylabel(r"$F_{UL}^{\sin(n\phi)} / F_{UU}$")
    ax.set_ylim(-0.2, 0.2)

    # x-limits from available points; otherwise a sensible default
    xs = np.concatenate([
        x1 if x1.size else np.array([]),
        x2 if x2.size else np.array([])
    ])
    if xs.size:
        lo = max(0.0, float(xs.min()) - 0.01)
        hi = float(xs.max()) + 0.01
        ax.set_xlim(lo, hi)
    else:
        ax.set_xlim(0.0, 0.40)

    ax.grid(True, linestyle="--", alpha=0.35)

    # Overlay constant-fit lines (thin dashed) in same colors
    xmin, xmax = ax.get_xlim()
    if np.isfinite(mu1):
        ax.plot([xmin, xmax], [mu1, mu1], linestyle="--", linewidth=1.0, color="tab:blue")
    if np.isfinite(mu2):
        ax.plot([xmin, xmax], [mu2, mu2], linestyle="--", linewidth=1.0, color="tab:orange")

    # Legend (two entries: the point sets)
    handles = [h for h in (h1, h2) if h is not None]
    if handles:
        ax.legend(handles=handles, frameon=True, fontsize=10,
                  title="Harmonics", title_fontsize=10, loc="best")
    else:
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                ha="center", va="center", fontsize=12, alpha=0.6)

# ─────────────────────────────────────────────────────────────────────
# New data (triples: {sinθγ, value, error})
# Only AULsinphi (n=1) and AULsin2phi (n=2) are required for the plots.
# ─────────────────────────────────────────────────────────────────────

# Low x_B & Low -t   (0.10<x_B<0.30, 0.05<-t<0.35)
enpiLowxBLowtGEchi2FitsAULsinphi = [
    [0.107516772, -0.064019517, 0.027852352],
    [0.158721085, -0.015667152, 0.006461520],
    [0.223208570,  0.001289905, 0.005487116],
    [0.284301854,  0.021711468, 0.009187490],
    [0.336696510,  0.020050283, 0.037320508],
]
enpiLowxBLowtGEchi2FitsAULsin2phi = [
    [0.107516772,  0.005803660, 0.055010545],
    [0.158721085, -0.005154474, 0.018514105],
    [0.223208570, -0.023305563, 0.010715054],
    [0.284301854, -0.013637025, 0.019766420],
    [0.336696510, -0.127297142, 0.073767144],
]

# Low x_B & High -t  (0.10<x_B<0.30, 0.85<-t<1.35)
enpiLowxBHightGEchi2FitsAULsinphi = [
    [0.105650048, -0.028064870, 0.076567872],
    [0.157680814,  0.038585539, 0.017023297],
    [0.224130532,  0.026030965, 0.012459477],
    [0.284491459,  0.063985998, 0.016559254],
    [0.336825566,  0.110784828, 0.000038842],
]
enpiLowxBHightGEchi2FitsAULsin2phi = [
    [0.105650048, -0.116516175, 0.151982281],
    [0.157680814,  0.024751149, 0.032143668],
    [0.224130532, -0.090728307, 0.029384592],
    [0.284491459, -0.098747896, 0.027101219],
    [0.336825566, -0.116687033, 0.000030434],
]

# High x_B & Low -t  (0.35<x_B<0.60, 0.05<-t<0.35)
enpiHighxBBLowtGEchi2FitsAULsinphi = [
    [0.000000000,  0.000000000, 0.010000782],  # sentinel -> drop
    [0.170104546,  0.047205526, 0.035515564],
    [0.232317770,  0.022489809, 0.015343618],
    [0.299502631,  0.034118607, 0.008742240],
    [0.361806661,  0.011211490, 0.007635931],
]
enpiHighxBBLowtGEchi2FitsAULsin2phi = [
    [0.000000000,  0.000000000, 0.010001391],  # sentinel -> drop
    [0.170104546,  0.073591162, 0.086565720],
    [0.232317770,  0.009229356, 0.033590381],
    [0.299502631,  0.002870764, 0.038150067],
    [0.361806661, -0.008715496, 0.012555314],
]

# High x_B & High -t (0.35<x_B<0.60, 0.85<-t<1.35)
enpiHighxBHightGEchi2FitsAULsinphi = [
    [0.000000000,  0.000000000, 0.010000782],  # sentinel -> drop
    [0.172592573, -0.068576063, 0.046779678],
    [0.232020736,  0.038101527, 0.013602812],
    [0.299690415,  0.039984526, 0.010228255],
    [0.359996179,  0.055971602, 0.006836161],
]
enpiHighxBHightGEchi2FitsAULsin2phi = [
    [0.000000000,  0.000000000, 0.010001391],  # sentinel -> drop
    [0.172592573, -0.193949995, 0.085291091],
    [0.232020736, -0.056253835, 0.028833524],
    [0.299690415, -0.102146191, 0.015516744],
    [0.359996179, -0.071449196, 0.013017852],
]

# ─────────────────────────────────────────────────────────────────────
# Prepare arrays (drop sentinel x<=0 rows)
# ─────────────────────────────────────────────────────────────────────
# Low x_B & Low -t
x_ll_1, y_ll_1, e_ll_1 = as_xyz(enpiLowxBLowtGEchi2FitsAULsinphi,  drop_x_le=0.0)
x_ll_2, y_ll_2, e_ll_2 = as_xyz(enpiLowxBLowtGEchi2FitsAULsin2phi, drop_x_le=0.0)

# Low x_B & High -t
x_lh_1, y_lh_1, e_lh_1 = as_xyz(enpiLowxBHightGEchi2FitsAULsinphi,  drop_x_le=0.0)
x_lh_2, y_lh_2, e_lh_2 = as_xyz(enpiLowxBHightGEchi2FitsAULsin2phi, drop_x_le=0.0)

# High x_B & Low -t
x_hl_1, y_hl_1, e_hl_1 = as_xyz(enpiHighxBBLowtGEchi2FitsAULsinphi,  drop_x_le=0.0)
x_hl_2, y_hl_2, e_hl_2 = as_xyz(enpiHighxBBLowtGEchi2FitsAULsin2phi, drop_x_le=0.0)

# High x_B & High -t
x_hh_1, y_hh_1, e_hh_1 = as_xyz(enpiHighxBHightGEchi2FitsAULsinphi,  drop_x_le=0.0)
x_hh_2, y_hh_2, e_hh_2 = as_xyz(enpiHighxBHightGEchi2FitsAULsin2phi, drop_x_le=0.0)

# ─────────────────────────────────────────────────────────────────────
# Plot
# ─────────────────────────────────────────────────────────────────────
outdir = os.path.join("output", "enpi+")
os.makedirs(outdir, exist_ok=True)

fig, axes = plt.subplots(2, 2, figsize=(12, 9))
fig.suptitle(r"$F_{UL}^{\sin(n\phi)}/F_{UU}$ vs $\sin\theta_{\gamma}$", fontsize=14, y=0.97)

# [0,0] Low x_B & Low -t
panel(
    axes[0, 0],
    x_ll_1, y_ll_1, e_ll_1,
    x_ll_2, y_ll_2, e_ll_2,
    title=r"Low $x_{B}$, Low $-t$:  $0.10<x_{B}<0.30$,  $0.05<-t<0.35$"
)

# [0,1] Low x_B & High -t
panel(
    axes[0, 1],
    x_lh_1, y_lh_1, e_lh_1,
    x_lh_2, y_lh_2, e_lh_2,
    title=r"Low $x_{B}$, High $-t$: $0.10<x_{B}<0.30$,  $0.85<-t<1.35$"
)

# [1,0] High x_B & Low -t
panel(
    axes[1, 0],
    x_hl_1, y_hl_1, e_hl_1,
    x_hl_2, y_hl_2, e_hl_2,
    title=r"High $x_{B}$, Low $-t$: $0.35<x_{B}<0.60$,  $0.05<-t<0.35$"
)

# [1,1] High x_B & High -t
panel(
    axes[1, 1],
    x_hh_1, y_hh_1, e_hh_1,
    x_hh_2, y_hh_2, e_hh_2,
    title=r"High $x_{B}$, High $-t$: $0.35<x_{B}<0.60$,  $0.85<-t<1.35$"
)

plt.tight_layout(rect=[0, 0.03, 1, 0.94], w_pad=2.0, h_pad=2.0)
outpath = os.path.join(outdir, "sinthetagamma_fits.pdf")
plt.savefig(outpath)
plt.close(fig)
print(f"Saved: {outpath}")