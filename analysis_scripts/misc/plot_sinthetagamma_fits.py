#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2x2 panels of F_UL^{sin(n phi)}/F_UU vs sin(theta_gamma) with uncertainties.

Top row:  low x_B  |  high x_B
Bottom:   low -t   |  high -t

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
    Weighted-mean fit (constant model).
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

    # Labels with χ²/ndf (escape LaTeX braces for Python formatting)
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

    # Legend (only two entries: the point sets)
    handles = [h for h in (h1, h2) if h is not None]
    if handles:
        ax.legend(handles=handles, frameon=True, fontsize=10,
                  title="Harmonics", title_fontsize=10, loc="best")
    else:
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                ha="center", va="center", fontsize=12, alpha=0.6)

# ─────────────────────────────────────────────────────────────────────
# Data you provided (triples: {sinθγ, value, error})
# ─────────────────────────────────────────────────────────────────────
# Low x_B
enpiLowxBGEchi2FitsAULsinphi = [
    [0.106224541, -0.028065860, 0.013605010],
    [0.156049384, -0.003632377, 0.004786986],
    [0.217968471, -0.000475566, 0.007246571],
    [0.270187215,  0.056693694, 0.021983002],
    [0.000000000,  0.000000000, 0.010000782],  # sentinel -> drop
]
enpiLowxBGEchi2FitsAULsin2phi = [
    [0.106224541,  0.015747976, 0.030275168],
    [0.156049384, -0.025927698, 0.011664101],
    [0.217968471, -0.036819170, 0.010764596],
    [0.270187215, -0.059681623, 0.023731383],
    [0.000000000,  0.000000000, 0.010001391],  # sentinel -> drop
]

# High x_B
enpiHighxBGEchi2FitsAULsinphi = [
    [0.000000000,  0.000000000, 0.010000782],  # sentinel -> drop
    [0.175687261, -0.016201462, 0.031354087],
    [0.232413857,  0.031997790, 0.009731937],
    [0.300164957,  0.040768159, 0.009676298],
    [0.359281080,  0.043694399, 0.005165615],
]
enpiHighxBGEchi2FitsAULsin2phi = [
    [0.000000000,  0.000000000, 0.010001391],  # sentinel -> drop
    [0.175687261, -0.177957475, 0.084246820],
    [0.232413857, -0.035707467, 0.025202947],
    [0.300164957, -0.078183037, 0.019622011],
    [0.359281080, -0.058298753, 0.010739176],
]

# Low -t
enpiLowtGEchi2FitsAULsinphi = [
    [0.107461948, -0.024570746, 0.030067875],
    [0.159045586, -0.008710291, 0.007323591],
    [0.224895070,  0.001769930, 0.005313844],
    [0.291905934,  0.029046547, 0.006820219],
    [0.356821364, -0.007561284, 0.008577452],
]
enpiLowtGEchi2FitsAULsin2phi = [
    [0.107461948,  0.055076753, 0.060566694],
    [0.159045586, -0.008341678, 0.016487063],
    [0.224895070, -0.018639517, 0.010753646],
    [0.291905934, -0.020406471, 0.015695833],
    [0.356821364, -0.020663126, 0.014489705],
]

# High -t  (note: name preserved as provided: "Hight")
enpiHightGEchi2FitsAULsinphi = [
    [0.105863431, -0.025144693, 0.053163791],
    [0.160573253, -0.023512803, 0.019664193],
    [0.228725090,  0.001328815, 0.020276061],
    [0.297207339,  0.027304978, 0.010367475],
    [0.359265763,  0.048359517, 0.007921233],
]
enpiHightGEchi2FitsAULsin2phi = [
    [0.105863431, -0.170791606, 0.143963579],
    [0.160573253, -0.075891112, 0.043966462],
    [0.228725090, -0.094003089, 0.025050109],
    [0.297207339, -0.059079167, 0.027004949],
    [0.359265763, -0.078686026, 0.012797081],
]

# ─────────────────────────────────────────────────────────────────────
# Prepare arrays (drop sentinel x<=0 rows)
# ─────────────────────────────────────────────────────────────────────
x_lowxb_1, y_lowxb_1, e_lowxb_1   = as_xyz(enpiLowxBGEchi2FitsAULsinphi,   drop_x_le=0.0)
x_lowxb_2, y_lowxb_2, e_lowxb_2   = as_xyz(enpiLowxBGEchi2FitsAULsin2phi,  drop_x_le=0.0)

x_highxb_1, y_highxb_1, e_highxb_1 = as_xyz(enpiHighxBGEchi2FitsAULsinphi,  drop_x_le=0.0)
x_highxb_2, y_highxb_2, e_highxb_2 = as_xyz(enpiHighxBGEchi2FitsAULsin2phi, drop_x_le=0.0)

x_lowt_1, y_lowt_1, e_lowt_1     = as_xyz(enpiLowtGEchi2FitsAULsinphi,   drop_x_le=0.0)
x_lowt_2, y_lowt_2, e_lowt_2     = as_xyz(enpiLowtGEchi2FitsAULsin2phi,  drop_x_le=0.0)

x_hight_1, y_hight_1, e_hight_1  = as_xyz(enpiHightGEchi2FitsAULsinphi,  drop_x_le=0.0)
x_hight_2, y_hight_2, e_hight_2  = as_xyz(enpiHightGEchi2FitsAULsin2phi, drop_x_le=0.0)

# ─────────────────────────────────────────────────────────────────────
# Plot
# ─────────────────────────────────────────────────────────────────────
outdir = os.path.join("output", "enpi+")
os.makedirs(outdir, exist_ok=True)

fig, axes = plt.subplots(2, 2, figsize=(12, 9))
fig.suptitle(r"$F_{UL}^{\sin(n\phi)}/F_{UU}$ vs $\sin\theta_{\gamma}$", fontsize=14, y=0.97)

# Titles include the selection ranges you provided
panel(axes[0, 0],
      x_lowxb_1, y_lowxb_1, e_lowxb_1,
      x_lowxb_2, y_lowxb_2, e_lowxb_2,
      title=r"Low $x_{B}$: $0.10<x_{B}<0.25$")

panel(axes[0, 1],
      x_highxb_1, y_highxb_1, e_highxb_1,
      x_highxb_2, y_highxb_2, e_highxb_2,
      title=r"High $x_{B}$: $0.45<x_{B}<0.60$")

panel(axes[1, 0],
      x_lowt_1, y_lowt_1, e_lowt_1,
      x_lowt_2, y_lowt_2, e_lowt_2,
      title=r"Low $-t$: $0.05<-t<0.25$")

panel(axes[1, 1],
      x_hight_1, y_hight_1, e_hight_1,
      x_hight_2, y_hight_2, e_hight_2,
      title=r"High $-t$: $0.80<-t<1.00$")

plt.tight_layout(rect=[0, 0.03, 1, 0.94], w_pad=2.0, h_pad=2.0)
outpath = os.path.join(outdir, "sinthetagamma_fits.pdf")
plt.savefig(outpath)
plt.close(fig)
print(f"Saved: {outpath}")