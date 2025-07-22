#!/usr/bin/env python3
"""
plot_ImH_vs_xi.py

Reads a FitH output file (fit_results_<TIMESTAMP>.txt), extracts
the original VGG and fitted ImH parameters, then plots Im H(ξ,−t)
vs ξ in a 2×3 grid of panels (one panel per −t).

Usage:
    python plot_ImH_vs_xi.py output/fit_results/fit_results_<TIMESTAMP>.txt
"""

import os
import sys
import re
import numpy as np
import matplotlib.pyplot as plt

# ─── Parse command‐line ────────────────────────────────────────────────────────
if len(sys.argv) != 2:
    print("Usage: python plot_ImH_vs_xi.py output/fit_results/fit_results_<TIMESTAMP>.txt")
    sys.exit(1)

fitfile = sys.argv[1]
m = re.search(r'fit_results_(\d{8}_\d{6})\.txt$', fitfile)
if not m:
    print("ERROR: Couldn't extract timestamp from filename:", fitfile)
    sys.exit(1)
timestamp = m.group(1)

# ─── Load fit results ─────────────────────────────────────────────────────────
def load_fit_results(fname):
    vals = None
    chi2ndf = None
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if line.startswith("# values"):
                vals = next(f).split()
            if line.startswith("# chi2"):
                parts = next(f).split()
                chi2ndf = float(parts[2])
    if vals is None or chi2ndf is None:
        raise RuntimeError("Failed to parse fit file")
    return np.array(list(map(float, vals))), chi2ndf

vals, chi2ndf = load_fit_results(fitfile)
renorm_fit, alpha0_fit, alpha1_fit, n_fit, b_fit, Mm2_fit, P_fit, renormReal_fit = vals

# ─── Original VGG defaults ───────────────────────────────────────────────────
renorm0, alpha0_0, alpha1_0 = 1.0, 0.43, 0.85
n0, b0, Mm2_0, P0 = 1.35, 0.4, 0.64, 1.0

# ─── ImH function ────────────────────────────────────────────────────────────
def ImH(xi, t, renorm, alpha0, alpha1, n_val, b_val, Mm2_val, P_val):
    r     = 0.9
    alpha = alpha0 + alpha1 * t
    pref  = np.pi * 5.0/9.0 * n_val * r / (1 + xi)
    xfac  = (2*xi/(1+xi))**(-alpha)
    yfac  = ((1 - xi)/(1+xi))**(b_val)
    tfac  = (1 - ((1 - xi)/(1+xi))*t/Mm2_val)**(-P_val)
    return renorm * pref * xfac * yfac * tfac * 2.0

def ImH_orig(xi, t):
    return ImH(xi, t, renorm0, alpha0_0, alpha1_0, n0, b0, Mm2_0, P0)

def ImH_fit(xi, t):
    return ImH(xi, t, renorm_fit, alpha0_fit, alpha1_fit, n_fit, b_fit, Mm2_fit, P_fit)

# ─── Plot style ───────────────────────────────────────────────────────────────
plt.style.use('classic')
plt.rcParams.update({
    'font.size':       14,
    'font.family':     'serif',
    'axes.facecolor':  'white',
    'axes.edgecolor':  'black',
    'legend.frameon':  True,
    'legend.framealpha': 1.0,
})

# ─── Prepare data ─────────────────────────────────────────────────────────────
t_values = np.linspace(0.1, 1.0, 6)   # six −t values equally spaced
xi_vals  = np.linspace(0, 0.5, 300)   # ξ range up to 0.5

orig_color = 'tab:blue'
fit_color  = 'tab:red'

# ─── Create 2×3 grid ─────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 3, figsize=(12, 6),
                         sharex=True, sharey=True)
axes = axes.flatten()
fig.subplots_adjust(wspace=0, hspace=0)

for ax, tt in zip(axes, t_values):
    neg_t = -tt
    y0 = ImH_orig(xi_vals, neg_t)
    y1 = ImH_fit( xi_vals, neg_t)

    ax.plot(xi_vals, y0, '-',  lw=2, color=orig_color,
            label="Original Parameters")
    ax.plot(xi_vals, y1, '--', lw=2, color=fit_color,
            label="RGA pass-1 BSA")

    # reposition the -t label a bit left and up
    ax.text(0.60, 0.25,
            rf'$-t = {tt:.2f}\,\mathrm{{GeV}}^2$',
            transform=ax.transAxes,
            ha='left', va='bottom',
            fontsize=12)

    ax.set_xlim(0, 0.5)
    ax.set_ylim(0, 12)

    ax.legend(loc='upper right', fontsize=10)

# ─── Global axis labels ───────────────────────────────────────────────────────
fig.text(0.5, 0.02,
         r'$\xi$',
         ha='center', va='center', fontsize=16)
fig.text(0.06, 0.5,
         r'$\mathrm{Im}\,H(\xi,\,-t)$',
         ha='center', va='center',
         rotation='vertical', fontsize=16)

# ─── Save ─────────────────────────────────────────────────────────────────────
outdir = 'output/plots'
os.makedirs(outdir, exist_ok=True)
outname = f'{outdir}/ImH_vs_xi_{timestamp}.pdf'
fig.savefig(outname, bbox_inches='tight')
print("Saved figure to", outname)