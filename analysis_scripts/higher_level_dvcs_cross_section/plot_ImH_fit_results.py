#!/usr/bin/env python3
"""
plot_ImH_vs_xi.py

Reads a FitH output file (fit_results_<TIMESTAMP>.txt), extracts
the original VGG and fitted ImH parameters, then plots Im H(ξ,−t)
vs ξ in a 2×2 grid of panels (one panel per −t).

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
                # parts: [chi2, ndf, chi2/ndf]
                chi2ndf = float(parts[2])
    if vals is None or chi2ndf is None:
        raise RuntimeError("Failed to parse fit file")
    return np.array(list(map(float, vals))), chi2ndf

vals, chi2ndf = load_fit_results(fitfile)
renorm_fit, alpha0_fit, alpha1_fit, n_fit, b_fit, Mm2_fit, P_fit, renormReal_fit = vals

# ─── Original VGG defaults ───────────────────────────────────────────────────
renorm0, alpha0_0, alpha1_0 = 1.0, 0.43, 0.85
n0,   b0,    Mm2_0,   P0    = 1.35, 0.4, 0.64, 1.0

# ─── ImH function ────────────────────────────────────────────────────────────
def ImH(xi, t, renorm, alpha0, alpha1, n_val, b_val, Mm2_val, P_val):
    r     = 0.9
    alpha = alpha0 + alpha1 * t
    n, b  = n_val, b_val
    Mm2, P= Mm2_val, P_val
    pref  = np.pi * 5.0/9.0 * n * r / (1 + xi)
    xfac  = (2*xi/(1+xi))**(-alpha)
    yfac  = ((1 - xi)/(1+xi))**(b)
    tfac  = (1 - ((1 - xi)/(1+xi))*t/Mm2)**(-P)
    return renorm * pref * xfac * yfac * tfac * 2.0

def ImH_orig(xi, t):
    return ImH(xi, t, renorm0, alpha0_0, alpha1_0, n0, b0, Mm2_0, P0)

def ImH_fit(xi, t):
    return ImH(xi, t, renorm_fit, alpha0_fit, alpha1_fit, n_fit, b_fit, Mm2_fit, P_fit)

# ─── Plot style ───────────────────────────────────────────────────────────────
plt.style.use('classic')
plt.rcParams.update({
    'font.size':        14,
    'font.family':      'serif',
    'axes.facecolor':   'white',
    'axes.edgecolor':   'black',
    'axes.grid':        True,
    'grid.linestyle':   '--',
    'grid.color':       '0.8',
    'legend.frameon':   True,
    'legend.framealpha':1.0,
})

# ─── Prepare data ─────────────────────────────────────────────────────────────
t_values = [0.1, 0.4, 0.7, 1.0]         # −t in GeV²
xi_vals  = np.linspace(0.02, 0.5, 300)   # ξ-range up to 0.5

# ─── Create 2×2 grid ─────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True, sharey=True)
axes = axes.flatten()

for ax, tt in zip(axes, t_values):
    neg_t = -tt
    y0 = ImH_orig(xi_vals, neg_t)
    y1 = ImH_fit( xi_vals, neg_t)
    ax.plot(xi_vals, y0, '-',  lw=2, color='C0', label='orig')
    ax.plot(xi_vals, y1, '--', lw=2, color='C1', label='fit')
    ax.set_title(rf'$-t = {tt:.1f}\ \mathrm{{GeV}}^2$')
    ax.set_xlim(0, 0.5)
    ax.set_xlabel(r'$\xi$')

axes[0].set_ylabel(r'$\mathrm{Im}\,H(\xi,\,-t)$')
axes[2].set_ylabel(r'$\mathrm{Im}\,H(\xi,\,-t)$')

# shared legend in top‐right panel
axes[1].legend(loc='upper right', fontsize=11)

# overall title
fig.suptitle(f'ImH vs ξ  (fit χ²/ndf = {chi2ndf:.2f})', y=1.02, fontsize=16)
plt.tight_layout()

# ─── Save ─────────────────────────────────────────────────────────────────────
outdir = 'output/plots'
os.makedirs(outdir, exist_ok=True)
outname = f'{outdir}/ImH_vs_xi_{timestamp}.pdf'
fig.savefig(outname, bbox_inches='tight')
print("Saved figure to", outname)