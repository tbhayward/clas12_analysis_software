#!/usr/bin/env python3
"""
plot_ImH_vs_xi.py

Reads a FitH output file (fit_results_<TIMESTAMP>.txt), extracts
the original VGG and fitted ImH parameters, then plots Im H(ξ,−t)
vs ξ for four fixed −t values in one panel.

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
    chi2 = ndf = chi2ndf = None
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if line.startswith("# values"):
                vals = next(f).split()
            if line.startswith("# chi2"):
                parts = next(f).split()
                chi2, ndf, chi2ndf = float(parts[0]), int(parts[1]), float(parts[2])
    if vals is None:
        raise RuntimeError("Could not find '# values' in fit file")
    vals = list(map(float, vals))
    return np.array(vals), chi2ndf

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
    return ImH(xi, t,
               renorm0, alpha0_0, alpha1_0,
               n0,      b0,      Mm2_0, P0)

def ImH_fit(xi, t):
    return ImH(xi, t,
               renorm_fit, alpha0_fit, alpha1_fit,
               n_fit,      b_fit,      Mm2_fit, P_fit)

# ─── Plot style ───────────────────────────────────────────────────────────────
plt.style.use('classic')
plt.rcParams.update({
    'font.size':       14,
    'font.family':     'serif',
    'axes.facecolor':  'white',
    'axes.edgecolor':  'black',
    'axes.grid':       True,
    'grid.linestyle':  '--',
    'grid.color':      '0.8',
    'legend.frameon':  True,
    'legend.framealpha': 1.0,
})

# ─── Prepare data ─────────────────────────────────────────────────────────────
# four −t values (GeV²) to compare:
t_abs = [0.1, 0.4, 0.7, 1.0]
colors = ['C0','C1','C2','C3']

# xi range:
xi_vals = np.linspace(0.02, 0.5, 200)

# ─── Single‐panel plot ────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8,6))

for tt, col in zip(t_abs, colors):
    t = -tt  # pass negative t into ImH
    y0 = ImH_orig( xi_vals, t)
    y1 = ImH_fit(  xi_vals, t)
    ax.plot( xi_vals, y0, '-',  lw=2, color=col,
             label=rf'orig, $-t={tt:.1f}$')
    ax.plot( xi_vals, y1, '--', lw=2, color=col,
             label=rf'fit,  $-t={tt:.1f}$')

ax.set_xlim(0, 0.6)
ax.set_ylim(0, None)        # auto‐scale top
ax.set_xlabel(r'$\xi$')
ax.set_ylabel(r'$\mathrm{Im}\,H(\xi,\,-t)$')
ax.set_title(f'ImH vs ξ (Fit χ²/ndf={chi2ndf:.2f})')
ax.legend(loc='upper right', fontsize=10)
plt.tight_layout()

# ─── Save ─────────────────────────────────────────────────────────────────────
outdir = 'output/plots'
os.makedirs(outdir, exist_ok=True)
outname = f'{outdir}/ImH_dependence_{timestamp}.pdf'
fig.savefig(outname)
print("Saved figure to", outname)