#!/usr/bin/env python3
"""
plot_ImH_fit_results.py

Usage:
    python plot_ImH_fit_results.py output/fit_results/fit_results_<TIMESTAMP>.txt

Reads the fitted parameters from your text file, then draws Im H(ξ,t)
for a 2×2 grid of ξ values, plotting Im H vs. −t for three fixed −t values.
Saves to output/plots/ImH_dependence_<TIMESTAMP>.pdf
"""

import os
import sys
import re
import numpy as np
import matplotlib.pyplot as plt

# ─── Parse command‐line ────────────────────────────────────────────────────────
if len(sys.argv) != 2:
    print("Usage: python plot_ImH_fit_results.py output/fit_results/fit_results_<TIMESTAMP>.txt")
    sys.exit(1)

fitfile = sys.argv[1]
m = re.search(r'fit_results_(\d{8}_\d{6})\.txt$', fitfile)
if not m:
    print("Couldn't extract timestamp from filename:", fitfile)
    sys.exit(1)
timestamp = m.group(1)

# ─── Load fit results ─────────────────────────────────────────────────────────
def load_fit(fname):
    with open(fname) as f:
        lines = [l.strip() for l in f if l.strip()]
    vals = errs = None
    chi2 = ndf = chi2ndf = None
    for i,l in enumerate(lines):
        if l.startswith("# values"):
            vals = list(map(float, lines[i+1].split()))
        if l.startswith("# errors"):
            errs = list(map(float, lines[i+1].split()))
        if l.startswith("# chi2"):
            parts = lines[i+1].split()
            chi2, ndf, chi2ndf = float(parts[0]), int(parts[1]), float(parts[2])
    if vals is None:
        raise RuntimeError("Could not find '# values' in fit file")
    return np.array(vals), np.array(errs), chi2, ndf, chi2ndf

vals, errs, chi2, ndf, chi2ndf = load_fit(fitfile)
renorm_fit, alpha0_fit, alpha1_fit, n_fit, b_fit, Mm2_fit, P_fit, renormReal_fit = vals

# ─── Original VGG defaults ───────────────────────────────────────────────────
renorm0, alpha0_0, alpha1_0 = 1.0, 0.43, 0.85
n0, b0, Mm2_0, P0           = 1.35, 0.4, 0.64, 1.0

# ─── ImH functions ────────────────────────────────────────────────────────────
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
    return ImH(xi, t, renorm_fit,
               alpha0_fit, alpha1_fit,
               n_fit, b_fit, Mm2_fit, P_fit)

# ─── Plot setup ───────────────────────────────────────────────────────────────
# try ggplot, fallback to default
try:
    plt.style.use('ggplot')
except OSError:
    pass

plt.rcParams.update({'font.size': 12, 'font.family': 'serif'})

# suggested ξ values in DVCS (ξ ~ xB/(2−xB) at moderate Q2) 
xi_vals = [0.05, 0.10, 0.20, 0.30]
t_vals  = [0.1, 0.4, 0.7]   # we'll plot at these −t

fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharey=True)
axes = axes.flatten()

colors = ['tab:blue', 'tab:orange', 'tab:green']

for ax, xi in zip(axes, xi_vals):
    for t, color in zip(t_vals, colors):
        # we only have three discrete t points, so plot them directly
        ts = np.array(t_vals)
        y0 = [ImH_orig(xi, -tt) for tt in ts]
        y1 = [ImH_fit( xi, -tt) for tt in ts]
        ax.plot(ts, y0, '-',  lw=2, color=color,
                label=f'orig, -t={t:.1f}')
        ax.plot(ts, y1, '--', lw=2, color=color,
                label=f'fit,  -t={t:.1f}')
    ax.set_title(f'ξ = {xi:.2f}')
    ax.set_xlabel('−t [GeV²]')
    ax.legend(fontsize=9)

axes[0].set_ylabel('ImH(ξ, t)')
axes[2].set_ylabel('ImH(ξ, t)')

plt.tight_layout()

# ─── Save ─────────────────────────────────────────────────────────────────────
outdir = 'output/plots'
os.makedirs(outdir, exist_ok=True)
outname = f'{outdir}/ImH_dependence_{timestamp}.pdf'
fig.savefig(outname)
print("Saved figure to", outname)