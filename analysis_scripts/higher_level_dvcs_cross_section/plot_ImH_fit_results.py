#!/usr/bin/env python3
"""
plot_ImH_vs_xi.py

Usage:
    python plot_ImH_vs_xi.py output/fit_results/fit_results_<TIMESTAMP>.txt

Reads the fitted parameters from your text file, then draws Im H(ξ,−t)
for a 2×3 grid of −t values, plotting Im H vs. ξ for six equally spaced
−t between 0.1 and 1.0 GeV². Saves to
output/plots/ImH_vs_xi_<TIMESTAMP>.pdf
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
    print("Couldn't extract timestamp from filename:", fitfile)
    sys.exit(1)
timestamp = m.group(1)

# ─── Load fit results ─────────────────────────────────────────────────────────
def parse_fit_results(fname):
    vals = None
    chi2 = ndf = chi2ndf = None
    with open(fname) as f:
        lines = [l.strip() for l in f if l.strip()]
    for i, l in enumerate(lines):
        if l.startswith("# values"):
            vals = list(map(float, lines[i+1].split()))
        if l.startswith("# chi2"):
            parts = lines[i+1].split()
            chi2, ndf, chi2ndf = float(parts[0]), int(parts[1]), float(parts[2])
    if vals is None:
        raise RuntimeError("Could not parse fit-values from file")
    return np.array(vals), chi2, ndf, chi2ndf

vals, chi2, ndf, chi2ndf = parse_fit_results(fitfile)
(renorm_fit, alpha0_fit, alpha1_fit,
 n_fit, b_fit, Mm2_fit, P_fit,
 renormR_fit) = vals

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
    return ImH(xi, t, renorm0,
               alpha0_0, alpha1_0,
               n0, b0, Mm2_0, P0)

def ImH_fit(xi, t):
    return ImH(xi, t, renorm_fit,
               alpha0_fit, alpha1_fit,
               n_fit, b_fit, Mm2_fit, P_fit)

# ─── Plot setup ───────────────────────────────────────────────────────────────
plt.style.use('classic')
plt.rcParams.update({
    'font.size': 14,
    'font.family': 'serif',
})

# six −t values equally spaced [0.1, 1.0]
t_vals = np.linspace(0.1, 1.0, 6)
# ξ range [0, 0.5]
xi = np.linspace(0, 0.5, 200)

fig, axes = plt.subplots(2, 3, figsize=(12, 8),
                         sharex=True, sharey=True)
axes = axes.flatten()

# styles
orig_style = {'color': 'tab:blue', 'linestyle': '-', 'linewidth': 2.5}
fit_style  = {'color': 'tab:red',  'linestyle': '--','linewidth': 2.5}

for idx, (ax, t) in enumerate(zip(axes, t_vals)):
    ax.plot(xi, ImH_orig(xi, -t), **orig_style)
    ax.plot(xi, ImH_fit(xi, -t),  **fit_style)

    # legend only in top-right subplot (idx==2)
    if idx == 2:
        ax.legend(["Original Parameters","RGA pass-1 BSA"],
                  loc='upper right', fontsize=10)

    # annotation of t at (0.65,0.70)
    ax.text(0.65, 0.70,
            rf"$-t = {t:.2f}\,\mathrm{{GeV}}^2$",
            transform=ax.transAxes,
            fontsize=12)

    # axes limits
    ax.set_xlim(0, 0.5)
    ax.set_ylim(0, 12)

# ─── Tweak tick‐labels ────────────────────────────────────────────────────────
# remove '0' on top-left y-axis
for lbl in axes[0].get_yticklabels():
    if lbl.get_text() == '0':
        lbl.set_visible(False)
# remove '0.0' on bottom-center & bottom-right x-axes
for ax in (axes[4], axes[5]):
    for lbl in ax.get_xticklabels():
        if lbl.get_text() == '0.0':
            lbl.set_visible(False)

# ─── Global axis labels & layout ───────────────────────────────────────────────
# shrink left margin, move y-label closer (x=0.06)
fig.subplots_adjust(left=0.10, right=0.98, bottom=0.08, top=0.97,
                    wspace=0, hspace=0)
fig.text(0.06, 0.5,
         r"$\mathrm{Im}\,H(\xi,\,-t)$",
         va='center', ha='center', rotation='vertical')
fig.text(0.5, 0.02,
         r"$\xi$",
         ha='center', va='center')

# ─── Save ─────────────────────────────────────────────────────────────────────
outdir = 'output/plots'
os.makedirs(outdir, exist_ok=True)
outname = f'{outdir}/ImH_vs_xi_{timestamp}.pdf'
fig.savefig(outname, bbox_inches='tight')
print("Saved figure to", outname)