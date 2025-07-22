#!/usr/bin/env python3
"""
plot_ImH_dependence.py

Usage:
    python plot_ImH_dependence.py output/fit_results/fit_results_<TIMESTAMP>.txt

Reads the fitted parameters from your text file, then makes two figures:
  1) Im H vs. ξ for six fixed −t between 0.1 and 1.0 GeV² (2×3 grid)
  2) Im H vs. −t for six fixed ξ between 0.05 and 0.50 (2×3 grid)

Saves to:
  output/plots/ImH_vs_xi_<TIMESTAMP>.pdf
  output/plots/ImH_vs_t_<TIMESTAMP>.pdf
"""

import os
import sys
import re

import numpy as np
import matplotlib.pyplot as plt

# ─── Parse command‐line ────────────────────────────────────────────────────────
if len(sys.argv) != 2:
    print("Usage: python plot_ImH_dependence.py output/fit_results/fit_results_<TIMESTAMP>.txt")
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
    pref  = np.pi * 5.0/9.0 * n_val * r / (1 + xi)
    xfac  = (2*xi/(1+xi))**(-alpha)
    yfac  = ((1 - xi)/(1+xi))**(b_val)
    tfac  = (1 - ((1 - xi)/(1+xi))*t/Mm2_val)**(-P_val)
    return renorm * pref * xfac * yfac * tfac * 2.0

def ImH_orig(xi, t):
    return ImH(xi, t, renorm0,
               alpha0_0, alpha1_0,
               n0, b0, Mm2_0, P0)

def ImH_fit(xi, t):
    return ImH(xi, t, renorm_fit,
               alpha0_fit, alpha1_fit,
               n_fit, b_fit, Mm2_fit, P_fit)

# ─── Common plot styling ───────────────────────────────────────────────────────
plt.style.use('classic')
plt.rcParams.update({
    'font.size': 14,
    'font.family': 'serif',
})

outdir = 'output/plots'
os.makedirs(outdir, exist_ok=True)

# ─── 1) ImH vs ξ for fixed −t ──────────────────────────────────────────────────
t_fixed = np.linspace(0.1, 1.0, 6)   # six −t values
xi = np.linspace(0, 0.5, 200)        # ξ range

fig1, axes1 = plt.subplots(2, 3, figsize=(12, 8),
                           sharex=True, sharey=True)
axes1 = axes1.flatten()

orig_style = {'color': 'tab:blue', 'linestyle': '-',  'linewidth': 2.5}
fit_style  = {'color': 'tab:red',  'linestyle': '--', 'linewidth': 2.5}

for idx, (ax, t) in enumerate(zip(axes1, t_fixed)):
    ax.plot(xi, ImH_orig(xi, -t), **orig_style)
    ax.plot(xi, ImH_fit(xi, -t),  **fit_style)
    if idx == 2:
        ax.legend(["Original Parameters","RGA pass-1 BSA"],
                  loc='upper right', fontsize=10)
    ax.text(0.65, 0.70,
            rf"$-t={t:.2f}\,\mathrm{{GeV}}^2$",
            transform=ax.transAxes, fontsize=12)
    ax.set_xlim(0, 0.5)
    ax.set_ylim(0, 12)

# hide zero tick labels as before
for lbl in axes1[0].get_yticklabels():
    if lbl.get_text() == '0':
        lbl.set_visible(False)
for ax in (axes1[4], axes1[5]):
    for lbl in ax.get_xticklabels():
        if lbl.get_text() == '0.0':
            lbl.set_visible(False)

fig1.subplots_adjust(left=0.10, right=0.98, bottom=0.08, top=0.97,
                     wspace=0, hspace=0)
fig1.text(0.06, 0.5, r"$\mathrm{Im}\,H(\xi,\,-t)$",
          va='center', ha='center', rotation='vertical')
fig1.text(0.5, 0.02, r"$\xi$", ha='center', va='center')

outname1 = f'{outdir}/ImH_vs_xi_{timestamp}.pdf'
fig1.savefig(outname1, bbox_inches='tight')
print("Saved:", outname1)

# ─── 2) ImH vs −t for fixed ξ ──────────────────────────────────────────────────
xi_fixed = np.linspace(0.05, 0.50, 6)  # six ξ values
t = np.linspace(0, 1.0, 200)           # continuous −t range

fig2, axes2 = plt.subplots(2, 3, figsize=(12, 8),
                           sharex=True, sharey=True)
axes2 = axes2.flatten()

for idx, (ax, xi0) in enumerate(zip(axes2, xi_fixed)):
    ax.plot(t, ImH_orig(xi0, -t), **orig_style)
    ax.plot(t, ImH_fit(xi0, -t),  **fit_style)
    if idx == 2:
        ax.legend(["Original Parameters","RGA pass-1 BSA"],
                  loc='upper right', fontsize=10)
    ax.text(0.65, 0.70,
            rf"$\xi = {xi0:.2f}$",
            transform=ax.transAxes, fontsize=12)
    ax.set_xlim(0, 1.0)
    ax.set_ylim(0, 12)

# hide zero tick labels similarly
for lbl in axes2[0].get_yticklabels():
    if lbl.get_text() == '0':
        lbl.set_visible(False)
for ax in (axes2[4], axes2[5]):
    for lbl in ax.get_xticklabels():
        if lbl.get_text() == '0.0':
            lbl.set_visible(False)

fig2.subplots_adjust(left=0.10, right=0.98, bottom=0.08, top=0.97,
                     wspace=0, hspace=0)
fig2.text(0.06, 0.5, r"$\mathrm{Im}\,H(\xi,\,-t)$",
          va='center', ha='center', rotation='vertical')
fig2.text(0.5, 0.02, r"$-t\,[\mathrm{GeV}^2]$",
          ha='center', va='center')

outname2 = f'{outdir}/ImH_vs_t_{timestamp}.pdf'
fig2.savefig(outname2, bbox_inches='tight')
print("Saved:", outname2)