#!/usr/bin/env python3
"""
plot_ImCFFs_fit_results.py

Usage:
    python plot_ImCFFs_fit_results.py output/fit_results/fit_results_<TIMESTAMP>.txt

Reads which CFFs were fit from the header of your results file, then for each
enabled Im CFF makes two figures:
  1) Im CFF vs. ξ for six fixed −t between 0.1 and 1.0 (GeV^{2}) (2×3 grid)
  2) Im CFF vs. −t for six fixed ξ between 0.05 and 0.50 (2×3 grid)

Saves to:
  output/plots/Im{CFF}_vs_xi_<TIMESTAMP>.pdf  
  output/plots/Im{CFF}_vs_t_<TIMESTAMP>.pdf
"""
import os
import sys
import re

import numpy as np
import matplotlib.pyplot as plt

# ─── Parse command‐line ────────────────────────────────────────────────────────
if len(sys.argv) != 2:
    print("Usage: python plot_ImCFFs_fit_results.py "
          "output/fit_results/fit_results_<TIMESTAMP>.txt")
    sys.exit(1)

fitfile = sys.argv[1]
m = re.search(r'fit_results_(\d{8}_\d{6})\.txt$', fitfile)
if not m:
    print("Couldn't extract timestamp from filename:", fitfile)
    sys.exit(1)
timestamp = m.group(1)

# ─── Load fit results & flags ─────────────────────────────────────────────────
def parse_fit_results(fname):
    params = {}
    with open(fname) as f:
        lines = [l.strip() for l in f if l.strip()]
    # find flags line: e.g. "H 1  Ht 1  E 0  Et 1"
    flag_line = next(l for l in lines if l.startswith("H "))
    toks = flag_line.split()
    flags = { toks[i]: int(toks[i+1]) for i in range(0, len(toks), 2) }

    # parameter names
    pnames = []
    for i,l in enumerate(lines):
        if l.startswith("# parameters"):
            pnames = lines[i].split()[2:]
            break
    # values
    vals = None
    chi2 = ndf = chi2ndf = None
    for i,l in enumerate(lines):
        if l.startswith("# values"):
            vals = list(map(float, lines[i+1].split()))
        if l.startswith("# chi2"):
            parts = lines[i+1].split()
            chi2, ndf, chi2ndf = float(parts[0]), int(parts[1]), float(parts[2])
    if vals is None:
        raise RuntimeError("Could not parse fit‐values from file")
    return flags, pnames, np.array(vals), chi2, ndf, chi2ndf

flags, pnames, vals, chi2, ndf, chi2ndf = parse_fit_results(fitfile)

# ─── Extract fit‐parameters for each CFF ───────────────────────────────────────
get_idx = lambda name: pnames.index(name) if name in pnames else None
renorm_fit = vals[get_idx("renormImag")]

fit_params = {}
for cff in ("H","Ht","E"):
    if flags[cff]:
        base = {"r":get_idx(f"r_{cff}"),
                "alpha0":get_idx(f"alpha0_{cff}"),
                "alpha1":get_idx(f"alpha1_{cff}"),
                "n":get_idx(f"n_{cff}"),
                "b":get_idx(f"b_{cff}"),
                "Mm2":get_idx(f"Mm2_{cff}"),
                "P":get_idx(f"P_{cff}")}
        fit_params[cff] = { key: vals[idx] for key,idx in base.items() }

if flags["Et"]:
    fit_params["Et"] = {}  # no shape params

# ─── VGG‐default parameters ────────────────────────────────────────────────────
defaults = {
    "H":  dict(r=0.9,   alpha0=0.43, alpha1=0.85, n=1.35, b=0.4,  Mm2=0.64, P=1.0, factor=2.0),
    "Ht": dict(r=7.0,   alpha0=0.43, alpha1=0.85, n=0.6,  b=2.0,  Mm2=0.8,  P=1.0, factor=0.4),
    "E":  dict(r=0.9,   alpha0=0.43, alpha1=0.85, n=1.35, b=0.4,  Mm2=0.64, P=1.0, factor=1.0),
    "Et": dict(factor=0.0)
}

# ─── Im CFF generator ──────────────────────────────────────────────────────────
def make_Im_func(cff, params, renorm):
    d = defaults[cff]
    def Im(xi, t):
        if cff=="Et":
            # always zero
            return np.zeros_like(xi if hasattr(xi,'__iter__') else t)
        r       = params.get("r",      d["r"])
        alpha0  = params.get("alpha0", d["alpha0"])
        alpha1  = params.get("alpha1", d["alpha1"])
        n_val   = params.get("n",      d["n"])
        b_val   = params.get("b",      d["b"])
        Mm2_val = params.get("Mm2",    d["Mm2"])
        P_val   = params.get("P",      d["P"])
        factor  = d["factor"]
        alpha = alpha0 + alpha1 * t
        pref  = np.pi * 5.0/9.0 * n_val * r / (1 + xi)
        xfac  = (2*xi/(1+xi))**(-alpha)
        yfac  = ((1 - xi)/(1+xi))**(b_val)
        tfac  = (1 - ((1 - xi)/(1+xi))*t/Mm2_val)**(-P_val)
        return renorm * pref * xfac * yfac * tfac * factor
    return Im

Im_funcs = {}
for cff in ("H","Ht","E","Et"):
    if flags[cff]:
        Im_funcs[cff] = make_Im_func(cff, fit_params.get(cff,{}), renorm_fit)

# ─── Common plot styling ───────────────────────────────────────────────────────
plt.style.use('classic')
plt.rcParams.update({
    'font.size': 14,
    'font.family': 'serif',
})

outdir = 'output/plots'
os.makedirs(outdir, exist_ok=True)

# ─── Fixed grids ──────────────────────────────────────────────────────────────
t_fixed  = np.linspace(0.1, 1.0, 6)
xi_range = np.linspace(0,   0.5, 200)
xi_fixed = np.linspace(0.05,0.50,6)
t_range  = np.linspace(0,   1.0, 200)

orig_style = {'color': 'tab:blue', 'linestyle': '-',  'linewidth': 2.5}
fit_style  = {'color': 'tab:red',  'linestyle': '--', 'linewidth': 2.5}
zero_line  = {'color': 'gray',     'linestyle': '--', 'linewidth': 1.0}

# LaTeX names for titles
latex_map = {
    "H":  r"H",
    "Ht": r"\tilde H",
    "E":  r"E",
    "Et": r"\tilde E"
}

# ─── Loop over each enabled CFF and make its two plots ─────────────────────────
for cff, Im_fit in Im_funcs.items():
    Im_orig = make_Im_func(cff, {}, 1.0)

    # 1) ImCFF vs ξ for fixed −t
    fig, axes = plt.subplots(2,3,figsize=(12,8),sharex=True,sharey=True)
    axes = axes.flatten()
    for i,(ax, tval) in enumerate(zip(axes, t_fixed)):
        ax.plot(xi_range, Im_orig(xi_range, -tval), **orig_style)
        ax.plot(xi_range, Im_fit(xi_range,  -tval), **fit_style)
        ax.axhline(0, **zero_line)
        ax.text(0.62, 0.72, rf"$-t={tval:.2f}\,\mathrm{{(GeV^2)}}$",
                transform=ax.transAxes, fontsize=12)
        ax.set_xlim(0,0.5)
        ax.set_ylim(-2,12)
        # remove the "-2" tick on top-left only
        if i==0:
            for lbl in ax.get_yticklabels():
                if lbl.get_text() == '-2':
                    lbl.set_visible(False)
    axes[2].legend(["Original Parameters","RGA pass-1 fit"],
                   loc='upper right', fontsize=10)
    fig.subplots_adjust(left=0.10,right=0.98,bottom=0.08,top=0.93,wspace=0,hspace=0)
    fig.suptitle(rf"$\mathrm{{Im}}\,{latex_map[cff]}(\xi,\,-t)$", y=0.97, fontsize=16)
    fig.text(0.06,0.5, fr"$\mathrm{{Im}}\,{latex_map[cff]}(\xi,\,-t)$",
             va='center',ha='center',rotation='vertical')
    fig.text(0.5,0.02, r"$\xi$", ha='center')
    out1 = f"{outdir}/Im{cff}_vs_xi_{timestamp}.pdf"
    fig.savefig(out1, bbox_inches='tight')
    print("Saved:", out1)
    plt.close(fig)

    # 2) ImCFF vs −t for fixed ξ
    fig, axes = plt.subplots(2,3,figsize=(12,8),sharex=True,sharey=True)
    axes = axes.flatten()
    for i,(ax, xi0) in enumerate(zip(axes, xi_fixed)):
        ax.plot(t_range, Im_orig(xi0, -t_range), **orig_style)
        ax.plot(t_range, Im_fit(xi0,  -t_range), **fit_style)
        ax.axhline(0, **zero_line)
        ax.text(0.62, 0.72, rf"$\xi={xi0:.2f}$",
                transform=ax.transAxes, fontsize=12)
        ax.set_xlim(0,1.0)
        ax.set_ylim(-2,12)
        # remove "0.0" ticks on bottom-middle (i==4) and bottom-right (i==5)
        if i in (4,5):
            for lbl in ax.get_xticklabels():
                if lbl.get_text() == '0.0':
                    lbl.set_visible(False)
    axes[2].legend(["Original Parameters","RGA pass-1 fit"],
                   loc='upper right', fontsize=10)
    fig.subplots_adjust(left=0.10,right=0.98,bottom=0.08,top=0.93,wspace=0,hspace=0)
    fig.suptitle(rf"$\mathrm{{Im}}\,{latex_map[cff]}(\xi,\,-t)$", y=0.97, fontsize=16)
    fig.text(0.06,0.5, fr"$\mathrm{{Im}}\,{latex_map[cff]}(\xi,\,-t)$",
             va='center',ha='center',rotation='vertical')
    fig.text(0.5,0.02, r"$-t\ (\mathrm{GeV^2})$", ha='center')
    out2 = f"{outdir}/Im{cff}_vs_t_{timestamp}.pdf"
    fig.savefig(out2, bbox_inches='tight')
    print("Saved:", out2)
    plt.close(fig)