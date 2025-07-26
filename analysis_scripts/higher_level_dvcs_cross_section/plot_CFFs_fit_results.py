#!/usr/bin/env python3
"""
plot_ImCFFs_fit_results.py

Usage:
    python plot_ImCFFs_fit_results.py output/fit_results/fit_results_<TIMESTAMP>.txt

Reads which CFFs were fit from the header of results file, then for each
enabled Im CFF makes two figures:
  1) Im CFF vs. ξ for six fixed −t between 0.1 and 1.0 (GeV²) (2×3 grid)
  2) Im CFF vs. −t for six fixed ξ between 0.05 and 0.50 (2×3 grid)

Includes uncertainty bands for fitted results using replica method (1σ).

Saves to:
  output/plots/Im{CFF}_vs_xi_<TIMESTAMP>.pdf  
  output/plots/Im{CFF}_vs_t_<TIMESTAMP>.pdf
"""
import os
import sys
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ─── Parse command-line ─────────────────────────────────────────────────────────
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

# ─── Load fit results & flags ───────────────────────────────────────────────────
def parse_fit_results(fname):
    with open(fname) as f:
        lines = [l.strip() for l in f if l.strip()]
    flag_line = next(l for l in lines if l.startswith("H "))
    toks = flag_line.split()
    flags = { toks[i]: int(toks[i+1]) for i in range(0, len(toks), 2) }
    for i,l in enumerate(lines):
        if l.startswith("# parameters"):
            pnames = l.split()[2:]
            break
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
    if vals is None or errs is None:
        raise RuntimeError("Could not parse fit-values/errors from file")
    return flags, pnames, np.array(vals), np.array(errs), chi2, ndf, chi2ndf

flags, pnames, vals, errs, chi2, ndf, chi2ndf = parse_fit_results(fitfile)
def get_idx(name):
    return pnames.index(name) if name in pnames else None

renorm_fit = vals[get_idx("renormImag")]
renorm_err = errs[get_idx("renormImag")]

# ─── Extract shape‐parameters only (no “r”) ─────────────────────────────────────
fit_params = {}
fit_errors = {}
for cff in ("H","Ht","E","Et"):
    if flags[cff]:
        keys = ["alpha0","alpha1","n","b","M2","P"]
        base = { k: get_idx(f"{k}_{cff}") for k in keys }
        fit_params[cff] = { k: vals[idx]  for k,idx in base.items() }
        fit_errors[cff] = { k: errs[idx]  for k,idx in base.items() }

# ─── Defaults including the original r and VGG correction factor ────────────────
defaults = {
    "H":  dict(r=0.9, alpha0=0.43, alpha1=0.85, n=1.35, b=0.4,  M2=0.64, P=1.0, corr=2.0),
    "Ht": dict(r=7.0, alpha0=0.43, alpha1=0.85, n=0.6,  b=2.0,  M2=0.8,  P=1.0, corr=0.4),
    "E":  dict(r=0.9, alpha0=0.43, alpha1=0.85, n=1.35, b=0.4,  M2=0.64, P=1.0, corr=1.0),
    "Et": dict(r=1.0, alpha0=0.0,  alpha1=0.0,  n=0.0,  b=0.0,  M2=0.0,  P=0.0, corr=1.0),
}

# ─── Build Im-CFF function ─────────────────────────────────────────────────────
def make_Im_func(cff, params, renorm):
    d = defaults[cff]
    def Im(xi, t):
        a0   = params.get("alpha0", d["alpha0"])
        a1   = params.get("alpha1", d["alpha1"])
        nval = params.get("n",      d["n"])
        bval = params.get("b",      d["b"])
        M2   = params.get("M2",     d["M2"])
        Pval = params.get("P",      d["P"])
        alpha = a0 + a1 * t
        pref = renorm * np.pi * 5.0/9.0 * nval * d["r"] / (1.0 + xi)
        xfac = (2*xi/(1.0+xi))**(-alpha)
        yfac = ((1.0 - xi)/(1.0+xi))**(bval)
        tfac = (1.0 - ((1.0 - xi)/(1.0+xi))*t/M2)**(-Pval)
        return pref * xfac * yfac * tfac * d["corr"]
    return Im

# ─── Replica‐band support (1σ = 16th–84th percentile) ─────────────────────────
def generate_replicas(central, errors, nrep=10000):
    reps = []
    for _ in range(nrep):
        d = {}
        for k,v in central.items():
            sigma = errors[k] / 1.96
            d[k] = np.random.normal(v, sigma)
        reps.append(d)
    return reps

def compute_uncertainty_band(cff, xi_vals, t_vals, nrep=10000):
    if cff not in fit_params:
        return None, None, None

    # draw replicas
    param_reps  = generate_replicas(fit_params[cff], fit_errors[cff], nrep)
    renorm_reps = np.random.normal(renorm_fit, renorm_err/1.96, nrep)

    # build all curves
    curves = np.empty((nrep, len(xi_vals) if np.ndim(xi_vals)>0 else len(t_vals)))
    for i in range(nrep):
        Im_rep = make_Im_func(cff, param_reps[i], renorm_reps[i])
        curves[i] = Im_rep(xi_vals, t_vals)

    # replace non‐finite → nan, then do nan‐percentiles
    curves = np.where(np.isfinite(curves), curves, np.nan)
    med = np.nanmedian(curves, axis=0)
    lo  = np.nanpercentile(curves, 16, axis=0)
    up  = np.nanpercentile(curves, 84, axis=0)
    return med, lo, up

# ─── Plot setup ────────────────────────────────────────────────────────────────
plt.style.use('classic')
plt.rcParams.update({'font.size':14,'font.family':'serif'})
outdir = 'output/plots'
os.makedirs(outdir, exist_ok=True)

xi_range  = np.linspace(0,0.5,200)
t_range   = np.linspace(0,1.0,200)
t_fixed   = np.linspace(0.1,1.0,6)
xi_fixed  = np.linspace(0.05,0.50,6)

orig_style = {'color':'tab:blue','linestyle':'-','linewidth':2.5}
fit_style  = {'color':'tab:red','linestyle':'--','linewidth':2.5}
band_style = {'color':'tab:red','alpha':0.2}
zero_line  = {'color':'gray','linestyle':'--','linewidth':1}

legend_elems = [
    Line2D([0],[0], color='tab:blue', linestyle='-', lw=2.5, label='Original model'),
    Line2D([0],[0], color='tab:red', linestyle='--', lw=2.5, label='Fit median'),
    Line2D([0],[0], color='tab:red', lw=6, alpha=0.2, label='1σ band'),
]

tex_map = {"H":"H", "Ht":r"\tilde H", "E":"E", "Et":r"\tilde E"}

# ─── Plot each enabled CFF ─────────────────────────────────────────────────────
for cff in ("H","Ht","E","Et"):
    if not flags[cff]:
        continue

    Im_orig = make_Im_func(cff, {}, 1.0)
    tex     = tex_map[cff]

    # — Im vs ξ at fixed t —
    fig, axes = plt.subplots(2,3, figsize=(12,8), sharex=True, sharey=True)
    axes = axes.flatten()
    fig.suptitle(rf"$\mathrm{{Im}}\,{tex}$", fontsize=16, y=0.95)

    for i,(ax,t0) in enumerate(zip(axes,t_fixed)):
        ax.plot(xi_range, Im_orig(xi_range,-t0), **orig_style)
        med,lo,up = compute_uncertainty_band(cff, xi_range, -t0)
        ax.plot(xi_range, med, **fit_style)
        ax.fill_between(xi_range, lo, up, **band_style)
        ax.axhline(0, **zero_line)

        ax.set_xlim(0,0.5)
        ax.set_ylim(-2,12)
        ax.set_xticks([0,0.1,0.2,0.3,0.4,0.5])
        ax.set_yticks([-2,0,3,6,9,12])

        if i%3==0:
            ax.set_ylabel(r"$\mathrm{Im}\,"+tex+r"(\xi,\,-t)$")
        else:
            ax.tick_params(labelleft=False)

        ax.set_xlabel(r"$\xi$")
        if i==0:
            yt=ax.get_yticklabels(); yt and yt[0].set_visible(False)
        if i in (4,5):
            for lbl in ax.get_xticklabels():
                if lbl.get_text() in ('0','0.0'): lbl.set_visible(False)

        ax.text(0.60,0.65, rf"$-t={t0:.2f}\,\mathrm{{GeV^2}}$",
                transform=ax.transAxes, fontsize=12)

    fig.subplots_adjust(left=0.08,right=0.98,bottom=0.08,top=0.92,
                        wspace=0.0,hspace=0.0)
    axes[2].legend(handles=legend_elems, loc='upper right', fontsize=10)
    fig.savefig(f"{outdir}/Im{cff}_vs_xi_{timestamp}.pdf", bbox_inches='tight')
    plt.close(fig)

    # — Im vs −t at fixed ξ —
    fig, axes = plt.subplots(2,3, figsize=(12,8), sharex=True, sharey=True)
    axes = axes.flatten()
    fig.suptitle(rf"$\mathrm{{Im}}\,{tex}$", fontsize=16, y=0.95)

    for i,(ax,x0) in enumerate(zip(axes,xi_fixed)):
        ax.plot(t_range, Im_orig(x0,-t_range), **orig_style)
        med,lo,up = compute_uncertainty_band(cff, x0, -t_range)
        ax.plot(t_range, med, **fit_style)
        ax.fill_between(t_range, lo, up, **band_style)
        ax.axhline(0, **zero_line)

        ax.set_xlim(0,1.0)
        ax.set_ylim(-2,12)
        ax.set_xticks([0,0.2,0.4,0.6,0.8,1.0])
        ax.set_yticks([-2,0,3,6,9,12])

        if i%3==0:
            ax.set_ylabel(r"$\mathrm{Im}\,"+tex+r"(\xi,\,-t)$")
        else:
            ax.tick_params(labelleft=False)

        ax.set_xlabel(r"$-t\;(\mathrm{GeV^2})$")
        if i==0:
            yt=ax.get_yticklabels(); yt and yt[0].set_visible(False)
        if i in (4,5):
            for lbl in ax.get_xticklabels():
                if lbl.get_text() in ('0','0.0'): lbl.set_visible(False)

        ax.text(0.60,0.65, rf"$\xi={x0:.2f}$",
                transform=ax.transAxes, fontsize=12)

    fig.subplots_adjust(left=0.08,right=0.98,bottom=0.08,top=0.92,
                        wspace=0.0,hspace=0.0)
    axes[2].legend(handles=legend_elems, loc='upper right', fontsize=10)
    fig.savefig(f"{outdir}/Im{cff}_vs_t_{timestamp}.pdf", bbox_inches='tight')
    plt.close(fig)