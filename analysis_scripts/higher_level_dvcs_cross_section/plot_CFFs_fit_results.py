#!/usr/bin/env python3
"""
plot_ImCFFs_fit_results.py

Usage:
    python plot_ImCFFs_fit_results.py output/fit_results/fit_results_<TIMESTAMP>.txt

Reads which CFFs were fit from the header of results file, then for each
enabled Im CFF makes two figures:
  1) Im CFF vs. ξ for six fixed −t between 0.1 and 1.0 (GeV²) (2×3 grid)
  2) Im CFF vs. −t for six fixed ξ between 0.05 and 0.50 (2×3 grid)

Includes uncertainty bands for fitted results using replica method (95% CI).

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

# ─── Extract fit-parameters ─────────────────────────────────────────────────────
fit_params = {}
fit_errors = {}
for cff in ("H","Ht","E","Et"):
    if flags[cff]:
        keys = ["r","alpha0","alpha1","n","b","Mm2","P"]
        base = { k: get_idx(f"{k}_{cff}") for k in keys }
        fit_params[cff] = { k: vals[idx] for k,idx in base.items() }
        fit_errors[cff] = { k: errs[idx] for k,idx in base.items() }

# ─── VGG-default params + factors ───────────────────────────────────────────────
defaults = {
    "H":  dict(r=0.9, alpha0=0.43, alpha1=0.85, n=1.35, b=0.4, Mm2=0.64, P=1.0, corr=2.0),
    "Ht": dict(r=7.0, alpha0=0.43, alpha1=0.85, n=0.6,  b=2.0, Mm2=0.8,  P=1.0, corr=0.4),
    "E":  dict(r=0.9, alpha0=0.43, alpha1=0.85, n=1.35, b=0.4, Mm2=0.64, P=1.0, corr=1.0),
    "Et": dict(r=1.0, alpha0=0.0,  alpha1=0.0,  n=0.0,  b=0.0, Mm2=0.0,  P=0.0, corr=1.0),
}

# ─── Build the Im‐CFF function exactly matching your C++ code ────────────────────
def make_Im_func(cff, params, renorm):
    d = defaults[cff]
    def Im(xi, t):
        a0   = params.get("alpha0", d["alpha0"])
        a1   = params.get("alpha1", d["alpha1"])
        nval = params.get("n",      d["n"])
        bval = params.get("b",      d["b"])
        Mm2  = params.get("Mm2",    d["Mm2"])
        Pval = params.get("P",      d["P"])
        r    = params.get("r",      d["r"])
        fac  = d["corr"]
        alpha = a0 + a1 * t
        pref = renorm * np.pi * 5.0/9.0 * nval * r / (1.0 + xi)
        xfac = (2*xi/(1.0+xi))**(-alpha)
        yfac = ((1.0 - xi)/(1.0+xi))**(bval)
        tfac = (1.0 - ((1.0 - xi)/(1.0+xi))*t/Mm2)**(-Pval)
        return pref * xfac * yfac * tfac * fac
    return Im

# ─── Replica‐band support ───────────────────────────────────────────────────────
def generate_replicas(central, errors, nrep=200):
    reps = []
    for _ in range(nrep):
        rdict = {}
        for k,v in central.items():
            sigma = errors[k] / 1.96
            rdict[k] = np.random.normal(v, sigma)
        reps.append(rdict)
    return reps

def compute_uncertainty_band(cff, xi_vals, t_vals, nrep=200):
    if cff not in fit_params:
        return None, None, None
    param_reps  = generate_replicas(fit_params[cff], fit_errors[cff], nrep)
    renorm_reps = np.random.normal(renorm_fit, renorm_err/1.96, nrep)
    curves = np.empty((nrep, len(xi_vals) if np.ndim(xi_vals)>0 else len(t_vals)))
    for i in range(nrep):
        Im_rep = make_Im_func(cff, param_reps[i], renorm_reps[i])
        # xi_vals vs t_vals must match DVCS sign convention: pass -t
        if np.ndim(t_vals) > 0:
            curves[i] = Im_rep(xi_vals, t_vals)
        else:
            curves[i] = Im_rep(xi_vals, t_vals)
    return (
        np.median(curves, axis=0),
        np.percentile(curves, 2.5, axis=0),
        np.percentile(curves, 97.5, axis=0),
    )

# ─── Set up plotting ──────────────────────────────────────────────────────────
plt.style.use('classic')
plt.rcParams.update({'font.size':14,'font.family':'serif'})

outdir = 'output/plots'
os.makedirs(outdir, exist_ok=True)

xi_range  = np.linspace(0,0.5,200)
t_range   = np.linspace(0,1.0,200)
t_fixed   = np.linspace(0.1,1.0,6)
xi_fixed  = np.linspace(0.05,0.50,6)

orig_style = {'color':'tab:blue','linestyle':'-','linewidth':2.5}
fit_style  = {'color':'tab:red', 'linestyle':'--','linewidth':2.5}
band_style = {'color':'tab:red', 'alpha':0.2}
zero_line  = {'color':'gray','linestyle':'--','linewidth':1}

legend_elems = [
    Line2D([0],[0], color='tab:blue', linestyle='-', lw=2.5, label='Original model'),
    Line2D([0],[0], color='tab:red',  linestyle='--', lw=2.5, label='Fit median'),
    plt.Rectangle((0,0),1,1,fc='tab:red',alpha=0.2,ec=None,label='95% CI'),
]

tex_map = {"H":"H", "Ht":r"\tilde H", "E":"E", "Et":r"\tilde E"}

# ─── Plot loop ────────────────────────────────────────────────────────────────
for cff in ("H","Ht","E","Et"):
    if not flags[cff]:
        continue

    Im_fit  = make_Im_func(cff, fit_params.get(cff, {}), renorm_fit)
    Im_orig = make_Im_func(cff, {}, 1.0)
    tex     = tex_map[cff]

    # — Im vs ξ @ fixed t —
    fig, axes = plt.subplots(2,3,figsize=(12,8), sharex=True, sharey=True)
    axes = axes.flatten()
    fig.suptitle(rf"$\mathrm{{Im}}\,{tex}$", fontsize=16, y=0.95)

    for ax, t0 in zip(axes, t_fixed):
        # pass -t0 to match C++ convention
        ax.plot(xi_range, Im_orig(xi_range, -t0), **orig_style)
        med, lo, up = compute_uncertainty_band(cff, xi_range, -t0)
        ax.plot(xi_range, med, **fit_style)
        ax.fill_between(xi_range, lo, up, **band_style)
        ax.axhline(0, **zero_line)
        ax.text(0.6, 0.75, rf"$-t={t0:.2f}$", transform=ax.transAxes)
        ax.set_xlim(0,0.5)
        ax.set_ylim(0,None)

    axes[2].legend(handles=legend_elems, loc='upper right', fontsize=10)
    fig.text(0.06,0.5, rf"$\mathrm{{Im}}\,{tex}(\xi,\,-t)$",
             va='center', ha='center', rotation='vertical')
    fig.text(0.5,0.02, r"$\xi$", ha='center')

    out1 = f"{outdir}/Im{cff}_vs_xi_{timestamp}.pdf"
    fig.savefig(out1, bbox_inches='tight')
    print("Saved:", out1)
    plt.close(fig)

    # — Im vs −t @ fixed ξ —
    fig, axes = plt.subplots(2,3,figsize=(12,8), sharex=True, sharey=True)
    axes = axes.flatten()
    fig.suptitle(rf"$\mathrm{{Im}}\,{tex}$", fontsize=16, y=0.95)

    for ax, x0 in zip(axes, xi_fixed):
        # again pass -t-range
        ax.plot(t_range, Im_orig(x0, -t_range), **orig_style)
        med, lo, up = compute_uncertainty_band(cff, x0, -t_range)
        ax.plot(t_range, med, **fit_style)
        ax.fill_between(t_range, lo, up, **band_style)
        ax.axhline(0, **zero_line)
        ax.text(0.6, 0.75, rf"$\xi={x0:.2f}$", transform=ax.transAxes)
        ax.set_xlim(0,1.0)
        ax.set_ylim(0,None)

    axes[2].legend(handles=legend_elems, loc='upper right', fontsize=10)
    fig.text(0.06,0.5, rf"$\mathrm{{Im}}\,{tex}(\xi,\,-t)$",
             va='center', ha='center', rotation='vertical')
    fig.text(0.5,0.02, r"$-t\,(\mathrm{GeV}^2)$", ha='center')

    out2 = f"{outdir}/Im{cff}_vs_t_{timestamp}.pdf"
    fig.savefig(out2, bbox_inches='tight')
    print("Saved:", out2)
    plt.close(fig)