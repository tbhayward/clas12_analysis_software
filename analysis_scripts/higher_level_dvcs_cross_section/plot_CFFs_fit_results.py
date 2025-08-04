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
    # flags line: e.g., "H 1 Ht 1 E 1 Et 0"
    flag_line = next((l for l in lines if re.match(r'^H\s+\d+', l)), None)
    if flag_line is None:
        raise RuntimeError("Could not find flags line (e.g., 'H 1 Ht 1 ...') in fit file")
    toks = flag_line.split()
    flags = {}
    for i in range(0, len(toks), 2):
        key = toks[i]
        try:
            val = int(toks[i+1])
        except:
            continue
        flags[key] = val

    pnames = []
    vals = errs = None
    chi2 = ndf = chi2ndf = None
    for i, l in enumerate(lines):
        if l.startswith("# parameters"):
            parts = l.split()
            pnames = parts[2:]
        elif l.startswith("# values"):
            if i + 1 < len(lines):
                vals = np.array(list(map(float, lines[i+1].split())))
        elif l.startswith("# errors"):
            if i + 1 < len(lines):
                errs = np.array(list(map(float, lines[i+1].split())))
        elif l.startswith("# chi2"):
            if i + 1 < len(lines):
                parts = lines[i+1].split()
                chi2 = float(parts[0])
                ndf = int(float(parts[1]))
                chi2ndf = float(parts[2])
    if vals is None or errs is None or not pnames:
        raise RuntimeError("Could not parse fit-values/errors/parameter names from file")
    return flags, pnames, vals, errs, chi2, ndf, chi2ndf

flags, pnames, vals, errs, chi2, ndf, chi2ndf = parse_fit_results(fitfile)

def get_idx(name):
    try:
        return pnames.index(name)
    except ValueError:
        return None

# renormImag is fixed in your C++ code to 1.0 unless you add it to the output
renorm_imag = 1.0
if get_idx("renormImag") is not None:
    renorm_imag = vals[get_idx("renormImag")]

# ─── Defaults from the C++ physically-motivated model ───────────────────────────
defaults = {
    "H":  dict(r=0.9,   n=1.25, alpha0=0.43, alpha1=0.85, b=0.4, M2=0.64, P=1.0),
    "Ht": dict(r=7.0,   n=0.6,  alpha0=0.43, alpha1=0.85, b=2.0, M2=0.8,  P=1.0),
    "E":  dict(r=0.9,   n=1.25, alpha0=0.43, alpha1=0.85, b=0.4, M2=0.64, P=1.0),
    "Et": dict(r=1.0,   n=0.6,  alpha0=0.0,  alpha1=0.0,  b=0.0, M2=0.0,  P=0.0),
}

# ─── Extract fit parameters safely ─────────────────────────────────────────────
fit_params = {}
fit_errors = {}
for cff in ("H", "Ht", "E", "Et"):
    if flags.get(cff, 0) != 1:
        continue
    # expected parameter names with suffix
    param_keys = ["r", "n", "alpha0", "alpha1", "b", "M2", "P"]
    central = {}
    error = {}
    for k in param_keys:
        name = f"{k}_{cff}"
        idx = get_idx(name)
        if idx is not None:
            central[k] = vals[idx]
            error[k] = errs[idx]
        else:
            # fallback to default
            central[k] = defaults[cff][k]
            error[k] = 0.0  # no uncertainty if not fitted
    fit_params[cff] = central
    fit_errors[cff] = error

# ─── Build Im-CFF function matching C++ physically-motivated ansatz ───────────
def make_Im_func(cff, params, renorm):
    d = defaults[cff]
    def Im(xi, t):
        # broadcast xi and t
        xi_arr = np.array(xi, copy=False)
        t_arr = np.array(t, copy=False)
        a0 = params.get("alpha0", d["alpha0"])
        a1 = params.get("alpha1", d["alpha1"])
        nval = params.get("n", d["n"])
        rval = params.get("r", d["r"])
        bval = params.get("b", d["b"])
        M2   = params.get("M2", d["M2"])
        Pval = params.get("P", d["P"])
        alpha = a0 + a1 * t_arr
        pref = renorm * (nval * rval) / (1.0 + xi_arr)
        xfac = (2 * xi_arr / (1.0 + xi_arr)) ** (-alpha)
        yfac = ((1.0 - xi_arr) / (1.0 + xi_arr)) ** (bval)
        # protect against division by zero in the third term if M2==0
        with np.errstate(divide='ignore', invalid='ignore'):
            tfac = (1.0 - ((1.0 - xi_arr) / (1.0 + xi_arr)) * t_arr / M2) ** (-Pval) if M2 != 0 else np.ones_like(xi_arr + t_arr)
        return pref * xfac * yfac * tfac
    return Im

# ─── Replica‐band support (1σ) ────────────────────────────────────────────────
def generate_replicas(central, errors, nrep=5000):
    reps = []
    for _ in range(nrep):
        d = {}
        for k, v in central.items():
            sigma = errors.get(k, 0.0)
            if sigma > 0:
                d[k] = np.random.normal(v, sigma)
            else:
                d[k] = v
        reps.append(d)
    return reps

def compute_uncertainty_band(cff, xi_vals, t_vals, nrep=5000):
    if cff not in fit_params:
        return None, None, None
    central = fit_params[cff]
    errors_dict = fit_errors[cff]
    # build replicas of shape parameters
    param_reps = generate_replicas(central, errors_dict, nrep)
    renorm_reps = np.full(nrep, renorm_imag)  # fixed here; adjust if you fit it

    # decide output length
    if np.ndim(xi_vals) > 0 and np.ndim(t_vals) == 0:
        N = len(xi_vals)
    elif np.ndim(t_vals) > 0 and np.ndim(xi_vals) == 0:
        N = len(t_vals)
    else:
        xi_arr = np.array(xi_vals)
        t_arr = np.array(t_vals)
        broadcast = np.broadcast(xi_arr, t_arr)
        N = broadcast.shape[0]

    curves = np.empty((nrep, N))
    for i in range(nrep):
        Im_rep = make_Im_func(cff, param_reps[i], renorm_reps[i])
        val = Im_rep(xi_vals, t_vals)
        curves[i] = np.array(val).reshape(-1)[:N]  # ensure shape matches

    # handle non-finite
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
    Line2D([0],[0], color='tab:blue', linestyle='-', lw=2.5, label='Default model'),
    Line2D([0],[0], color='tab:red', linestyle='--', lw=2.5, label='Fit median'),
    Line2D([0],[0], color='tab:red', lw=6, alpha=0.2, label='1σ band'),
]

tex_map = {"H":"H", "Ht":r"\tilde H", "E":"E", "Et":r"\tilde E"}

# ─── Plot each enabled CFF ─────────────────────────────────────────────────────
for cff in ("H","Ht","E","Et"):
    if not flags.get(cff, 0):
        continue

    # default (central) Im function and fitted one
    Im_default = make_Im_func(cff, defaults[cff], renorm_imag)
    Im_fit_central = make_Im_func(cff, fit_params.get(cff, defaults[cff]), renorm_imag)
    tex = tex_map[cff]

    # — Im vs ξ at fixed t —
    fig, axes = plt.subplots(2,3, figsize=(12,8), sharex=True, sharey=True)
    axes = axes.flatten()
    fig.suptitle(rf"$\mathrm{{Im}}\,{tex}$", fontsize=16, y=0.95)

    for i,(ax,t0) in enumerate(zip(axes,t_fixed)):
        ax.plot(xi_range, Im_default(xi_range, -t0), **orig_style)
        med, lo, up = compute_uncertainty_band(cff, xi_range, -t0)
        if med is not None:
            ax.plot(xi_range, med, **fit_style)
            ax.fill_between(xi_range, lo, up, **band_style)
        ax.axhline(0, **zero_line)

        ax.set_xlim(0,0.5)
        ax.set_ylim(-2,10)
        ax.set_xticks([0,0.1,0.2,0.3,0.4,0.5])
        ax.set_yticks([-2,0,2,4,6,8,10])

        if i%3==0:
            ax.set_ylabel(r"$\mathrm{Im}\,"+tex+r"(\xi,\,-t)$")
        else:
            ax.tick_params(labelleft=False)
        ax.set_xlabel(r"$\xi$")

        # hide the "-2" tick label on top-left
        if i == 0:
            for lbl in ax.get_yticklabels():
                if lbl.get_text() in ('-2','-2.0'):
                    lbl.set_visible(False)
        # avoid clutter from "0.0" on bottom center and right
        if i in (4,5):
            for lbl in ax.get_xticklabels():
                if lbl.get_text() in ('0','0.0'):
                    lbl.set_visible(False)

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
        ax.plot(t_range, Im_default(x0, -t_range), **orig_style)
        med, lo, up = compute_uncertainty_band(cff, x0, -t_range)
        if med is not None:
            ax.plot(t_range, med, **fit_style)
            ax.fill_between(t_range, lo, up, **band_style)
        ax.axhline(0, **zero_line)

        ax.set_xlim(0,1.0)
        ax.set_ylim(-2,10)
        ax.set_xticks([0,0.2,0.4,0.6,0.8,1.0])
        ax.set_yticks([-2,0,2,4,6,8,10])

        if i%3==0:
            ax.set_ylabel(r"$\mathrm{Im}\,"+tex+r"(\xi,\,-t)$")
        else:
            ax.tick_params(labelleft=False)
        ax.set_xlabel(r"$-t\;(\mathrm{GeV^2})$")

        if i == 0:
            for lbl in ax.get_yticklabels():
                if lbl.get_text() in ('-2','-2.0'):
                    lbl.set_visible(False)
        if i in (4,5):
            for lbl in ax.get_xticklabels():
                if lbl.get_text() in ('0','0.0'):
                    lbl.set_visible(False)

        ax.text(0.60,0.65, rf"$\xi={x0:.2f}$",
                transform=ax.transAxes, fontsize=12)

    fig.subplots_adjust(left=0.08,right=0.98,bottom=0.08,top=0.92,
                        wspace=0.0,hspace=0.0)
    axes[2].legend(handles=legend_elems, loc='upper right', fontsize=10)
    fig.savefig(f"{outdir}/Im{cff}_vs_t_{timestamp}.pdf", bbox_inches='tight')
    plt.close(fig)