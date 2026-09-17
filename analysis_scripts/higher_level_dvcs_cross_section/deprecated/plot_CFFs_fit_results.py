#!/usr/bin/env python3
"""
plot_ImCFFs_fit_results.py

Usage:
    python plot_CFFs_fit_results.py output/fit_results/fit_results_<TIMESTAMP>.txt
"""
import os
import sys
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ─── Parse command-line ───────────────────────────────────────────────────────
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

# ─── Small helpers ────────────────────────────────────────────────────────────
def build_six(minv, maxv):
    """Six values: 0.5*min, 20/40/60/80% in-range, 1.5*max."""
    if not np.isfinite(minv) or not np.isfinite(maxv) or maxv <= minv:
        return None
    return [0.5*minv,
            minv + 0.2*(maxv-minv),
            minv + 0.4*(maxv-minv),
            minv + 0.6*(maxv-minv),
            minv + 0.8*(maxv-minv),
            1.5*maxv]

def rng_str(lo, hi):
    """Plain numeric range string; no units, no math dollars."""
    if np.isfinite(lo) and np.isfinite(hi):
        return f"[{lo:.3f}, {hi:.3f}]"
    return "(unknown)"

# ─── Load fit results & flags ─────────────────────────────────────────────────
def parse_fit_results(fname):
    with open(fname) as f:
        lines = [l.rstrip("\n") for l in f]

    # flags line (e.g. "H 1 Ht 0 E 0 Et 0")
    flag_line = next((l.strip() for l in lines if re.match(r'^\s*H\s+\d+', l)), None)
    if flag_line is None:
        raise RuntimeError("Could not find flags line (e.g., 'H 1 Ht 1 ...')")
    toks = flag_line.split()
    flags = {}
    for i in range(0, len(toks), 2):
        key = toks[i]
        if i+1 < len(toks):
            try: flags[key] = int(toks[i+1])
            except: pass

    # parameter names / values / errors
    pnames, vals, errs = [], None, None
    chi2 = ndf = chi2ndf = None

    for i, l in enumerate(lines):
        if l.startswith("# parameters"):
            pnames = l.split()[2:]
        elif l.startswith("# values:"):
            if i + 1 < len(lines):
                vals = np.array([float(x) for x in lines[i+1].split()])
        elif l.startswith("# errors:"):
            if i + 1 < len(lines):
                errs = np.array([float(x) for x in lines[i+1].split()])
        elif l.startswith("# chi2"):
            if i + 1 < len(lines):
                parts = lines[i+1].split()
                chi2 = float(parts[0]); ndf = int(float(parts[1])); chi2ndf = float(parts[2])

    if vals is None or errs is None or not pnames:
        raise RuntimeError("Could not parse fit-values/errors/parameter names")

    # kinematic ranges (USED BINS block)
    xi_min = xi_max = mt_min = mt_max = np.nan
    for s in (l.strip() for l in lines):
        if s.startswith("xi_min"):
            parts = s.replace("  ", " ").split()
            try:
                xi_min = float(parts[1]) if parts[1] != "NA" else np.nan
                xi_max = float(parts[3]) if parts[3] != "NA" else np.nan
            except: pass
        elif s.startswith("-t_min"):
            parts = s.replace("  ", " ").split()
            try:
                mt_min = float(parts[1]) if parts[1] != "NA" else np.nan
                mt_max = float(parts[3]) if parts[3] != "NA" else np.nan
            except: pass

    return flags, pnames, vals, errs, chi2, ndf, chi2ndf, xi_min, xi_max, mt_min, mt_max

flags, pnames, vals, errs, chi2, ndf, chi2ndf, xi_min, xi_max, mt_min, mt_max = parse_fit_results(fitfile)

def get_idx(name):
    try: return pnames.index(name)
    except ValueError: return None

# renormImag is fixed 1.0 unless present
renorm_imag = 1.0
if get_idx("renormImag") is not None:
    renorm_imag = vals[get_idx("renormImag")]

# ─── Defaults for the “Default model” curve ───────────────────────────────────
defaults = {
    "H":  dict(r=0.9,   n=1.35, alpha0=0.43, alpha1=0.85, b=0.4, M2=0.64, P=1.0),
    "Ht": dict(r=7.0,   n=0.6,  alpha0=0.43, alpha1=0.85, b=2.0, M2=0.8,  P=1.0),
    "E":  dict(r=0.9,   n=1.35, alpha0=0.43, alpha1=0.85, b=0.4, M2=0.64, P=1.0),
    "Et": dict(r=1.0,   n=0.6,  alpha0=0.0,  alpha1=0.0,  b=0.0, M2=0.0,  P=0.0),
}

# ─── Extract fit parameters safely ────────────────────────────────────────────
fit_params, fit_errors = {}, {}
for cff in ("H", "Ht", "E", "Et"):
    if flags.get(cff, 0) != 1:
        continue
    keys = ["r","n","alpha0","alpha1","b","M2","P"]
    cent, err = {}, {}
    for k in keys:
        idx = get_idx(f"{k}_{cff}")
        if idx is not None:
            cent[k] = vals[idx]; err[k] = errs[idx]
        else:
            cent[k] = defaults[cff][k]; err[k] = 0.0
    fit_params[cff] = cent
    fit_errors[cff] = err

# ─── Im-CFF function (simple ansatz) ─────────────────────────────────────────
def make_Im_func(cff, params, renorm):
    d = defaults[cff]
    def Im(xi, t):
        xi_arr = np.array(xi, copy=False)
        t_arr  = np.array(t,  copy=False)
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
        with np.errstate(divide='ignore', invalid='ignore'):
            if M2 != 0:
                tfac = (1.0 - ((1.0 - xi_arr) / (1.0 + xi_arr)) * t_arr / M2) ** (-Pval)
            else:
                tfac = np.ones_like(xi_arr + t_arr)
        return pref * xfac * yfac * tfac
    return Im

# ─── Replica bands (1σ) ───────────────────────────────────────────────────────
def generate_replicas(central, errors, nrep=2000):
    reps = []
    for _ in range(nrep):
        d = {}
        for k, v in central.items():
            s = errors.get(k, 0.0)
            d[k] = np.random.normal(v, s) if s > 0 else v
        reps.append(d)
    return reps

def compute_uncertainty_band(cff, xi_vals, t_vals, nrep=2000):
    if cff not in fit_params: return None, None, None
    params = fit_params[cff]; errs = fit_errors[cff]
    reps = generate_replicas(params, errs, nrep)
    if np.ndim(xi_vals) > 0 and np.ndim(t_vals) == 0:
        N = len(xi_vals)
    elif np.ndim(t_vals) > 0 and np.ndim(xi_vals) == 0:
        N = len(t_vals)
    else:
        N = np.broadcast(np.array(xi_vals), np.array(t_vals)).shape[0]
    curves = np.empty((nrep, N))
    for i in range(nrep):
        Im_rep = make_Im_func(cff, reps[i], renorm_imag)
        curves[i] = np.array(Im_rep(xi_vals, t_vals)).reshape(-1)[:N]
    curves = np.where(np.isfinite(curves), curves, np.nan)
    med = np.nanmedian(curves, axis=0)
    lo  = np.nanpercentile(curves, 16, axis=0)
    up  = np.nanpercentile(curves, 84, axis=0)
    return med, lo, up

# ─── Dynamic kinematics from results file ─────────────────────────────────────
xi_ok = (np.isfinite(xi_min) and np.isfinite(xi_max) and xi_max > xi_min and xi_min > 0)
mt_ok = (np.isfinite(mt_min) and np.isfinite(mt_max) and mt_max > mt_min and mt_min >= 0)

if xi_ok:
    xi_lo_draw = max(1e-6, 0.5*xi_min)
    xi_hi_draw = min(1.0 - 1e-6, 1.5*xi_max)
else:
    xi_lo_draw, xi_hi_draw = 0.00, 0.50

if mt_ok:
    mt_lo_draw = max(0.0, 0.5*mt_min)
    mt_hi_draw = 1.5*mt_max
else:
    mt_lo_draw, mt_hi_draw = 0.00, 0.60

xi_range = np.linspace(xi_lo_draw, xi_hi_draw, 400)
t_range  = np.linspace(mt_lo_draw, mt_hi_draw, 400)

# six fixed panel values
t_fixed  = build_six(mt_min, mt_max) if mt_ok else [0.1,0.2,0.3,0.4,0.5,0.6]
xi_fixed = build_six(xi_min, xi_max) if xi_ok else [0.05,0.15,0.25,0.35,0.45,0.50]

# Titles: show USED-BIN ranges; keep units as LaTeX without extra dollars
xi_title_math = rng_str(xi_min, xi_max)
t_title_math  = (rng_str(mt_min, mt_max) + r"\,\mathrm{GeV}^2") if mt_ok else "(unknown)"

# ─── Plot setup ────────────────────────────────────────────────────────────────
plt.style.use('classic')
plt.rcParams.update({'font.size':14,'font.family':'serif'})
outdir = 'output/plots'
os.makedirs(outdir, exist_ok=True)

orig_style = {'color':'tab:blue','linestyle':'-','linewidth':2.5}
fit_style  = {'color':'tab:red','linestyle':'--','linewidth':2.5}
band_style = {'color':'tab:red','alpha':0.2}
zero_line  = {'color':'gray','linestyle':'--','linewidth':1}

legend_elems = [
    Line2D([0],[0], color='tab:blue', linestyle='-', lw=2.5, label='Default model'),
    Line2D([0],[0], color='tab:red',  linestyle='--', lw=2.5, label='Fit median'),
    Line2D([0],[0], color='tab:red',  lw=6, alpha=0.2, label='1σ band'),
]

# Common y-range and ticks
Y_MIN, Y_MAX = 0.0, 3.0
YTICKS_ALL  = [0, 1, 2, 3]
YTICKS_TL   = [1, 2, 3]     # top-left: hide 0 tick

tex_map = {"H":"H", "Ht":r"\tilde H", "E":"E", "Et":r"\tilde E"}

# ─── Plot each enabled CFF ────────────────────────────────────────────────────
for cff in ("H","Ht","E","Et"):
    if not flags.get(cff, 0):
        continue

    Im_default = make_Im_func(cff, defaults[cff], renorm_imag)
    tex = tex_map[cff]

    # — Im vs ξ at fixed -t —
    fig, axes = plt.subplots(2,3, figsize=(12,8), sharex=True, sharey=False)
    axes = axes.flatten()
    fig.suptitle(rf"$\mathrm{{Im}}\,{tex}$; "
                 rf"applicability: $\xi\in{xi_title_math}$, $-t\in{t_title_math}$",
                 fontsize=14, y=0.98)

    for i,(ax,mt0) in enumerate(zip(axes, t_fixed)):
        ax.plot(xi_range, Im_default(xi_range, -mt0), **orig_style)
        med, lo, up = compute_uncertainty_band(cff, xi_range, -mt0)
        if med is not None:
            ax.plot(xi_range, med, **fit_style)
            ax.fill_between(xi_range, lo, up, **band_style)
        ax.axhline(0, **zero_line)

        # y-lims & ticks
        ax.set_xlim(xi_range[0], xi_range[-1])
        ax.set_ylim(Y_MIN, Y_MAX)
        if i == 0:
            ax.set_yticks(YTICKS_TL)          # top-left: no 0 tick
        elif i % 3 == 0:
            ax.set_yticks(YTICKS_ALL)         # left column: full ticks
        else:
            ax.set_yticks(YTICKS_ALL)
            ax.tick_params(labelleft=False)   # hide labels on middle/right

        ax.set_xlabel(r"$\xi$")
        if i % 3 == 0:
            ax.set_ylabel(r"$\mathrm{Im}\,"+tex+r"(\xi,\,-t)$")

        # Move annotation lower so it doesn't overlap the curves
        ax.text(0.58, 0.60, rf"$-t={mt0:.3f}\,\mathrm{{GeV^2}}$",
                transform=ax.transAxes, fontsize=12)

    fig.subplots_adjust(left=0.08,right=0.98,bottom=0.08,top=0.90,
                        wspace=0.0,hspace=0.0)
    axes[2].legend(handles=legend_elems, loc='upper right', fontsize=10)
    fig.savefig(f"{outdir}/Im{cff}_vs_xi_{timestamp}.pdf", bbox_inches='tight')
    plt.close(fig)

    # — Im vs -t at fixed ξ —
    fig, axes = plt.subplots(2,3, figsize=(12,8), sharex=True, sharey=False)
    axes = axes.flatten()
    fig.suptitle(rf"$\mathrm{{Im}}\,{tex}$;  "
                 rf"applicability: $\xi\in{xi_title_math}$, $-t\in{t_title_math}$",
                 fontsize=14, y=0.98)

    for i,(ax,xi0) in enumerate(zip(axes, xi_fixed)):
        ax.plot(t_range, Im_default(xi0, -t_range), **orig_style)
        med, lo, up = compute_uncertainty_band(cff, xi0, -t_range)
        if med is not None:
            ax.plot(t_range, med, **fit_style)
            ax.fill_between(t_range, lo, up, **band_style)
        ax.axhline(0, **zero_line)

        # y-lims & ticks
        ax.set_xlim(t_range[0], t_range[-1])
        ax.set_ylim(Y_MIN, Y_MAX)
        if i == 0:
            ax.set_yticks(YTICKS_TL)          # top-left: no 0 tick
        elif i % 3 == 0:
            ax.set_yticks(YTICKS_ALL)         # left column: full ticks
        else:
            ax.set_yticks(YTICKS_ALL)
            ax.tick_params(labelleft=False)   # hide labels on middle/right

        ax.set_xlabel(r"$-t\;(\mathrm{GeV^2})$")
        if i % 3 == 0:
            ax.set_ylabel(r"$\mathrm{Im}\,"+tex+r"(\xi,\,-t)$")

        # Move annotation lower so it doesn't overlap the curves
        ax.text(0.62, 0.60, rf"$\xi={xi0:.3f}$",
                transform=ax.transAxes, fontsize=12)

    fig.subplots_adjust(left=0.08,right=0.98,bottom=0.08,top=0.90,
                        wspace=0.0,hspace=0.0)
    axes[2].legend(handles=legend_elems, loc='upper right', fontsize=10)
    fig.savefig(f"{outdir}/Im{cff}_vs_t_{timestamp}.pdf", bbox_inches='tight')
    plt.close(fig)