#!/usr/bin/env python3
"""
plot_CFFs_fit_RGA_RGK.py

Usage:
    python plot_CFFs_fit_RGA_RGK.py <fitfile_RGA.txt> <fitfile_RGK.txt>

Reads two CFF‐fit result files (RGA and RGK), makes for each:
  1) Im CFF vs. ξ for fixed −t (2×3 grid),
  2) Im CFF vs. −t for fixed ξ (2×3 grid).

First it will produce the "single‐file" plots (only RGK, labeled “RGK preliminary”),
then "comparison" plots overlaying RGA and RGK results (_2 suffix).

All plots go into output/plots/.
"""
import os
import sys
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from matplotlib.lines import Line2D

# ─── Parse command-line ─────────────────────────────────────────────────────────
if len(sys.argv) != 3:
    print("Usage: python plot_CFFs_fit_RGA_RGK.py "
          "<fitfile_RGA.txt> <fitfile_RGK.txt>")
    sys.exit(1)

fitfile_rga = sys.argv[1]
fitfile_rgk = sys.argv[2]

def extract_timestamp(fn):
    m = re.search(r'fit_results_(\d{8}_\d{6})\.txt$', fn)
    return m.group(1) if m else "unknown"

ts_rga = extract_timestamp(fitfile_rga)
ts_rgk = extract_timestamp(fitfile_rgk)

# ─── Common I/O & parsing ───────────────────────────────────────────────────────
def parse_fit_results(fname):
    with open(fname) as f:
        lines = [l.strip() for l in f if l.strip()]
    # flags
    flag_line = next(l for l in lines if l.startswith("H "))
    toks = flag_line.split()
    flags = { toks[i]: int(toks[i+1]) for i in range(0, len(toks), 2) }
    # parameter names
    pnames = []
    for i,l in enumerate(lines):
        if l.startswith("# parameters"):
            pnames = lines[i].split()[2:]
            break
    # values & errors
    vals = errs = None
    chi2 = ndf = chi2ndf = None
    for i,l in enumerate(lines):
        if l.startswith("# values"):
            vals = list(map(float, lines[i+1].split()))
        if l.startswith("# errors"):
            errs = list(map(float, lines[i+1].split()))
        if l.startswith("# chi2"):
            c, n, cndf = lines[i+1].split()
            chi2, ndf, chi2ndf = float(c), int(n), float(cndf)
    return flags, pnames, np.array(vals), np.array(errs), chi2, ndf, chi2ndf

flags_rga, pnames_rga, vals_rga, errs_rga, *_ = parse_fit_results(fitfile_rga)
flags_rgk, pnames_rgk, vals_rgk, errs_rgk, *_ = parse_fit_results(fitfile_rgk)

def get_idx(pnames, name):
    return pnames.index(name) if name in pnames else None

# ─── Build dicts for RGA & RGK ──────────────────────────────────────────────────
def build_param_dicts(pnames, vals, errs):
    ren_idx = get_idx(pnames, "renormImag")
    renorm = vals[ren_idx]
    renorm_err = errs[ren_idx]
    fps = {}
    fes = {}
    for cff in ("H","Ht","E","Et"):
        if cff_map[cff]["flag"]:
            keys = ["r","alpha0","alpha1","n","b","M2","P"]
            idxs = { k: get_idx(pnames, f"{k}_{cff}") for k in keys }
            fps[cff] = { k: vals[i] for k,i in idxs.items() }
            fes[cff] = { k: errs[i] for k,i in idxs.items() }
    return renorm, renorm_err, fps, fes

cff_map = {
    "H":  {"flag": None},
    "Ht": {"flag": None},
    "E":  {"flag": None},
    "Et": {"flag": None},
}
for c in cff_map:
    cff_map[c]["flag_rga"] = flags_rga[c]
    cff_map[c]["flag_rgk"] = flags_rgk[c]
    cff_map[c]["flag"]     = flags_rgk[c]  # we always plot RGK in single-file

renorm_rga, renerr_rga, fit_rga, err_rga = build_param_dicts(pnames_rga, vals_rga, errs_rga)
renorm_rgk, renerr_rgk, fit_rgk, err_rgk = build_param_dicts(pnames_rgk, vals_rgk, errs_rgk)

# ─── Default VGG shape & factor ────────────────────────────────────────────────
defaults = {
    "H":  dict(r=0.9,   alpha0=0.43, alpha1=0.85, n=1.35, b=0.4,  M2=0.64, P=1.0, factor=2.0),
    "Ht": dict(r=7.0,   alpha0=0.43, alpha1=0.85, n=0.6,  b=2.0,  M2=0.8,  P=1.0, factor=0.4),
    "E":  dict(r=0.9,   alpha0=0.43, alpha1=0.85, n=1.35, b=0.4,  M2=0.64, P=1.0, factor=1.0),
    "Et": dict(r=0.9,   alpha0=0.43, alpha1=0.85, n=0.0,  b=0.4,  M2=0.64, P=1.0, factor=1.0),
}

# ─── Im‐CFF factory ─────────────────────────────────────────────────────────────
def make_Im_func(cff, params, renorm):
    d = defaults[cff]
    def Im(xi, t):
        r     = params.get("r",      d["r"])
        a0    = params.get("alpha0", d["alpha0"])
        a1    = params.get("alpha1", d["alpha1"])
        nval  = params.get("n",      d["n"])
        bval  = params.get("b",      d["b"])
        M2    = params.get("M2",     d["M2"])
        Pval  = params.get("P",      d["P"])
        fac   = d["factor"]
        alpha = a0 + a1 * t
        pref  = np.pi * 5/9 * nval * r / (1 + xi)
        xfac  = (2*xi/(1+xi))**(-alpha)
        yfac  = ((1 - xi)/(1+xi))**(bval)
        tfac  = (1 - ((1 - xi)/(1+xi))*t/M2)**(-Pval)
        return renorm * pref * xfac * yfac * tfac * fac
    return Im

# ─── Replica‐band helper ────────────────────────────────────────────────────────
def generate_replicas(central, errors, n=1000):
    reps = []
    for _ in range(n):
        r = {}
        for k,v in central.items():
            sigma = errors[k]/1.96
            r[k] = np.random.normal(v, sigma)
        reps.append(r)
    return reps

def compute_band(cff, Im_rep, xi_vals, t_val):
    xi_a = np.atleast_1d(xi_vals)
    t_a  = np.atleast_1d(t_val)
    xi_b, t_b = np.broadcast_arrays(xi_a, -t_a)
    curves = []
    for params_rep, renorm_rep in zip(
        generate_replicas(fit_rgk[cff], err_rgk[cff], Nrep),
        np.random.normal(renorm_rgk, renerr_rgk/1.96, Nrep)
    ):
        Imf = make_Im_func(cff, params_rep, renorm_rep)
        curves.append(Imf(xi_b, t_b))
    arr = np.array(curves)
    return np.median(arr,0), np.percentile(arr,2.5,0), np.percentile(arr,97.5,0)

# ─── Plot settings ─────────────────────────────────────────────────────────────
plt.style.use("classic")
plt.rcParams.update({'font.size':14, 'font.family':'serif'})

outdir = "output/plots"
os.makedirs(outdir, exist_ok=True)

t_fixed  = np.linspace(0.1, 1.0, 6)
xi_fixed = np.linspace(0.05, 0.50, 6)
xi_range = np.linspace(0.0,  0.5, 200)
t_range  = np.linspace(0.0,  1.0, 200)
Nrep     = 500

orig_style = {'color':'tab:blue','linestyle':'-','linewidth':2.5}
fit_style  = {'color':'tab:red', 'linestyle':'--','linewidth':2.5}
comp_style = {'color':'tab:green','linestyle':':','linewidth':2.5}
band_style = {'color':'tab:red','alpha':0.2}
band2_style= {'color':'tab:green','alpha':0.2}
zero_line  = {'color':'gray','linestyle':'--','linewidth':1}

tex_map = {"H":"H","Ht":r"\tilde H","E":"E","Et":r"\tilde E"}

# ─── Single‐file “RGK only” plots ───────────────────────────────────────────────
for cff in ("H","Ht","E","Et"):
    if not flags_rgk[cff]: continue

    Im_rgk = make_Im_func(cff, fit_rgk[cff], renorm_rgk)
    Im_orig = make_Im_func(cff, {}, 1.0)
    tex = tex_map[cff]

    # Im vs ξ
    fig, axes = plt.subplots(2,3,figsize=(12,8), sharex=True, sharey=True)
    axes = axes.flatten()
    fig.suptitle(rf"$\mathrm{{Im}}\,{tex}$  (RGK preliminary)", fontsize=16, y=0.98)

    for ax, t in zip(axes, t_fixed):
        ax.plot(xi_range, Im_rgk(xi_range, -t), **fit_style)
        med, lo, up = compute_band(cff, Im_rgk, xi_range, t)
        ax.fill_between(xi_range, lo, up, **band_style)
        ax.axhline(0, **zero_line)
        ax.set_xlim(0,0.5); ax.set_ylim(-2,20)
        if ax in axes[:3]:
            ax.set_yticklabels(ax.get_yticks()[1:])  # drop the -2 label
    axes[2].legend(handles=[
        Line2D([0],[0], **fit_style, label="RGK preliminary")
    ], loc='upper right')
    fig.text(0.06,0.5, rf"$\mathrm{{Im}}\,{tex}(\xi,\,-t)$", va='center', ha='center', rotation='vertical')
    fig.text(0.5,0.02, r"$\xi$", ha='center')
    fig.subplots_adjust(left=0.1,right=0.98,bottom=0.08,top=0.92,wspace=0,hspace=0)
    fig.savefig(f"{outdir}/Im{cff}_vs_xi_{ts_rgk}.pdf", bbox_inches='tight')
    plt.close(fig)

    # Im vs −t
    fig, axes = plt.subplots(2,3,figsize=(12,8), sharex=True, sharey=True)
    axes = axes.flatten()
    fig.suptitle(rf"$\mathrm{{Im}}\,{tex}$  (RGK preliminary)", fontsize=16, y=0.98)

    for ax, xi0 in zip(axes, xi_fixed):
        ax.plot(t_range, Im_rgk(xi0, -t_range), **fit_style)
        med, lo, up = compute_band(cff, Im_rgk, xi0, t_range)
        ax.fill_between(t_range, lo, up, **band_style)
        ax.axhline(0, **zero_line)
        ax.set_xlim(0,1.0); ax.set_ylim(-2,20)
        if ax in axes[:3]:
            ax.set_yticklabels(ax.get_yticks()[1:])
    axes[2].legend(handles=[
        Line2D([0],[0], **fit_style, label="RGK preliminary")
    ], loc='upper right')
    fig.text(0.06,0.5, rf"$\mathrm{{Im}}\,{tex}(\xi,\,-t)$", va='center', ha='center', rotation='vertical')
    fig.text(0.5,0.02, r"$-t\,(\mathrm{GeV^2})$", ha='center')
    fig.subplots_adjust(left=0.1,right=0.98,bottom=0.08,top=0.92,wspace=0,hspace=0)
    fig.savefig(f"{outdir}/Im{cff}_vs_t_{ts_rgk}.pdf", bbox_inches='tight')
    plt.close(fig)

# ─── Comparison plots (RGA vs RGK) ─────────────────────────────────────────────
for cff in ("H","Ht","E","Et"):
    if not (flags_rga[cff] and flags_rgk[cff]): continue

    Im_rga = make_Im_func(cff, fit_rga[cff], renorm_rga)
    Im_rgk = make_Im_func(cff, fit_rgk[cff], renorm_rgk)
    tex = tex_map[cff]

    # Im vs ξ
    fig, axes = plt.subplots(2,3,figsize=(12,8), sharex=True, sharey=True)
    axes = axes.flatten()
    fig.suptitle(rf"$\mathrm{{Im}}\,{tex}$  (RGA vs RGK)", fontsize=16, y=0.98)

    for ax, t in zip(axes, t_fixed):
        ax.plot(xi_range, Im_rga(xi_range, -t), **orig_style)
        m1,l1,u1 = compute_band(cff, Im_rga, xi_range, t)
        ax.fill_between(xi_range, l1, u1, **{'color':'tab:blue','alpha':0.1})
        ax.plot(xi_range, Im_rgk(xi_range, -t), **fit_style)
        m2,l2,u2 = compute_band(cff, Im_rgk, xi_range, t)
        ax.fill_between(xi_range, l2, u2, **band_style)
        ax.axhline(0, **zero_line)
        ax.set_xlim(0,0.5); ax.set_ylim(-2,20)
        if ax in axes[:3]:
            ax.set_yticklabels(ax.get_yticks()[1:])
    axes[2].legend(handles=[
        Line2D([0],[0], **orig_style, label="RGA"),
        Line2D([0],[0], **fit_style, label="RGK preliminary")
    ], loc='upper right')
    fig.text(0.06,0.5, rf"$\mathrm{{Im}}\,{tex}(\xi,\,-t)$", va='center', ha='center', rotation='vertical')
    fig.text(0.5,0.02, r"$\xi$", ha='center')
    fig.subplots_adjust(left=0.1,right=0.98,bottom=0.08,top=0.92,wspace=0,hspace=0)
    fig.savefig(f"{outdir}/Im{cff}_vs_xi_{ts_rga}_{ts_rgk}_2.pdf", bbox_inches='tight')
    plt.close(fig)

    # Im vs −t
    fig, axes = plt.subplots(2,3,figsize=(12,8), sharex=True, sharey=True)
    axes = axes.flatten()
    fig.suptitle(rf"$\mathrm{{Im}}\,{tex}$  (RGA vs RGK)", fontsize=16, y=0.98)

    for ax, xi0 in zip(axes, xi_fixed):
        ax.plot(t_range, Im_rga(xi0, -t_range), **orig_style)
        m1,l1,u1 = compute_band(cff, Im_rga, xi0, t_range)
        ax.fill_between(t_range, l1, u1, **{'color':'tab:blue','alpha':0.1})
        ax.plot(t_range, Im_rgk(xi0, -t_range), **fit_style)
        m2,l2,u2 = compute_band(cff, Im_rgk, xi0, t_range)
        ax.fill_between(t_range, l2, u2, **band_style)
        ax.axhline(0, **zero_line)
        ax.set_xlim(0,1.0); ax.set_ylim(-2,20)
        if ax in axes[:3]:
            ax.set_yticklabels(ax.get_yticks()[1:])
    axes[2].legend(handles=[
        Line2D([0],[0], **orig_style, label="RGA"),
        Line2D([0],[0], **fit_style, label="RGK preliminary")
    ], loc='upper right')
    fig.text(0.06,0.5, rf"$\mathrm{{Im}}\,{tex}(\xi,\,-t)$", va='center', ha='center', rotation='vertical')
    fig.text(0.5,0.02, r"$-t\,(\mathrm{GeV^2})$", ha='center')
    fig.subplots_adjust(left=0.1,right=0.98,bottom=0.08,top=0.92,wspace=0,hspace=0)
    fig.savefig(f"{outdir}/Im{cff}_vs_t_{ts_rga}_{ts_rgk}_2.pdf", bbox_inches='tight')
    plt.close(fig)

print("Done plotting.")