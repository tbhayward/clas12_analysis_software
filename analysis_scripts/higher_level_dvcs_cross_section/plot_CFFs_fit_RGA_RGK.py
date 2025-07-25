#!/usr/bin/env python3
"""
plot_CFFs_fit_RGA_RGK.py

Usage:
    python plot_CFFs_fit_RGA_RGK.py fit_results_<TIMESTAMP1>.txt [fit_results_<TIMESTAMP2>.txt]

Generates Im CFF vs ξ and vs −t plots for the first file (RGK preliminary only),
and if a second file is provided, overplots that comparison as well.
Y‐axes run from −2 to 20, with the bottom tick suppressed on the top row.
"""
import os
import sys
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

# ─── Parse command‐line ─────────────────────────────────────────────────────────
if not (2 <= len(sys.argv) <= 3):
    print("Usage: python plot_CFFs_fit_RGA_RGK.py "
          "fit_results_<TIMESTAMP1>.txt [fit_results_<TIMESTAMP2>.txt]")
    sys.exit(1)

file1 = sys.argv[1]
file2 = sys.argv[2] if len(sys.argv) == 3 else None

def extract_timestamp(fn):
    m = re.search(r'fit_results_(\d{8}_\d{6})\.txt$', fn)
    if not m:
        raise ValueError(f"Cannot parse timestamp from {fn!r}")
    return m.group(1)

ts1 = extract_timestamp(file1)
ts2 = extract_timestamp(file2) if file2 else None

# ─── Load fit results & flags ───────────────────────────────────────────────────
def parse_fit_results(fname):
    with open(fname) as f:
        lines = [l.strip() for l in f if l.strip()]
    # flags line "H 1  Ht 1  E 0  Et 1"
    flag_line = next(l for l in lines if l.startswith("H "))
    toks = flag_line.split()
    flags = { toks[i]: int(toks[i+1]) for i in range(0, len(toks), 2) }
    # parameter names
    pnames = []
    for l in lines:
        if l.startswith("# parameters"):
            pnames = l.split()[2:]
            break
    # values & errors
    vals = errs = None
    for i,l in enumerate(lines):
        if l.startswith("# values"):
            vals = list(map(float, lines[i+1].split()))
        if l.startswith("# errors"):
            errs = list(map(float, lines[i+1].split()))
    if vals is None or errs is None:
        raise RuntimeError("Failed to parse values/errors")
    return flags, pnames, np.array(vals), np.array(errs)

flags1, pnames1, vals1, errs1 = parse_fit_results(file1)
if file2:
    flags2, pnames2, vals2, errs2 = parse_fit_results(file2)

def get_idx(pnames, name):
    return pnames.index(name) if name in pnames else None

# ─── Extract fit‐parameters ─────────────────────────────────────────────────────
renorm1 = vals1[get_idx(pnames1, "renormImag")]
renorm_err1 = errs1[get_idx(pnames1, "renormImag")]
fitp1 = {}
fite1 = {}
for cff in ("H","Ht","E","Et"):
    if flags1[cff]:
        keys = ["r","alpha0","alpha1","n","b","M2","P"]
        idxs = { k: get_idx(pnames1, f"{k}_{cff}") for k in keys }
        fitp1[cff] = { k: vals1[i] for k,i in idxs.items() }
        fite1[cff] = { k: errs1[i] for k,i in idxs.items() }

if file2:
    renorm2 = vals2[get_idx(pnames2, "renormImag")]
    renorm_err2 = errs2[get_idx(pnames2, "renormImag")]
    fitp2 = {}
    fite2 = {}
    for cff in ("H","Ht","E","Et"):
        if flags2[cff]:
            keys = ["r","alpha0","alpha1","n","b","M2","P"]
            idxs = { k: get_idx(pnames2, f"{k}_{cff}") for k in keys }
            fitp2[cff] = { k: vals2[i] for k,i in idxs.items() }
            fite2[cff] = { k: errs2[i] for k,i in idxs.items() }

# ─── VGG defaults + factors ─────────────────────────────────────────────────────
defaults = {
    "H":  dict(r=0.9, alpha0=0.43, alpha1=0.85, n=1.35, b=0.4,  M2=0.64, P=1.0, factor=2.0),
    "Ht": dict(r=7.0, alpha0=0.43, alpha1=0.85, n=0.6,  b=2.0,  M2=0.8,  P=1.0, factor=0.4),
    "E":  dict(r=0.9, alpha0=0.43, alpha1=0.85, n=1.35, b=0.4,  M2=0.64, P=1.0, factor=1.0),
    "Et": dict(r=0.9,alpha0=0.43, alpha1=0.85, n=0.0,  b=0.4,  M2=0.64, P=1.0, factor=1.0),
}

# ─── Im‐CFF builder ─────────────────────────────────────────────────────────────
def make_Im_func(cff, params, renorm):
    d = defaults[cff]
    def Im(xi, t):
        r    = params.get("r",      d["r"])
        a0   = params.get("alpha0", d["alpha0"])
        a1   = params.get("alpha1", d["alpha1"])
        nval = params.get("n",      d["n"])
        bval = params.get("b",      d["b"])
        M2   = params.get("M2",     d["M2"])
        Pval = params.get("P",      d["P"])
        fac  = d["factor"]
        alpha = a0 + a1 * t
        pref  = np.pi*5.0/9.0 * nval * r / (1+xi)
        xfac  = (2*xi/(1+xi))**(-alpha)
        yfac  = ((1 - xi)/(1+xi))**(bval)
        tfac  = (1 - ((1 - xi)/(1+xi))*t/M2)**(-Pval)
        return renorm * pref * xfac * yfac * tfac * fac
    return Im

# ─── Replica helpers ────────────────────────────────────────────────────────────
def generate_replica_params(central, errors, n=100):
    reps = []
    for _ in range(n):
        d = {}
        for k,v in central.items():
            sigma = errors[k]/1.96
            d[k] = np.random.normal(v, sigma)
        reps.append(d)
    return reps

def compute_uncertainty_band(cff, central, errors, renorm, renorm_err,
                             xi_vals, t_vals, n_reps=200):
    # build replicas
    param_reps  = generate_replica_params(central, errors, n_reps)
    renorm_reps = np.random.normal(renorm, renorm_err/1.96, n_reps)
    xi_a = np.atleast_1d(xi_vals)
    t_a  = np.atleast_1d(t_vals)
    xi_b, t_b = np.broadcast_arrays(xi_a, -t_a)
    curves = []
    for i in range(n_reps):
        Im_rep = make_Im_func(cff, param_reps[i], renorm_reps[i])
        curves.append(Im_rep(xi_b, t_b))
    A = np.array(curves)
    return np.median(A,0), np.percentile(A,2.5,0), np.percentile(A,97.5,0)

# ─── Plot setup ─────────────────────────────────────────────────────────────────
plt.style.use('classic')
plt.rcParams.update({'font.size':14,'font.family':'serif'})
outdir = 'output/plots'
os.makedirs(outdir, exist_ok=True)

t_fixed   = np.linspace(0.1,1.0,6)
xi_fixed  = np.linspace(0.05,0.50,6)
xi_range  = np.linspace(0.0,0.5,200)
t_range   = np.linspace(0.0,1.0,200)

fit_style  = {'color':'tab:red','linestyle':'--','linewidth':2.5}
band1      = {'color':'tab:red','alpha':0.2}
comp_style = {'color':'tab:green','linestyle':'--','linewidth':2.5}
band2      = {'color':'tab:green','alpha':0.2}
zero_line  = {'color':'gray','linestyle':'--','linewidth':1}

from matplotlib.lines import Line2D
tex_map = {"H":"H","Ht":r"\tilde H","E":"E","Et":r"\tilde E"}

for cff in ("H","Ht","E","Et"):
    if not flags1[cff]:
        continue

    # single‐file RGK preliminary
    Im1 = make_Im_func(cff, fitp1[cff], renorm1)

    # — Im vs ξ (fix t) —
    fig, axes = plt.subplots(2,3,figsize=(12,8), sharex=True, sharey=True)
    axes = axes.flatten()
    for ax, tval in zip(axes, t_fixed):
        m, lo, up = compute_uncertainty_band(
            cff, fitp1[cff], fite1[cff], renorm1, renorm_err1,
            xi_range, tval
        )
        ax.plot( xi_range, m, **fit_style )
        ax.fill_between( xi_range, lo, up, **band1 )
        ax.axhline(0, **zero_line)
        ax.set_xlim(0,0.5)
        ax.set_ylim(-2,20)
        ax.set_xlabel(r"$\xi$")
        ax.set_title(rf"$-t={tval:.2f}\,$GeV$^2$")
    # remove the "-2" tick on top row
    for ax in axes[:3]:
        yt = ax.get_yticks()[1:]
        ax.set_yticks(yt)
    axes[2].legend([Line2D([0],[0],**fit_style, label="RGK preliminary")],
                   loc='upper right')
    fig.tight_layout(rect=[0,0,1,0.95])
    fig.savefig(f"{outdir}/Im{cff}_vs_xi_{ts1}.pdf", bbox_inches='tight')
    plt.close(fig)

    # — Im vs −t (fix ξ) —
    fig, axes = plt.subplots(2,3,figsize=(12,8), sharex=True, sharey=True)
    axes = axes.flatten()
    for ax, xi0 in zip(axes, xi_fixed):
        m, lo, up = compute_uncertainty_band(
            cff, fitp1[cff], fite1[cff], renorm1, renorm_err1,
            xi0, t_range
        )
        ax.plot( t_range, m, **fit_style )
        ax.fill_between( t_range, lo, up, **band1 )
        ax.axhline(0, **zero_line)
        ax.set_xlim(0,1.0)
        ax.set_ylim(-2,20)
        ax.set_xlabel(r"$-t\ (\mathrm{GeV}^2)$")
        ax.set_title(rf"$\xi={xi0:.2f}$")
    for ax in axes[:3]:
        yt = ax.get_yticks()[1:]
        ax.set_yticks(yt)
    axes[2].legend([Line2D([0],[0],**fit_style, label="RGK preliminary")],
                   loc='upper right')
    fig.tight_layout(rect=[0,0,1,0.95])
    fig.savefig(f"{outdir}/Im{cff}_vs_t_{ts1}.pdf", bbox_inches='tight')
    plt.close(fig)

    # comparison if file2 provided
    if file2 and flags2[cff]:
        # — Im vs ξ comparison —
        fig, axes = plt.subplots(2,3,figsize=(12,8), sharex=True, sharey=True)
        axes = axes.flatten()
        for ax, tval in zip(axes, t_fixed):
            m1,lo1,up1 = compute_uncertainty_band(
                cff, fitp1[cff], fite1[cff], renorm1, renorm_err1,
                xi_range, tval
            )
            m2,lo2,up2 = compute_uncertainty_band(
                cff, fitp2[cff], fite2[cff], renorm2, renorm_err2,
                xi_range, tval
            )
            ax.plot( xi_range, m1, **fit_style )
            ax.fill_between( xi_range, lo1, up1, **band1 )
            ax.plot( xi_range, m2, **comp_style )
            ax.fill_between( xi_range, lo2, up2, **band2 )
            ax.axhline(0, **zero_line)
            ax.set_xlim(0,0.5)
            ax.set_ylim(-2,20)
            ax.set_xlabel(r"$\xi$")
            ax.set_title(rf"$-t={tval:.2f}\,$GeV$^2$")
        for ax in axes[:3]:
            yt = ax.get_yticks()[1:]
            ax.set_yticks(yt)
        axes[2].legend([
            Line2D([0],[0],**fit_style,  label="RGK preliminary"),
            Line2D([0],[0],**comp_style, label="Comparison")
        ], loc='upper right')
        fig.tight_layout(rect=[0,0,1,0.95])
        fig.savefig(f"{outdir}/Im{cff}_vs_xi_{ts1}_2.pdf", bbox_inches='tight')
        plt.close(fig)

        # — Im vs −t comparison —
        fig, axes = plt.subplots(2,3,figsize=(12,8), sharex=True, sharey=True)
        axes = axes.flatten()
        for ax, xi0 in zip(axes, xi_fixed):
            m1,lo1,up1 = compute_uncertainty_band(
                cff, fitp1[cff], fite1[cff], renorm1, renorm_err1,
                xi0, t_range
            )
            m2,lo2,up2 = compute_uncertainty_band(
                cff, fitp2[cff], fite2[cff], renorm2, renorm_err2,
                xi0, t_range
            )
            ax.plot( t_range, m1, **fit_style )
            ax.fill_between( t_range, lo1, up1, **band1 )
            ax.plot( t_range, m2, **comp_style )
            ax.fill_between( t_range, lo2, up2, **band2 )
            ax.axhline(0, **zero_line)
            ax.set_xlim(0,1.0)
            ax.set_ylim(-2,20)
            ax.set_xlabel(r"$-t\ (\mathrm{GeV}^2)$")
            ax.set_title(rf"$\xi={xi0:.2f}$")
        for ax in axes[:3]:
            yt = ax.get_yticks()[1:]
            ax.set_yticks(yt)
        axes[2].legend([
            Line2D([0],[0],**fit_style,  label="RGK preliminary"),
            Line2D([0],[0],**comp_style, label="Comparison")
        ], loc='upper right')
        fig.tight_layout(rect=[0,0,1,0.95])
        fig.savefig(f"{outdir}/Im{cff}_vs_t_{ts1}_2.pdf", bbox_inches='tight')
        plt.close(fig)

print("Done plotting.")