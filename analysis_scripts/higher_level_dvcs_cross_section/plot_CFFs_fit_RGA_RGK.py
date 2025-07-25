#!/usr/bin/env python3
"""
plot_CFFs_fit_RGA_RGK.py

Usage:
    python plot_CFFs_fit_RGA_RGK.py <fit_results1.txt> [<fit_results2.txt>]

Reads which CFFs were fit from the header of each results file, then for each
enabled Im CFF makes two sets of figures:
  1) “single” — Im CFF vs. ξ and vs. −t for the first fit only
  2) “comparison” — same plots overlaying fit1 (blue) and fit2 (red)

Always shows median curve and 95% CI bands. Y-axes run from 0 to 20 with ticks
at [0, 5, 10, 15, 20], but the top‐row plots omit the “0” tick so it doesn’t
overlap the bottom‐row 20.
"""
import os
import sys
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.stats import norm

# ─── Command‐line ────────────────────────────────────────────────────────────────
if not (3 <= len(sys.argv) <= 4):
    print("Usage: python plot_CFFs_fit_RGA_RGK.py <fit1.txt> [<fit2.txt>]")
    sys.exit(1)
fit1file = sys.argv[1]
fit2file = sys.argv[2] if len(sys.argv) == 3 else None

def parse_timestamp(fn):
    m = re.search(r'fit_results_(\d{8}_\d{6})\.txt$', fn)
    return m.group(1) if m else os.path.splitext(os.path.basename(fn))[0]

ts1 = parse_timestamp(fit1file)
ts2 = parse_timestamp(fit2file) if fit2file else None

# ─── Parser ─────────────────────────────────────────────────────────────────────
def parse_fit_results(fname):
    with open(fname) as f:
        lines = [L.strip() for L in f if L.strip()]
    # flags
    fl = next(L for L in lines if L.startswith("H "))
    toks = fl.split()
    flags = { toks[i]: int(toks[i+1]) for i in range(0, len(toks), 2) }
    # parameter names
    pnames = []
    for i,L in enumerate(lines):
        if L.startswith("# parameters"):
            pnames = lines[i].split()[2:]
            break
    # values & errors
    vals = errs = None
    chi2 = ndf = chi2ndf = None
    for i,L in enumerate(lines):
        if L.startswith("# values"):
            vals = list(map(float, lines[i+1].split()))
        if L.startswith("# errors"):
            errs = list(map(float, lines[i+1].split()))
        if L.startswith("# chi2"):
            a = lines[i+1].split()
            chi2,ndf,chi2ndf = float(a[0]), int(a[1]), float(a[2])
    if vals is None or errs is None:
        raise RuntimeError(f"Could not parse values/errors from {fname}")
    return flags, pnames, np.array(vals), np.array(errs), chi2, ndf, chi2ndf

# load both
flags1, pnames1, vals1, errs1, chi2_1, ndf1, _ = parse_fit_results(fit1file)
if fit2file:
    flags2, pnames2, vals2, errs2, chi2_2, ndf2, _ = parse_fit_results(fit2file)

# ─── Helpers ────────────────────────────────────────────────────────────────────
def get_idx(pnames, name):
    return pnames.index(name) if name in pnames else None

def build_fit_dicts(flags, pnames, vals, errs):
    renorm = vals[get_idx(pnames, "renormImag")]
    renorm_err = errs[get_idx(pnames, "renormImag")]
    fp = {}; fe = {}
    for cff in ("H","Ht","E","Et"):
        if flags[cff]:
            keys = ["r","alpha0","alpha1","n","b","M2","P"]
            idxs = { k: get_idx(pnames, f"{k}_{cff}") for k in keys }
            fp[cff] = { k: vals[i] for k,i in idxs.items() }
            fe[cff] = { k: errs[i] for k,i in idxs.items() }
    return renorm, renorm_err, fp, fe

ren1, ren1_err, fitp1, fite1 = build_fit_dicts(flags1, pnames1, vals1, errs1)
if fit2file:
    ren2, ren2_err, fitp2, fite2 = build_fit_dicts(flags2, pnames2, vals2, errs2)

# defaults + factor
defaults = {
    "H":  dict(r=0.9,   alpha0=0.43, alpha1=0.85, n=1.35, b=0.4,  M2=0.64, P=1.0, factor=2.0),
    "Ht": dict(r=7.0,   alpha0=0.43, alpha1=0.85, n=0.6,  b=2.0,  M2=0.8,  P=1.0, factor=0.4),
    "E":  dict(r=0.9,   alpha0=0.43, alpha1=0.85, n=1.35, b=0.4,  M2=0.64, P=1.0, factor=1.0),
    "Et": dict(r=0.9,   alpha0=0.43, alpha1=0.85, n=0.0,  b=0.4,  M2=0.64, P=1.0, factor=1.0),
}

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
        pref  = np.pi*5.0/9.0 * nval * r / (1 + xi)
        xfac  = (2*xi/(1+xi))**(-alpha)
        yfac  = ((1 - xi)/(1+xi))**(bval)
        tfac  = (1 - ((1 - xi)/(1+xi))*t/M2)**(-Pval)
        return renorm * pref * xfac * yfac * tfac * fac
    return Im

def generate_replica_params(central, errors, n=100):
    reps=[]
    for _ in range(n):
        rc={}
        for k,v in central.items():
            sigma = errors[k]/1.96
            rc[k] = np.random.normal(v, sigma)
        reps.append(rc)
    return reps

def compute_uncertainty_band(cff, xi_vals, t_vals, renorm, central, errors, n_reps=100):
    # always broadcast xi_vals, t_vals
    xi = np.atleast_1d(xi_vals)
    t  = np.atleast_1d(t_vals)
    xi_b, t_b = np.broadcast_arrays(xi, -t)
    reps = generate_replica_params(central[cff], errors[cff], n_reps)
    ren_reps = np.random.normal(renorm, renorm/np.sqrt(len(xi_b)), n_reps)
    allc=[]
    for i in range(n_reps):
        Imf = make_Im_func(cff, reps[i], ren_reps[i])
        allc.append(Imf(xi_b, t_b))
    A = np.array(allc)
    return np.median(A,0), np.percentile(A,2.5,0), np.percentile(A,97.5,0)

# ─── Plot setup ─────────────────────────────────────────────────────────────────
os.makedirs("output/plots", exist_ok=True)
t_fixed  = np.linspace(0.1,1.0,6)
xi_fixed = np.linspace(0.05,0.50,6)
xi_range = np.linspace(0.0,0.5,200)
t_range  = np.linspace(0.0,1.0,200)

plt.style.use('classic')
plt.rcParams.update({'font.size':14,'font.family':'serif'})

# styles
style1 = {'color':'tab:blue','linestyle':'-','linewidth':2.5}
band1  = {'color':'tab:blue','alpha':0.2}
style2 = {'color':'tab:red', 'linestyle':'--','linewidth':2.5}
band2  = {'color':'tab:red','alpha':0.2}
zero   = {'color':'gray','linestyle':'--','linewidth':1}

label1 = os.path.splitext(os.path.basename(fit1file))[0]
label2 = os.path.splitext(os.path.basename(fit2file))[0] if fit2file else None

# ─── Single‐file plots (fit1 only) ───────────────────────────────────────────────
for cff in fitp1:
    Im1 = make_Im_func(cff, fitp1[cff], ren1)
    tex = {"H":"H","Ht":r"\tilde H","E":"E","Et":r"\tilde E"}[cff]

    # Im vs ξ
    fig, axes = plt.subplots(2,3,figsize=(12,8), sharex=True, sharey=True)
    axes = axes.flatten()
    fig.suptitle(rf"$\mathrm{{Im}}\,{tex}$", fontsize=16, y=0.98)
    for ax, t in zip(axes, t_fixed):
        med, lo, up = compute_uncertainty_band(cff, xi_range, t, ren1, fitp1, fite1, 1000)
        ax.plot(xi_range, med,   **style1, label=label1)
        ax.fill_between(xi_range, lo, up, **band1)
        ax.axhline(0, **zero)
        ax.text(0.6,0.75, rf"$-t={t:.2f}$", transform=ax.transAxes)
        ax.set_xlim(0,0.5)
        ax.set_ylim(0,20)
        ax.set_yticks([0,5,10,15,20])
    # drop “0” tick in top row
    for ax in axes[:3]:
        ax.set_yticks([5,10,15,20])
    axes[2].legend(loc='upper right', fontsize=10)
    fig.text(0.06,0.5, rf"$\mathrm{{Im}}\,{tex}(\xi,\,-t)$",
             va='center', ha='center', rotation='vertical')
    fig.text(0.5,0.02, r"$\xi$", ha='center')
    fig.subplots_adjust(left=0.10,right=0.98,top=0.92,bottom=0.08,wspace=0,hspace=0)
    out = f"output/plots/{label1}_Im{cff}_vs_xi_{ts1}.pdf"
    fig.savefig(out, bbox_inches='tight'); print("Saved:", out)
    plt.close(fig)

    # Im vs −t
    fig, axes = plt.subplots(2,3,figsize=(12,8), sharex=True, sharey=True)
    axes = axes.flatten()
    fig.suptitle(rf"$\mathrm{{Im}}\,{tex}$", fontsize=16, y=0.98)
    for ax, xi0 in zip(axes, xi_fixed):
        med, lo, up = compute_uncertainty_band(cff, xi0, t_range, ren1, fitp1, fite1, 1000)
        ax.plot(t_range, med,   **style1)
        ax.fill_between(t_range, lo, up, **band1)
        ax.axhline(0, **zero)
        ax.text(0.6,0.75, rf"$\xi={xi0:.2f}$", transform=ax.transAxes)
        ax.set_xlim(0,1.0)
        ax.set_ylim(0,20)
        ax.set_yticks([0,5,10,15,20])
    for ax in axes[:3]:
        ax.set_yticks([5,10,15,20])
    axes[2].legend(loc='upper right', fontsize=10)
    fig.text(0.06,0.5, rf"$\mathrm{{Im}}\,{tex}(\xi,\,-t)$",
             va='center', ha='center', rotation='vertical')
    fig.text(0.5,0.02, r"$-t\,(\mathrm{GeV^2})$", ha='center')
    fig.subplots_adjust(left=0.10,right=0.98,top=0.92,bottom=0.08,wspace=0,hspace=0)
    out = f"output/plots/{label1}_Im{cff}_vs_t_{ts1}.pdf"
    fig.savefig(out, bbox_inches='tight'); print("Saved:", out)
    plt.close(fig)

# ─── Comparison plots (fit1 + fit2) ──────────────────────────────────────────────
if fit2file:
    for cff in fitp1:
        Im1 = make_Im_func(cff, fitp1[cff], ren1)
        Im2 = make_Im_func(cff, fitp2[cff], ren2)
        tex = {"H":"H","Ht":r"\tilde H","E":"E","Et":r"\tilde E"}[cff]

        # Im vs ξ
        fig, axes = plt.subplots(2,3,figsize=(12,8), sharex=True, sharey=True)
        axes = axes.flatten()
        fig.suptitle(rf"$\mathrm{{Im}}\,{tex}$", fontsize=16, y=0.98)
        for ax, t in zip(axes, t_fixed):
            m1, l1, u1 = compute_uncertainty_band(cff, xi_range, t, ren1, fitp1, fite1, 500)
            m2, l2, u2 = compute_uncertainty_band(cff, xi_range, t, ren2, fitp2, fite2, 500)
            ax.plot(xi_range, m1,   **style1, label=label1)
            ax.fill_between(xi_range, l1, u1, **band1)
            ax.plot(xi_range, m2,   **style2, label=label2)
            ax.fill_between(xi_range, l2, u2, **band2)
            ax.axhline(0, **zero)
            ax.text(0.6,0.75, rf"$-t={t:.2f}$", transform=ax.transAxes)
            ax.set_xlim(0,0.5)
            ax.set_ylim(0,20)
            ax.set_yticks([0,5,10,15,20])
        for ax in axes[:3]:
            ax.set_yticks([5,10,15,20])
        handles = [
            Line2D([0],[0],**style1, label=label1),
            Line2D([0],[0],**style2, label=label2),
            Patch(**band1, label=f"{label1} 95% CI"),
            Patch(**band2, label=f"{label2} 95% CI"),
        ]
        axes[2].legend(handles=handles, loc='upper right', fontsize=10)
        fig.text(0.06,0.5, rf"$\mathrm{{Im}}\,{tex}(\xi,\,-t)$",
                 va='center', ha='center', rotation='vertical')
        fig.text(0.5,0.02, r"$\xi$", ha='center')
        fig.subplots_adjust(left=0.10,right=0.98,top=0.92,bottom=0.08,wspace=0,hspace=0)
        out = f"output/plots/{label1}_vs_{label2}_Im{cff}_vs_xi_{ts1}_{ts2}.pdf"
        fig.savefig(out, bbox_inches='tight'); print("Saved:", out)
        plt.close(fig)

        # Im vs −t
        fig, axes = plt.subplots(2,3,figsize=(12,8), sharex=True, sharey=True)
        axes = axes.flatten()
        fig.suptitle(rf"$\mathrm{{Im}}\,{tex}$", fontsize=16, y=0.98)
        for ax, xi0 in zip(axes, xi_fixed):
            m1, l1, u1 = compute_uncertainty_band(cff, xi0, t_range, ren1, fitp1, fite1, 500)
            m2, l2, u2 = compute_uncertainty_band(cff, xi0, t_range, ren2, fitp2, fite2, 500)
            ax.plot(t_range, m1,   **style1)
            ax.fill_between(t_range, l1, u1, **band1)
            ax.plot(t_range, m2,   **style2)
            ax.fill_between(t_range, l2, u2, **band2)
            ax.axhline(0, **zero)
            ax.text(0.6,0.75, rf"$\xi={xi0:.2f}$", transform=ax.transAxes)
            ax.set_xlim(0,1.0)
            ax.set_ylim(0,20)
            ax.set_yticks([0,5,10,15,20])
        for ax in axes[:3]:
            ax.set_yticks([5,10,15,20])
        handles = [
            Line2D([0],[0],**style1, label=label1),
            Line2D([0],[0],**style2, label=label2),
            Patch(**band1, label=f"{label1} 95% CI"),
            Patch(**band2, label=f"{label2} 95% CI"),
        ]
        axes[2].legend(handles=handles, loc='upper right', fontsize=10)
        fig.text(0.06,0.5, rf"$\mathrm{{Im}}\,{tex}(\xi,\,-t)$",
                 va='center', ha='center', rotation='vertical')
        fig.text(0.5,0.02, r"$-t\,(\mathrm{GeV^2})$", ha='center')
        fig.subplots_adjust(left=0.10,right=0.98,top=0.92,bottom=0.08,wspace=0,hspace=0)
        out = f"output/plots/{label1}_vs_{label2}_Im{cff}_vs_t_{ts1}_{ts2}.pdf"
        fig.savefig(out, bbox_inches='tight'); print("Saved:", out)
        plt.close(fig)

print("All done.")