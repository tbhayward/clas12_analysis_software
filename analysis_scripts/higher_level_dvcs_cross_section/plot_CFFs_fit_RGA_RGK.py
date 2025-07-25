#!/usr/bin/env python3
"""
plot_ImCFFs_two_fits.py

Usage:
    python plot_ImCFFs_two_fits.py fit_results_<TS1>.txt fit_results_<TS2>.txt

1) From the *first* file, make
     ImCFF vs ξ  (2×3 grid)   →  only the fitted median+95%CI, label “RGK preliminary”
     ImCFF vs −t (2×3 grid)
   Files:  Im{CFF}_vs_xi_<TS1>.pdf  and Im{CFF}_vs_t_<TS1>.pdf

2) Make the *comparison* plots (same two grids) overplotting fit1 + fit2
   Files:  Im{CFF}_vs_xi_<TS1>_2.pdf  and Im{CFF}_vs_t_<TS1>_2.pdf
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
    print(__doc__.strip())
    sys.exit(1)

file1, file2 = sys.argv[1], sys.argv[2]
def extract_ts(fname):
    m = re.search(r'fit_results_(\d{8}_\d{6})\.txt$', fname)
    if not m:
        raise ValueError(f"Filename {fname!r} does not match expected pattern")
    return m.group(1)

ts1, ts2 = extract_ts(file1), extract_ts(file2)

# ─── Common parsing routine ─────────────────────────────────────────────────────
def parse_fit_results(fname):
    with open(fname) as f:
        lines = [l.strip() for l in f if l.strip()]
    # flags line
    flag_line = next(l for l in lines if l.startswith("H "))
    toks = flag_line.split()
    flags = { toks[i]: int(toks[i+1]) for i in range(0, len(toks), 2) }
    # parameter names
    for i,l in enumerate(lines):
        if l.startswith("# parameters"):
            pnames = lines[i].split()[2:]
            break
    else:
        raise RuntimeError("No '# parameters' line")
    # values/errors
    for i,l in enumerate(lines):
        if l.startswith("# values"):
            vals = np.array(list(map(float, lines[i+1].split())))
        if l.startswith("# errors"):
            errs = np.array(list(map(float, lines[i+1].split())))
        if l.startswith("# chi2"):
            parts = lines[i+1].split()
            chi2, ndf, chi2ndf = float(parts[0]), int(parts[1]), float(parts[2])
    # renormImag
    def idx(n): return pnames.index(n)
    renorm = vals[idx("renormImag")]
    renorm_err = errs[idx("renormImag")]
    # collect fit_params & errors
    fit_p, fit_e = {}, {}
    for cff in ("H","Ht","E","Et"):
        if flags[cff]:
            keys = ["r","alpha0","alpha1","n","b","M2","P"]
            ip = {}
            for k in keys:
                name = f"{k}_{cff}"
                ip[k] = vals[idx(name)], errs[idx(name)]
            fit_p[cff] = {k:iv for k,(iv,ie) in ip.items()}
            fit_e[cff] = {k:ie for k,(iv,ie) in ip.items()}
    return flags, fit_p, fit_e, renorm, renorm_err

flags1,  fit1, err1, ren1, ren1_err = parse_fit_results(file1)
flags2,  fit2, err2, ren2, ren2_err = parse_fit_results(file2)

# ─── VGG defaults & Im‐builder ──────────────────────────────────────────────────
defaults = {
    "H":  dict(r=0.9, alpha0=0.43, alpha1=0.85, n=1.35, b=0.4,  M2=0.64, P=1.0, factor=2.0),
    "Ht": dict(r=7.0, alpha0=0.43, alpha1=0.85, n=0.6,  b=2.0,  M2=0.8,  P=1.0, factor=0.4),
    "E":  dict(r=0.9, alpha0=0.43, alpha1=0.85, n=1.35, b=0.4,  M2=0.64, P=1.0, factor=1.0),
    "Et": dict(r=0.9, alpha0=0.43, alpha1=0.85, n=0.0,  b=0.4,  M2=0.64, P=1.0, factor=1.0),
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

        alpha = a0 + a1*t
        pref  = np.pi*5/9 * nval * r/(1+xi)
        xfac  = (2*xi/(1+xi))**(-alpha)
        yfac  = ((1-xi)/(1+xi))**(bval)
        tfac  = (1 - ((1-xi)/(1+xi))*t/M2)**(-Pval)
        return renorm * pref*xfac*yfac*tfac*fac
    return Im

def generate_replicas(central, errors, n=1000):
    reps=[]
    for _ in range(n):
        r={}
        for k,v in central.items():
            sigma=errors[k]/1.96
            r[k]=np.random.normal(v,sigma)
        reps.append(r)
    return reps

def compute_band(cff, Im_func, xi_vals, t_vals, nrep=1000):
    xi_a = np.atleast_1d(xi_vals)
    t_a  = np.atleast_1d(t_vals)
    xi_b, t_b = np.broadcast_arrays(xi_a, -t_a)
    arr = np.stack([ Im_func(xi_b, t_b) for _ in range(1) ],0)  # dummy
    # for bands we assume Im_func carries its own replica set
    # but here we just return np.nan band for single‐fit mode
    return Im_func(xi_b, t_b), None, None

def compute_band_rep(cff, central, errors, renorm, renorm_err, xi, t, nrep=500):
    xi_a = np.atleast_1d(xi)
    t_a  = np.atleast_1d(t)
    xi_b, t_b = np.broadcast_arrays(xi_a, -t_a)

    param_reps  = generate_replicas(central, errors, nrep)
    ren_reps    = np.random.normal(renorm, renorm_err/1.96, nrep)
    allC = np.zeros((nrep,)+xi_b.shape)
    for i in range(nrep):
        f = make_Im_func(cff, param_reps[i], ren_reps[i])
        allC[i] = f(xi_b, t_b)
    m = np.median(allC,0)
    lo = np.percentile(allC,2.5,0)
    hi = np.percentile(allC,97.5,0)
    return m, lo, hi

# ─── Plot settings ───────────────────────────────────────────────────────────────
outdir = "output/plots"
os.makedirs(outdir,exist_ok=True)

t_fixed   = np.linspace(0.1,1.0,6)
xi_range  = np.linspace(0.0,0.5,200)
xi_fixed  = np.linspace(0.05,0.50,6)
t_range   = np.linspace(0.0,1.0,200)

orig_style = {'color':'tab:blue','ls':'-','lw':2.5}
fit1_style= {'color':'tab:red',  'ls':'--','lw':2.5}
fit2_style= {'color':'tab:green','ls':'-.','lw':2.5}
band1     = {'color':'tab:red',   'alpha':0.2}
band2     = {'color':'tab:green', 'alpha':0.2}
zero_line = {'color':'gray','ls':'--','lw':1}

tex_map = {"H":"H","Ht":r"\tilde H","E":"E","Et":r"\tilde E"}
Nrep = 500

# ─── Make Im_funcs for each fit ─────────────────────────────────────────────────
Im1 = { c: make_Im_func(c, fit1[c], ren1) for c in fit1 }
Im2 = { c: make_Im_func(c, fit2[c], ren2) for c in fit2 }

# ─── Phase 1: ONLY fit1 (no original model) ────────────────────────────────────
for cff, f1 in Im1.items():
    tex = tex_map[cff]
    # Im vs ξ
    fig,axes = plt.subplots(2,3,figsize=(12,8),sharex=True,sharey=True)
    axes=axes.flatten()
    fig.suptitle(rf"$\mathrm{{Im}}\,{tex}$ — RGK preliminary",y=0.98,fontsize=16)
    for ax,t in zip(axes,t_fixed):
        m1,lo1,hi1 = compute_band_rep(cff, fit1[cff], err1[cff], ren1, ren1_err, xi_range, t, Nrep)
        ax.plot(xi_range,m1, **fit1_style, label="RGK preliminary")
        ax.fill_between(xi_range, lo1, hi1, **band1)
        ax.axhline(0,**zero_line)
        ax.set_xlim(0,0.5); ax.set_ylim(0,12)
        ax.text(0.6,0.7,rf"$-t={t:.2f}$",transform=ax.transAxes)
    axes[2].legend(loc='upper right')
    fig.text(0.06,0.5, rf"$\mathrm{{Im}}\,{tex}(\xi,\,-t)$",va='center',ha='center',rotation='vertical')
    fig.text(0.5,0.02,r"$\xi$",ha='center')
    fig.tight_layout(rect=[0,0,1,0.95])
    fig.savefig(f"{outdir}/Im{cff}_vs_xi_{ts1}.pdf",bbox_inches='tight')
    plt.close(fig)

    # Im vs −t
    fig,axes = plt.subplots(2,3,figsize=(12,8),sharex=True,sharey=True)
    axes=axes.flatten()
    fig.suptitle(rf"$\mathrm{{Im}}\,{tex}$ — RGK preliminary",y=0.98,fontsize=16)
    for ax,xi0 in zip(axes,xi_fixed):
        m1,lo1,hi1 = compute_band_rep(cff, fit1[cff], err1[cff], ren1, ren1_err, xi0, t_range, Nrep)
        ax.plot(t_range,m1, **fit1_style)
        ax.fill_between(t_range, lo1, hi1, **band1)
        ax.axhline(0,**zero_line)
        ax.set_xlim(0,1.0); ax.set_ylim((0,12) if xi0>0.2 else (-2,12))
        ax.text(0.6,0.7,rf"$\xi={xi0:.2f}$",transform=ax.transAxes)
    axes[2].legend(loc='upper right')
    fig.text(0.06,0.5, rf"$\mathrm{{Im}}\,{tex}(\xi,\,-t)$",va='center',ha='center',rotation='vertical')
    fig.text(0.5,0.02,r"$-t\,(\mathrm{GeV^2})$",ha='center')
    fig.tight_layout(rect=[0,0,1,0.95])
    fig.savefig(f"{outdir}/Im{cff}_vs_t_{ts1}.pdf",bbox_inches='tight')
    plt.close(fig)

# ─── Phase 2: compare fit1 vs fit2 ──────────────────────────────────────────────
for cff in sorted(set(Im1)|set(Im2)):
    tex = tex_map[cff]
    fig,axes = plt.subplots(2,3,figsize=(12,8),sharex=True,sharey=True)
    axes=axes.flatten()
    fig.suptitle(rf"$\mathrm{{Im}}\,{tex}$ — comparison",y=0.98,fontsize=16)

    for ax,t in zip(axes,t_fixed):
        if cff in Im1:
            m1,lo1,hi1 = compute_band_rep(cff, fit1[cff], err1[cff], ren1, ren1_err, xi_range, t, Nrep)
            ax.plot(xi_range,m1, **fit1_style, label="RGK preliminary")
            ax.fill_between(xi_range, lo1, hi1, **band1)
        if cff in Im2:
            m2,lo2,hi2 = compute_band_rep(cff, fit2[cff], err2[cff], ren2, ren2_err, xi_range, t, Nrep)
            ax.plot(xi_range,m2, **fit2_style, label=os.path.basename(file2))
            ax.fill_between(xi_range, lo2, hi2, **band2)

        ax.axhline(0,**zero_line)
        ax.set_xlim(0,0.5); ax.set_ylim(0,12)
        ax.text(0.6,0.7,rf"$-t={t:.2f}$",transform=ax.transAxes)

    axes[2].legend(loc='upper right')
    fig.text(0.06,0.5, rf"$\mathrm{{Im}}\,{tex}(\xi,\,-t)$",va='center',ha='center',rotation='vertical')
    fig.text(0.5,0.02,r"$\xi$",ha='center')
    fig.tight_layout(rect=[0,0,1,0.95])
    fig.savefig(f"{outdir}/Im{cff}_vs_xi_{ts1}_2.pdf",bbox_inches='tight')
    plt.close(fig)

    # vs −t
    fig,axes = plt.subplots(2,3,figsize=(12,8),sharex=True,sharey=True)
    axes=axes.flatten()
    fig.suptitle(rf"$\mathrm{{Im}}\,{tex}$ — comparison",y=0.98,fontsize=16)

    for ax,xi0 in zip(axes,xi_fixed):
        if cff in Im1:
            m1,lo1,hi1 = compute_band_rep(cff, fit1[cff], err1[cff], ren1, ren1_err, xi0, t_range, Nrep)
            ax.plot(t_range,m1, **fit1_style, label="RGK preliminary")
            ax.fill_between(t_range, lo1, hi1, **band1)
        if cff in Im2:
            m2,lo2,hi2 = compute_band_rep(cff, fit2[cff], err2[cff], ren2, ren2_err, xi0, t_range, Nrep)
            ax.plot(t_range,m2, **fit2_style, label=os.path.basename(file2))
            ax.fill_between(t_range, lo2, hi2, **band2)

        ax.axhline(0,**zero_line)
        ax.set_xlim(0,1.0); ax.set_ylim((0,12) if xi0>0.2 else (-2,12))
        ax.text(0.6,0.7,rf"$\xi={xi0:.2f}$",transform=ax.transAxes)

    axes[2].legend(loc='upper right')
    fig.text(0.06,0.5, rf"$\mathrm{{Im}}\,{tex}(\xi,\,-t)$",va='center',ha='center',rotation='vertical')
    fig.text(0.5,0.02,r"$-t\,(\mathrm{GeV^2})$",ha='center')
    fig.tight_layout(rect=[0,0,1,0.95])
    fig.savefig(f"{outdir}/Im{cff}_vs_t_{ts1}_2.pdf",bbox_inches='tight')
    plt.close(fig)

print("All done.")