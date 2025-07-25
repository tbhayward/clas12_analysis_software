#!/usr/bin/env python3
"""
plot_ImCFFs_two_fits.py

Usage:
    python plot_ImCFFs_two_fits.py fit_results_<TS1>.txt [fit_results_<TS2>.txt]

If only one file is provided, makes “phase 1” plots (just RGK preliminary fit).
If two files are provided, also makes “phase 2” comparison plots.
//
// Each run produces, for each enabled Im CFF (H, Ht, E, Et):
//   • Im{CFF}_vs_xi_<TS1>.pdf           — phase 1 xi‐scan (fit1 only)
//   • Im{CFF}_vs_t_<TS1>.pdf            — phase 1 t‐scan (fit1 only)
//   • Im{CFF}_vs_xi_<TS1>_2.pdf         — phase 2 xi‐scan (fit1 vs fit2)
//   • Im{CFF}_vs_t_<TS1>_2.pdf          — phase 2 t‐scan (fit1 vs fit2)
//
// All y-axes run from –2 to 20; on the top row of each grid the “–2” label is hidden.
"""
import os
import sys
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

# ─── parse args ─────────────────────────────────────────────────────────────────
if not (2 <= len(sys.argv) <= 3):
    print("Usage: python plot_ImCFFs_two_fits.py "
          "fit_results_<TS1>.txt [fit_results_<TS2>.txt]")
    sys.exit(1)
files = sys.argv[1:]
# extract timestamps from filenames
def extract_ts(fname):
    m = re.search(r'fit_results_(\d{8}_\d{6})\.txt$', fname)
    if not m:
        raise RuntimeError(f"Can't parse timestamp from '{fname}'")
    return m.group(1)
ts1 = extract_ts(files[0])
ts2 = extract_ts(files[1]) if len(files) == 2 else None

# ─── loader ─────────────────────────────────────────────────────────────────────
def parse_fit(fname):
    lines = [l.strip() for l in open(fname) if l.strip()]
    # flags
    flag_line = next(l for l in lines if l.startswith("H "))
    toks = flag_line.split()
    flags = { toks[i]: int(toks[i+1]) for i in range(0,len(toks),2) }
    # param names
    pnames = []
    for i,l in enumerate(lines):
        if l.startswith("# parameters"):
            pnames = l.split()[2:]
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
            c, n, cn = lines[i+1].split()
            chi2, ndf, chi2ndf = float(c), int(n), float(cn)
    vals = np.array(vals); errs = np.array(errs)
    return flags, pnames, vals, errs

flags1, pnames1, vals1, errs1 = parse_fit(files[0])
if ts2:
    flags2, pnames2, vals2, errs2 = parse_fit(files[1])

# ─── extract helpers ────────────────────────────────────────────────────────────
def get_idx(name, pnames):
    return pnames.index(name) if name in pnames else None

def unpack(flags, pnames, vals, errs):
    renorm = vals[get_idx("renormImag", pnames)]
    renorm_err = errs[get_idx("renormImag", pnames)]
    fp = {}; fe = {}
    for cff in ("H","Ht","E","Et"):
        if flags[cff]:
            ks = ["r","alpha0","alpha1","n","b","M2","P"]
            idxs = { k:get_idx(f"{k}_{cff}", pnames) for k in ks }
            fp[cff] = { k: vals[i] for k,i in idxs.items() }
            fe[cff] = { k: errs[i] for k,i in idxs.items() }
    return renorm, renorm_err, fp, fe

ren1, ren1_err, fp1, fe1 = unpack(flags1,pnames1,vals1,errs1)
if ts2:
    ren2, ren2_err, fp2, fe2 = unpack(flags2,pnames2,vals2,errs2)

# ─── defaults & factory ─────────────────────────────────────────────────────────
defaults = {
  "H":  dict(r=0.9, alpha0=0.43, alpha1=0.85, n=1.35, b=0.4, M2=0.64, P=1.0, factor=2.0),
  "Ht": dict(r=7.0, alpha0=0.43, alpha1=0.85, n=0.6,  b=2.0, M2=0.8,  P=1.0, factor=0.4),
  "E":  dict(r=0.9, alpha0=0.43, alpha1=0.85, n=1.35, b=0.4, M2=0.64, P=1.0, factor=1.0),
  "Et": dict(r=0.9, alpha0=0.43, alpha1=0.85, n=0.0,  b=0.4, M2=0.64, P=1.0, factor=1.0),
}
def make_Im_func(params, renorm):
    def Im(xi, t):
        d = defaults[cff]
        r    = params.get("r",      d["r"])
        a0   = params.get("alpha0", d["alpha0"])
        a1   = params.get("alpha1", d["alpha1"])
        nval = params.get("n",      d["n"])
        bval = params.get("b",      d["b"])
        M2   = params.get("M2",     d["M2"])
        Pval = params.get("P",      d["P"])
        fac  = d["factor"]
        alpha = a0 + a1*t
        pref  = np.pi*5/9.*nval*r/(1+xi)
        xfac  = (2*xi/(1+xi))**(-alpha)
        yfac  = ((1-xi)/(1+xi))**(bval)
        tfac  = (1 - ((1-xi)/(1+xi))*t/M2)**(-Pval)
        return renorm*pref*xfac*yfac*tfac*fac
    return Im

# ─── replicas & band ────────────────────────────────────────────────────────────
def generate_replicas(central, errors, n=200):
    reps=[]
    for _ in range(n):
        r={}
        for k,v in central.items():
            sigma = errors[k]/1.96
            r[k] = np.random.normal(v, sigma)
        reps.append(r)
    return reps

def compute_band(cff, xi_vals, t_vals, central, errors, renorm, renorm_err, N=200):
    # prepare xi,t grids
    xi = np.atleast_1d(xi_vals)
    t  = np.atleast_1d(t_vals)
    X, T = np.broadcast_arrays(xi, -t)
    reps = generate_replicas(central[cff], errors[cff], N)
    rens = np.random.normal(renorm, renorm_err/1.96, N)
    curves=[]
    for i in range(N):
        Im = make_Im_func(reps[i], rens[i])
        curves.append(Im(X,T))
    A = np.array(curves)
    return np.median(A,0), np.percentile(A,2.5,0), np.percentile(A,97.5,0)

# ─── plotting setup ─────────────────────────────────────────────────────────────
t_fixed  = np.linspace(0.1,1.0,6)
xi_range = np.linspace(0.0,0.5,200)
xi_fixed = np.linspace(0.05,0.5,6)
t_range  = np.linspace(0.0,1.0,200)
orig_style = {'color':'tab:blue','lw':2.5}
fit1_style = {'color':'tab:red','ls':'--','lw':2.5}
fit2_style = {'color':'tab:green','ls':'-.','lw':2.5}
band1 = {'color':'tab:red','alpha':0.2}
band2 = {'color':'tab:green','alpha':0.2}
zero = {'color':'gray','ls':'--','lw':1}
tex = {"H":"H","Ht":r"\tilde H","E":"E","Et":r"\tilde E"}
from matplotlib.lines import Line2D
legend1 = [ Line2D([0],[0],**orig_style, label='Original Model'),
            Line2D([0],[0],**fit1_style,  label='RGK preliminary'),
            plt.Rectangle((0,0),1,1,fc='tab:red',alpha=0.2,label='95% CI') ]
legend2 = [ Line2D([0],[0],**fit1_style, label='RGK preliminary'),
            Line2D([0],[0],**fit2_style, label=os.path.basename(files[1])),
            plt.Rectangle((0,0),1,1,fc='tab:red',alpha=0.2,label='95% CI 1'),
            plt.Rectangle((0,0),1,1,fc='tab:green',alpha=0.2,label='95% CI 2') ]

# ─── phase 1: fit1 only ──────────────────────────────────────────────────────────
for cff in fp1:
    Im1 = make_Im_func(fp1[cff], ren1)
    # -- vs xi --
    fig, axes = plt.subplots(2,3,figsize=(12,8),sharex=True,sharey=True)
    axes=axes.flatten()
    fig.suptitle(rf"$\mathrm{{Im}}\,{tex[cff]}$ — RGK preliminary",y=0.98,fontsize=16)
    for ax,t in zip(axes,t_fixed):
        m,lo,hi = compute_band(cff, xi_range, t,
                               fp1, fe1, ren1, ren1_err, N=500)
        ax.plot(xi_range,m, **fit1_style)
        ax.fill_between(xi_range,lo,hi, **band1)
        ax.axhline(0,**zero)
        ax.set_xlim(0,0.5); ax.set_ylim(-2,20)
    # hide “-2” on top row
    for ax in axes[:3]:
        yt = ax.get_yticklabels()
        if yt: yt[0].set_visible(False)
    axes[2].legend(handles=legend1,loc='upper right')
    fig.text(0.06,0.5,rf"$\mathrm{{Im}}\,{tex[cff]}(\xi,\,-t)$",
             va='center',ha='center',rotation='vertical')
    fig.text(0.5,0.02,r"$\xi$",ha='center')
    fig.tight_layout(rect=[0,0,1,0.95])
    fig.savefig(f"Im{cff}_vs_xi_{ts1}.pdf")
    plt.close(fig)

    # -- vs t --
    fig, axes = plt.subplots(2,3,figsize=(12,8),sharex=True,sharey=True)
    axes=axes.flatten()
    fig.suptitle(rf"$\mathrm{{Im}}\,{tex[cff]}$ — RGK preliminary",y=0.98,fontsize=16)
    for ax,xi0 in zip(axes,xi_fixed):
        m,lo,hi = compute_band(cff, xi0, t_range,
                               fp1, fe1, ren1, ren1_err, N=500)
        ax.plot(t_range,m, **fit1_style)
        ax.fill_between(t_range,lo,hi, **band1)
        ax.axhline(0,**zero)
        ax.set_xlim(0,1.0); ax.set_ylim(-2,20)
    for ax in axes[:3]:
        yt = ax.get_yticklabels()
        if yt: yt[0].set_visible(False)
    axes[2].legend(handles=legend1,loc='upper right')
    fig.text(0.06,0.5,rf"$\mathrm{{Im}}\,{tex[cff]}(\xi,\,-t)$",
             va='center',ha='center',rotation='vertical')
    fig.text(0.5,0.02,r"$-t\,(\mathrm{GeV^2})$",ha='center')
    fig.tight_layout(rect=[0,0,1,0.95])
    fig.savefig(f"Im{cff}_vs_t_{ts1}.pdf")
    plt.close(fig)

# ─── phase 2: comparison fit1 vs fit2 ────────────────────────────────────────────
if ts2:
    for cff in sorted(set(fp1)|set(fp2)):
        fig,axes = plt.subplots(2,3,figsize=(12,8),sharex=True,sharey=True)
        axes=axes.flatten()
        fig.suptitle(rf"$\mathrm{{Im}}\,{tex[cff]}$ — comparison",y=0.98,fontsize=16)

        # vs xi
        for ax,t in zip(axes,t_fixed):
            if cff in fp1:
                m1,lo1,hi1 = compute_band(cff, xi_range, t,
                                         fp1, fe1, ren1, ren1_err, N=500)
                ax.plot(xi_range,m1, **fit1_style)
                ax.fill_between(xi_range,lo1,hi1, **band1)
            if cff in fp2:
                m2,lo2,hi2 = compute_band(cff, xi_range, t,
                                         fp2, fe2, ren2, ren2_err, N=500)
                ax.plot(xi_range,m2, **fit2_style)
                ax.fill_between(xi_range,lo2,hi2, **band2)
            ax.axhline(0,**zero)
            ax.set_xlim(0,0.5); ax.set_ylim(-2,20)
        for ax in axes[:3]:
            yt = ax.get_yticklabels()
            if yt: yt[0].set_visible(False)
        axes[2].legend(handles=legend2,loc='upper right')
        fig.text(0.06,0.5,rf"$\mathrm{{Im}}\,{tex[cff]}(\xi,\,-t)$",
                 va='center',ha='center',rotation='vertical')
        fig.text(0.5,0.02,r"$\xi$",ha='center')
        fig.tight_layout(rect=[0,0,1,0.95])
        fig.savefig(f"Im{cff}_vs_xi_{ts1}_2.pdf")
        plt.close(fig)

        # vs t
        fig,axes = plt.subplots(2,3,figsize=(12,8),sharex=True,sharey=True)
        axes=axes.flatten()
        fig.suptitle(rf"$\mathrm{{Im}}\,{tex[cff]}$ — comparison",y=0.98,fontsize=16)
        for ax,xi0 in zip(axes,xi_fixed):
            if cff in fp1:
                m1,lo1,hi1 = compute_band(cff, xi0, t_range,
                                         fp1, fe1, ren1, ren1_err, N=500)
                ax.plot(t_range,m1, **fit1_style)
                ax.fill_between(t_range,lo1,hi1, **band1)
            if cff in fp2:
                m2,lo2,hi2 = compute_band(cff, xi0, t_range,
                                         fp2, fe2, ren2, ren2_err, N=500)
                ax.plot(t_range,m2, **fit2_style)
                ax.fill_between(t_range,lo2,hi2, **band2)
            ax.axhline(0,**zero)
            ax.set_xlim(0,1.0); ax.set_ylim(-2,20)
        for ax in axes[:3]:
            yt = ax.get_yticklabels()
            if yt: yt[0].set_visible(False)
        axes[2].legend(handles=legend2,loc='upper right')
        fig.text(0.06,0.5,rf"$\mathrm{{Im}}\,{tex[cff]}(\xi,\,-t)$",
                 va='center',ha='center',rotation='vertical')
        fig.text(0.5,0.02,r"$-t\,(\mathrm{GeV^2})$",ha='center')
        fig.tight_layout(rect=[0,0,1,0.95])
        fig.savefig(f"Im{cff}_vs_t_{ts1}_2.pdf")
        plt.close(fig)

print("Done.")