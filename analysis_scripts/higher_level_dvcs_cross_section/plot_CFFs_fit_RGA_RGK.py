#!/usr/bin/env python3
"""
plot_CFFs_fit_RGA_RGK.py

Usage:
    python plot_CFFs_fit_RGA_RGK.py <fit1.txt> [<fit2.txt>]

Produces “single” RGK‐only plots and (if given) “comparison” RGA vs RGK overlays,
with y∈[0,20], ticks at 0,5,10,15,20 (top row drops the 0), and no “0.0” on
bottom‐middle/right x‐axes. Uses N_REPLICAS=1000 for 95% CI bands.
"""
import os, sys, re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ─── Number of Monte Carlo replicas for uncertainty bands ───────────────────────
N_REPLICAS = 1000

# ─── CLI ────────────────────────────────────────────────────────────────────────
if not (2 <= len(sys.argv) <= 3):
    print("Usage: python plot_CFFs_fit_RGA_RGK.py <fit1.txt> [<fit2.txt>]")
    sys.exit(1)
fit1file = sys.argv[1]
fit2file = sys.argv[2] if len(sys.argv) == 3 else None

def parse_ts(fn):
    m = re.search(r'fit_results_(\d{8}_\d{6})\.txt$', fn)
    return m.group(1) if m else os.path.splitext(fn)[0]

ts1 = parse_ts(fit1file)
ts2 = parse_ts(fit2file) if fit2file else None

# ─── parse_fit_results ──────────────────────────────────────────────────────────
def parse_fit_results(fname):
    L = [l.strip() for l in open(fname) if l.strip()]
    fl = next(l for l in L if l.startswith("H "))
    toks = fl.split()
    flags = { toks[i]: int(toks[i+1]) for i in range(0, len(toks), 2) }
    pnames = []
    for i,l in enumerate(L):
        if l.startswith("# parameters"):
            pnames = L[i].split()[2:]
            break
    vals = errs = None
    chi2 = ndf = None
    for i,l in enumerate(L):
        if l.startswith("# values"):
            vals = list(map(float, L[i+1].split()))
        if l.startswith("# errors"):
            errs = list(map(float, L[i+1].split()))
        if l.startswith("# chi2"):
            parts = L[i+1].split()
            chi2, ndf = float(parts[0]), int(parts[1])
    return flags, pnames, np.array(vals), np.array(errs), chi2, ndf

flags1, pn1, vals1, errs1, chi21, ndf1 = parse_fit_results(fit1file)
if fit2file:
    flags2, pn2, vals2, errs2, chi22, ndf2 = parse_fit_results(fit2file)

def idx(pn, name):
    return pn.index(name)

def build(flags, pn, vals, errs):
    ren    = vals[idx(pn, "renormImag")]
    renerr = errs[idx(pn, "renormImag")]
    fp, fe = {}, {}
    for cff in ("H","Ht","E","Et"):
        if flags[cff]:
            keys = ["r","alpha0","alpha1","n","b","M2","P"]
            ids  = { k: idx(pn, f"{k}_{cff}") for k in keys }
            fp[cff] = { k: vals[i] for k,i in ids.items() }
            fe[cff] = { k: errs[i] for k,i in ids.items() }
    return ren, renerr, fp, fe

ren1, ren1err, fitp1, fite1 = build(flags1, pn1, vals1, errs1)
if fit2file:
    ren2, ren2err, fitp2, fite2 = build(flags2, pn2, vals2, errs2)

# ─── defaults + Im‐builder with safety guards ────────────────────────────────────
defaults = {
    "H":  dict(r=0.9,  alpha0=0.43, alpha1=0.85, n=1.35,  b=0.4, M2=0.64, P=1.0, factor=2.0),
    "Ht": dict(r=7.0,  alpha0=0.43, alpha1=0.85, n=0.6,   b=2.0, M2=0.8,  P=1.0, factor=0.4),
    "E":  dict(r=0.9,  alpha0=0.43, alpha1=0.85, n=1.35,  b=0.4, M2=0.64, P=1.0, factor=1.0),
    "Et": dict(r=0.9,  alpha0=0.43, alpha1=0.85, n=0.0,   b=0.4, M2=0.64, P=1.0, factor=1.0),
}

def make_Im_func(cff, pars, ren):
    d = defaults[cff]
    def Im(xi, t):
        # --- 1) xi_safe to avoid xi=0 division ---
        xi_safe = np.maximum(xi, 1e-6)
        r     = pars.get("r",      d["r"])
        a0    = pars.get("alpha0", d["alpha0"])
        a1    = pars.get("alpha1", d["alpha1"])
        nval  = pars.get("n",      d["n"])
        bval  = pars.get("b",      d["b"])
        M2    = pars.get("M2",     d["M2"])
        Pval  = pars.get("P",      d["P"])
        fac   = d["factor"]

        alpha = a0 + a1 * t
        pref  = np.pi * 5.0/9.0 * nval * r / (1 + xi_safe)

        with np.errstate(divide='ignore', invalid='ignore'):
            xfac = (2*xi_safe/(1+xi_safe))**(-alpha)
        yfac  = ((1 - xi_safe)/(1+xi_safe))**(bval)

        # --- 2) clip the base of the power to ≥0 before raising to non-int exponent ---
        base = 1 - ((1 - xi_safe)/(1+xi_safe)) * t / M2
        base_clipped = np.clip(base, 0.0, None)
        with np.errstate(invalid='ignore'):
            tfac = base_clipped**(-Pval)

        result = ren * pref * xfac * yfac * tfac * fac
        # --- 3) scrub any remaining NaNs → 0 ---
        return np.nan_to_num(result, nan=0.0, posinf=0.0, neginf=0.0)
    return Im

def gen_reps(cent, err, n=N_REPLICAS):
    reps = []
    for _ in range(n):
        d = {}
        for k,v in cent.items():
            d[k] = np.random.normal(v, err[k]/1.96)
        reps.append(d)
    return reps

def band(cff, xi_vals, t_vals, ren, cent, err, n=N_REPLICAS):
    xi   = np.atleast_1d(xi_vals)
    t    = np.atleast_1d(t_vals)
    xi_b, t_b = np.broadcast_arrays(xi, -t)
    reps      = gen_reps(cent[cff], err[cff], n)
    ren_reps  = np.random.normal(ren, ren/np.sqrt(len(xi_b)), n)
    allA = []
    for i in range(n):
        f = make_Im_func(cff, reps[i], ren_reps[i])
        allA.append(f(xi_b, t_b))
    M = np.array(allA)
    return np.median(M, axis=0), np.percentile(M, 2.5, axis=0), np.percentile(M, 97.5, axis=0)

# ─── plotting setup ─────────────────────────────────────────────────────────────
os.makedirs("output/plots", exist_ok=True)
t_fixed   = np.linspace(0.1, 1.0, 6)
xi_fixed  = np.linspace(0.05,0.50, 6)
xi_range  = np.linspace(0.0, 0.5, 200)
t_range   = np.linspace(0.0, 1.0, 200)

plt.style.use('classic')
plt.rcParams.update({'font.size':14,'font.family':'serif'})

style1 = {'color':'tab:red',   'ls':'--','lw':2.5}
band1  = {'color':'tab:red',   'alpha':0.2}
style2 = {'color':'tab:blue',  'ls':'-','lw':2.5}
band2  = {'color':'tab:blue',  'alpha':0.2}
zero   = {'color':'gray',      'ls':'--','lw':1}

lbl1   = "RGK preliminary"
lbl2   = "RGA pass-1"
tex_map = {"H":"H","Ht":r"\tilde H","E":"E","Et":r"\tilde E"}

# ─── Single RGK only ────────────────────────────────────────────────────────────
for cff in fitp1:
    tex = tex_map[cff]

    # Im vs ξ
    fig, axes = plt.subplots(2,3, figsize=(12,8), sharex=True, sharey=True)
    axes = axes.flatten()
    fig.suptitle(rf"$\mathrm{{Im}}\,{tex}$", fontsize=16, y=0.98)
    for ax, t in zip(axes, t_fixed):
        m, lo, up = band(cff, xi_range, t, ren1, fitp1, fite1)
        ax.plot( xi_range, m, **style1)
        ax.fill_between(xi_range, lo, up, **band1)
        ax.axhline(0, **zero)
        ax.text(0.6, 0.75, rf"$-t={t:.2f}$", transform=ax.transAxes)
        ax.set_xlim(0,0.5); ax.set_ylim(0,20); ax.set_yticks([0,5,10,15,20])
    # drop 0 tick in top row
    for ax in axes[:3]:
        ax.set_yticks([5,10,15,20])
    # hide "0.0" x label bottom-middle/right
    for ax in (axes[4], axes[5]):
        for lbl in ax.get_xticklabels():
            if lbl.get_text() in ('0','0.0'):
                lbl.set_visible(False)
    axes[2].legend(handles=[Line2D([0],[0],**style1,label=lbl1)],
                   loc='upper right', fontsize=10)
    fig.text(0.06,0.5, rf"$\mathrm{{Im}}\,{tex}(\xi,\,-t)$",
             va='center', ha='center', rotation='vertical')
    fig.text(0.5,0.02, r"$\xi$", ha='center')
    fig.subplots_adjust(left=0.10, right=0.98, top=0.92, bottom=0.08, wspace=0, hspace=0)
    out = f"output/plots/RGK_Im{cff}_vs_xi_{ts1}.pdf"
    fig.savefig(out, bbox_inches='tight')
    print("Saved:", out)
    plt.close(fig)

    # Im vs −t
    fig, axes = plt.subplots(2,3, figsize=(12,8), sharex=True, sharey=True)
    axes = axes.flatten()
    fig.suptitle(rf"$\mathrm{{Im}}\,{tex}$", fontsize=16, y=0.98)
    for ax, xi0 in zip(axes, xi_fixed):
        m, lo, up = band(cff, xi0, t_range, ren1, fitp1, fite1)
        ax.plot( t_range, m, **style1)
        ax.fill_between(t_range, lo, up, **band1)
        ax.axhline(0, **zero)
        ax.text(0.6, 0.75, rf"$\xi={xi0:.2f}$", transform=ax.transAxes)
        ax.set_xlim(0,1.0); ax.set_ylim(0,20); ax.set_yticks([0,5,10,15,20])
    for ax in axes[:3]:
        ax.set_yticks([5,10,15,20])
    for ax in (axes[4], axes[5]):
        for lbl in ax.get_xticklabels():
            if lbl.get_text() in ('0','0.0'): lbl.set_visible(False)
    axes[2].legend(handles=[Line2D([0],[0],**style1,label=lbl1)],
                   loc='upper right', fontsize=10)
    fig.text(0.06,0.5, rf"$\mathrm{{Im}}\,{tex}(\xi,\,-t)$",
             va='center', ha='center', rotation='vertical')
    fig.text(0.5,0.02, r"$-t\,(\mathrm{GeV^2})$", ha='center')
    fig.subplots_adjust(left=0.10, right=0.98, top=0.92, bottom=0.08, wspace=0, hspace=0)
    out = f"output/plots/RGK_Im{cff}_vs_t_{ts1}.pdf"
    fig.savefig(out, bbox_inches='tight')
    print("Saved:", out)
    plt.close(fig)

# ─── RGA vs RGK comparison ──────────────────────────────────────────────────────
if fit2file:
    for cff in fitp1:
        tex = tex_map[cff]

        # Im vs ξ
        fig, axes = plt.subplots(2,3, figsize=(12,8), sharex=True, sharey=True)
        axes = axes.flatten()
        fig.suptitle(rf"$\mathrm{{Im}}\,{tex}$", fontsize=16, y=0.98)
        for ax, t in zip(axes, t_fixed):
            m1, l1, u1 = band(cff, xi_range, t, ren1, fitp1, fite1)
            m2, l2, u2 = band(cff, xi_range, t, ren2, fitp2, fite2)
            ax.plot( xi_range, m2, **style2)
            ax.fill_between(xi_range, l2, u2, **band2)
            ax.plot( xi_range, m1, **style1)
            ax.fill_between(xi_range, l1, u1, **band1)
            ax.axhline(0, **zero)
            ax.text(0.6, 0.75, rf"$-t={t:.2f}$", transform=ax.transAxes)
            ax.set_xlim(0,0.5); ax.set_ylim(0,20); ax.set_yticks([0,5,10,15,20])
        for ax in axes[:3]:
            ax.set_yticks([5,10,15,20])
        for ax in (axes[4], axes[5]):
            for lbl in ax.get_xticklabels():
                if lbl.get_text() in ('0','0.0'): lbl.set_visible(False)
        axes[2].legend(handles=[
            Line2D([0],[0],**style2,label=lbl2),
            Line2D([0],[0],**style1,label=lbl1),
        ], loc='upper right', fontsize=10)
        fig.text(0.06,0.5, rf"$\mathrm{{Im}}\,{tex}(\xi,\,-t)$",
                 va='center', ha='center', rotation='vertical')
        fig.text(0.5,0.02, r"$\xi$", ha='center')
        fig.subplots_adjust(left=0.10, right=0.98, top=0.92, bottom=0.08, wspace=0, hspace=0)
        out = f"output/plots/RGA_vs_RGK_Im{cff}_vs_xi_{ts2}_{ts1}.pdf"
        fig.savefig(out, bbox_inches='tight')
        print("Saved:", out)
        plt.close(fig)

        # Im vs −t
        fig, axes = plt.subplots(2,3, figsize=(12,8), sharex=True, sharey=True)
        axes = axes.flatten()
        fig.suptitle(rf"$\mathrm{{Im}}\,{tex}$", fontsize=16, y=0.98)
        for ax, xi0 in zip(axes, xi_fixed):
            m1, l1, u1 = band(cff, xi0, t_range, ren1, fitp1, fite1)
            m2, l2, u2 = band(cff, xi0, t_range, ren2, fitp2, fite2)
            ax.plot( t_range, m2, **style2)
            ax.fill_between(t_range, l2, u2, **band2)
            ax.plot( t_range, m1, **style1)
            ax.fill_between(t_range, l1, u1, **band1)
            ax.axhline(0, **zero)
            ax.text(0.6, 0.75, rf"$\xi={xi0:.2f}$", transform=ax.transAxes)
            ax.set_xlim(0,1.0); ax.set_ylim(0,20); ax.set_yticks([0,5,10,15,20])
        for ax in axes[:3]:
            ax.set_yticks([5,10,15,20])
        for ax in (axes[4], axes[5]):
            for lbl in ax.get_xticklabels():
                if lbl.get_text() in ('0','0.0'): lbl.set_visible(False)
        axes[2].legend(handles=[
            Line2D([0],[0],**style2,label=lbl2),
            Line2D([0],[0],**style1,label=lbl1),
        ], loc='upper right', fontsize=10)
        fig.text(0.06,0.5, rf"$\mathrm{{Im}}\,{tex}(\xi,\,-t)$",
                 va='center', ha='center', rotation='vertical')
        fig.text(0.5,0.02, r"$-t\,(\mathrm{GeV^2})$", ha='center')
        fig.subplots_adjust(left=0.10, right=0.98, top=0.92, bottom=0.08, wspace=0, hspace=0)
        out = f"output/plots/RGA_vs_RGK_Im{cff}_vs_t_{ts2}_{ts1}.pdf"
        fig.savefig(out, bbox_inches='tight')
        print("Saved:", out)
        plt.close(fig)

print("All done.")