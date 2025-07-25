#!/usr/bin/env python3
"""
plot_CFFs_fit_results.py

Usage:
    python plot_CFFs_fit_results.py output/fit_results/fit_results_<TIMESTAMP>.txt

Reads which CFFs were fit from the header of results file, then for each
enabled Im CFF makes two figures:
  1) Im CFF vs. xi for six fixed −t between 0.1 and 1.0 (GeV^{2})
  2) Im CFF vs. −t for six fixed xi between 0.05 and 0.50

Includes uncertainty bands (95% CI) via replica method.
"""
import os
import sys
import re
import numpy as np
import matplotlib.pyplot as plt

# ─── Parse command-line ─────────────────────────────────────────────────────────
if len(sys.argv) != 2:
    print("Usage: python plot_CFFs_fit_results.py "
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

    pnames = []
    for i,l in enumerate(lines):
        if l.startswith("# parameters"):
            pnames = lines[i].split()[2:]
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

# ─── Extract fit-parameters ─────────────────────────────────────────────────────
def get_idx(name):
    return pnames.index(name) if name in pnames else None

renorm_fit = vals[get_idx("renormImag")]
renorm_err = errs[get_idx("renormImag")]

fit_params  = {}
fit_errors  = {}
for cff in ("H","Ht","E","Et"):
    if flags[cff]:
        ks = ["r","alpha0","alpha1","n","b","Mm2","P"]
        idxs = { k: get_idx(f"{k}_{cff}") for k in ks }
        fit_params[cff] = { k: vals[i] for k,i in idxs.items() }
        fit_errors[cff] = { k: errs[i] for k,i in idxs.items() }

# ─── VGG-default params + factors ───────────────────────────────────────────────
defaults = {
    "H":  dict(r=0.9,   alpha0=0.43, alpha1=0.85, n=1.35, b=0.4,  Mm2=0.64, P=1.0, factor=2.0),
    "Ht": dict(r=7.0,   alpha0=0.43, alpha1=0.85, n=0.6,  b=2.0,  Mm2=0.8,  P=1.0, factor=0.4),
    "E":  dict(r=0.9,   alpha0=0.43, alpha1=0.85, n=1.35, b=0.4,  Mm2=0.64, P=1.0, factor=1.0),
    "Et": dict(r=0.9,   alpha0=0.43, alpha1=0.85, n=0.0,  b=0.4,  Mm2=0.64, P=1.0, factor=1.0),
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
        Mm2  = params.get("Mm2",    d["Mm2"])
        Pval = params.get("P",      d["P"])
        fac  = d["factor"]

        alpha = a0 + a1 * t
        pref  = np.pi * 5.0/9.0 * nval * r / (1 + xi)
        xfac  = (2*xi/(1+xi))**(-alpha)
        yfac  = ((1 - xi)/(1+xi))**(bval)
        tfac  = (1 - ((1 - xi)/(1+xi))*t/Mm2)**(-Pval)
        return renorm * pref * xfac * yfac * tfac * fac
    return Im

# ─── Replica helpers ────────────────────────────────────────────────────────────
def generate_replica_params(central, errors, n=100):
    reps = []
    for _ in range(n):
        d = {}
        for k,v in central.items():
            sigma = errors[k] / 1.96
            d[k]   = np.random.normal(v, sigma)
        reps.append(d)
    return reps

# ─── Diagnostic compute_uncertainty_band ───────────────────────────────────────
def compute_uncertainty_band(cff, xi_vals, t_vals, n_reps=100):
    print(f"\nDBG: entering compute_uncertainty_band: cff={cff!r}, "
          f" xi_vals.shape={np.shape(xi_vals)}, t_vals.shape={np.shape(t_vals)}")

    if cff not in fit_params:
        return None, None, None

    param_reps  = generate_replica_params(fit_params[cff], fit_errors[cff], n_reps)
    renorm_reps = np.random.normal(renorm_fit, renorm_err/1.96, n_reps)

    xi_a = np.atleast_1d(xi_vals)
    t_a  = np.atleast_1d(t_vals)
    xi_b, t_b = np.broadcast_arrays(xi_a, -t_a)
    print(f"DBG: after broadcast xi_b.shape={xi_b.shape}, t_b.shape={t_b.shape}")

    all_curves = []
    for i in range(n_reps):
        Im_rep = make_Im_func(cff, param_reps[i], renorm_reps[i])
        try:
            curve = Im_rep(xi_b, t_b)
        except ValueError as err:
            print("!! ERROR inside Im_rep: xi_b.shape=", xi_b.shape,
                  " t_b.shape=", t_b.shape)
            raise
        all_curves.append(curve)

    all_curves = np.array(all_curves)
    return (
        np.median(all_curves, axis=0),
        np.percentile(all_curves,  2.5, axis=0),
        np.percentile(all_curves, 97.5, axis=0),
    )

# ─── Plot setup ─────────────────────────────────────────────────────────────────
Im_funcs = {
    cff: make_Im_func(cff, fit_params.get(cff, {}), renorm_fit)
    for cff in ("H","Ht","E","Et") if flags[cff]
}

plt.style.use('classic')
plt.rcParams.update({'font.size':14, 'font.family':'serif'})

outdir = 'output/plots'
os.makedirs(outdir, exist_ok=True)

t_fixed  = np.linspace(0.1, 1.0, 6)
xi_range = np.linspace(0.0, 0.5, 200)
xi_fixed = np.linspace(0.05,0.50, 6)
t_range  = np.linspace(0.0, 1.0, 200)

from matplotlib.lines import Line2D
legend_elems = [
    Line2D([0],[0], color='tab:blue', linestyle='-', lw=2.5, label='Original'),
    Line2D([0],[0], color='tab:red',  linestyle='--',lw=2.5, label='Fit'),
    plt.Rectangle((0,0),1,1, fc='tab:red',alpha=0.2, ec=None, label='95% CI'),
]

tex_map = {"H":"H", "Ht":r"\tilde H", "E":"E", "Et":r"\tilde E"}

for cff, Im_fit in Im_funcs.items():
    Im_orig = make_Im_func(cff, {}, 1.0)
    tex     = tex_map[cff]

    # — Im vs ξ —
    fig, axes = plt.subplots(2,3, figsize=(12,8), sharex=True, sharey=True)
    axes = axes.flatten()
    fig.suptitle(rf"$\mathrm{{Im}}\;{tex}$", fontsize=16, y=0.98)
    for ax, t in zip(axes, t_fixed):
        ax.plot(xi_range, Im_orig(xi_range, -t), color='tab:blue', lw=2.5)
        med, lo, up = compute_uncertainty_band(cff, xi_range, t, 500)
        if med is not None:
            ax.plot(xi_range, med, **{'color':'tab:red','ls':'--','lw':2.5})
            ax.fill_between(xi_range, lo, up, **{'color':'tab:red','alpha':0.2})
        ax.axhline(0, color='gray', ls='--')
        ax.text(0.6, 0.75, rf"$-t={t:.2f}$", transform=ax.transAxes)
        ax.set_xlim(0,0.5); ax.set_ylim(0,12)
    axes[2].legend(handles=legend_elems, loc='upper right')
    fig.text(0.06,0.5, rf"$\mathrm{{Im}}\;{tex}(\xi,\,-t)$",
             va='center', ha='center', rotation='vertical')
    fig.text(0.5,0.02, r"$\xi$", ha='center')
    fig.tight_layout(rect=[0,0,1,0.95])
    fig.savefig(f"{outdir}/Im{cff}_vs_xi_{timestamp}.pdf")
    plt.close(fig)

    # — Im vs −t —
    fig, axes = plt.subplots(2,3, figsize=(12,8), sharex=True, sharey=True)
    axes = axes.flatten()
    fig.suptitle(rf"$\mathrm{{Im}}\;{tex}$", fontsize=16, y=0.98)
    for ax, xi0 in zip(axes, xi_fixed):
        ax.plot(t_range, Im_orig(xi0, -t_range), color='tab:blue', lw=2.5)
        med, lo, up = compute_uncertainty_band(cff, xi0, t_range, 500)
        if med is not None:
            ax.plot(t_range, med, **{'color':'tab:red','ls':'--','lw':2.5})
            ax.fill_between(t_range, lo, up, **{'color':'tab:red','alpha':0.2})
        ax.axhline(0, color='gray', ls='--')
        ax.text(0.6, 0.75, rf"$\xi={xi0:.2f}$", transform=ax.transAxes)
        ax.set_xlim(0,1); ax.set_ylim((0,12) if xi0>0.2 else (-2,12))
    axes[2].legend(handles=legend_elems, loc='upper right')
    fig.text(0.06,0.5, rf"$\mathrm{{Im}}\;{tex}(\xi,\,-t)$",
             va='center', ha='center', rotation='vertical')
    fig.text(0.5,0.02, r"$-t$ (GeV$^2$)", ha='center')
    fig.tight_layout(rect=[0,0,1,0.95])
    fig.savefig(f"{outdir}/Im{cff}_vs_t_{timestamp}.pdf")
    plt.close(fig)

print("Done plotting.")
#test