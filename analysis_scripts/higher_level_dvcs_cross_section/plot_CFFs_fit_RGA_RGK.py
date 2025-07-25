#!/usr/bin/env python3
"""
plot_ImCFFs_fit_results.py

Usage:
    python plot_ImCFFs_fit_results.py fit_results_<TS1>.txt fit_results_<TS2>.txt

Reads two fit‐results files and for each enabled Im‐CFF:
  • Creates per‐file plots (only fitted result, no original model) with median and 95% CI,
    saved as Im{CFF}_vs_xi_<TS1>.pdf and Im{CFF}_vs_t_<TS1>.pdf in output/plots/.
  • Creates comparison plots overlaying both fits, saved as
    Im{CFF}_vs_xi_<TS1>_2.pdf and Im{CFF}_vs_t_<TS1>_2.pdf.
Y‐axes run from –2 to 20.
"""
import os, sys, re
import numpy as np
import matplotlib.pyplot as plt

# ――― Parse args ―――
if len(sys.argv) != 3:
    print("Usage: python plot_ImCFFs_fit_results.py <fit1.txt> <fit2.txt>")
    sys.exit(1)
file1, file2 = sys.argv[1], sys.argv[2]
def extract_ts(fname):
    m = re.search(r'fit_results_(\d{8}_\d{6})\.txt$', fname)
    if not m:
        raise RuntimeError(f"Can't parse timestamp from {fname}")
    return m.group(1)
ts1, ts2 = extract_ts(file1), extract_ts(file2)

# ――― Loader & parser ―――
def parse_fit(fname):
    with open(fname) as f:
        lines = [l.strip() for l in f if l.strip()]
    flag_line = next(l for l in lines if l.startswith("H "))
    toks = flag_line.split()
    flags = {toks[i]: int(toks[i+1]) for i in range(0, len(toks), 2)}
    # params
    pnames = []
    for i,l in enumerate(lines):
        if l.startswith("# parameters"):
            pnames = lines[i].split()[2:]; break
    vals = errs = None
    for i,l in enumerate(lines):
        if l.startswith("# values"):
            vals = list(map(float, lines[i+1].split()))
        if l.startswith("# errors"):
            errs = list(map(float, lines[i+1].split()))
    if vals is None or errs is None:
        raise RuntimeError("Bad fit file")
    return flags, pnames, np.array(vals), np.array(errs)

flags1, p1, v1, e1 = parse_fit(file1)
flags2, p2, v2, e2 = parse_fit(file2)

# ――― Extract blocks ―――
def build_dicts(pnames, vals, errs):
    renorm = vals[pnames.index("renormImag")]
    renorm_err = errs[pnames.index("renormImag")]
    fp, fe = {}, {}
    for cff in ("H","Ht","E","Et"):
        if flags1[cff]:
            keys = ["r","alpha0","alpha1","n","b","M2","P"]
            fp[cff] = {k: vals[pnames.index(f"{k}_{cff}")] for k in keys}
            fe[cff] = {k: errs[pnames.index(f"{k}_{cff}")] for k in keys}
    return renorm, renorm_err, fp, fe

r1, re1, fp1, fe1 = build_dicts(p1, v1, e1)
r2, re2, fp2, fe2 = build_dicts(p2, v2, e2)

# ――― Defaults ―――
defaults = {
    "H":  dict(r=0.9, alpha0=0.43, alpha1=0.85, n=1.35, b=0.4,  M2=0.64, P=1.0, factor=2.0),
    "Ht": dict(r=7.0, alpha0=0.43, alpha1=0.85, n=0.6,  b=2.0,  M2=0.8,  P=1.0, factor=0.4),
    "E":  dict(r=0.9, alpha0=0.43, alpha1=0.85, n=1.35, b=0.4,  M2=0.64, P=1.0, factor=1.0),
    "Et": dict(r=0.9, alpha0=0.43, alpha1=0.85, n=0.0,  b=0.4,  M2=0.64, P=1.0, factor=1.0),
}

def make_Im(cff, params, renorm):
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
        alpha = a0 + a1*t
        pref  = np.pi*5/9 * nval*r/(1+xi)
        xfac  = (2*xi/(1+xi))**(-alpha)
        yfac  = ((1-xi)/(1+xi))**(bval)
        tfac  = (1 - ((1-xi)/(1+xi))*t/M2)**(-Pval)
        return renorm * pref*xfac*yfac*tfac*fac
    return Im

# ――― Band calculator ―――
def compute_band(params, errors, renorm, xi, t, nrep=500):
    xi_a = np.atleast_1d(xi); t_a = np.atleast_1d(t)
    xi_b, t_b = np.broadcast_arrays(xi_a, -t_a)
    allc = []
    for _ in range(nrep):
        pr = {k: np.random.normal(params[k], errors[k]/1.96) for k in params}
        rr = np.random.normal(renorm, errors.get("renormImag", errors.get("renormImag",renorm))/1.96)
        f = make_Im(cff, pr, rr)
        allc.append(f(xi_b, t_b))
    allc = np.array(allc)
    return np.median(allc,0), np.percentile(allc,2.5,0), np.percentile(allc,97.5,0)

# ――― Plotting setup ―――
outdir = "output/plots"
os.makedirs(outdir, exist_ok=True)
t_fixed   = np.linspace(0.1,1.0,6)
xi_range  = np.linspace(0.0,0.5,200)
xi_fixed  = np.linspace(0.05,0.5,6)
t_range   = np.linspace(0.0,1.0,200)
import matplotlib.lines as mlines
for cff in ("H","Ht","E","Et"):
    if not flags1[cff]: continue
    Im1 = make_Im(cff, fp1[cff], r1)
    Im2 = make_Im(cff, fp2[cff], r2)

    # — single-file plots —
    for what, xs, ys, ts, tsign in [
        ("xi", xi_range, None, t_fixed, True),
        ("t",  t_range,  xi_fixed, None, False),
    ]:
        fig, axes = plt.subplots(2,3,figsize=(12,8), sharex=True, sharey=True)
        axes = axes.flatten()
        for ax,i in enumerate(axes):
            if what=="xi":
                tval = ts[i]
                med1,lo1,up1 = compute_band(fp1[cff], fe1[cff], r1, xi_range, tval)
                ax.plot(xi_range, med1, color='tab:red', linestyle='--', lw=2.5, label="RGK preliminary")
                ax.fill_between(xi_range, lo1, up1, color='tab:red', alpha=0.2)
                ax.set_xlabel(r"$\xi$"); ax.set_xlim(0,0.5)
                ax.set_title(rf"$-t={tval:.2f}$")
            else:
                xi0 = ys[i]
                med1,lo1,up1 = compute_band(fp1[cff], fe1[cff], r1, xi0, t_range)
                ax.plot(t_range, med1, color='tab:red', linestyle='--', lw=2.5, label="RGK preliminary")
                ax.fill_between(t_range, lo1, up1, color='tab:red', alpha=0.2)
                ax.set_xlabel(r"$-t\ (\mathrm{GeV}^2)$"); ax.set_xlim(0,1)
                ax.set_title(rf"$\xi={xi0:.2f}$")
            ax.axhline(0, color='gray', ls='--')
            ax.set_ylim(-2,20)
        # hide "-2" tick on top row
        for ax in axes[:3]:
            yts = ax.get_yticks()
            ax.set_yticks(yts[1:])  # drop the first
        axes[2].legend(loc='upper right'); fig.tight_layout()
        fn = os.path.join(outdir, f"Im{cff}_vs_{what}_{ts1}.pdf")
        fig.savefig(fn, bbox_inches='tight'); plt.close(fig)

    # — comparison plots —
    for what, xs, ys, ts, tsign in [
        ("xi", xi_range, None, t_fixed, True),
        ("t",  t_range,  xi_fixed, None, False),
    ]:
        fig, axes = plt.subplots(2,3,figsize=(12,8), sharex=True, sharey=True)
        axes = axes.flatten()
        for ax,i in enumerate(axes):
            if what=="xi":
                tval = ts[i]
                m1,l1,u1 = compute_band(fp1[cff], fe1[cff], r1, xi_range, tval)
                m2,l2,u2 = compute_band(fp2[cff], fe2[cff], r2, xi_range, tval)
                ax.plot(xi_range, m1, color='tab:red', linestyle='--', lw=2.5, label="RGK preliminary")
                ax.fill_between(xi_range, l1, u1, color='tab:red', alpha=0.2)
                ax.plot(xi_range, m2, color='tab:green', linestyle='--', lw=2.5, label="Comparison")
                ax.fill_between(xi_range, l2, u2, color='tab:green', alpha=0.2)
                ax.set_xlabel(r"$\xi$"); ax.set_xlim(0,0.5)
                ax.set_title(rf"$-t={tval:.2f}$")
            else:
                xi0 = ys[i]
                m1,l1,u1 = compute_band(fp1[cff], fe1[cff], r1, xi0, t_range)
                m2,l2,u2 = compute_band(fp2[cff], fe2[cff], r2, xi0, t_range)
                ax.plot(t_range, m1, color='tab:red', linestyle='--', lw=2.5, label="RGK preliminary")
                ax.fill_between(t_range, l1, u1, color='tab:red', alpha=0.2)
                ax.plot(t_range, m2, color='tab:green', linestyle='--', lw=2.5, label="Comparison")
                ax.fill_between(t_range, l2, u2, color='tab:green', alpha=0.2)
                ax.set_xlabel(r"$-t\ (\mathrm{GeV}^2)$"); ax.set_xlim(0,1)
                ax.set_title(rf"$\xi={xi0:.2f}$")
            ax.axhline(0, color='gray', ls='--')
            ax.set_ylim(-2,20)
        for ax in axes[:3]:
            yts = ax.get_yticks()
            ax.set_yticks(yts[1:])
        axes[2].legend(loc='upper right'); fig.tight_layout()
        fn = os.path.join(outdir, f"Im{cff}_vs_{what}_{ts1}_2.pdf")
        fig.savefig(fn, bbox_inches='tight'); plt.close(fig)