#!/usr/bin/env python3
"""
plot_ImCFFs_fit_results.py  (BSA per-bin: BKM full vs BH-form approximations)

Curves:
  • BKM ............. solid blue        (full dv.BSA())
  • BKM approx ...... dashed blue       (A sinφ / (1 + B cosφ), A from BKM, B from data or model)
  • RGK fit ......... solid orange      (same BH-form using fit params)

Usage:
  python plot_ImCFFs_fit_results.py output/fit_results/fit_results_<TIMESTAMP>.txt
         [--data PATH] [--B-source data|model] [--outdir output/plots] [--debug]
"""
import os, re, argparse
import numpy as np
import matplotlib.pyplot as plt
import ROOT  # PyROOT

# ──────────────────────────────────────────────────────────────────────────────
def _load_dvcs_once():
    if getattr(_load_dvcs_once, "_done", False):
        return
    ROOT.gInterpreter.Declare(r"""
#ifndef __DVCS_XSEC_BSA_PLOT_GUARD__
#define __DVCS_XSEC_BSA_PLOT_GUARD__
#include "DVCS_xsec.C"
#endif
""")
    _load_dvcs_once._done = True

def parse_fit_results(fname):
    with open(fname) as f:
        lines = [l.strip() for l in f if l.strip()]
    fl = next((l for l in lines
               if re.match(r'^\s*H\s+[01]\s+Ht\s+[01]\s+E\s+[01]\s+Et\s+[01]\s*$', l)), None)
    if fl is None: fl = next((l for l in lines if l.startswith("H ")), None)
    if fl is None: raise RuntimeError("No 'H ... Ht ... E ... Et ...' line in fit file")
    toks = fl.split()
    flags = {toks[i]: int(toks[i+1]) for i in range(0, len(toks), 2)}
    pnames, pvals = [], None
    for i, l in enumerate(lines):
        if l.startswith("# parameters"): pnames = l.split()[2:]
        elif l.startswith("# values") and i+1 < len(lines):
            pvals = list(map(float, lines[i+1].split()))
    if not pnames or pvals is None:
        raise RuntimeError("Couldn't parse parameter names/values from fit file")
    params = dict(zip(pnames, pvals))
    input_path = next((l.split(None,1)[1].strip() for l in lines if l.startswith("input")), None)
    m = re.search(r'fit_results_(\d{8}_\d{6})\.txt$', os.path.basename(fname))
    ts = (m.group(1) if m else "noTS")
    return flags, params, input_path, ts

def load_bins(datafile):
    bins, curr, prev_phi = [], {k: [] for k in ("phi","Q2","xB","t","Eb","A","sigA")}, None
    with open(datafile) as f:
        for line in f:
            if not line.strip() or line.lstrip().startswith("#"): continue
            phi, Q2, xB, t, Eb, A, sigA = map(float, line.split())
            if prev_phi is not None and phi < prev_phi:
                arr = {k: np.array(v, float) for k, v in curr.items()}
                arr["Q2m"], arr["xBm"], arr["tm"], arr["Ebm"] = (
                    float(arr["Q2"].mean()), float(arr["xB"].mean()),
                    float(arr["t"].mean()),  float(arr["Eb"].mean())
                )
                bins.append(arr); curr = {k: [] for k in curr}
            for k, v in zip(curr.keys(), (phi, Q2, xB, t, Eb, A, sigA)):
                curr[k].append(v)
            prev_phi = phi
    if curr["phi"]:
        arr = {k: np.array(v, float) for k, v in curr.items()}
        arr["Q2m"], arr["xBm"], arr["tm"], arr["Ebm"] = (
            float(arr["Q2"].mean()), float(arr["xB"].mean()),
            float(arr["t"].mean()),  float(arr["Eb"].mean())
        )
        bins.append(arr)
    return bins

def bkm_defaults():
    return {
        'renormImag': 1.0, 'renormReal': 1.0,
        'r_H': 0.9, 'n_H': 1.35, 'alpha0_H': 0.43, 'alpha1_H': 0.85, 'b_H': 0.4, 'M2_H': 0.64, 'P_H': 1.0,
        'r_Ht': 7.0, 'n_Ht': 0.6,  'alpha0_Ht': 0.43, 'alpha1_Ht': 0.85, 'b_Ht': 2.0, 'M2_Ht': 0.8,  'P_Ht': 1.0,
        'r_E': 0.9, 'n_E': 1.35, 'alpha0_E': 0.43, 'alpha1_E': 0.85, 'b_E': 0.4, 'M2_E': 0.64, 'P_E': 1.0,
        'r_Et': 1.0, 'n_Et': 0.0,  'alpha0_Et': 0.0,  'alpha1_Et': 0.0, 'b_Et': 0.0, 'M2_Et': 1.0, 'P_Et': 0.0,
    }

def push_globals(param_map, flags, *, label=""):
    # reset switches
    ROOT.gInterpreter.ProcessLine("hasH=0; hasHt=0; hasE=0; hasEt=0;")
    # Et OFF by design here
    flags_eff = dict(flags); flags_eff["Et"] = 0
    for cff in ("H","Ht","E","Et"):
        ROOT.gInterpreter.ProcessLine(f"has{cff} = {int(flags_eff.get(cff,0))};")
    # renorms
    ROOT.gInterpreter.ProcessLine(f"renormImag = {float(param_map.get('renormImag',1.0))};")
    ROOT.gInterpreter.ProcessLine(f"renormReal = {float(param_map.get('renormReal',1.0))};")
    # parameters
    def setp(name):
        if name in param_map:
            ROOT.gInterpreter.ProcessLine(f"{name} = {float(param_map[name])};")
    for cff in ("H","Ht","E","Et"):
        for k in ("r","n","alpha0","alpha1","b","M2","P"):
            setp(f"{k}_{cff}")
    for cff in ("H","Ht","E","Et"):
        for k in ("C0","MD2","lambda"):
            setp(f"{k}_{cff}")

# full dv.BSA curve
def full_bsa_curve(phi_deg, Q2, xB, t, Eb):
    out = np.empty_like(phi_deg, float)
    for i, ph in enumerate(phi_deg):
        dv = ROOT.BMK_DVCS(-1, +1, 0, Eb, xB, Q2, t, ph)
        out[i] = dv.BSA()
    return out

# A and B from model (BMK harmonics) – with BMK→Trento sign fix on B
def model_AB(Q2, xB, t, Eb):
    dv = ROOT.BMK_DVCS(-1, +1, 0, Eb, xB, Q2, t, 0.0)
    c0 = dv.c0_BH()
    if c0 == 0: return 0.0, 0.0
    A = dv.s1_I()/c0
    B = - dv.c1_BH()/c0   # sign flip: cos(phi_BMK) = -cos(phi_trento)
    return A, B

# Weighted fit of C + A sinφ/(1 + B cosφ) to get B from data (and C,A if needed)
def fit_B_from_data(phi_deg, y, dy):
    ph = np.deg2rad(phi_deg)
    s, c = np.sin(ph), np.cos(ph)
    w = 1.0/np.clip(dy, 1e-6, None)**2
    # grid scan over B to avoid SciPy; bounds keep denom away from poles
    Bgrid = np.linspace(-0.95, 0.95, 761)
    best = (np.inf, 0.0, 0.0, 0.0)  # chi2, B, C, A
    for B in Bgrid:
        denom = 1.0 + B*c
        X0 = np.ones_like(s)
        X1 = s/denom
        # weighted linear least squares for [C, A]
        W = np.diag(w)
        X = np.vstack([X0, X1]).T
        XT_W = X.T @ W
        try:
            theta = np.linalg.solve(XT_W @ X, XT_W @ y)
        except np.linalg.LinAlgError:
            continue
        Cfit, Afit = theta
        yhat = Cfit + Afit*X1
        chi2 = np.sum(((y - yhat)*np.sqrt(w))**2)
        if chi2 < best[0]:
            best = (chi2, B, Cfit, Afit)
    _, Bfit, Cfit, Afit = best
    return Bfit, Cfit, Afit

def bh_form_curve(phi_deg, A, B):
    ph = np.deg2rad(phi_deg)
    return A * np.sin(ph) / (1.0 + B * np.cos(ph))

# ──────────────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('fitfile')
    ap.add_argument('--data', default=None, help='Override data file; else use fit-file "input"')
    ap.add_argument('--B-source', choices=['data','model'], default='data',
                    help='Use per-bin B from data-fit (default) or from BH harmonics')
    ap.add_argument('--outdir', default='output/plots', help='Directory for PDFs')
    ap.add_argument('--debug', action='store_true')
    args = ap.parse_args()

    flags_fit, params_fit, input_path, ts = parse_fit_results(args.fitfile)
    datafile = args.data or input_path or 'imports/rgk_preliminary_bsa.txt'
    print(">> Using data file:", datafile)
    print(">> Fit flags:", flags_fit)
    print(">> Fit params:", params_fit)

    _load_dvcs_once()
    bins = load_bins(datafile)
    if not bins: raise RuntimeError(f"No φ-bins parsed from {datafile}")
    print(f">> Found {len(bins)} φ-bins")

    params_bkm = bkm_defaults()
    flags_bkm  = {'H':1,'Ht':1,'E':1,'Et':0}
    os.makedirs(args.outdir, exist_ok=True)

    for ibin, b in enumerate(bins, start=1):
        phi_data, A_data, dA = b["phi"], b["A"], b["sigA"]
        Q2m, xBm, tm, Ebm = b["Q2m"], b["xBm"], b["tm"], b["Ebm"]
        phi_grid = np.linspace(0.0, 360.0, 361)

        # fit B from data once per bin if requested
        B_from_data = None
        if args.B_source == 'data':
            B_from_data, Ctmp, Atmp = fit_B_from_data(phi_data, A_data, dA)

        # --- BKM full + approx
        push_globals(params_bkm, flags_bkm, label="BKM")
        bkm_full = full_bsa_curve(phi_grid, Q2m, xBm, tm, Ebm)
        A_bkm, B_bkm_model = model_AB(Q2m, xBm, tm, Ebm)
        B_bkm = B_from_data if B_from_data is not None else B_bkm_model
        bkm_approx = bh_form_curve(phi_grid, A_bkm, B_bkm)

        # --- RGK fit approx
        push_globals(params_fit, flags_fit, label="RGK fit")
        A_rgk, B_rgk_model = model_AB(Q2m, xBm, tm, Ebm)
        B_rgk = B_from_data if B_from_data is not None else B_rgk_model
        rgk_curve = bh_form_curve(phi_grid, A_rgk, B_rgk)

        if args.debug and ibin == 1:
            print(f"[BKM] A_model={A_bkm:+.6e}, B_model={B_bkm_model:+.6e}, "
                  f"B_data={B_from_data if B_from_data is not None else float('nan'):+.6e}")
            print(f"[RGK] A_model={A_rgk:+.6e}, B_model={B_rgk_model:+.6e}")

        # --- Plot
        fig, ax = plt.subplots(figsize=(8,5))
        ax.errorbar(phi_data, A_data, yerr=dA, fmt='o', ms=5, color='k', label='Data')
        # BKM full = solid blue
        ax.plot(phi_grid, bkm_full, '-',  lw=2.5, color='tab:blue',  label='BKM')
        # BKM approx = dashed blue
        ax.plot(phi_grid, bkm_approx, '--', lw=2.0, color='tab:blue',  label='BKM approx')
        # RGK fit = solid orange
        ax.plot(phi_grid, rgk_curve, '-',   lw=2.5, color='tab:orange', label='RGK fit')

        ax.set_xlim(0, 360)
        ax.set_xticks([0,60,120,180,240,300,360])
        ax.set_ylim(-0.6, 0.6)
        ax.set_xlabel(r'$\phi\;[\mathrm{deg}]$')
        ax.set_ylabel(r'$A_{LU}(\phi)$')
        ax.set_title(
            (r'$\langle Q^2\rangle={:.2f}\,\mathrm{{GeV}}^2,\;\langle x_B\rangle={:.3f},\;'
             r'\langle -t\rangle={:.3f}\,\mathrm{{GeV}}^2,\;\langle E_b\rangle={:.2f}\,\mathrm{{GeV}}$'
            ).format(Q2m, xBm, abs(tm), Ebm),
            pad=12
        )
        ax.legend(loc='upper right', frameon=True, edgecolor='k')
        plt.tight_layout()

        outname = os.path.join(
            args.outdir,
            f"BSA_bin{ibin:02d}_{ts}_Q2_{Q2m:.2f}_xB_{xBm:.3f}_mt_{abs(tm):.3f}.pdf"
        )
        fig.savefig(outname)
        print(f">> Saved bin {ibin} plot to {outname}")
        plt.close(fig)

if __name__ == "__main__":
    main()