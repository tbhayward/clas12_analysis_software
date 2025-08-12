#!/usr/bin/env python3
"""
plot_ImCFFs_fit_results.py   (repurposed)

Per φ-bin BSA plots showing:
  • BKM  ............ thin dashed blue  (full BSA from C++)
  • BKM approx ...... solid blue        (A sinφ / (1 + B cosφ), A=s1_I/c0_BH, B=c1_BH/c0_BH)
  • RGK fit ......... solid orange      (same BH-form, using fit parameters/flags)

Usage:
    python plot_ImCFFs_fit_results.py output/fit_results/fit_results_<TIMESTAMP>.txt
           [--data PATH] [--outdir output/plots] [--debug]

Notes:
  - Requires DVCS_xsec.C in CWD, with your latest Et fix (or Et disabled).
  - The fit file’s "input" line is used if --data isn’t provided.
"""

import os
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt
import ROOT  # PyROOT

# -----------------------------------------------------------------------------
# Load the C++ once (fresh guard to force current DVCS_xsec.C to compile)
# -----------------------------------------------------------------------------
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


# -----------------------------------------------------------------------------
# Fit-file parsing
# -----------------------------------------------------------------------------
def parse_fit_results(fname):
    with open(fname) as f:
        lines = [l.strip() for l in f if l.strip()]

    # flags line like "H 1 Ht 0 E 0 Et 0"
    flags_line = next((l for l in lines
                       if re.match(r'^\s*H\s+[01]\s+Ht\s+[01]\s+E\s+[01]\s+Et\s+[01]\s*$', l)), None)
    if flags_line is None:
        flags_line = next((l for l in lines if l.startswith("H ")), None)
    if flags_line is None:
        raise RuntimeError("No 'H ... Ht ... E ... Et ...' line in fit file")

    toks = flags_line.split()
    flags = {toks[i]: int(toks[i+1]) for i in range(0, len(toks), 2)}

    # parameter names and values
    pnames, pvals = [], None
    for i, l in enumerate(lines):
        if l.startswith("# parameters"):
            pnames = l.split()[2:]
        elif l.startswith("# values") and i+1 < len(lines):
            pvals = list(map(float, lines[i+1].split()))
    if not pnames or pvals is None:
        raise RuntimeError("Couldn't parse parameter names/values from fit file")
    param_map_fit = dict(zip(pnames, pvals))

    # optional input path
    input_path = next((l.split(None, 1)[1].strip()
                       for l in lines if l.startswith("input")), None)

    # timestamp convenience
    m = re.search(r'fit_results_(\d{8}_\d{6})\.txt$', os.path.basename(fname))
    ts = (m.group(1) if m else "noTS")
    return flags, param_map_fit, input_path, ts


# -----------------------------------------------------------------------------
# Data loading: split into φ-bins by wrap-around
# -----------------------------------------------------------------------------
def load_bins(datafile):
    bins, curr, prev_phi = [], {k: [] for k in ("phi","Q2","xB","t","Eb","A","sigA")}, None
    with open(datafile) as f:
        for line in f:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            phi, Q2, xB, t, Eb, A, sigA = map(float, line.split())
            if prev_phi is not None and phi < prev_phi:
                arr = {k: np.array(v, float) for k, v in curr.items()}
                arr["Q2m"], arr["xBm"], arr["tm"], arr["Ebm"] = (
                    float(arr["Q2"].mean()), float(arr["xB"].mean()),
                    float(arr["t"].mean()),  float(arr["Eb"].mean())
                )
                bins.append(arr)
                curr = {k: [] for k in curr}
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


# -----------------------------------------------------------------------------
# BKM default parameters (safe wrt Et)
# -----------------------------------------------------------------------------
def bkm_defaults():
    return {
        'renormImag': 1.0, 'renormReal': 1.0,
        # H
        'r_H': 0.9, 'n_H': 1.35, 'alpha0_H': 0.43, 'alpha1_H': 0.85, 'b_H': 0.4, 'M2_H': 0.64, 'P_H': 1.0,
        # Htilde
        'r_Ht': 7.0, 'n_Ht': 0.6,  'alpha0_Ht': 0.43, 'alpha1_Ht': 0.85, 'b_Ht': 2.0, 'M2_Ht': 0.8,  'P_Ht': 1.0,
        # E
        'r_E': 0.9, 'n_E': 1.35, 'alpha0_E': 0.43, 'alpha1_E': 0.85, 'b_E': 0.4, 'M2_E': 0.64, 'P_E': 1.0,
        # Et: keep harmless/off (n=0,P=0), nonzero M2 to avoid t/0 if on
        'r_Et': 1.0, 'n_Et': 0.0,  'alpha0_Et': 0.0,  'alpha1_Et': 0.0, 'b_Et': 0.0, 'M2_Et': 1.0, 'P_Et': 0.0,
    }


# -----------------------------------------------------------------------------
# Push globals into C++ (defensively disable Et)
# -----------------------------------------------------------------------------
def push_globals(param_map, flags, *, label=""):
    # reset switches
    ROOT.gInterpreter.ProcessLine("hasH=0; hasHt=0; hasE=0; hasEt=0;")

    # Et OFF no matter what (safe & per discussion)
    flags_eff = dict(flags)
    flags_eff["Et"] = 0
    for cff in ("H","Ht","E","Et"):
        ROOT.gInterpreter.ProcessLine(f"has{cff} = {int(flags_eff.get(cff,0))};")

    # renormalizations
    ROOT.gInterpreter.ProcessLine(f"renormImag = {float(param_map.get('renormImag',1.0))};")
    ROOT.gInterpreter.ProcessLine(f"renormReal = {float(param_map.get('renormReal',1.0))};")

    # parameter helper
    def setp(name):
        if name in param_map:
            ROOT.gInterpreter.ProcessLine(f"{name} = {float(param_map[name])};")

    for cff in ("H","Ht","E","Et"):
        for k in ("r","n","alpha0","alpha1","b","M2","P"):
            setp(f"{k}_{cff}")
    # optional subtraction constants if present
    for cff in ("H","Ht","E","Et"):
        for k in ("C0","MD2","lambda"):
            setp(f"{k}_{cff}")


# -----------------------------------------------------------------------------
# Physics helpers
# -----------------------------------------------------------------------------
def full_bsa_curve(phi_deg, Q2, xB, t, Eb):
    out = np.empty_like(phi_deg, float)
    for i, ph in enumerate(phi_deg):
        dv = ROOT.BMK_DVCS(-1, +1, 0, Eb, xB, Q2, t, ph)
        out[i] = dv.BSA()
    return out

def bh_form_params(Q2, xB, t, Eb):
    # coefficients do not depend on φ; evaluate once
    dv = ROOT.BMK_DVCS(-1, +1, 0, Eb, xB, Q2, t, 0.0)
    c0 = dv.c0_BH()
    if c0 == 0:
        return 0.0, 0.0
    A = dv.s1_I()/c0
    B = dv.c1_BH()/c0
    return A, B

def bh_form_curve(phi_deg, A, B):
    ph = np.deg2rad(phi_deg)
    return A * np.sin(ph) / (1.0 + B * np.cos(ph))


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('fitfile', help='fit_results_<TIMESTAMP>.txt')
    ap.add_argument('--data', default=None, help='Override BSA data file; else uses fit file "input"')
    ap.add_argument('--outdir', default='output/plots', help='Directory for PDFs')
    ap.add_argument('--debug', action='store_true', help='Print A,B for the first bin')
    args = ap.parse_args()

    # Parse fit results
    flags_fit, params_fit, input_path, ts = parse_fit_results(args.fitfile)

    # Choose data file
    datafile = args.data or input_path or 'imports/rgk_preliminary_bsa.txt'
    print(">> Using data file:", datafile)
    print(">> Fit flags:", flags_fit)
    print(">> Fit params:", params_fit)

    # Load C++
    _load_dvcs_once()

    # Load bins
    bins = load_bins(datafile)
    if not bins:
        raise RuntimeError(f"No φ-bins parsed from {datafile}")
    print(f">> Found {len(bins)} φ-bins")

    # BKM baseline: H,Ht,E on; Et off
    flags_bkm = {'H':1,'Ht':1,'E':1,'Et':0}
    params_bkm = bkm_defaults()

    os.makedirs(args.outdir, exist_ok=True)

    for ibin, b in enumerate(bins, start=1):
        phi_data, A_data, dA = b["phi"], b["A"], b["sigA"]
        Q2m, xBm, tm, Ebm = b["Q2m"], b["xBm"], b["tm"], b["Ebm"]

        # dense φ grid
        phi_grid = np.linspace(0.0, 360.0, 361)

        # --- BKM (full) + BKM (approx) ---
        push_globals(params_bkm, flags_bkm, label="BKM")
        bkm_full = full_bsa_curve(phi_grid, Q2m, xBm, tm, Ebm)
        A_bkm, B_bkm = bh_form_params(Q2m, xBm, tm, Ebm)
        bkm_approx = bh_form_curve(phi_grid, A_bkm, B_bkm)

        # --- RGK fit (approx only, as requested) ---
        push_globals(params_fit, flags_fit, label="RGK fit")
        A_rgk, B_rgk = bh_form_params(Q2m, xBm, tm, Ebm)
        rgk_curve = bh_form_curve(phi_grid, A_rgk, B_rgk)

        if args.debug and ibin == 1:
            print(f"[BKM]  A=s1_I/c0_BH={A_bkm:+.6e},  B=c1_BH/c0_BH={B_bkm:+.6e}")
            print(f"[RGK]  A=s1_I/c0_BH={A_rgk:+.6e},  B=c1_BH/c0_BH={B_rgk:+.6e}")

        # --- Plot ---
        fig, ax = plt.subplots(figsize=(8,5))
        # data
        ax.errorbar(phi_data, A_data, yerr=dA, fmt='o', ms=5, color='k', label='Data')

        # BKM approx (solid blue), BKM full (thin dashed blue)
        ax.plot(phi_grid, bkm_approx, '-',  lw=2.5, color='tab:blue',  label='BKM approx')
        ax.plot(phi_grid, bkm_full,   '--', lw=1.5, color='tab:blue',  label='BKM')

        # RGK fit (solid orange)
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