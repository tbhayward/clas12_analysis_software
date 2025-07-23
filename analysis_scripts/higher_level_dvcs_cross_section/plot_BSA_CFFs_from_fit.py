#!/usr/bin/env python3
"""
plot_BSA_CFFs_from_fit.py

Usage:
    python plot_BSA_CFFs_from_fit.py output/fit_results/fit_results_<TIMESTAMP>.txt [--CFFs 0|1]

    --CFFs 0 : Plot fit with only ImH turned on (default)
    --CFFs 1 : Plot dashed red line (fit ImH only) and green dot-dashed (fit all CFFs)
"""

import os
import sys
import re
import numpy as np
import matplotlib.pyplot as plt
import argparse
import ROOT

def parse_fit_results(fname):
    """Extract fit parameter names, values, and flags from your text file."""
    with open(fname) as f:
        lines = [l.strip() for l in f if l.strip()]
    flag_line = next(l for l in lines if l.startswith("H "))
    toks = flag_line.split()
    flags = {toks[i]: int(toks[i+1]) for i in range(0, len(toks), 2)}
    pnames = []
    for i, l in enumerate(lines):
        if l.startswith("# parameters"):
            pnames = lines[i].split()[2:]
            break
    vals = None
    for i, l in enumerate(lines):
        if l.startswith("# values"):
            vals = list(map(float, lines[i+1].split()))
    if vals is None:
        raise RuntimeError("Could not parse fit-values from file")
    return flags, pnames, np.array(vals)

def load_all_bins(datafile):
    bins = []
    curr = {k: [] for k in ("phi","Q2","xB","t","Eb","A","sigA")}
    prev_phi = None
    with open(datafile) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            phi, Q2, xB, t, Eb, A, sigA = map(float, line.split())
            if prev_phi is not None and phi < prev_phi:
                arr = {k: np.array(v) for k, v in curr.items()}
                arr["Q2m"] = arr["Q2"].mean()
                arr["xBm"] = arr["xB"].mean()
                arr["tm"]  = arr["t"].mean()
                arr["Ebm"] = arr["Eb"].mean()
                bins.append(arr)
                curr = {k: [] for k in curr}
            for k, v in zip(("phi","Q2","xB","t","Eb","A","sigA"),
                            (phi, Q2, xB, t, Eb, A, sigA)):
                curr[k].append(v)
            prev_phi = phi
    if curr["phi"]:
        arr = {k: np.array(v) for k, v in curr.items()}
        arr["Q2m"] = arr["Q2"].mean()
        arr["xBm"] = arr["xB"].mean()
        arr["tm"]  = arr["t"].mean()
        arr["Ebm"] = arr["Eb"].mean()
        bins.append(arr)
    return bins

def compute_bsa(phi_arr, Q2_arr, xB_arr, t_arr, Eb_arr, params, flags, tag=""):
    keys = ["renormImag", "r_H", "alpha0_H", "alpha1_H", "n_H", "b_H", "Mm2_H", "P_H",
            "r_Ht", "alpha0_Ht", "alpha1_Ht", "n_Ht", "b_Ht", "Mm2_Ht", "P_Ht",
            "r_E", "alpha0_E", "alpha1_E", "n_E", "b_E", "Mm2_E", "P_E",
            "r_Et", "alpha0_Et", "alpha1_Et", "n_Et", "b_Et", "Mm2_Et", "P_Et"]
    p = dict(zip(keys, params))
    ROOT.renormImag = p["renormImag"]
    ROOT.renormReal = 1.0  # Not used for BSA, safe default
    ROOT.hasH  = bool(flags.get("H",0))
    ROOT.hasHt = bool(flags.get("Ht",0))
    ROOT.hasE  = bool(flags.get("E",0))
    ROOT.hasEt = bool(flags.get("Et",0))
    if ROOT.hasH:
        ROOT.r_H      = p["r_H"]
        ROOT.alpha0_H = p["alpha0_H"]
        ROOT.alpha1_H = p["alpha1_H"]
        ROOT.n_H      = p["n_H"]
        ROOT.b_H      = p["b_H"]
        ROOT.Mm2_H    = p["Mm2_H"]
        ROOT.P_H      = p["P_H"]
    if ROOT.hasHt:
        ROOT.r_Ht      = p["r_Ht"]
        ROOT.alpha0_Ht = p["alpha0_Ht"]
        ROOT.alpha1_Ht = p["alpha1_Ht"]
        ROOT.n_Ht      = p["n_Ht"]
        ROOT.b_Ht      = p["b_Ht"]
        ROOT.Mm2_Ht    = p["Mm2_Ht"]
        ROOT.P_Ht      = p["P_Ht"]
    if ROOT.hasE:
        ROOT.r_E      = p["r_E"]
        ROOT.alpha0_E = p["alpha0_E"]
        ROOT.alpha1_E = p["alpha1_E"]
        ROOT.n_E      = p["n_E"]
        ROOT.b_E      = p["b_E"]
        ROOT.Mm2_E    = p["Mm2_E"]
        ROOT.P_E      = p["P_E"]
    if ROOT.hasEt:
        ROOT.r_Et      = p["r_Et"]
        ROOT.alpha0_Et = p["alpha0_Et"]
        ROOT.alpha1_Et = p["alpha1_Et"]
        ROOT.n_Et      = p["n_Et"]
        ROOT.b_Et      = p["b_Et"]
        ROOT.Mm2_Et    = p["Mm2_Et"]
        ROOT.P_Et      = p["P_Et"]
    bsas = []
    for i, (phi, Q2, xB, t, Eb) in enumerate(zip(
            phi_arr, Q2_arr, xB_arr, t_arr, Eb_arr)):
        dvcs = ROOT.BMK_DVCS(-1, 1, 0, Eb, xB, Q2, t, phi)
        mA   = dvcs.BSA()
        bsas.append(mA)
        if i < 3:
            print(f"[{tag}] φ={phi:6.1f}°, ξ={dvcs.xi:.3f}, t={t:.3f}, BSA={mA:.4f}")
    return np.array(bsas)

def main():
    parser = argparse.ArgumentParser(description="Plot BSA from fit results, with optional CFF controls.")
    parser.add_argument("fitfile", type=str,
                        help="fit results file (output/fit_results/fit_results_<TIMESTAMP>.txt)")
    parser.add_argument("--CFFs", type=int, default=0,
                        help="0: show ImH-only fit (default). 1: also plot model with all fitted CFFs.")
    args = parser.parse_args()

    m_ts = re.search(r'fit_results_(\d{8}_\d{6})\.txt$', args.fitfile)
    if not m_ts:
        print("ERROR: can't extract timestamp from", args.fitfile)
        sys.exit(1)
    timestamp = m_ts.group(1)
    flags, pnames, vals = parse_fit_results(args.fitfile)
    print(">> Flags:", flags)
    print(">> Fitted parameters:", dict(zip(pnames, vals)))

    ROOT.gSystem.Load('./DVCS_xsec_C.so')
    datafile = 'imports/rga_prl_bsa.txt'
    bins = load_all_bins(datafile)
    print(f">> Found {len(bins)} φ-bins")
    outdir = 'output/plots'
    os.makedirs(outdir, exist_ok=True)

    # VGG default (original model), all CFFs on
    defaults = [
        1.0, 0.43, 0.85, 1.35, 0.4, 0.64, 1.0, 1.0,   # H
        7.0, 0.43, 0.85, 0.6, 2.0, 0.8, 1.0,          # Ht
        0.9, 0.43, 0.85, 1.35, 0.4, 0.64, 1.0,        # E
        0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0             # Et
    ]
    all_on_flags = {'H':1, 'Ht':1, 'E':1, 'Et':1}
    imh_flags = {'H': 1, 'Ht': 0, 'E': 0, 'Et': 0}

    for idx, b in enumerate(bins, start=1):
        phi_data = b["phi"]
        As, sigAs = b["A"], b["sigA"]
        phi_grid = np.linspace(0, 360, 100)
        Q2g = np.full_like(phi_grid, b["Q2m"])
        xBg = np.full_like(phi_grid, b["xBm"])
        tg  = np.full_like(phi_grid, b["tm"])
        Ebg = np.full_like(phi_grid, b["Ebm"])

        # Always plot original model (all CFFs on, VGG default)
        bsas_orig = compute_bsa(phi_grid, Q2g, xBg, tg, Ebg, defaults, all_on_flags, tag=f"bin{idx}-origVGG")

        # Dashed red: ImH-only, fit value
        bsas_fitH = None
        if args.CFFs in [0,1]:
            fitH = vals.copy()
            # Zero out non-H CFFs for fitH
            for key in ["r_Ht", "alpha0_Ht", "alpha1_Ht", "n_Ht", "b_Ht", "Mm2_Ht", "P_Ht",
                        "r_E", "alpha0_E", "alpha1_E", "n_E", "b_E", "Mm2_E", "P_E",
                        "r_Et", "alpha0_Et", "alpha1_Et", "n_Et", "b_Et", "Mm2_Et", "P_Et"]:
                if key in pnames:
                    fitH[pnames.index(key)] = 0.0
            bsas_fitH = compute_bsa(phi_grid, Q2g, xBg, tg, Ebg, fitH, imh_flags, tag=f"bin{idx}-fitH")

        # Green dot-dashed: all CFFs on, fit values
        bsas_fitAll = None
        if args.CFFs == 1:
            all_flags = {k: int(flags.get(k,0)) for k in ['H','Ht','E','Et']}
            bsas_fitAll = compute_bsa(phi_grid, Q2g, xBg, tg, Ebg, vals, all_flags, tag=f"bin{idx}-fitAll")

        fig, ax = plt.subplots(figsize=(8,5))
        ax.errorbar(phi_data, As, yerr=sigAs, fmt='o', ms=5, color='k', label='Data')
        ax.plot(phi_grid, bsas_orig, '-', lw=2, color='tab:blue', label='Original Model (VGG)')
        if bsas_fitH is not None:
            ax.plot(phi_grid, bsas_fitH, '--', lw=2, color='tab:red', label='Fit ImH only')
        if bsas_fitAll is not None:
            ax.plot(phi_grid, bsas_fitAll, '-.', lw=2, color='tab:green', label='Fit all CFFs')

        ax.set_xlim(0, 360)
        ax.set_xticks([0, 60, 120, 180, 240, 300, 360])
        ax.set_ylim(-0.6, 0.6)
        ax.set_xlabel(r'$\phi\;[\mathrm{deg}]$')
        ax.set_ylabel(r'$A_{LU}(\phi)$')
        Q2m, xBm, tm = b["Q2m"], b["xBm"], b["tm"]
        ax.set_title(
            (r'$\langle Q^2\rangle={:.2f}\,\mathrm{{GeV}}^2,\;'
             r'\langle x_B\rangle={:.3f},\;\langle -t\rangle={:.3f}\,\mathrm{{GeV}}^2$'
            ).format(Q2m, xBm, -tm),
            pad=12
        )
        ax.legend(loc='upper right', frameon=True, edgecolor='k')
        plt.tight_layout()
        fname = (f'{outdir}/BSA_bin{idx:02d}_'
                 f'{timestamp}_Q2_{Q2m:.2f}_xB_{xBm:.3f}_t_{abs(tm):.3f}.pdf')
        fig.savefig(fname)
        print(f">> Saved bin {idx} plot to {fname}")
        plt.close(fig)

if __name__ == '__main__':
    main()