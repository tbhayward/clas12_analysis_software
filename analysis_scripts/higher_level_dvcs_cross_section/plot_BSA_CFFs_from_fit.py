#!/usr/bin/env python3
"""
plot_BSA_CFFs_from_fit.py

Usage:
    python plot_BSA_CFFs_from_fit.py output/fit_results/fit_results_<TIMESTAMP>.txt --CFFs [0|1]
"""
import os
import sys
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt
import ROOT

def parse_fit_results(fname):
    """Extract the flags, parameter names, and fit values from your text file."""
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
            break
    if vals is None or not pnames:
        raise RuntimeError("Couldn't parse fit file")
    return flags, pnames, vals

def load_all_bins(datafile):
    """Read imports/rga_prl_bsa.txt, split into φ-bins whenever φ wraps around."""
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
                arr["Q2m"], arr["xBm"], arr["tm"], arr["Ebm"] = (
                    arr["Q2"].mean(), arr["xB"].mean(),
                    arr["t"].mean(),  arr["Eb"].mean()
                )
                bins.append(arr)
                curr = {k: [] for k in curr}
            for k, v in zip(curr.keys(), (phi, Q2, xB, t, Eb, A, sigA)):
                curr[k].append(v)
            prev_phi = phi
    if curr["phi"]:
        arr = {k: np.array(v) for k, v in curr.items()}
        arr["Q2m"], arr["xBm"], arr["tm"], arr["Ebm"] = (
            arr["Q2"].mean(), arr["xB"].mean(),
            arr["t"].mean(),  arr["Eb"].mean()
        )
        bins.append(arr)
    return bins

# VGG defaults for all parameters
orig_defaults = {
    'renormImag':1.0, 'renormReal':1.0,
    'r_H':0.9,'alpha0_H':0.43,'alpha1_H':0.85,'n_H':1.35,'b_H':0.4,'Mm2_H':0.64,'P_H':1.0,
    'r_Ht':7.0,'alpha0_Ht':0.43,'alpha1_Ht':0.85,'n_Ht':0.6,'b_Ht':2.0,'Mm2_Ht':0.8,'P_Ht':1.0,
    'r_E':0.9,'alpha0_E':0.43,'alpha1_E':0.85,'n_E':1.35,'b_E':0.4,'Mm2_E':0.64,'P_E':1.0,
    'r_Et':1.0,'alpha0_Et':0.0,'alpha1_Et':0.0,'n_Et':0.0,'b_Et':0.0,'Mm2_Et':0.0,'P_Et':0.0,
}

def compute_bsa(phi_arr, Q2_arr, xB_arr, t_arr, Eb_arr, param_map, flags, tag=""):
    """
    Compile-set the C++ globals via ProcessLine(), then build BMK_DVCS and return dvcs.BSA().
    """
    # renormalizations
    ROOT.gInterpreter.ProcessLine(
        f"renormImag = {param_map.get('renormImag', orig_defaults['renormImag'])};"
    )
    ROOT.gInterpreter.ProcessLine(
        f"renormReal = {param_map.get('renormReal', orig_defaults['renormReal'])};"
    )

    # switch each CFF on/off & set its parameters
    for cff in ("H","Ht","E","Et"):
        ROOT.gInterpreter.ProcessLine(f"has{cff} = {int(flags[cff])};")
        if flags[cff]:
            for k in ("r", "alpha0", "alpha1", "n", "b", "Mm2", "P"):
                key = f"{k}_{cff}"
                v = param_map.get(key, orig_defaults[key])
                ROOT.gInterpreter.ProcessLine(f"{key} = {v};")

    bsas = []
    for i, (phi, Q2, xB, t, Eb) in enumerate(zip(
            phi_arr, Q2_arr, xB_arr, t_arr, Eb_arr)):
        dvcs = ROOT.BMK_DVCS(-1,1,0,Eb,xB,Q2,t,phi)
        mA   = dvcs.BSA()
        bsas.append(mA)
        if i<5:
            print(f"[{tag}] φ={phi:6.1f}°, ξ={dvcs.xi:.3f}, t={t:.3f}, BSA={mA:.4f}")
    return np.array(bsas)

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('fitfile', help='Fit results file')
    parser.add_argument('--CFFs', type=int, choices=[0,1], default=0,
                        help='0: use only ImH for fitted model, 1: show both ImH-only and full CFF models')
    args = parser.parse_args()

    m = re.search(r'fit_results_(\d{8}_\d{6})\.txt$', args.fitfile)
    if not m:
        print("ERROR: can't extract timestamp"); sys.exit(1)
    timestamp = m.group(1)

    flags, pnames, vals = parse_fit_results(args.fitfile)
    param_map = dict(zip(pnames, vals))
    print(">> Flags:", flags)
    print(">> Fitted parameters:", param_map)

    # Instead of loading a .so, compile the entire macro into Cling
    ROOT.gInterpreter.ProcessLine('#include "DVCS_xsec.C"')

    bins = load_all_bins('imports/rga_prl_bsa.txt')
    print(f">> Found {len(bins)} φ-bins")

    os.makedirs('output/plots', exist_ok=True)

    for idx, b in enumerate(bins, start=1):
        phi_data, As, sigAs = b["phi"], b["A"], b["sigA"]
        phi_grid = np.linspace(0,360,100)
        Q2g = np.full_like(phi_grid, b["Q2m"])
        xBg = np.full_like(phi_grid, b["xBm"])
        tg  = np.full_like(phi_grid, b["tm"])
        Ebg = np.full_like(phi_grid, b["Ebm"])

        # Compute original model with all CFFs using VGG defaults
        bsas_orig = compute_bsa(phi_grid, Q2g, xBg, tg, Ebg,
                                orig_defaults, flags, tag=f"bin{idx}-orig")

        # Compute fitted model with flags from fit
        bsas_fit = compute_bsa(phi_grid, Q2g, xBg, tg, Ebg,
                               param_map, flags, tag=f"bin{idx}-fit")

        fig, ax = plt.subplots(figsize=(8,5))
        ax.errorbar(phi_data, As, yerr=sigAs, fmt='o', ms=5,
                    color='k', label='Data')
        ax.plot(phi_grid, bsas_orig, '-',  lw=2,
                color='tab:blue', label='Original Model')

        if args.CFFs == 0:
            # Only show fitted model with ImH
            ax.plot(phi_grid, bsas_fit,  '--', lw=2,
                    color='tab:red',  label='Fitted Model (ImH only)')
        else:
            # Compute fitted model with only ImH (turn off other CFFs)
            flags_imh_only = flags.copy()
            flags_imh_only['Ht'] = 0
            flags_imh_only['E'] = 0
            flags_imh_only['Et'] = 0
            bsas_fit_imh = compute_bsa(phi_grid, Q2g, xBg, tg, Ebg,
                                       param_map, flags_imh_only, tag=f"bin{idx}-fit-ImH")
            
            # Plot both models
            ax.plot(phi_grid, bsas_fit_imh, '--', lw=2,
                    color='tab:red',  label='Fitted Model (ImH only)')
            ax.plot(phi_grid, bsas_fit,  '-.', lw=2,
                    color='tab:green',  label='Fitted Model (all CFFs)')

        ax.set_xlim(0,360)
        ax.set_xticks([0,60,120,180,240,300,360])
        ax.set_ylim(-0.6,0.6)
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

        fname = (f'output/plots/BSA_bin{idx:02d}_'
                 f'{timestamp}_Q2_{Q2m:.2f}_xB_{xBm:.3f}_t_{abs(tm):.3f}.pdf')
        fig.savefig(fname)
        print(f">> Saved bin {idx} plot to {fname}")
        plt.close(fig)

if __name__=='__main__':
    main()