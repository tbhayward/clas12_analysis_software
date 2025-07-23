#!/usr/bin/env python3
"""
plot_BSA_CFFs_from_fit.py

Usage:
    python plot_BSA_CFFs_from_fit.py [--CFFs {0,1}] output/fit_results/fit_results_<TIMESTAMP>.txt

Options:
    --CFFs 0   Show only ImH‐only original prediction
    --CFFs 1   Show both ImH‐only fitted (red dashed) and all-CFFs fitted (green dash-dot)
"""
import os
import sys
import re
import argparse

import numpy as np
import matplotlib.pyplot as plt
import ROOT

# VGG defaults for Im‐parts + renorm
orig_defaults = {
    'renormImag':1.0, 'renormReal':1.0,
    # H
    'r_H':0.9,'alpha0_H':0.43,'alpha1_H':0.85,
    'n_H':1.35,'b_H':0.4,'Mm2_H':0.64,'P_H':1.0,
    # Ht
    'r_Ht':7.0,'alpha0_Ht':0.43,'alpha1_Ht':0.85,
    'n_Ht':0.6,'b_Ht':2.0,'Mm2_Ht':0.8,'P_Ht':1.0,
    # E
    'r_E':0.9,'alpha0_E':0.43,'alpha1_E':0.85,
    'n_E':1.35,'b_E':0.4,'Mm2_E':0.64,'P_E':1.0,
    # Et
    'r_Et':1.0,'alpha0_Et':0.0,'alpha1_Et':0.0,
    'n_Et':0.0,'b_Et':0.0,'Mm2_Et':0.0,'P_Et':0.0,
}
shape_keys = ["r","alpha0","alpha1","n","b","Mm2","P"]

def parse_fit_results(fname):
    """Read flags, parameter names, and fit values."""
    with open(fname) as f:
        lines = [l.strip() for l in f if l.strip()]
    # flags line: e.g. "H 1  Ht 1  E 0  Et 1"
    flag_line = next(l for l in lines if l.startswith("H "))
    toks = flag_line.split()
    flags = {toks[i]: int(toks[i+1]) for i in range(0, len(toks), 2)}
    # parameter names
    for l in lines:
        if l.startswith("# parameters"):
            pnames = l.split()[2:]
            break
    # parameter values
    for i, l in enumerate(lines):
        if l.startswith("# values"):
            vals = list(map(float, lines[i+1].split()))
            break
    return flags, pnames, vals

def load_all_bins(datafile):
    """Split imports/rga_prl_bsa.txt into φ‐bins."""
    bins = []
    curr = {k: [] for k in ("phi","Q2","xB","t","Eb","A","sigA")}
    prev_phi = None
    with open(datafile) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            phi, Q2, xB, t, Eb, A, sigA = map(float, line.split())
            if prev_phi is not None and phi < prev_phi:
                arr = {k: np.array(v) for k,v in curr.items()}
                arr["Q2m"] = arr["Q2"].mean()
                arr["xBm"] = arr["xB"].mean()
                arr["tm"]  = arr["t"].mean()
                arr["Ebm"] = arr["Eb"].mean()
                bins.append(arr)
                curr = {k: [] for k in curr}
            for k,v in zip(curr.keys(), (phi,Q2,xB,t,Eb,A,sigA)):
                curr[k].append(v)
            prev_phi = phi
    if curr["phi"]:
        arr = {k: np.array(v) for k,v in curr.items()}
        arr["Q2m"] = arr["Q2"].mean()
        arr["xBm"] = arr["xB"].mean()
        arr["tm"]  = arr["t"].mean()
        arr["Ebm"] = arr["Eb"].mean()
        bins.append(arr)
    return bins

def compute_bsa(phi_arr, Q2_arr, xB_arr, t_arr, Eb_arr,
                param_map, flags, tag=""):
    """
    Set the C++ globals via Cling, then build BMK_DVCS & return dvcs.BSA().
    """
    interp = ROOT.gInterpreter
    # renormalizations
    interp.ProcessLine(f"renormImag = {param_map.get('renormImag',orig_defaults['renormImag'])};")
    interp.ProcessLine(f"renormReal = {param_map.get('renormReal',orig_defaults['renormReal'])};")
    # each CFF on/off + its shape‐params
    for cff in ("H","Ht","E","Et"):
        interp.ProcessLine(f"has{cff} = {int(flags[cff])};")
        if flags[cff]:
            for k in shape_keys:
                key = f"{k}_{cff}"
                v = param_map.get(key, orig_defaults[key])
                interp.ProcessLine(f"{key} = {v};")
    # now loop
    out = []
    for i, (phi,Q2,xB,t,Eb) in enumerate(zip(phi_arr,Q2_arr,xB_arr,t_arr,Eb_arr)):
        dvcs = ROOT.BMK_DVCS(-1,1,0,Eb,xB,Q2,t,phi)
        val  = dvcs.BSA()
        out.append(val)
        if i<5:
            print(f"[{tag}] φ={phi:.1f}°, ξ={dvcs.xi:.3f}, t={t:.3f}, BSA={val:.4f}")
    return np.array(out)

def main():
    p = argparse.ArgumentParser(
        description="Plot BSA vs φ for either ImH‐only or all‐CFF fits"
    )
    p.add_argument("fitfile",
                   help="output/fit_results/fit_results_<TIMESTAMP>.txt")
    p.add_argument("--CFFs", type=int, choices=[0,1], default=1,
                   help="0: ImH-only prediction; 1: ImH-only fit + all-CFF fit")
    args = p.parse_args()

    # parse fit file
    flags, pnames, vals = parse_fit_results(args.fitfile)
    fit_map = dict(zip(pnames, vals))
    print(">> Flags:", flags)
    print(">> Fitted parameters:", fit_map)

    # load the full macro into Cling so all globals exist
    ROOT.gInterpreter.ProcessLine('#include "DVCS_xsec.C"')

    bins = load_all_bins("imports/rga_prl_bsa.txt")
    print(f">> Found {len(bins)} φ-bins")

    os.makedirs("output/plots", exist_ok=True)

    # prepare two flag‐sets for H‐only vs all
    flags_all    = flags
    flags_H_only = flags.copy()
    flags_H_only.update({'Ht':0,'E':0,'Et':0})

    for idx, b in enumerate(bins, start=1):
        phi_data, As, sigAs = b["phi"], b["A"], b["sigA"]
        phi_grid = np.linspace(0,360,100)
        Q2g = np.full_like(phi_grid, b["Q2m"])
        xBg = np.full_like(phi_grid, b["xBm"])
        tg  = np.full_like(phi_grid, b["tm"])
        Ebg = np.full_like(phi_grid, b["Ebm"])

        # always compute the ImH-only **original** prediction, for CFFs=0
        bsas_orig_H = compute_bsa(phi_grid, Q2g, xBg, tg, Ebg,
                                  orig_defaults, flags_H_only, tag=f"bin{idx}-origH")

        # if CFFs=1, also compute the two fitted curves
        if args.CFFs == 1:
            bsas_fit_H = compute_bsa(phi_grid, Q2g, xBg, tg, Ebg,
                                     fit_map, flags_H_only, tag=f"bin{idx}-fitH")
            bsas_fit_all = compute_bsa(phi_grid, Q2g, xBg, tg, Ebg,
                                       fit_map, flags_all,  tag=f"bin{idx}-fitAll")

        # now plot
        fig, ax = plt.subplots(figsize=(8,5))
        ax.errorbar(phi_data, As, yerr=sigAs, fmt='o', ms=5,
                    color='k', label='Data')

        if args.CFFs == 0:
            ax.plot(phi_grid, bsas_orig_H, '-', lw=2,
                    color='tab:blue', label='ImH-only₀ model')
        else:
            ax.plot(phi_grid, bsas_fit_H, '--', lw=2,
                    color='tab:red',    label='ImH-only fit')
            ax.plot(phi_grid, bsas_fit_all, '-.', lw=2,
                    color='tab:green',  label='All-CFFs fit')

        ax.set_xlim(0,360)
        ax.set_xticks([0,60,120,180,240,300,360])
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

        outname = (f"output/plots/BSA_bin{idx:02d}_"
                   f"{re.search(r'_(\d{{8}}_\d{{6}})\\.txt$', args.fitfile).group(1)}"
                   f"_Q2_{Q2m:.2f}_xB_{xBm:.3f}_t_{abs(tm):.3f}.pdf")
        fig.savefig(outname)
        print(f">> Saved bin {idx} → {outname}")
        plt.close(fig)

if __name__ == "__main__":
    main()