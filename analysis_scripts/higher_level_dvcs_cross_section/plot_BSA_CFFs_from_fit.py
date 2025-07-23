#!/usr/bin/env python3
"""
plot_BSA_CFFs_from_fit.py

Usage:
    python plot_BSA_CFFs_from_fit.py output/fit_results/fit_results_<TIMESTAMP>.txt [--CFFs 0|1]

Options:
    --CFFs 0 : Plot "Fitted Model" using only ImH on (all others off)
    --CFFs 1 : Plot "Fitted Model" using all fitted CFFs that are enabled in the fit file

This script plots BSA data and models for each φ-bin using your fit parameters.
"""
import os
import sys
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt
import ROOT

# ------------------- Parse Command-Line Arguments ----------------------
parser = argparse.ArgumentParser(description="Plot BSA with fit CFFs.")
parser.add_argument('fitfile', help="fit results file (output/fit_results/fit_results_*.txt)")
parser.add_argument('--CFFs', type=int, default=1,
                    help="0: plot fitted model with only ImH; 1: plot fitted model with all fitted CFFs (default=1)")
args = parser.parse_args()

# ------------------- Parse Fit File ------------------------------------
def parse_fit_file(fname):
    with open(fname) as f:
        lines = [l.strip() for l in f if l.strip()]
    # Get flags
    flag_line = next(l for l in lines if l.startswith("H "))
    toks = flag_line.split()
    flags = {toks[i]: int(toks[i+1]) for i in range(0, len(toks), 2)}
    # Parameter names
    pnames = []
    for l in lines:
        if l.startswith("# parameters"):
            pnames = l.split()[2:]
            break
    # Values
    vals = None
    for i, l in enumerate(lines):
        if l.startswith("# values"):
            vals = list(map(float, lines[i+1].split()))
            break
    if vals is None:
        raise RuntimeError("No '# values' block in fit file!")
    params = dict(zip(pnames, vals))
    return flags, params

flags, fit_params = parse_fit_file(args.fitfile)
print(">> Flags:", flags)
print(">> Fitted parameters:", fit_params)

# ------------------- Load Data and Bin Info ----------------------------
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

# ------------------- Default VGG Model Parameters ----------------------
defaults = {
    "renormImag": 1.0,
    "r_H": 0.9, "alpha0_H": 0.43, "alpha1_H": 0.85, "n_H": 1.35, "b_H": 0.4, "Mm2_H": 0.64, "P_H": 1.0,
    "r_Ht": 7.0, "alpha0_Ht": 0.43, "alpha1_Ht": 0.85, "n_Ht": 0.6, "b_Ht": 2.0, "Mm2_Ht": 0.8, "P_Ht": 1.0,
    "r_E": 0.9, "alpha0_E": 0.43, "alpha1_E": 0.85, "n_E": 1.35, "b_E": 0.4, "Mm2_E": 0.64, "P_E": 1.0,
    "r_Et": 1.0, "alpha0_Et": 0.43, "alpha1_Et": 0.85, "n_Et": 1.35, "b_Et": 0.4, "Mm2_Et": 0.64, "P_Et": 1.0,
    "renormReal": 1.0
}

# Helper: CFF name blocks and corresponding flag
CFF_blocks = [
    ("H",   ["r_H", "alpha0_H", "alpha1_H", "n_H", "b_H", "Mm2_H", "P_H"]),
    ("Ht",  ["r_Ht", "alpha0_Ht", "alpha1_Ht", "n_Ht", "b_Ht", "Mm2_Ht", "P_Ht"]),
    ("E",   ["r_E", "alpha0_E", "alpha1_E", "n_E", "b_E", "Mm2_E", "P_E"]),
    ("Et",  ["r_Et", "alpha0_Et", "alpha1_Et", "n_Et", "b_Et", "Mm2_Et", "P_Et"])
]

# ------------------- PyROOT CFF Setter Utility -------------------------
def set_CFFs_in_ROOT(params, active_flags):
    # Set renormImag, renormReal first
    if hasattr(ROOT, "renormImag"): ROOT.renormImag = params.get("renormImag", 1.0)
    if hasattr(ROOT, "renormReal"): ROOT.renormReal = params.get("renormReal", 1.0)
    # Set CFFs and their on/off flags
    for blockname, blockvars in CFF_blocks:
        cff_flag = bool(active_flags.get(blockname, 0))
        if hasattr(ROOT, f"has{blockname}"): setattr(ROOT, f"has{blockname}", cff_flag)
        if cff_flag:
            for v in blockvars:
                if hasattr(ROOT, v):
                    try:
                        setattr(ROOT, v, params[v])
                    except KeyError:
                        pass # safe to skip missing

# ------------------- Compute BSA Utility -------------------------------
def compute_bsa(phi_arr, Q2_arr, xB_arr, t_arr, Eb_arr, params, flags, tag=""):
    set_CFFs_in_ROOT(params, flags)
    bsas = []
    for i, (phi, Q2, xB, t, Eb) in enumerate(zip(phi_arr, Q2_arr, xB_arr, t_arr, Eb_arr)):
        dvcs = ROOT.BMK_DVCS(-1, 1, 0, Eb, xB, Q2, t, phi)
        mA   = dvcs.BSA()
        bsas.append(mA)
        if i < 3:
            print(f"[{tag}] φ={phi:6.1f}°, ξ={dvcs.xi:.3f}, t={t:.3f}, BSA={mA:.4f}")
    return np.array(bsas)

# ------------------- Main Plotting Routine -----------------------------
def main():
    m = re.search(r'_(\d{8}_\d{6})\.txt$', args.fitfile)
    if not m:
        print("ERROR: can't extract timestamp from", args.fitfile)
        sys.exit(1)
    timestamp = m.group(1)
    ROOT.gSystem.Load('./DVCS_xsec_C.so')

    datafile = 'imports/rga_prl_bsa.txt'
    bins = load_all_bins(datafile)
    print(f">> Found {len(bins)} φ-bins")

    outdir = 'output/plots'
    os.makedirs(outdir, exist_ok=True)

    # For original model, always turn on all CFFs
    all_on_flags = {cff: 1 for cff, _ in CFF_blocks}
    defaults_used = {k: defaults[k] for k in defaults if k in defaults}

    # For ImH only case (fitted), turn on only ImH, set all other flags off
    onlyH_flags = {"H": 1, "Ht": 0, "E": 0, "Et": 0}

    for idx, b in enumerate(bins, start=1):
        phi_data = b["phi"]
        As, sigAs = b["A"], b["sigA"]
        phi_grid = np.linspace(0, 360, 100)
        Q2g = np.full_like(phi_grid, b["Q2m"])
        xBg = np.full_like(phi_grid, b["xBm"])
        tg  = np.full_like(phi_grid, b["tm"])
        Ebg = np.full_like(phi_grid, b["Ebm"])

        # (1) Original model: all VGG defaults, all CFFs on
        bsas_orig = compute_bsa(phi_grid, Q2g, xBg, tg, Ebg, defaults, all_on_flags, tag=f"bin{idx}-origVGG")

        # (2) Fitted model(s)
        # Choose flags for fit model:
        if args.CFFs == 0:
            # Fitted line: only ImH, all others off
            fit_flags = onlyH_flags
            # Fit parameters: renormImag + H only, other params fallback to VGG defaults
            fitH_params = {k: fit_params[k] if k in fit_params else defaults[k] for k in defaults}
            for cff, block in CFF_blocks:
                if cff != "H":
                    for v in block:
                        fitH_params[v] = defaults[v]
            bsas_fit = compute_bsa(phi_grid, Q2g, xBg, tg, Ebg, fitH_params, fit_flags, tag=f"bin{idx}-fitH")
            # Fitted model with all CFFs (for overplot): turn on all CFFs that were fitted
            all_fit_flags = {cff: int(flags[cff]) for cff, _ in CFF_blocks}
            bsas_fit_all = compute_bsa(phi_grid, Q2g, xBg, tg, Ebg, fit_params, all_fit_flags, tag=f"bin{idx}-fitALL")
        else:
            # Only plot the model with all fitted CFFs
            all_fit_flags = {cff: int(flags[cff]) for cff, _ in CFF_blocks}
            bsas_fit = compute_bsa(phi_grid, Q2g, xBg, tg, Ebg, fit_params, all_fit_flags, tag=f"bin{idx}-fitALL")
            bsas_fit_all = None

        # Plot
        fig, ax = plt.subplots(figsize=(8,5))
        ax.errorbar(phi_data, As, yerr=sigAs, fmt='o', ms=5, color='k', label='Data')
        ax.plot(phi_grid, bsas_orig, '-', lw=2, color='tab:blue', label='Original Model')
        if args.CFFs == 0:
            ax.plot(phi_grid, bsas_fit, '--', lw=2, color='tab:red', label='Fit: Only ImH')
            ax.plot(phi_grid, bsas_fit_all, '-.', lw=2, color='tab:green', label='Fit: All CFFs')
        else:
            ax.plot(phi_grid, bsas_fit, '--', lw=2, color='tab:red', label='Fit: All CFFs')
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