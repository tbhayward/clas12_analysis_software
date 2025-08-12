#!/usr/bin/env python3
"""
plot_BSA_CFFs_from_fit.py

Usage:
    python plot_BSA_CFFs_from_fit.py output/fit_results/fit_results_<TIMESTAMP>.txt --CFFs [0|1]

Reads the fit-results text file produced by fit_CFFs.C (now includes an 'input' line),
loads that BSA file, and plots per-φ bin:
  - data with error bars
  - "original/default" model curve
  - fitted model curve(s): ImH-only (default) or ImH-only + full CFFs (if --CFFs 1)

Assumes input BSA files have columns:
# phi(deg) q2(GeV2) xb t(GeV2) Eb(GeV) A sigA
with t being the physical (negative) Mandelstam t (so -t > 0).
"""
import os
import sys
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt
import ROOT

# --------------------------------------------------------------------------------------
# Parse fit results (flags, parameter names/values, and 'input' BSA path)
# --------------------------------------------------------------------------------------
def parse_fit_results(fname):
    with open(fname) as f:
        lines = [l.rstrip("\n") for l in f if l.strip()]

    # flags line (e.g. "H 1 Ht 0 E 0 Et 0")
    flag_line = next((l for l in lines if re.match(r'^\s*H\s+\d+', l)), None)
    if flag_line is None:
        raise RuntimeError("Couldn't find flags line in fit file")
    toks = flag_line.split()
    flags = {}
    for i in range(0, len(toks), 2):
        key = toks[i]
        try:
            flags[key] = int(toks[i+1])
        except Exception:
            pass

    # parameter names / values / errors
    pnames = []
    vals = None
    for i, l in enumerate(lines):
        if l.startswith("# parameters"):
            pnames = l.split()[2:]
        elif l.startswith("# values:"):
            if i + 1 < len(lines):
                vals = [float(x) for x in lines[i+1].split()]

    if not pnames or vals is None:
        raise RuntimeError("Couldn't parse parameter names/values from fit file")

    # input path line: "input       <path>"
    input_line = next((l for l in lines if l.startswith("input")), None)
    input_path = None
    if input_line is not None:
        parts = input_line.split()
        # e.g., ["input", "imports/rgk_preliminary_bsa.txt"]
        if len(parts) >= 2:
            input_path = parts[1]

    return flags, pnames, vals, input_path

# --------------------------------------------------------------------------------------
# Read BSA file and split into φ-bins (whenever φ wraps around)
# --------------------------------------------------------------------------------------
def load_all_bins(datafile):
    """
    Read a BSA text file with columns:
    phi  Q2  xB  t  Eb  A  sigA
    Split into bins when phi wraps (decreases). 't' is negative; we keep it that way.
    """
    bins = []
    curr = {k: [] for k in ("phi","Q2","xB","t","Eb","A","sigA")}
    prev_phi = None
    with open(datafile) as f:
        for line in f:
            if (not line.strip()) or line.lstrip().startswith("#"):
                continue
            phi, Q2, xB, t, Eb, A, sigA = map(float, line.split())
            if prev_phi is not None and phi < prev_phi:
                arr = {k: np.array(v) for k, v in curr.items()}
                arr["Q2m"] = arr["Q2"].mean()
                arr["xBm"] = arr["xB"].mean()
                arr["tm"]  = arr["t"].mean()   # negative (physical t)
                arr["Ebm"] = arr["Eb"].mean()
                bins.append(arr)
                curr = {k: [] for k in curr}
            for k, v in zip(curr.keys(), (phi, Q2, xB, t, Eb, A, sigA)):
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

# --------------------------------------------------------------------------------------
# Default model parameters (match C++ globals; include r_* and correct M2_* names)
# --------------------------------------------------------------------------------------
orig_defaults = {
    'renormImag': 1.0,
    'renormReal': 1.0,

    # H
    'r_H': 0.9, 'n_H': 1.25, 'alpha0_H': 0.43, 'alpha1_H': 0.85,
    'b_H': 0.4, 'M2_H': 0.64, 'P_H': 1.0,

    # Htilde
    'r_Ht': 7.0, 'n_Ht': 0.6, 'alpha0_Ht': 0.43, 'alpha1_Ht': 0.85,
    'b_Ht': 2.0, 'M2_Ht': 0.8, 'P_Ht': 1.0,

    # E
    'r_E': 0.9, 'n_E': 1.25, 'alpha0_E': 0.43, 'alpha1_E': 0.85,
    'b_E': 0.4, 'M2_E': 0.64, 'P_E': 1.0,

    # Etilde
    'r_Et': 1.0, 'n_Et': 0.6, 'alpha0_Et': 0.0, 'alpha1_Et': 0.0,
    'b_Et': 0.0, 'M2_Et': 0.0, 'P_Et': 0.0,
}

# --------------------------------------------------------------------------------------
# Compute BSA by setting C++ globals then calling BMK_DVCS
# --------------------------------------------------------------------------------------
def compute_bsa(phi_arr, Q2_arr, xB_arr, t_arr, Eb_arr, param_map, flags, tag=""):
    """
    Sets ROOT globals then returns BSA(phi) from BMK_DVCS.
    Expects 't_arr' to be negative (physical t). BMK_DVCS signature: (.., Eb, xB, Q2, t, phi)
    """
    # renormalizations
    ROOT.gInterpreter.ProcessLine(f"renormImag = {param_map.get('renormImag', orig_defaults['renormImag'])};")
    ROOT.gInterpreter.ProcessLine(f"renormReal = {param_map.get('renormReal', orig_defaults['renormReal'])};")

    # which CFFs are active
    for cff in ("H","Ht","E","Et"):
        ROOT.gInterpreter.ProcessLine(f"has{cff} = {int(flags.get(cff,0))};")

    # set parameters for active CFFs
    for cff in ("H","Ht","E","Et"):
        if not flags.get(cff,0):
            continue
        # include r, n, alpha0, alpha1, b, M2, P
        for k in ("r","n","alpha0","alpha1","b","M2","P"):
            key = f"{k}_{cff}"
            val = param_map.get(key, orig_defaults.get(key))
            if val is None:
                continue
            ROOT.gInterpreter.ProcessLine(f"{key} = {val};")

    bsas = []
    for i, (phi, Q2, xB, t, Eb) in enumerate(zip(phi_arr, Q2_arr, xB_arr, t_arr, Eb_arr)):
        dvcs = ROOT.BMK_DVCS(-1, 1, 0, Eb, xB, Q2, t, phi)
        mA   = dvcs.BSA()
        bsas.append(mA)
        if i < 3:
            try:
                xi = dvcs.xi
            except Exception:
                xi = np.nan
            print(f"[{tag}] φ={phi:6.1f}°, ξ={xi:.3f}, t={t:.3f}, BSA={mA:.4f}")
    return np.array(bsas)

# --------------------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('fitfile', help='Fit results file from fit_CFFs.C')
    parser.add_argument('--CFFs', type=int, choices=[0,1], default=0,
                        help='0: show ImH-only fit; 1: show ImH-only AND full-CFF fit')
    parser.add_argument('--override-data', default=None,
                        help='Optional path to BSA data file to override the one in fitfile')
    args = parser.parse_args()

    m = re.search(r'fit_results_(\d{8}_\d{6})\.txt$', args.fitfile)
    if not m:
        print("ERROR: cannot extract timestamp from fitfile name"); sys.exit(1)
    timestamp = m.group(1)

    flags, pnames, vals, input_path = parse_fit_results(args.fitfile)
    param_map = dict(zip(pnames, vals))

    # Resolve which BSA file to load
    datafile = args.override-data if args.override_data else input_path
    if not datafile:
        # fallback (old default)
        datafile = 'imports/rga_prl_bsa.txt'
        print("WARNING: no 'input' line found in fitfile; falling back to", datafile)
    if not os.path.isfile(datafile):
        print("ERROR: data file not found:", datafile)
        sys.exit(1)

    print(">> Flags:", flags)
    print(">> Fitted parameters (subset shown):", {k: param_map[k] for k in pnames[:min(10,len(pnames))]}, "...")
    print(">> Using BSA data file:", datafile)

    # Make sure DVCS_xsec.C is loaded
    ROOT.gInterpreter.ProcessLine('#include "DVCS_xsec.C"')

    bins = load_all_bins(datafile)
    print(f">> Found {len(bins)} φ-bins")
    if not bins:
        print("ERROR: no bins parsed from data file (check format)"); sys.exit(1)

    os.makedirs('output/plots', exist_ok=True)

    for idx, b in enumerate(bins, start=1):
        phi_data, As, sigAs = b["phi"], b["A"], b["sigA"]

        # grid at same kinematics (mean per bin)
        phi_grid = np.linspace(0, 360, 181)
        Q2g = np.full_like(phi_grid, b["Q2m"])
        xBg = np.full_like(phi_grid, b["xBm"])
        tg  = np.full_like(phi_grid, b["tm"])   # negative t
        Ebg = np.full_like(phi_grid, b["Ebm"])

        # Original/default model
        bsas_orig = compute_bsa(phi_grid, Q2g, xBg, tg, Ebg,
                                orig_defaults, flags, tag=f"bin{idx}-orig")

        # Fitted model (all params present in fit file)
        bsas_fit_all = compute_bsa(phi_grid, Q2g, xBg, tg, Ebg,
                                   param_map, flags, tag=f"bin{idx}-fit-all")

        fig, ax = plt.subplots(figsize=(8,5))
        ax.errorbar(phi_data, As, yerr=sigAs, fmt='o', ms=5, color='k', label='Data')
        ax.plot(phi_grid, bsas_orig, '-',  lw=2, color='tab:blue',  label='Default model')

        if args.CFFs == 0:
            # Show only ImH fit (others off)
            flags_imh_only = flags.copy()
            for other in ("Ht","E","Et"):
                flags_imh_only[other] = 0
            bsas_fit_imh = compute_bsa(phi_grid, Q2g, xBg, tg, Ebg,
                                       param_map, flags_imh_only, tag=f"bin{idx}-fit-ImH")
            ax.plot(phi_grid, bsas_fit_imh, '--', lw=2, color='tab:red', label='Fitted (ImH)')
        else:
            # Show ImH-only and full
            flags_imh_only = flags.copy()
            for other in ("Ht","E","Et"):
                flags_imh_only[other] = 0
            bsas_fit_imh = compute_bsa(phi_grid, Q2g, xBg, tg, Ebg,
                                       param_map, flags_imh_only, tag=f"bin{idx}-fit-ImH")
            ax.plot(phi_grid, bsas_fit_imh, '--', lw=2, color='tab:red',   label='Fitted (ImH)')
            ax.plot(phi_grid, bsas_fit_all, '-.', lw=2, color='tab:green', label='Fitted (all CFFs)')

        ax.set_xlim(0, 360)
        ax.set_xticks([0,60,120,180,240,300,360])
        ax.set_ylim(-0.6, 0.6)
        ax.set_xlabel(r'$\phi\;[\mathrm{deg}]$')
        ax.set_ylabel(r'$A_{LU}(\phi)$')

        Q2m, xBm, tm = b["Q2m"], b["xBm"], b["tm"]  # tm is negative
        ax.set_title(
            (r'$\langle Q^2\rangle={:.2f}\,\mathrm{{GeV}}^2,\;'
             r'\langle x_B\rangle={:.3f},\;\langle -t\rangle={:.3f}\,\mathrm{{GeV}}^2$').format(
                Q2m, xBm, -tm),
            pad=12
        )

        ax.legend(loc='upper right', frameon=True, edgecolor='k')
        plt.tight_layout()

        fname = (f'output/plots/BSA_bin{idx:02d}_'
                 f'{timestamp}_Q2_{Q2m:.2f}_xB_{xBm:.3f}_mt_{(-tm):.3f}.pdf')
        fig.savefig(fname)
        print(f">> Saved bin {idx} plot to {fname}")
        plt.close(fig)

if __name__ == '__main__':
    main()