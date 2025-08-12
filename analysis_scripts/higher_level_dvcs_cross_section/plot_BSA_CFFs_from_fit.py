#!/usr/bin/env python3
"""
plot_BSA_CFFs_from_fit.py

Usage:
    python plot_BSA_CFFs_from_fit.py output/fit_results/fit_results_<TIMESTAMP>.txt --CFFs [0|1]

What this does
--------------
• Parses the fit-results text file produced by fit_CFFs.C
  - Reads which CFFs were enabled (H, Ht, E, Et)
  - Reads the fitted parameter names and values
  - Reads the original BSA input file path from the "input" line
• Loads DVCS_xsec.C into ROOT exactly once (guarded) and pushes the flags/parameters
• Reconstructs BSA(φ) per φ-bin from the original BSA input file:
  - Plots data points with error bars
  - Overlays "original-model" BSA and "fitted" BSA curves
  - If --CFFs=1, also shows the ImH-only fitted curve separately from the "all enabled CFFs" fit
• Saves one PDF per φ-bin into output/plots/
"""

import os
import sys
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt

import ROOT  # PyROOT (cling)

# -----------------------------------------------------------------------------
# One-time, guarded load of DVCS_xsec.C
# -----------------------------------------------------------------------------
def _load_dvcs_once():
    """Load DVCS_xsec.C into the ROOT interpreter exactly once, guarded."""
    if getattr(_load_dvcs_once, "_done", False):
        return
    ROOT.gInterpreter.Declare(r"""
#ifndef DVCS_XSEC_ONCE_GUARD
#define DVCS_XSEC_ONCE_GUARD
#include "DVCS_xsec.C"
#endif
""")
    _load_dvcs_once._done = True


# -----------------------------------------------------------------------------
# Fit-file parsing
# -----------------------------------------------------------------------------
def parse_fit_results(fname):
    """
    Extract:
      - flags dict like {'H':1,'Ht':0,'E':0,'Et':0}
      - list of parameter names
      - list of parameter values
      - input data file path from 'input       <path>' line (if present)
      - strategy (int) if present
    """
    with open(fname) as f:
        lines = [l.strip() for l in f if l.strip()]

    # Flags line, e.g. "H 1 Ht 0 E 0 Et 0"
    flag_line = None
    for l in lines:
        if re.match(r'^\s*H\s+[01]\s+Ht\s+[01]\s+E\s+[01]\s+Et\s+[01]\s*$', l):
            flag_line = l
            break
    if flag_line is None:
        # fallback: first line starting with "H "
        flag_line = next((l for l in lines if l.startswith("H ")), None)
    if flag_line is None:
        raise RuntimeError("Could not find 'H ... Ht ... E ... Et ...' flags line in fit file")

    toks = flag_line.split()
    flags = {}
    for i in range(0, len(toks), 2):
        key = toks[i]
        val = int(toks[i + 1])
        flags[key] = val

    # Parameter names, values, (errors not needed here)
    pnames = []
    pvals = None
    for i, l in enumerate(lines):
        if l.startswith("# parameters"):
            pnames = l.split()[2:]
        elif l.startswith("# values"):
            if i + 1 < len(lines):
                pvals = list(map(float, lines[i + 1].split()))

    if not pnames or pvals is None:
        raise RuntimeError("Couldn't parse parameter names/values from fit file")

    # Optional: input data file path (from fit_CFFs.C output 'input       <path>')
    input_file = None
    for l in lines:
        if l.startswith("input"):
            # everything after 'input' (split once) is the path
            parts = l.split(None, 1)
            if len(parts) == 2:
                input_file = parts[1].strip()
            break

    # Optional: strategy
    strategy = None
    for l in lines:
        if l.startswith("strategy"):
            parts = l.split()
            if len(parts) >= 2:
                try:
                    strategy = int(parts[1])
                except Exception:
                    pass
            break

    return flags, pnames, pvals, input_file, strategy


# -----------------------------------------------------------------------------
# Read BSA input file and split into φ-bins (wrapping when φ decreases)
# -----------------------------------------------------------------------------
def load_all_bins(datafile):
    """
    Read a BSA data file with columns:
      phi  Q2  xB  t  Eb  A  sigA
    Split into bins whenever phi wraps around (phi decreases).
    Returns a list of dicts with numpy arrays and per-bin means.
    """
    bins = []
    curr = {k: [] for k in ("phi", "Q2", "xB", "t", "Eb", "A", "sigA")}
    prev_phi = None
    with open(datafile) as f:
        for line in f:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            phi, Q2, xB, t, Eb, A, sigA = map(float, line.split())
            if prev_phi is not None and phi < prev_phi:
                arr = {k: np.array(v) for k, v in curr.items()}
                arr["Q2m"], arr["xBm"], arr["tm"], arr["Ebm"] = (
                    arr["Q2"].mean(), arr["xB"].mean(), arr["t"].mean(), arr["Eb"].mean()
                )
                bins.append(arr)
                curr = {k: [] for k in curr}
            for k, v in zip(curr.keys(), (phi, Q2, xB, t, Eb, A, sigA)):
                curr[k].append(v)
            prev_phi = phi
    if curr["phi"]:
        arr = {k: np.array(v) for k, v in curr.items()}
        arr["Q2m"], arr["xBm"], arr["tm"], arr["Ebm"] = (
            arr["Q2"].mean(), arr["xB"].mean(), arr["t"].mean(), arr["Eb"].mean()
        )
        bins.append(arr)
    return bins


# -----------------------------------------------------------------------------
# "Original" model defaults (match DVCS_xsec.C global variable names)
# -----------------------------------------------------------------------------
def original_model_defaults():
    """
    Return a dict of baseline parameters keyed exactly as in DVCS_xsec.C:
      r_H, n_H, alpha0_H, alpha1_H, b_H, M2_H, P_H,
      r_Ht, n_Ht, alpha0_Ht, ...
      r_E, n_E, ...
      r_Et, n_Et, ...
      renormImag, renormReal
    """
    return {
        'renormImag': 1.0,
        'renormReal': 1.0,
        # H
        'r_H': 0.9, 'n_H': 1.25, 'alpha0_H': 0.43, 'alpha1_H': 0.85, 'b_H': 0.4, 'M2_H': 0.64, 'P_H': 1.0,
        # Htilde
        'r_Ht': 7.0, 'n_Ht': 0.6, 'alpha0_Ht': 0.43, 'alpha1_Ht': 0.85, 'b_Ht': 2.0, 'M2_Ht': 0.8, 'P_Ht': 1.0,
        # E
        'r_E': 0.9, 'n_E': 1.25, 'alpha0_E': 0.43, 'alpha1_E': 0.85, 'b_E': 0.4, 'M2_E': 0.64, 'P_E': 1.0,
        # Etilde
        'r_Et': 1.0, 'n_Et': 0.6, 'alpha0_Et': 0.0, 'alpha1_Et': 0.0, 'b_Et': 0.0, 'M2_Et': 0.0, 'P_Et': 0.0,
    }


# -----------------------------------------------------------------------------
# Push flags and parameters into the already-loaded DVCS_xsec.C (no re-includes)
# -----------------------------------------------------------------------------
def push_globals(param_map, flags, *, strategy1_like=True):
    """
    Set global flags/parameters inside DVCS_xsec.C via ROOT.gInterpreter.ProcessLine.
    Only assigns keys present in param_map (so missing keys keep their DVCS defaults).
    """
    # Flags first
    ROOT.gInterpreter.ProcessLine("hasH=0; hasHt=0; hasE=0; hasEt=0;")
    for cff in ("H", "Ht", "E", "Et"):
        ROOT.gInterpreter.ProcessLine(f"has{cff} = {int(flags.get(cff, 0))};")

    # Renormalizations (BSA uses only renormImag; renormReal used in Re-parts)
    renI = float(param_map.get('renormImag', 1.0))
    ROOT.gInterpreter.ProcessLine(f"renormImag = {renI};")

    # If you want to mimic "strategy 1" semantics (Im-only) you can set renormReal=1.0 safely;
    # BSA does not depend on Re-parts in these formulas (it uses Im pieces in s1_I).
    renR = float(param_map.get('renormReal', 1.0))
    ROOT.gInterpreter.ProcessLine(f"renormReal = {renR};")

    # Helper to set one variable if present
    def set_if_present(name):
        if name in param_map:
            ROOT.gInterpreter.ProcessLine(f"{name} = {float(param_map[name])};")

    # Imaginary-ansatz parameters per CFF
    for cff in ("H", "Ht", "E", "Et"):
        if not int(flags.get(cff, 0)):
            continue
        for k in ("r", "n", "alpha0", "alpha1", "b", "M2", "P"):
            set_if_present(f"{k}_{cff}")

    # Optional real-part subtraction constants (used only if present)
    for cff in ("H", "Ht", "E", "Et"):
        for k in ("C0", "MD2", "lambda"):
            set_if_present(f"{k}_{cff}")


# -----------------------------------------------------------------------------
# Compute BSA curve for arrays of kinematics
# -----------------------------------------------------------------------------
def compute_bsa(phi_arr, Q2_arr, xB_arr, t_arr, Eb_arr, *, debug_tag=None):
    """
    Build BMK_DVCS and return dvcs.BSA() for each point.
    Assumes DVCS_xsec.C is already loaded and globals are set.
    """
    out = np.empty_like(phi_arr, dtype=float)
    for i, (phi, Q2, xB, t, Eb) in enumerate(zip(phi_arr, Q2_arr, xB_arr, t_arr, Eb_arr)):
        dvcs = ROOT.BMK_DVCS(-1, 1, 0, Eb, xB, Q2, t, phi)  # constructor forces t<=0 internally
        mA = dvcs.BSA()
        out[i] = mA
        if debug_tag and i < 3:
            # quick peek at a few points
            ROOT.gInterpreter.ProcessLine("/* debug */")
            print(f"[{debug_tag}] φ={phi:6.1f}°, Q2={Q2:.3f}, xB={xB:.3f}, -t={abs(t):.3f} -> BSA={mA:.6f}")
    return out


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('fitfile', help='Fit results file (output/fit_results/fit_results_<TIMESTAMP>.txt)')
    parser.add_argument('--CFFs', type=int, choices=[0, 1], default=0,
                        help='0: show only ImH fit (or whatever flags were used), '
                             '1: compare ImH-only vs all enabled CFFs from the fit')
    parser.add_argument('--data', default=None,
                        help='Override BSA input data file (otherwise taken from "input" line in fit file).')
    args = parser.parse_args()

    # Parse fit file
    flags, pnames, pvals, input_path, strategy = parse_fit_results(args.fitfile)
    param_map_fit = dict(zip(pnames, pvals))

    # If no input path in fit file, allow override or fall back to your usual file
    if args.data is not None:
        datafile = args.data
    elif input_path is not None and os.path.exists(input_path):
        datafile = input_path
    else:
        # fallback (adjust if your repo uses a different path)
        datafile = 'imports/rgk_preliminary_bsa.txt'

    print(">> Flags from fit:", flags)
    print(">> Using data file:", datafile)
    print(">> Fitted parameters:", param_map_fit)

    # Load DVCS_xsec.C once
    _load_dvcs_once()

    # Load data bins
    bins = load_all_bins(datafile)
    print(f">> Found {len(bins)} φ-bins")

    os.makedirs('output/plots', exist_ok=True)

    # Prepare an "original model" parameter map using DVCS defaults
    param_map_orig = original_model_defaults()

    # Plot per φ-bin
    for idx, b in enumerate(bins, start=1):
        phi_data, As, sigAs = b["phi"], b["A"], b["sigA"]

        # Grid to draw smooth model curves
        phi_grid = np.linspace(0.0, 360.0, 361)
        Q2g = np.full_like(phi_grid, b["Q2m"])
        xBg = np.full_like(phi_grid, b["xBm"])
        tg  = np.full_like(phi_grid, b["tm"])     # note: DVCS_xsec.C enforces t<=0 internally
        Ebg = np.full_like(phi_grid, b["Ebm"])

        # --- Original model curve ---
        push_globals(param_map_orig, flags, strategy1_like=True)
        bsas_orig = compute_bsa(phi_grid, Q2g, xBg, tg, Ebg, debug_tag=f"bin{idx}-orig")

        # --- Fitted model curves ---
        # Full "as fit" configuration
        push_globals(param_map_fit, flags, strategy1_like=(strategy == 1 or strategy is None))
        bsas_fit_all = compute_bsa(phi_grid, Q2g, xBg, tg, Ebg, debug_tag=f"bin{idx}-fitAll")

        # Optionally, ImH-only comparison: turn off Ht/E/Et regardless of fit flags
        bsas_fit_imh = None
        if args.CFFs == 1:
            flags_imh_only = dict(flags)
            for other in ("Ht", "E", "Et"):
                flags_imh_only[other] = 0
            push_globals(param_map_fit, flags_imh_only, strategy1_like=(strategy == 1 or strategy is None))
            bsas_fit_imh = compute_bsa(phi_grid, Q2g, xBg, tg, Ebg, debug_tag=f"bin{idx}-fitImH")

        # --- Plot ---
        fig, ax = plt.subplots(figsize=(8, 5))

        # Data with errors
        ax.errorbar(phi_data, As, yerr=sigAs, fmt='o', ms=5, color='k', label='Data')

        # Original & fitted
        ax.plot(phi_grid, bsas_orig, '-', lw=2, color='tab:blue', label='Original model')
        ax.plot(phi_grid, bsas_fit_all, '--', lw=2, color='tab:red',
                label=('Fitted model (enabled CFFs)' if args.CFFs == 1 else 'Fitted model'))

        if bsas_fit_imh is not None:
            ax.plot(phi_grid, bsas_fit_imh, '-.', lw=2, color='tab:green', label='Fitted model (ImH only)')

        ax.set_xlim(0, 360)
        ax.set_xticks([0, 60, 120, 180, 240, 300, 360])
        ax.set_ylim(-0.6, 0.6)
        ax.set_xlabel(r'$\phi\;[\mathrm{deg}]$')
        ax.set_ylabel(r'$A_{LU}(\phi)$')

        Q2m, xBm, tm = b["Q2m"], b["xBm"], b["tm"]
        ax.set_title(
            (r'$\langle Q^2\rangle={:.2f}\,\mathrm{{GeV}}^2,\;'
             r'\langle x_B\rangle={:.3f},\;\langle -t\rangle={:.3f}\,\mathrm{{GeV}}^2$').format(
                Q2m, xBm, abs(tm)
            ),
            pad=12
        )

        ax.legend(loc='upper right', frameon=True, edgecolor='k')
        plt.tight_layout()

        # Save
        # Use a timestamp extracted from the fitfile name if present
        m = re.search(r'fit_results_(\d{8}_\d{6})\.txt$', os.path.basename(args.fitfile))
        ts = (m.group(1) if m else "noTS")
        fname = (f'output/plots/BSA_bin{idx:02d}_{ts}_Q2_{Q2m:.2f}_xB_{xBm:.3f}_mt_{abs(tm):.3f}.pdf')
        fig.savefig(fname)
        print(f">> Saved bin {idx} plot to {fname}")
        plt.close(fig)


if __name__ == '__main__':
    main()