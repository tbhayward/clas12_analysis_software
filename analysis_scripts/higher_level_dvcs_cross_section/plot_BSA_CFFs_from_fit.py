#!/usr/bin/env python3
"""
plot_BSA_CFFs_from_fit.py

Usage:
    python plot_BSA_CFFs_from_fit.py output/fit_results/fit_results_<TIMESTAMP>.txt [--CFFs 0|1] [--data PATH] [--debug]

Notes
-----
• The "φ=0°,1°,2°" lines are just a sanity print of the first 3 points on a dense
  φ-grid (0..360° step 1°). The full curve is computed.

• If BSA=0 everywhere, the usual causes are:
  (a) s1_I=0 (e.g. K=0 near t≈t_min, or all Im CFFs effectively zero/disabled),
  (b) renorms/flags/params not applied into the loaded C++,
  (c) something pathological in kinematics (y=0, etc).

• Pass --debug to print internals (xi, t_min, K, BH2, DVCS2, BHDVCS, s1_I, GetImH)
  for a few φ points of the first bin.
"""

import os
import sys
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt
import ROOT  # PyROOT

# -----------------------------------------------------------------------------
# One-time, guarded load of DVCS_xsec.C
# -----------------------------------------------------------------------------
def _load_dvcs_once():
    """Load DVCS_xsec.C exactly once into cling without redefinitions."""
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
# Parse the fit-results text file produced by the fitter
# -----------------------------------------------------------------------------
def parse_fit_results(fname):
    """
    Returns:
      flags: dict {'H':1,'Ht':0,'E':0,'Et':0}
      pnames: list of parameter names
      pvals:  list of corresponding values
      input_file: path in line "input <path>" (if present), else None
      strategy: int if present, else None
    """
    with open(fname) as f:
        lines = [l.strip() for l in f if l.strip()]

    # flags
    flag_line = None
    for l in lines:
        if re.match(r'^\s*H\s+[01]\s+Ht\s+[01]\s+E\s+[01]\s+Et\s+[01]\s*$', l):
            flag_line = l
            break
    if flag_line is None:
        flag_line = next((l for l in lines if l.startswith("H ")), None)
    if flag_line is None:
        raise RuntimeError("Could not find 'H ... Ht ... E ... Et ...' flags line in fit file")

    toks = flag_line.split()
    flags = {}
    for i in range(0, len(toks), 2):
        key = toks[i]
        val = int(toks[i + 1])
        flags[key] = val

    # parameters & values
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

    # optional input file path
    input_file = None
    for l in lines:
        if l.startswith("input"):
            parts = l.split(None, 1)
            if len(parts) == 2:
                input_file = parts[1].strip()
            break

    # optional strategy
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
# Read a BSA data file and split into φ-bins (wrap when φ decreases)
# -----------------------------------------------------------------------------
def load_all_bins(datafile):
    """
    Expect columns: phi  Q2  xB  t  Eb  A  sigA
    Splits into bins when phi wraps around (decreases).
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
# "Original" model defaults (match DVCS_xsec.C variable names)
# -----------------------------------------------------------------------------
def original_model_defaults():
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
# Push flags and parameters into the already-loaded DVCS_xsec.C
# -----------------------------------------------------------------------------
def push_globals(param_map, flags):
    """
    Sets the global switches and parameters in the C++:
      hasH, hasHt, hasE, hasEt; renormImag/Real; and any of
      r_*, n_*, alpha0_*, alpha1_*, b_*, M2_*, P_* present in param_map.
    """
    # switches
    ROOT.gInterpreter.ProcessLine("hasH=0; hasHt=0; hasE=0; hasEt=0;")
    for cff in ("H", "Ht", "E", "Et"):
        ROOT.gInterpreter.ProcessLine(f"has{cff} = {int(flags.get(cff, 0))};")

    # renormalizations
    renI = float(param_map.get('renormImag', 1.0))
    renR = float(param_map.get('renormReal', 1.0))
    ROOT.gInterpreter.ProcessLine(f"renormImag = {renI};")
    ROOT.gInterpreter.ProcessLine(f"renormReal = {renR};")

    # helper
    def set_if_present(name):
        if name in param_map:
            ROOT.gInterpreter.ProcessLine(f"{name} = {float(param_map[name])};")

    # imaginary-ansatz parameters
    for cff in ("H", "Ht", "E", "Et"):
        if not int(flags.get(cff, 0)):
            continue
        for k in ("r", "n", "alpha0", "alpha1", "b", "M2", "P"):
            set_if_present(f"{k}_{cff}")

    # optional real-part subtraction constants (only if present)
    for cff in ("H", "Ht", "E", "Et"):
        for k in ("C0", "MD2", "lambda"):
            set_if_present(f"{k}_{cff}")


# -----------------------------------------------------------------------------
# Compute BSA curve; optional deep debug print for first few points
# -----------------------------------------------------------------------------
def compute_bsa(phi_arr, Q2_arr, xB_arr, t_arr, Eb_arr, *, debug=False, debug_tag=None, max_debug_points=3):
    out = np.empty_like(phi_arr, dtype=float)
    for i, (phi, Q2, xB, t, Eb) in enumerate(zip(phi_arr, Q2_arr, xB_arr, t_arr, Eb_arr)):
        dvcs = ROOT.BMK_DVCS(-1, 1, 0, Eb, xB, Q2, t, phi)
        mA = dvcs.BSA()
        out[i] = mA

        if debug and i < max_debug_points:
            # Internals that matter for ALU
            xi = dvcs.xi
            tmin = dvcs.t_min
            K = dvcs.K
            bh2 = dvcs.BH2()
            dv2 = dvcs.DVCS2()
            inter = dvcs.BHDVCS()
            s1i = dvcs.s1_I()  # should carry L_beam
            # evaluate ImH at (xi,t) to check it's not zero
            ROOT.gInterpreter.ProcessLine("double __tmpImH = GetImH(%g, %g);" % (xi, dvcs.t))
            # fetch it back by printing (cling easiest way)
            ROOT.gInterpreter.ProcessLine('std::cout.setf(std::ios::scientific);')
            ROOT.gInterpreter.ProcessLine('std::cout<<"[samp] ImH="<<__tmpImH<<std::endl;')
            print(f"[{debug_tag}] φ={phi:6.1f}°, Q2={Q2:.3f}, xB={xB:.3f}, -t={abs(dvcs.t):.3f} | "
                  f"xi={xi:.4f}, t_min={tmin:.4f}, K={K:.6e}, "
                  f"BH2={bh2:.6e}, DVCS2={dv2:.6e}, I={inter:.6e}, s1_I={s1i:.6e} | BSA={mA:.6e}")

    return out


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('fitfile', help='output/fit_results/fit_results_<TIMESTAMP>.txt')
    parser.add_argument('--CFFs', type=int, choices=[0, 1], default=0,
                        help='0: show only fitted curve using enabled CFFs; '
                             '1: also show ImH-only fitted curve')
    parser.add_argument('--data', default=None,
                        help='Override BSA input file (else use "input" line from fit file or fallback).')
    parser.add_argument('--debug', action='store_true', help='Print internal diagnostics for first bin.')
    args = parser.parse_args()

    # Parse fit file
    flags, pnames, pvals, input_path, strategy = parse_fit_results(args.fitfile)
    param_map_fit = dict(zip(pnames, pvals))

    # Decide data file
    if args.data:
        datafile = args.data
    elif input_path and os.path.exists(input_path):
        datafile = input_path
    else:
        datafile = 'imports/rgk_preliminary_bsa.txt'

    print(">> Flags from fit:", flags)
    print(">> Using data file:", datafile)
    print(">> Fitted parameters:", param_map_fit)

    # Load C++ once
    _load_dvcs_once()

    # Load data bins
    bins = load_all_bins(datafile)
    print(f">> Found {len(bins)} φ-bins")

    os.makedirs('output/plots', exist_ok=True)

    # Prepare defaults for the "original model" curve
    param_map_orig = original_model_defaults()

    # Optional: quick dump of C++ state after pushing globals
    def dump_cpp_state(label):
        ROOT.gInterpreter.ProcessLine(rf'std::cout<<"[{label}] hasH="<<hasH<<" hasHt="<<hasHt<<" hasE="<<hasE<<" hasEt="<<hasEt'
                                      r'<<" renormImag="<<renormImag<<" renormReal="<<renormReal<<std::endl;')

    # Plot per bin
    for idx, b in enumerate(bins, start=1):
        phi_data, As, sigAs = b["phi"], b["A"], b["sigA"]

        # Dense grid for model curves
        phi_grid = np.linspace(0.0, 360.0, 361)  # 0..360 step 1°
        Q2g = np.full_like(phi_grid, b["Q2m"])
        xBg = np.full_like(phi_grid, b["xBm"])
        tg  = np.full_like(phi_grid, b["tm"])   # constructor enforces t<=0
        Ebg = np.full_like(phi_grid, b["Ebm"])

        # --- Original model ---
        push_globals(param_map_orig, flags)
        if args.debug and idx == 1:
            dump_cpp_state("orig")
        bsas_orig = compute_bsa(phi_grid, Q2g, xBg, tg, Ebg,
                                debug=(args.debug and idx == 1),
                                debug_tag=f"bin{idx}-orig")

        # --- Fitted model (enabled CFFs per fit flags) ---
        push_globals(param_map_fit, flags)
        if args.debug and idx == 1:
            dump_cpp_state("fitAll")
        bsas_fit_all = compute_bsa(phi_grid, Q2g, xBg, tg, Ebg,
                                   debug=(args.debug and idx == 1),
                                   debug_tag=f"bin{idx}-fitAll")

        # Optional ImH-only overlay
        bsas_fit_imh = None
        if args.CFFs == 1:
            flags_imh_only = dict(flags)
            for other in ("Ht", "E", "Et"):
                flags_imh_only[other] = 0
            push_globals(param_map_fit, flags_imh_only)
            if args.debug and idx == 1:
                dump_cpp_state("fitImH")
            bsas_fit_imh = compute_bsa(phi_grid, Q2g, xBg, tg, Ebg,
                                       debug=False,  # already printed above
                                       debug_tag=f"bin{idx}-fitImH")

        # --- Plot ---
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.errorbar(phi_data, As, yerr=sigAs, fmt='o', ms=5, color='k', label='Data')
        ax.plot(phi_grid, bsas_orig, '-',  lw=2, color='tab:blue', label='Original model')
        ax.plot(phi_grid, bsas_fit_all, '--', lw=2, color='tab:red',
                label=('Fitted (enabled CFFs)' if args.CFFs == 1 else 'Fitted'))

        if bsas_fit_imh is not None:
            ax.plot(phi_grid, bsas_fit_imh, '-.', lw=2, color='tab:green', label='Fitted (ImH only)')

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
        m = re.search(r'fit_results_(\d{8}_\d{6})\.txt$', os.path.basename(args.fitfile))
        ts = (m.group(1) if m else "noTS")
        fname = (f'output/plots/BSA_bin{idx:02d}_{ts}_Q2_{Q2m:.2f}_xB_{xBm:.3f}_mt_{abs(tm):.3f}.pdf')
        fig.savefig(fname)
        print(f">> Saved bin {idx} plot to {fname}")
        plt.close(fig)

        # For debugging speed, you can break after first bin by uncommenting:
        # if args.debug: break


if __name__ == '__main__':
    main()