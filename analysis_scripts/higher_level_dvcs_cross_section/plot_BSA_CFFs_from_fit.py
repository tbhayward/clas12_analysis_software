#!/usr/bin/env python3
"""
plot_BSA_CFFs_from_fit.py

Usage:
    python plot_BSA_CFFs_from_fit.py output/fit_results/fit_results_<TIMESTAMP>.txt --CFFs [0|1] [--override-data PATH] [--debug]
"""
import os
import sys
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt
import ROOT

def parse_fit_results(fname):
    with open(fname) as f:
        lines = [l.rstrip("\n") for l in f if l.strip()]

    # flags: "H 1 Ht 0 E 0 Et 0"
    flag_line = next((l for l in lines if re.match(r'^\s*H\s+\d+', l)), None)
    if flag_line is None:
        raise RuntimeError("Couldn't find flags line in fit file")
    toks = flag_line.split()
    flags = {}
    for i in range(0, len(toks), 2):
        k = toks[i]
        try:
            flags[k] = int(toks[i+1])
        except Exception:
            pass

    pnames, vals, errs = [], None, None
    for i, l in enumerate(lines):
        if l.startswith("# parameters"):
            pnames = l.split()[2:]
        elif l.startswith("# values:"):
            if i+1 < len(lines):
                vals = [float(x) for x in lines[i+1].split()]
        elif l.startswith("# errors:"):
            if i+1 < len(lines):
                errs = [float(x) for x in lines[i+1].split()]

    if not pnames or vals is None:
        raise RuntimeError("Couldn't parse parameters/values")

    input_path = None
    input_line = next((l for l in lines if l.startswith("input")), None)
    if input_line:
        parts = input_line.split()
        if len(parts) >= 2:
            input_path = parts[1]

    return flags, pnames, vals, errs, input_path

def load_all_bins(datafile):
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
                arr["tm"]  = arr["t"].mean()   # keep negative t
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

# Default model values (match DVCS_xsec.C globals)
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

def push_globals(param_map, flags, *, strategy1_like=True):
    """
    Push C++ globals into the ROOT JIT:
      - reset has* to 0, then set from flags
      - set renormImag=1.0; renormReal=0 if strategy-1-like (no Re fit)
      - push Im-ansatz parameters (r, n, alpha0, alpha1, b, M2, P)
    """
    # ensure file is loaded
    ROOT.gInterpreter.ProcessLine('#include "DVCS_xsec.C"')

    # Reset all switches first, then enable according to flags
    ROOT.gInterpreter.ProcessLine("hasH=0; hasHt=0; hasE=0; hasEt=0;")
    for cff in ("H","Ht","E","Et"):
        v = int(flags.get(cff, 0))
        ROOT.gInterpreter.ProcessLine(f"has{cff} = {v};")

    # Renormalizations
    renI = float(param_map.get('renormImag', orig_defaults['renormImag']))
    ROOT.gInterpreter.ProcessLine(f"renormImag = {renI};")
    if strategy1_like:
        ROOT.gInterpreter.ProcessLine("renormReal = 0.0;")  # no Re in BSA-only fits
    else:
        renR = float(param_map.get('renormReal', orig_defaults['renormReal']))
        ROOT.gInterpreter.ProcessLine(f"renormReal = {renR};")

    # Push Im-parameters for any enabled CFF
    for cff in ("H","Ht","E","Et"):
        if not int(flags.get(cff,0)):
            continue
        for k in ("r","n","alpha0","alpha1","b","M2","P"):
            key = f"{k}_{cff}"
            val = param_map.get(key, orig_defaults.get(key))
            if val is None:
                continue
            ROOT.gInterpreter.ProcessLine(f"{key} = {float(val)};")

def compute_bsa(phi_arr, Q2_arr, xB_arr, t_arr, Eb_arr, param_map, flags, *,
                strategy1_like=True, debug=False, tag=""):
    """
    Build BMK_DVCS objects & return BSA(φ). t_arr should be negative.
    """
    push_globals(param_map, flags, strategy1_like=strategy1_like)

    bsas = []
    for i, (phi, Q2, xB, t, Eb) in enumerate(zip(phi_arr, Q2_arr, xB_arr, t_arr, Eb_arr)):
        dvcs = ROOT.BMK_DVCS(-1, 1, 0, Eb, xB, Q2, t, phi)  # will flip helicity inside BSA()
        mA   = dvcs.BSA()
        bsas.append(mA)

        if debug and i < 3:
            # probe s1_I (helicity-odd piece) vs c0_BH just to confirm it's alive
            # set L_beam = +1 and sample the coefficients
            dvcs.L_beam = 1
            s1 = dvcs.s1_I()
            c0 = dvcs.c0_BH()
            try:
                xi = dvcs.xi
            except Exception:
                xi = np.nan
            print(f"[{tag}] φ={phi:6.1f}°, xi={xi:.3f}, t={t:.3f}  BSA={mA:.5f}  s1_I={s1:.3e}  c0_BH={c0:.3e}")

    return np.array(bsas)

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('fitfile', help='fit results file from fit_CFFs.C')
    ap.add_argument('--CFFs', type=int, choices=[0,1], default=0,
                    help='0: show ImH-only fit; 1: also show full flagged CFF fit')
    ap.add_argument('--override-data', dest='override_data', default=None,
                    help='override BSA data file path')
    ap.add_argument('--debug', action='store_true', help='print a few s1_I/c0_BH probes')
    args = ap.parse_args()

    m = re.search(r'fit_results_(\d{8}_\d{6})\.txt$', args.fitfile)
    if not m:
        print("ERROR: cannot extract timestamp from fitfile name"); sys.exit(1)
    timestamp = m.group(1)

    flags, pnames, vals, errs, input_path = parse_fit_results(args.fitfile)
    param_map = dict(zip(pnames, vals))

    datafile = args.override_data if args.override_data else input_path
    if not datafile:
        datafile = 'imports/rga_prl_bsa.txt'
        print("WARNING: no 'input' path in fit file; falling back to", datafile)
    if not os.path.isfile(datafile):
        print("ERROR: data file not found:", datafile); sys.exit(1)

    print(">> Flags:", flags)
    print(">> Using BSA data file:", datafile)

    ROOT.gInterpreter.ProcessLine('#include "DVCS_xsec.C"')

    bins = load_all_bins(datafile)
    print(f">> Found {len(bins)} φ-bins")
    if not bins:
        print("ERROR: no bins parsed from data file"); sys.exit(1)

    os.makedirs('output/plots', exist_ok=True)

    # strategy-1-like if the fit file did not include any Re-part parameters
    strategy1_like = not any(n in param_map for n in ('renormReal','C0_H','MD2_H','lambda_H'))

    for idx, b in enumerate(bins, start=1):
        phi_data, As, sigAs = b["phi"], b["A"], b["sigA"]

        # grid at mean kinematics
        phi_grid = np.linspace(0, 360, 181)
        Q2g = np.full_like(phi_grid, b["Q2m"])
        xBg = np.full_like(phi_grid, b["xBm"])
        tg  = np.full_like(phi_grid, b["tm"])   # negative
        Ebg = np.full_like(phi_grid, b["Ebm"])

        # default model, using the same active CFFs as the fit flags
        bsas_orig = compute_bsa(phi_grid, Q2g, xBg, tg, Ebg,
                                orig_defaults, flags,
                                strategy1_like=strategy1_like,
                                debug=args.debug, tag=f"bin{idx}-orig")

        # fitted model: full flagged set
        bsas_fit_all = compute_bsa(phi_grid, Q2g, xBg, tg, Ebg,
                                   param_map, flags,
                                   strategy1_like=strategy1_like,
                                   debug=args.debug, tag=f"bin{idx}-fit-all")

        fig, ax = plt.subplots(figsize=(8,5))
        ax.errorbar(phi_data, As, yerr=sigAs, fmt='o', ms=5, color='k', label='Data')
        ax.plot(phi_grid, bsas_orig, '-',  lw=2, color='tab:blue',  label='Default model')

        if args.CFFs == 0:
            # ImH-only curve (others off)
            flags_imh = flags.copy()
            for o in ("Ht","E","Et"):
                flags_imh[o] = 0
            bsas_fit_imh = compute_bsa(phi_grid, Q2g, xBg, tg, Ebg,
                                       param_map, flags_imh,
                                       strategy1_like=strategy1_like,
                                       debug=args.debug, tag=f"bin{idx}-fit-ImH")
            ax.plot(phi_grid, bsas_fit_imh, '--', lw=2, color='tab:red', label='Fitted (ImH)')
        else:
            # show both ImH-only and full flagged
            flags_imh = flags.copy()
            for o in ("Ht","E","Et"):
                flags_imh[o] = 0
            bsas_fit_imh = compute_bsa(phi_grid, Q2g, xBg, tg, Ebg,
                                       param_map, flags_imh,
                                       strategy1_like=strategy1_like,
                                       debug=args.debug, tag=f"bin{idx}-fit-ImH")
            ax.plot(phi_grid, bsas_fit_imh, '--', lw=2, color='tab:red',   label='Fitted (ImH)')
            ax.plot(phi_grid, bsas_fit_all, '-.', lw=2, color='tab:green', label='Fitted (all flagged)')

        ax.set_xlim(0, 360)
        ax.set_xticks([0,60,120,180,240,300,360])
        ax.set_ylim(-0.6, 0.6)
        ax.set_xlabel(r'$\phi\;[\mathrm{deg}]$')
        ax.set_ylabel(r'$A_{LU}(\phi)$')

        Q2m, xBm, tm = b["Q2m"], b["xBm"], b["tm"]
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