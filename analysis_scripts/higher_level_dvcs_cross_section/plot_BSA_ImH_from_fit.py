#!/usr/bin/env python3
"""
plot_BSA_ImH_from_fit.py

Usage:
    python plot_BSA_ImH_from_fit.py output/fit_results/fit_results_<TIMESTAMP>.txt

Reads the first φ‐bin from imports/rga_prl_bsa.txt, computes ⟨Q²⟩, ⟨x_B⟩, ⟨−t⟩,
and plots:
  • Data (black circles)
  • Original DVCS_xsec.C defaults (solid blue)
  • Fitted parameters (dashed red)

Saves to output/plots/BSA_bin1_<TIMESTAMP>.pdf
"""

import os
import sys
import re
import numpy as np
import matplotlib.pyplot as plt
import ROOT

def parse_fit_results(fname):
    """Pull the 8 fit parameters (values only) from your text file."""
    with open(fname) as f:
        lines = [l.strip() for l in f if l.strip()]
    vals = None
    for i, L in enumerate(lines):
        if L.startswith("# values"):
            vals = list(map(float, lines[i+1].split()))
            break
    if vals is None:
        raise RuntimeError("Couldn't find '# values' in fit file")
    return vals  # [renormImag, alpha0, alpha1, n_val, b_val, Mm2_val, P_val, renormReal]

def load_first_bin(datafile):
    """Read imports/rga_prl_bsa.txt and return only the first φ-bin."""
    phis = []; Q2s = []; xBs = []; ts = []; Ebs = []; As = []; sigAs = []
    prev_phi = -1.0
    with open(datafile) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"): 
                continue
            phi, Q2, xB, t, Eb, A, sigA = map(float, line.split())
            # φ has wrapped → new bin
            if phis and phi < prev_phi:
                break
            phis.append(phi)
            Q2s.append(Q2)
            xBs.append(xB)
            ts.append(t)
            Ebs.append(Eb)
            As.append(A)
            sigAs.append(sigA)
            prev_phi = phi
    return (np.array(phis), np.array(Q2s), np.array(xBs),
            np.array(ts),  np.array(Ebs),
            np.array(As),  np.array(sigAs))

def compute_bsa(phi_arr, Q2_arr, xB_arr, t_arr, Eb_arr, params):
    """
    Given arrays of kinematics and an 8‐element params list,
    assign to C++ globals and call BMK_DVCS(...).BSA() for each φ.
    """
    # unpack and assign
    ROOT.renormImag = params[0]
    ROOT.alpha0     = params[1]
    ROOT.alpha1     = params[2]
    ROOT.n_val      = params[3]
    ROOT.b_val      = params[4]
    ROOT.Mm2_val    = params[5]
    ROOT.P_val      = params[6]
    ROOT.renormReal = params[7]

    out = []
    for phi, Q2, xB, t, Eb in zip(phi_arr, Q2_arr, xB_arr, t_arr, Eb_arr):
        dvcs = ROOT.BMK_DVCS(-1, 1, 0, Eb, xB, Q2, t, phi)
        out.append(dvcs.BSA())
    return np.array(out)

def main():
    # ─── Command‐line ──────────────────────────────────────────────────────────
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)
    fitfile = sys.argv[1]
    m = re.search(r'fit_results_(\d{8}_\d{6})\.txt$', fitfile)
    if not m:
        print("Couldn't extract timestamp from filename:", fitfile)
        sys.exit(1)
    timestamp = m.group(1)

    # ─── Parse fit parameters ─────────────────────────────────────────────────
    params_fit = parse_fit_results(fitfile)

    # ─── Load & compile DVCS_xsec.C via its shared object ─────────────────────
    # (Adjust path if needed; this assumes you've built DVCS_xsec_C.so in cwd)
    ROOT.gSystem.Load('./DVCS_xsec_C.so')

    # ─── Load the first φ‐bin from data ────────────────────────────────────────
    datafile = 'imports/rga_prl_bsa.txt'
    phi, Q2s, xBs, ts, Ebs, As, sigAs = load_first_bin(datafile)

    # ─── Compute three curves: data, original, fit ────────────────────────────
    # original defaults  = [1.0,0.43,0.85,1.35,0.4,0.64,1.0,1.0]
    defaults = [1.0, 0.43, 0.85, 1.35, 0.4, 0.64, 1.0, 1.0]
    bsas_orig = compute_bsa(phi, Q2s, xBs, ts, Ebs, defaults)
    bsas_fit  = compute_bsa(phi, Q2s, xBs, ts, Ebs, params_fit)

    # ─── Average kinematics for the title ─────────────────────────────────────
    Q2m = Q2s.mean()
    xBm = xBs.mean()
    tm  = ts.mean()

    # ─── Plot ─────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(8,5))
    ax.errorbar(phi, As,  yerr=sigAs, fmt='o', ms=5, color='k', label='Data')
    ax.plot(   phi, bsas_orig, '-',  lw=1.5, color='tab:blue',
               label='Original Parameters')
    ax.plot(   phi, bsas_fit,  '--', lw=1.5, color='tab:red',
               label='RGA pass-1 BSA')

    ax.set_xlim(0, 360)
    ax.set_xticks([0, 60, 120, 180, 240, 300, 360])
    ax.set_ylim(-0.6, 0.6)
    ax.set_xlabel(r'$\phi\;[\mathrm{deg}]$')
    ax.set_ylabel(r'$A_{LU}(\phi)$')

    # escaped braces around \mathrm{GeV}
    ax.set_title(
        (r'$\langle Q^2\rangle={:.2f}\,\mathrm{{GeV}}^2,\;'
         r'\langle x_B\rangle={:.3f},\;\langle -t\rangle={:.3f}\,\mathrm{{GeV}}^2$'
        ).format(Q2m, xBm, -tm),
        pad=12
    )

    ax.legend(loc='upper right', frameon=False)

    # ─── Save ─────────────────────────────────────────────────────────────────
    outdir = 'output/plots'
    os.makedirs(outdir, exist_ok=True)
    outfile = f'{outdir}/BSA_bin1_{timestamp}.pdf'
    fig.savefig(outfile)
    print("Saved figure to", outfile)

if __name__ == '__main__':
    main()