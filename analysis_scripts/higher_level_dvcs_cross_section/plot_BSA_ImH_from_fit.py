#!/usr/bin/env python3
"""
plot_BSA_ImH_from_fit.py

Usage:
    python plot_BSA_ImH_from_fit.py output/fit_results/fit_results_<TIMESTAMP>.txt
"""

import os, sys, re
import numpy as np
import matplotlib.pyplot as plt
import ROOT

def parse_fit_results(fname):
    """Extract the eight fit values from your text file."""
    with open(fname) as f:
        lines = [l.strip() for l in f if l.strip()]
    for i, L in enumerate(lines):
        if L.startswith("# values"):
            return list(map(float, lines[i+1].split()))
    raise RuntimeError("Couldn't find '# values' in fit file")

def load_first_bin(datafile):
    """
    Read imports/rga_prl_bsa.txt until phi wraps around.
    Returns arrays of φ, Q2, xB, t, Eb, A, sigA for the first bin.
    """
    phis, Q2s, xBs, ts, Ebs, As, sigAs = [], [], [], [], [], [], []
    prev_phi = -1.0
    with open(datafile) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            phi, Q2, xB, t, Eb, A, sigA = map(float, line.split())
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
    return (np.array(phis),
            np.array(Q2s),
            np.array(xBs),
            np.array(ts),
            np.array(Ebs),
            np.array(As),
            np.array(sigAs))

def compute_bsa(phi_arr, Q2_arr, xB_arr, t_arr, Eb_arr, params, tag=""):
    """
    For each i, build BMK_DVCS(-1,1,0,Eb,xB,Q2,t,phi) and return dvcs.BSA().
    Prints the first five for quick debug.
    """
    (rI, a0, a1, nv, bv, m2, Pv, rR) = params
    # assign C globals
    ROOT.renormImag = rI
    ROOT.alpha0     = a0
    ROOT.alpha1     = a1
    ROOT.n_val      = nv
    ROOT.b_val      = bv
    ROOT.Mm2_val    = m2
    ROOT.P_val      = Pv
    ROOT.renormReal = rR
    # switch on H only
    ROOT.hasH  = True
    ROOT.hasHt = False
    ROOT.hasE  = False
    ROOT.hasEt = False

    bsas = []
    for i, (phi, Q2, xB, t, Eb) in enumerate(zip(
            phi_arr, Q2_arr, xB_arr, t_arr, Eb_arr)):
        dvcs = ROOT.BMK_DVCS(-1, 1, 0, Eb, xB, Q2, t, phi)
        mA   = dvcs.BSA()
        bsas.append(mA)
        if i < 5:
            print(f"[{tag}] φ={phi:6.1f}°, ξ={dvcs.xi:.3f}, t={t:.3f}, BSA={mA:.4f}")
    return np.array(bsas)

def main():
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)

    fitfile = sys.argv[1]
    m = re.search(r'fit_results_(\d{8}_\d{6})\.txt$', fitfile)
    if not m:
        print("ERROR: can't extract timestamp from", fitfile)
        sys.exit(1)
    timestamp = m.group(1)

    # 1) parse fit
    params_fit = parse_fit_results(fitfile)
    print(">> Fitted parameters:", params_fit)

    # 2) load C++ DVCS library
    ROOT.gSystem.Load('./DVCS_xsec_C.so')

    # 3) read first φ‐bin
    datafile = 'imports/rga_prl_bsa.txt'
    phi_data, Q2s, xBs, ts, Ebs, As, sigAs = load_first_bin(datafile)
    print(f">> Loaded {len(phi_data)} points in first φ‐bin")

    # 4) build 100‐point φ grid
    phi_grid = np.linspace(0, 360, 100)
    Q2m, xBm, tm, Ebm = Q2s.mean(), xBs.mean(), ts.mean(), Ebs.mean()
    Q2g = np.full_like(phi_grid, Q2m)
    xBg = np.full_like(phi_grid, xBm)
    tg  = np.full_like(phi_grid, tm)
    Ebg = np.full_like(phi_grid, Ebm)

    # 5) compute BSA on grid
    defaults       = [1.0, 0.43, 0.85, 1.35, 0.4, 0.64, 1.0, 1.0]
    bsas_orig_grid = compute_bsa(phi_grid, Q2g, xBg, tg, Ebg,
                                 defaults, tag="orig-grid")
    bsas_fit_grid  = compute_bsa(phi_grid, Q2g, xBg, tg, Ebg,
                                 params_fit, tag="fit-grid")

    # 6) plot
    fig, ax = plt.subplots(figsize=(8, 5))
    # data
    ax.errorbar(phi_data, As, yerr=sigAs, fmt='o', ms=5,
                color='k', label='Data')
    # model lines
    ax.plot(phi_grid, bsas_orig_grid, '-',  lw=2,
            color='tab:blue', label='Original Model')
    ax.plot(phi_grid, bsas_fit_grid,  '--', lw=2,
            color='tab:red',  label='Fitted Model')

    # formatting
    ax.set_xlim(0, 360)
    ax.set_xticks([0, 60, 120, 180, 240, 300, 360])
    ax.set_ylim(-0.6, 0.6)
    ax.set_xlabel(r'$\phi\;[\mathrm{deg}]$')
    ax.set_ylabel(r'$A_{LU}(\phi)$')

    ax.set_title(
        (r'$\langle Q^2\rangle={:.2f}\,\mathrm{{GeV}}^2,\;\langle x_B\rangle={:.3f},\;\langle -t\rangle={:.3f}\,\mathrm{{GeV}}^2$'
        ).format(Q2m, xBm, -tm),
        pad=12
    )

    ax.legend(loc='upper right', frameon=True, edgecolor='k')
    plt.tight_layout()

    outdir = 'output/plots'
    os.makedirs(outdir, exist_ok=True)
    outname = f'{outdir}/BSA_bin1_{timestamp}.pdf'
    fig.savefig(outname)
    print(">> Saved", outname)

if __name__ == '__main__':
    main()