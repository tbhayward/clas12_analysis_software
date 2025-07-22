#!/usr/bin/env python3
"""
plot_BSA_ImH_from_fit.py

Usage:
    python plot_BSA_ImH_from_fit.py output/fit_results/fit_results_<TIMESTAMP>.txt
"""
import os
import sys
import re
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

def load_all_bins(datafile):
    """
    Read imports/rga_prl_bsa.txt, split into bins when φ wraps around.
    Returns a list of dicts with arrays and their means.
    """
    bins = []
    curr = {k: [] for k in ("phi","Q2","xB","t","Eb","A","sigA")}
    prev_phi = None

    with open(datafile) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            phi, Q2, xB, t, Eb, A, sigA = map(float, line.split())
            # new bin if phi decreases
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

    # last bin
    if curr["phi"]:
        arr = {k: np.array(v) for k, v in curr.items()}
        arr["Q2m"] = arr["Q2"].mean()
        arr["xBm"] = arr["xB"].mean()
        arr["tm"]  = arr["t"].mean()
        arr["Ebm"] = arr["Eb"].mean()
        bins.append(arr)

    return bins

def compute_bsa(phi_arr, Q2_arr, xB_arr, t_arr, Eb_arr, params, tag=""):
    """
    Build BMK_DVCS(-1,1,0,Eb,xB,Q2,t,phi) for each point and return dvcs.BSA().
    Prints the first five for debug.
    """
    (rI, a0, a1, nv, bv, m2, Pv, rR) = params

    # set C++ globals
    ROOT.renormImag = rI
    ROOT.alpha0     = a0
    ROOT.alpha1     = a1
    ROOT.n_val      = nv
    ROOT.b_val      = bv
    ROOT.Mm2_val    = m2
    ROOT.P_val      = Pv
    ROOT.renormReal = rR

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

    # 1) parse fit parameters
    params_fit = parse_fit_results(fitfile)
    print(">> Fitted parameters:", params_fit)

    # 2) load C++ DVCS library
    ROOT.gSystem.Load('./DVCS_xsec_C.so')

    # 3) load all φ-bins
    datafile = 'imports/rga_prl_bsa.txt'
    bins = load_all_bins(datafile)
    print(f">> Found {len(bins)} φ-bins")

    # ensure output directory
    outdir = 'output/plots'
    os.makedirs(outdir, exist_ok=True)

    # 4) loop through each bin and make one canvas
    defaults = [1.0, 0.43, 0.85, 1.35, 0.4, 0.64, 1.0, 1.0]
    for idx, b in enumerate(bins, start=1):
        phi_data = b["phi"]
        As, sigAs = b["A"], b["sigA"]

        # build 100-point φ grid
        phi_grid = np.linspace(0, 360, 100)
        Q2g = np.full_like(phi_grid, b["Q2m"])
        xBg = np.full_like(phi_grid, b["xBm"])
        tg  = np.full_like(phi_grid, b["tm"])
        Ebg = np.full_like(phi_grid, b["Ebm"])

        # compute models
        bsas_orig = compute_bsa(phi_grid, Q2g, xBg, tg, Ebg,
                                defaults, tag=f"bin{idx}-orig")
        bsas_fit  = compute_bsa(phi_grid, Q2g, xBg, tg, Ebg,
                                params_fit, tag=f"bin{idx}-fit")

        # plot
        fig, ax = plt.subplots(figsize=(8,5))
        ax.errorbar(phi_data, As, yerr=sigAs, fmt='o', ms=5,
                    color='k', label='Data')
        ax.plot(phi_grid, bsas_orig, '-',  lw=2,
                color='tab:blue', label='Original Model')
        ax.plot(phi_grid, bsas_fit,  '--', lw=2,
                color='tab:red',  label='Fitted Model')

        ax.set_xlim(0, 360)
        ax.set_xticks([0, 60, 120, 180, 240, 300, 360])
        ax.set_ylim(-0.6, 0.6)

        ax.set_xlabel(r'$\phi\;[\mathrm{deg}]$')
        ax.set_ylabel(r'$A_{LU}(\phi)$')

        # title with mean Q2, xB, -t
        Q2m, xBm, tm = b["Q2m"], b["xBm"], b["tm"]
        ax.set_title(
            (r'$\langle Q^2\rangle={:.2f}\,\mathrm{{GeV}}^2,\;'
             r'\langle x_B\rangle={:.3f},\;\langle -t\rangle={:.3f}\,\mathrm{{GeV}}^2$'
            ).format(Q2m, xBm, -tm),
            pad=12
        )

        ax.legend(loc='upper right', frameon=True, edgecolor='k')
        plt.tight_layout()

        # filename with bin index, timestamp, Q2, xB, t
        fname = (f'{outdir}/BSA_bin{idx:02d}_'
                 f'{timestamp}_Q2_{Q2m:.2f}_xB_{xBm:.3f}_t_{abs(tm):.3f}.pdf')
        fig.savefig(fname)
        print(f">> Saved bin {idx} plot to {fname}")
        plt.close(fig)

if __name__ == '__main__':
    main()