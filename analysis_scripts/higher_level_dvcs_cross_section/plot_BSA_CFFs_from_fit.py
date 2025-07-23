#!/usr/bin/env python3
"""
plot_BSA_CFFs_from_fits.py

Usage:
    python plot_BSA_CFFs_from_fits.py output/fit_results/fit_results_<TIMESTAMP>.txt
"""
import os
import sys
import re
import numpy as np
import matplotlib.pyplot as plt
import ROOT

def parse_fit_results(fname):
    """
    Extract the fit flags, parameter names, and values from your text file.
    """
    with open(fname) as f:
        lines = [l.strip() for l in f if l.strip()]
    # flags line, e.g. "H 1  Ht 1  E 0  Et 1"
    flag_line = next(l for l in lines if l.startswith("H "))
    toks = flag_line.split()
    flags = { toks[i]: int(toks[i+1]) for i in range(0, len(toks), 2) }
    # parameter names
    pnames = []
    for i, l in enumerate(lines):
        if l.startswith("# parameters"):
            pnames = lines[i].split()[2:]
            break
    # values
    vals = None
    for i, l in enumerate(lines):
        if l.startswith("# values"):
            vals = list(map(float, lines[i+1].split()))
            break
    if vals is None or not pnames:
        raise RuntimeError("Couldn't parse fit file")
    return flags, pnames, vals

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
            if prev_phi is not None and phi < prev_phi:
                arr = {k: np.array(v) for k, v in curr.items()}
                arr["Q2m"] = arr["Q2"].mean()
                arr["xBm"] = arr["xB"].mean()
                arr["tm"]  = arr["t"].mean()
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

# default ("original") parameter values:
orig_defaults = {
    'renormImag': 1.0,
    'renormReal': 1.0,
    'r_H': 0.9,   'alpha0_H': 0.43, 'alpha1_H': 0.85, 'n_H': 1.35, 'b_H': 0.4,  'Mm2_H': 0.64, 'P_H': 1.0,
    'r_Ht': 7.0,  'alpha0_Ht':0.43, 'alpha1_Ht':0.85,'n_Ht':0.6,  'b_Ht':2.0,  'Mm2_Ht':0.8,  'P_Ht':1.0,
    'r_E': 0.9,   'alpha0_E': 0.43, 'alpha1_E': 0.85,'n_E': 1.35, 'b_E': 0.4,  'Mm2_E': 0.64, 'P_E': 1.0,
    'r_Et':1.0,   'alpha0_Et':0.0,  'alpha1_Et':0.0, 'n_Et': 0.0,  'b_Et': 0.0,  'Mm2_Et': 0.0,  'P_Et': 0.0,
}
shape_keys = ["r","alpha0","alpha1","n","b","Mm2","P"]

def compute_bsa(phi_arr, Q2_arr, xB_arr, t_arr, Eb_arr, param_map, flags, tag=""):
    """
    Build BMK_DVCS(-1,1,0,Eb,xB,Q2,t,phi) for each point and return dvcs.BSA().
    Uses param_map for both original and fitted parameters.
    """
    # set global renormalizations
    ROOT.renormImag = param_map.get('renormImag', orig_defaults['renormImag'])
    ROOT.renormReal = param_map.get('renormReal', orig_defaults['renormReal'])
    # set which CFFs are active and their parameters
    for cff in ("H","Ht","E","Et"):
        setattr(ROOT, f"has{cff}", bool(flags[cff]))
        if flags[cff]:
            for k in shape_keys:
                key = f"{k}_{cff}"
                val = param_map.get(key, orig_defaults[key])
                setattr(ROOT, key, val)

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

    # 1) parse fit parameters and flags
    flags, pnames, vals = parse_fit_results(fitfile)
    param_map = dict(zip(pnames, vals))
    print(">> Flags:", flags)
    print(">> Fitted parameters:", param_map)

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
    for idx, b in enumerate(bins, start=1):
        phi_data = b["phi"]
        As, sigAs = b["A"], b["sigA"]

        # build 100-point φ grid
        phi_grid = np.linspace(0, 360, 100)
        Q2g = np.full_like(phi_grid, b["Q2m"])
        xBg = np.full_like(phi_grid, b["xBm"])
        tg  = np.full_like(phi_grid, b["tm"])
        Ebg = np.full_like(phi_grid, b["Ebm"])

        # compute original & fitted BSAs
        bsas_orig = compute_bsa(phi_grid, Q2g, xBg, tg, Ebg,
                                orig_defaults, flags, tag=f"bin{idx}-orig")
        bsas_fit  = compute_bsa(phi_grid, Q2g, xBg, tg, Ebg,
                                param_map,    flags, tag=f"bin{idx}-fit")

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

        fname = (f'{outdir}/BSA_bin{idx:02d}_'
                 f'{timestamp}_Q2_{Q2m:.2f}_xB_{xBm:.3f}_t_{abs(tm):.3f}.pdf')
        fig.savefig(fname)
        print(f">> Saved bin {idx} plot to {fname}")
        plt.close(fig)

if __name__ == '__main__':
    main()