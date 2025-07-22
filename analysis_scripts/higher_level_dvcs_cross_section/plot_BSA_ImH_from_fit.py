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

def load_all_bins(datafile):
    """
    Read imports/rga_prl_bsa.txt, split into bins when φ wraps around.
    Returns a list of dicts with arrays and their means.
    """
    bins = []
    current = {k: [] for k in ("phi","Q2","xB","t","Eb","A","sigA")}
    prev_phi = None

    with open(datafile) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            phi, Q2, xB, t, Eb, A, sigA = map(float, line.split())
            # new bin if φ decreases
            if prev_phi is not None and phi < prev_phi:
                arr = {k: np.array(v) for k, v in current.items()}
                arr["Q2m"] = arr["Q2"].mean()
                arr["xBm"] = arr["xB"].mean()
                arr["tm"]  = arr["t"].mean()
                arr["Ebm"] = arr["Eb"].mean()
                bins.append(arr)
                current = {k: [] for k in current}
            # append to current
            current["phi"].append(phi)
            current["Q2"].append(Q2)
            current["xB"].append(xB)
            current["t"].append(t)
            current["Eb"].append(Eb)
            current["A"].append(A)
            current["sigA"].append(sigA)
            prev_phi = phi

    # last bin
    if current["phi"]:
        arr = {k: np.array(v) for k, v in current.items()}
        arr["Q2m"] = arr["Q2"].mean()
        arr["xBm"] = arr["xB"].mean()
        arr["tm"]  = arr["t"].mean()
        arr["Ebm"] = arr["Eb"].mean()
        bins.append(arr)

    return bins

def compute_bsa(phi_arr, Q2_arr, xB_arr, t_arr, Eb_arr, params, tag=""):
    """
    For each i, build BMK_DVCS(-1,1,0,Eb,xB,Q2,t,phi) and return dvcs.BSA().
    Prints the first five for quick debug.
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

def cluster_bins_by_t(bins, tol=0.01):
    """
    Group bins by their mean t within tol. Returns list of
    {'center':t_center, 'indices':[bin_idx,...]}.
    """
    clusters = []
    for idx, b in enumerate(bins):
        t0 = b["tm"]
        placed = False
        for cl in clusters:
            if abs(t0 - cl["center"]) <= tol:
                cl["indices"].append(idx)
                # update cluster center
                cl["center"] = np.mean([bins[i]["tm"] for i in cl["indices"]])
                placed = True
                break
        if not placed:
            clusters.append({"center": t0, "indices": [idx]})
    return clusters

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

    # 3) load all φ‐bins
    datafile = 'imports/rga_prl_bsa.txt'
    bins = load_all_bins(datafile)
    print(f">> Found {len(bins)} φ‐bins")

    # 4) cluster by t
    clusters = cluster_bins_by_t(bins, tol=0.01)
    print(f">> Grouped into {len(clusters)} |t|-clusters")

    # 5) loop over t-clusters
    for cl in clusters:
        t_center = cl["center"]
        idxs     = cl["indices"]

        # unique sorted grid points
        Q2_vals = sorted({round(bins[i]["Q2m"], 3) for i in idxs})
        xB_vals = sorted({round(bins[i]["xBm"], 3) for i in idxs})
        nrows, ncols = len(Q2_vals), len(xB_vals)
        print(f"Cluster |t|≈{abs(t_center):.3f} → grid {nrows}×{ncols}")

        # prepare figure & axes
        fig, axes = plt.subplots(nrows, ncols,
                                 figsize=(4*ncols, 3*nrows),
                                 sharex=True, sharey=True)
        # unify axes array shape
        if nrows == 1 and ncols == 1:
            axes = np.array([[axes]])
        elif nrows == 1:
            axes = axes[np.newaxis, :]
        elif ncols == 1:
            axes = axes[:, np.newaxis]

        # top title with mean -t
        fig.suptitle(rf'$\langle -t\rangle = {abs(t_center):.3f}\,\mathrm{{GeV}}^2$',
                     y=0.98, fontsize=14)

        # for each bin in this cluster
        for i in idxs:
            b = bins[i]
            q2m = round(b["Q2m"], 3)
            xbm = round(b["xBm"], 3)
            row = Q2_vals.index(q2m)
            col = xB_vals.index(xbm)
            ax  = axes[row, col]

            # data
            ax.errorbar(b["phi"], b["A"], yerr=b["sigA"],
                        fmt='o', ms=4, color='k')

            # model on fine φ grid
            phi_grid = np.linspace(0, 360, 100)
            Q2g = np.full_like(phi_grid, b["Q2m"])
            xBg = np.full_like(phi_grid, b["xBm"])
            tg  = np.full_like(phi_grid, b["tm"])
            Ebg = np.full_like(phi_grid, b["Ebm"])

            defaults     = [1.0, 0.43, 0.85, 1.35, 0.4, 0.64, 1.0, 1.0]
            bsas_orig    = compute_bsa(phi_grid, Q2g, xBg, tg, Ebg,
                                       defaults, tag="orig-grid")
            bsas_fitted  = compute_bsa(phi_grid, Q2g, xBg, tg, Ebg,
                                       params_fit, tag="fit-grid")

            ax.plot(phi_grid, bsas_orig,   '-',  lw=2, color='tab:blue')
            ax.plot(phi_grid, bsas_fitted, '--', lw=2, color='tab:red')

            # subplot title
            ax.set_title(rf'$Q^2={b["Q2m"]:.2f},\,x_B={b["xBm"]:.3f}$',
                         pad=6, fontsize=10)

        # label only bottom row & left column
        for col in range(ncols):
            ax = axes[-1, col]
            ax.set_xlabel(r'$\phi\;[\mathrm{deg}]$')
            ax.set_xticks([0, 60, 120, 180, 240, 300, 360])
        for row in range(nrows):
            ax = axes[row, 0]
            ax.set_ylabel(r'$A_{LU}(\phi)$')
            ax.set_ylim(-0.6, 0.6)

        # one legend in the top-right subplot
        leg_ax = axes[0, -1]
        leg_ax.legend(["Data", "Original Model", "Fitted Model"],
                      loc='upper right', frameon=True, edgecolor='k')

        plt.tight_layout(rect=[0, 0, 1, 0.96])

        # save
        outdir = 'output/plots'
        os.makedirs(outdir, exist_ok=True)
        tfmt = f"{abs(t_center):.3f}".replace('.', 'p')
        outname = f'{outdir}/BSA_t_{tfmt}_{timestamp}.pdf'
        fig.savefig(outname)
        print(">> Saved", outname)
        plt.close(fig)


if __name__ == '__main__':
    main()