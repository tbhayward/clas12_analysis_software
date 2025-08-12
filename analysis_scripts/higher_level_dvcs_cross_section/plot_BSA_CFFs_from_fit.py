#!/usr/bin/env python3
"""
plot_BSA_default.py

Plot A_LU(φ) in each (Q2, xB, t, Eb) bin from a BSA text file using the *default* model
in DVCS_xsec.C (all CFFs ON, default parameters, renormImag=renormReal=1).

Usage:
    python plot_BSA_default.py --data imports/rgk_preliminary_bsa.txt [--outdir output/plots] [--bh-limit] [--debug]
"""

import os
import argparse
import time
import numpy as np
import matplotlib.pyplot as plt
import ROOT

# --------------------------------------------------------------------------------------
# ROOT / C++ glue
#   Guard the include so DVCS_xsec.C isn't pulled in multiple times (cling remembers macros)
# --------------------------------------------------------------------------------------
ROOT.gInterpreter.Declare(r"""
#ifndef __DVCS_XSEC_ONCE__
#define __DVCS_XSEC_ONCE__
#include "DVCS_xsec.C"
#endif
""")

def set_default_model(force_all_cffs=True):
    """
    Push the 'default' model into the C++ globals:
      - renormImag=renormReal=1
      - hasH/Ht/E/Et = 1 (if force_all_cffs) so default curve is non-zero
      - (re)set parameter values to the DVCS_xsec.C defaults (explicitly, to be bulletproof)
    """
    # Renormalizations
    ROOT.gInterpreter.ProcessLine("renormImag = 1.0;")
    ROOT.gInterpreter.ProcessLine("renormReal = 1.0;")

    # CFF flags
    if force_all_cffs:
        ROOT.gInterpreter.ProcessLine("hasH  = 1;")
        ROOT.gInterpreter.ProcessLine("hasHt = 1;")
        ROOT.gInterpreter.ProcessLine("hasE  = 1;")
        ROOT.gInterpreter.ProcessLine("hasEt = 1;")

    # Explicitly mirror DVCS_xsec.C defaults (harmless if already set)
    ROOT.gInterpreter.ProcessLine("r_H=0.9; n_H=1.25; alpha0_H=0.43; alpha1_H=0.85; b_H=0.4; M2_H=0.64; P_H=1.0;")
    ROOT.gInterpreter.ProcessLine("r_Ht=7.0; n_Ht=0.6; alpha0_Ht=0.43; alpha1_Ht=0.85; b_Ht=2.0; M2_Ht=0.8; P_Ht=1.0;")
    ROOT.gInterpreter.ProcessLine("r_E=0.9; n_E=1.25; alpha0_E=0.43; alpha1_E=0.85; b_E=0.4; M2_E=0.64; P_E=1.0;")
    ROOT.gInterpreter.ProcessLine("r_Et=1.0; n_Et=0.6; alpha0_Et=0.0; alpha1_Et=0.0; b_Et=0.0; M2_Et=0.0; P_Et=0.0;")


# --------------------------------------------------------------------------------------
# Data loading & binning (bins defined by φ wrap-around)
# --------------------------------------------------------------------------------------
def load_phi_bins(datafile):
    """
    Read a BSA text file with columns:
        phi  Q2  xB  t  Eb  A  sigA
    and split into bins each time φ decreases (wrap-around).
    Returns list of dicts with arrays and per-bin means (Q2m,xBm,tm,Ebm).
    """
    bins = []
    curr = {k: [] for k in ("phi","Q2","xB","t","Eb","A","sigA")}
    prev_phi = None
    with open(datafile) as f:
        for line in f:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            phi, Q2, xB, t, Eb, A, sigA = map(float, line.split())
            if prev_phi is not None and phi < prev_phi:
                arr = {k: np.array(v, dtype=float) for k, v in curr.items()}
                arr["Q2m"], arr["xBm"], arr["tm"], arr["Ebm"] = (
                    float(arr["Q2"].mean()), float(arr["xB"].mean()),
                    float(arr["t"].mean()),  float(arr["Eb"].mean())
                )
                bins.append(arr)
                curr = {k: [] for k in curr}
            for k, v in zip(curr.keys(), (phi, Q2, xB, t, Eb, A, sigA)):
                curr[k].append(v)
            prev_phi = phi
    if curr["phi"]:
        arr = {k: np.array(v, dtype=float) for k, v in curr.items()}
        arr["Q2m"], arr["xBm"], arr["tm"], arr["Ebm"] = (
            float(arr["Q2"].mean()), float(arr["xB"].mean()),
            float(arr["t"].mean()),  float(arr["Eb"].mean())
        )
        bins.append(arr)
    return bins


# --------------------------------------------------------------------------------------
# Physics helpers
# --------------------------------------------------------------------------------------
def bsa_curve(phi_deg, Q2, xB, t, Eb):
    """
    Compute full BSA(φ) from BMK_DVCS::BSA() for arrays of φ (degrees),
    holding (Q2,xB,t,Eb) fixed for the bin.
    """
    out = np.empty_like(phi_deg, dtype=float)
    for i, ph in enumerate(phi_deg):
        dv = ROOT.BMK_DVCS(-1, +1, 0, Eb, xB, Q2, t, ph)  # ctor takes φ in degrees
        out[i] = dv.BSA()
    return out


def bh_limit_sinusoid(Q2, xB, t, Eb):
    """
    Return the BH-dominant approximation A(φ) ≈ [s1_I / c0_BH] * sin(φ).
    Note s1_I and c0_BH don't depend on φ (for unpolarized target), so we evaluate once.
    """
    dv = ROOT.BMK_DVCS(-1, +1, 0, Eb, xB, Q2, t, 0.0)
    amp = dv.s1_I() / dv.c0_BH() if dv.c0_BH() != 0 else 0.0
    def A(phi_deg):
        return amp * np.sin(np.deg2rad(phi_deg))
    return A


# --------------------------------------------------------------------------------------
# Plotting
# --------------------------------------------------------------------------------------
def plot_bins(bins, outdir, show_bh_limit=False, debug=False):
    os.makedirs(outdir, exist_ok=True)
    ts = time.strftime("%Y%m%d_%H%M%S")

    for ibin, b in enumerate(bins, start=1):
        phi_data, A_data, dA_data = b["phi"], b["A"], b["sigA"]
        Q2m, xBm, tm, Ebm = b["Q2m"], b["xBm"], b["tm"], b["Ebm"]

        # Dense φ grid
        phi_grid = np.linspace(0.0, 360.0, 361)

        # Full model curve
        A_full = bsa_curve(phi_grid, Q2m, xBm, tm, Ebm)

        if debug and ibin == 1:
            # Quick sanity: print at φ where sin is large
            for probe in (30.0, 90.0, 150.0):
                dv = ROOT.BMK_DVCS(-1, +1, 0, Ebm, xBm, Q2m, tm, probe)
                print(f"[bin{ibin:02d} dbg] φ={probe:6.1f}°  BSA={dv.BSA():+.6e}  "
                      f"BH2={dv.BH2():.3e}  DVCS2={dv.DVCS2():.3e}  I={dv.BHDVCS():.3e}  s1_I={dv.s1_I():+.6e}")

        # Plot
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.errorbar(phi_data, A_data, yerr=dA_data, fmt='o', ms=5, color='k', label='Data')

        ax.plot(phi_grid, A_full, '-', lw=2, label='Default model (full BSA)')

        if show_bh_limit:
            A_bh = bh_limit_sinusoid(Q2m, xBm, tm, Ebm)
            ax.plot(phi_grid, A_bh(phi_grid), '--', lw=2, label=r'BH limit: $\frac{s_1^I}{c_0^{BH}}\sin\phi$')

        ax.set_xlim(0, 360)
        ax.set_xticks([0, 60, 120, 180, 240, 300, 360])
        ax.set_ylim(-0.6, 0.6)
        ax.set_xlabel(r'$\phi\;\mathrm{[deg]}$')
        ax.set_ylabel(r'$A_{LU}(\phi)$')

        ax.set_title(
            r'$\langle Q^2\rangle={:.2f}\,\mathrm{{GeV}}^2,\;\langle x_B\rangle={:.3f},\;\langle -t\rangle={:.3f}\,\mathrm{{GeV}}^2$'
            .format(Q2m, xBm, -tm),
            pad=10
        )
        ax.legend(loc='upper right', frameon=True)
        plt.tight_layout()

        fname = os.path.join(outdir, f"BSA_default_bin{ibin:02d}_{ts}_Q2_{Q2m:.2f}_xB_{xBm:.3f}_t_{abs(tm):.3f}.pdf")
        fig.savefig(fname)
        print(f">> Saved bin {ibin} plot to {fname}")
        plt.close(fig)


# --------------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default="imports/rgk_preliminary_bsa.txt",
                    help="BSA text file (columns: phi Q2 xB t Eb A sigA)")
    ap.add_argument("--outdir", default="output/plots", help="Output directory for PDFs")
    ap.add_argument("--bh-limit", action="store_true", help="Overlay BH-limit sinusoid A ~ (s1_I/c0_BH) sin φ")
    ap.add_argument("--debug", action="store_true", help="Print a few diagnostic values for the first bin")
    args = ap.parse_args()

    set_default_model(force_all_cffs=True)
    bins = load_phi_bins(args.data)
    if not bins:
        raise RuntimeError(f"No bins parsed from {args.data}")

    print(f">> Loaded {len(bins)} φ-bins from {args.data}")
    plot_bins(bins, args.outdir, show_bh_limit=args.bh_limit, debug=args.debug)

if __name__ == "__main__":
    main()