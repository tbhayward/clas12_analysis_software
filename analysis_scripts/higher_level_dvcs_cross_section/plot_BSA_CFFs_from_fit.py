#!/usr/bin/env python3
"""
plot_BSA_CFFs_from_fit.py

Plots A_LU(φ) per (Q2, xB, t, Eb) bin using the BH-denominator approximation:

  • BKM (solid blue):
        A_LU(φ) ≈ [ (s1_I / c0_BH) * sin φ ] / [ 1 + (c1_BH/c0_BH) * cos φ ]

  • RGK fit (solid orange):
        Same formula, but s1_I is evaluated after pushing the ImCFF
        parameters/flags from the provided fit_results_*.txt file.

Usage:
  python plot_BSA_CFFs_from_fit.py output/fit_results/fit_results_<TIMESTAMP>.txt
         [--data PATH]
         [--outdir output/plots]
         [--debug]
"""

import os, re, argparse
import numpy as np
import matplotlib.pyplot as plt
import ROOT  # PyROOT


# ──────────────────────────── C++ loader ────────────────────────────
def _load_dvcs_once():
    """Load DVCS_xsec.C exactly once into cling."""
    if getattr(_load_dvcs_once, "_done", False):
        return
    ROOT.gInterpreter.Declare(r"""
#ifndef __DVCS_XSEC_V2__
#define __DVCS_XSEC_V2__
#include "DVCS_xsec.C"
#endif
""")
    _load_dvcs_once._done = True


# ─────────────────────── fit-file parsing ───────────────────────────
def parse_fit_results(fname):
    with open(fname) as f:
        lines = [l.strip() for l in f if l.strip()]

    # flags line like: "H 1 Ht 0 E 1 Et 0"
    flags_line = next((l for l in lines
                       if re.match(r'^\s*H\s+[01]\s+Ht\s+[01]\s+E\s+[01]\s+Et\s+[01]\s*$', l)), None)
    if flags_line is None:
        flags_line = next((l for l in lines if l.startswith("H ")), None)
    if flags_line is None:
        raise RuntimeError("No 'H ... Ht ... E ... Et ...' flags line found in fit file")
    toks = flags_line.split()
    flags = {toks[i]: int(toks[i+1]) for i in range(0, len(toks), 2)}

    # parameter names/values
    pnames, pvals = [], None
    for i, l in enumerate(lines):
        if l.startswith("# parameters"):
            pnames = l.split()[2:]
        elif l.startswith("# values") and i+1 < len(lines):
            pvals = list(map(float, lines[i+1].split()))
    if not pnames or pvals is None:
        raise RuntimeError("Couldn't parse parameter names/values from fit file")
    param_map_fit = dict(zip(pnames, pvals))

    # optional input path for BSA data
    input_path = next((l.split(None, 1)[1].strip()
                       for l in lines if l.startswith("input")), None)

    return flags, param_map_fit, input_path


# ───────────────────── data loader (φ-binned) ───────────────────────
def load_bins(datafile):
    """
    Expect columns: phi  Q2  xB  t  Eb  A  sigA
    Split into bins when φ decreases (wrap).
    """
    bins, curr, prev_phi = [], {k: [] for k in ("phi","Q2","xB","t","Eb","A","sigA")}, None
    with open(datafile) as f:
        for line in f:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            phi, Q2, xB, t, Eb, A, sigA = map(float, line.split())
            if prev_phi is not None and phi < prev_phi:
                arr = {k: np.array(v, float) for k, v in curr.items()}
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
        arr = {k: np.array(v, float) for k, v in curr.items()}
        arr["Q2m"], arr["xBm"], arr["tm"], arr["Ebm"] = (
            float(arr["Q2"].mean()), float(arr["xB"].mean()),
            float(arr["t"].mean()),  float(arr["Eb"].mean())
        )
        bins.append(arr)
    return bins


# ─────────────── defaults for the BKM “baseline” ────────────────
def bkm_defaults():
    return {
        'renormImag': 1.0, 'renormReal': 1.0,
        # H
        'r_H': 0.9, 'n_H': 1.35, 'alpha0_H': 0.43, 'alpha1_H': 0.85, 'b_H': 0.4, 'M2_H': 0.64, 'P_H': 1.0,
        # Htilde
        'r_Ht': 7.0, 'n_Ht': 0.6, 'alpha0_Ht': 0.43, 'alpha1_Ht': 0.85, 'b_Ht': 2.0, 'M2_Ht': 0.8, 'P_Ht': 1.0,
        # E
        'r_E': 0.9, 'n_E': 1.35, 'alpha0_E': 0.43, 'alpha1_E': 0.85, 'b_E': 0.4, 'M2_E': 0.64, 'P_E': 1.0,
        # Et — keep harmless
        'r_Et': 1.0, 'n_Et': 0.0, 'alpha0_Et': 0.0, 'alpha1_Et': 0.0, 'b_Et': 0.0, 'M2_Et': 1.0, 'P_Et': 0.0,
    }


# ───────── push globals into C++ (and set switches) ───────────────
def push_globals(param_map, flags, *, safety_disable_Et=True, label=""):
    ROOT.gInterpreter.ProcessLine("hasH=0; hasHt=0; hasE=0; hasEt=0;")

    # Defensive Et handling (avoid t/0)
    et_flag = int(flags.get("Et", 0))
    if safety_disable_Et and et_flag:
        M2 = float(param_map.get("M2_Et", 0.0))
        P  = float(param_map.get("P_Et", 0.0))
        n  = float(param_map.get("n_Et", 0.0))
        if abs(M2) < 1e-12 and abs(P) > 1e-12:
            print(f"[warn] Disabling Et for {label}: M2_Et≈0 with P_Et≠0.")
            et_flag = 0
        if n == 0.0:
            et_flag = 0

    for cff, on in (("H",flags.get("H",0)), ("Ht",flags.get("Ht",0)),
                    ("E",flags.get("E",0)),  ("Et",et_flag)):
        ROOT.gInterpreter.ProcessLine(f"has{cff} = {int(on)};")

    ROOT.gInterpreter.ProcessLine(f"renormImag = {float(param_map.get('renormImag',1.0))};")
    ROOT.gInterpreter.ProcessLine(f"renormReal = {float(param_map.get('renormReal',1.0))};")

    def setp(name):
        if name in param_map:
            ROOT.gInterpreter.ProcessLine(f"{name} = {float(param_map[name])};")

    for cff in ("H","Ht","E","Et"):
        for k in ("r","n","alpha0","alpha1","b","M2","P"):
            setp(f"{k}_{cff}")

    for cff in ("H","Ht","E","Et"):
        for k in ("C0","MD2","lambda"):
            setp(f"{k}_{cff}")


# ─────────────── helpers to compute A & B ────────────────
def get_A_B(Q2, xB, t, Eb):
    """
    Return A = s1_I/c0_BH and B = c1_BH/c0_BH at fixed kinematics.
    (s1_I carries the ImCFF dependence; c0_BH,c1_BH are BH-only.)
    """
    dv = ROOT.BMK_DVCS(-1, +1, 0, Eb, xB, Q2, t, 0.0)
    c0 = dv.c0_BH()
    c1 = dv.c1_BH()
    s1 = dv.s1_I()
    if c0 == 0:
        return 0.0, 0.0
    return s1 / c0, c1 / c0


def approx_curve(phi_deg, A, B):
    """A_LU(φ) ≈ [A sinφ] / [1 + B cosφ]."""
    ph = np.deg2rad(phi_deg)
    denom = 1.0 + B * np.cos(ph)
    # avoid accidental blow-ups
    eps = 1e-14
    denom = np.where(np.abs(denom) < eps, np.sign(denom) * eps, denom)
    return (A * np.sin(ph)) / denom


# ───────────────────────────── main ───────────────────────────────
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('fitfile', help='fit_results_<TIMESTAMP>.txt from your fitter')
    ap.add_argument('--data', default=None, help='Override BSA data file; else use "input" in fit file or fallback')
    ap.add_argument('--outdir', default='output/plots', help='Output directory for PDFs')
    ap.add_argument('--debug', action='store_true', help='Print internal diagnostics for first bin')
    args = ap.parse_args()

    # Parse fit file
    flags_fit, param_map_fit, input_path = parse_fit_results(args.fitfile)

    # Decide data file
    datafile = (args.data or input_path or 'imports/rgk_preliminary_bsa.txt')
    print(">> Using data file:", datafile)
    print(">> Fit flags:", flags_fit)

    # Load C++
    _load_dvcs_once()

    # Load bins
    bins = load_bins(datafile)
    print(f">> Found {len(bins)} φ-bins")

    os.makedirs(args.outdir, exist_ok=True)

    # For filename timestamp
    m = re.search(r'fit_results_(\d{8}_\d{6})\.txt$', os.path.basename(args.fitfile))
    ts = (m.group(1) if m else "noTS")

    # BKM baseline flags: H,Ht,E on; Et off
    flags_bkm = {'H':1,'Ht':1,'E':1,'Et':0}
    params_bkm = bkm_defaults()

    # Loop bins
    for ibin, b in enumerate(bins, start=1):
        phi_data, A_data, dA = b["phi"], b["A"], b["sigA"]
        Q2m, xBm, tm, Ebm = b["Q2m"], b["xBm"], b["tm"], b["Ebm"]

        # Dense grid
        phi_grid = np.linspace(0.0, 360.0, 361)

        # --- BKM baseline (A,B from defaults) ---
        push_globals(params_bkm, flags_bkm, label="BKM")
        if args.debug and ibin == 1:
            ROOT.gInterpreter.ProcessLine(r'std::cout<<"[BKM] hasH="<<hasH<<" hasHt="<<hasHt<<" hasE="<<hasE<<" hasEt="<<hasEt<<std::endl;')
        A_bkm, B_bkm = get_A_B(Q2m, xBm, tm, Ebm)
        bkm_curve = approx_curve(phi_grid, A_bkm, B_bkm)

        # --- RGK fit (A from fit ImCFFs; B still BH-only) ---
        push_globals(param_map_fit, flags_fit, label="RGK-fit")
        if args.debug and ibin == 1:
            ROOT.gInterpreter.ProcessLine(r'std::cout<<"[RGK] hasH="<<hasH<<" hasHt="<<hasHt<<" hasE="<<hasE<<" hasEt="<<hasEt<<std::endl;')
        A_rgk, B_rgk = get_A_B(Q2m, xBm, tm, Ebm)  # c0,c1 (BH) same; s1_I reflects fit
        rgk_curve = approx_curve(phi_grid, A_rgk, B_rgk)

        # --- Plot ---
        fig, ax = plt.subplots(figsize=(8,5))
        ax.errorbar(phi_data, A_data, yerr=dA, fmt='o', ms=5, color='k', label='Data')

        ax.plot(phi_grid, bkm_curve, '-',  lw=2, color='tab:blue',   label='BKM')
        ax.plot(phi_grid, rgk_curve, '-',  lw=2, color='tab:orange', label='RGK fit')

        ax.set_xlim(0, 360)
        ax.set_xticks([0,60,120,180,240,300,360])
        ax.set_ylim(-0.6, 0.6)
        ax.set_xlabel(r'$\phi\;(\mathrm{deg})$')
        ax.set_ylabel(r'$A_{LU}(\phi)$')
        ax.set_title(
            (r'$\langle Q^2\rangle={:.2f}\,\mathrm{{GeV}}^2,\;\langle x_B\rangle={:.3f},\;'
             r'\langle -t\rangle={:.3f}\,\mathrm{{GeV}}^2,\;\langle E_b\rangle={:.2f}\,\mathrm{{GeV}}$'
            ).format(Q2m, xBm, abs(tm), Ebm),
            pad=12
        )
        ax.legend(loc='upper right', frameon=True, edgecolor='k')
        plt.tight_layout()

        outname = os.path.join(
            args.outdir,
            f"BSA_bin{ibin:02d}_{ts}_Q2_{Q2m:.2f}_xB_{xBm:.3f}_mt_{abs(tm)::.3f}.pdf"
        )
        fig.savefig(outname)
        print(f">> Saved bin {ibin} plot to {outname}")
        plt.close(fig)


if __name__ == "__main__":
    main()