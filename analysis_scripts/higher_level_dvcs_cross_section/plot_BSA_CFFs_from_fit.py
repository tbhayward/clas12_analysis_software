#!/usr/bin/env python3
"""
plot_BSA_CFFs_from_fit.py

Plots A_LU(φ) per (Q2, xB, t, Eb) bin using the BH-only approximation:
  • BKM:     A1*sinφ / (1 + B1*cosφ), with  A1 = (s1_I/c0_BH) * [xB/y * (1+eps^2)^2]
             evaluated using the BKM default CFF parameters (H,H~ ,E on; Et off).
  • RGK fit: same formula, but A1 uses s1_I computed with your fitted Im-CFFs.
             (B1 = c1_BH/c0_BH in both curves; BH coefficients depend only on kinematics.)

Usage:
  python plot_BSA_CFFs_from_fit.py output/fit_results/fit_results_<TIMESTAMP>.txt
         [--data PATH] [--outdir output/plots] [--debug]
"""

import os, re, argparse
import numpy as np
import matplotlib.pyplot as plt
import ROOT  # PyROOT

# ──────────────────────────── C++ loader ────────────────────────────
def _load_dvcs_once():
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

    # optional input path
    input_path = next((l.split(None, 1)[1].strip()
                       for l in lines if l.startswith("input")), None)

    return flags, param_map_fit, input_path

# ───────────────────── data loader (φ-binned) ───────────────────────
def load_bins(datafile):
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
        # Et — keep harmless (imaginary effectively off)
        'r_Et': 1.0, 'n_Et': 0.0, 'alpha0_Et': 0.0, 'alpha1_Et': 0.0, 'b_Et': 0.0, 'M2_Et': 1.0, 'P_Et': 0.0,
    }

# ───────── push globals into C++ (and set switches) ───────────────
def push_globals(param_map, flags, *, safety_disable_Et=True, label=""):
    ROOT.gInterpreter.ProcessLine("hasH=0; hasHt=0; hasE=0; hasEt=0;")

    # Defensive Et handling
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

# ─────────────── helpers: BH coeffs, A1B1 curve, normalization ───────────────
def bh_coeffs_and_s1I(Q2, xB, t, Eb):
    """Return (s1_I, c0_BH, c1_BH) at fixed kinematics (φ-independent)."""
    dv = ROOT.BMK_DVCS(-1, +1, 0, Eb, xB, Q2, t, 0.0)
    s1I = dv.s1_I()
    c0  = dv.c0_BH()
    c1  = dv.c1_BH()
    return s1I, c0, c1

def approx_curve_A1B1(phi_deg, A1, B1):
    """A(φ) = A1 * sinφ / (1 + B1 cosφ)."""
    ph = np.deg2rad(phi_deg)
    denom = 1.0 + B1*np.cos(ph)
    eps = 1e-14
    denom = np.where(np.abs(denom) < eps, np.sign(denom)*eps, denom)
    return A1 * np.sin(ph) / denom

def bh_norm_ratio(Q2, xB, Eb):
    """
    Overall factor converting (s1_I / c0_BH) into the BH-only A1 for A_LU:
        scale = (xB / y) * (1 + eps^2)^2
      with y = Q2 / (2 M xB Eb), eps = 2 xB M / sqrt(Q2).
    """
    M = 0.93827
    y = Q2 / (2.0 * M * xB * Eb)
    eps2 = (2.0 * xB * M / np.sqrt(Q2))**2
    return (xB / y) * (1.0 + eps2)**2

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

    # Data file
    datafile = (args.data or input_path or 'imports/rgk_preliminary_bsa.txt')
    print(">> Using data file:", datafile)
    print(">> Fit flags:", flags_fit)

    # Load C++
    _load_dvcs_once()

    # Load bins
    bins = load_bins(datafile)
    print(f">> Found {len(bins)} φ-bins")

    os.makedirs(args.outdir, exist_ok=True)

    # BKM baseline flags: H,Ht,E on; Et off
    flags_bkm = {'H':1,'Ht':1,'E':1,'Et':0}
    params_bkm = bkm_defaults()

    # Timestamp for filenames
    m = re.search(r'fit_results_(\d{8}_\d{6})\.txt$', os.path.basename(args.fitfile))
    ts = (m.group(1) if m else "noTS")

    # Loop bins
    for ibin, b in enumerate(bins, start=1):
        phi_data, A_data, dA = b["phi"], b["A"], b["sigA"]
        Q2m, xBm, tm, Ebm = b["Q2m"], b["xBm"], b["tm"], b["Ebm"]

        # Dense grid
        phi_grid = np.linspace(0.0, 360.0, 361)

        # Common kinematic scale for A1
        scale = bh_norm_ratio(Q2m, xBm, Ebm)
        if args.debug and ibin == 1:
            M = 0.93827
            y = Q2m / (2.0 * M * xBm * Ebm)
            eps2 = (2.0 * xBm * M / np.sqrt(Q2m))**2
            print(f"[dbg] scale=(xB/y)*(1+eps^2)^2 = {scale:.6e}  (y={y:.6e}, eps2={eps2:.6e})")

        # --- BKM A1,B1 ---
        push_globals(params_bkm, flags_bkm, label="BKM")
        s1I_bkm, c0_bh, c1_bh = bh_coeffs_and_s1I(Q2m, xBm, tm, Ebm)
        A1_bkm = (s1I_bkm / c0_bh) * scale
        B1_bh  = c1_bh / c0_bh
        bkm_curve = approx_curve_A1B1(phi_grid, A1_bkm, B1_bh)

        # --- RGK fit A1 (same B1 from BH) ---
        push_globals(param_map_fit, flags_fit, label="RGK-fit")
        s1I_rgk, c0_bh_chk, c1_bh_chk = bh_coeffs_and_s1I(Q2m, xBm, tm, Ebm)
        if args.debug and ibin == 1:
            print(f"[dbg] c0_BH={c0_bh:.6e}  c1_BH={c1_bh:.6e}  "
                  f"(RGK pass: c0={c0_bh_chk:.6e}, c1={c1_bh_chk:.6e})")
            print(f"[dbg] A1_BKM={(s1I_bkm/c0_bh):.6e}*scale -> {A1_bkm:.6e}; "
                  f"A1_RGK={(s1I_rgk/c0_bh_chk):.6e}*scale -> {(s1I_rgk/c0_bh_chk)*scale:.6e}")
        # use BH from either pass (should be identical)
        A1_rgk = (s1I_rgk / c0_bh_chk) * scale
        B1_used = c1_bh_chk / c0_bh_chk
        rgk_curve = approx_curve_A1B1(phi_grid, A1_rgk, B1_used)

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
            f"BSA_bin{ibin:02d}_{ts}_Q2_{Q2m:.2f}_xB_{xBm:.3f}_mt_{abs(tm):.3f}.pdf"
        )
        fig.savefig(outname)
        print(f">> Saved bin {ibin} plot to {outname}")
        plt.close(fig)

if __name__ == "__main__":
    main()