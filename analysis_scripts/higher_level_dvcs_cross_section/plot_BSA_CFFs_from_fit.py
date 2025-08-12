#!/usr/bin/env python3
"""
plot_BSA_CFFs_from_fit.py

Plot A_LU(φ) per (Q2, xB, t, Eb) bin for:
  • BKM baseline ("BKM")
  • Your fitted model ("RGK 6 & 7 GeV Fit")
Optional overlays:
  • ImH-only curve from the fit (--CFFs 1)
  • BH-dominant sinusoid A ≈ (s1_I/c0_BH) sin φ (--bh-limit)

Usage:
  python plot_BSA_CFFs_from_fit.py output/fit_results/fit_results_<TIMESTAMP>.txt
         [--data PATH] [--CFFs 0|1] [--bh-limit] [--outdir output/plots] [--debug]

Notes:
  - Requires DVCS_xsec.C in your CWD (with your Et fix or Et disabled).
  - This loader uses a NEW include guard so cling actually compiles the current file.
"""

import os, re, argparse
import numpy as np
import matplotlib.pyplot as plt
import ROOT  # PyROOT

# -----------------------------------------------------------------------------
# C++ loader — use a fresh macro guard to ensure the current DVCS_xsec.C is compiled
# -----------------------------------------------------------------------------
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


# -----------------------------------------------------------------------------
# Fit-file parsing
# -----------------------------------------------------------------------------
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

    # optional input and strategy
    input_path = next((l.split(None, 1)[1].strip()
                       for l in lines if l.startswith("input")), None)
    strategy = None
    sline = next((l for l in lines if l.startswith("strategy")), None)
    if sline:
        try: strategy = int(sline.split()[1])
        except: pass

    return flags, param_map_fit, input_path, strategy


# -----------------------------------------------------------------------------
# Data loading (bins defined by φ wrap-around)
# -----------------------------------------------------------------------------
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


# -----------------------------------------------------------------------------
# BKM "original" baseline defaults (safe wrt Et)
# -----------------------------------------------------------------------------
def bkm_defaults():
    return {
        'renormImag': 1.0, 'renormReal': 1.0,
        # H
        'r_H': 0.9, 'n_H': 1.25, 'alpha0_H': 0.43, 'alpha1_H': 0.85, 'b_H': 0.4, 'M2_H': 0.64, 'P_H': 1.0,
        # Htilde
        'r_Ht': 7.0, 'n_Ht': 0.6, 'alpha0_Ht': 0.43, 'alpha1_Ht': 0.85, 'b_Ht': 2.0, 'M2_Ht': 0.8, 'P_Ht': 1.0,
        # E
        'r_E': 0.9, 'n_E': 1.25, 'alpha0_E': 0.43, 'alpha1_E': 0.85, 'b_E': 0.4, 'M2_E': 0.64, 'P_E': 1.0,
        # Et — keep harmless: n_Et=0, P_Et=0, M2_Et nonzero to avoid accidental t/0 evaluation
        'r_Et': 1.0, 'n_Et': 0.0, 'alpha0_Et': 0.0, 'alpha1_Et': 0.0, 'b_Et': 0.0, 'M2_Et': 1.0, 'P_Et': 0.0,
    }


# -----------------------------------------------------------------------------
# Push globals into the C++ world (with a little safety for Et)
# -----------------------------------------------------------------------------
def push_globals(param_map, flags, *, safety_disable_Et=True, label=""):
    # Beam-target CFF switches
    # (We *don’t* trust previous state.)
    ROOT.gInterpreter.ProcessLine("hasH=0; hasHt=0; hasE=0; hasEt=0;")

    # Defensive Et handling — only enable if not obviously unsafe
    et_flag = int(flags.get("Et", 0))
    if safety_disable_Et and et_flag:
        M2 = float(param_map.get("M2_Et", 0.0))
        P  = float(param_map.get("P_Et", 0.0))
        n  = float(param_map.get("n_Et", 0.0))
        if abs(M2) < 1e-12 and abs(P) > 1e-12:
            print(f"[warn] Disabling Et for {label}: M2_Et≈0 with P_Et≠0 would cause t/0.")
            et_flag = 0
        if n == 0.0:
            # No imaginary amplitude anyway — harmless to keep off.
            et_flag = 0

    for cff, on in (("H",flags.get("H",0)), ("Ht",flags.get("Ht",0)),
                    ("E",flags.get("E",0)),  ("Et",et_flag)):
        ROOT.gInterpreter.ProcessLine(f"has{cff} = {int(on)};")

    # Renormalizations
    ROOT.gInterpreter.ProcessLine(f"renormImag = {float(param_map.get('renormImag',1.0))};")
    ROOT.gInterpreter.ProcessLine(f"renormReal = {float(param_map.get('renormReal',1.0))};")

    # Set ansatz parameters that exist in the map
    def setp(name):
        if name in param_map:
            ROOT.gInterpreter.ProcessLine(f"{name} = {float(param_map[name])};")

    for cff in ("H","Ht","E","Et"):
        for k in ("r","n","alpha0","alpha1","b","M2","P"):
            setp(f"{k}_{cff}")

    # Optional real-part subtraction constants (only if present)
    for cff in ("H","Ht","E","Et"):
        for k in ("C0","MD2","lambda"):
            setp(f"{k}_{cff}")


# -----------------------------------------------------------------------------
# Physics helpers
# -----------------------------------------------------------------------------
def bsa_curve(phi_deg, Q2, xB, t, Eb, *, debug=False, tag="", max_dbg_pts=3):
    out = np.empty_like(phi_deg, float)
    for i, ph in enumerate(phi_deg):
        dv = ROOT.BMK_DVCS(-1, +1, 0, Eb, xB, Q2, t, ph)  # ctor takes φ in degrees
        out[i] = dv.BSA()
        if debug and i < max_dbg_pts:
            print(f"[{tag}] φ={ph:6.1f}°  BSA={out[i]:+.6e}  "
                  f"BH2={dv.BH2():.3e}  DVCS2={dv.DVCS2():.3e}  I={dv.BHDVCS():.3e}  s1_I={dv.s1_I():+.6e}")
    return out

def bh_limit_sinusoid(Q2, xB, t, Eb):
    dv = ROOT.BMK_DVCS(-1, +1, 0, Eb, xB, Q2, t, 0.0)
    amp = dv.s1_I()/dv.c0_BH() if dv.c0_BH()!=0 else 0.0
    return lambda phi_deg: amp * np.sin(np.deg2rad(phi_deg))


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('fitfile', help='fit_results_<TIMESTAMP>.txt from your fitter')
    ap.add_argument('--data', default=None, help='Override BSA data file; else use "input" line in fit file or fallback')
    ap.add_argument('--CFFs', type=int, choices=[0,1], default=0,
                    help='If 1, also show ImH-only curve from the fit')
    ap.add_argument('--bh-limit', action='store_true',
                    help='Overlay BH-limit: A ≈ (s1_I/c0_BH) sin φ')
    ap.add_argument('--outdir', default='output/plots', help='Output directory for PDFs')
    ap.add_argument('--debug', action='store_true', help='Print internal diagnostics for first bin')
    args = ap.parse_args()

    # Parse fit file
    flags_fit, param_map_fit, input_path, _strategy = parse_fit_results(args.fitfile)

    # Decide data file
    datafile = (args.data or input_path or 'imports/rgk_preliminary_bsa.txt')
    print(">> Using data file:", datafile)
    print(">> Fit flags:", flags_fit)
    print(">> Fit params:", param_map_fit)

    # Load C++
    _load_dvcs_once()

    # Load bins
    bins = load_bins(datafile)
    if not bins:
        raise RuntimeError(f"No φ-bins parsed from {datafile}")
    print(f">> Found {len(bins)} φ-bins")

    os.makedirs(args.outdir, exist_ok=True)

    # BKM baseline flags: H,Ht,E on; Et off (safe default)
    flags_bkm = {'H':1,'Ht':1,'E':1,'Et':0}
    params_bkm = bkm_defaults()

    # For filename timestamp
    m = re.search(r'fit_results_(\d{8}_\d{6})\.txt$', os.path.basename(args.fitfile))
    ts = (m.group(1) if m else "noTS")

    # Loop bins
    for ibin, b in enumerate(bins, start=1):
        phi_data, A_data, dA = b["phi"], b["A"], b["sigA"]
        Q2m, xBm, tm, Ebm = b["Q2m"], b["xBm"], b["tm"], b["Ebm"]

        # Dense grid
        phi_grid = np.linspace(0.0, 360.0, 361)
        Q2g = np.full_like(phi_grid, Q2m)
        xBg = np.full_like(phi_grid, xBm)
        tg  = np.full_like(phi_grid, tm)
        Ebg = np.full_like(phi_grid, Ebm)

        # --- BKM baseline ---
        push_globals(params_bkm, flags_bkm, label="BKM")
        if args.debug and ibin == 1:
            ROOT.gInterpreter.ProcessLine(r'std::cout<<"[BKM] hasH="<<hasH<<" hasHt="<<hasHt<<" hasE="<<hasE<<" hasEt="<<hasEt<<std::endl;')
        bkm_curve = bsa_curve(phi_grid, Q2m, xBm, tm, Ebm,
                              debug=(args.debug and ibin==1), tag=f"bin{ibin}-BKM")

        # --- RGK 6 & 7 GeV Fit (respect fit flags) ---
        push_globals(param_map_fit, flags_fit, label="RGK-fit")
        if args.debug and ibin == 1:
            ROOT.gInterpreter.ProcessLine(r'std::cout<<"[RGK] hasH="<<hasH<<" hasHt="<<hasHt<<" hasE="<<hasE<<" hasEt="<<hasEt<<std::endl;')
        rgk_curve = bsa_curve(phi_grid, Q2m, xBm, tm, Ebm,
                              debug=(args.debug and ibin==1), tag=f"bin{ibin}-RGK")

        # Optional ImH-only overlay using fit params
        imh_curve = None
        if args.CFFs == 1:
            flags_imh = dict(H=flags_fit.get("H",1), Ht=0, E=0, Et=0)
            push_globals(param_map_fit, flags_imh, label="RGK-ImH")
            imh_curve = bsa_curve(phi_grid, Q2m, xBm, tm, Ebm)

        # Optional BH-limit
        A_bh = None
        if args.bh_limit:
            A_bh = bh_limit_sinusoid(Q2m, xBm, tm, Ebm)

        # --- Plot ---
        fig, ax = plt.subplots(figsize=(8,5))
        ax.errorbar(phi_data, A_data, yerr=dA, fmt='o', ms=5, color='k', label='Data')

        ax.plot(phi_grid, bkm_curve, '-',  lw=2, label='BKM')
        ax.plot(phi_grid, rgk_curve, '--', lw=2, label='RGK 6 & 7 GeV Fit')

        if imh_curve is not None:
            ax.plot(phi_grid, imh_curve, '-.', lw=2, label='RGK (ImH only)')

        if A_bh is not None:
            ax.plot(phi_grid, A_bh(phi_grid), ':', lw=2, label=r'BH limit: $s_1^I/c_0^{BH}\,\sin\phi$')

        ax.set_xlim(0, 360)
        ax.set_xticks([0,60,120,180,240,300,360])
        ax.set_ylim(-0.6, 0.6)
        ax.set_xlabel(r'$\phi\;[\mathrm{deg}]$')
        ax.set_ylabel(r'$A_{LU}(\phi)$')
        ax.set_title(
            (r'$\langle Q^2\rangle={:.2f}\,\mathrm{{GeV}}^2,\;\langle x_B\rangle={:.3f},\;'
             r'\langle -t\rangle={:.3f}\,\mathrm{{GeV}}^2,\;\langle E_b\rangle={:.2f}\,\mathrm{{GeV}}$'
            ).format(Q2m, xBm, abs(tm), Ebm),
            pad=12
        )
        ax.legend(loc='upper right', frameon=True, edgecolor='k')
        plt.tight_layout()

        outname = os.path.join(args.outdir,
            f"BSA_bin{ibin:02d}_{ts}_Q2_{Q2m:.2f}_xB_{xBm:.3f}_mt_{abs(tm):.3f}.pdf")
        fig.savefig(outname)
        print(f">> Saved bin {ibin} plot to {outname}")
        plt.close(fig)


if __name__ == "__main__":
    main()