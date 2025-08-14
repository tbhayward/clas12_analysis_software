#!/usr/bin/env python3
"""
plot_8.5GeV_prediction.py
--------------------------------
Predict A_LU(φ) = [A*sinφ] / [1 + B cosφ + C cos2φ] at a chosen beam energy
(default 8.5 GeV), using the SAME (Q2, xB, t) binning as an existing data file.

Curves shown per bin:
  • BKM (ALU):  A = (R * s1_I^BKM) / c0_BH,    with R=(xB/y)*(1+eps^2)^2
  • RGK (ALU):  A = (s1_I^fit)    / c0_BH      (R absorbed in Im–CFF fit)

There is NO data plotted here; we only draw the two theory predictions using
the bin-mean kinematics from the data file but replacing Eb → Ebeam (e.g. 8.5 GeV).

Usage:
  python plot_8.5GeV_prediction.py output/fit_results/fit_results_<TS>.txt
         [--data PATH]                         # if omitted, taken from "input" line of fit file
         [--approx-order {1,2}]                # 1: use c0−c1 cosφ ; 2: also include c2 cos2φ
         [--Ebeam 8.5]                         # GeV
         [--outdir output/plots_8p5GeV]
         [--debug]
"""

import os, re, argparse
import numpy as np
import matplotlib.pyplot as plt
import ROOT  # PyROOT


# ─────────────────────────── C++ loader ────────────────────────────
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
def _parse_range_line(line, key_a, key_b):
    toks = line.split()
    out = {}
    try:
        if len(toks) >= 4:
            va = toks[1]; vb = toks[3]
            out[key_a] = float(va) if va.upper() != "NA" else None
            out[key_b] = float(vb) if vb.upper() != "NA" else None
    except ValueError:
        pass
    return out

def parse_fit_results(fname):
    with open(fname) as f:
        lines = [l.strip() for l in f if l.strip()]

    # flags like: "H 1 Ht 0 E 0 Et 0"
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

    # original data path (optional)
    input_path = next((l.split(None, 1)[1].strip()
                       for l in lines if l.startswith("input")), None)

    # xi and -t ranges
    ranges = {'xi_min':None, 'xi_max':None, 'mt_min':None, 'mt_max':None}
    for l in lines:
        if l.startswith("xi_min"):
            ranges.update(_parse_range_line(l, "xi_min", "xi_max"))
        elif l.startswith("-t_min"):
            ranges.update(_parse_range_line(l, "mt_min", "mt_max"))

    return flags, param_map_fit, input_path, ranges


# ───────────────────── data loader (φ-binned) ───────────────────────
def load_bins(datafile):
    """
    Read φ-sorted points and break into bins at φ wraparound.
    We only need the mean (Q2,xB,t,Eb) per bin; φ array is for plotting.
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


# ───────────── helpers: BH coeffs, s1_I, R ─────────────
def bh_coeffs(Q2, xB, t, Eb):
    """Return BH c0,c1,c2 (independent of φ) for given kinematics & beam energy."""
    dv = ROOT.BMK_DVCS(-1, +1, 0, Eb, xB, Q2, t, 0.0)
    return float(dv.c0_BH()), float(dv.c1_BH()), float(dv.c2_BH())

def s1I_and_K(Q2, xB, t, Eb):
    """Return (s1_I, R, y, eps2) at the given kinematics & beam energy, with current globals."""
    dv = ROOT.BMK_DVCS(-1, +1, 0, Eb, xB, Q2, t, 0.0)
    s1I  = float(dv.s1_I())
    y    = float(dv.y)
    eps2 = float(dv.eps2)
    R    = float((xB / y) * (1.0 + eps2)**2)
    return s1I, R, y, eps2


# ─────────────── ALU curve builder (order 1 or 2) ────────────────
def alu_curve(phi_deg, A1, B1, C1=0.0, order=1):
    ph = np.deg2rad(phi_deg)
    if order == 2:
        denom = 1.0 + B1*np.cos(ph) + C1*np.cos(2*ph)
    else:
        denom = 1.0 + B1*np.cos(ph)
    eps = 1e-14
    denom = np.where(np.abs(denom) < eps, np.sign(denom)*eps, denom)
    return A1 * np.sin(ph) / denom


# ───────────────────────────── main ───────────────────────────────
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('fitfile', help='fit_results_<TIMESTAMP>.txt from your fitter')
    ap.add_argument('--data', default=None, help='Path to an existing BSA φ-binned file to define bins')
    ap.add_argument('--approx-order', type=int, choices=[1,2], default=1,
                    help='Denominator order for theory: 1→use c0−c1 cosφ ; 2→also include c2 cos2φ')
    ap.add_argument('--Ebeam', type=float, default=8.5, help='Prediction beam energy in GeV (default 8.5)')
    ap.add_argument('--outdir', default='output/plots_8p5GeV', help='Output directory for PDFs')
    ap.add_argument('--debug', action='store_true', help='Print detailed diagnostics for the first bin')
    args = ap.parse_args()

    # Parse fit results and ranges
    flags_fit, param_map_fit, input_path, ranges = parse_fit_results(args.fitfile)

    # Use provided datafile, else the one recorded in the fit file, else a default
    datafile = (args.data or input_path or 'imports/rgk_rga_combined.txt')
    print(">> Using bin-defining data file:", datafile)
    print(">> Prediction beam energy (GeV):", args.Ebeam)
    print(">> Fit flags:", flags_fit)

    # Load C++
    _load_dvcs_once()

    # Load bins (we only use means <Q2>,<xB>,<t>; Eb from file is ignored)
    bins = load_bins(datafile)
    print(f">> Found {len(bins)} φ-bins")

    os.makedirs(args.outdir, exist_ok=True)

    # Baseline BKM input (Im H, H~, E all on; Et off)
    flags_bkm = {'H':1,'Ht':1,'E':1,'Et':0}
    params_bkm = bkm_defaults()

    # Timestamp tag from fitfile name (if present)
    m = re.search(r'fit_results_(\d{8}_\d{6})\.txt$', os.path.basename(args.fitfile))
    ts = (m.group(1) if m else "noTS")

    # Loop over bins
    for ibin, b in enumerate(bins, start=1):
        phi_grid = np.linspace(0.0, 360.0, 361)

        Q2m, xBm, tm = b["Q2m"], b["xBm"], b["tm"]
        Eb_pred = float(args.Ebeam)

        xi_m = xBm / (2.0 - xBm)
        mt_m = abs(tm)

        # BH coefficients at the PREDICTION beam energy
        c0_bh, c1_bh, c2_bh = bh_coeffs(Q2m, xBm, tm, Eb_pred)
        B_bh = -(c1_bh / c0_bh)
        C_bh =  (c2_bh / c0_bh) if args.approx_order == 2 else 0.0

        # --- BKM curve (baseline) ---
        push_globals(params_bkm, flags_bkm, label="BKM")
        s1I_bkm, R_bkm, *_ = s1I_and_K(Q2m, xBm, tm, Eb_pred)
        A1_bkm = (R_bkm * s1I_bkm) / c0_bh
        bkm_curve = alu_curve(phi_grid, A1_bkm, B_bh, C_bh, order=args.approx_order)

        # --- RGK fit curve (your Im–CFFs & flags) ---
        push_globals(param_map_fit, flags_fit, label="RGK-fit")
        s1I_rgk, R_rgk, *_ = s1I_and_K(Q2m, xBm, tm, Eb_pred)
        A1_rgk = (s1I_rgk) / c0_bh
        rgk_curve = alu_curve(phi_grid, A1_rgk, B_bh, C_bh, order=args.approx_order)

        # Range flag
        outside = False
        if ranges.get('xi_min') is not None and ranges.get('xi_max') is not None:
            if not (ranges['xi_min'] <= xi_m <= ranges['xi_max']):
                outside = True
        if ranges.get('mt_min') is not None and ranges.get('mt_max') is not None:
            if not (ranges['mt_min'] <= mt_m <= ranges['mt_max']):
                outside = True

        # Optional diagnostics for the first bin
        if args.debug and ibin == 1:
            print(f"\n=== Diagnostics for bin {ibin:02d} (prediction) ===")
            print(f"  <Q2>= {Q2m:6.3f} GeV^2, <xB>= {xBm:0.3f}, <xi>= {xi_m:0.3f}, <-t>= {mt_m:0.3f} GeV^2,  Eb(pred)={Eb_pred:0.2f} GeV")
            print(f"  BH (pred): B = -c1/c0 = {B_bh:.6f},  C = c2/c0 = {C_bh:.6f}")
            print(f"  R_BKM={R_bkm:.5f}, R_RGK={R_rgk:.5f}, s1I(BKM)={s1I_bkm:.6e}, s1I(RGK)={s1I_rgk:.6e}")
            print(f"  A1_BKM={A1_bkm:.6e},  A1_RGK={A1_rgk:.6e}")
            # restore RGK globals if needed later
            push_globals(param_map_fit, flags_fit, label="RGK-fit-restore")

        # --- Plot (theory only) ---
        fig, ax = plt.subplots(figsize=(8,5))
        ax.plot(phi_grid, bkm_curve, '-',  lw=2, color='tab:blue',
                label=f'BKM (ALU): A={A1_bkm:.3f}')
        ax.plot(phi_grid, rgk_curve, '-',  lw=2, color='tab:orange',
                label=f'RGK fit (ALU): A={A1_rgk:.3f}')

        ax.set_xlim(0, 360)
        ax.set_xticks([0,60,120,180,240,300,360])
        ax.set_ylim(-0.6, 0.6)
        ax.set_xlabel(r'$\phi\;(\mathrm{deg})$')
        ax.set_ylabel(r'$A_{LU}(\phi)$')
        ax.set_title(
            (r'Prediction @ $E_b={:.1f}\,$GeV  |  '
             r'$\langle Q^2\rangle={:.2f}\,\mathrm{{GeV}}^2,\;'
             r'\langle x_B\rangle={:.3f},\;\langle \xi\rangle={:.3f},\;'
             r'\langle -t\rangle={:.3f}\,\mathrm{{GeV}}^2$'
            ).format(Eb_pred, Q2m, xBm, xi_m, mt_m),
            pad=12
        )

        ax.legend(loc='upper right', frameon=True, edgecolor='k')

        if outside:
            ax.text(0.02, 0.05, "bin outside of fit range",
                    transform=ax.transAxes, color='crimson',
                    fontsize=11, fontweight='bold')

        plt.tight_layout()

        outname = os.path.join(
            args.outdir,
            f"ALU_pred8p5_bin{ibin:02d}_{ts}_Q2_{Q2m:.2f}_xB_{xBm:.3f}_xi_{xi_m:.3f}_mt_{mt_m:.3f}_ord{args.approx_order}.pdf"
        )
        fig.savefig(outname)
        print(f">> Saved prediction for bin {ibin} to {outname}")
        plt.close(fig)


if __name__ == "__main__":
    main()