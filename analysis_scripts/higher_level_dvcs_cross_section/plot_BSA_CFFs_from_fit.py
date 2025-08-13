#!/usr/bin/env python3
"""
plot_BSA_CFFs_from_fit.py  —  ALU(φ) with BH-only denominator for theory,
and (A,B,C) fitted to data.

Theory curves:
  • BKM (ALU):  A = (R * s1_I^BKM) / c0_BH      with R=(xB/y)*(1+eps^2)^2
  • RGK (ALU):  A = (s1_I^fit)   / c0_BH        (R absorbed in the Im–CFF fit)

Data curve:
  • Fit A,B,(C) in  A*sinφ / (1 + B cosφ [+ C cos2φ])  by weighted least squares.

Usage:
  python plot_BSA_CFFs_from_fit.py output/fit_results/fit_results_<TS>.txt
         [--data PATH]
         [--approx-order {1,2}]    # 1: fit A,B ; 2: fit A,B,C
         [--outdir output/plots]
         [--debug]
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

    flags_line = next((l for l in lines
                       if re.match(r'^\s*H\s+[01]\s+Ht\s+[01]\s+E\s+[01]\s+Et\s+[01]\s*$', l)), None)
    if flags_line is None:
        flags_line = next((l for l in lines if l.startswith("H ")), None)
    if flags_line is None:
        raise RuntimeError("No 'H ... Ht ... E ... Et ...' flags line found in fit file")
    toks = flags_line.split()
    flags = {toks[i]: int(toks[i+1]) for i in range(0, len(toks), 2)}

    pnames, pvals = [], None
    for i, l in enumerate(lines):
        if l.startswith("# parameters"):
            pnames = l.split()[2:]
        elif l.startswith("# values") and i+1 < len(lines):
            pvals = list(map(float, lines[i+1].split()))
    if not pnames or pvals is None:
        raise RuntimeError("Couldn't parse parameter names/values from fit file")
    param_map_fit = dict(zip(pnames, pvals))

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
        # Et — keep harmless
        'r_Et': 1.0, 'n_Et': 0.0, 'alpha0_Et': 0.0, 'alpha1_Et': 0.0, 'b_Et': 0.0, 'M2_Et': 1.0, 'P_Et': 0.0,
    }


# ───────── push globals into C++ (and set switches) ───────────────
def push_globals(param_map, flags, *, safety_disable_Et=True, label=""):
    ROOT.gInterpreter.ProcessLine("hasH=0; hasHt=0; hasE=0; hasEt=0;")

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


# ───────────── helpers: BH coeffs, s1_I, K, components ─────────────
def bh_coeffs(Q2, xB, t, Eb):
    dv = ROOT.BMK_DVCS(-1, +1, 0, Eb, xB, Q2, t, 0.0)
    return float(dv.c0_BH()), float(dv.c1_BH()), float(dv.c2_BH())

def s1I_and_K(Q2, xB, t, Eb):
    dv = ROOT.BMK_DVCS(-1, +1, 0, Eb, xB, Q2, t, 0.0)
    s1I  = float(dv.s1_I())
    y    = float(dv.y)
    eps2 = float(dv.eps2)
    R    = float((xB / y) * (1.0 + eps2)**2)
    return s1I, R, y, eps2

def s1I_components(Q2, xB, t, Eb, *, param_map, flags):
    out = {'H':0.0, 'Ht':0.0, 'E':0.0}
    for cff in ('H','Ht','E'):
        if not int(flags.get(cff,0)):
            out[cff] = 0.0
            continue
        fl = {'H':0,'Ht':0,'E':0,'Et':0}
        fl[cff] = 1
        push_globals(param_map, fl, label=f"comp-{cff}")
        dv = ROOT.BMK_DVCS(-1, +1, 0, Eb, xB, Q2, t, 0.0)
        out[cff] = float(dv.s1_I())
    return out


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


# ─────────── weighted fit of (A,B[,C]) to data using ROOT ──────────
def fit_data_ABC(phi_deg, y, yerr, *, order, A0, B0, C0):
    n = len(phi_deg)
    gr = ROOT.TGraphErrors(n)
    for i in range(n):
        gr.SetPoint(i, float(phi_deg[i]), float(y[i]))
        gr.SetPointError(i, 0.0, float(yerr[i]))

    if order == 2:
        formula = ("[0]*sin(x*TMath::Pi()/180.)"
                   "/(1 + [1]*cos(x*TMath::Pi()/180.) + [2]*cos(2*x*TMath::Pi()/180.))")
    else:
        formula = ("[0]*sin(x*TMath::Pi()/180.)"
                   "/(1 + [1]*cos(x*TMath::Pi()/180.))")
    f = ROOT.TF1("f_data_fit", formula, 0, 360)
    f.SetParameters(float(A0), float(B0), float(C0))
    if order == 1:
        f.FixParameter(2, 0.0)  # no cos2φ term

    gr.Fit(f, "Q0R")  # quiet, no draw, respect range

    A  = float(f.GetParameter(0)); dA = float(f.GetParError(0))
    B  = float(f.GetParameter(1)); dB = float(f.GetParError(1))
    C  = float(f.GetParameter(2)) if order == 2 else 0.0
    dC = float(f.GetParError(2)) if order == 2 else 0.0
    chi2 = float(f.GetChisquare()); ndf = int(f.GetNDF())
    return dict(A=A, dA=dA, B=B, dB=dB, C=C, dC=dC, chi2=chi2, ndf=ndf)


# ───────────────────────────── main ───────────────────────────────
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('fitfile', help='fit_results_<TIMESTAMP>.txt from your fitter')
    ap.add_argument('--data', default=None, help='Override BSA data file; else use "input" in fit file or fallback')
    ap.add_argument('--approx-order', type=int, choices=[1,2], default=1,
                    help='Denominator order for theory AND for data fit: 1→A,B ; 2→A,B,C')
    ap.add_argument('--outdir', default='output/plots', help='Output directory for PDFs')
    ap.add_argument('--debug', action='store_true', help='Print detailed diagnostics for the first bin')
    args = ap.parse_args()

    # Parse fit file and choose data file
    flags_fit, param_map_fit, input_path = parse_fit_results(args.fitfile)
    datafile = (args.data or input_path or 'imports/rgk_preliminary_bsa.txt')
    print(">> Using data file:", datafile)
    print(">> Fit flags:", flags_fit)

    # Load C++
    _load_dvcs_once()

    # Load φ-bins
    bins = load_bins(datafile)
    print(f">> Found {len(bins)} φ-bins")

    os.makedirs(args.outdir, exist_ok=True)

    # BKM baseline flags: H,Ht,E on; Et off (safe)
    flags_bkm = {'H':1,'Ht':1,'E':1,'Et':0}
    params_bkm = bkm_defaults()

    # Timestamp tag from fitfile name
    m = re.search(r'fit_results_(\d{8}_\d{6})\.txt$', os.path.basename(args.fitfile))
    ts = (m.group(1) if m else "noTS")

    # Loop over bins
    for ibin, b in enumerate(bins, start=1):
        phi_data, A_data, dA = b["phi"], b["A"], b["sigA"]
        Q2m, xBm, tm, Ebm = b["Q2m"], b["xBm"], b["tm"], b["Ebm"]

        # Dense φ for smooth curves
        phi_grid = np.linspace(0.0, 360.0, 361)

        # --- BH coefficients (for theory curves) ---
        c0_bh, c1_bh, c2_bh = bh_coeffs(Q2m, xBm, tm, Ebm)
        B_bh = -(c1_bh / c0_bh)
        C_bh =  (c2_bh / c0_bh) if args.approx_order == 2 else 0.0

        # --- BKM theory (ALU) ---
        push_globals(params_bkm, flags_bkm, label="BKM")
        s1I_bkm, R_bkm, *_ = s1I_and_K(Q2m, xBm, tm, Ebm)
        A1_bkm = (R_bkm * s1I_bkm) / c0_bh
        bkm_curve = alu_curve(phi_grid, A1_bkm, B_bh, C_bh, order=args.approx_order)

        # --- RGK theory (ALU) ---
        push_globals(param_map_fit, flags_fit, label="RGK-fit")
        s1I_rgk, R_rgk, *_ = s1I_and_K(Q2m, xBm, tm, Ebm)
        A1_rgk = (s1I_rgk) / c0_bh
        rgk_curve = alu_curve(phi_grid, A1_rgk, B_bh, C_bh, order=args.approx_order)

        # --- Data fit: free (A,B[,C]) with BH values as starting guesses ---
        A0_guess = float(np.nan_to_num(np.average(A_data, weights=1/np.maximum(dA,1e-9)**2)))
        fit_res = fit_data_ABC(
            phi_data, A_data, dA,
            order=args.approx_order,
            A0=A0_guess if np.isfinite(A0_guess) else A1_bkm,
            B0=B_bh, C0=C_bh
        )
        data_curve = alu_curve(phi_grid, fit_res["A"], fit_res["B"], fit_res["C"], order=args.approx_order)

        # ── Diagnostics (first bin) ──
        if args.debug and ibin == 1:
            print(f"\n=== Diagnostics for bin {ibin:02d} ===")
            print(f"  <Q2>= {Q2m:6.3f} GeV^2, <xB>= {xBm:0.3f}, <-t>= {abs(tm):0.3f} GeV^2, <Eb>= {Ebm:0.2f} GeV")
            print("  (1) K-factor & s1_I ratio:")
            ratio = (s1I_rgk / s1I_bkm) if abs(s1I_bkm) > 0 else float('nan')
            print(f"      R_BKM = {R_bkm:.5f},  R_RGK = {R_rgk:.5f}")
            print(f"      s1_I(BKM) = {s1I_bkm:.6e},  s1_I(RGK) = {s1I_rgk:.6e},  ratio RGK/BKM = {ratio:.6f}")
            print(f"      A1_BKM = {A1_bkm:.6e},  A1_RGK = {A1_rgk:.6e}")
            print("  (2) BH B,C used for theory:")
            print(f"      B_BH = -c1/c0 = {B_bh:.6f},  C_BH = c2/c0 = {C_bh:.6f}")
            comp_bkm = s1I_components(Q2m, xBm, tm, Ebm, param_map=params_bkm,    flags=flags_bkm)
            comp_rgk = s1I_components(Q2m, xBm, tm, Ebm, param_map=param_map_fit, flags=flags_fit)
            tot_bkm = comp_bkm['H'] + comp_bkm['Ht'] + comp_bkm['E']
            tot_rgk = comp_rgk['H'] + comp_rgk['Ht'] + comp_rgk['E']
            print("  (3) s1_I components (Im H / Im H~ / Im E):")
            print(f"      BKM: H={comp_bkm['H']:.6e}, H~={comp_bkm['Ht']:.6e}, E={comp_bkm['E']:.6e}  -> total={tot_bkm:.6e}")
            print(f"      RGK: H={comp_rgk['H']:.6e}, H~={comp_rgk['Ht']:.6e}, E={comp_rgk['E']:.6e}  -> total={tot_rgk:.6e}")
            print("  (4) Data fit results (free A,B[,C]):")
            print(f"      A = {fit_res['A']:.6f} ± {fit_res['dA']:.6f}")
            print(f"      B = {fit_res['B']:.6f} ± {fit_res['dB']:.6f}")
            if args.approx_order == 2:
                print(f"      C = {fit_res['C']:.6f} ± {fit_res['dC']:.6f}")
            print(f"      χ²/ndf = {fit_res['chi2']:.2f}/{fit_res['ndf']} = "
                  f"{(fit_res['chi2']/fit_res['ndf']) if fit_res['ndf']>0 else float('nan'):.2f}")
            # restore RGK globals after diagnostics component toggles
            push_globals(param_map_fit, flags_fit, label="RGK-fit-restore")

        # --- Plot ---
        fig, ax = plt.subplots(figsize=(8,5))
        ax.errorbar(phi_data, A_data, yerr=dA, fmt='o', ms=5, color='k', label='Data')

        ax.plot(phi_grid, bkm_curve, '-',  lw=2, color='tab:blue',
                label=f'BKM (ALU): A={A1_bkm:.3f}')
        ax.plot(phi_grid, rgk_curve, '-',  lw=2, color='tab:orange',
                label=f'RGK fit (ALU): A={A1_rgk:.3f}')
        ax.plot(phi_grid, data_curve, '--', lw=2, color='tab:green',
                label=(f'Data fit: A={fit_res["A"]:.3f}±{fit_res["dA"]:.3f}'))

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
            f"ALU_bin{ibin:02d}_{ts}_Q2_{Q2m:.2f}_xB_{xBm:.3f}_mt_{abs(tm):.3f}_ord{args.approx_order}.pdf"
        )
        fig.savefig(outname)
        print(f">> Saved bin {ibin} plot to {outname}")
        plt.close(fig)


if __name__ == "__main__":
    main()