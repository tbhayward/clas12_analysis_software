#!/usr/bin/env python3
"""
plot_BSA_CFFs_from_fit.py

Plots A_LU(φ) per (Q2, xB, t, Eb) bin for:
  • BKM (full) — solid blue
  • RGK fit for Im CFF(s) — solid orange:
        A(φ) ≈ [s1_I · sinφ] / [c0_BH + c1_BH cosφ (+ c2_BH cos2φ if --approx-order=2)]
        (s1_I from your fitted Im-CFFs; denominator is BH-only)
  • BKM with Im CFF(s) from fit — dashed green:
        --green-mode=num-only (default): numerator from fitted Im-CFFs, denominator from default BKM
        --green-mode=full: both numerator & denominator from the hybrid params (Im from fit)

Usage:
  python plot_BSA_CFFs_from_fit.py output/fit_results/fit_results_<TIMESTAMP>.txt
         [--data PATH]
         [--approx-order {1,2}]
         [--B-source {bh,data}]
         [--green-mode {num-only,full}]
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
        'r_H': 0.9, 'n_H': 1.35, 'alpha0_H': 0.43, 'alpha1_H': 0.85, 'b_H': 0.4, 'M2_H': 0.64, 'P_H': 1.0,
        'r_Ht': 7.0, 'n_Ht': 0.6,  'alpha0_Ht': 0.43, 'alpha1_Ht': 0.85, 'b_Ht': 2.0, 'M2_Ht': 0.8,  'P_Ht': 1.0,
        'r_E': 0.9, 'n_E': 1.35, 'alpha0_E': 0.43, 'alpha1_E': 0.85, 'b_E': 0.4, 'M2_E': 0.64, 'P_E': 1.0,
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


# ─────────────── helpers to compute curves/coeffs ────────────────
def bsa_full_curve(phi_deg, Q2, xB, t, Eb, *, debug=False, tag="", max_dbg_pts=3):
    out = np.empty_like(phi_deg, float)
    for i, ph in enumerate(phi_deg):
        dv = ROOT.BMK_DVCS(-1, +1, 0, Eb, xB, Q2, t, ph)
        out[i] = dv.BSA()
        if debug and i < max_dbg_pts:
            print(f"[{tag}] φ={ph:6.1f}°  BSA={out[i]:+.6e}  "
                  f"BH2={dv.BH2():.3e}  DVCS2={dv.DVCS2():.3e}  I={dv.BHDVCS():.3e}  s1_I={dv.s1_I():+.6e}")
    return out

def bh_coeffs_and_s1I(Q2, xB, t, Eb):
    dv = ROOT.BMK_DVCS(-1, +1, 0, Eb, xB, Q2, t, 0.0)
    return dv.s1_I(), dv.c0_BH(), dv.c1_BH(), dv.c2_BH()

def xs_num_den_over_phi(param_map, flags, phi_deg, Q2, xB, t, Eb):
    """Return arrays: numerator=(σ+−σ−)/2 and denominator=(σ++σ−)/2 for given parameters."""
    push_globals(param_map, flags, label="xs_num_den")
    num = np.empty_like(phi_deg, float)
    den = np.empty_like(phi_deg, float)
    for i, ph in enumerate(phi_deg):
        dv_plus  = ROOT.BMK_DVCS(-1, +1, 0, Eb, xB, Q2, t, float(ph))
        dv_minus = ROOT.BMK_DVCS(-1, -1, 0, Eb, xB, Q2, t, float(ph))
        xs_p = dv_plus.CrossSection()
        xs_m = dv_minus.CrossSection()
        num[i] = 0.5*(xs_p - xs_m)
        den[i] = 0.5*(xs_p + xs_m)
    return num, den


def fit_B_from_data(phi_data, A_data, dA):
    import numpy.linalg as LA
    ph = np.deg2rad(phi_data)
    w  = 1.0 / np.maximum(dA, 1e-8)**2
    X = np.vstack([np.sin(ph), np.sin(ph)*np.cos(ph)]).T
    W = np.diag(w)
    y = A_data
    beta = LA.pinv(X.T @ W @ X) @ (X.T @ W @ y)
    a, b = beta
    B1 = b / max(a, 1e-12)
    A1 = a
    yfit = (a*np.sin(ph) + b*np.sin(ph)*np.cos(ph))
    chi2 = float(((y - yfit)**2 * w).sum())
    ndf  = max(len(y) - 2, 1)
    return B1, A1, chi2/ndf


def build_imfit_param_map(fit_param_map):
    """Hybrid map: start with BKM defaults, override ONLY Im-CFF ansatz params present in the fit."""
    pm = bkm_defaults().copy()
    for cff in ("H","Ht","E","Et"):
        for k in ("r","n","alpha0","alpha1","b","M2","P"):
            key = f"{k}_{cff}"
            if key in fit_param_map:
                pm[key] = float(fit_param_map[key])
    return pm


def latex_im_label(flags_fit):
    parts = []
    if flags_fit.get("H",0):  parts.append(r"$\mathrm{Im}\,\mathcal{H}$")
    if flags_fit.get("Ht",0): parts.append(r"$\mathrm{Im}\,\widetilde{\mathcal{H}}$")
    if flags_fit.get("E",0):  parts.append(r"$\mathrm{Im}\,\mathcal{E}$")
    if flags_fit.get("Et",0): parts.append(r"$\mathrm{Im}\,\widetilde{\mathcal{E}}$")
    return (parts[0] if len(parts)==1 else (r",\;".join(parts) if parts else r"Im CFF(s)"))


# ───────────────────────────── main ───────────────────────────────
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('fitfile', help='fit_results_<TIMESTAMP>.txt from your fitter')
    ap.add_argument('--data', default=None, help='Override BSA data file; else use "input" in fit file or fallback')
    ap.add_argument('--approx-order', type=int, choices=[1,2], default=2,
                    help='Use c0+c1 cosφ (1) or c0+c1 cosφ+c2 cos2φ (2) in the BH denominator (orange curve)')
    ap.add_argument('--B-source', choices=['bh','data'], default='bh',
                    help='For orange: use BH Fourier coeffs for B (default) or fit B1 from data (order=1 only)')
    ap.add_argument('--green-mode', choices=['num-only','full'], default='num-only',
                    help='Green curve: "num-only" (default) uses numerator from Im-fit over BKM default denominator; "full" replaces Im everywhere.')
    ap.add_argument('--outdir', default='output/plots', help='Output directory for PDFs')
    ap.add_argument('--debug', action='store_true', help='Print internal diagnostics for first bin')
    args = ap.parse_args()

    flags_fit, param_map_fit, input_path = parse_fit_results(args.fitfile)
    datafile = (args.data or input_path or 'imports/rgk_preliminary_bsa.txt')
    print(">> Using data file:", datafile)
    print(">> Fit flags:", flags_fit)

    _load_dvcs_once()
    bins = load_bins(datafile)
    print(f">> Found {len(bins)} φ-bins")
    os.makedirs(args.outdir, exist_ok=True)

    flags_bkm = {'H':1,'Ht':1,'E':1,'Et':0}
    params_bkm = bkm_defaults()
    params_imfit = build_imfit_param_map(param_map_fit)
    imlabel = latex_im_label(flags_fit)

    m = re.search(r'fit_results_(\d{8}_\d{6})\.txt$', os.path.basename(args.fitfile))
    ts = (m.group(1) if m else "noTS")

    for ibin, b in enumerate(bins, start=1):
        phi_data, A_data, dA = b["phi"], b["A"], b["sigA"]
        Q2m, xBm, tm, Ebm = b["Q2m"], b["xBm"], b["tm"], b["Ebm"]

        phi_grid = np.linspace(0.0, 360.0, 361)
        ph = np.deg2rad(phi_grid)

        # ── Blue: BKM full BSA
        push_globals(params_bkm, flags_bkm, label="BKM")
        bkm_full = bsa_full_curve(phi_grid, Q2m, xBm, tm, Ebm,
                                  debug=(args.debug and ibin==1), tag=f"bin{ibin}-BKM")

        # ── Orange: your previous approx using s1_I and BH-only denominator
        push_globals(param_map_fit, flags_fit, label="RGK-fit")
        s1I_rgk, c0_rgk, c1_rgk, c2_rgk = bh_coeffs_and_s1I(Q2m, xBm, tm, Ebm)
        if args.approx_order == 1 and args.B_source == 'data':
            # quick per-bin B1 fit
            import numpy.linalg as LA
            pp = np.deg2rad(phi_data)
            w  = 1.0/np.maximum(dA,1e-8)**2
            X = np.vstack([np.sin(pp), np.sin(pp)*np.cos(pp)]).T
            W = np.diag(w); y=A_data
            a,bcoef = (LA.pinv(X.T@W@X) @ (X.T@W@y))
            B1_used, B2_used = bcoef/max(a,1e-12), 0.0
        else:
            B1_used = (c1_rgk/c0_rgk) if c0_rgk!=0 else 0.0
            B2_used = (c2_rgk/c0_rgk) if (args.approx_order==2 and c0_rgk!=0) else 0.0
        denom_bh = c0_rgk*(1.0 + B1_used*np.cos(ph) + (B2_used*np.cos(2*ph) if args.approx_order==2 else 0.0))
        denom_bh = np.where(np.abs(denom_bh)<1e-14, np.sign(denom_bh)*1e-14, denom_bh)
        rgk_curve = (s1I_rgk * np.sin(ph)) / denom_bh

        # ── Green: consistent BSA from σ±
        # Denominator_DEFAULT from BKM
        _, den_bkm = xs_num_den_over_phi(params_bkm, flags_bkm, phi_grid, Q2m, xBm, tm, Ebm)

        # Numerator/Denominator with Im from fit
        num_imfit, den_imfit = xs_num_den_over_phi(params_imfit, flags_bkm, phi_grid, Q2m, xBm, tm, Ebm)

        if args.green_mode == 'full':
            den_for_green = den_imfit
        else:  # num-only
            den_for_green = den_bkm

        # protect against zeros
        den_for_green = np.where(np.abs(den_for_green)<1e-20,
                                 np.sign(den_for_green)*1e-20, den_for_green)
        green_curve = num_imfit / den_for_green

        # ── Plot
        fig, ax = plt.subplots(figsize=(8,5))
        ax.errorbar(phi_data, A_data, yerr=dA, fmt='o', ms=5, color='k', label='Data')
        ax.plot(phi_grid, bkm_full, '-',  lw=2, color='tab:blue',   label='BKM')
        ax.plot(phi_grid, rgk_curve, '-', lw=2, color='tab:orange', label=f'RGK fit for {imlabel}')
        ax.plot(phi_grid, green_curve, '--', lw=2, color='tab:green',
                label=(f'BKM with {imlabel} from fit'
                       + (' (full)' if args.green_mode=='full' else ' (num-only)')))

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