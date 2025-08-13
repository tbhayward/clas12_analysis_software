#!/usr/bin/env python3
"""
plot_BSA_CFFs_from_fit.py

Plots A_LU(φ) per (Q2, xB, t, Eb) bin for:
  • BKM (full) — solid blue
  • RGK fit for Im CFF(s) — solid orange:
        A(φ) ≈ [s1_I · sinφ] / [c0_BH + c1_BH cosφ (+ c2_BH cos2φ if --approx-order=2)]
        where s1_I uses your fitted ImCFFs; denominator is BH-only.
  • BKM with Im CFF(s) from fit — dashed green:
        Full BSA like BKM, but replace any available *imaginary-part* CFF parameters
        (r, n, alpha0, alpha1, b, M2, P) with the fitted values from your results file.
        All else (real parts, subtraction constants) stays at BKM defaults.

Usage:
  python plot_BSA_CFFs_from_fit.py output/fit_results/fit_results_<TIMESTAMP>.txt
         [--data PATH]
         [--approx-order {1,2}]    # denom: c0+c1 cosφ (1) or c0+c1 cosφ+c2 cos2φ (2)
         [--B-source {bh,data}]    # use BH harmonics (default) or fit B1 from data (only for order=1)
         [--outdir output/plots]
         [--debug]

Notes
-----
• The orange RGK curve depends only on ImCFFs via s1_I; the denominator is pure BH.
• The green curve shows what the *full* BKM prediction would look like if you swap in your fitted
  ImCFF parameters (for any CFFs present in the fit file), leaving the rest untouched.
• We keep Ē disabled as a safe default to avoid M2=0 / P≠0 pathologies.
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

    # optional input
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
    """Return (s1_I, c0, c1, c2) for fixed kinematics (independent of φ)."""
    dv = ROOT.BMK_DVCS(-1, +1, 0, Eb, xB, Q2, t, 0.0)
    s1I = dv.s1_I()
    c0  = dv.c0_BH()
    c1  = dv.c1_BH()
    c2  = dv.c2_BH()
    return s1I, c0, c1, c2

def fit_B_from_data(phi_data, A_data, dA):
    """Fit A_LU(φ) = [A1 sinφ] / [1 + B1 cosφ]; return B1, A1, χ²/ndf."""
    import numpy.linalg as LA
    ph = np.deg2rad(phi_data)
    w  = 1.0 / np.maximum(dA, 1e-8)**2
    # linearize: A*(1 + B cosφ) ≈ A1 sinφ  →  A ≈ a sinφ + b sinφ cosφ ; then B ≈ b/a
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


# ─────────── build ImCFF-substitution map for the green curve ───────────
def build_imfit_param_map(fit_param_map):
    """
    Start from BKM defaults and *only* overwrite imaginary-part ansatz params
    (r, n, alpha0, alpha1, b, M2, P) for any CFFs present in the fit file.
    Real-part parameters (renormReal, C0, MD2, lambda) are left at defaults.
    """
    pm = bkm_defaults().copy()
    for cff in ("H", "Ht", "E", "Et"):
        for k in ("r","n","alpha0","alpha1","b","M2","P"):
            key = f"{k}_{cff}"
            if key in fit_param_map:
                pm[key] = float(fit_param_map[key])
    # Keep renormImag/Real as defaults (1.0). Do NOT import C0/MD2/lambda.
    return pm


# ───────────────────────────── labeling helpers ─────────────────────────────
def latex_im_label(flags_fit):
    parts = []
    if flags_fit.get("H",0):  parts.append(r"$\mathrm{Im}\,\mathcal{H}$")
    if flags_fit.get("Ht",0): parts.append(r"$\mathrm{Im}\,\widetilde{\mathcal{H}}$")
    if flags_fit.get("E",0):  parts.append(r"$\mathrm{Im}\,\mathcal{E}$")
    if flags_fit.get("Et",0): parts.append(r"$\mathrm{Im}\,\widetilde{\mathcal{E}}$")
    if not parts:
        return r"Im CFF(s)"  # fallback
    if len(parts) == 1:
        return parts[0]
    return r",\;".join(parts)


# ───────────────────────────── main ───────────────────────────────
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('fitfile', help='fit_results_<TIMESTAMP>.txt from your fitter')
    ap.add_argument('--data', default=None, help='Override BSA data file; else use "input" in fit file or fallback')
    ap.add_argument('--approx-order', type=int, choices=[1,2], default=2,
                    help='Use c0+c1 cosφ (1) or c0+c1 cosφ+c2 cos2φ (2) in the BH denominator')
    ap.add_argument('--B-source', choices=['bh','data'], default='bh',
                    help='Use BH Fourier coeffs for B (default) or fit B1 from data (only valid for order=1)')
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

    # BKM baseline flags: H,Ht,E on; Et off (safe)
    flags_bkm = {'H':1,'Ht':1,'E':1,'Et':0}
    params_bkm = bkm_defaults()

    # For filename timestamp
    m = re.search(r'fit_results_(\d{8}_\d{6})\.txt$', os.path.basename(args.fitfile))
    ts = (m.group(1) if m else "noTS")

    # Pretty label for which ImCFFs were fitted
    imlabel = latex_im_label(flags_fit)

    # Loop bins
    for ibin, b in enumerate(bins, start=1):
        phi_data, A_data, dA = b["phi"], b["A"], b["sigA"]
        Q2m, xBm, tm, Ebm = b["Q2m"], b["xBm"], b["tm"], b["Ebm"]

        # Dense grid
        phi_grid = np.linspace(0.0, 360.0, 361)

        # --- BKM baseline (full) ---
        push_globals(params_bkm, flags_bkm, label="BKM")
        if args.debug and ibin == 1:
            ROOT.gInterpreter.ProcessLine(r'std::cout<<"[BKM] hasH="<<hasH<<" hasHt="<<hasHt<<" hasE="<<hasE<<" hasEt="<<hasEt<<std::endl;')
        bkm_full = bsa_full_curve(phi_grid, Q2m, xBm, tm, Ebm,
                                  debug=(args.debug and ibin==1), tag=f"bin{ibin}-BKM")

        # --- BKM full with ImCFFs substituted from fit (green dashed) ---
        params_bkm_imfit = build_imfit_param_map(param_map_fit)
        push_globals(params_bkm_imfit, flags_bkm, label="BKM+ImFit")
        if args.debug and ibin == 1:
            ROOT.gInterpreter.ProcessLine(r'std::cout<<"[BKM+ImFit] hasH="<<hasH<<" hasHt="<<hasHt<<" hasE="<<hasE<<" hasEt="<<hasEt<<std::endl;')
        bkm_with_imfit = bsa_full_curve(phi_grid, Q2m, xBm, tm, Ebm,
                                        debug=False, tag=f"bin{ibin}-BKM-ImFit")

        # --- RGK fit approx (use your fit params; BH denominator choice per args) ---
        push_globals(param_map_fit, flags_fit, label="RGK-fit")
        if args.debug and ibin == 1:
            ROOT.gInterpreter.ProcessLine(r'std::cout<<"[RGK] hasH="<<hasH<<" hasHt="<<hasHt<<" hasE="<<hasE<<" hasEt="<<hasEt<<std::endl;')

        s1I_rgk, c0_rgk, c1_rgk, c2_rgk = bh_coeffs_and_s1I(Q2m, xBm, tm, Ebm)

        # Choose B from BH (default) or fit-from-data (only for order=1)
        if args.approx_order == 1 and args.B_source == 'data':
            B1_fit, _, _ = fit_B_from_data(phi_data, A_data, dA)
            B1_used, B2_used = B1_fit, 0.0
        else:
            B1_used = c1_rgk / c0_rgk if c0_rgk != 0 else 0.0
            B2_used = (c2_rgk / c0_rgk) if (args.approx_order == 2 and c0_rgk != 0) else 0.0

        ph = np.deg2rad(phi_grid)
        denom = c0_rgk * (1.0 + B1_used*np.cos(ph) + (B2_used*np.cos(2*ph) if args.approx_order==2 else 0.0))
        eps = 1e-14
        denom = np.where(np.abs(denom) < eps, np.sign(denom)*eps, denom)
        rgk_curve = (s1I_rgk * np.sin(ph)) / denom

        # --- Plot ---
        fig, ax = plt.subplots(figsize=(8,5))
        ax.errorbar(phi_data, A_data, yerr=dA, fmt='o', ms=5, color='k', label='Data')
        ax.plot(phi_grid, bkm_full, '-',   lw=2, color='tab:blue',   label='BKM')
        ax.plot(phi_grid, rgk_curve, '-',  lw=2, color='tab:orange', label=f'RGK fit for {imlabel}')
        ax.plot(phi_grid, bkm_with_imfit, '--', lw=2, color='tab:green', label=f'BKM with {imlabel} from fit')

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