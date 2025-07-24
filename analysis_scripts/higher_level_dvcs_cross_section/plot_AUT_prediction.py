#!/usr/bin/env python3
"""
plot_BSA_and_AUT_from_fit.py

Usage:
    python plot_BSA_and_AUT_from_fit.py output/fit_results/fit_results_<TIMESTAMP>.txt

Reads which CFFs were fit from the header of the results file, then:
  - For each φ-bin in the BSA data, makes a 1×2 canvas:
      • Left: beam-spin asymmetry A_LU data, original model, and fitted model
      • Right: target-spin asymmetry A_UT predictions:
          – E-only (solid green)
          – H+Ht+E median (thick dark-green)
          – 95% CI band around H+Ht+E (filled dark-green)
Saves to:
  output/plots/BSA_AUT_bin{BIN}_{TIMESTAMP}.pdf
"""
import os
import sys
import re
import numpy as np
import matplotlib.pyplot as plt
import ROOT

# ─── Parse command-line ─────────────────────────────────────────────────────────
if len(sys.argv) != 2:
    print(__doc__)
    sys.exit(1)
fitfile = sys.argv[1]
m = re.search(r'fit_results_(\d{8}_\d{6})\.txt$', fitfile)
if not m:
    print("Error: can't extract timestamp from filename:", fitfile)
    sys.exit(1)
timestamp = m.group(1)

# ─── Parse fit results ──────────────────────────────────────────────────────────
def parse_fit_results(fname):
    with open(fname) as f:
        lines = [l.strip() for l in f if l.strip()]
    # flags line: "H 1  Ht 1  E 1  Et 0"
    flag_line = next(l for l in lines if l.startswith("H "))
    toks = flag_line.split()
    flags = {toks[i]: int(toks[i+1]) for i in range(0, len(toks), 2)}
    # parameter names
    pnames = []
    for i,l in enumerate(lines):
        if l.startswith("# parameters"):
            pnames = lines[i].split()[2:]
            break
    # values and errors
    vals = errs = None
    for i,l in enumerate(lines):
        if l.startswith("# values"):
            vals = list(map(float, lines[i+1].split()))
        if l.startswith("# errors"):
            errs = list(map(float, lines[i+1].split()))
    if vals is None or errs is None:
        raise RuntimeError("Couldn't parse fit-values/errors from file")
    return flags, pnames, np.array(vals), np.array(errs)

flags, pnames, vals, errs = parse_fit_results(fitfile)
param_map = dict(zip(pnames, vals))
error_map = dict(zip(pnames, errs))

# ─── VGG defaults for original model ────────────────────────────────────────────
orig_defaults = {
    'renormImag': 1.0, 'renormReal': 1.0,
    'r_H':0.9,   'alpha0_H':0.43, 'alpha1_H':0.85, 'n_H':1.35, 'b_H':0.4,  'Mm2_H':0.64, 'P_H':1.0,
    'r_Ht':7.0,  'alpha0_Ht':0.43,'alpha1_Ht':0.85,'n_Ht':0.6,  'b_Ht':2.0,  'Mm2_Ht':0.8,  'P_Ht':1.0,
    'r_E':0.9,   'alpha0_E':0.43, 'alpha1_E':0.85, 'n_E':1.35, 'b_E':0.4,  'Mm2_E':0.64, 'P_E':1.0,
    'r_Et':1.0,  'alpha0_Et':0.0, 'alpha1_Et':0.0, 'n_Et':0.0,  'b_Et':0.0,  'Mm2_Et':0.0,  'P_Et':0.0,
}

# ─── Load BSA φ-bins ─────────────────────────────────────────────────────────────
def load_all_bins(datafile):
    bins = []
    curr = {k: [] for k in ("phi","Q2","xB","t","Eb","A","sigA")}
    prev_phi = None
    with open(datafile) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            phi, Q2, xB, t, Eb, A, sigA = map(float, line.split())
            # detect wrap-around for φ
            if prev_phi is not None and phi < prev_phi:
                arr = {k: np.array(v) for k,v in curr.items()}
                arr["Q2m"], arr["xBm"], arr["tm"], arr["Ebm"] = (
                    arr["Q2"].mean(), arr["xB"].mean(),
                    arr["t"].mean(),  arr["Eb"].mean()
                )
                bins.append(arr)
                curr = {k: [] for k in curr}
            for k,v in zip(curr.keys(), (phi,Q2,xB,t,Eb,A,sigA)):
                curr[k].append(v)
            prev_phi = phi
    if curr["phi"]:
        arr = {k: np.array(v) for k,v in curr.items()}
        arr["Q2m"], arr["xBm"], arr["tm"], arr["Ebm"] = (
            arr["Q2"].mean(), arr["xB"].mean(),
            arr["t"].mean(),  arr["Eb"].mean()
        )
        bins.append(arr)
    return bins

bins = load_all_bins('imports/rga_prl_bsa.txt')

# ─── Compile the DVCS code into Cling ────────────────────────────────────────────
ROOT.gInterpreter.ProcessLine('#include "DVCS_xsec.C"')

# ─── Compute asymmetry helper ───────────────────────────────────────────────────
def compute_asymmetry(phi_arr, Q2_arr, xB_arr, t_arr, Eb_arr,
                      param_map, flags, kind="BSA"):
    """
    kind == "BSA" -> dvcs.BSA()
    kind == "AUT" -> dvcs.TTSAx()
    """
    # set renormalizations
    ROOT.gInterpreter.ProcessLine(f"renormImag = {param_map.get('renormImag', orig_defaults['renormImag'])};")
    ROOT.gInterpreter.ProcessLine(f"renormReal = {param_map.get('renormReal', orig_defaults['renormReal'])};")
    # switch CFFs and set parameters
    for cff in ("H","Ht","E","Et"):
        ROOT.gInterpreter.ProcessLine(f"has{cff} = {int(flags[cff])};")
        if flags[cff]:
            for k in ("r","alpha0","alpha1","n","b","Mm2","P"):
                key = f"{k}_{cff}"
                val = param_map.get(key, orig_defaults[key])
                ROOT.gInterpreter.ProcessLine(f"{key} = {val};")
    # compute
    out = []
    for phi,Q2,xB,t,Eb in zip(phi_arr, Q2_arr, xB_arr, t_arr, Eb_arr):
        if kind == "BSA":
            dvcs = ROOT.BMK_DVCS(-1, 1, 0, Eb, xB, Q2, t, phi)
            out.append(dvcs.BSA())
        else:
            # Transverse target asymmetry in x-direction
            dvcs = ROOT.BMK_DVCS(-1, 0, 1, Eb, xB, Q2, t, phi)
            out.append(dvcs.TTSAx())
    return np.array(out)

# ─── Replica parameter maps for uncertainty band ────────────────────────────────
def make_replica_maps(param_map, error_map, n_reps=200):
    reps = []
    keys = list(param_map.keys())
    for _ in range(n_reps):
        rp = {}
        for k in keys:
            mu = param_map[k]
            sigma = error_map[k]
            rp[k] = np.random.normal(mu, sigma)
        reps.append(rp)
    return reps

replica_maps = make_replica_maps(param_map, error_map, n_reps=200)

# ─── Plot loop ────────────────────────────────────────────────────────────────
os.makedirs('output/plots', exist_ok=True)
phi_grid = np.linspace(0, 360, 200)

for idx, b in enumerate(bins, start=1):
    # data & fixed grids
    phi_data, As, sigAs = b["phi"], b["A"], b["sigA"]
    Q2g = np.full_like(phi_grid, b["Q2m"])
    xBg = np.full_like(phi_grid, b["xBm"])
    tg  = np.full_like(phi_grid, b["tm"])
    Ebg = np.full_like(phi_grid, b["Ebm"])
    
    # 1×2 canvas
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14,5), sharey=True)
    
    # ── Left: BSA data + models ────────────────────────────────────────────────
    # Original model (VGG defaults)
    bsas_orig = compute_asymmetry(phi_grid, Q2g, xBg, tg, Ebg,
                                  orig_defaults, flags, kind="BSA")
    # Fitted model (all flags from fit)
    bsas_fit  = compute_asymmetry(phi_grid, Q2g, xBg, tg, Ebg,
                                  param_map, flags, kind="BSA")
    ax1.errorbar(phi_data, As, yerr=sigAs, fmt='o', ms=5,
                 color='k', label='Data')
    ax1.plot(phi_grid, bsas_orig, '-', lw=2,
             color='tab:blue', label='Original Model')
    ax1.plot(phi_grid, bsas_fit,  '--', lw=2,
             color='tab:red',  label='Fitted Model')
    ax1.set(xlim=(0,360), xticks=np.arange(0,361,60),
            ylim=(-0.6,0.6), xlabel=r'$\phi\,[°]$', ylabel=r'$A_{LU}(\phi)$')
    ax1.set_title("Beam-spin asymmetry")
    ax1.legend(loc='upper right', frameon=True)
    
    # ── Right: AUT predictions ────────────────────────────────────────────────
    # (a) E-only
    flags_E_only = flags.copy()
    flags_E_only.update({"H":0, "Ht":0, "E":1, "Et":0})
    aut_E_only = compute_asymmetry(phi_grid, Q2g, xBg, tg, Ebg,
                                   param_map, flags_E_only, kind="AUT")
    # (b) H+Ht+E central
    aut_central = compute_asymmetry(phi_grid, Q2g, xBg, tg, Ebg,
                                    param_map, flags, kind="AUT")
    # (c) uncertainty band from replicas
    all_aut = np.array([
        compute_asymmetry(phi_grid, Q2g, xBg, tg, Ebg, rp, flags, kind="AUT")
        for rp in replica_maps
    ])
    aut_low  = np.percentile(all_aut, 2.5, axis=0)
    aut_high = np.percentile(all_aut, 97.5, axis=0)
    
    ax2.plot(phi_grid, aut_E_only, '-', lw=1.5,
             color='tab:green', label='E-only')
    ax2.plot(phi_grid, aut_central, '-', lw=2.0,
             color='darkgreen', label='H+Ht+E median')
    ax2.fill_between(phi_grid, aut_low, aut_high,
                     color='darkgreen', alpha=0.3,
                     label='95% CI (H+Ht+E)')
    
    ax2.set(xlim=(0,360), xticks=np.arange(0,361,60),
            ylim=(-0.6,0.6), xlabel=r'$\phi\,[°]$')
    ax2.set_title("Target-spin asymmetry $A_{UT}$")
    ax2.legend(loc='upper right', frameon=True)
    
    plt.tight_layout()
    outname = (f'output/plots/BSA_AUT_bin{idx:02d}_{timestamp}'
               f'_Q2_{b["Q2m"]:.2f}_xB_{b["xBm"]:.3f}_t_{abs(b["tm"]):.3f}.pdf')
    fig.savefig(outname, bbox_inches='tight')
    print(f"Saved: {outname}")
    plt.close(fig)