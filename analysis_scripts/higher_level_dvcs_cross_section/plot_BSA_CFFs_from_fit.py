#!/usr/bin/env python3
"""
plot_BSA_AUT_from_fit.py

Usage:
    python plot_BSA_AUT_from_fit.py output/fit_results/fit_results_<TIMESTAMP>.txt

Reads the fitted CFF parameters and their uncertainties, loads all BSA data,
splits it into φ-bins, and for each bin makes a 1×2 figure:
  - Left: data + original & fitted BSA predictions
  - Right: two solid‐green lines:
      • E-only AUT prediction 
      • H+E AUT prediction (median) with a 95% CI band

Saves to:
  output/plots/BSA_AUT_bin{BIN:02d}_{TIMESTAMP}.pdf
"""
import os
import sys
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt
import ROOT

# ─── Parse command-line ─────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('fitfile', help='Fit results file')
args = parser.parse_args()

m = re.search(r'fit_results_(\d{8}_\d{6})\.txt$', args.fitfile)
if not m:
    print("ERROR: can't extract timestamp from filename")
    sys.exit(1)
timestamp = m.group(1)

# ─── Load fit results & flags ───────────────────────────────────────────────────
def parse_fit_results(fname):
    with open(fname) as f:
        lines = [l.strip() for l in f if l.strip()]
    flag_line = next(l for l in lines if l.startswith("H "))
    toks = flag_line.split()
    flags = { toks[i]: int(toks[i+1]) for i in range(0, len(toks), 2) }
    pnames = []
    for i,l in enumerate(lines):
        if l.startswith("# parameters"):
            pnames = lines[i].split()[2:]
            break
    vals = errs = None
    for i,l in enumerate(lines):
        if l.startswith("# values"):
            vals = list(map(float, lines[i+1].split()))
        if l.startswith("# errors"):
            errs = list(map(float, lines[i+1].split()))
    if vals is None or errs is None:
        raise RuntimeError("Could not parse fit-values/errors from file")
    return flags, pnames, np.array(vals), np.array(errs)

flags, pnames, vals, errs = parse_fit_results(args.fitfile)

def get_idx(name):
    return pnames.index(name) if name in pnames else None

# build central & error dicts
central = {}
errors  = {}
# renormImag
central["renormImag"] = vals[get_idx("renormImag")]
errors ["renormImag"] = errs[get_idx("renormImag")]
# shape parameters
for cff in ("H","Ht","E","Et"):
    if flags[cff]:
        for k in ("r","alpha0","alpha1","n","b","Mm2","P"):
            key = f"{k}_{cff}"
            idx = get_idx(key)
            central[key] = vals[idx]
            errors [key] = errs[idx]

# ─── Replica generation ─────────────────────────────────────────────────────────
def generate_replicas(central_params, param_errors, n=100):
    reps = []
    for _ in range(n):
        pm = {}
        for k,v in central_params.items():
            sigma = param_errors[k]/1.96
            pm[k] = np.random.normal(v, sigma)
        reps.append(pm)
    return reps

N_REP = 100
replica_params = generate_replicas(central, errors, n=N_REP)

# ─── Load & bin BSA data ─────────────────────────────────────────────────────────
def load_all_bins(fname):
    bins = []
    curr = {k: [] for k in ("phi","Q2","xB","t","Eb","A","sigA")}
    prev = None
    with open(fname) as f:
        for line in f:
            if not line.strip() or line.startswith("#"): continue
            φ,Q2,xB,t,Eb,A,σA = map(float, line.split())
            if prev is not None and φ < prev:
                arr = {k: np.array(v) for k,v in curr.items()}
                arr.update({
                    "Q2m": arr["Q2"].mean(),
                    "xBm": arr["xB"].mean(),
                    "tm":  arr["t"].mean(),
                    "Ebm": arr["Eb"].mean()
                })
                bins.append(arr)
                curr = {k: [] for k in curr}
            for k,v in zip(curr, (φ,Q2,xB,t,Eb,A,σA)):
                curr[k].append(v)
            prev = φ
    if curr["phi"]:
        arr = {k: np.array(v) for k,v in curr.items()}
        arr.update({
            "Q2m": arr["Q2"].mean(),
            "xBm": arr["xB"].mean(),
            "tm":  arr["t"].mean(),
            "Ebm": arr["Eb"].mean()
        })
        bins.append(arr)
    return bins

bins = load_all_bins("imports/rga_prl_bsa.txt")

# ─── Prepare ROOT DVCS code ──────────────────────────────────────────────────────
ROOT.gInterpreter.ProcessLine('#include "DVCS_xsec.C"')

def compute_asymmetry(phi_arr, Q2_arr, xB_arr, t_arr, Eb_arr,
                      param_map, flags_map, asym="BSA"):
    # renormalizations
    ROOT.gInterpreter.ProcessLine(f"renormImag = {param_map.get('renormImag',1.0)};")
    ROOT.gInterpreter.ProcessLine("renormReal = 1.0;")
    # flags & params
    for cff in ("H","Ht","E","Et"):
        ROOT.gInterpreter.ProcessLine(f"has{cff} = {int(flags_map[cff])};")
        if flags_map[cff]:
            for k in ("r","alpha0","alpha1","n","b","Mm2","P"):
                ROOT.gInterpreter.ProcessLine(
                    f"{k}_{cff} = {param_map[k+'_'+cff]};"
                )
    out = []
    for φ,Q2,xB,t,Eb in zip(phi_arr, Q2_arr, xB_arr, t_arr, Eb_arr):
        dvcs = ROOT.BMK_DVCS(-1, 0, 0, Eb, xB, Q2, t, φ)
        out.append(dvcs.BSA() if asym=="BSA" else dvcs.TTSAx())
    return np.array(out)

phi_grid = np.linspace(0,360,200)
os.makedirs("output/plots", exist_ok=True)

for ibin, b in enumerate(bins, start=1):
    # data
    φ_dat, A_dat, σA = b["phi"], b["A"], b["sigA"]
    # kinematic grids
    Q2g = np.full_like(phi_grid, b["Q2m"])
    xBg = np.full_like(phi_grid, b["xBm"])
    tg  = np.full_like(phi_grid, b["tm"])
    Ebg = np.full_like(phi_grid, b["Ebm"])

    # 1) BSA panel
    # original
    orig = central.copy()
    orig["renormImag"] = 1.0
    bsas_orig = compute_asymmetry(phi_grid, Q2g, xBg, tg, Ebg,
                                  orig, flags, asym="BSA")
    # fit central
    fit_map   = dict(zip(pnames, vals))
    bsas_fit  = compute_asymmetry(phi_grid, Q2g, xBg, tg, Ebg,
                                  fit_map, flags, asym="BSA")

    # 2) AUT panel
    # E-only
    flags_E = {"H":0,"Ht":0,"E":1,"Et":0}
    aut_E = compute_asymmetry(phi_grid, Q2g, xBg, tg, Ebg,
                              fit_map, flags_E, asym="AUT")
    # H+E
    flags_HE = {"H":1,"Ht":0,"E":1,"Et":0}
    aut_med = compute_asymmetry(phi_grid, Q2g, xBg, tg, Ebg,
                                fit_map, flags_HE, asym="AUT")
    all_aut = np.array([
        compute_asymmetry(phi_grid, Q2g, xBg, tg, Ebg, rp, flags_HE, asym="AUT")
        for rp in replica_params
    ])
    aut_lo = np.percentile(all_aut, 2.5,  axis=0)
    aut_hi = np.percentile(all_aut, 97.5, axis=0)

    # draw
    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(12,5), sharex=True)
    # suptitle with kinematics
    fig.suptitle(
        (r'$\langle Q^2\rangle={:.2f}\,\mathrm{{GeV}}^2,\;\langle x_B\rangle={:.3f},\;\langle -t\rangle={:.3f}\,\mathrm{{GeV}}^2$'
         ).format(b["Q2m"], b["xBm"], -b["tm"]),
        fontsize=14, y=1.02
    )

    # left: BSA
    ax1.errorbar(φ_dat, A_dat, yerr=σA, fmt='o', ms=5, color='k',
                 label='Data')
    ax1.plot(phi_grid, bsas_orig, '-',  lw=2, color='tab:blue',
             label='Original')
    ax1.plot(phi_grid, bsas_fit,  '--', lw=2, color='tab:red',
             label='Fit')
    ax1.set(xlim=(0,360), xticks=np.arange(0,361,60),
            ylim=(-0.6,0.6), xlabel=r'$\phi\,[°]$', ylabel=r'$A_{LU}$')
    ax1.set_title("Beam-spin asymmetry")
    ax1.legend(loc='upper right', frameon=True)

    # right: AUT
    ax2.plot(phi_grid, aut_E,   '-', lw=1.5,
             color='tab:green', label='AUT, E-only')
    ax2.plot(phi_grid, aut_med, '-', lw=2.0,
             color='darkgreen', label='AUT, H+E median')
    ax2.fill_between(phi_grid, aut_lo, aut_hi,
                     color='darkgreen', alpha=0.3,
                     label='95% CI (H+E)')
    ax2.set(xlim=(0,360), xticks=np.arange(0,361,60),
            ylim=(-0.6,0.6), xlabel=r'$\phi\,[°]$', ylabel=r'$A_{UT}$')
    ax2.set_title("Target-spin asymmetry prediction")
    ax2.legend(loc='upper right', frameon=True)

    plt.tight_layout()
    outfn = f"output/plots/BSA_AUT_bin{ibin:02d}_{timestamp}.pdf"
    fig.savefig(outfn, bbox_inches='tight')
    print("Saved:", outfn)
    plt.close(fig)