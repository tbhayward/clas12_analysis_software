#!/usr/bin/env python3
"""
plot_BSA_AUT_from_fit.py

Usage:
    python plot_BSA_AUT_from_fit.py output/fit_results/fit_results_<TIMESTAMP>.txt

Reads the fitted CFF parameters and their uncertainties, loads all BSA data,
splits it into φ-bins, and for each bin makes a 1×2 figure:
  - Left: data + original & fitted BSA predictions (as before)
  - Right: 95% CI prediction for the target-spin asymmetry AUT vs φ

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
    # flags line: "H 1  Ht 0  E 1  Et 0"
    flag_line = next(l for l in lines if l.startswith("H "))
    toks = flag_line.split()
    flags = { toks[i]: int(toks[i+1]) for i in range(0, len(toks), 2) }
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
        raise RuntimeError("Could not parse fit-values/errors from file")
    return flags, pnames, np.array(vals), np.array(errs)

flags, pnames, vals, errs = parse_fit_results(args.fitfile)

def get_idx(name):
    return pnames.index(name) if name in pnames else None

# flatten central params & errors for replicas
central = {}
errors  = {}
# renormImag
ri_idx = get_idx("renormImag")
central["renormImag"] = vals[ri_idx]
errors["renormImag"]  = errs[ri_idx]
# shape parameters
for cff in ("H","Ht","E","Et"):
    if flags[cff]:
        for k in ("r","alpha0","alpha1","n","b","Mm2","P"):
            key = f"{k}_{cff}"
            idx = get_idx(key)
            central[key] = vals[idx]
            errors[key]  = errs[idx]

# ─── Replica generation ─────────────────────────────────────────────────────────
def generate_replicas(central_params, param_errors, n=100):
    """Return list of param_map dicts sampled from 95% CI errors."""
    reps = []
    for _ in range(n):
        pm = {}
        for k, v in central_params.items():
            sigma = param_errors[k] / 1.96
            pm[k] = np.random.normal(v, sigma)
        reps.append(pm)
    return reps

# ─── Load BSA data & split into φ-bins ────────────────────────────────────────────
def load_all_bins(fname):
    bins = []
    curr = {k: [] for k in ("phi","Q2","xB","t","Eb","A","sigA")}
    prev = None
    with open(fname) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            phi,Q2,xB,t,Eb,A,sigA = map(float, line.split())
            if prev is not None and phi < prev:
                arr = {k: np.array(v) for k,v in curr.items()}
                arr.update({
                    "Q2m": arr["Q2"].mean(),
                    "xBm": arr["xB"].mean(),
                    "tm":  arr["t"].mean(),
                    "Ebm": arr["Eb"].mean()
                })
                bins.append(arr)
                curr = {k: [] for k in curr}
            for k,v in zip(curr, (phi,Q2,xB,t,Eb,A,sigA)):
                curr[k].append(v)
            prev = phi
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
                      param_map, flags, asym="BSA"):
    """
    Set globals in ROOT, then loop over kinematics and return
    an array of asymmetries ('BSA' or 'AUT' via TTSAx).
    """
    # renormalizations
    ROOT.gInterpreter.ProcessLine(f"renormImag = {param_map.get('renormImag',1.0)};")
    # real part left at default=1.0 (not fitted for Im-only fits)
    ROOT.gInterpreter.ProcessLine(f"renormReal = 1.0;")
    # flags & parameters
    for cff in ("H","Ht","E","Et"):
        ROOT.gInterpreter.ProcessLine(f"has{cff} = {int(flags[cff])};")
        if flags[cff]:
            for k in ("r","alpha0","alpha1","n","b","Mm2","P"):
                key = f"{k}_{cff}"
                ROOT.gInterpreter.ProcessLine(f"{key} = {param_map.get(key)};")
    # compute
    vals = []
    for phi, Q2, xB, t, Eb in zip(phi_arr, Q2_arr, xB_arr, t_arr, Eb_arr):
        dvcs = ROOT.BMK_DVCS(-1, 0, 0, Eb, xB, Q2, t, phi)
        if asym=="BSA":
            vals.append(dvcs.BSA())
        elif asym=="AUT":
            vals.append(dvcs.TTSAx())
        else:
            raise ValueError("Unknown asymmetry type")
    return np.array(vals)

# ─── Main plotting loop ─────────────────────────────────────────────────────────
os.makedirs("output/plots", exist_ok=True)

# prepare replicas for AUT
N_REP = 100
replica_params = generate_replicas(central, errors, n=N_REP)

phi_grid = np.linspace(0, 360, 200)

for ibin, b in enumerate(bins, start=1):
    # data
    phi_data, As, sigAs = b["phi"], b["A"], b["sigA"]
    # kinematic grids
    Q2g = np.full_like(phi_grid, b["Q2m"])
    xBg = np.full_like(phi_grid, b["xBm"])
    tg  = np.full_like(phi_grid, b["tm"])
    Ebg = np.full_like(phi_grid, b["Ebm"])

    # left panel: BSA
    # original = defaults, flags
    orig_defaults = {**{k: central[k] for k in central}, **{k: 1.0 for k in ("renormImag",)}}
    bsas_orig = compute_asymmetry(phi_grid, Q2g, xBg, tg, Ebg,
                                  orig_defaults, flags, asym="BSA")
    # fitted central
    fitted = dict(zip(pnames, vals))
    bsas_fit  = compute_asymmetry(phi_grid, Q2g, xBg, tg, Ebg,
                                  fitted, flags, asym="BSA")

    # right panel: AUT prediction
    # central
    aut_central = compute_asymmetry(phi_grid, Q2g, xBg, tg, Ebg,
                                    fitted, flags, asym="AUT")
    # replicas
    all_aut = np.array([
        compute_asymmetry(phi_grid, Q2g, xBg, tg, Ebg, rp, flags, asym="AUT")
        for rp in replica_params
    ])
    aut_med   = np.median(all_aut, axis=0)
    aut_low   = np.percentile(all_aut, 2.5, axis=0)
    aut_high  = np.percentile(all_aut, 97.5, axis=0)

    # make figure
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharex=True)
    # left: BSA
    ax = axes[0]
    ax.errorbar(phi_data, As, yerr=sigAs, fmt='o', ms=5, color='k', label='Data')
    ax.plot(phi_grid, bsas_orig, '-',  lw=2, color='tab:blue', label='Original Model')
    ax.plot(phi_grid, bsas_fit,  '--', lw=2, color='tab:red',  label='Fitted Model')
    ax.set_xlim(0,360)
    ax.set_xticks([0,60,120,180,240,300,360])
    ax.set_ylim(-0.6, 0.6)
    ax.set_xlabel(r'$\phi\;[\mathrm{deg}]$')
    ax.set_ylabel(r'$A_{LU}(\phi)$')
    ax.set_title("Beam‐spin asymmetry")
    ax.legend(loc='upper right', frameon=True)

    # right: AUT
    ax = axes[1]
    ax.fill_between(phi_grid, aut_low, aut_high,
                    color='tab:green', alpha=0.3, label='95% CI')
    ax.plot(phi_grid, aut_med, '-', lw=2, color='tab:green', label='Median Prediction')
    ax.set_xlim(0,360)
    ax.set_xticks([0,60,120,180,240,300,360])
    ax.set_ylim(-0.6, 0.6)
    ax.set_xlabel(r'$\phi\;[\mathrm{deg}]$')
    ax.set_ylabel(r'$A_{UT}(\phi)$')
    ax.set_title("Target‐spin asymmetry (prediction)")
    ax.legend(loc='upper right', frameon=True)

    plt.tight_layout()
    outname = f"output/plots/BSA_AUT_bin{ibin:02d}_{timestamp}.pdf"
    fig.savefig(outname, bbox_inches='tight')
    print(f"Saved: {outname}")
    plt.close(fig)