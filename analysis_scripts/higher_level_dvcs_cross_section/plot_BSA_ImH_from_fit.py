#!/usr/bin/env python3
"""
compare_bsa_bin1.py

Compare the first φ‐scan bin of BSA data to the original and fitted BMK_DVCS predictions.

Usage:
    python compare_bsa_bin1.py <bsa_data.txt> <fit_results.txt>
"""

import os
import sys
import re
import numpy as np
import matplotlib.pyplot as plt

# ─── PyROOT setup ─────────────────────────────────────────────────────────────
import ROOT
ROOT.gROOT.SetBatch(True)
# Load and compile your DVCS_xsec.C (must be in cwd or adjust path)
if not ROOT.gSystem.Load("DVCS_xsec_C") >= 0:
    ROOT.gROOT.ProcessLine('.L DVCS_xsec.C+')

# ─── parse command-line ────────────────────────────────────────────────────────
if len(sys.argv) != 3:
    print(__doc__)
    sys.exit(1)

bsa_file, fit_file = sys.argv[1], sys.argv[2]

m = re.search(r'fit_results_(\d{8}_\d{6})\.txt$', fit_file)
if not m:
    print("ERROR: can't extract timestamp from fit file name")
    sys.exit(1)
timestamp = m.group(1)

# ─── read first bin of BSA data ────────────────────────────────────────────────
phis = []
As   = []
errs = []
with open(bsa_file) as f:
    prev_phi = None
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'): continue
        phi, Q2, xB, t, Eb, A, sigA = map(float, line.split())
        if prev_phi is None:
            prev_phi = phi
        elif phi < prev_phi:
            # φ wrapped around → new bin starts
            break
        phis.append(phi)
        As.append(A)
        errs.append(sigA)
        prev_phi = phi

phis = np.array(phis)
As   = np.array(As)
errs = np.array(errs)

# ─── parse your fit results ───────────────────────────────────────────────────
vals = None
with open(fit_file) as f:
    lines = [L.strip() for L in f if L.strip()]
for i,L in enumerate(lines):
    if L.startswith("# values"):
        vals = list(map(float, lines[i+1].split()))
if vals is None:
    raise RuntimeError("Couldn't find '# values' in fit results")
renormI, alpha0, alpha1, n_val, b_val, Mm2_val, P_val, renormR = vals

# ─── set global model parameters ──────────────────────────────────────────────
ROOT.renormImag = renormI
ROOT.alpha0     = alpha0
ROOT.alpha1     = alpha1
ROOT.n_val      = n_val
ROOT.b_val      = b_val
ROOT.Mm2_val    = Mm2_val
ROOT.P_val      = P_val
ROOT.renormReal = renormR

# also store the originals for convenience
orig = {
    'renormImag': 1.0,
    'alpha0':     0.43,
    'alpha1':     0.85,
    'n_val':      1.35,
    'b_val':      0.4,
    'Mm2_val':    0.64,
    'P_val':      1.0,
}

# ─── helper to compute BSA(φ) via PyROOT ────────────────────────────────────────
def compute_bsa(phi_arr, params):
    # set the global params in ROOT
    for name, val in params.items():
        setattr(ROOT, name, val)
    out = []
    for φ in phi_arr:
        dvcs = ROOT.BMK_DVCS(-1, 1, 0,
                             Eb, xB, Q2, t, φ)
        out.append(dvcs.BSA())
    return np.array(out)

# But we need the kinematics (Q2, xB, t, Eb) for this first bin → just take from first line
# Re-read that one:
with open(bsa_file) as f:
    for line in f:
        if line.strip().startswith('#'): continue
        phi0, Q2, xB, t, Eb, *_ = map(float, line.split())
        break

# ─── compute model curves ─────────────────────────────────────────────────────
bsa_orig = compute_bsa(phis, orig)
bsa_fit  = compute_bsa(phis, {
    'renormImag': renormI,
    'alpha0':     alpha0,
    'alpha1':     alpha1,
    'n_val':      n_val,
    'b_val':      b_val,
    'Mm2_val':    Mm2_val,
    'P_val':      P_val,
})

# ─── plot ─────────────────────────────────────────────────────────────────────
plt.style.use('classic')
plt.rcParams.update({
    'font.size': 14,
    'font.family': 'serif',
})

fig, ax = plt.subplots(figsize=(8,5))

ax.errorbar(phis, As, yerr=errs,
            fmt='o', color='k', label='Data')

ax.plot(phis, bsa_orig, '-', color='tab:blue', lw=2,
        label='Original parameters')
ax.plot(phis, bsa_fit,  '--', color='tab:red',  lw=2,
        label='RGA pass-1 fit')

ax.set_xlim(0, 360)
ax.set_xlabel(r'$\phi\;[\mathrm{deg}]$')
ax.set_ylabel(r'$A_{LU}(\phi)$')
ax.legend(loc='upper right', frameon=False)

outdir = 'output/plots'
os.makedirs(outdir, exist_ok=True)
outname = f"{outdir}/BSA_bin1_{timestamp}.pdf"
fig.tight_layout()
fig.savefig(outname)
print("Saved plot to", outname)