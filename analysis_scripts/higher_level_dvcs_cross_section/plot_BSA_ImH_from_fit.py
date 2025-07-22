#!/usr/bin/env python3
"""
compare_bsa_bin1.py

Compare the first φ‐scan bin of BSA data (from imports/rga_prl_bsa.txt)
to the original and fitted BMK_DVCS predictions.

Usage:
    python compare_bsa_bin1.py <fit_results.txt>
"""
import os
import sys
import re
import numpy as np
import matplotlib.pyplot as plt

# ─── PyROOT setup ─────────────────────────────────────────────────────────────
import ROOT
ROOT.gROOT.SetBatch(True)
# compile+load your DVCS_xsec.C once and for all
ROOT.gROOT.ProcessLine('.L DVCS_xsec.C+')

# ─── command‐line ─────────────────────────────────────────────────────────────
if len(sys.argv) != 2:
    print(__doc__)
    sys.exit(1)

fit_file = sys.argv[1]
# hard-coded BSA data path:
bsa_file = 'imports/rga_prl_bsa.txt'

m = re.search(r'fit_results_(\d{8}_\d{6})\.txt$', fit_file)
if not m:
    print("ERROR: can't extract timestamp from fit file name")
    sys.exit(1)
timestamp = m.group(1)

# ─── read first φ‐scan bin from BSA data ───────────────────────────────────────
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
            # φ wrapped → next bin
            break
        phis.append(phi)
        As.append(A)
        errs.append(sigA)
        prev_phi = phi

phis = np.array(phis); As = np.array(As); errs = np.array(errs)

# ─── grab the kinematics from the first line ─────────────────────────────────
with open(bsa_file) as f:
    for line in f:
        if line.strip().startswith('#'): continue
        phi0, Q2, xB, t, Eb, *_ = map(float, line.split())
        break

# ─── parse the fit results ───────────────────────────────────────────────────
vals = None
with open(fit_file) as f:
    lines = [L.strip() for L in f if L.strip()]
for i,L in enumerate(lines):
    if L.startswith("# values"):
        vals = list(map(float, lines[i+1].split()))
if vals is None:
    raise RuntimeError("Couldn't find '# values' in fit results")
renormI, alpha0, alpha1, n_val, b_val, Mm2_val, P_val, renormR = vals

# ─── set globals in ROOT for your fit and for the original defaults ──────────
# fitted
ROOT.renormImag = renormI
ROOT.alpha0     = alpha0
ROOT.alpha1     = alpha1
ROOT.n_val      = n_val
ROOT.b_val      = b_val
ROOT.Mm2_val    = Mm2_val
ROOT.P_val      = P_val
ROOT.renormReal = renormR

# original VGG defaults
orig = {
    'renormImag': 1.0,
    'alpha0':     0.43,
    'alpha1':     0.85,
    'n_val':      1.35,
    'b_val':      0.4,
    'Mm2_val':    0.64,
    'P_val':      1.0,
}

# ─── helper: compute BSA(φ) calling the C++ class via PyROOT ─────────────────
def compute_bsa(phi_arr, params):
    for name,val in params.items():
        setattr(ROOT, name, val)
    out = []
    for φ in phi_arr:
        dvcs = ROOT.BMK_DVCS(-1, 1, 0,
                             Eb, xB, Q2, t, φ)
        out.append(dvcs.BSA())
    return np.array(out)

# ─── build the three curves ───────────────────────────────────────────────────
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