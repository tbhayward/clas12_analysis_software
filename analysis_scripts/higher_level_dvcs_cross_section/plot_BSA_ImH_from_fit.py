#!/usr/bin/env python3
"""
compare_bsa_bin1.py

Compare the first φ‐scan bin of BSA data (from imports/rga_prl_bsa.txt)
to the original and fitted BMK_DVCS predictions.

Usage:
    python compare_bsa_bin1.py output/fit_results/fit_results_<TIMESTAMP>.txt
"""
import os, sys, re
import numpy as np
import matplotlib.pyplot as plt

# ─── PyROOT setup ─────────────────────────────────────────────────────────────
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.ProcessLine('.L DVCS_xsec.C+')

# ─── parse args ───────────────────────────────────────────────────────────────
if len(sys.argv) != 2:
    print(__doc__); sys.exit(1)
fit_file = sys.argv[1]
m = re.search(r'fit_results_(\d{8}_\d{6})\.txt$', fit_file)
if not m:
    print("ERROR: can't extract timestamp from fit filename"); sys.exit(1)
timestamp = m.group(1)

# ─── load the first φ‐scan bin of BSA data ─────────────────────────────────────
bsa_file = 'imports/rga_prl_bsa.txt'
phis = []; As = []; errs = []
Q2s = []; xBs = []; ts = []; Ebs = []
with open(bsa_file) as f:
    prev_phi = None
    for line in f:
        if not line.strip() or line.startswith('#'): continue
        phi, Q2, xB, t, Eb, A, sigA = map(float, line.split())
        if prev_phi is None:
            prev_phi = phi
        elif phi < prev_phi:
            # wrapped → new bin
            break
        phis.append(phi); As.append(A); errs.append(sigA)
        Q2s.append(Q2); xBs.append(xB); ts.append(t); Ebs.append(Eb)
        prev_phi = phi

phis = np.array(phis); As = np.array(As); errs = np.array(errs)
Q2m  = np.mean(Q2s)
xBm  = np.mean(xBs)
tm   = np.mean(ts)      # note ts are negative; we'll plot −t

# ─── parse the fit results ───────────────────────────────────────────────────
vals = None
with open(fit_file) as f:
    lines = [L.strip() for L in f if L.strip()]
for i,L in enumerate(lines):
    if L.startswith("# values"):
        vals = list(map(float, lines[i+1].split()))
if vals is None:
    raise RuntimeError("Couldn't find '# values' in fit file")
renormI, alpha0, alpha1, n_val, b_val, Mm2_val, P_val, renormR = vals

# ─── push globals into ROOT for both “orig” and “fit” cases ────────────────
# enable only the H‐term
ROOT.hasH  = True
ROOT.hasHt = False
ROOT.hasE  = False
ROOT.hasEt = False

# original VGG defaults
orig = dict(
    renormImag=1.0, alpha0=0.43, alpha1=0.85,
    n_val=1.35,    b_val=0.4,
    Mm2_val=0.64,  P_val=1.0,
)

# fitted
fit  = dict(
    renormImag=renormI, alpha0=alpha0, alpha1=alpha1,
    n_val=n_val,        b_val=b_val,
    Mm2_val=Mm2_val,    P_val=P_val,
)

# ─── helper to compute BSA(φ) via PyROOT BMK_DVCS ──────────────────────────
def compute_bsa(phi_arr, params):
    for name,val in params.items():
        setattr(ROOT, name, val)
    out = []
    for φ in phi_arr:
        dvcs = ROOT.BMK_DVCS(-1, 1, 0,
                             np.mean(Ebs),  # beam energy
                             xBm, Q2m, tm, φ)
        out.append(dvcs.BSA())
    return np.array(out)

bsa_orig = compute_bsa(phis, orig)
bsa_fit  = compute_bsa(phis, fit)

# ─── make the plot ───────────────────────────────────────────────────────────
plt.style.use('classic')
plt.rcParams.update({'font.size':14, 'font.family':'serif'})

fig, ax = plt.subplots(figsize=(8,5))
ax.errorbar(phis, As, yerr=errs, fmt='o', color='k', label='Data')
ax.plot(phis, bsa_orig, '-', color='tab:blue', lw=2, label='Original parameters')
ax.plot(phis, bsa_fit,  '--', color='tab:red',  lw=2, label='RGA pass-1 fit')

# axis ranges & ticks
ax.set_xlim(0, 360)
ax.set_xticks([0,60,120,180,240,300,360])
ax.set_ylim(-0.6, 0.6)
ax.set_xlabel(r'$\phi\;[\mathrm{deg}]$')
ax.set_ylabel(r'$A_{LU}(\phi)$')

# title with averaged kinematics
ax.set_title(
    (r'$\langle Q^2\rangle={:.2f}\,\mathrm{{GeV}^2},\;\langle x_B\rangle={:.3f},\;\langle -t\rangle={:.3f}\,\mathrm{{GeV}^2}$'
    ).format(Q2m, xBm, -tm),
    pad=12
)

ax.legend(loc='upper right', frameon=False)
fig.tight_layout()

outdir = 'output/plots'
os.makedirs(outdir, exist_ok=True)
outname = f"{outdir}/BSA_bin1_{timestamp}.pdf"
fig.savefig(outname)
print("Saved plot to", outname)