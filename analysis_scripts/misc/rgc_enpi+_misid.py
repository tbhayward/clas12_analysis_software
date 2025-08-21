#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Misidentification vs -t for four xB slices (MC).

For each xB slice and -t bin:
  fraction_pi_to_e = N(matching_e_pid == -211) / N(total)
  fraction_pi_to_k = N(matching_p1_pid == -321) / N(total)

Cuts:
  0.81 < Mx2 < 1.00
  x in slice range
  -t in one of the six bins defined below
  (No fiducial-status cuts.)

Y-axis is a FRACTION with limits [0, 0.01].

Usage:
  python misid_vs_t.py <input.root>

Output:
  output/enpi+/misidentification.pdf
  (Also prints exact values per bin to terminal.)
"""

import sys
import os
import numpy as np

try:
    import uproot
except Exception:
    print("[ERROR] Please install uproot (e.g. pip install uproot)")
    sys.exit(1)

import matplotlib.pyplot as plt

# ─────────────────────────────────────────────────────────────────────
# Config
# ─────────────────────────────────────────────────────────────────────
TREE_NAME = "PhysicsEvents"

# xB slices
X_SLICES = {
    "Low":     (0.10, 0.25),
    "MidLow":  (0.25, 0.35),
    "MidHigh": (0.35, 0.45),
    "High":    (0.45, 0.60),
}
SLICE_ORDER = ["Low", "MidLow", "MidHigh", "High"]

# -t bin edges in GeV^2 (ascending) → bins [edge[i], edge[i+1])
T_EDGES_POS = np.array([0.05, 0.25, 0.45, 0.65, 0.85, 1.05, 1.25], dtype=float)
T_CENTERS   = 0.5 * (T_EDGES_POS[:-1] + T_EDGES_POS[1:])
NTBINS      = len(T_EDGES_POS) - 1

OUT_DIR  = os.path.join("output", "enpi+")
OUT_PATH = os.path.join(OUT_DIR, "misidentification.pdf")

# Plot styling (two series)
STYLE_E = dict(marker="o",  color="tab:blue",  label=r"$\pi^{-}\!\to e^{-}$")
STYLE_K = dict(marker="^",  color="tab:orange",label=r"$\pi^{-}\!\to K^{-}$")
MSIZE   = 36

# ─────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────
def slice_index_from_x(x):
    """Return array of indices in {0..3} for xB slice, else -1."""
    idx = np.full(x.shape, -1, dtype=np.int16)
    for k, name in enumerate(SLICE_ORDER):
        xa, xb = X_SLICES[name]
        m = (x > xa) & (x < xb)
        idx[m] = k
    return idx

def tbin_index_from_tpos(tpos):
    """Return -t bin index in {0..NTBINS-1}, else -1 if outside range."""
    bix = np.digitize(tpos, T_EDGES_POS, right=False) - 1
    bad = (bix < 0) | (bix >= NTBINS)
    bix[bad] = -1
    return bix.astype(np.int16)

# ─────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────
def main():
    if len(sys.argv) != 2:
        print("Usage: python misid_vs_t.py <input.root>")
        sys.exit(1)

    infile = sys.argv[1]
    if not os.path.isfile(infile):
        print(f"[ERROR] File not found: {infile}")
        sys.exit(1)

    branches = [
        "x", "t", "Mx2",
        "matching_e_pid", "matching_p1_pid",
    ]

    totals = np.zeros((4, NTBINS), dtype=np.int64)
    bad_e  = np.zeros((4, NTBINS), dtype=np.int64)  # pi- -> e-
    bad_k  = np.zeros((4, NTBINS), dtype=np.int64)  # pi- -> K-

    step = "200 MB"
    try:
        itr = uproot.iterate({infile: TREE_NAME}, branches, step_size=step, library="np")
    except Exception as e:
        print(f"[ERROR] Failed to open '{infile}' or tree '{TREE_NAME}': {e}")
        sys.exit(1)

    for arrays in itr:
        x   = arrays["x"]
        t   = arrays["t"]
        mx2 = arrays["Mx2"]
        ep  = arrays["matching_e_pid"]
        p1  = arrays["matching_p1_pid"]

        # Exclusivity window only
        base = (mx2 > 0.81) & (mx2 < 1.00)
        if not np.any(base):
            continue

        tpos = -t[base]
        xs   = x[base]
        e_is_pi_minus = (ep[base] != 11)   # π− → e− (electron candidate is actually π−)
        p1_is_k_minus = (p1[base] != 211)   # π− → K− (hadron candidate flagged as K−)

        sidx = slice_index_from_x(xs)           # 0..3 or -1
        tbix = tbin_index_from_tpos(tpos)       # 0..5 or -1
        valid = (sidx >= 0) & (tbix >= 0)
        if not np.any(valid):
            continue

        sidx = sidx[valid]
        tbix = tbix[valid]

        np.add.at(totals, (sidx, tbix), 1)

        ee = e_is_pi_minus[valid]
        kk = p1_is_k_minus[valid]
        if np.any(ee):
            np.add.at(bad_e, (sidx[ee], tbix[ee]), 1)
        if np.any(kk):
            np.add.at(bad_k, (sidx[kk], tbix[kk]), 1)

    # Fractions (not percent)
    with np.errstate(divide="ignore", invalid="ignore"):
        frac_e = np.where(totals > 0, bad_e / totals, np.nan)
        frac_k = np.where(totals > 0, bad_k / totals, np.nan)

    # ── Print exact values to terminal ───────────────────────────────
    print("\n=== Misidentification fractions by xB slice and -t bin ===")
    for si, sname in enumerate(SLICE_ORDER):
        xa, xb = X_SLICES[sname]
        print(f"\nSlice {sname}  ( {xa:.2f} < x_B < {xb:.2f} )")
        print("bin  [-t_min, -t_max)   N_total   N(pi->e)  frac(pi->e)   N(pi->K)  frac(pi->K)")
        for bi in range(NTBINS):
            tmin, tmax = T_EDGES_POS[bi], T_EDGES_POS[bi+1]
            nt  = int(totals[si, bi])
            ne  = int(bad_e[si, bi])
            nk  = int(bad_k[si, bi])
            fe  = float(frac_e[si, bi]) if nt > 0 else float("nan")
            fk  = float(frac_k[si, bi]) if nt > 0 else float("nan")
            print(f"{bi:>2d}   [{tmin:5.2f}, {tmax:5.2f})   {nt:7d}   {ne:8d}   {fe:11.6f}   {nk:8d}   {fk:11.6f}")

    # ── Plotting ─────────────────────────────────────────────────────
    os.makedirs(OUT_DIR, exist_ok=True)
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.6), sharey=True)
    fig.suptitle(
        r"Misidentification vs $-t$"
        "\n" r"$0.81<M_{x}^{2}<1.00$  (no fiducial cuts)",
        fontsize=13, y=0.98
    )

    for k, name in enumerate(SLICE_ORDER):
        ax = axes[k]
        y_e = frac_e[k, :]
        y_k = frac_k[k, :]
        m_e = np.isfinite(y_e)
        m_k = np.isfinite(y_k)

        if np.any(m_e):
            ax.scatter(T_CENTERS[m_e], y_e[m_e], s=MSIZE, **STYLE_E)
        if np.any(m_k):
            ax.scatter(T_CENTERS[m_k], y_k[m_k], s=MSIZE, **STYLE_K)

        ax.set_xlim(T_EDGES_POS[0], T_EDGES_POS[-1])
        ax.set_ylim(0.0, 0.01)  # fraction scale 0..1%
        ax.set_xlabel(r"$-t\ (\mathrm{GeV}^{2})$")
        if k == 0:
            ax.set_ylabel("misID fraction")
        ax.grid(True, linestyle="--", alpha=0.35)
        xa, xb = X_SLICES[name]
        ax.set_title(rf"${xa:.2f} < x_{{B}} < {xb:.2f}$", fontsize=11)
        ax.legend(loc="upper left", frameon=True, edgecolor="black", fontsize=10)

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(OUT_PATH, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved: {OUT_PATH}")

if __name__ == "__main__":
    main()