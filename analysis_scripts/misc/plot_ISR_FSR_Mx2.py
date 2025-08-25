#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare baseline vs simulated ISR+FSR for three observables on a 1x3 canvas:
  [Left]   xB (branch: 'x'), show migration INTO [0.10, 0.25]
  [Middle] -t' (from branch 'tprime'; plot -t'), migration INTO [0.05, 0.25]
  [Right]  Mx2 (branch: 'Mx2'), migration INTO (0.81, 1.00)

"Migration" means: events that were OUTSIDE the target range in the baseline file,
but INSIDE the target range in the ISR+FSR file (matched by 'evnnum').
Those baseline values are drawn as a thin dashed gray histogram.

Usage:
  python isr_fsr_compare_1x3.py <baseline.root> <isr_fsr.root>

Output:
  output/enpi+/ISR_FSR_Mx2.pdf
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

try:
    import uproot
except Exception:
    print("[ERROR] Please install uproot (e.g. pip install uproot)")
    sys.exit(1)

TREE = "PhysicsEvents"
OUT_DIR = os.path.join("output", "enpi+")
OUT_PATH = os.path.join(OUT_DIR, "ISR_FSR_Mx2.pdf")

# Binning / ranges
XB_RANGE     = (0.00, 1.00)
XB_BINS      = 50
XB_TARGET    = (0.10, 0.25)

TP_RANGE     = (0.00, 1.30)   # for -t'
TP_BINS      = 50
TP_TARGET    = (0.05, 0.25)   # in -t'

MX2_RANGE    = (0.00, 1.30)
MX2_BINS     = 50
MX2_TARGET   = (0.81, 1.00)

def load_branches(root_path, branches):
    """Load selected branches from TREE as numpy arrays."""
    try:
        with uproot.open(root_path) as f:
            t = f[TREE]
            arrs = t.arrays(branches, library="np")
        return {k: np.asarray(arrs[k]) for k in branches}
    except Exception as e:
        print(f"[ERROR] Failed to read '{root_path}' ({TREE}): {e}")
        sys.exit(2)

def make_event_map(evnums, xb, tprime, mx2):
    """
    Build a map: evnnum -> (xb, minus_tprime, mx2).
    If duplicate evnnums exist, the last one wins.
    """
    ev = np.asarray(evnums)
    # Ensure 1D and integer-like keys
    try:
        keys = ev.astype(np.int64, copy=False)
    except Exception:
        keys = ev  # fallback (should still be hashable)
    minus_tp = -np.asarray(tprime, dtype=float)
    xb = np.asarray(xb, dtype=float)
    mx2 = np.asarray(mx2, dtype=float)

    d = {}
    for k, xv, mtp, m2 in zip(keys, xb, minus_tp, mx2):
        d[k] = (xv, mtp, m2)
    return d

def migrated_baseline_values(basemap, evnums_rad, xb_rad, tprime_rad, mx2_rad,
                             which="xb", target=(0,1)):
    """
    Identify baseline values for events that were OUTSIDE the target range in baseline
    but INSIDE in the ISR+FSR sample, matched by evnnum.

    which: "xb", "tprime" (interpreted as -t'), or "mx2"
    Returns: np.array of the *baseline* observable values for migrated events.
    """
    keys_rad = np.asarray(evnums_rad).astype(np.int64, copy=False)
    minus_tp_rad = -np.asarray(tprime_rad, dtype=float)
    xb_rad = np.asarray(xb_rad, dtype=float)
    mx2_rad = np.asarray(mx2_rad, dtype=float)

    lo, hi = target
    acc = []

    for k, xr, mtr, m2r in zip(keys_rad, xb_rad, minus_tp_rad, mx2_rad):
        base_tuple = basemap.get(k)
        if base_tuple is None:
            continue
        xb_b, mtp_b, mx2_b = base_tuple

        if which == "xb":
            in_rad  = (xr  >= lo) and (xr  <= hi)
            in_base = (xb_b >= lo) and (xb_b <= hi)
            if in_rad and (not in_base):
                acc.append(xb_b)
        elif which == "tprime":
            in_rad  = (mtr  >= lo) and (mtr  <= hi)  # using -t'
            in_base = (mtp_b >= lo) and (mtp_b <= hi)
            if in_rad and (not in_base):
                acc.append(mtp_b)
        elif which == "mx2":
            in_rad  = (m2r  > lo) and (m2r  < hi)    # match your open interval
            in_base = (mx2_b > lo) and (mx2_b < hi)
            if in_rad and (not in_base):
                acc.append(mx2_b)
        else:
            raise ValueError("which must be one of {'xb','tprime','mx2'}")

    return np.array(acc, dtype=float)

def main():
    if len(sys.argv) != 3:
        print("Usage: python isr_fsr_compare_1x3.py <baseline.root> <isr_fsr.root>")
        sys.exit(1)

    base_path, rad_path = sys.argv[1], sys.argv[2]
    if not os.path.isfile(base_path):
        print(f"[ERROR] Baseline file not found: {base_path}")
        sys.exit(1)
    if not os.path.isfile(rad_path):
        print(f"[ERROR] ISR+FSR file not found: {rad_path}")
        sys.exit(1)

    # Load minimal branches
    branches = ["evnnum", "x", "tprime", "Mx2"]
    base = load_branches(base_path, branches)
    rad  = load_branches(rad_path,  branches)

    # Build baseline map for matching
    base_map = make_event_map(base["evnnum"], base["x"], base["tprime"], base["Mx2"])

    # Prepare arrays to histogram
    xb_base = np.asarray(base["x"], dtype=float)
    xb_rad  = np.asarray(rad["x"], dtype=float)

    mtp_base = -np.asarray(base["tprime"], dtype=float)  # -t'
    mtp_rad  = -np.asarray(rad["tprime"], dtype=float)

    mx2_base = np.asarray(base["Mx2"], dtype=float)
    mx2_rad  = np.asarray(rad["Mx2"], dtype=float)

    # Migration selections (baseline values of events that 'moved in' under ISR+FSR)
    xb_migr = migrated_baseline_values(base_map, rad["evnnum"], rad["x"], rad["tprime"], rad["Mx2"],
                                       which="xb",     target=XB_TARGET)
    tp_migr = migrated_baseline_values(base_map, rad["evnnum"], rad["x"], rad["tprime"], rad["Mx2"],
                                       which="tprime", target=TP_TARGET)
    m2_migr = migrated_baseline_values(base_map, rad["evnnum"], rad["x"], rad["tprime"], rad["Mx2"],
                                       which="mx2",    target=MX2_TARGET)

    # Plotting
    os.makedirs(OUT_DIR, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6))
    ax_xb, ax_tp, ax_m2 = axes

    # Shared style
    base_style = dict(histtype="step", linewidth=1.6, color="tab:blue",   label="baseline")
    rad_style  = dict(histtype="step", linewidth=1.6, color="tab:orange", label="sim. ISR+FSR")
    mig_style  = dict(histtype="step", linewidth=1.0, color="gray",       linestyle="--",
                      label="migrated-in (see text)")

    # Left: xB
    ax_xb.hist(xb_base, bins=XB_BINS, range=XB_RANGE, **base_style)
    ax_xb.hist(xb_rad,  bins=XB_BINS, range=XB_RANGE, **rad_style)
    if xb_migr.size:
        ax_xb.hist(xb_migr, bins=XB_BINS, range=XB_RANGE, **mig_style)
    ax_xb.set_xlabel(r"$x_{B}$")
    ax_xb.set_ylabel("counts")
    ax_xb.set_title(r"Low $x_{B}$ bin: $0.10\le x_{B}\le0.25$")
    ax_xb.legend(loc="upper right", frameon=True, edgecolor="black")

    # Middle: -t'
    ax_tp.hist(mtp_base, bins=TP_BINS, range=TP_RANGE, **base_style)
    ax_tp.hist(mtp_rad,  bins=TP_BINS, range=TP_RANGE, **rad_style)
    if tp_migr.size:
        ax_tp.hist(tp_migr, bins=TP_BINS, range=TP_RANGE, **mig_style)
    ax_tp.set_xlabel(r"$-t'\ (\mathrm{GeV}^{2})$")
    ax_tp.set_ylabel("counts")
    ax_tp.set_title(r"Low $-t'$ bin: $0.05\le -t' \le 0.25$")
    ax_tp.legend(loc="upper right", frameon=True, edgecolor="black")

    # Right: Mx2
    ax_m2.hist(mx2_base, bins=MX2_BINS, range=MX2_RANGE, **base_style)
    ax_m2.hist(mx2_rad,  bins=MX2_BINS, range=MX2_RANGE, **rad_style)
    if m2_migr.size:
        ax_m2.hist(m2_migr, bins=MX2_BINS, range=MX2_RANGE, **mig_style)
    ax_m2.set_xlabel(r"$M_{x}^{2}\ (\mathrm{GeV}^{2})$")
    ax_m2.set_ylabel("counts")
    ax_m2.set_title(r"$0.81<M_{x}^{2}<1.00$ migration")
    ax_m2.legend(loc="upper right", frameon=True, edgecolor="black")

    fig.tight_layout()
    fig.savefig(OUT_PATH, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {OUT_PATH}")

if __name__ == "__main__":
    main()