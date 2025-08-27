#!/usr/bin/env python3
"""
Build ISR-proxy angle kernels from e'γ open_angle.

- Input: ROOT file with tree "PhysicsEvents" containing:
    * open_angle  [deg]
    * p_p         [GeV]  (photon momentum, treat as E_gamma)

- Output: prints to stdout Java-ready arrays for 5 Egamma bins:
    [0,1), [1,2), [2,3), [3,4), [4,10] GeV

Usage:
    python make_isr_theta_histos.py \
        --root /volatile/clas12/thayward/egamma/rga_fa18_inb_egamma_short.root \
        --tree PhysicsEvents \
        --theta-max 8.0 \
        --nbins 100 \
        --mode counts        # or pmf or cdf

Notes:
- `counts` prints int[] counts per bin (copy/paste into Java).
- `pmf` prints normalized probability mass over bins (sums to 1).
- `cdf` prints cumulative sum of the PMF (last element = 1).
"""

import argparse
import numpy as np

try:
    import uproot
except ImportError as e:
    print("ERROR: This script requires the 'uproot' package. Try: pip install uproot")
    raise

# ----------------------------------------
def java_array(name, arr, dtype="double"):
    """
    Format a numpy array as a Java array initializer string.
    dtype: "double" or "int"
    """
    if dtype == "int":
        body = ", ".join(str(int(x)) for x in arr.tolist())
        return f"int[] {name} = new int[]{{{body}}};"
    else:
        # use repr with sufficient precision, avoid scientific for small numbers
        body = ", ".join(f"{float(x):.16g}" for x in arr.tolist())
        return f"double[] {name} = new double[]{{{body}}};"
# ----------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="/volatile/clas12/thayward/egamma/rga_fa18_inb_egamma_short.root",
                    help="Path to input ROOT file")
    ap.add_argument("--tree", default="PhysicsEvents", help="TTree name")
    ap.add_argument("--theta-max", type=float, default=8.0, help="max open_angle [deg]")
    ap.add_argument("--nbins", type=int, default=100, help="number of theta bins")
    ap.add_argument("--mode", choices=["counts", "pmf", "cdf"], default="counts",
                    help="what to print for each Egamma bin")
    args = ap.parse_args()

    # Define Egamma bin edges (GeV) and labels for printing
    e_edges = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 10.0], dtype=float)
    e_labels = [f"E{int(e_edges[i])}_{int(e_edges[i+1])}" for i in range(len(e_edges)-1)]

    # Theta binning in degrees (common for all)
    th_edges_deg = np.linspace(0.0, args.theta_max, args.nbins + 1)
    th_width_deg = th_edges_deg[1] - th_edges_deg[0]

    # Load variables
    with uproot.open(args.root) as f:
        if args.tree not in f:
            raise RuntimeError(f"Tree '{args.tree}' not found in file '{args.root}'")
        #endif
        T = f[args.tree]
        arrs = T.arrays(["open_angle", "p_p"], library="np")
        theta_deg = arrs["open_angle"].astype(float)
        egamma   = arrs["p_p"].astype(float)

    # Basic selection: valid angles within [0, theta_max]
    sel = np.isfinite(theta_deg) & np.isfinite(egamma) & (theta_deg >= 0.0) & (theta_deg <= args.theta_max)
    theta_deg = theta_deg[sel]
    egamma    = egamma[sel]

    # Print the common theta-edges so you can hard-code once
    print("// Common theta bin edges (degrees):")
    print(java_array("theta_edges_deg", th_edges_deg, dtype="double"))
    print("")

    # For each Egamma bin, fill histogram
    for i in range(len(e_edges) - 1):
        lo = e_edges[i]
        hi = e_edges[i+1]

        # Last bin inclusive on high edge
        if i == len(e_edges) - 2:
            mask = (egamma >= lo) & (egamma <= hi)
        else:
            mask = (egamma >= lo) & (egamma < hi)
        #endif

        th_bin = theta_deg[mask]
        counts, _ = np.histogram(th_bin, bins=th_edges_deg)

        # Prepare outputs
        label = e_labels[i]  # e.g., "E0_1"
        if args.mode == "counts":
            print(f"// Egamma in [{lo:.1f}, {hi:.1f}{']' if i==len(e_edges)-2 else ')'} GeV], N = {int(mask.sum())}")
            print(java_array(f"h_counts_{label}", counts, dtype="int"))
            print("")
        elif args.mode == "pmf":
            total = counts.sum()
            pmf = counts.astype(float) / total if total > 0 else counts.astype(float)
            # Ensure exact normalization for safety
            if pmf.sum() > 0:
                pmf = pmf / pmf.sum()
            #endif
            print(f"// Egamma in [{lo:.1f}, {hi:.1f}{']' if i==len(e_edges)-2 else ')'} GeV], N = {int(mask.sum())}")
            print(java_array(f"h_pmf_{label}", pmf, dtype="double"))
            print("")
        else:  # cdf
            total = counts.sum()
            pmf = counts.astype(float) / total if total > 0 else counts.astype(float)
            if pmf.sum() > 0:
                pmf = pmf / pmf.sum()
            #endif
            cdf = np.cumsum(pmf)
            # Clip numerical noise
            cdf[-1] = 1.0 if cdf.size > 0 else 0.0
            print(f"// Egamma in [{lo:.1f}, {hi:.1f}{']' if i==len(e_edges)-2 else ')'} GeV], N = {int(mask.sum())}")
            print(java_array(f"h_cdf_{label}", cdf, dtype="double"))
            print("")
        #endfor
# ----------------------------------------

if __name__ == "__main__":
    main()
# endif