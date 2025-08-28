#!/usr/bin/env python3
"""
Build ISR-proxy angle kernels from e'γ open_angle, with a small-angle patch.

- Input TTree ("PhysicsEvents") branches:
    * open_angle  (deg)  -- e'–γ opening angle
    * p_p         (GeV)  -- photon momentum ≈ E_gamma

- Outputs:
    1) Prints to stdout Java-ready arrays (per Eγ bin) for:
         - 'counts' : int[] counts per θ bin (patched)
         - 'pmf'    : double[] normalized probability mass (patched)
         - 'cdf'    : double[] cumulative (last element = 1 if N>0, patched)
       using common θ-bin edges.

    2) Saves validation plots per Eγ bin to --outdir (default: output/validation_inverseISR/):
         - {label}_overlay_full.pdf/png : 0–8° measured vs patched + fit
         - {label}_overlay_zoom.pdf/png : 0–1° zoom overlay
         - {label}_ratio_zoom.pdf/png   : 0–1° ratio (patched/measured)
         - {label}_cdf_compare.pdf/png  : CDF before vs after

Small-angle patching:
    * Fit density y(θ) = A * (θ + θ0)^(-p) on a "trust" region (default 1–8°).
    * θ0 chosen by a short grid search (default 0.02,0.05,0.10 deg) to minimize
      weighted SSE in log–log space.
    * Predict sub-degree (θ < --patch-below, default 1°) bin contents via bin-wise integrals
      of the fitted function; optionally enforce "only raise" (default True).
    * Above the patch boundary, measured bins are left unchanged.

Usage:
  python determine_egamma_open_angle_CDF.py \
      --root /volatile/clas12/thayward/egamma/rga_fa18_inb_egamma.root \
      --tree PhysicsEvents \
      --theta-max 8.0 \
      --nbins 100 \
      --mode cdf \
      --outdir output/validation_inverseISR

Notes:
  - Doubles printed with four decimals.
  - Last Eγ bin is inclusive on the high edge.
"""

import argparse
import os
from math import exp, log
import numpy as np
import matplotlib.pyplot as plt

try:
    import uproot
except ImportError:
    print("ERROR: This script requires the 'uproot' package. Try: pip install uproot")
    raise

# ------------------------------ Java formatting helpers ------------------------------
def java_array(name, arr, dtype="double"):
    """
    Format a numpy array as a Java array initializer string.
    dtype: "double" or "int"
    Doubles are formatted to four decimals.
    """
    if dtype == "int":
        body = ", ".join(str(int(x)) for x in np.asarray(arr).tolist())
        return f"int[] {name} = new int[]{{{body}}};"
    else:
        body = ", ".join(f"{float(x):.4f}" for x in np.asarray(arr).tolist())
        return f"double[] {name} = new double[]{{{body}}};"

def edge_label(v):
    """Make a Java-friendly label from a float (e.g., 0.25 -> '0p25', 3.0 -> '3')."""
    s = f"{v:.2f}".rstrip('0').rstrip('.')
    return s.replace('.', 'p')

# ------------------------------ Fitting / patching -----------------------------------
def _log_polyfit(x, y, w=None):
    """
    Weighted polyfit of log(y) vs log(x). Returns (A, p) for y = A * x^(-p),
    where A = exp(intercept), p = -slope.
    Only finite positive (x,y) contribute.
    """
    x = np.asarray(x); y = np.asarray(y)
    m = (x > 0) & (y > 0) & np.isfinite(x) & np.isfinite(y)
    if w is None:
        w = np.ones_like(x)
    else:
        w = np.asarray(w)
    m &= np.isfinite(w) & (w > 0)
    if np.count_nonzero(m) < 2:
        return None  # not enough points
    lx = np.log(x[m])
    ly = np.log(y[m])
    ww = np.sqrt(w[m])  # np.polyfit expects weights ~ 1/sigma; here we want ~sqrt(counts)
    try:
        coeffs = np.polyfit(lx, ly, deg=1, w=ww)
        slope, intercept = coeffs[0], coeffs[1]
        A = exp(intercept)
        p = -slope
        return A, p
    except Exception:
        return None

def _integral_powerlaw(A, p, t0, tL, tH):
    """∫_{tL}^{tH} A (t + t0)^(-p) dt, safe for p≈1."""
    if p == 1.0:
        return A * np.log((tH + t0) / (tL + t0))
    else:
        q = 1.0 - p
        return A / q * ((tH + t0)**q - (tL + t0)**q)

def fit_and_patch(theta_edges, counts, patch_below_deg=1.0,
                  fit_min=1.0, fit_max=8.0,
                  theta0_grid=(0.02, 0.05, 0.10),
                  only_raise=True):
    """
    Fit density in [fit_min, fit_max] and patch bins with upper edge <= patch_below_deg.
    Returns: counts_patched, fit_info (dict).
    """
    counts = np.asarray(counts, dtype=float)
    edges = np.asarray(theta_edges, dtype=float)
    centers = 0.5*(edges[:-1] + edges[1:])
    dtheta = np.diff(edges)
    assert counts.size == dtheta.size

    # Density (counts per degree) for fitting
    with np.errstate(divide="ignore", invalid="ignore"):
        density = counts / dtheta
    # Select trust region bins
    m_fit = (centers >= fit_min) & (centers <= fit_max) & (density > 0)
    if np.count_nonzero(m_fit) < 5:
        # Too few bins to fit; return unmodified
        return counts.copy(), {
            "A": np.nan, "p": np.nan, "theta0": np.nan,
            "chi2": np.nan, "ndf": 0, "used_fit": False
        }

    x_fit = centers[m_fit]
    y_fit = density[m_fit]
    w_fit = counts[m_fit].clip(min=1.0)  # Poisson weights ~ counts

    best = None
    for t0 in theta0_grid:
        params = _log_polyfit(x_fit + t0, y_fit, w_fit)
        if params is None:
            continue
        A, p = params
        # Predict density on fit region and compute weighted SSE in linear space
        y_pred = A * np.power(x_fit + t0, -p)
        # Variance of density ~ counts / dtheta^2; use chi2-like metric
        sigma2 = (counts[m_fit] / (dtheta[m_fit]**2)).clip(min=1e-12)
        chi2 = np.sum((y_fit - y_pred)**2 / sigma2)
        ndf = max(1, np.count_nonzero(m_fit) - 2)
        score = chi2 / ndf
        if (best is None) or (score < best["score"]):
            best = {"A": A, "p": p, "theta0": t0, "chi2": chi2, "ndf": ndf, "score": score}

    # Guardrails / fallback
    if best is None or not np.isfinite(best["A"]) or not np.isfinite(best["p"]) or best["A"] <= 0 or best["p"] <= 0:
        # Conservative fallback: p=1, choose A by area matching over fit region
        t0 = 0.05
        area_data = np.sum(counts[m_fit])
        # area_model = ∫ A/(θ+t0)^1 over fit bins; solve A so areas match
        area_unit = 0.0
        for tL, tH in zip(edges[:-1][m_fit], edges[1:][m_fit]):
            area_unit += _integral_powerlaw(1.0, 1.0, t0, tL, tH)
        A = (area_unit > 0) and (area_data / area_unit) or 1.0
        best = {"A": A, "p": 1.0, "theta0": t0, "chi2": np.nan, "ndf": max(1, np.count_nonzero(m_fit)-1), "score": np.nan}

    # Patch sub-degree bins via integrals
    counts_patched = counts.copy()
    m_patch = edges[1:] <= patch_below_deg + 1e-12
    if np.any(m_patch):
        A, p, t0 = best["A"], best["p"], best["theta0"]
        pred = np.zeros_like(counts_patched)
        for i, (tL, tH) in enumerate(zip(edges[:-1], edges[1:])):
            if not m_patch[i]:
                continue
            pred[i] = _integral_powerlaw(A, p, t0, tL, tH)
        if only_raise:
            counts_patched[m_patch] = np.maximum(counts_patched[m_patch], pred[m_patch])
        else:
            counts_patched[m_patch] = pred[m_patch]

    best["used_fit"] = True
    return counts_patched, best

# ------------------------------ Plotting helpers -------------------------------------
def _safe_mkdir(d):
    os.makedirs(d, exist_ok=True)

def plot_validation(bin_label, e_range_str,
                    theta_edges, counts_raw, counts_patched,
                    fit_info, outdir, save_png=True):
    """
    Make four plots:
      - overlay_full (0–8°)
      - overlay_zoom  (0–1°)
      - ratio_zoom    (patched/measured in 0–1°)
      - cdf_compare   (0–8° CDF)
    """
    _safe_mkdir(outdir)
    edges = theta_edges
    centers = 0.5*(edges[:-1] + edges[1:])
    dtheta = np.diff(edges)

    # Fit curve for display
    A, p, t0 = fit_info["A"], fit_info["p"], fit_info["theta0"]
    have_fit = np.isfinite(A) and np.isfinite(p) and np.isfinite(t0) and fit_info.get("used_fit", False)
    theta_dense = np.linspace(edges[0], edges[-1], 800)
    curve_dense = None
    if have_fit:
        curve_dense = A * np.power(theta_dense + t0, -p)

    # ---------- overlay_full ----------
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.step(edges[:-1], counts_raw, where="post", alpha=0.7, label="measured (raw)")
    ax.step(edges[:-1], counts_patched, where="post", alpha=0.9, label="patched", linewidth=1.5)
    if have_fit:
        ax2 = ax.twinx()
        ax2.plot(theta_dense, curve_dense, linestyle='--', alpha=0.7, label="fit density")
        ax2.set_ylabel("f(θ) (counts/deg)", rotation=270, labelpad=15)
        ax2.grid(False)
        # combine legends
        lines, labels = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax.legend(lines + lines2, labels + labels2, loc="upper right")
    else:
        ax.legend(loc="upper right")
    ax.set_title(f"open_angle {e_range_str}  (full range)")
    ax.set_xlabel(r"$\theta_{e'\gamma}$ (deg)")
    ax.set_ylabel("counts")
    txt = (f"fit: A={A:.3g}, p={p:.3g}, θ0={t0:.3g}°\n"
           f"χ²/ndf={fit_info['chi2']:.2f}/{fit_info['ndf']} = "
           f"{(fit_info['chi2']/fit_info['ndf']):.2f}" if np.isfinite(fit_info["chi2"]) else
           f"fit: A={A:.3g}, p={p:.3g}, θ0={t0:.3g}°")
    ax.text(0.02, 0.98, txt, transform=ax.transAxes, va="top", ha="left",
            fontsize=9, bbox=dict(boxstyle="round", facecolor="w", alpha=0.7))
    ax.grid(alpha=0.25)
    base = os.path.join(outdir, f"{bin_label}_overlay_full")
    fig.tight_layout()
    fig.savefig(base + ".pdf")
    if save_png: fig.savefig(base + ".png", dpi=160)
    plt.close(fig)

    # ---------- overlay_zoom (0–1°) ----------
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.step(edges[:-1], counts_raw, where="post", alpha=0.7, label="measured (raw)")
    ax.step(edges[:-1], counts_patched, where="post", alpha=0.9, label="patched", linewidth=1.5)
    if have_fit:
        ax2 = ax.twinx()
        ax2.plot(theta_dense, curve_dense, linestyle='--', alpha=0.7, label="fit density")
        ax2.set_ylabel("f(θ) (counts/deg)", rotation=270, labelpad=15)
        ax2.grid(False)
        lines, labels = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax.legend(lines + lines2, labels + labels2, loc="upper right")
    else:
        ax.legend(loc="upper right")
    ax.set_xlim(0.0, 1.0)
    ax.set_title(f"open_angle {e_range_str}  (zoom: 0–1°)")
    ax.set_xlabel(r"$\theta_{e'\gamma}$ (deg)")
    ax.set_ylabel("counts")
    ax.grid(alpha=0.25)
    base = os.path.join(outdir, f"{bin_label}_overlay_zoom")
    fig.tight_layout()
    fig.savefig(base + ".pdf")
    if save_png: fig.savefig(base + ".png", dpi=160)
    plt.close(fig)

    # ---------- ratio_zoom (patched / measured) ----------
    fig, ax = plt.subplots(figsize=(7, 3.6))
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(counts_raw > 0, counts_patched / counts_raw, np.nan)
    ax.step(edges[:-1], ratio, where="post")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(bottom=0)
    ax.set_title(f"patched / measured  (zoom: 0–1°)  {e_range_str}")
    ax.set_xlabel(r"$\theta_{e'\gamma}$ (deg)")
    ax.set_ylabel("ratio")
    ax.grid(alpha=0.25)
    base = os.path.join(outdir, f"{bin_label}_ratio_zoom")
    fig.tight_layout()
    fig.savefig(base + ".pdf")
    if save_png: fig.savefig(base + ".png", dpi=160)
    plt.close(fig)

    # ---------- CDF compare ----------
    def _cdf_from_counts(cnt):
        tot = float(np.sum(cnt))
        if tot <= 0:
            return np.zeros_like(cnt, dtype=float)
        pmf = cnt / tot
        cdf = np.cumsum(pmf)
        cdf = np.clip(cdf, 0.0, 1.0)
        if cdf.size: cdf[-1] = 1.0
        return cdf

    cdf_raw = _cdf_from_counts(counts_raw)
    cdf_pat = _cdf_from_counts(counts_patched)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.step(edges[:-1], cdf_raw, where="post", label="CDF (raw)")
    ax.step(edges[:-1], cdf_pat, where="post", label="CDF (patched)")
    ax.set_xlim(0.0, edges[-1])
    ax.set_ylim(0.0, 1.0)
    ax.set_title(f"CDF compare  {e_range_str}")
    ax.set_xlabel(r"$\theta_{e'\gamma}$ (deg)")
    ax.set_ylabel("CDF")
    ax.legend(loc="lower right")
    ax.grid(alpha=0.25)
    base = os.path.join(outdir, f"{bin_label}_cdf_compare")
    fig.tight_layout()
    fig.savefig(base + ".pdf")
    if save_png: fig.savefig(base + ".png", dpi=160)
    plt.close(fig)

# ------------------------------ Main -----------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, help="Path to input ROOT file")
    ap.add_argument("--tree", default="PhysicsEvents", help="TTree name")
    ap.add_argument("--theta-max", type=float, default=8.0, help="max open_angle (deg)")
    ap.add_argument("--nbins", type=int, default=100, help="number of theta bins")
    ap.add_argument("--mode", choices=["counts", "pmf", "cdf"], default="counts",
                    help="what to print for each Egamma bin (patched)")
    ap.add_argument("--outdir", default="output/validation_inverseISR", help="directory for validation plots")

    # Patching / fitting controls
    ap.add_argument("--patch-below", type=float, default=1.0,
                    help="replace bins with upper edge <= this (deg) using fit")
    ap.add_argument("--fit-min", type=float, default=1.0, help="fit lower bound (deg)")
    ap.add_argument("--fit-max", type=float, default=8.0, help="fit upper bound (deg)")
    ap.add_argument("--theta0-grid", type=str, default="0.02,0.05,0.10",
                    help="comma-separated θ0 (deg) values to try (e.g. '0.02,0.05,0.10')")
    ap.add_argument("--allow-lower", action="store_true",
                    help="if set, patch may lower counts below measured in 0–patch-below; default is only-raise")
    ap.add_argument("--no-plots", action="store_true", help="disable saving validation plots")

    args = ap.parse_args()

    theta0_grid = tuple(float(x.strip()) for x in args.theta0_grid.split(","))

    # E_gamma bins (GeV); last bin inclusive
    e_edges = np.array(
        [0.0, 0.25, 0.50, 0.75, 1.00, 1.50, 2.00, 2.50, 3.0, 3.5, 4.0, 4.5, 5.0, 10.0],
        dtype=float
    )
    e_labels = [f"E{edge_label(e_edges[i])}_{edge_label(e_edges[i+1])}"
                for i in range(len(e_edges)-1)]

    # Common theta bins (deg)
    th_edges_deg = np.linspace(0.0, args.theta_max, args.nbins + 1)

    # Load variables
    with uproot.open(args.root) as f:
        if args.tree not in f:
            raise RuntimeError(f"Tree '{args.tree}' not found in file '{args.root}'")
        T = f[args.tree]
        arrs = T.arrays(["open_angle", "p_p"], library="np")
        theta_deg = np.asarray(arrs["open_angle"], dtype=float)
        egamma    = np.asarray(arrs["p_p"], dtype=float)

    # Basic selection: valid and within [0, theta_max]
    sel = np.isfinite(theta_deg) & np.isfinite(egamma) & (theta_deg >= 0.0) & (theta_deg <= args.theta_max)
    theta_deg = theta_deg[sel]
    egamma    = egamma[sel]

    # Print the common theta-edges once
    print("// Common theta bin edges (degrees):")
    print(java_array("theta_edges_deg", th_edges_deg, dtype="double"))
    print("")

    # Prepare output directory for plots
    if not args.no_plots:
        os.makedirs(args.outdir, exist_ok=True)

    # Helper to emit arrays from counts
    def emit_arrays(label, counts, lo, hi, inclusive):
        range_str = f"[{lo:.4f}, {hi:.4f}{']' if inclusive else ')'} GeV]"
        N = int(np.sum(counts))
        total = float(np.sum(counts))
        if args.mode == "counts":
            print(f"// Egamma in {range_str}, N = {N}")
            print(java_array(f"h_counts_{label}", counts.astype(int), dtype="int"))
            print("")
        elif args.mode == "pmf":
            if total > 0:
                pmf = counts / total
                pmf = np.clip(pmf, 0.0, 1.0)
                pmf = pmf / pmf.sum()
            else:
                pmf = np.zeros_like(counts, dtype=float)
            print(f"// Egamma in {range_str}, N = {N}")
            print(java_array(f"h_pmf_{label}", pmf, dtype="double"))
            print("")
        else:  # cdf
            if total > 0:
                pmf = counts / total
                pmf = np.clip(pmf, 0.0, 1.0)
                pmf = pmf / pmf.sum()
                cdf = np.cumsum(pmf)
                cdf = np.clip(cdf, 0.0, 1.0)
                cdf[-1] = 1.0
            else:
                cdf = np.zeros_like(counts, dtype=float)
            print(f"// Egamma in {range_str}, N = {N}")
            print(java_array(f"h_cdf_{label}", cdf, dtype="double"))
            print("")

    # Loop over Eγ bins
    last = len(e_edges) - 1
    for i in range(last):
        lo, hi = e_edges[i], e_edges[i+1]
        inclusive = (i == last - 1)
        if inclusive:
            mask = (egamma >= lo) & (egamma <= hi)
        else:
            mask = (egamma >= lo) & (egamma < hi)

        th_bin = theta_deg[mask]
        counts_raw, _ = np.histogram(th_bin, bins=th_edges_deg)

        # Patch small angles
        counts_patched, fit_info = fit_and_patch(
            th_edges_deg, counts_raw,
            patch_below_deg=args.patch_below,
            fit_min=args.fit_min, fit_max=args.fit_max,
            theta0_grid=theta0_grid,
            only_raise=(not args.allow_lower)
        )

        # Plots
        if not args.no_plots:
            plot_validation(
                bin_label=e_labels[i],
                e_range_str=f"[{lo:.2f}, {hi:.2f}{']' if inclusive else ')'} GeV]",
                theta_edges=th_edges_deg,
                counts_raw=counts_raw,
                counts_patched=counts_patched,
                fit_info=fit_info,
                outdir=args.outdir
            )

        # Emit arrays (PATCHED histogram drives PMF/CDF)
        emit_arrays(e_labels[i], counts_patched, lo, hi, inclusive)

# -------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()