#!/usr/bin/env python3
"""
backsolve_photon_efficiency_corrections.py

Exploratory diagnostic only.

Infer constant FD and FT photon cross-section correction multipliers from the
pass-2 cross-section CSV by asking what two constants best bring the measured
cross sections into agreement with an external reference.

The fitted model is

    sigma_i(corrected)
      = sigma_i(raw) * [ f_FD,i * C_FD + f_FT,i * C_FT ],

where the topology yields are formed from the four unpolarized 10.6-GeV
periods (Fa18 Inb/Out and Sp18 Inb/Out):

    f_FD,i = (Y_FD,FD + Y_CD,FD) / (Y_FD,FD + Y_CD,FD + Y_CD,FT)
    f_FT,i = Y_CD,FT             / (Y_FD,FD + Y_CD,FD + Y_CD,FT).

Sp19 Inb is not included in these mixture weights because it is the
10.2-GeV period.

C_FD and C_FT are CROSS-SECTION multipliers = epsilon_MC / epsilon_data.

The corresponding efficiency ratios are therefore

    epsilon_data / epsilon_MC = 1 / C.

Two reference modes are supported:

  1. KM15:
       use the pass-2 KM15 values saved by compare_world_dvcs_cross_sections.py.

  2. Lee:
       use exact same-bin Lee 2026 matches saved by the same comparison script.
       Lee is transported from the Lee mean kinematics to the Hayward mean
       kinematics with the saved KM15 transport factor before fitting.

This is NOT a production efficiency extraction.  It is deliberately a
"what correction would the cross sections prefer?" diagnostic.

In addition to the two-parameter global fit, the script now reports the direct
pointwise correction

    C_i = sigma_reference / sigma_Hayward_raw

for FD-pure (default f_FT <= 0.05), FT-pure (default f_FT >= 0.95), mixed,
and all points.  The median and central 68% interval are unweighted so that a
few tiny-error points cannot dictate the answer.  Where both measurements have
point uncertainties (notably Lee), a one-parameter weighted scale using the
combined Hayward+reference uncertainty is also reported.

Example:

  python external_scripts/backsolve_photon_efficiency_corrections.py \
      --pass2-file output/csvs/dvcs_pass2_analysis.csv \
      --comparison-dir ../output/world_dvcs_cross_section_comparison_photon_eff_prelim \
      --outdir ../output/photon_efficiency_backsolve

"""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd


XS_COL = "normed cross sections, ep->epg, exp, 10.6 GeV, unpol"


def parse_args():
    p = argparse.ArgumentParser(
        description="Back-solve constant FD/FT photon correction factors from pass-2 cross sections."
    )
    p.add_argument(
        "--pass2-file",
        default="output/csvs/dvcs_pass2_analysis.csv",
        help="UNCORRECTED pass-2 cross-section CSV.",
    )
    p.add_argument(
        "--comparison-dir",
        default="../output/world_dvcs_cross_section_comparison_photon_eff_prelim",
        help="Output directory from compare_world_dvcs_cross_sections.py.",
    )
    p.add_argument(
        "--outdir",
        default="../output/photon_efficiency_backsolve",
        help="Directory for diagnostic CSV outputs.",
    )
    p.add_argument(
        "--reference",
        choices=("both", "km15", "lee"),
        default="both",
        help="Reference used to back-solve the correction.",
    )
    p.add_argument(
        "--bh-min",
        type=float,
        default=0.0,
        help=(
            "Optional minimum BH/KM15 fraction for the KM15 fit. "
            "Use e.g. 0.8 or 0.9 to emphasize BH-dominated points."
        ),
    )
    p.add_argument(
        "--phi-edge-max",
        type=float,
        default=None,
        help=(
            "Optional |phi_edge| maximum in degrees for KM15, where "
            "phi_edge=min(phi,360-phi). Example: 40 selects phi near 0/360."
        ),
    )
    p.add_argument(
        "--no-point-weights",
        action="store_true",
        help="Fit all points with equal weight instead of point-by-point uncertainties.",
    )
    p.add_argument(
        "--pure-threshold",
        type=float,
        default=0.05,
        help=(
            "Topology-purity threshold. FD-pure means f_FT <= threshold and "
            "FT-pure means f_FT >= 1-threshold. Default: 0.05."
        ),
    )
    p.add_argument(
        "--bootstrap",
        type=int,
        default=1000,
        help="Number of row-bootstrap replicas for diagnostic factor uncertainties; 0 disables.",
    )
    p.add_argument(
        "--seed",
        type=int,
        default=12345,
        help="Random seed for bootstrap.",
    )
    return p.parse_args()


def parse_tuple_first(x):
    """Return the first numeric component of scalar or '(value,...)' CSV cells."""
    if x is None:
        return np.nan
    if isinstance(x, (int, float, np.integer, np.floating)):
        return float(x) if np.isfinite(x) else np.nan

    s = str(x).strip()
    if not s or s.lower() in {"nan", "none"}:
        return np.nan

    if s.startswith("(") and s.endswith(")"):
        s = s[1:-1]
    #endif

    first = s.split(",")[0].strip()
    try:
        return float(first)
    except Exception:
        return np.nan
    #endtry


def parse_tuple_second(x):
    """Return the second numeric component of '(value,error,...)' cells."""
    if x is None:
        return np.nan

    s = str(x).strip()
    if s.startswith("(") and s.endswith(")"):
        s = s[1:-1]
    #endif

    parts = [q.strip() for q in s.split(",")]
    if len(parts) < 2:
        return np.nan
    #endif

    try:
        return float(parts[1])
    except Exception:
        return np.nan
    #endtry


def norm_name(s):
    return re.sub(r"\s+", " ", str(s).strip().lower())


def resolve_comparison_file(comparison_dir, filename):
    """
    Find a named comparison-output file robustly.

    Normally files live in:
        <comparison_dir>/tables/<filename>

    But earlier runs in this analysis accidentally used a nested cache/
    directory as --outdir.  To make this diagnostic insensitive to that
    bookkeeping, first try the canonical location and then search recursively.
    """
    comparison_dir = Path(comparison_dir)

    preferred = comparison_dir / "tables" / filename
    if preferred.exists():
        return preferred
    #endif

    matches = list(comparison_dir.rglob(filename))
    if len(matches) == 1:
        print(f"[comparison file] using nested path: {matches[0]}")
        return matches[0]
    #endif

    if len(matches) > 1:
        # Prefer the shallowest path, then newest modification time.
        matches = sorted(
            matches,
            key=lambda p: (len(p.relative_to(comparison_dir).parts), -p.stat().st_mtime),
        )
        print("[comparison file] multiple matches found:")
        for p in matches:
            print(f"  {p}")
        #endfor
        print(f"[comparison file] choosing: {matches[0]}")
        return matches[0]
    #endif

    # Give a useful directory diagnostic before failing.
    nearby = []
    if comparison_dir.exists():
        for p in comparison_dir.rglob("*.csv"):
            nearby.append(p)
            if len(nearby) >= 40:
                break
            #endif
        #endfor
    #endif

    msg = [
        f"Could not find '{filename}' beneath comparison directory:",
        f"  {comparison_dir}",
    ]
    if nearby:
        msg.append("CSV files found beneath that directory:")
        msg.extend(f"  {p}" for p in nearby)
    else:
        msg.append("No CSV files were found beneath that directory.")
    #endif

    raise FileNotFoundError("\n".join(msg))



TENP6_PERIODS = ("Fa18 Inb", "Fa18 Out", "Sp18 Inb", "Sp18 Out")


def find_topology_yield_columns(df):
    """
    Locate the UNPOLARIZED normalized raw-yield columns for the four 10.6-GeV
    periods separately, for:
      (FD, FD), (CD, FD), (CD, FT)

    The CSV does not contain a single "10.6 GeV" topology-yield column.  The
    combined 10.6-GeV cross section is built from Fa18 Inb, Fa18 Out,
    Sp18 Inb, and Sp18 Out, while Sp19 Inb is the 10.2-GeV period.

    For the sole purpose of constructing the FD/FT topology MIXTURE, this
    script sums the normalized raw yields over the four 10.6-GeV periods.
    A common factor of four would cancel in the topology fractions, so this is
    equivalent to taking the four-period mean normalized raw yield.

    This is a diagnostic mixture estimate, not a redefinition of the published
    10.6-GeV cross-section combination.
    """
    candidates = [
        c for c in df.columns
        if "normalized raw yield" in norm_name(c)
        and "ep->epg" in norm_name(c)
        and norm_name(c).endswith(", unpol")
    ]

    topology_text = {
        "FD_FD": "(fd,fd)",
        "CD_FD": "(cd,fd)",
        "CD_FT": "(cd,ft)",
    }

    out = {key: [] for key in topology_text}

    for key, topo in topology_text.items():
        for period in TENP6_PERIODS:
            period_norm = norm_name(period).replace(" ", "")
            matches = []

            for c in candidates:
                n = norm_name(c).replace(" ", "")
                if topo in n and period_norm in n:
                    matches.append(c)
                #endif
            #endfor

            if len(matches) != 1:
                raise RuntimeError(
                    "Could not uniquely locate the requested normalized raw-yield "
                    f"column for topology {key}, period '{period}'.\n"
                    f"Matches: {matches}\n"
                    "Available unpolarized normalized raw-yield columns were:\n  "
                    + "\n  ".join(candidates)
                )
            #endif

            out[key].append(matches[0])
        #endfor
    #endfor

    return out


def build_pass2_dataframe(path):
    df = pd.read_csv(path, low_memory=False).copy()

    if XS_COL not in df.columns:
        raise RuntimeError(f"Required cross-section column not found:\n  {XS_COL}")
    #endif

    ycols = find_topology_yield_columns(df)

    df["_xs_raw"] = df[XS_COL].map(parse_tuple_first)
    df["_xs_stat_raw"] = df[XS_COL].map(parse_tuple_second)

    # Build each topology's effective 10.6-GeV normalized raw yield by summing
    # the four same-energy periods.  Using the mean instead would give exactly
    # the same FD/FT fractions.
    for key, cols in ycols.items():
        pieces = []
        for col in cols:
            pieces.append(df[col].map(parse_tuple_first))
        #endfor

        tmp = pd.concat(pieces, axis=1)
        n_valid = tmp.notna().sum(axis=1)

        # min_count=1 prevents "all missing" from silently becoming zero.
        summed = tmp.sum(axis=1, min_count=1)

        # Require at least one contributing 10.6-GeV period.  Normally all four
        # are present for populated bins; sparse bins can legitimately have fewer.
        df[f"_yield_{key}"] = np.where(n_valid > 0, summed, np.nan)
        df[f"_nperiod_{key}"] = n_valid
    #endfor

    df["_yield_FD"] = (
        df["_yield_FD_FD"].fillna(0.0)
        + df["_yield_CD_FD"].fillna(0.0)
    )
    df["_yield_FT"] = df["_yield_CD_FT"].fillna(0.0)
    df["_yield_total"] = df["_yield_FD"] + df["_yield_FT"]

    good_y = np.isfinite(df["_yield_total"]) & (df["_yield_total"] > 0.0)
    df["_f_FD"] = np.where(
        good_y,
        df["_yield_FD"] / df["_yield_total"],
        np.nan,
    )
    df["_f_FT"] = np.where(
        good_y,
        df["_yield_FT"] / df["_yield_total"],
        np.nan,
    )

    print("[columns] 10.6-GeV unpolarized topology yields:")
    for key, cols in ycols.items():
        print(f"  {key}:")
        for period, col in zip(TENP6_PERIODS, cols):
            print(f"    {period:9s}: {col}")
        #endfor
    #endfor

    finite_mix = np.isfinite(df["_f_FT"])
    if finite_mix.any():
        print(
            "[topology mixture] FT fraction: "
            f"min={df.loc[finite_mix, '_f_FT'].min():.4f}, "
            f"median={df.loc[finite_mix, '_f_FT'].median():.4f}, "
            f"max={df.loc[finite_mix, '_f_FT'].max():.4f}"
        )
    #endif

    return df



def solve_linear(df, target_col, unc_col, use_weights=True):
    """
    Weighted linear least squares for:
      target = xs_raw * fFD * CFD + xs_raw * fFT * CFT
    """
    x1 = df["_xs_raw"].to_numpy(float) * df["_f_FD"].to_numpy(float)
    x2 = df["_xs_raw"].to_numpy(float) * df["_f_FT"].to_numpy(float)
    y = df[target_col].to_numpy(float)

    good = np.isfinite(x1) & np.isfinite(x2) & np.isfinite(y) & (y > 0)

    if use_weights:
        u = df[unc_col].to_numpy(float)
        good &= np.isfinite(u) & (u > 0)
        w = np.zeros_like(y)
        w[good] = 1.0 / (u[good] ** 2)
    else:
        w = np.ones_like(y)
    #endif

    X = np.column_stack([x1[good], x2[good]])
    yy = y[good]
    ww = w[good]

    if X.shape[0] < 3:
        raise RuntimeError("Too few valid points for a two-parameter FD/FT fit.")
    #endif

    sw = np.sqrt(ww)
    Xw = X * sw[:, None]
    yw = yy * sw

    beta, _, _, _ = np.linalg.lstsq(Xw, yw, rcond=None)
    c_fd, c_ft = map(float, beta)

    pred = X @ beta
    resid = yy - pred

    if use_weights:
        chi2 = float(np.sum((resid * sw) ** 2))
    else:
        # diagnostic only: use unweighted squared fractional residual scale
        chi2 = float(np.sum(resid ** 2))
    #endif

    n = len(yy)
    ndof = max(n - 2, 1)

    # Formal covariance for the specified point uncertainties.
    xtwx = Xw.T @ Xw
    cov = np.linalg.pinv(xtwx)
    if not use_weights:
        cov *= chi2 / ndof
    #endif

    err_fd = float(math.sqrt(max(cov[0, 0], 0.0)))
    err_ft = float(math.sqrt(max(cov[1, 1], 0.0)))
    corr = (
        float(cov[0, 1] / math.sqrt(cov[0, 0] * cov[1, 1]))
        if cov[0, 0] > 0 and cov[1, 1] > 0
        else np.nan
    )

    result = {
        "N": n,
        "C_FD": c_fd,
        "C_FT": c_ft,
        "formal_err_C_FD": err_fd,
        "formal_err_C_FT": err_ft,
        "corr_CFD_CFT": corr,
        "eps_data_over_mc_FD": 1.0 / c_fd if c_fd > 0 else np.nan,
        "eps_data_over_mc_FT": 1.0 / c_ft if c_ft > 0 else np.nan,
        "chi2": chi2,
        "ndof": ndof,
        "chi2_per_dof": chi2 / ndof,
    }

    idx = df.index.to_numpy()[good]
    return result, idx, pred, resid



def solve_linear_combined_uncertainty(
    df,
    target_col,
    target_unc_col,
    hayward_unc_col,
    max_iter=100,
    tol=1.0e-10,
):
    """
    Two-parameter FD/FT fit using BOTH reference and Hayward point uncertainties.

    The model is
        target_i = xs_i * (fFD_i*C_FD + fFT_i*C_FT)

    and the variance used at each iteration is
        var_i = target_unc_i^2 + (Ceff_i * hayward_unc_i)^2.

    Because Ceff depends on the fitted parameters, the weights are updated
    iteratively until the two correction factors stop changing.
    """
    xraw = df["_xs_raw"].to_numpy(float)
    ffd = df["_f_FD"].to_numpy(float)
    fft = df["_f_FT"].to_numpy(float)
    y = df[target_col].to_numpy(float)
    uy = df[target_unc_col].to_numpy(float)
    ux = df[hayward_unc_col].to_numpy(float)

    good = (
        np.isfinite(xraw)
        & np.isfinite(ffd)
        & np.isfinite(fft)
        & np.isfinite(y)
        & np.isfinite(uy)
        & np.isfinite(ux)
        & (xraw > 0)
        & (y > 0)
        & (uy > 0)
        & (ux >= 0)
    )

    X = np.column_stack(
        [
            xraw[good] * ffd[good],
            xraw[good] * fft[good],
        ]
    )
    yy = y[good]
    uyy = uy[good]
    uxx = ux[good]
    ffd_g = ffd[good]
    fft_g = fft[good]

    if X.shape[0] < 3:
        raise RuntimeError(
            "Too few valid points for combined-uncertainty FD/FT fit."
        )
    #endif

    # Start from the reference-only weighted solution.
    start_df = df.loc[df.index.to_numpy()[good]].copy()
    start, _, _, _ = solve_linear(
        start_df,
        target_col,
        target_unc_col,
        use_weights=True,
    )
    beta = np.array([start["C_FD"], start["C_FT"]], dtype=float)

    converged = False
    niter = 0
    for niter in range(1, max_iter + 1):
        ceff = ffd_g * beta[0] + fft_g * beta[1]
        var = uyy**2 + (ceff * uxx)**2
        valid_var = np.isfinite(var) & (var > 0)

        Xv = X[valid_var]
        yv = yy[valid_var]
        vv = var[valid_var]

        sw = 1.0 / np.sqrt(vv)
        Xw = Xv * sw[:, None]
        yw = yv * sw

        beta_new, _, _, _ = np.linalg.lstsq(Xw, yw, rcond=None)

        if np.max(np.abs(beta_new - beta)) < tol:
            beta = beta_new
            converged = True
            break
        #endif

        beta = beta_new
    #endfor

    ceff = ffd_g * beta[0] + fft_g * beta[1]
    var = uyy**2 + (ceff * uxx)**2
    valid_var = np.isfinite(var) & (var > 0)

    Xv = X[valid_var]
    yv = yy[valid_var]
    vv = var[valid_var]
    pred = Xv @ beta
    resid = yv - pred

    chi2 = float(np.sum(resid**2 / vv))
    n = len(yv)
    ndof = max(n - 2, 1)

    sw = 1.0 / np.sqrt(vv)
    Xw = Xv * sw[:, None]
    cov = np.linalg.pinv(Xw.T @ Xw)

    c_fd, c_ft = map(float, beta)
    efd = float(math.sqrt(max(cov[0, 0], 0.0)))
    eft = float(math.sqrt(max(cov[1, 1], 0.0)))
    corr = (
        float(cov[0, 1] / math.sqrt(cov[0, 0] * cov[1, 1]))
        if cov[0, 0] > 0 and cov[1, 1] > 0
        else np.nan
    )

    original_idx = df.index.to_numpy()[good][valid_var]

    return {
        "N": n,
        "C_FD": c_fd,
        "C_FT": c_ft,
        "formal_err_C_FD": efd,
        "formal_err_C_FT": eft,
        "corr_CFD_CFT": corr,
        "eps_data_over_mc_FD": 1.0 / c_fd if c_fd > 0 else np.nan,
        "eps_data_over_mc_FT": 1.0 / c_ft if c_ft > 0 else np.nan,
        "chi2": chi2,
        "ndof": ndof,
        "chi2_per_dof": chi2 / ndof,
        "iterations": niter,
        "converged": converged,
    }, original_idx, pred, resid


def solve_one_scale_combined_uncertainty(
    x,
    y,
    ux,
    uy,
    max_iter=100,
    tol=1.0e-12,
):
    """
    Fit y = C*x using uncertainties on BOTH x and y.

    The effective variance is
        uy^2 + (C*ux)^2,
    so the weight is updated iteratively.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ux = np.asarray(ux, float)
    uy = np.asarray(uy, float)

    good = (
        np.isfinite(x)
        & np.isfinite(y)
        & np.isfinite(ux)
        & np.isfinite(uy)
        & (x > 0)
        & (y > 0)
        & (ux >= 0)
        & (uy > 0)
    )
    x = x[good]
    y = y[good]
    ux = ux[good]
    uy = uy[good]

    if len(x) < 2:
        return {
            "weighted_C": np.nan,
            "weighted_C_err": np.nan,
            "weighted_chi2_per_dof": np.nan,
            "weighted_N": len(x),
        }
    #endif

    # Median ratio is a stable initial value.
    C = float(np.median(y / x))

    for _ in range(max_iter):
        var = uy**2 + (C * ux)**2
        w = 1.0 / var
        denom = float(np.sum(w * x * x))
        if denom <= 0:
            break
        #endif
        Cnew = float(np.sum(w * x * y) / denom)
        if abs(Cnew - C) < tol:
            C = Cnew
            break
        #endif
        C = Cnew
    #endfor

    var = uy**2 + (C * ux)**2
    w = 1.0 / var
    denom = float(np.sum(w * x * x))
    Cerr = math.sqrt(1.0 / denom) if denom > 0 else np.nan
    chi2 = float(np.sum((y - C * x) ** 2 / var))
    ndof = max(len(x) - 1, 1)

    return {
        "weighted_C": C,
        "weighted_C_err": Cerr,
        "weighted_chi2_per_dof": chi2 / ndof,
        "weighted_N": len(x),
    }


def summarize_pointwise_corrections(
    work,
    label,
    pure_threshold,
    outdir,
):
    """
    Report the direct pointwise correction distribution

        C_i = sigma_reference / sigma_Hayward_raw

    separately for FD-pure, FT-pure, and mixed topology populations.

    The median and central 68% interval are intentionally unweighted: they
    answer the simple diagnostic question "what scale does a typical point
    prefer?" without allowing a small number of tiny-error points to dominate.
    """
    if not (0.0 <= pure_threshold < 0.5):
        raise ValueError("--pure-threshold must satisfy 0 <= threshold < 0.5")
    #endif

    d = work.copy()

    d["_C_point"] = d["_target"] / d["_xs_raw"]
    d["_eps_data_over_mc_point"] = np.where(
        d["_C_point"] > 0,
        1.0 / d["_C_point"],
        np.nan,
    )

    d["_topology_class"] = "mixed"
    d.loc[d["_f_FT"] <= pure_threshold, "_topology_class"] = "FD-pure"
    d.loc[d["_f_FT"] >= 1.0 - pure_threshold, "_topology_class"] = "FT-pure"

    rows = []

    for cls in ("FD-pure", "FT-pure", "mixed", "all"):
        if cls == "all":
            q = d.copy()
        else:
            q = d[d["_topology_class"] == cls].copy()
        #endif

        good = (
            np.isfinite(q["_C_point"])
            & (q["_C_point"] > 0)
            & np.isfinite(q["_xs_raw"])
            & (q["_xs_raw"] > 0)
            & np.isfinite(q["_target"])
            & (q["_target"] > 0)
        )
        q = q[good].copy()

        if len(q) == 0:
            continue
        #endif

        c = q["_C_point"].to_numpy(float)
        eps = 1.0 / c

        c16, c50, c84 = np.percentile(c, [16, 50, 84])
        e16, e50, e84 = np.percentile(eps, [16, 50, 84])

        row = {
            "reference": label,
            "topology_class": cls,
            "pure_threshold": pure_threshold,
            "N": len(q),
            "C_median": float(c50),
            "C_mean": float(np.mean(c)),
            "C_p16": float(c16),
            "C_p84": float(c84),
            "eps_data_over_mc_median": float(e50),
            "eps_data_over_mc_mean": float(np.mean(eps)),
            "eps_data_over_mc_p16": float(e16),
            "eps_data_over_mc_p84": float(e84),
            "median_f_FT": float(np.median(q["_f_FT"])),
        }

        # If both Hayward and reference uncertainties are present, also report
        # a one-parameter weighted scale using BOTH. This is particularly useful
        # for the Lee exact-same-bin comparison.
        if "_hayward_unc" in q.columns and "_target_unc" in q.columns:
            weighted = solve_one_scale_combined_uncertainty(
                q["_xs_raw"].to_numpy(float),
                q["_target"].to_numpy(float),
                q["_hayward_unc"].to_numpy(float),
                q["_target_unc"].to_numpy(float),
            )
            row.update(weighted)

            wc = weighted["weighted_C"]
            wce = weighted["weighted_C_err"]
            if np.isfinite(wc) and wc > 0:
                row["weighted_eps_data_over_mc"] = 1.0 / wc
                row["weighted_eps_data_over_mc_err"] = (
                    wce / (wc * wc) if np.isfinite(wce) else np.nan
                )
            else:
                row["weighted_eps_data_over_mc"] = np.nan
                row["weighted_eps_data_over_mc_err"] = np.nan
            #endif
        #endif

        rows.append(row)
    #endfor

    summary = pd.DataFrame(rows)
    summary.to_csv(
        outdir / f"{label}_pointwise_correction_summary.csv",
        index=False,
    )

    keep = [
        c for c in [
            "_source_row_int",
            "_xs_raw",
            "_hayward_unc",
            "_target",
            "_target_unc",
            "_f_FD",
            "_f_FT",
            "_topology_class",
            "_C_point",
            "_eps_data_over_mc_point",
            "xB",
            "Q2",
            "t_abs",
            "phi_deg",
            "g_theta",
            "p_theta",
            "angle_region_2d",
            "photon_angle_region",
            "proton_angle_region",
            "hayward_point_id",
            "lee_point_id",
        ]
        if c in d.columns
    ]
    d[keep].to_csv(
        outdir / f"{label}_pointwise_corrections.csv",
        index=False,
    )

    return summary


def pretty_print_pointwise(summary):
    print("")
    print("Direct pointwise correction distributions:")
    print("  C_i = sigma_reference / sigma_Hayward_raw")
    print("  (median and central 68% interval are unweighted)")
    print("")

    for _, r in summary.iterrows():
        cls = str(r["topology_class"])
        print(
            f"  {cls:7s} N={int(r['N']):4d}: "
            f"C median={r['C_median']:.4f} "
            f"[{r['C_p16']:.4f}, {r['C_p84']:.4f}]"
        )
        print(
            " " * 20
            + f"eps_data/eps_MC median={r['eps_data_over_mc_median']:.4f} "
            f"[{r['eps_data_over_mc_p16']:.4f}, "
            f"{r['eps_data_over_mc_p84']:.4f}]"
        )
        if "weighted_C" in r.index and np.isfinite(r.get("weighted_C", np.nan)):
            print(
                " " * 20
                + f"combined-unc weighted C={r['weighted_C']:.4f}"
                + (
                    f" +/- {r['weighted_C_err']:.4f}"
                    if np.isfinite(r.get("weighted_C_err", np.nan))
                    else ""
                )
                + f", chi2/dof={r['weighted_chi2_per_dof']:.3f}"
            )
        #endif
    #endfor

def bootstrap_factors(df, target_col, unc_col, use_weights, nrep, seed):
    if nrep <= 0:
        return {}
    #endif

    work = df[
        np.isfinite(df["_xs_raw"])
        & np.isfinite(df["_f_FD"])
        & np.isfinite(df["_f_FT"])
        & np.isfinite(df[target_col])
    ].copy()

    if use_weights:
        work = work[np.isfinite(work[unc_col]) & (work[unc_col] > 0)].copy()
    #endif

    if len(work) < 3:
        return {}
    #endif

    rng = np.random.default_rng(seed)
    vals = []

    for _ in range(nrep):
        take = rng.integers(0, len(work), size=len(work))
        sample = work.iloc[take].reset_index(drop=True)
        try:
            r, _, _, _ = solve_linear(sample, target_col, unc_col, use_weights)
            if np.isfinite(r["C_FD"]) and np.isfinite(r["C_FT"]):
                vals.append((r["C_FD"], r["C_FT"]))
            #endif
        except Exception:
            pass
        #endtry
    #endfor

    if not vals:
        return {}
    #endif

    a = np.asarray(vals, float)
    out = {"bootstrap_valid": len(a)}

    for j, name in enumerate(("C_FD", "C_FT")):
        q16, q50, q84 = np.percentile(a[:, j], [16, 50, 84])
        out[f"bootstrap_{name}_median"] = float(q50)
        out[f"bootstrap_{name}_p16"] = float(q16)
        out[f"bootstrap_{name}_p84"] = float(q84)
    #endfor

    # Convert each replica to epsilon_data / epsilon_MC.
    inv = np.where(a > 0, 1.0 / a, np.nan)
    for j, name in enumerate(("FD", "FT")):
        z = inv[:, j]
        z = z[np.isfinite(z)]
        if len(z):
            q16, q50, q84 = np.percentile(z, [16, 50, 84])
            out[f"bootstrap_eps_data_over_mc_{name}_median"] = float(q50)
            out[f"bootstrap_eps_data_over_mc_{name}_p16"] = float(q16)
            out[f"bootstrap_eps_data_over_mc_{name}_p84"] = float(q84)
        #endif
    #endfor

    return out


def make_km15_fit(pass2, comparison_dir, args):
    model_file = resolve_comparison_file(
        comparison_dir,
        "canonical_world_data_with_models.csv",
    )

    ref = pd.read_csv(model_file, low_memory=False)
    ref = ref[ref["dataset"].astype(str).eq("pass2")].copy()

    # compare_world_dvcs_cross_sections.py saves source_row as the row index of
    # the supplied pass-2 CSV after parsing.  Use it to attach model values.
    ref["_source_row_int"] = pd.to_numeric(ref["source_row"], errors="coerce")
    ref = ref[np.isfinite(ref["_source_row_int"])].copy()
    ref["_source_row_int"] = ref["_source_row_int"].astype(int)

    work = pass2.copy()
    work["_source_row_int"] = np.arange(len(work), dtype=int)
    keep_cols = [
        "_source_row_int",
        "km15_native",
        "bh_native",
        "bh_fraction_km15",
        "point_unc_abs",
        "phi_deg",
        "xB",
        "Q2",
        "t_abs",
        "g_theta",
        "p_theta",
    ]
    keep_cols = [c for c in keep_cols if c in ref.columns]
    work = work.merge(ref[keep_cols], on="_source_row_int", how="inner")

    work["_target"] = pd.to_numeric(work["km15_native"], errors="coerce")

    # KM15 is treated as the reference curve without an assigned model error.
    # point_unc_abs is the Hayward point uncertainty saved in the canonical
    # comparison table, so use it as the measurement uncertainty.
    work["_hayward_unc"] = pd.to_numeric(
        work["point_unc_abs"],
        errors="coerce",
    )
    work["_target_unc"] = work["_hayward_unc"]

    # The uncertainty to minimize should track the measured point uncertainty.
    # Because sigma_corrected = C_eff * sigma_raw, the exact uncertainty would
    # also scale with C_eff.  For this quick inverse diagnostic we keep the
    # published/raw point uncertainty fixed; the unweighted solution can be
    # requested as a robustness check.
    if args.bh_min > 0:
        work = work[pd.to_numeric(work["bh_fraction_km15"], errors="coerce") >= args.bh_min]
    #endif

    if args.phi_edge_max is not None:
        phi = np.mod(pd.to_numeric(work["phi_deg"], errors="coerce"), 360.0)
        phi_edge = np.minimum(phi, 360.0 - phi)
        work = work[phi_edge <= args.phi_edge_max]
    #endif

    label = "KM15"
    if args.bh_min > 0:
        label += f"_BHfrac_ge_{args.bh_min:g}"
    #endif
    if args.phi_edge_max is not None:
        label += f"_phiEdge_le_{args.phi_edge_max:g}"
    #endif

    return label, work



def make_lee_fit(pass2, comparison_dir):
    lee_file = resolve_comparison_file(
        comparison_dir,
        "lee_hayward_exact_same_bin_points.csv",
    )
    canonical_file = resolve_comparison_file(
        comparison_dir,
        "canonical_world_data_with_models.csv",
    )

    ref = pd.read_csv(lee_file, low_memory=False).copy()
    canonical = pd.read_csv(canonical_file, low_memory=False)
    canonical = canonical[canonical["dataset"].astype(str).eq("pass2")].copy()

    # The exact-same-bin table identifies the Hayward point by point_id
    # (for example "pass2:5").  canonical_world_data_with_models.csv contains
    # both point_id and source_row, so use that explicit mapping rather than
    # trying to recover the row by floating-point kinematic matching.
    canonical["_source_row_int"] = pd.to_numeric(
        canonical["source_row"],
        errors="coerce",
    )
    canonical = canonical[np.isfinite(canonical["_source_row_int"])].copy()
    canonical["_source_row_int"] = canonical["_source_row_int"].astype(int)

    point_map = canonical[
        [
            "point_id",
            "_source_row_int",
            "xB",
            "Q2",
            "t_abs",
            "phi_deg",
            "point_unc_abs",
        ]
    ].drop_duplicates("point_id")

    work_ref = ref.merge(
        point_map,
        left_on="hayward_point_id",
        right_on="point_id",
        how="inner",
        validate="many_to_one",
    )

    work_ref["_hayward_unc"] = pd.to_numeric(
        work_ref["point_unc_abs"],
        errors="coerce",
    )

    if len(work_ref) == 0:
        raise RuntimeError(
            "Could not map Lee exact-same-bin points onto pass-2 source rows."
        )
    #endif

    # Transport the Lee measurement from the Lee mean kinematics to the
    # Hayward mean kinematics with the already-saved KM15 transport factor.
    work_ref["_target"] = (
        pd.to_numeric(work_ref["lee_xs"], errors="coerce")
        * pd.to_numeric(
            work_ref["km15_lee_to_hayward_mean_transport_factor"],
            errors="coerce",
        )
    )

    # Apply the same multiplicative transport to Lee's point uncertainty.
    work_ref["_target_unc"] = (
        pd.to_numeric(work_ref["lee_point_unc"], errors="coerce")
        * np.abs(
            pd.to_numeric(
                work_ref["km15_lee_to_hayward_mean_transport_factor"],
                errors="coerce",
            )
        )
    )

    keep_cols = [
        "_source_row_int",
        "_target",
        "_target_unc",
        "_hayward_unc",
        "lee_xs",
        "lee_point_unc",
        "km15_lee_to_hayward_mean_transport_factor",
        "angle_region_2d",
        "photon_angle_region",
        "proton_angle_region",
        "hayward_point_id",
        "lee_point_id",
    ]
    keep_cols = [c for c in keep_cols if c in work_ref.columns]

    mapdf = work_ref[keep_cols].drop_duplicates("_source_row_int")

    work = pass2.copy()
    work["_source_row_int"] = np.arange(len(work), dtype=int)
    work = work.merge(mapdf, on="_source_row_int", how="inner")

    print(
        f"[Lee mapping] exact same-bin rows mapped to original pass-2 CSV: "
        f"{len(work)}"
    )

    return "Lee2026_exact_same_bin_transported", work


def run_one_fit(label, work, args, outdir):
    use_weights = not args.no_point_weights

    result, idx, pred, resid = solve_linear(
        work,
        target_col="_target",
        unc_col="_target_unc",
        use_weights=use_weights,
    )

    boot = bootstrap_factors(
        work,
        target_col="_target",
        unc_col="_target_unc",
        use_weights=use_weights,
        nrep=args.bootstrap,
        seed=args.seed,
    )
    result.update(boot)
    result["reference"] = label
    result["weighted"] = use_weights
    result["fit_type"] = "original_global"

    # Save per-point diagnostic table for the original global fit.
    detail = work.loc[idx].copy()
    detail["_fit_target"] = detail["_target"]
    detail["_fit_prediction"] = pred
    detail["_fit_residual"] = resid
    detail["_fit_Ceff"] = (
        detail["_f_FD"] * result["C_FD"]
        + detail["_f_FT"] * result["C_FT"]
    )
    detail["_fit_eps_data_over_mc_eff"] = np.where(
        detail["_fit_Ceff"] > 0,
        1.0 / detail["_fit_Ceff"],
        np.nan,
    )
    detail.to_csv(outdir / f"{label}_points.csv", index=False)

    extra_results = []

    # For Lee, perform a second global fit with BOTH Hayward and Lee point
    # uncertainties. The original result is retained for backward comparison.
    if (
        label.startswith("Lee")
        and "_hayward_unc" in work.columns
        and not args.no_point_weights
    ):
        comb, cidx, cpred, cresid = solve_linear_combined_uncertainty(
            work,
            target_col="_target",
            target_unc_col="_target_unc",
            hayward_unc_col="_hayward_unc",
        )
        comb["reference"] = label
        comb["weighted"] = True
        comb["fit_type"] = "global_combined_hayward_plus_reference_unc"
        extra_results.append(comb)

        cdetail = work.loc[cidx].copy()
        cdetail["_fit_target"] = cdetail["_target"]
        cdetail["_fit_prediction"] = cpred
        cdetail["_fit_residual"] = cresid
        cdetail["_fit_Ceff"] = (
            cdetail["_f_FD"] * comb["C_FD"]
            + cdetail["_f_FT"] * comb["C_FT"]
        )
        cdetail.to_csv(
            outdir / f"{label}_combined_uncertainty_points.csv",
            index=False,
        )
    #endif

    pointwise = summarize_pointwise_corrections(
        work,
        label,
        args.pure_threshold,
        outdir,
    )

    return result, extra_results, pointwise


def pretty_print(r):
    print("\n" + "=" * 78)
    print(f"Reference: {r['reference']}")
    if "fit_type" in r:
        print(f"Fit type : {r['fit_type']}")
    #endif
    print(f"N points : {r['N']}")
    print("-" * 78)
    print("Cross-section multipliers C = epsilon_MC / epsilon_data")
    print(
        f"  FD: C_FD = {r['C_FD']:.5f}  "
        f"(formal {r['formal_err_C_FD']:.5f})"
    )
    print(
        f"  FT: C_FT = {r['C_FT']:.5f}  "
        f"(formal {r['formal_err_C_FT']:.5f})"
    )
    print("")
    print("Equivalent efficiency ratios epsilon_data / epsilon_MC")
    print(f"  FD: {r['eps_data_over_mc_FD']:.5f}")
    print(f"  FT: {r['eps_data_over_mc_FT']:.5f}")

    if "bootstrap_C_FD_median" in r:
        print("")
        print("Row-bootstrap central 68% intervals")
        print(
            f"  C_FD = {r['bootstrap_C_FD_median']:.5f} "
            f"[{r['bootstrap_C_FD_p16']:.5f}, {r['bootstrap_C_FD_p84']:.5f}]"
        )
        print(
            f"  C_FT = {r['bootstrap_C_FT_median']:.5f} "
            f"[{r['bootstrap_C_FT_p16']:.5f}, {r['bootstrap_C_FT_p84']:.5f}]"
        )
        print(
            "  eps_data/eps_MC FD = "
            f"{r['bootstrap_eps_data_over_mc_FD_median']:.5f} "
            f"[{r['bootstrap_eps_data_over_mc_FD_p16']:.5f}, "
            f"{r['bootstrap_eps_data_over_mc_FD_p84']:.5f}]"
        )
        print(
            "  eps_data/eps_MC FT = "
            f"{r['bootstrap_eps_data_over_mc_FT_median']:.5f} "
            f"[{r['bootstrap_eps_data_over_mc_FT_p16']:.5f}, "
            f"{r['bootstrap_eps_data_over_mc_FT_p84']:.5f}]"
        )
    #endif

    print("")
    print(f"FD/FT factor correlation = {r['corr_CFD_CFT']:.4f}")
    print(f"chi2/ndof = {r['chi2']:.2f}/{r['ndof']} = {r['chi2_per_dof']:.3f}")
    if "iterations" in r:
        print(
            f"combined-uncertainty iterations = {r['iterations']}, "
            f"converged = {r['converged']}"
        )
    #endif
    print("=" * 78)


def main():
    args = parse_args()

    pass2_file = Path(args.pass2_file).expanduser().resolve()
    comparison_dir = Path(args.comparison_dir).expanduser().resolve()
    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"[input] pass2       : {pass2_file}")
    print(f"[input] comparison  : {comparison_dir}")
    print(f"[output]             : {outdir}")

    pass2 = build_pass2_dataframe(pass2_file)

    results = []

    if args.reference in ("both", "km15"):
        label, work = make_km15_fit(pass2, comparison_dir, args)
        r, extras, pointwise = run_one_fit(label, work, args, outdir)
        results.append(r)
        results.extend(extras)
        pretty_print(r)
        for er in extras:
            pretty_print(er)
        #endfor
        pretty_print_pointwise(pointwise)
    #endif

    if args.reference in ("both", "lee"):
        label, work = make_lee_fit(pass2, comparison_dir)
        r, extras, pointwise = run_one_fit(label, work, args, outdir)
        results.append(r)
        results.extend(extras)
        pretty_print(r)
        for er in extras:
            pretty_print(er)
        #endfor
        pretty_print_pointwise(pointwise)
    #endif

    summary = pd.DataFrame(results)
    summary.to_csv(outdir / "backsolved_fd_ft_corrections.csv", index=False)

    print("\nCurrent photon-efficiency study values for comparison:")
    print("  FD high-E: epsilon_data/epsilon_MC = 0.70240")
    print("             -> C_FD = epsilon_MC/epsilon_data = 1.42369")
    print("  FT high-E: epsilon_data/epsilon_MC = 0.35061")
    print("             -> C_FT = epsilon_MC/epsilon_data = 2.85220")
    print("")
    print("IMPORTANT: these fitted factors are reference-dependent cross-section")
    print("diagnostics, not independent detector-efficiency measurements.")
    print(f"\n[wrote] {outdir / 'backsolved_fd_ft_corrections.csv'}")
    print("[wrote] per-reference *_pointwise_correction_summary.csv")
    print("[wrote] per-reference *_pointwise_corrections.csv")


if __name__ == "__main__":
    main()
