#!/usr/bin/env python3
"""
diagnose_dvcs_sector_acceptance.py

CSV-level diagnostic suite for CLAS12 DVCS sector extractions.

Purpose
-------
Diagnose why sector-restricted extractions have very different numbers of usable
bins before attributing the effect to detector inefficiency or applying any
production correction.

The script compares:
  * inclusive parent
  * photon FD sector selections (CD-FD parent and CD-FD-S1...S6)
  * electron FD sector selections (eS1...eS6; missing files are allowed)
  * proton FD sector selections (pFD-S1...S6)
  * proton CD sector selections (pCD-S1...S3)

For every period and sector it explicitly examines:
  N_gen, N_rec, acceptance, DATA signal yield, and extracted cross section.

It produces:
  1. failure-reason tables
  2. generated/reconstructed MC diagnostics
  3. acceptance and sector/parent acceptance plots
  4. DATA and MC sector-fraction comparisons
  5. sector/parent cross-section plots
  6. valid/invalid-bin maps in analysis kinematics
  7. focused Sp18 Out pathology plots
  8. CSV files listing every classified bin

No production corrections are applied.

Example
-------
python diagnose_dvcs_sector_acceptance.py \
    --input-dir output/csvs \
    --output-dir output/sector_acceptance_diagnostics

If the archive has been unpacked into ./csvs:
python diagnose_dvcs_sector_acceptance.py --input-dir csvs
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


PERIODS = ["Fa18 Inb", "Fa18 Out", "Sp19 Inb", "Sp18 Inb", "Sp18 Out"]

KINEMATICS = {
    "xB": {
        "avg": "xBavg, {period}",
        "min": "xBmin",
        "max": "xBmax",
        "label": r"$x_B$",
    },
    "Q2": {
        "avg": "Q2avg, {period}",
        "min": "Q2min",
        "max": "Q2max",
        "label": r"$Q^2$ (GeV$^2$)",
    },
    "t": {
        "avg": "t_abs_avg, {period}",
        "min": "t_abs_min",
        "max": "t_abs_max",
        "label": r"$|t|$ (GeV$^2$)",
    },
    "phi": {
        "avg": "phiavg, {period}",
        "min": "phimin",
        "max": "phimax",
        "label": r"$\phi$ (deg)",
    },
}

GEOMETRY = {
    "e_theta": (r"$\theta_e$ (deg)", "e_theta, {period}"),
    "p_theta": (r"$\theta_p$ (deg)", "p_theta, {period}"),
    "p_phi": (r"$\phi_p$ (deg)", "p_phi, {period}"),
    "g_theta": (r"$\theta_\gamma$ (deg)", "g_theta, {period}"),
}

GROUP_SPECS = {
    "photon": {
        "parent": "dvcs_pass2_CD-FD.csv",
        "members": {f"S{i}": f"dvcs_pass2_CD-FD-S{i}.csv" for i in range(1, 7)},
        "title": "Photon FD sector",
    },
    "electron": {
        "parent": "dvcs_pass2_INCLUSIVE.csv",
        "members": {f"S{i}": f"dvcs_pass2_eS{i}.csv" for i in range(1, 7)},
        "title": "Electron FD sector",
    },
    "proton_FD": {
        "parent": "dvcs_pass2_INCLUSIVE.csv",
        "members": {f"S{i}": f"dvcs_pass2_pFD-S{i}.csv" for i in range(1, 7)},
        "title": "Proton FD sector",
    },
    "proton_CD": {
        "parent": "dvcs_pass2_INCLUSIVE.csv",
        "members": {f"S{i}": f"dvcs_pass2_pCD-S{i}.csv" for i in range(1, 4)},
        "title": "Proton CD sector",
    },
}

FAILURE_ORDER = [
    "parent_xs_invalid",
    "generated_mc_invalid",
    "generated_mc_zero",
    "reconstructed_mc_invalid",
    "reconstructed_mc_zero",
    "acceptance_invalid",
    "acceptance_nonpositive",
    "signal_invalid",
    "signal_nonpositive",
    "sector_xs_invalid",
    "sector_xs_nonpositive",
    "usable",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Diagnose DVCS sector acceptance and missing-bin structure."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("output/csvs"),
        help="Directory containing dvcs_pass2_*.csv files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output/sector_acceptance_diagnostics"),
        help="Output directory.",
    )
    parser.add_argument(
        "--periods",
        nargs="+",
        default=PERIODS,
        choices=PERIODS,
        help="Periods to analyze.",
    )
    parser.add_argument(
        "--focus-period",
        default="Sp18 Out",
        choices=PERIODS,
        help="Period used for detailed pathology maps.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=180,
        help="PNG resolution.",
    )
    parser.add_argument(
        "--min-positive",
        type=float,
        default=0.0,
        help="Strict lower bound used when testing positive yields/acceptance/XS.",
    )
    return parser.parse_args()


def sanitize(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("_")


def read_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, low_memory=False)
    if "bin index" not in df.columns:
        raise RuntimeError(f"{path} has no 'bin index' column.")
    # endif
    df = df.copy()
    df["bin index"] = pd.to_numeric(df["bin index"], errors="coerce")
    return df


def find_column(df: pd.DataFrame, candidates: Iterable[str]) -> Optional[str]:
    for col in candidates:
        if col in df.columns:
            return col
        # endif
    # endfor
    return None


def columns_for_period(df: pd.DataFrame, period: str) -> Dict[str, Optional[str]]:
    return {
        "generated": find_column(
            df,
            [f"generated yield, ep->epg, mc, {period}"],
        ),
        "reconstructed": find_column(
            df,
            [f"reconstructed yield, ep->epg, mc, {period}"],
        ),
        "acceptance": find_column(
            df,
            [f"acceptance, {period}"],
        ),
        "signal": find_column(
            df,
            [f"signal yield, ep->epg, exp, {period}, unpol"],
        ),
        "xs": find_column(
            df,
            [f"cross sections, ep->epg, exp, {period}, unpol"],
        ),
    }


def require_period_columns(df: pd.DataFrame, period: str, context: str) -> Dict[str, str]:
    cols = columns_for_period(df, period)
    missing = [key for key, col in cols.items() if col is None]
    if missing:
        raise RuntimeError(
            f"{context}, {period}: missing exact expected columns: {', '.join(missing)}"
        )
    # endif
    return {key: str(col) for key, col in cols.items()}


def topology_columns(df: pd.DataFrame, period: str) -> Dict[str, Optional[str]]:
    return {
        "FD_FD": find_column(
            df, [f"reconstructed yield, ep->epg, (FD, FD), mc, {period}"]
        ),
        "CD_FD": find_column(
            df, [f"reconstructed yield, ep->epg, (CD, FD), mc, {period}"]
        ),
        "CD_FT": find_column(
            df, [f"reconstructed yield, ep->epg, (CD, FT), mc, {period}"]
        ),
    }


def numeric(series: pd.Series) -> pd.Series:
    """
    Convert a DVCS CSV quantity to its central numerical value.

    The production CSV commonly stores quantities as tuple-like strings,
    e.g. "(12.34,0.56)" or "(12.34,0.56,...)", rather than plain scalars.
    pd.to_numeric() alone therefore turns valid physics entries into NaN.
    This parser extracts the first tuple component and also accepts ordinary
    numeric/scalar strings.
    """
    if pd.api.types.is_numeric_dtype(series):
        return pd.to_numeric(series, errors="coerce")
    # endif

    s = series.astype("string").str.strip()

    # ROOT/CSV tuple quantities are serialized with an additional literal
    # quote character, e.g. '"(83221,288.48,0)"'.  Search for the first
    # numeric token anywhere in the field rather than requiring the number
    # to follow '(' immediately.
    first = s.str.extract(
        r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)",
        expand=False,
    )
    return pd.to_numeric(first, errors="coerce")


def finite_positive(series: pd.Series, threshold: float = 0.0) -> pd.Series:
    x = numeric(series)
    return np.isfinite(x) & (x > threshold)


def safe_ratio(num: pd.Series, den: pd.Series) -> pd.Series:
    n = numeric(num)
    d = numeric(den)
    out = pd.Series(np.nan, index=n.index, dtype=float)
    mask = np.isfinite(n) & np.isfinite(d) & (d != 0)
    out.loc[mask] = n.loc[mask] / d.loc[mask]
    return out


def merge_parent_sector(parent: pd.DataFrame, sector: pd.DataFrame) -> pd.DataFrame:
    common = [c for c in ["bin index", "Bin Name"] if c in parent.columns and c in sector.columns]
    if "bin index" not in common:
        raise RuntimeError("Cannot merge without bin index.")
    # endif

    keep_parent = ["bin index"]
    if "Bin Name" in parent.columns:
        keep_parent.append("Bin Name")
    # endif

    # Keep all sector columns, but parent physics columns receive _parent suffix.
    merged = sector.merge(
        parent,
        on="bin index",
        how="outer",
        suffixes=("_sector", "_parent"),
        indicator=True,
    )
    return merged


def get_merged_column(
    merged: pd.DataFrame,
    original: Optional[str],
    side: str,
) -> Optional[str]:
    if original is None:
        return None
    # endif
    suffixed = f"{original}_{side}"
    if suffixed in merged.columns:
        return suffixed
    # endif
    if original in merged.columns:
        return original
    # endif
    return None


def classify_bins(
    parent: pd.DataFrame,
    sector: pd.DataFrame,
    period: str,
    threshold: float,
) -> pd.DataFrame:
    pc = columns_for_period(parent, period)
    sc = columns_for_period(sector, period)

    merged = merge_parent_sector(parent, sector)

    p_xs_col = get_merged_column(merged, pc["xs"], "parent")
    s_gen_col = get_merged_column(merged, sc["generated"], "sector")
    s_rec_col = get_merged_column(merged, sc["reconstructed"], "sector")
    s_acc_col = get_merged_column(merged, sc["acceptance"], "sector")
    s_sig_col = get_merged_column(merged, sc["signal"], "sector")
    s_xs_col = get_merged_column(merged, sc["xs"], "sector")

    required = {
        "parent xs": p_xs_col,
        "sector generated MC": s_gen_col,
        "sector reconstructed MC": s_rec_col,
        "sector acceptance": s_acc_col,
        "sector signal": s_sig_col,
        "sector xs": s_xs_col,
    }
    missing = [name for name, col in required.items() if col is None]
    if missing:
        raise RuntimeError(
            f"{period}: required columns absent: {', '.join(missing)}"
        )
    # endif

    p_xs = numeric(merged[p_xs_col])
    gen = numeric(merged[s_gen_col])
    rec = numeric(merged[s_rec_col])
    acc = numeric(merged[s_acc_col])
    sig = numeric(merged[s_sig_col])
    xs = numeric(merged[s_xs_col])

    reason = pd.Series("usable", index=merged.index, dtype=object)

    # Classification is deliberately sequential. A bin is assigned to the
    # earliest failed stage so that counts are mutually exclusive.
    masks = [
        ("parent_xs_invalid", ~np.isfinite(p_xs)),
        ("generated_mc_invalid", np.isfinite(p_xs) & ~np.isfinite(gen)),
        ("generated_mc_zero", np.isfinite(p_xs) & np.isfinite(gen) & (gen <= threshold)),
        (
            "reconstructed_mc_invalid",
            np.isfinite(p_xs) & (gen > threshold) & ~np.isfinite(rec),
        ),
        (
            "reconstructed_mc_zero",
            np.isfinite(p_xs) & (gen > threshold) & np.isfinite(rec) & (rec <= threshold),
        ),
        (
            "acceptance_invalid",
            np.isfinite(p_xs) & (gen > threshold) & (rec > threshold) & ~np.isfinite(acc),
        ),
        (
            "acceptance_nonpositive",
            np.isfinite(p_xs) & (gen > threshold) & (rec > threshold)
            & np.isfinite(acc) & (acc <= threshold),
        ),
        (
            "signal_invalid",
            np.isfinite(p_xs) & (gen > threshold) & (rec > threshold)
            & (acc > threshold) & ~np.isfinite(sig),
        ),
        (
            "signal_nonpositive",
            np.isfinite(p_xs) & (gen > threshold) & (rec > threshold)
            & (acc > threshold) & np.isfinite(sig) & (sig <= threshold),
        ),
        (
            "sector_xs_invalid",
            np.isfinite(p_xs) & (gen > threshold) & (rec > threshold)
            & (acc > threshold) & (sig > threshold) & ~np.isfinite(xs),
        ),
        (
            "sector_xs_nonpositive",
            np.isfinite(p_xs) & (gen > threshold) & (rec > threshold)
            & (acc > threshold) & (sig > threshold)
            & np.isfinite(xs) & (xs <= threshold),
        ),
    ]

    assigned = pd.Series(False, index=merged.index)
    for label, mask in masks:
        use = mask & ~assigned
        reason.loc[use] = label
        assigned |= use
    # endfor

    out = pd.DataFrame(
        {
            "bin index": merged["bin index"],
            "merge_status": merged["_merge"].astype(str),
            "failure_reason": reason,
            "parent_xs": p_xs,
            "sector_generated_mc": gen,
            "sector_reconstructed_mc": rec,
            "sector_acceptance": acc,
            "sector_signal_data": sig,
            "sector_xs": xs,
        }
    )

    # Add kinematics/geometry, preferring sector averages and then parent.
    for key, spec in KINEMATICS.items():
        original = spec["avg"].format(period=period)
        scol = get_merged_column(merged, original, "sector")
        pcol = get_merged_column(merged, original, "parent")
        if scol is not None:
            out[key] = numeric(merged[scol])
        elif pcol is not None:
            out[key] = numeric(merged[pcol])
        else:
            out[key] = np.nan
        # endif
    # endfor

    for key, (_, template) in GEOMETRY.items():
        original = template.format(period=period)
        scol = get_merged_column(merged, original, "sector")
        pcol = get_merged_column(merged, original, "parent")
        if scol is not None:
            out[key] = numeric(merged[scol])
        elif pcol is not None:
            out[key] = numeric(merged[pcol])
        else:
            out[key] = np.nan
        # endif
    # endfor

    return out


def setup_axis(ax: plt.Axes, xlabel: str, ylabel: str, title: str) -> None:
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(alpha=0.25)


def binned_median(
    x: pd.Series,
    y: pd.Series,
    nbins: int = 12,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    x = numeric(x)
    y = numeric(y)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask].to_numpy()
    y = y[mask].to_numpy()

    if len(x) < 3:
        return np.array([]), np.array([]), np.array([]), np.array([])
    # endif

    xmin, xmax = np.nanmin(x), np.nanmax(x)
    if not np.isfinite(xmin) or not np.isfinite(xmax) or xmin == xmax:
        return np.array([]), np.array([]), np.array([]), np.array([])
    # endif

    edges = np.linspace(xmin, xmax, nbins + 1)
    centers, medians, lo, hi = [], [], [], []

    for left, right in zip(edges[:-1], edges[1:]):
        m = (x >= left) & (x < right)
        if right == edges[-1]:
            m = (x >= left) & (x <= right)
        # endif
        if np.count_nonzero(m) < 2:
            continue
        # endif
        vals = y[m]
        centers.append(0.5 * (left + right))
        medians.append(np.nanmedian(vals))
        lo.append(np.nanpercentile(vals, 16))
        hi.append(np.nanpercentile(vals, 84))
    # endfor

    return (
        np.asarray(centers),
        np.asarray(medians),
        np.asarray(lo),
        np.asarray(hi),
    )


def plot_acceptance_by_kinematic(
    datasets: Dict[str, pd.DataFrame],
    period: str,
    group_title: str,
    outdir: Path,
    dpi: int,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.ravel()

    for ax, (kin, spec) in zip(axes, KINEMATICS.items()):
        for label, df in datasets.items():
            c = columns_for_period(df, period)
            if c["acceptance"] is None:
                continue
            # endif
            xcol = spec["avg"].format(period=period)
            if xcol not in df.columns:
                continue
            # endif
            x = numeric(df[xcol])
            y = numeric(df[c["acceptance"]])
            xc, med, lo, hi = binned_median(x, y)
            if len(xc):
                ax.plot(xc, med, marker="o", linewidth=1.2, markersize=3, label=label)
            # endif
        # endfor
        setup_axis(ax, spec["label"], "Acceptance", f"{group_title}: {period}")
    # endfor

    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=min(6, len(labels)))
    # endif
    fig.suptitle(f"{group_title} acceptance vs kinematics — {period}", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(outdir / f"acceptance_vs_kinematics_{sanitize(period)}.png", dpi=dpi)
    plt.close(fig)


def plot_ratio_by_kinematic(
    parent: pd.DataFrame,
    datasets: Dict[str, pd.DataFrame],
    period: str,
    quantity: str,
    ylabel: str,
    group_title: str,
    outdir: Path,
    dpi: int,
) -> None:
    pc = columns_for_period(parent, period)
    pcol = pc[quantity]
    if pcol is None:
        return
    # endif

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.ravel()

    for ax, (kin, spec) in zip(axes, KINEMATICS.items()):
        for label, sector in datasets.items():
            sc = columns_for_period(sector, period)
            scol = sc[quantity]
            if scol is None:
                continue
            # endif

            kincol = spec["avg"].format(period=period)
            if kincol not in sector.columns:
                continue
            # endif

            m = sector[["bin index", kincol, scol]].merge(
                parent[["bin index", pcol]],
                on="bin index",
                how="inner",
                suffixes=("_sector", "_parent"),
            )
            ratio = safe_ratio(m[f"{scol}_sector"], m[f"{pcol}_parent"])
            xc, med, lo, hi = binned_median(m[kincol], ratio)
            if len(xc):
                ax.plot(xc, med, marker="o", linewidth=1.2, markersize=3, label=label)
            # endif
        # endfor
        ax.axhline(1.0, linewidth=1.0, linestyle="--")
        setup_axis(ax, spec["label"], ylabel, f"{group_title}: {period}")
    # endfor

    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=min(6, len(labels)))
    # endif
    fig.suptitle(f"{group_title} {ylabel} vs kinematics — {period}", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(
        outdir / f"{sanitize(quantity)}_sector_over_parent_{sanitize(period)}.png",
        dpi=dpi,
    )
    plt.close(fig)


def plot_mc_acceptance_closure(
    parent: pd.DataFrame,
    datasets: Dict[str, pd.DataFrame],
    period: str,
    group_title: str,
    outdir: Path,
    dpi: int,
) -> None:
    pc = columns_for_period(parent, period)
    if pc["acceptance"] is None or pc["generated"] is None or pc["reconstructed"] is None:
        return
    # endif

    rows = []
    for label, sector in datasets.items():
        sc = columns_for_period(sector, period)
        needed = [sc["acceptance"], sc["generated"], sc["reconstructed"]]
        if any(x is None for x in needed):
            continue
        # endif

        cols_s = ["bin index", sc["acceptance"], sc["generated"], sc["reconstructed"]]
        cols_p = ["bin index", pc["acceptance"], pc["generated"], pc["reconstructed"]]
        m = sector[cols_s].merge(
            parent[cols_p],
            on="bin index",
            how="inner",
            suffixes=("_sector", "_parent"),
        )

        acc_ratio = safe_ratio(
            m[f"{sc['acceptance']}_sector"],
            m[f"{pc['acceptance']}_parent"],
        )
        rec_ratio = safe_ratio(
            m[f"{sc['reconstructed']}_sector"],
            m[f"{pc['reconstructed']}_parent"],
        )
        gen_ratio = safe_ratio(
            m[f"{sc['generated']}_sector"],
            m[f"{pc['generated']}_parent"],
        )
        expected_acc_ratio = safe_ratio(rec_ratio, gen_ratio)

        good = (
            np.isfinite(acc_ratio)
            & np.isfinite(rec_ratio)
            & np.isfinite(gen_ratio)
            & np.isfinite(expected_acc_ratio)
        )
        for a, r, g, e in zip(
            acc_ratio[good],
            rec_ratio[good],
            gen_ratio[good],
            expected_acc_ratio[good],
        ):
            rows.append(
                {
                    "sector": label,
                    "acceptance_ratio": a,
                    "reconstructed_ratio": r,
                    "generated_ratio": g,
                    "rec_over_gen_ratio": e,
                }
            )
        # endfor
    # endfor

    if not rows:
        return
    # endif

    d = pd.DataFrame(rows)

    fig, ax = plt.subplots(figsize=(8, 7))
    for label, sub in d.groupby("sector"):
        ax.scatter(
            sub["rec_over_gen_ratio"],
            sub["acceptance_ratio"],
            s=10,
            alpha=0.35,
            label=label,
        )
    # endfor

    finite = np.isfinite(d["rec_over_gen_ratio"]) & np.isfinite(d["acceptance_ratio"])
    vals = np.concatenate(
        [
            d.loc[finite, "rec_over_gen_ratio"].to_numpy(),
            d.loc[finite, "acceptance_ratio"].to_numpy(),
        ]
    )
    if len(vals):
        lo = max(0.0, np.nanpercentile(vals, 1))
        hi = np.nanpercentile(vals, 99)
        if np.isfinite(lo) and np.isfinite(hi) and hi > lo:
            ax.plot([lo, hi], [lo, hi], linestyle="--", linewidth=1.0)
            ax.set_xlim(lo, hi)
            ax.set_ylim(lo, hi)
        # endif
    # endif

    setup_axis(
        ax,
        r"$(N_{\rm rec,s}/N_{\rm rec,p})/(N_{\rm gen,s}/N_{\rm gen,p})$",
        r"$A_s/A_p$",
        f"Acceptance bookkeeping closure — {period}",
    )
    ax.legend(ncol=2)
    fig.tight_layout()
    fig.savefig(outdir / f"acceptance_mc_closure_{sanitize(period)}.png", dpi=dpi)
    plt.close(fig)

    summary = (
        d.groupby("sector")
        .agg(
            n=("acceptance_ratio", "size"),
            median_acceptance_ratio=("acceptance_ratio", "median"),
            median_reconstructed_ratio=("reconstructed_ratio", "median"),
            median_generated_ratio=("generated_ratio", "median"),
            median_rec_over_gen_ratio=("rec_over_gen_ratio", "median"),
        )
        .reset_index()
    )
    summary.to_csv(outdir / f"acceptance_mc_closure_{sanitize(period)}.csv", index=False)


def plot_data_mc_sector_fraction(
    parent: pd.DataFrame,
    datasets: Dict[str, pd.DataFrame],
    period: str,
    group_title: str,
    outdir: Path,
    dpi: int,
) -> None:
    pc = columns_for_period(parent, period)
    if pc["signal"] is None or pc["reconstructed"] is None:
        return
    # endif

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    for label, sector in datasets.items():
        sc = columns_for_period(sector, period)
        if sc["signal"] is None or sc["reconstructed"] is None:
            continue
        # endif

        m = sector[
            ["bin index", sc["signal"], sc["reconstructed"], f"t_abs_avg, {period}"]
        ].merge(
            parent[["bin index", pc["signal"], pc["reconstructed"]]],
            on="bin index",
            how="inner",
            suffixes=("_sector", "_parent"),
        )

        data_frac = safe_ratio(
            m[f"{sc['signal']}_sector"],
            m[f"{pc['signal']}_parent"],
        )
        mc_frac = safe_ratio(
            m[f"{sc['reconstructed']}_sector"],
            m[f"{pc['reconstructed']}_parent"],
        )
        double_ratio = safe_ratio(data_frac, mc_frac)

        x = numeric(m[f"t_abs_avg, {period}"])
        xc, med, _, _ = binned_median(x, data_frac)
        if len(xc):
            axes[0].plot(xc, med, marker="o", markersize=3, label=f"{label} DATA")
        # endif
        xc, med, _, _ = binned_median(x, mc_frac)
        if len(xc):
            axes[0].plot(xc, med, marker="x", markersize=3, linestyle="--", label=f"{label} MC")
        # endif

        xc, med, _, _ = binned_median(x, double_ratio)
        if len(xc):
            axes[1].plot(xc, med, marker="o", markersize=3, label=label)
        # endif
    # endfor

    setup_axis(
        axes[0],
        r"$|t|$ (GeV$^2$)",
        "Sector / parent fraction",
        f"{group_title}: DATA and reconstructed MC",
    )
    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        axes[0].legend(ncol=2, fontsize=8)
    # endif

    axes[1].axhline(1.0, linestyle="--", linewidth=1.0)
    setup_axis(
        axes[1],
        r"$|t|$ (GeV$^2$)",
        "(DATA sector fraction) / (MC sector fraction)",
        f"{group_title}: DATA/MC population ratio",
    )
    handles, labels = axes[1].get_legend_handles_labels()
    if handles:
        axes[1].legend(ncol=2, fontsize=8)
    # endif

    fig.suptitle(f"{period}", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(outdir / f"data_mc_sector_fraction_vs_t_{sanitize(period)}.png", dpi=dpi)
    plt.close(fig)


def plot_failure_counts(
    classifications: Dict[str, pd.DataFrame],
    period: str,
    group_title: str,
    outdir: Path,
    dpi: int,
) -> pd.DataFrame:
    rows = []
    for sector, df in classifications.items():
        counts = df["failure_reason"].value_counts()
        for reason in FAILURE_ORDER:
            rows.append(
                {
                    "period": period,
                    "sector": sector,
                    "failure_reason": reason,
                    "count": int(counts.get(reason, 0)),
                }
            )
        # endfor
    # endfor

    summary = pd.DataFrame(rows)
    pivot = summary.pivot(index="sector", columns="failure_reason", values="count").fillna(0)

    # Plot the physically most useful mutually-exclusive categories.
    plot_cols = [
        c for c in FAILURE_ORDER
        if c in pivot.columns and c != "parent_xs_invalid"
    ]

    fig, ax = plt.subplots(figsize=(12, 6))
    bottom = np.zeros(len(pivot))
    x = np.arange(len(pivot))

    for reason in plot_cols:
        vals = pivot[reason].to_numpy()
        ax.bar(x, vals, bottom=bottom, label=reason)
        bottom += vals
    # endfor

    ax.set_xticks(x)
    ax.set_xticklabels(pivot.index)
    ax.set_ylabel("Number of parent-valid analysis bins")
    ax.set_xlabel("Sector")
    ax.set_title(f"{group_title}: sequential bin classification — {period}")
    ax.legend(fontsize=8, ncol=2)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(outdir / f"failure_counts_{sanitize(period)}.png", dpi=dpi)
    plt.close(fig)

    return summary


def plot_failure_maps(
    classified: pd.DataFrame,
    sector: str,
    period: str,
    group_title: str,
    outdir: Path,
    dpi: int,
) -> None:
    # Show where usable vs non-usable bins lie in three complementary planes.
    pairs = [
        ("xB", "t", r"$x_B$", r"$|t|$ (GeV$^2$)"),
        ("Q2", "t", r"$Q^2$ (GeV$^2$)", r"$|t|$ (GeV$^2$)"),
        ("phi", "t", r"$\phi$ (deg)", r"$|t|$ (GeV$^2$)"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.2))

    usable = classified["failure_reason"] == "usable"
    parent_valid = classified["failure_reason"] != "parent_xs_invalid"

    for ax, (xkey, ykey, xlabel, ylabel) in zip(axes, pairs):
        bad = parent_valid & ~usable
        good = parent_valid & usable

        ax.scatter(
            classified.loc[bad, xkey],
            classified.loc[bad, ykey],
            s=13,
            alpha=0.45,
            marker="x",
            label="Not usable",
        )
        ax.scatter(
            classified.loc[good, xkey],
            classified.loc[good, ykey],
            s=10,
            alpha=0.35,
            marker="o",
            label="Usable",
        )
        setup_axis(ax, xlabel, ylabel, f"{sector}")
    # endfor

    axes[0].legend()
    fig.suptitle(f"{group_title} {sector}: usable-bin geography — {period}", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(
        outdir / f"failure_map_{sanitize(sector)}_{sanitize(period)}.png",
        dpi=dpi,
    )
    plt.close(fig)


def plot_failure_reason_kinematics(
    classified: pd.DataFrame,
    sector: str,
    period: str,
    group_title: str,
    outdir: Path,
    dpi: int,
) -> None:
    reasons = [
        r for r in FAILURE_ORDER
        if r != "parent_xs_invalid" and (classified["failure_reason"] == r).any()
    ]
    if not reasons:
        return
    # endif

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.ravel()

    for ax, (kin, spec) in zip(axes, KINEMATICS.items()):
        for reason in reasons:
            vals = numeric(classified.loc[classified["failure_reason"] == reason, kin])
            vals = vals[np.isfinite(vals)]
            if len(vals) < 2:
                continue
            # endif
            ax.hist(vals, bins=20, histtype="step", linewidth=1.3, label=reason)
        # endfor
        setup_axis(ax, spec["label"], "Bins", f"{sector}")
    # endfor

    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=3, fontsize=8)
    # endif
    fig.suptitle(f"{group_title} {sector}: failure reasons in kinematics — {period}", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(
        outdir / f"failure_reason_kinematics_{sanitize(sector)}_{sanitize(period)}.png",
        dpi=dpi,
    )
    plt.close(fig)


def plot_generated_reconstructed_acceptance(
    datasets: Dict[str, pd.DataFrame],
    period: str,
    group_title: str,
    outdir: Path,
    dpi: int,
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.3))

    for label, df in datasets.items():
        c = columns_for_period(df, period)
        tcol = f"t_abs_avg, {period}"
        if tcol not in df.columns:
            continue
        # endif

        for ax, key, ylabel in [
            (axes[0], "generated", r"Generated MC yield"),
            (axes[1], "reconstructed", r"Reconstructed MC yield"),
            (axes[2], "acceptance", r"Acceptance"),
        ]:
            col = c[key]
            if col is None:
                continue
            # endif
            xc, med, _, _ = binned_median(df[tcol], df[col])
            if len(xc):
                ax.plot(xc, med, marker="o", markersize=3, label=label)
            # endif
        # endfor
    # endfor

    for ax, ylabel in zip(
        axes,
        ["Generated MC yield", "Reconstructed MC yield", "Acceptance"],
    ):
        setup_axis(ax, r"$|t|$ (GeV$^2$)", ylabel, group_title)
    # endfor

    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        axes[0].legend(ncol=2, fontsize=8)
    # endif
    fig.suptitle(f"MC population and acceptance — {period}", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(
        outdir / f"generated_reconstructed_acceptance_vs_t_{sanitize(period)}.png",
        dpi=dpi,
    )
    plt.close(fig)


def plot_focus_geometry(
    classifications: Dict[str, pd.DataFrame],
    period: str,
    group_title: str,
    outdir: Path,
    dpi: int,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.ravel()

    for ax, (key, (xlabel, _)) in zip(axes, GEOMETRY.items()):
        for sector, df in classifications.items():
            bad = (
                (df["failure_reason"] != "usable")
                & (df["failure_reason"] != "parent_xs_invalid")
            )
            vals = numeric(df.loc[bad, key])
            vals = vals[np.isfinite(vals)]
            if len(vals) < 3:
                continue
            # endif
            ax.hist(vals, bins=20, histtype="step", density=True, linewidth=1.2, label=sector)
        # endfor
        setup_axis(ax, xlabel, "Normalized bin density", "Non-usable parent-valid bins")
    # endfor

    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=min(6, len(labels)))
    # endif
    fig.suptitle(f"{group_title}: geometry of failed bins — {period}", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(outdir / f"failed_bin_geometry_{sanitize(period)}.png", dpi=dpi)
    plt.close(fig)



def write_numerical_audit(
    parent: pd.DataFrame,
    datasets: Dict[str, pd.DataFrame],
    period: str,
    group_name: str,
    outdir: Path,
) -> pd.DataFrame:
    """
    Write one row per sector with explicit counts at every stage of the chain.
    Counts are evaluated on bins with a finite, positive parent cross section.
    """
    pc = require_period_columns(parent, period, f"{group_name} parent")
    rows = []

    for sector_name, sector in datasets.items():
        sc = require_period_columns(sector, period, f"{group_name} {sector_name}")

        keep_p = ["bin index", pc["xs"]]
        keep_s = [
            "bin index",
            sc["generated"],
            sc["reconstructed"],
            sc["acceptance"],
            sc["signal"],
            sc["xs"],
        ]
        # Rename before merging.  Pandas only applies merge suffixes to
        # overlapping column names.  Here parent contributes only its XS
        # column, so generated/reconstructed/acceptance/signal from the sector
        # would otherwise retain their unsuffixed names.
        sector_map = {
            sc["generated"]: "sector_generated",
            sc["reconstructed"]: "sector_reconstructed",
            sc["acceptance"]: "sector_acceptance",
            sc["signal"]: "sector_signal",
            sc["xs"]: "sector_xs",
        }
        parent_map = {pc["xs"]: "parent_xs"}

        m = sector[keep_s].rename(columns=sector_map).merge(
            parent[keep_p].rename(columns=parent_map),
            on="bin index",
            how="outer",
        )

        p_xs = numeric(m["parent_xs"])
        gen = numeric(m["sector_generated"])
        rec = numeric(m["sector_reconstructed"])
        acc = numeric(m["sector_acceptance"])
        sig = numeric(m["sector_signal"])
        xs = numeric(m["sector_xs"])

        base = np.isfinite(p_xs) & (p_xs > 0)
        gen_ok = base & np.isfinite(gen) & (gen > 0)
        rec_ok = gen_ok & np.isfinite(rec) & (rec > 0)
        acc_ok = rec_ok & np.isfinite(acc) & (acc > 0)
        sig_ok = acc_ok & np.isfinite(sig) & (sig > 0)
        xs_ok = sig_ok & np.isfinite(xs) & (xs > 0)

        rows.append(
            {
                "group": group_name,
                "period": period,
                "sector": sector_name,
                "parent_positive_xs_bins": int(base.sum()),
                "sector_positive_generated_mc": int(gen_ok.sum()),
                "lost_at_generated_mc": int(base.sum() - gen_ok.sum()),
                "sector_positive_reconstructed_mc": int(rec_ok.sum()),
                "lost_at_reconstructed_mc": int(gen_ok.sum() - rec_ok.sum()),
                "sector_positive_acceptance": int(acc_ok.sum()),
                "lost_at_acceptance": int(rec_ok.sum() - acc_ok.sum()),
                "sector_positive_data_signal": int(sig_ok.sum()),
                "lost_at_data_signal": int(acc_ok.sum() - sig_ok.sum()),
                "sector_positive_xs": int(xs_ok.sum()),
                "lost_at_cross_section": int(sig_ok.sum() - xs_ok.sum()),
            }
        )
    # endfor

    audit = pd.DataFrame(rows)
    audit.to_csv(outdir / f"numerical_chain_audit_{sanitize(period)}.csv", index=False)
    return audit


def plot_chain_survival(
    audit: pd.DataFrame,
    period: str,
    group_title: str,
    outdir: Path,
    dpi: int,
) -> None:
    stages = [
        ("parent_positive_xs_bins", "Parent XS"),
        ("sector_positive_generated_mc", "Ngen"),
        ("sector_positive_reconstructed_mc", "Nrec"),
        ("sector_positive_acceptance", "Acceptance"),
        ("sector_positive_data_signal", "DATA signal"),
        ("sector_positive_xs", "Sector XS"),
    ]

    fig, ax = plt.subplots(figsize=(10, 6))
    x = np.arange(len(stages))
    plotted = False

    for _, row in audit.iterrows():
        y = [row[col] for col, _ in stages]
        ax.plot(x, y, marker="o", linewidth=1.5, label=row["sector"])
        plotted = True
    # endfor

    ax.set_xticks(x)
    ax.set_xticklabels([label for _, label in stages])
    ax.set_ylabel("Number of bins")
    ax.set_xlabel("Sequential analysis stage")
    ax.set_title(f"{group_title}: bin survival chain — {period}")
    ax.grid(alpha=0.25)
    if plotted:
        ax.legend(ncol=2)
    # endif
    fig.tight_layout()
    fig.savefig(outdir / f"bin_survival_chain_{sanitize(period)}.png", dpi=dpi)
    plt.close(fig)


def plot_exact_stage_ratios(
    parent: pd.DataFrame,
    datasets: Dict[str, pd.DataFrame],
    period: str,
    group_name: str,
    group_title: str,
    outdir: Path,
    dpi: int,
) -> None:
    """
    Plot exact sector/parent ratios for Ngen, Nrec, acceptance, signal, and XS
    versus |t|. No assumption is made that sector and parent Ngen are equal.
    """
    pc = require_period_columns(parent, period, f"{group_name} parent")
    stages = [
        ("generated", r"$N_{\rm gen,s}/N_{\rm gen,p}$"),
        ("reconstructed", r"$N_{\rm rec,s}/N_{\rm rec,p}$"),
        ("acceptance", r"$A_s/A_p$"),
        ("signal", r"$N_{\rm sig,s}^{DATA}/N_{\rm sig,p}^{DATA}$"),
        ("xs", r"$\sigma_s/\sigma_p$"),
    ]

    fig, axes = plt.subplots(3, 2, figsize=(14, 13))
    axes = axes.ravel()

    for iax, (key, ylabel) in enumerate(stages):
        ax = axes[iax]
        plotted = False

        for sector_name, sector in datasets.items():
            sc = require_period_columns(sector, period, f"{group_name} {sector_name}")
            tcol = f"t_abs_avg, {period}"
            if tcol not in sector.columns:
                raise RuntimeError(f"{group_name} {sector_name}, {period}: missing {tcol}")
            # endif

            m = sector[["bin index", tcol, sc[key]]].merge(
                parent[["bin index", pc[key]]],
                on="bin index",
                how="inner",
                suffixes=("_sector", "_parent"),
            )
            ratio = safe_ratio(
                m[f"{sc[key]}_sector"],
                m[f"{pc[key]}_parent"],
            )
            xc, med, _, _ = binned_median(numeric(m[tcol]), ratio)
            if len(xc):
                ax.plot(xc, med, marker="o", markersize=3, linewidth=1.2, label=sector_name)
                plotted = True
            # endif
        # endfor

        ax.axhline(1.0, linestyle="--", linewidth=1.0)
        setup_axis(ax, r"$|t|$ (GeV$^2$)", ylabel, "")
        if plotted:
            ax.legend(ncol=2, fontsize=8)
        else:
            raise RuntimeError(
                f"{group_name}, {period}: no plottable values for exact {key} ratio"
            )
        # endif
    # endfor

    axes[-1].axis("off")
    fig.suptitle(f"{group_title}: exact sector/parent analysis chain — {period}", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(outdir / f"exact_stage_ratios_vs_t_{sanitize(period)}.png", dpi=dpi)
    plt.close(fig)


def plot_absolute_mc_diagnostics(
    datasets: Dict[str, pd.DataFrame],
    period: str,
    group_name: str,
    group_title: str,
    outdir: Path,
    dpi: int,
) -> None:
    """
    Show the absolute Ngen, Nrec, and acceptance values versus |t|.
    Also show Nrec versus Ngen bin-by-bin to expose empty/low-stat MC bins.
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    axes = axes.ravel()
    plotted = [False, False, False, False]

    for sector_name, df in datasets.items():
        c = require_period_columns(df, period, f"{group_name} {sector_name}")
        tcol = f"t_abs_avg, {period}"
        if tcol not in df.columns:
            raise RuntimeError(f"{group_name} {sector_name}, {period}: missing {tcol}")
        # endif

        x = numeric(df[tcol])
        gen = numeric(df[c["generated"]])
        rec = numeric(df[c["reconstructed"]])
        acc = numeric(df[c["acceptance"]])

        for i, (vals, ylabel) in enumerate(
            [
                (gen, "Generated MC yield"),
                (rec, "Reconstructed MC yield"),
                (acc, "Acceptance"),
            ]
        ):
            xc, med, _, _ = binned_median(x, vals)
            if len(xc):
                axes[i].plot(xc, med, marker="o", markersize=3, linewidth=1.2, label=sector_name)
                plotted[i] = True
            # endif
        # endfor

        good = np.isfinite(gen) & np.isfinite(rec) & (gen >= 0) & (rec >= 0)
        if good.any():
            axes[3].scatter(gen[good], rec[good], s=8, alpha=0.25, label=sector_name)
            plotted[3] = True
        # endif
    # endfor

    setup_axis(axes[0], r"$|t|$ (GeV$^2$)", "Generated MC yield", "")
    setup_axis(axes[1], r"$|t|$ (GeV$^2$)", "Reconstructed MC yield", "")
    setup_axis(axes[2], r"$|t|$ (GeV$^2$)", "Acceptance", "")
    setup_axis(axes[3], "Generated MC yield", "Reconstructed MC yield", "")

    for i, ax in enumerate(axes):
        if not plotted[i]:
            raise RuntimeError(
                f"{group_name}, {period}: absolute MC diagnostic panel {i} has no values"
            )
        # endif
        ax.legend(ncol=2, fontsize=8)
    # endfor

    fig.suptitle(f"{group_title}: absolute MC population and acceptance — {period}", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(outdir / f"absolute_mc_acceptance_{sanitize(period)}.png", dpi=dpi)
    plt.close(fig)


def plot_topology_fractions(
    datasets: Dict[str, pd.DataFrame],
    period: str,
    group_name: str,
    group_title: str,
    outdir: Path,
    dpi: int,
) -> None:
    """
    Decompose reconstructed epg MC into (FD,FD), (CD,FD), and (CD,FT)
    topology fractions versus |t| for each sector-restricted extraction.
    """
    topology_labels = {
        "FD_FD": "(FD, FD)",
        "CD_FD": "(CD, FD)",
        "CD_FT": "(CD, FT)",
    }

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.5))
    plotted = [False, False, False]

    for sector_name, df in datasets.items():
        tc = topology_columns(df, period)
        if any(col is None for col in tc.values()):
            missing = [k for k, col in tc.items() if col is None]
            raise RuntimeError(
                f"{group_name} {sector_name}, {period}: missing topology columns {missing}"
            )
        # endif

        tcol = f"t_abs_avg, {period}"
        if tcol not in df.columns:
            raise RuntimeError(f"{group_name} {sector_name}, {period}: missing {tcol}")
        # endif

        vals = {key: numeric(df[col]) for key, col in tc.items()}
        total = vals["FD_FD"] + vals["CD_FD"] + vals["CD_FT"]

        for i, key in enumerate(["FD_FD", "CD_FD", "CD_FT"]):
            frac = safe_ratio(vals[key], total)
            xc, med, _, _ = binned_median(numeric(df[tcol]), frac)
            if len(xc):
                axes[i].plot(xc, med, marker="o", markersize=3, linewidth=1.2, label=sector_name)
                plotted[i] = True
            # endif
        # endfor
    # endfor

    for i, key in enumerate(["FD_FD", "CD_FD", "CD_FT"]):
        setup_axis(
            axes[i],
            r"$|t|$ (GeV$^2$)",
            "Fraction of reconstructed epγ MC",
            topology_labels[key],
        )
        axes[i].set_ylim(0.0, 1.0)
        if plotted[i]:
            axes[i].legend(ncol=2, fontsize=8)
        else:
            raise RuntimeError(
                f"{group_name}, {period}: no plottable topology fraction for {key}"
            )
        # endif
    # endfor

    fig.suptitle(f"{group_title}: reconstructed-MC topology composition — {period}", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(outdir / f"topology_fractions_vs_t_{sanitize(period)}.png", dpi=dpi)
    plt.close(fig)


def plot_failure_reason_maps_colored(
    classified: pd.DataFrame,
    sector: str,
    period: str,
    group_title: str,
    outdir: Path,
    dpi: int,
) -> None:
    """
    Four phase-space maps with points separated by the actual sequential
    failure stage. This is intentionally categorical rather than binary.
    """
    pairs = [
        ("xB", "Q2", r"$x_B$", r"$Q^2$ (GeV$^2$)"),
        ("xB", "t", r"$x_B$", r"$|t|$ (GeV$^2$)"),
        ("Q2", "t", r"$Q^2$ (GeV$^2$)", r"$|t|$ (GeV$^2$)"),
        ("phi", "t", r"$\phi$ (deg)", r"$|t|$ (GeV$^2$)"),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    axes = axes.ravel()

    reasons = [
        r for r in FAILURE_ORDER
        if r != "parent_xs_invalid" and (classified["failure_reason"] == r).any()
    ]

    for ax, (xkey, ykey, xlabel, ylabel) in zip(axes, pairs):
        for reason in reasons:
            sub = classified[classified["failure_reason"] == reason]
            ax.scatter(
                numeric(sub[xkey]),
                numeric(sub[ykey]),
                s=14,
                alpha=0.55,
                label=reason,
            )
        # endfor
        setup_axis(ax, xlabel, ylabel, "")
    # endfor

    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=3, fontsize=8)
    # endif
    fig.suptitle(
        f"{group_title} {sector}: exact failure-stage geography — {period}",
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(
        outdir / f"failure_stage_maps_{sanitize(sector)}_{sanitize(period)}.png",
        dpi=dpi,
    )
    plt.close(fig)

def write_readme(output_dir: Path) -> None:
    text = """DVCS sector acceptance diagnostics

Interpretation rules
====================
0. Production quantities stored as tuple-like strings are parsed using their
   first component as the central value. This is required for these CSVs;
   plain pd.to_numeric() would incorrectly convert them to NaN.
1. A bin is classified sequentially. The first failed stage wins.
2. "generated_mc_zero" means the CSV explicitly contains a finite generated MC
   yield <= the configured threshold. It does NOT mean the script inferred the
   cause from the acceptance.
3. "reconstructed_mc_zero" requires finite positive generated MC and an
   explicitly finite reconstructed MC yield <= threshold.
4. "acceptance_invalid/nonpositive" is only used after both generated and
   reconstructed MC passed.
5. "signal_invalid/nonpositive" is only used after MC and acceptance passed.
6. "sector_xs_invalid/nonpositive" is only used after all upstream quantities
   passed.
7. "usable" means every preceding quantity and the sector cross section are
   finite and positive.

The acceptance_mc_closure plots test the bookkeeping identity

    A_s/A_parent =
      (Nrec_s/Nrec_parent) / (Ngen_s/Ngen_parent)

bin by bin.  This is preferable to assuming that the generated denominator is
the same in sector and parent extractions.

The script does not apply any correction and does not modify production CSVs.
"""
    (output_dir / "README.txt").write_text(text)


def main() -> None:
    args = parse_args()
    input_dir = args.input_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    write_readme(output_dir)

    print(f"[input]  {input_dir}")
    print(f"[output] {output_dir}")

    all_failure_rows = []
    all_bin_rows = []
    inventory_rows = []

    for group_name, spec in GROUP_SPECS.items():
        parent_path = input_dir / spec["parent"]
        if not parent_path.exists():
            print(f"[skip] {group_name}: missing parent {parent_path.name}")
            continue
        # endif

        parent = read_csv(parent_path)
        members: Dict[str, pd.DataFrame] = {}

        for label, filename in spec["members"].items():
            path = input_dir / filename
            inventory_rows.append(
                {
                    "group": group_name,
                    "sector": label,
                    "filename": filename,
                    "exists": path.exists(),
                }
            )
            if not path.exists():
                print(f"[missing] {group_name} {label}: {filename}")
                continue
            # endif
            members[label] = read_csv(path)
        # endfor

        if not members:
            continue
        # endif

        group_dir = output_dir / group_name
        group_dir.mkdir(parents=True, exist_ok=True)

        for period in args.periods:
            period_dir = group_dir / sanitize(period)
            period_dir.mkdir(parents=True, exist_ok=True)

            print(f"[analyze] {group_name:10s} {period}")

            classifications: Dict[str, pd.DataFrame] = {}
            for label, sector in members.items():
                try:
                    cdf = classify_bins(
                        parent,
                        sector,
                        period,
                        threshold=args.min_positive,
                    )
                except RuntimeError as exc:
                    print(f"  [classification skipped] {label}: {exc}")
                    continue
                # endif

                cdf.insert(0, "group", group_name)
                cdf.insert(1, "sector", label)
                cdf.insert(2, "period", period)
                classifications[label] = cdf

                cdf.to_csv(
                    period_dir / f"bin_classification_{sanitize(label)}.csv",
                    index=False,
                )
                all_bin_rows.append(cdf)

                if period == args.focus_period:
                    plot_failure_maps(
                        cdf,
                        label,
                        period,
                        spec["title"],
                        period_dir,
                        args.dpi,
                    )
                    plot_failure_reason_kinematics(
                        cdf,
                        label,
                        period,
                        spec["title"],
                        period_dir,
                        args.dpi,
                    )
                # endif
            # endfor

            if classifications:
                audit = write_numerical_audit(
                    parent,
                    members,
                    period,
                    group_name,
                    period_dir,
                )
                plot_chain_survival(
                    audit,
                    period,
                    spec["title"],
                    period_dir,
                    args.dpi,
                )
                plot_exact_stage_ratios(
                    parent,
                    members,
                    period,
                    group_name,
                    spec["title"],
                    period_dir,
                    args.dpi,
                )
                plot_absolute_mc_diagnostics(
                    members,
                    period,
                    group_name,
                    spec["title"],
                    period_dir,
                    args.dpi,
                )
                plot_topology_fractions(
                    members,
                    period,
                    group_name,
                    spec["title"],
                    period_dir,
                    args.dpi,
                )

                if period == args.focus_period:
                    for label, cdf in classifications.items():
                        plot_failure_reason_maps_colored(
                            cdf,
                            label,
                            period,
                            spec["title"],
                            period_dir,
                            args.dpi,
                        )
                    # endfor
                # endif

                failure_summary = plot_failure_counts(
                    classifications,
                    period,
                    spec["title"],
                    period_dir,
                    args.dpi,
                )
                failure_summary.insert(0, "group", group_name)
                all_failure_rows.append(failure_summary)

                plot_acceptance_by_kinematic(
                    members,
                    period,
                    spec["title"],
                    period_dir,
                    args.dpi,
                )
                plot_ratio_by_kinematic(
                    parent,
                    members,
                    period,
                    "acceptance",
                    r"$A_{\rm sector}/A_{\rm parent}$",
                    spec["title"],
                    period_dir,
                    args.dpi,
                )
                plot_ratio_by_kinematic(
                    parent,
                    members,
                    period,
                    "xs",
                    r"$\sigma_{\rm sector}/\sigma_{\rm parent}$",
                    spec["title"],
                    period_dir,
                    args.dpi,
                )
                plot_generated_reconstructed_acceptance(
                    members,
                    period,
                    spec["title"],
                    period_dir,
                    args.dpi,
                )
                plot_mc_acceptance_closure(
                    parent,
                    members,
                    period,
                    spec["title"],
                    period_dir,
                    args.dpi,
                )
                plot_data_mc_sector_fraction(
                    parent,
                    members,
                    period,
                    spec["title"],
                    period_dir,
                    args.dpi,
                )

                if period == args.focus_period:
                    plot_focus_geometry(
                        classifications,
                        period,
                        spec["title"],
                        period_dir,
                        args.dpi,
                    )
                # endif
            # endif
        # endfor
    # endfor

    pd.DataFrame(inventory_rows).to_csv(
        output_dir / "file_inventory.csv",
        index=False,
    )

    if all_failure_rows:
        pd.concat(all_failure_rows, ignore_index=True).to_csv(
            output_dir / "failure_summary_all.csv",
            index=False,
        )
    # endif

    if all_bin_rows:
        pd.concat(all_bin_rows, ignore_index=True).to_csv(
            output_dir / "bin_classification_all.csv",
            index=False,
        )
    # endif

    print("")
    print("Done.")
    print(f"Main summary: {output_dir / 'failure_summary_all.csv'}")
    print(f"All bin classifications: {output_dir / 'bin_classification_all.csv'}")
    print(f"Plots are under: {output_dir}/<group>/<period>/")


if __name__ == "__main__":
    main()
