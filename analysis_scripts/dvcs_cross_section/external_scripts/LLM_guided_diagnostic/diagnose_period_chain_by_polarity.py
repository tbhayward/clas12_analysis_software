#!/usr/bin/env python3
"""
diagnose_period_chain_by_polarity.py

Purpose
-------
Diagnose where RGA DVCS period-to-period discrepancies enter the correction
chain without treating any single period as truth.

The script reads one pass-2 CSV, usually all_data.csv, and compares each period
to a same-polarity group reference:

  inbending:  Fa18 Inb, Sp18 Inb, Sp19 Inb
  outbending: Fa18 Out, Sp18 Out

For the final post-acceptance stages, the heatmap also prints a second
parenthesized line: the ratio to an all-five-period reference. This tests
whether the final corrected quantities are consistent across both torus
polarities/configurations.

Raw yields are compared only after luminosity normalization. Acceptance-corrected
and cross-section stages are compared directly through their physical values.

Outputs
-------
  period_chain_by_polarity_summary.csv
  period_chain_worst_bins.csv
  period_chain_by_polarity_heatmap.png
  period_chain_by_phi.png
  period_chain_by_xB.png
  period_chain_by_Q2.png
  period_chain_by_t_abs.png

Example
-------
  python diagnose_period_chain_by_polarity.py \
      --csv all_data.csv \
      --outdir output/period_chain_by_polarity
"""

import argparse
import csv
import math
import os
import re
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm


PERIODS = ["Fa18 Inb", "Sp18 Inb", "Sp19 Inb", "Fa18 Out", "Sp18 Out"]
GROUPS = {
    "inbending": ["Fa18 Inb", "Sp18 Inb", "Sp19 Inb"],
    "outbending": ["Fa18 Out", "Sp18 Out"],
}
TOPOLOGIES = ["(FD, FD)", "(CD, FD)", "(CD, FT)"]

STAGE_ORDER = [
    # DVCS candidate chain
    "raw_epg_over_lumi",
    "current_eff_epg",
    "normalized_raw_epg",

    # Exclusive pi0 control channel chain
    "raw_eppi0_over_lumi",
    "current_eff_eppi0",
    "normalized_raw_eppi0",

    # Background subtraction and final correction chain
    "contamination_ratio",
    "signal_over_lumi",
    "acceptance",
    "acceptance_corrected_over_lumi",
    "cross_section",
    "normed_cross_section",
]

STAGE_LABELS = {
    "raw_epg_over_lumi": "raw epγ / L",
    "current_eff_epg": "curr. eff. epγ",
    "normalized_raw_epg": "eff. corr. raw epγ / L",
    "raw_eppi0_over_lumi": "raw epπ0 / L",
    "current_eff_eppi0": "curr. eff. epπ0",
    "normalized_raw_eppi0": "eff. corr. raw epπ0 / L",
    "contamination_ratio": "π0 contam.",
    "signal_over_lumi": "signal / L",
    "acceptance": "acceptance",
    "acceptance_corrected_over_lumi": "acc.-corr. / L",
    "cross_section": "cross section",
    "normed_cross_section": "normed xsec",
}

ALL_PERIOD_REFERENCE_STAGES = {
    "acceptance_corrected_over_lumi",
    "cross_section",
    "normed_cross_section",
}


# ----------------------------- tuple/cell helpers -----------------------------

def parse_tuple3(cell):
    """Parse cells like '(value, stat, sys)' or a plain number."""
    if cell is None:
        return (np.nan, np.nan, np.nan)

    s = str(cell).strip()
    if not s:
        return (np.nan, np.nan, np.nan)

    # Strip repeated quote layers.
    changed = True
    while changed and len(s) >= 2:
        changed = False
        if (s[0] == s[-1]) and s[0] in ['"', "'"]:
            s = s[1:-1].strip()
            changed = True

    if len(s) >= 2 and s[0] == "(" and s[-1] == ")":
        s = s[1:-1].strip()

    if not s:
        return (np.nan, np.nan, np.nan)

    parts = [p.strip() for p in s.split(",")]
    vals = []
    for p in parts[:3]:
        try:
            vals.append(float(p))
        except ValueError:
            vals.append(np.nan)

    while len(vals) < 3:
        vals.append(np.nan)

    return tuple(vals[:3])


def get_triple(row, col):
    if col not in row.index:
        return (np.nan, np.nan, np.nan)
    return parse_tuple3(row[col])


def quad_sum(values):
    values = [v for v in values if np.isfinite(v)]
    if not values:
        return np.nan
    return float(math.sqrt(sum(v * v for v in values)))


def add_topology_triples(row, prefix, reaction, period):
    """Sum topology-specific yield columns and combine errors in quadrature."""
    vals, stats, syss = [], [], []

    for topo in TOPOLOGIES:
        col = f"{prefix}, {reaction}, {topo}, exp, {period}, unpol"
        v, st, sy = get_triple(row, col)
        if np.isfinite(v):
            vals.append(v)
        if np.isfinite(st):
            stats.append(st)
        if np.isfinite(sy):
            syss.append(sy)

    if not vals:
        return (np.nan, np.nan, np.nan)

    return (float(sum(vals)), quad_sum(stats), quad_sum(syss))


def div_triple_by_scalar(triple, denom):
    v, st, sy = triple
    if not np.isfinite(denom) or denom <= 0.0:
        return (np.nan, np.nan, np.nan)
    return (
        v / denom if np.isfinite(v) else np.nan,
        st / denom if np.isfinite(st) else np.nan,
        sy / denom if np.isfinite(sy) else np.nan,
    )


def lumi_for_period(row, period):
    col = f"integrated luminosity, {period} (nC)"
    v, _, _ = get_triple(row, col)
    return v


# ----------------------------- stage extraction ------------------------------

def stage_triple(row, period, stage):
    lumi = lumi_for_period(row, period)

    if stage == "raw_epg_over_lumi":
        y = add_topology_triples(row, "raw yield", "ep->epg", period)
        return div_triple_by_scalar(y, lumi)

    if stage == "raw_eppi0_over_lumi":
        y = add_topology_triples(row, "raw yield", "ep->eppi0", period)
        return div_triple_by_scalar(y, lumi)

    if stage == "current_eff_epg":
        return get_triple(row, f"current efficiency factor, ep->epg, exp, {period}")

    if stage == "normalized_raw_epg":
        y = add_topology_triples(row, "normalized raw yield", "ep->epg", period)
        return div_triple_by_scalar(y, lumi)

    if stage == "current_eff_eppi0":
        return get_triple(row, f"current efficiency factor, ep->eppi0, exp, {period}")

    if stage == "normalized_raw_eppi0":
        y = add_topology_triples(row, "normalized raw yield", "ep->eppi0", period)
        return div_triple_by_scalar(y, lumi)

    if stage == "contamination_ratio":
        return get_triple(row, f"contamination ratio, {period}")

    if stage == "signal_over_lumi":
        y = get_triple(row, f"signal yield, ep->epg, exp, {period}, unpol")
        return div_triple_by_scalar(y, lumi)

    if stage == "acceptance":
        return get_triple(row, f"acceptance, {period}")

    if stage == "acceptance_corrected_over_lumi":
        y = get_triple(row, f"acceptance corrected yield, ep->epg, exp, {period}, unpol")
        return div_triple_by_scalar(y, lumi)

    if stage == "cross_section":
        return get_triple(row, f"cross sections, ep->epg, exp, {period}, unpol")

    if stage == "normed_cross_section":
        return get_triple(row, f"normed cross sections, ep->epg, exp, {period}, unpol")

    raise ValueError(f"Unknown stage: {stage}")


def weighted_group_reference(values, errors, use_weighted=True):
    """Return weighted mean if possible; otherwise return unweighted median."""
    vals = np.asarray(values, dtype=float)
    errs = np.asarray(errors, dtype=float)

    good = np.isfinite(vals) & (vals > 0.0)
    if not np.any(good):
        return np.nan

    vals_good = vals[good]
    errs_good = errs[good]

    if use_weighted:
        wgood = np.isfinite(errs_good) & (errs_good > 0.0)
        if np.any(wgood):
            weights = 1.0 / (errs_good[wgood] ** 2)
            return float(np.sum(weights * vals_good[wgood]) / np.sum(weights))

    return float(np.nanmedian(vals_good))


def finite_positive(x):
    return np.isfinite(x) and x > 0.0


# ----------------------------- diagnostics table -----------------------------

def build_long_table(df, use_weighted=True):
    records = []

    for i, row in df.iterrows():
        bin_index = row.get("bin index", i)

        # Build an all-five-period reference for each stage. This is mainly used
        # for the final post-acceptance stages, but computing it generically keeps
        # the bookkeeping simple and transparent.
        all_refs = {}
        all_triples_by_stage = {}

        for stage in STAGE_ORDER:
            vals = []
            errs = []
            triples = {}

            for period in PERIODS:
                tr = stage_triple(row, period, stage)
                triples[period] = tr
                vals.append(tr[0])
                errs.append(tr[1])

            all_refs[stage] = weighted_group_reference(vals, errs, use_weighted=use_weighted)
            all_triples_by_stage[stage] = triples

        for group, periods in GROUPS.items():
            for stage in STAGE_ORDER:
                vals = []
                errs = []
                triples = {}

                for period in periods:
                    tr = all_triples_by_stage[stage][period]
                    triples[period] = tr
                    vals.append(tr[0])
                    errs.append(tr[1])

                ref = weighted_group_reference(vals, errs, use_weighted=use_weighted)
                all_ref = all_refs.get(stage, np.nan)

                if not finite_positive(ref):
                    continue

                for period in periods:
                    value, stat, sys = triples[period]
                    if not finite_positive(value):
                        continue

                    ratio_to_group = value / ref
                    ratio_to_all = value / all_ref if finite_positive(all_ref) else np.nan

                    records.append({
                        "row": i,
                        "bin index": bin_index,
                        "group": group,
                        "period": period,
                        "stage": stage,
                        "stage label": STAGE_LABELS[stage],
                        "value": value,
                        "stat": stat,
                        "sys": sys,
                        "group_reference": ref,
                        "all_period_reference": all_ref,
                        "ratio_to_group": ratio_to_group,
                        "ratio_to_all": ratio_to_all,
                        "abs_frac_dev": abs(ratio_to_group - 1.0),
                        "abs_frac_dev_all": abs(ratio_to_all - 1.0) if np.isfinite(ratio_to_all) else np.nan,
                        "xB index": row.get("xB index", np.nan),
                        "Q2 index": row.get("Q2 index", np.nan),
                        "t_abs index": row.get("t_abs index", np.nan),
                        "phi index": row.get("phi index", np.nan),
                        "xBavg": row.get("xBavg", np.nan),
                        "Q2avg": row.get("Q2avg", np.nan),
                        "t_abs_avg": row.get("t_abs_avg", np.nan),
                        "phiavg": row.get("phiavg", np.nan),
                    })

    return pd.DataFrame.from_records(records)

def robust_median_and_unc(values):
    r = np.asarray(values, dtype=float)
    r = r[np.isfinite(r) & (r > 0.0)]

    if len(r) == 0:
        return (np.nan, np.nan, np.nan, np.nan, 0)

    median = float(np.nanmedian(r))
    p16 = float(np.nanpercentile(r, 16))
    p84 = float(np.nanpercentile(r, 84))

    # Robust approximate statistical uncertainty on the median ratio.
    # For a Gaussian-like distribution, SE(median) ~= 1.253 * sigma / sqrt(N).
    # We estimate sigma robustly as half the central 68% interval.
    robust_sigma = 0.5 * (p84 - p16)
    median_stat_unc = float(1.253 * robust_sigma / math.sqrt(len(r)))

    return (median, median_stat_unc, p16, p84, int(len(r)))


def summarize(long_df):
    rows = []
    grouped = long_df.groupby(["group", "period", "stage", "stage label"], sort=False)

    for (group, period, stage, stage_label), g in grouped:
        r_group = g["ratio_to_group"].to_numpy(dtype=float)
        median, median_unc, p16, p84, n = robust_median_and_unc(r_group)
        if n == 0:
            continue

        r_all = g["ratio_to_all"].to_numpy(dtype=float) if "ratio_to_all" in g.columns else np.array([])
        all_median, all_unc, all_p16, all_p84, all_n = robust_median_and_unc(r_all)

        rows.append({
            "group": group,
            "period": period,
            "stage": stage,
            "stage label": stage_label,
            "valid bins": n,
            "median ratio": median,
            "median ratio stat unc": median_unc,
            "p16": p16,
            "p84": p84,
            "median abs frac dev": float(np.nanmedian(np.abs(r_group[np.isfinite(r_group)] - 1.0))),
            "frac outside 10pct": float(np.mean(np.abs(r_group[np.isfinite(r_group)] - 1.0) > 0.10)),
            "frac outside 20pct": float(np.mean(np.abs(r_group[np.isfinite(r_group)] - 1.0) > 0.20)),
            "all-period valid bins": all_n,
            "median ratio to all periods": all_median,
            "median ratio to all periods stat unc": all_unc,
            "all-period p16": all_p16,
            "all-period p84": all_p84,
        })

    return pd.DataFrame(rows)


# ---------------------------------- plots ------------------------------------

def save_heatmap(summary_df, outpath):
    pivot = summary_df.pivot_table(
        index=["group", "period"],
        columns="stage label",
        values="median ratio",
        aggfunc="first",
    )
    unc_pivot = summary_df.pivot_table(
        index=["group", "period"],
        columns="stage label",
        values="median ratio stat unc",
        aggfunc="first",
    )
    all_pivot = summary_df.pivot_table(
        index=["group", "period"],
        columns="stage label",
        values="median ratio to all periods",
        aggfunc="first",
    )
    all_unc_pivot = summary_df.pivot_table(
        index=["group", "period"],
        columns="stage label",
        values="median ratio to all periods stat unc",
        aggfunc="first",
    )

    cols = [STAGE_LABELS[s] for s in STAGE_ORDER if STAGE_LABELS[s] in pivot.columns]
    pivot = pivot[cols]
    unc_pivot = unc_pivot[cols]
    all_pivot = all_pivot[cols]
    all_unc_pivot = all_unc_pivot[cols]

    final_stage_labels = {STAGE_LABELS[s] for s in ALL_PERIOD_REFERENCE_STAGES}

    fig, ax = plt.subplots(figsize=(16.5, 6.5))
    mat = pivot.to_numpy(dtype=float)
    unc = unc_pivot.to_numpy(dtype=float)
    all_mat = all_pivot.to_numpy(dtype=float)
    all_unc = all_unc_pivot.to_numpy(dtype=float)

    # Use both same-polarity and all-period final-stage ratios when setting the
    # symmetric color range, so the parenthesized all-period diagnostic is not
    # visually outside the stated color scale.
    finite_values = [mat[np.isfinite(mat)]]
    if all_mat.size > 0:
        finite_values.append(all_mat[np.isfinite(all_mat)])
    finite = np.concatenate([v for v in finite_values if v.size > 0]) if any(v.size > 0 for v in finite_values) else np.array([])

    if finite.size > 0:
        max_abs_dev = max(0.50, float(np.nanmax(np.abs(finite - 1.0))))
    else:
        max_abs_dev = 0.50
    vmin = 1.0 - max_abs_dev
    vmax = 1.0 + max_abs_dev

    # Very light blue is good/near unity; red is bad on either side of unity.
    # The color scale is symmetric around one, so equally large low/high deviations
    # have equal visual severity.
    cmap = LinearSegmentedColormap.from_list(
        "red_lightblue_red",
        [(0.0, "#b2182b"), (0.5, "#dbeeff"), (1.0, "#b2182b")],
    )
    norm = TwoSlopeNorm(vmin=vmin, vcenter=1.0, vmax=vmax)

    im = ax.imshow(mat, aspect="auto", cmap=cmap, norm=norm)
    ax.set_xticks(np.arange(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=35, ha="right")
    ax.set_yticks(np.arange(len(pivot.index)))
    ax.set_yticklabels([f"{g}: {p}" for g, p in pivot.index])
    ax.set_title("Median ratio to same-polarity group reference")

    for iy in range(mat.shape[0]):
        for ix in range(mat.shape[1]):
            col_label = pivot.columns[ix]

            if np.isfinite(mat[iy, ix]):
                if np.isfinite(unc[iy, ix]):
                    label = f"{mat[iy, ix]:.2f}±{unc[iy, ix]:.2f}"
                else:
                    label = f"{mat[iy, ix]:.2f}"

                # For the final post-acceptance stages, add a second line in
                # parentheses showing the ratio to an all-five-period reference.
                # This tests whether the fully corrected quantities remain
                # consistent once both torus polarities/configurations are pooled.
                if col_label in final_stage_labels and np.isfinite(all_mat[iy, ix]):
                    if np.isfinite(all_unc[iy, ix]):
                        label += f"\n({all_mat[iy, ix]:.2f}±{all_unc[iy, ix]:.2f})"
                    else:
                        label += f"\n({all_mat[iy, ix]:.2f})"

                ax.text(ix, iy, label, ha="center", va="center", fontsize=6.5, color="black")

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("same-polarity ratio; light blue = close to 1, red = far from 1")

    ax.text(
        0.0,
        -0.22,
        "For final three columns, parenthesized second line = ratio to all-five-period reference.",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
    )

    fig.tight_layout()
    fig.savefig(outpath, dpi=180)
    plt.close(fig)

def save_projection(long_df, var, label, outpath):
    df = long_df.copy()
    df = df[np.isfinite(pd.to_numeric(df[var], errors="coerce"))]
    if df.empty:
        return

    stages_to_plot = [
        "raw_epg_over_lumi",
        "signal_over_lumi",
        "acceptance",
        "acceptance_corrected_over_lumi",
        "normed_cross_section",
    ]

    nrows = len(GROUPS)
    fig, axes = plt.subplots(nrows=nrows, ncols=1, figsize=(12, 7.0), sharex=True)
    if nrows == 1:
        axes = [axes]

    for ax, (group, periods) in zip(axes, GROUPS.items()):
        for period in periods:
            for stage in stages_to_plot:
                sub = df[(df["group"] == group) & (df["period"] == period) & (df["stage"] == stage)]
                if sub.empty:
                    continue
                proj = sub.groupby(var)["ratio_to_group"].median().reset_index()
                ax.plot(proj[var], proj["ratio_to_group"], marker="o", linewidth=1.0,
                        label=f"{period}: {STAGE_LABELS[stage]}")

        ax.axhline(1.0, linewidth=1.0)
        ax.set_ylabel("ratio")
        ax.set_title(group)
        ax.set_ylim(0.0, 2.0)
        ax.grid(True, alpha=0.3)

    axes[-1].set_xlabel(label)
    axes[0].legend(fontsize=7, ncol=2, loc="upper right")
    fig.suptitle(f"Correction-chain ratios vs {label}")
    fig.tight_layout()
    fig.savefig(outpath, dpi=180)
    plt.close(fig)


def save_worst_bins(long_df, outpath, n=200):
    cols = [
        "group", "period", "stage", "stage label", "ratio_to_group", "ratio_to_all",
        "abs_frac_dev", "abs_frac_dev_all", "bin index", "xB index", "Q2 index", "t_abs index", "phi index",
        "xBavg", "Q2avg", "t_abs_avg", "phiavg",
    ]
    worst = long_df.sort_values("abs_frac_dev", ascending=False).head(n)
    worst[cols].to_csv(outpath, index=False)


# ----------------------------------- main ------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Polarity-separated correction-chain diagnostics for RGA DVCS CSVs.")
    parser.add_argument("--csv", required=True, help="Input pass-2 CSV, usually all_data.csv")
    parser.add_argument("--outdir", default="output/period_chain_by_polarity", help="Output directory")
    parser.add_argument("--reference", choices=["weighted", "median"], default="weighted",
                        help="Same-polarity group reference convention")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    df = pd.read_csv(args.csv, dtype=str, keep_default_na=False)

    use_weighted = (args.reference == "weighted")
    long_df = build_long_table(df, use_weighted=use_weighted)
    summary_df = summarize(long_df)

    long_csv = os.path.join(args.outdir, "period_chain_by_polarity_long.csv")
    summary_csv = os.path.join(args.outdir, "period_chain_by_polarity_summary.csv")
    worst_csv = os.path.join(args.outdir, "period_chain_worst_bins.csv")

    long_df.to_csv(long_csv, index=False)
    summary_df.to_csv(summary_csv, index=False)
    save_worst_bins(long_df, worst_csv)

    save_heatmap(summary_df, os.path.join(args.outdir, "period_chain_by_polarity_heatmap.png"))

    projections = [
        ("phi index", "phi bin index", "period_chain_by_phi.png"),
        ("xB index", "xB bin index", "period_chain_by_xB.png"),
        ("Q2 index", "Q2 bin index", "period_chain_by_Q2.png"),
        ("t_abs index", "|t| bin index", "period_chain_by_t_abs.png"),
    ]

    for var, label, fname in projections:
        if var in long_df.columns:
            save_projection(long_df, var, label, os.path.join(args.outdir, fname))

    print(f"[diagnostics] wrote: {summary_csv}")
    print(f"[diagnostics] wrote: {worst_csv}")
    print(f"[diagnostics] wrote plots to: {args.outdir}")


if __name__ == "__main__":
    main()
