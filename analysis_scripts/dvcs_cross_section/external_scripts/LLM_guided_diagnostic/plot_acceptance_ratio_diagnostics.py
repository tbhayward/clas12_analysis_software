#!/usr/bin/env python3
"""
plot_acceptance_ratio_diagnostics_by_topology.py

Make 2x4 diagnostic canvases showing where the MC acceptance differs between
Fa18 and Sp18 periods.

If one CSV is supplied:
  - produce one plot for that CSV.

If four CSVs are supplied:
  - assume order:
      all_data.csv FD-FD.csv CD-FD.csv CD-FT.csv
  - produce one plot for each selection.

Main plot:
  - Fa18 Inb / Sp18 Inb acceptance
  - Fa18 Out / Sp18 Out acceptance

Each subplot shows the ratio of the mean acceptance versus one variable:
  x_B, Q^2, |t|, phi, theta_e, theta_p, theta_gamma

The eighth pad is intentionally left blank.

Outputs:
  - acceptance_ratio_2x4.png                       for one input CSV
  - acceptance_ratio_2x4_<selection>.png           for four input CSVs
  - acceptance_ratio_diagnostics.csv
  - acceptance_ratio_worst_bins.csv
  - acceptance_ratio_summary.txt
"""

import argparse
import math
import os
import re
from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


PERIOD_PAIRS = [
    ("Fa18 Inb", "Sp18 Inb", "Fa18 Inb / Sp18 Inb"),
    ("Fa18 Out", "Sp18 Out", "Fa18 Out / Sp18 Out"),
]

PERIODS = ["Fa18 Inb", "Sp18 Inb", "Fa18 Out", "Sp18 Out"]

FOUR_CSV_SELECTION_NAMES = ["all_data", "FD-FD", "CD-FD", "CD-FT"]

VARIABLES = [
    {
        "key": "xB",
        "title": r"$x_B$",
        "xlabel": r"$x_B$",
        "type": "native",
        "min_col": "xBmin",
        "max_col": "xBmax",
    },
    {
        "key": "Q2",
        "title": r"$Q^2$",
        "xlabel": r"$Q^2$ (GeV$^2$)",
        "type": "native",
        "min_col": "Q2min",
        "max_col": "Q2max",
    },
    {
        "key": "t",
        "title": r"$|t|$",
        "xlabel": r"$|t|$ (GeV$^2$)",
        "type": "native",
        "min_col": "t_abs_min",
        "max_col": "t_abs_max",
    },
    {
        "key": "phi",
        "title": r"$\phi$",
        "xlabel": r"$\phi$ (deg)",
        "type": "native",
        "min_col": "phimin",
        "max_col": "phimax",
    },
    {
        "key": "e_theta",
        "title": r"$\theta_e$",
        "xlabel": r"$\theta_e$ (deg)",
        "type": "period_avg",
        "prefix": "e_theta",
    },
    {
        "key": "p_theta",
        "title": r"$\theta_p$",
        "xlabel": r"$\theta_p$ (deg)",
        "type": "period_avg",
        "prefix": "p_theta",
    },
    {
        "key": "g_theta",
        "title": r"$\theta_\gamma$",
        "xlabel": r"$\theta_\gamma$ (deg)",
        "type": "period_avg",
        "prefix": "g_theta",
    },
]


@dataclass
class BinResult:
    selection: str
    variable: str
    bin_label: str
    x_center: float
    x_low: float
    x_high: float
    pair_label: str
    numerator_period: str
    denominator_period: str
    n_rows: int
    mean_acc_num: float
    mean_acc_den: float
    sem_acc_num: float
    sem_acc_den: float
    ratio: float
    ratio_err: float
    mean_gen_num: float
    mean_gen_den: float
    mean_reco_num: float
    mean_reco_den: float
    mean_reco_cc_num: float
    mean_reco_cc_den: float


def strip_cell(cell) -> str:
    if pd.isna(cell):
        return ""
    #endif

    s = str(cell).strip()

    changed = True
    while changed and len(s) >= 2:
        changed = False
        if (s[0] == '"' and s[-1] == '"') or (s[0] == "'" and s[-1] == "'"):
            s = s[1:-1].strip()
            changed = True
        #endif
    #endwhile

    return s


def parse_tuple_value(cell) -> Tuple[float, float, float]:
    s = strip_cell(cell)

    if not s:
        return (np.nan, np.nan, np.nan)
    #endif

    if s.startswith("(") and s.endswith(")"):
        s = s[1:-1].strip()
    #endif

    parts = [p.strip() for p in s.split(",")]

    try:
        value = float(parts[0]) if len(parts) >= 1 and parts[0] != "" else np.nan
        stat = float(parts[1]) if len(parts) >= 2 and parts[1] != "" else 0.0
        sys = float(parts[2]) if len(parts) >= 3 and parts[2] != "" else 0.0
    except ValueError:
        return (np.nan, np.nan, np.nan)
    #endtry

    return (value, stat, sys)


def parse_tuple_series(series: pd.Series) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    vals = []
    stats = []
    syss = []

    for cell in series:
        v, s, y = parse_tuple_value(cell)
        vals.append(v)
        stats.append(s)
        syss.append(y)
    #endfor

    return np.array(vals, dtype=float), np.array(stats, dtype=float), np.array(syss, dtype=float)


def numeric_col(df: pd.DataFrame, col: str) -> np.ndarray:
    if col not in df.columns:
        return np.full(len(df), np.nan)
    #endif

    return pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)


def value_col(df: pd.DataFrame, col: str) -> np.ndarray:
    if col not in df.columns:
        return np.full(len(df), np.nan)
    #endif

    vals, _, _ = parse_tuple_series(df[col])
    return vals


def first_existing_value_col(df: pd.DataFrame, candidates: List[str]) -> np.ndarray:
    for col in candidates:
        if col in df.columns:
            return value_col(df, col)
        #endif
    #endfor

    return np.full(len(df), np.nan)


def safe_mean(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]

    if len(x) == 0:
        return np.nan
    #endif

    return float(np.mean(x))


def safe_sem(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]

    if len(x) <= 1:
        return 0.0
    #endif

    return float(np.std(x, ddof=1) / math.sqrt(len(x)))


def ratio_and_error(num: float, den: float, num_err: float, den_err: float) -> Tuple[float, float]:
    if not np.isfinite(num) or not np.isfinite(den) or den <= 0.0:
        return (np.nan, np.nan)
    #endif

    ratio = num / den

    rel2 = 0.0

    if num > 0.0 and np.isfinite(num_err):
        rel2 += (num_err / num) ** 2
    #endif

    if den > 0.0 and np.isfinite(den_err):
        rel2 += (den_err / den) ** 2
    #endif

    return (ratio, abs(ratio) * math.sqrt(rel2))


def get_period_average_x(df: pd.DataFrame, prefix: str, p1: str, p2: str) -> np.ndarray:
    c1 = f"{prefix}, {p1}"
    c2 = f"{prefix}, {p2}"

    x1 = numeric_col(df, c1)
    x2 = numeric_col(df, c2)

    both = np.vstack([x1, x2])

    with np.errstate(invalid="ignore"):
        x = np.nanmean(both, axis=0)
    #endwith

    return x


def build_bins_for_variable(
    df: pd.DataFrame,
    variable: Dict,
    p1: str,
    p2: str,
    theta_bins: int,
    phi_combine: int,
) -> List[Tuple[str, float, float, np.ndarray]]:
    if variable["type"] == "native":
        lo = numeric_col(df, variable["min_col"])
        hi = numeric_col(df, variable["max_col"])

        valid = np.isfinite(lo) & np.isfinite(hi) & (hi > lo)
        keys = {}

        for i in np.where(valid)[0]:
            key = (round(float(lo[i]), 10), round(float(hi[i]), 10))

            if key not in keys:
                keys[key] = []
            #endif

            keys[key].append(i)
        #endfor

        sorted_keys = sorted(keys.keys(), key=lambda k: 0.5 * (k[0] + k[1]))

        if variable["key"] == "phi" and phi_combine > 1:
            grouped_keys = []

            for istart in range(0, len(sorted_keys), phi_combine):
                grouped_keys.append(sorted_keys[istart:istart + phi_combine])
            #endfor
        else:
            grouped_keys = [[key] for key in sorted_keys]
        #endif

        bins = []

        for group in grouped_keys:
            idx_list = []

            for key in group:
                idx_list.extend(keys[key])
            #endfor

            idx = np.array(idx_list, dtype=int)
            mask = np.zeros(len(df), dtype=bool)
            mask[idx] = True

            x_low = float(group[0][0])
            x_high = float(group[-1][1])
            label = f"{x_low:.4g}-{x_high:.4g}"

            bins.append((label, x_low, x_high, mask))
        #endfor

        return bins
    #endif

    x = get_period_average_x(df, variable["prefix"], p1, p2)
    valid_x = x[np.isfinite(x)]

    if len(valid_x) == 0:
        return []
    #endif

    q = np.linspace(0.0, 1.0, theta_bins + 1)
    edges = np.quantile(valid_x, q)
    edges = np.unique(edges)

    if len(edges) < 2:
        return []
    #endif

    bins = []

    for ibin in range(len(edges) - 1):
        lo = float(edges[ibin])
        hi = float(edges[ibin + 1])

        if ibin == len(edges) - 2:
            mask = np.isfinite(x) & (x >= lo) & (x <= hi)
        else:
            mask = np.isfinite(x) & (x >= lo) & (x < hi)
        #endif

        if np.sum(mask) == 0:
            continue
        #endif

        label = f"{lo:.3g}-{hi:.3g}"
        bins.append((label, lo, hi, mask))
    #endfor

    return bins


def build_bin_results(
    df: pd.DataFrame,
    selection: str,
    theta_bins: int,
    min_rows: int,
    phi_combine: int,
) -> List[BinResult]:
    results = []

    cache = {}

    for period in PERIODS:
        cache[("acceptance", period)] = value_col(df, f"acceptance, {period}")
        cache[("generated", period)] = value_col(df, f"generated yield, ep->epg, mc, {period}")
        cache[("reco", period)] = value_col(df, f"reconstructed yield, ep->epg, mc, {period}")

        cache[("reco_cc", period)] = first_existing_value_col(
            df,
            [
                f"reconstructed current corrected yield, ep->epg, mc, {period}",
                f"current corrected reconstructed yield, ep->epg, mc, {period}",
            ],
        )
    #endfor

    for variable in VARIABLES:
        for p_num, p_den, pair_label in PERIOD_PAIRS:
            bins = build_bins_for_variable(df, variable, p_num, p_den, theta_bins, phi_combine)

            acc_num_all = cache[("acceptance", p_num)]
            acc_den_all = cache[("acceptance", p_den)]

            gen_num_all = cache[("generated", p_num)]
            gen_den_all = cache[("generated", p_den)]

            reco_num_all = cache[("reco", p_num)]
            reco_den_all = cache[("reco", p_den)]

            reco_cc_num_all = cache[("reco_cc", p_num)]
            reco_cc_den_all = cache[("reco_cc", p_den)]

            for bin_label, x_low, x_high, base_mask in bins:
                valid = (
                    base_mask
                    & np.isfinite(acc_num_all)
                    & np.isfinite(acc_den_all)
                    & (acc_num_all > 0.0)
                    & (acc_den_all > 0.0)
                )

                n_rows = int(np.sum(valid))

                if n_rows < min_rows:
                    continue
                #endif

                mean_acc_num = safe_mean(acc_num_all[valid])
                mean_acc_den = safe_mean(acc_den_all[valid])
                sem_acc_num = safe_sem(acc_num_all[valid])
                sem_acc_den = safe_sem(acc_den_all[valid])

                ratio, ratio_err = ratio_and_error(
                    mean_acc_num,
                    mean_acc_den,
                    sem_acc_num,
                    sem_acc_den,
                )

                x_center = 0.5 * (x_low + x_high)

                results.append(
                    BinResult(
                        selection=selection,
                        variable=variable["key"],
                        bin_label=bin_label,
                        x_center=x_center,
                        x_low=x_low,
                        x_high=x_high,
                        pair_label=pair_label,
                        numerator_period=p_num,
                        denominator_period=p_den,
                        n_rows=n_rows,
                        mean_acc_num=mean_acc_num,
                        mean_acc_den=mean_acc_den,
                        sem_acc_num=sem_acc_num,
                        sem_acc_den=sem_acc_den,
                        ratio=ratio,
                        ratio_err=ratio_err,
                        mean_gen_num=safe_mean(gen_num_all[valid]),
                        mean_gen_den=safe_mean(gen_den_all[valid]),
                        mean_reco_num=safe_mean(reco_num_all[valid]),
                        mean_reco_den=safe_mean(reco_den_all[valid]),
                        mean_reco_cc_num=safe_mean(reco_cc_num_all[valid]),
                        mean_reco_cc_den=safe_mean(reco_cc_den_all[valid]),
                    )
                )
            #endfor
        #endfor
    #endfor

    return results


def results_to_dataframe(results: List[BinResult]) -> pd.DataFrame:
    rows = []

    for r in results:
        d = r.__dict__.copy()

        d["abs_log_ratio"] = abs(math.log(r.ratio)) if np.isfinite(r.ratio) and r.ratio > 0.0 else np.nan
        d["percent_deviation"] = 100.0 * (r.ratio - 1.0) if np.isfinite(r.ratio) else np.nan

        d["gen_ratio"] = (
            r.mean_gen_num / r.mean_gen_den
            if np.isfinite(r.mean_gen_den) and r.mean_gen_den > 0.0
            else np.nan
        )

        d["reco_ratio"] = (
            r.mean_reco_num / r.mean_reco_den
            if np.isfinite(r.mean_reco_den) and r.mean_reco_den > 0.0
            else np.nan
        )

        d["reco_cc_ratio"] = (
            r.mean_reco_cc_num / r.mean_reco_cc_den
            if np.isfinite(r.mean_reco_cc_den) and r.mean_reco_cc_den > 0.0
            else np.nan
        )

        rows.append(d)
    #endfor

    return pd.DataFrame(rows)


def make_plot(results_df: pd.DataFrame, out_png: str, selection: str) -> None:
    fig, axes = plt.subplots(2, 4, figsize=(20, 9), constrained_layout=True)
    axes_flat = axes.ravel()

    pair_styles = {
        "Fa18 Inb / Sp18 Inb": {
            "marker": "o",
            "linestyle": "-",
            "label": "Fa18 Inb / Sp18 Inb",
        },
        "Fa18 Out / Sp18 Out": {
            "marker": "s",
            "linestyle": "--",
            "label": "Fa18 Out / Sp18 Out",
        },
    }

    for i, variable in enumerate(VARIABLES):
        ax = axes_flat[i]
        sub = results_df[
            (results_df["selection"] == selection)
            & (results_df["variable"] == variable["key"])
        ].copy()

        for pair_label in [p[2] for p in PERIOD_PAIRS]:
            s = sub[sub["pair_label"] == pair_label].sort_values("x_center")

            if s.empty:
                continue
            #endif

            st = pair_styles[pair_label]

            ax.errorbar(
                s["x_center"].to_numpy(),
                s["ratio"].to_numpy(),
                yerr=s["ratio_err"].to_numpy(),
                marker=st["marker"],
                linestyle=st["linestyle"],
                linewidth=1.8,
                markersize=5.5,
                capsize=2.5,
                label=st["label"],
            )
        #endfor

        ax.axhline(1.0, linewidth=1.2, linestyle=":")
        ax.set_title(variable["title"], fontsize=14)
        ax.set_xlabel(variable["xlabel"], fontsize=12)
        ax.set_ylabel("mean acceptance ratio", fontsize=12)
        ax.grid(True, alpha=0.35)
        ax.tick_params(axis="both", labelsize=10)
        ax.set_ylim(0.6, 2.0)
        ax.legend(fontsize=9, frameon=True, loc="best")
    #endfor

    axes_flat[7].axis("off")
    axes_flat[7].text(
        0.04,
        0.84,
        "Acceptance diagnostic",
        fontsize=16,
        weight="bold",
        transform=axes_flat[7].transAxes,
    )
    axes_flat[7].text(
        0.04,
        0.70,
        "Plotted quantity:\n"
        r"$\langle A_{\mathrm{Fa18}}\rangle / \langle A_{\mathrm{Sp18}}\rangle$"
        "\n\n"
        "Each point averages the MC acceptance\n"
        "over all analysis bins in that x-axis bin.\n\n"
        "Use the output CSV to inspect generated\n"
        "and reconstructed MC yield ratios.",
        fontsize=12,
        transform=axes_flat[7].transAxes,
        va="top",
    )

    fig.suptitle(
        rf"Acceptance ratio diagnostics: {selection}, Fa18/Sp18 inbending and outbending",
        fontsize=18,
    )

    fig.savefig(out_png, dpi=200)
    plt.close(fig)


def print_console_summary(results_df: pd.DataFrame, worst_df: pd.DataFrame) -> None:
    print("\n[acceptance diagnostics] Median ratio summary by selection, variable, and pair")
    print("---------------------------------------------------------------------------")

    for selection in sorted(results_df["selection"].dropna().unique()):
        print(f"\nSelection: {selection}")

        for variable in [v["key"] for v in VARIABLES]:
            for pair_label in [p[2] for p in PERIOD_PAIRS]:
                s = results_df[
                    (results_df["selection"] == selection)
                    & (results_df["variable"] == variable)
                    & (results_df["pair_label"] == pair_label)
                ]

                if s.empty:
                    continue
                #endif

                med = float(np.nanmedian(s["ratio"]))
                p16 = float(np.nanpercentile(s["ratio"], 16))
                p84 = float(np.nanpercentile(s["ratio"], 84))

                print(
                    f"{variable:8s}  {pair_label:22s}  "
                    f"median={med:7.3f}  p16={p16:7.3f}  p84={p84:7.3f}"
                )
            #endfor
        #endfor
    #endfor

    print("\n[acceptance diagnostics] Largest absolute log-ratio bins")
    print("---------------------------------------------------------")

    for _, row in worst_df.head(12).iterrows():
        print(
            f"{row['selection']:9s}  {row['variable']:8s}  {row['pair_label']:22s}  "
            f"bin={str(row['bin_label']):>15s}  ratio={row['ratio']:.3f}  "
            f"gen={row['gen_ratio']:.3f}  reco={row['reco_ratio']:.3f}  "
            f"reco_cc={row['reco_cc_ratio']:.3f}  n={int(row['n_rows'])}"
        )
    #endfor


def write_summary_text(path: str, results_df: pd.DataFrame) -> None:
    with open(path, "w") as f:
        for selection in sorted(results_df["selection"].dropna().unique()):
            f.write(f"Selection: {selection}\n")
            f.write("=" * (11 + len(selection)) + "\n")

            for pair_label in [p[2] for p in PERIOD_PAIRS]:
                f.write(f"\n{pair_label}\n")
                f.write("-" * len(pair_label) + "\n")

                for variable in [v["key"] for v in VARIABLES]:
                    s = results_df[
                        (results_df["selection"] == selection)
                        & (results_df["variable"] == variable)
                        & (results_df["pair_label"] == pair_label)
                    ]

                    if s.empty:
                        continue
                    #endif

                    f.write(
                        f"{variable:10s} "
                        f"acceptance median={np.nanmedian(s['ratio']):.4g}, "
                        f"p16={np.nanpercentile(s['ratio'], 16):.4g}, "
                        f"p84={np.nanpercentile(s['ratio'], 84):.4g}, "
                        f"gen median={np.nanmedian(s['gen_ratio']):.4g}, "
                        f"reco median={np.nanmedian(s['reco_ratio']):.4g}, "
                        f"reco_cc median={np.nanmedian(s['reco_cc_ratio']):.4g}\n"
                    )
                #endfor
            #endfor

            f.write("\n\n")
        #endfor
    #endwith


def infer_selection_names(csv_paths: List[str]) -> List[str]:
    if len(csv_paths) == 1:
        return [os.path.splitext(os.path.basename(csv_paths[0]))[0]]
    #endif

    if len(csv_paths) == 4:
        return FOUR_CSV_SELECTION_NAMES
    #endif

    raise RuntimeError(
        "Expected either exactly 1 CSV or exactly 4 CSVs. "
        "For four CSVs, order must be: all_data.csv FD-FD.csv CD-FD.csv CD-FT.csv"
    )


def safe_selection_name(selection: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", selection)


def check_required_columns(df: pd.DataFrame) -> List[str]:
    missing = []

    for period in PERIODS:
        c = f"acceptance, {period}"
        if c not in df.columns:
            missing.append(c)
        #endif
    #endfor

    for variable in VARIABLES:
        if variable["type"] == "native":
            for c in [variable["min_col"], variable["max_col"]]:
                if c not in df.columns:
                    missing.append(c)
                #endif
            #endfor
        else:
            for period in PERIODS:
                c = f"{variable['prefix']}, {period}"
                if c not in df.columns:
                    missing.append(c)
                #endif
            #endfor
        #endif
    #endfor

    return missing


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Plot Fa18/Sp18 MC acceptance ratios versus kinematic variables. "
            "Accepts either one CSV or four CSVs ordered as all_data, FD-FD, CD-FD, CD-FT."
        )
    )
    parser.add_argument(
        "--csv",
        nargs="+",
        required=True,
        help=(
            "Input CSV(s). Either one CSV or four CSVs. "
            "If four: all_data.csv FD-FD.csv CD-FD.csv CD-FT.csv"
        ),
    )
    parser.add_argument(
        "--outdir",
        default="output/acceptance_ratio_diagnostics",
        help="Output directory.",
    )
    parser.add_argument(
        "--theta-bins",
        type=int,
        default=10,
        help="Number of quantile bins for theta variables.",
    )
    parser.add_argument(
        "--phi-combine",
        type=int,
        default=3,
        help="Number of adjacent native phi bins to combine. Default 3 gives 8 phi points for 24 native bins.",
    )
    parser.add_argument(
        "--min-rows",
        type=int,
        default=2,
        help="Minimum rows required per plotted bin.",
    )

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    try:
        selection_names = infer_selection_names(args.csv)
    except RuntimeError as e:
        print(f"[acceptance diagnostics] ERROR: {e}")
        raise SystemExit(1)
    #endtry

    all_results = []

    for selection, csv_path in zip(selection_names, args.csv):
        df = pd.read_csv(csv_path)

        missing = check_required_columns(df)

        if missing:
            print(f"[acceptance diagnostics] ERROR: {csv_path} is missing required columns:")
            for c in missing:
                print(f"  {c}")
            #endfor
            raise SystemExit(1)
        #endif

        results = build_bin_results(
            df=df,
            selection=selection,
            theta_bins=args.theta_bins,
            min_rows=args.min_rows,
            phi_combine=max(1, args.phi_combine),
        )

        results_df = results_to_dataframe(results)

        if results_df.empty:
            print(f"[acceptance diagnostics] WARNING: no valid results for {selection}")
            continue
        #endif

        all_results.append(results_df)

        if len(args.csv) == 1:
            out_png = os.path.join(args.outdir, "acceptance_ratio_2x4.png")
        else:
            out_png = os.path.join(
                args.outdir,
                f"acceptance_ratio_2x4_{safe_selection_name(selection)}.png",
            )
        #endif

        make_plot(results_df, out_png, selection)
        print(f"[acceptance diagnostics] Wrote plot: {out_png}")
    #endfor

    if not all_results:
        print("[acceptance diagnostics] ERROR: no valid acceptance-ratio bins were produced.")
        raise SystemExit(1)
    #endif

    combined_df = pd.concat(all_results, ignore_index=True)

    diag_csv = os.path.join(args.outdir, "acceptance_ratio_diagnostics.csv")
    worst_csv = os.path.join(args.outdir, "acceptance_ratio_worst_bins.csv")
    summary_txt = os.path.join(args.outdir, "acceptance_ratio_summary.txt")

    combined_df.to_csv(diag_csv, index=False)

    worst_df = combined_df.sort_values("abs_log_ratio", ascending=False)
    worst_df.to_csv(worst_csv, index=False)

    write_summary_text(summary_txt, combined_df)

    print_console_summary(combined_df, worst_df)

    print("\n[acceptance diagnostics] Wrote:")
    print(f"  {diag_csv}")
    print(f"  {worst_csv}")
    print(f"  {summary_txt}")


if __name__ == "__main__":
    main()