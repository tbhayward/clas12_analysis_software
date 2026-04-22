#!/usr/bin/env python3

import os
import re
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.lines import Line2D

INPUT_FILE = "/u/home/thayward/clas12_analysis_software/analysis_scripts/asymmetry_extraction/output/results/asymmetries_rgc_combined_inb_NH3_epi+_2_timeStamp_04_22_130920.txt"
OUTPUT_FILE = "output/enpi+_sector_dependence.pdf"

ASYMMETRY_CONFIGS = [
    {
        "key": "ALUsinphi",
        "ylabel": r"$F_{LU}^{\sin\phi}/F_{UU}$",
        "ylim_top": (-0.2, 0.2),
        "ylim_pull": (-5.0, 5.0),
    },
    {
        "key": "AULsinphi",
        "ylabel": r"$F_{UL}^{\sin\phi}/F_{UU}$",
        "ylim_top": (-0.2, 0.2),
        "ylim_pull": (-5.0, 5.0),
    },
    {
        "key": "AULsin2phi",
        "ylabel": r"$F_{UL}^{\sin2\phi}/F_{UU}$",
        "ylim_top": (-0.2, 0.2),
        "ylim_pull": (-5.0, 5.0),
    },
    {
        "key": "ALL",
        "ylabel": r"$F_{LL}/F_{UU}$",
        "ylim_top": (-0.2, 0.8),
        "ylim_pull": (-5.0, 5.0),
    },
    {
        "key": "ALLcosphi",
        "ylabel": r"$F_{LL}^{\cos\phi}/F_{UU}$",
        "ylim_top": (-0.2, 0.8),
        "ylim_pull": (-5.0, 5.0),
    },
]

SECTOR_STYLES = {
    1: {"label": "Sector 1", "marker": "o", "linestyle": "-", "mfc": "none"},
    2: {"label": "Sector 2", "marker": "s", "linestyle": "-", "mfc": "none"},
    3: {"label": "Sector 3", "marker": "^", "linestyle": "-", "mfc": "none"},
    4: {"label": "Sector 4", "marker": "D", "linestyle": "-", "mfc": "none"},
    5: {"label": "Sector 5", "marker": "v", "linestyle": "-", "mfc": "none"},
    6: {"label": "Sector 6", "marker": "P", "linestyle": "-", "mfc": "none"},
}


def ensure_output_dir(path):
    outdir = os.path.dirname(path)
    if outdir != "":
        os.makedirs(outdir, exist_ok=True)
    #endif
#endfor


def parse_triplets(triplet_block):
    triplet_pattern = re.compile(
        r"\{\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*,\s*"
        r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*,\s*"
        r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*\}"
    )

    matches = triplet_pattern.findall(triplet_block)

    triplets = []
    for match in matches:
        xb = float(match[0])
        val = float(match[1])
        err = float(match[2])
        triplets.append((xb, val, err))
    #endfor

    return triplets
#endfor


def parse_input_file(filename):
    with open(filename, "r") as f:
        text = f.read()
    #endif

    pattern = re.compile(
        r"enpi_sector(\d+)GEchi2Fits([A-Za-z0-9_]+)\s*=\s*\{(.*?)\};",
        re.DOTALL
    )

    data = {}
    for match in pattern.finditer(text):
        sector = int(match.group(1))
        asym_name = match.group(2)
        triplet_block = match.group(3)
        triplets = parse_triplets(triplet_block)

        if sector not in data:
            data[sector] = {}
        #endif

        data[sector][asym_name] = triplets
    #endfor

    return data
#endfor


def weighted_mean(values, errors):
    values = np.asarray(values, dtype=float)
    errors = np.asarray(errors, dtype=float)

    good = np.isfinite(values) & np.isfinite(errors) & (errors > 0.0)
    if not np.any(good):
        return np.nan, np.nan
    #endif

    weights = 1.0 / (errors[good] * errors[good])
    mean = np.sum(weights * values[good]) / np.sum(weights)
    mean_err = math.sqrt(1.0 / np.sum(weights))

    return mean, mean_err
#endfor


def build_asymmetry_arrays(parsed_data, asym_key):
    sectors = sorted(parsed_data.keys())
    sector_arrays = {}

    n_points = None

    for sector in sectors:
        if asym_key not in parsed_data[sector]:
            raise RuntimeError("Missing asymmetry {} for sector {}".format(asym_key, sector))
        #endif

        triplets = parsed_data[sector][asym_key]

        xvals = np.array([t[0] for t in triplets], dtype=float)
        yvals = np.array([t[1] for t in triplets], dtype=float)
        evals = np.array([t[2] for t in triplets], dtype=float)

        if n_points is None:
            n_points = len(xvals)
        else:
            if len(xvals) != n_points:
                raise RuntimeError(
                    "Inconsistent number of points across sectors for asymmetry {}: "
                    "sector {} has {} points but expected {}".format(asym_key, sector, len(xvals), n_points)
                )
            #endif
        #endif

        sector_arrays[sector] = {
            "x": xvals,
            "y": yvals,
            "err": evals,
        }
    #endfor

    x_plot = np.zeros(n_points, dtype=float)
    means = np.zeros(n_points, dtype=float)
    mean_errs = np.zeros(n_points, dtype=float)
    rms_per_bin = np.zeros(n_points, dtype=float)

    pulls = {}
    for sector in sectors:
        pulls[sector] = np.zeros(n_points, dtype=float)
    #endfor

    for i in range(n_points):
        xvals_i = []
        vals_i = []
        errs_i = []

        for sector in sectors:
            xvals_i.append(sector_arrays[sector]["x"][i])
            vals_i.append(sector_arrays[sector]["y"][i])
            errs_i.append(sector_arrays[sector]["err"][i])
        #endfor

        xvals_i = np.array(xvals_i, dtype=float)
        vals_i = np.array(vals_i, dtype=float)
        errs_i = np.array(errs_i, dtype=float)

        x_plot[i] = np.mean(xvals_i)

        mean_i, mean_err_i = weighted_mean(vals_i, errs_i)
        means[i] = mean_i
        mean_errs[i] = mean_err_i

        residuals = vals_i - mean_i
        rms_per_bin[i] = math.sqrt(np.mean(residuals * residuals))

        for isector, sector in enumerate(sectors):
            err = errs_i[isector]
            if err > 0.0:
                pulls[sector][i] = (vals_i[isector] - mean_i) / err
            else:
                pulls[sector][i] = np.nan
            #endif
        #endfor
    #endfor

    avg_rms = float(np.mean(rms_per_bin))

    return {
        "x_plot": x_plot,
        "sector_arrays": sector_arrays,
        "means": means,
        "mean_errs": mean_errs,
        "pulls": pulls,
        "rms_per_bin": rms_per_bin,
        "avg_rms": avg_rms,
    }
#endfor


def make_plot(all_results):
    fig = plt.figure(figsize=(18, 11))

    outer = GridSpec(
        2,
        3,
        figure=fig,
        left=0.06,
        right=0.98,
        bottom=0.07,
        top=0.95,
        wspace=0.28,
        hspace=0.24
    )

    legend_handles = []
    legend_labels = []

    for idx, cfg in enumerate(ASYMMETRY_CONFIGS):
        row = idx // 3
        col = idx % 3

        inner = GridSpecFromSubplotSpec(
            2,
            1,
            subplot_spec=outer[row, col],
            height_ratios=[3.0, 1.4],
            hspace=0.08
        )

        ax_top = fig.add_subplot(inner[0])
        ax_bot = fig.add_subplot(inner[1], sharex=ax_top)

        result = all_results[cfg["key"]]
        x_mean = result["x_plot"]

        for sector in sorted(result["sector_arrays"].keys()):
            style = SECTOR_STYLES[sector]
            x = result["sector_arrays"][sector]["x"]
            y = result["sector_arrays"][sector]["y"]
            yerr = result["sector_arrays"][sector]["err"]
            pull = result["pulls"][sector]

            eb = ax_top.errorbar(
                x,
                y,
                yerr=yerr,
                marker=style["marker"],
                linestyle=style["linestyle"],
                linewidth=1.2,
                markersize=5.5,
                capsize=2.5,
                markerfacecolor=style["mfc"],
                label=style["label"]
            )

            ax_bot.plot(
                x_mean,
                pull,
                marker=style["marker"],
                linestyle=style["linestyle"],
                linewidth=1.0,
                markersize=5.0,
                markerfacecolor=style["mfc"],
                color=eb[0].get_color()
            )

            if idx == 0:
                legend_handles.append(
                    Line2D(
                        [0],
                        [0],
                        marker=style["marker"],
                        linestyle=style["linestyle"],
                        linewidth=1.2,
                        markersize=6.0,
                        markerfacecolor=style["mfc"],
                        color=eb[0].get_color()
                    )
                )
                legend_labels.append(style["label"])
            #endif
        #endfor

        ax_top.errorbar(
            x_mean,
            result["means"],
            yerr=result["mean_errs"],
            marker="o",
            linestyle="--",
            linewidth=1.6,
            markersize=4.5,
            color="black",
            capsize=2.5,
            label="Weighted mean"
        )

        ax_top.set_xlim(0.0, 0.6)
        ax_top.set_ylim(cfg["ylim_top"])
        ax_top.set_ylabel(cfg["ylabel"], fontsize=13)
        ax_top.grid(True, alpha=0.3)
        ax_top.tick_params(axis="x", labelbottom=False)

        ax_bot.axhline(0.0, color="gray", linestyle="--", linewidth=1.0)
        ax_bot.set_xlim(0.0, 0.6)
        ax_bot.set_ylim(cfg["ylim_pull"])
        ax_bot.set_xlabel(r"$x_{B}$", fontsize=13)
        ax_bot.set_ylabel("Pull", fontsize=12)
        ax_bot.grid(True, alpha=0.3)

        ax_top.text(
            0.03,
            0.93,
            "Avg RMS = {:.4f}".format(result["avg_rms"]),
            transform=ax_top.transAxes,
            ha="left",
            va="top",
            fontsize=10,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.85, edgecolor="0.7")
        )
    #endfor

    legend_ax = fig.add_subplot(outer[1, 2])
    legend_ax.axis("off")

    legend_handles_with_mean = list(legend_handles)
    legend_labels_with_mean = list(legend_labels)
    legend_handles_with_mean.append(
        Line2D([0], [0], marker="o", linestyle="--", linewidth=1.6, markersize=5.0, color="black")
    )
    legend_labels_with_mean.append("Weighted mean")

    legend_ax.text(
        0.5,
        0.82,
        "Sector dependence",
        ha="center",
        va="center",
        fontsize=16
    )

    legend_ax.legend(
        legend_handles_with_mean,
        legend_labels_with_mean,
        loc="center",
        frameon=True,
        fontsize=13,
        ncol=1
    )

    fig.savefig(OUTPUT_FILE)
    plt.close(fig)
#endfor


def main():
    ensure_output_dir(OUTPUT_FILE)

    parsed_data = parse_input_file(INPUT_FILE)

    all_results = {}
    for cfg in ASYMMETRY_CONFIGS:
        all_results[cfg["key"]] = build_asymmetry_arrays(parsed_data, cfg["key"])
    #endfor

    make_plot(all_results)

    print("")
    print("Average RMS of sector values with respect to the weighted mean:")
    overall_rms_list = []

    for cfg in ASYMMETRY_CONFIGS:
        avg_rms = all_results[cfg["key"]]["avg_rms"]
        overall_rms_list.append(avg_rms)
        print("  {:<12s} : {:.6f}".format(cfg["key"], avg_rms))
    #endfor

    overall_avg = float(np.mean(np.array(overall_rms_list, dtype=float)))

    print("")
    print("Overall average across all five asymmetries: {:.6f}".format(overall_avg))
    print("")
    print("Saved plot to {}".format(OUTPUT_FILE))
#endfor


if __name__ == "__main__":
    main()
#endif