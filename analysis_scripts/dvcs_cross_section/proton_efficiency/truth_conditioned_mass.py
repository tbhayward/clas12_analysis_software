#!/usr/bin/env python3
"""Inspect reconstructed missing mass conditional on a truth-near proton tag.

This reads the histogram cache made by mass_shape_diagnostics.py. It does not
fit a signal, subtract a background, or calculate a proton efficiency.
"""

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def fraction(hist, mask):
    total = int(hist.sum())
    return float(hist[mask].sum() / total) if total else float("nan")


def histogram_median(hist, centers):
    total = int(hist.sum())
    return float(centers[np.searchsorted(np.cumsum(hist), (total + 1) // 2)]) if total else float("nan")


def read_cache(path):
    with np.load(path, allow_pickle=False) as source:
        keys = ("m2_edges", "theta_edges", "p_edges", "data_all", "mc_all",
                "data_matched", "mc_matched", "mc_truth_valid_near",
                "mc_truth_valid_near_aligned", "mc_truth_valid_low",
                "mc_truth_valid_high", "mc_truth_invalid_or_nonfinite")
        missing = set(keys) - set(source.files)
        if missing:
            raise ValueError(f"Cache lacks required histograms: {sorted(missing)}")
        cache = {key: source[key].copy() for key in keys}
    shape = (len(cache["m2_edges"]) - 1, len(cache["theta_edges"]) - 1,
             len(cache["p_edges"]) - 1)
    for key in keys[3:]:
        if cache[key].shape != shape or np.any(cache[key] < 0):
            raise ValueError(f"Invalid histogram {key}: expected shape {shape}")
    if np.any(cache["data_matched"] > cache["data_all"]):
        raise ValueError("Data matched counts exceed all tags")
    if np.any(cache["mc_matched"] > cache["mc_all"]):
        raise ValueError("MC matched counts exceed all tags")
    categories = sum(cache[f"mc_truth_{key}"] for key in
                     ("valid_low", "valid_near", "valid_high", "invalid_or_nonfinite"))
    if not np.array_equal(categories, cache["mc_all"]):
        raise ValueError("MC truth categories do not sum to all MC tags")
    if np.any(cache["mc_truth_valid_near_aligned"] > cache["mc_truth_valid_near"]):
        raise ValueError("Aligned counts exceed truth-near counts")
    return cache


def write_metrics(cache, output):
    edges = cache["m2_edges"]
    centers = (edges[:-1] + edges[1:]) / 2
    low = (centers >= .43) & (centers < .73)
    peak = (centers >= .73) & (centers < 1.03)
    columns = ("theta_low_deg", "theta_high_deg", "p_low_GeV", "p_high_GeV",
               "data_tags", "mc_tags", "mc_truth_near_tags", "mc_truth_near_aligned_tags",
               "data_matched_fraction", "mc_matched_fraction",
               "data_reco_m2_median", "mc_reco_m2_median", "mc_truth_near_reco_m2_median",
               "data_reco_peak_fraction", "mc_reco_peak_fraction",
               "mc_truth_near_reco_peak_fraction",
               "mc_truth_near_fraction_of_all", "mc_truth_near_fraction_of_reco_low_shoulder")
    with (output / "4_slice_metrics.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        for j in range(len(cache["theta_edges"]) - 1):
            for k in range(len(cache["p_edges"]) - 1):
                data = cache["data_all"][:, j, k]
                mc = cache["mc_all"][:, j, k]
                near = cache["mc_truth_valid_near"][:, j, k]
                aligned = cache["mc_truth_valid_near_aligned"][:, j, k]
                data_matched = cache["data_matched"][:, j, k]
                mc_matched = cache["mc_matched"][:, j, k]
                writer.writerow((cache["theta_edges"][j], cache["theta_edges"][j + 1],
                                 cache["p_edges"][k], cache["p_edges"][k + 1],
                                 int(data.sum()), int(mc.sum()), int(near.sum()),
                                 int(aligned.sum()), data_matched.sum() / data.sum()
                                 if data.sum() else float("nan"),
                                 mc_matched.sum() / mc.sum() if mc.sum() else float("nan"),
                                 histogram_median(data, centers), histogram_median(mc, centers),
                                 histogram_median(near, centers), fraction(data, peak),
                                 fraction(mc, peak), fraction(near, peak),
                                 near.sum() / mc.sum() if mc.sum() else float("nan"),
                                 near[low].sum() / mc[low].sum() if mc[low].sum() else float("nan")))


def plot_shapes(cache, output, theta_index):
    edges = cache["m2_edges"]
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex=True)
    series = (("data_all", "Data: all tags", "tab:blue"),
              ("mc_all", "MC: all tags", "tab:orange"),
              ("mc_truth_valid_near", "MC: truth-near proton", "tab:green"))
    for k, ax in enumerate(axes.flat):
        for key, label, color in series:
            hist = cache[key][:, theta_index, k]
            total = int(hist.sum())
            if total:
                ax.stairs(hist / total, edges, color=color, linewidth=1.7,
                          label=f"{label} (N={total:,})")
        ax.axvspan(.73, 1.03, color="gray", alpha=.11, label="Nominal proton window")
        ax.axvline(.9382720813**2, color="black", linestyle=":", linewidth=1)
        ax.set_title(f"Predicted p = {cache['p_edges'][k]:.2f}–{cache['p_edges'][k+1]:.2f} GeV")
        ax.set_ylabel("Fraction / 0.01 GeV²")
        ax.set_xlabel("Reconstructed missing mass² [GeV²]")
        ax.set_xlim(-.2, 1.5)
        ax.grid(alpha=.18)
        ax.legend(fontsize=8)
    low, high = cache["theta_edges"][theta_index:theta_index + 2]
    fig.suptitle(f"Predicted proton θ = {low:g}–{high:g}°: truth-conditioned mass", fontsize=15)
    fig.tight_layout()
    fig.savefig(output / f"{theta_index-3}_truth_near_shapes_theta_{low:g}_{high:g}.png", dpi=150)
    plt.close(fig)


def plot_summary(cache, output):
    edges = cache["m2_edges"]
    centers = (edges[:-1] + edges[1:]) / 2
    low = (centers >= .43) & (centers < .73)
    peak = (centers >= .73) & (centers < 1.03)
    all_mc = cache["mc_all"]
    near = cache["mc_truth_valid_near"]
    data = cache["data_all"]
    metrics = (
        ("Data: reconstructed proton-window fraction", data[peak].sum(axis=0), data.sum(axis=0)),
        ("MC: reconstructed proton-window fraction", all_mc[peak].sum(axis=0), all_mc.sum(axis=0)),
        ("Truth-near MC: reconstructed proton-window fraction", near[peak].sum(axis=0), near.sum(axis=0)),
        ("Truth-near fraction of MC low-mass shoulder", near[low].sum(axis=0), all_mc[low].sum(axis=0)),
    )
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharex=True, sharey=True)
    for ax, (title, numerator, denominator) in zip(axes.flat, metrics):
        values = np.divide(numerator, denominator, out=np.full(denominator.shape, np.nan),
                           where=denominator >= 100)
        image = ax.pcolormesh(cache["p_edges"], cache["theta_edges"], values,
                              cmap="viridis", vmin=0, vmax=1, shading="flat")
        for j in range(values.shape[0]):
            for k in range(values.shape[1]):
                if np.isfinite(values[j, k]):
                    ax.text((cache["p_edges"][k] + cache["p_edges"][k+1])/2,
                            (cache["theta_edges"][j] + cache["theta_edges"][j+1])/2,
                            f"{values[j,k]:.2f}", ha="center", va="center", fontsize=8,
                            color="white" if values[j,k] < .5 else "black")
        ax.set_title(title)
        ax.set_xlabel("Predicted proton momentum [GeV]")
        ax.set_ylabel("Predicted proton θ [degrees]")
        fig.colorbar(image, ax=ax, label="Fraction")
    fig.suptitle("Saved-mass-range tags; blank cells have fewer than 100 denominator tags")
    fig.tight_layout()
    fig.savefig(output / "0_mass_window_and_truth_summary.png", dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", required=True, type=Path,
                        help="9_mass_shape_histograms.npz from mass_shape_diagnostics.py")
    parser.add_argument("--output", type=Path,
                        help="Default: truth_conditioned_v1 beside the cache")
    args = parser.parse_args()
    cache = read_cache(args.cache)
    output = args.output or args.cache.parent / "truth_conditioned_v1"
    output.mkdir(parents=True, exist_ok=True)
    plot_summary(cache, output)
    for j, low in enumerate(cache["theta_edges"][:-1]):
        if low in (40, 50, 60):
            plot_shapes(cache, output, j)
    write_metrics(cache, output)
    print(f"[OUTPUT] {output}")
    print("[NOTE] Truth-near means truth missing mass² within proton mass² ±0.15 GeV².")
    print("[NOTE] Shapes are normalized independently; these are diagnostics, not efficiencies.")


if __name__ == "__main__":
    main()
