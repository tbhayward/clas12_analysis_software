#!/usr/bin/env python3
"""First-pass, fully integrated e pi+ pi- tag-and-probe study.

Run from any directory. Input ROOT parts must contain the campaign's Tags tree.
The window and matching thresholds are exploratory, not calibrated selections.
"""

import argparse
from array import array
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import json
import math
from pathlib import Path
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares
from scipy.special import ndtr

FIELDS = ("missing_M2", "missing_theta", "missing_p", "best_selected_rec_index",
          "best_selected_match_angle", "best_selected_relative_delta_p")
PROTON_MASS_GEV = 0.9382720813
P_MAX_GEV = math.sqrt(1.0 + 1.0 / (4.0 * PROTON_MASS_GEV**2))
EDGES = {
    "missing_M2": np.linspace(-1.0, 3.0, 161),
    "missing_theta": np.linspace(0.0, 130.0, 131),
    "missing_p": np.linspace(0.0, P_MAX_GEV, 115),
    "best_selected_match_angle": np.linspace(0.0, 30.0, 121),
    "best_selected_relative_delta_p": np.linspace(-1.0, 1.0, 121),
}
LABELS = {
    "missing_M2": r"$M_X^2(e\pi^+\pi^-)$ (GeV$^2$)",
    "missing_theta": r"Predicted proton $\theta$ (degrees)",
    "missing_p": r"Predicted proton momentum (GeV)",
    "best_selected_match_angle": "Best selected proton opening angle (degrees)",
    "best_selected_relative_delta_p": r"$(p_{rec}-p_{miss})/p_{miss}$",
}
THETA_BINS = np.array([0, 10, 20, 30, 40, 50, 60, 70, 180], dtype=float)
MOMENTUM_BINS = np.array([0, .2, .4, .6, .8, 1., P_MAX_GEV], dtype=float)
GRID_STAGES = ("window", "matched", "peak", "peak_matched", "side_left",
               "side_left_matched", "side_right", "side_right_matched")
FIT_EDGES = np.linspace(0.43, 1.33, 91)
# The broad regions are diagnostics; no fit result is used as a production weight.
FIT_REGIONS = (
    ("visible", 0.0, 70.0, 0.2, P_MAX_GEV),
    ("forward_0_40", 0.0, 40.0, 0.2, P_MAX_GEV),
    ("cd_40_70_low_p", 40.0, 70.0, 0.2, 0.7),
    ("cd_40_70_high_p", 40.0, 70.0, 0.7, P_MAX_GEV),
    ("high_angle_50_70", 50.0, 70.0, 0.2, P_MAX_GEV),
)


def empty_result():
    return {"input": 0, "all": 0, "window": 0, "selected": 0, "matched": 0,
            "grid": {stage: np.zeros((len(THETA_BINS) - 1, len(MOMENTUM_BINS) - 1),
                                    dtype=np.int64) for stage in GRID_STAGES},
            "fit_hist": {region[0]: {stage: np.zeros(len(FIT_EDGES) - 1, dtype=np.int64)
                                     for stage in ("all", "matched")}
                         for region in FIT_REGIONS},
            "hist": {name: {stage: np.zeros(len(edges) - 1, dtype=np.int64)
                            for stage in (("all", "window", "matched", "matched_all")
                                          if name == "missing_M2" else ("all", "window", "matched"))}
                     for name, edges in EDGES.items()}}


def process_file(task):
    filename, lo, hi, max_angle, max_rel_dp, chunk_size, entry_stop = task
    import ROOT

    ROOT.TH1.AddDirectory(False)
    result = empty_result()
    started = time.perf_counter()
    source = ROOT.TFile.Open(filename, "READ")
    if not source or source.IsZombie():
        raise RuntimeError(f"Cannot open ROOT file: {filename}")
    opened = time.perf_counter()
    try:
        tree = source.Get("Tags")
        if not tree:
            raise RuntimeError(f"Missing Tags tree: {filename}")
        missing = set(FIELDS) - {branch.GetName() for branch in tree.GetListOfBranches()}
        if missing:
            raise RuntimeError(f"{filename}: missing Tags branches {sorted(missing)}")
        tree_ready = time.perf_counter()
        input_frame = ROOT.RDataFrame(tree)
        if entry_stop is not None:
            input_frame = input_frame.Range(0, min(entry_stop, int(tree.GetEntries())))
        frame = input_frame.Filter(f"std::isfinite(missing_p) && missing_p >= 0.0 "
                                   f"&& missing_p < {P_MAX_GEV:.17g}")
        window = frame.Filter(f"std::isfinite(missing_M2) && missing_M2 >= {lo:.17g} "
                              f"&& missing_M2 <= {hi:.17g}")
        selected = window.Filter("best_selected_rec_index >= 0")
        match_cut = ("std::isfinite(best_selected_match_angle) && "
                     f"best_selected_match_angle < {max_angle:.17g}")
        if max_rel_dp is not None:
            match_cut += (" && std::isfinite(best_selected_relative_delta_p) && "
                          f"std::abs(best_selected_relative_delta_p) < {max_rel_dp:.17g}")
        matched = selected.Filter(match_cut)
        stages = {"all": frame, "window": window, "matched": matched}
        stages["matched_all"] = frame.Filter("best_selected_rec_index >= 0").Filter(match_cut)
        peak_half_width = 0.15
        peak_center = 0.9382720813 ** 2
        regions = {"peak": (peak_center - peak_half_width, peak_center + peak_half_width),
                   "side_left": (peak_center - 3 * peak_half_width,
                                 peak_center - peak_half_width),
                   "side_right": (peak_center + peak_half_width,
                                  peak_center + 3 * peak_half_width)}
        for name, (low, high) in regions.items():
            region = frame.Filter(f"std::isfinite(missing_M2) && missing_M2 >= {low:.17g} "
                                  f"&& missing_M2 < {high:.17g}")
            stages[name] = region
            stages[name + "_matched"] = region.Filter("best_selected_rec_index >= 0").Filter(match_cut)
        counts = {"input": input_frame.Count(), "all": frame.Count(), "window": window.Count(),
                  "selected": selected.Count(), "matched": matched.Count()}
        histograms = {}
        for name, edges in EDGES.items():
            histogram_stages = (("all", "window", "matched", "matched_all")
                                if name == "missing_M2" else ("all", "window", "matched"))
            for stage in histogram_stages:
                node = stages[stage]
                model = (f"h_{name}_{stage}", "", len(edges) - 1,
                         float(edges[0]), float(edges[-1]))
                histograms[name, stage] = node.Filter(f"std::isfinite({name})").Histo1D(model, name)
        x_bins = array("d", THETA_BINS)
        y_bins = array("d", MOMENTUM_BINS)
        grid_histograms = {}
        for stage in GRID_STAGES:
            model = ROOT.RDF.TH2DModel(f"grid_{stage}", "", len(x_bins) - 1, x_bins,
                                       len(y_bins) - 1, y_bins)
            grid_histograms[stage] = stages[stage].Filter(
                "std::isfinite(missing_theta) && std::isfinite(missing_p)").Histo2D(
                    model, "missing_theta", "missing_p")
        fit_histograms = {}
        for region_name, theta_low, theta_high, p_low, p_high in FIT_REGIONS:
            region_cut = (f"std::isfinite(missing_theta) && missing_theta >= {theta_low:.17g} "
                          f"&& missing_theta < {theta_high:.17g} && missing_p >= {p_low:.17g} "
                          f"&& missing_p < {p_high:.17g} && std::isfinite(missing_M2)")
            for stage in ("all", "matched"):
                model = (f"fit_{region_name}_{stage}", "", len(FIT_EDGES) - 1,
                         float(FIT_EDGES[0]), float(FIT_EDGES[-1]))
                # The displayed matched-rate histogram is mass-window selected;
                # a mass-spectrum fit must instead use matched tags before that cut.
                fit_source = stages["matched_all"] if stage == "matched" else stages["all"]
                fit_histograms[region_name, stage] = fit_source.Filter(region_cut).Histo1D(
                    model, "missing_M2")
        ROOT.RDF.RunGraphs(list(counts.values()) + list(histograms.values()) +
                           list(grid_histograms.values()) + list(fit_histograms.values()))
        for stage, handle in counts.items():
            result[stage] = int(handle.GetValue())
        for (name, stage), handle in histograms.items():
            hist = handle.GetValue()
            result["hist"][name][stage] = np.fromiter(
                (int(hist.GetBinContent(i)) for i in range(1, hist.GetNbinsX() + 1)),
                dtype=np.int64, count=hist.GetNbinsX())
        for stage, handle in grid_histograms.items():
            hist = handle.GetValue()
            result["grid"][stage] = np.array(
                [[int(hist.GetBinContent(i + 1, j + 1)) for j in range(len(y_bins) - 1)]
                 for i in range(len(x_bins) - 1)], dtype=np.int64)
        for (region_name, stage), handle in fit_histograms.items():
            hist = handle.GetValue()
            result["fit_hist"][region_name][stage] = np.fromiter(
                (int(hist.GetBinContent(i)) for i in range(1, hist.GetNbinsX() + 1)),
                dtype=np.int64, count=hist.GetNbinsX())
    finally:
        source.Close()
    finished = time.perf_counter()
    result["timing"] = {"open_s": opened - started,
                        "tree_lookup_s": tree_ready - opened,
                        "read_and_hist_s": finished - tree_ready,
                        "total_s": finished - started}
    return result


def add_result(total, part):
    for key in ("input", "all", "window", "selected", "matched"):
        total[key] += part[key]
    for name in EDGES:
        for stage in total["hist"][name]:
            total["hist"][name][stage] += part["hist"][name][stage]
    for stage in GRID_STAGES:
        total["grid"][stage] += part["grid"][stage]
    for region_name, _, _, _, _ in FIT_REGIONS:
        for stage in ("all", "matched"):
            total["fit_hist"][region_name][stage] += part["fit_hist"][region_name][stage]


def plot_counts(path, name, histograms, logy=False, title=None):
    edges = EDGES[name]
    centers = (edges[1:] + edges[:-1]) / 2
    fig, ax = plt.subplots(figsize=(9, 6))
    for label, values in histograms:
        ax.step(centers, values, where="mid", label=label)
    ax.set(xlabel=LABELS[name], ylabel="Tags / bin", title=title)
    if logy:
        ax.set_yscale("log")
        ax.set_ylim(bottom=0.8)
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def ratio(numerator, denominator):
    value = np.divide(numerator, denominator, out=np.full_like(numerator, np.nan, dtype=float),
                      where=denominator > 0)
    binomial = (denominator > 0) & (value >= 0) & (value <= 1)
    variance = np.divide(value * (1 - value), denominator,
                         out=np.full_like(value, np.nan), where=binomial)
    error = np.sqrt(variance)
    return value, error


def plot_ratio(path, name, data, mc, title):
    edges = EDGES[name]
    centers = (edges[1:] + edges[:-1]) / 2
    fig, ax = plt.subplots(figsize=(9, 6))
    for label, result in (("Data", data), ("MC", mc)):
        eff, err = ratio(result["hist"][name]["matched"], result["hist"][name]["window"])
        ax.errorbar(centers, eff, yerr=err, fmt=".", markersize=3, capsize=0,
                    label=label)
    ax.set(xlabel=LABELS[name], ylabel="Matched / missing-mass-window tags", title=title)
    ax.set_ylim(0, 1.05)
    ax.set_xlim((0.2, P_MAX_GEV) if name == "missing_p" else (0, 70))
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_mass_overview(path, results):
    edges = EDGES["missing_M2"]
    centers = (edges[1:] + edges[:-1]) / 2
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharex=True)
    for ax, stage, title in zip(axes, ("all", "matched_all"),
                                ("All tags", "Matched protons")):
        for mode in ("data", "mc"):
            values = results[mode]["hist"]["missing_M2"][stage]
            ax.step(centers, values * 1e6 / results[mode]["input"], where="mid",
                    label=mode.upper())
        ax.axvspan(0.9382720813**2 - .15, 0.9382720813**2 + .15,
                   color="gray", alpha=.2, label="Peak" if stage == "all" else None)
        ax.axvspan(0.9382720813**2 - .45, 0.9382720813**2 - .15,
                   color="orange", alpha=.12)
        ax.axvspan(0.9382720813**2 + .15, 0.9382720813**2 + .45,
                   color="orange", alpha=.12)
        ax.set(xlim=(0.3, 1.5), xlabel=LABELS["missing_M2"],
               ylabel="Tags per million input tags / bin", title=title)
        ax.grid(alpha=.2)
        ax.legend()
    fig.suptitle("Peak (gray) and equal-total-width sidebands (orange)")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_match_angle_scan(path, results, chosen_angle):
    edges = EDGES["best_selected_match_angle"]
    fig, ax = plt.subplots(figsize=(9, 6))
    for mode in ("data", "mc"):
        selected = results[mode]["selected"]
        histogram = results[mode]["hist"]["best_selected_match_angle"]["window"]
        survival = np.cumsum(histogram) / selected if selected else np.zeros(len(histogram))
        ax.step(edges[1:], survival, where="post", label=mode.upper())
    ax.axvline(chosen_angle, color="black", linestyle="--", label=f"Current {chosen_angle:g}° cut")
    ax.set(xlim=(0, 30), ylim=(0, 1), xlabel="Best selected proton opening angle cut (degrees)",
           ylabel="Fraction of selected protons retained", title="Opening-angle cut scan in mass window")
    ax.grid(alpha=.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def fit_mass_peak(counts):
    """Exploratory Gaussian signal plus positive quadratic Bernstein background."""
    observed = counts.astype(float)
    centers = (FIT_EDGES[:-1] + FIT_EDGES[1:]) / 2
    coordinate = (centers - FIT_EDGES[0]) / (FIT_EDGES[-1] - FIT_EDGES[0])
    background_basis = np.stack(((1 - coordinate)**2,
                                 2 * coordinate * (1 - coordinate), coordinate**2))
    total = float(observed.sum())
    if total < 200:
        return {"status": "insufficient_counts", "entries": int(total)}
    scale = max(total / len(observed), 1.0)
    weights = np.sqrt(np.maximum(observed, 1.0))

    def components(parameters):
        yield_scaled, mean, width, *background_scaled = parameters
        signal_shape = np.diff(ndtr((FIT_EDGES - mean) / width))
        signal_shape /= signal_shape.sum()
        signal = scale * yield_scaled * signal_shape
        background = scale * np.asarray(background_scaled) @ background_basis
        return signal, background

    edge_left = float(np.median(observed[:10])) / scale
    edge_right = float(np.median(observed[-10:])) / scale
    peak_mean = float(np.clip(centers[np.argmax(observed)], .825, .955))
    initial_yield = max((total - (edge_left + edge_right) * scale * len(observed) / 2) / scale,
                        1.0)
    lower = [0.0, .82, .025, 0.0, 0.0, 0.0]
    upper = [2.0 * total / scale, .96, .18, 5.0 * observed.max() / scale,
             5.0 * observed.max() / scale, 5.0 * observed.max() / scale]
    solutions = []
    for start_width in (.045, .075, .115):
        initial = [initial_yield, peak_mean, start_width, edge_left,
                   (edge_left + edge_right) / 2, edge_right]
        solution = least_squares(
            lambda parameters: (sum(components(parameters)) - observed) / weights,
            initial, bounds=(lower, upper), x_scale="jac", max_nfev=1500)
        if solution.success and np.all(np.isfinite(solution.x)):
            solutions.append(solution)
    if not solutions:
        return {"status": "optimizer_failed", "entries": int(total)}
    best = min(solutions, key=lambda solution: float(np.dot(solution.fun, solution.fun)))
    signal, background = components(best.x)
    yield_value = float(signal.sum())
    chi2 = float(np.dot(best.fun, best.fun))
    ndf = len(observed) - len(best.x)
    chi2_ndf = chi2 / ndf
    quality = ("review_shape" if chi2_ndf > 3.0 else
               "weak_signal" if yield_value < max(100.0, 5.0 * math.sqrt(total)) else "ok")
    return {"status": quality, "entries": int(total), "signal_yield": yield_value,
            "mean_GeV2": float(best.x[1]), "width_GeV2": float(best.x[2]),
            "chi2": chi2, "ndf": ndf, "chi2_ndf": chi2_ndf,
            "signal_curve": signal.tolist(), "background_curve": background.tolist(),
            "total_curve": (signal + background).tolist()}


def plot_fit_region(path, region, results, fitted):
    name, theta_low, theta_high, p_low, p_high = region
    centers = (FIT_EDGES[:-1] + FIT_EDGES[1:]) / 2
    fig, axes = plt.subplots(4, 2, figsize=(13, 11), sharex=True,
                             gridspec_kw={"height_ratios": [3, 1, 3, 1]})
    for column, mode in enumerate(("data", "mc")):
        for stage, row in (("all", 0), ("matched", 2)):
            observed = results[mode]["fit_hist"][name][stage]
            fit = fitted[mode][stage]
            main, residual = axes[row, column], axes[row + 1, column]
            main.errorbar(centers, observed, yerr=np.sqrt(observed), fmt=".",
                          markersize=3, color="black", label="Observed")
            if "total_curve" in fit:
                model = np.asarray(fit["total_curve"])
                main.plot(centers, model, color="tab:blue", label="Signal + background")
                main.plot(centers, fit["background_curve"], color="tab:orange",
                          linestyle="--", label="Background")
                residual.plot(centers, (observed - model) / np.sqrt(np.maximum(observed, 1)),
                              ".", markersize=3, color="black")
                detail = (f"S={fit['signal_yield']:.0f}, σ={fit['width_GeV2']:.3f}, "
                          f"χ²/ndf={fit['chi2_ndf']:.1f}; {fit['status']}")
            else:
                detail = fit["status"]
            main.set_title(f"{mode.upper()} {stage}: {detail}", fontsize=10)
            main.set_ylabel("Tags / bin")
            main.grid(alpha=.2)
            residual.axhline(0, color="tab:blue", linewidth=1)
            residual.set_ylabel("Pull")
            residual.grid(alpha=.2)
            if row == 0 and column == 0:
                main.legend(fontsize=8)
        axes[3, column].set_xlabel(LABELS["missing_M2"])
        axes[3, column].set_xlim(FIT_EDGES[0], FIT_EDGES[-1])
    fig.suptitle(f"{name}: {theta_low:g}–{theta_high:g}°; {p_low:.2f}–{p_high:.2f} GeV")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def save_fit_products(output, period, results):
    rows = []
    details = {"period": period, "model": "Gaussian + positive quadratic Bernstein background",
               "fit_M2_GeV2": [float(FIT_EDGES[0]), float(FIT_EDGES[-1])],
               "notes": "Independent all/matched fits; exploratory ratios have no fit uncertainty yet.",
               "regions": {}}
    for index, region in enumerate(FIT_REGIONS, 9):
        name, theta_low, theta_high, p_low, p_high = region
        fitted = {mode: {stage: fit_mass_peak(results[mode]["fit_hist"][name][stage])
                         for stage in ("all", "matched")} for mode in ("data", "mc")}
        plot_fit_region(output / f"{index}_fit_{name}.png", region, results, fitted)
        rates = {}
        for mode in ("data", "mc"):
            all_fit, matched_fit = fitted[mode]["all"], fitted[mode]["matched"]
            if "signal_yield" in all_fit and "signal_yield" in matched_fit and \
                    all_fit["signal_yield"] > 0 and \
                    0 <= matched_fit["signal_yield"] <= all_fit["signal_yield"]:
                rates[mode] = matched_fit["signal_yield"] / all_fit["signal_yield"]
            else:
                rates[mode] = None
        scale = (rates["data"] / rates["mc"]
                 if rates["data"] is not None and rates["mc"] else None)
        quality = all(fitted[mode][stage]["status"] == "ok"
                      for mode in ("data", "mc") for stage in ("all", "matched"))
        details["regions"][name] = {
            "theta_deg": [theta_low, theta_high], "predicted_p_GeV": [p_low, p_high],
            "fits": {mode: {stage: {key: value for key, value in fitted[mode][stage].items()
                                    if not key.endswith("_curve")}
                            for stage in ("all", "matched")} for mode in ("data", "mc")},
            "rate": rates, "data_mc_ratio": scale, "all_four_fits_ok": quality}
        rows.append((name, theta_low, theta_high, p_low, p_high,
                     rates["data"], rates["mc"], scale, int(quality),
                     *(fitted[mode][stage].get("signal_yield")
                       for mode in ("data", "mc") for stage in ("all", "matched")),
                     *(fitted[mode][stage].get("chi2_ndf")
                       for mode in ("data", "mc") for stage in ("all", "matched"))))
        print(f"[fit] {name}: data/MC={scale if scale is not None else 'unavailable'} "
              f"all_four_fits_ok={quality}", flush=True)
    with (output / "14_fit_regions.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("region", "theta_low_deg", "theta_high_deg", "p_low_GeV", "p_high_GeV",
                         "data_rate", "mc_rate", "data_mc_ratio", "all_four_fits_ok",
                         "data_all_signal", "data_matched_signal", "mc_all_signal",
                         "mc_matched_signal", "data_all_chi2_ndf", "data_matched_chi2_ndf",
                         "mc_all_chi2_ndf", "mc_matched_chi2_ndf"))
        writer.writerows(rows)
    (output / "15_fit_regions.json").write_text(json.dumps(details, indent=2) + "\n")
    return details


def save_grid(path, values, title, color_label, vmin=None, vmax=None, log=False):
    from matplotlib.colors import LogNorm

    fig, ax = plt.subplots(figsize=(11, 7))
    image = ax.pcolormesh(THETA_BINS, MOMENTUM_BINS, values.T,
                         cmap="viridis", shading="flat",
                         norm=LogNorm(vmin=1) if log else None,
                         vmin=None if log else vmin, vmax=None if log else vmax)
    fig.colorbar(image, ax=ax, label=color_label)
    ax.set(xlabel=LABELS["missing_theta"], ylabel=LABELS["missing_p"], title=title)
    ax.set_xlim(0, 70)
    ax.set_ylim(0.2, P_MAX_GEV)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def save_grid_products(output, period, results):
    rates = {}
    corrected = {}
    for mode in ("data", "mc"):
        grid = results[mode]["grid"]
        rates[mode], _ = ratio(grid["matched"], grid["window"])
        corrected[mode] = {
            "denominator": grid["peak"].astype(float) -
            0.5 * (grid["side_left"] + grid["side_right"]),
            "numerator": grid["peak_matched"].astype(float) -
            0.5 * (grid["side_left_matched"] + grid["side_right_matched"]),
        }
        corrected[mode]["rate"], _ = ratio(corrected[mode]["numerator"],
                                           corrected[mode]["denominator"])
        corrected[mode]["rate"] = np.where(
            (corrected[mode]["numerator"] > 0) &
            (corrected[mode]["numerator"] <= corrected[mode]["denominator"]),
            corrected[mode]["rate"], np.nan)
    data, mc = results["data"]["grid"], results["mc"]["grid"]
    valid = ((data["window"] >= 100) & (mc["window"] >= 100) &
             (data["matched"] >= 20) & (mc["matched"] >= 20))
    raw_scale = np.divide(rates["data"], rates["mc"],
                          out=np.full_like(rates["data"], np.nan),
                          where=valid & (rates["mc"] > 0))
    save_grid(output / "5_raw_theta_p_data_mc.png", raw_scale,
              f"{period}: raw window data / MC rate",
              "Data / MC", vmin=0.5, vmax=2.0)
    corrected_valid = (valid &
                       (corrected["data"]["denominator"] >= 100) &
                       (corrected["mc"]["denominator"] >= 100) &
                       (corrected["data"]["numerator"] >= 20) &
                       (corrected["mc"]["numerator"] >= 20) &
                       np.isfinite(corrected["data"]["rate"]) &
                       np.isfinite(corrected["mc"]["rate"]))
    corrected_scale = np.divide(corrected["data"]["rate"], corrected["mc"]["rate"],
                                out=np.full_like(raw_scale, np.nan),
                                where=corrected_valid & np.isfinite(corrected["mc"]["rate"]) &
                                      (corrected["mc"]["rate"] > 0))
    save_grid(output / "6_sideband_theta_p_data_mc.png", corrected_scale,
              f"{period}: sideband-subtracted data / MC rate", "Data / MC",
              vmin=0.5, vmax=2.0)
    with (output / "7_theta_p_cells.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("theta_low_deg", "theta_high_deg", "p_low_GeV", "p_high_GeV",
                         "data_window", "data_matched", "data_rate", "mc_window",
                         "mc_matched", "mc_rate", "raw_data_mc_ratio", "raw_plot_quality",
                         "data_peak", "data_left", "data_right", "data_peak_matched",
                         "data_left_matched", "data_right_matched", "data_signal_denominator",
                         "data_signal_numerator", "data_sideband_rate", "mc_peak", "mc_left",
                         "mc_right", "mc_peak_matched", "mc_left_matched", "mc_right_matched",
                         "mc_signal_denominator", "mc_signal_numerator", "mc_sideband_rate",
                         "sideband_data_mc_ratio", "sideband_plot_quality"))
        for i in range(len(THETA_BINS) - 1):
            for j in range(len(MOMENTUM_BINS) - 1):
                writer.writerow((THETA_BINS[i], THETA_BINS[i + 1],
                                 MOMENTUM_BINS[j], MOMENTUM_BINS[j + 1],
                                 data["window"][i, j], data["matched"][i, j],
                                 rates["data"][i, j], mc["window"][i, j],
                                 mc["matched"][i, j], rates["mc"][i, j],
                                 raw_scale[i, j], int(valid[i, j]),
                                 *(data[key][i, j] for key in ("peak", "side_left", "side_right",
                                   "peak_matched", "side_left_matched", "side_right_matched")),
                                 corrected["data"]["denominator"][i, j],
                                 corrected["data"]["numerator"][i, j],
                                 corrected["data"]["rate"][i, j],
                                 *(mc[key][i, j] for key in ("peak", "side_left", "side_right",
                                   "peak_matched", "side_left_matched", "side_right_matched")),
                                 corrected["mc"]["denominator"][i, j],
                                 corrected["mc"]["numerator"][i, j],
                                 corrected["mc"]["rate"][i, j],
                                 corrected_scale[i, j], int(corrected_valid[i, j])))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, default=Path("output/proton_efficiency"))
    parser.add_argument("--period", default="fa18_inb")
    parser.add_argument("--workers", type=int, default=10)
    parser.add_argument("--mass-half-width", type=float, default=0.30)
    parser.add_argument("--match-angle", type=float, default=10.0)
    parser.add_argument("--max-relative-dp", type=float, default=None,
                        help="Optional absolute relative momentum difference cut")
    parser.add_argument("--chunk-size", type=int, default=100000,
                        help="Legacy option; native ROOT scans each file in one pass")
    parser.add_argument("--fast-mode", action="store_true",
                        help="Preview the first 10000 tags from one completed file per mode")
    parser.add_argument("--max-tags-per-mode", type=int, default=0,
                        help="Read at most this many Tags rows each for data and MC; 0 means all")
    args = parser.parse_args()
    if not 1 <= args.workers <= 10 or args.chunk_size < 1:
        parser.error("workers must be 1..10 and chunk-size must be positive")
    if args.mass_half_width <= 0 or args.match_angle <= 0:
        parser.error("mass-half-width and match-angle must be positive")
    if args.max_relative_dp is not None and args.max_relative_dp <= 0:
        parser.error("max-relative-dp must be positive")
    if args.max_tags_per_mode < 0:
        parser.error("max-tags-per-mode must be nonnegative")
    if args.fast_mode and args.max_tags_per_mode:
        parser.error("--fast-mode and --max-tags-per-mode cannot be combined")
    center = 0.9382720813 ** 2
    lo, hi = center - args.mass_half_width, center + args.mass_half_width
    study_name = (("fast_preview" if args.fast_mode else
                   f"limited_{args.max_tags_per_mode}" if args.max_tags_per_mode else "integrated")
                  + "_fit_v2_dvcs")
    output = args.output_root / args.period / study_name
    output.mkdir(parents=True, exist_ok=True)
    results = {}
    files_by_mode = {}
    for mode in ("data", "mc"):
        discover_start = time.perf_counter()
        directory = args.input_root / args.period / mode / "parts"
        files = sorted(path for path in directory.glob("*.root")
                       if Path(str(path)[:-5] + ".ok.json").is_file())
        if not files:
            raise SystemExit(f"No completed ROOT parts for {args.period}/{mode} in {directory}")
        print(f"[{mode}] discovered {len(files)} completed files in "
              f"{time.perf_counter() - discover_start:.1f}s", flush=True)
        if args.fast_mode:
            chosen = next((path for path in files if json.loads(
                Path(str(path)[:-5] + ".ok.json").read_text()).get("selected_tags", 0) >= 10000), None)
            if chosen is None:
                raise SystemExit(f"No completed {args.period}/{mode} part has 10000 tags")
            files = [chosen]
            print(f"[{mode}] FAST PREVIEW: reading first 10000 tags from {files[0].name}",
                  flush=True)
        planned = []
        remaining = args.max_tags_per_mode
        for path in files:
            if args.fast_mode:
                entry_stop = 10000
            elif args.max_tags_per_mode:
                available = int(json.loads(Path(str(path)[:-5] + ".ok.json").read_text())
                                .get("selected_tags", 0))
                if available <= 0:
                    continue
                entry_stop = min(available, remaining)
                remaining -= entry_stop
            else:
                entry_stop = None
            planned.append((str(path), lo, hi, args.match_angle, args.max_relative_dp,
                            args.chunk_size, entry_stop))
            if args.fast_mode or (args.max_tags_per_mode and remaining == 0):
                break
        if args.max_tags_per_mode and remaining:
            raise SystemExit(f"{args.period}/{mode} has fewer than "
                             f"{args.max_tags_per_mode} completed tags; short by {remaining}")
        tasks = planned
        files_by_mode[mode] = len(tasks)
        print(f"[{mode}] planned {len(tasks)} files for "
              f"{args.max_tags_per_mode or 'all available'} tags", flush=True)
        total = empty_result()
        if args.fast_mode:
            print(f"[{mode}] opening and reading {files[0].name}", flush=True)
            part = process_file(tasks[0])
            add_result(total, part)
            print(f"[{mode}] timing: open={part['timing']['open_s']:.1f}s "
                  f"tree_lookup={part['timing']['tree_lookup_s']:.1f}s "
                  f"read+hist={part['timing']['read_and_hist_s']:.1f}s", flush=True)
        else:
            with ProcessPoolExecutor(max_workers=min(args.workers, len(tasks))) as pool:
                futures = {pool.submit(process_file, task): task[0] for task in tasks}
                for index, future in enumerate(as_completed(futures), 1):
                    try:
                        part = future.result()
                        add_result(total, part)
                    except Exception as error:
                        raise RuntimeError(f"Failed reading {futures[future]}: {error}") from error
                    if index % 10 == 0 or index == len(tasks):
                        print(f"[{mode}] {index}/{len(tasks)} files ({100 * index / len(tasks):.1f}%) "
                              f"latest open={part['timing']['open_s']:.1f}s "
                              f"tree_lookup={part['timing']['tree_lookup_s']:.1f}s "
                              f"read+hist={part['timing']['read_and_hist_s']:.1f}s", flush=True)
        results[mode] = total
        print(f"[{mode}] input_tags={total['input']} DVCS_momentum_tags={total['all']} "
              f"window={total['window']} "
              f"selected={total['selected']} matched={total['matched']}", flush=True)
    comparison_start = time.perf_counter()
    plot_mass_overview(output / "1_missing_mass_regions.png", results)
    plot_match_angle_scan(output / "2_match_angle_scan.png", results, args.match_angle)
    plot_ratio(output / "3_momentum_rate.png", "missing_p",
               results["data"], results["mc"], f"{args.period}: raw window rate by momentum")
    plot_ratio(output / "4_theta_rate.png", "missing_theta",
               results["data"], results["mc"], f"{args.period}: raw window rate by angle")
    save_grid_products(output, args.period, results)
    fit_details = save_fit_products(output, args.period, results)
    print(f"[comparison] plotted in {time.perf_counter() - comparison_start:.1f}s", flush=True)
    summary = {"period": args.period, "input_root": str(args.input_root),
               "fast_mode": args.fast_mode,
               "predicted_proton_p_max_GeV": P_MAX_GEV,
               "p_max_definition": "exclusive proton recoil at |t|=1 GeV^2",
               "sideband_peak_M2_GeV2": [0.9382720813**2 - .15, 0.9382720813**2 + .15],
               "sideband_left_M2_GeV2": [0.9382720813**2 - .45, 0.9382720813**2 - .15],
               "sideband_right_M2_GeV2": [0.9382720813**2 + .15, 0.9382720813**2 + .45],
               "max_tags_per_mode": args.max_tags_per_mode,
               "fit_model": fit_details["model"],
               "fit_region_summary": "15_fit_regions.json",
               "files": files_by_mode, "mass_window_GeV2": [lo, hi],
               "match_angle_deg": args.match_angle, "max_relative_dp": args.max_relative_dp,
               "counts": {mode: {key: result[key] for key in ("input", "all", "window", "selected", "matched")}
                          for mode, result in results.items()}}
    for mode in ("data", "mc"):
        n, d = results[mode]["matched"], results[mode]["window"]
        summary["counts"][mode]["provisional_rate"] = n / d if d else None
    a = summary["counts"]["data"]["provisional_rate"]
    b = summary["counts"]["mc"]["provisional_rate"]
    summary["provisional_data_mc_ratio"] = a / b if a is not None and b else None
    summary["sideband_integrated"] = {}
    for mode, result in results.items():
        grid = result["grid"]
        denominator = int(grid["peak"].sum()) - 0.5 * int(
            (grid["side_left"] + grid["side_right"]).sum())
        numerator = int(grid["peak_matched"].sum()) - 0.5 * int(
            (grid["side_left_matched"] + grid["side_right_matched"]).sum())
        summary["sideband_integrated"][mode] = {
            "signal_denominator": denominator, "signal_numerator": numerator,
            "provisional_rate": numerator / denominator if denominator > 0 else None}
    (output / "8_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"[OUTPUT] {output}", flush=True)


if __name__ == "__main__":
    main()
