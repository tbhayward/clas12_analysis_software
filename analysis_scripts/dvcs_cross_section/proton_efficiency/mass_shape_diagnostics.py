#!/usr/bin/env python3
"""Inspect e pi+ pi- missing-mass shapes without assuming a signal fit model.

Read completed campaign ROOT parts once, cache the binned spectra, and compare
all, proton-matched, and nonmatched tags in fine proton-angle/momentum cells.
The outputs are diagnostics, not proton-efficiency corrections.
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
from matplotlib.colors import LogNorm
import numpy as np

PROTON_MASS_GEV = 0.9382720813
P_MAX_GEV = math.sqrt(1.0 + 1.0 / (4.0 * PROTON_MASS_GEV**2))
M2_EDGES = np.linspace(-0.2, 1.8, 201)
THETA_EDGES = np.array([0, 10, 20, 30, 40, 50, 60, 70], dtype=float)
P_EDGES = np.array([0.2, 0.5, 0.7, 0.9, P_MAX_GEV], dtype=float)
MODES = ("data", "mc")
STAGES = ("all", "matched", "not_matched")
SHAPE = (len(M2_EDGES) - 1, len(THETA_EDGES) - 1, len(P_EDGES) - 1)
REQUIRED_BRANCHES = {"missing_M2", "missing_theta", "missing_p",
                     "best_selected_rec_index", "best_selected_match_angle",
                     "best_selected_relative_delta_p"}
MC_TRUTH_BRANCHES = {"truth_tag_match_valid", "truth_missing_M2",
                     "truth_best_index", "truth_best_match_angle"}
TRUTH_STAGES = ("valid_low", "valid_near", "valid_high", "valid_near_aligned")


def read_part(task):
    filename, entry_stop, match_angle, max_relative_dp, mode = task
    import ROOT

    ROOT.TH1.AddDirectory(False)
    source = ROOT.TFile.Open(filename, "READ")
    if not source or source.IsZombie():
        raise RuntimeError(f"Cannot open ROOT file: {filename}")
    try:
        tree = source.Get("Tags")
        if not tree:
            raise RuntimeError(f"Missing Tags tree: {filename}")
        branches = {branch.GetName() for branch in tree.GetListOfBranches()}
        missing = (REQUIRED_BRANCHES | (MC_TRUTH_BRANCHES if mode == "mc" else set())) - branches
        if missing:
            raise RuntimeError(f"Missing Tags branches in {filename}: {sorted(missing)}")
        frame = ROOT.RDataFrame(tree)
        if entry_stop is not None:
            frame = frame.Range(0, min(entry_stop, int(tree.GetEntries())))
        valid = frame.Filter(
            f"std::isfinite(missing_M2) && std::isfinite(missing_theta) && "
            f"std::isfinite(missing_p) && missing_theta >= {THETA_EDGES[0]:.17g} && "
            f"missing_theta < {THETA_EDGES[-1]:.17g} && "
            f"missing_p >= {P_EDGES[0]:.17g} && missing_p < {P_EDGES[-1]:.17g}")
        matched_cut = ("best_selected_rec_index >= 0 && "
                       "std::isfinite(best_selected_match_angle) && "
                       f"best_selected_match_angle < {match_angle:.17g}")
        if max_relative_dp is not None:
            matched_cut += (" && std::isfinite(best_selected_relative_delta_p) && "
                            f"std::abs(best_selected_relative_delta_p) < {max_relative_dp:.17g}")
        matched = valid.Filter(matched_cut)
        m2_bins, theta_bins, p_bins = (array("d", edges) for edges in
                                       (M2_EDGES, THETA_EDGES, P_EDGES))
        model = lambda name: ROOT.RDF.TH3DModel(
            name, "", len(M2_EDGES) - 1, m2_bins,
            len(THETA_EDGES) - 1, theta_bins,
            len(P_EDGES) - 1, p_bins)
        input_count = frame.Count()
        valid_count = valid.Count()
        matched_count = matched.Count()
        all_hist = valid.Histo3D(model("shape_all"), "missing_M2", "missing_theta", "missing_p")
        matched_hist = matched.Histo3D(
            model("shape_matched"), "missing_M2", "missing_theta", "missing_p")
        truth_histograms = {}
        if mode == "mc":
            truth_valid = valid.Filter("truth_tag_match_valid > 0.5 && "
                                       "std::isfinite(truth_missing_M2)")
            center = PROTON_MASS_GEV**2
            low = truth_valid.Filter(f"truth_missing_M2 < {center - .15:.17g}")
            near = truth_valid.Filter(f"truth_missing_M2 >= {center - .15:.17g} && "
                                      f"truth_missing_M2 < {center + .15:.17g}")
            high = truth_valid.Filter(f"truth_missing_M2 >= {center + .15:.17g}")
            aligned = near.Filter("truth_best_index >= 0 && "
                                  "std::isfinite(truth_best_match_angle) && "
                                  f"truth_best_match_angle < {match_angle:.17g}")
            for name, node in (("valid_low", low), ("valid_near", near),
                               ("valid_high", high), ("valid_near_aligned", aligned)):
                truth_histograms[name] = node.Histo3D(
                    model(f"truth_{name}"), "missing_M2", "missing_theta", "missing_p")
        ROOT.RDF.RunGraphs([input_count, valid_count, matched_count, all_hist, matched_hist] +
                           list(truth_histograms.values()))

        def values(handle):
            hist = handle.GetValue()
            return np.array(
                [[[int(hist.GetBinContent(i + 1, j + 1, k + 1))
                   for k in range(SHAPE[2])]
                  for j in range(SHAPE[1])]
                 for i in range(SHAPE[0])], dtype=np.int64)

        all_values, matched_values = values(all_hist), values(matched_hist)
        if np.any(matched_values > all_values):
            raise RuntimeError(f"Matched tags exceed all tags in {filename}")
        truth_values = {name: values(handle) for name, handle in truth_histograms.items()}
        return {"input": int(input_count.GetValue()), "selected": int(valid_count.GetValue()),
                "matched_count": int(matched_count.GetValue()), "all": all_values,
                "matched": matched_values, "truth": truth_values,
                "branches": sorted(branches)}
    finally:
        source.Close()


def completed_files(directory):
    return sorted(path for path in directory.glob("*.root")
                  if path.with_suffix(".ok.json").is_file())


def scan_mode(input_root, period, mode, workers, max_tags, match_angle, max_relative_dp):
    files = completed_files(input_root / period / mode / "parts")
    if not files:
        raise RuntimeError(f"No completed ROOT parts for {period}/{mode}")
    print(f"[{mode}] discovered {len(files)} completed files", flush=True)
    tasks = []
    remaining = max_tags
    for path in files:
        if max_tags:
            metadata = json.loads(path.with_suffix(".ok.json").read_text())
            available = int(metadata.get("selected_tags", 0))
            if available <= 0:
                continue
            stop = min(available, remaining)
            remaining -= stop
        else:
            stop = None
        tasks.append((str(path), stop, match_angle, max_relative_dp, mode))
        if max_tags and remaining == 0:
            break
    if max_tags and remaining:
        raise RuntimeError(f"{period}/{mode} has fewer than {max_tags} completed tags")
    total = {"input": 0, "selected": 0, "matched_count": 0,
             "all": np.zeros(SHAPE, dtype=np.int64),
             "matched": np.zeros(SHAPE, dtype=np.int64), "branches": set(),
             "files": len(tasks)}
    if mode == "mc":
        total["truth"] = {name: np.zeros(SHAPE, dtype=np.int64) for name in TRUTH_STAGES}
    print(f"[{mode}] scanning {len(tasks)} files with {workers} workers", flush=True)
    with ProcessPoolExecutor(max_workers=min(workers, len(tasks))) as pool:
        futures = {pool.submit(read_part, task): task[0] for task in tasks}
        for index, future in enumerate(as_completed(futures), 1):
            try:
                part = future.result()
            except Exception as error:
                raise RuntimeError(f"Failed reading {futures[future]}: {error}") from error
            for key in ("input", "selected", "matched_count"):
                total[key] += part[key]
            total["all"] += part["all"]
            total["matched"] += part["matched"]
            if mode == "mc":
                for name in TRUTH_STAGES:
                    total["truth"][name] += part["truth"][name]
            total["branches"].update(part["branches"])
            if index % 10 == 0 or index == len(tasks):
                print(f"[{mode}] {index}/{len(tasks)} files ({100 * index / len(tasks):.1f}%)",
                      flush=True)
    total["not_matched"] = total["all"] - total["matched"]
    if mode == "mc":
        classified = sum(total["truth"][name] for name in TRUTH_STAGES[:3])
        total["truth"]["invalid_or_nonfinite"] = total["all"] - classified
        if np.any(total["truth"]["invalid_or_nonfinite"] < 0):
            raise RuntimeError("MC truth categories exceed all tags")
    print(f"[{mode}] input_tags={total['input']} selected_tags={total['selected']} "
          f"matched_tags={total['matched_count']}", flush=True)
    return total


def cache_path(output):
    return output / "9_mass_shape_histograms.npz"


def save_cache(output, results, match_angle, max_relative_dp):
    arrays = {"m2_edges": M2_EDGES, "theta_edges": THETA_EDGES, "p_edges": P_EDGES,
              "match_angle": np.array(match_angle),
              "max_relative_dp": np.array(np.nan if max_relative_dp is None else max_relative_dp)}
    for mode in MODES:
        for stage in STAGES:
            arrays[f"{mode}_{stage}"] = results[mode][stage]
        for key in ("input", "selected", "matched_count", "files"):
            arrays[f"{mode}_{key}_count"] = np.array(results[mode][key])
    for name, values in results["mc"]["truth"].items():
        arrays[f"mc_truth_{name}"] = values
    np.savez_compressed(cache_path(output), **arrays)


def load_cache(path):
    with np.load(path, allow_pickle=False) as cache:
        for key, expected in (("m2_edges", M2_EDGES), ("theta_edges", THETA_EDGES),
                              ("p_edges", P_EDGES)):
            if not np.array_equal(cache[key], expected):
                raise RuntimeError(f"Cache has incompatible {key}: {path}")
        results = {}
        for mode in MODES:
            result = {stage: cache[f"{mode}_{stage}"].astype(np.int64) for stage in STAGES}
            for stage in STAGES:
                if result[stage].shape != SHAPE:
                    raise RuntimeError(f"Cache has incompatible {mode}/{stage} shape")
            for key in ("input", "selected", "matched_count", "files"):
                result[key] = int(cache[f"{mode}_{key}_count"])
            result["branches"] = set()
            if mode == "mc":
                result["truth"] = {name: cache[f"mc_truth_{name}"].astype(np.int64)
                                   for name in (*TRUTH_STAGES, "invalid_or_nonfinite")}
            results[mode] = result
        return results, float(cache["match_angle"]), float(cache["max_relative_dp"])


def shape_line(ax, result, stage, theta_index, p_indices, label, color):
    histogram = result[stage][:, theta_index, p_indices].sum(axis=1)
    count = int(histogram.sum())
    if count:
        ax.step(M2_EDGES[:-1], histogram / count, where="post", color=color,
                label=f"{label} (N={count:,})")
    return count


def plot_high_angle(output, results, stage, number):
    # Three high-angle slices, each shown below and above 0.7 GeV.
    fig, axes = plt.subplots(3, 2, figsize=(13, 10), sharex=True)
    for row, theta_index in enumerate((4, 5, 6)):
        for col, p_indices in enumerate(((0, 1), (2, 3))):
            ax = axes[row, col]
            for mode, color in (("data", "tab:blue"), ("mc", "tab:orange")):
                shape_line(ax, results[mode], stage, theta_index, p_indices,
                           mode.upper(), color)
            p_low, p_high = P_EDGES[p_indices[0]], P_EDGES[p_indices[-1] + 1]
            ax.set_title(f"{THETA_EDGES[theta_index]:g}–{THETA_EDGES[theta_index+1]:g}°; "
                         f"p={p_low:.2f}–{p_high:.2f} GeV")
            ax.set_ylabel("Fraction / 0.01 GeV²")
            ax.grid(alpha=.2)
            ax.legend(fontsize=8)
    for ax in axes[-1]:
        ax.set_xlabel(r"$M_X^2(e\pi^+\pi^-)$ [GeV$^2$]")
        ax.set_xlim(M2_EDGES[0], M2_EDGES[-1])
    fig.suptitle(f"{stage.replace('_', ' ')} tags: data and MC mass shapes", fontsize=15)
    fig.tight_layout()
    fig.savefig(output / f"{number}_{stage}_high_angle_mass.png", dpi=150)
    plt.close(fig)


def plot_mass_angle_maps(output, results, p_index, number):
    fig, axes = plt.subplots(3, 2, figsize=(13, 12), sharex=True, sharey=True)
    for row, stage in enumerate(STAGES):
        for col, mode in enumerate(MODES):
            ax = axes[row, col]
            counts = results[mode][stage][:, :, p_index]
            totals = counts.sum(axis=0)
            fractions = np.divide(counts, totals[None, :],
                                  out=np.full(counts.shape, np.nan), where=totals[None, :] >= 100)
            image = ax.pcolormesh(M2_EDGES, THETA_EDGES, fractions.T,
                                  cmap="viridis", norm=LogNorm(vmin=1e-4, vmax=1e-1),
                                  shading="flat")
            ax.set_title(f"{mode.upper()} {stage.replace('_', ' ')}; N={int(totals.sum()):,}")
            ax.set_ylabel("Predicted proton θ [degrees]")
            fig.colorbar(image, ax=ax, label="Fraction / mass bin")
    for ax in axes[-1]:
        ax.set_xlabel(r"$M_X^2(e\pi^+\pi^-)$ [GeV$^2$]")
        ax.set_xlim(M2_EDGES[0], M2_EDGES[-1])
    fig.suptitle(f"p={P_EDGES[p_index]:.2f}–{P_EDGES[p_index+1]:.2f} GeV; "
                 "white angle rows have fewer than 100 tags")
    fig.tight_layout()
    fig.savefig(output / f"{number}_mass_angle_p{p_index+1}.png", dpi=150)
    plt.close(fig)


def plot_truth_high_angle(output, results):
    mc = results["mc"]
    centers = (M2_EDGES[:-1] + M2_EDGES[1:]) / 2
    fig, axes = plt.subplots(3, 2, figsize=(13, 10), sharex=True)
    categories = ("valid_low", "valid_near", "valid_high", "invalid_or_nonfinite")
    labels = ("Truth M² below proton window", "Truth M² near proton",
              "Truth M² above proton window", "No matched charged truth tag")
    for row, theta_index in enumerate((4, 5, 6)):
        for col, p_indices in enumerate(((0, 1), (2, 3))):
            ax = axes[row, col]
            total = mc["all"][:, theta_index, p_indices].sum(axis=1)
            total_count = int(total.sum())
            if total_count:
                components = [mc["truth"][name][:, theta_index, p_indices].sum(axis=1) /
                              total_count for name in categories]
                ax.stackplot(centers, *components, labels=labels, alpha=.7)
                ax.step(M2_EDGES[:-1], total / total_count, where="post",
                        color="black", linewidth=1, label="All MC tags")
                aligned = (mc["truth"]["valid_near_aligned"][:, theta_index, p_indices]
                           .sum(axis=1) / total_count)
                ax.step(M2_EDGES[:-1], aligned, where="post", color="red",
                        linestyle="--", linewidth=1, label="Truth-near, proton aligned")
            p_low, p_high = P_EDGES[p_indices[0]], P_EDGES[p_indices[-1] + 1]
            ax.set_title(f"{THETA_EDGES[theta_index]:g}–{THETA_EDGES[theta_index+1]:g}°; "
                         f"p={p_low:.2f}–{p_high:.2f} GeV; N={total_count:,}")
            ax.set_ylabel("Fraction of all MC tags / bin")
            ax.grid(alpha=.2)
    for ax in axes[-1]:
        ax.set_xlabel(r"Reconstructed $M_X^2(e\pi^+\pi^-)$ [GeV$^2$]")
        ax.set_xlim(M2_EDGES[0], M2_EDGES[-1])
    handles, legend_labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, legend_labels, loc="lower center", ncol=3, fontsize=9)
    fig.suptitle("MC truth composition of reconstructed missing-mass shapes")
    fig.tight_layout(rect=(0, .07, 1, .97))
    fig.savefig(output / "8_mc_truth_high_angle_mass.png", dpi=150)
    plt.close(fig)


def write_cell_counts(output, results):
    centers = (M2_EDGES[:-1] + M2_EDGES[1:]) / 2
    ranges = ((.43, .73), (.73, 1.03), (1.03, 1.33))
    with (output / "10_slice_counts.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("mode", "stage", "theta_low_deg", "theta_high_deg",
                         "p_low_GeV", "p_high_GeV", "tags_in_saved_mass_range",
                         "low_shoulder_M2_0.43_0.73", "peak_M2_0.73_1.03",
                         "high_side_M2_1.03_1.33"))
        for mode in MODES:
            for stage in STAGES:
                for j in range(len(THETA_EDGES) - 1):
                    for k in range(len(P_EDGES) - 1):
                        values = results[mode][stage][:, j, k]
                        writer.writerow((mode, stage, THETA_EDGES[j], THETA_EDGES[j + 1],
                                         P_EDGES[k], P_EDGES[k + 1], int(values.sum()),
                                         *(int(values[(centers >= low) & (centers < high)].sum())
                                           for low, high in ranges)))
    with (output / "11_mc_truth_slice_counts.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("truth_category", "theta_low_deg", "theta_high_deg",
                         "p_low_GeV", "p_high_GeV", "tags_in_saved_mass_range",
                         "reco_low_shoulder_M2_0.43_0.73", "reco_peak_M2_0.73_1.03",
                         "reco_high_side_M2_1.03_1.33"))
        for category, histogram in results["mc"]["truth"].items():
            for j in range(len(THETA_EDGES) - 1):
                for k in range(len(P_EDGES) - 1):
                    values = histogram[:, j, k]
                    writer.writerow((category, THETA_EDGES[j], THETA_EDGES[j + 1],
                                     P_EDGES[k], P_EDGES[k + 1], int(values.sum()),
                                     *(int(values[(centers >= low) & (centers < high)].sum())
                                       for low, high in ranges)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path)
    parser.add_argument("--from-cache", type=Path,
                        help="Replot a saved 9_mass_shape_histograms.npz without ROOT input")
    parser.add_argument("--output-root", type=Path, default=Path("output/proton_efficiency"))
    parser.add_argument("--period", default="fa18_inb")
    parser.add_argument("--workers", type=int, default=10)
    parser.add_argument("--max-tags-per-mode", type=int, default=0,
                        help="0 reads all completed Tags rows")
    parser.add_argument("--match-angle", type=float, default=10.0)
    parser.add_argument("--max-relative-dp", type=float)
    args = parser.parse_args()
    if args.from_cache is None and args.input_root is None:
        parser.error("--input-root is required unless --from-cache is given")
    if not 1 <= args.workers <= 10 or args.max_tags_per_mode < 0:
        parser.error("workers must be 1..10 and max-tags-per-mode must be nonnegative")
    if args.match_angle <= 0 or (args.max_relative_dp is not None and args.max_relative_dp <= 0):
        parser.error("matching thresholds must be positive")
    suffix = (f"limited_{args.max_tags_per_mode}_" if args.max_tags_per_mode else "")
    output = args.output_root / args.period / f"{suffix}mass_shape_v1"
    output.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    if args.from_cache is None:
        results = {mode: scan_mode(args.input_root, args.period, mode, args.workers,
                                   args.max_tags_per_mode, args.match_angle,
                                   args.max_relative_dp) for mode in MODES}
        save_cache(output, results, args.match_angle, args.max_relative_dp)
        match_angle, max_relative_dp = args.match_angle, args.max_relative_dp
    else:
        results, match_angle, cached_dp = load_cache(args.from_cache)
        max_relative_dp = None if math.isnan(cached_dp) else cached_dp
        save_cache(output, results, match_angle, max_relative_dp)
        print(f"[CACHE] loaded {args.from_cache}", flush=True)
    for stage, number in zip(STAGES, (1, 2, 3)):
        plot_high_angle(output, results, stage, number)
    for p_index in range(len(P_EDGES) - 1):
        plot_mass_angle_maps(output, results, p_index, 4 + p_index)
    plot_truth_high_angle(output, results)
    write_cell_counts(output, results)
    summary = {"period": args.period, "input_root": str(args.input_root) if args.input_root else None,
               "from_cache": str(args.from_cache) if args.from_cache else None,
               "mass_edges_GeV2": [float(M2_EDGES[0]), float(M2_EDGES[-1]), len(M2_EDGES) - 1],
               "theta_edges_deg": THETA_EDGES.tolist(), "predicted_p_edges_GeV": P_EDGES.tolist(),
               "match_angle_deg": match_angle, "max_relative_dp": max_relative_dp,
               "cache": cache_path(output).name,
               "counts": {mode: {key: results[mode][key] for key in
                                  ("files", "input", "selected", "matched_count")}
                          for mode in MODES},
               "mc_truth_branches": sorted(name for name in results["mc"]["branches"]
                                           if name.startswith("truth_")),
               "mc_truth_categories": {
                   name: int(values.sum()) for name, values in results["mc"]["truth"].items()},
               "mc_truth_definitions": {
                   "valid": "Reconstructed e, pi+, pi- each matched via MC::RecMatch to correct MC::Particle PID",
                   "low_near_high": "truth_missing_M2 split at proton mass squared ±0.15 GeV²",
                   "aligned": "truth-near with a generated proton within match_angle_deg of truth missing direction",
                   "invalid_or_nonfinite": "No complete matched charged truth tag or nonfinite truth missing mass"},
               "notes": "No signal model or efficiency correction is fitted. "
                        "not_matched is the complement of matched within the visible tags."}
    (output / "12_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"[OUTPUT] {output} elapsed={time.perf_counter() - started:.1f}s", flush=True)


if __name__ == "__main__":
    main()
