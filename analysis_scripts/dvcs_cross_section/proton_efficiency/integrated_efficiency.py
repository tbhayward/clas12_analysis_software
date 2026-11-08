#!/usr/bin/env python3
"""First-pass, fully integrated e pi+ pi- tag-and-probe study.

Run from any directory. Input ROOT parts must contain the campaign's Tags tree.
The window and matching thresholds are exploratory, not calibrated selections.
"""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot

FIELDS = ("missing_M2", "missing_theta", "missing_p", "best_selected_rec_index",
          "best_selected_match_angle", "best_selected_relative_delta_p")
EDGES = {
    "missing_M2": np.linspace(-1.0, 3.0, 161),
    "missing_theta": np.linspace(0.0, 130.0, 131),
    "missing_p": np.linspace(0.0, 1.5, 121),
    "best_selected_match_angle": np.linspace(0.0, 30.0, 121),
    "best_selected_relative_delta_p": np.linspace(-1.0, 1.0, 121),
}
LABELS = {
    "missing_M2": r"$M_X^2(e\pi^+\pi^-)$ [GeV$^2$]",
    "missing_theta": r"Predicted proton $\theta$ [degrees]",
    "missing_p": r"Predicted proton momentum [GeV]",
    "best_selected_match_angle": "Best selected proton opening angle [degrees]",
    "best_selected_relative_delta_p": r"$(p_{rec}-p_{miss})/p_{miss}$",
}


def empty_result():
    return {"all": 0, "window": 0, "selected": 0, "matched": 0,
            "hist": {name: {stage: np.zeros(len(edges) - 1, dtype=np.int64)
                            for stage in ("all", "window", "matched")}
                     for name, edges in EDGES.items()}}


def process_file(task):
    filename, lo, hi, max_angle, max_rel_dp, chunk_size, entry_stop = task
    result = empty_result()
    with uproot.open(filename) as source:
        tree = source["Tags"]
        missing = set(FIELDS) - set(tree.keys())
        if missing:
            raise RuntimeError(f"{filename}: missing Tags branches {sorted(missing)}")
        for arrays in tree.iterate(FIELDS, step_size=chunk_size, entry_stop=entry_stop,
                                   library="np"):
            mass = arrays["missing_M2"]
            window = np.isfinite(mass) & (mass >= lo) & (mass <= hi)
            selected = window & (arrays["best_selected_rec_index"] >= 0)
            angle = arrays["best_selected_match_angle"]
            rel_dp = arrays["best_selected_relative_delta_p"]
            matched = selected & np.isfinite(angle) & (angle < max_angle)
            if max_rel_dp is not None:
                matched &= np.isfinite(rel_dp) & (np.abs(rel_dp) < max_rel_dp)
            result["all"] += len(mass)
            result["window"] += int(window.sum())
            result["selected"] += int(selected.sum())
            result["matched"] += int(matched.sum())
            for name, edges in EDGES.items():
                values = arrays[name]
                finite = np.isfinite(values)
                for stage, mask in (("all", finite), ("window", finite & window),
                                    ("matched", finite & matched)):
                    result["hist"][name][stage] += np.histogram(values[mask], edges)[0]
    return result


def add_result(total, part):
    for key in ("all", "window", "selected", "matched"):
        total[key] += part[key]
    for name in EDGES:
        for stage in ("all", "window", "matched"):
            total["hist"][name][stage] += part["hist"][name][stage]


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
    error = np.sqrt(np.divide(value * (1 - value), denominator,
                              out=np.full_like(value, np.nan), where=denominator > 0))
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
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


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
    parser.add_argument("--chunk-size", type=int, default=100000)
    parser.add_argument("--fast-mode", action="store_true",
                        help="Preview the first 10000 tags from one completed file per mode")
    args = parser.parse_args()
    if not 1 <= args.workers <= 10 or args.chunk_size < 1:
        parser.error("workers must be 1..10 and chunk-size must be positive")
    if args.mass_half_width <= 0 or args.match_angle <= 0:
        parser.error("mass-half-width and match-angle must be positive")
    if args.max_relative_dp is not None and args.max_relative_dp <= 0:
        parser.error("max-relative-dp must be positive")
    center = 0.9382720813 ** 2
    lo, hi = center - args.mass_half_width, center + args.mass_half_width
    output = args.output_root / args.period / ("fast_preview" if args.fast_mode else "integrated")
    output.mkdir(parents=True, exist_ok=True)
    results = {}
    files_by_mode = {}
    for mode in ("data", "mc"):
        directory = args.input_root / args.period / mode / "parts"
        files = sorted(path for path in directory.glob("*.root")
                       if Path(str(path)[:-5] + ".ok.json").is_file())
        if not files:
            raise SystemExit(f"No completed ROOT parts for {args.period}/{mode} in {directory}")
        if args.fast_mode:
            chosen = next((path for path in files if json.loads(
                Path(str(path)[:-5] + ".ok.json").read_text()).get("selected_tags", 0) >= 10000), None)
            if chosen is None:
                raise SystemExit(f"No completed {args.period}/{mode} part has 10000 tags")
            files = [chosen]
            print(f"[{mode}] FAST PREVIEW: reading first 10000 tags from {files[0].name}",
                  flush=True)
        files_by_mode[mode] = len(files)
        tasks = [(str(path), lo, hi, args.match_angle, args.max_relative_dp,
                  args.chunk_size, 10000 if args.fast_mode else None)
                 for path in files]
        total = empty_result()
        with ProcessPoolExecutor(max_workers=min(args.workers, len(tasks))) as pool:
            futures = {pool.submit(process_file, task): task[0] for task in tasks}
            for index, future in enumerate(as_completed(futures), 1):
                try:
                    add_result(total, future.result())
                except Exception as error:
                    raise RuntimeError(f"Failed reading {futures[future]}: {error}") from error
                if args.fast_mode or index % 10 == 0 or index == len(tasks):
                    print(f"[{mode}] {index}/{len(tasks)} files ({100 * index / len(tasks):.1f}%)", flush=True)
        results[mode] = total
        print(f"[{mode}] tags={total['all']} window={total['window']} "
              f"selected={total['selected']} matched={total['matched']}", flush=True)
        for name in EDGES:
            plot_counts(output / f"{mode}_{name}.png", name,
                        [("All tags", total["hist"][name]["all"]),
                         ("Mass window", total["hist"][name]["window"]),
                         ("Matched", total["hist"][name]["matched"])],
                        logy=True, title=f"{args.period} {mode}: {name}")
    for name in ("missing_theta", "missing_p"):
        plot_ratio(output / f"data_vs_mc_{name}_rate.png", name,
                   results["data"], results["mc"], f"{args.period}: provisional tag-and-probe rate")
    summary = {"period": args.period, "input_root": str(args.input_root),
               "fast_mode": args.fast_mode,
               "files": files_by_mode, "mass_window_GeV2": [lo, hi],
               "match_angle_deg": args.match_angle, "max_relative_dp": args.max_relative_dp,
               "counts": {mode: {key: result[key] for key in ("all", "window", "selected", "matched")}
                          for mode, result in results.items()}}
    for mode in ("data", "mc"):
        n, d = results[mode]["matched"], results[mode]["window"]
        summary["counts"][mode]["provisional_rate"] = n / d if d else None
    a = summary["counts"]["data"]["provisional_rate"]
    b = summary["counts"]["mc"]["provisional_rate"]
    summary["provisional_data_mc_ratio"] = a / b if a is not None and b else None
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"[OUTPUT] {output}", flush=True)


if __name__ == "__main__":
    main()
