#!/usr/bin/env python3
"""
plot_external_isr_diagnostics_v2.py

Produce an analysis-note-oriented diagnostic package for the additional
external-ISR transformation applied to RGC exclusive e n pi+ data.

The script compares, event by event,

    *_ISR_mom_corrections.root

with

    *_ISR_externalISR_mom_corrections.root

using the `external_radiation_entry` branch written by apply_external_isr.py.
By default it discovers the three NH3 period pairs (Su22, Fa22, Sp23) in the
paper_versions directory. Other targets can be selected with --target.

Outputs
-------
  01_beam_energy_evolution.png
      Nominal beam energy, internal-ISR effective energy, and
      internal+external-ISR effective energy.

  02_external_photon_energy.png
      External photon-energy distribution, with useful tail fractions.

  03_absolute_kinematic_migrations.png
      Delta Q2, Delta W, Delta xB, Delta t, Delta Mx2, and Delta DepA.

  04_fractional_kinematic_migrations.png
      Fractional changes in Q2, W, xB, |t|, Mx2 (scaled), and DepA.

  05_kinematic_correlations.png
      Before-versus-after correlations for Q2, xB, W, and -t.

  06_period_comparison.png
      Period-by-period external photon energy and fractional migration checks.

  07_survival_summary.png
      Input/output survival, endpoint-cap counts, and analysis-cut migration.

  external_isr_summary.json
      Machine-readable statistics suitable for analysis-note prose.

  external_isr_summary.csv
      Compact period-by-period table.

Typical usage
-------------
  python3 plot_external_isr_diagnostics_v2.py

Explicit directory and target:
  python3 plot_external_isr_diagnostics_v2.py \
      --input-dir /work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/paper_versions \
      --target NH3

Representative-period note plots only (while still summarizing all periods):
  python3 plot_external_isr_diagnostics_v2.py --representative-period fa22

Notes
-----
* Angles in the input tree are assumed to be radians, consistently with the
  transformation script. This diagnostic package does not alter them.
* The default "analysis cuts" are the common DIS cuts Q2 > 1 GeV^2, W > 2 GeV,
  and y < 0.8. They are configurable from the command line.
* The event survival fraction is taken from the provenance JSON when present,
  because events rejected as nonphysical do not appear in the output tree.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import matplotlib.pyplot as plt
import numpy as np
import uproot


DEFAULT_INPUT_DIR = Path(
    "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/paper_versions"
)
DEFAULT_TREE = "PhysicsEvents"
DEFAULT_OUTPUT_SUBDIR = "external_isr_note_diagnostics"
PERIOD_ORDER = ("su22", "fa22", "sp23")
PERIOD_LABELS = {"su22": "Su22", "fa22": "Fa22", "sp23": "Sp23"}
PERIOD_BEAM_ENERGIES = {"su22": 10.5473, "fa22": 10.5563, "sp23": 10.5593}


@dataclass(frozen=True)
class FilePair:
    period: str
    target: str
    before: Path
    after: Path
    provenance: Path | None


@dataclass
class PairData:
    pair: FilePair
    before: dict[str, np.ndarray]
    after: dict[str, np.ndarray]
    provenance: dict[str, Any]
    input_entries: int
    output_entries: int
    removed_entries: int
    endpoint_caps: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create external-ISR validation plots and numerical summaries."
    )
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--tree-name", default=DEFAULT_TREE)
    parser.add_argument(
        "--target",
        default="NH3",
        choices=("NH3", "C", "CH2", "He", "ET"),
    )
    parser.add_argument(
        "--periods",
        nargs="+",
        default=list(PERIOD_ORDER),
        choices=PERIOD_ORDER,
    )
    parser.add_argument(
        "--representative-period",
        default="fa22",
        choices=PERIOD_ORDER,
        help="Period used for the detailed note-ready migration panels.",
    )
    parser.add_argument("--q2-min", type=float, default=1.0)
    parser.add_argument("--w-min", type=float, default=2.0)
    parser.add_argument("--y-max", type=float, default=0.8)
    parser.add_argument(
        "--max-events-per-period",
        type=int,
        default=0,
        help=(
            "Optional deterministic plotting cap per period. Zero keeps all "
            "surviving events. Summary counts always use full provenance counts."
        ),
    )
    parser.add_argument("--dpi", type=int, default=180)
    return parser.parse_args()


def period_from_name(path: Path) -> str | None:
    lower = path.name.lower()
    for period in PERIOD_ORDER:
        if f"rgc_{period}_" in lower:
            return period
    return None


def discover_pairs(input_dir: Path, target: str, periods: Iterable[str]) -> list[FilePair]:
    wanted = set(periods)
    pattern = f"rgc_*_inb_{target}_epi+_ISR_externalISR_mom_corrections.root"
    after_files = sorted(input_dir.glob(pattern))
    pairs: list[FilePair] = []

    for after in after_files:
        period = period_from_name(after)
        if period is None or period not in wanted:
            continue
        before_name = after.name.replace(
            "_ISR_externalISR_mom_corrections.root",
            "_ISR_mom_corrections.root",
        )
        before = after.with_name(before_name)
        if not before.exists():
            raise FileNotFoundError(
                f"Missing corresponding internal-ISR file for {after.name}: {before}"
            )
        provenance = after.with_suffix(".provenance.json")
        pairs.append(
            FilePair(
                period=period,
                target=target,
                before=before,
                after=after,
                provenance=provenance if provenance.exists() else None,
            )
        )

    found = {pair.period for pair in pairs}
    missing = sorted(wanted - found)
    if missing:
        raise FileNotFoundError(
            "Could not find external-ISR pairs for: " + ", ".join(missing)
        )
    return sorted(pairs, key=lambda p: PERIOD_ORDER.index(p.period))


def resolve_tree(
    root_file: uproot.ReadOnlyDirectory,
    requested: str,
    path: Path,
    provenance_name: str | None = None,
) -> uproot.behaviors.TTree.TTree:
    """Resolve a TTree robustly, including files written with a cycle in its name."""
    candidates_to_try: list[str] = []
    for name in (requested, provenance_name):
        if not name:
            continue
        candidates_to_try.extend([name, str(name).split(";")[0]])

    seen: set[str] = set()
    for name in candidates_to_try:
        if name in seen:
            continue
        seen.add(name)
        try:
            obj = root_file[name]
        except (KeyError, uproot.KeyInFileError):
            continue
        if hasattr(obj, "arrays") and hasattr(obj, "num_entries"):
            return obj

    tree_like: list[tuple[str, Any]] = []
    for key in root_file.keys():
        try:
            obj = root_file[key]
        except Exception:
            continue
        if hasattr(obj, "arrays") and hasattr(obj, "num_entries"):
            tree_like.append((str(key), obj))

    if len(tree_like) == 1:
        return tree_like[0][1]

    available = [name for name, _ in tree_like]
    raise KeyError(
        f"Tree {requested!r} not found unambiguously in {path}. "
        f"Tree-like objects: {available}"
    )


def require_branches(tree: uproot.behaviors.TTree.TTree, names: Iterable[str], path: Path) -> None:
    available = set(tree.keys())
    missing = [name for name in names if name not in available]
    if missing:
        raise KeyError(f"{path} is missing required branches: {', '.join(missing)}")


def safe_take_source(
    source_tree: uproot.behaviors.TTree.TTree,
    entries: np.ndarray,
    branches: list[str],
) -> dict[str, np.ndarray]:
    """Read only the contiguous source range needed for sorted absolute entries."""
    if entries.size == 0:
        return {name: np.array([], dtype=float) for name in branches}
    if np.any(np.diff(entries) < 0):
        raise RuntimeError("external_radiation_entry is not monotonically increasing.")
    start = int(entries[0])
    stop = int(entries[-1]) + 1
    block = source_tree.arrays(branches, entry_start=start, entry_stop=stop, library="np")
    offsets = entries - start
    return {name: np.asarray(block[name])[offsets] for name in branches}


def deterministic_subsample(data: PairData, maximum: int) -> PairData:
    if maximum <= 0 or data.output_entries <= maximum:
        return data
    n = len(next(iter(data.after.values())))
    if n <= maximum:
        return data
    idx = np.linspace(0, n - 1, maximum, dtype=np.int64)
    before = {name: values[idx] for name, values in data.before.items()}
    after = {name: values[idx] for name, values in data.after.items()}
    return PairData(
        pair=data.pair,
        before=before,
        after=after,
        provenance=data.provenance,
        input_entries=data.input_entries,
        output_entries=data.output_entries,
        removed_entries=data.removed_entries,
        endpoint_caps=data.endpoint_caps,
    )


def load_pair(pair: FilePair, tree_name: str, max_events: int) -> PairData:
    before_branches = ["Q2", "W", "x", "y", "t", "Mx2", "DepA"]
    after_branches = before_branches + [
        "external_radiation_entry",
        "Egamma_external",
        "Egamma_total",
        "effective_beam_energy_externalISR",
        "external_radiator_thickness",
    ]

    provenance: dict[str, Any] = {}
    if pair.provenance is not None:
        provenance = json.loads(pair.provenance.read_text(encoding="utf-8"))

    with uproot.open(pair.after) as f_after:
        after_tree = resolve_tree(
            f_after,
            tree_name,
            pair.after,
            provenance_name=provenance.get("output_tree"),
        )
        require_branches(after_tree, after_branches, pair.after)
        after = {
            name: np.asarray(values)
            for name, values in after_tree.arrays(after_branches, library="np").items()
        }
        output_entries = int(after_tree.num_entries)

    entries = np.asarray(after.pop("external_radiation_entry"), dtype=np.int64)

    with uproot.open(pair.before) as f_before:
        before_tree = resolve_tree(
            f_before,
            tree_name,
            pair.before,
            provenance_name=provenance.get("input_tree"),
        )
        require_branches(before_tree, before_branches + ["runnum"], pair.before)
        before = safe_take_source(before_tree, entries, before_branches + ["runnum"])
        input_entries_tree = int(before_tree.num_entries)

    input_entries = int(provenance.get("input_entries", input_entries_tree))
    output_entries_meta = int(provenance.get("output_entries", output_entries))
    if output_entries_meta != output_entries:
        raise RuntimeError(
            f"Output-entry mismatch for {pair.after}: tree={output_entries}, "
            f"provenance={output_entries_meta}"
        )
    removed = int(
        provenance.get(
            "removed_nonphysical_recalculated_events",
            input_entries - output_entries,
        )
    )
    endpoint_caps = int(
        provenance.get("external_energy_endpoint_capped_events", 0)
    )

    nominal = np.full(output_entries, PERIOD_BEAM_ENERGIES[pair.period], dtype=float)
    # Egamma_total = Egamma_internal + Egamma_external.
    internal_egamma = np.asarray(after["Egamma_total"], dtype=float) - np.asarray(
        after["Egamma_external"], dtype=float
    )
    after["beam_nominal"] = nominal
    after["beam_internal"] = nominal - internal_egamma

    loaded = PairData(
        pair=pair,
        before=before,
        after=after,
        provenance=provenance,
        input_entries=input_entries,
        output_entries=output_entries,
        removed_entries=removed,
        endpoint_caps=endpoint_caps,
    )
    return deterministic_subsample(loaded, max_events)


def finite(values: np.ndarray) -> np.ndarray:
    return np.asarray(values, dtype=float)[np.isfinite(values)]


def percentile_limits(values: np.ndarray, low: float = 0.2, high: float = 99.8) -> tuple[float, float]:
    arr = finite(values)
    if arr.size == 0:
        return -1.0, 1.0
    lo, hi = np.percentile(arr, [low, high])
    if not np.isfinite(lo) or not np.isfinite(hi) or lo == hi:
        center = float(np.nanmean(arr)) if arr.size else 0.0
        width = max(abs(center) * 0.1, 1.0)
        return center - width, center + width
    pad = 0.04 * (hi - lo)
    return float(lo - pad), float(hi + pad)


def symmetric_limits(values: np.ndarray, percentile: float = 99.5) -> tuple[float, float]:
    arr = np.abs(finite(values))
    if arr.size == 0:
        return -1.0, 1.0
    edge = float(np.percentile(arr, percentile))
    if edge <= 0.0 or not np.isfinite(edge):
        edge = 1.0
    return -1.05 * edge, 1.05 * edge


def hist_step(ax: plt.Axes, arrays: list[np.ndarray], labels: list[str], bins: np.ndarray | int, **kwargs: Any) -> None:
    for values, label in zip(arrays, labels):
        ax.hist(finite(values), bins=bins, histtype="step", linewidth=1.4, label=label, **kwargs)


def analysis_mask(values: dict[str, np.ndarray], q2_min: float, w_min: float, y_max: float) -> np.ndarray:
    return (
        np.isfinite(values["Q2"])
        & np.isfinite(values["W"])
        & np.isfinite(values["y"])
        & (values["Q2"] > q2_min)
        & (values["W"] > w_min)
        & (values["y"] < y_max)
    )


def delta(data: PairData, name: str) -> np.ndarray:
    return np.asarray(data.after[name], dtype=float) - np.asarray(data.before[name], dtype=float)


def fractional_delta(data: PairData, name: str, absolute_denominator: bool = False) -> np.ndarray:
    before = np.asarray(data.before[name], dtype=float)
    denom = np.abs(before) if absolute_denominator else before
    out = np.full_like(before, np.nan, dtype=float)
    valid = np.isfinite(before) & np.isfinite(data.after[name]) & (np.abs(denom) > 1e-12)
    out[valid] = (np.asarray(data.after[name], dtype=float)[valid] - before[valid]) / denom[valid]
    return out


def combined_arrays(data_list: list[PairData], source: str, name: str) -> np.ndarray:
    return np.concatenate([np.asarray(getattr(item, source)[name]) for item in data_list])


def annotate_stats(ax: plt.Axes, values: np.ndarray, unit: str = "") -> None:
    arr = finite(values)
    if arr.size == 0:
        return
    text = (
        f"Mean = {np.mean(arr):.4g}{unit}\n"
        f"Median = {np.median(arr):.4g}{unit}\n"
        f"RMS = {np.std(arr):.4g}{unit}"
    )
    ax.text(
        0.98,
        0.96,
        text,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=9,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85},
    )


def plot_beam_energy(data_list: list[PairData], output: Path, dpi: int) -> None:
    nominal = combined_arrays(data_list, "after", "beam_nominal")
    internal = combined_arrays(data_list, "after", "beam_internal")
    external = combined_arrays(data_list, "after", "effective_beam_energy_externalISR")

    fig, ax = plt.subplots(figsize=(8.2, 5.6))
    lo = max(0.0, min(np.percentile(internal, 0.1), np.percentile(external, 0.1)) - 0.1)
    hi = max(nominal) + 0.03
    bins = np.linspace(lo, hi, 170)
    hist_step(
        ax,
        [nominal, internal, external],
        ["Nominal beam", "After internal ISR", "After internal + external ISR"],
        bins,
        density=True,
    )
    ax.set_yscale("log")
    ax.set_xlabel("Effective incident-electron energy (GeV)")
    ax.set_ylabel("Normalized event density")
    ax.set_title(f"{data_list[0].pair.target}: beam-energy evolution")
    ax.legend()
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(output, dpi=dpi)
    plt.close(fig)


def plot_external_energy(data_list: list[PairData], output: Path, dpi: int) -> None:
    values = combined_arrays(data_list, "after", "Egamma_external")
    arr = finite(values)
    upper = max(float(np.percentile(arr, 99.8)), 0.05)
    bins = np.linspace(0.0, upper, 150)

    fig, ax = plt.subplots(figsize=(8.2, 5.6))
    ax.hist(arr, bins=bins, histtype="step", linewidth=1.5)
    ax.set_yscale("log")
    ax.set_xlabel(r"Additional external-radiation energy $E_{\gamma}^{\rm ext}$ (GeV)")
    ax.set_ylabel("Events")
    ax.set_title(f"{data_list[0].pair.target}: sampled external ISR")
    annotate_stats(ax, arr, " GeV")
    thresholds = (0.05, 0.10, 0.25, 0.50)
    lines = [f"P(Eγext > {int(t*1000)} MeV) = {np.mean(arr > t):.3%}" for t in thresholds]
    ax.text(
        0.98,
        0.60,
        "\n".join(lines),
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=9,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85},
    )
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(output, dpi=dpi)
    plt.close(fig)


def plot_absolute_migrations(data: PairData, output: Path, dpi: int) -> None:
    specs = [
        ("Q2", r"$\Delta Q^2$ (GeV$^2$)"),
        ("W", r"$\Delta W$ (GeV)"),
        ("x", r"$\Delta x_B$"),
        ("t", r"$\Delta t$ (GeV$^2$)"),
        ("Mx2", r"$\Delta M_X^2$ (GeV$^2$)"),
        ("DepA", r"$\Delta D_A$"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(14.2, 8.2))
    for ax, (name, xlabel) in zip(axes.flat, specs):
        values = delta(data, name)
        lo, hi = symmetric_limits(values)
        ax.hist(finite(values), bins=np.linspace(lo, hi, 120), histtype="step", linewidth=1.3)
        ax.axvline(0.0, linestyle="--", linewidth=1.0)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Events")
        ax.grid(alpha=0.25)
        annotate_stats(ax, values)
    fig.suptitle(
        f"{PERIOD_LABELS[data.pair.period]} {data.pair.target}: event-by-event kinematic migration",
        y=0.995,
    )
    fig.tight_layout()
    fig.savefig(output, dpi=dpi)
    plt.close(fig)


def plot_fractional_migrations(data: PairData, output: Path, dpi: int) -> None:
    values = {
        "Q2": 100.0 * fractional_delta(data, "Q2"),
        "W": 100.0 * fractional_delta(data, "W"),
        "x": 100.0 * fractional_delta(data, "x"),
        "t": 100.0 * fractional_delta(data, "t", absolute_denominator=True),
        # Mx2 can cross zero, so display Delta Mx2 rather than a singular ratio.
        "Mx2": delta(data, "Mx2"),
        "DepA": 100.0 * fractional_delta(data, "DepA"),
    }
    specs = [
        ("Q2", r"$100\,\Delta Q^2/Q^2$ (%)"),
        ("W", r"$100\,\Delta W/W$ (%)"),
        ("x", r"$100\,\Delta x_B/x_B$ (%)"),
        ("t", r"$100\,\Delta t/|t|$ (%)"),
        ("Mx2", r"$\Delta M_X^2$ (GeV$^2$)"),
        ("DepA", r"$100\,\Delta D_A/D_A$ (%)"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(14.2, 8.2))
    for ax, (name, xlabel) in zip(axes.flat, specs):
        arr = values[name]
        lo, hi = symmetric_limits(arr)
        ax.hist(finite(arr), bins=np.linspace(lo, hi, 120), histtype="step", linewidth=1.3)
        ax.axvline(0.0, linestyle="--", linewidth=1.0)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Events")
        ax.grid(alpha=0.25)
        annotate_stats(ax, arr)
    fig.suptitle(
        f"{PERIOD_LABELS[data.pair.period]} {data.pair.target}: relative kinematic migration",
        y=0.995,
    )
    fig.tight_layout()
    fig.savefig(output, dpi=dpi)
    plt.close(fig)


def plot_correlations(data: PairData, output: Path, dpi: int) -> None:
    specs = [
        ("Q2", r"$Q^2$ (GeV$^2$)"),
        ("x", r"$x_B$"),
        ("W", r"$W$ (GeV)"),
        ("t", r"$t$ (GeV$^2$)"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 9.4))
    for ax, (name, label) in zip(axes.flat, specs):
        x = np.asarray(data.before[name], dtype=float)
        y = np.asarray(data.after[name], dtype=float)
        valid = np.isfinite(x) & np.isfinite(y)
        lo, hi = percentile_limits(np.concatenate([x[valid], y[valid]]), 0.2, 99.8)
        h = ax.hist2d(x[valid], y[valid], bins=110, range=((lo, hi), (lo, hi)), norm="log")
        ax.plot([lo, hi], [lo, hi], linestyle="--", linewidth=1.0)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_xlabel(f"Internal ISR: {label}")
        ax.set_ylabel(f"Internal + external ISR: {label}")
        fig.colorbar(h[3], ax=ax, label="Events")
    fig.suptitle(
        f"{PERIOD_LABELS[data.pair.period]} {data.pair.target}: before-versus-after correlations",
        y=0.995,
    )
    fig.tight_layout()
    fig.savefig(output, dpi=dpi)
    plt.close(fig)


def plot_period_comparison(data_list: list[PairData], output: Path, dpi: int) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8))
    labels = [PERIOD_LABELS[d.pair.period] for d in data_list]

    ext = [d.after["Egamma_external"] for d in data_list]
    max_ext = max(np.percentile(finite(v), 99.5) for v in ext)
    hist_step(axes[0], ext, labels, np.linspace(0.0, max(max_ext, 0.05), 120), density=True)
    axes[0].set_yscale("log")
    axes[0].set_xlabel(r"$E_{\gamma}^{\rm ext}$ (GeV)")
    axes[0].set_ylabel("Normalized event density")

    dq2 = [100.0 * fractional_delta(d, "Q2") for d in data_list]
    edge_q2 = max(np.percentile(np.abs(finite(v)), 99.0) for v in dq2)
    hist_step(axes[1], dq2, labels, np.linspace(-edge_q2, edge_q2, 120), density=True)
    axes[1].set_xlabel(r"$100\,\Delta Q^2/Q^2$ (%)")
    axes[1].set_ylabel("Normalized event density")

    dx = [100.0 * fractional_delta(d, "x") for d in data_list]
    edge_x = max(np.percentile(np.abs(finite(v)), 99.0) for v in dx)
    hist_step(axes[2], dx, labels, np.linspace(-edge_x, edge_x, 120), density=True)
    axes[2].set_xlabel(r"$100\,\Delta x_B/x_B$ (%)")
    axes[2].set_ylabel("Normalized event density")

    for ax in axes:
        ax.axvline(0.0, linestyle="--", linewidth=1.0)
        ax.grid(alpha=0.25)
        ax.legend()
    fig.suptitle(f"{data_list[0].pair.target}: period stability of the external-ISR response")
    fig.tight_layout()
    fig.savefig(output, dpi=dpi)
    plt.close(fig)


def cut_migration_counts(data: PairData, args: argparse.Namespace) -> dict[str, int]:
    before_mask = analysis_mask(data.before, args.q2_min, args.w_min, args.y_max)
    after_mask = analysis_mask(data.after, args.q2_min, args.w_min, args.y_max)
    return {
        "before_pass": int(np.count_nonzero(before_mask)),
        "after_pass": int(np.count_nonzero(after_mask)),
        "pass_to_fail": int(np.count_nonzero(before_mask & ~after_mask)),
        "fail_to_pass": int(np.count_nonzero(~before_mask & after_mask)),
        "pass_both": int(np.count_nonzero(before_mask & after_mask)),
    }


def plot_survival(data_list: list[PairData], output: Path, args: argparse.Namespace) -> None:
    labels = [PERIOD_LABELS[d.pair.period] for d in data_list]
    input_counts = np.array([d.input_entries for d in data_list], dtype=float)
    output_counts = np.array([d.output_entries for d in data_list], dtype=float)
    survival = 100.0 * output_counts / input_counts
    cut_counts = [cut_migration_counts(d, args) for d in data_list]
    before_cut = np.array([item["before_pass"] for item in cut_counts], dtype=float)
    after_cut = np.array([item["after_pass"] for item in cut_counts], dtype=float)
    cut_ratio = np.divide(100.0 * after_cut, before_cut, out=np.full_like(after_cut, np.nan), where=before_cut > 0)

    x = np.arange(len(data_list))
    width = 0.36
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.0))

    axes[0].bar(x - width / 2, input_counts, width, label="Input ISR entries")
    axes[0].bar(x + width / 2, output_counts, width, label="Surviving external-ISR entries")
    axes[0].set_xticks(x, labels)
    axes[0].set_ylabel("Events")
    axes[0].set_title("Tree-level survival")
    axes[0].legend()
    axes[0].grid(axis="y", alpha=0.25)
    for i, value in enumerate(survival):
        axes[0].text(i, output_counts[i], f"{value:.3f}%", ha="center", va="bottom", fontsize=9)

    axes[1].bar(x - width / 2, before_cut, width, label="Pass before")
    axes[1].bar(x + width / 2, after_cut, width, label="Pass after")
    axes[1].set_xticks(x, labels)
    axes[1].set_ylabel("Events")
    axes[1].set_title(
        rf"DIS-cut migration: $Q^2>{args.q2_min:g}$ GeV$^2$, $W>{args.w_min:g}$ GeV, $y<{args.y_max:g}$"
    )
    axes[1].legend()
    axes[1].grid(axis="y", alpha=0.25)
    for i, value in enumerate(cut_ratio):
        axes[1].text(i, after_cut[i], f"{value:.3f}%", ha="center", va="bottom", fontsize=9)

    fig.suptitle(f"{data_list[0].pair.target}: external-ISR event survival and cut migration")
    fig.tight_layout()
    fig.savefig(output, dpi=args.dpi)
    plt.close(fig)


def summarize_array(values: np.ndarray) -> dict[str, float]:
    arr = finite(values)
    if arr.size == 0:
        return {key: float("nan") for key in ("mean", "std", "median", "p05", "p95", "p99", "min", "max")}
    return {
        "mean": float(np.mean(arr)),
        "std": float(np.std(arr)),
        "median": float(np.median(arr)),
        "p05": float(np.percentile(arr, 5)),
        "p95": float(np.percentile(arr, 95)),
        "p99": float(np.percentile(arr, 99)),
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
    }


def period_summary(data: PairData, args: argparse.Namespace) -> dict[str, Any]:
    ext = np.asarray(data.after["Egamma_external"], dtype=float)
    d_q2 = 100.0 * fractional_delta(data, "Q2")
    d_x = 100.0 * fractional_delta(data, "x")
    d_w = 100.0 * fractional_delta(data, "W")
    d_t = 100.0 * fractional_delta(data, "t", absolute_denominator=True)
    d_mx2 = delta(data, "Mx2")
    d_depa = 100.0 * fractional_delta(data, "DepA")
    cuts = cut_migration_counts(data, args)

    return {
        "period": data.pair.period,
        "period_label": PERIOD_LABELS[data.pair.period],
        "target": data.pair.target,
        "input_file": str(data.pair.before),
        "output_file": str(data.pair.after),
        "input_entries": data.input_entries,
        "output_entries": data.output_entries,
        "removed_nonphysical_entries": data.removed_entries,
        "survival_fraction": data.output_entries / data.input_entries,
        "endpoint_capped_entries": data.endpoint_caps,
        "external_photon_energy_GeV": summarize_array(ext),
        "external_photon_tail_fractions": {
            "above_50_MeV": float(np.mean(ext > 0.050)),
            "above_100_MeV": float(np.mean(ext > 0.100)),
            "above_250_MeV": float(np.mean(ext > 0.250)),
            "above_500_MeV": float(np.mean(ext > 0.500)),
        },
        "fractional_migrations_percent": {
            "Q2": summarize_array(d_q2),
            "W": summarize_array(d_w),
            "xB": summarize_array(d_x),
            "t_over_abs_t": summarize_array(d_t),
            "DepA": summarize_array(d_depa),
        },
        "absolute_migrations": {
            "Mx2_GeV2": summarize_array(d_mx2),
        },
        "large_fractional_migration_fractions": {
            "abs_delta_Q2_over_Q2_above_1pct": float(np.mean(np.abs(d_q2) > 1.0)),
            "abs_delta_Q2_over_Q2_above_5pct": float(np.mean(np.abs(d_q2) > 5.0)),
            "abs_delta_Q2_over_Q2_above_10pct": float(np.mean(np.abs(d_q2) > 10.0)),
            "abs_delta_xB_over_xB_above_1pct": float(np.mean(np.abs(d_x) > 1.0)),
            "abs_delta_xB_over_xB_above_5pct": float(np.mean(np.abs(d_x) > 5.0)),
            "abs_delta_xB_over_xB_above_10pct": float(np.mean(np.abs(d_x) > 10.0)),
        },
        "analysis_cut_counts": cuts,
        "analysis_cut_after_over_before": (
            cuts["after_pass"] / cuts["before_pass"] if cuts["before_pass"] else None
        ),
    }


def write_summary(data_list: list[PairData], output_dir: Path, args: argparse.Namespace) -> None:
    summaries = [period_summary(data, args) for data in data_list]

    all_ext = combined_arrays(data_list, "after", "Egamma_external")
    combined = {
        "target": args.target,
        "periods": list(args.periods),
        "analysis_cuts": {
            "Q2_min_GeV2": args.q2_min,
            "W_min_GeV": args.w_min,
            "y_max": args.y_max,
        },
        "total_input_entries": int(sum(d.input_entries for d in data_list)),
        "total_output_entries": int(sum(d.output_entries for d in data_list)),
        "total_removed_nonphysical_entries": int(sum(d.removed_entries for d in data_list)),
        "total_endpoint_capped_entries": int(sum(d.endpoint_caps for d in data_list)),
        "overall_survival_fraction": (
            sum(d.output_entries for d in data_list) / sum(d.input_entries for d in data_list)
        ),
        "external_photon_energy_GeV": summarize_array(all_ext),
        "external_photon_tail_fractions": {
            "above_50_MeV": float(np.mean(all_ext > 0.050)),
            "above_100_MeV": float(np.mean(all_ext > 0.100)),
            "above_250_MeV": float(np.mean(all_ext > 0.250)),
            "above_500_MeV": float(np.mean(all_ext > 0.500)),
        },
    }

    payload = {
        "schema_version": 1,
        "script": Path(__file__).name,
        "combined": combined,
        "period_summaries": summaries,
    }
    (output_dir / "external_isr_summary.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )

    csv_path = output_dir / "external_isr_summary.csv"
    fieldnames = [
        "period",
        "target",
        "input_entries",
        "output_entries",
        "removed_entries",
        "survival_percent",
        "endpoint_caps",
        "mean_Egamma_external_MeV",
        "median_Egamma_external_MeV",
        "fraction_above_100_MeV_percent",
        "mean_delta_Q2_over_Q2_percent",
        "rms_delta_Q2_over_Q2_percent",
        "mean_delta_xB_over_xB_percent",
        "rms_delta_xB_over_xB_percent",
        "before_DIS_cut_entries",
        "after_DIS_cut_entries",
        "pass_to_fail",
        "fail_to_pass",
    ]
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for item in summaries:
            writer.writerow({
                "period": item["period_label"],
                "target": item["target"],
                "input_entries": item["input_entries"],
                "output_entries": item["output_entries"],
                "removed_entries": item["removed_nonphysical_entries"],
                "survival_percent": 100.0 * item["survival_fraction"],
                "endpoint_caps": item["endpoint_capped_entries"],
                "mean_Egamma_external_MeV": 1000.0 * item["external_photon_energy_GeV"]["mean"],
                "median_Egamma_external_MeV": 1000.0 * item["external_photon_energy_GeV"]["median"],
                "fraction_above_100_MeV_percent": 100.0 * item["external_photon_tail_fractions"]["above_100_MeV"],
                "mean_delta_Q2_over_Q2_percent": item["fractional_migrations_percent"]["Q2"]["mean"],
                "rms_delta_Q2_over_Q2_percent": item["fractional_migrations_percent"]["Q2"]["std"],
                "mean_delta_xB_over_xB_percent": item["fractional_migrations_percent"]["xB"]["mean"],
                "rms_delta_xB_over_xB_percent": item["fractional_migrations_percent"]["xB"]["std"],
                "before_DIS_cut_entries": item["analysis_cut_counts"]["before_pass"],
                "after_DIS_cut_entries": item["analysis_cut_counts"]["after_pass"],
                "pass_to_fail": item["analysis_cut_counts"]["pass_to_fail"],
                "fail_to_pass": item["analysis_cut_counts"]["fail_to_pass"],
            })


def main() -> int:
    args = parse_args()
    input_dir = args.input_dir.expanduser().resolve()
    if not input_dir.is_dir():
        raise NotADirectoryError(f"Input directory does not exist: {input_dir}")

    output_dir = (
        args.output_dir.expanduser().resolve()
        if args.output_dir
        else input_dir / DEFAULT_OUTPUT_SUBDIR / args.target
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    pairs = discover_pairs(input_dir, args.target, args.periods)
    data_list = [
        load_pair(pair, args.tree_name, args.max_events_per_period)
        for pair in pairs
    ]
    by_period = {item.pair.period: item for item in data_list}
    if args.representative_period not in by_period:
        raise RuntimeError(
            f"Representative period {args.representative_period} was not loaded."
        )
    representative = by_period[args.representative_period]

    plot_beam_energy(data_list, output_dir / "01_beam_energy_evolution.png", args.dpi)
    plot_external_energy(data_list, output_dir / "02_external_photon_energy.png", args.dpi)
    plot_absolute_migrations(
        representative,
        output_dir / "03_absolute_kinematic_migrations.png",
        args.dpi,
    )
    plot_fractional_migrations(
        representative,
        output_dir / "04_fractional_kinematic_migrations.png",
        args.dpi,
    )
    plot_correlations(
        representative,
        output_dir / "05_kinematic_correlations.png",
        args.dpi,
    )
    plot_period_comparison(data_list, output_dir / "06_period_comparison.png", args.dpi)
    plot_survival(data_list, output_dir / "07_survival_summary.png", args)
    write_summary(data_list, output_dir, args)

    print("\nExternal-ISR diagnostic package complete")
    print("=" * 72)
    print(f"Target:                 {args.target}")
    print(f"Periods:                {', '.join(args.periods)}")
    print(f"Representative period: {args.representative_period}")
    print(f"Output directory:       {output_dir}")
    print("Created:")
    for path in sorted(output_dir.iterdir()):
        print(f"  {path.name}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"\nFATAL ERROR:\n{exc}", file=sys.stderr)
        raise SystemExit(1)
