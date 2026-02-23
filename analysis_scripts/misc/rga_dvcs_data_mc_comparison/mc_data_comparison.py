#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
calibration_vertex_vz_compare.py

Compare DATA vs DVCSGEN MC for particle vertex z (particle_vz), for:
  - electron (pid = 11)
  - proton   (pid = 2212)

Photon plots are removed.

Canvas layout (2x3):
  Top row:    Sp18 Inb | Sp18 Out | (empty)
  Bottom row: Fa18 Inb | Fa18 Out | Sp19 Inb

DATA histogram: black
MC histogram:   red
Y scale: log

Cuts (THIS VERSION):
  - We compute the Highest Density Interval (HDI): the shortest interval containing
    a user-specified fraction of events.
  - The coverage fraction is provided at runtime, e.g.:
      ./calibration_vertex_vz_compare.py --coverage 0.99

Parallel processing:
  - At most 5 workers.

X-axis / histogram range:
  - Fixed to [-20, 20] (cm).

Output:
  output/vz_electron.png
  output/vz_proton.png

Dependencies:
  pip install uproot numpy matplotlib
"""

import os
import sys
import argparse
import traceback
from dataclasses import dataclass
from typing import Dict, Tuple

import numpy as np
import uproot
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor, as_completed


# -----------------------------------------------------------------------------
# USER SETTINGS (STATIC)
# -----------------------------------------------------------------------------

TREE_NAME = "PhysicsEvents"

# Histogram range and binning
VZ_MIN = -20.0
VZ_MAX = 20.0
N_BINS = 200

# Hard limit on workers
MAX_WORKERS = 5

OUTDIR = "output"

DATA_FILES = {
    "Sp18 Inb": "/volatile/clas12/thayward/dvcs_temp/calibration_data_rga_sp18_inb.root",
    "Sp18 Out": "/volatile/clas12/thayward/dvcs_temp/calibration_data_rga_sp18_out.root",
    "Fa18 Inb": "/volatile/clas12/thayward/dvcs_temp/calibration_data_rga_fa18_inb.root",
    "Fa18 Out": "/volatile/clas12/thayward/dvcs_temp/calibration_data_rga_fa18_out.root",
    "Sp19 Inb": "/volatile/clas12/thayward/dvcs_temp/calibration_data_rga_sp19_inb.root",
}

MC_FILES = {
    "Sp18 Inb": "/volatile/clas12/thayward/dvcs_temp/calibration_dvcsgen_rga_sp18_inb.root",
    "Sp18 Out": "/volatile/clas12/thayward/dvcs_temp/calibration_dvcsgen_rga_sp18_out.root",
    "Fa18 Inb": "/volatile/clas12/thayward/dvcs_temp/calibration_dvcsgen_rga_fa18_inb.root",
    "Fa18 Out": "/volatile/clas12/thayward/dvcs_temp/calibration_dvcsgen_rga_fa18_out.root",
    "Sp19 Inb": "/volatile/clas12/thayward/dvcs_temp/calibration_dvcsgen_rga_sp19_inb.root",
}

PANEL_POS = {
    "Sp18 Inb": (0, 0),
    "Sp18 Out": (0, 1),
    "Fa18 Inb": (1, 0),
    "Fa18 Out": (1, 1),
    "Sp19 Inb": (1, 2),
}

# -----------------------------------------------------------------------------
# LEGACY PASS2 VERTEX CUTS (COMMENTED OUT; kept for easy revert)
# -----------------------------------------------------------------------------
# VTX_CUTS_POS = {
#     "Sp18 Inb": (-7.8790, 1.5150),
#     "Sp18 Out": (-6.6667, 2.7273),
#     "Fa18 Inb": (-8.4850, 0.6060),
#     "Fa18 Out": (-6.9700, 1.8180),
#     "Sp19 Inb": (-8.4850, 0.6060),
# }
#
# VTX_CUTS_NEG = {
#     "Sp18 Inb": (-6.0606, 1.8182),
#     "Sp18 Out": (-7.2730, 0.9091),
#     "Fa18 Inb": (-6.3640, 1.5150),
#     "Fa18 Out": (-7.8790, 0.3030),
#     "Sp19 Inb": (-6.3640, 1.5150),
# }
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# INTERNALS
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class Hist1D:
    counts: np.ndarray          # normalized to unit integral over [VZ_MIN, VZ_MAX]
    edges: np.ndarray
    centers: np.ndarray
    n_selected: int

    peak_index: int
    mode_vz: float
    mode_count: int
    mode_frac: float

    cut_low: float
    cut_high: float
    frac_in_cut: float
    n_in_cut: int
    hdi_width: float


def fatal(msg: str) -> None:
    sys.stderr.write(f"\nFATAL: {msg}\n")
    sys.exit(1)
    #endif


def ensure_outdir(path: str) -> None:
    if not os.path.isdir(path):
        os.makedirs(path, exist_ok=True)
    #endif


def validate_inputs(coverage: float) -> None:
    if VZ_MAX <= VZ_MIN:
        fatal("Invalid VZ range.")
    #endif
    if N_BINS <= 0:
        fatal("Invalid N_BINS.")
    #endif
    if MAX_WORKERS <= 0:
        fatal("Invalid MAX_WORKERS.")
    #endif
    if coverage <= 0.0 or coverage >= 1.0:
        fatal(f"--coverage must be in (0,1). Got {coverage}")
    #endif

    for label, path in DATA_FILES.items():
        if not os.path.isfile(path):
            fatal(f"Missing DATA file for '{label}': {path}")
        #endif
    #endfor
    for label, path in MC_FILES.items():
        if not os.path.isfile(path):
            fatal(f"Missing MC file for '{label}': {path}")
        #endif
    #endfor

    if set(DATA_FILES.keys()) != set(MC_FILES.keys()):
        fatal("DATA_FILES and MC_FILES keys differ.")
    #endif

    for k in DATA_FILES.keys():
        if k not in PANEL_POS:
            fatal(f"PANEL_POS missing period '{k}'")
        #endif
    #endfor


def normalize_to_integral(counts: np.ndarray) -> np.ndarray:
    s = float(np.sum(counts))
    if s <= 0.0:
        return np.zeros_like(counts, dtype=np.float64)
    #endif
    return counts.astype(np.float64) / s


def compute_mode_from_raw(counts_raw: np.ndarray, centers: np.ndarray, n_selected: int) -> Tuple[int, float, int, float]:
    if counts_raw.size == 0:
        return 0, 0.0, 0, 0.0
    #endif
    peak_index = int(np.argmax(counts_raw))
    mode_count = int(counts_raw[peak_index])
    mode_vz = float(centers[peak_index])
    mode_frac = float(mode_count) / float(n_selected) if n_selected > 0 else 0.0
    return peak_index, mode_vz, mode_count, mode_frac


def hdi_interval(sorted_vals: np.ndarray, coverage: float) -> Tuple[float, float]:
    """
    Highest Density Interval (HDI): shortest interval containing `coverage` fraction of points.

    Input must be sorted.
    Returns (low, high).
    """
    n = int(sorted_vals.size)
    if n <= 0:
        fatal("HDI: empty array.")
    #endif

    m = int(np.floor(coverage * n))
    if m < 1:
        fatal("HDI: coverage too small for sample size.")
    #endif
    if m >= n:
        return float(sorted_vals[0]), float(sorted_vals[-1])
    #endif

    widths = sorted_vals[m - 1 :] - sorted_vals[: n - (m - 1)]
    i_best = int(np.argmin(widths))
    low = float(sorted_vals[i_best])
    high = float(sorted_vals[i_best + m - 1])

    if not (low < high):
        fatal(f"HDI produced non-ordered interval: low={low}, high={high}")
    #endif

    return low, high


def compute_hist_for_file(args: Tuple[str, str, int, str, float]) -> Tuple[str, str, int, Hist1D]:
    period_label, root_path, pid, tree_name, coverage = args

    try:
        with uproot.open(root_path) as f:
            if tree_name not in f:
                raise KeyError(f"TTree '{tree_name}' not found. Keys: {list(f.keys())}")
            #endif
            tree = f[tree_name]

            required = ["particle_pid", "particle_vz"]
            for br in required:
                if br not in tree.keys():
                    raise KeyError(f"Branch '{br}' not found in tree '{tree_name}'.")
                #endif
            #endfor

            arr = tree.arrays(required, library="np")
            pids = arr["particle_pid"]
            vz = arr["particle_vz"]

            mask = (pids == pid) & np.isfinite(vz)
            vz_sel = vz[mask].astype(np.float64)

            counts_raw, edges = np.histogram(vz_sel, bins=N_BINS, range=(VZ_MIN, VZ_MAX))
            counts_norm = normalize_to_integral(counts_raw)
            centers = 0.5 * (edges[:-1] + edges[1:])

            peak_index, mode_vz, mode_count, mode_frac = compute_mode_from_raw(
                counts_raw=counts_raw, centers=centers, n_selected=int(vz_sel.size)
            )

            if vz_sel.size <= 0:
                cut_low = float("nan")
                cut_high = float("nan")
                n_in_cut = 0
                frac_in_cut = 0.0
                hdi_width = float("nan")
            else:
                vz_sorted = np.sort(vz_sel)
                cut_low, cut_high = hdi_interval(vz_sorted, coverage)
                hdi_width = float(cut_high - cut_low)

                in_cut = (vz_sel >= cut_low) & (vz_sel <= cut_high)
                n_in_cut = int(np.count_nonzero(in_cut))
                frac_in_cut = float(n_in_cut) / float(vz_sel.size)
            #endif

            h = Hist1D(
                counts=counts_norm.astype(np.float64),
                edges=edges.astype(np.float64),
                centers=centers.astype(np.float64),
                n_selected=int(vz_sel.size),
                peak_index=int(peak_index),
                mode_vz=float(mode_vz),
                mode_count=int(mode_count),
                mode_frac=float(mode_frac),
                cut_low=float(cut_low),
                cut_high=float(cut_high),
                frac_in_cut=float(frac_in_cut),
                n_in_cut=int(n_in_cut),
                hdi_width=float(hdi_width),
            )
            return (period_label, root_path, pid, h)

    except Exception as e:
        tb = traceback.format_exc()
        raise RuntimeError(
            f"Failed processing file:\n"
            f"  period={period_label}\n"
            f"  path={root_path}\n"
            f"  pid={pid}\n"
            f"  tree={tree_name}\n"
            f"  error={repr(e)}\n\n{tb}"
        )


def run_parallel_hists(file_map: Dict[str, str], pid: int, coverage: float) -> Dict[str, Hist1D]:
    items = sorted(file_map.items(), key=lambda kv: kv[0])
    tasks = [(label, path, pid, TREE_NAME, coverage) for (label, path) in items]

    n_workers = min(MAX_WORKERS, len(tasks))
    if n_workers < 1:
        fatal("No tasks.")
    #endif

    out: Dict[str, Hist1D] = {}
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        futures = [ex.submit(compute_hist_for_file, t) for t in tasks]
        for fut in as_completed(futures):
            period_label, root_path, pid_ret, hist = fut.result()
            out[period_label] = hist
        #endfor
    #endwith

    for period_label in file_map.keys():
        if period_label not in out:
            fatal(f"Missing result for {period_label}")
        #endif
    #endfor

    return out


def plot_2x3_canvas(
    title: str,
    data_hists: Dict[str, Hist1D],
    mc_hists: Dict[str, Hist1D],
    outpath: str,
    coverage: float,
) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(16, 9), sharex=True, sharey=False)
    fig.suptitle(title, fontsize=18)

    for r in range(2):
        for c in range(3):
            axes[r, c].axis("off")
        #endfor
    #endfor

    for period_label, (r, c) in PANEL_POS.items():
        ax = axes[r, c]
        ax.axis("on")

        dh = data_hists[period_label]
        mh = mc_hists[period_label]

        data_label = f"data: {100.0*dh.frac_in_cut:.2f}% in-cut (N={dh.n_selected})"
        mc_label = f"mc: {100.0*mh.frac_in_cut:.2f}% in-cut (N={mh.n_selected})"

        ax.step(dh.centers, dh.counts, where="mid", color="black", linewidth=1.2, label=data_label)
        ax.step(mh.centers, mh.counts, where="mid", color="red", linewidth=1.2, label=mc_label)

        ax.axvline(dh.cut_low, color="black", linestyle="--", linewidth=1.0)
        ax.axvline(dh.cut_high, color="black", linestyle="--", linewidth=1.0)
        ax.axvline(mh.cut_low, color="red", linestyle="--", linewidth=1.0)
        ax.axvline(mh.cut_high, color="red", linestyle="--", linewidth=1.0)

        ax.set_title(period_label, fontsize=13)
        ax.set_xlim(VZ_MIN, VZ_MAX)
        ax.set_xlabel("particle_vz (cm)", fontsize=12)
        ax.set_ylabel("normalized counts", fontsize=12)
        ax.grid(True, alpha=0.25)
        ax.legend(loc="upper right", fontsize=9, frameon=True)

        ax.set_yscale("log")
        positive_vals = np.concatenate([dh.counts[dh.counts > 0.0], mh.counts[mh.counts > 0.0]])
        y_min = float(np.min(positive_vals)) if positive_vals.size > 0 else 1e-12
        ax.set_ylim(bottom=y_min)

        ax.text(
            0.98,
            0.06,
            f"HDI {100.0*coverage:.1f}%:\n"
            f"data: ({dh.cut_low:.3f}, {dh.cut_high:.3f})\n"
            f"mc:   ({mh.cut_low:.3f}, {mh.cut_high:.3f})",
            transform=ax.transAxes,
            fontsize=9,
            verticalalignment="bottom",
            horizontalalignment="right",
        )

    #endfor

    axes[0, 2].axis("off")
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.95])
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def print_summary(pid: int, pid_label: str, data_hists: Dict[str, Hist1D], mc_hists: Dict[str, Hist1D], coverage: float) -> None:
    print("")
    print("------------------------------------------------------------")
    print(f"SUMMARY: {pid_label} (pid={pid})")
    print(f"  Range=({VZ_MIN:.1f},{VZ_MAX:.1f}) (cm), bins={N_BINS}, HDI coverage={coverage:.5f}")
    print("------------------------------------------------------------")
    periods = sorted(PANEL_POS.keys(), key=lambda k: (PANEL_POS[k][0], PANEL_POS[k][1]))
    for period_label in periods:
        d = data_hists[period_label]
        m = mc_hists[period_label]
        print(f"{period_label}:")
        print(f"  data: N={d.n_selected}  in_cut={100.0*d.frac_in_cut:.3f}%  cut=({d.cut_low:.4f},{d.cut_high:.4f})  width={d.hdi_width:.4f}  mode_vz={d.mode_vz:.4f}  mode_frac={100.0*d.mode_frac:.3f}%")
        print(f"  mc:   N={m.n_selected}  in_cut={100.0*m.frac_in_cut:.3f}%  cut=({m.cut_low:.4f},{m.cut_high:.4f})  width={m.hdi_width:.4f}  mode_vz={m.mode_vz:.4f}  mode_frac={100.0*m.mode_frac:.3f}%")
    #endfor


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Vertex z calibration plots with HDI-based central coverage cuts (data vs mc)."
    )
    p.add_argument(
        "--coverage",
        type=float,
        required=True,
        help="Central coverage fraction for HDI window (e.g. 0.99). Must be in (0,1).",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    coverage = float(args.coverage)

    validate_inputs(coverage)
    ensure_outdir(OUTDIR)

    particles = [
        (11, "electron", "vz_electron.png"),
        (2212, "proton", "vz_proton.png"),
    ]

    for pid, pid_label, fname in particles:
        data_hists = run_parallel_hists(DATA_FILES, pid, coverage)
        mc_hists = run_parallel_hists(MC_FILES, pid, coverage)

        print_summary(pid, pid_label, data_hists, mc_hists, coverage)

        outpath = os.path.join(OUTDIR, fname)
        title = f"Vertex z comparison: {pid_label} (pid={pid}) [unit-normalized, log-y, HDI {100.0*coverage:.1f}%]"
        plot_2x3_canvas(title, data_hists, mc_hists, outpath, coverage)

        print(f"Wrote: {outpath}")
    #endfor


if __name__ == "__main__":
    main()
#endif