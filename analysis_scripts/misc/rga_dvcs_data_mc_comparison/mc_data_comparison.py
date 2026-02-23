#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
calibration_vertex_vz_compare.py

Make 2x3 canvases comparing DATA vs DVCSGEN MC for particle vertex z (particle_vz),
separately for:
  - electron (pid = 11)
  - proton   (pid = 2212)

NOTE:
  - Photon plots are REMOVED (you said you don't cut on photon vertex).

Layout (2x3):
  Top row:    Sp18 Inb | Sp18 Out | (empty)
  Bottom row: Fa18 Inb | Fa18 Out | Sp19 Inb

DATA histogram: black
MC histogram:   red

Histogram normalization:
  - Each histogram is normalized to its integral (unit area) over the plotted range.

Parallel processing:
  - Up to 5 workers HARD LIMIT.
  - We process ALL data files in parallel, then ALL mc files in parallel.

Features:
  1) Draw vertical cut lines for your pass2-derived vertex cuts (per RGA period).
     - electron (pid=11): charge<0 window from your Java (per period).
     - proton   (pid=2212): charge>0 window from your Java (per period).
  2) Legend includes the percentage of selected events inside the cut window.
  3) Print out the MODE bin entry info (per period, for data and mc).
  4) Cut description text is in the bottom-right of each subplot.
  5) Y-axis is LOG scale.
  6) NEW: Extend x-axis (and histogram range) to [-14, 14] (cm).
  7) NEW: For each histogram, compute the points (bin centers) on the LEFT and RIGHT
          side of the peak (mode bin) whose normalized bin content is closest to 2%.
          These are printed to the console for both DATA and MC.

Interpretation of "closest to 2%":
  - We use the normalized histogram bin content (unit area normalization),
    so "2%" corresponds to y = 0.02.
  - We find, among bins strictly left of the peak bin, the bin whose normalized
    content minimizes |y - 0.02|. Same on the right side.

Output:
  output/vz_electron.png
  output/vz_proton.png

Dependencies:
  pip install uproot numpy matplotlib
"""

import os
import sys
import traceback
from dataclasses import dataclass
from typing import Dict, Tuple

import numpy as np
import uproot
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor, as_completed


# -----------------------------------------------------------------------------
# USER SETTINGS
# -----------------------------------------------------------------------------

TREE_NAME = "PhysicsEvents"  # set this to the TTree name in your ROOT files

# Histogram binning for particle_vz
VZ_MIN = -14.0
VZ_MAX = 14.0
N_BINS = 200

# Target normalized level for "2%"
TARGET_LEVEL = 0.02

# Hard limit on workers
MAX_WORKERS = 5

# Output directory
OUTDIR = "output"

# File map: period key -> ROOT path
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

# Panel placement in the 2x3 canvas
# (row, col) with row=0 top, row=1 bottom; col=0..2
PANEL_POS = {
    "Sp18 Inb": (0, 0),
    "Sp18 Out": (0, 1),
    # (0,2) intentionally empty
    "Fa18 Inb": (1, 0),
    "Fa18 Out": (1, 1),
    "Sp19 Inb": (1, 2),
}

# -----------------------------------------------------------------------------
# PASS2 DERIVED VERTEX CUTS (mapped to the 5 RGA periods)
# -----------------------------------------------------------------------------

VTX_CUTS_POS = {
    "Sp18 Inb": (-7.8790, 1.5150),
    "Sp18 Out": (-6.6667, 2.7273),
    "Fa18 Inb": (-8.4850, 0.6060),
    "Fa18 Out": (-6.9700, 1.8180),
    "Sp19 Inb": (-8.4850, 0.6060),
}

VTX_CUTS_NEG = {
    "Sp18 Inb": (-6.0606, 1.8182),
    "Sp18 Out": (-7.2730, 0.9091),
    "Fa18 Inb": (-6.3640, 1.5150),
    "Fa18 Out": (-7.8790, 0.3030),
    "Sp19 Inb": (-6.3640, 1.5150),
}


# -----------------------------------------------------------------------------
# INTERNALS
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class Hist1D:
    counts: np.ndarray          # normalized counts (unit area)
    edges: np.ndarray
    n_selected: int
    n_in_cut: int
    frac_in_cut: float
    cut_low: float
    cut_high: float

    # Peak (mode) info from RAW counts
    peak_index: int
    mode_vz: float
    mode_count: int
    mode_frac: float

    # Closest-to-2% points, left and right of peak (using normalized counts)
    left_2pct_vz: float
    left_2pct_y: float
    right_2pct_vz: float
    right_2pct_y: float


def fatal(msg: str) -> None:
    sys.stderr.write(f"\nFATAL: {msg}\n")
    sys.exit(1)
    #endif


def ensure_outdir(path: str) -> None:
    if not os.path.isdir(path):
        os.makedirs(path, exist_ok=True)
    #endif


def validate_inputs() -> None:
    if VZ_MAX <= VZ_MIN:
        fatal(f"Invalid VZ range: VZ_MIN={VZ_MIN} must be < VZ_MAX={VZ_MAX}")
    #endif
    if N_BINS <= 0:
        fatal(f"Invalid N_BINS={N_BINS}, must be > 0")
    #endif
    if MAX_WORKERS <= 0:
        fatal(f"Invalid MAX_WORKERS={MAX_WORKERS}, must be > 0")
    #endif
    if TARGET_LEVEL <= 0.0:
        fatal(f"Invalid TARGET_LEVEL={TARGET_LEVEL}, must be > 0")
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

    data_keys = set(DATA_FILES.keys())
    mc_keys = set(MC_FILES.keys())
    if data_keys != mc_keys:
        missing_in_mc = sorted(list(data_keys - mc_keys))
        missing_in_data = sorted(list(mc_keys - data_keys))
        fatal(
            "DATA_FILES and MC_FILES period keys do not match.\n"
            f"  Missing in MC: {missing_in_mc}\n"
            f"  Missing in DATA: {missing_in_data}"
        )
    #endif

    for k in DATA_FILES.keys():
        if k not in PANEL_POS:
            fatal(f"PANEL_POS missing placement for period '{k}'")
        #endif
        if k not in VTX_CUTS_POS or k not in VTX_CUTS_NEG:
            fatal(f"Vertex cut maps missing period '{k}'")
        #endif
    #endfor

    if (0, 2) in PANEL_POS.values():
        fatal("PANEL_POS should NOT assign any period to (0,2); that panel is reserved as empty.")
    #endif


def normalize_to_integral(counts: np.ndarray) -> np.ndarray:
    integral = float(np.sum(counts))
    if integral <= 0.0:
        return np.zeros_like(counts, dtype=np.float64)
    #endif
    return counts.astype(np.float64) / integral


def cut_window_for_pid_and_period(pid: int, period_label: str) -> Tuple[float, float]:
    if pid == 11:
        return VTX_CUTS_NEG[period_label]
    elif pid == 2212:
        return VTX_CUTS_POS[period_label]
    #endif
    return (-9.0, 2.0)


def compute_mode_from_hist(counts_raw: np.ndarray, edges: np.ndarray, n_selected: int) -> Tuple[int, float, int, float]:
    """
    Mode defined as the bin with maximum raw count.
    Returns (peak_index, mode_vz_center, mode_count, mode_frac).
    """
    if counts_raw.size == 0:
        return 0, 0.0, 0, 0.0
    #endif

    peak_index = int(np.argmax(counts_raw))
    mode_count = int(counts_raw[peak_index])

    centers = 0.5 * (edges[:-1] + edges[1:])
    mode_vz = float(centers[peak_index])

    if n_selected > 0:
        mode_frac = float(mode_count) / float(n_selected)
    else:
        mode_frac = 0.0
    #endif

    return peak_index, mode_vz, mode_count, mode_frac


def compute_closest_to_level_left_right(
    counts_norm: np.ndarray,
    edges: np.ndarray,
    peak_index: int,
    level: float,
) -> Tuple[float, float, float, float]:
    """
    Find the bin (center) on the LEFT of peak_index with normalized content closest to 'level',
    and similarly on the RIGHT.

    Returns:
      (left_vz, left_y, right_vz, right_y)

    If left or right side has no bins (degenerate), we return (nan, nan, nan, nan).
    """
    centers = 0.5 * (edges[:-1] + edges[1:])

    # Left side: indices [0, peak_index-1]
    if peak_index <= 0:
        left_vz = float("nan")
        left_y = float("nan")
    else:
        left_slice = counts_norm[:peak_index]
        left_diff = np.abs(left_slice - level)
        i_left = int(np.argmin(left_diff))
        left_vz = float(centers[i_left])
        left_y = float(left_slice[i_left])
    #endif

    # Right side: indices [peak_index+1, end)
    if peak_index >= (counts_norm.size - 1):
        right_vz = float("nan")
        right_y = float("nan")
    else:
        right_slice = counts_norm[peak_index + 1 :]
        right_diff = np.abs(right_slice - level)
        i_rel = int(np.argmin(right_diff))
        i_right = peak_index + 1 + i_rel
        right_vz = float(centers[i_right])
        right_y = float(counts_norm[i_right])
    #endif

    return left_vz, left_y, right_vz, right_y


def compute_hist_for_file(args: Tuple[str, str, int, str]) -> Tuple[str, str, int, Hist1D]:
    period_label, root_path, pid, tree_name = args

    try:
        with uproot.open(root_path) as f:
            if tree_name not in f:
                raise KeyError(
                    f"TTree '{tree_name}' not found in file. Available keys: {list(f.keys())}"
                )
            #endif

            tree = f[tree_name]

            required = ["particle_pid", "particle_vz"]
            for br in required:
                if br not in tree.keys():
                    raise KeyError(
                        f"Branch '{br}' not found in TTree '{tree_name}'. "
                        f"Available branches (first 50): {list(tree.keys())[:50]}"
                    )
                #endif
            #endfor

            arrays = tree.arrays(required, library="np")
            pids = arrays["particle_pid"]
            vz = arrays["particle_vz"]

            mask = (pids == pid) & np.isfinite(vz)
            vz_sel = vz[mask]

            counts_raw, edges = np.histogram(vz_sel, bins=N_BINS, range=(VZ_MIN, VZ_MAX))
            counts_norm = normalize_to_integral(counts_raw)

            cut_low, cut_high = cut_window_for_pid_and_period(pid, period_label)

            if vz_sel.size > 0:
                in_cut = (vz_sel > cut_low) & (vz_sel < cut_high)
                n_in_cut = int(np.count_nonzero(in_cut))
                frac_in_cut = float(n_in_cut) / float(vz_sel.size)
            else:
                n_in_cut = 0
                frac_in_cut = 0.0
            #endif

            peak_index, mode_vz, mode_count, mode_frac = compute_mode_from_hist(counts_raw, edges, int(vz_sel.size))

            left_vz, left_y, right_vz, right_y = compute_closest_to_level_left_right(
                counts_norm=counts_norm,
                edges=edges,
                peak_index=peak_index,
                level=TARGET_LEVEL,
            )

            h = Hist1D(
                counts=counts_norm.astype(np.float64),
                edges=edges.astype(np.float64),
                n_selected=int(vz_sel.size),
                n_in_cut=n_in_cut,
                frac_in_cut=frac_in_cut,
                cut_low=float(cut_low),
                cut_high=float(cut_high),
                peak_index=int(peak_index),
                mode_vz=float(mode_vz),
                mode_count=int(mode_count),
                mode_frac=float(mode_frac),
                left_2pct_vz=float(left_vz),
                left_2pct_y=float(left_y),
                right_2pct_vz=float(right_vz),
                right_2pct_y=float(right_y),
            )
            return (period_label, root_path, pid, h)

    except Exception as e:
        tb = traceback.format_exc()
        raise RuntimeError(
            f"Failed processing file:\n"
            f"  period = {period_label}\n"
            f"  path   = {root_path}\n"
            f"  pid    = {pid}\n"
            f"  tree   = {tree_name}\n"
            f"  error  = {repr(e)}\n\n{tb}"
        )


def run_parallel_hists(file_map: Dict[str, str], pid: int, tree_name: str) -> Dict[str, Hist1D]:
    items = sorted(file_map.items(), key=lambda kv: kv[0])
    tasks = [(label, path, pid, tree_name) for (label, path) in items]

    n_workers = min(MAX_WORKERS, len(tasks))
    if n_workers < 1:
        fatal("No tasks to run in run_parallel_hists()")
    #endif

    out: Dict[str, Hist1D] = {}

    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        futures = [ex.submit(compute_hist_for_file, t) for t in tasks]
        for fut in as_completed(futures):
            period_label, root_path, pid_ret, hist = fut.result()
            if period_label in out:
                fatal(f"Duplicate histogram result for period '{period_label}'")
            #endif
            out[period_label] = hist
        #endfor
    #endwith

    for period_label in file_map.keys():
        if period_label not in out:
            fatal(f"Missing histogram result for period '{period_label}' (pid={pid})")
        #endif
    #endfor

    return out


def plot_2x3_canvas(
    title: str,
    data_hists: Dict[str, Hist1D],
    mc_hists: Dict[str, Hist1D],
    outpath: str,
) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(16, 9), sharex=True, sharey=False)
    fig.suptitle(title, fontsize=18)

    for r in range(2):
        for c in range(3):
            axes[r, c].axis("off")
        #endfor
    #endfor

    any_period = next(iter(data_hists.keys()))
    edges = data_hists[any_period].edges
    centers = 0.5 * (edges[:-1] + edges[1:])

    for period_label, (r, c) in PANEL_POS.items():
        ax = axes[r, c]
        ax.axis("on")

        dh = data_hists[period_label]
        mh = mc_hists[period_label]

        if dh.edges.shape != mh.edges.shape or not np.allclose(dh.edges, mh.edges):
            fatal(f"Histogram edges mismatch for period '{period_label}'")
        #endif

        cut_low = dh.cut_low
        cut_high = dh.cut_high

        data_label = f"data: {100.0 * dh.frac_in_cut:.2f}% in-cut (N={dh.n_selected})"
        mc_label = f"mc: {100.0 * mh.frac_in_cut:.2f}% in-cut (N={mh.n_selected})"

        ax.step(centers, dh.counts, where="mid", color="black", linewidth=1.2, label=data_label)
        ax.step(centers, mh.counts, where="mid", color="red", linewidth=1.2, label=mc_label)

        ax.axvline(cut_low, color="black", linestyle="--", linewidth=1.0)
        ax.axvline(cut_high, color="black", linestyle="--", linewidth=1.0)

        ax.set_title(period_label, fontsize=13)
        ax.set_xlim(VZ_MIN, VZ_MAX)
        ax.set_xlabel("particle_vz (cm)", fontsize=12)
        ax.set_ylabel("normalized counts", fontsize=12)
        ax.grid(True, alpha=0.25)
        ax.legend(loc="upper right", fontsize=9, frameon=True)

        # Log y
        ax.set_yscale("log")
        positive_vals = np.concatenate([dh.counts[dh.counts > 0.0], mh.counts[mh.counts > 0.0]])
        if positive_vals.size > 0:
            y_min = float(np.min(positive_vals))
        else:
            y_min = 1e-12
        #endif
        ax.set_ylim(bottom=y_min)

        # Cut description bottom-right
        ax.text(
            0.98,
            0.06,
            f"cut: ({cut_low:.3f}, {cut_high:.3f}) (cm)",
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


def print_mode_and_2pct_summary(pid: int, pid_label: str, data_hists: Dict[str, Hist1D], mc_hists: Dict[str, Hist1D]) -> None:
    print("")
    print("------------------------------------------------------------")
    print(f"SUMMARY: {pid_label} (pid={pid})")
    print(f"  Histogram range = ({VZ_MIN:.1f}, {VZ_MAX:.1f}) (cm), N_BINS={N_BINS}, target level = {TARGET_LEVEL:.3f}")
    print("  mode is peak bin (max RAW count). 2% points are closest bins to y=0.02 on each side of the peak.")
    print("------------------------------------------------------------")
    periods = sorted(PANEL_POS.keys(), key=lambda k: (PANEL_POS[k][0], PANEL_POS[k][1]))
    for period_label in periods:
        dh = data_hists[period_label]
        mh = mc_hists[period_label]

        print(f"{period_label}:")
        print(f"  cut_window = ({dh.cut_low:.4f}, {dh.cut_high:.4f}) (cm)")
        print(f"  data: N={dh.n_selected}  in_cut={100.0*dh.frac_in_cut:.3f}%  mode_vz={dh.mode_vz:.4f} (cm)  mode_count={dh.mode_count}  mode_frac={100.0*dh.mode_frac:.3f}%")
        print(f"        left@2pct:  vz={dh.left_2pct_vz:.4f} (cm)  y={dh.left_2pct_y:.6f}   right@2pct: vz={dh.right_2pct_vz:.4f} (cm)  y={dh.right_2pct_y:.6f}")
        print(f"  mc:   N={mh.n_selected}  in_cut={100.0*mh.frac_in_cut:.3f}%  mode_vz={mh.mode_vz:.4f} (cm)  mode_count={mh.mode_count}  mode_frac={100.0*mh.mode_frac:.3f}%")
        print(f"        left@2pct:  vz={mh.left_2pct_vz:.4f} (cm)  y={mh.left_2pct_y:.6f}   right@2pct: vz={mh.right_2pct_vz:.4f} (cm)  y={mh.right_2pct_y:.6f}")
    #endfor


def main() -> None:
    validate_inputs()
    ensure_outdir(OUTDIR)

    particles = [
        (11, "electron", "vz_electron.png"),
        (2212, "proton", "vz_proton.png"),
    ]

    for pid, pid_label, fname in particles:
        data_hists = run_parallel_hists(DATA_FILES, pid, TREE_NAME)
        mc_hists = run_parallel_hists(MC_FILES, pid, TREE_NAME)

        print_mode_and_2pct_summary(pid, pid_label, data_hists, mc_hists)

        outpath = os.path.join(OUTDIR, fname)
        title = f"Vertex z comparison: {pid_label} (pid={pid}) [unit-normalized, log-y]"
        plot_2x3_canvas(title, data_hists, mc_hists, outpath)

        print(f"Wrote: {outpath}")
    #endfor


if __name__ == "__main__":
    main()
#endif