#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
calibration_vertex_vz_compare.py

Make 2x3 canvases comparing DATA vs DVCSGEN MC for particle vertex z (particle_vz),
separately for:
  - electron (pid = 11)
  - proton   (pid = 2212)
  - photon   (pid = 22)

Layout (2x3):
  Top row:    Sp18 Inb | Sp18 Out | (empty)
  Bottom row: Fa18 Inb | Fa18 Out | Sp19 Inb

DATA histogram: black
MC histogram:   red

IMPORTANT UPDATE:
  - Histograms are normalized to their integral (each period + dataset separately).

Parallel processing:
  - Up to 5 workers HARD LIMIT.
  - We process ALL data files in parallel, then ALL mc files in parallel.

Output:
  output/vz_electron.png
  output/vz_proton.png
  output/vz_photon.png

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
VZ_MIN = -20.0
VZ_MAX = 20.0
N_BINS = 200

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
# INTERNALS
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class Hist1D:
    counts: np.ndarray
    edges: np.ndarray
    n_selected: int


def fatal(msg: str) -> None:
    sys.stderr.write(f"\nFATAL: {msg}\n")
    sys.exit(1)


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

    # Ensure files exist
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

    # Ensure period keys match exactly
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

    # Ensure PANEL_POS covers all periods
    for k in DATA_FILES.keys():
        if k not in PANEL_POS:
            fatal(f"PANEL_POS missing placement for period '{k}'")
        #endif
    #endfor
    if (0, 2) in PANEL_POS.values():
        fatal("PANEL_POS should NOT assign any period to (0,2); that panel is reserved as empty.")
    #endif


def normalize_to_integral(counts: np.ndarray) -> np.ndarray:
    """
    Normalize histogram counts to unit integral.
    If integral is zero, returns an array of zeros (fail-safe, deterministic).
    """
    integral = float(np.sum(counts))
    if integral <= 0.0:
        return np.zeros_like(counts, dtype=np.float64)
    #endif
    return counts.astype(np.float64) / integral


def compute_hist_for_file(args: Tuple[str, str, int, str]) -> Tuple[str, str, int, Hist1D]:
    """
    Worker function (must be top-level for multiprocessing pickling).

    args = (period_label, root_path, pid, tree_name)

    Returns:
      (period_label, root_path, pid, Hist1D)

    NOTE:
      - The returned histogram counts are normalized to unit integral.
      - n_selected is still the raw number of selected entries (pre-normalization),
        so legends can show the statistics.
    """
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

            h = Hist1D(
                counts=counts_norm.astype(np.float64),
                edges=edges.astype(np.float64),
                n_selected=int(vz_sel.size),
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
    """
    Compute normalized histograms for all periods in file_map for a given pid.
    Parallelized with a hard max of MAX_WORKERS.

    Returns:
      dict period_label -> Hist1D
    """
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
    """
    Make a 2x3 matplotlib figure with overlays (DATA vs MC) for each period.
    Counts are already normalized to unit integral.
    """
    fig, axes = plt.subplots(2, 3, figsize=(16, 9), sharex=True, sharey=False)
    fig.suptitle(title, fontsize=18)

    # Turn off all axes first
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

        ax.step(centers, dh.counts, where="mid", color="black", linewidth=1.2, label=f"data (N={dh.n_selected})")
        ax.step(centers, mh.counts, where="mid", color="red", linewidth=1.2, label=f"mc (N={mh.n_selected})")

        ax.set_title(period_label, fontsize=13)
        ax.set_xlim(VZ_MIN, VZ_MAX)
        ax.set_xlabel("particle_vz (cm)", fontsize=12)
        ax.set_ylabel("normalized counts", fontsize=12)
        ax.grid(True, alpha=0.25)

        ax.legend(loc="upper right", fontsize=10, frameon=True)

        # Optional: ensure y starts at 0 for normalized plots
        ax.set_ylim(bottom=0.0)

    #endfor

    # Keep top-right empty
    axes[0, 2].axis("off")

    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.95])
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def main() -> None:
    validate_inputs()
    ensure_outdir(OUTDIR)

    particles = [
        (11, "electron", "vz_electron.png"),
        (2212, "proton", "vz_proton.png"),
        (22, "photon", "vz_photon.png"),
    ]

    for pid, pid_label, fname in particles:
        data_hists = run_parallel_hists(DATA_FILES, pid, TREE_NAME)
        mc_hists = run_parallel_hists(MC_FILES, pid, TREE_NAME)

        outpath = os.path.join(OUTDIR, fname)
        title = f"Vertex z comparison: {pid_label} (pid={pid}) [unit-normalized]"
        plot_2x3_canvas(title, data_hists, mc_hists, outpath)

        print(f"Wrote: {outpath}")
    #endfor


if __name__ == "__main__":
    main()
#endif