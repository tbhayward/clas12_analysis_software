#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
calibration_cal_energy_1_compare.py

Compare DATA vs DVCSGEN MC for calorimeter deposited energy in layer 1:
  - cal_energy_1

Selections:
  - electrons only: particle_pid == 11

Canvas layout (2x3):
  Top row:    Sp18 Inb | Sp18 Out | (empty)
  Bottom row: Fa18 Inb | Fa18 Out | Sp19 Inb

DATA histogram: black
MC histogram:   red

Parallel processing:
  - At most 5 workers HARD LIMIT.
  - We process ALL data files in parallel, then ALL mc files in parallel.

Output:
  output/cal_energy_1_electron.png

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
# STATIC SETTINGS
# -----------------------------------------------------------------------------

TREE_NAME = "PhysicsEvents"

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
    # (0,2) intentionally empty
    "Fa18 Inb": (1, 0),
    "Fa18 Out": (1, 1),
    "Sp19 Inb": (1, 2),
}

PID_ELECTRON = 11


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
    #endif


def ensure_outdir(path: str) -> None:
    if not os.path.isdir(path):
        os.makedirs(path, exist_ok=True)
    #endif


def normalize_to_integral(counts: np.ndarray) -> np.ndarray:
    s = float(np.sum(counts))
    if s <= 0.0:
        return np.zeros_like(counts, dtype=np.float64)
    #endif
    return counts.astype(np.float64) / s


def validate_inputs(args: argparse.Namespace) -> None:
    if args.emax <= args.emin:
        fatal(f"Invalid energy range: emin={args.emin} must be < emax={args.emax}")
    #endif
    if args.bins <= 0:
        fatal(f"Invalid --bins={args.bins}, must be > 0")
    #endif
    if args.max_workers <= 0 or args.max_workers > 5:
        fatal(f"--max-workers must be in [1,5]. Got {args.max_workers}")
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


def compute_hist_for_file(args: Tuple[str, str, float, float, int, bool]) -> Tuple[str, Hist1D]:
    """
    Worker:
      - Select electrons (particle_pid == 11)
      - Read cal_energy_1
      - Histogram in [emin, emax] with nbins
      - Optionally normalize to unit integral
    """
    period_label, root_path, emin, emax, nbins, do_norm = args

    try:
        with uproot.open(root_path) as f:
            if TREE_NAME not in f:
                raise KeyError(f"TTree '{TREE_NAME}' not found. Keys: {list(f.keys())}")
            #endif

            tree = f[TREE_NAME]
            required = ["particle_pid", "cal_energy_1"]
            for br in required:
                if br not in tree.keys():
                    raise KeyError(f"Branch '{br}' not found in tree '{TREE_NAME}'.")
                #endif
            #endfor

            arr = tree.arrays(required, library="np")
            pid = arr["particle_pid"]
            e1 = arr["cal_energy_1"]

            mask = (pid == PID_ELECTRON) & np.isfinite(e1)
            e_sel = e1[mask].astype(np.float64)

            counts, edges = np.histogram(e_sel, bins=nbins, range=(emin, emax))
            if do_norm:
                counts = normalize_to_integral(counts)
            else:
                counts = counts.astype(np.float64)
            #endif

            h = Hist1D(counts=counts, edges=edges.astype(np.float64), n_selected=int(e_sel.size))
            return period_label, h

    except Exception as e:
        tb = traceback.format_exc()
        raise RuntimeError(
            f"Failed processing file:\n"
            f"  period={period_label}\n"
            f"  path={root_path}\n"
            f"  error={repr(e)}\n\n{tb}"
        )


def run_parallel_hists(file_map: Dict[str, str], emin: float, emax: float, nbins: int, do_norm: bool, max_workers: int) -> Dict[str, Hist1D]:
    items = sorted(file_map.items(), key=lambda kv: kv[0])
    tasks = [(label, path, emin, emax, nbins, do_norm) for (label, path) in items]

    n_workers = min(max_workers, len(tasks))
    if n_workers < 1:
        fatal("No tasks.")
    #endif

    out: Dict[str, Hist1D] = {}
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        futures = [ex.submit(compute_hist_for_file, t) for t in tasks]
        for fut in as_completed(futures):
            period_label, hist = fut.result()
            out[period_label] = hist
        #endfor
    #endwith

    for period_label in file_map.keys():
        if period_label not in out:
            fatal(f"Missing histogram result for '{period_label}'")
        #endif
    #endfor

    return out


def plot_2x3_canvas(
    title: str,
    data_hists: Dict[str, Hist1D],
    mc_hists: Dict[str, Hist1D],
    outpath: str,
    emin: float,
    emax: float,
    do_norm: bool,
    logy: bool,
) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(16, 9), sharex=True, sharey=False)
    fig.suptitle(title, fontsize=18)

    # Disable all, then enable only used
    for r in range(2):
        for c in range(3):
            axes[r, c].axis("off")
        #endfor
    #endfor

    # Common bin centers
    any_period = next(iter(data_hists.keys()))
    edges = data_hists[any_period].edges
    centers = 0.5 * (edges[:-1] + edges[1:])

    y_label = "normalized counts" if do_norm else "counts"

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
        ax.set_xlim(emin, emax)
        ax.set_xlabel("cal_energy_1 (GeV)", fontsize=12)
        ax.set_ylabel(y_label, fontsize=12)
        ax.grid(True, alpha=0.25)
        ax.legend(loc="upper right", fontsize=10, frameon=True)

        if logy:
            ax.set_yscale("log")
            positive_vals = np.concatenate([dh.counts[dh.counts > 0.0], mh.counts[mh.counts > 0.0]])
            y_min = float(np.min(positive_vals)) if positive_vals.size > 0 else 1e-12
            ax.set_ylim(bottom=y_min)
        else:
            ax.set_ylim(bottom=0.0)
        #endif

    #endfor

    # Keep top-right empty
    axes[0, 2].axis("off")

    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.95])
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compare cal_energy_1 between data and mc for electrons.")
    p.add_argument("--emin", type=float, default=0.0, help="Histogram minimum for cal_energy_1 (GeV).")
    p.add_argument("--emax", type=float, default=6.0, help="Histogram maximum for cal_energy_1 (GeV).")
    p.add_argument("--bins", type=int, default=240, help="Number of histogram bins.")
    p.add_argument("--normalize", action="store_true", help="Normalize each histogram to unit integral.")
    p.add_argument("--logy", action="store_true", help="Use log scale on y axis.")
    p.add_argument("--max-workers", type=int, default=5, help="Max parallel workers (hard limit 5).")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    validate_inputs(args)
    ensure_outdir(OUTDIR)

    # Data first, then MC (as you’ve been doing)
    data_hists = run_parallel_hists(DATA_FILES, args.emin, args.emax, args.bins, args.normalize, args.max_workers)
    mc_hists = run_parallel_hists(MC_FILES, args.emin, args.emax, args.bins, args.normalize, args.max_workers)

    outpath = os.path.join(OUTDIR, "cal_energy_1_electron.png")
    title = "Calorimeter energy deposit comparison: electrons (pid=11), cal_energy_1"
    if args.normalize:
        title += " [unit-normalized]"
    if args.logy:
        title += " [log-y]"
    #endif

    plot_2x3_canvas(
        title=title,
        data_hists=data_hists,
        mc_hists=mc_hists,
        outpath=outpath,
        emin=args.emin,
        emax=args.emax,
        do_norm=args.normalize,
        logy=args.logy,
    )

    print(f"Wrote: {outpath}")
    print("Done.")


if __name__ == "__main__":
    main()
#endif