#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
calibration_sampling_fraction_vs_p_compare.py

Compare DATA vs DVCSGEN MC for the MEAN sampling fraction as a function of momentum p,
for electrons only (particle_pid == 11).

Sampling fraction definition:
  SF = (cal_energy_1 + cal_energy_4 + cal_energy_7) / p

We compute, in bins of p, the mean SF:
  <SF>(p-bin) = (1/N) * sum_i SF_i   for events in that p bin

Canvas layout (2x3):
  Top row:    Sp18 Inb | Sp18 Out | (empty)
  Bottom row: Fa18 Inb | Fa18 Out | Sp19 Inb

DATA: black points
MC:   red points

Parallel processing:
  - At most 5 workers HARD LIMIT.
  - We process ALL data files in parallel, then ALL mc files in parallel.

p range:
  - 2 to 8 (GeV)

Output:
  output/sampling_fraction_vs_p_electron.png

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
OUTDIR = "output"

# HARD LIMIT on workers
MAX_WORKERS_HARD = 5

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
class MeanProfile:
    p_centers: np.ndarray   # bin centers
    mean_sf: np.ndarray     # mean SF in each p bin (NaN where empty)
    err_sf: np.ndarray      # standard error on mean (NaN where empty)
    n_in_bin: np.ndarray    # counts per bin (int)


def fatal(msg: str) -> None:
    sys.stderr.write(f"\nFATAL: {msg}\n")
    sys.exit(1)
    #endif


def ensure_outdir(path: str) -> None:
    if not os.path.isdir(path):
        os.makedirs(path, exist_ok=True)
    #endif


def validate_inputs(args: argparse.Namespace) -> None:
    if args.pmax <= args.pmin:
        fatal(f"Invalid p range: pmin={args.pmin} must be < pmax={args.pmax}")
    #endif
    if args.pbins <= 0:
        fatal(f"Invalid --pbins={args.pbins}, must be > 0")
    #endif
    if args.max_workers <= 0 or args.max_workers > MAX_WORKERS_HARD:
        fatal(f"--max-workers must be in [1,{MAX_WORKERS_HARD}]. Got {args.max_workers}")
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


def compute_profile_for_file(args: Tuple[str, str, float, float, int]) -> Tuple[str, MeanProfile]:
    """
    Worker:
      - Select electrons (particle_pid == 11)
      - Require finite p, cal_energy_1/4/7 and p>0
      - Compute SF = (e1+e4+e7)/p
      - Bin in p and compute mean SF and standard error on mean in each bin

    Returns:
      (period_label, MeanProfile)
    """
    period_label, root_path, pmin, pmax, pbins = args

    try:
        with uproot.open(root_path) as f:
            if TREE_NAME not in f:
                raise KeyError(f"TTree '{TREE_NAME}' not found. Keys: {list(f.keys())}")
            #endif
            tree = f[TREE_NAME]

            required = ["particle_pid", "p", "cal_energy_1", "cal_energy_4", "cal_energy_7"]
            for br in required:
                if br not in tree.keys():
                    raise KeyError(f"Branch '{br}' not found in tree '{TREE_NAME}'.")
                #endif
            #endfor

            arr = tree.arrays(required, library="np")
            pid = arr["particle_pid"]
            p = arr["p"]
            e1 = arr["cal_energy_1"]
            e4 = arr["cal_energy_4"]
            e7 = arr["cal_energy_7"]

            mask = (
                (pid == PID_ELECTRON)
                & np.isfinite(p)
                & np.isfinite(e1)
                & np.isfinite(e4)
                & np.isfinite(e7)
                & (p > 0.0)
            )
            p_sel = p[mask].astype(np.float64)
            sf_sel = (e1[mask].astype(np.float64) + e4[mask].astype(np.float64) + e7[mask].astype(np.float64)) / p_sel

            # p binning
            edges = np.linspace(pmin, pmax, pbins + 1, dtype=np.float64)
            centers = 0.5 * (edges[:-1] + edges[1:])

            # Bin indices: 0..pbins-1, ignore under/overflow by mask
            bin_idx = np.digitize(p_sel, edges) - 1
            inrange = (bin_idx >= 0) & (bin_idx < pbins)
            bin_idx = bin_idx[inrange]
            sf_sel = sf_sel[inrange]

            # Aggregate sums and sums of squares for mean and std error
            n = np.zeros(pbins, dtype=np.int64)
            s1 = np.zeros(pbins, dtype=np.float64)
            s2 = np.zeros(pbins, dtype=np.float64)

            # Vectorized accumulation
            np.add.at(n, bin_idx, 1)
            np.add.at(s1, bin_idx, sf_sel)
            np.add.at(s2, bin_idx, sf_sel * sf_sel)

            mean = np.full(pbins, np.nan, dtype=np.float64)
            err = np.full(pbins, np.nan, dtype=np.float64)

            nonzero = n > 0
            mean[nonzero] = s1[nonzero] / n[nonzero]

            # Sample variance estimate in each bin (unbiased) then SEM
            # var = (sum(x^2) - sum(x)^2 / n) / (n-1)
            gt1 = n > 1
            var = np.zeros(pbins, dtype=np.float64)
            var[gt1] = (s2[gt1] - (s1[gt1] * s1[gt1]) / n[gt1]) / (n[gt1] - 1.0)
            var[var < 0.0] = 0.0  # numerical guard
            err[gt1] = np.sqrt(var[gt1] / n[gt1])
            # bins with n==1 keep err=NaN

            prof = MeanProfile(p_centers=centers, mean_sf=mean, err_sf=err, n_in_bin=n)
            return period_label, prof

    except Exception as e:
        tb = traceback.format_exc()
        raise RuntimeError(
            f"Failed processing file:\n"
            f"  period={period_label}\n"
            f"  path={root_path}\n"
            f"  error={repr(e)}\n\n{tb}"
        )


def run_parallel_profiles(file_map: Dict[str, str], pmin: float, pmax: float, pbins: int, max_workers: int) -> Dict[str, MeanProfile]:
    items = sorted(file_map.items(), key=lambda kv: kv[0])
    tasks = [(label, path, pmin, pmax, pbins) for (label, path) in items]

    n_workers = min(max_workers, len(tasks))
    if n_workers < 1:
        fatal("No tasks.")
    #endif

    out: Dict[str, MeanProfile] = {}
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        futures = [ex.submit(compute_profile_for_file, t) for t in tasks]
        for fut in as_completed(futures):
            period_label, prof = fut.result()
            out[period_label] = prof
        #endfor
    #endwith

    for period_label in file_map.keys():
        if period_label not in out:
            fatal(f"Missing profile for '{period_label}'")
        #endif
    #endfor

    return out


def plot_2x3_canvas(
    title: str,
    data_prof: Dict[str, MeanProfile],
    mc_prof: Dict[str, MeanProfile],
    outpath: str,
    pmin: float,
    pmax: float,
) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(16, 9), sharex=True, sharey=False)
    fig.suptitle(title, fontsize=18)

    # Disable all, then enable only used
    for r in range(2):
        for c in range(3):
            axes[r, c].axis("off")
        #endfor
    #endfor

    for period_label, (r, c) in PANEL_POS.items():
        ax = axes[r, c]
        ax.axis("on")

        dp = data_prof[period_label]
        mp = mc_prof[period_label]

        # Keep only finite mean points
        dmask = np.isfinite(dp.mean_sf)
        mmask = np.isfinite(mp.mean_sf)

        ax.errorbar(
            dp.p_centers[dmask],
            dp.mean_sf[dmask],
            yerr=dp.err_sf[dmask],
            fmt="o",
            color="black",
            markersize=3,
            linewidth=1.0,
            label="data",
        )
        ax.errorbar(
            mp.p_centers[mmask],
            mp.mean_sf[mmask],
            yerr=mp.err_sf[mmask],
            fmt="o",
            color="red",
            markersize=3,
            linewidth=1.0,
            label="mc",
        )

        ax.set_title(period_label, fontsize=13)
        ax.set_xlim(pmin, pmax)
        ax.set_xlabel("p (GeV)", fontsize=12)
        ax.set_ylabel("mean sampling fraction", fontsize=12)
        ax.grid(True, alpha=0.25)
        ax.legend(loc="upper right", fontsize=10, frameon=True)

    #endfor

    axes[0, 2].axis("off")

    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.95])
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Mean sampling fraction vs p (data vs mc), electrons only.")
    p.add_argument("--pmin", type=float, default=2.0, help="Minimum p (GeV).")
    p.add_argument("--pmax", type=float, default=8.0, help="Maximum p (GeV).")
    p.add_argument("--pbins", type=int, default=60, help="Number of p bins.")
    p.add_argument("--max-workers", type=int, default=5, help="Max parallel workers (hard limit 5).")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    validate_inputs(args)
    ensure_outdir(OUTDIR)

    data_prof = run_parallel_profiles(DATA_FILES, args.pmin, args.pmax, args.pbins, args.max_workers)
    mc_prof = run_parallel_profiles(MC_FILES, args.pmin, args.pmax, args.pbins, args.max_workers)

    outpath = os.path.join(OUTDIR, "sampling_fraction_vs_p_electron.png")
    title = "Mean sampling fraction vs p: electrons (pid=11), (cal_energy_1+4+7)/p"
    plot_2x3_canvas(
        title=title,
        data_prof=data_prof,
        mc_prof=mc_prof,
        outpath=outpath,
        pmin=args.pmin,
        pmax=args.pmax,
    )

    print(f"Wrote: {outpath}")
    print("Done.")


if __name__ == "__main__":
    main()
#endif