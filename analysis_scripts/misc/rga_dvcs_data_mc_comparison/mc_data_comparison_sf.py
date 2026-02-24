#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
calibration_sampling_fraction_vs_p_compare.py

Make a 2x3 canvas comparing DATA vs DVCSGEN MC for the MEAN sampling fraction
as a function of momentum p, for ELECTRONS only (particle_pid == 11).

Sampling fraction:
  sf = (cal_energy_1 + cal_energy_4 + cal_energy_7) / p

What is plotted:
  - In p-bins, plot mean(sf) with standard error on the mean (SEM)
  - DATA in black, MC in red

Cuts and survival rates (this version):
  A) Sector- and period-dependent polynomial band cut (from your Java), event-by-event.
     We compute pass fraction = N_pass / N_used for that cut.
  B) Additional simple cut: sf >= 0.19 (requested), also event-by-event.
     We compute pass fraction = N_pass_simple / N_used for that cut.

Overlay:
  - Dashed blue lines show a representative (median-over-sectors) lower/upper band for the
    polynomial cut (DATA period definition).
  - (Optional visual) we draw a horizontal dashed blue line at sf = 0.19 to show the simple cut.

Axes:
  - p from 2 to 8 (GeV)
  - y fixed to [0.1, 0.35]

Parallel processing:
  - At most 5 workers HARD LIMIT.
  - We process all DATA in parallel, then all MC in parallel.

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
MAX_WORKERS_HARD = 5

P_MIN_DEFAULT = 2.0
P_MAX_DEFAULT = 8.0

Y_MIN = 0.1
Y_MAX = 0.35

PID_ELECTRON = 11

# Additional simple SF cut (requested)
SF_SIMPLE_MIN = 0.19

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


# -----------------------------------------------------------------------------
# CUT COEFFICIENTS FROM YOUR JAVA (RGA only + MC block)
# Each is 6 sectors x 3 coefficients for quadratic: a0 + a1*p + a2*p^2
# -----------------------------------------------------------------------------

def coeffs_rga_fa18_inb() -> Tuple[np.ndarray, np.ndarray]:
    lo = np.zeros((6, 3), dtype=np.float64)
    hi = np.zeros((6, 3), dtype=np.float64)

    lo[0] = [0.182257, 0.007442, -0.000758]
    hi[0] = [0.304421, -0.002123, -0.000286]
    lo[1] = [0.179615, 0.009837, -0.000981]
    hi[1] = [0.306626, -0.001156, -0.000418]
    lo[2] = [0.179768, 0.011258, -0.001113]
    hi[2] = [0.304496, 0.000808, -0.000712]
    lo[3] = [0.179226, 0.010001, -0.000851]
    hi[3] = [0.310082, -0.002499, -0.000111]
    lo[4] = [0.174914, 0.011152, -0.001008]
    hi[4] = [0.315217, -0.006281, 0.000256]
    lo[5] = [0.182785, 0.008533, -0.000850]
    hi[5] = [0.307112, -0.000586, -0.000559]

    return lo, hi
    #endif


def coeffs_rga_fa18_out() -> Tuple[np.ndarray, np.ndarray]:
    lo = np.zeros((6, 3), dtype=np.float64)
    hi = np.zeros((6, 3), dtype=np.float64)

    lo[0] = [0.186279, 0.006740, -0.000434]
    hi[0] = [0.303829, -0.003178, 0.000029]
    lo[1] = [0.186446, 0.006058, -0.000374]
    hi[1] = [0.298393, -0.001543, -0.000090]
    lo[2] = [0.186023, 0.008322, -0.000555]
    hi[2] = [0.300515, -0.000776, -0.000302]
    lo[3] = [0.186014, 0.006690, -0.000387]
    hi[3] = [0.299888, -0.001181, -0.000234]
    lo[4] = [0.187277, 0.004536, -0.000115]
    hi[4] = [0.297246, -0.002602, 0.000075]
    lo[5] = [0.181060, 0.008314, -0.000519]
    hi[5] = [0.297379, -0.000577, -0.000283]

    return lo, hi
    #endif


def coeffs_rga_sp19_inb() -> Tuple[np.ndarray, np.ndarray]:
    lo = np.zeros((6, 3), dtype=np.float64)
    hi = np.zeros((6, 3), dtype=np.float64)

    lo[0] = [0.166035, 0.017046, -0.001806]
    hi[0] = [0.323043, -0.010838, 0.000676]
    lo[1] = [0.187470, 0.004908, -0.000521]
    hi[1] = [0.296504, 0.004220, -0.001003]
    lo[2] = [0.178773, 0.012186, -0.001253]
    hi[2] = [0.302329, 0.001279, -0.000800]
    lo[3] = [0.176627, 0.011507, -0.001057]
    hi[3] = [0.293717, 0.003929, -0.000917]
    lo[4] = [0.173851, 0.012135, -0.001183]
    hi[4] = [0.300202, 0.001455, -0.000786]
    lo[5] = [0.183544, 0.007833, -0.000789]
    hi[5] = [0.303698, 0.000896, -0.000796]

    return lo, hi
    #endif


def coeffs_rga_sp18_inb() -> Tuple[np.ndarray, np.ndarray]:
    lo = np.zeros((6, 3), dtype=np.float64)
    hi = np.zeros((6, 3), dtype=np.float64)

    lo[0] = [0.187361, 0.002859, -0.000360]
    hi[0] = [0.291145, 0.010578, -0.001185]
    lo[1] = [0.184836, 0.005616, -0.000512]
    hi[1] = [0.303975, 0.003104, -0.000616]
    lo[2] = [0.176353, 0.013075, -0.001351]
    hi[2] = [0.310923, -0.001850, -0.000212]
    lo[3] = [0.189493, 0.002891, -0.000338]
    hi[3] = [0.326210, -0.010605, 0.001433]
    lo[4] = [0.167201, 0.015198, -0.001506]
    hi[4] = [0.325967, -0.009674, 0.000799]
    lo[5] = [0.178730, 0.010864, -0.001048]
    hi[5] = [0.311884, -0.002032, -0.000310]

    return lo, hi
    #endif


def coeffs_rga_sp18_out() -> Tuple[np.ndarray, np.ndarray]:
    lo = np.zeros((6, 3), dtype=np.float64)
    hi = np.zeros((6, 3), dtype=np.float64)

    lo[0] = [0.181450, 0.007710, -0.000542]
    hi[0] = [0.297929, -0.003157, 0.000260]
    lo[1] = [0.189908, 0.004218, -0.000162]
    hi[1] = [0.295094, -0.001258, -0.000010]
    lo[2] = [0.184224, 0.008849, -0.000674]
    hi[2] = [0.297369, -0.000138, -0.000264]
    lo[3] = [0.181792, 0.008900, -0.000568]
    hi[3] = [0.294249, 0.001436, -0.000410]
    lo[4] = [0.183149, 0.006206, -0.000264]
    hi[4] = [0.301745, -0.004616, 0.000288]
    lo[5] = [0.170606, 0.013561, -0.001025]
    hi[5] = [0.294579, 0.001291, -0.000489]

    return lo, hi
    #endif


def coeffs_mc_run11() -> Tuple[np.ndarray, np.ndarray]:
    lo = np.zeros((6, 3), dtype=np.float64)
    hi = np.zeros((6, 3), dtype=np.float64)

    lo[0] = [0.182342, 0.010612, -0.000989]
    hi[0] = [0.316417, -0.007887, 0.000497]
    lo[1] = [0.174139, 0.015019, -0.001926]
    hi[1] = [0.319061, -0.009477, 0.000800]
    lo[2] = [0.168963, 0.017387, -0.002206]
    hi[2] = [0.318202, -0.008970, 0.000741]
    lo[3] = [0.176921, 0.012012, -0.001191]
    hi[3] = [0.315149, -0.008103, 0.000559]
    lo[4] = [0.177276, 0.013250, -0.001440]
    hi[4] = [0.316444, -0.008119, 0.000571]
    lo[5] = [0.180294, 0.011205, -0.001183]
    hi[5] = [0.317384, -0.008433, 0.000605]

    return lo, hi
    #endif


def coeffs_for_period(period_label: str, dataset: str) -> Tuple[np.ndarray, np.ndarray, bool]:
    if dataset == "mc":
        lo, hi = coeffs_mc_run11()
        return lo, hi, True
    #endif

    if period_label == "Fa18 Inb":
        lo, hi = coeffs_rga_fa18_inb()
        return lo, hi, True
    elif period_label == "Fa18 Out":
        lo, hi = coeffs_rga_fa18_out()
        return lo, hi, True
    elif period_label == "Sp19 Inb":
        lo, hi = coeffs_rga_sp19_inb()
        return lo, hi, True
    elif period_label == "Sp18 Inb":
        lo, hi = coeffs_rga_sp18_inb()
        return lo, hi, True
    elif period_label == "Sp18 Out":
        lo, hi = coeffs_rga_sp18_out()
        return lo, hi, True
    #endif

    # Java default: return sf > 0.19
    lo = np.zeros((6, 3), dtype=np.float64)
    hi = np.zeros((6, 3), dtype=np.float64)
    lo[:, 0] = SF_SIMPLE_MIN
    lo[:, 1] = 0.0
    lo[:, 2] = 0.0
    return lo, hi, False
    #endif


def eval_quad(coeffs: np.ndarray, p: np.ndarray) -> np.ndarray:
    a0 = coeffs[:, 0:1]
    a1 = coeffs[:, 1:2]
    a2 = coeffs[:, 2:3]
    return a0 + a1 * p[None, :] + a2 * (p[None, :] * p[None, :])
    #endif


# -----------------------------------------------------------------------------
# COMPUTATION
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class MeanProfile:
    p_centers: np.ndarray
    mean_sf: np.ndarray
    err_sf: np.ndarray
    n_in_bin: np.ndarray

    pass_fraction_poly: float
    pass_n_poly: int

    pass_fraction_simple: float
    pass_n_simple: int

    n_used: int


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
        fatal("Invalid p range.")
    #endif
    if args.pbins <= 0:
        fatal("Invalid --pbins.")
    #endif
    if args.max_workers <= 0 or args.max_workers > MAX_WORKERS_HARD:
        fatal(f"--max-workers must be in [1,{MAX_WORKERS_HARD}]")
    #endif
    if set(DATA_FILES.keys()) != set(MC_FILES.keys()):
        fatal("DATA_FILES and MC_FILES keys differ.")
    #endif

    for label, path in DATA_FILES.items():
        if not os.path.isfile(path):
            fatal(f"Missing DATA file '{label}': {path}")
        #endif
    #endfor
    for label, path in MC_FILES.items():
        if not os.path.isfile(path):
            fatal(f"Missing MC file '{label}': {path}")
        #endif
    #endfor


def compute_profile_for_file(args: Tuple[str, str, float, float, int, str]) -> Tuple[str, MeanProfile]:
    period_label, root_path, pmin, pmax, pbins, dataset = args

    try:
        with uproot.open(root_path) as f:
            if TREE_NAME not in f:
                raise KeyError(f"TTree '{TREE_NAME}' not found. Keys: {list(f.keys())}")
            #endif
            tree = f[TREE_NAME]

            required = ["particle_pid", "p", "cal_energy_1", "cal_energy_4", "cal_energy_7", "cal_sector"]
            for br in required:
                if br not in tree.keys():
                    raise KeyError(
                        f"Required branch '{br}' not found in tree '{TREE_NAME}'. "
                        "This script requires 'cal_sector' (1..6) to apply sector-dependent SF cuts."
                    )
                #endif
            #endfor

            arr = tree.arrays(required, library="np")
            pid = arr["particle_pid"]
            p = arr["p"]
            e1 = arr["cal_energy_1"]
            e4 = arr["cal_energy_4"]
            e7 = arr["cal_energy_7"]
            sector = arr["cal_sector"]

            base = (
                (pid == PID_ELECTRON)
                & np.isfinite(p)
                & np.isfinite(e1)
                & np.isfinite(e4)
                & np.isfinite(e7)
                & np.isfinite(sector)
                & (p > 0.0)
                & (e1 >= 0.0)
                & (e4 >= 0.0)
                & (e7 >= 0.0)
            )

            p_sel = p[base].astype(np.float64)
            sf_sel = (e1[base].astype(np.float64) + e4[base].astype(np.float64) + e7[base].astype(np.float64)) / p_sel
            sec_sel = sector[base].astype(np.int64)

            sec_ok = (sec_sel >= 1) & (sec_sel <= 6)
            p_sel = p_sel[sec_ok]
            sf_sel = sf_sel[sec_ok]
            sec_sel = sec_sel[sec_ok]

            ok = np.isfinite(sf_sel) & (sf_sel >= 0.0)
            p_sel = p_sel[ok]
            sf_sel = sf_sel[ok]
            sec_sel = sec_sel[ok]

            n_used = int(sf_sel.size)

            # Simple cut: sf >= 0.19
            pass_simple = (sf_sel >= SF_SIMPLE_MIN)
            pass_n_simple = int(np.count_nonzero(pass_simple))
            pass_fraction_simple = float(pass_n_simple) / float(n_used) if n_used > 0 else 0.0

            # Polynomial band cut (sector-dependent)
            lo_coeffs, hi_coeffs, has_upper = coeffs_for_period(period_label, dataset)
            idx = sec_sel - 1  # 0..5

            lo = lo_coeffs[idx, 0] + lo_coeffs[idx, 1] * p_sel + lo_coeffs[idx, 2] * p_sel * p_sel
            if has_upper:
                hi = hi_coeffs[idx, 0] + hi_coeffs[idx, 1] * p_sel + hi_coeffs[idx, 2] * p_sel * p_sel
                pass_poly = (sf_sel > lo) & (sf_sel < hi)
            else:
                pass_poly = (sf_sel > lo)
            #endif

            pass_n_poly = int(np.count_nonzero(pass_poly))
            pass_fraction_poly = float(pass_n_poly) / float(n_used) if n_used > 0 else 0.0

            # Mean(sf) vs p (NO additional cuts; it is the raw mean SF)
            edges = np.linspace(pmin, pmax, pbins + 1, dtype=np.float64)
            centers = 0.5 * (edges[:-1] + edges[1:])

            bin_idx = np.digitize(p_sel, edges) - 1
            inrange = (bin_idx >= 0) & (bin_idx < pbins)
            bin_idx = bin_idx[inrange]
            sf_use = sf_sel[inrange]

            n = np.zeros(pbins, dtype=np.int64)
            s1 = np.zeros(pbins, dtype=np.float64)
            s2 = np.zeros(pbins, dtype=np.float64)

            np.add.at(n, bin_idx, 1)
            np.add.at(s1, bin_idx, sf_use)
            np.add.at(s2, bin_idx, sf_use * sf_use)

            mean = np.full(pbins, np.nan, dtype=np.float64)
            err = np.full(pbins, np.nan, dtype=np.float64)

            nonzero = n > 0
            mean[nonzero] = s1[nonzero] / n[nonzero]

            gt1 = n > 1
            var = np.zeros(pbins, dtype=np.float64)
            var[gt1] = (s2[gt1] - (s1[gt1] * s1[gt1]) / n[gt1]) / (n[gt1] - 1.0)
            var[var < 0.0] = 0.0
            err[gt1] = np.sqrt(var[gt1] / n[gt1])

            prof = MeanProfile(
                p_centers=centers,
                mean_sf=mean,
                err_sf=err,
                n_in_bin=n,
                pass_fraction_poly=pass_fraction_poly,
                pass_n_poly=pass_n_poly,
                pass_fraction_simple=pass_fraction_simple,
                pass_n_simple=pass_n_simple,
                n_used=n_used,
            )
            return period_label, prof

    except Exception as e:
        tb = traceback.format_exc()
        raise RuntimeError(
            f"Failed processing file:\n"
            f"  period={period_label}\n"
            f"  path={root_path}\n"
            f"  dataset={dataset}\n"
            f"  error={repr(e)}\n\n{tb}"
        )


def run_parallel_profiles(file_map: Dict[str, str], pmin: float, pmax: float, pbins: int, dataset: str, max_workers: int) -> Dict[str, MeanProfile]:
    items = sorted(file_map.items(), key=lambda kv: kv[0])
    tasks = [(label, path, pmin, pmax, pbins, dataset) for (label, path) in items]

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
            fatal(f"Missing profile for '{period_label}' ({dataset})")
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
    fig, axes = plt.subplots(2, 3, figsize=(16, 9), sharex=True, sharey=True)
    fig.suptitle(title, fontsize=18)

    for r in range(2):
        for c in range(3):
            axes[r, c].axis("off")
        #endfor
    #endfor

    p_grid = np.linspace(pmin, pmax, 400, dtype=np.float64)

    for period_label, (r, c) in PANEL_POS.items():
        ax = axes[r, c]
        ax.axis("on")

        dp = data_prof[period_label]
        mp = mc_prof[period_label]

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
            label=f"data poly={100.0*dp.pass_fraction_poly:.2f}%  sf>=0.19={100.0*dp.pass_fraction_simple:.2f}%",
        )
        ax.errorbar(
            mp.p_centers[mmask],
            mp.mean_sf[mmask],
            yerr=mp.err_sf[mmask],
            fmt="o",
            color="red",
            markersize=3,
            linewidth=1.0,
            label=f"mc poly={100.0*mp.pass_fraction_poly:.2f}%  sf>=0.19={100.0*mp.pass_fraction_simple:.2f}%",
        )

        # Polynomial cut curves (representative median-over-sectors for DATA)
        lo_d, hi_d, has_upper_d = coeffs_for_period(period_label, "data")
        lo_all = eval_quad(lo_d, p_grid)
        lo_line = np.median(lo_all, axis=0)

        hi_line = None
        if has_upper_d:
            hi_all = eval_quad(hi_d, p_grid)
            hi_line = np.median(hi_all, axis=0)
        #endif

        ax.plot(p_grid, lo_line, linestyle="--", color="blue", linewidth=1.2, label="SF cut band + sf>=0.19")
        if hi_line is not None:
            ax.plot(p_grid, hi_line, linestyle="--", color="blue", linewidth=1.2)
        #endif

        # Simple cut line (sf >= 0.19)
        ax.axhline(SF_SIMPLE_MIN, linestyle="--", color="blue", linewidth=1.0)

        ax.set_title(period_label, fontsize=13)
        ax.set_xlim(pmin, pmax)
        ax.set_ylim(Y_MIN, Y_MAX)
        ax.set_xlabel("p (GeV)", fontsize=12)
        ax.set_ylabel("mean sampling fraction", fontsize=12)
        ax.grid(True, alpha=0.25)
        ax.legend(loc="upper right", fontsize=8, frameon=True)

    #endfor

    axes[0, 2].axis("off")
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.95])
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Mean sampling fraction vs p (data vs mc), electrons only, with SF cuts.")
    p.add_argument("--pmin", type=float, default=P_MIN_DEFAULT, help="Minimum p (GeV).")
    p.add_argument("--pmax", type=float, default=P_MAX_DEFAULT, help="Maximum p (GeV).")
    p.add_argument("--pbins", type=int, default=60, help="Number of p bins.")
    p.add_argument("--max-workers", type=int, default=5, help="Max parallel workers (hard limit 5).")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    validate_inputs(args)

    if not os.path.isdir(OUTDIR):
        os.makedirs(OUTDIR, exist_ok=True)
    #endif

    data_prof = run_parallel_profiles(DATA_FILES, args.pmin, args.pmax, args.pbins, "data", args.max_workers)
    mc_prof = run_parallel_profiles(MC_FILES, args.pmin, args.pmax, args.pbins, "mc", args.max_workers)

    outpath = os.path.join(OUTDIR, "sampling_fraction_vs_p_electron.png")
    title = "Mean sampling fraction vs p: electrons (pid=11), (cal_energy_1+4+7)/p"
    plot_2x3_canvas(title, data_prof, mc_prof, outpath, args.pmin, args.pmax)

    print("")
    print("Sampling fraction survival rates (denominator = n_used after basic validity + sector 1..6):")
    periods = sorted(PANEL_POS.keys(), key=lambda k: (PANEL_POS[k][0], PANEL_POS[k][1]))
    for period_label in periods:
        dp = data_prof[period_label]
        mp = mc_prof[period_label]
        print(
            f"{period_label}: "
            f"data poly={100.0*dp.pass_fraction_poly:.3f}% (N={dp.pass_n_poly}/{dp.n_used})  "
            f"data sf>=0.19={100.0*dp.pass_fraction_simple:.3f}% (N={dp.pass_n_simple}/{dp.n_used})  ||  "
            f"mc poly={100.0*mp.pass_fraction_poly:.3f}% (N={mp.pass_n_poly}/{mp.n_used})  "
            f"mc sf>=0.19={100.0*mp.pass_fraction_simple:.3f}% (N={mp.pass_n_simple}/{mp.n_used})"
        )
    #endfor

    print(f"Wrote: {outpath}")
    print("Done.")


if __name__ == "__main__":
    main()
#endif