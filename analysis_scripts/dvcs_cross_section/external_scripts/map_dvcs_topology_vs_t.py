#!/usr/bin/env python3
"""
map_dvcs_topology_vs_t.py

Standalone event-level diagnostic for the CLAS12 pass-2 DVCS analysis.

Purpose
-------
Map the strong exclusive-kinematics correlation between |t|, recoil-proton
kinematics, and real-photon kinematics before assigning the low-|t| pass-2
deficit to either proton or photon reconstruction.

For each run period and for the combined sample, the script reports versus |t|:

  * proton region fractions: FD / CD
  * photon region fractions: FT / FD
  * joint proton-photon topology fractions:
        FDp-FDg, FDp-FTg, CDp-FDg, CDp-FTg
  * distributions and summary statistics for:
        proton p, theta, phi
        photon p, theta, phi

The script reads the same pass-2 DVCS DATA ROOT trees listed in load_trees.cpp.
It deliberately does NOT apply the Krishna Neupane efficiency correction or
any photon-efficiency correction. This is a population/kinematic map.

Selection
---------
By default it applies only the analysis-wide baseline cuts that can be
reproduced unambiguously from the event tree:

    |t| < 1 GeV^2
    open_angle_ep2 > 5 deg

It does NOT apply the obsolete special Sp18-Out sector-quality exclusions.

It also does NOT reproduce the topology/period-dependent optimized exclusivity
quantile cuts from the full C++ analysis. This is intentional: the first
diagnostic asks where the candidate DVCS population lives in detector and
kinematic space. Use --no-baseline-cuts if an entirely preselection-level map
is desired.

Requirements
------------
Python 3, numpy, pandas, matplotlib, uproot

Example
-------
python3 external_scripts/map_dvcs_topology_vs_t.py

python3 external_scripts/map_dvcs_topology_vs_t.py \
    --output-dir output/topology_vs_t_map
"""

import argparse
import math
import os
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", f"/tmp/matplotlib-{os.getuid()}")
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)

import matplotlib
matplotlib.use("Agg", force=True)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot


TREE_NAME = "PhysicsEvents"

DEFAULT_FILES = {
    "Sp18 Inb": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp18_inb_epgamma.root",
    "Sp18 Out": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp18_out_epgamma.root",
    "Fa18 Inb": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_inb_epgamma.root",
    "Fa18 Out": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_out_epgamma.root",
    "Sp19 Inb": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp19_inb_epgamma.root",
}

PERIOD_ORDER = ["Sp18 Inb", "Sp18 Out", "Fa18 Inb", "Fa18 Out", "Sp19 Inb"]

REQUIRED_BRANCHES = [
    "t1",
    "detector1",
    "detector2",
    "p1_p",
    "p1_theta",
    "p1_phi",
    "p2_p",
    "p2_theta",
    "p2_phi",
]

OPTIONAL_BRANCHES = [
    "open_angle_ep2",
    "e_phi",
    "runnum",
]

RAD2DEG = 180.0 / math.pi

T_EDGES_DEFAULT = [0.0, 0.20, 0.30, 0.40, 1.00]
T_LABELS_DEFAULT = [
    "|t| < 0.20",
    "0.20 <= |t| < 0.30",
    "0.30 <= |t| < 0.40",
    "0.40 <= |t| < 1.00",
]

VARIABLES = {
    "p_p": {
        "label": r"$p_p$ (GeV)",
        "title": "proton momentum",
        "bins": np.linspace(0.0, 3.5, 71),
    },
    "p_theta_deg": {
        "label": r"$\theta_p$ (deg)",
        "title": "proton polar angle",
        "bins": np.linspace(0.0, 80.0, 81),
    },
    "p_phi_deg": {
        "label": r"$\phi_p$ (deg)",
        "title": "proton azimuth",
        "bins": np.linspace(0.0, 360.0, 73),
    },
    "g_p": {
        "label": r"$p_\gamma$ (GeV)",
        "title": "photon momentum",
        "bins": np.linspace(0.0, 10.5, 71),
    },
    "g_theta_deg": {
        "label": r"$\theta_\gamma$ (deg)",
        "title": "photon polar angle",
        "bins": np.linspace(0.0, 45.0, 91),
    },
    "g_phi_deg": {
        "label": r"$\phi_\gamma$ (deg)",
        "title": "photon azimuth",
        "bins": np.linspace(0.0, 360.0, 73),
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Map DVCS proton/photon detector topology and kinematics versus |t|."
    )
    parser.add_argument(
        "--output-dir",
        default="output/topology_vs_t_map",
        help="Output directory.",
    )
    parser.add_argument(
        "--period",
        action="append",
        choices=PERIOD_ORDER,
        help="Run only selected period(s). May be supplied more than once.",
    )
    parser.add_argument(
        "--t-edges",
        nargs="+",
        type=float,
        default=T_EDGES_DEFAULT,
        help="|t| bin edges in GeV^2. Default: 0 0.20 0.30 0.40 1.00",
    )
    parser.add_argument(
        "--step-size",
        default="250 MB",
        help="uproot iterate step size. Default: 250 MB",
    )
    parser.add_argument(
        "--no-baseline-cuts",
        action="store_true",
        help="Do not apply |t|<1 and open_angle_ep2>5 baseline cuts.",
    )
    parser.add_argument(
        "--max-events",
        type=int,
        default=None,
        help="Optional per-period event cap for quick tests.",
    )
    return parser.parse_args()


def wrap_deg(phi_rad: np.ndarray) -> np.ndarray:
    return np.mod(np.asarray(phi_rad, dtype=float) * RAD2DEG, 360.0)


def fd_sector(phi_rad: np.ndarray) -> np.ndarray:
    """
    Match the six 60-degree CLAS12 FD sectors:
      S1: 330-30, S2: 30-90, ..., S6: 270-330.
    """
    phi = wrap_deg(phi_rad)
    return ((np.floor((phi + 30.0) / 60.0).astype(int)) % 6) + 1


def t_bin_labels(edges: List[float]) -> List[str]:
    labels = []
    for i in range(len(edges) - 1):
        lo, hi = edges[i], edges[i + 1]
        if i == 0 and abs(lo) < 1e-12:
            labels.append(f"|t| < {hi:.2f}")
        else:
            labels.append(f"{lo:.2f} <= |t| < {hi:.2f}")
    #endfor
    return labels


def topology_label(det1: np.ndarray, det2: np.ndarray) -> np.ndarray:
    out = np.full(len(det1), "other", dtype=object)
    out[(det1 == 1) & (det2 == 1)] = "FDp-FDg"
    out[(det1 == 1) & (det2 == 0)] = "FDp-FTg"
    out[(det1 == 2) & (det2 == 1)] = "CDp-FDg"
    out[(det1 == 2) & (det2 == 0)] = "CDp-FTg"
    return out


def validate_tree(path: str) -> List[str]:
    if not Path(path).exists():
        raise FileNotFoundError(path)

    with uproot.open(path) as root_file:
        if TREE_NAME not in root_file:
            raise RuntimeError(f"{path}: missing tree '{TREE_NAME}'")
        tree = root_file[TREE_NAME]
        available = set(tree.keys())
    #endif

    missing = [name for name in REQUIRED_BRANCHES if name not in available]
    if missing:
        raise RuntimeError(f"{path}: missing required branches: {', '.join(missing)}")
    #endif

    return [name for name in OPTIONAL_BRANCHES if name in available]


def apply_selection(
    arrays: Dict[str, np.ndarray],
    period: str,
    apply_baseline: bool,
) -> np.ndarray:
    n = len(arrays["t1"])
    keep = np.ones(n, dtype=bool)

    tabs = np.abs(np.asarray(arrays["t1"], dtype=float))
    keep &= np.isfinite(tabs)

    if apply_baseline:
        keep &= tabs < 1.0

        if "open_angle_ep2" in arrays:
            open_angle = np.asarray(arrays["open_angle_ep2"], dtype=float)
            keep &= np.isfinite(open_angle) & (open_angle > 5.0)
        #endif
    #endif

    return keep


def chunk_to_frame(
    arrays: Dict[str, np.ndarray],
    period: str,
    keep: np.ndarray,
    edges: List[float],
    labels: List[str],
) -> pd.DataFrame:
    det1 = np.asarray(arrays["detector1"], dtype=int)[keep]
    det2 = np.asarray(arrays["detector2"], dtype=int)[keep]
    tabs = np.abs(np.asarray(arrays["t1"], dtype=float)[keep])

    p_p = np.asarray(arrays["p1_p"], dtype=float)[keep]
    p_th = np.asarray(arrays["p1_theta"], dtype=float)[keep] * RAD2DEG
    p_ph = wrap_deg(np.asarray(arrays["p1_phi"], dtype=float)[keep])

    g_p = np.asarray(arrays["p2_p"], dtype=float)[keep]
    g_th = np.asarray(arrays["p2_theta"], dtype=float)[keep] * RAD2DEG
    g_ph = wrap_deg(np.asarray(arrays["p2_phi"], dtype=float)[keep])

    ibin = np.searchsorted(np.asarray(edges), tabs, side="right") - 1
    valid_tbin = (ibin >= 0) & (ibin < len(labels))

    frame = pd.DataFrame({
        "period": period,
        "tabs": tabs[valid_tbin],
        "t_bin": np.asarray(labels, dtype=object)[ibin[valid_tbin]],
        "detector1": det1[valid_tbin],
        "detector2": det2[valid_tbin],
        "proton_region": np.where(det1[valid_tbin] == 1, "FD",
                           np.where(det1[valid_tbin] == 2, "CD", "other")),
        "photon_region": np.where(det2[valid_tbin] == 0, "FT",
                           np.where(det2[valid_tbin] == 1, "FD", "other")),
        "topology": topology_label(det1[valid_tbin], det2[valid_tbin]),
        "p_p": p_p[valid_tbin],
        "p_theta_deg": p_th[valid_tbin],
        "p_phi_deg": p_ph[valid_tbin],
        "g_p": g_p[valid_tbin],
        "g_theta_deg": g_th[valid_tbin],
        "g_phi_deg": g_ph[valid_tbin],
    })
    return frame


def load_period(
    period: str,
    path: str,
    edges: List[float],
    labels: List[str],
    step_size: str,
    apply_baseline: bool,
    max_events: int | None,
) -> pd.DataFrame:
    optional = validate_tree(path)
    branches = REQUIRED_BRANCHES + optional

    frames = []
    n_read = 0
    n_kept = 0

    print(f"[map] {period}: {path}")

    source = f"{path}:{TREE_NAME}"
    for arrays_ak in uproot.iterate(
        source,
        expressions=branches,
        step_size=step_size,
        library="np",
    ):
        arrays = dict(arrays_ak)

        if max_events is not None:
            remaining = max_events - n_read
            if remaining <= 0:
                break
            #endif

            n_chunk = len(arrays["t1"])
            if n_chunk > remaining:
                arrays = {key: value[:remaining] for key, value in arrays.items()}
            #endif
        #endif

        n_here = len(arrays["t1"])
        n_read += n_here

        keep = apply_selection(
            arrays,
            period,
            apply_baseline=apply_baseline,
        )
        n_kept += int(np.count_nonzero(keep))

        frame = chunk_to_frame(arrays, period, keep, edges, labels)
        if not frame.empty:
            frames.append(frame)
        #endif

        if max_events is not None and n_read >= max_events:
            break
        #endif
    #endfor

    print(f"[map] {period}: read {n_read:,}; selected {n_kept:,}")

    if not frames:
        return pd.DataFrame()
    #endif

    return pd.concat(frames, ignore_index=True)


def fraction_table(
    df: pd.DataFrame,
    labels: List[str],
    category_col: str,
    categories: List[str],
) -> pd.DataFrame:
    rows = []
    samples = [(period, df[df["period"] == period]) for period in PERIOD_ORDER]
    samples.append(("Combined", df))

    for sample_name, sample in samples:
        if sample.empty:
            continue
        #endif

        for t_label in labels:
            sub = sample[sample["t_bin"] == t_label]
            n = len(sub)
            row = {
                "sample": sample_name,
                "t_bin": t_label,
                "N": n,
                "mean_abs_t_GeV2": sub["tabs"].mean() if n else np.nan,
            }

            for category in categories:
                count = int(np.count_nonzero(sub[category_col].to_numpy() == category))
                row[f"N_{category}"] = count
                row[f"fraction_{category}"] = count / n if n else np.nan
            #endfor

            rows.append(row)
        #endfor
    #endfor

    return pd.DataFrame(rows)


def kinematic_summary(df: pd.DataFrame, labels: List[str]) -> pd.DataFrame:
    rows = []
    samples = [(period, df[df["period"] == period]) for period in PERIOD_ORDER]
    samples.append(("Combined", df))

    for sample_name, sample in samples:
        if sample.empty:
            continue
        #endif

        for t_label in labels:
            sub = sample[sample["t_bin"] == t_label]
            if sub.empty:
                continue
            #endif

            for variable in VARIABLES:
                values = sub[variable].to_numpy(dtype=float)
                values = values[np.isfinite(values)]
                if len(values) == 0:
                    continue
                #endif

                q16, q50, q84 = np.percentile(values, [16.0, 50.0, 84.0])
                rows.append({
                    "sample": sample_name,
                    "t_bin": t_label,
                    "variable": variable,
                    "N": len(values),
                    "mean": float(np.mean(values)),
                    "median": float(q50),
                    "q16": float(q16),
                    "q84": float(q84),
                })
            #endfor
        #endfor
    #endfor

    return pd.DataFrame(rows)


def plot_fraction_vs_t(
    table: pd.DataFrame,
    edges: List[float],
    categories: List[str],
    output: Path,
    title: str,
) -> None:
    combined = table[table["sample"] == "Combined"].copy()
    if combined.empty:
        return
    #endif

    x = combined["mean_abs_t_GeV2"].to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(8.5, 6.0))
    for category in categories:
        y = combined[f"fraction_{category}"].to_numpy(dtype=float)
        ax.plot(x, y, marker="o", linewidth=1.5, label=category)
    #endfor

    for edge in edges[1:-1]:
        ax.axvline(edge, linewidth=0.8, alpha=0.35)
    #endfor

    ax.set_xlabel(r"$|t|$ (GeV$^2$)")
    ax.set_ylabel("event fraction")
    ax.set_ylim(0.0, 1.0)
    ax.set_title(title)
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_period_topology_fractions(
    table: pd.DataFrame,
    edges: List[float],
    categories: List[str],
    output: Path,
) -> None:
    periods = [p for p in PERIOD_ORDER if p in set(table["sample"])]
    if not periods:
        return
    #endif

    fig, axes = plt.subplots(
        len(periods),
        1,
        figsize=(8.5, 2.5 * len(periods)),
        sharex=True,
        sharey=True,
    )
    if len(periods) == 1:
        axes = [axes]
    #endif

    for ax, period in zip(axes, periods):
        sub = table[table["sample"] == period]
        x = sub["mean_abs_t_GeV2"].to_numpy(dtype=float)

        for category in categories:
            y = sub[f"fraction_{category}"].to_numpy(dtype=float)
            ax.plot(x, y, marker="o", linewidth=1.2, label=category)
        #endfor

        for edge in edges[1:-1]:
            ax.axvline(edge, linewidth=0.7, alpha=0.25)
        #endfor

        ax.set_ylabel("fraction")
        ax.set_ylim(0.0, 1.0)
        ax.set_title(period)
        ax.grid(alpha=0.20)
    #endfor

    axes[0].legend(ncol=4, frameon=False, fontsize=9)
    axes[-1].set_xlabel(r"$|t|$ (GeV$^2$)")
    fig.suptitle("Joint proton-photon topology fractions versus |t|", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_variable_distributions(
    df: pd.DataFrame,
    labels: List[str],
    variable: str,
    output: Path,
) -> None:
    cfg = VARIABLES[variable]
    values_all = df[variable].to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(8.5, 6.0))

    for t_label in labels:
        values = df.loc[df["t_bin"] == t_label, variable].to_numpy(dtype=float)
        values = values[np.isfinite(values)]
        if len(values) == 0:
            continue
        #endif

        hist, bins = np.histogram(values, bins=cfg["bins"], density=False)
        width = np.diff(bins)
        area = np.sum(hist * width)
        if area <= 0:
            continue
        #endif

        density = hist / area
        centers = 0.5 * (bins[:-1] + bins[1:])
        ax.step(centers, density, where="mid", linewidth=1.4, label=f"{t_label} (N={len(values):,})")
    #endfor

    ax.set_xlabel(cfg["label"])
    ax.set_ylabel("normalized event density")
    ax.set_title(f"{cfg['title']} by |t| interval — combined periods")
    ax.grid(alpha=0.20)
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_joint_theta_map(df: pd.DataFrame, output: Path) -> None:
    finite = (
        np.isfinite(df["p_theta_deg"].to_numpy(dtype=float))
        & np.isfinite(df["g_theta_deg"].to_numpy(dtype=float))
    )
    sub = df.loc[finite]
    if sub.empty:
        return
    #endif

    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    h = ax.hist2d(
        sub["p_theta_deg"],
        sub["g_theta_deg"],
        bins=[np.linspace(0, 80, 81), np.linspace(0, 45, 91)],
        norm=matplotlib.colors.LogNorm(),
    )
    fig.colorbar(h[3], ax=ax, label="events")
    ax.axvline(37.0, linewidth=1.0, linestyle="--", label=r"Neupane FD/CD boundary ($37^\circ$)")
    ax.axhline(5.5, linewidth=1.0, linestyle="--", label=r"FT/FD photon boundary ($5.5^\circ$)")
    ax.set_xlabel(r"$\theta_p$ (deg)")
    ax.set_ylabel(r"$\theta_\gamma$ (deg)")
    ax.set_title(r"Exclusive proton-photon angular correlation")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def write_readme(
    output_dir: Path,
    args: argparse.Namespace,
    periods: List[str],
    labels: List[str],
) -> None:
    text = f"""DVCS topology-versus-|t| diagnostic

Periods:
  {", ".join(periods)}

|t| intervals:
  {", ".join(labels)}

Selection:
  baseline cuts: {"OFF" if args.no_baseline_cuts else "ON: |t| < 1 GeV^2 and open_angle_ep2 > 5 deg"}
  obsolete Sp18-Out special sector-quality cuts: NOT APPLIED

Important:
  This script maps event populations and reconstructed kinematics only.
  It does not apply the Krishna Neupane proton-efficiency correction.
  It does not apply a photon-efficiency correction.
  It does not reproduce the full topology/period-dependent optimized
  exclusivity-cut machinery from the C++ production extraction.

CSV outputs:
  proton_region_fractions_vs_t.csv
  photon_region_fractions_vs_t.csv
  joint_topology_fractions_vs_t.csv
  kinematic_summary_vs_t.csv

PNG outputs:
  proton_region_fractions_vs_t.png
  photon_region_fractions_vs_t.png
  joint_topology_fractions_vs_t.png
  joint_topology_fractions_by_period.png
  proton_momentum_by_t.png
  proton_theta_by_t.png
  proton_phi_by_t.png
  photon_momentum_by_t.png
  photon_theta_by_t.png
  photon_phi_by_t.png
  proton_photon_theta_correlation.png
"""
    (output_dir / "README.txt").write_text(text)


def main() -> None:
    args = parse_args()

    edges = list(args.t_edges)
    if len(edges) < 2 or any(b <= a for a, b in zip(edges[:-1], edges[1:])):
        raise ValueError("--t-edges must be strictly increasing.")
    #endif

    labels = t_bin_labels(edges)
    periods = args.period if args.period else PERIOD_ORDER

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    frames = []
    for period in periods:
        frame = load_period(
            period=period,
            path=DEFAULT_FILES[period],
            edges=edges,
            labels=labels,
            step_size=args.step_size,
            apply_baseline=not args.no_baseline_cuts,
            max_events=args.max_events,
        )
        if not frame.empty:
            frames.append(frame)
        #endif
    #endfor

    if not frames:
        raise RuntimeError("No events survived / no input frames were produced.")
    #endif

    df = pd.concat(frames, ignore_index=True)
    print(f"[map] total selected events in requested |t| range: {len(df):,}")

    proton_table = fraction_table(
        df, labels, "proton_region", ["FD", "CD", "other"]
    )
    photon_table = fraction_table(
        df, labels, "photon_region", ["FT", "FD", "other"]
    )
    topology_categories = ["FDp-FDg", "FDp-FTg", "CDp-FDg", "CDp-FTg", "other"]
    topology_table = fraction_table(
        df, labels, "topology", topology_categories
    )
    kin_table = kinematic_summary(df, labels)

    proton_table.to_csv(output_dir / "proton_region_fractions_vs_t.csv", index=False)
    photon_table.to_csv(output_dir / "photon_region_fractions_vs_t.csv", index=False)
    topology_table.to_csv(output_dir / "joint_topology_fractions_vs_t.csv", index=False)
    kin_table.to_csv(output_dir / "kinematic_summary_vs_t.csv", index=False)

    plot_fraction_vs_t(
        proton_table,
        edges,
        ["FD", "CD"],
        output_dir / "proton_region_fractions_vs_t.png",
        "Proton detector-region fractions versus |t|",
    )
    plot_fraction_vs_t(
        photon_table,
        edges,
        ["FT", "FD"],
        output_dir / "photon_region_fractions_vs_t.png",
        "Photon detector-region fractions versus |t|",
    )
    plot_fraction_vs_t(
        topology_table,
        edges,
        topology_categories[:-1],
        output_dir / "joint_topology_fractions_vs_t.png",
        "Joint proton-photon topology fractions versus |t|",
    )
    plot_period_topology_fractions(
        topology_table,
        edges,
        topology_categories[:-1],
        output_dir / "joint_topology_fractions_by_period.png",
    )

    output_names = {
        "p_p": "proton_momentum_by_t.png",
        "p_theta_deg": "proton_theta_by_t.png",
        "p_phi_deg": "proton_phi_by_t.png",
        "g_p": "photon_momentum_by_t.png",
        "g_theta_deg": "photon_theta_by_t.png",
        "g_phi_deg": "photon_phi_by_t.png",
    }

    for variable, filename in output_names.items():
        plot_variable_distributions(
            df,
            labels,
            variable,
            output_dir / filename,
        )
    #endfor

    plot_joint_theta_map(
        df,
        output_dir / "proton_photon_theta_correlation.png",
    )

    write_readme(output_dir, args, periods, labels)

    print(f"[map] wrote outputs to {output_dir}")
    print("[map] key first-look file: joint_topology_fractions_vs_t.csv")
    print("[map] key first-look plot: joint_topology_fractions_by_period.png")


if __name__ == "__main__":
    main()
