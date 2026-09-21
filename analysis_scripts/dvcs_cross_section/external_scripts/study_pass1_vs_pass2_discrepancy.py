#!/usr/bin/env python3
"""
study_pass1_vs_pass2_discrepancy.py

Standalone pass-1/pass-2 discrepancy and detector-topology diagnostic for CLAS12 DVCS.

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
python3 external_scripts/study_pass1_vs_pass2_discrepancy.py

python3 external_scripts/study_pass1_vs_pass2_discrepancy.py
"""

import argparse
import math
import os
from pathlib import Path
from typing import Dict, Iterable, List, Tuple
import re

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
        description="Study the pass-1/pass-2 DVCS discrepancy and map detector topology versus |t|."
    )
    parser.add_argument(
        "--output-dir",
        default=str(Path(__file__).resolve().parent / "output" / "study_pass1_vs_pass2_discrepancy"),
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
    parser.add_argument("--csv-dir", default=str(Path(__file__).resolve().parent.parent / "output" / "csvs"), help="Directory containing existing pass-2 CSV outputs.")
    parser.add_argument("--imports-dir", default=str(Path(__file__).resolve().parent.parent / "imports"), help="Directory containing pass-1 imports.")
    parser.add_argument("--pass1-authoritative", default=str(Path(__file__).resolve().parent.parent / "imports" / "clasdb_E214M1.txt"), help="Released pass-1 cross-section table.")
    parser.add_argument("--inclusive-csv", default=None, help="Explicit inclusive pass-2 CSV; otherwise auto-discover.")
    parser.add_argument("--cd-fd-csv", default=None, help="Explicit CD-FD pass-2 CSV; otherwise auto-discover.")
    parser.add_argument("--cd-ft-csv", default=None, help="Explicit CD-FT pass-2 CSV; otherwise auto-discover.")
    parser.add_argument("--fd-fd-csv", default=None, help="Explicit FD-FD pass-2 CSV; otherwise auto-discover.")
    parser.add_argument("--match-dx", type=float, default=0.015, help="Maximum |Delta xB| for released pass-1 matching.")
    parser.add_argument("--match-dq2", type=float, default=0.35, help="Maximum |Delta Q2| (GeV^2) for released pass-1 matching.")
    parser.add_argument("--match-dt", type=float, default=0.03, help="Maximum |Delta |t|| (GeV^2) for released pass-1 matching.")
    parser.add_argument("--match-dphi", type=float, default=8.0, help="Maximum circular |Delta phi| (deg) for released pass-1 matching.")
    parser.add_argument("--min-bin-events", type=float, default=50.0,
                        help="Require reconstructed DVCS-MC yield strictly greater than this value for a bin to enter the primary comparison; DATA yield is unrestricted. Default: 50.")
    parser.add_argument("--valerii-map", default=None,
                        help="Optional FD Valerii reproduction CSV (valerii_fd_results.csv; legacy cell-validity schema also accepted). Invalid/out-of-range cells receive no correction and are reported as uncovered.")
    parser.add_argument("--skip-event-map", action="store_true", help="Skip the ROOT-tree topology/kinematic map and run only the existing-CSV comparison.")
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



# -----------------------------------------------------------------------------
# Pass-1 versus pass-2 cross-section comparison
# -----------------------------------------------------------------------------

ANALYSIS_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CSV_DIR = ANALYSIS_ROOT / "output" / "csvs"
DEFAULT_IMPORTS_DIR = ANALYSIS_ROOT / "imports"


def parse_numeric(series: pd.Series) -> pd.Series:
    """Parse both ordinary numbers and the analysis CSV object-string format."""
    s = series.astype("string").str.strip()
    first = s.str.extract(
        r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)",
        expand=False,
    )
    return pd.to_numeric(first, errors="coerce")


def normalized_name(name: str) -> str:
    return "".join(ch.lower() for ch in str(name) if ch.isalnum())


def find_column(df: pd.DataFrame, candidates: List[str], contains: List[str] | None = None) -> str:
    by_norm = {normalized_name(c): c for c in df.columns}
    for candidate in candidates:
        key = normalized_name(candidate)
        if key in by_norm:
            return by_norm[key]
        #endif
    #endfor

    if contains:
        tokens = [normalized_name(x) for x in contains]
        matches = [
            c for c in df.columns
            if all(token in normalized_name(c) for token in tokens)
        ]
        if len(matches) == 1:
            return matches[0]
        #endif
        if len(matches) > 1:
            raise RuntimeError(f"Ambiguous column search {contains}: {matches[:12]}")
        #endif
    #endif

    raise RuntimeError(
        "Could not identify required column. Tried: " + ", ".join(candidates)
    )


def find_kinematic_columns(df: pd.DataFrame) -> Dict[str, str]:
    # Canonical combined-pass2 kinematic columns. Do not guess among the
    # period-specific averages: the cross section used below is the combined
    # 10.6-GeV result, so use the matching combined bin averages.
    columns = {
        "xB": "xBavg, 10.6 GeV",
        "Q2": "Q2avg, 10.6 GeV",
        "t": "t_abs_avg, 10.6 GeV",
        "phi": "phiavg, 10.6 GeV",
    }
    missing = [name for name in columns.values() if name not in df.columns]
    if missing:
        raise RuntimeError(
            "Missing expected canonical pass-2 kinematic column(s): "
            + ", ".join(missing)
        )
    #endif
    return columns


def combined_xs_column(df: pd.DataFrame) -> str:
    column = "cross sections, ep->epg, exp, 10.6 GeV, unpol"
    if column not in df.columns:
        raise RuntimeError(f"Missing expected canonical pass-2 cross-section column: {column}")
    #endif
    return column


PASS2_PERIODS = ["Fa18 Inb", "Fa18 Out", "Sp18 Inb", "Sp18 Out", "Sp19 Inb"]
STAT_THRESHOLDS = [0, 5, 10, 20, 30, 50, 100, 200, 500]


def sum_exact_period_columns(raw: pd.DataFrame, template: str, quantity: str) -> pd.Series:
    """Sum an exact per-period analysis quantity, preserving NaN when no period contributes."""
    values = []
    missing = []
    for period in PASS2_PERIODS:
        column = template.format(period=period)
        if column not in raw.columns:
            missing.append(column)
            continue
        #endif
        values.append(parse_numeric(raw[column]))
    #endfor

    if missing:
        raise RuntimeError(
            f"Missing expected pass-2 {quantity} column(s) in statistics study: "
            + "; ".join(missing)
        )
    #endif
    frame = pd.concat(values, axis=1)
    return frame.sum(axis=1, min_count=1)


def extract_pass2_bin_statistics(raw: pd.DataFrame) -> pd.DataFrame:
    """Extract the DATA and MC populations entering each combined pass-2 bin."""
    n_data = sum_exact_period_columns(
        raw,
        "signal yield, ep->epg, exp, {period}, unpol",
        "signal-yield",
    )
    n_rec = sum_exact_period_columns(
        raw,
        "reconstructed yield, ep->epg, mc, {period}",
        "reconstructed-MC",
    )
    n_gen = sum_exact_period_columns(
        raw,
        "generated yield, ep->epg, mc, {period}",
        "generated-MC",
    )

    out = pd.DataFrame({
        "N_data_signal": n_data,
        "N_mc_rec": n_rec,
        "N_mc_gen": n_gen,
    })
    out["acceptance_counts"] = out["N_mc_rec"] / out["N_mc_gen"]
    # Binomial counting approximation.  This is a diagnostic of MC-counting
    # stability, not a replacement for the production acceptance uncertainty.
    a = out["acceptance_counts"].to_numpy(float)
    ngen = out["N_mc_gen"].to_numpy(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        sigma_a = np.sqrt(np.clip(a * (1.0 - a), 0.0, None) / ngen)
        frac = sigma_a / a
    #endwith
    out["mc_acceptance_frac_stat"] = frac
    return out


def read_pass2_csv(path: Path, label: str) -> pd.DataFrame:
    raw = pd.read_csv(path, low_memory=False)
    kin = find_kinematic_columns(raw)
    xs_col = combined_xs_column(raw)

    stats = extract_pass2_bin_statistics(raw)
    out = pd.DataFrame({
        "xB": parse_numeric(raw[kin["xB"]]),
        "Q2": parse_numeric(raw[kin["Q2"]]),
        "t": parse_numeric(raw[kin["t"]]).abs(),
        "phi": np.mod(parse_numeric(raw[kin["phi"]]), 360.0),
        "xs_pass2": parse_numeric(raw[xs_col]),
        "N_data_signal": stats["N_data_signal"],
        "N_mc_rec": stats["N_mc_rec"],
        "N_mc_gen": stats["N_mc_gen"],
        "acceptance_counts": stats["acceptance_counts"],
        "mc_acceptance_frac_stat": stats["mc_acceptance_frac_stat"],
    })
    out["topology"] = label
    out["source_file"] = str(path)
    # Production-quality acceptance requirement: reconstructed MC strictly > 50.
    # DATA yield is deliberately unrestricted; its statistical uncertainty carries
    # the information from low-statistics DATA bins.
    out["valid_statistics_50"] = (
        np.isfinite(out["N_mc_rec"]) & (out["N_mc_rec"] > 50.0)
    )
    out = out.replace([np.inf, -np.inf], np.nan)
    out = out.dropna(subset=["xB", "Q2", "t", "phi", "xs_pass2"])
    out = out[out["xs_pass2"] > 0].reset_index(drop=True)
    print(f"[pass1/pass2] {label}: loaded {len(out):,} positive pass-2 cross sections from {path.name}")
    return out


def discover_csv(csv_dir: Path, include_tokens: List[str], exclude_tokens: List[str] | None = None) -> Path | None:
    files = sorted(csv_dir.glob("*.csv"))
    include = [normalized_name(x) for x in include_tokens]
    exclude = [normalized_name(x) for x in (exclude_tokens or [])]
    matches = []
    for path in files:
        key = normalized_name(path.stem)
        if all(token in key for token in include) and not any(token in key for token in exclude):
            matches.append(path)
        #endif
    #endfor

    if len(matches) == 1:
        return matches[0]
    #endif
    if len(matches) > 1:
        print(f"[pass1/pass2] WARNING: multiple files match {include_tokens}: {[p.name for p in matches]}")
        return matches[-1]
    #endif
    return None


def resolve_pass2_files(args, csv_dir: Path) -> Dict[str, Path]:
    """Resolve the four parent pass-2 CSVs using their exact production filenames."""
    exact = {
        "Inclusive": csv_dir / "dvcs_pass2_INCLUSIVE.csv",
        "CD-FD": csv_dir / "dvcs_pass2_CD-FD.csv",
        "CD-FT": csv_dir / "dvcs_pass2_CD-FT.csv",
        "FD-FD": csv_dir / "dvcs_pass2_FD-FD.csv",
    }

    resolved = {}
    for label, path in exact.items():
        if not path.exists():
            raise FileNotFoundError(
                f"Required {label} CSV not found at exact path: {path}"
            )
        #endif
        resolved[label] = path
        print(f"[pass1/pass2] {label}: {path.name}")
    #endfor

    return resolved


def read_authoritative_pass1(path: Path) -> pd.DataFrame:
    """
    Read imports/clasdb_E214M1.txt.

    The file has human-readable metadata/header lines followed by numerical
    rows with exactly:
        bin  xB  Q2  t  phi  xs  stat  syst

    Do not ask pandas to infer the mixed-format header.  Keep only lines
    containing exactly eight numeric fields.
    """
    numeric = re.compile(
        r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$"
    )
    rows = []

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.strip().split()
            if len(fields) != 8:
                continue
            #endif
            if not all(numeric.fullmatch(field) for field in fields):
                continue
            #endif
            rows.append([float(field) for field in fields])
        #endfor
    #endwith

    if not rows:
        raise RuntimeError(
            f"No 8-column numerical rows found in authoritative pass-1 file: {path}"
        )
    #endif

    df = pd.DataFrame(
        rows,
        columns=[
            "bin",
            "xB",
            "Q2",
            "t",
            "phi",
            "xs_pass1",
            "stat_pass1",
            "syst_pass1",
        ],
    )

    df["bin"] = df["bin"].astype(int)
    df["t"] = df["t"].abs()
    df["phi"] = np.mod(df["phi"], 360.0)

    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna(
        subset=["xB", "Q2", "t", "phi", "xs_pass1", "stat_pass1", "syst_pass1"]
    )
    df = df[df["xs_pass1"] > 0].reset_index(drop=True)
    # Unique identifier for an individual released (xB,Q2,t,phi) point.
    # The CLAS database `bin` value is shared by all phi points in a 3-D bin,
    # so it must NOT be used by itself for topology closure matching.
    df["pass1_row_id"] = np.arange(len(df), dtype=int)

    print(
        f"[pass1/pass2] authoritative pass-1: loaded {len(df):,} positive points "
        f"from {path.name}"
    )
    return df


def phi_distance(a: np.ndarray, b: float) -> np.ndarray:
    d = np.abs(a - b)
    return np.minimum(d, 360.0 - d)


def match_to_authoritative_pass1(
    p2: pd.DataFrame,
    p1: pd.DataFrame,
    dx: float,
    dq2: float,
    dt: float,
    dphi: float,
) -> pd.DataFrame:
    """Nearest 4-D kinematic match, with explicit maximum differences in each coordinate."""
    p1x = p1["xB"].to_numpy(float)
    p1q = p1["Q2"].to_numpy(float)
    p1t = p1["t"].to_numpy(float)
    p1p = p1["phi"].to_numpy(float)

    rows = []
    for row in p2.itertuples(index=False):
        ddx = np.abs(p1x - row.xB)
        ddq = np.abs(p1q - row.Q2)
        ddt = np.abs(p1t - row.t)
        ddp = phi_distance(p1p, row.phi)
        allowed = (ddx <= dx) & (ddq <= dq2) & (ddt <= dt) & (ddp <= dphi)
        if not np.any(allowed):
            continue
        #endif

        idxs = np.flatnonzero(allowed)
        metric = (
            (ddx[idxs] / dx) ** 2
            + (ddq[idxs] / dq2) ** 2
            + (ddt[idxs] / dt) ** 2
            + (ddp[idxs] / dphi) ** 2
        )
        j = idxs[int(np.argmin(metric))]
        ref = p1.iloc[j]
        rows.append({
            "topology": row.topology,
            "xB_pass2": row.xB,
            "Q2_pass2": row.Q2,
            "t_pass2": row.t,
            "phi_pass2": row.phi,
            "xs_pass2": row.xs_pass2,
            "pass1_row_id": int(ref["pass1_row_id"]),
            "pass1_bin": ref["bin"],
            "xB_pass1": ref["xB"],
            "Q2_pass1": ref["Q2"],
            "t_pass1": ref["t"],
            "phi_pass1": ref["phi"],
            "xs_pass1": ref["xs_pass1"],
            "stat_pass1": ref["stat_pass1"],
            "syst_pass1": ref["syst_pass1"],
            "dxB": ddx[j],
            "dQ2": ddq[j],
            "dt": ddt[j],
            "dphi": ddp[j],
            "pass2_over_pass1": row.xs_pass2 / ref["xs_pass1"],
            "N_data_signal": row.N_data_signal,
            "N_mc_rec": row.N_mc_rec,
            "N_mc_gen": row.N_mc_gen,
            "acceptance_counts": row.acceptance_counts,
            "mc_acceptance_frac_stat": row.mc_acceptance_frac_stat,
            "valid_statistics_50": bool(row.valid_statistics_50),
        })
    #endfor
    return pd.DataFrame(rows)


def discrepancy_summary(matched: pd.DataFrame, edges: List[float], labels: List[str]) -> pd.DataFrame:
    rows = []
    for topology in matched["topology"].drop_duplicates():
        top = matched[matched["topology"] == topology]
        for i, t_label in enumerate(labels):
            lo, hi = edges[i], edges[i + 1]
            sub = top[(top["t_pass2"] >= lo) & (top["t_pass2"] < hi)]
            ratio = sub["pass2_over_pass1"].to_numpy(float)
            ratio = ratio[np.isfinite(ratio) & (ratio > 0)]
            if len(ratio) == 0:
                continue
            #endif
            q16, med, q84 = np.percentile(ratio, [16, 50, 84])
            rows.append({
                "topology": topology,
                "t_bin": t_label,
                "N": len(ratio),
                "mean_t_GeV2": sub["t_pass2"].mean(),
                "median_pass2_over_pass1": med,
                "q16_pass2_over_pass1": q16,
                "q84_pass2_over_pass1": q84,
                "mean_pass2_over_pass1": np.mean(ratio),
            })
        #endfor
    #endfor
    return pd.DataFrame(rows)


def threshold_scan(
    matched: pd.DataFrame,
    edges: List[float],
    labels: List[str],
    variable: str,
    thresholds: List[int],
) -> pd.DataFrame:
    """Scan a minimum pass-2 DATA or reconstructed-MC population requirement."""
    rows = []
    for topology in matched["topology"].drop_duplicates():
        top = matched[matched["topology"] == topology]
        for i, t_label in enumerate(labels):
            lo, hi = edges[i], edges[i + 1]
            tsub = top[(top["t_pass2"] >= lo) & (top["t_pass2"] < hi)]
            for threshold in thresholds:
                sub = tsub[tsub[variable] >= threshold]
                ratio = sub["pass2_over_pass1"].to_numpy(float)
                ratio = ratio[np.isfinite(ratio) & (ratio > 0)]
                if len(ratio) == 0:
                    rows.append({
                        "topology": topology, "t_bin": t_label,
                        "threshold_variable": variable, "minimum": threshold,
                        "N": 0, "median": np.nan, "q16": np.nan, "q84": np.nan,
                    })
                    continue
                #endif
                q16, med, q84 = np.percentile(ratio, [16, 50, 84])
                rows.append({
                    "topology": topology, "t_bin": t_label,
                    "threshold_variable": variable, "minimum": threshold,
                    "N": len(ratio), "median": med, "q16": q16, "q84": q84,
                })
            #endfor
        #endfor
    #endfor
    return pd.DataFrame(rows)


def binned_ratio_vs_statistics(matched: pd.DataFrame, variable: str) -> pd.DataFrame:
    """Robust pass-2/pass-1 summaries in logarithmic DATA/MC-statistics bins."""
    positive = matched[np.isfinite(matched[variable]) & (matched[variable] > 0)].copy()
    if positive.empty:
        return pd.DataFrame()
    #endif
    edges = np.array([1, 3, 5, 10, 20, 30, 50, 100, 200, 500, 1000, 3000, 10000, np.inf], float)
    rows = []
    for topology in positive["topology"].drop_duplicates():
        top = positive[positive["topology"] == topology]
        for lo, hi in zip(edges[:-1], edges[1:]):
            sub = top[(top[variable] >= lo) & (top[variable] < hi)]
            ratio = sub["pass2_over_pass1"].to_numpy(float)
            ratio = ratio[np.isfinite(ratio) & (ratio > 0)]
            if len(ratio) == 0:
                continue
            #endif
            q16, med, q84 = np.percentile(ratio, [16, 50, 84])
            x = np.median(sub[variable].to_numpy(float))
            rows.append({
                "topology": topology, "variable": variable,
                "bin_low": lo, "bin_high": hi, "x_median": x, "N": len(ratio),
                "median": med, "q16": q16, "q84": q84,
            })
        #endfor
    #endfor
    return pd.DataFrame(rows)


def plot_ratio_vs_statistics(table: pd.DataFrame, variable: str, output: Path) -> None:
    if table.empty:
        return
    #endif
    fig, ax = plt.subplots(figsize=(8.5, 6.0))
    for topology in table["topology"].drop_duplicates():
        sub = table[table["topology"] == topology]
        x = sub["x_median"].to_numpy(float)
        y = sub["median"].to_numpy(float)
        lo = y - sub["q16"].to_numpy(float)
        hi = sub["q84"].to_numpy(float) - y
        ax.errorbar(x, y, yerr=np.vstack([lo, hi]), marker="o", capsize=2, linewidth=1.0, label=topology)
    #endfor
    ax.axhline(1.0, linewidth=1.0, linestyle="--")
    ax.set_xscale("log")
    ax.set_ylim(0.0, 1.6)
    ax.set_xlabel("pass-2 signal yield" if variable == "N_data_signal" else "pass-2 reconstructed MC events")
    ax.set_ylabel(r"median $\sigma_{\rm pass2}/\sigma_{\rm pass1}$")
    ax.set_title("Pass-2 / pass-1 stability versus bin statistics")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_threshold_scan(scan: pd.DataFrame, variable: str, output: Path) -> None:
    if scan.empty:
        return
    #endif
    topologies = [x for x in ["Inclusive", "CD-FT", "CD-FD", "FD-FD"] if x in set(scan["topology"])]
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.0), sharex=True, sharey=True)
    for ax, topology in zip(axes.flat, topologies):
        top = scan[scan["topology"] == topology]
        for t_bin in top["t_bin"].drop_duplicates():
            sub = top[top["t_bin"] == t_bin]
            ax.plot(sub["minimum"], sub["median"], marker="o", linewidth=1.2, label=t_bin)
        #endfor
        ax.axhline(1.0, linewidth=1.0, linestyle="--")
        ax.set_title(topology)
        ax.set_ylim(0.0, 1.5)
        ax.grid(alpha=0.25)
    #endfor
    for ax in axes[-1, :]:
        ax.set_xlabel("minimum signal yield" if variable == "N_data_signal" else "minimum reconstructed MC events")
    #endfor
    for ax in axes[:, 0]:
        ax.set_ylabel(r"median $\sigma_{\rm pass2}/\sigma_{\rm pass1}$")
    #endfor
    handles, leglabels = axes.flat[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, leglabels, loc="upper center", ncol=min(4, len(handles)), frameon=False)
    #endif
    fig.suptitle("Threshold stability of pass-2 / pass-1", y=0.985)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(output, dpi=180)
    plt.close(fig)


def run_statistics_stability_study(matched: pd.DataFrame, comparison_dir: Path, edges: List[float], labels: List[str]) -> None:
    stats_dir = comparison_dir / "statistics_stability"
    stats_dir.mkdir(parents=True, exist_ok=True)

    columns = [
        "topology", "pass1_row_id", "pass1_bin", "t_pass2", "pass2_over_pass1",
        "N_data_signal", "N_mc_rec", "N_mc_gen", "acceptance_counts", "mc_acceptance_frac_stat",
    ]
    matched[columns].to_csv(stats_dir / "matched_bin_statistics.csv", index=False)

    for variable, stem in [("N_data_signal", "data_signal"), ("N_mc_rec", "mc_reconstructed")]:
        binned = binned_ratio_vs_statistics(matched, variable)
        binned.to_csv(stats_dir / f"pass2_over_pass1_vs_{stem}_statistics.csv", index=False)
        plot_ratio_vs_statistics(binned, variable, stats_dir / f"pass2_over_pass1_vs_{stem}_statistics.png")

        scan = threshold_scan(matched, edges, labels, variable, STAT_THRESHOLDS)
        scan.to_csv(stats_dir / f"threshold_scan_{stem}.csv", index=False)
        plot_threshold_scan(scan, variable, stats_dir / f"threshold_scan_{stem}.png")
    #endfor

    print(f"[statistics] wrote DATA/MC stability study to {stats_dir}")
    print("[statistics] no minimum-statistics cut is imposed; threshold scans are diagnostic only")


def direct_topology_closure(matched: pd.DataFrame, a: str, b: str) -> pd.DataFrame:
    """Direct pass-2 topology ratio for the same released pass-1 (xB,Q2,t,phi) point."""
    aa = matched[matched["topology"] == a].copy()
    bb = matched[matched["topology"] == b].copy()
    if aa.empty or bb.empty:
        return pd.DataFrame()
    #endif

    # `pass1_bin` alone is NOT unique: it contains many phi points.  Use the
    # unique released-row identifier assigned by read_authoritative_pass1().
    sort_cols = ["pass1_row_id", "dt", "dQ2", "dxB", "dphi"]
    aa = aa.sort_values(sort_cols).drop_duplicates("pass1_row_id")
    bb = bb.sort_values(sort_cols).drop_duplicates("pass1_row_id")

    keep_a = [
        "pass1_row_id", "pass1_bin", "xB_pass1", "Q2_pass1", "t_pass1",
        "phi_pass1", "xs_pass1", "xs_pass2", "pass2_over_pass1"
    ]
    keep_b = ["pass1_row_id", "xs_pass2", "pass2_over_pass1"]
    out = aa[keep_a].merge(
        bb[keep_b], on="pass1_row_id", suffixes=(f"_{a}", f"_{b}")
    )
    out[f"xs_{a}_over_{b}"] = out[f"xs_pass2_{a}"] / out[f"xs_pass2_{b}"]
    return out

def plot_pass1_ratios(summary: pd.DataFrame, output: Path) -> None:
    if summary.empty:
        return
    #endif
    fig, ax = plt.subplots(figsize=(8.5, 6.0))
    for topology in summary["topology"].drop_duplicates():
        sub = summary[summary["topology"] == topology]
        x = sub["mean_t_GeV2"].to_numpy(float)
        y = sub["median_pass2_over_pass1"].to_numpy(float)
        lo = y - sub["q16_pass2_over_pass1"].to_numpy(float)
        hi = sub["q84_pass2_over_pass1"].to_numpy(float) - y
        ax.errorbar(x, y, yerr=np.vstack([lo, hi]), marker="o", capsize=3, linewidth=1.2, label=topology)
    #endfor
    ax.axhline(1.0, linewidth=1.0, linestyle="--")
    ax.set_xlabel(r"$|t|$ (GeV$^2$)")
    ax.set_ylabel(r"median $\sigma_{\rm pass2}/\sigma_{\rm pass1}$")
    ax.set_title("Pass-2 / released pass-1 cross section by detector topology")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_direct_closure(closure: pd.DataFrame, a: str, b: str, output: Path) -> None:
    if closure.empty:
        return
    #endif
    ycol = f"xs_{a}_over_{b}"
    fig, ax = plt.subplots(figsize=(8.5, 6.0))
    ax.scatter(closure["t_pass1"], closure[ycol], s=18, alpha=0.55)
    ax.axhline(1.0, linewidth=1.0, linestyle="--")
    ax.set_xlabel(r"$|t|$ (GeV$^2$)")
    ax.set_ylabel(rf"$\sigma_{{{a}}}/\sigma_{{{b}}}$")
    ax.set_title(f"Direct pass-2 topology closure: {a} versus {b}")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def run_pass1_pass2_study(args: argparse.Namespace, output_dir: Path, edges: List[float], labels: List[str]) -> None:
    comparison_dir = output_dir / "pass1_vs_pass2_topology"
    comparison_dir.mkdir(parents=True, exist_ok=True)

    pass1_path = Path(args.pass1_authoritative)
    if not pass1_path.exists():
        print(f"[pass1/pass2] WARNING: authoritative pass-1 file not found: {pass1_path}; skipping CSV comparison.")
        return
    #endif

    csv_dir = Path(args.csv_dir)
    files = resolve_pass2_files(args, csv_dir)
    if not files:
        print("[pass1/pass2] WARNING: no pass-2 CSVs resolved; skipping CSV comparison.")
        return
    #endif

    p1 = read_authoritative_pass1(pass1_path)
    matched_frames = []
    for label, path in files.items():
        p2 = read_pass2_csv(path, label)
        matched = match_to_authoritative_pass1(
            p2, p1,
            dx=args.match_dx,
            dq2=args.match_dq2,
            dt=args.match_dt,
            dphi=args.match_dphi,
        )
        print(f"[pass1/pass2] {label}: matched {len(matched):,}/{len(p2):,} pass-2 points to released pass-1")
        if not matched.empty:
            matched_frames.append(matched)
        #endif
    #endfor

    if not matched_frames:
        print("[pass1/pass2] WARNING: no pass-1/pass-2 matches survived.")
        return
    #endif

    all_matched = pd.concat(matched_frames, ignore_index=True)
    all_matched.to_csv(comparison_dir / "pass2_vs_authoritative_pass1_matched_points_all_statistics.csv", index=False)
    run_statistics_stability_study(all_matched, comparison_dir, edges, labels)

    primary = all_matched[
        np.isfinite(all_matched["N_mc_rec"])
        & (all_matched["N_mc_rec"] > args.min_bin_events)
    ].copy()
    primary.to_csv(comparison_dir / "pass2_vs_authoritative_pass1_matched_points.csv", index=False)
    print(
        f"[pass1/pass2] primary validity: reconstructed MC > {args.min_bin_events:g}; "
        f"no DATA-yield threshold; "
        f"kept {len(primary):,}/{len(all_matched):,} matched topology-points"
    )
    summary = discrepancy_summary(primary, edges, labels)
    summary.to_csv(comparison_dir / "pass2_over_pass1_summary_vs_t.csv", index=False)
    plot_pass1_ratios(summary, comparison_dir / "pass2_over_pass1_by_topology_vs_t.png")

    pairs = [("CD-FD", "CD-FT"), ("FD-FD", "CD-FD"), ("Inclusive", "CD-FD"), ("Inclusive", "CD-FT")]
    closure_summaries = []
    for a, b in pairs:
        closure = direct_topology_closure(primary, a, b)
        if closure.empty:
            continue
        #endif
        safe = f"{a}_over_{b}".replace("-", "_")
        closure.to_csv(comparison_dir / f"direct_{safe}_matched_bins.csv", index=False)
        plot_direct_closure(closure, a, b, comparison_dir / f"direct_{safe}_vs_t.png")
        ratio = closure[f"xs_{a}_over_{b}"].to_numpy(float)
        finite = np.isfinite(ratio) & (ratio > 0)
        q16, med, q84 = np.percentile(ratio[finite], [16, 50, 84])
        closure_summaries.append({"ratio": f"{a}/{b}", "N": int(np.count_nonzero(finite)), "median": med, "q16": q16, "q84": q84})
    #endfor
    pd.DataFrame(closure_summaries).to_csv(comparison_dir / "direct_topology_closure_summary.csv", index=False)

    print(f"[pass1/pass2] wrote comparison outputs to {comparison_dir}")
    print("[pass1/pass2] key plot: pass2_over_pass1_by_topology_vs_t.png")
    print("[pass1/pass2] key direct closure: direct_CD_FD_over_CD_FT_vs_t.png")



def neupane_efficiency_ratio(p_gev: np.ndarray, theta_deg: np.ndarray, phi_deg: np.ndarray) -> np.ndarray:
    """Exact production parameterization from total_counts.cpp."""
    p = np.asarray(p_gev, float)
    th = np.asarray(theta_deg, float)
    ph = np.mod(np.asarray(phi_deg, float), 360.0)
    out = np.ones_like(p)

    fd_coeff = np.asarray([
        [0.04437, -0.14271, 1.03439], [0.00490, 0.00554, 0.91770],
        [0.03671, -0.11680, 1.03002], [0.01863, -0.07756, 1.02308],
        [0.04915, -0.17173, 1.11768], [0.01077, -0.01328, 0.96242],
    ])
    cd_coeff = np.asarray([
        [0.20052, -0.79964, 1.38699], [0.16842, -0.64970, 1.31246],
        [0.18845, -0.75824, 1.41677],
    ])
    finite = np.isfinite(p) & np.isfinite(th) & np.isfinite(ph)
    fd = finite & (th < 37.0)
    cd = finite & ~fd
    for mask, coeff, width, plo, phi in [(fd, fd_coeff, 60.0, 0.4, 4.0), (cd, cd_coeff, 120.0, 0.5, 2.2)]:
        if not np.any(mask):
            continue
        #endif
        sec = np.floor(ph[mask] / width).astype(int)
        sec = np.clip(sec, 0, len(coeff) - 1)
        pe = np.clip(p[mask], plo, phi)
        c = coeff[sec]
        ratio = c[:, 0] * pe * pe + c[:, 1] * pe + c[:, 2]
        good = np.isfinite(ratio) & (ratio > 0.20) & (ratio < 1.80)
        vals = np.ones_like(ratio)
        vals[good] = ratio[good]
        out[mask] = vals
    #endfor
    return out


def load_valerii_cells(path: Path) -> pd.DataFrame:
    """Load either the current mode-6 reproduction or the legacy Valerii-map schema.

    Internal ``correction`` is always the cross-section multiplier epsilon_MC/epsilon_DATA.
    The current mode-6 ``corr_2s`` column is epsilon_DATA/epsilon_MC, so it is inverted.
    """
    v = pd.read_csv(path)
    current = {"p_lo", "p_hi", "theta_lo", "theta_hi", "phi_lo", "phi_hi",
               "data_fit_valid", "mc_fit_valid", "corr_2s"}
    legacy = {"p_low_GeV", "p_high_GeV", "theta_low_deg", "theta_high_deg",
              "phi_low_deg", "phi_high_deg", "joint_valid", "data_eff", "mc_eff"}

    if current.issubset(v.columns):
        v = v.copy()
        v["p_low_GeV"] = pd.to_numeric(v["p_lo"], errors="coerce")
        v["p_high_GeV"] = pd.to_numeric(v["p_hi"], errors="coerce")
        v["theta_low_deg"] = pd.to_numeric(v["theta_lo"], errors="coerce")
        v["theta_high_deg"] = pd.to_numeric(v["theta_hi"], errors="coerce")
        v["phi_low_deg"] = pd.to_numeric(v["phi_lo"], errors="coerce")
        v["phi_high_deg"] = pd.to_numeric(v["phi_hi"], errors="coerce")
        de = pd.to_numeric(v["data_eff_2s"], errors="coerce") if "data_eff_2s" in v else np.nan
        me = pd.to_numeric(v["mc_eff_2s"], errors="coerce") if "mc_eff_2s" in v else np.nan
        ratio = pd.to_numeric(v["corr_2s"], errors="coerce")
        v["joint_valid"] = (pd.to_numeric(v["data_fit_valid"], errors="coerce") == 1) & (pd.to_numeric(v["mc_fit_valid"], errors="coerce") == 1)
        v["data_eff"] = de
        v["mc_eff"] = me
        v["eff_data_over_mc"] = ratio
        v["correction"] = 1.0 / ratio
        v["partial_rel_unc"] = pd.to_numeric(v.get("partial_rel_unc", np.nan), errors="coerce")
    elif legacy.issubset(v.columns):
        v = v.copy()
        v["data_eff"] = pd.to_numeric(v["data_eff"], errors="coerce")
        v["mc_eff"] = pd.to_numeric(v["mc_eff"], errors="coerce")
        v["eff_data_over_mc"] = v["data_eff"] / v["mc_eff"]
        v["correction"] = v["mc_eff"] / v["data_eff"]
        if "partial_rel_unc" not in v:
            v["partial_rel_unc"] = np.nan
        #endif
    else:
        raise RuntimeError("Unrecognized Valerii map schema. Expected mode-6 valerii_fd_results.csv or legacy cell-validity CSV.")
    #endif
    return v

def fold_valerii_on_events(df: pd.DataFrame, cells: pd.DataFrame) -> pd.DataFrame:
    """Fold the FD map through reconstructed Fa18 DVCS photons with no extrapolation."""
    work = df[df["period"].isin(["Fa18 Inb", "Fa18 Out"])].copy()
    work["valerii_correction"] = 1.0
    work["valerii_eff_data_over_mc"] = np.nan
    work["valerii_partial_rel_unc"] = np.nan
    work["valerii_covered"] = False
    work["valerii_coverage_reason"] = np.where(work["photon_region"] == "FD", "invalid_or_uncovered_cell", "not_FD")

    wphi = np.where(work["g_phi_deg"].to_numpy(float) >= 330.0,
                    work["g_phi_deg"].to_numpy(float) - 360.0,
                    work["g_phi_deg"].to_numpy(float))
    gp = work["g_p"].to_numpy(float)
    gt = work["g_theta_deg"].to_numpy(float)
    is_fd = work["photon_region"].to_numpy() == "FD"
    corr = np.ones(len(work), float)
    effratio = np.full(len(work), np.nan, float)
    relunc = np.full(len(work), np.nan, float)
    covered = np.zeros(len(work), bool)

    pmin = float(pd.to_numeric(cells["p_low_GeV"], errors="coerce").min())
    pmax = float(pd.to_numeric(cells["p_high_GeV"], errors="coerce").max())
    reasons = work["valerii_coverage_reason"].to_numpy(object)
    reasons[is_fd & (gp < pmin)] = "p_below_map"
    reasons[is_fd & (gp >= pmax)] = "p_above_map"

    for row in cells.itertuples(index=False):
        joint = str(row.joint_valid).strip().lower() in {"1", "1.0", "true"}
        c = float(row.correction)
        if not joint or not np.isfinite(c) or c <= 0:
            continue
        #endif
        mask = (is_fd & (gp >= float(row.p_low_GeV)) & (gp < float(row.p_high_GeV))
                & (gt >= float(row.theta_low_deg)) & (gt < float(row.theta_high_deg))
                & (wphi >= float(row.phi_low_deg)) & (wphi < float(row.phi_high_deg)))
        corr[mask] = c
        effratio[mask] = float(row.eff_data_over_mc)
        try:
            relunc[mask] = float(row.partial_rel_unc)
        except (TypeError, ValueError):
            pass
        #endtry
        covered[mask] = True
        reasons[mask] = "covered"
    #endfor
    work["valerii_correction"] = corr
    work["valerii_eff_data_over_mc"] = effratio
    work["valerii_partial_rel_unc"] = relunc
    work["valerii_covered"] = covered
    work["valerii_coverage_reason"] = reasons
    return work

def correction_summary_from_events(df: pd.DataFrame, labels: List[str]) -> pd.DataFrame:
    work = df.copy()
    ratio = neupane_efficiency_ratio(work["p_p"].to_numpy(float), work["p_theta_deg"].to_numpy(float), work["p_phi_deg"].to_numpy(float))
    work["neupane_weight"] = 1.0 / ratio
    work["neupane_p_clamped_low"] = ((work["proton_region"] == "CD") & (work["p_p"] < 0.5)) | ((work["proton_region"] == "FD") & (work["p_p"] < 0.4))
    rows = []
    for sample, sub in [("Combined", work)] + [(p, work[work["period"] == p]) for p in PERIOD_ORDER]:
        for tlabel in labels:
            z = sub[sub["t_bin"] == tlabel]
            if z.empty:
                continue
            #endif
            rows.append({"sample": sample, "t_bin": tlabel, "N": len(z),
                         "mean_neupane_yield_multiplier": z["neupane_weight"].mean(),
                         "median_neupane_yield_multiplier": z["neupane_weight"].median(),
                         "fraction_neupane_low_p_clamped": z["neupane_p_clamped_low"].mean()})
        #endfor
    #endfor
    return pd.DataFrame(rows)


def summarize_sangbaek_pass2(repo_root: Path, edges: List[float], labels: List[str]) -> pd.DataFrame:
    path = repo_root / "output" / "data_mc_normalization" / "production_fit_test" / "eppi0_predicted_cross_section_shift.csv"
    if not path.exists():
        print(f"[correction-folding] pass-2 Sangbaek-style prediction not found: {path}")
        print("[correction-folding] run the existing eppi0 normalization validation that writes production_fit_test; not fabricating a replacement.")
        return pd.DataFrame()
    #endif
    x = pd.read_csv(path)
    x["t_center"] = 0.5 * (x["t_abs_min"] + x["t_abs_max"])
    x["t_bin"] = pd.cut(x["t_center"], bins=edges, labels=labels, right=False)
    rows=[]
    for (period,tbin), z in x.groupby(["period","t_bin"], observed=True):
        good=z[np.isfinite(z["sigma_theta_p_over_raw"]) & (z["rec_raw"]>0)]
        if good.empty: continue
        #endif
        w=good["rec_raw"].to_numpy(float); c=good["sigma_theta_p_over_raw"].to_numpy(float)
        rows.append({"period":period,"t_bin":str(tbin),"N_bins":len(good),"rec_mc_weighted_xs_multiplier":np.average(c,weights=w),"median_xs_multiplier":np.median(c)})
    #endfor
    return pd.DataFrame(rows)


def run_correction_folding(df: pd.DataFrame, args: argparse.Namespace, output_dir: Path, edges: List[float], labels: List[str]) -> None:
    out = output_dir / "correction_folding"
    out.mkdir(parents=True, exist_ok=True)
    neu = correction_summary_from_events(df, labels)
    neu.to_csv(out / "krishna_neupane_correction_vs_t.csv", index=False)

    repo_root = Path(__file__).resolve().parent.parent
    sang = summarize_sangbaek_pass2(repo_root, edges, labels)
    if not sang.empty:
        sang.to_csv(out / "pass2_sangbaek_style_correction_vs_t.csv", index=False)
    #endif

    if args.valerii_map:
        vp = Path(args.valerii_map)
        cells = load_valerii_cells(vp)
        vf = fold_valerii_on_events(df, cells)
        rows=[]
        for (period,tbin,topo), z in vf.groupby(["period","t_bin","topology"], observed=True):
            fd = z[z["photon_region"] == "FD"]
            cov = fd[fd["valerii_covered"]]
            rows.append({"period":period,"t_bin":tbin,"topology":topo,"N_events":len(z),"N_FD_photon":len(fd),
                         "N_FD_covered":len(cov),
                         "fraction_all_events_map_covered":z["valerii_covered"].mean(),
                         "fraction_FD_events_map_covered":fd["valerii_covered"].mean() if len(fd) else np.nan,
                         "fraction_FD_p_above_map":(fd["valerii_coverage_reason"] == "p_above_map").mean() if len(fd) else np.nan,
                         "fraction_FD_invalid_or_uncovered":(fd["valerii_coverage_reason"] == "invalid_or_uncovered_cell").mean() if len(fd) else np.nan,
                         "mean_multiplier_all_events_no_extrapolation":z["valerii_correction"].mean(),
                         "mean_multiplier_covered_FD":cov["valerii_correction"].mean() if len(cov) else np.nan,
                         "median_multiplier_covered_FD":cov["valerii_correction"].median() if len(cov) else np.nan,
                         "mean_eff_data_over_mc_covered_FD":cov["valerii_eff_data_over_mc"].mean() if len(cov) else np.nan,
                         "fraction_covered_FD_partial_rel_unc_gt30pct":(cov["valerii_partial_rel_unc"] > 0.30).mean() if len(cov) else np.nan})
        #endfor
        pd.DataFrame(rows).to_csv(out / "valerii_fa18_fold_vs_t_topology.csv", index=False)
        vf[["period","t_bin","topology","tabs","g_p","g_theta_deg","g_phi_deg","valerii_covered","valerii_coverage_reason","valerii_eff_data_over_mc","valerii_partial_rel_unc","valerii_correction"]].to_csv(out / "valerii_fa18_event_fold.csv", index=False)
        print(f"[correction-folding] Valerii map folded from {vp}")
    else:
        print("[correction-folding] --valerii-map not supplied; Krishna Neupane and available pass-2 Sangbaek-style outputs were still evaluated.")
    #endif
    print(f"[correction-folding] wrote outputs to {out}")

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

    run_pass1_pass2_study(args, output_dir, edges, labels)

    if args.skip_event_map:
        print(f"[study] CSV-only study complete: {output_dir}")
        return
    #endif

    topology_output_dir = output_dir / "topology_vs_t_map"
    topology_output_dir.mkdir(parents=True, exist_ok=True)

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

    run_correction_folding(df, args, output_dir, edges, labels)

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

    proton_table.to_csv(topology_output_dir / "proton_region_fractions_vs_t.csv", index=False)
    photon_table.to_csv(topology_output_dir / "photon_region_fractions_vs_t.csv", index=False)
    topology_table.to_csv(topology_output_dir / "joint_topology_fractions_vs_t.csv", index=False)
    kin_table.to_csv(topology_output_dir / "kinematic_summary_vs_t.csv", index=False)

    plot_fraction_vs_t(
        proton_table,
        edges,
        ["FD", "CD"],
        topology_output_dir / "proton_region_fractions_vs_t.png",
        "Proton detector-region fractions versus |t|",
    )
    plot_fraction_vs_t(
        photon_table,
        edges,
        ["FT", "FD"],
        topology_output_dir / "photon_region_fractions_vs_t.png",
        "Photon detector-region fractions versus |t|",
    )
    plot_fraction_vs_t(
        topology_table,
        edges,
        topology_categories[:-1],
        topology_output_dir / "joint_topology_fractions_vs_t.png",
        "Joint proton-photon topology fractions versus |t|",
    )
    plot_period_topology_fractions(
        topology_table,
        edges,
        topology_categories[:-1],
        topology_output_dir / "joint_topology_fractions_by_period.png",
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
            topology_output_dir / filename,
        )
    #endfor

    plot_joint_theta_map(
        df,
        topology_output_dir / "proton_photon_theta_correlation.png",
    )

    write_readme(topology_output_dir, args, periods, labels)

    print(f"[map] wrote outputs to {topology_output_dir}")
    print("[map] key first-look file: joint_topology_fractions_vs_t.csv")
    print("[map] key first-look plot: joint_topology_fractions_by_period.png")


if __name__ == "__main__":
    main()
