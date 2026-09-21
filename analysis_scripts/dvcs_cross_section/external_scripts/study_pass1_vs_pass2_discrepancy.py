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
import re
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
        description="Study the pass-1/pass-2 DVCS discrepancy and map detector topology versus |t|."
    )
    parser.add_argument(
        "--output-dir",
        default="output/study_pass1_vs_pass2_discrepancy",
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
    return {
        "xB": find_column(df, ["xB", "x_B", "x", "<xB>", "<x_B>"], ["xb"]),
        "Q2": find_column(df, ["Q2", "Q^2", "<Q2>", "q2"], ["q2"]),
        "t": find_column(df, ["t", "-t", "|t|", "<t>"], None),
        "phi": find_column(df, ["phi", "phi_deg", "phi (deg)", "<phi>"], ["phi"]),
    }


def combined_xs_column(df: pd.DataFrame) -> str:
    preferred = "cross sections, ep->epg, exp, 10.6 GeV, unpol"
    if preferred in df.columns:
        return preferred
    #endif

    matches = [
        c for c in df.columns
        if "crosssections" in normalized_name(c)
        and "epepg" in normalized_name(c)
        and "exp" in normalized_name(c)
        and "unpol" in normalized_name(c)
        and "106gev" in normalized_name(c)
    ]
    if len(matches) == 1:
        return matches[0]
    #endif

    raise RuntimeError(
        "Could not uniquely identify combined 10.6-GeV unpolarized DVCS cross-section column. "
        f"Candidates found: {matches[:12]}"
    )


def read_pass2_csv(path: Path, label: str) -> pd.DataFrame:
    raw = pd.read_csv(path, low_memory=False)
    kin = find_kinematic_columns(raw)
    xs_col = combined_xs_column(raw)

    out = pd.DataFrame({
        "xB": parse_numeric(raw[kin["xB"]]),
        "Q2": parse_numeric(raw[kin["Q2"]]),
        "t": parse_numeric(raw[kin["t"]]).abs(),
        "phi": np.mod(parse_numeric(raw[kin["phi"]]), 360.0),
        "xs_pass2": parse_numeric(raw[xs_col]),
    })
    out["topology"] = label
    out["source_file"] = str(path)
    out = out.replace([np.inf, -np.inf], np.nan)
    out = out.dropna(subset=["xB", "Q2", "t", "phi", "xs_pass2"])
    out = out[out["xs_pass2"] > 0].reset_index(drop=True)
    print(f"[pass1/pass2] {label}: loaded {len(out):,} positive pass-2 cross sections from {path.name}")
    return out


def discover_csv(csv_dir: Path, include_tokens: List[str], exclude_tokens: List[str] | None = None) -> Path | None:
    """
    Discover a pass-2 CSV without assuming one exact filename convention.

    The analysis suite has used names containing forms such as CD-FD, CD_FD,
    CDFD, etc.  Compare both a separator-preserving lowercase name and a
    fully normalized alphanumeric name.
    """
    files = sorted(csv_dir.glob("*.csv"))
    include = [normalized_name(x) for x in include_tokens]
    exclude = [normalized_name(x) for x in (exclude_tokens or [])]
    matches = []

    for path in files:
        raw_key = path.stem.lower()
        norm_key = normalized_name(path.stem)

        include_ok = all(token in norm_key for token in include)
        exclude_hit = any(token in norm_key for token in exclude)

        if include_ok and not exclude_hit:
            matches.append(path)
        #endif
    #endfor

    if len(matches) == 1:
        return matches[0]
    #endif
    if len(matches) > 1:
        print(f"[pass1/pass2] WARNING: multiple files match {include_tokens}: {[p.name for p in matches]}")
        # Prefer the shortest basename: parent topology files are normally
        # shorter than sector-resolved derivatives.
        matches = sorted(matches, key=lambda p: (len(p.stem), p.name))
        return matches[0]
    #endif
    return None

def resolve_pass2_files(args: argparse.Namespace) -> Dict[str, Path]:
    csv_dir = Path(args.csv_dir)
    explicit = {
        "Inclusive": args.inclusive_csv,
        "CD-FD": args.cd_fd_csv,
        "CD-FT": args.cd_ft_csv,
        "FD-FD": args.fd_fd_csv,
    }
    expected_names = {
        "Inclusive": "dvcs_pass2_INCLUSIVE.csv",
        "CD-FD": "dvcs_pass2_CD-FD.csv",
        "CD-FT": "dvcs_pass2_CD-FT.csv",
        "FD-FD": "dvcs_pass2_FD-FD.csv",
    }

    resolved = {}
    for label in ["Inclusive", "CD-FD", "CD-FT", "FD-FD"]:
        if explicit[label]:
            path = Path(explicit[label])
        else:
            expected = csv_dir / expected_names[label]
            path = expected if expected.exists() else None
        #endif

        if path is None or not path.exists():
            print(f"[pass1/pass2] WARNING: no {label} CSV found; skipping it.")
            if label != "Inclusive":
                available = [p.name for p in sorted(csv_dir.glob("*.csv"))]
                print(f"[pass1/pass2] Available CSVs in {csv_dir}:")
                for fname in available:
                    print(f"    {fname}")
                #endfor
            #endif
            continue
        #endif
        resolved[label] = path
    #endfor
    return resolved


def read_authoritative_pass1(path: Path) -> pd.DataFrame:
    """
    Read the released CLAS database E214M1 text file.

    The file begins with title/author/header/unit lines followed by
    whitespace-separated rows:
      bin x Q2 t_average phi xs stat syst

    Parse only lines containing exactly eight numeric fields.  This avoids
    treating the human-readable header as a variable-width CSV.
    """
    number = re.compile(
        r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$"
    )
    rows = []

    with path.open("r", errors="replace") as handle:
        for line in handle:
            fields = line.strip().split()
            if len(fields) != 8:
                continue
            #endif
            if not all(number.match(field) for field in fields):
                continue
            #endif
            rows.append([float(field) for field in fields])
        #endfor
    #endwith

    if not rows:
        raise RuntimeError(f"{path}: found no eight-column numeric pass-1 rows")
    #endif

    df = pd.DataFrame(
        rows,
        columns=["bin", "xB", "Q2", "t", "phi", "xs_pass1", "stat_pass1", "syst_pass1"],
    )
    df["bin"] = df["bin"].astype(int)
    df["t"] = df["t"].abs()
    df["phi"] = np.mod(df["phi"], 360.0)
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna(subset=["xB", "Q2", "t", "phi", "xs_pass1"])
    df = df[df["xs_pass1"] > 0].reset_index(drop=True)

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


def direct_topology_closure(matched: pd.DataFrame, a: str, b: str) -> pd.DataFrame:
    """Direct pass-2 topology ratio in bins sharing the same released pass-1 point."""
    aa = matched[matched["topology"] == a].copy()
    bb = matched[matched["topology"] == b].copy()
    if aa.empty or bb.empty:
        return pd.DataFrame()
    #endif

    aa = aa.sort_values(["pass1_bin", "dt", "dQ2", "dxB", "dphi"]).drop_duplicates("pass1_bin")
    bb = bb.sort_values(["pass1_bin", "dt", "dQ2", "dxB", "dphi"]).drop_duplicates("pass1_bin")
    keep_a = ["pass1_bin", "xB_pass1", "Q2_pass1", "t_pass1", "phi_pass1", "xs_pass1", "xs_pass2", "pass2_over_pass1"]
    keep_b = ["pass1_bin", "xs_pass2", "pass2_over_pass1"]
    out = aa[keep_a].merge(bb[keep_b], on="pass1_bin", suffixes=(f"_{a}", f"_{b}"))
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

    files = resolve_pass2_files(args)
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
    all_matched.to_csv(comparison_dir / "pass2_vs_authoritative_pass1_matched_points.csv", index=False)
    summary = discrepancy_summary(all_matched, edges, labels)
    summary.to_csv(comparison_dir / "pass2_over_pass1_summary_vs_t.csv", index=False)
    plot_pass1_ratios(summary, comparison_dir / "pass2_over_pass1_by_topology_vs_t.png")

    pairs = [("CD-FD", "CD-FT"), ("FD-FD", "CD-FD"), ("Inclusive", "CD-FD"), ("Inclusive", "CD-FT")]
    closure_summaries = []
    for a, b in pairs:
        closure = direct_topology_closure(all_matched, a, b)
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
