#!/usr/bin/env python3
"""
derive_photon_efficiency_scale_factors_concept.py

Fast, standalone proof-of-concept script for the RGA photon-efficiency study.

Purpose
-------
Make only the first-stage epgamma shape-comparison canvases, with essentially
no analysis machinery beyond the minimum needed to compare the samples.

For each of the five RGA periods and separately for FT and FD photons, make a
3x8 canvas comparing:

    1. Mx2
    2. Mx2_1
    3. Mx2_2
    4. Emiss2
    5. E_tag (= p2_p, because p2 is always the photon)
    6. Delta_phi (tree branch, displayed as a residual about 180 degrees)
    7. pTmiss
    8. theta_gamma_gamma

Samples:
    * nSidis data
    * exclusive-pi0 AAO MC reconstructed as epgamma
    * BH/DVCS DVCSgen MC

CLASDIS is intentionally not used anywhere in this script.

The first row uses ONLY

    angle(e, gamma) > 5 degrees,

apart from assigning the reconstructed photon to FT or FD by its polar angle.

The second row cumulatively adds

    Mx2_1 < 0.15 GeV^2.

The third row cumulatively adds

    Emiss2 > 1.0 GeV.

The FT/FD angular assignment is:

    FT: 2.4 <= theta_gamma < 5.0 degrees
    FD: 6.0 <= theta_gamma < 35.0 degrees

There are deliberately no missing-mass, missing-energy, pTmiss, Delta_phi,
electron-momentum, tag-energy, or other exclusivity cuts.

Efficiency design
-----------------
This script is intentionally optimized for speed:

    * ROOT files are streamed in chunks with uproot.iterate().
    * Events are never accumulated in large in-memory feature tables.
    * Histograms are filled immediately and the chunk is discarded.
    * A given ROOT sample is read only once; FT and FD histograms are filled
      simultaneously from the same chunk.
    * Run periods are processed concurrently with ProcessPoolExecutor.
    * ROOT files within one period are still read sequentially to avoid nested
      I/O contention and excessive memory pressure.
    * Only the 12 branches actually needed here are read.
    * No eppi0 files, event association, fitting, bootstrapping, CLASDIS,
      Stage-II, Stage-III, or grand diagnostics are touched.

Examples
--------
Run all five periods concurrently with up to 4 million events per ROOT sample:

    python derive_photon_efficiency_scale_factors_concept_v2.py --workers 5

Run one period:

    python derive_photon_efficiency_scale_factors_concept.py \
        --period fa18_inb

Use all entries:

    python derive_photon_efficiency_scale_factors_concept.py \
        --max-entries 0

Outputs are written under:

    output/photon_efficiency_concept/stage1_shape_comparison/
"""

from __future__ import annotations

import argparse
import math
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Optional, Sequence, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot


# =============================================================================
# Configuration
# =============================================================================

TREE_DEFAULT = "PhysicsEvents"

THETA_EGAMMA_MIN_DEG = 5.0

# Keep the same explicit photon angular regions used in the current study.
FT_THETA_MIN_DEG = 2.4
FT_THETA_MAX_DEG = 5.0
FD_THETA_MIN_DEG = 6.0
FD_THETA_MAX_DEG = 35.0

# Only these branches are read from the very large epgamma ROOT files.
BRANCHES: Tuple[str, ...] = (
    "Mx2",
    "Mx2_1",
    "Mx2_2",
    "Emiss2",
    "Delta_phi",
    "pTmiss",
    "theta_gamma_gamma",
    "e_theta",
    "e_phi",
    "p2_p",
    "p2_theta",
    "p2_phi",
)

# Plot ranges are intentionally centralized here so that the concept study is
# trivial to iterate on.
#
# tuple format:
#   (branch/key, title, x label, x_min, x_max, n_bins)
PLOT_SPECS = (
    (
        "Mx2",
        r"$MM^2(ep\gamma X)$",
        r"$MM^2(ep\gamma X)$ (GeV$^2$)",
        -0.05,
        0.10,
        100,
    ),
    (
        "Mx2_1",
        r"$MM^2(epX)$",
        r"$MM^2(epX)$ (GeV$^2$)",
        -0.30,
        0.70,
        140,
    ),
    (
        "Mx2_2",
        r"$MM^2(e\gamma X)$",
        r"$MM^2(e\gamma X)$ (GeV$^2$)",
        0.00,
        5.00,
        140,
    ),
    (
        "Emiss2",
        r"$E_{\rm miss}(ep\gamma X)$",
        r"$E_{\rm miss}(ep\gamma X)$ (GeV)",
        -1.0,
        5.0,
        120,
    ),
    (
        "E_tag",
        r"$E_{\rm tag}$",
        r"$E_{\rm tag}$ (GeV)",
        0.0,
        9.5,
        120,
    ),
    (
        "Delta_phi_residual_deg",
        r"$\Delta\phi(p,\gamma)$",
        r"$\Delta\phi(p,\gamma)$ residual from $180^\circ$ (deg)",
        -25.0,
        25.0,
        120,
    ),
    (
        "pTmiss",
        r"$p_{T,\rm miss}$",
        r"$p_{T,\rm miss}$ (GeV)",
        0.0,
        1.60,
        160,
    ),
    (
        "theta_gamma_gamma",
        r"$\theta_{\gamma\gamma}$",
        r"$\theta_{\gamma\gamma}$ (deg)",
        0.0,
        6.0,
        120,
    ),
)


@dataclass(frozen=True)
class Period:
    key: str
    label: str
    data: str
    pi0_mc: str
    dvcs_mc: str


_BASE = (
    "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/"
    "dvcsgen_files_greater_than_0.40GeV"
)

_DATA_BASE = (
    "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/"
    "dvcs/efficiency_study"
)

PERIODS: Tuple[Period, ...] = (
    Period(
        "fa18_inb",
        "Fa18 Inb",
        f"{_DATA_BASE}/nSidis_rga_fa18_inb_epgamma.root",
        f"{_BASE}/bkg_rga_fa18_inb_epgamma_0.40GeV.root",
        f"{_BASE}/dvcsgen_rga_fa18_inb_epgamma_0.40GeV.root",
    ),
    Period(
        "fa18_out",
        "Fa18 Out",
        f"{_DATA_BASE}/nSidis_rga_fa18_out_epgamma.root",
        f"{_BASE}/bkg_rga_fa18_out_epgamma_0.40GeV.root",
        f"{_BASE}/dvcsgen_rga_fa18_out_epgamma_0.40GeV.root",
    ),
    Period(
        "sp18_inb",
        "Sp18 Inb",
        f"{_DATA_BASE}/nSidis_rga_sp18_inb_epgamma.root",
        f"{_BASE}/bkg_rga_sp18_inb_epgamma_0.40GeV.root",
        f"{_BASE}/dvcsgen_rga_sp18_inb_epgamma_0.40GeV.root",
    ),
    Period(
        "sp18_out",
        "Sp18 Out",
        f"{_DATA_BASE}/nSidis_rga_sp18_out_epgamma.root",
        f"{_BASE}/bkg_rga_sp18_out_epgamma_0.40GeV.root",
        f"{_BASE}/dvcsgen_rga_sp18_out_epgamma_0.40GeV.root",
    ),
    Period(
        "sp19_inb",
        "Sp19 Inb",
        f"{_DATA_BASE}/nSidis_rga_sp19_inb_epgamma.root",
        f"{_BASE}/bkg_rga_sp19_inb_epgamma_0.40GeV.root",
        f"{_BASE}/dvcsgen_rga_sp19_inb_epgamma_0.40GeV.root",
    ),
)


SAMPLES = (
    ("data", "data", "black"),
    ("pi0", r"exclusive $\pi^0$ (AAO)", "tab:red"),
    ("dvcs", "BH/DVCS (DVCSgen)", "tab:blue"),
)


# =============================================================================
# Small utilities
# =============================================================================

def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Fast standalone Stage-I photon-efficiency shape comparison: "
            "data vs AAO pi0 vs DVCSgen, separately in FT and FD."
        )
    )
    parser.add_argument(
        "--period",
        action="append",
        choices=[p.key for p in PERIODS],
        help=(
            "Run only this period. May be supplied more than once. "
            "Default: all five RGA periods."
        ),
    )
    parser.add_argument(
        "--tree",
        default=TREE_DEFAULT,
        help=f"ROOT tree name. Default: {TREE_DEFAULT}.",
    )
    parser.add_argument(
        "--max-entries",
        type=int,
        default=4_000_000,
        help=(
            "Maximum entries read from each ROOT sample. "
            "Use 0 for the entire tree. Default: 4000000."
        ),
    )
    parser.add_argument(
        "--step-size",
        type=int,
        default=500_000,
        help="Entries per uproot chunk. Default: 500000.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=5,
        help=(
            "Maximum number of run periods processed concurrently. "
            "Useful range is 1-5 because there are five RGA periods. "
            "Default: 5."
        ),
    )
    parser.add_argument(
        "--output",
        default="output/photon_efficiency_concept/stage1_shape_comparison",
        help=(
            "Output directory. Default: "
            "output/photon_efficiency_concept/stage1_shape_comparison"
        ),
    )
    return parser.parse_args()


def get_tree(root_file: uproot.ReadOnlyDirectory, requested: str):
    if requested in root_file:
        return root_file[requested], requested
    #endif

    # Be forgiving if the requested name is absent, but do not silently choose
    # an arbitrary non-tree object.
    for key, classname in root_file.classnames().items():
        if classname.startswith("TTree") or "RNTuple" in classname:
            clean = key.split(";")[0]
            return root_file[key], clean
        #endif
    #endfor

    raise RuntimeError(
        f"No TTree/RNTuple found in ROOT file; requested '{requested}'."
    )


def preflight_file(path: str, tree_name: str) -> Tuple[str, int]:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(path)
    #endif

    with uproot.open(path) as root_file:
        tree, found_name = get_tree(root_file, tree_name)
        available = set(tree.keys())
        missing = [branch for branch in BRANCHES if branch not in available]
        if missing:
            raise RuntimeError(
                f"{path}: tree '{found_name}' is missing required branches: "
                + ", ".join(missing)
            )
        #endif
        return found_name, int(tree.num_entries)
    #endwith


def infer_angle_unit(
    e_theta: np.ndarray,
    e_phi: np.ndarray,
    g_theta: np.ndarray,
    g_phi: np.ndarray,
) -> str:
    # The production files have historically contained angular variables in
    # either radians or degrees depending on the processing path. Infer the
    # convention locally and robustly.
    probe = np.concatenate(
        [
            np.abs(np.asarray(e_theta, dtype=float)[:10000]),
            np.abs(np.asarray(e_phi, dtype=float)[:10000]),
            np.abs(np.asarray(g_theta, dtype=float)[:10000]),
            np.abs(np.asarray(g_phi, dtype=float)[:10000]),
        ]
    )
    probe = probe[np.isfinite(probe)]
    if probe.size == 0:
        return "rad"
    #endif
    return "deg" if float(np.nanpercentile(probe, 99.0)) > 7.0 else "rad"


def opening_angle_deg(
    theta1: np.ndarray,
    phi1: np.ndarray,
    theta2: np.ndarray,
    phi2: np.ndarray,
    unit: str,
) -> np.ndarray:
    if unit == "deg":
        theta1 = np.deg2rad(theta1)
        phi1 = np.deg2rad(phi1)
        theta2 = np.deg2rad(theta2)
        phi2 = np.deg2rad(phi2)
    else:
        theta1 = np.asarray(theta1, dtype=float)
        phi1 = np.asarray(phi1, dtype=float)
        theta2 = np.asarray(theta2, dtype=float)
        phi2 = np.asarray(phi2, dtype=float)
    #endif

    cos_angle = (
        np.cos(theta1) * np.cos(theta2)
        + np.sin(theta1) * np.sin(theta2) * np.cos(phi1 - phi2)
    )
    return np.rad2deg(np.arccos(np.clip(cos_angle, -1.0, 1.0)))


def photon_theta_deg(values: np.ndarray, unit: str) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    return values if unit == "deg" else np.rad2deg(values)


def delta_phi_residual_deg(values: np.ndarray) -> np.ndarray:
    """
    Use the stored ROOT-tree Delta_phi branch; do NOT reconstruct Delta_phi
    from particle azimuths here.

    For plotting only, express that stored coplanarity angle as its signed
    residual about 180 degrees, so a perfectly back-to-back p-gamma pair is
    displayed at 0 degrees. The current processing stores Delta_phi in radians,
    centered near pi. A defensive degree fallback is included for portability.
    """
    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]
    if finite.size and float(np.nanpercentile(np.abs(finite[:10000]), 99.0)) > 7.0:
        residual = (values - 180.0 + 180.0) % 360.0 - 180.0
        return residual
    #endif

    residual_rad = (values - math.pi + math.pi) % (2.0 * math.pi) - math.pi
    return np.rad2deg(residual_rad)


def empty_histograms() -> Dict[str, Dict[str, Dict[str, np.ndarray]]]:
    out: Dict[str, Dict[str, Dict[str, np.ndarray]]] = {}
    for region in ("FT", "FD"):
        out[region] = {}
        for row in ("minimal", "mx2_1_cut", "emiss_cut"):
            out[region][row] = {}
            for key, _title, _xlabel, x_min, x_max, n_bins in PLOT_SPECS:
                out[region][row][key] = np.zeros(n_bins, dtype=np.int64)
            #endfor
        #endfor
    #endfor
    return out


def histogram_edges() -> Dict[str, np.ndarray]:
    return {
        key: np.linspace(x_min, x_max, n_bins + 1)
        for key, _title, _xlabel, x_min, x_max, n_bins in PLOT_SPECS
    }


# =============================================================================
# Streaming histogrammer
# =============================================================================

def stream_sample_histograms(
    path: str,
    sample_label: str,
    tree_name: str,
    max_entries: int,
    step_size: int,
) -> Tuple[Dict[str, Dict[str, Dict[str, np.ndarray]]], Dict[str, int]]:
    """
    Read one ROOT sample exactly once and fill both FT and FD histograms.

    Nothing event-level survives beyond the current chunk.
    """
    hist = empty_histograms()
    edges = histogram_edges()
    counts = {
        "read": 0,
        "opening_angle": 0,
        "FT": 0,
        "FD": 0,
        "FT_mx2_1_cut": 0,
        "FD_mx2_1_cut": 0,
        "FT_emiss_cut": 0,
        "FD_emiss_cut": 0,
    }

    with uproot.open(path) as root_file:
        tree, found_tree = get_tree(root_file, tree_name)
        total = int(tree.num_entries)
        entry_stop = total if max_entries <= 0 else min(total, max_entries)

        log(
            f"{sample_label}: tree '{found_tree}', streaming "
            f"{entry_stop:,}/{total:,} entries in {step_size:,}-entry chunks."
        )

        for arrays in tree.iterate(
            expressions=list(BRANCHES),
            entry_start=0,
            entry_stop=entry_stop,
            step_size=step_size,
            library="np",
        ):
            n = len(arrays["Mx2"])
            counts["read"] += n

            angle_unit = infer_angle_unit(
                arrays["e_theta"],
                arrays["e_phi"],
                arrays["p2_theta"],
                arrays["p2_phi"],
            )

            theta_egamma = opening_angle_deg(
                arrays["e_theta"],
                arrays["e_phi"],
                arrays["p2_theta"],
                arrays["p2_phi"],
                angle_unit,
            )
            photon_theta = photon_theta_deg(arrays["p2_theta"], angle_unit)

            # This is intentionally the only event-selection cut.
            base = np.isfinite(theta_egamma) & (theta_egamma > THETA_EGAMMA_MIN_DEG)
            counts["opening_angle"] += int(np.count_nonzero(base))

            region_masks = {
                "FT": (
                    base
                    & np.isfinite(photon_theta)
                    & (photon_theta >= FT_THETA_MIN_DEG)
                    & (photon_theta < FT_THETA_MAX_DEG)
                ),
                "FD": (
                    base
                    & np.isfinite(photon_theta)
                    & (photon_theta >= FD_THETA_MIN_DEG)
                    & (photon_theta < FD_THETA_MAX_DEG)
                ),
            }

            values = {
                "Mx2": np.asarray(arrays["Mx2"], dtype=float),
                "Mx2_1": np.asarray(arrays["Mx2_1"], dtype=float),
                "Mx2_2": np.asarray(arrays["Mx2_2"], dtype=float),
                "Emiss2": np.asarray(arrays["Emiss2"], dtype=float),
                # For a photon, E = |p|, and p2 is always the photon in these
                # processed epgamma trees.
                "E_tag": np.asarray(arrays["p2_p"], dtype=float),
                "Delta_phi_residual_deg": delta_phi_residual_deg(
                    arrays["Delta_phi"]
                ),
                "pTmiss": np.asarray(arrays["pTmiss"], dtype=float),
                # kinematic_variables.java explicitly stores this in degrees.
                "theta_gamma_gamma": np.asarray(
                    arrays["theta_gamma_gamma"], dtype=float
                ),
            }

            for region, region_mask in region_masks.items():
                counts[region] += int(np.count_nonzero(region_mask))

                # Row 1: only Angle(e,gamma) > 5 deg plus FT/FD assignment.
                mx2_1_mask = (
                    region_mask
                    & np.isfinite(values["Mx2_1"])
                    & (values["Mx2_1"] < 0.15)
                )

                # Row 3 is cumulative: Row 2 plus Emiss2 > 1 GeV.
                emiss_mask = (
                    mx2_1_mask
                    & np.isfinite(values["Emiss2"])
                    & (values["Emiss2"] > 1.0)
                )

                row_masks = {
                    "minimal": region_mask,
                    "mx2_1_cut": mx2_1_mask,
                    "emiss_cut": emiss_mask,
                }

                counts[f"{region}_mx2_1_cut"] += int(
                    np.count_nonzero(mx2_1_mask)
                )
                counts[f"{region}_emiss_cut"] += int(
                    np.count_nonzero(emiss_mask)
                )

                for row_name, row_mask in row_masks.items():
                    for key in values:
                        v = values[key]
                        mask = row_mask & np.isfinite(v)
                        if np.any(mask):
                            h, _ = np.histogram(v[mask], bins=edges[key])
                            hist[region][row_name][key] += h.astype(
                                np.int64,
                                copy=False,
                            )
                        #endif
                    #endfor
                #endfor
            #endfor

            if counts["read"] == entry_stop or counts["read"] % max(step_size, 1_000_000) == 0:
                frac = 100.0 * counts["read"] / entry_stop if entry_stop else 100.0
                log(
                    f"{sample_label}: {counts['read']:,}/{entry_stop:,} "
                    f"({frac:.1f}%)"
                )
            #endif
        #endfor
    #endwith

    log(
        f"{sample_label}: retained after Angle(e,gamma)>5 deg: "
        f"FT={counts['FT']:,}, FD={counts['FD']:,}; "
        f"after additionally Mx2_1<0.15 GeV^2: "
        f"FT={counts['FT_mx2_1_cut']:,}, "
        f"FD={counts['FD_mx2_1_cut']:,}; "
        f"after additionally Emiss2>1 GeV: "
        f"FT={counts['FT_emiss_cut']:,}, "
        f"FD={counts['FD_emiss_cut']:,}."
    )
    return hist, counts


# =============================================================================
# Plotting
# =============================================================================

def normalized_hist(counts: np.ndarray) -> np.ndarray:
    counts = np.asarray(counts, dtype=float)
    total = float(np.sum(counts))
    if total <= 0.0:
        return np.zeros_like(counts)
    #endif
    return counts / total


def make_canvas(
    period: Period,
    region: str,
    sample_hists: Dict[str, Dict[str, Dict[str, Dict[str, np.ndarray]]]],
    sample_counts: Dict[str, Dict[str, int]],
    output_dir: Path,
) -> Path:
    fig, axes = plt.subplots(3, 8, figsize=(27.5, 11.8))
    edges_by_key = histogram_edges()

    sample_meta = {key: (label, color) for key, label, color in SAMPLES}
    rows = (
        (
            "minimal",
            r"Row 1: only $\mathrm{Angle}(e,\gamma)>5^\circ$",
        ),
        (
            "mx2_1_cut",
            r"Row 2: additionally $Mx2_1<0.15$ GeV$^2$",
        ),
        (
            "emiss_cut",
            r"Row 3: additionally $E_{\rm miss}>1$ GeV",
        ),
    )

    legend_handles = []
    legend_labels = []

    for irow, (row_key, row_label) in enumerate(rows):
        for icol, spec in enumerate(PLOT_SPECS):
            ax = axes[irow, icol]
            key, title, xlabel, _x_min, _x_max, _n_bins = spec
            edges = edges_by_key[key]
            centers = 0.5 * (edges[:-1] + edges[1:])

            for sample_key, _label, _color in SAMPLES:
                label, color = sample_meta[sample_key]
                y = normalized_hist(
                    sample_hists[sample_key][region][row_key][key]
                )
                line = ax.step(
                    centers,
                    y,
                    where="mid",
                    linewidth=1.5,
                    color=color,
                    label=label,
                )[0]

                if irow == 0 and icol == 0:
                    legend_handles.append(line)
                    legend_labels.append(label)
                #endif
            #endfor

            if irow == 0:
                ax.set_title(title, fontsize=11)
            #endif

            ax.set_xlabel(xlabel, fontsize=9.2)
            if icol == 0:
                ax.set_ylabel(row_label + "\nunit-normalized", fontsize=9.2)
            else:
                ax.set_ylabel("unit-normalized", fontsize=9.2)
            #endif

            ax.set_xlim(edges[0], edges[-1])
            ax.tick_params(axis="both", labelsize=7.8)
            ax.grid(alpha=0.18)

            # Show the cumulative selection thresholds on the relevant panels.
            if key == "Mx2_1":
                ax.axvline(
                    0.15,
                    linestyle="--",
                    linewidth=1.0,
                    color="0.35",
                )
            #endif
            if key == "Emiss2":
                ax.axvline(
                    1.0,
                    linestyle="--",
                    linewidth=1.0,
                    color="0.35",
                )
            #endif
        #endfor
    #endfor

    region_text = (
        rf"FT: ${FT_THETA_MIN_DEG:.1f}^\circ\leq\theta_\gamma"
        rf"<{FT_THETA_MAX_DEG:.1f}^\circ$"
        if region == "FT"
        else
        rf"FD: ${FD_THETA_MIN_DEG:.1f}^\circ\leq\theta_\gamma"
        rf"<{FD_THETA_MAX_DEG:.1f}^\circ$"
    )

    row1_counts = ", ".join(
        f"{sample_meta[key][0]}={sample_counts[key][region]:,}"
        for key, _label, _color in SAMPLES
    )
    row2_counts = ", ".join(
        f"{sample_meta[key][0]}="
        f"{sample_counts[key][region + '_mx2_1_cut']:,}"
        for key, _label, _color in SAMPLES
    )
    row3_counts = ", ".join(
        f"{sample_meta[key][0]}="
        f"{sample_counts[key][region + '_emiss_cut']:,}"
        for key, _label, _color in SAMPLES
    )

    fig.suptitle(
        f"{period.label}: {region} epgamma shape comparison\n"
        rf"Row 1: $\mathrm{{Angle}}(e,\gamma)>"
        rf"{THETA_EGAMMA_MIN_DEG:.0f}^\circ$; "
        + region_text
        + r"; Row 2: $+\,Mx2_1<0.15$ GeV$^2$"
        + r"; Row 3: $+\,E_{\rm miss}>1$ GeV"
        + "\n"
        + f"Counts — row 1: {row1_counts}    |    "
        + f"row 2: {row2_counts}    |    row 3: {row3_counts}",
        fontsize=10.8,
        y=0.992,
    )

    # One boxed canvas-level legend centered directly beneath the title.
    fig.legend(
        legend_handles,
        legend_labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.905),
        ncol=3,
        fontsize=9.5,
        frameon=True,
        fancybox=False,
        edgecolor="black",
    )

    fig.subplots_adjust(
        left=0.048,
        right=0.995,
        bottom=0.075,
        top=0.82,
        wspace=0.36,
        hspace=0.46,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    out = output_dir / f"canvas_shape_comparison_{period.key}_{region.lower()}.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


# =============================================================================
# Period driver
# =============================================================================

def process_period(
    period: Period,
    tree_name: str,
    max_entries: int,
    step_size: int,
    output_dir: Path,
) -> None:
    t0 = time.perf_counter()
    log(f"{period.label}: starting.")

    paths = {
        "data": period.data,
        "pi0": period.pi0_mc,
        "dvcs": period.dvcs_mc,
    }

    # Fail early, before reading millions of events.
    for sample_key, path in paths.items():
        found_tree, total = preflight_file(path, tree_name)
        log(
            f"{period.label}: preflight {sample_key}: "
            f"{Path(path).name}, tree '{found_tree}', {total:,} entries."
        )
    #endfor

    sample_hists: Dict[
        str,
        Dict[str, Dict[str, Dict[str, np.ndarray]]],
    ] = {}
    sample_counts: Dict[str, Dict[str, int]] = {}

    for sample_key, label, _color in SAMPLES:
        hists, counts = stream_sample_histograms(
            paths[sample_key],
            f"{period.label} {label}",
            tree_name,
            max_entries,
            step_size,
        )
        sample_hists[sample_key] = hists
        sample_counts[sample_key] = counts
    #endfor

    for region in ("FT", "FD"):
        out = make_canvas(
            period,
            region,
            sample_hists,
            sample_counts,
            output_dir,
        )
        log(f"{period.label}: wrote {out}")
    #endfor

    log(
        f"{period.label}: complete in "
        f"{time.perf_counter() - t0:.1f} s."
    )


def process_period_worker(
    period: Period,
    tree_name: str,
    max_entries: int,
    step_size: int,
    output_dir_str: str,
) -> str:
    """
    Picklable process-level wrapper for one run period.

    Each worker streams that period's data, AAO, and DVCSgen files
    sequentially. FT and FD are filled simultaneously within each file,
    avoiding nested ROOT I/O parallelism.
    """
    process_period(
        period,
        tree_name,
        max_entries,
        step_size,
        Path(output_dir_str),
    )
    return period.key


def main() -> int:
    args = parse_args()

    if args.max_entries < 0:
        raise ValueError("--max-entries must be >= 0.")
    #endif
    if args.step_size <= 0:
        raise ValueError("--step-size must be > 0.")
    #endif
    if args.workers <= 0:
        raise ValueError("--workers must be > 0.")
    #endif

    selected_keys = set(args.period or [p.key for p in PERIODS])
    selected_periods = [p for p in PERIODS if p.key in selected_keys]

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    log(
        "Standalone photon-efficiency concept study: "
        "data + AAO pi0 + DVCSgen only; CLASDIS disabled."
    )
    log(
        "ONLY event-selection cut: "
        f"Angle(e,gamma) > {THETA_EGAMMA_MIN_DEG:.1f} deg. "
        f"Photon regions: FT {FT_THETA_MIN_DEG:.1f}-{FT_THETA_MAX_DEG:.1f} deg, "
        f"FD {FD_THETA_MIN_DEG:.1f}-{FD_THETA_MAX_DEG:.1f} deg."
    )

    n_workers = min(int(args.workers), len(selected_periods))
    log(
        f"Period-level parallelism: {n_workers} process(es) for "
        f"{len(selected_periods)} selected period(s)."
    )

    if n_workers == 1:
        for period in selected_periods:
            process_period(
                period,
                args.tree,
                args.max_entries,
                args.step_size,
                output_dir,
            )
        #endfor
    else:
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            future_to_period = {
                executor.submit(
                    process_period_worker,
                    period,
                    args.tree,
                    args.max_entries,
                    args.step_size,
                    str(output_dir),
                ): period
                for period in selected_periods
            }

            for future in as_completed(future_to_period):
                period = future_to_period[future]
                try:
                    completed_key = future.result()
                    log(f"{period.label}: worker completed ({completed_key}).")
                except Exception as exc:
                    raise RuntimeError(
                        f"Parallel concept-study processing failed for "
                        f"{period.label}: {exc}"
                    ) from exc
                #endtry
            #endfor
        #endwith
    #endif

    log(f"Done. Outputs are in {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
