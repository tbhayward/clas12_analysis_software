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
    * The purity/shape stage reads only the compact epgamma branch set.
    * The final efficiency stage then makes ONE additional epgamma pass and
      reads the companion eppi0 data/AAO files needed for tag-probe matching.
    * CLASDIS and the old broad development machinery remain disabled.
    * A lightweight AAO-vs-DVCS template-fit closure test is run after every
      cumulative optimizer step to quantify whether the surviving templates
      remain distinguishable enough to fit their mixture fraction.

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

By default all shape-comparison panels use a logarithmic y-axis. Add

    --linear

for a linear y-axis.

The script runs the greedy rectangular-cut optimizer BEFORE making the
shape-comparison canvases. The optimizer is NEVER combined across periods:
each of the five run periods is optimized independently, and FT/FD are
optimized independently within each period.

Before optimization, a loose exclusivity preselection is imposed:

    Mx2_1 < 0.30 GeV^2

by default. The real-data interval

    0.20 <= Mx2_1 < 0.30 GeV^2

is used as an empirical nonexclusive/SIDIS-like sideband. Because Mx2_1
defines both the loose preselection and the sideband, Mx2_1 itself is excluded
from the greedy threshold scan. The optimizer therefore seeks reconstructed
cuts that simultaneously retain AAO pi0, suppress DVCSgen, and suppress this
data sideband.

After the optimizer selects its period- and detector-specific cut sequence,
the shape-comparison canvas is constructed with one cumulative row per selected
cut. The same selected cuts are then shown on data, AAO, and DVCSgen.

The optimizer is run independently in four predicted-probe-energy bins:

    0.4 <= E_probe < 1.0 GeV
    1.0 <= E_probe < 2.0 GeV
    2.0 <= E_probe < 4.0 GeV
    4.0 <= E_probe < 9.5 GeV

with E_probe identified here with the stored Emiss2 quantity used by the
missing-photon construction. Because E_probe now defines the optimization bin,
Emiss2 is removed completely from the optimizer variables: the optimizer may
not sculpt the probe-energy spectrum inside or across the requested bins.

Shape canvases are written under:

    output/photon_efficiency_concept/stage1_shape_comparison/

The script intentionally keeps machine-readable output compact.  For each
period it writes ONE combined CSV containing both optimizer and closure results:

    output/photon_efficiency_concept/summary/
        concept_summary_<period>.csv

No per-bin TXT files and no per-bin optimizer/closure CSV files are produced.

The main diagnostics are visual:

    output/photon_efficiency_concept/rectangular_optimizer/
        rectangular_optimizer_<period>.png

    output/photon_efficiency_concept/template_fit_closure/
        template_fit_closure_<period>.png
        template_fit_response_<period>.png

The first closure canvas shows RMS error in the fitted pi0 fraction versus cut
step. The response canvas shows fitted versus injected pi0 fraction, with the
ideal y=x relation, so loss of template identifiability is visible directly.

A final real-data composition-fit stage then chooses one closure-safe cumulative
cut step independently in every period x detector x E_probe bin and fits the
surviving nSidis data with THREE templates:

    data = f_pi0 * AAO_morphed
         + f_DVCS * DVCSgen
         + f_SB * data_sideband,

with non-negative fractions constrained to sum to one.

The fit target uses Mx2_1 below the sideband lower edge. The empirical
nonexclusive template comes from the disjoint configured real-data Mx2_1
sideband and receives exactly the same cumulative optimizer cuts. Mx2_1 and
Emiss2/E_probe are excluded as real-data fit discriminators.

Only AAO is allowed to morph: a bounded centroid shift plus additional
Gaussian smearing in units of fixed fit-bin widths. DVCSgen and the empirical
sideband template remain fixed. The fit plots report the deviance before and
after AAO morphing so the benefit of the detector-resolution nuisance model is
immediately visible.

A dedicated sideband-transfer validation is also performed after the operating
cut step is selected. The real data are divided into three adjacent Mx2_1
slices,

    0.15 <= Mx2_1 < 0.20 GeV^2
    0.20 <= Mx2_1 < 0.25 GeV^2
    0.25 <= Mx2_1 < 0.30 GeV^2,

and the SAME cumulative optimizer cuts are applied to all three. Their
unit-normalized shapes are compared in the exact variable used by the
composition fit. Pairwise Jensen-Shannon divergences are reported on the plot.
The lowest slice is deliberately treated as a signal-adjacent diagnostic, not
as an independent background template; the key background-transfer test is the stability of the eta-safe
0.15-0.19 and 0.19-0.23 GeV^2 slices.

The empirical sideband is used adaptively.  A three-template fit is attempted
only when, after the selected operating cuts, both eta-safe sidebands have
adequate statistics and the near/far sideband shapes are sufficiently stable.
Otherwise the code falls back explicitly to an AAO + DVCSgen fit, retaining the
same bounded AAO-only morphing.  The fallback is labeled on the plots and in
the single per-period summary CSV rather than silently pretending that the
sideband fraction is known.

TAG-FD / predicted-probe-sector diagnostics are now added WITHOUT refitting the physics composition
independently in each sector. The optimizer cuts, operating step, fit
discriminator, fitted composition fractions, and AAO morph nuisance parameters
remain those obtained from the period-level FD_all fit. The same quantities
are then applied independently to sectors 1-6 using predicted-probe-sector-specific data/MC
template shapes. This exposes detector-sector response differences without
allowing a sector-by-sector composition fit to absorb them.

The sector diagnostic files are kept in one dedicated hierarchy:

    output/photon_efficiency_concept/fd_sector_diagnostics/<period>/
        sector_template_fit_<Eprobe-bin>.png
        sector_overview_<period>.png
        sector_summary_<period>.csv

There are four sector-fit canvases plus one overview and one compact CSV per
period; no per-sector files are generated.

The visual outputs per period are:

    output/photon_efficiency_concept/data_template_fit/
        data_template_fit_<period>.png
        pi0_fraction_summary_<period>.png

The first is a 4x2 canvas showing data, fitted total, the two fitted template
components, the chosen cut step/discriminator, f_pi0 +/- statistical error,
and deviance/ndf. The second summarizes the fitted pi0 fraction versus E_probe
for FD and FT.


Final tag-probe efficiency extraction
-------------------------------------
After the period-level tag-purity fits are complete, the script performs the
actual photon-efficiency extraction.

The tag/probe roles are deliberately distinct:

  * reconstructed p2 is the TAG photon;
  * tag FT/FD classification therefore comes from reconstructed p2_theta;
  * the missing PARTNER/PROBE is predicted from

        p_probe^pred = p_beam - p_e - p_p - p_gamma,tag;

  * predicted probe theta/phi determine the detector efficiency bin.

The data denominator is a purity-weighted expected-pi0 count.  In one probe
efficiency bin b,

    N_expected^data(b)
      = sum_{tag category t} f_pi0(t,E_probe) N_tag^data(t,E_probe,b).

The fitted f_pi0 is NOT allowed to depend on predicted probe sector.

The AAO denominator is simply the number of selected reconstructed epgamma tags
because the AAO sample is exclusive pi0 MC.

The numerator asks whether the predicted partner is present in the reconstructed
eppi0 sample:

  * DATA: epgamma/eppi0 association first requires exact (runnum,evnum), then
    checks the reconstructed e/p parent consistency and the reconstructed-side
    pi0/tag remainder.
  * AAO MC: runnum/evnum are deliberately NOT used (MC runnum is always 11 and
    evnum is not a reliable cross-file identifier).  The two trees are matched
    exclusively through reconstructed electron/proton Cartesian momenta using
    a six-dimensional cKDTree, as in the established photon-efficiency study.

The reconstructed partner is obtained from the eppi0 tree as

    P_probe^reco = P_pi0^reco - k_tag^reco,

because these eppi0 trees store the reconstructed pi0 four-vector rather than
the two daughter-photon four-vectors.

A clean association requires the reconstructed pi0 mass window, an
approximately massless positive-energy remainder, and the photon threshold.
The final probe match uses Delta p_x, Delta p_y, Delta p_z windows derived from
the AAO MC itself.  The MC residual center and robust sigma are determined
separately for predicted FT and FD probes; the same MC-derived windows are then
applied to data and MC.

Efficiency binning:

  * predicted FT: one FT bin in every E_probe bin;
  * predicted FD, E_probe < 4 GeV: sectors 1-6 independently;
  * predicted FD, 4 <= E_probe < 9.5 GeV: all six sectors combined (FD_all).

All five run periods remain independent.

Outputs are intentionally compact and plot-driven:

    output/photon_efficiency_concept/efficiency/<period>/
        matching_resolution_<period>.png
        association_cutflow_<period>.png
        selected_event_overlap_<period>.png
        matching_resolution_<period>.png
        association_cutflow_<period>.png
        efficiency_counts_<period>.png
        efficiency_summary_<period>.png
        efficiency_summary_<period>.csv

There is exactly one machine-readable efficiency CSV per period; no per-bin or
per-sector CSV/TXT files are created.

For every cumulative cut step, the closure test automatically chooses the
single surviving reconstructed variable with the largest AAO-vs-DVCS
Jensen-Shannon divergence, builds independent training and held-out templates,
creates pseudo-data at known pi0 fractions, and refits those fractions.
"""

from __future__ import annotations

import argparse
import csv
import math
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot
from scipy.spatial import cKDTree


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


# Particle masses used in the explicit four-vector tag-probe construction.
M_E_GEV = 0.00051099895
M_P_GEV = 0.93827208816

# The highest FD energy bin is intentionally NOT split by sector because the
# present statistics are too sparse after the purity/exclusivity selections.
HIGH_ENERGY_FD_COMBINE_MIN_GEV = 4.0


# Predicted missing-photon energy bins used for the purity optimization.
# These are intentionally kept fixed for this concept iteration.
#
# tuple format:
#   (key, lower edge [GeV], upper edge [GeV])
E_PROBE_BINS: Tuple[Tuple[str, float, float], ...] = (
    ("0p4_1p0", 0.4, 1.0),
    ("1p0_2p0", 1.0, 2.0),
    ("2p0_4p0", 2.0, 4.0),
    ("4p0_9p5", 4.0, 9.5),
)


# Adjacent Mx2_1 slices used only to validate whether the empirical
# nonexclusive sideband shape is stable as one approaches the signal region.
# The composition fit is deliberately kept below the eta region.
#
#   signal:       Mx2_1 < 0.15 GeV^2
#   near SB:      0.15 <= Mx2_1 < 0.19 GeV^2  (fit template)
#   farther SB:   0.19 <= Mx2_1 < 0.23 GeV^2  (transfer validation)
#
# The eta mass-squared is ~0.300 GeV^2, so neither sideband approaches the
# exclusive eta peak.  A separate eta-control region is retained only as a
# diagnostic so possible eta contamination remains visible rather than being
# silently discarded.
SIGNAL_MX2_1_MAX = 0.15
NEAR_SIDEBAND_MIN = 0.15
NEAR_SIDEBAND_MAX = 0.19
FAR_SIDEBAND_MIN = 0.19
FAR_SIDEBAND_MAX = 0.23

ETA_MASS_GEV = 0.547862
ETA_MASS2_GEV2 = ETA_MASS_GEV * ETA_MASS_GEV
PI0_MASS_GEV = 0.134977
PI0_MASS2_GEV2 = PI0_MASS_GEV * PI0_MASS_GEV

# Used only for the transfer/stability plot.
SIDEBAND_VALIDATION_SLICES: Tuple[Tuple[str, float, float], ...] = (
    ("signal_edge", 0.11, 0.15),
    ("near", NEAR_SIDEBAND_MIN, NEAR_SIDEBAND_MAX),
    ("far", FAR_SIDEBAND_MIN, FAR_SIDEBAND_MAX),
)

# Only these branches are read from the very large epgamma ROOT files.
BRANCHES: Tuple[str, ...] = (
    "Mx2",
    "Mx2_1",
    "Mx2_2",
    "Emiss2",
    "Delta_phi",
    "pTmiss",
    "theta_gamma_gamma",
    "e_p",
    "e_theta",
    "e_phi",
    "p1_p",
    "p1_theta",
    "p1_phi",
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
    eppi0_data: str
    eppi0_pi0_mc: str


_BASE = (
    "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/"
    "dvcsgen_files_greater_than_0.40GeV"
)

_DATA_BASE = (
    "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/"
    "dvcs/efficiency_study"
)

BEAM_ENERGY_GEV = {
    "fa18_inb": 10.604,
    "fa18_out": 10.604,
    "sp18_inb": 10.594,
    "sp18_out": 10.594,
    "sp19_inb": 10.200,
}


PERIODS: Tuple[Period, ...] = (
    Period(
        "fa18_inb",
        "Fa18 Inb",
        f"{_DATA_BASE}/nSidis_rga_fa18_inb_epgamma.root",
        f"{_BASE}/bkg_rga_fa18_inb_epgamma_0.40GeV.root",
        f"{_BASE}/dvcsgen_rga_fa18_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/"
        "efficiency_study/"
        "nSidis_rga_fa18_inb_eppi0.root",
        f"{_BASE}/aaogen_rga_fa18_inb_eppi0_0.40GeV.root",
    ),
    Period(
        "fa18_out",
        "Fa18 Out",
        f"{_DATA_BASE}/nSidis_rga_fa18_out_epgamma.root",
        f"{_BASE}/bkg_rga_fa18_out_epgamma_0.40GeV.root",
        f"{_BASE}/dvcsgen_rga_fa18_out_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/"
        "efficiency_study/"
        "nSidis_rga_fa18_out_eppi0.root",
        f"{_BASE}/aaogen_rga_fa18_out_eppi0_0.40GeV.root",
    ),
    Period(
        "sp18_inb",
        "Sp18 Inb",
        f"{_DATA_BASE}/nSidis_rga_sp18_inb_epgamma.root",
        f"{_BASE}/bkg_rga_sp18_inb_epgamma_0.40GeV.root",
        f"{_BASE}/dvcsgen_rga_sp18_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/"
        "efficiency_study/"
        "nSidis_rga_sp18_inb_eppi0.root",
        f"{_BASE}/aaogen_rga_sp18_inb_eppi0_0.40GeV.root",
    ),
    Period(
        "sp18_out",
        "Sp18 Out",
        f"{_DATA_BASE}/nSidis_rga_sp18_out_epgamma.root",
        f"{_BASE}/bkg_rga_sp18_out_epgamma_0.40GeV.root",
        f"{_BASE}/dvcsgen_rga_sp18_out_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/"
        "efficiency_study/"
        "nSidis_rga_sp18_out_eppi0.root",
        f"{_BASE}/aaogen_rga_sp18_out_eppi0_0.40GeV.root",
    ),
    Period(
        "sp19_inb",
        "Sp19 Inb",
        f"{_DATA_BASE}/nSidis_rga_sp19_inb_epgamma.root",
        f"{_BASE}/bkg_rga_sp19_inb_epgamma_0.40GeV.root",
        f"{_BASE}/dvcsgen_rga_sp19_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/"
        "efficiency_study/"
        "nSidis_rga_sp19_inb_eppi0.root",
        f"{_BASE}/aaogen_rga_sp19_inb_eppi0_0.40GeV.root",
    ),
)


SAMPLES = (
    ("data", "data", "black"),
    ("pi0", r"exclusive $\pi^0$ (AAO)", "tab:red"),
    ("dvcs", "BH/DVCS (DVCSgen)", "tab:blue"),
)


# Reconstructed variables available to the first rectangular-cut optimizer.
# The optimizer is run separately for every period and independently in FT/FD.
OPTIMIZER_FEATURES: Tuple[str, ...] = (
    "Mx2",
    "Mx2_1",
    "Mx2_2",
    "Emiss2",
    "E_tag",
    "Delta_phi_residual_deg",
    "pTmiss",
    "theta_gamma_gamma",
)

# Mx2_1 defines the loose exclusivity preselection/data sideband, while
# Emiss2 defines the E_probe bin. Both are deliberately excluded from the
# greedy threshold scan.
OPTIMIZER_SCAN_FEATURES: Tuple[str, ...] = tuple(
    feature
    for feature in OPTIMIZER_FEATURES
    if feature not in ("Mx2_1", "Emiss2")
)

# Every remaining optimizer variable may use either one-sided direction.
# Emiss2/E_probe is not present here because it defines the energy bin itself.
OPTIMIZER_ALLOWED_DIRECTIONS = {
    feature: ("lt", "gt")
    for feature in OPTIMIZER_SCAN_FEATURES
}


# Mx2_1 defines signal versus sideband and Emiss2 defines the E_probe bin,
# therefore the real-data composition fit uses only these remaining variables.
DATA_FIT_FEATURES: Tuple[str, ...] = tuple(OPTIMIZER_SCAN_FEATURES)

# Extra event metadata retained in the same reservoir. It is never scanned by
# the optimizer and exists only for downstream detector diagnostics.
EVENT_METADATA: Tuple[str, ...] = (
    "tag_theta_deg",
    "pred_probe_sector",
    "pred_probe_theta_deg",
    "pred_probe_phi_deg",
)
EVENT_COLUMNS: Tuple[str, ...] = OPTIMIZER_FEATURES + EVENT_METADATA


# Additional branches needed only during the final tag-probe efficiency pass.
EFF_EPG_REQUIRED: Tuple[str, ...] = (
    "e_p", "e_theta", "e_phi",
    "p1_p", "p1_theta", "p1_phi",
    "p2_p", "p2_theta", "p2_phi",
    "Mx2", "Mx2_1", "Mx2_2",
    "Emiss2", "Delta_phi", "pTmiss", "theta_gamma_gamma",
)

EFF_EPG_DATA_IDS: Tuple[str, ...] = ("runnum", "evnum")

EFF_EPPIO_REQUIRED: Tuple[str, ...] = (
    "e_p", "e_theta", "e_phi",
    "p1_p", "p1_theta", "p1_phi",
    "p2_p", "p2_theta", "p2_phi",
    "Mh_gammagamma",
)

EFF_EPPIO_DATA_IDS: Tuple[str, ...] = ("runnum", "evnum")


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
        "--coverage-only",
        action="store_true",
        help=(
            "Run only the raw full-file nSidis epgamma/eppi0 run-event "
            "coverage audit for the requested periods, write the coverage "
            "plots, and exit."
        ),
    )
    parser.add_argument(
        "--coverage-min-eppi0-overlap",
        type=float,
        default=0.50,
        help=(
            "Minimum fraction of unique eppi0 (runnum,evnum) keys that must "
            "also occur in the corresponding epgamma file before the full "
            "efficiency stage is allowed to proceed. Default: 0.50."
        ),
    )
    parser.add_argument(
        "--skip-efficiency",
        action="store_true",
        help=(
            "Run the purity/optimizer/template study only and skip the final "
            "epgamma<->eppi0 tag-probe efficiency extraction."
        ),
    )
    parser.add_argument(
        "--efficiency-eppi0-max-entries",
        type=int,
        default=0,
        help=(
            "Maximum entries loaded from each companion eppi0 file for the "
            "efficiency stage. 0 means the full eppi0 tree (recommended even "
            "when --max-entries limits epgamma, so cross-file matching is not "
            "artificially truncated). Default: 0."
        ),
    )
    parser.add_argument(
        "--parent-component-tol",
        type=float,
        default=0.002,
        help=(
            "Maximum absolute reconstructed e/p Cartesian-component mismatch "
            "(GeV) used for AAO epgamma<->eppi0 matching and as the data-parent "
            "consistency check. Default: 0.002."
        ),
    )
    parser.add_argument(
        "--parent-distance-max",
        type=float,
        default=2.0,
        help=(
            "Maximum six-dimensional cKDTree distance after scaling each e/p "
            "component by --parent-component-tol. Used only for MC. Default: 2."
        ),
    )
    parser.add_argument(
        "--kdtree-workers",
        type=int,
        default=1,
        help=(
            "Threads used inside each MC cKDTree query. Keep at 1 when periods "
            "are processed in parallel to avoid oversubscription. Default: 1."
        ),
    )
    parser.add_argument(
        "--assoc-mgg-min",
        type=float,
        default=0.11,
        help="Lower reconstructed pi0 mass for a clean tag association. Default: 0.11 GeV.",
    )
    parser.add_argument(
        "--assoc-mgg-max",
        type=float,
        default=0.16,
        help="Upper reconstructed pi0 mass for a clean tag association. Default: 0.16 GeV.",
    )
    parser.add_argument(
        "--assoc-remainder-mass2-max",
        type=float,
        default=1.0e-3,
        help=(
            "Require |(P_pi0-k_tag)^2| below this value for a clean "
            "reconstructed-side association. Default: 1e-3 GeV^2."
        ),
    )
    parser.add_argument(
        "--assoc-probe-energy-min",
        type=float,
        default=0.40,
        help="Minimum reconstructed companion-photon energy. Default: 0.40 GeV.",
    )
    parser.add_argument(
        "--probe-match-nsigma",
        type=float,
        default=3.0,
        help=(
            "Final |Delta p_i-center_i| < N sigma_i probe-match window, with "
            "centers/sigmas derived from AAO. Default: 3."
        ),
    )
    parser.add_argument(
        "--probe-match-min-resolution-events",
        type=int,
        default=200,
        help=(
            "Minimum clean AAO pairs required to derive a separate FT or FD "
            "Delta-p resolution. Otherwise that region falls back to the "
            "period-wide AAO residual model. Default: 200."
        ),
    )
    parser.add_argument(
        "--linear",
        action="store_true",
        help=(
            "Use a linear y-axis for the shape-comparison canvases. "
            "Default: logarithmic y-axis."
        ),
    )
    parser.add_argument(
        "--no-optimize-cuts",
        action="store_true",
        help=(
            "Disable the iterative rectangular-cut optimizer. "
            "Default: run it independently for every selected period and "
            "separately in FT and FD."
        ),
    )
    parser.add_argument(
        "--optimizer-max-events",
        type=int,
        default=200_000,
        help=(
            "Maximum events retained in memory per sample, period, detector region, "
            "and E_probe bin for cut optimization. A uniform "
            "priority reservoir is used while streaming. Default: 200000."
        ),
    )
    parser.add_argument(
        "--optimizer-steps",
        type=int,
        default=4,
        help="Maximum number of greedy rectangular cuts. Default: 4.",
    )
    parser.add_argument(
        "--optimizer-quantiles",
        type=int,
        default=80,
        help=(
            "Number of quantile positions scanned per variable and direction "
            "at each optimizer step. Default: 80."
        ),
    )
    parser.add_argument(
        "--optimizer-min-step-pi0-eff",
        type=float,
        default=0.75,
        help=(
            "Each newly added cut must retain at least this fraction of the "
            "AAO events surviving the previous step. Default: 0.75."
        ),
    )
    parser.add_argument(
        "--optimizer-min-total-pi0-eff",
        type=float,
        default=0.20,
        help=(
            "Never accept a cut that reduces cumulative AAO retention below "
            "this fraction of the baseline. Default: 0.20."
        ),
    )
    parser.add_argument(
        "--optimizer-min-improvement",
        type=float,
        default=1.02,
        help=(
            "Stop if the best new cut improves the current separation score "
            "by less than this multiplicative factor. Default: 1.02."
        ),
    )
    parser.add_argument(
        "--optimizer-mx2-1-preselection-max",
        type=float,
        default=SIGNAL_MX2_1_MAX,
        help=(
            "Signal-region exclusivity preselection Mx2_1 < value in GeV^2, "
            "applied before AAO/DVCS optimization. Default: 0.15."
        ),
    )
    parser.add_argument(
        "--optimizer-data-sideband-min",
        type=float,
        default=NEAR_SIDEBAND_MIN,
        help=(
            "Lower Mx2_1 edge of the empirical nonexclusive background proxy. "
            "Default: 0.15 GeV^2."
        ),
    )
    parser.add_argument(
        "--optimizer-data-sideband-max",
        type=float,
        default=FAR_SIDEBAND_MAX,
        help=(
            "Upper Mx2_1 edge of the empirical background proxy used by the "
            "optimizer. Default: 0.23 GeV^2, deliberately below eta."
        ),
    )
    parser.add_argument(
        "--optimizer-data-sideband-weight",
        type=float,
        default=1.0,
        help=(
            "Relative weight lambda of data-sideband retention in the "
            "optimizer score denominator. Default: 1.0."
        ),
    )
    parser.add_argument(
        "--closure-toys",
        type=int,
        default=100,
        help=(
            "Pseudoexperiments per true pi0 fraction and cumulative optimizer "
            "step for template-fit closure. Default: 100."
        ),
    )
    parser.add_argument(
        "--closure-events",
        type=int,
        default=5000,
        help=(
            "Reference pseudo-data event count used in each closure toy. "
            "This fixes the statistical scale for comparing identifiability "
            "between cut steps. Default: 5000."
        ),
    )
    parser.add_argument(
        "--closure-bins",
        type=int,
        default=24,
        help=(
            "Maximum quantile bins for the automatically selected 1D closure "
            "discriminator. Default: 24."
        ),
    )
    parser.add_argument(
        "--data-fit-max-closure-bias",
        type=float,
        default=0.01,
        help=(
            "Maximum allowed absolute closure bias at injected f_pi0=0.5 "
            "when selecting the real-data fit operating step. Default: 0.01."
        ),
    )
    parser.add_argument(
        "--data-fit-max-closure-rms",
        type=float,
        default=0.03,
        help=(
            "Absolute ceiling on closure RMS when selecting the real-data fit "
            "operating step. Default: 0.03."
        ),
    )
    parser.add_argument(
        "--data-fit-background-near-best",
        type=float,
        default=1.10,
        help=(
            "Among closure-safe cut steps, choose the earliest step whose "
            "DVCS+data-sideband retention is within this factor of the best "
            "available value. Default: 1.10."
        ),
    )
    parser.add_argument(
        "--data-fit-morph-shift-max-bins",
        type=float,
        default=1.5,
        help=(
            "Maximum absolute AAO centroid shift in fixed fit-bin widths. "
            "Default: 1.5."
        ),
    )
    parser.add_argument(
        "--data-fit-morph-smear-max-bins",
        type=float,
        default=1.5,
        help=(
            "Maximum extra Gaussian AAO smearing in fixed fit-bin widths. "
            "Default: 1.5."
        ),
    )
    parser.add_argument(
        "--data-fit-morph-shift-steps",
        type=int,
        default=13,
        help="Number of bounded AAO shift values scanned. Default: 13.",
    )
    parser.add_argument(
        "--data-fit-morph-smear-steps",
        type=int,
        default=9,
        help="Number of bounded AAO smearing values scanned. Default: 9.",
    )
    parser.add_argument(
        "--data-fit-sideband-min-events",
        type=int,
        default=100,
        help=(
            "Minimum in-range events required independently in the eta-safe "
            "near and far sidebands before the empirical sideband is allowed "
            "as a third fit template. Default: 100."
        ),
    )
    parser.add_argument(
        "--data-fit-sideband-max-js",
        type=float,
        default=0.05,
        help=(
            "Maximum allowed Jensen-Shannon divergence between the eta-safe "
            "near (0.15-0.19) and far (0.19-0.23) sideband shapes before "
            "using the near sideband as a transported third template. "
            "Default: 0.05."
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
) -> dict:
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



def _angles_to_radians(
    theta: np.ndarray,
    phi: np.ndarray,
    unit: str,
) -> Tuple[np.ndarray, np.ndarray]:
    """Convert stored polar/azimuthal angles to radians."""
    theta = np.asarray(theta, dtype=float)
    phi = np.asarray(phi, dtype=float)

    if unit == "deg":
        return np.deg2rad(theta), np.deg2rad(phi)
    #endif

    return theta, phi


def predicted_probe_direction_deg(
    arrays: Dict[str, np.ndarray],
    angle_unit: str,
    beam_energy_gev: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Predict the missing probe-photon direction from the reconstructed
    e + p + tag-gamma system.

    With the target proton at rest, the missing three-momentum is

        p_probe^pred = p_beam - p_e - p_p - p_gamma,tag.

    The incoming electron beam is taken along +z with |p_beam| ~= E_beam.
    Only the direction is needed here, so the target mass does not enter.

    Returns:
        theta_pred_deg, phi_pred_deg

    with phi wrapped to [0, 360).
    """
    e_theta, e_phi = _angles_to_radians(
        arrays["e_theta"],
        arrays["e_phi"],
        angle_unit,
    )
    p_theta, p_phi = _angles_to_radians(
        arrays["p1_theta"],
        arrays["p1_phi"],
        angle_unit,
    )
    g_theta, g_phi = _angles_to_radians(
        arrays["p2_theta"],
        arrays["p2_phi"],
        angle_unit,
    )

    e_p = np.asarray(arrays["e_p"], dtype=float)
    p_p = np.asarray(arrays["p1_p"], dtype=float)
    g_p = np.asarray(arrays["p2_p"], dtype=float)

    e_px = e_p * np.sin(e_theta) * np.cos(e_phi)
    e_py = e_p * np.sin(e_theta) * np.sin(e_phi)
    e_pz = e_p * np.cos(e_theta)

    p_px = p_p * np.sin(p_theta) * np.cos(p_phi)
    p_py = p_p * np.sin(p_theta) * np.sin(p_phi)
    p_pz = p_p * np.cos(p_theta)

    g_px = g_p * np.sin(g_theta) * np.cos(g_phi)
    g_py = g_p * np.sin(g_theta) * np.sin(g_phi)
    g_pz = g_p * np.cos(g_theta)

    pred_px = -(e_px + p_px + g_px)
    pred_py = -(e_py + p_py + g_py)
    pred_pz = (
        float(beam_energy_gev)
        - e_pz
        - p_pz
        - g_pz
    )

    pred_pt = np.hypot(pred_px, pred_py)
    pred_p = np.sqrt(
        pred_px * pred_px
        + pred_py * pred_py
        + pred_pz * pred_pz
    )

    theta_deg = np.rad2deg(
        np.arctan2(pred_pt, pred_pz)
    )
    phi_deg = np.mod(
        np.rad2deg(np.arctan2(pred_py, pred_px)),
        360.0,
    )

    valid = (
        np.isfinite(pred_p)
        & (pred_p > 1.0e-9)
        & np.isfinite(theta_deg)
        & np.isfinite(phi_deg)
    )

    theta_deg = np.where(valid, theta_deg, np.nan)
    phi_deg = np.where(valid, phi_deg, np.nan)

    return theta_deg, phi_deg


def fd_sector_from_phi_deg(phi_deg: np.ndarray) -> np.ndarray:
    """
    CLAS12 FD sector assigned from the PREDICTED probe-photon azimuth.

      S1: 330-360 or   0-30 deg
      S2:  30-90 deg
      S3:  90-150 deg
      S4: 150-210 deg
      S5: 210-270 deg
      S6: 270-330 deg

    Non-finite entries receive sector 0.
    """
    phi = np.asarray(phi_deg, dtype=float)
    wrapped = np.mod(phi, 360.0)
    sector = np.zeros(phi.shape, dtype=np.float32)

    finite = np.isfinite(wrapped)
    sector[
        finite
        & ((wrapped >= 330.0) | (wrapped < 30.0))
    ] = 1.0
    sector[
        finite
        & (wrapped >= 30.0)
        & (wrapped < 90.0)
    ] = 2.0
    sector[
        finite
        & (wrapped >= 90.0)
        & (wrapped < 150.0)
    ] = 3.0
    sector[
        finite
        & (wrapped >= 150.0)
        & (wrapped < 210.0)
    ] = 4.0
    sector[
        finite
        & (wrapped >= 210.0)
        & (wrapped < 270.0)
    ] = 5.0
    sector[
        finite
        & (wrapped >= 270.0)
        & (wrapped < 330.0)
    ] = 6.0

    return sector


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
# Rectangular-cut optimizer support
# =============================================================================

class PriorityReservoir:
    """
    Uniform fixed-size event sample maintained with random priorities.

    Every accepted event receives an independent U(0,1) priority. We retain
    only the K smallest priorities seen so far. This is a vectorized uniform
    reservoir across all streamed chunks and avoids keeping the full event
    sample in memory.
    """

    def __init__(self, capacity: int, n_features: int, seed: int) -> None:
        self.capacity = int(capacity)
        self.n_features = int(n_features)
        self.rng = np.random.default_rng(int(seed))
        self.priorities = np.empty(0, dtype=np.float32)
        self.values = np.empty((0, n_features), dtype=np.float32)

    def update(self, values: np.ndarray) -> None:
        if self.capacity <= 0 or values.size == 0:
            return
        #endif

        values = np.asarray(values, dtype=np.float32)
        finite = np.all(np.isfinite(values), axis=1)
        values = values[finite]
        if values.shape[0] == 0:
            return
        #endif

        priorities = self.rng.random(values.shape[0]).astype(np.float32)

        if self.values.shape[0] == 0:
            combined_values = values
            combined_priorities = priorities
        else:
            combined_values = np.concatenate((self.values, values), axis=0)
            combined_priorities = np.concatenate(
                (self.priorities, priorities),
                axis=0,
            )
        #endif

        if combined_values.shape[0] <= self.capacity:
            self.values = combined_values
            self.priorities = combined_priorities
            return
        #endif

        keep = np.argpartition(
            combined_priorities,
            self.capacity - 1,
        )[: self.capacity]
        self.values = combined_values[keep]
        self.priorities = combined_priorities[keep]

    def array(self) -> np.ndarray:
        return np.asarray(self.values, dtype=np.float32)


def optimizer_matrix(values: Dict[str, np.ndarray], mask: np.ndarray) -> np.ndarray:
    """
    Build the reservoir event matrix. Physics variables come first, followed by
    detector metadata that the optimizer never scans.
    """
    if not np.any(mask):
        return np.empty((0, len(EVENT_COLUMNS)), dtype=np.float32)
    #endif

    return np.column_stack(
        [
            np.asarray(values[key][mask], dtype=np.float32)
            for key in EVENT_COLUMNS
        ]
    )


def _candidate_thresholds(
    pi_values: np.ndarray,
    dvcs_values: np.ndarray,
    n_quantiles: int,
) -> np.ndarray:
    """
    Quantile-based thresholds avoid wasting scan points in empty phase space.
    """
    combined = np.concatenate((pi_values, dvcs_values))
    combined = combined[np.isfinite(combined)]

    if combined.size < 10:
        return np.empty(0, dtype=float)
    #endif

    q = np.linspace(0.02, 0.98, max(3, int(n_quantiles)))
    thresholds = np.quantile(combined, q)
    return np.unique(np.asarray(thresholds, dtype=float))


def _candidate_thresholds_three(
    pi_values: np.ndarray,
    dvcs_values: np.ndarray,
    side_values: np.ndarray,
    n_quantiles: int,
) -> np.ndarray:
    """Quantile thresholds from surviving AAO, DVCSgen, and data sideband."""
    combined = np.concatenate((pi_values, dvcs_values, side_values))
    combined = combined[np.isfinite(combined)]

    if combined.size < 10:
        return np.empty(0, dtype=float)
    #endif

    q = np.linspace(0.02, 0.98, max(3, int(n_quantiles)))
    return np.unique(np.asarray(np.quantile(combined, q), dtype=float))


def optimize_rectangular_cuts(
    pi0_events: np.ndarray,
    dvcs_events: np.ndarray,
    data_sideband_events: np.ndarray,
    max_steps: int,
    n_quantiles: int,
    min_step_pi0_eff: float,
    min_total_pi0_eff: float,
    min_improvement: float,
    data_sideband_weight: float,
) -> list:
    """
    Greedy one-sided optimizer balancing AAO retention against two backgrounds.

    At each step it scans one-sided rectangular thresholds over
    OPTIMIZER_SCAN_FEATURES and maximizes

        F = eps_pi0 / sqrt(
            eps_DVCS + lambda * eps_data_sideband + 1e-9
        )

    using cumulative efficiencies relative to the corresponding baseline
    populations.

    Mx2_1 is excluded from the scan because it defines both the loose
    exclusivity preselection and the real-data sideband.

    Emiss2 is also excluded from the scan because it defines the E_probe bin.
    Thus the optimizer cannot improve purity by discarding part of the
    requested photon-energy interval.
    """
    if (
        pi0_events.shape[0] == 0
        or dvcs_events.shape[0] == 0
        or data_sideband_events.shape[0] == 0
    ):
        return []
    #endif

    pi_mask = np.ones(pi0_events.shape[0], dtype=bool)
    dvcs_mask = np.ones(dvcs_events.shape[0], dtype=bool)
    side_mask = np.ones(data_sideband_events.shape[0], dtype=bool)

    n_pi0_0 = int(pi0_events.shape[0])
    n_dvcs_0 = int(dvcs_events.shape[0])
    n_side_0 = int(data_sideband_events.shape[0])

    current_score = 1.0 / math.sqrt(1.0 + float(data_sideband_weight))
    used_features = set()
    results = []

    for step in range(1, int(max_steps) + 1):
        n_pi_before = int(np.count_nonzero(pi_mask))
        n_dvcs_before = int(np.count_nonzero(dvcs_mask))
        n_side_before = int(np.count_nonzero(side_mask))

        if min(n_pi_before, n_dvcs_before, n_side_before) == 0:
            break
        #endif

        best = None

        for feature in OPTIMIZER_SCAN_FEATURES:
            if feature in used_features:
                continue
            #endif

            ifeature = FEATURE_INDEX[feature]
            thresholds = _candidate_thresholds_three(
                pi0_events[pi_mask, ifeature],
                dvcs_events[dvcs_mask, ifeature],
                data_sideband_events[side_mask, ifeature],
                n_quantiles,
            )

            for threshold in thresholds:
                for direction in OPTIMIZER_ALLOWED_DIRECTIONS[feature]:
                    if direction == "lt":
                        pi_local = pi0_events[:, ifeature] < threshold
                        dvcs_local = dvcs_events[:, ifeature] < threshold
                        side_local = data_sideband_events[:, ifeature] < threshold
                    else:
                        pi_local = pi0_events[:, ifeature] > threshold
                        dvcs_local = dvcs_events[:, ifeature] > threshold
                        side_local = data_sideband_events[:, ifeature] > threshold
                    #endif

                    pi_test = pi_mask & pi_local
                    dvcs_test = dvcs_mask & dvcs_local
                    side_test = side_mask & side_local

                    n_pi = int(np.count_nonzero(pi_test))
                    n_dvcs = int(np.count_nonzero(dvcs_test))
                    n_side = int(np.count_nonzero(side_test))

                    step_pi_eff = n_pi / n_pi_before
                    step_dvcs_eff = n_dvcs / n_dvcs_before
                    step_side_eff = n_side / n_side_before

                    total_pi_eff = n_pi / n_pi0_0
                    total_dvcs_eff = n_dvcs / n_dvcs_0
                    total_side_eff = n_side / n_side_0

                    if step_pi_eff < min_step_pi0_eff:
                        continue
                    #endif
                    if total_pi_eff < min_total_pi0_eff:
                        continue
                    #endif

                    score = total_pi_eff / math.sqrt(
                        total_dvcs_eff
                        + float(data_sideband_weight) * total_side_eff
                        + 1.0e-9
                    )

                    if best is None or score > best["score"]:
                        best = {
                            "step": step,
                            "feature": feature,
                            "direction": direction,
                            "threshold": float(threshold),
                            "score": float(score),
                            "step_pi0_eff": float(step_pi_eff),
                            "step_dvcs_eff": float(step_dvcs_eff),
                            "step_data_sideband_eff": float(step_side_eff),
                            "total_pi0_eff": float(total_pi_eff),
                            "total_dvcs_eff": float(total_dvcs_eff),
                            "total_data_sideband_eff": float(total_side_eff),
                            "n_pi0": n_pi,
                            "n_dvcs": n_dvcs,
                            "n_data_sideband": n_side,
                            "pi_mask": pi_test,
                            "dvcs_mask": dvcs_test,
                            "side_mask": side_test,
                        }
                    #endif
                #endfor
            #endfor
        #endfor

        if best is None:
            break
        #endif

        improvement = best["score"] / max(current_score, 1.0e-12)
        if improvement < min_improvement:
            break
        #endif

        results.append(
            {
                key: value
                for key, value in best.items()
                if key not in ("pi_mask", "dvcs_mask", "side_mask")
            }
        )

        pi_mask = best["pi_mask"]
        dvcs_mask = best["dvcs_mask"]
        side_mask = best["side_mask"]
        current_score = best["score"]
        used_features.add(best["feature"])
    #endfor

    return results


def cut_expression(result: dict) -> str:
    op = "<" if result["direction"] == "lt" else ">"
    return f"{result['feature']} {op} {result['threshold']:.6g}"



# =============================================================================
# Template-fit closure / identifiability study
# =============================================================================

CLOSURE_TRUE_FRACTIONS: Tuple[float, ...] = (0.20, 0.50, 0.80)


def jensen_shannon_divergence(
    p: np.ndarray,
    q: np.ndarray,
) -> float:
    """Jensen-Shannon divergence in natural-log units."""
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)

    p = p / np.sum(p)
    q = q / np.sum(q)
    m = 0.5 * (p + q)

    return float(
        0.5 * np.sum(p * np.log(p / m))
        + 0.5 * np.sum(q * np.log(q / m))
    )


def closure_histogram_probabilities(
    values: np.ndarray,
    edges: np.ndarray,
    pseudocount: float = 0.5,
) -> np.ndarray:
    """
    Histogram probability vector with a small Jeffreys-like pseudocount so
    finite MC statistics never produce impossible zero-probability bins.
    """
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    counts, _ = np.histogram(values, bins=edges)
    probs = counts.astype(float) + float(pseudocount)
    return probs / np.sum(probs)


def closure_quantile_edges(
    pi_values: np.ndarray,
    dvcs_values: np.ndarray,
    n_bins: int,
) -> np.ndarray:
    """
    Quantile bins from the combined training templates.

    Quantile binning gives approximately useful occupancy everywhere and is
    more stable than a fixed-width histogram for long-tailed variables.
    """
    combined = np.concatenate((pi_values, dvcs_values))
    combined = combined[np.isfinite(combined)]

    if combined.size < 20:
        return np.empty(0, dtype=float)
    #endif

    edges = np.quantile(
        combined,
        np.linspace(0.0, 1.0, max(4, int(n_bins)) + 1),
    )
    edges = np.unique(np.asarray(edges, dtype=float))

    if edges.size < 5:
        return np.empty(0, dtype=float)
    #endif

    # Include the extrema robustly in np.histogram.
    edges[0] = np.nextafter(edges[0], -np.inf)
    edges[-1] = np.nextafter(edges[-1], np.inf)
    return edges


def fit_two_template_fraction(
    counts: np.ndarray,
    p_pi0: np.ndarray,
    p_dvcs: np.ndarray,
) -> Tuple[float, float]:
    """
    Maximum-likelihood fit of

        p_i(f) = f p_pi0,i + (1-f) p_dvcs,i

    to multinomial-binned pseudo-data.

    The score equation is monotonic in f, so a bounded bisection gives the
    exact one-parameter maximum without scipy.  The returned uncertainty is
    the local inverse-curvature estimate.
    """
    counts = np.asarray(counts, dtype=float)
    p_pi0 = np.asarray(p_pi0, dtype=float)
    p_dvcs = np.asarray(p_dvcs, dtype=float)
    delta = p_pi0 - p_dvcs

    def gradient(f_value: float) -> float:
        mixture = p_dvcs + f_value * delta
        mixture = np.clip(mixture, 1.0e-15, None)
        return float(np.sum(counts * delta / mixture))

    g0 = gradient(0.0)
    g1 = gradient(1.0)

    if g0 <= 0.0:
        f_hat = 0.0
    elif g1 >= 0.0:
        f_hat = 1.0
    else:
        lo = 0.0
        hi = 1.0
        for _ in range(50):
            mid = 0.5 * (lo + hi)
            if gradient(mid) > 0.0:
                lo = mid
            else:
                hi = mid
            #endif
        #endfor
        f_hat = 0.5 * (lo + hi)
    #endif

    mixture = p_dvcs + f_hat * delta
    mixture = np.clip(mixture, 1.0e-15, None)
    information = float(
        np.sum(counts * (delta * delta) / (mixture * mixture))
    )
    sigma = (
        1.0 / math.sqrt(information)
        if information > 0.0
        else float("inf")
    )

    return float(f_hat), float(sigma)


def build_cumulative_optimizer_masks(
    events: np.ndarray,
    optimizer_results: list,
) -> list:
    """Baseline plus one cumulative mask after each optimizer-selected cut."""
    masks = [np.ones(events.shape[0], dtype=bool)]
    current = masks[0]

    for result in optimizer_results:
        current = apply_optimizer_cut(events, current, result)
        masks.append(current.copy())
    #endfor

    return masks


def run_template_fit_closure(
    pi0_events: np.ndarray,
    dvcs_events: np.ndarray,
    optimizer_results: list,
    n_toys: int,
    n_pseudodata_events: int,
    n_bins: int,
    seed: int,
) -> list:
    """
    Test AAO/DVCS mixture identifiability at every cumulative cut step.

    Procedure at each step:
      1. Apply all cuts through that step.
      2. Split surviving AAO and DVCS events into independent template-training
         and held-out pseudo-data halves.
      3. Among the allowed reconstructed variables, find the single 1D
         variable with the largest AAO-vs-DVCS Jensen-Shannon divergence in the
         training templates.
      4. Build pseudo-data from the held-out shapes at true pi0 fractions
         0.20, 0.50, and 0.80.
      5. Fit each pseudo-data sample with the independent training templates.

    Large fit RMS/bias after later cuts means that purity optimization has made
    the residual AAO and DVCS templates too degenerate for a reliable fraction
    extraction.
    """
    if (
        pi0_events.shape[0] < 40
        or dvcs_events.shape[0] < 40
        or n_toys <= 0
        or n_pseudodata_events <= 0
    ):
        return []
    #endif

    rng = np.random.default_rng(int(seed))

    # Split once before cut application so the training/held-out samples remain
    # statistically independent at every cumulative optimizer step.
    pi_perm = rng.permutation(pi0_events.shape[0])
    dvcs_perm = rng.permutation(dvcs_events.shape[0])

    pi_train_flag = np.zeros(pi0_events.shape[0], dtype=bool)
    dvcs_train_flag = np.zeros(dvcs_events.shape[0], dtype=bool)
    pi_train_flag[pi_perm[: pi_perm.size // 2]] = True
    dvcs_train_flag[dvcs_perm[: dvcs_perm.size // 2]] = True

    pi_test_flag = ~pi_train_flag
    dvcs_test_flag = ~dvcs_train_flag

    pi_masks = build_cumulative_optimizer_masks(
        pi0_events,
        optimizer_results,
    )
    dvcs_masks = build_cumulative_optimizer_masks(
        dvcs_events,
        optimizer_results,
    )

    closure_results = []

    for step_index, (pi_mask, dvcs_mask) in enumerate(
        zip(pi_masks, dvcs_masks)
    ):
        pi_train = pi_mask & pi_train_flag
        pi_test = pi_mask & pi_test_flag
        dvcs_train = dvcs_mask & dvcs_train_flag
        dvcs_test = dvcs_mask & dvcs_test_flag

        n_pi_train = int(np.count_nonzero(pi_train))
        n_pi_test = int(np.count_nonzero(pi_test))
        n_dvcs_train = int(np.count_nonzero(dvcs_train))
        n_dvcs_test = int(np.count_nonzero(dvcs_test))

        if min(n_pi_train, n_pi_test, n_dvcs_train, n_dvcs_test) < 20:
            continue
        #endif

        best_feature = None
        best_edges = None
        best_p_pi_train = None
        best_p_dvcs_train = None
        best_js = -1.0

        # Emiss2 and Mx2_1 are still valid closure diagnostics even though they
        # are not optimizer variables.  The closure study asks only whether
        # any surviving reconstructed 1D shape retains normalization leverage.
        for feature in OPTIMIZER_FEATURES:
            ifeature = FEATURE_INDEX[feature]

            pi_values = pi0_events[pi_train, ifeature]
            dvcs_values = dvcs_events[dvcs_train, ifeature]

            edges = closure_quantile_edges(
                pi_values,
                dvcs_values,
                n_bins,
            )
            if edges.size == 0:
                continue
            #endif

            p_pi_train = closure_histogram_probabilities(
                pi_values,
                edges,
            )
            p_dvcs_train = closure_histogram_probabilities(
                dvcs_values,
                edges,
            )

            js = jensen_shannon_divergence(
                p_pi_train,
                p_dvcs_train,
            )

            if js > best_js:
                best_js = js
                best_feature = feature
                best_edges = edges
                best_p_pi_train = p_pi_train
                best_p_dvcs_train = p_dvcs_train
            #endif
        #endfor

        if best_feature is None:
            continue
        #endif

        ifeature = FEATURE_INDEX[best_feature]
        p_pi_test = closure_histogram_probabilities(
            pi0_events[pi_test, ifeature],
            best_edges,
        )
        p_dvcs_test = closure_histogram_probabilities(
            dvcs_events[dvcs_test, ifeature],
            best_edges,
        )

        for true_fraction in CLOSURE_TRUE_FRACTIONS:
            truth_probability = (
                true_fraction * p_pi_test
                + (1.0 - true_fraction) * p_dvcs_test
            )
            truth_probability /= np.sum(truth_probability)

            fitted = np.empty(n_toys, dtype=float)
            fitted_sigma = np.empty(n_toys, dtype=float)

            for itoy in range(n_toys):
                toy_counts = rng.multinomial(
                    int(n_pseudodata_events),
                    truth_probability,
                )
                fitted[itoy], fitted_sigma[itoy] = (
                    fit_two_template_fraction(
                        toy_counts,
                        best_p_pi_train,
                        best_p_dvcs_train,
                    )
                )
            #endfor

            delta = fitted - float(true_fraction)

            closure_results.append(
                {
                    "step": int(step_index),
                    "feature": str(best_feature),
                    "js_divergence": float(best_js),
                    "true_fraction": float(true_fraction),
                    "mean_fit": float(np.mean(fitted)),
                    "bias": float(np.mean(delta)),
                    "rms_error": float(np.sqrt(np.mean(delta * delta))),
                    "std_fit": float(np.std(fitted, ddof=1)),
                    "mean_sigma": float(np.mean(fitted_sigma)),
                    "boundary_fraction": float(
                        np.mean(
                            (fitted <= 1.0e-6)
                            | (fitted >= 1.0 - 1.0e-6)
                        )
                    ),
                    "n_pi0_train": n_pi_train,
                    "n_pi0_test": n_pi_test,
                    "n_dvcs_train": n_dvcs_train,
                    "n_dvcs_test": n_dvcs_test,
                    "n_toys": int(n_toys),
                    "n_pseudodata_events": int(n_pseudodata_events),
                }
            )
        #endfor
    #endfor

    return closure_results


def write_template_fit_closure_outputs(
    period: Period,
    region: str,
    energy_bin_key: str,
    energy_min: float,
    energy_max: float,
    closure_results: list,
    output_dir: Path,
) -> None:
    """Write detailed closure results for one period/region/E_probe bin."""
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = (
        f"template_fit_closure_{period.key}_{region.lower()}_"
        f"eprobe_{energy_bin_key}"
    )

    csv_path = output_dir / f"{stem}.csv"
    txt_path = output_dir / f"{stem}.txt"

    fields = (
        "step",
        "feature",
        "js_divergence",
        "true_fraction",
        "mean_fit",
        "bias",
        "rms_error",
        "std_fit",
        "mean_sigma",
        "boundary_fraction",
        "n_pi0_train",
        "n_pi0_test",
        "n_dvcs_train",
        "n_dvcs_test",
        "n_toys",
        "n_pseudodata_events",
    )

    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for result in closure_results:
            writer.writerow({field: result[field] for field in fields})
        #endfor
    #endwith

    lines = [
        (
            f"{period.label} {region}: template-fit closure, "
            f"{energy_min:g} <= E_probe < {energy_max:g} GeV"
        ),
        "=" * 88,
        (
            "At each cumulative optimizer step, the single reconstructed "
            "variable with maximum AAO-vs-DVCS Jensen-Shannon divergence is "
            "used for an independent train/held-out two-template closure fit."
        ),
        "",
    ]

    if not closure_results:
        lines.append("Insufficient statistics for closure.")
    else:
        steps = sorted({int(result["step"]) for result in closure_results})
        for step in steps:
            step_rows = [
                result
                for result in closure_results
                if int(result["step"]) == step
            ]
            if not step_rows:
                continue
            #endif

            lines.append(
                f"Step {step}: best feature={step_rows[0]['feature']}, "
                f"JS={step_rows[0]['js_divergence']:.5f}"
            )

            for result in step_rows:
                lines.append(
                    f"  f_true={result['true_fraction']:.2f}: "
                    f"<f_fit>={result['mean_fit']:.4f}, "
                    f"bias={result['bias']:+.4f}, "
                    f"RMS={result['rms_error']:.4f}, "
                    f"<sigma_fit>={result['mean_sigma']:.4f}, "
                    f"boundary={100.0*result['boundary_fraction']:.1f}%"
                )
            #endfor
            lines.append("")
        #endfor
    #endif

    txt_path.write_text("\n".join(lines) + "\n")


def make_combined_closure_canvas(
    period: Period,
    optimizer_summary: dict,
    output_dir: Path,
) -> Path:
    """
    One 4x2 closure canvas per period.

    Each panel shows RMS(f_fit-f_true) versus cumulative cut step for the three
    injected pi0 fractions.  Growth of the RMS after later cuts is the direct
    warning that template identifiability is being lost.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(
        len(E_PROBE_BINS),
        2,
        figsize=(13.8, 4.0 * len(E_PROBE_BINS) + 1.0),
        squeeze=False,
    )

    legend_handles = []
    legend_labels = []

    for irow, (bin_key, e_min, e_max) in enumerate(E_PROBE_BINS):
        bin_summary = optimizer_summary.get("energy_bins", {}).get(
            bin_key,
            {},
        )

        for icol, region in enumerate(("FD", "FT")):
            ax = axes[irow, icol]
            region_summary = bin_summary.get("regions", {}).get(region)
            closure_results = (
                region_summary.get("closure_results", [])
                if region_summary is not None
                else []
            )

            if not closure_results:
                ax.text(
                    0.5,
                    0.5,
                    "Insufficient closure statistics",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                )
                ax.set_title(
                    f"{region}: {e_min:g} <= E_probe < {e_max:g} GeV"
                )
                continue
            #endif

            for true_fraction in CLOSURE_TRUE_FRACTIONS:
                rows = [
                    result
                    for result in closure_results
                    if abs(
                        result["true_fraction"] - true_fraction
                    ) < 1.0e-9
                ]
                rows.sort(key=lambda result: result["step"])

                line, = ax.plot(
                    [result["step"] for result in rows],
                    [result["rms_error"] for result in rows],
                    marker="o",
                    label=rf"$f_{{\pi^0}}^{{true}}={true_fraction:.1f}$",
                )

                if irow == 0 and icol == 0:
                    legend_handles.append(line)
                    legend_labels.append(
                        rf"$f_{{\pi^0}}^{{true}}={true_fraction:.1f}$"
                    )
                #endif
            #endfor

            # Annotate which 1D shape carried the most separation at each step.
            step_rows = {}
            for result in closure_results:
                step_rows.setdefault(
                    int(result["step"]),
                    result,
                )
            #endfor

            annotation = "\n".join(
                (
                    f"{step}: {row['feature']} "
                    f"(JS={row['js_divergence']:.3f})"
                )
                for step, row in sorted(step_rows.items())
            )
            ax.text(
                0.98,
                0.97,
                annotation,
                transform=ax.transAxes,
                fontsize=7.4,
                va="top",
                ha="right",
                bbox={
                    "facecolor": "white",
                    "alpha": 0.82,
                    "edgecolor": "0.6",
                },
            )

            ax.set_xlabel("Cumulative cut step")
            ax.set_ylabel(r"closure RMS of $f_{\pi^0}$")
            ax.set_ylim(bottom=0.0)
            ax.grid(alpha=0.25)
            ax.set_title(
                f"{region}: {e_min:g} <= E_probe < {e_max:g} GeV"
            )
        #endfor
    #endfor

    if legend_handles:
        fig.legend(
            legend_handles,
            legend_labels,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.975),
            ncol=3,
            frameon=True,
        )
    #endif

    fig.suptitle(
        f"{period.label}: AAO/DVCS template-fit closure by cut step",
        fontsize=13,
        y=0.998,
    )
    fig.subplots_adjust(
        left=0.08,
        right=0.985,
        bottom=0.055,
        top=0.93,
        wspace=0.18,
        hspace=0.42,
    )

    out = output_dir / f"template_fit_closure_{period.key}.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out




# =============================================================================
# Real-data two-template composition fit
# =============================================================================

def closure_row_for_step(
    closure_results: list,
    step: int,
    true_fraction: float = 0.50,
) -> Optional[dict]:
    """Return the closure result nearest the requested injected fraction."""
    rows = [
        result
        for result in closure_results
        if int(result["step"]) == int(step)
    ]
    if not rows:
        return None
    #endif

    return min(
        rows,
        key=lambda result: abs(
            float(result["true_fraction"]) - float(true_fraction)
        ),
    )


def optimizer_retention_at_step(
    optimizer_results: list,
    step: int,
) -> Tuple[float, float, float]:
    """Cumulative AAO, DVCSgen, and data-sideband retention at one step."""
    if int(step) <= 0:
        return 1.0, 1.0, 1.0
    #endif

    for result in optimizer_results:
        if int(result["step"]) == int(step):
            return (
                float(result["total_pi0_eff"]),
                float(result["total_dvcs_eff"]),
                float(result["total_data_sideband_eff"]),
            )
        #endif
    #endfor

    return float("nan"), float("nan"), float("nan")


def choose_data_fit_operating_step(
    optimizer_results: list,
    closure_results: list,
    sideband_weight: float,
    max_abs_bias: float,
    max_rms: float,
    near_best_factor: float,
) -> dict:
    """
    Pick a cut step that rejects backgrounds without sacrificing closure.

    Closure safety is evaluated at injected f_pi0=0.5. Relative to the
    baseline closure RMS, the allowed RMS is

        min(max_rms, max(1.5 * RMS_0, RMS_0 + 0.005)).

    A step must also satisfy |closure bias| <= max_abs_bias.

    Among closure-safe steps, compute

        B = eps_DVCS + lambda * eps_data_sideband.

    The chosen operating point is the EARLIEST step whose B is within
    near_best_factor of the minimum closure-safe B. This avoids adding a later
    cut for a negligible purity gain.
    """
    baseline = closure_row_for_step(
        closure_results,
        0,
        true_fraction=0.50,
    )

    if baseline is None:
        return {
            "step": 0,
            "status": "fallback_no_baseline_closure",
            "closure_rms_limit": float(max_rms),
            "background_metric": 1.0 + float(sideband_weight),
        }
    #endif

    baseline_rms = float(baseline["rms_error"])
    rms_limit = min(
        float(max_rms),
        max(1.5 * baseline_rms, baseline_rms + 0.005),
    )

    max_step = len(optimizer_results)
    candidates = []

    for step in range(max_step + 1):
        closure = closure_row_for_step(
            closure_results,
            step,
            true_fraction=0.50,
        )
        if closure is None:
            continue
        #endif

        if abs(float(closure["bias"])) > float(max_abs_bias):
            continue
        #endif
        if float(closure["rms_error"]) > rms_limit:
            continue
        #endif

        _pi_eff, dvcs_eff, side_eff = optimizer_retention_at_step(
            optimizer_results,
            step,
        )
        if not np.isfinite(dvcs_eff) or not np.isfinite(side_eff):
            continue
        #endif

        background_metric = (
            float(dvcs_eff)
            + float(sideband_weight) * float(side_eff)
        )

        candidates.append(
            {
                "step": int(step),
                "closure": closure,
                "background_metric": float(background_metric),
            }
        )
    #endfor

    if not candidates:
        return {
            "step": 0,
            "status": "fallback_no_closure_safe_step",
            "closure_rms_limit": float(rms_limit),
            "background_metric": 1.0 + float(sideband_weight),
        }
    #endif

    best_background = min(
        candidate["background_metric"]
        for candidate in candidates
    )
    allowed_background = (
        float(near_best_factor) * best_background
        + 1.0e-15
    )

    near_best = [
        candidate
        for candidate in candidates
        if candidate["background_metric"] <= allowed_background
    ]
    chosen = min(near_best, key=lambda candidate: candidate["step"])

    return {
        "step": int(chosen["step"]),
        "status": "closure_safe_near_best_background",
        "closure_rms_limit": float(rms_limit),
        "background_metric": float(chosen["background_metric"]),
    }


def best_template_feature_for_events(
    pi0_events: np.ndarray,
    dvcs_events: np.ndarray,
    n_bins: int,
) -> Tuple[Optional[str], Optional[np.ndarray], float]:
    """
    Find the surviving 1D reconstructed variable with maximum JS divergence.

    This is primarily a fallback if closure information for the chosen step is
    unavailable; normally the real-data fit uses the discriminator selected by
    the closure study itself.
    """
    best_feature = None
    best_edges = None
    best_js = -1.0

    for feature in DATA_FIT_FEATURES:
        ifeature = FEATURE_INDEX[feature]
        pi_values = pi0_events[:, ifeature]
        dvcs_values = dvcs_events[:, ifeature]

        edges = closure_quantile_edges(
            pi_values,
            dvcs_values,
            n_bins,
        )
        if edges.size == 0:
            continue
        #endif

        p_pi = closure_histogram_probabilities(
            pi_values,
            edges,
        )
        p_dvcs = closure_histogram_probabilities(
            dvcs_values,
            edges,
        )
        js = jensen_shannon_divergence(p_pi, p_dvcs)

        if js > best_js:
            best_feature = feature
            best_edges = edges
            best_js = float(js)
        #endif
    #endfor

    return best_feature, best_edges, best_js



def morph_template_probability(
    probability: np.ndarray,
    shift_bins: float,
    smear_sigma_bins: float,
) -> np.ndarray:
    """Bounded AAO shift plus additional Gaussian smearing."""
    probability = np.asarray(probability, dtype=float)
    x = np.arange(probability.size, dtype=float)

    shifted = np.interp(
        x - float(shift_bins),
        x,
        probability,
        left=0.0,
        right=0.0,
    )

    sigma = float(smear_sigma_bins)
    if sigma > 1.0e-9:
        half_width = max(2, int(math.ceil(4.0 * sigma)))
        kernel_x = np.arange(-half_width, half_width + 1, dtype=float)
        kernel = np.exp(-0.5 * (kernel_x / sigma) ** 2)
        kernel /= np.sum(kernel)
        shifted = np.convolve(shifted, kernel, mode="same")
    #endif

    shifted = np.clip(shifted, 1.0e-15, None)
    shifted /= np.sum(shifted)
    return shifted


def simplex_fraction_grid(
    n_steps: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Triangular grid for f_pi0 + f_DVCS + f_SB = 1."""
    values = np.linspace(0.0, 1.0, max(8, int(n_steps)) + 1)
    f_pi0 = []
    f_dvcs = []
    f_side = []

    for pi_value in values:
        for side_value in values:
            if pi_value + side_value > 1.0 + 1.0e-12:
                continue
            #endif
            f_pi0.append(pi_value)
            f_side.append(side_value)
            f_dvcs.append(1.0 - pi_value - side_value)
        #endfor
    #endfor

    return (
        np.asarray(f_pi0, dtype=float),
        np.asarray(f_dvcs, dtype=float),
        np.asarray(f_side, dtype=float),
    )


def optimize_three_template_fractions(
    counts: np.ndarray,
    p_pi0: np.ndarray,
    p_dvcs: np.ndarray,
    p_sideband: np.ndarray,
    coarse_steps: int = 40,
) -> dict:
    """Fast, scipy-free, bounded three-template multinomial fit."""
    counts = np.asarray(counts, dtype=float)

    def evaluate(
        f_pi0: np.ndarray,
        f_dvcs: np.ndarray,
        f_side: np.ndarray,
    ) -> Tuple[int, np.ndarray]:
        mixtures = (
            f_pi0[:, None] * p_pi0[None, :]
            + f_dvcs[:, None] * p_dvcs[None, :]
            + f_side[:, None] * p_sideband[None, :]
        )
        mixtures = np.clip(mixtures, 1.0e-15, None)
        ll = np.sum(counts[None, :] * np.log(mixtures), axis=1)
        return int(np.argmax(ll)), ll

    f_pi0, f_dvcs, f_side = simplex_fraction_grid(coarse_steps)
    best_idx, _ = evaluate(f_pi0, f_dvcs, f_side)

    best_pi0 = float(f_pi0[best_idx])
    best_side = float(f_side[best_idx])
    coarse_spacing = 1.0 / max(8, int(coarse_steps))
    fine_spacing = coarse_spacing / 8.0

    pi_values = np.arange(
        max(0.0, best_pi0 - 2.0 * coarse_spacing),
        min(1.0, best_pi0 + 2.0 * coarse_spacing) + 0.5 * fine_spacing,
        fine_spacing,
    )
    side_values = np.arange(
        max(0.0, best_side - 2.0 * coarse_spacing),
        min(1.0, best_side + 2.0 * coarse_spacing) + 0.5 * fine_spacing,
        fine_spacing,
    )

    fine_pi0 = []
    fine_dvcs = []
    fine_side = []

    for pi_value in pi_values:
        for side_value in side_values:
            if pi_value + side_value > 1.0 + 1.0e-12:
                continue
            #endif
            fine_pi0.append(pi_value)
            fine_side.append(side_value)
            fine_dvcs.append(1.0 - pi_value - side_value)
        #endfor
    #endfor

    fine_pi0 = np.asarray(fine_pi0, dtype=float)
    fine_dvcs = np.asarray(fine_dvcs, dtype=float)
    fine_side = np.asarray(fine_side, dtype=float)

    best_idx, ll = evaluate(fine_pi0, fine_dvcs, fine_side)

    return {
        "fraction_pi0": float(fine_pi0[best_idx]),
        "fraction_dvcs": float(fine_dvcs[best_idx]),
        "fraction_sideband": float(fine_side[best_idx]),
        "log_likelihood": float(ll[best_idx]),
    }


def fit_three_templates_with_aao_morph(
    counts: np.ndarray,
    p_pi0_nominal: np.ndarray,
    p_dvcs: np.ndarray,
    p_sideband: np.ndarray,
    shift_max_bins: float,
    smear_max_bins: float,
    shift_steps: int,
    smear_steps: int,
) -> dict:
    """Profile bounded AAO-only morphing together with three fractions."""
    nominal = optimize_three_template_fractions(
        counts,
        p_pi0_nominal,
        p_dvcs,
        p_sideband,
    )

    best = None
    for shift_bins in np.linspace(
        -float(shift_max_bins),
        float(shift_max_bins),
        max(3, int(shift_steps)),
    ):
        for smear_bins in np.linspace(
            0.0,
            float(smear_max_bins),
            max(2, int(smear_steps)),
        ):
            p_pi0 = morph_template_probability(
                p_pi0_nominal,
                shift_bins,
                smear_bins,
            )
            fit = optimize_three_template_fractions(
                counts,
                p_pi0,
                p_dvcs,
                p_sideband,
            )
            if best is None or fit["log_likelihood"] > best["log_likelihood"]:
                best = dict(fit)
                best["shift_bins"] = float(shift_bins)
                best["smear_sigma_bins"] = float(smear_bins)
                best["p_pi0_morphed"] = p_pi0
            #endif
        #endfor
    #endfor

    best["nominal_fit"] = nominal
    return best


def approximate_fraction_uncertainty_three_template(
    counts: np.ndarray,
    p_pi0: np.ndarray,
    p_dvcs: np.ndarray,
    p_sideband: np.ndarray,
    fractions: Tuple[float, float, float],
) -> float:
    """
    Local curvature estimate of sigma(f_pi0) at fixed AAO morph, with the
    sideband fraction profiled. Morph-parameter uncertainty is not included.
    """
    f_pi0, f_dvcs, f_side = fractions
    mixture = (
        f_pi0 * p_pi0
        + f_dvcs * p_dvcs
        + f_side * p_sideband
    )
    mixture = np.clip(mixture, 1.0e-15, None)

    d_pi0 = p_pi0 - p_dvcs
    d_side = p_sideband - p_dvcs

    info = np.array(
        [
            [
                np.sum(counts * d_pi0 * d_pi0 / (mixture * mixture)),
                np.sum(counts * d_pi0 * d_side / (mixture * mixture)),
            ],
            [
                np.sum(counts * d_pi0 * d_side / (mixture * mixture)),
                np.sum(counts * d_side * d_side / (mixture * mixture)),
            ],
        ],
        dtype=float,
    )

    try:
        covariance = np.linalg.inv(info)
        return math.sqrt(max(float(covariance[0, 0]), 0.0))
    except np.linalg.LinAlgError:
        return float("inf")
    #endtry


def poisson_deviance(
    observed: np.ndarray,
    expected: np.ndarray,
) -> float:
    """Poisson deviance for a binned goodness-of-fit diagnostic."""
    observed = np.asarray(observed, dtype=float)
    expected = np.clip(np.asarray(expected, dtype=float), 1.0e-15, None)

    positive = observed > 0.0
    terms = np.array(expected - observed, copy=True)
    terms[positive] += (
        observed[positive]
        * np.log(observed[positive] / expected[positive])
    )
    return float(2.0 * np.sum(terms))



def choose_real_data_fit_feature(
    pi0_events: np.ndarray,
    dvcs_events: np.ndarray,
    optimizer_results: list,
    chosen_step: int,
    closure_results: list,
    n_bins: int,
) -> Optional[str]:
    """
    Resolve the real-data discriminator once, before deciding whether the
    empirical sideband is stable enough to enter the fit.

    Prefer the closure-selected discriminator at the chosen operating step.
    If that variable is not allowed for the real-data fit, fall back to the
    strongest AAO/DVCS discriminator among DATA_FIT_FEATURES.
    """
    closure = closure_row_for_step(
        closure_results,
        chosen_step,
        true_fraction=0.50,
    )
    if closure is not None:
        candidate = str(closure["feature"])
        if candidate in DATA_FIT_FEATURES:
            return candidate
        #endif
    #endif

    pi_masks = build_cumulative_optimizer_masks(
        pi0_events,
        optimizer_results,
    )
    dvcs_masks = build_cumulative_optimizer_masks(
        dvcs_events,
        optimizer_results,
    )
    step = min(
        int(chosen_step),
        len(pi_masks) - 1,
        len(dvcs_masks) - 1,
    )

    pi_selected = pi0_events[pi_masks[step]]
    dvcs_selected = dvcs_events[dvcs_masks[step]]

    feature, _edges, _js = best_template_feature_for_events(
        pi_selected,
        dvcs_selected,
        n_bins,
    )
    return feature


def fit_real_data_two_templates_morphed(
    data_events: np.ndarray,
    pi0_events: np.ndarray,
    dvcs_events: np.ndarray,
    optimizer_results: list,
    chosen_step: int,
    feature: str,
    n_bins: int,
    morph_shift_max_bins: float,
    morph_smear_max_bins: float,
    morph_shift_steps: int,
    morph_smear_steps: int,
    fallback_reason: str,
) -> dict:
    """
    Adaptive fallback fit: signal-like data = morphed AAO + fixed DVCSgen.

    The AAO nuisance prescription is identical to the three-template case.
    """
    pi_masks = build_cumulative_optimizer_masks(
        pi0_events,
        optimizer_results,
    )
    dvcs_masks = build_cumulative_optimizer_masks(
        dvcs_events,
        optimizer_results,
    )
    data_masks = build_cumulative_optimizer_masks(
        data_events,
        optimizer_results,
    )

    chosen_step = min(
        int(chosen_step),
        len(pi_masks) - 1,
        len(dvcs_masks) - 1,
        len(data_masks) - 1,
    )

    pi_selected = pi0_events[pi_masks[chosen_step]]
    dvcs_selected = dvcs_events[dvcs_masks[chosen_step]]
    data_selected = data_events[data_masks[chosen_step]]

    if min(
        pi_selected.shape[0],
        dvcs_selected.shape[0],
        data_selected.shape[0],
    ) < 20:
        return {
            "success": False,
            "status": "insufficient_selected_statistics",
            "step": int(chosen_step),
            "fit_mode": "two_template_fallback",
            "fallback_reason": fallback_reason,
        }
    #endif

    edges = data_fit_edges_for_feature(feature, n_bins=n_bins)
    if edges.size < 5:
        return {
            "success": False,
            "status": "no_valid_fit_range",
            "step": int(chosen_step),
            "fit_mode": "two_template_fallback",
            "fallback_reason": fallback_reason,
        }
    #endif

    ifeature = FEATURE_INDEX[feature]
    fit_min = float(edges[0])
    fit_max = float(edges[-1])

    def extract(events: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        values = np.asarray(events[:, ifeature], dtype=float)
        mask = (
            np.isfinite(values)
            & (values >= fit_min)
            & (values < fit_max)
        )
        return values, mask

    pi_values, pi_mask = extract(pi_selected)
    dvcs_values, dvcs_mask = extract(dvcs_selected)
    data_values, data_mask = extract(data_selected)

    if min(
        np.count_nonzero(pi_mask),
        np.count_nonzero(dvcs_mask),
    ) < 20:
        return {
            "success": False,
            "status": "insufficient_template_statistics_in_fit_range",
            "step": int(chosen_step),
            "fit_mode": "two_template_fallback",
            "fallback_reason": fallback_reason,
        }
    #endif

    p_pi0_nominal = closure_histogram_probabilities(
        pi_values[pi_mask],
        edges,
    )
    p_dvcs = closure_histogram_probabilities(
        dvcs_values[dvcs_mask],
        edges,
    )
    data_counts, _ = np.histogram(
        data_values[data_mask],
        bins=edges,
    )

    if np.sum(data_counts) < 20:
        return {
            "success": False,
            "status": "insufficient_data_in_fit_histogram",
            "step": int(chosen_step),
            "fit_mode": "two_template_fallback",
            "fallback_reason": fallback_reason,
        }
    #endif

    # Unmorphed reference.
    f_nominal, _sigma_nominal = fit_two_template_fraction(
        data_counts,
        p_pi0_nominal,
        p_dvcs,
    )
    nominal_probability = (
        f_nominal * p_pi0_nominal
        + (1.0 - f_nominal) * p_dvcs
    )

    best = None
    for shift_bins in np.linspace(
        -float(morph_shift_max_bins),
        float(morph_shift_max_bins),
        max(3, int(morph_shift_steps)),
    ):
        for smear_bins in np.linspace(
            0.0,
            float(morph_smear_max_bins),
            max(2, int(morph_smear_steps)),
        ):
            p_pi0 = morph_template_probability(
                p_pi0_nominal,
                shift_bins,
                smear_bins,
            )
            f_pi0, sigma_pi0 = fit_two_template_fraction(
                data_counts,
                p_pi0,
                p_dvcs,
            )
            probability = (
                f_pi0 * p_pi0
                + (1.0 - f_pi0) * p_dvcs
            )
            probability = np.clip(probability, 1.0e-15, None)
            log_likelihood = float(
                np.sum(data_counts * np.log(probability))
            )

            if (
                best is None
                or log_likelihood > best["log_likelihood"]
            ):
                best = {
                    "fraction_pi0": float(f_pi0),
                    "fraction_pi0_stat": float(sigma_pi0),
                    "shift_bins": float(shift_bins),
                    "smear_sigma_bins": float(smear_bins),
                    "p_pi0": p_pi0,
                    "probability": probability,
                    "log_likelihood": log_likelihood,
                }
            #endif
        #endfor
    #endfor

    n_data = int(np.sum(data_counts))
    f_pi0 = best["fraction_pi0"]
    f_dvcs = 1.0 - f_pi0
    expected_total = n_data * best["probability"]
    expected_pi0 = n_data * f_pi0 * best["p_pi0"]
    expected_dvcs = n_data * f_dvcs * p_dvcs
    expected_sideband = np.zeros_like(expected_total)

    expected_nominal = n_data * nominal_probability
    deviance = poisson_deviance(data_counts, expected_total)
    nominal_deviance = poisson_deviance(
        data_counts,
        expected_nominal,
    )

    # f_pi0 + AAO shift + AAO smear.
    ndf = max(int(data_counts.size) - 3, 1)
    bin_width = float(edges[1] - edges[0])

    return {
        "success": True,
        "status": "ok",
        "fit_mode": "two_template_fallback",
        "fallback_reason": str(fallback_reason),
        "step": int(chosen_step),
        "feature": str(feature),
        "fraction_pi0": float(f_pi0),
        "fraction_dvcs": float(f_dvcs),
        "fraction_sideband": 0.0,
        "fraction_pi0_stat": float(best["fraction_pi0_stat"]),
        "morph_shift_bins": float(best["shift_bins"]),
        "morph_smear_sigma_bins": float(best["smear_sigma_bins"]),
        "morph_shift_x": float(best["shift_bins"] * bin_width),
        "morph_smear_sigma_x": float(
            best["smear_sigma_bins"] * bin_width
        ),
        "deviance": float(deviance),
        "nominal_deviance": float(nominal_deviance),
        "ndf": int(ndf),
        "deviance_per_ndf": float(deviance / ndf),
        "nominal_deviance_per_ndf": float(
            nominal_deviance / ndf
        ),
        "n_data": n_data,
        "n_data_selected_total": int(data_selected.shape[0]),
        "n_pi0_template": int(np.count_nonzero(pi_mask)),
        "n_dvcs_template": int(np.count_nonzero(dvcs_mask)),
        "n_sideband_template": 0,
        "fit_range_min": fit_min,
        "fit_range_max": fit_max,
        "data_in_range_fraction": float(
            n_data / max(int(data_selected.shape[0]), 1)
        ),
        "edges": np.asarray(edges, dtype=float),
        "data_counts": np.asarray(data_counts, dtype=float),
        "expected_total": np.asarray(expected_total, dtype=float),
        "expected_nominal": np.asarray(expected_nominal, dtype=float),
        "expected_pi0": np.asarray(expected_pi0, dtype=float),
        "expected_dvcs": np.asarray(expected_dvcs, dtype=float),
        "expected_sideband": np.asarray(
            expected_sideband,
            dtype=float,
        ),
    }


def fit_real_data_three_templates(
    data_events: np.ndarray,
    pi0_events: np.ndarray,
    dvcs_events: np.ndarray,
    sideband_events: np.ndarray,
    optimizer_results: list,
    chosen_step: int,
    closure_results: list,
    n_bins: int,
    forced_feature: Optional[str],
    morph_shift_max_bins: float,
    morph_smear_max_bins: float,
    morph_shift_steps: int,
    morph_smear_steps: int,
) -> dict:
    """Fit signal-like data as morphed-AAO + DVCSgen + empirical sideband."""
    pi_masks = build_cumulative_optimizer_masks(pi0_events, optimizer_results)
    dvcs_masks = build_cumulative_optimizer_masks(dvcs_events, optimizer_results)
    data_masks = build_cumulative_optimizer_masks(data_events, optimizer_results)
    side_masks = build_cumulative_optimizer_masks(sideband_events, optimizer_results)

    chosen_step = min(
        int(chosen_step),
        len(pi_masks) - 1,
        len(dvcs_masks) - 1,
        len(data_masks) - 1,
        len(side_masks) - 1,
    )

    pi_selected = pi0_events[pi_masks[chosen_step]]
    dvcs_selected = dvcs_events[dvcs_masks[chosen_step]]
    data_selected = data_events[data_masks[chosen_step]]
    side_selected = sideband_events[side_masks[chosen_step]]

    if min(
        pi_selected.shape[0],
        dvcs_selected.shape[0],
        data_selected.shape[0],
        side_selected.shape[0],
    ) < 20:
        return {
            "success": False,
            "status": "insufficient_selected_statistics",
            "step": int(chosen_step),
        }
    #endif

    feature = (
        str(forced_feature)
        if forced_feature in DATA_FIT_FEATURES
        else None
    )

    if feature is None:
        closure = closure_row_for_step(
            closure_results,
            chosen_step,
            true_fraction=0.50,
        )
        if closure is not None:
            candidate = str(closure["feature"])
            if candidate in DATA_FIT_FEATURES:
                feature = candidate
            #endif
        #endif
    #endif

    if feature is None:
        feature, _unused_edges, _unused_js = best_template_feature_for_events(
            pi_selected,
            dvcs_selected,
            n_bins,
        )
    #endif

    if feature is None:
        return {
            "success": False,
            "status": "no_valid_discriminator",
            "step": int(chosen_step),
        }
    #endif

    edges = data_fit_edges_for_feature(feature, n_bins=n_bins)
    if edges.size < 5:
        return {
            "success": False,
            "status": "no_valid_fit_range",
            "step": int(chosen_step),
        }
    #endif

    ifeature = FEATURE_INDEX[feature]
    fit_min = float(edges[0])
    fit_max = float(edges[-1])

    def selected_values(events: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        values = np.asarray(events[:, ifeature], dtype=float)
        mask = (
            np.isfinite(values)
            & (values >= fit_min)
            & (values < fit_max)
        )
        return values, mask

    pi_values, pi_mask = selected_values(pi_selected)
    dvcs_values, dvcs_mask = selected_values(dvcs_selected)
    data_values, data_mask = selected_values(data_selected)
    side_values, side_mask = selected_values(side_selected)

    if min(
        np.count_nonzero(pi_mask),
        np.count_nonzero(dvcs_mask),
        np.count_nonzero(side_mask),
    ) < 20:
        return {
            "success": False,
            "status": "insufficient_template_statistics_in_fit_range",
            "step": int(chosen_step),
        }
    #endif

    p_pi0_nominal = closure_histogram_probabilities(pi_values[pi_mask], edges)
    p_dvcs = closure_histogram_probabilities(dvcs_values[dvcs_mask], edges)
    p_sideband = closure_histogram_probabilities(side_values[side_mask], edges)

    data_counts, _ = np.histogram(data_values[data_mask], bins=edges)
    if np.sum(data_counts) < 20:
        return {
            "success": False,
            "status": "insufficient_data_in_fit_histogram",
            "step": int(chosen_step),
        }
    #endif

    profiled = fit_three_templates_with_aao_morph(
        data_counts,
        p_pi0_nominal,
        p_dvcs,
        p_sideband,
        shift_max_bins=morph_shift_max_bins,
        smear_max_bins=morph_smear_max_bins,
        shift_steps=morph_shift_steps,
        smear_steps=morph_smear_steps,
    )

    f_pi0 = float(profiled["fraction_pi0"])
    f_dvcs = float(profiled["fraction_dvcs"])
    f_side = float(profiled["fraction_sideband"])
    p_pi0 = np.asarray(profiled["p_pi0_morphed"], dtype=float)

    n_data = int(np.sum(data_counts))
    total_probability = (
        f_pi0 * p_pi0
        + f_dvcs * p_dvcs
        + f_side * p_sideband
    )
    expected_total = n_data * total_probability
    expected_pi0 = n_data * f_pi0 * p_pi0
    expected_dvcs = n_data * f_dvcs * p_dvcs
    expected_sideband = n_data * f_side * p_sideband

    nominal_fit = profiled["nominal_fit"]
    nominal_probability = (
        nominal_fit["fraction_pi0"] * p_pi0_nominal
        + nominal_fit["fraction_dvcs"] * p_dvcs
        + nominal_fit["fraction_sideband"] * p_sideband
    )
    expected_nominal = n_data * nominal_probability

    deviance = poisson_deviance(data_counts, expected_total)
    nominal_deviance = poisson_deviance(data_counts, expected_nominal)

    # 2 independent fractions + AAO shift + AAO smear.
    ndf = max(int(data_counts.size) - 4, 1)

    sigma_pi0 = approximate_fraction_uncertainty_three_template(
        data_counts,
        p_pi0,
        p_dvcs,
        p_sideband,
        (f_pi0, f_dvcs, f_side),
    )

    bin_width = float(edges[1] - edges[0])

    return {
        "success": True,
        "status": "ok",
        "fit_mode": "three_template",
        "fallback_reason": "",
        "step": int(chosen_step),
        "feature": str(feature),
        "fraction_pi0": f_pi0,
        "fraction_dvcs": f_dvcs,
        "fraction_sideband": f_side,
        "fraction_pi0_stat": float(sigma_pi0),
        "morph_shift_bins": float(profiled["shift_bins"]),
        "morph_smear_sigma_bins": float(profiled["smear_sigma_bins"]),
        "morph_shift_x": float(profiled["shift_bins"] * bin_width),
        "morph_smear_sigma_x": float(
            profiled["smear_sigma_bins"] * bin_width
        ),
        "deviance": float(deviance),
        "nominal_deviance": float(nominal_deviance),
        "ndf": int(ndf),
        "deviance_per_ndf": float(deviance / ndf),
        "nominal_deviance_per_ndf": float(nominal_deviance / ndf),
        "n_data": n_data,
        "n_data_selected_total": int(data_selected.shape[0]),
        "n_pi0_template": int(np.count_nonzero(pi_mask)),
        "n_dvcs_template": int(np.count_nonzero(dvcs_mask)),
        "n_sideband_template": int(np.count_nonzero(side_mask)),
        "fit_range_min": fit_min,
        "fit_range_max": fit_max,
        "data_in_range_fraction": float(
            n_data / max(int(data_selected.shape[0]), 1)
        ),
        "edges": np.asarray(edges, dtype=float),
        "data_counts": np.asarray(data_counts, dtype=float),
        "expected_total": np.asarray(expected_total, dtype=float),
        "expected_nominal": np.asarray(expected_nominal, dtype=float),
        "expected_pi0": np.asarray(expected_pi0, dtype=float),
        "expected_dvcs": np.asarray(expected_dvcs, dtype=float),
        "expected_sideband": np.asarray(expected_sideband, dtype=float),
    }



def data_fit_edges_for_feature(
    feature: str,
    n_bins: int = 30,
) -> np.ndarray:
    """
    Fixed-width, physically bounded histogram edges for the REAL-DATA fit.

    Unlike the closure study, which deliberately uses quantile bins for stable
    pseudoexperiment statistics, the real-data composition fit uses the
    established plotting range for the selected discriminator.  Events outside
    this range are explicitly excluded rather than being absorbed into a huge
    terminal quantile bin.
    """
    for key, _title, _xlabel, x_min, x_max, _plot_bins in PLOT_SPECS:
        if key == feature:
            return np.linspace(
                float(x_min),
                float(x_max),
                int(n_bins) + 1,
            )
        #endif
    #endfor

    return np.empty(0, dtype=float)


def feature_plot_label(feature: str) -> str:
    """Readable x-axis label for a fit discriminator."""
    labels = {
        "Mx2": r"$MM^2(ep\gamma X)$ (GeV$^2$)",
        "Mx2_1": r"$MM^2(epX)$ (GeV$^2$)",
        "Mx2_2": r"$MM^2(e\gamma X)$ (GeV$^2$)",
        "Emiss2": r"$E_{\rm probe}$ (GeV)",
        "E_tag": r"$E_{\rm tag}$ (GeV)",
        "Delta_phi_residual_deg": (
            r"$\Delta\phi(p,\gamma)-180^\circ$ (deg)"
        ),
        "pTmiss": r"$p_{T,\rm miss}$ (GeV)",
        "theta_gamma_gamma": r"$\theta_{\gamma\gamma}$ (deg)",
    }
    return labels.get(feature, feature)



# =============================================================================
# TAG-FD / predicted-probe-sector diagnostics
# =============================================================================

def fd_sector_mask(events: np.ndarray, sector: int) -> np.ndarray:
    """Mask one PREDICTED probe-photon FD sector in a reservoir array."""
    values = np.asarray(
        events[:, EVENT_INDEX["pred_probe_sector"]],
        dtype=float,
    )
    return np.isfinite(values) & (np.rint(values).astype(int) == int(sector))


def build_fixed_composition_sector_result(
    data_all: np.ndarray,
    pi0_all: np.ndarray,
    dvcs_all: np.ndarray,
    optimizer_results: list,
    global_fit: dict,
    sector: int,
    n_bins: int,
) -> dict:
    """
    Sector detector diagnostic using the GLOBAL FD composition and morph.

    No fraction or nuisance parameter is refitted by sector. Only the
    sector-specific reconstructed template shapes change.
    """
    if not global_fit or not global_fit.get("success", False):
        return {
            "success": False,
            "status": "global_fd_fit_unavailable",
            "sector": int(sector),
        }
    #endif

    feature = str(global_fit["feature"])
    step = int(global_fit["step"])
    mx2_idx = FEATURE_INDEX["Mx2_1"]

    def sector_signal(events: np.ndarray) -> np.ndarray:
        return events[
            fd_sector_mask(events, sector)
            & np.isfinite(events[:, mx2_idx])
            & (events[:, mx2_idx] < SIGNAL_MX2_1_MAX)
        ]

    signal_data = sector_signal(data_all)
    signal_pi0 = sector_signal(pi0_all)
    signal_dvcs = sector_signal(dvcs_all)
    near_sideband = data_all[
        fd_sector_mask(data_all, sector)
        & np.isfinite(data_all[:, mx2_idx])
        & (data_all[:, mx2_idx] >= NEAR_SIDEBAND_MIN)
        & (data_all[:, mx2_idx] < NEAR_SIDEBAND_MAX)
    ]

    def apply_step(events: np.ndarray) -> np.ndarray:
        masks = build_cumulative_optimizer_masks(events, optimizer_results)
        return events[masks[min(step, len(masks) - 1)]]

    signal_data = apply_step(signal_data)
    signal_pi0 = apply_step(signal_pi0)
    signal_dvcs = apply_step(signal_dvcs)
    near_sideband = apply_step(near_sideband)

    edges = data_fit_edges_for_feature(feature, n_bins=n_bins)
    if edges.size < 5:
        return {
            "success": False,
            "status": "invalid_fit_edges",
            "sector": int(sector),
        }
    #endif

    ifeature = FEATURE_INDEX[feature]
    fit_min = float(edges[0])
    fit_max = float(edges[-1])

    def template_probability(events: np.ndarray) -> Tuple[np.ndarray, int]:
        values = np.asarray(events[:, ifeature], dtype=float)
        in_range = (
            np.isfinite(values)
            & (values >= fit_min)
            & (values < fit_max)
        )
        counts, _ = np.histogram(values[in_range], bins=edges)
        probs = counts.astype(float) + 0.5
        probs /= np.sum(probs)
        return probs, int(np.sum(counts))

    p_pi0_nominal, n_pi0 = template_probability(signal_pi0)
    p_dvcs, n_dvcs = template_probability(signal_dvcs)

    data_values = np.asarray(signal_data[:, ifeature], dtype=float)
    data_in_range = (
        np.isfinite(data_values)
        & (data_values >= fit_min)
        & (data_values < fit_max)
    )
    data_counts, _ = np.histogram(data_values[data_in_range], bins=edges)
    n_data = int(np.sum(data_counts))

    use_sideband = str(global_fit.get("fit_mode", "")) == "three_template"
    if use_sideband:
        p_sideband, n_sideband = template_probability(near_sideband)
    else:
        p_sideband = np.zeros_like(p_pi0_nominal)
        n_sideband = 0
    #endif

    required = [n_data, n_pi0, n_dvcs]
    if use_sideband:
        required.append(n_sideband)
    #endif
    if min(required) < 20:
        return {
            "success": False,
            "status": "insufficient_sector_statistics",
            "sector": int(sector),
            "n_data": n_data,
            "n_pi0_template": n_pi0,
            "n_dvcs_template": n_dvcs,
            "n_sideband_template": n_sideband,
        }
    #endif

    p_pi0 = morph_template_probability(
        p_pi0_nominal,
        float(global_fit.get("morph_shift_bins", 0.0)),
        float(global_fit.get("morph_smear_sigma_bins", 0.0)),
    )

    f_pi0 = float(global_fit["fraction_pi0"])
    f_dvcs = float(global_fit["fraction_dvcs"])
    f_side = float(global_fit["fraction_sideband"]) if use_sideband else 0.0

    total_probability = (
        f_pi0 * p_pi0
        + f_dvcs * p_dvcs
        + f_side * p_sideband
    )
    total_probability = np.clip(total_probability, 1.0e-15, None)
    total_probability /= np.sum(total_probability)

    expected_total = n_data * total_probability
    expected_pi0 = n_data * f_pi0 * p_pi0
    expected_dvcs = n_data * f_dvcs * p_dvcs
    expected_sideband = n_data * f_side * p_sideband

    pull = (
        data_counts - expected_total
    ) / np.sqrt(np.maximum(expected_total, 1.0))
    deviance = poisson_deviance(data_counts, expected_total)

    return {
        "success": True,
        "status": "ok",
        "sector": int(sector),
        "fit_mode": str(global_fit.get("fit_mode", "")),
        "feature": feature,
        "step": step,
        "fraction_pi0": f_pi0,
        "fraction_dvcs": f_dvcs,
        "fraction_sideband": f_side,
        "n_data": n_data,
        "n_pi0_template": n_pi0,
        "n_dvcs_template": n_dvcs,
        "n_sideband_template": n_sideband,
        "deviance": float(deviance),
        "deviance_per_bin": float(deviance / max(data_counts.size, 1)),
        "rms_pull": float(np.sqrt(np.mean(pull * pull))),
        "max_abs_pull": float(np.max(np.abs(pull))),
        "edges": np.asarray(edges, dtype=float),
        "data_counts": np.asarray(data_counts, dtype=float),
        "expected_total": np.asarray(expected_total, dtype=float),
        "expected_pi0": np.asarray(expected_pi0, dtype=float),
        "expected_dvcs": np.asarray(expected_dvcs, dtype=float),
        "expected_sideband": np.asarray(expected_sideband, dtype=float),
        "pull": np.asarray(pull, dtype=float),
    }


def build_fd_sector_diagnostics(
    optimizer_events: dict,
    optimizer_summary: dict,
    n_bins: int,
) -> dict:
    """Build sectors 1-6 using the global FD fit in each E_probe bin."""
    result = {"energy_bins": {}}

    for bin_key, e_min, e_max in E_PROBE_BINS:
        fd_summary = (
            optimizer_summary.get("energy_bins", {})
            .get(bin_key, {})
            .get("regions", {})
            .get("FD")
        )
        global_fit = fd_summary.get("data_fit") if fd_summary else None
        optimizer_results = fd_summary.get("results", []) if fd_summary else []

        sector_results = {}
        for sector in range(1, 7):
            sector_results[sector] = build_fixed_composition_sector_result(
                optimizer_events["data"]["FD"][bin_key],
                optimizer_events["pi0"]["FD"][bin_key],
                optimizer_events["dvcs"]["FD"][bin_key],
                optimizer_results,
                global_fit,
                sector,
                n_bins,
            )
        #endfor

        result["energy_bins"][bin_key] = {
            "energy_min": e_min,
            "energy_max": e_max,
            "global_fit": global_fit,
            "sectors": sector_results,
        }
    #endfor

    return result


def make_fd_sector_template_fit_canvases(
    period: Period,
    sector_diagnostics: dict,
    output_dir: Path,
) -> list:
    """One 3x2 sector canvas per E_probe bin; no per-sector files."""
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = []

    for bin_key, e_min, e_max in E_PROBE_BINS:
        bin_result = sector_diagnostics["energy_bins"][bin_key]
        global_fit = bin_result.get("global_fit")

        fig = plt.figure(figsize=(15.0, 14.2))
        outer = fig.add_gridspec(
            3, 2,
            hspace=0.36,
            wspace=0.20,
            left=0.07,
            right=0.985,
            bottom=0.055,
            top=0.91,
        )
        legend_handles = []
        legend_labels = []

        for sector in range(1, 7):
            row = (sector - 1) // 2
            col = (sector - 1) % 2
            inner = outer[row, col].subgridspec(
                2, 1,
                height_ratios=(3.2, 1.0),
                hspace=0.06,
            )
            ax = fig.add_subplot(inner[0, 0])
            ax_pull = fig.add_subplot(inner[1, 0], sharex=ax)
            result = bin_result["sectors"][sector]

            if not result.get("success", False):
                ax.text(
                    0.5, 0.5,
                    (
                        f"predicted probe sector {sector}\n"
                        f"{result.get('status', 'unavailable')}\n"
                        f"Ndata={result.get('n_data', 0):,}"
                    ),
                    ha="center", va="center", transform=ax.transAxes,
                )
                ax.set_title(f"FD sector {sector}")
                ax_pull.axis("off")
                continue
            #endif

            edges = np.asarray(result["edges"], dtype=float)
            centers = 0.5 * (edges[:-1] + edges[1:])
            widths = edges[1:] - edges[:-1]
            counts = np.asarray(result["data_counts"], dtype=float)

            data_handle = ax.errorbar(
                centers, counts,
                yerr=np.sqrt(np.maximum(counts, 1.0)),
                xerr=0.5 * widths,
                fmt="o", markersize=3.0, linewidth=0.9,
                label="data",
            )
            total_handle, = ax.step(
                centers, result["expected_total"],
                where="mid", linewidth=1.6,
                label="fixed global composition",
            )
            pi_handle, = ax.step(
                centers, result["expected_pi0"],
                where="mid", linewidth=1.2, linestyle="--",
                label=r"$\pi^0$",
            )
            dvcs_handle, = ax.step(
                centers, result["expected_dvcs"],
                where="mid", linewidth=1.2, linestyle=":",
                label="DVCS",
            )

            side_handle = None
            if result["fit_mode"] == "three_template":
                side_handle, = ax.step(
                    centers, result["expected_sideband"],
                    where="mid", linewidth=1.1, linestyle="-.",
                    label="sideband",
                )
            #endif

            if sector == 1:
                legend_handles = [
                    data_handle, total_handle, pi_handle, dvcs_handle
                ]
                legend_labels = [
                    "data", "fixed global composition", r"$\pi^0$", "DVCS"
                ]
                if side_handle is not None:
                    legend_handles.append(side_handle)
                    legend_labels.append("sideband")
                #endif
            #endif

            ax.set_title(
                f"predicted probe sector {sector} | {feature_plot_label(result['feature'])}",
                fontsize=9.5,
            )
            ax.set_ylabel("Events / bin")
            ax.grid(alpha=0.18)

            annotation = (
                f"global {result['fit_mode']}\n"
                rf"$f_{{\pi^0}}={result['fraction_pi0']:.3f}$, "
                rf"$f_{{DVCS}}={result['fraction_dvcs']:.3f}$"
                + (
                    rf", $f_{{SB}}={result['fraction_sideband']:.3f}$"
                    if result["fit_mode"] == "three_template" else ""
                )
                + "\n"
                + (
                    f"D/Nbin={result['deviance_per_bin']:.2f}; "
                    f"RMS pull={result['rms_pull']:.2f}; N={result['n_data']:,}"
                )
            )
            ax.text(
                0.98, 0.96, annotation,
                transform=ax.transAxes,
                ha="right", va="top", fontsize=7.4,
                bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "0.6"},
            )

            ax_pull.axhline(0.0, linestyle="--", linewidth=0.9)
            ax_pull.axhline(3.0, linestyle=":", linewidth=0.7)
            ax_pull.axhline(-3.0, linestyle=":", linewidth=0.7)
            ax_pull.plot(centers, result["pull"], "o", markersize=2.7)
            ax_pull.set_ylabel("pull", fontsize=8.5)
            ax_pull.set_xlabel(
                feature_plot_label(result["feature"]), fontsize=8.8
            )
            ax_pull.grid(alpha=0.16)
            plt.setp(ax.get_xticklabels(), visible=False)
        #endfor

        if legend_handles:
            fig.legend(
                legend_handles, legend_labels,
                loc="upper center", bbox_to_anchor=(0.5, 0.955),
                ncol=len(legend_handles), frameon=True,
            )
        #endif

        fit_description = "global FD fit unavailable"
        if global_fit and global_fit.get("success", False):
            fit_description = (
                f"global fractions fixed: "
                f"f_pi0={global_fit['fraction_pi0']:.3f}, "
                f"f_DVCS={global_fit['fraction_dvcs']:.3f}, "
                f"f_SB={global_fit['fraction_sideband']:.3f}"
            )
        #endif

        fig.suptitle(
            (
                f"{period.label}: TAG-FD / predicted-probe-sector diagnostics, "
                f"{e_min:g} <= E_probe < {e_max:g} GeV\n"
                + fit_description
            ),
            fontsize=12.5, y=0.995,
        )

        out = output_dir / f"sector_template_fit_{bin_key}.png"
        fig.savefig(out, dpi=180)
        plt.close(fig)
        outputs.append(out)
    #endfor

    return outputs


def make_fd_sector_overview_canvas(
    period: Period,
    sector_diagnostics: dict,
    output_dir: Path,
) -> Path:
    """Sector x E_probe overview of residual quality."""
    output_dir.mkdir(parents=True, exist_ok=True)
    shape = (6, len(E_PROBE_BINS))
    deviance = np.full(shape, np.nan)
    rms_pull = np.full(shape, np.nan)

    for ibin, (bin_key, _emin, _emax) in enumerate(E_PROBE_BINS):
        for sector in range(1, 7):
            result = sector_diagnostics["energy_bins"][bin_key]["sectors"][sector]
            if result.get("success", False):
                deviance[sector - 1, ibin] = result["deviance_per_bin"]
                rms_pull[sector - 1, ibin] = result["rms_pull"]
            #endif
        #endfor
    #endfor

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.0), constrained_layout=True)
    xlabels = [f"{emin:g}-{emax:g}" for _key, emin, emax in E_PROBE_BINS]
    ylabels = [f"S{s}" for s in range(1, 7)]

    for ax, matrix, title, cbar_label in (
        (axes[0], deviance, "Poisson deviance / histogram bin", "D / Nbin"),
        (axes[1], rms_pull, "RMS Pearson pull", "RMS pull"),
    ):
        image = ax.imshow(matrix, aspect="auto", origin="upper")
        ax.set_xticks(range(len(xlabels)))
        ax.set_xticklabels(xlabels)
        ax.set_yticks(range(6))
        ax.set_yticklabels(ylabels)
        ax.set_xlabel(r"$E_{\rm probe}$ bin (GeV)")
        ax.set_ylabel("predicted probe FD sector")
        ax.set_title(title)
        for isector in range(6):
            for ibin in range(len(E_PROBE_BINS)):
                value = matrix[isector, ibin]
                if np.isfinite(value):
                    ax.text(
                        ibin, isector, f"{value:.1f}",
                        ha="center", va="center", fontsize=8,
                    )
                #endif
            #endfor
        #endfor
        cbar = fig.colorbar(image, ax=ax)
        cbar.set_label(cbar_label)
    #endfor

    fig.suptitle(
        (
            f"{period.label}: FD sector residual overview\n"
            "global FD composition/cuts/morph fixed in every sector"
        ),
        fontsize=12.5,
    )
    out = output_dir / f"sector_overview_{period.key}.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def write_fd_sector_summary_csv(
    period: Period,
    sector_diagnostics: dict,
    output_dir: Path,
) -> Path:
    """One compact FD-sector summary CSV per period."""
    output_dir.mkdir(parents=True, exist_ok=True)
    out = output_dir / f"sector_summary_{period.key}.csv"

    fields = (
        "period", "energy_bin", "energy_min", "energy_max", "sector",
        "status", "global_fit_mode", "feature", "step",
        "global_fraction_pi0", "global_fraction_dvcs",
        "global_fraction_sideband", "n_data", "n_pi0_template",
        "n_dvcs_template", "n_sideband_template", "deviance",
        "deviance_per_bin", "rms_pull", "max_abs_pull",
    )

    rows = []
    for bin_key, e_min, e_max in E_PROBE_BINS:
        for sector in range(1, 7):
            result = sector_diagnostics["energy_bins"][bin_key]["sectors"][sector]
            rows.append({
                "period": period.key,
                "energy_bin": bin_key,
                "energy_min": e_min,
                "energy_max": e_max,
                "sector": sector,
                "status": result.get("status", ""),
                "global_fit_mode": result.get("fit_mode", ""),
                "feature": result.get("feature", ""),
                "step": result.get("step", ""),
                "global_fraction_pi0": result.get("fraction_pi0", ""),
                "global_fraction_dvcs": result.get("fraction_dvcs", ""),
                "global_fraction_sideband": result.get("fraction_sideband", ""),
                "n_data": result.get("n_data", ""),
                "n_pi0_template": result.get("n_pi0_template", ""),
                "n_dvcs_template": result.get("n_dvcs_template", ""),
                "n_sideband_template": result.get("n_sideband_template", ""),
                "deviance": result.get("deviance", ""),
                "deviance_per_bin": result.get("deviance_per_bin", ""),
                "rms_pull": result.get("rms_pull", ""),
                "max_abs_pull": result.get("max_abs_pull", ""),
            })
        #endfor
    #endfor

    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    #endwith
    return out


def make_combined_data_template_fit_canvas(
    period: Period,
    optimizer_summary: dict,
    output_dir: Path,
) -> Path:
    """
    One real-data template-fit canvas per period.

    Rows are E_probe bins and columns are FD/FT. Each logical panel is built
    with a NESTED GridSpec:
      * top: data + fitted AAO + fitted DVCSgen + total;
      * bottom: Pearson pull (data-fit)/sqrt(fit).

    The nested layout keeps the fit and pull visually grouped while giving a
    much larger gap between successive E_probe rows. This prevents the pull
    x-axis label from colliding with the title of the next fit panel.

    The real-data fit uses fixed-width bins over the established physical
    plotting range of the selected discriminator. Events outside that range
    are excluded explicitly, never folded into a giant edge bin.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    n_energy = len(E_PROBE_BINS)

    fig = plt.figure(
        figsize=(14.8, 6.4 * n_energy + 1.5),
    )

    # Outer grid controls spacing BETWEEN distinct E_probe regions.
    outer = fig.add_gridspec(
        n_energy,
        2,
        hspace=0.52,
        wspace=0.20,
        left=0.075,
        right=0.985,
        bottom=0.045,
        top=0.945,
    )

    legend_handles = []
    legend_labels = []

    for irow, (bin_key, e_min, e_max) in enumerate(E_PROBE_BINS):
        bin_summary = optimizer_summary.get("energy_bins", {}).get(
            bin_key,
            {},
        )

        for icol, region in enumerate(("FD", "FT")):
            # Inner grid keeps main panel and pull panel close together.
            inner = outer[irow, icol].subgridspec(
                2,
                1,
                height_ratios=(3.5, 1.0),
                hspace=0.06,
            )

            ax = fig.add_subplot(inner[0, 0])
            ax_pull = fig.add_subplot(
                inner[1, 0],
                sharex=ax,
            )

            region_summary = bin_summary.get("regions", {}).get(region)
            fit = (
                region_summary.get("data_fit")
                if region_summary is not None
                else None
            )

            if not fit or not fit.get("success", False):
                status = (
                    fit.get("status", "no fit")
                    if fit
                    else "no fit"
                )
                ax.text(
                    0.5,
                    0.5,
                    f"Real-data fit unavailable\n{status}",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                )
                ax.set_title(
                    (
                        f"{region}: {e_min:g} <= E_probe < {e_max:g} GeV\n"
                        "fit unavailable"
                    ),
                    fontsize=10.5,
                    pad=8,
                )
                ax_pull.axis("off")
                continue
            #endif

            edges = np.asarray(fit["edges"], dtype=float)
            centers = 0.5 * (edges[:-1] + edges[1:])
            widths = edges[1:] - edges[:-1]
            counts = np.asarray(fit["data_counts"], dtype=float)
            expected = np.asarray(fit["expected_total"], dtype=float)

            data_handle = ax.errorbar(
                centers,
                counts,
                yerr=np.sqrt(np.maximum(counts, 1.0)),
                xerr=0.5 * widths,
                fmt="o",
                markersize=3.4,
                linewidth=1.0,
                label="data",
            )
            total_handle, = ax.step(
                centers,
                expected,
                where="mid",
                linewidth=1.8,
                label="AAO + DVCS fit",
            )
            pi_handle, = ax.step(
                centers,
                fit["expected_pi0"],
                where="mid",
                linewidth=1.3,
                linestyle="--",
                label=r"fitted $\pi^0$ component",
            )
            dvcs_handle, = ax.step(
                centers,
                fit["expected_dvcs"],
                where="mid",
                linewidth=1.3,
                linestyle=":",
                label="fitted DVCS component",
            )
            side_handle, = ax.step(
                centers,
                fit["expected_sideband"],
                where="mid",
                linewidth=1.2,
                linestyle="-.",
                label="fitted data-sideband component",
            )

            if irow == 0 and icol == 0:
                legend_handles = [
                    data_handle,
                    total_handle,
                    pi_handle,
                    dvcs_handle,
                    side_handle,
                ]
                legend_labels = [
                    "data",
                    "AAO + DVCS + sideband fit",
                    r"fitted morphed $\pi^0$ component",
                    "fitted DVCS component",
                    "fitted data-sideband component",
                ]
            #endif

            # Put the discriminator on its own second title line. This is much
            # easier to read than appending it to the E_probe title.
            ax.set_title(
                (
                    f"{region}: {e_min:g} <= E_probe < {e_max:g} GeV\n"
                    f"fit variable: {feature_plot_label(fit['feature'])}"
                ),
                pad=8,
                fontsize=10.2,
            )
            ax.set_ylabel("Events / bin")
            ax.grid(alpha=0.20)

            pull = (
                counts - expected
            ) / np.sqrt(np.maximum(expected, 1.0))

            ax_pull.axhline(
                0.0,
                linewidth=1.0,
                linestyle="--",
            )
            ax_pull.axhline(
                3.0,
                linewidth=0.8,
                linestyle=":",
            )
            ax_pull.axhline(
                -3.0,
                linewidth=0.8,
                linestyle=":",
            )
            ax_pull.plot(
                centers,
                pull,
                "o",
                markersize=3.0,
            )
            ax_pull.set_ylabel(
                "pull",
                fontsize=9.0,
                labelpad=4,
            )

            # The variable is already explicit in the panel title. Keep the
            # x-axis label compact and only on the pull panel.
            ax_pull.set_xlabel(
                feature_plot_label(fit["feature"]),
                labelpad=5,
                fontsize=9.2,
            )
            ax_pull.grid(alpha=0.18)

            mode_label = (
                "3-template"
                if fit.get("fit_mode") == "three_template"
                else "2-template fallback"
            )
            annotation = (
                f"selected cut step {fit['step']}; {mode_label}\n"
                rf"$f_{{\pi^0}}={fit['fraction_pi0']:.3f}$, "
                rf"$f_{{DVCS}}={fit['fraction_dvcs']:.3f}$, "
                rf"$f_{{SB}}={fit['fraction_sideband']:.3f}$"
                + "\n"
                + (
                    f"AAO morph: shift={fit['morph_shift_bins']:+.2f} bins, "
                    f"smear={fit['morph_smear_sigma_bins']:.2f} bins"
                )
                + "\n"
                + (
                    f"D/ndf: nominal {fit['nominal_deviance_per_ndf']:.2f} "
                    f"-> morphed {fit['deviance_per_ndf']:.2f}"
                )
                + f", N={fit['n_data']:,}\n"
                + (
                    f"fit range: {fit['fit_range_min']:.3g} to "
                    f"{fit['fit_range_max']:.3g}; "
                    f"in range={100.0*fit['data_in_range_fraction']:.1f}%"
                )
                + (
                    (
                        "\nSB fallback: "
                        + str(fit.get("fallback_reason", ""))
                    )
                    if fit.get("fit_mode") == "two_template_fallback"
                    else ""
                )
            )
            ax.text(
                0.98,
                0.96,
                annotation,
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=7.6,
                bbox={
                    "facecolor": "white",
                    "alpha": 0.86,
                    "edgecolor": "0.6",
                },
            )

            plt.setp(ax.get_xticklabels(), visible=False)
        #endfor
    #endfor

    if legend_handles:
        fig.legend(
            legend_handles,
            legend_labels,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.982),
            ncol=5,
            frameon=True,
        )
    #endif

    fig.suptitle(
        f"{period.label}: real-data AAO + DVCS + sideband template fits",
        fontsize=13,
        y=0.999,
    )

    out = output_dir / f"data_template_fit_{period.key}.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out



def make_pi0_fraction_summary_canvas(
    period: Period,
    optimizer_summary: dict,
    output_dir: Path,
) -> Path:
    """Fitted AAO/DVCS/sideband composition versus E_probe for FD and FT."""
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(11.8, 4.9),
        sharey=True,
    )

    for ax, region in zip(axes, ("FD", "FT")):
        x = []
        xerr_low = []
        xerr_high = []
        f_pi0 = []
        f_dvcs = []
        f_side = []
        fit_modes = []

        for bin_key, e_min, e_max in E_PROBE_BINS:
            region_summary = (
                optimizer_summary.get("energy_bins", {})
                .get(bin_key, {})
                .get("regions", {})
                .get(region)
            )
            if region_summary is None:
                continue
            #endif

            fit = region_summary.get("data_fit")
            if not fit or not fit.get("success", False):
                continue
            #endif

            center = 0.5 * (e_min + e_max)
            x.append(center)
            xerr_low.append(center - e_min)
            xerr_high.append(e_max - center)
            f_pi0.append(float(fit["fraction_pi0"]))
            f_dvcs.append(float(fit["fraction_dvcs"]))
            f_side.append(float(fit["fraction_sideband"]))
            fit_modes.append(str(fit.get("fit_mode", "unknown")))
        #endfor

        if x:
            xerr = np.vstack((xerr_low, xerr_high))
            ax.errorbar(x, f_pi0, xerr=xerr, fmt="o-", label=r"$f_{\pi^0}$")
            ax.errorbar(x, f_dvcs, xerr=xerr, fmt="s--", label=r"$f_{\rm DVCS}$")
            ax.errorbar(x, f_side, xerr=xerr, fmt="^-.", label=r"$f_{\rm sideband}$")

            for xi, yi, mode in zip(x, f_pi0, fit_modes):
                ax.annotate(
                    "3T" if mode == "three_template" else "2T",
                    (xi, yi),
                    xytext=(4, 5),
                    textcoords="offset points",
                    fontsize=7.5,
                )
            #endfor
        #endif

        ax.set_xlim(0.35, 9.6)
        ax.set_ylim(0.0, 1.0)
        ax.set_xlabel(r"$E_{\rm probe}$ (GeV)")
        ax.set_title(region)
        ax.grid(alpha=0.25)
        ax.legend(frameon=True)
    #endfor

    axes[0].set_ylabel("Fitted composition fraction")
    fig.suptitle(
        f"{period.label}: preliminary three-template composition",
        fontsize=13,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))

    out = output_dir / f"composition_fraction_summary_{period.key}.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out



# =============================================================================
# Empirical sideband-transfer validation
# =============================================================================

def build_sideband_stability_result(
    data_all: np.ndarray,
    optimizer_results: list,
    chosen_step: int,
    feature: str,
    n_bins: int,
) -> dict:
    """
    Compare adjacent Mx2_1 slices after the exact cumulative cuts used for the
    real-data composition fit.

    The returned histograms are unit normalized and use the same fixed physical
    fit range as the selected real-data discriminator.
    """
    if feature not in DATA_FIT_FEATURES:
        return {
            "success": False,
            "status": "invalid_fit_feature",
        }
    #endif

    edges = data_fit_edges_for_feature(feature, n_bins=n_bins)
    if edges.size < 5:
        return {
            "success": False,
            "status": "invalid_fit_edges",
        }
    #endif

    mx2_idx = FEATURE_INDEX["Mx2_1"]
    ifeature = FEATURE_INDEX[feature]

    slice_events = {}
    slice_counts = {}

    for slice_key, low, high in SIDEBAND_VALIDATION_SLICES:
        mask = (
            np.isfinite(data_all[:, mx2_idx])
            & (data_all[:, mx2_idx] >= low)
            & (data_all[:, mx2_idx] < high)
        )
        events = data_all[mask]

        masks = build_cumulative_optimizer_masks(
            events,
            optimizer_results,
        )
        step = min(int(chosen_step), len(masks) - 1)
        selected = events[masks[step]]

        values = np.asarray(selected[:, ifeature], dtype=float)
        in_range = (
            np.isfinite(values)
            & (values >= edges[0])
            & (values < edges[-1])
        )
        values = values[in_range]

        counts, _ = np.histogram(values, bins=edges)
        probs = counts.astype(float) + 0.5
        probs /= np.sum(probs)

        slice_events[slice_key] = probs
        slice_counts[slice_key] = int(np.sum(counts))
    #endfor

    if min(slice_counts.values()) < 20:
        return {
            "success": False,
            "status": "insufficient_slice_statistics",
            "edges": edges,
            "slice_counts": slice_counts,
        }
    #endif

    p_signal_edge = slice_events["signal_edge"]
    p_near = slice_events["near"]
    p_far = slice_events["far"]

    return {
        "success": True,
        "status": "ok",
        "feature": feature,
        "step": int(chosen_step),
        "edges": np.asarray(edges, dtype=float),
        "probabilities": slice_events,
        "slice_counts": slice_counts,
        "js_signal_edge_near": float(
            jensen_shannon_divergence(p_signal_edge, p_near)
        ),
        "js_near_far": float(
            jensen_shannon_divergence(p_near, p_far)
        ),
        "js_signal_edge_far": float(
            jensen_shannon_divergence(p_signal_edge, p_far)
        ),
    }



def decide_sideband_template_usage(
    stability: dict,
    min_events: int,
    max_js: float,
) -> Tuple[bool, str]:
    """
    Decide whether the near eta-safe sideband is trustworthy enough to enter
    the composition fit as a third template.
    """
    if not stability.get("success", False):
        return (
            False,
            f"sideband_validation_{stability.get('status', 'failed')}",
        )
    #endif

    counts = stability.get("slice_counts", {})
    n_near = int(counts.get("near", 0))
    n_far = int(counts.get("far", 0))

    if n_near < int(min_events) or n_far < int(min_events):
        return (
            False,
            (
                f"sideband_stats_near={n_near}_far={n_far}"
                f"_below_{int(min_events)}"
            ),
        )
    #endif

    js = float(stability.get("js_near_far", float("inf")))
    if not np.isfinite(js) or js > float(max_js):
        return (
            False,
            f"sideband_js={js:.4f}_above_{float(max_js):.4f}",
        )
    #endif

    return (
        True,
        (
            f"sideband_accepted_near={n_near}_far={n_far}"
            f"_js={js:.4f}"
        ),
    )


def make_sideband_stability_canvas(
    period: Period,
    optimizer_summary: dict,
    output_dir: Path,
) -> Path:
    """
    One 4x2 sideband-transfer canvas per period.

    The three displayed real-data control slices are:
      signal edge : 0.11 <= Mx2_1 < 0.15
      near SB     : 0.15 <= Mx2_1 < 0.19
      far SB      : 0.19 <= Mx2_1 < 0.23

    The near sideband is the actual empirical third template used in the
    composition fit.  The farther sideband tests transfer stability while
    staying well below m_eta^2 ~= 0.300 GeV^2.  The signal-edge slice is only
    a diagnostic because it can already contain AAO/DVCS signal.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    n_energy = len(E_PROBE_BINS)

    fig = plt.figure(
        figsize=(14.8, 6.2 * n_energy + 1.4),
    )
    outer = fig.add_gridspec(
        n_energy,
        2,
        hspace=0.50,
        wspace=0.20,
        left=0.075,
        right=0.985,
        bottom=0.045,
        top=0.945,
    )

    legend_handles = []
    legend_labels = []

    labels = {
        "signal_edge": r"signal edge: $0.11\leq Mx2_1<0.15$",
        "near": r"near SB: $0.15\leq Mx2_1<0.19$",
        "far": r"far SB: $0.19\leq Mx2_1<0.23$",
    }
    linestyles = {
        "signal_edge": "-",
        "near": "--",
        "far": ":",
    }

    for irow, (bin_key, e_min, e_max) in enumerate(E_PROBE_BINS):
        bin_summary = optimizer_summary.get("energy_bins", {}).get(
            bin_key,
            {},
        )

        for icol, region in enumerate(("FD", "FT")):
            inner = outer[irow, icol].subgridspec(
                2,
                1,
                height_ratios=(3.3, 1.0),
                hspace=0.06,
            )
            ax = fig.add_subplot(inner[0, 0])
            ax_ratio = fig.add_subplot(
                inner[1, 0],
                sharex=ax,
            )

            region_summary = bin_summary.get("regions", {}).get(region)
            stability = (
                region_summary.get("sideband_stability")
                if region_summary is not None
                else None
            )

            if not stability or not stability.get("success", False):
                status = (
                    stability.get("status", "no result")
                    if stability
                    else "no result"
                )
                ax.text(
                    0.5,
                    0.5,
                    f"Sideband stability unavailable\n{status}",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                )
                ax.set_title(
                    f"{region}: {e_min:g} <= E_probe < {e_max:g} GeV",
                    fontsize=10.5,
                )
                ax_ratio.axis("off")
                continue
            #endif

            edges = np.asarray(stability["edges"], dtype=float)
            centers = 0.5 * (edges[:-1] + edges[1:])
            p_near = np.asarray(
                stability["probabilities"]["near"],
                dtype=float,
            )

            for slice_key in ("signal_edge", "near", "far"):
                probs = np.asarray(
                    stability["probabilities"][slice_key],
                    dtype=float,
                )

                line, = ax.step(
                    centers,
                    probs,
                    where="mid",
                    linewidth=1.5,
                    linestyle=linestyles[slice_key],
                    label=labels[slice_key],
                )

                if irow == 0 and icol == 0:
                    legend_handles.append(line)
                    legend_labels.append(labels[slice_key])
                #endif

                if slice_key != "near":
                    ratio = np.divide(
                        probs,
                        p_near,
                        out=np.full_like(probs, np.nan),
                        where=p_near > 1.0e-12,
                    )
                    ax_ratio.step(
                        centers,
                        ratio,
                        where="mid",
                        linewidth=1.2,
                        linestyle=linestyles[slice_key],
                    )
                #endif
            #endfor

            ax.set_ylabel("unit-normalized")
            ax.grid(alpha=0.20)
            ax.set_title(
                (
                    f"{region}: {e_min:g} <= E_probe < {e_max:g} GeV\n"
                    f"fit variable: {feature_plot_label(stability['feature'])}"
                ),
                fontsize=10.2,
                pad=8,
            )

            annotation = (
                rf"$JS_{{near,far}}={stability['js_near_far']:.3f}$"
                + "\n"
                + (
                    rf"$JS_{{signal\ edge,near}}="
                    rf"{stability['js_signal_edge_near']:.3f}$"
                )
                + "\n"
                + (
                    "N=("
                    f"{stability['slice_counts']['signal_edge']:,}, "
                    f"{stability['slice_counts']['near']:,}, "
                    f"{stability['slice_counts']['far']:,})"
                )
            )
            ax.text(
                0.98,
                0.96,
                annotation,
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=8.0,
                bbox={
                    "facecolor": "white",
                    "alpha": 0.86,
                    "edgecolor": "0.6",
                },
            )

            ax_ratio.axhline(
                1.0,
                linestyle="--",
                linewidth=1.0,
            )
            ax_ratio.set_ylabel("ratio\nto near SB", fontsize=8.5)
            ax_ratio.set_xlabel(
                feature_plot_label(stability["feature"]),
                fontsize=9.2,
                labelpad=5,
            )
            ax_ratio.set_ylim(0.0, 2.0)
            ax_ratio.grid(alpha=0.18)
            plt.setp(ax.get_xticklabels(), visible=False)
        #endfor
    #endfor

    if legend_handles:
        fig.legend(
            legend_handles,
            legend_labels,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.982),
            ncol=3,
            frameon=True,
        )
    #endif

    fig.suptitle(
        f"{period.label}: eta-safe Mx2_1 sideband-shape stability",
        fontsize=13,
        y=0.999,
    )

    out = output_dir / f"sideband_stability_{period.key}.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out



def make_mx2_1_control_region_canvas(
    period: Period,
    optimizer_summary: dict,
    optimizer_events: dict,
    output_dir: Path,
) -> Path:
    """
    Show where the signal/sideband definitions lie relative to the pi0 and eta
    missing-mass-squared locations after the selected operating cuts.

    This is diagnostic only; the eta region is never used as a fit template.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(
        len(E_PROBE_BINS),
        2,
        figsize=(13.8, 3.6 * len(E_PROBE_BINS) + 1.2),
        squeeze=False,
    )

    bins = np.linspace(-0.05, 0.40, 91)
    centers = 0.5 * (bins[:-1] + bins[1:])

    for irow, (bin_key, e_min, e_max) in enumerate(E_PROBE_BINS):
        for icol, region in enumerate(("FD", "FT")):
            ax = axes[irow, icol]

            region_summary = (
                optimizer_summary.get("energy_bins", {})
                .get(bin_key, {})
                .get("regions", {})
                .get(region)
            )
            if region_summary is None:
                continue
            #endif

            data_all = optimizer_events["data"][region][bin_key]
            optimizer_results = region_summary.get("results", [])
            data_fit = region_summary.get("data_fit", {})
            chosen_step = int(data_fit.get("step", 0))

            masks = build_cumulative_optimizer_masks(
                data_all,
                optimizer_results,
            )
            chosen_step = min(chosen_step, len(masks) - 1)
            selected = data_all[masks[chosen_step]]

            values = np.asarray(
                selected[:, FEATURE_INDEX["Mx2_1"]],
                dtype=float,
            )
            values = values[np.isfinite(values)]

            counts, _ = np.histogram(values, bins=bins)
            y = counts.astype(float)
            if np.sum(y) > 0.0:
                y /= np.sum(y)
            #endif

            ax.step(
                centers,
                y,
                where="mid",
                linewidth=1.4,
                label="data",
            )

            for x_value, linestyle in (
                (SIGNAL_MX2_1_MAX, "--"),
                (NEAR_SIDEBAND_MAX, "--"),
                (FAR_SIDEBAND_MAX, "--"),
            ):
                ax.axvline(
                    x_value,
                    linestyle=linestyle,
                    linewidth=1.0,
                )
            #endfor

            ax.axvline(
                PI0_MASS2_GEV2,
                linestyle=":",
                linewidth=1.2,
            )
            ax.axvline(
                ETA_MASS2_GEV2,
                linestyle=":",
                linewidth=1.2,
            )

            ax.text(
                PI0_MASS2_GEV2,
                0.97,
                r"$m_{\pi^0}^2$",
                transform=ax.get_xaxis_transform(),
                rotation=90,
                va="top",
                ha="right",
                fontsize=8,
            )
            ax.text(
                ETA_MASS2_GEV2,
                0.97,
                r"$m_{\eta}^2$",
                transform=ax.get_xaxis_transform(),
                rotation=90,
                va="top",
                ha="right",
                fontsize=8,
            )

            ax.axvspan(
                -0.05,
                SIGNAL_MX2_1_MAX,
                alpha=0.08,
            )
            ax.axvspan(
                NEAR_SIDEBAND_MIN,
                NEAR_SIDEBAND_MAX,
                alpha=0.08,
            )
            ax.axvspan(
                FAR_SIDEBAND_MIN,
                FAR_SIDEBAND_MAX,
                alpha=0.08,
            )

            ax.set_xlim(-0.05, 0.40)
            ax.set_xlabel(r"$Mx2_1=MM^2(epX)$ (GeV$^2$)")
            ax.set_ylabel("unit-normalized")
            ax.set_title(
                f"{region}: {e_min:g} <= E_probe < {e_max:g} GeV; "
                f"step {chosen_step}",
                fontsize=9.5,
            )
            ax.grid(alpha=0.18)
        #endfor
    #endfor

    fig.suptitle(
        (
            f"{period.label}: Mx2_1 control regions "
            "(eta shown but excluded from background template)"
        ),
        fontsize=12.5,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.965))

    out = output_dir / f"mx2_1_control_regions_{period.key}.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def write_period_summary_csv(
    period: Period,
    optimizer_summary: dict,
    output_dir: Path,
) -> Path:
    """
    Write ONE compact CSV per period containing both optimizer and closure
    information for all E_probe bins and both detector regions.

    This replaces the previous many-per-bin CSV/TXT output scheme.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    out = output_dir / f"concept_summary_{period.key}.csv"

    fields = (
        "record_type",
        "period",
        "region",
        "energy_bin",
        "energy_min",
        "energy_max",
        "step",
        "feature",
        "direction",
        "threshold",
        "step_pi0_eff",
        "step_dvcs_eff",
        "step_data_sideband_eff",
        "total_pi0_eff",
        "total_dvcs_eff",
        "total_data_sideband_eff",
        "score",
        "js_divergence",
        "true_fraction",
        "mean_fit",
        "bias",
        "rms_error",
        "std_fit",
        "mean_sigma",
        "boundary_fraction",
        "n_pi0_baseline",
        "n_dvcs_baseline",
        "n_data_sideband_baseline",
        "operating_step",
        "operating_step_status",
        "data_fit_feature",
        "data_fit_mode",
        "data_fit_fallback_reason",
        "data_fit_sideband_accepted",
        "data_fit_sideband_decision",
        "data_fit_fraction_pi0",
        "data_fit_fraction_pi0_stat",
        "data_fit_fraction_dvcs",
        "data_fit_fraction_sideband",
        "data_fit_morph_shift_bins",
        "data_fit_morph_smear_sigma_bins",
        "data_fit_nominal_deviance_per_ndf",
        "data_fit_deviance",
        "data_fit_ndf",
        "data_fit_deviance_per_ndf",
        "data_fit_n_data",
        "sideband_js_signal_edge_near",
        "sideband_js_near_far",
        "sideband_js_signal_edge_far",
        "sideband_n_signal_edge",
        "sideband_n_near",
        "sideband_n_far",
    )

    rows = []

    for bin_key, e_min, e_max in E_PROBE_BINS:
        bin_summary = optimizer_summary.get("energy_bins", {}).get(
            bin_key,
            {},
        )

        for region in ("FT", "FD"):
            region_summary = bin_summary.get("regions", {}).get(region)
            if region_summary is None:
                continue
            #endif

            common = {
                "period": period.key,
                "region": region,
                "energy_bin": bin_key,
                "energy_min": e_min,
                "energy_max": e_max,
                "n_pi0_baseline": region_summary.get(
                    "n_pi0_baseline",
                    "",
                ),
                "n_dvcs_baseline": region_summary.get(
                    "n_dvcs_baseline",
                    "",
                ),
                "n_data_sideband_baseline": region_summary.get(
                    "n_data_sideband_baseline",
                    "",
                ),
            }

            for result in region_summary.get("results", []):
                row = {field: "" for field in fields}
                row.update(common)
                row.update(
                    {
                        "record_type": "optimizer",
                        "step": result["step"],
                        "feature": result["feature"],
                        "direction": result["direction"],
                        "threshold": result["threshold"],
                        "step_pi0_eff": result["step_pi0_eff"],
                        "step_dvcs_eff": result["step_dvcs_eff"],
                        "step_data_sideband_eff": result[
                            "step_data_sideband_eff"
                        ],
                        "total_pi0_eff": result["total_pi0_eff"],
                        "total_dvcs_eff": result["total_dvcs_eff"],
                        "total_data_sideband_eff": result[
                            "total_data_sideband_eff"
                        ],
                        "score": result["score"],
                    }
                )
                rows.append(row)
            #endfor

            for result in region_summary.get("closure_results", []):
                row = {field: "" for field in fields}
                row.update(common)
                row.update(
                    {
                        "record_type": "closure",
                        "step": result["step"],
                        "feature": result["feature"],
                        "js_divergence": result["js_divergence"],
                        "true_fraction": result["true_fraction"],
                        "mean_fit": result["mean_fit"],
                        "bias": result["bias"],
                        "rms_error": result["rms_error"],
                        "std_fit": result["std_fit"],
                        "mean_sigma": result["mean_sigma"],
                        "boundary_fraction": result[
                            "boundary_fraction"
                        ],
                    }
                )
                rows.append(row)
            #endfor

            data_fit = region_summary.get("data_fit")
            if data_fit is not None:
                row = {field: "" for field in fields}
                row.update(common)
                row.update(
                    {
                        "record_type": "data_fit",
                        "operating_step": data_fit.get("step", ""),
                        "operating_step_status": data_fit.get(
                            "operating_step_status",
                            data_fit.get("status", ""),
                        ),
                        "data_fit_feature": data_fit.get("feature", ""),
                        "data_fit_mode": data_fit.get("fit_mode", ""),
                        "data_fit_fallback_reason": data_fit.get(
                            "fallback_reason",
                            "",
                        ),
                        "data_fit_sideband_accepted": data_fit.get(
                            "sideband_accepted",
                            "",
                        ),
                        "data_fit_sideband_decision": data_fit.get(
                            "sideband_decision",
                            "",
                        ),
                        "data_fit_fraction_pi0": data_fit.get(
                            "fraction_pi0",
                            "",
                        ),
                        "data_fit_fraction_pi0_stat": data_fit.get(
                            "fraction_pi0_stat",
                            "",
                        ),
                        "data_fit_fraction_dvcs": data_fit.get(
                            "fraction_dvcs",
                            "",
                        ),
                        "data_fit_fraction_sideband": data_fit.get(
                            "fraction_sideband",
                            "",
                        ),
                        "data_fit_morph_shift_bins": data_fit.get(
                            "morph_shift_bins",
                            "",
                        ),
                        "data_fit_morph_smear_sigma_bins": data_fit.get(
                            "morph_smear_sigma_bins",
                            "",
                        ),
                        "data_fit_nominal_deviance_per_ndf": data_fit.get(
                            "nominal_deviance_per_ndf",
                            "",
                        ),
                        "data_fit_deviance": data_fit.get("deviance", ""),
                        "data_fit_ndf": data_fit.get("ndf", ""),
                        "data_fit_deviance_per_ndf": data_fit.get(
                            "deviance_per_ndf",
                            "",
                        ),
                        "data_fit_n_data": data_fit.get("n_data", ""),
                    }
                )
                stability = region_summary.get(
                    "sideband_stability",
                    {},
                )
                if stability.get("success", False):
                    row.update(
                        {
                            "sideband_js_signal_edge_near": stability[
                                "js_signal_edge_near"
                            ],
                            "sideband_js_near_far": stability[
                                "js_near_far"
                            ],
                            "sideband_js_signal_edge_far": stability[
                                "js_signal_edge_far"
                            ],
                            "sideband_n_signal_edge": stability[
                                "slice_counts"
                            ]["signal_edge"],
                            "sideband_n_near": stability[
                                "slice_counts"
                            ]["near"],
                            "sideband_n_far": stability[
                                "slice_counts"
                            ]["far"],
                        }
                    )
                #endif
                rows.append(row)
            #endif
        #endfor
    #endfor

    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    #endwith

    return out


def make_combined_closure_response_canvas(
    period: Period,
    optimizer_summary: dict,
    output_dir: Path,
) -> Path:
    """
    One 4x2 fitted-vs-injected pi0-fraction canvas per period.

    Rows are E_probe bins and columns are FD/FT.  Each line is a cumulative
    optimizer step.  The dashed diagonal is perfect closure.  This directly
    visualizes whether the AAO/DVCS templates still contain enough shape
    information to determine their relative normalization.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(
        len(E_PROBE_BINS),
        2,
        figsize=(13.8, 4.0 * len(E_PROBE_BINS) + 1.0),
        squeeze=False,
    )

    legend_handles = []
    legend_labels = []

    for irow, (bin_key, e_min, e_max) in enumerate(E_PROBE_BINS):
        bin_summary = optimizer_summary.get("energy_bins", {}).get(
            bin_key,
            {},
        )

        for icol, region in enumerate(("FD", "FT")):
            ax = axes[irow, icol]
            region_summary = bin_summary.get("regions", {}).get(region)

            closure_results = (
                region_summary.get("closure_results", [])
                if region_summary is not None
                else []
            )

            ax.plot(
                [0.0, 1.0],
                [0.0, 1.0],
                linestyle="--",
                linewidth=1.0,
                label="ideal",
            )

            if not closure_results:
                ax.text(
                    0.5,
                    0.5,
                    "Insufficient closure statistics",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                )
            else:
                steps = sorted(
                    {
                        int(result["step"])
                        for result in closure_results
                    }
                )

                for step in steps:
                    rows = [
                        result
                        for result in closure_results
                        if int(result["step"]) == step
                    ]
                    rows.sort(
                        key=lambda result: result["true_fraction"]
                    )

                    line, = ax.plot(
                        [row["true_fraction"] for row in rows],
                        [row["mean_fit"] for row in rows],
                        marker="o",
                        linewidth=1.3,
                        label=f"step {step}",
                    )

                    if irow == 0 and icol == 0:
                        legend_handles.append(line)
                        legend_labels.append(f"step {step}")
                    #endif
                #endfor

                # Also annotate the best remaining 1D discriminator at the
                # final cut step and its JS divergence.
                final_step = max(steps)
                final_rows = [
                    row
                    for row in closure_results
                    if int(row["step"]) == final_step
                ]
                if final_rows:
                    ref = final_rows[0]
                    ax.text(
                        0.03,
                        0.97,
                        (
                            f"final best: {ref['feature']}\n"
                            f"JS={ref['js_divergence']:.3f}"
                        ),
                        transform=ax.transAxes,
                        va="top",
                        ha="left",
                        fontsize=8.0,
                        bbox={
                            "facecolor": "white",
                            "alpha": 0.82,
                            "edgecolor": "0.6",
                        },
                    )
                #endif
            #endif

            ax.set_xlim(0.0, 1.0)
            ax.set_ylim(0.0, 1.0)
            ax.set_xlabel(r"Injected $f_{\pi^0}$")
            ax.set_ylabel(r"Mean fitted $f_{\pi^0}$")
            ax.set_title(
                f"{region}: {e_min:g} <= E_probe < {e_max:g} GeV"
            )
            ax.grid(alpha=0.25)
        #endfor
    #endfor

    if legend_handles:
        fig.legend(
            legend_handles,
            legend_labels,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.975),
            ncol=min(5, len(legend_handles)),
            frameon=True,
        )
    #endif

    fig.suptitle(
        f"{period.label}: template-fit response by cumulative cut step",
        fontsize=13,
        y=0.998,
    )
    fig.subplots_adjust(
        left=0.08,
        right=0.985,
        bottom=0.055,
        top=0.93,
        wspace=0.18,
        hspace=0.42,
    )

    out = output_dir / f"template_fit_response_{period.key}.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def write_optimizer_outputs(
    period: Period,
    region: str,
    energy_bin_key: str,
    energy_min: float,
    energy_max: float,
    pi0_events: np.ndarray,
    dvcs_events: np.ndarray,
    data_sideband_events: np.ndarray,
    results: list,
    output_dir: Path,
) -> None:
    """
    Write detailed per-region CSV/TXT optimizer results.

    Plotting is handled separately so each period gets one combined FD+FT
    progression canvas rather than separate FD and FT PNG files.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = (
        f"rectangular_optimizer_{period.key}_{region.lower()}_"
        f"eprobe_{energy_bin_key}"
    )

    csv_path = output_dir / f"{stem}.csv"
    txt_path = output_dir / f"{stem}.txt"

    fields = (
        "step",
        "feature",
        "direction",
        "threshold",
        "step_pi0_eff",
        "step_dvcs_eff",
        "step_data_sideband_eff",
        "total_pi0_eff",
        "total_dvcs_eff",
        "total_data_sideband_eff",
        "score",
        "n_pi0",
        "n_dvcs",
        "n_data_sideband",
    )

    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for result in results:
            writer.writerow({key: result[key] for key in fields})
        #endfor
    #endwith

    lines = [
        (
            f"{period.label} {region} rectangular-cut optimizer: "
            f"{energy_min:g} <= E_probe < {energy_max:g} GeV"
        ),
        "=" * 72,
        (
            "Baseline: "
            f"AAO={pi0_events.shape[0]:,}, "
            f"DVCSgen={dvcs_events.shape[0]:,}, "
            f"data sideband={data_sideband_events.shape[0]:,}"
        ),
        (
            "Objective: maximize cumulative eps_pi0 / "
            "sqrt(eps_DVCS + lambda*eps_data_sideband), with configured "
            "AAO-retention constraints."
        ),
        "",
    ]

    if not results:
        lines.append("No accepted optimizer step.")
    else:
        for result in results:
            lines.extend(
                [
                    f"Step {result['step']}: {cut_expression(result)}",
                    (
                        f"  step retention: "
                        f"AAO={100.0*result['step_pi0_eff']:.2f}%  "
                        f"DVCSgen={100.0*result['step_dvcs_eff']:.2f}%  "
                        f"data-SB={100.0*result['step_data_sideband_eff']:.2f}%"
                    ),
                    (
                        f"  cumulative:     "
                        f"AAO={100.0*result['total_pi0_eff']:.2f}%  "
                        f"DVCSgen={100.0*result['total_dvcs_eff']:.2f}%  "
                        f"data-SB={100.0*result['total_data_sideband_eff']:.2f}%  "
                        f"F={result['score']:.4g}"
                    ),
                    "",
                ]
            )
        #endfor
    #endif

    txt_path.write_text("\n".join(lines) + "\n")


def make_combined_optimizer_canvas(
    period: Period,
    optimizer_summary: dict,
    output_dir: Path,
) -> Path:
    """
    One optimizer progression canvas per period.

    Rows are E_probe bins; columns are FD and FT. This keeps all eight
    period-specific optimizations in one compact diagnostic.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(
        len(E_PROBE_BINS),
        2,
        figsize=(13.8, 4.0 * len(E_PROBE_BINS) + 1.0),
        sharey=True,
        squeeze=False,
    )

    legend_handles = []
    legend_labels = []

    for irow, (bin_key, e_min, e_max) in enumerate(E_PROBE_BINS):
        bin_summary = optimizer_summary.get("energy_bins", {}).get(bin_key, {})

        for icol, region in enumerate(("FD", "FT")):
            ax = axes[irow, icol]
            region_summary = bin_summary.get("regions", {}).get(region)

            if region_summary is None:
                ax.text(
                    0.5,
                    0.5,
                    "No optimizer result",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                )
                ax.set_title(
                    f"{region}: {e_min:g}-{e_max:g} GeV"
                )
                continue
            #endif

            results = region_summary.get("results", [])
            x = [0] + [int(result["step"]) for result in results]
            pi_eff = [1.0] + [
                result["total_pi0_eff"]
                for result in results
            ]
            dvcs_eff = [1.0] + [
                result["total_dvcs_eff"]
                for result in results
            ]
            side_eff = [1.0] + [
                result["total_data_sideband_eff"]
                for result in results
            ]

            line_pi, = ax.plot(
                x,
                pi_eff,
                marker="o",
                label=r"AAO $\pi^0$",
            )
            line_dvcs, = ax.plot(
                x,
                dvcs_eff,
                marker="o",
                label="DVCSgen",
            )
            line_side, = ax.plot(
                x,
                side_eff,
                marker="o",
                label="data sideband",
            )

            if irow == 0 and icol == 0:
                legend_handles = [line_pi, line_dvcs, line_side]
                legend_labels = [
                    r"AAO $\pi^0$",
                    "DVCSgen",
                    "data sideband",
                ]
            #endif

            ax.set_yscale("log")
            ax.set_ylim(1.0e-4, 1.1)
            ax.set_xticks(x)
            ax.set_xlabel("Greedy cut step")
            ax.set_title(
                f"{region}: {e_min:g} <= E_probe < {e_max:g} GeV\n"
                f"AAO={region_summary['n_pi0_baseline']:,}, "
                f"DVCS={region_summary['n_dvcs_baseline']:,}, "
                f"data-SB={region_summary['n_data_sideband_baseline']:,}",
                fontsize=10.0,
            )
            ax.grid(alpha=0.25)

            if results:
                annotation = "\n".join(
                    f"{r['step']}. {cut_expression(r)}"
                    for r in results
                )
                final = results[-1]
                annotation += (
                    "\n"
                    f"final: AAO {100.0*final['total_pi0_eff']:.1f}%, "
                    f"DVCS {100.0*final['total_dvcs_eff']:.2f}%, "
                    f"data-SB "
                    f"{100.0*final['total_data_sideband_eff']:.2f}%"
                )
                ax.text(
                    0.02,
                    0.03,
                    annotation,
                    transform=ax.transAxes,
                    fontsize=7.8,
                    va="bottom",
                    ha="left",
                    bbox={
                        "facecolor": "white",
                        "alpha": 0.88,
                        "edgecolor": "0.5",
                    },
                )
            #endif
        #endfor
    #endfor

    for ax in axes[:, 0]:
        ax.set_ylabel("Cumulative retained fraction")
    #endfor

    if legend_handles:
        fig.legend(
            legend_handles,
            legend_labels,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.975),
            ncol=3,
            frameon=True,
        )
    #endif

    fig.suptitle(
        f"{period.label}: E_probe-binned rectangular-cut optimization",
        fontsize=13,
        y=0.998,
    )
    fig.subplots_adjust(
        left=0.08,
        right=0.985,
        bottom=0.055,
        top=0.93,
        wspace=0.15,
        hspace=0.42,
    )

    out = output_dir / f"rectangular_optimizer_{period.key}.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


# =============================================================================
# Streaming histogrammer
# =============================================================================

def stream_sample_histograms(
    path: str,
    sample_label: str,
    tree_name: str,
    max_entries: int,
    step_size: int,
    beam_energy_gev: float,
    collect_optimizer: bool = False,
    optimizer_max_events: int = 0,
    optimizer_seed: int = 12345,
) -> Tuple[
    Dict[str, Dict[str, Dict[str, np.ndarray]]],
    Dict[str, int],
    Dict[str, Dict[str, np.ndarray]],
]:
    """
    Read one ROOT sample exactly once and fill tag-FT and tag-FD histograms.

    FT/FD here always refers to the reconstructed tag photon p2.  Predicted
    probe theta/phi/sector are retained independently as event metadata.

    When optimizer collection is enabled, maintain an independent bounded
    uniform reservoir for every detector region AND every E_probe bin. This
    avoids undersampling a sparse energy bin when another energy region is much
    more populated, while still requiring only one ROOT-file pass.
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

    region_seed_offset = {"FT": 11, "FD": 29}
    reservoirs = {
        region: {
            bin_key: PriorityReservoir(
                optimizer_max_events,
                len(EVENT_COLUMNS),
                optimizer_seed
                + region_seed_offset[region]
                + 1000 * ibin,
            )
            for ibin, (bin_key, _emin, _emax) in enumerate(E_PROBE_BINS)
        }
        for region in ("FT", "FD")
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
            # Reconstructed p2 is the TAG photon.  Its theta defines the
            # FT/FD composition category.
            tag_theta = photon_theta_deg(
                arrays["p2_theta"],
                angle_unit,
            )

            # The missing momentum predicts the PARTNER/PROBE photon.  Its
            # direction is retained separately for the eventual efficiency
            # binning and FD-sector diagnostics.
            pred_probe_theta, pred_probe_phi = (
                predicted_probe_direction_deg(
                    arrays,
                    angle_unit,
                    beam_energy_gev,
                )
            )
            pred_probe_sector = fd_sector_from_phi_deg(
                pred_probe_phi
            )

            base = np.isfinite(theta_egamma) & (
                theta_egamma > THETA_EGAMMA_MIN_DEG
            )
            counts["opening_angle"] += int(np.count_nonzero(base))

            # Composition categories are defined by the RECONSTRUCTED TAG
            # photon p2.  Predicted probe direction is deliberately NOT used
            # here; it is a separate downstream efficiency coordinate.
            region_masks = {
                "FT": (
                    base
                    & np.isfinite(tag_theta)
                    & (tag_theta >= FT_THETA_MIN_DEG)
                    & (tag_theta < FT_THETA_MAX_DEG)
                ),
                "FD": (
                    base
                    & np.isfinite(tag_theta)
                    & (tag_theta >= FD_THETA_MIN_DEG)
                    & (tag_theta < FD_THETA_MAX_DEG)
                ),
            }

            values = {
                "Mx2": np.asarray(arrays["Mx2"], dtype=float),
                "Mx2_1": np.asarray(arrays["Mx2_1"], dtype=float),
                "Mx2_2": np.asarray(arrays["Mx2_2"], dtype=float),
                "Emiss2": np.asarray(arrays["Emiss2"], dtype=float),
                "E_tag": np.asarray(arrays["p2_p"], dtype=float),
                "Delta_phi_residual_deg": delta_phi_residual_deg(
                    arrays["Delta_phi"]
                ),
                "pTmiss": np.asarray(arrays["pTmiss"], dtype=float),
                "theta_gamma_gamma": np.asarray(
                    arrays["theta_gamma_gamma"],
                    dtype=float,
                ),
                "tag_theta_deg": tag_theta,
                "pred_probe_sector": pred_probe_sector,
                "pred_probe_theta_deg": pred_probe_theta,
                "pred_probe_phi_deg": pred_probe_phi,
            }

            e_probe = values["Emiss2"]

            for region, region_mask in region_masks.items():
                counts[region] += int(np.count_nonzero(region_mask))

                if collect_optimizer:
                    for bin_key, e_min, e_max in E_PROBE_BINS:
                        energy_mask = (
                            region_mask
                            & np.isfinite(e_probe)
                            & (e_probe >= e_min)
                            & (e_probe < e_max)
                        )
                        reservoirs[region][bin_key].update(
                            optimizer_matrix(values, energy_mask)
                        )
                    #endfor
                #endif

                # Retain the old fixed diagnostic histograms only for the
                # --no-optimize-cuts fallback path.
                mx2_1_mask = (
                    region_mask
                    & np.isfinite(values["Mx2_1"])
                    & (values["Mx2_1"] < 0.15)
                )

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
                    for key in OPTIMIZER_FEATURES:
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

            if (
                counts["read"] == entry_stop
                or counts["read"] % max(step_size, 1_000_000) == 0
            ):
                frac = (
                    100.0 * counts["read"] / entry_stop
                    if entry_stop
                    else 100.0
                )
                log(
                    f"{sample_label}: {counts['read']:,}/{entry_stop:,} "
                    f"({frac:.1f}%)"
                )
            #endif
        #endfor
    #endwith

    if collect_optimizer:
        energy_counts = []
        for region in ("FT", "FD"):
            pieces = []
            for bin_key, e_min, e_max in E_PROBE_BINS:
                n_bin = reservoirs[region][bin_key].array().shape[0]
                pieces.append(
                    f"{e_min:g}-{e_max:g} GeV={n_bin:,}"
                )
            #endfor
            energy_counts.append(f"{region}: " + ", ".join(pieces))
        #endfor
        log(
            f"{sample_label}: optimizer reservoirs after "
            f"Angle(e,gamma)>{THETA_EGAMMA_MIN_DEG:g} deg -> "
            + " | ".join(energy_counts)
        )
    else:
        log(
            f"{sample_label}: retained after Angle(e,gamma)>5 deg: "
            f"FT={counts['FT']:,}, FD={counts['FD']:,}."
        )
    #endif

    optimizer_events = {
        region: {
            bin_key: reservoirs[region][bin_key].array()
            for bin_key, _e_min, _e_max in E_PROBE_BINS
        }
        for region in ("FT", "FD")
    }
    return hist, counts, optimizer_events


# =============================================================================
# Dynamic post-optimizer shape-comparison support
# =============================================================================

FEATURE_INDEX = {
    feature: i
    for i, feature in enumerate(OPTIMIZER_FEATURES)
}


EVENT_INDEX = {
    column: i
    for i, column in enumerate(EVENT_COLUMNS)
}


def apply_optimizer_cut(
    events: np.ndarray,
    current_mask: np.ndarray,
    result: dict,
) -> np.ndarray:
    """Apply one selected rectangular cut cumulatively to one sample."""
    feature = result["feature"]
    threshold = float(result["threshold"])
    values = events[:, FEATURE_INDEX[feature]]

    if result["direction"] == "lt":
        return current_mask & np.isfinite(values) & (values < threshold)
    #endif

    return current_mask & np.isfinite(values) & (values > threshold)


def optimizer_row_label(result: dict) -> str:
    """Compact human-readable label for one optimizer-selected cut."""
    op = "<" if result["direction"] == "lt" else ">"
    feature = result["feature"]

    pretty = {
        "Mx2": r"$Mx2$",
        "Mx2_1": r"$Mx2_1$",
        "Mx2_2": r"$Mx2_2$",
        "Emiss2": r"$E_{\rm miss}$",
        "E_tag": r"$E_{\rm tag}$",
        "Delta_phi_residual_deg": r"$\Delta\phi(p,\gamma)$",
        "pTmiss": r"$p_{T,\rm miss}$",
        "theta_gamma_gamma": r"$\theta_{\gamma\gamma}$",
    }[feature]

    return rf"$+\,{pretty}{op}{result['threshold']:.3g}$"


def make_optimizer_cut_shape_canvas(
    period: Period,
    region: str,
    energy_bin_key: str,
    energy_min: float,
    energy_max: float,
    sample_events: Dict[str, Dict[str, Dict[str, np.ndarray]]],
    optimizer_results: list,
    output_dir: Path,
    linear_y: bool,
    mx2_1_preselection_max: float,
) -> Path:
    """
    Build one post-optimizer shape canvas for one period, detector region, and
    E_probe bin.

    Row 1 is the energy-bin baseline plus the loose Mx2_1 preselection. Each
    subsequent row cumulatively applies the optimizer-selected cut sequence.
    """
    n_rows = 1 + len(optimizer_results)
    n_cols = len(PLOT_SPECS)

    fig_height = max(4.2, 3.25 * n_rows + 1.8)
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(27.5, fig_height),
        squeeze=False,
    )

    sample_meta = {key: (label, color) for key, label, color in SAMPLES}
    edges_by_key = histogram_edges()

    masks_by_sample = {}
    counts_by_sample = {}

    for sample_key, _label, _color in SAMPLES:
        events = sample_events[sample_key][region][energy_bin_key]
        mx2_1_values = events[:, FEATURE_INDEX["Mx2_1"]]
        baseline_mask = (
            np.isfinite(mx2_1_values)
            & (mx2_1_values < mx2_1_preselection_max)
        )
        masks = [baseline_mask]

        current = masks[0]
        for result in optimizer_results:
            current = apply_optimizer_cut(events, current, result)
            masks.append(current.copy())
        #endfor

        masks_by_sample[sample_key] = masks
        counts_by_sample[sample_key] = [
            int(np.count_nonzero(mask))
            for mask in masks
        ]
    #endfor

    legend_handles = []
    legend_labels = []

    row_labels = [
        (
            rf"${energy_min:g}\leq E_{{\rm probe}}<{energy_max:g}$ GeV"
            + "\n"
            + rf"$Mx2_1<{mx2_1_preselection_max:.3g}$ GeV$^2$"
        )
    ]
    row_labels.extend(
        [
            f"step {i}: {optimizer_row_label(result)}"
            for i, result in enumerate(optimizer_results, start=1)
        ]
    )

    for irow in range(n_rows):
        for icol, spec in enumerate(PLOT_SPECS):
            ax = axes[irow, icol]
            key, title, xlabel, _xmin, _xmax, _nbins = spec
            edges = edges_by_key[key]
            centers = 0.5 * (edges[:-1] + edges[1:])
            ifeature = FEATURE_INDEX[key]

            for sample_key, _label, _color in SAMPLES:
                label, color = sample_meta[sample_key]
                events = sample_events[sample_key][region][energy_bin_key]
                mask = masks_by_sample[sample_key][irow]

                values = events[:, ifeature]
                finite = mask & np.isfinite(values)
                h, _ = np.histogram(values[finite], bins=edges)

                y = normalized_hist(h)
                line = ax.step(
                    centers,
                    y,
                    where="mid",
                    linewidth=1.45,
                    color=color,
                    label=label,
                )[0]

                if irow == 0 and icol == 0:
                    legend_handles.append(line)
                    legend_labels.append(label)
                #endif
            #endfor

            if irow == 0:
                ax.set_title(title, fontsize=10.5)
            #endif

            ax.set_xlabel(xlabel, fontsize=8.8)
            if icol == 0:
                ax.set_ylabel(
                    row_labels[irow] + "\nunit-normalized",
                    fontsize=8.8,
                )
            else:
                ax.set_ylabel("unit-normalized", fontsize=8.8)
            #endif

            ax.set_xlim(edges[0], edges[-1])

            if not linear_y:
                ax.set_yscale("log")
                ax.set_ylim(1.0e-3, 1.0)
            #endif

            ax.tick_params(axis="both", labelsize=7.4)
            ax.grid(alpha=0.18)

            if key == "Mx2_1":
                ax.axvline(
                    mx2_1_preselection_max,
                    linestyle="--",
                    linewidth=0.9,
                    color="0.35",
                )
            #endif

            # Emiss2 is the E_probe binning variable, not an optimizer cut.
            if key == "Emiss2":
                ax.axvline(
                    energy_min,
                    linestyle=":",
                    linewidth=1.0,
                    color="0.35",
                )
                ax.axvline(
                    energy_max,
                    linestyle=":",
                    linewidth=1.0,
                    color="0.35",
                )
            #endif

            for applied_result in optimizer_results[:irow]:
                if applied_result["feature"] != key:
                    continue
                #endif
                ax.axvline(
                    float(applied_result["threshold"]),
                    linestyle="--",
                    linewidth=0.9,
                    color="0.35",
                )
            #endfor
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

    final_count_text = ", ".join(
        f"{sample_meta[key][0]}={counts_by_sample[key][-1]:,}"
        for key, _label, _color in SAMPLES
    )

    cut_sequence = (
        " -> ".join(cut_expression(result) for result in optimizer_results)
        if optimizer_results
        else "no optimizer cut accepted"
    )

    fig.suptitle(
        f"{period.label}: {region}, "
        rf"${energy_min:g}\leq E_{{\rm probe}}<{energy_max:g}$ GeV\n"
        rf"baseline: $\mathrm{{Angle}}(e,\gamma)>"
        rf"{THETA_EGAMMA_MIN_DEG:.0f}^\circ$; "
        + rf"$Mx2_1<{mx2_1_preselection_max:.3g}$ GeV$^2$; "
        + region_text
        + "\n"
        + f"optimizer sequence: {cut_sequence}"
        + "\n"
        + f"final sampled counts: {final_count_text}",
        fontsize=10.5,
        y=0.995,
    )

    fig.legend(
        legend_handles,
        legend_labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.915),
        ncol=3,
        fontsize=9.2,
        frameon=True,
        fancybox=False,
        edgecolor="black",
    )

    fig.subplots_adjust(
        left=0.047,
        right=0.995,
        bottom=0.055,
        top=0.83,
        wspace=0.35,
        hspace=0.48,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    out = (
        output_dir
        / (
            f"canvas_shape_comparison_{period.key}_{region.lower()}_"
            f"eprobe_{energy_bin_key}_optimized.png"
        )
    )
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out



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
    linear_y: bool,
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

            if not linear_y:
                # The histograms are unit-normalized. A small positive floor
                # keeps empty bins from causing log-scale warnings while not
                # distorting any nonzero bin content.
                ax.set_yscale("log")
                ax.set_ylim(1.0e-3, 1.0)
            #endif

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
        + ("\nLinear y scale" if linear_y else "\nLogarithmic y scale")
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
# Raw nSidis epgamma <-> eppi0 run/event coverage audit
# =============================================================================

def load_eppi0_coverage_reference(
    path: str,
    tree_name: str,
    step_size: int,
    label: str,
) -> dict:
    """
    Load the SMALLER eppi0 run/event set once.

    The coverage question we actually need to answer is whether each eppi0
    event exists somewhere in the much larger nSidis epgamma file.  We
    therefore sort only the eppi0 keys and never construct/sort the full
    epgamma key array.
    """
    found_tree, total = preflight_custom_branches(
        path,
        tree_name,
        ("runnum", "evnum"),
    )
    log(
        f"{label}: coverage reference reading full eppi0 tree "
        f"({total:,} entries) from '{found_tree}'."
    )

    key_parts: List[np.ndarray] = []
    run_parts: List[np.ndarray] = []

    for arrays in iterate_efficiency_tree_arrays(
        path,
        found_tree,
        ("runnum", "evnum"),
        0,
        step_size,
    ):
        run = np.asarray(arrays["runnum"], dtype=np.int64)
        ev = np.asarray(arrays["evnum"], dtype=np.int64)
        key_parts.append(packed_event_keys_eff(run, ev))
        run_parts.append(run)
    #endfor

    if not key_parts:
        return {
            "keys": np.asarray([], dtype=np.uint64),
            "runs": np.asarray([], dtype=np.int64),
            "n_rows": 0,
        }
    #endif

    keys = np.concatenate(key_parts)
    runs = np.concatenate(run_parts)
    n_rows = int(len(keys))

    # Sorting ~10^6 eppi0 keys is cheap.  This is the ONLY global key sort in
    # the coverage audit.
    order = np.argsort(keys, kind="mergesort")
    keys = keys[order]
    runs = runs[order]

    unique = np.ones(len(keys), dtype=bool)
    if len(keys) > 1:
        unique[1:] = keys[1:] != keys[:-1]
    #endif

    return {
        "keys": keys[unique],
        "runs": runs[unique],
        "n_rows": n_rows,
    }


def raw_run_event_coverage(
    period: Period,
    tree_name: str,
    step_size: int,
) -> dict:
    """
    Fast FULL-file exact-event coverage audit.

    Algorithm:
      1. Load/sort/uniquify only the relatively small eppi0 key set.
      2. Stream the huge epgamma tree chunk-by-chunk.
      3. For each epgamma key, binary-search the sorted eppi0 keys.
      4. Mark the corresponding eppi0 event as found.
      5. Stop early if every eppi0 event has been found.

    We never concatenate or sort the tens-of-millions-entry epgamma tree.
    File ordering is irrelevant.
    """
    reference = load_eppi0_coverage_reference(
        period.eppi0_data,
        tree_name,
        step_size,
        f"{period.label} eppi0",
    )
    pi_keys = np.asarray(reference["keys"], dtype=np.uint64)
    pi_runs = np.asarray(reference["runs"], dtype=np.int64)

    if len(pi_keys) == 0:
        return {
            "period": period.key,
            "rows": [],
            "n_epgamma_rows_scanned": 0,
            "n_epgamma_rows_total": 0,
            "n_eppi0_rows": int(reference["n_rows"]),
            "n_eppi0_unique": 0,
            "n_common_unique": 0,
            "common_over_eppi0": float("nan"),
            "n_epgamma_runs_seen": 0,
            "n_eppi0_runs": 0,
            "n_common_runs": 0,
            "epgamma_only_runs": [],
            "eppi0_only_runs": [],
            "stopped_early": False,
        }
    #endif

    pi_run_values, pi_run_counts = np.unique(
        pi_runs,
        return_counts=True,
    )
    pi_count_map = dict(
        zip(pi_run_values.tolist(), pi_run_counts.tolist())
    )

    matched_pi = np.zeros(len(pi_keys), dtype=bool)

    epg_found_tree, epg_total = preflight_custom_branches(
        period.data,
        tree_name,
        ("runnum", "evnum"),
    )
    log(
        f"{period.label} epgamma: fast coverage scan of full tree "
        f"({epg_total:,} rows) from '{epg_found_tree}'; "
        f"searching against {len(pi_keys):,} unique eppi0 events."
    )

    # Only run-level epgamma row counts are accumulated.  There is no reason
    # to globally uniquify the enormous epgamma sample for the eppi0-coverage
    # question.
    epg_rows_by_run: Dict[int, int] = {}
    epg_runs_seen: set = set()
    rows_scanned = 0
    stopped_early = False
    t_scan = time.perf_counter()

    for arrays in iterate_efficiency_tree_arrays(
        period.data,
        epg_found_tree,
        ("runnum", "evnum"),
        0,
        step_size,
    ):
        run = np.asarray(arrays["runnum"], dtype=np.int64)
        ev = np.asarray(arrays["evnum"], dtype=np.int64)
        keys = packed_event_keys_eff(run, ev)
        rows_scanned += len(keys)

        chunk_runs, chunk_counts = np.unique(
            run,
            return_counts=True,
        )
        for run_value, count in zip(chunk_runs, chunk_counts):
            run_i = int(run_value)
            epg_runs_seen.add(run_i)
            epg_rows_by_run[run_i] = (
                epg_rows_by_run.get(run_i, 0) + int(count)
            )
        #endfor

        # Binary search the sorted SMALL reference array.
        pos = np.searchsorted(pi_keys, keys, side="left")
        in_bounds = pos < len(pi_keys)
        if np.any(in_bounds):
            epg_idx = np.flatnonzero(in_bounds)
            candidate_pos = pos[in_bounds]
            exact = pi_keys[candidate_pos] == keys[epg_idx]
            if np.any(exact):
                matched_pi[candidate_pos[exact]] = True
            #endif
        #endif

        if np.all(matched_pi):
            stopped_early = True
            log(
                f"{period.label}: all {len(pi_keys):,} unique eppi0 events "
                f"found after scanning {rows_scanned:,}/{epg_total:,} "
                f"epgamma rows; stopping coverage scan early."
            )
            break
        #endif

        # Lightweight progress every ~10M scanned rows.
        if rows_scanned % 10_000_000 < len(keys):
            elapsed = max(time.perf_counter() - t_scan, 1.0e-9)
            found = int(np.count_nonzero(matched_pi))
            log(
                f"{period.label}: coverage scan {rows_scanned:,}/"
                f"{epg_total:,} epgamma rows; "
                f"found {found:,}/{len(pi_keys):,} eppi0 events "
                f"({100.0 * found / len(pi_keys):.2f}%); "
                f"{rows_scanned / elapsed / 1.0e6:.2f} M rows/s."
            )
        #endif
    #endfor

    common_keys = pi_keys[matched_pi]
    common_runs = pi_runs[matched_pi]
    common_run_values, common_run_counts = np.unique(
        common_runs,
        return_counts=True,
    )
    common_count_map = dict(
        zip(common_run_values.tolist(), common_run_counts.tolist())
    )

    all_runs = np.union1d(
        np.asarray(sorted(epg_runs_seen), dtype=np.int64),
        pi_run_values,
    )

    rows = []
    for run_value in all_runs:
        run_i = int(run_value)
        n_epg_rows = int(epg_rows_by_run.get(run_i, 0))
        n_pi = int(pi_count_map.get(run_i, 0))
        n_common = int(common_count_map.get(run_i, 0))

        rows.append(
            {
                "run": run_i,
                "epgamma_rows_scanned": n_epg_rows,
                "eppi0_unique": n_pi,
                "common_unique": n_common,
                "common_over_eppi0": (
                    n_common / n_pi if n_pi > 0 else float("nan")
                ),
            }
        )
    #endfor

    n_common = int(np.count_nonzero(matched_pi))
    common_run_set = set(common_run_values.astype(int).tolist())
    epg_run_set = set(int(x) for x in epg_runs_seen)
    pi_run_set = set(pi_run_values.astype(int).tolist())

    return {
        "period": period.key,
        "rows": rows,
        "n_epgamma_rows_scanned": int(rows_scanned),
        "n_epgamma_rows_total": int(epg_total),
        "n_eppi0_rows": int(reference["n_rows"]),
        "n_eppi0_unique": int(len(pi_keys)),
        "n_common_unique": n_common,
        "common_over_eppi0": (
            n_common / len(pi_keys) if len(pi_keys) else float("nan")
        ),
        "n_epgamma_runs_seen": int(len(epg_run_set)),
        "n_eppi0_runs": int(len(pi_run_set)),
        "n_common_runs": int(len(common_run_set)),
        "epgamma_only_runs": sorted(epg_run_set - pi_run_set),
        "eppi0_only_runs": sorted(pi_run_set - epg_run_set),
        "stopped_early": bool(stopped_early),
    }


def make_raw_coverage_canvas(
    period: Period,
    coverage: dict,
    output_dir: Path,
) -> Path:
    """
    Compact full-file coverage diagnostic.

    The relevant physics/production metric is common/eppi0: essentially every
    event reconstructed in eppi0 should occur somewhere in the more inclusive
    nSidis epgamma tree.  We deliberately do not globally unique/sort epgamma.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = coverage["rows"]

    runs = np.asarray([r["run"] for r in rows], dtype=int)
    epg_rows = np.asarray(
        [r["epgamma_rows_scanned"] for r in rows],
        dtype=float,
    )
    pi = np.asarray(
        [r["eppi0_unique"] for r in rows],
        dtype=float,
    )
    common = np.asarray(
        [r["common_unique"] for r in rows],
        dtype=float,
    )
    frac_pi = np.asarray(
        [r["common_over_eppi0"] for r in rows],
        dtype=float,
    )

    fig = plt.figure(figsize=(14.5, 9.0))
    gs = fig.add_gridspec(
        3,
        1,
        height_ratios=(2.5, 2.1, 1.25),
        hspace=0.30,
        left=0.075,
        right=0.985,
        bottom=0.07,
        top=0.92,
    )

    ax_counts = fig.add_subplot(gs[0, 0])
    ax_frac = fig.add_subplot(gs[1, 0], sharex=ax_counts)
    ax_text = fig.add_subplot(gs[2, 0])

    ax_counts.step(
        runs,
        epg_rows,
        where="mid",
        linewidth=1.2,
        label="epgamma rows scanned / run",
    )
    ax_counts.step(
        runs,
        pi,
        where="mid",
        linewidth=1.3,
        label="eppi0 unique events / run",
    )
    ax_counts.step(
        runs,
        common,
        where="mid",
        linewidth=1.3,
        linestyle="--",
        label="eppi0 events found in epgamma",
    )
    ax_counts.set_ylabel("Events / run")
    ax_counts.set_yscale("log")
    ax_counts.grid(alpha=0.18)
    ax_counts.legend(frameon=True)

    ax_frac.plot(
        runs,
        frac_pi,
        "o",
        markersize=3.2,
        label="common / eppi0",
    )
    ax_frac.axhline(1.0, linestyle="--", linewidth=0.9)
    ax_frac.set_ylim(-0.03, 1.05)
    ax_frac.set_ylabel("eppi0 coverage in epgamma")
    ax_frac.set_xlabel("Run number")
    ax_frac.grid(alpha=0.18)
    ax_frac.legend(frameon=True)

    ax_text.axis("off")
    epg_only = coverage["epgamma_only_runs"]
    pi_only = coverage["eppi0_only_runs"]
    scan_text = (
        f"{coverage['n_epgamma_rows_scanned']:,}/"
        f"{coverage['n_epgamma_rows_total']:,}"
    )
    lines = [
        (
            f"epgamma rows scanned: {scan_text}"
            + (" (early stop: all eppi0 found)" if coverage["stopped_early"] else "")
        ),
        (
            f"eppi0: {coverage['n_eppi0_rows']:,} rows -> "
            f"{coverage['n_eppi0_unique']:,} unique events, "
            f"{coverage['n_eppi0_runs']} runs"
        ),
        (
            f"found in epgamma: {coverage['n_common_unique']:,} unique eppi0 "
            f"events in {coverage['n_common_runs']} runs"
        ),
        (
            f"common/eppi0 = {coverage['common_over_eppi0']:.3%}"
        ),
        (
            "epgamma-only runs seen: "
            + (
                ", ".join(map(str, epg_only[:30]))
                + (" ..." if len(epg_only) > 30 else "")
                if epg_only else "none"
            )
        ),
        (
            "eppi0-only runs: "
            + (
                ", ".join(map(str, pi_only[:30]))
                + (" ..." if len(pi_only) > 30 else "")
                if pi_only else "none"
            )
        ),
    ]
    ax_text.text(
        0.01,
        0.98,
        "\n".join(lines),
        ha="left",
        va="top",
        fontsize=10,
        family="monospace",
        transform=ax_text.transAxes,
    )

    fig.suptitle(
        (
            f"{period.label}: fast raw nSidis eppi0 coverage in epgamma\n"
            "FULL files; exact (runnum,evnum); no physics cuts; no global epgamma sort"
        ),
        fontsize=13,
    )

    out = output_dir / f"event_coverage_{period.key}.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out



def run_raw_coverage_audit(
    period: Period,
    tree_name: str,
    step_size: int,
    output_root: Path,
    min_eppi0_overlap: float,
    fail_on_pathology: bool,
) -> dict:
    coverage = raw_run_event_coverage(period, tree_name, step_size)
    outdir = output_root / "coverage" / period.key
    plot = make_raw_coverage_canvas(period, coverage, outdir)

    log(
        (
            f"{period.label}: RAW COVERAGE: "
            f"epgamma rows scanned={coverage['n_epgamma_rows_scanned']:,}/"
            f"{coverage['n_epgamma_rows_total']:,}; "
            f"eppi0={coverage['n_eppi0_unique']:,} unique events/"
            f"{coverage['n_eppi0_runs']} runs; "
            f"found={coverage['n_common_unique']:,}/"
            f"{coverage['n_common_runs']} runs; "
            f"common/eppi0={coverage['common_over_eppi0']:.3%}; "
            f"early_stop={coverage['stopped_early']}. "
            f"Plot: {plot}"
        )
    )

    pathological = (
        coverage["n_common_runs"] == 0
        or coverage["n_common_unique"] == 0
        or not np.isfinite(coverage["common_over_eppi0"])
        or coverage["common_over_eppi0"] < min_eppi0_overlap
    )

    if pathological and fail_on_pathology:
        raise RuntimeError(
            "\n".join(
                [
                    f"{period.label}: PATHOLOGICAL raw data-file coverage.",
                    (
                        f"  common/eppi0 = "
                        f"{coverage['common_over_eppi0']:.3%}; "
                        f"required >= {min_eppi0_overlap:.3%}"
                    ),
                    (
                        f"  common runs = {coverage['n_common_runs']} "
                        f"(epgamma runs seen {coverage['n_epgamma_runs_seen']}, "
                        f"eppi0 {coverage['n_eppi0_runs']})"
                    ),
                    (
                        "  epgamma-only runs seen: "
                        + (
                            ", ".join(
                                map(str, coverage["epgamma_only_runs"][:40])
                            )
                            if coverage["epgamma_only_runs"]
                            else "none"
                        )
                    ),
                    (
                        "  eppi0-only runs: "
                        + (
                            ", ".join(
                                map(str, coverage["eppi0_only_runs"][:40])
                            )
                            if coverage["eppi0_only_runs"]
                            else "none"
                        )
                    ),
                    (
                        "Efficiency extraction is aborted before matching. "
                        "Verify that both ROOT files are the corresponding "
                        "complete nSidis productions."
                    ),
                ]
            )
        )
    #endif
    return coverage


# =============================================================================
# Final tag-probe efficiency extraction
# =============================================================================

@dataclass
class EfficiencyEPPi0Store:
    electron_p3: np.ndarray
    proton_p3: np.ndarray
    pi0_p3: np.ndarray
    pi0_mass: np.ndarray
    runnum: Optional[np.ndarray]
    evnum: Optional[np.ndarray]
    angle_unit: str


@dataclass
class EfficiencyResolution:
    center: np.ndarray
    sigma: np.ndarray
    n: int
    source: str


def spherical_to_cartesian(
    p: np.ndarray,
    theta: np.ndarray,
    phi: np.ndarray,
    unit: str,
) -> np.ndarray:
    """Vectorized spherical -> Cartesian momentum."""
    theta_r, phi_r = _angles_to_radians(theta, phi, unit)
    p = np.asarray(p, dtype=float)
    st = np.sin(theta_r)
    return np.column_stack(
        (
            p * st * np.cos(phi_r),
            p * st * np.sin(phi_r),
            p * np.cos(theta_r),
        )
    )


def efficiency_predicted_probe(
    arrays: Dict[str, np.ndarray],
    unit: str,
    beam_energy_gev: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Return (pred_E, pred_p3, pred_theta_deg) from exact reconstructed
    four-vector missing kinematics.
    """
    e3 = spherical_to_cartesian(
        arrays["e_p"], arrays["e_theta"], arrays["e_phi"], unit
    )
    p3 = spherical_to_cartesian(
        arrays["p1_p"], arrays["p1_theta"], arrays["p1_phi"], unit
    )
    tag3 = spherical_to_cartesian(
        arrays["p2_p"], arrays["p2_theta"], arrays["p2_phi"], unit
    )

    e_p = np.linalg.norm(e3, axis=1)
    p_p = np.linalg.norm(p3, axis=1)
    tag_e = np.linalg.norm(tag3, axis=1)
    e_e = np.sqrt(e_p * e_p + M_E_GEV * M_E_GEV)
    p_e = np.sqrt(p_p * p_p + M_P_GEV * M_P_GEV)

    beam_p = math.sqrt(
        max(0.0, beam_energy_gev * beam_energy_gev - M_E_GEV * M_E_GEV)
    )
    initial_p3 = np.asarray([0.0, 0.0, beam_p], dtype=float)

    pred_e = beam_energy_gev + M_P_GEV - e_e - p_e - tag_e
    pred3 = initial_p3[None, :] - e3 - p3 - tag3
    pred_p = np.linalg.norm(pred3, axis=1)

    with np.errstate(invalid="ignore", divide="ignore"):
        pred_theta = np.degrees(
            np.arccos(
                np.clip(
                    pred3[:, 2] / pred_p,
                    -1.0,
                    1.0,
                )
            )
        )
    #endwith
    return pred_e, pred3, pred_theta


def efficiency_epgamma_feature_matrix(
    arrays: Dict[str, np.ndarray],
    beam_energy_gev: float,
) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    """
    Build exactly the same purity-cut variables as the streaming concept stage,
    plus explicit tag/probe vectors used by the efficiency association.
    """
    unit = infer_angle_unit(
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
        unit,
    )
    tag_theta = photon_theta_deg(arrays["p2_theta"], unit)

    pred_e, pred3, pred_theta = efficiency_predicted_probe(
        arrays,
        unit,
        beam_energy_gev,
    )
    pred_phi = np.mod(
        np.degrees(np.arctan2(pred3[:, 1], pred3[:, 0])),
        360.0,
    )
    pred_sector = fd_sector_from_phi_deg(pred_phi)

    values = {
        "Mx2": np.asarray(arrays["Mx2"], dtype=float),
        "Mx2_1": np.asarray(arrays["Mx2_1"], dtype=float),
        "Mx2_2": np.asarray(arrays["Mx2_2"], dtype=float),
        "Emiss2": np.asarray(arrays["Emiss2"], dtype=float),
        "E_tag": np.asarray(arrays["p2_p"], dtype=float),
        "Delta_phi_residual_deg": delta_phi_residual_deg(
            arrays["Delta_phi"]
        ),
        "pTmiss": np.asarray(arrays["pTmiss"], dtype=float),
        "theta_gamma_gamma": np.asarray(
            arrays["theta_gamma_gamma"], dtype=float
        ),
        "tag_theta_deg": tag_theta,
        "pred_probe_sector": pred_sector,
        "pred_probe_theta_deg": pred_theta,
        "pred_probe_phi_deg": pred_phi,
    }
    matrix = np.column_stack(
        [
            np.asarray(values[key], dtype=np.float32)
            for key in EVENT_COLUMNS
        ]
    )

    extra = {
        "unit": unit,
        "theta_egamma_deg": theta_egamma,
        "pred_probe_energy_4vec": pred_e,
        "pred_probe_p3": pred3,
        "electron_p3": spherical_to_cartesian(
            arrays["e_p"], arrays["e_theta"], arrays["e_phi"], unit
        ),
        "proton_p3": spherical_to_cartesian(
            arrays["p1_p"], arrays["p1_theta"], arrays["p1_phi"], unit
        ),
        "tag_p3": spherical_to_cartesian(
            arrays["p2_p"], arrays["p2_theta"], arrays["p2_phi"], unit
        ),
        "tag_theta_deg": tag_theta,
        "pred_probe_theta_deg": pred_theta,
        "pred_probe_phi_deg": pred_phi,
        "pred_probe_sector": pred_sector,
    }
    return matrix, extra


def efficiency_probe_bin(
    energy_bin_key: str,
    energy_min: float,
    pred_theta_deg: np.ndarray,
    pred_sector: np.ndarray,
) -> np.ndarray:
    """
    Vector string labels for the downstream efficiency coordinate.

    FT is kept as one bin.  FD uses sectors 1-6 below 4 GeV and FD_all in the
    highest energy bin.
    """
    theta = np.asarray(pred_theta_deg, dtype=float)
    sector = np.rint(np.asarray(pred_sector, dtype=float)).astype(int)
    labels = np.full(theta.shape, "", dtype="<U12")

    ft = (
        np.isfinite(theta)
        & (theta >= FT_THETA_MIN_DEG)
        & (theta < FT_THETA_MAX_DEG)
    )
    fd = (
        np.isfinite(theta)
        & (theta >= FD_THETA_MIN_DEG)
        & (theta < FD_THETA_MAX_DEG)
    )
    labels[ft] = "FT"

    if energy_min >= HIGH_ENERGY_FD_COMBINE_MIN_GEV:
        labels[fd] = "FD_all"
    else:
        for s in range(1, 7):
            labels[fd & (sector == s)] = f"FD_S{s}"
        #endfor
    #endif
    return labels


def efficiency_bin_order(
    energy_min: float,
) -> Tuple[str, ...]:
    if energy_min >= HIGH_ENERGY_FD_COMBINE_MIN_GEV:
        return ("FT", "FD_all")
    #endif
    return ("FT", "FD_S1", "FD_S2", "FD_S3", "FD_S4", "FD_S5", "FD_S6")


def preflight_custom_branches(
    path: str,
    tree_name: str,
    required: Sequence[str],
) -> Tuple[str, int]:
    """Preflight a ROOT tree against an explicit branch list."""
    if not Path(path).exists():
        raise FileNotFoundError(
            f"Required ROOT file does not exist: {path}\n"
            "For the data efficiency numerator, this must be the nSidis "
            "eppi0 companion production corresponding to the nSidis epgamma "
            "denominator file."
        )
    #endif
    with uproot.open(path) as root_file:
        tree, found = get_tree(root_file, tree_name)
        available = set(tree.keys())
        missing = [name for name in required if name not in available]
        if missing:
            raise RuntimeError(
                f"{path}: tree '{found}' missing required efficiency branches: "
                + ", ".join(missing)
            )
        #endif
        return found, int(tree.num_entries)
    #endwith


def iterate_efficiency_tree_arrays(
    path: str,
    tree_name: str,
    expressions: Sequence[str],
    max_entries: int,
    step_size: int,
):
    """Yield bounded NumPy chunks from one ROOT tree."""
    with uproot.open(path) as root_file:
        tree, _ = get_tree(root_file, tree_name)
        total = int(tree.num_entries)
        entry_stop = total if max_entries <= 0 else min(total, int(max_entries))
        for arrays in tree.iterate(
            expressions=list(expressions),
            entry_start=0,
            entry_stop=entry_stop,
            step_size=step_size,
            library="np",
        ):
            yield arrays
        #endfor
    #endwith


def load_eppi0_efficiency_store(
    path: str,
    tree_name: str,
    max_entries: int,
    step_size: int,
    require_event_ids: bool,
    label: str,
) -> EfficiencyEPPi0Store:
    """Load only the compact reconstructed eppi0 information needed to match."""
    required = list(EFF_EPPIO_REQUIRED)
    if require_event_ids:
        required.extend(EFF_EPPIO_DATA_IDS)
    #endif

    found_tree, total = preflight_custom_branches(
        path, tree_name, required
    )
    n_read = total if max_entries <= 0 else min(total, max_entries)
    log(
        f"{label}: loading eppi0 match store: {n_read:,}/{total:,} entries "
        f"from '{found_tree}'."
    )

    pieces: Dict[str, List[np.ndarray]] = {name: [] for name in required}
    for arrays in iterate_efficiency_tree_arrays(
        path,
        found_tree,
        required,
        max_entries,
        step_size,
    ):
        for name in required:
            pieces[name].append(np.asarray(arrays[name]))
        #endfor
    #endfor

    arr = {
        name: np.concatenate(values) if values else np.asarray([])
        for name, values in pieces.items()
    }
    unit = infer_angle_unit(
        arr["e_theta"], arr["e_phi"], arr["p2_theta"], arr["p2_phi"]
    )

    store = EfficiencyEPPi0Store(
        electron_p3=spherical_to_cartesian(
            arr["e_p"], arr["e_theta"], arr["e_phi"], unit
        ).astype(np.float32, copy=False),
        proton_p3=spherical_to_cartesian(
            arr["p1_p"], arr["p1_theta"], arr["p1_phi"], unit
        ).astype(np.float32, copy=False),
        pi0_p3=spherical_to_cartesian(
            arr["p2_p"], arr["p2_theta"], arr["p2_phi"], unit
        ).astype(np.float32, copy=False),
        pi0_mass=np.asarray(arr["Mh_gammagamma"], dtype=np.float32),
        runnum=(
            np.asarray(arr["runnum"], dtype=np.int64)
            if require_event_ids else None
        ),
        evnum=(
            np.asarray(arr["evnum"], dtype=np.int64)
            if require_event_ids else None
        ),
        angle_unit=unit,
    )
    return store


def packed_event_keys_eff(
    runnum: np.ndarray,
    evnum: np.ndarray,
) -> np.ndarray:
    run = np.asarray(runnum, dtype=np.int64)
    ev = np.asarray(evnum, dtype=np.int64)
    if np.any(run < 0) or np.any(ev < 0):
        raise ValueError("Negative data runnum/evnum encountered.")
    #endif
    return (
        (run.astype(np.uint64) << np.uint64(32))
        | ev.astype(np.uint64)
    )


def clean_reconstructed_probe_from_pairs(
    electron_p3: np.ndarray,
    proton_p3: np.ndarray,
    tag_p3: np.ndarray,
    pred_probe_p3: np.ndarray,
    eppi0: EfficiencyEPPi0Store,
    eppi0_index: np.ndarray,
    mgg_min: float,
    mgg_max: float,
    remainder_mass2_max: float,
    reco_probe_energy_min: float,
) -> Dict[str, np.ndarray]:
    """
    Reconstruct P_probe = P_pi0 - k_tag for already-associated parent pairs.
    """
    j = np.asarray(eppi0_index, dtype=np.int64)
    pi3 = np.asarray(eppi0.pi0_p3[j], dtype=float)
    pi_mass = np.asarray(eppi0.pi0_mass[j], dtype=float)
    tag3 = np.asarray(tag_p3, dtype=float)

    pi_p = np.linalg.norm(pi3, axis=1)
    tag_e = np.linalg.norm(tag3, axis=1)
    pi_e = np.sqrt(np.maximum(0.0, pi_p * pi_p + pi_mass * pi_mass))

    reco_e = pi_e - tag_e
    reco3 = pi3 - tag3
    reco_p = np.linalg.norm(reco3, axis=1)
    reco_m2 = reco_e * reco_e - reco_p * reco_p
    delta_p3 = reco3 - np.asarray(pred_probe_p3, dtype=float)

    finite = (
        np.isfinite(pi_mass)
        & np.isfinite(reco_e)
        & np.isfinite(reco_p)
        & np.all(np.isfinite(delta_p3), axis=1)
    )
    positive = finite & (reco_e > 0.0) & (reco_p > 0.0)
    mass_window = (
        positive
        & (pi_mass >= mgg_min)
        & (pi_mass <= mgg_max)
    )
    mass_shell = (
        mass_window
        & np.isfinite(reco_m2)
        & (np.abs(reco_m2) <= remainder_mass2_max)
    )
    threshold = mass_shell & (reco_e >= reco_probe_energy_min)

    return {
        "pi0_mass": pi_mass,
        "reco_probe_energy": reco_e,
        "reco_probe_mass2": reco_m2,
        "delta_p3": delta_p3,
        "positive": positive,
        "mass_window": mass_window,
        "mass_shell": mass_shell,
        "threshold": threshold,
    }


def build_mc_parent_index(
    eppi0: EfficiencyEPPi0Store,
    component_tolerance: float,
) -> Tuple[cKDTree, np.ndarray]:
    good = (
        np.all(np.isfinite(eppi0.electron_p3), axis=1)
        & np.all(np.isfinite(eppi0.proton_p3), axis=1)
    )
    indices = np.flatnonzero(good)
    search = np.empty((len(indices), 6), dtype=np.float32)
    search[:, :3] = (
        eppi0.electron_p3[indices] / component_tolerance
    ).astype(np.float32, copy=False)
    search[:, 3:] = (
        eppi0.proton_p3[indices] / component_tolerance
    ).astype(np.float32, copy=False)
    tree = cKDTree(search, compact_nodes=True, balanced_tree=True)
    return tree, indices


def match_mc_parent_chunk(
    electron_p3: np.ndarray,
    proton_p3: np.ndarray,
    tree: cKDTree,
    tree_indices: np.ndarray,
    eppi0: EfficiencyEPPi0Store,
    component_tolerance: float,
    distance_max: float,
    workers: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Match MC ONLY through reconstructed electron/proton kinematics.

    No runnum or evnum enters this function.
    """
    n = len(electron_p3)
    if n == 0:
        return (
            np.asarray([], dtype=np.int64),
            np.asarray([], dtype=np.int64),
            np.asarray([], dtype=float),
        )
    #endif

    query = np.empty((n, 6), dtype=np.float32)
    query[:, :3] = (
        np.asarray(electron_p3, dtype=np.float32) / component_tolerance
    )
    query[:, 3:] = (
        np.asarray(proton_p3, dtype=np.float32) / component_tolerance
    )
    distance, local = tree.query(
        query,
        k=1,
        workers=max(1, int(workers)),
    )
    candidate = tree_indices[local]

    de = np.abs(
        np.asarray(electron_p3, dtype=float)
        - np.asarray(eppi0.electron_p3[candidate], dtype=float)
    )
    dp = np.abs(
        np.asarray(proton_p3, dtype=float)
        - np.asarray(eppi0.proton_p3[candidate], dtype=float)
    )
    max_component = np.maximum(np.max(de, axis=1), np.max(dp, axis=1))
    accepted = (
        np.isfinite(distance)
        & (distance <= distance_max)
        & (max_component <= component_tolerance)
    )
    epg_index = np.flatnonzero(accepted)
    return epg_index, candidate[accepted], max_component[accepted]


def prepare_data_event_index(
    eppi0: EfficiencyEPPi0Store,
) -> Tuple[np.ndarray, np.ndarray]:
    if eppi0.runnum is None or eppi0.evnum is None:
        raise ValueError("Data eppi0 store has no run/event identifiers.")
    #endif
    keys = packed_event_keys_eff(eppi0.runnum, eppi0.evnum)
    order = np.argsort(keys, kind="mergesort")
    return keys[order], order


def match_data_parent_chunk(
    runnum: np.ndarray,
    evnum: np.ndarray,
    electron_p3: np.ndarray,
    proton_p3: np.ndarray,
    sorted_eppi0_keys: np.ndarray,
    sorted_eppi0_order: np.ndarray,
    eppi0: EfficiencyEPPi0Store,
    component_tolerance: float,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, int]]:
    """
    DATA association:
      1) exact runnum/evnum event identity,
      2) reconstructed e/p parent consistency.

    The best same-event candidate is the one with the smallest six-component
    e/p distance. This remains independent of the probe-momentum residual.
    """
    keys = packed_event_keys_eff(runnum, evnum)
    left = np.searchsorted(sorted_eppi0_keys, keys, side="left")
    right = np.searchsorted(sorted_eppi0_keys, keys, side="right")
    counts = right - left
    same_event = counts > 0

    total_pairs = int(np.sum(counts))
    if total_pairs == 0:
        return (
            np.asarray([], dtype=np.int64),
            np.asarray([], dtype=np.int64),
            {
                "selected_tags": int(len(keys)),
                "same_event": 0,
                "same_event_pairs": 0,
                "parent_consistent": 0,
            },
        )
    #endif

    epg_rep = np.repeat(np.arange(len(keys), dtype=np.int64), counts)
    left_rep = np.repeat(left.astype(np.int64), counts)
    prefix = np.repeat((np.cumsum(counts) - counts), counts)
    offsets = np.arange(total_pairs, dtype=np.int64) - prefix
    pi_idx = sorted_eppi0_order[left_rep + offsets]

    de = np.abs(
        np.asarray(electron_p3[epg_rep], dtype=float)
        - np.asarray(eppi0.electron_p3[pi_idx], dtype=float)
    )
    dp = np.abs(
        np.asarray(proton_p3[epg_rep], dtype=float)
        - np.asarray(eppi0.proton_p3[pi_idx], dtype=float)
    )
    max_component = np.maximum(np.max(de, axis=1), np.max(dp, axis=1))
    good_parent = np.isfinite(max_component) & (
        max_component <= component_tolerance
    )

    eligible = np.flatnonzero(good_parent)
    if len(eligible) == 0:
        return (
            np.asarray([], dtype=np.int64),
            np.asarray([], dtype=np.int64),
            {
                "selected_tags": int(len(keys)),
                "same_event": int(np.count_nonzero(same_event)),
                "same_event_pairs": total_pairs,
                "parent_consistent": 0,
            },
        )
    #endif

    # One candidate per epgamma tag, chosen only from parent consistency.
    score = max_component[eligible]
    g = epg_rep[eligible]
    order = np.lexsort((score, g))
    eligible = eligible[order]
    g = epg_rep[eligible]
    first = np.ones(len(eligible), dtype=bool)
    if len(eligible) > 1:
        first[1:] = g[1:] != g[:-1]
    #endif
    chosen = eligible[first]

    return (
        epg_rep[chosen],
        pi_idx[chosen],
        {
            "selected_tags": int(len(keys)),
            "same_event": int(np.count_nonzero(same_event)),
            "same_event_pairs": total_pairs,
            "parent_consistent": int(len(chosen)),
        },
    )


def robust_resolution_model(
    delta_p3: np.ndarray,
    min_events: int,
    source: str,
) -> EfficiencyResolution:
    values = np.asarray(delta_p3, dtype=float)
    values = values[np.all(np.isfinite(values), axis=1)]
    if len(values) == 0:
        return EfficiencyResolution(
            center=np.zeros(3, dtype=float),
            sigma=np.full(3, np.nan, dtype=float),
            n=0,
            source=source,
        )
    #endif

    center = np.nanmedian(values, axis=0)
    q_lo = np.nanpercentile(values, 15.865, axis=0)
    q_hi = np.nanpercentile(values, 84.135, axis=0)
    sigma = 0.5 * (q_hi - q_lo)

    # Guard against pathological zero-width quantization.
    sigma = np.maximum(sigma, 1.0e-4)
    return EfficiencyResolution(
        center=np.asarray(center, dtype=float),
        sigma=np.asarray(sigma, dtype=float),
        n=int(len(values)),
        source=source,
    )


def apply_resolution_match(
    delta_p3: np.ndarray,
    model: EfficiencyResolution,
    nsigma: float,
) -> np.ndarray:
    delta = np.asarray(delta_p3, dtype=float)
    if (
        model.n <= 0
        or not np.all(np.isfinite(model.center))
        or not np.all(np.isfinite(model.sigma))
    ):
        return np.zeros(len(delta), dtype=bool)
    #endif
    z = np.abs(
        (delta - model.center[None, :])
        / model.sigma[None, :]
    )
    return np.all(np.isfinite(z) & (z <= nsigma), axis=1)


def initialize_efficiency_accumulator() -> Dict[str, dict]:
    out: Dict[str, dict] = {}
    for key, emin, emax in E_PROBE_BINS:
        for probe_bin in efficiency_bin_order(emin):
            out[f"{key}:{probe_bin}"] = {
                "energy_bin": key,
                "energy_min": float(emin),
                "energy_max": float(emax),
                "probe_bin": probe_bin,
                "data_raw_tags": 0,
                "data_expected_pi0": 0.0,
                "data_expected_variance": 0.0,
                "data_matched": 0,
                "data_tagFT_raw": 0,
                "data_tagFD_raw": 0,
                "data_tagFT_fpi0": float("nan"),
                "data_tagFT_fpi0_stat": float("nan"),
                "data_tagFD_fpi0": float("nan"),
                "data_tagFD_fpi0_stat": float("nan"),
                "mc_denominator": 0,
                "mc_matched": 0,
                "mc_tagFT_raw": 0,
                "mc_tagFD_raw": 0,
            }
        #endfor
    #endfor
    return out


def efficiency_selected_indices(
    matrix: np.ndarray,
    theta_egamma_deg: np.ndarray,
    tag_theta_deg: np.ndarray,
    region: str,
    energy_min: float,
    energy_max: float,
    optimizer_results: list,
    operating_step: int,
) -> np.ndarray:
    """Apply exactly the purity-fit signal selection to one epgamma chunk."""
    emiss_idx = FEATURE_INDEX["Emiss2"]
    mx2_idx = FEATURE_INDEX["Mx2_1"]

    base = (
        np.isfinite(theta_egamma_deg)
        & (theta_egamma_deg > THETA_EGAMMA_MIN_DEG)
        & np.isfinite(matrix[:, emiss_idx])
        & (matrix[:, emiss_idx] >= energy_min)
        & (matrix[:, emiss_idx] < energy_max)
        & np.isfinite(matrix[:, mx2_idx])
        & (matrix[:, mx2_idx] < SIGNAL_MX2_1_MAX)
    )

    if region == "FT":
        base &= (
            np.isfinite(tag_theta_deg)
            & (tag_theta_deg >= FT_THETA_MIN_DEG)
            & (tag_theta_deg < FT_THETA_MAX_DEG)
        )
    else:
        base &= (
            np.isfinite(tag_theta_deg)
            & (tag_theta_deg >= FD_THETA_MIN_DEG)
            & (tag_theta_deg < FD_THETA_MAX_DEG)
        )
    #endif

    idx = np.flatnonzero(base)
    if len(idx) == 0:
        return idx
    #endif

    subset = matrix[idx]
    masks = build_cumulative_optimizer_masks(
        subset,
        optimizer_results,
    )
    step = min(int(operating_step), len(masks) - 1)
    return idx[masks[step]]


def period_efficiency_fit_config(
    optimizer_summary: dict,
    bin_key: str,
    tag_region: str,
) -> Optional[dict]:
    region_summary = (
        optimizer_summary.get("energy_bins", {})
        .get(bin_key, {})
        .get("regions", {})
        .get(tag_region)
    )
    if region_summary is None:
        return None
    #endif
    fit = region_summary.get("data_fit")
    if not fit or not fit.get("success", False):
        return None
    #endif
    return {
        "fit": fit,
        "optimizer_results": region_summary.get("results", []),
        "step": int(fit["step"]),
    }


def efficiency_resolution_region(
    pred_theta_deg: np.ndarray,
) -> np.ndarray:
    theta = np.asarray(pred_theta_deg, dtype=float)
    region = np.full(theta.shape, "", dtype="<U3")
    region[
        np.isfinite(theta)
        & (theta >= FT_THETA_MIN_DEG)
        & (theta < FT_THETA_MAX_DEG)
    ] = "FT"
    region[
        np.isfinite(theta)
        & (theta >= FD_THETA_MIN_DEG)
        & (theta < FD_THETA_MAX_DEG)
    ] = "FD"
    return region


def accumulate_selected_tags(
    accumulator: Dict[str, dict],
    bin_key: str,
    energy_min: float,
    probe_labels: np.ndarray,
    tag_region: str,
    is_data: bool,
    f_pi0: float = 1.0,
    f_pi0_stat: float = 0.0,
) -> None:
    for probe_bin in efficiency_bin_order(energy_min):
        count = int(np.count_nonzero(probe_labels == probe_bin))
        if count == 0:
            continue
        #endif
        row = accumulator[f"{bin_key}:{probe_bin}"]
        if is_data:
            row["data_raw_tags"] += count
            row["data_expected_pi0"] += f_pi0 * count
            count_key = (
                "data_tagFT_raw" if tag_region == "FT" else "data_tagFD_raw"
            )
            row[count_key] += count
            prefix = "data_tagFT" if tag_region == "FT" else "data_tagFD"
            row[f"{prefix}_fpi0"] = float(f_pi0)
            row[f"{prefix}_fpi0_stat"] = float(f_pi0_stat)
        else:
            row["mc_denominator"] += count
            row[
                "mc_tagFT_raw" if tag_region == "FT" else "mc_tagFD_raw"
            ] += count
        #endif
    #endfor


def accumulate_matched_tags(
    accumulator: Dict[str, dict],
    bin_key: str,
    energy_min: float,
    probe_labels: np.ndarray,
    matched_mask: np.ndarray,
    is_data: bool,
) -> None:
    for probe_bin in efficiency_bin_order(energy_min):
        count = int(
            np.count_nonzero(
                (probe_labels == probe_bin)
                & np.asarray(matched_mask, dtype=bool)
            )
        )
        if count == 0:
            continue
        #endif
        row = accumulator[f"{bin_key}:{probe_bin}"]
        if is_data:
            row["data_matched"] += count
        else:
            row["mc_matched"] += count
        #endif
    #endfor


def process_mc_efficiency_pass(
    period: Period,
    optimizer_summary: dict,
    tree_name: str,
    max_entries: int,
    step_size: int,
    eppi0: EfficiencyEPPi0Store,
    parent_component_tol: float,
    parent_distance_max: float,
    kdtree_workers: int,
    assoc_mgg_min: float,
    assoc_mgg_max: float,
    assoc_remainder_mass2_max: float,
    assoc_probe_energy_min: float,
    accumulator: Dict[str, dict],
) -> Dict[str, object]:
    """
    One additional AAO epgamma pass:
      denominator counting + e/p cKDTree parent association + clean residuals.
    """
    tree, tree_indices = build_mc_parent_index(
        eppi0,
        parent_component_tol,
    )

    found_tree, total = preflight_custom_branches(
        period.pi0_mc,
        tree_name,
        EFF_EPG_REQUIRED,
    )
    residual_parts = {"FT": [], "FD": [], "ALL": []}
    clean_records: List[dict] = []
    counters = {
        "selected_tags": 0,
        "parent_matched": 0,
        "positive": 0,
        "mass_window": 0,
        "mass_shell": 0,
        "threshold": 0,
    }
    matched_parent_indices_all: List[np.ndarray] = []

    for arrays in iterate_efficiency_tree_arrays(
        period.pi0_mc,
        found_tree,
        EFF_EPG_REQUIRED,
        max_entries,
        step_size,
    ):
        matrix, extra = efficiency_epgamma_feature_matrix(
            arrays,
            BEAM_ENERGY_GEV[period.key],
        )

        for bin_key, emin, emax in E_PROBE_BINS:
            for tag_region in ("FT", "FD"):
                config = period_efficiency_fit_config(
                    optimizer_summary,
                    bin_key,
                    tag_region,
                )
                if config is None:
                    continue
                #endif

                selected = efficiency_selected_indices(
                    matrix,
                    extra["theta_egamma_deg"],
                    extra["tag_theta_deg"],
                    tag_region,
                    emin,
                    emax,
                    config["optimizer_results"],
                    config["step"],
                )
                if len(selected) == 0:
                    continue
                #endif

                probe_labels = efficiency_probe_bin(
                    bin_key,
                    emin,
                    extra["pred_probe_theta_deg"][selected],
                    matrix[selected, EVENT_INDEX["pred_probe_sector"]],
                )
                valid_probe = probe_labels != ""
                selected = selected[valid_probe]
                probe_labels = probe_labels[valid_probe]
                if len(selected) == 0:
                    continue
                #endif

                counters["selected_tags"] += len(selected)
                accumulate_selected_tags(
                    accumulator,
                    bin_key,
                    emin,
                    probe_labels,
                    tag_region,
                    is_data=False,
                )

                local_epg, pi_idx, _parent_delta = match_mc_parent_chunk(
                    extra["electron_p3"][selected],
                    extra["proton_p3"][selected],
                    tree,
                    tree_indices,
                    eppi0,
                    parent_component_tol,
                    parent_distance_max,
                    kdtree_workers,
                )
                if len(local_epg) == 0:
                    continue
                #endif
                counters["parent_matched"] += len(local_epg)
                matched_parent_indices_all.append(pi_idx)

                selected_m = selected[local_epg]
                probe_labels_m = probe_labels[local_epg]

                reco = clean_reconstructed_probe_from_pairs(
                    extra["electron_p3"][selected_m],
                    extra["proton_p3"][selected_m],
                    extra["tag_p3"][selected_m],
                    extra["pred_probe_p3"][selected_m],
                    eppi0,
                    pi_idx,
                    assoc_mgg_min,
                    assoc_mgg_max,
                    assoc_remainder_mass2_max,
                    assoc_probe_energy_min,
                )

                for stage in ("positive", "mass_window", "mass_shell", "threshold"):
                    counters[stage] += int(np.count_nonzero(reco[stage]))
                #endfor

                clean = np.asarray(reco["threshold"], dtype=bool)
                if not np.any(clean):
                    continue
                #endif

                delta = np.asarray(reco["delta_p3"][clean], dtype=float)
                theta = extra["pred_probe_theta_deg"][selected_m][clean]
                res_region = efficiency_resolution_region(theta)
                residual_parts["ALL"].append(delta)
                for rr in ("FT", "FD"):
                    if np.any(res_region == rr):
                        residual_parts[rr].append(delta[res_region == rr])
                    #endif
                #endfor

                clean_records.append(
                    {
                        "bin_key": bin_key,
                        "energy_min": emin,
                        "probe_labels": probe_labels_m[clean],
                        "resolution_region": res_region,
                        "delta_p3": delta,
                    }
                )
            #endfor
        #endfor
    #endfor

    residuals = {}
    for region, parts in residual_parts.items():
        residuals[region] = (
            np.concatenate(parts, axis=0)
            if parts else np.empty((0, 3), dtype=float)
        )
    #endfor

    # Directed-tag duplication audit: how often multiple selected epgamma tags
    # found the same AAO eppi0 parent.
    duplicate_fraction = float("nan")
    if matched_parent_indices_all:
        all_idx = np.concatenate(matched_parent_indices_all)
        if len(all_idx):
            _, multiplicity = np.unique(all_idx, return_counts=True)
            duplicate_fraction = float(
                np.sum(multiplicity[multiplicity > 1])
                / max(np.sum(multiplicity), 1)
            )
        #endif
    #endif

    return {
        "residuals": residuals,
        "clean_records": clean_records,
        "counters": counters,
        "directed_tag_duplicate_fraction": duplicate_fraction,
        "tree_entries_total": int(total),
    }


def derive_efficiency_resolutions(
    mc_pass: Dict[str, object],
    min_events: int,
) -> Dict[str, EfficiencyResolution]:
    all_model = robust_resolution_model(
        mc_pass["residuals"]["ALL"],
        min_events,
        "period-wide AAO",
    )
    models = {}
    for region in ("FT", "FD"):
        values = mc_pass["residuals"][region]
        if len(values) >= min_events:
            models[region] = robust_resolution_model(
                values,
                min_events,
                f"{region}-specific AAO",
            )
        else:
            models[region] = EfficiencyResolution(
                center=all_model.center.copy(),
                sigma=all_model.sigma.copy(),
                n=int(len(values)),
                source=f"period-wide AAO fallback ({len(values)} local)",
            )
        #endif
    #endfor
    return models


def finalize_mc_numerator(
    accumulator: Dict[str, dict],
    mc_pass: Dict[str, object],
    models: Dict[str, EfficiencyResolution],
    nsigma: float,
) -> Dict[str, int]:
    matched_total = 0
    clean_total = 0

    for record in mc_pass["clean_records"]:
        delta = record["delta_p3"]
        region = record["resolution_region"]
        matched = np.zeros(len(delta), dtype=bool)
        for rr in ("FT", "FD"):
            mask = region == rr
            if np.any(mask):
                matched[mask] = apply_resolution_match(
                    delta[mask],
                    models[rr],
                    nsigma,
                )
            #endif
        #endfor

        clean_total += len(delta)
        matched_total += int(np.count_nonzero(matched))
        accumulate_matched_tags(
            accumulator,
            record["bin_key"],
            record["energy_min"],
            record["probe_labels"],
            matched,
            is_data=False,
        )
    #endfor
    return {
        "clean_associations": int(clean_total),
        "probe_matched": int(matched_total),
    }


def process_data_efficiency_pass(
    period: Period,
    optimizer_summary: dict,
    tree_name: str,
    max_entries: int,
    step_size: int,
    eppi0: EfficiencyEPPi0Store,
    parent_component_tol: float,
    assoc_mgg_min: float,
    assoc_mgg_max: float,
    assoc_remainder_mass2_max: float,
    assoc_probe_energy_min: float,
    models: Dict[str, EfficiencyResolution],
    nsigma: float,
    accumulator: Dict[str, dict],
) -> Dict[str, object]:
    """
    One additional nSidis epgamma pass:
      purity-weighted denominator + exact-event/e-p association + probe match.
    """
    sorted_keys, sorted_order = prepare_data_event_index(eppi0)

    required = list(EFF_EPG_REQUIRED) + list(EFF_EPG_DATA_IDS)
    found_tree, total = preflight_custom_branches(
        period.data,
        tree_name,
        required,
    )
    counters = {
        "selected_tags": 0,
        "same_event": 0,
        "same_event_pairs": 0,
        "parent_consistent": 0,
        "positive": 0,
        "mass_window": 0,
        "mass_shell": 0,
        "threshold": 0,
        "probe_matched": 0,
        "purity_fit_unavailable_tags": 0,
    }
    residual_parts = {"FT": [], "FD": []}
    selected_event_key_parts: List[np.ndarray] = []

    # Keep only compact packed event keys for the selected-tag diagnostic.
    # There are just 4 E_probe bins x 2 tag regions, so this remains small
    # compared with the ROOT data and does not require another file pass.
    selected_keys_by_category: Dict[Tuple[str, str], List[np.ndarray]] = {
        (bin_key, tag_region): []
        for bin_key, _emin, _emax in E_PROBE_BINS
        for tag_region in ("FT", "FD")
    }

    for arrays in iterate_efficiency_tree_arrays(
        period.data,
        found_tree,
        required,
        max_entries,
        step_size,
    ):
        matrix, extra = efficiency_epgamma_feature_matrix(
            arrays,
            BEAM_ENERGY_GEV[period.key],
        )

        for bin_key, emin, emax in E_PROBE_BINS:
            for tag_region in ("FT", "FD"):
                config = period_efficiency_fit_config(
                    optimizer_summary,
                    bin_key,
                    tag_region,
                )
                if config is None:
                    continue
                #endif
                fit = config["fit"]

                selected = efficiency_selected_indices(
                    matrix,
                    extra["theta_egamma_deg"],
                    extra["tag_theta_deg"],
                    tag_region,
                    emin,
                    emax,
                    config["optimizer_results"],
                    config["step"],
                )
                if len(selected) == 0:
                    continue
                #endif

                probe_labels = efficiency_probe_bin(
                    bin_key,
                    emin,
                    extra["pred_probe_theta_deg"][selected],
                    matrix[selected, EVENT_INDEX["pred_probe_sector"]],
                )
                valid_probe = probe_labels != ""
                selected = selected[valid_probe]
                probe_labels = probe_labels[valid_probe]
                if len(selected) == 0:
                    continue
                #endif

                counters["selected_tags"] += len(selected)
                selected_keys = packed_event_keys_eff(
                    arrays["runnum"][selected],
                    arrays["evnum"][selected],
                )
                selected_event_key_parts.append(selected_keys)
                selected_keys_by_category[(bin_key, tag_region)].append(
                    selected_keys
                )

                accumulate_selected_tags(
                    accumulator,
                    bin_key,
                    emin,
                    probe_labels,
                    tag_region,
                    is_data=True,
                    f_pi0=float(fit["fraction_pi0"]),
                    f_pi0_stat=float(fit.get("fraction_pi0_stat", 0.0)),
                )

                local_epg, pi_idx, assoc_counts = match_data_parent_chunk(
                    arrays["runnum"][selected],
                    arrays["evnum"][selected],
                    extra["electron_p3"][selected],
                    extra["proton_p3"][selected],
                    sorted_keys,
                    sorted_order,
                    eppi0,
                    parent_component_tol,
                )
                counters["same_event"] += assoc_counts["same_event"]
                counters["same_event_pairs"] += assoc_counts["same_event_pairs"]
                counters["parent_consistent"] += assoc_counts["parent_consistent"]

                if len(local_epg) == 0:
                    continue
                #endif

                selected_m = selected[local_epg]
                probe_labels_m = probe_labels[local_epg]

                reco = clean_reconstructed_probe_from_pairs(
                    extra["electron_p3"][selected_m],
                    extra["proton_p3"][selected_m],
                    extra["tag_p3"][selected_m],
                    extra["pred_probe_p3"][selected_m],
                    eppi0,
                    pi_idx,
                    assoc_mgg_min,
                    assoc_mgg_max,
                    assoc_remainder_mass2_max,
                    assoc_probe_energy_min,
                )
                for stage in ("positive", "mass_window", "mass_shell", "threshold"):
                    counters[stage] += int(np.count_nonzero(reco[stage]))
                #endfor

                clean = np.asarray(reco["threshold"], dtype=bool)
                if not np.any(clean):
                    continue
                #endif

                delta = np.asarray(reco["delta_p3"][clean], dtype=float)
                theta = extra["pred_probe_theta_deg"][selected_m][clean]
                rr = efficiency_resolution_region(theta)
                probe_labels_c = probe_labels_m[clean]

                matched = np.zeros(len(delta), dtype=bool)
                for region_name in ("FT", "FD"):
                    mask = rr == region_name
                    if np.any(mask):
                        matched[mask] = apply_resolution_match(
                            delta[mask],
                            models[region_name],
                            nsigma,
                        )
                        # Keep a bounded diagnostic sample.
                        residual_parts[region_name].append(
                            delta[mask][:50_000]
                        )
                    #endif
                #endfor

                counters["probe_matched"] += int(np.count_nonzero(matched))
                accumulate_matched_tags(
                    accumulator,
                    bin_key,
                    emin,
                    probe_labels_c,
                    matched,
                    is_data=True,
                )
            #endfor
        #endfor
    #endfor

    # ------------------------------------------------------------------
    # Selected-tag event-overlap diagnostic.
    #
    # This is intentionally evaluated at the EXACT point where the purity
    # selection has finished, but BEFORE e/p parent checks or reconstructed
    # pi0/probe requirements.  It therefore answers the narrow question:
    #
    #   "Do the physical events selected by the epgamma purity machinery
    #    actually exist in the companion eppi0 file?"
    #
    # Row ordering cannot matter because this is a packed-key set lookup.
    # ------------------------------------------------------------------
    eppi0_unique_keys = np.unique(sorted_keys)
    selected_overlap_summary = []

    for bin_key, emin, emax in E_PROBE_BINS:
        for tag_region in ("FT", "FD"):
            parts = selected_keys_by_category[(bin_key, tag_region)]
            if parts:
                selected_rows_keys = np.concatenate(parts)
            else:
                selected_rows_keys = np.asarray([], dtype=np.uint64)
            #endif

            n_rows = int(len(selected_rows_keys))
            if n_rows:
                unique_keys, multiplicity = np.unique(
                    selected_rows_keys,
                    return_counts=True,
                )
                positions = np.searchsorted(
                    eppi0_unique_keys,
                    unique_keys,
                    side="left",
                )
                in_bounds = positions < len(eppi0_unique_keys)
                in_eppi0 = np.zeros(len(unique_keys), dtype=bool)
                if np.any(in_bounds):
                    valid_idx = np.flatnonzero(in_bounds)
                    in_eppi0[valid_idx] = (
                        eppi0_unique_keys[positions[valid_idx]]
                        == unique_keys[valid_idx]
                    )
                #endif

                n_unique = int(len(unique_keys))
                n_overlap_unique = int(np.count_nonzero(in_eppi0))
                n_overlap_rows = int(np.sum(multiplicity[in_eppi0]))
                mean_multiplicity = float(
                    n_rows / max(n_unique, 1)
                )
                repeated_event_fraction = float(
                    np.count_nonzero(multiplicity > 1)
                    / max(n_unique, 1)
                )
                max_multiplicity = int(np.max(multiplicity))
            else:
                n_unique = 0
                n_overlap_unique = 0
                n_overlap_rows = 0
                mean_multiplicity = float("nan")
                repeated_event_fraction = float("nan")
                max_multiplicity = 0
            #endif

            config = period_efficiency_fit_config(
                optimizer_summary,
                bin_key,
                tag_region,
            )
            fit = config["fit"] if config is not None else None
            f_pi0 = (
                float(fit["fraction_pi0"])
                if fit is not None
                else float("nan")
            )

            selected_overlap_summary.append(
                {
                    "energy_bin": bin_key,
                    "energy_min": float(emin),
                    "energy_max": float(emax),
                    "tag_region": tag_region,
                    "selected_rows": n_rows,
                    "selected_unique_events": n_unique,
                    "overlap_unique_events": n_overlap_unique,
                    "overlap_rows": n_overlap_rows,
                    "overlap_unique_fraction": (
                        n_overlap_unique / n_unique
                        if n_unique > 0
                        else float("nan")
                    ),
                    "overlap_row_fraction": (
                        n_overlap_rows / n_rows
                        if n_rows > 0
                        else float("nan")
                    ),
                    "mean_rows_per_unique_event": mean_multiplicity,
                    "repeated_event_fraction": repeated_event_fraction,
                    "max_rows_per_event": max_multiplicity,
                    "f_pi0_fit": f_pi0,
                }
            )
        #endfor
    #endfor

    duplicate_fraction = float("nan")
    if selected_event_key_parts:
        keys = np.concatenate(selected_event_key_parts)
        if len(keys):
            _, multiplicity = np.unique(keys, return_counts=True)
            duplicate_fraction = float(
                np.sum(multiplicity[multiplicity > 1])
                / max(np.sum(multiplicity), 1)
            )
        #endif
    #endif

    residuals = {}
    for rr in ("FT", "FD"):
        residuals[rr] = (
            np.concatenate(residual_parts[rr], axis=0)
            if residual_parts[rr]
            else np.empty((0, 3), dtype=float)
        )
    #endfor

    return {
        "counters": counters,
        "residuals": residuals,
        "selected_overlap_summary": selected_overlap_summary,
        "directed_tag_event_multiplicity_fraction": duplicate_fraction,
        "tree_entries_total": int(total),
    }


def calculate_efficiency_results(
    accumulator: Dict[str, dict],
) -> List[dict]:
    rows: List[dict] = []

    for key, row in accumulator.items():
        data_den = float(row["data_expected_pi0"])
        data_num = int(row["data_matched"])
        mc_den = int(row["mc_denominator"])
        mc_num = int(row["mc_matched"])

        data_eff = data_num / data_den if data_den > 0.0 else float("nan")
        mc_eff = mc_num / mc_den if mc_den > 0 else float("nan")

        # Keep the common f_pi0 uncertainty correlated across ROOT chunks.
        data_den_variance = 0.0
        for prefix in ("data_tagFT", "data_tagFD"):
            n_tag = int(row[f"{prefix}_raw"])
            f = float(row[f"{prefix}_fpi0"])
            sf = float(row[f"{prefix}_fpi0_stat"])
            if n_tag > 0 and np.isfinite(f):
                data_den_variance += f * f * n_tag
                if np.isfinite(sf):
                    data_den_variance += (n_tag * sf) ** 2
                #endif
            #endif
        #endfor
        row["data_expected_variance"] = data_den_variance
        data_den_sigma = math.sqrt(max(data_den_variance, 0.0))
        if data_num > 0 and data_den > 0:
            data_eff_sigma = data_eff * math.sqrt(
                1.0 / data_num
                + (data_den_sigma / data_den) ** 2
            )
        else:
            data_eff_sigma = float("nan")
        #endif

        if mc_den > 0 and np.isfinite(mc_eff):
            mc_eff_sigma = math.sqrt(
                max(mc_eff * (1.0 - mc_eff), 0.0) / mc_den
            )
        else:
            mc_eff_sigma = float("nan")
        #endif

        ratio = (
            data_eff / mc_eff
            if np.isfinite(data_eff)
            and np.isfinite(mc_eff)
            and mc_eff > 0.0
            else float("nan")
        )
        if (
            np.isfinite(ratio)
            and data_eff > 0.0
            and mc_eff > 0.0
            and np.isfinite(data_eff_sigma)
            and np.isfinite(mc_eff_sigma)
        ):
            ratio_sigma = ratio * math.sqrt(
                (data_eff_sigma / data_eff) ** 2
                + (mc_eff_sigma / mc_eff) ** 2
            )
        else:
            ratio_sigma = float("nan")
        #endif

        rows.append(
            {
                **row,
                "data_expected_sigma": data_den_sigma,
                "data_efficiency": data_eff,
                "data_efficiency_sigma_preliminary": data_eff_sigma,
                "mc_efficiency": mc_eff,
                "mc_efficiency_sigma": mc_eff_sigma,
                "data_mc_ratio": ratio,
                "data_mc_ratio_sigma_preliminary": ratio_sigma,
            }
        )
    #endfor

    order_map = {key: i for i, (key, _a, _b) in enumerate(E_PROBE_BINS)}
    rows.sort(
        key=lambda r: (
            order_map[r["energy_bin"]],
            efficiency_bin_order(r["energy_min"]).index(r["probe_bin"]),
        )
    )
    return rows



def make_selected_event_overlap_canvas(
    period: Period,
    data_pass: Dict[str, object],
    output_dir: Path,
) -> Path:
    """
    Diagnose selected epgamma -> eppi0 event overlap BEFORE parent/probe cuts.

    Each panel corresponds to one E_probe bin and one reconstructed TAG region.
    The bar is the fraction of selected UNIQUE physical events that exist in
    the eppi0 file.  The fitted f_pi0 is overlaid only as context; it is NOT an
    expected equality because eppi0 additionally requires reconstruction of the
    partner photon.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    summary = data_pass.get("selected_overlap_summary", [])

    lookup = {
        (row["energy_bin"], row["tag_region"]): row
        for row in summary
    }

    fig, axes = plt.subplots(
        len(E_PROBE_BINS),
        2,
        figsize=(12.5, 3.1 * len(E_PROBE_BINS) + 1.0),
        squeeze=False,
    )

    for irow, (bin_key, emin, emax) in enumerate(E_PROBE_BINS):
        for icol, tag_region in enumerate(("FT", "FD")):
            ax = axes[irow, icol]
            row = lookup.get((bin_key, tag_region))

            if row is None or row["selected_unique_events"] <= 0:
                ax.text(
                    0.5,
                    0.5,
                    "No selected events",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                )
                ax.set_title(
                    f"tag {tag_region}: {emin:g} <= E_probe < {emax:g} GeV"
                )
                ax.set_ylim(0.0, 1.05)
                continue
            #endif

            overlap = float(row["overlap_unique_fraction"])
            f_pi0 = float(row["f_pi0_fit"])

            ax.bar(
                [0],
                [overlap],
                width=0.55,
                label="selected unique events also in eppi0",
            )
            if np.isfinite(f_pi0):
                ax.axhline(
                    f_pi0,
                    linestyle="--",
                    linewidth=1.2,
                    label=r"fitted $f_{\pi^0}$ (context only)",
                )
            #endif

            ax.set_xlim(-0.75, 0.75)
            ax.set_ylim(0.0, 1.05)
            ax.set_xticks([0])
            ax.set_xticklabels(["event overlap"])
            ax.set_ylabel("fraction")
            ax.grid(axis="y", alpha=0.18)
            ax.set_title(
                f"tag {tag_region}: {emin:g} <= E_probe < {emax:g} GeV",
                fontsize=10,
            )

            annotation = (
                f"selected rows = {row['selected_rows']:,}\n"
                f"unique events = {row['selected_unique_events']:,}\n"
                f"unique in eppi0 = {row['overlap_unique_events']:,}\n"
                f"overlap = {100.0*overlap:.2f}%\n"
                f"rows/event = {row['mean_rows_per_unique_event']:.3f}\n"
                f"events with >1 row = "
                f"{100.0*row['repeated_event_fraction']:.2f}%\n"
                f"max rows/event = {row['max_rows_per_event']}"
            )
            ax.text(
                0.98,
                0.96,
                annotation,
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=8.1,
                bbox={
                    "facecolor": "white",
                    "alpha": 0.88,
                    "edgecolor": "0.6",
                },
            )

            if irow == 0 and icol == 0:
                ax.legend(
                    loc="lower left",
                    frameon=True,
                    fontsize=8,
                )
            #endif
        #endfor
    #endfor

    fig.suptitle(
        (
            f"{period.label}: selected epgamma event overlap with eppi0\n"
            "exact (runnum,evnum), BEFORE e/p-parent, pi0, or probe-matching cuts"
        ),
        fontsize=13,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.955))

    out = output_dir / f"selected_event_overlap_{period.key}.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def make_matching_resolution_canvas(
    period: Period,
    mc_pass: Dict[str, object],
    data_pass: Dict[str, object],
    models: Dict[str, EfficiencyResolution],
    nsigma: float,
    output_dir: Path,
) -> Path:
    """AAO-derived Delta-p matching windows with data overlaid diagnostically."""
    output_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(
        2,
        3,
        figsize=(14.5, 8.2),
        squeeze=False,
    )

    labels = (r"$\Delta p_x$ (GeV)", r"$\Delta p_y$ (GeV)", r"$\Delta p_z$ (GeV)")
    for irow, region in enumerate(("FT", "FD")):
        mc = np.asarray(mc_pass["residuals"][region], dtype=float)
        data = np.asarray(data_pass["residuals"][region], dtype=float)
        model = models[region]

        for icomp in range(3):
            ax = axes[irow, icomp]
            if len(mc):
                lo, hi = np.nanpercentile(mc[:, icomp], [0.5, 99.5])
                span = max(hi - lo, 0.02)
                lo -= 0.1 * span
                hi += 0.1 * span
                ax.hist(
                    mc[:, icomp],
                    bins=80,
                    range=(lo, hi),
                    histtype="step",
                    density=True,
                    linewidth=1.5,
                    label="AAO clean association",
                )
                if len(data):
                    ax.hist(
                        data[:, icomp],
                        bins=80,
                        range=(lo, hi),
                        histtype="step",
                        density=True,
                        linewidth=1.2,
                        label="data clean association",
                    )
                #endif
                center = model.center[icomp]
                sigma = model.sigma[icomp]
                ax.axvline(center, linestyle="--", linewidth=1.0)
                ax.axvline(
                    center - nsigma * sigma,
                    linestyle=":",
                    linewidth=1.0,
                )
                ax.axvline(
                    center + nsigma * sigma,
                    linestyle=":",
                    linewidth=1.0,
                )
            #endif
            ax.set_xlabel(labels[icomp])
            ax.set_ylabel("density")
            ax.grid(alpha=0.18)
            ax.set_title(
                (
                    f"{region}: center={model.center[icomp]:+.4f}, "
                    f"sigma={model.sigma[icomp]:.4f} GeV\n"
                    f"N_AAO={model.n:,}; {model.source}"
                ),
                fontsize=9.0,
            )
        #endfor
    #endfor

    handles, labels_legend = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(
            handles,
            labels_legend,
            loc="upper center",
            ncol=2,
            frameon=True,
        )
    #endif
    fig.suptitle(
        (
            f"{period.label}: reconstructed-probe momentum matching\n"
            f"final window = MC center +/- {nsigma:g} robust sigma"
        ),
        fontsize=13,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.93))
    out = output_dir / f"matching_resolution_{period.key}.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def make_association_cutflow_canvas(
    period: Period,
    mc_pass: Dict[str, object],
    mc_final: Dict[str, int],
    data_pass: Dict[str, object],
    output_dir: Path,
) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)

    mc = mc_pass["counters"]
    da = data_pass["counters"]

    stages = (
        "selected tags",
        "parent/event",
        "positive remainder",
        r"$M_{\gamma\gamma}$",
        "massless remainder",
        r"$E_{\rm probe}^{reco}>0.4$",
        r"$\Delta\vec p$ match",
    )
    mc_values = np.asarray(
        [
            mc["selected_tags"],
            mc["parent_matched"],
            mc["positive"],
            mc["mass_window"],
            mc["mass_shell"],
            mc["threshold"],
            mc_final["probe_matched"],
        ],
        dtype=float,
    )
    data_values = np.asarray(
        [
            da["selected_tags"],
            da["parent_consistent"],
            da["positive"],
            da["mass_window"],
            da["mass_shell"],
            da["threshold"],
            da["probe_matched"],
        ],
        dtype=float,
    )

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.0))
    x = np.arange(len(stages))

    for ax, vals, title in (
        (axes[0], data_values, "data"),
        (axes[1], mc_values, "AAO MC"),
    ):
        denom = vals[0] if vals[0] > 0 else 1.0
        ax.bar(x, vals / denom)
        ax.set_xticks(x)
        ax.set_xticklabels(stages, rotation=32, ha="right")
        ax.set_ylim(0.0, 1.05)
        ax.set_ylabel("fraction of selected tags")
        ax.set_title(title)
        ax.grid(axis="y", alpha=0.2)
        for xi, value in zip(x, vals):
            ax.text(
                xi,
                min(value / denom + 0.025, 1.02),
                f"{int(value):,}",
                ha="center",
                va="bottom",
                fontsize=7.5,
                rotation=90,
            )
        #endfor
    #endfor

    fig.suptitle(
        (
            f"{period.label}: tag-probe association cut flow\n"
            "data uses exact run/event + e/p parent; MC uses e/p cKDTree only"
        ),
        fontsize=12.5,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.93))
    out = output_dir / f"association_cutflow_{period.key}.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def make_efficiency_counts_canvas(
    period: Period,
    rows: Sequence[dict],
    output_dir: Path,
) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(
        len(E_PROBE_BINS),
        2,
        figsize=(13.0, 3.3 * len(E_PROBE_BINS) + 0.8),
        squeeze=False,
    )

    for irow, (bin_key, emin, emax) in enumerate(E_PROBE_BINS):
        subset = [r for r in rows if r["energy_bin"] == bin_key]
        labels = [r["probe_bin"] for r in subset]
        x = np.arange(len(labels))

        ax = axes[irow, 0]
        ax.bar(x - 0.18, [r["data_expected_pi0"] for r in subset], width=0.36, label="expected pi0 tags")
        ax.bar(x + 0.18, [r["data_matched"] for r in subset], width=0.36, label="matched probes")
        ax.set_title(f"Data: {emin:g} <= E_probe < {emax:g} GeV")
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.set_ylabel("events")
        ax.grid(axis="y", alpha=0.18)

        ax = axes[irow, 1]
        ax.bar(x - 0.18, [r["mc_denominator"] for r in subset], width=0.36, label="AAO tags")
        ax.bar(x + 0.18, [r["mc_matched"] for r in subset], width=0.36, label="matched probes")
        ax.set_title(f"AAO MC: {emin:g} <= E_probe < {emax:g} GeV")
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.set_ylabel("events")
        ax.grid(axis="y", alpha=0.18)

        if irow == 0:
            axes[irow, 0].legend(frameon=True)
            axes[irow, 1].legend(frameon=True)
        #endif
    #endfor

    fig.suptitle(
        f"{period.label}: efficiency numerator/denominator statistics",
        fontsize=13,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.965))
    out = output_dir / f"efficiency_counts_{period.key}.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def make_efficiency_summary_canvas(
    period: Period,
    rows: Sequence[dict],
    output_dir: Path,
) -> Path:
    """
    Four rows (E_probe bins) x three columns:
      epsilon_data, epsilon_MC, epsilon_data/epsilon_MC.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(
        len(E_PROBE_BINS),
        3,
        figsize=(14.5, 3.2 * len(E_PROBE_BINS) + 0.8),
        squeeze=False,
    )

    for irow, (bin_key, emin, emax) in enumerate(E_PROBE_BINS):
        subset = [r for r in rows if r["energy_bin"] == bin_key]
        labels = [r["probe_bin"] for r in subset]
        x = np.arange(len(labels))

        panels = (
            (
                "data_efficiency",
                "data_efficiency_sigma_preliminary",
                r"$\epsilon_{\rm data}$",
                (0.0, 1.2),
            ),
            (
                "mc_efficiency",
                "mc_efficiency_sigma",
                r"$\epsilon_{\rm MC}$",
                (0.0, 1.05),
            ),
            (
                "data_mc_ratio",
                "data_mc_ratio_sigma_preliminary",
                r"$\epsilon_{\rm data}/\epsilon_{\rm MC}$",
                (0.5, 1.5),
            ),
        )

        for icol, (ykey, ekey, ylabel, ylim) in enumerate(panels):
            ax = axes[irow, icol]
            y = np.asarray([r[ykey] for r in subset], dtype=float)
            yerr = np.asarray([r[ekey] for r in subset], dtype=float)
            finite = np.isfinite(y)
            ax.errorbar(
                x[finite],
                y[finite],
                yerr=np.where(
                    np.isfinite(yerr[finite]),
                    yerr[finite],
                    0.0,
                ),
                fmt="o",
                capsize=2,
            )
            if icol == 2:
                ax.axhline(1.0, linestyle="--", linewidth=1.0)
            #endif
            ax.set_xticks(x)
            ax.set_xticklabels(labels)
            ax.set_ylim(*ylim)
            ax.set_ylabel(ylabel)
            ax.grid(alpha=0.2)
            ax.set_title(
                f"{emin:g} <= E_probe < {emax:g} GeV"
                if icol == 1 else ""
            )
        #endfor
    #endfor

    fig.suptitle(
        (
            f"{period.label}: preliminary photon reconstruction efficiency\n"
            "FD sectors separate below 4 GeV; 4-9.5 GeV FD combined"
        ),
        fontsize=13,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.955))
    out = output_dir / f"efficiency_summary_{period.key}.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def write_efficiency_summary_csv(
    period: Period,
    rows: Sequence[dict],
    mc_pass: Dict[str, object],
    data_pass: Dict[str, object],
    models: Dict[str, EfficiencyResolution],
    nsigma: float,
    output_dir: Path,
) -> Path:
    """Exactly one compact machine-readable efficiency table per period."""
    output_dir.mkdir(parents=True, exist_ok=True)
    out = output_dir / f"efficiency_summary_{period.key}.csv"

    fields = (
        "period",
        "energy_bin",
        "energy_min",
        "energy_max",
        "probe_bin",
        "data_raw_tags",
        "data_tagFT_raw",
        "data_tagFD_raw",
        "data_expected_pi0",
        "data_expected_sigma",
        "data_matched",
        "data_efficiency",
        "data_efficiency_sigma_preliminary",
        "mc_denominator",
        "mc_tagFT_raw",
        "mc_tagFD_raw",
        "mc_matched",
        "mc_efficiency",
        "mc_efficiency_sigma",
        "data_mc_ratio",
        "data_mc_ratio_sigma_preliminary",
        "match_nsigma",
        "ft_dp_center_x",
        "ft_dp_center_y",
        "ft_dp_center_z",
        "ft_dp_sigma_x",
        "ft_dp_sigma_y",
        "ft_dp_sigma_z",
        "fd_dp_center_x",
        "fd_dp_center_y",
        "fd_dp_center_z",
        "fd_dp_sigma_x",
        "fd_dp_sigma_y",
        "fd_dp_sigma_z",
    )

    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    **{key: row.get(key, "") for key in fields},
                    "period": period.key,
                    "match_nsigma": nsigma,
                    "ft_dp_center_x": models["FT"].center[0],
                    "ft_dp_center_y": models["FT"].center[1],
                    "ft_dp_center_z": models["FT"].center[2],
                    "ft_dp_sigma_x": models["FT"].sigma[0],
                    "ft_dp_sigma_y": models["FT"].sigma[1],
                    "ft_dp_sigma_z": models["FT"].sigma[2],
                    "fd_dp_center_x": models["FD"].center[0],
                    "fd_dp_center_y": models["FD"].center[1],
                    "fd_dp_center_z": models["FD"].center[2],
                    "fd_dp_sigma_x": models["FD"].sigma[0],
                    "fd_dp_sigma_y": models["FD"].sigma[1],
                    "fd_dp_sigma_z": models["FD"].sigma[2],
                }
            )
        #endfor
    #endwith
    return out


def run_efficiency_extraction(
    period: Period,
    optimizer_summary: dict,
    tree_name: str,
    max_entries: int,
    eppi0_max_entries: int,
    step_size: int,
    parent_component_tol: float,
    parent_distance_max: float,
    kdtree_workers: int,
    assoc_mgg_min: float,
    assoc_mgg_max: float,
    assoc_remainder_mass2_max: float,
    assoc_probe_energy_min: float,
    probe_match_nsigma: float,
    min_resolution_events: int,
    coverage_min_eppi0_overlap: float,
    output_root: Path,
) -> dict:
    """End-to-end efficiency stage for one period."""
    t0 = time.perf_counter()
    outdir = output_root / "efficiency" / period.key
    outdir.mkdir(parents=True, exist_ok=True)

    # First establish that the two DATA skims really cover the same event
    # population. This is a FULL-tree set intersection, independent of
    # --max-entries, row ordering, and all physics cuts.
    coverage = run_raw_coverage_audit(
        period,
        tree_name,
        step_size,
        output_root,
        min_eppi0_overlap=coverage_min_eppi0_overlap,
        fail_on_pathology=True,
    )

    # Preflight first so a missing companion file fails before expensive work.
    for label, path, req in (
        (
            "data eppi0",
            period.eppi0_data,
            EFF_EPPIO_REQUIRED + EFF_EPPIO_DATA_IDS,
        ),
        (
            "AAO eppi0",
            period.eppi0_pi0_mc,
            EFF_EPPIO_REQUIRED,
        ),
    ):
        found, total = preflight_custom_branches(
            path,
            tree_name,
            req,
        )
        log(
            f"{period.label}: efficiency preflight {label}: "
            f"{Path(path).name}, tree '{found}', {total:,} entries."
        )
    #endfor

    # Full eppi0 stores are recommended even in a limited epgamma concept run;
    # otherwise the two independently truncated files need not overlap.
    data_eppi0 = load_eppi0_efficiency_store(
        period.eppi0_data,
        tree_name,
        eppi0_max_entries,
        step_size,
        require_event_ids=True,
        label=f"{period.label} data",
    )
    mc_eppi0 = load_eppi0_efficiency_store(
        period.eppi0_pi0_mc,
        tree_name,
        eppi0_max_entries,
        step_size,
        require_event_ids=False,
        label=f"{period.label} AAO",
    )

    accumulator = initialize_efficiency_accumulator()

    log(
        f"{period.label}: efficiency MC pass — runnum/evnum intentionally ignored; "
        "matching AAO epgamma<->eppi0 by reconstructed e/p only."
    )
    mc_pass = process_mc_efficiency_pass(
        period,
        optimizer_summary,
        tree_name,
        max_entries,
        step_size,
        mc_eppi0,
        parent_component_tol,
        parent_distance_max,
        kdtree_workers,
        assoc_mgg_min,
        assoc_mgg_max,
        assoc_remainder_mass2_max,
        assoc_probe_energy_min,
        accumulator,
    )

    models = derive_efficiency_resolutions(
        mc_pass,
        min_resolution_events,
    )
    mc_final = finalize_mc_numerator(
        accumulator,
        mc_pass,
        models,
        probe_match_nsigma,
    )

    log(
        f"{period.label}: efficiency data pass — exact (runnum,evnum) + "
        "reconstructed e/p parent consistency."
    )
    data_pass = process_data_efficiency_pass(
        period,
        optimizer_summary,
        tree_name,
        max_entries,
        step_size,
        data_eppi0,
        parent_component_tol,
        assoc_mgg_min,
        assoc_mgg_max,
        assoc_remainder_mass2_max,
        assoc_probe_energy_min,
        models,
        probe_match_nsigma,
        accumulator,
    )

    rows = calculate_efficiency_results(accumulator)

    selected_overlap_plot = make_selected_event_overlap_canvas(
        period,
        data_pass,
        outdir,
    )
    matching_plot = make_matching_resolution_canvas(
        period,
        mc_pass,
        data_pass,
        models,
        probe_match_nsigma,
        outdir,
    )
    cutflow_plot = make_association_cutflow_canvas(
        period,
        mc_pass,
        mc_final,
        data_pass,
        outdir,
    )
    counts_plot = make_efficiency_counts_canvas(
        period,
        rows,
        outdir,
    )
    summary_plot = make_efficiency_summary_canvas(
        period,
        rows,
        outdir,
    )
    summary_csv = write_efficiency_summary_csv(
        period,
        rows,
        mc_pass,
        data_pass,
        models,
        probe_match_nsigma,
        outdir,
    )

    log(
        f"{period.label}: efficiency outputs -> {outdir}; "
        f"plots: {selected_overlap_plot.name}, {matching_plot.name}, "
        f"{cutflow_plot.name}, {counts_plot.name}, {summary_plot.name}; "
        f"one CSV: {summary_csv.name}."
    )
    for overlap_row in data_pass.get("selected_overlap_summary", []):
        log(
            f"{period.label}: SELECTED EVENT OVERLAP "
            f"{overlap_row['energy_bin']} tag-{overlap_row['tag_region']}: "
            f"{overlap_row['overlap_unique_events']:,}/"
            f"{overlap_row['selected_unique_events']:,} unique events "
            f"({100.0*overlap_row['overlap_unique_fraction']:.2f}%) are in "
            f"eppi0; {overlap_row['selected_rows']:,} selected rows, "
            f"{overlap_row['mean_rows_per_unique_event']:.3f} rows/event, "
            f"f_pi0_fit={overlap_row['f_pi0_fit']:.3f}."
        )
    #endfor

    log(
        f"{period.label}: directed-tag audit: "
        f"data selected-event multiplicity fraction="
        f"{data_pass['directed_tag_event_multiplicity_fraction']:.3f}; "
        f"AAO repeated-eppi0-parent fraction="
        f"{mc_pass['directed_tag_duplicate_fraction']:.3f}. "
        "These are diagnostics; directed tags are retained rather than silently "
        "deduplicated."
    )

    return {
        "rows": rows,
        "data_counters": data_pass["counters"],
        "selected_overlap_summary": data_pass.get(
            "selected_overlap_summary",
            [],
        ),
        "mc_counters": mc_pass["counters"],
        "mc_final": mc_final,
        "resolution_ft": {
            "center": models["FT"].center.tolist(),
            "sigma": models["FT"].sigma.tolist(),
            "n": models["FT"].n,
            "source": models["FT"].source,
        },
        "resolution_fd": {
            "center": models["FD"].center.tolist(),
            "sigma": models["FD"].sigma.tolist(),
            "n": models["FD"].n,
            "source": models["FD"].source,
        },
        "coverage": {
            "n_epgamma_rows_scanned": coverage["n_epgamma_rows_scanned"],
            "n_epgamma_rows_total": coverage["n_epgamma_rows_total"],
            "n_eppi0_unique": coverage["n_eppi0_unique"],
            "n_common_unique": coverage["n_common_unique"],
            "common_over_eppi0": coverage["common_over_eppi0"],
            "n_epgamma_runs_seen": coverage["n_epgamma_runs_seen"],
            "n_eppi0_runs": coverage["n_eppi0_runs"],
            "n_common_runs": coverage["n_common_runs"],
            "stopped_early": coverage["stopped_early"],
        },
        "elapsed_s": float(time.perf_counter() - t0),
    }


# =============================================================================
# Period driver
# =============================================================================

def process_period(
    period: Period,
    tree_name: str,
    max_entries: int,
    step_size: int,
    output_dir: Path,
    linear_y: bool,
    run_optimizer: bool,
    optimizer_max_events: int,
    optimizer_steps: int,
    optimizer_quantiles: int,
    optimizer_min_step_pi0_eff: float,
    optimizer_min_total_pi0_eff: float,
    optimizer_min_improvement: float,
    optimizer_mx2_1_preselection_max: float,
    optimizer_data_sideband_min: float,
    optimizer_data_sideband_max: float,
    optimizer_data_sideband_weight: float,
    closure_toys: int,
    closure_events: int,
    closure_bins: int,
    data_fit_max_closure_bias: float,
    data_fit_max_closure_rms: float,
    data_fit_background_near_best: float,
    data_fit_morph_shift_max_bins: float,
    data_fit_morph_smear_max_bins: float,
    data_fit_morph_shift_steps: int,
    data_fit_morph_smear_steps: int,
    data_fit_sideband_min_events: int,
    data_fit_sideband_max_js: float,
    run_efficiency: bool,
    efficiency_eppi0_max_entries: int,
    parent_component_tol: float,
    parent_distance_max: float,
    kdtree_workers: int,
    assoc_mgg_min: float,
    assoc_mgg_max: float,
    assoc_remainder_mass2_max: float,
    assoc_probe_energy_min: float,
    probe_match_nsigma: float,
    probe_match_min_resolution_events: int,
    coverage_min_eppi0_overlap: float,
) -> dict:
    t0 = time.perf_counter()
    log(f"{period.label}: starting.")

    optimizer_summary = {
        "period_key": period.key,
        "period_label": period.label,
        "regions": {},
    }

    paths = {
        "data": period.data,
        "pi0": period.pi0_mc,
        "dvcs": period.dvcs_mc,
    }

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
    optimizer_events: Dict[str, Dict[str, np.ndarray]] = {}

    period_index = [p.key for p in PERIODS].index(period.key)

    for sample_index, (sample_key, label, _color) in enumerate(SAMPLES):
        # In optimizer mode, retain a bounded uniform event reservoir for
        # all three samples. AAO and DVCSgen determine the cuts; the data
        # reservoir is used only afterward to visualize those same cuts.
        collect_optimizer = run_optimizer

        hists, counts, opt_events = stream_sample_histograms(
            paths[sample_key],
            f"{period.label} {label}",
            tree_name,
            max_entries,
            step_size,
            BEAM_ENERGY_GEV[period.key],
            collect_optimizer=collect_optimizer,
            optimizer_max_events=optimizer_max_events,
            optimizer_seed=(
                100_000
                + 10_000 * period_index
                + 100 * sample_index
            ),
        )
        sample_hists[sample_key] = hists
        sample_counts[sample_key] = counts

        if collect_optimizer:
            optimizer_events[sample_key] = opt_events
        #endif
    #endfor

    if run_optimizer:
        optimizer_outdir = output_dir.parent / "rectangular_optimizer"

        optimizer_summary["energy_bins"] = {}

        for bin_key, energy_min, energy_max in E_PROBE_BINS:
            optimizer_summary["energy_bins"][bin_key] = {
                "energy_min": float(energy_min),
                "energy_max": float(energy_max),
                "regions": {},
            }

            for region in ("FT", "FD"):
                pi0_all = optimizer_events["pi0"][region][bin_key]
                dvcs_all = optimizer_events["dvcs"][region][bin_key]
                data_all = optimizer_events["data"][region][bin_key]

                mx2_idx = FEATURE_INDEX["Mx2_1"]

                pi0_pre = (
                    np.isfinite(pi0_all[:, mx2_idx])
                    & (
                        pi0_all[:, mx2_idx]
                        < optimizer_mx2_1_preselection_max
                    )
                )
                dvcs_pre = (
                    np.isfinite(dvcs_all[:, mx2_idx])
                    & (
                        dvcs_all[:, mx2_idx]
                        < optimizer_mx2_1_preselection_max
                    )
                )
                sideband_mask = (
                    np.isfinite(data_all[:, mx2_idx])
                    & (
                        data_all[:, mx2_idx]
                        >= optimizer_data_sideband_min
                    )
                    & (
                        data_all[:, mx2_idx]
                        < optimizer_data_sideband_max
                    )
                )

                pi0_events = pi0_all[pi0_pre]
                dvcs_events = dvcs_all[dvcs_pre]
                data_sideband_events = data_all[sideband_mask]

                log(
                    f"{period.label} {region}, "
                    f"{energy_min:g}<=E_probe<{energy_max:g} GeV: "
                    f"after Mx2_1<"
                    f"{optimizer_mx2_1_preselection_max:.3g} GeV^2 -> "
                    f"AAO={pi0_events.shape[0]:,}, "
                    f"DVCSgen={dvcs_events.shape[0]:,}; "
                    f"data sideband "
                    f"{optimizer_data_sideband_min:.3g}<=Mx2_1<"
                    f"{optimizer_data_sideband_max:.3g} GeV^2 -> "
                    f"{data_sideband_events.shape[0]:,}."
                )

                results = optimize_rectangular_cuts(
                    pi0_events,
                    dvcs_events,
                    data_sideband_events,
                    max_steps=optimizer_steps,
                    n_quantiles=optimizer_quantiles,
                    min_step_pi0_eff=optimizer_min_step_pi0_eff,
                    min_total_pi0_eff=optimizer_min_total_pi0_eff,
                    min_improvement=optimizer_min_improvement,
                    data_sideband_weight=optimizer_data_sideband_weight,
                )

                closure_results = run_template_fit_closure(
                    pi0_events,
                    dvcs_events,
                    results,
                    n_toys=closure_toys,
                    n_pseudodata_events=closure_events,
                    n_bins=closure_bins,
                    seed=(
                        700_000
                        + 10_000 * period_index
                        + 1000 * [key for key, _lo, _hi in E_PROBE_BINS].index(
                            bin_key
                        )
                        + (0 if region == "FT" else 100)
                    ),
                )

                optimizer_summary["energy_bins"][bin_key]["regions"][region] = {
                    "n_pi0_baseline": int(pi0_events.shape[0]),
                    "n_dvcs_baseline": int(dvcs_events.shape[0]),
                    "n_data_sideband_baseline": int(
                        data_sideband_events.shape[0]
                    ),
                    "results": results,
                    "closure_results": closure_results,
                }
            #endfor
        #endfor

        # --------------------------------------------------------------
        # Real-data AAO + DVCSgen composition fit.
        # --------------------------------------------------------------
        for bin_key, energy_min, energy_max in E_PROBE_BINS:
            for region in ("FT", "FD"):
                region_summary = (
                    optimizer_summary["energy_bins"][bin_key]
                    ["regions"][region]
                )
                optimizer_results = region_summary.get("results", [])
                closure_results = region_summary.get(
                    "closure_results",
                    [],
                )

                operating = choose_data_fit_operating_step(
                    optimizer_results,
                    closure_results,
                    sideband_weight=optimizer_data_sideband_weight,
                    max_abs_bias=data_fit_max_closure_bias,
                    max_rms=data_fit_max_closure_rms,
                    near_best_factor=data_fit_background_near_best,
                )

                mx2_idx = FEATURE_INDEX["Mx2_1"]

                pi0_all = optimizer_events["pi0"][region][bin_key]
                dvcs_all = optimizer_events["dvcs"][region][bin_key]
                data_all = optimizer_events["data"][region][bin_key]

                # The real-data fit target and empirical third template are
                # disjoint and eta-safe by construction.
                signal_max = SIGNAL_MX2_1_MAX

                pi0_events = pi0_all[
                    np.isfinite(pi0_all[:, mx2_idx])
                    & (pi0_all[:, mx2_idx] < signal_max)
                ]
                dvcs_events = dvcs_all[
                    np.isfinite(dvcs_all[:, mx2_idx])
                    & (dvcs_all[:, mx2_idx] < signal_max)
                ]
                data_events = data_all[
                    np.isfinite(data_all[:, mx2_idx])
                    & (data_all[:, mx2_idx] < signal_max)
                ]

                # ONLY the closest eta-safe sideband is transported into the
                # signal region as the empirical third template.
                sideband_events = data_all[
                    np.isfinite(data_all[:, mx2_idx])
                    & (
                        data_all[:, mx2_idx]
                        >= NEAR_SIDEBAND_MIN
                    )
                    & (
                        data_all[:, mx2_idx]
                        < NEAR_SIDEBAND_MAX
                    )
                ]

                fit_feature = choose_real_data_fit_feature(
                    pi0_events,
                    dvcs_events,
                    optimizer_results,
                    chosen_step=int(operating["step"]),
                    closure_results=closure_results,
                    n_bins=closure_bins,
                )

                if fit_feature is None:
                    sideband_stability = {
                        "success": False,
                        "status": "no_valid_fit_feature",
                    }
                    use_sideband = False
                    sideband_decision = "no_valid_fit_feature"
                else:
                    sideband_stability = build_sideband_stability_result(
                        data_all,
                        optimizer_results,
                        chosen_step=int(operating["step"]),
                        feature=str(fit_feature),
                        n_bins=closure_bins,
                    )
                    use_sideband, sideband_decision = (
                        decide_sideband_template_usage(
                            sideband_stability,
                            min_events=data_fit_sideband_min_events,
                            max_js=data_fit_sideband_max_js,
                        )
                    )
                #endif

                if use_sideband:
                    data_fit = fit_real_data_three_templates(
                        data_events,
                        pi0_events,
                        dvcs_events,
                        sideband_events,
                        optimizer_results,
                        chosen_step=int(operating["step"]),
                        closure_results=closure_results,
                        n_bins=closure_bins,
                        forced_feature=fit_feature,
                        morph_shift_max_bins=data_fit_morph_shift_max_bins,
                        morph_smear_max_bins=data_fit_morph_smear_max_bins,
                        morph_shift_steps=data_fit_morph_shift_steps,
                        morph_smear_steps=data_fit_morph_smear_steps,
                    )

                    # If the accepted sideband still fails technically inside
                    # the three-template fit, fall back rather than lose the
                    # bin entirely.
                    if not data_fit.get("success", False):
                        fallback_reason = (
                            "three_template_failed_"
                            + data_fit.get("status", "unknown")
                        )
                        data_fit = fit_real_data_two_templates_morphed(
                            data_events,
                            pi0_events,
                            dvcs_events,
                            optimizer_results,
                            chosen_step=int(operating["step"]),
                            feature=str(fit_feature),
                            n_bins=closure_bins,
                            morph_shift_max_bins=data_fit_morph_shift_max_bins,
                            morph_smear_max_bins=data_fit_morph_smear_max_bins,
                            morph_shift_steps=data_fit_morph_shift_steps,
                            morph_smear_steps=data_fit_morph_smear_steps,
                            fallback_reason=fallback_reason,
                        )
                    #endif
                else:
                    if fit_feature is None:
                        data_fit = {
                            "success": False,
                            "status": "no_valid_fit_feature",
                            "fit_mode": "unavailable",
                            "fallback_reason": sideband_decision,
                            "step": int(operating["step"]),
                        }
                    else:
                        data_fit = fit_real_data_two_templates_morphed(
                            data_events,
                            pi0_events,
                            dvcs_events,
                            optimizer_results,
                            chosen_step=int(operating["step"]),
                            feature=str(fit_feature),
                            n_bins=closure_bins,
                            morph_shift_max_bins=data_fit_morph_shift_max_bins,
                            morph_smear_max_bins=data_fit_morph_smear_max_bins,
                            morph_shift_steps=data_fit_morph_shift_steps,
                            morph_smear_steps=data_fit_morph_smear_steps,
                            fallback_reason=sideband_decision,
                        )
                    #endif
                #endif

                data_fit["operating_step_status"] = operating["status"]
                data_fit["closure_rms_limit"] = operating[
                    "closure_rms_limit"
                ]
                data_fit["background_metric"] = operating[
                    "background_metric"
                ]
                data_fit["sideband_decision"] = sideband_decision
                data_fit["sideband_accepted"] = bool(use_sideband)

                region_summary["data_fit"] = data_fit
                region_summary["sideband_stability"] = (
                    sideband_stability
                )
                if data_fit.get("success", False):
                    log(
                        f"{period.label} {region} "
                        f"{energy_min:g}-{energy_max:g} GeV: "
                        f"real-data fit step {data_fit['step']} "
                        f"({data_fit['feature']}), "
                        f"mode={data_fit.get('fit_mode', 'unknown')}, "
                        f"f_pi0={data_fit['fraction_pi0']:.4f}, "
                        f"f_DVCS={data_fit['fraction_dvcs']:.4f}, "
                        f"f_SB={data_fit['fraction_sideband']:.4f}, "
                        f"AAO morph=({data_fit['morph_shift_bins']:+.2f}, "
                        f"{data_fit['morph_smear_sigma_bins']:.2f}) bins, "
                        f"D/ndf={data_fit['nominal_deviance_per_ndf']:.2f}"
                        f"->{data_fit['deviance_per_ndf']:.2f}."
                    )
                else:
                    log(
                        f"{period.label} {region} "
                        f"{energy_min:g}-{energy_max:g} GeV: "
                        f"real-data fit unavailable "
                        f"({data_fit.get('status', 'unknown')})."
                    )
                #endif
            #endfor
        #endfor

        combined_optimizer_plot = make_combined_optimizer_canvas(
            period,
            optimizer_summary,
            optimizer_outdir,
        )
        log(
            f"{period.label}: wrote E_probe-binned FD+FT optimizer canvas "
            f"{combined_optimizer_plot}"
        )

        combined_closure_plot = make_combined_closure_canvas(
            period,
            optimizer_summary,
            output_dir.parent / "template_fit_closure",
        )
        log(
            f"{period.label}: wrote combined template-fit closure canvas "
            f"{combined_closure_plot}"
        )

        combined_response_plot = (
            make_combined_closure_response_canvas(
                period,
                optimizer_summary,
                output_dir.parent / "template_fit_closure",
            )
        )
        log(
            f"{period.label}: wrote fitted-vs-injected closure response "
            f"canvas {combined_response_plot}"
        )

        data_fit_outdir = output_dir.parent / "data_template_fit"

        data_fit_canvas = make_combined_data_template_fit_canvas(
            period,
            optimizer_summary,
            data_fit_outdir,
        )
        log(
            f"{period.label}: wrote combined real-data template-fit canvas "
            f"{data_fit_canvas}"
        )

        fraction_summary_canvas = make_pi0_fraction_summary_canvas(
            period,
            optimizer_summary,
            data_fit_outdir,
        )
        log(
            f"{period.label}: wrote composition-fraction summary canvas "
            f"{fraction_summary_canvas}"
        )

        sideband_stability_canvas = make_sideband_stability_canvas(
            period,
            optimizer_summary,
            data_fit_outdir,
        )
        log(
            f"{period.label}: wrote eta-safe sideband-stability canvas "
            f"{sideband_stability_canvas}"
        )

        mx2_control_canvas = make_mx2_1_control_region_canvas(
            period,
            optimizer_summary,
            optimizer_events,
            data_fit_outdir,
        )
        log(
            f"{period.label}: wrote Mx2_1 pi0/sideband/eta control canvas "
            f"{mx2_control_canvas}"
        )

        # TAG-FD / predicted-probe-sector diagnostics: deliberately keep the period-level FD
        # composition fit fixed and only expose sector-dependent shapes.
        sector_diagnostics = build_fd_sector_diagnostics(
            optimizer_events,
            optimizer_summary,
            n_bins=closure_bins,
        )
        sector_outdir = (
            output_dir.parent
            / "fd_sector_diagnostics"
            / period.key
        )

        sector_fit_canvases = make_fd_sector_template_fit_canvases(
            period,
            sector_diagnostics,
            sector_outdir,
        )
        log(
            f"{period.label}: wrote {len(sector_fit_canvases)} FD sector "
            f"canvases under {sector_outdir}."
        )

        sector_overview = make_fd_sector_overview_canvas(
            period,
            sector_diagnostics,
            sector_outdir,
        )
        log(f"{period.label}: wrote FD sector overview {sector_overview}.")

        sector_csv = write_fd_sector_summary_csv(
            period,
            sector_diagnostics,
            sector_outdir,
        )
        log(f"{period.label}: wrote FD sector summary {sector_csv}.")

        summary_csv = write_period_summary_csv(
            period,
            optimizer_summary,
            output_dir.parent / "summary",
        )
        log(
            f"{period.label}: wrote single combined optimizer/closure/data-fit CSV "
            f"{summary_csv}"
        )

        # Make one cumulative shape canvas for every period x detector x
        # E_probe bin, using exactly the cuts selected in that bin.
        for bin_key, energy_min, energy_max in E_PROBE_BINS:
            for region in ("FT", "FD"):
                region_results = (
                    optimizer_summary["energy_bins"][bin_key]
                    ["regions"][region]["results"]
                )
                shape_out = make_optimizer_cut_shape_canvas(
                    period,
                    region,
                    bin_key,
                    energy_min,
                    energy_max,
                    optimizer_events,
                    region_results,
                    output_dir,
                    linear_y,
                    optimizer_mx2_1_preselection_max,
                )
                log(
                    f"{period.label}: wrote {region} "
                    f"{energy_min:g}-{energy_max:g} GeV "
                    f"post-optimizer shape canvas {shape_out}"
                )
            #endfor
        #endfor
    else:
        # If optimization is explicitly disabled, preserve the old fixed
        # three-row unbinned concept canvases.
        for region in ("FT", "FD"):
            out = make_canvas(
                period,
                region,
                sample_hists,
                sample_counts,
                output_dir,
                linear_y,
            )
            log(f"{period.label}: wrote {out}")
        #endfor
    #endif

    if run_efficiency and run_optimizer:
        efficiency_summary = run_efficiency_extraction(
            period,
            optimizer_summary,
            tree_name,
            max_entries,
            efficiency_eppi0_max_entries,
            step_size,
            parent_component_tol,
            parent_distance_max,
            kdtree_workers,
            assoc_mgg_min,
            assoc_mgg_max,
            assoc_remainder_mass2_max,
            assoc_probe_energy_min,
            probe_match_nsigma,
            probe_match_min_resolution_events,
            coverage_min_eppi0_overlap,
            output_dir.parent,
        )
        optimizer_summary["efficiency"] = efficiency_summary
    elif run_efficiency and not run_optimizer:
        log(
            f"{period.label}: efficiency extraction skipped because the "
            "purity optimizer/data-fit stage is disabled."
        )
    #endif

    log(
        f"{period.label}: complete in "
        f"{time.perf_counter() - t0:.1f} s."
    )

    return optimizer_summary


def process_period_worker(
    period: Period,
    tree_name: str,
    max_entries: int,
    step_size: int,
    output_dir_str: str,
    linear_y: bool,
    run_optimizer: bool,
    optimizer_max_events: int,
    optimizer_steps: int,
    optimizer_quantiles: int,
    optimizer_min_step_pi0_eff: float,
    optimizer_min_total_pi0_eff: float,
    optimizer_min_improvement: float,
    optimizer_mx2_1_preselection_max: float,
    optimizer_data_sideband_min: float,
    optimizer_data_sideband_max: float,
    optimizer_data_sideband_weight: float,
    closure_toys: int,
    closure_events: int,
    closure_bins: int,
    data_fit_max_closure_bias: float,
    data_fit_max_closure_rms: float,
    data_fit_background_near_best: float,
    data_fit_morph_shift_max_bins: float,
    data_fit_morph_smear_max_bins: float,
    data_fit_morph_shift_steps: int,
    data_fit_morph_smear_steps: int,
    data_fit_sideband_min_events: int,
    data_fit_sideband_max_js: float,
    run_efficiency: bool,
    efficiency_eppi0_max_entries: int,
    parent_component_tol: float,
    parent_distance_max: float,
    kdtree_workers: int,
    assoc_mgg_min: float,
    assoc_mgg_max: float,
    assoc_remainder_mass2_max: float,
    assoc_probe_energy_min: float,
    probe_match_nsigma: float,
    probe_match_min_resolution_events: int,
    coverage_min_eppi0_overlap: float,
) -> str:
    """
    Picklable process-level wrapper for one run period.

    Each worker streams that period's data, AAO, and DVCSgen files
    sequentially. FT and FD are filled simultaneously within each file,
    avoiding nested ROOT I/O parallelism.
    """
    return process_period(
        period,
        tree_name,
        max_entries,
        step_size,
        Path(output_dir_str),
        linear_y,
        run_optimizer,
        optimizer_max_events,
        optimizer_steps,
        optimizer_quantiles,
        optimizer_min_step_pi0_eff,
        optimizer_min_total_pi0_eff,
        optimizer_min_improvement,
        optimizer_mx2_1_preselection_max,
        optimizer_data_sideband_min,
        optimizer_data_sideband_max,
        optimizer_data_sideband_weight,
        closure_toys,
        closure_events,
        closure_bins,
        data_fit_max_closure_bias,
        data_fit_max_closure_rms,
        data_fit_background_near_best,
        data_fit_morph_shift_max_bins,
        data_fit_morph_smear_max_bins,
        data_fit_morph_shift_steps,
        data_fit_morph_smear_steps,
        data_fit_sideband_min_events,
        data_fit_sideband_max_js,
        run_efficiency,
        efficiency_eppi0_max_entries,
        parent_component_tol,
        parent_distance_max,
        kdtree_workers,
        assoc_mgg_min,
        assoc_mgg_max,
        assoc_remainder_mass2_max,
        assoc_probe_energy_min,
        probe_match_nsigma,
        probe_match_min_resolution_events,
        coverage_min_eppi0_overlap,
    )



def print_optimizer_summary(
    summaries_by_period: Dict[str, dict],
    selected_periods: Sequence[Period],
) -> None:
    """
    Print all optimizer results once, serially, after parallel processing ends.
    Results are grouped by period, then E_probe bin, then detector region.
    """
    print("", flush=True)
    print("=" * 100, flush=True)
    print("E_PROBE-BINNED RECTANGULAR CUT OPTIMIZER SUMMARY", flush=True)
    print("=" * 100, flush=True)

    for period in selected_periods:
        summary = summaries_by_period.get(period.key)
        if summary is None:
            continue
        #endif

        print(f"\n{period.label}", flush=True)
        print("-" * len(period.label), flush=True)

        energy_summaries = summary.get("energy_bins", {})
        if not energy_summaries:
            print(
                "  Optimizer disabled or no optimizer result.",
                flush=True,
            )
            continue
        #endif

        for bin_key, energy_min, energy_max in E_PROBE_BINS:
            bin_summary = energy_summaries.get(bin_key)
            if bin_summary is None:
                continue
            #endif

            print(
                f"\n  {energy_min:g} <= E_probe < {energy_max:g} GeV",
                flush=True,
            )

            for region in ("FT", "FD"):
                region_summary = bin_summary.get("regions", {}).get(region)
                if region_summary is None:
                    continue
                #endif

                print(
                    f"    {region}: baseline "
                    f"AAO={region_summary['n_pi0_baseline']:,}, "
                    f"DVCSgen={region_summary['n_dvcs_baseline']:,}, "
                    f"data-SB="
                    f"{region_summary['n_data_sideband_baseline']:,}",
                    flush=True,
                )

                results = region_summary.get("results", [])
                if not results:
                    print(
                        "      No cut passed the configured "
                        "retention/improvement constraints.",
                        flush=True,
                    )
                    continue
                #endif

                for result in results:
                    print(
                        f"      {result['step']}. "
                        f"{cut_expression(result)}"
                        f"  | step AAO="
                        f"{100.0*result['step_pi0_eff']:.2f}%"
                        f", DVCS="
                        f"{100.0*result['step_dvcs_eff']:.2f}%"
                        f", data-SB="
                        f"{100.0*result['step_data_sideband_eff']:.2f}%"
                        f"  | cumulative AAO="
                        f"{100.0*result['total_pi0_eff']:.2f}%"
                        f", DVCS="
                        f"{100.0*result['total_dvcs_eff']:.2f}%"
                        f", data-SB="
                        f"{100.0*result['total_data_sideband_eff']:.2f}%",
                        flush=True,
                    )
                #endfor

                final = results[-1]
                print(
                    f"      FINAL: AAO="
                    f"{100.0*final['total_pi0_eff']:.2f}%, "
                    f"DVCSgen="
                    f"{100.0*final['total_dvcs_eff']:.2f}%, "
                    f"data-SB="
                    f"{100.0*final['total_data_sideband_eff']:.2f}%, "
                    f"F={final['score']:.4g}",
                    flush=True,
                )

                closure_results = region_summary.get(
                    "closure_results",
                    [],
                )
                final_closure = [
                    result
                    for result in closure_results
                    if int(result["step"]) == len(results)
                    and abs(result["true_fraction"] - 0.50) < 1.0e-9
                ]
                if final_closure:
                    closure = final_closure[0]
                    print(
                        f"      CLOSURE @ final step, f_pi0=0.50: "
                        f"best={closure['feature']}, "
                        f"JS={closure['js_divergence']:.4f}, "
                        f"<f_fit>={closure['mean_fit']:.4f}, "
                        f"bias={closure['bias']:+.4f}, "
                        f"RMS={closure['rms_error']:.4f}",
                        flush=True,
                    )
                #endif

                data_fit = region_summary.get("data_fit")
                if data_fit and data_fit.get("success", False):
                    print(
                        f"      DATA FIT: step {data_fit['step']}, "
                        f"{data_fit['feature']}, "
                        f"mode={data_fit.get('fit_mode', 'unknown')}, "
                        f"f_pi0={data_fit['fraction_pi0']:.4f}, "
                        f"f_DVCS={data_fit['fraction_dvcs']:.4f}, "
                        f"f_SB={data_fit['fraction_sideband']:.4f}, "
                        f"AAO morph=({data_fit['morph_shift_bins']:+.2f},"
                        f"{data_fit['morph_smear_sigma_bins']:.2f}) bins, "
                        f"D/ndf={data_fit['nominal_deviance_per_ndf']:.2f}"
                        f"->{data_fit['deviance_per_ndf']:.2f}",
                        flush=True,
                    )
                #endif
                stability = region_summary.get(
                    "sideband_stability",
                    {},
                )
                if stability.get("success", False):
                    print(
                        f"      SIDEBAND SHAPE: JS(0.15-0.19,0.19-0.23)="
                        f"{stability['js_near_far']:.4f}; "
                        f"JS(signal-edge 0.11-0.15, near SB)="
                        f"{stability['js_signal_edge_near']:.4f}",
                        flush=True,
                    )
                #endif
            #endfor
        #endfor
    #endfor

    print("\n" + "=" * 100, flush=True)


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
    if args.optimizer_max_events <= 0:
        raise ValueError("--optimizer-max-events must be > 0.")
    #endif
    if args.optimizer_steps < 0:
        raise ValueError("--optimizer-steps must be >= 0.")
    #endif
    if args.optimizer_quantiles < 3:
        raise ValueError("--optimizer-quantiles must be >= 3.")
    #endif
    if not 0.0 < args.optimizer_min_step_pi0_eff <= 1.0:
        raise ValueError("--optimizer-min-step-pi0-eff must be in (0,1].")
    #endif
    if not 0.0 < args.optimizer_min_total_pi0_eff <= 1.0:
        raise ValueError("--optimizer-min-total-pi0-eff must be in (0,1].")
    #endif
    if args.optimizer_min_improvement < 1.0:
        raise ValueError("--optimizer-min-improvement must be >= 1.")
    #endif
    if (
        args.optimizer_data_sideband_min
        >= args.optimizer_data_sideband_max
    ):
        raise ValueError(
            "--optimizer-data-sideband-min must be < "
            "--optimizer-data-sideband-max."
        )
    #endif
    if args.optimizer_data_sideband_weight < 0.0:
        raise ValueError("--optimizer-data-sideband-weight must be >= 0.")
    #endif
    if args.closure_toys <= 0:
        raise ValueError("--closure-toys must be > 0.")
    #endif
    if args.closure_events <= 0:
        raise ValueError("--closure-events must be > 0.")
    #endif
    if args.closure_bins < 4:
        raise ValueError("--closure-bins must be >= 4.")
    #endif
    if args.data_fit_max_closure_bias < 0.0:
        raise ValueError("--data-fit-max-closure-bias must be >= 0.")
    #endif
    if args.data_fit_max_closure_rms <= 0.0:
        raise ValueError("--data-fit-max-closure-rms must be > 0.")
    #endif
    if args.data_fit_background_near_best < 1.0:
        raise ValueError(
            "--data-fit-background-near-best must be >= 1."
        )
    #endif
    if args.data_fit_morph_shift_max_bins < 0.0:
        raise ValueError("--data-fit-morph-shift-max-bins must be >= 0.")
    #endif
    if args.data_fit_morph_smear_max_bins < 0.0:
        raise ValueError("--data-fit-morph-smear-max-bins must be >= 0.")
    #endif
    if args.data_fit_morph_shift_steps < 3:
        raise ValueError("--data-fit-morph-shift-steps must be >= 3.")
    #endif
    if args.data_fit_morph_smear_steps < 2:
        raise ValueError("--data-fit-morph-smear-steps must be >= 2.")
    #endif
    if args.data_fit_sideband_min_events < 20:
        raise ValueError("--data-fit-sideband-min-events must be >= 20.")
    #endif
    if args.data_fit_sideband_max_js < 0.0:
        raise ValueError("--data-fit-sideband-max-js must be >= 0.")
    #endif
    if not (0.0 <= args.coverage_min_eppi0_overlap <= 1.0):
        raise ValueError(
            "--coverage-min-eppi0-overlap must be between 0 and 1."
        )
    #endif
    if args.efficiency_eppi0_max_entries < 0:
        raise ValueError("--efficiency-eppi0-max-entries must be >= 0.")
    #endif
    if args.parent_component_tol <= 0.0:
        raise ValueError("--parent-component-tol must be > 0.")
    #endif
    if args.parent_distance_max <= 0.0:
        raise ValueError("--parent-distance-max must be > 0.")
    #endif
    if args.kdtree_workers <= 0:
        raise ValueError("--kdtree-workers must be > 0.")
    #endif
    if not args.assoc_mgg_min < args.assoc_mgg_max:
        raise ValueError("--assoc-mgg-min must be < --assoc-mgg-max.")
    #endif
    if args.assoc_remainder_mass2_max <= 0.0:
        raise ValueError("--assoc-remainder-mass2-max must be > 0.")
    #endif
    if args.assoc_probe_energy_min < 0.0:
        raise ValueError("--assoc-probe-energy-min must be >= 0.")
    #endif
    if args.probe_match_nsigma <= 0.0:
        raise ValueError("--probe-match-nsigma must be > 0.")
    #endif
    if args.probe_match_min_resolution_events < 20:
        raise ValueError("--probe-match-min-resolution-events must be >= 20.")
    #endif

    selected_keys = set(args.period or [p.key for p in PERIODS])
    selected_periods = [p for p in PERIODS if p.key in selected_keys]

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.coverage_only:
        log(
            "Coverage-only mode: scanning FULL nSidis epgamma/eppi0 trees. "
            "--max-entries is intentionally ignored."
        )
        failures = []
        for period in selected_periods:
            try:
                run_raw_coverage_audit(
                    period,
                    args.tree,
                    args.step_size,
                    output_dir,
                    min_eppi0_overlap=args.coverage_min_eppi0_overlap,
                    fail_on_pathology=False,
                )
            except Exception as exc:
                failures.append((period.label, str(exc)))
            #endtry
        #endfor

        if failures:
            print("\nCoverage-audit failures:", flush=True)
            for label, message in failures:
                print(f"  {label}: {message}", flush=True)
            #endfor
            return 1
        #endif

        print(
            "\nCoverage-only audit complete. "
            f"See {output_dir / 'coverage'}/<period>/.",
            flush=True,
        )
        return 0
    #endif

    log(
        "Standalone photon-efficiency concept study: "
        "data + AAO pi0 + DVCSgen only; CLASDIS disabled."
    )
    log(
        "Shape-comparison y-axis: "
        + ("linear (--linear override)." if args.linear else "logarithmic (default).")
    )
    if args.no_optimize_cuts:
        log("Rectangular-cut optimizer: disabled.")
    else:
        log(
            "Rectangular-cut optimizer: enabled independently for every "
            "period and separately in FT/FD; "
            f"reservoir <= {args.optimizer_max_events:,} events/class/region, "
            f"max steps={args.optimizer_steps}; "
            f"loose Mx2_1<{args.optimizer_mx2_1_preselection_max:.3g} GeV^2; "
            f"data sideband {args.optimizer_data_sideband_min:.3g}<=Mx2_1<"
            f"{args.optimizer_data_sideband_max:.3g} GeV^2; "
            f"sideband weight={args.optimizer_data_sideband_weight:.3g}; "
            "E_probe (= Emiss2) bins: "
            "0.4-1, 1-2, 2-4, 4-9.5 GeV; "
            "Emiss2 excluded from optimizer thresholds; "
            f"closure={args.closure_toys} toys x "
            f"{args.closure_events} events at f_pi0=0.2,0.5,0.8; "
            "real-data AAO+DVCS fit enabled with closure-safe operating-step "
            f"selection (|bias|<={args.data_fit_max_closure_bias:g}, "
            f"RMS<={args.data_fit_max_closure_rms:g}); "
            "real-data fit uses target Mx2_1<0.15 plus AAO + DVCSgen + "
            "near empirical sideband 0.15-0.19 GeV^2, with "
            f"AAO-only |shift|<={args.data_fit_morph_shift_max_bins:g} "
            f"and smear<={args.data_fit_morph_smear_max_bins:g} bins; "
            "third-template sideband accepted only when near/far each have "
            f">={args.data_fit_sideband_min_events} events and "
            f"JS<={args.data_fit_sideband_max_js:g}, otherwise AAO+DVCS "
            "fallback is used; TAG-FD / predicted-probe-sector diagnostics keep the global FD "
            "composition/cuts/discriminator/morph fixed in sectors 1-6 and "
            "do not refit f_pi0 by sector. FT/FD composition categories "
            "come from reconstructed tag p2_theta; only the downstream probe "
            "theta/phi/sector comes from p_beam - p_e - p_p - p_gamma,tag."
        )
    #endif
    log(
        "Raw data-file coverage audit: enabled before efficiency extraction. "
        "FULL eppi0 run-event keys are searched in streamed epgamma chunks, independent of "
        "--max-entries and row ordering; "
        f"required common/eppi0 >= {args.coverage_min_eppi0_overlap:.0%}."
    )
    if args.skip_efficiency:
        log("Final tag-probe efficiency extraction: disabled (--skip-efficiency).")
    else:
        log(
            "Final tag-probe efficiency extraction: enabled. "
            "Tag FT/FD from reconstructed p2; predicted probe determines FT/FD "
            "efficiency bin and FD sector. FD sectors 1-6 below 4 GeV; "
            "4-9.5 GeV FD combined. DATA association = exact run/event + e/p; "
            "AAO association = e/p cKDTree only (MC run/ev never used); "
            f"probe match = AAO-derived +/-{args.probe_match_nsigma:g} sigma "
            "windows in Delta px,py,pz."
        )
    #endif

    log(
        "Common pre-binning selection: "
        f"Angle(e,gamma) > {THETA_EGAMMA_MIN_DEG:.1f} deg. "
        f"Photon regions: FT {FT_THETA_MIN_DEG:.1f}-{FT_THETA_MAX_DEG:.1f} deg, "
        f"FD {FD_THETA_MIN_DEG:.1f}-{FD_THETA_MAX_DEG:.1f} deg."
    )

    n_workers = min(int(args.workers), len(selected_periods))
    log(
        f"Period-level parallelism: {n_workers} process(es) for "
        f"{len(selected_periods)} selected period(s)."
    )

    summaries_by_period: Dict[str, dict] = {}

    if n_workers == 1:
        for period in selected_periods:
            summaries_by_period[period.key] = process_period(
                period,
                args.tree,
                args.max_entries,
                args.step_size,
                output_dir,
                args.linear,
                not args.no_optimize_cuts,
                args.optimizer_max_events,
                args.optimizer_steps,
                args.optimizer_quantiles,
                args.optimizer_min_step_pi0_eff,
                args.optimizer_min_total_pi0_eff,
                args.optimizer_min_improvement,
                args.optimizer_mx2_1_preselection_max,
                args.optimizer_data_sideband_min,
                args.optimizer_data_sideband_max,
                args.optimizer_data_sideband_weight,
                args.closure_toys,
                args.closure_events,
                args.closure_bins,
                args.data_fit_max_closure_bias,
                args.data_fit_max_closure_rms,
                args.data_fit_background_near_best,
                args.data_fit_morph_shift_max_bins,
                args.data_fit_morph_smear_max_bins,
                args.data_fit_morph_shift_steps,
                args.data_fit_morph_smear_steps,
                args.data_fit_sideband_min_events,
                args.data_fit_sideband_max_js,
                not args.skip_efficiency,
                args.efficiency_eppi0_max_entries,
                args.parent_component_tol,
                args.parent_distance_max,
                args.kdtree_workers,
                args.assoc_mgg_min,
                args.assoc_mgg_max,
                args.assoc_remainder_mass2_max,
                args.assoc_probe_energy_min,
                args.probe_match_nsigma,
                args.probe_match_min_resolution_events,
                args.coverage_min_eppi0_overlap,
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
                    args.linear,
                    not args.no_optimize_cuts,
                    args.optimizer_max_events,
                    args.optimizer_steps,
                    args.optimizer_quantiles,
                    args.optimizer_min_step_pi0_eff,
                    args.optimizer_min_total_pi0_eff,
                    args.optimizer_min_improvement,
                    args.optimizer_mx2_1_preselection_max,
                    args.optimizer_data_sideband_min,
                    args.optimizer_data_sideband_max,
                    args.optimizer_data_sideband_weight,
                    args.closure_toys,
                    args.closure_events,
                    args.closure_bins,
                    args.data_fit_max_closure_bias,
                    args.data_fit_max_closure_rms,
                    args.data_fit_background_near_best,
                    args.data_fit_morph_shift_max_bins,
                    args.data_fit_morph_smear_max_bins,
                    args.data_fit_morph_shift_steps,
                    args.data_fit_morph_smear_steps,
                    args.data_fit_sideband_min_events,
                    args.data_fit_sideband_max_js,
                    not args.skip_efficiency,
                    args.efficiency_eppi0_max_entries,
                    args.parent_component_tol,
                    args.parent_distance_max,
                    args.kdtree_workers,
                    args.assoc_mgg_min,
                    args.assoc_mgg_max,
                    args.assoc_remainder_mass2_max,
                    args.assoc_probe_energy_min,
                    args.probe_match_nsigma,
                    args.probe_match_min_resolution_events,
                    args.coverage_min_eppi0_overlap,
                ): period
                for period in selected_periods
            }

            for future in as_completed(future_to_period):
                period = future_to_period[future]
                try:
                    summaries_by_period[period.key] = future.result()
                except Exception as exc:
                    raise RuntimeError(
                        f"Parallel concept-study processing failed for "
                        f"{period.label}: {exc}"
                    ) from exc
                #endtry
            #endfor
        #endwith
    #endif

    if not args.no_optimize_cuts:
        print_optimizer_summary(
            summaries_by_period,
            selected_periods,
        )
    #endif

    log(f"Done. Outputs are in {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
