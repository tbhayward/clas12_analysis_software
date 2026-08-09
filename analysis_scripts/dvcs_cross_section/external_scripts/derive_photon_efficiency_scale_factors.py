#!/usr/bin/env python3
"""
derive_photon_efficiency_scale_factors_v78_simple_summary_json_fix.py

Photon-efficiency study with separate exact-one-photon and exact-two-photon
data categories. Events with three or more reconstructed-photon entries are
rejected.

Data-side changes in this revision:

  * remove the Emiss2 > 1 GeV lower cut;
  * remove the Mx2_2 > 2 GeV^2 lower cut;
  * remove Mx2_2 from the data shape comparison and data template fit;
  * retain the complete Emiss2 spectrum, including the BH/DVCS peak near zero;
  * define the after-cuts row using only pTmiss > 0.5 GeV, in addition to the
    open_angle_ep2 > 5 deg requirement common to both rows.

The four fitted data projections are Delta_phi, theta_gamma_gamma, pTmiss, and
Emiss2. Each is fitted independently, followed by a shared pi0-fraction fit
with variable-specific template morphing. Mx2_2 remains available to the
established optional MC-efficiency machinery elsewhere in this script.

Template-fit plotting conventions
---------------------------------

The advanced data-template canvases use the same visual conventions as
plot_exclusivity_data_dvcs_pi0_mc.py:

  * black points: e'p'gamma data;
  * thin dashed blue: raw DVCS MC shape;
  * thin dashed red: raw epi0 MC shape;
  * solid blue: fitted DVCS component;
  * solid red: fitted epi0 component;
  * dashed green: total two-template fit;
  * one shared legend centered above the panels;
  * per-panel fit-quality annotation in the upper-right corner;
  * lowercase 'events / bin' and 'data / fit' axis labels.



Shape-comparison kinematic diagnostics
--------------------------------------

The one-photon and two-photon before/after shape canvases additionally show
p1_p, p1_theta, p2_p, and p2_theta from the native epgamma tree branches.
These four distributions are diagnostic only and do not enter the pi0-fraction
template fits. Their default plotting ranges are 0--5 GeV, 0--70 deg,
0--10 GeV, and 0--40 deg, respectively.

Shape-comparison diagnostics
----------------------------

In addition to the four data-fit variables, the before/after shape-comparison
canvases include p1_p, p1_theta, p2_p, and p2_theta. These particle-kinematic
distributions are diagnostic only and do not enter the fitted pi0 fraction.

Plotting conventions
--------------------

Both the shape-comparison and template-fit canvases follow the visual
conventions used by plot_exclusivity_data_dvcs_pi0_mc.py: black data, blue
DVCS MC, red epi0 MC, dashed green total fit, thin dashed raw-template
outlines, solid fitted components, a shared frame-free legend, and matching
axis-label and horizontal-grid conventions.

Kinematic diagnostic canvases
-----------------------------

The four additional particle-kinematic diagnostics are written on separate
before/after canvases rather than appended to the fitted-variable canvases:

    p1_p, p1_theta, p2_p, p2_theta.

The p1_theta and p2_theta ROOT branches are stored in radians. They are
converted entry-by-entry to degrees before histogramming. Their configured
histogram limits and axis labels are therefore interpreted in degrees.

For every run period, separate exact-one-photon and exact-two-photon canvases
are produced for:

  * the four fitted shape variables;
  * the four particle-kinematic diagnostics.

Multiplicity reset after corrected epgamma production
------------------------------------------------------

The epgamma files were regenerated after correcting the photon loop so that
the leading reconstructed photon is once again included. This revision removes
the temporary hard rejection of events represented by three or more photon
entries.

No event is rejected solely because of photon-entry multiplicity.

The data diagnostic stage now keeps four independent categories:

  * exactly one photon entry;
  * exactly two photon entries;
  * three or more photon entries;
  * all multiplicities combined.

Multiplicity histograms and CSV/JSON summaries still report exact
multiplicities 1 through 6 and >6. The template fits remain restricted to the
exact-one and exact-two categories for now; the >=3 and inclusive categories
are diagnostic only.

The MC-efficiency opportunity stage also no longer requires multiplicity <3.
Every entry with a valid event identity is allowed to proceed to the ordinary
kinematic and topology selections.

Revision: restored MC stage and corrected multiplicity summaries
---------------------------------------------------------------

The full AAOGEN MC photon-efficiency stage is again enabled by default and
writes under output/photon_efficiency_study/mc/. It may be disabled explicitly
with --skip-mc-efficiency. The legacy --run-mc-efficiency option remains
accepted for command-line compatibility.

The temporary p1_p, p1_theta, p2_p, and p2_theta shape-comparison canvases
have been removed. The data shape-comparison stage now uses only:

    Delta_phi, theta_gamma_gamma, pTmiss, Emiss2.

A multiplicity-summary bug has also been corrected. The previous version
passed the array of already-unique event signatures into a helper that called
np.unique() again. Consequently every event appeared to have multiplicity one
in the summary plots even though the exact-one, exact-two, and >=3 masks were
correct. The summary is now constructed directly from the multiplicities found
during the identity prepass.

Revision: restore multiplicity cut and updated diagnostic ranges
---------------------------------------------------------------

This revision restores the hard event-level requirement that only events with
one or two reconstructed-photon entries are retained. Events with three or
more reconstructed-photon entries are rejected consistently in the data-shape
stage and the MC-efficiency opportunity stage.

The data shape-comparison ranges are now:

    theta_gamma_gamma: 0 to 10 deg
    pTmiss:            0 to 1.25 GeV
    Emiss2:            0 to 4 GeV

No pTmiss lower cut is applied in the after-selection row.

The four particle-kinematic diagnostics are restored:

    p1_p, p1_theta, p2_p, p2_theta

and are placed in the bottom row of the same shape-comparison canvas. The top
row contains the four exclusivity-shape variables.

Revision: three independent reconstructed-photon template cases
---------------------------------------------------------------

Shape comparisons and template fits are produced separately for:

  1. exact-one-photon events;
  2. exact-two-photon events, selecting the more energetic p2 entry;
  3. exact-two-photon events, selecting the less energetic p2 entry.

The two entries in every exact-two event are ranked by p2_p. Global ROOT entry
index is used as a deterministic tie breaker. Each exact-two event contributes
exactly once to each energy-ranked case. Multiplicity >=3 remains rejected.

Revision: algorithmic runtime optimizations
-------------------------------------------

This revision preserves the physics selections, three photon-template cases,
fit models, plotting products, and MC-efficiency outputs while reducing the
dominant runtime costs.

Data-stage changes:

  * exact-two event ranking is vectorized after one stable sort rather than
    repeatedly scanning the complete inverse-index array for every event;
  * the multiplicity prepass returns one dense uint8 category code per ROOT
    entry, so the histogram pass uses direct slices instead of repeated
    structured-signature construction and np.isin membership tests;
  * the duplicate after-selection histograms are no longer filled a second
    time when no additional cut is active; the serialized after-selection
    products are copied from the same finalized counts;
  * per-file timing is recorded for the ranking prepass and histogram pass.

MC-stage changes:

  * read_opportunities() can skip construction and exact-deduplication of the
    unused photon-candidate catalog;
  * the MC-efficiency processor uses that fast path;
  * read_epgg() applies entry_stop directly in uproot.iterate();
  * scalar opening-angle calculations replace one-element NumPy allocations
    in the truth-partner matching hot loop.

Fit-stage changes:

  * --fast-fits provides an explicit diagnostic mode with fewer starts,
    looser optimizer tolerances, and fewer coordinate iterations;
  * invariant fit-context arrays, bin widths, and two-component template
    partitions are prepared once and reused by objective evaluations;
  * stage timings are written into the data and MC metadata.

Revision: category-specific ranges and exact-two diphoton-mass diagnostics
-------------------------------------------------------------------------

The one-photon case now uses the restricted ranges:

    p1_p:              0 to 1.25 GeV
    theta_gamma_gamma: 0 to 3 deg
    pTmiss:            0 to 0.7 GeV
    Emiss2:            0 to 2 GeV

The two exact-two-photon cases retain the broader fit-variable ranges:

    theta_gamma_gamma: 0 to 10 deg
    pTmiss:            0 to 1.25 GeV
    Emiss2:            0 to 4 GeV

For exact-two events, the bottom row of the shape-comparison canvas no longer
contains p1_p or p1_theta. It instead contains the reconstructed diphoton mass,
followed by p2_p and p2_theta. The diphoton mass is calculated directly from
the two reconstructed photon four-vectors in each exact-two event:

    M_gg^2 = 2 E1 E2 (1 - cos(theta_12)).

The same event-level diphoton mass is attached to both the more-energetic and
less-energetic exact-two entries. The shape-comparison plot marks the nominal
pi0 window 0.10 < M_gg < 0.17 GeV and reports the corresponding data fraction.

Revision: p2_phi branch-resolution fix
--------------------------------------

The exact-two diphoton-mass prepass requires p2_phi in addition to p2_p and
p2_theta. This revision adds p2_phi as a direct logical branch alias and
includes it in the branch-resolution request made by the data shape processor.

Revision: exact-two pi0 mass cut and before/after canvases
----------------------------------------------------------

The nominal reconstructed pi0 mass window is now:

    0.11 < M_gamma_gamma < 0.16 GeV.

For exact-one-photon events, the before and after populations are identical:
no diphoton-mass requirement can be applied.

For both exact-two-photon categories, the after population requires the event
diphoton mass to lie inside the nominal pi0 window. The more-energetic and
less-energetic entries from a given exact-two event therefore receive the same
event-level pass/fail decision.

The shape-comparison canvases are now 2 rows by 5 columns:

    top row:    before the exact-two M_gamma_gamma cut;
    bottom row: after the exact-two M_gamma_gamma cut.

The five displayed variables are:

    Delta_phi, theta_gamma_gamma, pTmiss, Emiss2, and
    either p1_p for the exact-one case or M_gamma_gamma for exact-two cases.

The p2_p and p2_theta shape diagnostics are removed completely.

Revision: exact-two internal branch-resolution fix
--------------------------------------------------

Removing p2_p and p2_theta from the displayed shape diagnostics must not remove
them from the exact-two event machinery. The photon-ranking and diphoton-mass
prepass still requires p2_p, p2_theta, and p2_phi. This revision explicitly
resolves all three internal branches even though only M_gamma_gamma is displayed
from that information.

Revision: exclusivity-style coherent one-photon template morphing
-----------------------------------------------------------------

The exact-one-photon template fits now deliberately follow the validated
plot_exclusivity_data_dvcs_pi0_mc.py strategy as closely as possible.

For exact-one events, both DVCSGEN and AAOGEN receive independent nuisance
parameters, but each complete template is morphed coherently:

    Delta_phi:          additive shift + Gaussian broadening;
    theta_gamma_gamma:  log(theta+epsilon) shift + Gaussian broadening;
    pTmiss:             log(pTmiss+epsilon) shift + Gaussian broadening;
    Emiss2:             additive shift + Gaussian broadening.

The former two-component theta_gamma_gamma model and affine/asymmetric pTmiss
and Emiss2 models are therefore not used in the exact-one fit. This prevents a
single DVCS template from being split into independently movable sub-peaks.

The coherent transforms use the same integrated-Gaussian transport approach as
the validated exclusivity fitter, including the same nuisance-prior philosophy:
additive shifts/smears use --fit-shift-prior-bins and
--fit-smear-prior-bins, while positive/log variables use fixed transformed-space
prior widths of 0.20 for shifts and 0.40 for broadenings.

Exact-two events retain the established category-specific modeling because the
0.11 < M_gamma_gamma < 0.16 GeV selection already supplies the primary
two-photon discrimination.

Emiss2 is restricted to 0--2 GeV for every category.

Independent per-variable fractions remain diagnostics. The shared fraction is
still determined separately for each run period; no cross-period fraction
constraint is imposed.

Revision: category-aware model-callsite consistency fix
-------------------------------------------------------

variable_model_spec() became category dependent in the exclusivity-style
one-photon rewrite. This revision updates every remaining call site so the
multiplicity category is always explicit. Startup/preflight and JSON metadata
now record model definitions separately for the one-photon and exact-two
categories.

The fraction-consistency and flattened fit-row builders were also made
explicitly category-aware rather than relying on the backward-compatible
DATA_SHAPE_VARIABLES alias.

Revision: shared-fit category propagation fix
---------------------------------------------

AdvancedSharedFit now stores multiplicity_key explicitly. The category-aware
fit plotting, fraction-consistency plotting, and flattened fit-row machinery
introduced in the previous revision can therefore retrieve the photon category
directly from the returned fit summary.

Revision: raw-template primary diagnostic versus morphed comparison
-------------------------------------------------------------------

For each fit category the script now evaluates a strict raw-template,
fraction-only decomposition in addition to the constrained morphed-template
fit. This is especially important for the exact-one-photon sample, where the
raw DVCSGEN and AAOGEN shapes visibly bracket the data well.

The raw model has exactly one shared physics parameter:

    D_i = N_i [(1-f_pi0) T_DVCS,i^raw + f_pi0 T_pi0,i^raw].

No shifts, smearing, affine transformations, component splitting, or other
shape nuisance parameters enter the raw-template result.

For every variable the script records:
    * the independent raw-template fraction;
    * the independent morphed-template fraction;
    * the shared raw-template fraction;
    * the shared morphed-template fraction;
    * Delta f_morph = f_morphed - f_raw;
    * raw and morphed Poisson deviances.

The fit canvases overlay both total models and show both data/model ratios.
The fraction-consistency plots compare raw and morphed independent fractions
directly. Existing morphed-fit fields are retained for backward compatibility.

The raw-template result is labeled as the primary diagnostic for the exact-one
sample; the morphed result is a model-sensitivity comparison.

Revision: temporary --simple data-side efficiency override
-----------------------------------------------------------

A temporary --simple mode is available as a complete override of the normal
DATA section only. It does not remove or alter the existing data-template or
MC-efficiency implementations.

In --simple mode the epgamma DATA tree is split entry-by-entry at 2 GeV:

    low-energy tags:      p2_p < 2 GeV;
    high-energy partners: p2_p > 2 GeV.

Every low-energy photon is provisionally assumed to be the lower-energy photon
from a pi0 decay whose other photon should lie above 2 GeV. The expected
high-energy partner direction is reconstructed from the beam, target, electron,
proton, and low-energy photon four-vectors using the same reconstruct_probe()
helper as the established MC-efficiency study.

Each low-energy tag is assigned to FT or to one of the six FD sectors from the
predicted partner direction. A tag is counted as reconstructed when the same
(runnum, evnum) occurs anywhere in the >2 GeV photon set. No DVCS/pi0 template
decomposition, invariant-mass requirement, angular partner requirement, or
other data-side fit is used in this temporary mode.

The simple efficiency in each predicted detector region is therefore:

    epsilon = N(low-E tags whose event has a >2 GeV photon)
              / N(all low-E tags in that predicted region).

Results are printed to the console for FT and FD sectors 1--6, separately by
run period and for all selected periods combined. They are also written under:

    output/photon_efficiency_study/data/simple/

The normal data analysis remains the default when --simple is absent. The
separate MC-efficiency stage is unaffected and still runs unless --no-mc /
--skip-mc-efficiency is supplied.

Revision: --simple CSV-writer fix
---------------------------------

The initial --simple implementation called write_rows_csv(), but the main
script has no function with that name. This revision adds a small local CSV
writer dedicated to the temporary simple-data audit and updates the simple
stage to use it.

No physics logic or efficiency definitions are changed.

Revision: --simple JSON-writer fix
----------------------------------

The temporary --simple stage previously passed default=json_default to
json.dump(), but no json_default() helper exists in this script. This revision
adds a dedicated serializer for the simple-mode audit payload.

No physics, matching, efficiency, detector-classification, or normal-analysis
logic is changed.

Revision: --simple same-event pi0 mass requirement
-------------------------------------------------

The temporary --simple DATA override now requires the low-energy photon and at
least one >2 GeV photon from the same exact (runnum, evnum) to form a
reconstructed pi0 candidate in the nominal mass window:

    0.11 < M_gamma_gamma < 0.16 GeV.

The denominator and predicted FT/FD-region assignment are unchanged. If an
event contains multiple >2 GeV photons, all are tested and the low-energy tag
is counted once if any candidate passes the diphoton-mass requirement.

Revision: --simple reconstructed-data versus reconstructed-AAOGEN comparison
----------------------------------------------------------------------------

When --simple is enabled, the script now runs the same reconstructed-photon
algorithm on the AAOGEN epgamma sample in addition to data.

For both samples:
  * every reconstructed photon with p2_p < 2 GeV is a denominator opportunity;
  * the expected partner direction is reconstructed from beam+target-e-p-low_gamma;
  * the opportunity is assigned to FT or FD sector 1--6 from that direction;
  * success requires a reconstructed photon with p2_p > 2 GeV from the same
    reconstructed event and 0.11 < M_gamma_gamma < 0.16 GeV.

Data event identity is exact (runnum, evnum). AAOGEN epgamma files do not have
usable run/event identifiers, so reconstructed MC entries are grouped with the
same established quantized electron+proton kinematic identity used elsewhere
in this script for MC multiplicity reconstruction.

This simple reconstructed-AAOGEN diagnostic uses no truth matching and does
not replace the established AAOGEN MC-efficiency stage. At the very end of a
--simple run, the script prints a compact data/simple-MC comparison for FT and
aggregate FD, period by period and combined.

Revision: --simple FT/FD data-vs-AAOGEN diagnostic package
----------------------------------------------------------

The temporary --simple mode now records and plots reconstructed-level
diagnostics for data and reconstructed AAOGEN using the identical simple
selection.

A new top-level parent directory is created:

    output/photon_efficiency_study/simple/

with per-period diagnostics in:

    output/photon_efficiency_study/simple/diagnostics/

Two 2x3 canvases are written for every run period. Rows are predicted FT and
aggregate predicted FD. The first canvas contains:

    1. low-energy tag photon momentum/energy, p_gamma < 2 GeV;
    2. predicted missing-partner theta_X;
    3. reconstructed M_gamma_gamma for every same-event low/high pair before
       the pi0 mass requirement.

The second canvas contains:

    1. low-energy tag photon theta;
    2. same-event high-energy photon energy before the pi0 mass requirement;
    3. high-energy photon energy for pairs satisfying
       0.11 < M_gamma_gamma < 0.16 GeV.

All data and reconstructed-AAOGEN distributions are independently unit
normalized on the diagnostic plots. Predicted-theta panels mark the FT and FD
acceptance boundaries, and M_gamma_gamma panels mark the nominal pi0 window.

The diagnostic bookkeeping is accumulated during the existing --simple pass;
the ROOT files are not reread for plotting. Existing simple numerical outputs,
the reconstructed-data/reconstructed-AAOGEN comparison, the full data-template
analysis, and the truth-informed AAOGEN MC-efficiency stage are all preserved.

Revision: --simple predicted-theta diagnostic bookkeeping fix
--------------------------------------------------------------

The v76 diagnostic package declared low_probe_theta_deg_chunks in the simple
DATA processor but failed to append the reconstructed probe theta values to
that list. Low-energy opportunities therefore existed while the diagnostic
theta list remained empty, causing np.concatenate() to fail.

This revision explicitly computes probe_theta_deg once in the data processor,
uses it for detector classification, and appends it alongside probe energy.
Both simple DATA and reconstructed-AAOGEN processors also receive guarded
diagnostic concatenation and explicit array-length consistency checks.

No event selection, pi0 mass logic, matching, efficiency definition, detector
boundaries, template analysis, or truth-informed MC-efficiency logic changes.

Revision: --simple top-level JSON serialization fix
---------------------------------------------------

The simple FT/FD diagnostic histograms contain NumPy arrays. The dedicated
simple-data and reconstructed-AAOGEN JSON outputs already serialize those with
simple_json_default(), but the final top-level study-summary JSON did not.
This revision applies the same serializer to that final summary write.

No physics or analysis logic changes.

"""


from __future__ import annotations

import argparse
import concurrent.futures
import csv
import json
import math
import os
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", f"/tmp/matplotlib-{os.getuid()}")
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import minimize, minimize_scalar
from scipy.special import ndtr
from scipy.signal import find_peaks

try:
    import uproot
except ImportError as exc:
    raise SystemExit("This script requires uproot: python3 -m pip install uproot") from exc
#endif



TREE_NAME = "PhysicsEvents"
DEFAULT_OUTPUT_DIR = "output/photon_efficiency_study"
DEFAULT_STEP_SIZE = "200 MB"
MAX_WORKERS = 8

PROTON_MASS_GEV = 0.9382720813
ELECTRON_MASS_GEV = 0.00051099895


@dataclass(frozen=True)
class PeriodConfig:
    key: str
    label: str
    beam_energy_GeV: float
    epg_data: str
    epgg_data: str
    dvcs_mc: str
    pi0_epg_mc: str
    pi0_epgg_mc: str


PERIODS: Tuple[PeriodConfig, ...] = (
    PeriodConfig(
        "fa18_inb", "Fa18 Inb", 10.604,
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/rga_fa18_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_fa18_inb_eppi0.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/dvcsgen_rga_fa18_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/bkg_rga_fa18_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_fa18_inb_50nA_10604MeV.root",
    ),
    PeriodConfig(
        "fa18_out", "Fa18 Out", 10.604,
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/rga_fa18_out_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_fa18_out_eppi0.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/dvcsgen_rga_fa18_out_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/bkg_rga_fa18_out_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_fa18_out_50nA_10604MeV.root",
    ),
    PeriodConfig(
        "sp18_inb", "Sp18 Inb", 10.594,
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/rga_sp18_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_sp18_inb_eppi0.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/dvcsgen_rga_sp18_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/bkg_rga_sp18_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_sp18_inb_50nA_10594MeV.root",
    ),
    PeriodConfig(
        "sp18_out", "Sp18 Out", 10.594,
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/rga_sp18_out_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_sp18_out_eppi0.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/dvcsgen_rga_sp18_out_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/bkg_rga_sp18_out_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_sp18_out_45nA_10594MeV.root",
    ),
    PeriodConfig(
        "sp19_inb", "Sp19 Inb", 10.200,
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/rga_sp19_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_sp19_inb_eppi0.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/dvcsgen_rga_sp19_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/bkg_rga_sp19_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_sp19_inb_50nA_10200MeV.root",
    ),
)


ALIASES: Mapping[str, Tuple[str, ...]] = {
    "runnum": ("runnum", "run", "run_number"),
    "eventnum": ("evnum", "eventnum", "event", "event_number"),
    "e_p": ("e_p", "p_e", "electron_p"),
    "e_theta": ("e_theta", "theta_e", "electron_theta"),
    "e_phi": ("e_phi", "phi_e", "electron_phi"),
    "p_p": ("p1_p", "p_p", "proton_p"),
    "p_theta": ("p1_theta", "p_theta", "proton_theta"),
    "p_phi": ("p1_phi", "p_phi", "proton_phi"),
    "p1_p": ("p1_p",),
    "p1_theta": ("p1_theta",),
    "p2_p": ("p2_p",),
    "p2_theta": ("p2_theta",),
    "p2_phi": ("p2_phi",),
    "tag_E": ("p2_p", "g1_p", "g1_E", "gamma1_p", "photon1_p"),
    "tag_theta": ("p2_theta", "g1_theta", "gamma1_theta", "photon1_theta"),
    "tag_phi": ("p2_phi", "g1_phi", "gamma1_phi", "photon1_phi"),
    "tag_detector": ("detector2", "g1_detector", "gamma1_detector"),
    "proton_detector": ("detector1", "p_detector", "proton_detector"),
    "fiducial_status": ("fiducial_status",),
    "Q2": ("Q2", "q2"),
    "W": ("W", "w"),
    "y": ("y", "inelasticity"),
    "z": ("z",),
    "t1": ("t1", "t", "minus_t"),
    "open_angle_ep2": ("open_angle_ep2", "open_angle_eg", "open_angle_ephoton"),
    "Delta_phi": ("Delta_phi", "delta_phi", "dphi"),
    "theta_cm": ("theta", "theta2", "theta_2"),
    "theta_gamma_gamma": ("theta_gamma_gamma", "theta_pi0_pi0"),
    "pTmiss": ("pTmiss", "ptmiss", "pT_miss"),
    "Emiss2": ("Emiss2", "E_miss2", "Emiss_sq"),
    "Mx2": ("Mx2", "Mx2_epg", "Mx2_epgamma", "Mx2_eppi0", "Mx2_epi0"),
    "Mx2_1": ("Mx2_1", "Mx2_ep", "Mx2_x1", "Mx2_proton", "Mx2_p"),
    "Mx2_2": ("Mx2_2", "Mx2_egamma", "Mx2_gamma", "Mx2_pi0", "Mx2_x2"),
    # Native epgammagamma reconstruction.
    "pi0_p": ("p2_p", "pi0_p"),
    "pi0_theta": ("p2_theta", "pi0_theta"),
    "pi0_phi": ("p2_phi", "pi0_phi"),
    "pi0_mass": ("Mh_gammagamma",),
    "open1": ("open_angle_egamma1",),
    "open2": ("open_angle_egamma2",),
    "trento1": ("gamma_phi1",),
    "trento2": ("gamma_phi2",),
    "native_gamma_detector1": ("detector_gamma1",),
    "native_gamma_detector2": ("detector_gamma2",),
}


@dataclass(frozen=True)
class FitVariable:
    key: str
    label: str
    bins: int
    low: float
    high: float


FIT_VARIABLES: Tuple[FitVariable, ...] = (
    FitVariable(
        "Delta_phi",
        r"$\Delta\phi$ (rad)",
        120,
        0.5 * math.pi,
        1.5 * math.pi,
    ),
    FitVariable(
        "theta_cm",
        r"$\theta_{p\gamma}^{\rm CM}$ (rad)",
        120,
        1.0,
        math.pi,
    ),
    FitVariable(
        "theta_gamma_gamma",
        r"$\theta_{\gamma\gamma}$ (deg)",
        120,
        0.0,
        40.0,
    ),
    FitVariable(
        "pTmiss",
        r"$p_T^{\rm miss}$ (GeV)",
        125,
        0.0,
        2.0,
    ),
    FitVariable(
        "Emiss2",
        r"$E_{\rm miss}$ (GeV)",
        120,
        0.0,
        9.0,
    ),
    FitVariable(
        "Mx2",
        r"$M_x^2$ (GeV$^2$)",
        100,
        -0.03,
        0.03,
    ),
    FitVariable(
        "Mx2_2",
        r"$M_{x2}^2$ (GeV$^2$)",
        160,
        0.0,
        16.0,
    ),
)


ONE_PHOTON_FIT_VARIABLES: Tuple[FitVariable, ...] = (
    FitVariable(
        "Delta_phi",
        r"$\Delta\phi$ (rad)",
        120,
        0.5 * math.pi,
        1.5 * math.pi,
    ),
    FitVariable(
        "theta_gamma_gamma",
        r"$\theta_{\gamma\gamma}$ (deg)",
        120,
        0.0,
        3.0,
    ),
    FitVariable(
        "pTmiss",
        r"$p_T^{\mathrm{miss}}$ (GeV)",
        120,
        0.0,
        0.7,
    ),
    FitVariable(
        "Emiss2",
        r"$E_{\mathrm{miss}}$ (GeV)",
        120,
        0.0,
        2.0,
    ),
)

TWO_PHOTON_FIT_VARIABLES: Tuple[FitVariable, ...] = (
    FitVariable(
        "Delta_phi",
        r"$\Delta\phi$ (rad)",
        120,
        0.5 * math.pi,
        1.5 * math.pi,
    ),
    FitVariable(
        "theta_gamma_gamma",
        r"$\theta_{\gamma\gamma}$ (deg)",
        120,
        0.0,
        10.0,
    ),
    FitVariable(
        "pTmiss",
        r"$p_T^{\mathrm{miss}}$ (GeV)",
        120,
        0.0,
        1.25,
    ),
    FitVariable(
        "Emiss2",
        r"$E_{\mathrm{miss}}$ (GeV)",
        120,
        0.0,
        2.0,
    ),
)

ONE_PHOTON_DISPLAY_VARIABLES: Tuple[FitVariable, ...] = (
    *ONE_PHOTON_FIT_VARIABLES,
    FitVariable(
        "p1_p",
        r"$p_1$ momentum (GeV)",
        120,
        0.0,
        1.25,
    ),
)

TWO_PHOTON_DISPLAY_VARIABLES: Tuple[FitVariable, ...] = (
    *TWO_PHOTON_FIT_VARIABLES,
    FitVariable(
        "diphoton_mass",
        r"$M_{\gamma\gamma}$ (GeV)",
        120,
        0.0,
        0.30,
    ),
)

# Backward-compatible aliases used by fit-model metadata.
DATA_SHAPE_VARIABLES: Tuple[FitVariable, ...] = ONE_PHOTON_FIT_VARIABLES
KINEMATIC_SHAPE_VARIABLES: Tuple[FitVariable, ...] = (
    ONE_PHOTON_DISPLAY_VARIABLES[4:],
)
SHAPE_COMPARISON_VARIABLES: Tuple[FitVariable, ...] = (
    *ONE_PHOTON_DISPLAY_VARIABLES,
)

PI0_MASS_WINDOW: Tuple[float, float] = (0.11, 0.16)


def fit_variables_for_category(
    category_key: str,
) -> Tuple[FitVariable, ...]:
    if category_key == "one_photon":
        return ONE_PHOTON_FIT_VARIABLES
    # endif
    if category_key in (
        "two_photon_more_energetic",
        "two_photon_less_energetic",
    ):
        return TWO_PHOTON_FIT_VARIABLES
    # endif
    raise KeyError(f"Unknown photon category: {category_key}")


def display_variables_for_category(
    category_key: str,
) -> Tuple[FitVariable, ...]:
    if category_key == "one_photon":
        return ONE_PHOTON_DISPLAY_VARIABLES
    # endif
    if category_key in (
        "two_photon_more_energetic",
        "two_photon_less_energetic",
    ):
        return TWO_PHOTON_DISPLAY_VARIABLES
    # endif
    raise KeyError(f"Unknown photon category: {category_key}")


def shape_variables_for_category(
    category_key: str,
) -> Tuple[FitVariable, ...]:
    return display_variables_for_category(category_key)


def raw_shape_branch_keys() -> Tuple[str, ...]:
    """Return the union of non-derived branches needed by any category."""
    keys: List[str] = []
    for variable in (
        *ONE_PHOTON_DISPLAY_VARIABLES,
        *TWO_PHOTON_DISPLAY_VARIABLES,
    ):
        if variable.key == "diphoton_mass":
            continue
        # endif
        if variable.key not in keys:
            keys.append(variable.key)
        # endif
    # endfor
    return tuple(keys)

# Shape-comparison diagnostics include the four fit projections above plus
# native proton/photon kinematics stored in the epgamma trees. These additional
# variables do not enter either the independent or shared pi0-fraction fits.


FRACTION_DRIVERS: Tuple[str, ...] = (
    "theta_gamma_gamma",
)

TOPOLOGY_INFO: Mapping[str, Tuple[str, int, int]] = {
    "CD_FT": ("(CD, FT)", 2, 0),
    "CD_FD": ("(CD, FD)", 2, 1),
    "FD_FD": ("(FD, FD)", 1, 1),
}


@dataclass
class OpportunityRecords:
    runnum: np.ndarray
    eventnum: np.ndarray
    e_p: np.ndarray
    e_theta: np.ndarray
    e_phi: np.ndarray
    p_p: np.ndarray
    p_theta: np.ndarray
    p_phi: np.ndarray
    tag_E: np.ndarray
    tag_theta: np.ndarray
    tag_phi: np.ndarray
    tag_detector: np.ndarray
    proton_detector: np.ndarray
    probe_E: np.ndarray
    probe_p: np.ndarray
    probe_theta: np.ndarray
    probe_phi: np.ndarray
    probe_m2: np.ndarray
    probe_E_minus_p: np.ndarray
    probe_detector: np.ndarray
    probe_sector: np.ndarray
    topology_code: np.ndarray
    Delta_phi: np.ndarray
    theta_cm: np.ndarray
    theta_gamma_gamma: np.ndarray
    pTmiss: np.ndarray
    Emiss2: np.ndarray
    Mx2: np.ndarray
    Mx2_2: np.ndarray

    def size(self) -> int:
        return int(self.tag_E.size)


@dataclass
class PhotonCandidateRecords:
    runnum: np.ndarray
    eventnum: np.ndarray
    e_p: np.ndarray
    e_theta: np.ndarray
    e_phi: np.ndarray
    p_p: np.ndarray
    p_theta: np.ndarray
    p_phi: np.ndarray
    photon_E: np.ndarray
    photon_theta: np.ndarray
    photon_phi: np.ndarray
    photon_detector: np.ndarray

    def size(self) -> int:
        return int(self.photon_E.size)


@dataclass
class IdentityCatalog:
    runnum: np.ndarray
    eventnum: np.ndarray
    e_p: np.ndarray
    e_theta: np.ndarray
    e_phi: np.ndarray
    p_p: np.ndarray
    p_theta: np.ndarray
    p_phi: np.ndarray

    def size(self) -> int:
        return int(self.e_p.size)


@dataclass
class EpggRecords:
    runnum: np.ndarray
    eventnum: np.ndarray
    e_p: np.ndarray
    e_theta: np.ndarray
    e_phi: np.ndarray
    p_p: np.ndarray
    p_theta: np.ndarray
    p_phi: np.ndarray
    g1_E: np.ndarray
    g1_theta: np.ndarray
    g1_phi: np.ndarray
    g2_E: np.ndarray
    g2_theta: np.ndarray
    g2_phi: np.ndarray
    solution_valid: np.ndarray
    solution_closure: np.ndarray
    detector1: np.ndarray
    detector2: np.ndarray

    def size(self) -> int:
        return int(self.g1_E.shape[0])


@dataclass
class PartnerMatchResult:
    matched: np.ndarray
    catalog_confirmed: np.ndarray
    partner_index: np.ndarray
    probe_angle_deg: np.ndarray
    probe_relative_E: np.ndarray
    pair_mass_GeV: np.ndarray
    score: np.ndarray
    partner_E: np.ndarray
    tag_angle_deg: np.ndarray
    tag_relative_E: np.ndarray
    solution_closure: np.ndarray
    summary: Dict[str, object]


@dataclass
class EfficiencyRow:
    period: str
    period_label: str
    detector: str
    sector: int
    observed_data: float
    expected_data: float
    expected_data_err: float
    found_data: float
    raw_efficiency_data: float
    efficiency_data: float
    efficiency_data_err: float
    expected_mc: float
    found_mc: float
    raw_efficiency_mc: float
    efficiency_mc: float
    efficiency_mc_err: float
    raw_scale_factor: float
    scale_factor: float
    scale_factor_err: float
    population_valid: bool
    population_failure_reasons: str
    fit_success: bool
    fit_quality_warning: bool
    fit_warning_reasons: str
    fit_fraction_pi0: float
    fit_fraction_pi0_err: float
    fit_deviance: float
    fit_ndf: int
    fit_reduced_deviance: float




def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def elapsed_seconds(start_time: float) -> float:
    """Return elapsed wall-clock seconds from a perf_counter timestamp."""
    return float(time.perf_counter() - start_time)


def maybe_log_timing(
    args: argparse.Namespace,
    label: str,
    seconds: float,
) -> None:
    """Emit a timing line only when explicitly requested."""
    if getattr(args, "print_stage_timings", False):
        log(f"TIMING {label}: {seconds:.3f} s")
    # endif




def selected_periods(args: argparse.Namespace) -> List[PeriodConfig]:
    if not args.period:
        return list(PERIODS)
    # endif
    requested = set(args.period)
    return [period for period in PERIODS if period.key in requested]


def require_tree(path: str) -> Tuple[int, List[str]]:
    if not Path(path).is_file():
        raise FileNotFoundError(f"Missing input ROOT file: {path}")
    # endif
    with uproot.open(path) as root_file:
        if TREE_NAME not in root_file:
            raise KeyError(f"Tree '{TREE_NAME}' is absent from {path}")
        # endif
        tree = root_file[TREE_NAME]
        return int(tree.num_entries), [str(key) for key in tree.keys()]
    # endwith


def resolve_branches(
    path: str,
    logical_names: Sequence[str],
    optional: Sequence[str] = (),
) -> Dict[str, Optional[str]]:
    _, keys = require_tree(path)
    available = set(keys)
    optional_set = set(optional)
    resolved: Dict[str, Optional[str]] = {}
    missing: List[str] = []

    for logical in logical_names:
        aliases = ALIASES.get(logical, (logical,))
        branch = next((candidate for candidate in aliases if candidate in available), None)
        resolved[logical] = branch
        if branch is None and logical not in optional_set:
            missing.append(f"{logical}: tried {', '.join(aliases)}")
        # endif
    # endfor

    if missing:
        raise KeyError(
            f"Missing required branches in {path}:\n  " + "\n  ".join(missing)
        )
    # endif
    return resolved


def finite_array(
    arrays: Mapping[str, np.ndarray],
    branch: Optional[str],
    default: float = math.nan,
    dtype=np.float64,
) -> np.ndarray:
    n = len(next(iter(arrays.values())))
    if branch is None:
        return np.full(n, default, dtype=dtype)
    # endif
    return np.asarray(arrays[branch], dtype=dtype)


def spherical_to_cartesian(
    momentum: np.ndarray,
    theta: np.ndarray,
    phi: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    sin_theta = np.sin(theta)
    return (
        momentum * sin_theta * np.cos(phi),
        momentum * sin_theta * np.sin(phi),
        momentum * np.cos(theta),
    )


def massive_four_vector(
    momentum: np.ndarray,
    theta: np.ndarray,
    phi: np.ndarray,
    mass: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    px, py, pz = spherical_to_cartesian(momentum, theta, phi)
    energy = np.sqrt(np.maximum(momentum * momentum + mass * mass, 0.0))
    return energy, px, py, pz


def photon_four_vector(
    energy: np.ndarray,
    theta: np.ndarray,
    phi: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    px, py, pz = spherical_to_cartesian(energy, theta, phi)
    return energy, px, py, pz


def reconstruct_probe(
    beam_energy: float,
    e_p: np.ndarray,
    e_theta: np.ndarray,
    e_phi: np.ndarray,
    p_p: np.ndarray,
    p_theta: np.ndarray,
    p_phi: np.ndarray,
    tag_E: np.ndarray,
    tag_theta: np.ndarray,
    tag_phi: np.ndarray,
) -> Dict[str, np.ndarray]:
    e_E, e_px, e_py, e_pz = massive_four_vector(
        e_p, e_theta, e_phi, ELECTRON_MASS_GEV
    )
    p_E, p_px, p_py, p_pz = massive_four_vector(
        p_p, p_theta, p_phi, PROTON_MASS_GEV
    )
    g_E, g_px, g_py, g_pz = photon_four_vector(tag_E, tag_theta, tag_phi)

    probe_E = beam_energy + PROTON_MASS_GEV - e_E - p_E - g_E
    probe_px = -e_px - p_px - g_px
    probe_py = -e_py - p_py - g_py
    probe_pz = beam_energy - e_pz - p_pz - g_pz
    probe_p = np.sqrt(np.maximum(probe_px**2 + probe_py**2 + probe_pz**2, 0.0))
    probe_pt = np.sqrt(np.maximum(probe_px**2 + probe_py**2, 0.0))
    probe_theta = np.arctan2(probe_pt, probe_pz)
    probe_phi = np.mod(np.arctan2(probe_py, probe_px), 2.0 * math.pi)
    return {
        "E": probe_E,
        "p": probe_p,
        "theta": probe_theta,
        "phi": probe_phi,
        "m2": probe_E**2 - probe_p**2,
        "E_minus_p": probe_E - probe_p,
    }


def fd_sector_from_phi(phi_rad: np.ndarray) -> np.ndarray:
    phi_deg = np.mod(np.degrees(phi_rad), 360.0)
    sector = np.full(phi_deg.shape, -1, dtype=np.int16)
    sector[(phi_deg >= 330.0) | (phi_deg < 30.0)] = 1
    sector[(phi_deg >= 30.0) & (phi_deg < 90.0)] = 2
    sector[(phi_deg >= 90.0) & (phi_deg < 150.0)] = 3
    sector[(phi_deg >= 150.0) & (phi_deg < 210.0)] = 4
    sector[(phi_deg >= 210.0) & (phi_deg < 270.0)] = 5
    sector[(phi_deg >= 270.0) & (phi_deg < 330.0)] = 6
    return sector


def classify_probe(
    theta_deg: np.ndarray,
    phi_rad: np.ndarray,
    args: argparse.Namespace,
) -> Tuple[np.ndarray, np.ndarray]:
    detector = np.full(theta_deg.shape, -1, dtype=np.int16)
    detector[
        (theta_deg >= args.ft_theta_min)
        & (theta_deg < args.ft_theta_max)
    ] = 0
    detector[
        (theta_deg >= args.fd_theta_min)
        & (theta_deg < args.fd_theta_max)
    ] = 1
    sector = fd_sector_from_phi(phi_rad)
    sector[detector != 1] = 0
    return detector, sector


def topology_code(
    proton_detector: np.ndarray,
    tag_detector: np.ndarray,
) -> np.ndarray:
    code = np.full(proton_detector.shape, -1, dtype=np.int8)
    code[(proton_detector == 2) & (tag_detector == 0)] = 0  # CD_FT
    code[(proton_detector == 2) & (tag_detector == 1)] = 1  # CD_FD
    code[(proton_detector == 1) & (tag_detector == 1)] = 2  # FD_FD
    return code


def topology_key(code: int) -> str:
    return {0: "CD_FT", 1: "CD_FD", 2: "FD_FD"}[int(code)]


def opening_angle_deg(
    theta1: np.ndarray,
    phi1: np.ndarray,
    theta2: np.ndarray,
    phi2: np.ndarray,
) -> np.ndarray:
    cosine = (
        np.sin(theta1) * np.sin(theta2) * np.cos(phi1 - phi2)
        + np.cos(theta1) * np.cos(theta2)
    )
    return np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))


def opening_angle_deg_scalar(
    theta1: float,
    phi1: float,
    theta2: float,
    phi2: float,
) -> float:
    """Allocation-free scalar opening angle in degrees."""
    cosine = (
        math.cos(theta1) * math.cos(theta2)
        + math.sin(theta1)
        * math.sin(theta2)
        * math.cos(phi1 - phi2)
    )
    return math.degrees(
        math.acos(min(1.0, max(-1.0, cosine)))
    )


def concat(parts: Sequence[np.ndarray], dtype=np.float64) -> np.ndarray:
    if not parts:
        return np.asarray([], dtype=dtype)
    # endif
    return np.concatenate(parts).astype(dtype, copy=False)


def subset_opportunities(records: OpportunityRecords, mask: np.ndarray) -> OpportunityRecords:
    return OpportunityRecords(
        **{
            name: np.asarray(getattr(records, name))[mask]
            for name in records.__dataclass_fields__
        }
    )


def exact_duplicate_keep_mask(
    arrays: Sequence[np.ndarray],
    decimals: int = 10,
) -> np.ndarray:
    if not arrays or len(arrays[0]) == 0:
        return np.ones(0, dtype=bool)
    # endif
    rounded = [
        np.round(np.asarray(values, dtype=np.float64), decimals)
        for values in arrays
    ]
    matrix = np.column_stack(rounded)
    _, first = np.unique(matrix, axis=0, return_index=True)
    keep = np.zeros(matrix.shape[0], dtype=bool)
    keep[np.sort(first)] = True
    return keep



CUT_STAGE_ORDER: Tuple[str, ...] = (
    "tree_entries",
    "photon_multiplicity_lt_3",
    "finite_kinematics",
    "Q2",
    "W",
    "y",
    "z_tag",
    "minus_t1_tag",
    "open_angle_ep2_tag",
    "fiducial_status",
    "tag_energy",
    "probe_finite",
    "probe_energy_min",
    "probe_energy_max",
    "probe_m2",
    "probe_E_minus_p",
    "probe_acceptance",
    "supported_topology",
    "selected_before_deduplication",
    "selected_after_deduplication",
)


def cut_stage_definitions(
    resolved: Mapping[str, Optional[str]],
    args: argparse.Namespace,
) -> Dict[str, Dict[str, object]]:
    """Return human-readable definitions for every sequential opportunity cut."""
    return {
        "tree_entries": {
            "label": "Input ROOT-tree entries",
            "condition": "all entries",
            "branch": None,
            "applied": True,
            "scope": "input",
        },
        "photon_multiplicity_lt_3": {
            "label": "Event photon-entry multiplicity",
            "condition": (
                "retain only events with one or two reconstructed-photon "
                "entries; reject multiplicity >= 3"
            ),
            "branch": (
                "runnum,evnum for data; rounded reconstructed e/p "
                "kinematics for MC"
            ),
            "applied": True,
            "scope": "hard event-level cut",
        },
        "finite_kinematics": {
            "label": "Finite positive reconstructed e, p, and tag photon",
            "condition": (
                "finite e/p/gamma momentum and angles; "
                "p_e > 0, p_p > 0, E_tag > 0"
            ),
            "branch": "e_p,e_theta,e_phi,p1_p,p1_theta,p1_phi,p2_p,p2_theta,p2_phi",
            "applied": True,
            "scope": "common event/tag",
        },
        "Q2": {
            "label": "DIS virtuality",
            "condition": f"Q2 > {args.Q2_min:g} GeV^2",
            "branch": resolved.get("Q2"),
            "applied": resolved.get("Q2") is not None,
            "scope": "common event",
        },
        "W": {
            "label": "DIS invariant mass",
            "condition": f"W > {args.W_min:g} GeV",
            "branch": resolved.get("W"),
            "applied": resolved.get("W") is not None,
            "scope": "common event",
        },
        "y": {
            "label": "DIS inelasticity",
            "condition": f"y < {args.y_max:g}",
            "branch": resolved.get("y"),
            "applied": resolved.get("y") is not None,
            "scope": "common event",
        },
        "z_tag": {
            "label": "Tag-photon z requirement",
            "condition": (
                f"NOT APPLIED: z(tag photon) > {args.z_min:g} "
                "was removed from the nominal tag-and-probe selection"
            ),
            "branch": resolved.get("z"),
            "applied": False,
            "scope": "REMOVED — tag-dependent DVCS cut",
        },
        "minus_t1_tag": {
            "label": "Tag-photon -t1 requirement",
            "condition": (
                f"NOT APPLIED: -t1(tag photon) < {args.minus_t_max:g} GeV^2 "
                "was removed from the nominal tag-and-probe selection"
            ),
            "branch": resolved.get("t1"),
            "applied": False,
            "scope": "REMOVED — tag-dependent DVCS cut",
        },
        "open_angle_ep2_tag": {
            "label": "Electron–tag-photon opening angle",
            "condition": (
                f"NOT APPLIED: open_angle(e, tag photon) > "
                f"{args.open_angle_min_deg:g} deg was removed from the "
                "nominal tag-and-probe selection"
            ),
            "branch": resolved.get("open_angle_ep2"),
            "applied": False,
            "scope": "REMOVED — tag-dependent DVCS cut",
        },
        "fiducial_status": {
            "label": "Combined fiducial status",
            "condition": "fiducial_status == 111",
            "branch": resolved.get("fiducial_status"),
            "applied": (
                bool(args.require_fiducial_status_111)
                and resolved.get("fiducial_status") is not None
            ),
            "scope": "optional reconstructed-particle fiducial",
        },
        "tag_energy": {
            "label": "Reconstructed tag-photon energy",
            "condition": f"E_tag >= {args.tag_E_min:g} GeV",
            "branch": resolved.get("tag_E"),
            "applied": True,
            "scope": "tag",
        },
        "probe_finite": {
            "label": "Finite physical missing-photon four-vector",
            "condition": (
                "finite E_X, |p_X|, theta_X, phi_X, m_X^2, E_X-|p_X|; "
                "E_X > 0 and |p_X| > 0"
            ),
            "branch": "derived from beam + target - e - p - tag",
            "applied": True,
            "scope": "predicted probe",
        },
        "probe_energy_min": {
            "label": "Predicted probe minimum energy",
            "condition": f"E_X >= {args.probe_E_min:g} GeV",
            "branch": "derived E_X",
            "applied": True,
            "scope": "predicted probe",
        },
        "probe_energy_max": {
            "label": "Predicted probe maximum energy",
            "condition": f"E_X < {args.probe_E_max:g} GeV",
            "branch": "derived E_X",
            "applied": True,
            "scope": "predicted probe",
        },
        "probe_m2": {
            "label": "Predicted probe mass-squared consistency",
            "condition": f"|m_X^2| < {args.probe_m2_abs_max:g} GeV^2",
            "branch": "derived m_X^2",
            "applied": True,
            "scope": "predicted probe",
        },
        "probe_E_minus_p": {
            "label": "Predicted probe massless-energy consistency",
            "condition": (
                f"|E_X - |p_X|| < "
                f"{args.probe_E_minus_p_abs_max:g} GeV"
            ),
            "branch": "derived E_X-|p_X|",
            "applied": True,
            "scope": "predicted probe",
        },
        "probe_acceptance": {
            "label": "Predicted probe points into FT or FD",
            "condition": (
                f"FT: {args.ft_theta_min:g} <= theta_X < "
                f"{args.ft_theta_max:g} deg; FD: "
                f"{args.fd_theta_min:g} <= theta_X < "
                f"{args.fd_theta_max:g} deg"
            ),
            "branch": "derived theta_X and phi_X",
            "applied": True,
            "scope": "predicted probe detector",
        },
        "supported_topology": {
            "label": "Supported reconstructed proton/tag topology",
            "condition": "proton and tag detector codes map to a supported topology",
            "branch": "detector1, detector2",
            "applied": True,
            "scope": "bookkeeping only; proton topologies later combined",
        },
        "selected_before_deduplication": {
            "label": "Selected directed opportunities before MC deduplication",
            "condition": "all preceding requirements",
            "branch": None,
            "applied": True,
            "scope": "selected sample",
        },
        "selected_after_deduplication": {
            "label": "Selected directed opportunities after MC deduplication",
            "condition": "exact duplicate removal for MC; unchanged for data",
            "branch": None,
            "applied": True,
            "scope": "final fit denominator population",
        },
    }


def stage_mask_with_optional_branch(
    previous: np.ndarray,
    arrays: Mapping[str, np.ndarray],
    branch: Optional[str],
    predicate,
) -> np.ndarray:
    """Apply one branch cut sequentially; absent optional branches leave the mask unchanged."""
    if branch is None:
        return previous.copy()
    # endif
    values = finite_array(arrays, branch)
    return previous & np.isfinite(values) & predicate(values)


def finalize_cutflow_report(cutflow: Dict[str, object]) -> Dict[str, object]:
    """Add incremental and total survival rates to the raw cumulative counters."""
    definitions = cutflow["stage_definitions"]
    counts = cutflow["stage_counts"]
    rows: List[Dict[str, object]] = []
    total = int(counts.get("tree_entries", 0))
    previous_count = total

    for index, stage in enumerate(CUT_STAGE_ORDER):
        count = int(counts.get(stage, 0))
        definition = definitions[stage]
        applied = bool(definition["applied"])
        reference = total if index == 0 else previous_count
        incremental = (
            float(count / reference)
            if reference > 0
            else math.nan
        )
        total_survival = (
            float(count / total)
            if total > 0
            else math.nan
        )
        rows.append(
            {
                "order": index,
                "stage": stage,
                "label": definition["label"],
                "condition": definition["condition"],
                "branch": definition["branch"],
                "scope": definition["scope"],
                "applied": applied,
                "count": count,
                "rejected_at_stage": (
                    0 if index == 0 else max(previous_count - count, 0)
                ),
                "incremental_survival": incremental,
                "incremental_rejection": (
                    1.0 - incremental
                    if math.isfinite(incremental)
                    else math.nan
                ),
                "total_survival": total_survival,
            }
        )
        previous_count = count
    # endfor

    cutflow["stages"] = rows
    cutflow["largest_incremental_loss"] = max(
        rows[1:],
        key=lambda row: (
            row["incremental_rejection"]
            if math.isfinite(row["incremental_rejection"])
            else -1.0
        ),
    ) if len(rows) > 1 else None
    return cutflow


def common_global_mask(
    arrays: Mapping[str, np.ndarray],
    resolved: Mapping[str, Optional[str]],
    args: argparse.Namespace,
) -> np.ndarray:
    n = len(next(iter(arrays.values())))
    mask = np.ones(n, dtype=bool)

    def lower(key: str, threshold: float) -> None:
        nonlocal mask
        branch = resolved.get(key)
        if branch is not None:
            values = finite_array(arrays, branch)
            mask &= np.isfinite(values) & (values > threshold)
        # endif

    def upper(key: str, threshold: float) -> None:
        nonlocal mask
        branch = resolved.get(key)
        if branch is not None:
            values = finite_array(arrays, branch)
            mask &= np.isfinite(values) & (values < threshold)
        # endif

    lower("Q2", args.Q2_min)
    lower("W", args.W_min)
    upper("y", args.y_max)
    lower("z", args.z_min)

    # Tag-dependent z, t1, and electron-tag opening-angle cuts are
    # intentionally omitted from the nominal opportunity definition.

    if args.require_fiducial_status_111 and resolved.get("fiducial_status") is not None:
        status = finite_array(arrays, resolved["fiducial_status"])
        mask &= np.isfinite(status) & (status == 111)
    # endif
    return mask


def empty_opportunity_store() -> Dict[str, List[np.ndarray]]:
    return {name: [] for name in OpportunityRecords.__dataclass_fields__}


def append_opportunity_store(
    store: Dict[str, List[np.ndarray]],
    values: Mapping[str, np.ndarray],
    mask: np.ndarray,
) -> None:
    for name in store:
        store[name].append(np.asarray(values[name])[mask])
    # endfor


def finalize_opportunity_store(store: Mapping[str, Sequence[np.ndarray]]) -> OpportunityRecords:
    integer_fields = {
        "runnum", "eventnum", "tag_detector", "proton_detector",
        "probe_detector", "probe_sector", "topology_code",
    }
    payload = {}
    for name, parts in store.items():
        dtype = np.int64 if name in {"runnum", "eventnum"} else (
            np.int16 if name in integer_fields else np.float64
        )
        payload[name] = concat(parts, dtype=dtype)
    # endfor
    return OpportunityRecords(**payload)


def empty_candidate_store() -> Dict[str, List[np.ndarray]]:
    return {name: [] for name in PhotonCandidateRecords.__dataclass_fields__}


def append_candidate_store(
    store: Dict[str, List[np.ndarray]],
    values: Mapping[str, np.ndarray],
    mask: np.ndarray,
) -> None:
    for name in store:
        store[name].append(np.asarray(values[name])[mask])
    # endfor


def finalize_candidate_store(
    store: Mapping[str, Sequence[np.ndarray]],
) -> PhotonCandidateRecords:
    integer_fields = {"runnum", "eventnum", "photon_detector"}
    payload = {}
    for name, parts in store.items():
        dtype = np.int64 if name in {"runnum", "eventnum"} else (
            np.int16 if name in integer_fields else np.float64
        )
        payload[name] = concat(parts, dtype=dtype)
    # endfor
    return PhotonCandidateRecords(**payload)


def subset_candidates(
    records: PhotonCandidateRecords,
    mask: np.ndarray,
) -> PhotonCandidateRecords:
    return PhotonCandidateRecords(
        **{
            name: np.asarray(getattr(records, name))[mask]
            for name in records.__dataclass_fields__
        }
    )


def read_opportunities(
    path: str,
    beam_energy: float,
    role: str,
    args: argparse.Namespace,
    deduplicate: bool,
    collect_candidates: bool = True,
) -> Tuple[OpportunityRecords, PhotonCandidateRecords, Dict[str, object]]:
    logical = (
        "runnum", "eventnum",
        "e_p", "e_theta", "e_phi",
        "p_p", "p_theta", "p_phi",
        "tag_E", "tag_theta", "tag_phi",
        "tag_detector", "proton_detector",
        "fiducial_status", "Q2", "W", "y", "z", "t1", "open_angle_ep2",
        "Delta_phi", "theta_cm", "theta_gamma_gamma", "pTmiss",
        "Emiss2", "Mx2", "Mx2_1", "Mx2_2",
    )
    optional = (
        "runnum", "eventnum", "tag_detector", "proton_detector",
        "fiducial_status", "Q2", "W", "y", "z", "t1", "open_angle_ep2",
        "theta_cm", "Mx2_1",
    )
    resolved = resolve_branches(path, logical, optional=optional)
    expressions = sorted({branch for branch in resolved.values() if branch is not None})

    tree_entries, _available_keys = require_tree(path)
    entry_stop = (
        min(tree_entries, int(args.max_events))
        if args.max_events is not None
        else tree_entries
    )
    mc_identity_resolved = {
        logical_name: resolved[logical_name]
        for logical_name in identity_logical_names("mc")
    }
    (
        _all_mc_signatures,
        one_mc_signatures,
        two_mc_signatures,
        _three_plus_mc_signatures,
        multiplicity_audit,
    ) = scan_event_multiplicity(
        path,
        "mc",
        mc_identity_resolved,
        args,
        entry_stop,
    )

    store = empty_opportunity_store()
    candidate_store = (
        empty_candidate_store()
        if collect_candidates
        else None
    )
    stage_definitions = cut_stage_definitions(resolved, args)
    cutflow: Dict[str, object] = {
        "role": role,
        "path": path,
        "stage_definitions": stage_definitions,
        "stage_counts": {stage: 0 for stage in CUT_STAGE_ORDER},
        "detector_counts": {
            stage: {"FT": 0, "FD": 0}
            for stage in (
                "probe_finite",
                "probe_energy_min",
                "probe_energy_max",
                "probe_m2",
                "probe_E_minus_p",
                "probe_acceptance",
                "supported_topology",
                "selected_before_deduplication",
                "selected_after_deduplication",
            )
        },
        "partner_candidates_finite_and_E_above_threshold": 0,
        "multiplicity_audit": {
            "condition": (
                "all photon-entry multiplicities are retained; multiplicity "
                "is reported for diagnostics only"
            ),
            **multiplicity_audit,
        },
    }
    stage_counts = cutflow["stage_counts"]

    seen = 0
    log(f"Reading {role} from {path}")
    for arrays in uproot.iterate(
        f"{path}:{TREE_NAME}",
        expressions=expressions,
        step_size=args.step_size,
        entry_stop=entry_stop,
        library="np",
    ):
        n = len(next(iter(arrays.values())))
        if args.max_events is not None and seen >= args.max_events:
            break
        # endif
        if args.max_events is not None and seen + n > args.max_events:
            n = args.max_events - seen
            arrays = {key: value[:n] for key, value in arrays.items()}
        # endif
        seen += n
        stage_counts["tree_entries"] += int(n)

        def arr(key: str, default: float = math.nan, dtype=np.float64) -> np.ndarray:
            return finite_array(arrays, resolved.get(key), default=default, dtype=dtype)

        e_p, e_theta, e_phi = arr("e_p"), arr("e_theta"), arr("e_phi")
        p_p, p_theta, p_phi = arr("p_p"), arr("p_theta"), arr("p_phi")
        tag_E = arr("tag_E")
        tag_theta = arr("tag_theta")
        tag_phi = np.mod(arr("tag_phi"), 2.0 * math.pi)
        tag_detector = arr("tag_detector", default=-1, dtype=np.int16)
        proton_detector = arr("proton_detector", default=-1, dtype=np.int16)

        identity_arrays = {
            "e_p": e_p,
            "e_theta": e_theta,
            "e_phi": e_phi,
            "p_p": p_p,
            "p_theta": p_theta,
            "p_phi": p_phi,
        }
        (
            one_multiplicity_keep,
            two_multiplicity_keep,
            _three_plus_multiplicity_keep,
            valid_identity,
        ) = multiplicity_entry_masks(
            identity_arrays,
            "mc",
            one_mc_signatures,
            two_mc_signatures,
            _three_plus_mc_signatures,
            args,
        )
        multiplicity_keep = (
            one_multiplicity_keep | two_multiplicity_keep
        )
        stage_counts["photon_multiplicity_lt_3"] += int(
            np.count_nonzero(multiplicity_keep)
        )

        finite = multiplicity_keep & (
            np.isfinite(e_p) & np.isfinite(e_theta) & np.isfinite(e_phi)
            & np.isfinite(p_p) & np.isfinite(p_theta) & np.isfinite(p_phi)
            & np.isfinite(tag_E) & np.isfinite(tag_theta) & np.isfinite(tag_phi)
            & (e_p > 0.0) & (p_p > 0.0) & (tag_E > 0.0)
        )
        stage_counts["finite_kinematics"] += int(np.count_nonzero(finite))

        mask_Q2 = stage_mask_with_optional_branch(
            finite, arrays, resolved.get("Q2"),
            lambda values: values > args.Q2_min,
        )
        stage_counts["Q2"] += int(np.count_nonzero(mask_Q2))

        mask_W = stage_mask_with_optional_branch(
            mask_Q2, arrays, resolved.get("W"),
            lambda values: values > args.W_min,
        )
        stage_counts["W"] += int(np.count_nonzero(mask_W))

        mask_y = stage_mask_with_optional_branch(
            mask_W, arrays, resolved.get("y"),
            lambda values: values < args.y_max,
        )
        stage_counts["y"] += int(np.count_nonzero(mask_y))

        # These three DVCS-style cuts are tag dependent. The reconstructed
        # tag may be the low-energy pi0 daughter, whereas the efficiency probe
        # is the predicted high-energy photon. They are therefore retained as
        # explicit audit stages but intentionally do not alter the mask.
        mask_z = mask_y.copy()
        stage_counts["z_tag"] += int(np.count_nonzero(mask_z))

        mask_t1 = mask_z.copy()
        stage_counts["minus_t1_tag"] += int(np.count_nonzero(mask_t1))

        mask_open = mask_t1.copy()
        stage_counts["open_angle_ep2_tag"] += int(np.count_nonzero(mask_open))

        if (
            args.require_fiducial_status_111
            and resolved.get("fiducial_status") is not None
        ):
            status = finite_array(arrays, resolved["fiducial_status"])
            mask_fiducial = (
                mask_open
                & np.isfinite(status)
                & (status == 111)
            )
        else:
            mask_fiducial = mask_open.copy()
        # endif
        stage_counts["fiducial_status"] += int(
            np.count_nonzero(mask_fiducial)
        )

        tag_mask = mask_fiducial & (tag_E >= args.tag_E_min)
        stage_counts["tag_energy"] += int(np.count_nonzero(tag_mask))

        # The partner pool intentionally uses every finite reconstructed photon
        # above threshold, independent of the tag-specific opportunity cuts.
        candidate_mask = finite & (tag_E >= args.tag_E_min)
        candidate_values = {
            "runnum": arr("runnum", default=-1, dtype=np.int64),
            "eventnum": arr("eventnum", default=-1, dtype=np.int64),
            "e_p": e_p, "e_theta": e_theta, "e_phi": e_phi,
            "p_p": p_p, "p_theta": p_theta, "p_phi": p_phi,
            "photon_E": tag_E,
            "photon_theta": tag_theta,
            "photon_phi": tag_phi,
            "photon_detector": tag_detector,
        }
        if collect_candidates:
            append_candidate_store(
                candidate_store,
                candidate_values,
                candidate_mask,
            )
            cutflow[
                "partner_candidates_finite_and_E_above_threshold"
            ] += int(np.count_nonzero(candidate_mask))
        # endif

        probe = reconstruct_probe(
            beam_energy,
            e_p, e_theta, e_phi,
            p_p, p_theta, p_phi,
            tag_E, tag_theta, tag_phi,
        )
        probe_theta_deg = np.degrees(probe["theta"])
        probe_detector, probe_sector = classify_probe(
            probe_theta_deg, probe["phi"], args
        )

        def accumulate_detector_stage(stage: str, mask: np.ndarray) -> None:
            cutflow["detector_counts"][stage]["FT"] += int(
                np.count_nonzero(mask & (probe_detector == 0))
            )
            cutflow["detector_counts"][stage]["FD"] += int(
                np.count_nonzero(mask & (probe_detector == 1))
            )
        # enddef

        probe_finite = (
            tag_mask
            & np.isfinite(probe["E"]) & np.isfinite(probe["p"])
            & np.isfinite(probe["theta"]) & np.isfinite(probe["phi"])
            & np.isfinite(probe["m2"]) & np.isfinite(probe["E_minus_p"])
            & (probe["E"] > 0.0) & (probe["p"] > 0.0)
        )
        stage_counts["probe_finite"] += int(np.count_nonzero(probe_finite))
        accumulate_detector_stage("probe_finite", probe_finite)

        probe_E_min = probe_finite & (probe["E"] >= args.probe_E_min)
        stage_counts["probe_energy_min"] += int(
            np.count_nonzero(probe_E_min)
        )
        accumulate_detector_stage("probe_energy_min", probe_E_min)

        probe_E_window = probe_E_min & (probe["E"] < args.probe_E_max)
        stage_counts["probe_energy_max"] += int(
            np.count_nonzero(probe_E_window)
        )
        accumulate_detector_stage("probe_energy_max", probe_E_window)

        probe_m2 = (
            probe_E_window
            & (np.abs(probe["m2"]) < args.probe_m2_abs_max)
        )
        stage_counts["probe_m2"] += int(np.count_nonzero(probe_m2))
        accumulate_detector_stage("probe_m2", probe_m2)

        photon_like = (
            probe_m2
            & (
                np.abs(probe["E_minus_p"])
                < args.probe_E_minus_p_abs_max
            )
        )
        stage_counts["probe_E_minus_p"] += int(
            np.count_nonzero(photon_like)
        )
        accumulate_detector_stage("probe_E_minus_p", photon_like)

        accepted = photon_like & (probe_detector >= 0)
        stage_counts["probe_acceptance"] += int(np.count_nonzero(accepted))
        accumulate_detector_stage("probe_acceptance", accepted)

        topo = topology_code(proton_detector, tag_detector)
        selected = accepted & (topo >= 0)
        stage_counts["supported_topology"] += int(np.count_nonzero(selected))
        stage_counts["selected_before_deduplication"] += int(
            np.count_nonzero(selected)
        )
        accumulate_detector_stage("supported_topology", selected)
        accumulate_detector_stage("selected_before_deduplication", selected)

        values = {
            "runnum": arr("runnum", default=-1, dtype=np.int64),
            "eventnum": arr("eventnum", default=-1, dtype=np.int64),
            "e_p": e_p, "e_theta": e_theta, "e_phi": e_phi,
            "p_p": p_p, "p_theta": p_theta, "p_phi": p_phi,
            "tag_E": tag_E, "tag_theta": tag_theta, "tag_phi": tag_phi,
            "tag_detector": tag_detector,
            "proton_detector": proton_detector,
            "probe_E": probe["E"], "probe_p": probe["p"],
            "probe_theta": probe["theta"], "probe_phi": probe["phi"],
            "probe_m2": probe["m2"],
            "probe_E_minus_p": probe["E_minus_p"],
            "probe_detector": probe_detector,
            "probe_sector": probe_sector,
            "topology_code": topo,
            "Delta_phi": arr("Delta_phi"),
            "theta_cm": arr("theta_cm"),
            "theta_gamma_gamma": arr("theta_gamma_gamma"),
            "pTmiss": arr("pTmiss"),
            "Emiss2": arr("Emiss2"),
            "Mx2": arr("Mx2"),
            "Mx2_2": arr("Mx2_2"),
        }
        append_opportunity_store(store, values, selected)
    # endfor

    records = finalize_opportunity_store(store)

    if collect_candidates:
        candidates = finalize_candidate_store(candidate_store)
        if deduplicate and candidates.size() > 0:
            candidate_keep = exact_duplicate_keep_mask(
                [
                    candidates.e_p,
                    candidates.e_theta,
                    candidates.e_phi,
                    candidates.p_p,
                    candidates.p_theta,
                    candidates.p_phi,
                    candidates.photon_E,
                    candidates.photon_theta,
                    candidates.photon_phi,
                ],
                decimals=args.mc_signature_decimals,
            )
            candidates = subset_candidates(
                candidates,
                candidate_keep,
            )
        # endif
    else:
        candidates = finalize_candidate_store(
            empty_candidate_store()
        )
    # endif
    if deduplicate and records.size() > 0:
        keep = exact_duplicate_keep_mask(
            [
                records.e_p, records.e_theta, records.e_phi,
                records.p_p, records.p_theta, records.p_phi,
                records.tag_E, records.tag_theta, records.tag_phi,
            ],
            decimals=args.mc_signature_decimals,
        )
        records = subset_opportunities(records, keep)
    # endif
    stage_counts["selected_after_deduplication"] = records.size()
    cutflow["detector_counts"]["selected_after_deduplication"]["FT"] = int(
        np.count_nonzero(records.probe_detector == 0)
    )
    cutflow["detector_counts"]["selected_after_deduplication"]["FD"] = int(
        np.count_nonzero(records.probe_detector == 1)
    )
    cutflow["branch_mapping"] = resolved
    cutflow["partner_candidate_records"] = candidates.size()
    cutflow = finalize_cutflow_report(cutflow)
    return records, candidates, cutflow


def unit_vector(theta: float, phi: float) -> np.ndarray:
    return np.asarray(
        [
            math.sin(theta) * math.cos(phi),
            math.sin(theta) * math.sin(phi),
            math.cos(theta),
        ],
        dtype=float,
    )

def q_trento_basis(
    beam_energy: float,
    electron_p: float,
    electron_theta: float,
    electron_phi: float,
) -> Optional[Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
    electron_direction = unit_vector(electron_theta, electron_phi)
    electron_vector = electron_p * electron_direction
    q_vector = np.asarray([0.0, 0.0, beam_energy]) - electron_vector
    q_magnitude = float(np.linalg.norm(q_vector))
    if q_magnitude <= 1.0e-12:
        return None
    # endif
    q_hat = q_vector / q_magnitude

    electron_transverse = (
        electron_direction
        - float(np.dot(electron_direction, q_hat)) * q_hat
    )
    transverse_magnitude = float(np.linalg.norm(electron_transverse))
    if transverse_magnitude <= 1.0e-12:
        return None
    # endif

    x_hat = electron_transverse / transverse_magnitude
    y_hat = np.cross(q_hat, x_hat)
    y_hat /= float(np.linalg.norm(y_hat))
    return q_hat, x_hat, y_hat, electron_direction

def solve_q_polar_angles(
    electron_direction: np.ndarray,
    q_hat: np.ndarray,
    x_hat: np.ndarray,
    trento_phi: float,
    opening_angle: float,
) -> List[float]:
    """
    Solve for the photon polar angle about q.

    With n_gamma = sin(theta_q)(cos(phi)x + sin(phi)y) + cos(theta_q)q,
    the stored electron-photon opening angle gives

        A sin(theta_q) + B cos(theta_q) = cos(alpha),

    where A = (e_hat dot x_hat) cos(phi) and B = e_hat dot q_hat.
    """
    A = float(np.dot(electron_direction, x_hat)) * math.cos(trento_phi)
    B = float(np.dot(electron_direction, q_hat))
    C = math.cos(opening_angle)
    radius = math.hypot(A, B)
    if radius <= 1.0e-12:
        return []
    # endif

    ratio = C / radius
    if ratio < -1.0 - 1.0e-9 or ratio > 1.0 + 1.0e-9:
        return []
    # endif
    ratio = float(np.clip(ratio, -1.0, 1.0))

    phase = math.atan2(A, B)
    offset = math.acos(ratio)
    candidates: List[float] = []
    for raw in (phase + offset, phase - offset):
        theta = raw % (2.0 * math.pi)
        if theta > math.pi:
            theta = 2.0 * math.pi - theta
        # endif
        if (
            0.0 <= theta <= math.pi
            and all(abs(theta - existing) > 1.0e-9 for existing in candidates)
        ):
            candidates.append(theta)
        # endif
    # endfor
    return candidates

def lab_direction_from_trento(
    theta_q: float,
    trento_phi: float,
    q_hat: np.ndarray,
    x_hat: np.ndarray,
    y_hat: np.ndarray,
) -> np.ndarray:
    direction = (
        math.sin(theta_q)
        * (
            math.cos(trento_phi) * x_hat
            + math.sin(trento_phi) * y_hat
        )
        + math.cos(theta_q) * q_hat
    )
    direction /= float(np.linalg.norm(direction))
    return direction

def direction_to_angles(direction: np.ndarray) -> Tuple[float, float]:
    theta = math.acos(float(np.clip(direction[2], -1.0, 1.0)))
    phi = math.atan2(float(direction[1]), float(direction[0]))
    return theta, float(np.mod(phi, 2.0 * math.pi))

def reconstruct_daughter_solutions(
    beam_energy: float,
    electron_p: float,
    electron_theta: float,
    electron_phi: float,
    pi0_p: float,
    pi0_theta: float,
    pi0_phi: float,
    pi0_mass: float,
    opening1: float,
    opening2: float,
    trento1: float,
    trento2: float,
    args: argparse.Namespace,
) -> List[Tuple[float, float, float, float, float, float, float]]:
    basis = q_trento_basis(
        beam_energy,
        electron_p,
        electron_theta,
        electron_phi,
    )
    if basis is None:
        return []
    # endif
    q_hat, x_hat, y_hat, electron_direction = basis

    theta1_candidates = solve_q_polar_angles(
        electron_direction,
        q_hat,
        x_hat,
        trento1,
        opening1,
    )
    theta2_candidates = solve_q_polar_angles(
        electron_direction,
        q_hat,
        x_hat,
        trento2,
        opening2,
    )
    if not theta1_candidates or not theta2_candidates:
        return []
    # endif

    pi0_direction = unit_vector(pi0_theta, pi0_phi)
    pi0_vector = pi0_p * pi0_direction
    pi0_energy = math.sqrt(max(pi0_p**2 + pi0_mass**2, 0.0))
    target = np.asarray([pi0_energy, *pi0_vector], dtype=float)
    normalization = max(pi0_energy, 1.0e-12)

    solutions = []
    for theta1_q in theta1_candidates:
        direction1 = lab_direction_from_trento(
            theta1_q, trento1, q_hat, x_hat, y_hat
        )
        theta1_lab, phi1_lab = direction_to_angles(direction1)

        for theta2_q in theta2_candidates:
            direction2 = lab_direction_from_trento(
                theta2_q, trento2, q_hat, x_hat, y_hat
            )
            theta2_lab, phi2_lab = direction_to_angles(direction2)

            response = np.asarray(
                [
                    [1.0, 1.0],
                    [direction1[0], direction2[0]],
                    [direction1[1], direction2[1]],
                    [direction1[2], direction2[2]],
                ],
                dtype=float,
            )
            energies, _, rank, _ = np.linalg.lstsq(
                response,
                target,
                rcond=None,
            )
            if (
                rank < 2
                or energies[0] <= args.native_minimum_photon_E
                or energies[1] <= args.native_minimum_photon_E
            ):
                continue
            # endif

            closure = float(
                np.linalg.norm(response @ energies - target)
                / normalization
            )
            if closure > args.native_closure_max:
                continue
            # endif

            solutions.append(
                (
                    closure,
                    float(energies[0]),
                    theta1_lab,
                    phi1_lab,
                    float(energies[1]),
                    theta2_lab,
                    phi2_lab,
                )
            )
        # endfor
    # endfor

    solutions.sort(key=lambda item: item[0])
    return solutions[:4]

def read_epgg(
    path: str,
    beam_energy: float,
    mode: str,
    args: argparse.Namespace,
) -> Tuple[EpggRecords, Dict[str, int]]:
    logical = (
        "runnum", "eventnum",
        "e_p", "e_theta", "e_phi",
        "p_p", "p_theta", "p_phi",
        "pi0_p", "pi0_theta", "pi0_phi", "pi0_mass",
        "open1", "open2", "trento1", "trento2",
        "native_gamma_detector1", "native_gamma_detector2",
    )
    optional_identity = (
        ("runnum", "eventnum")
        if "AAOGEN" in mode.upper()
        else ()
    )
    resolved = resolve_branches(path, logical, optional=optional_identity)
    expressions = sorted(
        {branch for branch in resolved.values() if branch is not None}
    )

    scalar_store: Dict[str, List[np.ndarray]] = {
        name: [] for name in (
            "runnum", "eventnum",
            "e_p", "e_theta", "e_phi",
            "p_p", "p_theta", "p_phi",
            "detector1", "detector2",
        )
    }
    solution_store: Dict[str, List[np.ndarray]] = {
        name: [] for name in (
            "g1_E", "g1_theta", "g1_phi",
            "g2_E", "g2_theta", "g2_phi",
            "solution_valid", "solution_closure",
        )
    }
    counts = {
        "tree_entries": 0,
        "events_with_at_least_one_solution": 0,
        "valid_solution_count": 0,
    }

    tree_entries, _available_keys = require_tree(path)
    entry_stop = (
        min(tree_entries, int(args.max_events))
        if args.max_events is not None
        else tree_entries
    )

    seen = 0
    log(f"Reading {mode} epgammagamma records from {path}")
    for arrays in uproot.iterate(
        f"{path}:{TREE_NAME}",
        expressions=expressions,
        step_size=args.step_size,
        entry_stop=entry_stop,
        library="np",
    ):
        n = len(next(iter(arrays.values())))
        seen += n
        counts["tree_entries"] += int(n)

        def arr(name: str, default=math.nan, dtype=np.float64):
            return finite_array(arrays, resolved.get(name), default, dtype)

        e_p = arr("e_p")
        e_theta = arr("e_theta")
        e_phi = np.mod(arr("e_phi"), 2.0 * math.pi)
        p_p = arr("p_p")
        p_theta = arr("p_theta")
        p_phi = np.mod(arr("p_phi"), 2.0 * math.pi)
        pi0_p = arr("pi0_p")
        pi0_theta = arr("pi0_theta")
        pi0_phi = np.mod(arr("pi0_phi"), 2.0 * math.pi)
        pi0_mass = arr("pi0_mass")
        open1 = np.deg2rad(arr("open1"))
        open2 = np.deg2rad(arr("open2"))
        trento1 = np.mod(arr("trento1"), 2.0 * math.pi)
        trento2 = np.mod(arr("trento2"), 2.0 * math.pi)

        max_solutions = 4
        g1_E = np.full((n, max_solutions), np.nan)
        g1_theta = np.full((n, max_solutions), np.nan)
        g1_phi = np.full((n, max_solutions), np.nan)
        g2_E = np.full((n, max_solutions), np.nan)
        g2_theta = np.full((n, max_solutions), np.nan)
        g2_phi = np.full((n, max_solutions), np.nan)
        valid = np.zeros((n, max_solutions), dtype=bool)
        closure = np.full((n, max_solutions), np.nan)

        for index in range(n):
            values = (
                e_p[index], e_theta[index], e_phi[index],
                pi0_p[index], pi0_theta[index], pi0_phi[index],
                pi0_mass[index], open1[index], open2[index],
                trento1[index], trento2[index],
            )
            if not all(math.isfinite(float(value)) for value in values):
                continue
            # endif

            solutions = reconstruct_daughter_solutions(
                beam_energy,
                float(e_p[index]),
                float(e_theta[index]),
                float(e_phi[index]),
                float(pi0_p[index]),
                float(pi0_theta[index]),
                float(pi0_phi[index]),
                float(pi0_mass[index]),
                float(open1[index]),
                float(open2[index]),
                float(trento1[index]),
                float(trento2[index]),
                args,
            )
            for solution_index, solution in enumerate(solutions):
                (
                    closure_value,
                    energy1, theta1, phi1,
                    energy2, theta2, phi2,
                ) = solution
                closure[index, solution_index] = closure_value
                g1_E[index, solution_index] = energy1
                g1_theta[index, solution_index] = theta1
                g1_phi[index, solution_index] = phi1
                g2_E[index, solution_index] = energy2
                g2_theta[index, solution_index] = theta2
                g2_phi[index, solution_index] = phi2
                valid[index, solution_index] = True
            # endfor
        # endfor

        counts["events_with_at_least_one_solution"] += int(
            np.count_nonzero(np.any(valid, axis=1))
        )
        counts["valid_solution_count"] += int(np.count_nonzero(valid))

        scalar_values = {
            "runnum": arr("runnum", -1, np.int64),
            "eventnum": arr("eventnum", -1, np.int64),
            "e_p": e_p,
            "e_theta": e_theta,
            "e_phi": e_phi,
            "p_p": p_p,
            "p_theta": p_theta,
            "p_phi": p_phi,
            "detector1": arr("native_gamma_detector1", -1, np.int16),
            "detector2": arr("native_gamma_detector2", -1, np.int16),
        }
        for name in scalar_store:
            scalar_store[name].append(np.asarray(scalar_values[name]))
        # endfor

        solution_values = {
            "g1_E": g1_E,
            "g1_theta": g1_theta,
            "g1_phi": g1_phi,
            "g2_E": g2_E,
            "g2_theta": g2_theta,
            "g2_phi": g2_phi,
            "solution_valid": valid,
            "solution_closure": closure,
        }
        for name in solution_store:
            solution_store[name].append(solution_values[name])
        # endfor
    # endfor

    payload = {}
    for name, parts in scalar_store.items():
        if parts:
            payload[name] = np.concatenate(parts)
        else:
            dtype = (
                np.int64
                if name in {"runnum", "eventnum"}
                else np.int16
                if name in {"detector1", "detector2"}
                else np.float64
            )
            payload[name] = np.asarray([], dtype=dtype)
        # endif
    # endfor
    for name, parts in solution_store.items():
        if parts:
            payload[name] = np.concatenate(parts, axis=0)
        else:
            dtype = bool if name == "solution_valid" else np.float64
            payload[name] = np.empty((0, 4), dtype=dtype)
        # endif
    # endfor

    return EpggRecords(**payload), counts

def read_identity_catalog(
    path: str,
    role: str,
    mode: str,
    args: argparse.Namespace,
) -> Tuple[IdentityCatalog, Dict[str, object]]:
    logical = (
        "runnum", "eventnum",
        "e_p", "e_theta", "e_phi",
        "p_p", "p_theta", "p_phi",
    )
    optional = ("runnum", "eventnum") if mode == "mc" else ()
    resolved = resolve_branches(path, logical, optional=optional)
    expressions = sorted({branch for branch in resolved.values() if branch is not None})
    store = {name: [] for name in IdentityCatalog.__dataclass_fields__}
    seen = 0
    log(f"Reading {role} identity catalog from {path}")

    for arrays in uproot.iterate(
        f"{path}:{TREE_NAME}",
        expressions=expressions,
        step_size=args.step_size,
        library="np",
    ):
        n = len(next(iter(arrays.values())))
        if args.max_events is not None and seen >= args.max_events:
            break
        # endif
        if args.max_events is not None and seen + n > args.max_events:
            n = args.max_events - seen
            arrays = {key: value[:n] for key, value in arrays.items()}
        # endif
        seen += n

        for name in store:
            branch = resolved.get(name)
            default = -1 if name in {"runnum", "eventnum"} else math.nan
            dtype = np.int64 if name in {"runnum", "eventnum"} else np.float64
            store[name].append(finite_array(arrays, branch, default=default, dtype=dtype))
        # endfor
    # endfor

    catalog = IdentityCatalog(
        runnum=concat(store["runnum"], np.int64),
        eventnum=concat(store["eventnum"], np.int64),
        e_p=concat(store["e_p"]),
        e_theta=concat(store["e_theta"]),
        e_phi=concat(store["e_phi"]),
        p_p=concat(store["p_p"]),
        p_theta=concat(store["p_theta"]),
        p_phi=concat(store["p_phi"]),
    )

    keys = identity_keys_data(catalog) if mode == "data" else identity_keys_mc(
        catalog, args.mc_signature_decimals
    )
    unique_count = int(np.unique(structured_keys(keys)).size) if catalog.size() else 0
    diagnostics = {
        "role": role,
        "mode": mode,
        "path": path,
        "tree_entries_read": int(seen),
        "catalog_records": catalog.size(),
        "unique_identity_keys": unique_count,
    }
    return catalog, diagnostics


def identity_keys_data(records) -> np.ndarray:
    return np.column_stack(
        (
            np.asarray(records.runnum, dtype=np.int64),
            np.asarray(records.eventnum, dtype=np.int64),
        )
    )


def identity_keys_mc(records, decimals: int) -> np.ndarray:
    values = (
        records.e_p, records.e_theta, records.e_phi,
        records.p_p, records.p_theta, records.p_phi,
    )
    return np.column_stack(
        [np.round(np.asarray(value, dtype=np.float64), decimals) for value in values]
    )


def identity_keys_epgg_data(records: EpggRecords) -> np.ndarray:
    return np.column_stack((records.runnum.astype(np.int64), records.eventnum.astype(np.int64)))


def identity_keys_epgg_mc(records: EpggRecords, decimals: int) -> np.ndarray:
    return np.column_stack([
        np.round(records.e_p, decimals),
        np.round(records.e_theta, decimals),
        np.round(records.e_phi, decimals),
        np.round(records.p_p, decimals),
        np.round(records.p_theta, decimals),
        np.round(records.p_phi, decimals),
    ])


def structured_keys(matrix: np.ndarray) -> np.ndarray:
    matrix = np.ascontiguousarray(matrix)
    dtype = np.dtype([(f"f{index}", matrix.dtype) for index in range(matrix.shape[1])])
    return matrix.view(dtype).reshape(-1)


def group_indices(keys: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    structured = structured_keys(keys)
    _, inverse, counts = np.unique(structured, return_inverse=True, return_counts=True)
    order = np.argsort(inverse, kind="stable")
    offsets = np.empty(counts.size + 1, dtype=np.int64)
    offsets[0] = 0
    np.cumsum(counts, out=offsets[1:])
    return inverse, order, offsets


def photon_pair_mass(
    E1: float,
    theta1: float,
    phi1: float,
    E2: float,
    theta2: float,
    phi2: float,
) -> float:
    cosine = (
        math.sin(theta1) * math.sin(theta2) * math.cos(phi1 - phi2)
        + math.cos(theta1) * math.cos(theta2)
    )
    mass2 = 2.0 * E1 * E2 * (1.0 - max(-1.0, min(1.0, cosine)))
    return math.sqrt(max(mass2, 0.0))


def match_truth_partners(
    opportunities: OpportunityRecords,
    epgg: EpggRecords,
    mode: str,
    args: argparse.Namespace,
) -> PartnerMatchResult:
    n = opportunities.size()
    matched = np.zeros(n, dtype=bool)
    catalog_confirmed = np.zeros(n, dtype=bool)
    partner_index = np.full(n, -1, dtype=np.int64)
    probe_angle = np.full(n, np.nan)
    probe_relative_E = np.full(n, np.nan)
    pair_mass = np.full(n, np.nan)
    score = np.full(n, np.nan)
    partner_E = np.full(n, np.nan)
    tag_angle = np.full(n, np.nan)
    tag_relative_E = np.full(n, np.nan)
    closure = np.full(n, np.nan)

    opportunity_keys = identity_keys_data(opportunities) if mode == "data" else identity_keys_mc(opportunities, args.mc_signature_decimals)
    native_keys = identity_keys_epgg_data(epgg) if mode == "data" else identity_keys_epgg_mc(epgg, args.mc_signature_decimals)

    native_structured = structured_keys(native_keys)
    unique, inverse, counts = np.unique(native_structured, return_inverse=True, return_counts=True)
    order = np.argsort(inverse, kind="stable")
    offsets = np.empty(counts.size + 1, dtype=np.int64)
    offsets[0] = 0
    np.cumsum(counts, out=offsets[1:])
    group_by_key = {unique[i].tobytes(): i for i in range(unique.size)}

    identity_matches = valid_opportunities = tag_matches = partner_threshold = probe_passes = 0
    records_tested = solutions_tested = 0

    for oi, key in enumerate(structured_keys(opportunity_keys)):
        group = group_by_key.get(key.tobytes())
        if group is None:
            continue
        # endif
        catalog_confirmed[oi] = True
        identity_matches += 1
        best = None
        any_valid = False
        members = order[offsets[group]:offsets[group + 1]]
        for ni in members:
            records_tested += 1
            for si in range(epgg.solution_valid.shape[1]):
                if not epgg.solution_valid[ni, si]:
                    continue
                # endif
                any_valid = True
                solutions_tested += 1
                assignments = (
                    (epgg.g1_E[ni,si], epgg.g1_theta[ni,si], epgg.g1_phi[ni,si],
                     epgg.g2_E[ni,si], epgg.g2_theta[ni,si], epgg.g2_phi[ni,si]),
                    (epgg.g2_E[ni,si], epgg.g2_theta[ni,si], epgg.g2_phi[ni,si],
                     epgg.g1_E[ni,si], epgg.g1_theta[ni,si], epgg.g1_phi[ni,si]),
                )
                for tE,tth,tph,pE,pth,pph in assignments:
                    ta = opening_angle_deg_scalar(
                        float(opportunities.tag_theta[oi]),
                        float(opportunities.tag_phi[oi]),
                        float(tth),
                        float(tph),
                    )
                    tr = abs(opportunities.tag_E[oi]-tE)/max(tE,1e-12)
                    ts = (ta/max(args.tag_match_angle_max_deg,1e-12))**2 + (tr/max(args.tag_match_relative_E_max,1e-12))**2
                    if best is None or ts < best[0]:
                        pa = opening_angle_deg_scalar(
                            float(opportunities.probe_theta[oi]),
                            float(opportunities.probe_phi[oi]),
                            float(pth),
                            float(pph),
                        )
                        pr = abs(opportunities.probe_E[oi]-pE)/max(pE,1e-12)
                        pm = photon_pair_mass(float(opportunities.tag_E[oi]),float(opportunities.tag_theta[oi]),float(opportunities.tag_phi[oi]),float(pE),float(pth),float(pph))
                        best=(ts,ta,tr,pa,pr,pm,float(pE),int(ni),float(epgg.solution_closure[ni,si]))
                    # endif
                # endfor
            # endfor
        # endfor
        if any_valid:
            valid_opportunities += 1
        # endif
        if best is None:
            continue
        # endif
        ts,ta,tr,pa,pr,pm,pE,ni,cl = best
        partner_index[oi]=ni; probe_angle[oi]=pa; probe_relative_E[oi]=pr
        pair_mass[oi]=pm; score[oi]=ts; partner_E[oi]=pE
        tag_angle[oi]=ta; tag_relative_E[oi]=tr; closure[oi]=cl
        if not (ta < args.tag_match_angle_max_deg and tr < args.tag_match_relative_E_max):
            continue
        # endif
        tag_matches += 1
        if pE < args.found_probe_E_min:
            continue
        # endif
        partner_threshold += 1
        probe_ok = pa < args.probe_match_angle_max_deg and pr < args.probe_match_relative_E_max
        probe_passes += int(probe_ok)
        if args.require_probe_residual_match and not probe_ok:
            continue
        # endif
        matched[oi]=True
    # endfor

    summary={
        "mode":mode,"opportunities":n,"identity_matches":identity_matches,
        "valid_solution_opportunities":valid_opportunities,"tag_matched":tag_matches,
        "partner_above_threshold":partner_threshold,"probe_residual_passes":probe_passes,
        "matched_opportunities":int(np.count_nonzero(matched)),
        "epgammagamma_records_tested":records_tested,"daughter_solutions_tested":solutions_tested,
        "fraction_identity_matched":identity_matches/n if n else math.nan,
        "fraction_matched":float(np.count_nonzero(matched))/n if n else math.nan,
        "conditional_valid_solution_given_identity":valid_opportunities/identity_matches if identity_matches else math.nan,
        "conditional_tag_match_given_valid_solution":tag_matches/valid_opportunities if valid_opportunities else math.nan,
        "require_probe_residual_match":bool(args.require_probe_residual_match),
    }
    return PartnerMatchResult(
        matched,catalog_confirmed,partner_index,probe_angle,probe_relative_E,pair_mass,score,
        partner_E,tag_angle,tag_relative_E,closure,summary
    )



def category_mask(
    records: OpportunityRecords,
    detector: str,
    sector: int,
) -> np.ndarray:
    if detector == "FT":
        return records.probe_detector == 0
    # endif
    return (records.probe_detector == 1) & (records.probe_sector == sector)


def variable_values(records: OpportunityRecords, key: str) -> np.ndarray:
    return np.asarray(getattr(records, key), dtype=np.float64)


def histograms_for_mask(
    records: OpportunityRecords,
    mask: np.ndarray,
) -> Dict[str, np.ndarray]:
    histograms: Dict[str, np.ndarray] = {}
    for variable in FIT_VARIABLES:
        values = variable_values(records, variable.key)[mask]
        values = values[np.isfinite(values)]
        histograms[variable.key], _ = np.histogram(
            values,
            bins=variable.bins,
            range=(variable.low, variable.high),
        )
    # endfor
    return histograms







def plot_expected_probe_diagnostics(
    path: Path,
    period_label: str,
    samples: Mapping[str, OpportunityRecords],
) -> None:
    variables = (
        ("tag_E", np.linspace(0.4, 9.5, 92), r"$E_{\rm tag}$ (GeV)"),
        ("probe_E", np.linspace(2.0, 9.5, 76), r"$E_{\rm probe}^{\rm pred}$ (GeV)"),
        ("probe_theta", np.linspace(math.radians(2.5), math.radians(35.0), 80), r"$\theta_{\rm probe}^{\rm pred}$ (rad)"),
        ("probe_phi", np.linspace(0.0, 2.0 * math.pi, 73), r"$\phi_{\rm probe}^{\rm pred}$ (rad)"),
        ("probe_m2", np.linspace(-0.1, 0.1, 101), r"$m_X^2$ (GeV$^2$)"),
        ("probe_E_minus_p", np.linspace(-0.1, 0.1, 101), r"$E_X-|\vec p_X|$ (GeV)"),
    )
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    for axis, (name, bins, label) in zip(axes.flat, variables):
        for sample, records in samples.items():
            values = np.asarray(getattr(records, name), dtype=float)
            values = values[np.isfinite(values)]
            if values.size == 0:
                continue
            # endif
            axis.hist(
                values,
                bins=bins,
                weights=np.ones(values.size) / values.size,
                histtype="step",
                linewidth=1.3,
                label=f"{sample} ({values.size:,})",
            )
        # endfor
        axis.set_xlabel(label)
        axis.set_ylabel("Fraction / bin")
        axis.grid(alpha=0.25)
        axis.legend(fontsize=7)
    # endfor
    fig.suptitle(f"{period_label}: selected expected-probe opportunities")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=180)
    plt.close(fig)

def plot_matching_residuals(
    path: Path,
    period_label: str,
    data_match: PartnerMatchResult,
    mc_match: PartnerMatchResult,
) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    panels = (
        ("tag_angle_deg", "Tag angular residual (deg)", 0.0, 5.0),
        ("tag_relative_E", "Tag relative energy residual", 0.0, 0.5),
        ("probe_angle_deg", "Predicted-X to partner angle (deg)", 0.0, 20.0),
        ("probe_relative_E", "Predicted-X to partner relative energy", 0.0, 1.5),
        ("partner_E", "Reconstructed partner energy (GeV)", 0.0, 10.0),
        ("solution_closure", "Native daughter closure", 0.0, 0.10),
    )
    for axis,(field,label,low,high) in zip(axes.flat,panels):
        for sample,result in (("Data",data_match),("AAOGEN MC",mc_match)):
            values=np.asarray(getattr(result,field),dtype=float)
            values=values[np.isfinite(values)]
            if values.size:
                hist,edges=np.histogram(values,bins=np.linspace(low,high,101))
                if hist.sum()>0: hist=hist/hist.sum()
                centers=0.5*(edges[:-1]+edges[1:])
                axis.step(centers,hist,where="mid",label=sample)
            # endif
        # endfor
        axis.set_xlabel(label); axis.set_ylabel("Normalized entries")
        axis.set_xlim(low,high); axis.grid(alpha=0.25); axis.legend(fontsize=8)
    # endfor
    fig.suptitle(f"{period_label}: native truth-partner matching")
    fig.tight_layout(rect=(0,0,1,0.95)); fig.savefig(path,dpi=180); plt.close(fig)

def plot_cutflow(path: Path, title: str, cutflow: Mapping[str, object]) -> None:
    rows = cutflow.get("stages", [])
    if not rows:
        return
    # endif
    labels = [str(row["stage"]) for row in rows]
    counts = np.asarray([float(row["count"]) for row in rows])
    total_survival = np.asarray(
        [float(row["total_survival"]) for row in rows]
    )
    incremental = np.asarray(
        [float(row["incremental_survival"]) for row in rows]
    )

    fig, axes = plt.subplots(2, 1, figsize=(14, 10), sharex=True)
    x = np.arange(len(rows))

    axes[0].step(x, np.maximum(counts, 0.5), where="mid", linewidth=1.5)
    axes[0].scatter(x, np.maximum(counts, 0.5), s=25)
    axes[0].set_yscale("log")
    axes[0].set_ylabel("Cumulative entries")
    axes[0].grid(alpha=0.25)

    axes[1].plot(x, 100.0 * total_survival, marker="o", label="Total survival")
    axes[1].plot(x, 100.0 * incremental, marker="s", label="Incremental survival")
    axes[1].set_ylabel("Survival (%)")
    axes[1].set_ylim(0.0, 105.0)
    axes[1].grid(alpha=0.25)
    axes[1].legend()
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels, rotation=45, ha="right")

    suspect_stages = {"z_tag", "minus_t1_tag", "open_angle_ep2_tag"}
    for axis in axes:
        for index, stage in enumerate(labels):
            if stage in suspect_stages:
                axis.axvspan(index - 0.45, index + 0.45, alpha=0.10)
            # endif
        # endfor
    # endfor

    largest = cutflow.get("largest_incremental_loss")
    subtitle = ""
    if isinstance(largest, Mapping):
        subtitle = (
            f"Largest incremental loss: {largest['stage']} — "
            f"{100.0 * float(largest['incremental_rejection']):.2f}% rejected"
        )
    # endif
    fig.suptitle(f"{title}\n{subtitle}")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=180)
    plt.close(fig)

def write_cutflow_csv(path: Path, cutflow: Mapping[str, object]) -> None:
    """Write the finalized sequential cutflow table to CSV."""
    rows = list(cutflow.get("stages", []))
    if not rows:
        return
    # endif

    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "order",
        "stage",
        "label",
        "condition",
        "branch",
        "scope",
        "applied",
        "count",
        "rejected_at_stage",
        "incremental_survival",
        "incremental_rejection",
        "total_survival",
    ]

    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    field: row.get(field)
                    for field in fieldnames
                }
            )

def write_cutflow_text(path: Path, cutflow: Mapping[str, object]) -> None:
    lines = [
        f"Role: {cutflow.get('role')}",
        f"Path: {cutflow.get('path')}",
        "",
        (
            "Columns: cumulative count | incremental survival | total survival "
            "| rejected at stage | cut description"
        ),
        "",
    ]
    for row in cutflow.get("stages", []):
        branch = row.get("branch")
        branch_text = f"; branch={branch}" if branch else ""
        applied_text = "applied" if row.get("applied") else "NOT APPLIED (branch/option absent)"
        lines.append(
            f"{int(row['order']):2d} {row['stage']:<34s} "
            f"{int(row['count']):>12,d} | "
            f"{100.0 * float(row['incremental_survival']):>8.3f}% | "
            f"{100.0 * float(row['total_survival']):>8.3f}% | "
            f"rejected={int(row['rejected_at_stage']):>10,d} | "
            f"{row['condition']} [{applied_text}; scope={row['scope']}"
            f"{branch_text}]"
        )
    # endfor

    lines.extend(
        [
            "",
            "Predicted-probe detector populations:",
            json.dumps(cutflow.get("detector_counts", {}), indent=2),
            "",
            "NOTE: z_tag, minus_t1_tag, and open_angle_ep2_tag are intentionally "
            "NOT APPLIED. They remain in the stage table with 100% incremental "
            "survival so the removed selection is explicit in the provenance.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")

def emit_cutflow_diagnostics(
    period_dir: Path,
    period_label: str,
    named_cutflows: Mapping[str, Mapping[str, object]],
) -> None:
    audit_dir = period_dir / "cutflow_diagnostics"
    audit_dir.mkdir(parents=True, exist_ok=True)

    for sample_name, cutflow in named_cutflows.items():
        safe = sample_name.replace(" ", "_").replace("/", "_")
        write_cutflow_csv(audit_dir / f"{safe}_cutflow.csv", cutflow)
        write_cutflow_text(audit_dir / f"{safe}_cutflow.txt", cutflow)
        plot_cutflow(
            audit_dir / f"{safe}_cutflow.png",
            f"{period_label}: {sample_name} opportunity cutflow",
            cutflow,
        )

        largest = cutflow.get("largest_incremental_loss")
        if isinstance(largest, Mapping):
            log(
                f"{period_label} {sample_name} cutflow: "
                f"{int(largest['count']):,} survive {largest['stage']}; "
                f"incremental survival="
                f"{100.0 * float(largest['incremental_survival']):.3f}%, "
                f"total survival="
                f"{100.0 * float(largest['total_survival']):.3f}%"
            )
        # endif

        detector_counts = cutflow.get("detector_counts", {})
        for row in cutflow.get("stages", []):
            detector_suffix = ""
            if row["stage"] in detector_counts:
                stage_detector_counts = detector_counts[row["stage"]]
                detector_suffix = (
                    f"; FT={int(stage_detector_counts.get('FT', 0)):,}, "
                    f"FD={int(stage_detector_counts.get('FD', 0)):,}"
                )
            # endif
            log(
                f"{period_label} {sample_name}: "
                f"{row['stage']}: {int(row['count']):,} "
                f"({100.0 * float(row['incremental_survival']):.3f}% of prior, "
                f"{100.0 * float(row['total_survival']):.3f}% of input)"
                f"{detector_suffix} — {row['condition']}"
            )


@dataclass
class FailOnlyEfficiencyRow:
    period: str
    period_label: str
    detector: str
    sector: int
    data_pass_count: int
    data_fail_count: int
    pi0_mc_pass_count: int
    pi0_mc_fail_count: int
    dvcs_mc_count: int
    fit_fraction_pi0_fail: float
    fit_fraction_pi0_fail_stat_err: float
    fit_fraction_pi0_fail_model_spread: float
    fitted_pi0_fail: float
    fitted_pi0_fail_stat_err: float
    fitted_pi0_fail_model_spread: float
    efficiency_data: float
    efficiency_data_stat_err: float
    efficiency_data_model_err: float
    efficiency_mc: float
    efficiency_mc_err: float
    scale_factor: float
    scale_factor_stat_err: float
    scale_factor_model_err: float
    fit_success: bool
    fit_message: str
    fit_deviance: float
    fit_ndf: int
    fit_reduced_deviance: float
    fit_quality_warning: bool
    per_variable_fractions: str

@dataclass
class FractionFit:
    success: bool
    message: str
    fraction: float
    stat_error: float
    nll: float
    deviance: float
    ndf: int
    reduced_deviance: float
    per_variable: Dict[str, Dict[str, float]]
    payload: Dict[str, Dict[str, object]]
    warnings: List[str]

FAIL_FIT_DRIVERS: Tuple[str,...] = FRACTION_DRIVERS


@dataclass
class MCEfficiencyRow:
    period: str
    period_label: str
    topology: str
    opportunities: int
    reconstructed: int
    missed: int
    efficiency: float
    efficiency_error: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Produce before-cuts and after-cuts normalized epgamma shape "
            "comparisons, with the completed AAOGEN MC efficiency available "
            "as an optional stage."
        )
    )
    parser.add_argument(
        "--period",
        action="append",
        choices=[period.key for period in PERIODS],
        help="Restrict processing to one or more run periods.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=8,
        help=f"Number of period workers; hard maximum is {MAX_WORKERS}.",
    )
    parser.add_argument("--step-size", default=DEFAULT_STEP_SIZE)
    parser.add_argument(
        "--max-events",
        type=int,
        default=None,
        help="Optional debugging cap per input ROOT tree.",
    )
    parser.add_argument(
        "--output-dir",
        default=DEFAULT_OUTPUT_DIR,
        help="Common output root containing top-level data/ and mc/ folders.",
    )
    parser.add_argument(
        "--shape-log-y",
        action="store_true",
        help="Also write logarithmic-y normalized shape canvases.",
    )
    parser.add_argument(
        "--simple",
        action="store_true",
        help=(
            "Temporarily replace the normal DATA shape/template section with "
            "a direct <2 GeV tag / >2 GeV same-(runnum,evnum) partner "
            "efficiency study. The separate MC-efficiency stage is unchanged."
        ),
    )
    parser.add_argument(
        "--fast-fits",
        action="store_true",
        help=(
            "Use diagnostic-speed template-fit settings: fewer multistarts, "
            "looser optimizer tolerances, and fewer coordinate iterations. "
            "The default retains the full production fit settings."
        ),
    )
    parser.add_argument(
        "--print-stage-timings",
        action="store_true",
        help="Print detailed per-file and per-stage runtime measurements.",
    )

    parser.add_argument(
        "--data-pTmiss-min",
        type=float,
        default=0.5,
        help=(
            "Deprecated compatibility option. No pTmiss lower cut is applied "
            "in the current data-template workflow."
        ),
    )

    parser.add_argument(
        "--no-mc",
        "--skip-mc-efficiency",
        dest="skip_mc_efficiency",
        action="store_true",
        help=(
            "Skip the full AAOGEN MC photon-efficiency stage while running "
            "the data diagnostics and fits. The MC study runs by default."
        ),
    )
    parser.add_argument(
        "--run-mc-efficiency",
        action="store_true",
        help=(
            "Deprecated compatibility option. The full MC-efficiency stage "
            "now runs by default."
        ),
    )

    # Options retained by the optional AAOGEN MC efficiency stage.
    parser.add_argument("--tag-E-min", type=float, default=0.40)
    parser.add_argument("--probe-E-min", type=float, default=2.00)
    parser.add_argument("--probe-E-max", type=float, default=9.50)
    parser.add_argument("--probe-m2-abs-max", type=float, default=0.10)
    parser.add_argument(
        "--probe-E-minus-p-abs-max",
        type=float,
        default=0.10,
    )
    parser.add_argument("--ft-theta-min", type=float, default=2.5)
    parser.add_argument("--ft-theta-max", type=float, default=5.0)
    parser.add_argument("--fd-theta-min", type=float, default=5.0)
    parser.add_argument("--fd-theta-max", type=float, default=35.0)
    parser.add_argument("--Q2-min", type=float, default=1.0)
    parser.add_argument("--W-min", type=float, default=2.0)
    parser.add_argument("--y-max", type=float, default=0.8)
    parser.add_argument("--z-min", type=float, default=0.65)
    parser.add_argument("--minus-t-max", type=float, default=1.0)
    parser.add_argument("--open-angle-min-deg", type=float, default=5.0)
    parser.add_argument(
        "--require-fiducial-status-111",
        action="store_true",
    )
    parser.add_argument("--tag-match-angle-max-deg", type=float, default=3.0)
    parser.add_argument(
        "--tag-match-relative-E-max",
        type=float,
        default=0.35,
    )
    parser.add_argument(
        "--probe-match-angle-max-deg",
        type=float,
        default=3.0,
    )
    parser.add_argument(
        "--probe-match-relative-E-max",
        type=float,
        default=0.35,
    )
    parser.add_argument("--mc-signature-decimals", type=int, default=10)
    parser.add_argument("--native-closure-max", type=float, default=0.10)
    parser.add_argument(
        "--native-minimum-photon-E",
        type=float,
        default=0.20,
    )
    parser.add_argument("--found-probe-E-min", type=float, default=2.00)
    parser.add_argument(
        "--require-probe-residual-match",
        action="store_true",
    )
    parser.add_argument("--energy-bins", type=int, default=15)
    parser.add_argument("--theta-bins", type=int, default=13)
    parser.add_argument("--phi-bins", type=int, default=36)
    parser.add_argument(
        "--minimum-bin-opportunities",
        type=int,
        default=20,
    )
    parser.add_argument(
        "--skip-data-template-fits",
        action="store_true",
        help="Skip the advanced DVCSGEN+AAOGEN data-template fits.",
    )
    parser.add_argument(
        "--allow-data-fit-failures",
        action="store_true",
        help=(
            "Serialize failed categories and continue instead of aborting "
            "the complete run."
        ),
    )
    parser.add_argument(
        "--fit-min-counts",
        type=int,
        default=100,
        help="Minimum in-range data count required in every projection.",
    )
    parser.add_argument(
        "--fit-max-shift-bins",
        type=float,
        default=12.0,
        help="Maximum absolute component shift in transformed histogram bins.",
    )
    parser.add_argument(
        "--fit-max-smear-bins",
        type=float,
        default=20.0,
        help="Maximum additional Gaussian broadening in bins.",
    )
    parser.add_argument(
        "--fit-max-log-stretch",
        type=float,
        default=0.45,
        help=(
            "Maximum absolute logarithm of the affine scale factor. "
            "A value of 0.45 permits scales from exp(-0.45) to exp(0.45)."
        ),
    )
    parser.add_argument(
        "--fit-shift-prior-bins",
        type=float,
        default=4.0,
        help="Gaussian prior width for component shifts in bins.",
    )
    parser.add_argument(
        "--fit-smear-prior-bins",
        type=float,
        default=8.0,
        help="Half-Gaussian prior width for additional broadenings in bins.",
    )
    parser.add_argument(
        "--fit-log-stretch-prior",
        type=float,
        default=0.18,
        help="Gaussian prior width for the logarithmic affine scale.",
    )
    parser.add_argument(
        "--fit-component-weight-prior",
        type=float,
        default=0.55,
        help=(
            "Gaussian prior width for the component-weight logit change in "
            "the two-component Delta_phi and theta_gamma_gamma models."
        ),
    )
    parser.add_argument(
        "--fit-mean-order-penalty",
        type=float,
        default=1.0e6,
        help=(
            "Quadratic penalty coefficient when the morphed DVCSGEN mean "
            "exceeds the morphed AAOGEN mean for Emiss2 or Mx2_2."
        ),
    )
    parser.add_argument(
        "--fit-coordinate-iterations",
        type=int,
        default=5,
        help="Maximum shared-fraction coordinate-minimization iterations.",
    )
    parser.add_argument(
        "--fit-independent-multistarts",
        type=int,
        default=4,
        help="Number of independent-fraction starting points per variable.",
    )
    parser.add_argument(
        "--disable-fit-nuisance-penalties",
        action="store_true",
        help="Disable all template-morph nuisance penalties.",
    )
    parser.add_argument(
        "--fit-dpi",
        type=int,
        default=180,
        help="Resolution of template-fit diagnostic plots.",
    )
    return parser.parse_args()


def binomial_efficiency(
    reconstructed: int,
    opportunities: int,
) -> Tuple[float, float]:
    if opportunities <= 0:
        return math.nan, math.nan
    # endif
    if reconstructed < 0 or reconstructed > opportunities:
        raise ValueError(
            "Invalid binomial population: "
            f"reconstructed={reconstructed}, opportunities={opportunities}"
        )
    # endif

    efficiency = reconstructed / opportunities
    uncertainty = math.sqrt(
        max(efficiency * (1.0 - efficiency) / opportunities, 0.0)
    )
    return efficiency, uncertainty


def topology_mask(
    records: OpportunityRecords,
    topology: str,
) -> np.ndarray:
    if topology == "FT":
        return records.probe_detector == 0
    # endif
    if topology == "FD":
        return records.probe_detector == 1
    # endif
    raise ValueError(f"Unsupported topology: {topology}")


def efficiency_row(
    period: PeriodConfig,
    topology: str,
    mask: np.ndarray,
    matched: np.ndarray,
) -> MCEfficiencyRow:
    opportunities = int(np.count_nonzero(mask))
    reconstructed = int(np.count_nonzero(mask & matched))
    efficiency, uncertainty = binomial_efficiency(
        reconstructed,
        opportunities,
    )
    return MCEfficiencyRow(
        period=period.key,
        period_label=period.label,
        topology=topology,
        opportunities=opportunities,
        reconstructed=reconstructed,
        missed=opportunities - reconstructed,
        efficiency=efficiency,
        efficiency_error=uncertainty,
    )


def binned_efficiency(
    values: np.ndarray,
    selected: np.ndarray,
    matched: np.ndarray,
    edges: np.ndarray,
) -> List[Dict[str, object]]:
    values = np.asarray(values, dtype=float)
    selected = np.asarray(selected, dtype=bool)
    matched = np.asarray(matched, dtype=bool)

    finite = selected & np.isfinite(values)
    denominator, _ = np.histogram(values[finite], bins=edges)
    numerator, _ = np.histogram(values[finite & matched], bins=edges)

    rows: List[Dict[str, object]] = []
    for index in range(edges.size - 1):
        efficiency, uncertainty = binomial_efficiency(
            int(numerator[index]),
            int(denominator[index]),
        )
        rows.append(
            {
                "bin_index": index,
                "bin_low": float(edges[index]),
                "bin_high": float(edges[index + 1]),
                "bin_center": float(
                    0.5 * (edges[index] + edges[index + 1])
                ),
                "opportunities": int(denominator[index]),
                "reconstructed": int(numerator[index]),
                "missed": int(denominator[index] - numerator[index]),
                "efficiency": efficiency,
                "efficiency_error": uncertainty,
            }
        )
    # endfor
    return rows


def write_dict_rows(
    path: Path,
    rows: Sequence[Mapping[str, object]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    # endif

    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0].keys()),
        )
        writer.writeheader()
        writer.writerows(rows)
    # endwith


def plot_matching_residuals_mc_only(
    path: Path,
    period_label: str,
    match: PartnerMatchResult,
) -> None:
    panels = (
        ("tag_angle_deg", "Tag angular residual (deg)", 0.0, 5.0),
        (
            "tag_relative_E",
            "Tag relative energy residual",
            0.0,
            0.5,
        ),
        (
            "probe_angle_deg",
            "Predicted probe to partner angle (deg)",
            0.0,
            20.0,
        ),
        (
            "probe_relative_E",
            "Predicted probe to partner relative energy",
            0.0,
            1.5,
        ),
        (
            "partner_E",
            "Reconstructed partner energy (GeV)",
            0.0,
            10.0,
        ),
        (
            "solution_closure",
            "Native daughter-solution closure",
            0.0,
            0.10,
        ),
    )

    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    for axis, (field, label, low, high) in zip(axes.flat, panels):
        values = np.asarray(getattr(match, field), dtype=float)
        values = values[np.isfinite(values)]
        counts, edges = np.histogram(
            values,
            bins=np.linspace(low, high, 101),
        )
        if counts.sum() > 0:
            counts = counts / counts.sum()
        # endif
        centers = 0.5 * (edges[:-1] + edges[1:])
        axis.step(centers, counts, where="mid")
        axis.set_xlabel(label)
        axis.set_ylabel("Normalized entries")
        axis.set_xlim(low, high)
        axis.grid(alpha=0.25)
    # endfor

    fig.suptitle(f"{period_label}: AAOGEN truth-partner matching")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_efficiency_dependence(
    path: Path,
    title: str,
    variable_label: str,
    topology_rows: Mapping[str, Sequence[Mapping[str, object]]],
    minimum_opportunities: int,
) -> None:
    fig, axis = plt.subplots(figsize=(11, 7))

    for topology, rows in topology_rows.items():
        centers = np.asarray(
            [row["bin_center"] for row in rows],
            dtype=float,
        )
        efficiencies = np.asarray(
            [row["efficiency"] for row in rows],
            dtype=float,
        )
        uncertainties = np.asarray(
            [row["efficiency_error"] for row in rows],
            dtype=float,
        )
        support = np.asarray(
            [row["opportunities"] for row in rows],
            dtype=int,
        )
        valid = (
            np.isfinite(efficiencies)
            & np.isfinite(uncertainties)
            & (support >= minimum_opportunities)
        )
        axis.errorbar(
            centers[valid],
            efficiencies[valid],
            yerr=uncertainties[valid],
            fmt="o-",
            capsize=2,
            linewidth=1.0,
            markersize=4,
            label=topology,
        )
    # endfor

    axis.set_xlabel(variable_label)
    axis.set_ylabel(r"$\epsilon_{\mathrm{MC}}$")
    axis.set_ylim(0.0, 1.05)
    axis.grid(alpha=0.25)
    axis.legend()
    axis.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_fd_sector_diagnostics(
    path: Path,
    period_label: str,
    records: OpportunityRecords,
    matched: np.ndarray,
) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for sector in range(1, 7):
        mask = (
            (records.probe_detector == 1)
            & (records.probe_sector == sector)
        )
        opportunities = int(np.count_nonzero(mask))
        reconstructed = int(np.count_nonzero(mask & matched))
        efficiency, uncertainty = binomial_efficiency(
            reconstructed,
            opportunities,
        )
        rows.append(
            {
                "sector": sector,
                "opportunities": opportunities,
                "reconstructed": reconstructed,
                "missed": opportunities - reconstructed,
                "efficiency": efficiency,
                "efficiency_error": uncertainty,
            }
        )
    # endfor

    x = np.arange(1, 7)
    efficiency = np.asarray(
        [row["efficiency"] for row in rows],
        dtype=float,
    )
    uncertainty = np.asarray(
        [row["efficiency_error"] for row in rows],
        dtype=float,
    )

    fig, axis = plt.subplots(figsize=(10, 6))
    axis.errorbar(
        x,
        efficiency,
        yerr=uncertainty,
        fmt="o",
        capsize=3,
    )
    axis.set_xticks(x)
    axis.set_xlabel("Predicted FD photon sector")
    axis.set_ylabel(r"$\epsilon_{\mathrm{MC}}$")
    axis.set_ylim(0.0, 1.05)
    axis.grid(alpha=0.25)
    axis.set_title(
        f"{period_label}: FD-sector diagnostic photon efficiency"
    )
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return rows



RAW_SHAPE_SAMPLE_INFO: Tuple[Tuple[str, str, str], ...] = (
    ("data", r"$e'p'\gamma$ data", "epg_data"),
    ("dvcsgen", "DVCS MC", "dvcs_mc"),
    ("aaogen", r"$e\pi^0$ MC as DVCS", "pi0_epg_mc"),
)


def resolve_branches_from_keys(
    path: str,
    available_keys: Sequence[str],
    logical_names: Sequence[str],
) -> Dict[str, str]:
    """
    Resolve logical branches from an already opened tree branch list.

    This avoids reopening every ROOT file inside resolve_branches().
    """
    available = set(available_keys)
    resolved: Dict[str, str] = {}
    missing: List[str] = []

    for logical_name in logical_names:
        aliases = ALIASES.get(logical_name, (logical_name,))
        branch = next(
            (candidate for candidate in aliases if candidate in available),
            None,
        )
        if branch is None:
            missing.append(
                f"{logical_name}: tried {', '.join(aliases)}"
            )
        else:
            resolved[logical_name] = branch
        # endif
    # endfor

    if missing:
        raise KeyError(
            f"Missing required branches in {path}:\n  "
            + "\n  ".join(missing)
        )
    # endif
    return resolved


def empty_stage_histograms(
    variables: Sequence[FitVariable],
) -> Dict[str, Dict[str, object]]:
    histograms: Dict[str, Dict[str, object]] = {}
    for variable in variables:
        edges = np.linspace(
            variable.low,
            variable.high,
            variable.bins + 1,
        )
        histograms[variable.key] = {
            "edges": edges,
            "counts": np.zeros(variable.bins, dtype=np.int64),
            "selected_entries": 0,
            "in_range_entries": 0,
            "underflow_entries": 0,
            "overflow_entries": 0,
        }
    # endfor
    return histograms


def update_stage_histograms(
    stage_histograms: Dict[str, Dict[str, object]],
    arrays_by_logical_name: Mapping[str, np.ndarray],
    mask: np.ndarray,
    variables: Sequence[FitVariable],
) -> None:
    """Fill the requested category-specific shape histograms."""
    radians_to_degrees = 180.0 / math.pi

    for variable in variables:
        selected = np.asarray(
            arrays_by_logical_name[variable.key],
            dtype=float,
        )[mask]

        if variable.key in ("p1_theta", "p2_theta"):
            selected = selected * radians_to_degrees
        # endif

        histogram = stage_histograms[variable.key]
        edges = np.asarray(histogram["edges"], dtype=float)
        counts, _ = np.histogram(selected, bins=edges)

        histogram["counts"] += counts.astype(np.int64)
        histogram["selected_entries"] += int(selected.size)
        histogram["in_range_entries"] += int(np.sum(counts))
        histogram["underflow_entries"] += int(
            np.count_nonzero(selected < variable.low)
        )
        histogram["overflow_entries"] += int(
            np.count_nonzero(selected >= variable.high)
        )
    # endfor
    # endfor
    # endfor
    # endfor
    # endfor


def finalize_stage_histograms(
    stage_histograms: Dict[str, Dict[str, object]],
    variables: Sequence[FitVariable],
) -> Dict[str, Dict[str, object]]:
    finalized: Dict[str, Dict[str, object]] = {}

    for variable in variables:
        histogram = stage_histograms[variable.key]
        counts = np.asarray(histogram["counts"], dtype=np.int64)
        in_range = int(histogram["in_range_entries"])
        selected = int(histogram["selected_entries"])

        finalized[variable.key] = {
            "edges": np.asarray(histogram["edges"], dtype=float),
            "counts": counts,
            "unit_area": (
                counts.astype(float) / in_range
                if in_range > 0
                else np.zeros_like(counts, dtype=float)
            ),
            "selected_entries": selected,
            "in_range_entries": in_range,
            "underflow_entries": int(
                histogram["underflow_entries"]
            ),
            "overflow_entries": int(
                histogram["overflow_entries"]
            ),
            "in_range_fraction_of_selected": (
                in_range / selected if selected > 0 else math.nan
            ),
            "underflow_fraction_of_selected": (
                int(histogram["underflow_entries"]) / selected
                if selected > 0
                else math.nan
            ),
            "overflow_fraction_of_selected": (
                int(histogram["overflow_entries"]) / selected
                if selected > 0
                else math.nan
            ),
        }
    # endfor

    return finalized


def uint64_mix(values: np.ndarray) -> np.ndarray:
    """
    SplitMix64 finalizer used to build a deterministic 128-bit MC signature.
    """
    values = np.asarray(values, dtype=np.uint64)
    values = values ^ (values >> np.uint64(30))
    values = values * np.uint64(0xBF58476D1CE4E5B9)
    values = values ^ (values >> np.uint64(27))
    values = values * np.uint64(0x94D049BB133111EB)
    values = values ^ (values >> np.uint64(31))
    return values


def data_identity_signature(
    arrays_by_logical_name: Mapping[str, np.ndarray],
) -> Tuple[np.ndarray, np.ndarray]:
    runnum = np.asarray(
        arrays_by_logical_name["runnum"],
        dtype=float,
    )
    eventnum = np.asarray(
        arrays_by_logical_name["eventnum"],
        dtype=float,
    )
    valid = np.isfinite(runnum) & np.isfinite(eventnum)

    signature = np.empty(
        int(np.count_nonzero(valid)),
        dtype=[("runnum", "<i8"), ("eventnum", "<i8")],
    )
    signature["runnum"] = np.rint(runnum[valid]).astype(np.int64)
    signature["eventnum"] = np.rint(eventnum[valid]).astype(np.int64)
    return signature, valid


def mc_identity_signature(
    arrays_by_logical_name: Mapping[str, np.ndarray],
    decimals: int,
) -> Tuple[np.ndarray, np.ndarray]:
    logical_names = (
        "e_p",
        "e_theta",
        "e_phi",
        "p_p",
        "p_theta",
        "p_phi",
    )
    columns = [
        np.asarray(arrays_by_logical_name[name], dtype=float)
        for name in logical_names
    ]

    valid = np.ones(columns[0].size, dtype=bool)
    for column in columns:
        valid &= np.isfinite(column)
    # endfor

    scale = float(10 ** int(decimals))
    quantized = [
        np.rint(column[valid] * scale).astype(np.int64)
        for column in columns
    ]

    count = int(np.count_nonzero(valid))
    h1 = np.full(count, np.uint64(0x243F6A8885A308D3))
    h2 = np.full(count, np.uint64(0x13198A2E03707344))

    for index, column in enumerate(quantized):
        unsigned = column.view(np.uint64)
        mixed_a = uint64_mix(
            unsigned
            + np.uint64(
                (0x9E3779B97F4A7C15 * (index + 1))
                & 0xFFFFFFFFFFFFFFFF
            )
        )
        mixed_b = uint64_mix(
            unsigned
            ^ np.uint64(
                (0xD1B54A32D192ED03 * (index + 1))
                & 0xFFFFFFFFFFFFFFFF
            )
        )
        h1 = uint64_mix(h1 ^ mixed_a)
        h2 = uint64_mix(h2 + mixed_b)
    # endfor

    signature = np.empty(
        count,
        dtype=[("hash1", "<u8"), ("hash2", "<u8")],
    )
    signature["hash1"] = h1
    signature["hash2"] = h2
    return signature, valid


def photon_multiplicity_summary(
    multiplicity_audit: Mapping[str, object],
    period: PeriodConfig,
    sample_key: str,
    sample_label: str,
    args: argparse.Namespace,
) -> Dict[str, object]:
    """
    Serialize the exact event multiplicities computed during the identity
    prepass.

    The input already contains counts derived from the entry-level event
    signatures. It must not be passed through np.unique() a second time.
    """
    return {
        "period": period.key,
        "period_label": period.label,
        "sample": sample_key,
        "sample_label": sample_label,
        "identity_method": (
            "exact runnum + evnum"
            if sample_key == "data"
            else (
                "128-bit hash of rounded reconstructed electron/proton "
                f"kinematics ({args.mc_signature_decimals} decimals)"
            )
        ),
        "tree_entries_read": int(
            multiplicity_audit["entries_scanned"]
        ),
        "valid_identity_entries": int(
            multiplicity_audit["valid_identity_entries"]
        ),
        "unique_events": int(
            multiplicity_audit["unique_events"]
        ),
        "events_with_1_photon": int(
            multiplicity_audit["one_photon_events"]
        ),
        "events_with_2_photons": int(
            multiplicity_audit["two_photon_events"]
        ),
        "events_with_3_photons": int(
            multiplicity_audit["three_photon_events"]
        ),
        "events_with_4_photons": int(
            multiplicity_audit["four_photon_events"]
        ),
        "events_with_5_photons": int(
            multiplicity_audit["five_photon_events"]
        ),
        "events_with_6_photons": int(
            multiplicity_audit["six_photon_events"]
        ),
        "events_with_more_than_6_photons": int(
            multiplicity_audit["above_six_photon_events"]
        ),
        "maximum_photons_in_event": int(
            multiplicity_audit["maximum_multiplicity"]
        ),
        "photon_entries_accounted_for": int(
            multiplicity_audit["photon_entries_accounted_for"]
        ),
    }



def identity_logical_names(sample_key: str) -> Tuple[str, ...]:
    if sample_key == "data":
        return ("runnum", "eventnum")
    # endif
    return (
        "e_p",
        "e_theta",
        "e_phi",
        "p_p",
        "p_theta",
        "p_phi",
    )


def build_identity_signature(
    arrays_by_logical_name: Mapping[str, np.ndarray],
    sample_key: str,
    args: argparse.Namespace,
) -> Tuple[np.ndarray, np.ndarray]:
    if sample_key == "data":
        return data_identity_signature(arrays_by_logical_name)
    # endif
    return mc_identity_signature(
        arrays_by_logical_name,
        args.mc_signature_decimals,
    )


def scan_event_multiplicity(
    path: str,
    sample_key: str,
    resolved: Mapping[str, str],
    args: argparse.Namespace,
    entry_stop: int,
) -> Tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    Dict[str, int],
]:
    """
    Scan identity branches and classify events by photon-entry multiplicity.

    This function is diagnostic only. It does not define an event rejection.

    Returns:
      all unique event signatures,
      exact-one-photon event signatures,
      exact-two-photon event signatures,
      three-or-more-photon event signatures,
      entry/event accounting.
    """
    logical_names = identity_logical_names(sample_key)
    expressions = sorted(
        {resolved[logical_name] for logical_name in logical_names}
    )
    signature_chunks: List[np.ndarray] = []
    valid_identity_entries = 0
    entries_scanned = 0

    for arrays in uproot.iterate(
        f"{path}:{TREE_NAME}",
        expressions=expressions,
        step_size=args.step_size,
        entry_stop=entry_stop,
        library="np",
    ):
        chunk_size = len(next(iter(arrays.values()))) if arrays else 0
        entries_scanned += chunk_size

        logical_arrays = {
            logical_name: np.asarray(
                arrays[resolved[logical_name]],
                dtype=float,
            )
            for logical_name in logical_names
        }
        signatures, valid_identity = build_identity_signature(
            logical_arrays,
            sample_key,
            args,
        )
        valid_identity_entries += int(np.count_nonzero(valid_identity))
        signature_chunks.append(signatures)
    # endfor

    if signature_chunks:
        all_entry_signatures = np.concatenate(signature_chunks)
    else:
        if sample_key == "data":
            all_entry_signatures = np.empty(
                0,
                dtype=[("runnum", "<i8"), ("eventnum", "<i8")],
            )
        else:
            all_entry_signatures = np.empty(
                0,
                dtype=[("hash1", "<u8"), ("hash2", "<u8")],
            )
        # endif
    # endif

    if all_entry_signatures.size > 0:
        unique_signatures, multiplicities = np.unique(
            all_entry_signatures,
            return_counts=True,
        )
        one_photon_signatures = unique_signatures[multiplicities == 1]
        two_photon_signatures = unique_signatures[multiplicities == 2]
        three_plus_signatures = unique_signatures[multiplicities >= 3]
        three_plus_entries = int(
            np.sum(multiplicities[multiplicities >= 3])
        )
    else:
        unique_signatures = all_entry_signatures.copy()
        one_photon_signatures = all_entry_signatures.copy()
        two_photon_signatures = all_entry_signatures.copy()
        three_plus_signatures = all_entry_signatures.copy()
        multiplicities = np.empty(0, dtype=np.int64)
        three_plus_entries = 0
    # endif

    exact_event_counts = {
        multiplicity: int(
            np.count_nonzero(multiplicities == multiplicity)
        )
        for multiplicity in range(1, 7)
    }
    above_six_events = int(
        np.count_nonzero(multiplicities > 6)
    )
    accounted_entries = int(np.sum(multiplicities))

    accounting = {
        "entries_scanned": int(entries_scanned),
        "valid_identity_entries": int(valid_identity_entries),
        "unique_events": int(unique_signatures.size),
        "one_photon_events": exact_event_counts[1],
        "two_photon_events": exact_event_counts[2],
        "three_photon_events": exact_event_counts[3],
        "four_photon_events": exact_event_counts[4],
        "five_photon_events": exact_event_counts[5],
        "six_photon_events": exact_event_counts[6],
        "above_six_photon_events": above_six_events,
        "three_plus_photon_events": int(
            np.count_nonzero(multiplicities >= 3)
        ),
        "three_plus_photon_entries": three_plus_entries,
        "photon_entries_accounted_for": accounted_entries,
        "all_multiplicities_retained": False,
        "events_rejected_by_multiplicity": int(np.count_nonzero(multiplicities >= 3)),
        "entries_rejected_by_multiplicity": three_plus_entries,
        "maximum_multiplicity": (
            int(np.max(multiplicities))
            if multiplicities.size > 0
            else 0
        ),
    }
    return (
        unique_signatures,
        one_photon_signatures,
        two_photon_signatures,
        three_plus_signatures,
        accounting,
    )


def scan_shape_multiplicity_and_energy_ranking(
    path: str,
    sample_key: str,
    resolved: Mapping[str, str],
    args: argparse.Namespace,
    entry_stop: int,
) -> Tuple[np.ndarray, Dict[str, int], Dict[str, float]]:
    """
    Compute event multiplicities and one dense category code per ROOT entry.

    Category codes:
      0 = rejected or invalid identity;
      1 = exact-one-photon entry;
      2 = more energetic entry in an exact-two event;
      3 = less energetic entry in an exact-two event.

    The implementation is O(N log N), dominated by one stable lexicographic
    sort. It does not scan the complete inverse-index array once per event.
    """
    start_time = time.perf_counter()

    identity_names = list(identity_logical_names(sample_key))
    logical_names = list(
        dict.fromkeys(
            identity_names
            + ["p2_p", "p2_theta", "p2_phi"]
        )
    )
    expressions = sorted(
        {resolved[name] for name in logical_names}
    )

    signature_chunks: List[np.ndarray] = []
    energy_chunks: List[np.ndarray] = []
    theta_chunks: List[np.ndarray] = []
    phi_chunks: List[np.ndarray] = []
    ordinal_chunks: List[np.ndarray] = []
    valid_identity_entries = 0
    entries_scanned = 0

    io_start = time.perf_counter()
    for arrays in uproot.iterate(
        f"{path}:{TREE_NAME}",
        expressions=expressions,
        step_size=args.step_size,
        entry_stop=entry_stop,
        library="np",
    ):
        chunk_size = len(next(iter(arrays.values()))) if arrays else 0
        chunk_start = entries_scanned
        entries_scanned += chunk_size

        logical_arrays = {
            name: np.asarray(arrays[resolved[name]], dtype=float)
            for name in logical_names
        }
        signatures, valid_identity = build_identity_signature(
            logical_arrays,
            sample_key,
            args,
        )
        valid_identity_entries += int(np.count_nonzero(valid_identity))

        signature_chunks.append(signatures)
        energy_chunks.append(
            np.asarray(
                logical_arrays["p2_p"][valid_identity],
                dtype=float,
            )
        )
        theta_chunks.append(
            np.asarray(
                logical_arrays["p2_theta"][valid_identity],
                dtype=float,
            )
        )
        phi_chunks.append(
            np.asarray(
                logical_arrays["p2_phi"][valid_identity],
                dtype=float,
            )
        )
        ordinal_chunks.append(
            np.arange(
                chunk_start,
                chunk_start + chunk_size,
                dtype=np.int64,
            )[valid_identity]
        )
    # endfor
    io_seconds = elapsed_seconds(io_start)

    if signature_chunks:
        entry_signatures = np.concatenate(signature_chunks)
        entry_energies = np.concatenate(energy_chunks)
        entry_thetas = np.concatenate(theta_chunks)
        entry_phis = np.concatenate(phi_chunks)
        entry_ordinals = np.concatenate(ordinal_chunks)
    else:
        signature_dtype = (
            [("runnum", "<i8"), ("eventnum", "<i8")]
            if sample_key == "data"
            else [("hash1", "<u8"), ("hash2", "<u8")]
        )
        entry_signatures = np.empty(0, dtype=signature_dtype)
        entry_energies = np.empty(0, dtype=float)
        entry_thetas = np.empty(0, dtype=float)
        entry_phis = np.empty(0, dtype=float)
        entry_ordinals = np.empty(0, dtype=np.int64)
    # endif

    grouping_start = time.perf_counter()
    entry_category = np.zeros(entries_scanned, dtype=np.uint8)
    entry_diphoton_mass = np.full(
        entries_scanned,
        np.nan,
        dtype=np.float64,
    )

    if entry_signatures.size:
        unique_signatures, inverse, multiplicities = np.unique(
            entry_signatures,
            return_inverse=True,
            return_counts=True,
        )

        # Stable sort by group and then by original entry ordinal. Every group
        # becomes one contiguous segment.
        order = np.lexsort((entry_ordinals, inverse))
        sorted_groups = inverse[order]
        sorted_ordinals = entry_ordinals[order]
        sorted_energies = entry_energies[order]
        sorted_thetas = entry_thetas[order]
        sorted_phis = entry_phis[order]

        group_starts = np.empty(multiplicities.size, dtype=np.int64)
        group_starts[0] = 0
        if multiplicities.size > 1:
            np.cumsum(
                multiplicities[:-1],
                out=group_starts[1:],
            )
        # endif

        one_groups = np.flatnonzero(multiplicities == 1)
        two_groups = np.flatnonzero(multiplicities == 2)

        if one_groups.size:
            one_positions = group_starts[one_groups]
            one_ordinals = sorted_ordinals[one_positions]
            entry_category[one_ordinals] = 1
        # endif

        equal_energy_ties = 0
        nonfinite_pairs = 0

        if two_groups.size:
            first_positions = group_starts[two_groups]
            second_positions = first_positions + 1

            first_ordinals = sorted_ordinals[first_positions]
            second_ordinals = sorted_ordinals[second_positions]
            first_energies = sorted_energies[first_positions]
            second_energies = sorted_energies[second_positions]

            first_finite = np.isfinite(first_energies)
            second_finite = np.isfinite(second_energies)
            both_finite = first_finite & second_finite
            unequal = both_finite & (
                first_energies != second_energies
            )

            first_is_high = np.zeros(two_groups.size, dtype=bool)
            first_is_high[unequal] = (
                first_energies[unequal]
                > second_energies[unequal]
            )

            only_first_finite = first_finite & ~second_finite
            only_second_finite = second_finite & ~first_finite
            first_is_high[only_first_finite] = True
            first_is_high[only_second_finite] = False

            unresolved = ~unequal & ~only_first_finite & ~only_second_finite
            # Deterministic entry-order tie break for equal or doubly
            # non-finite energies.
            first_is_high[unresolved] = (
                first_ordinals[unresolved]
                < second_ordinals[unresolved]
            )

            equal_energy_ties = int(
                np.count_nonzero(
                    both_finite
                    & (first_energies == second_energies)
                )
            )
            nonfinite_pairs = int(
                np.count_nonzero(~both_finite)
            )

            high_ordinals = np.where(
                first_is_high,
                first_ordinals,
                second_ordinals,
            )
            low_ordinals = np.where(
                first_is_high,
                second_ordinals,
                first_ordinals,
            )

            entry_category[high_ordinals] = 2
            entry_category[low_ordinals] = 3

            first_thetas = sorted_thetas[first_positions]
            second_thetas = sorted_thetas[second_positions]
            first_phis = sorted_phis[first_positions]
            second_phis = sorted_phis[second_positions]

            opening_cosine = (
                np.cos(first_thetas) * np.cos(second_thetas)
                + np.sin(first_thetas)
                * np.sin(second_thetas)
                * np.cos(first_phis - second_phis)
            )
            opening_cosine = np.clip(
                opening_cosine,
                -1.0,
                1.0,
            )
            mass_squared = (
                2.0
                * first_energies
                * second_energies
                * (1.0 - opening_cosine)
            )
            valid_mass = (
                np.isfinite(mass_squared)
                & (mass_squared >= 0.0)
            )
            pair_mass = np.full(
                two_groups.size,
                np.nan,
                dtype=np.float64,
            )
            pair_mass[valid_mass] = np.sqrt(
                mass_squared[valid_mass]
            )
            entry_diphoton_mass[high_ordinals] = pair_mass
            entry_diphoton_mass[low_ordinals] = pair_mass
        # endif
    else:
        unique_signatures = entry_signatures.copy()
        multiplicities = np.empty(0, dtype=np.int64)
        equal_energy_ties = 0
        nonfinite_pairs = 0
    # endif

    grouping_seconds = elapsed_seconds(grouping_start)

    exact_counts = {
        value: int(np.count_nonzero(multiplicities == value))
        for value in range(1, 7)
    }
    three_plus_entries = int(
        np.sum(multiplicities[multiplicities >= 3])
    )

    accounting = {
        "entries_scanned": int(entries_scanned),
        "valid_identity_entries": int(valid_identity_entries),
        "unique_events": int(unique_signatures.size),
        "one_photon_events": exact_counts[1],
        "two_photon_events": exact_counts[2],
        "three_photon_events": exact_counts[3],
        "four_photon_events": exact_counts[4],
        "five_photon_events": exact_counts[5],
        "six_photon_events": exact_counts[6],
        "above_six_photon_events": int(
            np.count_nonzero(multiplicities > 6)
        ),
        "three_plus_photon_events": int(
            np.count_nonzero(multiplicities >= 3)
        ),
        "three_plus_photon_entries": three_plus_entries,
        "photon_entries_accounted_for": int(np.sum(multiplicities)),
        "all_multiplicities_retained": False,
        "events_rejected_by_multiplicity": int(
            np.count_nonzero(multiplicities >= 3)
        ),
        "entries_rejected_by_multiplicity": three_plus_entries,
        "maximum_multiplicity": (
            int(np.max(multiplicities))
            if multiplicities.size
            else 0
        ),
        "two_photon_more_energetic_entries": int(
            np.count_nonzero(entry_category == 2)
        ),
        "two_photon_less_energetic_entries": int(
            np.count_nonzero(entry_category == 3)
        ),
        "two_photon_equal_energy_ties": int(equal_energy_ties),
        "two_photon_nonfinite_energy_pairs": int(nonfinite_pairs),
        "two_photon_ranking_branch": "p2_p",
    }

    timings = {
        "multiplicity_prepass_total_seconds": elapsed_seconds(
            start_time
        ),
        "multiplicity_prepass_root_io_seconds": io_seconds,
        "multiplicity_grouping_seconds": grouping_seconds,
    }
    maybe_log_timing(
        args,
        f"{Path(path).name} multiplicity/ranking prepass",
        timings["multiplicity_prepass_total_seconds"],
    )

    return (
        entry_category,
        entry_diphoton_mass,
        accounting,
        timings,
    )


def multiplicity_entry_masks(
    arrays_by_logical_name: Mapping[str, np.ndarray],
    sample_key: str,
    one_photon_signatures: np.ndarray,
    two_photon_signatures: np.ndarray,
    three_plus_signatures: np.ndarray,
    args: argparse.Namespace,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Return full-length masks for exact-one, exact-two, >=3, and all valid
    event identities. No multiplicity category is rejected.
    """
    signatures, valid_identity = build_identity_signature(
        arrays_by_logical_name,
        sample_key,
        args,
    )

    one_valid = np.zeros(signatures.size, dtype=bool)
    two_valid = np.zeros(signatures.size, dtype=bool)
    three_plus_valid = np.zeros(signatures.size, dtype=bool)

    if signatures.size > 0:
        if one_photon_signatures.size > 0:
            one_valid = np.isin(signatures, one_photon_signatures)
        # endif
        if two_photon_signatures.size > 0:
            two_valid = np.isin(signatures, two_photon_signatures)
        # endif
        if three_plus_signatures.size > 0:
            three_plus_valid = np.isin(
                signatures,
                three_plus_signatures,
            )
        # endif
    # endif

    one_full = np.zeros(valid_identity.size, dtype=bool)
    two_full = np.zeros(valid_identity.size, dtype=bool)
    three_plus_full = np.zeros(valid_identity.size, dtype=bool)

    one_full[valid_identity] = one_valid
    two_full[valid_identity] = two_valid
    three_plus_full[valid_identity] = three_plus_valid

    return (
        one_full,
        two_full,
        three_plus_full,
        valid_identity,
    )



def process_data_shape_sample(
    period: PeriodConfig,
    sample_key: str,
    sample_label: str,
    period_attribute: str,
    args_dict: Mapping[str, object],
) -> Tuple[str, str, Dict[str, object]]:
    """
    Fill category-specific before/after populations.

    Exact-one:
        after == before.

    Exact-two:
        after additionally requires
        PI0_MASS_WINDOW[0] < M_gamma_gamma < PI0_MASS_WINDOW[1].
    """
    total_start = time.perf_counter()
    args = argparse.Namespace(**args_dict)
    path = getattr(period, period_attribute)
    total_entries, available_keys = require_tree(path)

    shape_names = [
        *raw_shape_branch_keys(),
        # Required internally for exact-two photon energy ranking and Mgg.
        # These need not be displayed to still be required inputs.
        "p2_p",
        "p2_theta",
        "p2_phi",
        "open_angle_ep2",
    ]
    logical_names = list(
        dict.fromkeys(
            shape_names + list(identity_logical_names(sample_key))
        )
    )
    resolved = resolve_branches_from_keys(
        path,
        available_keys,
        logical_names,
    )
    expressions = sorted(set(resolved.values()))

    entry_stop = (
        min(total_entries, int(args.max_events))
        if args.max_events is not None
        else total_entries
    )

    (
        entry_category,
        entry_diphoton_mass,
        multiplicity_audit,
        prepass_timings,
    ) = scan_shape_multiplicity_and_energy_ranking(
        path,
        sample_key,
        resolved,
        args,
        entry_stop,
    )

    multiplicity = photon_multiplicity_summary(
        multiplicity_audit,
        period,
        sample_key,
        sample_label,
        args,
    )

    category_specs = (
        ("one_photon", 1),
        ("two_photon_more_energetic", 2),
        ("two_photon_less_energetic", 3),
    )
    category_variables = {
        key: display_variables_for_category(key)
        for key, _code in category_specs
    }

    before_histograms = {
        key: empty_stage_histograms(category_variables[key])
        for key, _code in category_specs
    }
    after_histograms = {
        key: empty_stage_histograms(category_variables[key])
        for key, _code in category_specs
    }

    cutflow_counts = {
        key: {
            "finite_required_branches": 0,
            "open_angle_ep2_gt_5_deg": 0,
            "after_pi0_mass_window": 0,
        }
        for key, _code in category_specs
    }

    pi0_mass_diagnostics = {
        key: {
            "window_low_GeV": PI0_MASS_WINDOW[0],
            "window_high_GeV": PI0_MASS_WINDOW[1],
            "before_entries": 0,
            "after_entries": 0,
            "fraction_in_window": math.nan,
        }
        for key in (
            "two_photon_more_energetic",
            "two_photon_less_energetic",
        )
    }

    entries_read = 0
    histogram_start = time.perf_counter()

    log(
        f"{period.label} {sample_label}: one="
        f"{multiplicity_audit['one_photon_events']:,}, two="
        f"{multiplicity_audit['two_photon_events']:,}, rejected >=3="
        f"{multiplicity_audit['three_plus_photon_events']:,}; "
        f"exact-two after cut: "
        f"{PI0_MASS_WINDOW[0]:.2f}<Mgg<"
        f"{PI0_MASS_WINDOW[1]:.2f} GeV"
    )

    for arrays in uproot.iterate(
        f"{path}:{TREE_NAME}",
        expressions=expressions,
        step_size=args.step_size,
        entry_stop=entry_stop,
        library="np",
    ):
        chunk_size = len(next(iter(arrays.values()))) if arrays else 0
        chunk_start = entries_read
        chunk_stop = chunk_start + chunk_size
        chunk_category = entry_category[chunk_start:chunk_stop]
        entries_read = chunk_stop

        logical_arrays = {
            name: np.asarray(arrays[resolved[name]], dtype=float)
            for name in logical_names
        }
        logical_arrays["diphoton_mass"] = (
            entry_diphoton_mass[chunk_start:chunk_stop]
        )

        open_angle_finite = np.isfinite(
            logical_arrays["open_angle_ep2"]
        )
        open_angle_mask = (
            logical_arrays["open_angle_ep2"] > 5.0
        )

        diphoton_mass = logical_arrays["diphoton_mass"]
        pi0_window_mask = (
            np.isfinite(diphoton_mass)
            & (diphoton_mass > PI0_MASS_WINDOW[0])
            & (diphoton_mass < PI0_MASS_WINDOW[1])
        )

        for category_key, category_code in category_specs:
            variables = category_variables[category_key]
            identity_mask = chunk_category == category_code

            common_finite = open_angle_finite.copy()
            for variable in variables:
                common_finite &= np.isfinite(
                    logical_arrays[variable.key]
                )
            # endfor

            finite_mask = identity_mask & common_finite
            cutflow_counts[category_key][
                "finite_required_branches"
            ] += int(np.count_nonzero(finite_mask))

            before_mask = finite_mask & open_angle_mask
            before_count = int(np.count_nonzero(before_mask))
            cutflow_counts[category_key][
                "open_angle_ep2_gt_5_deg"
            ] += before_count

            if category_key == "one_photon":
                after_mask = before_mask
            else:
                after_mask = before_mask & pi0_window_mask
            # endif

            after_count = int(np.count_nonzero(after_mask))
            cutflow_counts[category_key][
                "after_pi0_mass_window"
            ] += after_count

            update_stage_histograms(
                before_histograms[category_key],
                logical_arrays,
                before_mask,
                variables,
            )
            update_stage_histograms(
                after_histograms[category_key],
                logical_arrays,
                after_mask,
                variables,
            )

            if category_key in pi0_mass_diagnostics:
                diagnostic = pi0_mass_diagnostics[category_key]
                diagnostic["before_entries"] += before_count
                diagnostic["after_entries"] += after_count
            # endif
        # endfor
    # endfor

    histogram_seconds = elapsed_seconds(histogram_start)

    if entries_read != multiplicity_audit["entries_scanned"]:
        raise RuntimeError(
            f"{period.label} {sample_label}: prepass/shape-pass entry "
            f"count mismatch ({multiplicity_audit['entries_scanned']:,} "
            f"versus {entries_read:,})."
        )
    # endif

    expected_two = int(multiplicity_audit["two_photon_events"])
    if (
        np.count_nonzero(entry_category == 2) != expected_two
        or np.count_nonzero(entry_category == 3) != expected_two
    ):
        raise RuntimeError(
            f"{period.label} {sample_label}: dense exact-two category "
            "counts do not match exact-two event count."
        )
    # endif

    for diagnostic in pi0_mass_diagnostics.values():
        before_entries = int(diagnostic["before_entries"])
        after_entries = int(diagnostic["after_entries"])
        diagnostic["fraction_in_window"] = (
            after_entries / before_entries
            if before_entries > 0
            else math.nan
        )
    # endfor

    finalized_before = {
        key: finalize_stage_histograms(
            before_histograms[key],
            category_variables[key],
        )
        for key, _code in category_specs
    }
    finalized_after = {
        key: finalize_stage_histograms(
            after_histograms[key],
            category_variables[key],
        )
        for key, _code in category_specs
    }

    total_seconds = elapsed_seconds(total_start)
    timings = {
        **prepass_timings,
        "shape_histogram_pass_seconds": histogram_seconds,
        "shape_sample_total_seconds": total_seconds,
    }
    maybe_log_timing(
        args,
        f"{period.label} {sample_label} histogram pass",
        histogram_seconds,
    )
    maybe_log_timing(
        args,
        f"{period.label} {sample_label} total",
        total_seconds,
    )

    payload = {
        "period": period.key,
        "period_label": period.label,
        "sample": sample_key,
        "sample_label": sample_label,
        "path": path,
        "tree_entries": total_entries,
        "entries_read": entries_read,
        "resolved_branches": resolved,
        "multiplicity_audit": {
            "condition": (
                "retain exact-one events; split exact-two events into "
                "more-/less-energetic entries; reject multiplicity >=3"
            ),
            **multiplicity_audit,
        },
        "cutflow": {
            "tree_entries_read": entries_read,
            "valid_identity_entries": int(
                multiplicity_audit["valid_identity_entries"]
            ),
            **cutflow_counts,
        },
        "photon_multiplicity": multiplicity,
        "pi0_mass_diagnostics": pi0_mass_diagnostics,
        "timings": timings,
    }

    for category_key, _code in category_specs:
        payload[category_key] = {
            "before_histograms": finalized_before[category_key],
            "after_histograms": finalized_after[category_key],
        }
    # endfor

    return period.key, sample_key, payload


def plot_before_after_shape_canvas(
    path: Path,
    period_label: str,
    multiplicity_label: str,
    sample_payloads: Mapping[str, Mapping[str, object]],
    multiplicity_key: str,
    variables: Sequence[FitVariable],
    canvas_description: str,
    args: argparse.Namespace,
) -> None:
    """
    Write a 2x5 category-specific shape-comparison canvas.

    Top row:
        before the exact-two M_gamma_gamma requirement.

    Bottom row:
        after the exact-two M_gamma_gamma requirement.

    Exact-one events have identical before and after rows.
    """
    path.parent.mkdir(parents=True, exist_ok=True)

    sample_colors = {
        "data": "black",
        "dvcsgen": "tab:blue",
        "aaogen": "tab:red",
    }

    display_variables = display_variables_for_category(
        multiplicity_key
    )
    if len(display_variables) != 5:
        raise RuntimeError(
            f"{multiplicity_key}: expected exactly five display variables."
        )
    # endif

    fig, axes = plt.subplots(
        2,
        5,
        figsize=(23.0, 9.4),
        squeeze=False,
    )

    row_specs = (
        ("before_histograms", "before"),
        ("after_histograms", "after"),
    )

    for row_index, (stage_name, row_label) in enumerate(row_specs):
        for column_index, variable in enumerate(display_variables):
            axis = axes[row_index, column_index]

            for sample_key, sample_label, _ in RAW_SHAPE_SAMPLE_INFO:
                histogram = sample_payloads[sample_key][
                    multiplicity_key
                ][stage_name][variable.key]
                edges = np.asarray(
                    histogram["edges"],
                    dtype=float,
                )
                normalized = np.asarray(
                    histogram["unit_area"],
                    dtype=float,
                )

                label = (
                    f"{sample_label} "
                    f"({histogram['in_range_entries']:,})"
                )

                if (
                    variable.key == "diphoton_mass"
                    and sample_key == "data"
                    and row_index == 0
                ):
                    diagnostic = sample_payloads[sample_key][
                        "pi0_mass_diagnostics"
                    ][multiplicity_key]
                    fraction = float(
                        diagnostic["fraction_in_window"]
                    )
                    if np.isfinite(fraction):
                        label += (
                            f"; in window={100.0 * fraction:.1f}%"
                        )
                    # endif
                # endif

                axis.stairs(
                    normalized,
                    edges,
                    color=sample_colors[sample_key],
                    linewidth=(
                        1.45 if sample_key == "data" else 1.25
                    ),
                    linestyle="-",
                    label=label,
                    zorder=(
                        4 if sample_key == "data" else
                        3 if sample_key == "dvcsgen" else 2
                    ),
                )
            # endfor

            if variable.key == "diphoton_mass":
                axis.axvspan(
                    PI0_MASS_WINDOW[0],
                    PI0_MASS_WINDOW[1],
                    alpha=0.15,
                    label=(
                        r"$\pi^0$ mass requirement "
                        f"({PI0_MASS_WINDOW[0]:.2f}-"
                        f"{PI0_MASS_WINDOW[1]:.2f} GeV)"
                    ),
                )
                axis.axvline(
                    0.1349768,
                    linewidth=1.0,
                    linestyle="--",
                )
            # endif

            axis.set_xlim(variable.low, variable.high)
            axis.set_xlabel(variable.label)
            axis.set_ylabel("fraction / bin")
            axis.set_title(
                f"{row_label}: {variable.label}"
            )
            axis.grid(axis="y", alpha=0.25)

            if args.shape_log_y:
                positive_values: List[float] = []
                for sample_key, _, _ in RAW_SHAPE_SAMPLE_INFO:
                    values = np.asarray(
                        sample_payloads[sample_key][
                            multiplicity_key
                        ][stage_name][variable.key]["unit_area"],
                        dtype=float,
                    )
                    positive_values.extend(
                        values[values > 0.0].tolist()
                    )
                # endfor
                if positive_values:
                    axis.set_yscale("log")
                    axis.set_ylim(
                        bottom=max(
                            0.5 * min(positive_values),
                            1.0e-8,
                        )
                    )
                # endif
            else:
                axis.set_ylim(bottom=0.0)
            # endif
        # endfor
    # endfor

    handles, labels = axes[0, 0].get_legend_handles_labels()
    if multiplicity_key != "one_photon":
        mass_handles, mass_labels = axes[0, 4].get_legend_handles_labels()
        handles = [*handles, *mass_handles]
        labels = [*labels, *mass_labels]
    # endif

    if multiplicity_key == "one_photon":
        selection_text = (
            r"open_angle_ep2 $>5^\circ$; no additional after cut; "
            "rows are intentionally identical"
        )
    else:
        selection_text = (
            r"open_angle_ep2 $>5^\circ$; bottom row additionally requires "
            f"{PI0_MASS_WINDOW[0]:.2f}"
            r"$<M_{\gamma\gamma}<$"
            f"{PI0_MASS_WINDOW[1]:.2f} GeV"
        )
    # endif

    fig.suptitle(
        f"Unit-area epgamma shape comparisons: "
        f"{period_label}, {multiplicity_label}\n"
        f"{selection_text}; each sample independently normalized",
        fontsize=15,
        y=0.992,
    )
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.88),
        ncol=min(5, max(1, len(labels))),
        frameon=False,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.81))
    fig.savefig(
        path,
        dpi=180,
        bbox_inches="tight",
    )
    plt.close(fig)


def serializable_histograms(
    histograms: Mapping[str, Mapping[str, object]],
) -> Dict[str, Dict[str, object]]:
    output: Dict[str, Dict[str, object]] = {}
    for variable_key, histogram in histograms.items():
        output[variable_key] = {
            key: (
                value.tolist()
                if isinstance(value, np.ndarray)
                else value
            )
            for key, value in histogram.items()
        }
    # endfor
    return output



def multiplicity_categories() -> Tuple[Tuple[str, str], ...]:
    return (
        ("1", "events_with_1_photon"),
        ("2", "events_with_2_photons"),
        ("3", "events_with_3_photons"),
        ("4", "events_with_4_photons"),
        ("5", "events_with_5_photons"),
        ("6", "events_with_6_photons"),
        (">6", "events_with_more_than_6_photons"),
    )


def plot_period_photon_multiplicity(
    path: Path,
    period_label: str,
    rows: Sequence[Mapping[str, object]],
) -> None:
    categories = multiplicity_categories()
    labels = [label for label, _ in categories]
    x = np.arange(len(categories), dtype=float)

    sample_rows = {
        str(row["sample"]): row
        for row in rows
    }
    offsets = {
        "data": -0.24,
        "dvcsgen": 0.0,
        "aaogen": 0.24,
    }
    width = 0.22

    fig, axes = plt.subplots(2, 1, figsize=(12, 10), sharex=True)

    for sample_key, sample_label, _ in RAW_SHAPE_SAMPLE_INFO:
        if sample_key not in sample_rows:
            raise RuntimeError(
                f"{period_label}: missing multiplicity row for "
                f"{sample_key}."
            )
        # endif

        row = sample_rows[sample_key]
        totals = np.asarray(
            [int(row[field_name]) for _, field_name in categories],
            dtype=float,
        )
        unique_events = int(row["unique_events"])
        percentages = (
            100.0 * totals / unique_events
            if unique_events > 0
            else np.zeros_like(totals)
        )

        axes[0].bar(
            x + offsets[sample_key],
            totals,
            width=width,
            label=sample_label,
        )
        axes[1].bar(
            x + offsets[sample_key],
            percentages,
            width=width,
            label=sample_label,
        )
    # endfor

    axes[0].set_ylabel("Unique events")
    axes[0].set_title("Total event count")
    axes[0].grid(axis="y", alpha=0.25)
    axes[0].legend()

    axes[1].set_ylabel("Percentage of unique events (%)")
    axes[1].set_xlabel("Number of reconstructed-photon entries")
    axes[1].set_title("Percentage of unique events")
    axes[1].set_xticks(x, labels)
    axes[1].set_ylim(bottom=0.0)
    axes[1].grid(axis="y", alpha=0.25)
    axes[1].legend()

    fig.suptitle(
        f"{period_label}: epgamma photon-entry multiplicity"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_all_periods_photon_multiplicity(
    path: Path,
    rows: Sequence[Mapping[str, object]],
) -> None:
    categories = multiplicity_categories()
    labels = [label for label, _ in categories]
    x = np.arange(len(categories), dtype=float)

    offsets = {
        "data": -0.24,
        "dvcsgen": 0.0,
        "aaogen": 0.24,
    }
    width = 0.22

    selected_periods = [
        period
        for period in PERIODS
        if any(str(row["period"]) == period.key for row in rows)
    ]

    fig, axes = plt.subplots(
        len(selected_periods),
        2,
        figsize=(18, 4.2 * len(selected_periods)),
        squeeze=False,
    )

    for period_index, period in enumerate(selected_periods):
        period_rows = {
            str(row["sample"]): row
            for row in rows
            if str(row["period"]) == period.key
        }

        total_axis = axes[period_index, 0]
        percent_axis = axes[period_index, 1]

        for sample_key, sample_label, _ in RAW_SHAPE_SAMPLE_INFO:
            if sample_key not in period_rows:
                raise RuntimeError(
                    f"{period.label}: missing multiplicity row for "
                    f"{sample_key}."
                )
            # endif

            row = period_rows[sample_key]
            totals = np.asarray(
                [int(row[field_name]) for _, field_name in categories],
                dtype=float,
            )
            unique_events = int(row["unique_events"])
            percentages = (
                100.0 * totals / unique_events
                if unique_events > 0
                else np.zeros_like(totals)
            )

            total_axis.bar(
                x + offsets[sample_key],
                totals,
                width=width,
                label=sample_label,
            )
            percent_axis.bar(
                x + offsets[sample_key],
                percentages,
                width=width,
                label=sample_label,
            )
        # endfor

        total_axis.set_title(f"{period.label}: total events")
        total_axis.set_ylabel("Unique events")
        total_axis.set_xticks(x, labels)
        total_axis.grid(axis="y", alpha=0.25)

        percent_axis.set_title(f"{period.label}: percentage")
        percent_axis.set_ylabel("Percentage of unique events (%)")
        percent_axis.set_xticks(x, labels)
        percent_axis.set_ylim(bottom=0.0)
        percent_axis.grid(axis="y", alpha=0.25)

        if period_index == 0:
            total_axis.legend(fontsize=8)
            percent_axis.legend(fontsize=8)
        # endif
    # endfor

    axes[-1, 0].set_xlabel("Number of reconstructed-photon entries")
    axes[-1, 1].set_xlabel("Number of reconstructed-photon entries")

    fig.suptitle(
        "epgamma photon-entry multiplicity by run period"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(path, dpi=180)
    plt.close(fig)



def write_photon_multiplicity_outputs(
    data_root: Path,
    rows: Sequence[Mapping[str, object]],
) -> None:
    ordered_fields = [
        "period",
        "period_label",
        "sample",
        "sample_label",
        "identity_method",
        "tree_entries_read",
        "valid_identity_entries",
        "unique_events",
        "events_with_1_photon",
        "events_with_2_photons",
        "events_with_3_photons",
        "events_with_4_photons",
        "events_with_5_photons",
        "events_with_6_photons",
        "events_with_more_than_6_photons",
        "maximum_photons_in_event",
        "photon_entries_accounted_for",
    ]

    with open(
        data_root / "photon_multiplicity_summary.csv",
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=ordered_fields)
        writer.writeheader()
        writer.writerows(rows)
    # endwith

    with open(
        data_root / "photon_multiplicity_summary.json",
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(list(rows), handle, indent=2)
    # endwith

    plot_dir = data_root / "photon_multiplicity"
    plot_dir.mkdir(parents=True, exist_ok=True)

    for period in PERIODS:
        period_rows = [
            row
            for row in rows
            if str(row["period"]) == period.key
        ]
        if not period_rows:
            continue
        # endif

        plot_period_photon_multiplicity(
            plot_dir / f"{period.key}_photon_multiplicity.png",
            period.label,
            period_rows,
        )
    # endfor

    plot_all_periods_photon_multiplicity(
        plot_dir / "all_periods_photon_multiplicity.png",
        rows,
    )
    # endwith



@dataclass
class AdvancedVariableFit:
    success: bool
    message: str
    variable: str
    model_name: str
    f_pi0: float
    f_pi0_err: float
    independent_f_pi0: float
    independent_f_pi0_err: float
    deviance: float
    ndf: int
    nuisance_count: int
    model_counts: np.ndarray
    dvcs_component_counts: np.ndarray
    pi0_component_counts: np.ndarray
    raw_dvcs_component_counts: np.ndarray
    raw_pi0_component_counts: np.ndarray
    morphed_dvcs_probability: np.ndarray
    morphed_pi0_probability: np.ndarray
    raw_dvcs_probability: np.ndarray
    raw_pi0_probability: np.ndarray
    dvcs_nuisance: Dict[str, float]
    pi0_nuisance: Dict[str, float]
    divider_dvcs: float
    divider_pi0: float
    fit_total: float
    data_total: float
    raw_shared_f_pi0: float
    raw_shared_f_pi0_err: float
    raw_independent_f_pi0: float
    raw_independent_f_pi0_err: float
    raw_shared_deviance: float
    raw_independent_deviance: float
    raw_shared_model_counts: np.ndarray
    raw_shared_dvcs_component_counts: np.ndarray
    raw_shared_pi0_component_counts: np.ndarray


@dataclass
class AdvancedSharedFit:
    success: bool
    message: str
    f_pi0: float
    f_pi0_err: float
    deviance: float
    ndf: int
    variable_results: Dict[str, AdvancedVariableFit]
    coordinate_iterations: int
    multiplicity_key: str
    raw_f_pi0: float
    raw_f_pi0_err: float
    raw_deviance: float
    raw_ndf: int


def normalized_probability(counts: np.ndarray) -> np.ndarray:
    values = np.asarray(counts, dtype=float)
    values = np.clip(values, 0.0, None)
    total = float(np.sum(values))
    if total <= 0.0:
        return np.zeros_like(values, dtype=float)
    # endif
    return values / total


def poisson_deviance(
    observed: np.ndarray,
    expected: np.ndarray,
) -> float:
    observed = np.asarray(observed, dtype=float)
    expected = np.clip(
        np.asarray(expected, dtype=float),
        1.0e-12,
        None,
    )
    term = expected - observed
    positive = observed > 0.0
    term[positive] += (
        observed[positive]
        * np.log(observed[positive] / expected[positive])
    )
    return float(2.0 * np.sum(term))


def weighted_mean(
    probability: np.ndarray,
    centers: np.ndarray,
) -> float:
    probability = normalized_probability(probability)
    if np.sum(probability) <= 0.0:
        return math.nan
    # endif
    return float(np.sum(probability * centers))


def deposit_probability(
    source_probability: np.ndarray,
    mapped_centers: np.ndarray,
    output_centers: np.ndarray,
) -> np.ndarray:
    """
    Deposit probability linearly on the output-bin centers.

    Probability mapped outside the configured plotting range is accumulated
    into the nearest edge bin. This preserves template normalization.
    """
    source = normalized_probability(source_probability)
    output = np.zeros_like(output_centers, dtype=float)

    if source.size == 0:
        return output
    # endif

    indices = np.searchsorted(output_centers, mapped_centers)
    indices = np.clip(indices, 1, output_centers.size - 1)
    left = indices - 1
    right = indices

    x_left = output_centers[left]
    x_right = output_centers[right]
    denominator = np.maximum(x_right - x_left, 1.0e-12)
    right_weight = np.clip(
        (mapped_centers - x_left) / denominator,
        0.0,
        1.0,
    )
    left_weight = 1.0 - right_weight

    below = mapped_centers <= output_centers[0]
    above = mapped_centers >= output_centers[-1]
    left[below] = 0
    right[below] = 0
    left_weight[below] = 1.0
    right_weight[below] = 0.0
    left[above] = output_centers.size - 1
    right[above] = output_centers.size - 1
    left_weight[above] = 1.0
    right_weight[above] = 0.0

    np.add.at(output, left, source * left_weight)
    np.add.at(output, right, source * right_weight)
    return normalized_probability(output)


def asymmetric_gaussian_broadening(
    probability: np.ndarray,
    sigma_left_bins: float,
    sigma_right_bins: float,
) -> np.ndarray:
    probability = normalized_probability(probability)
    if probability.size == 0:
        return probability
    # endif

    mode_index = int(np.argmax(probability))
    left_component = probability.copy()
    right_component = probability.copy()
    left_component[mode_index + 1:] = 0.0
    right_component[:mode_index + 1] = 0.0

    if sigma_left_bins > 1.0e-8:
        left_component = gaussian_filter1d(
            left_component,
            sigma=float(sigma_left_bins),
            mode="constant",
            cval=0.0,
            truncate=4.0,
        )
    # endif
    if sigma_right_bins > 1.0e-8:
        right_component = gaussian_filter1d(
            right_component,
            sigma=float(sigma_right_bins),
            mode="constant",
            cval=0.0,
            truncate=4.0,
        )
    # endif

    return normalized_probability(left_component + right_component)


def affine_morph_probability(
    raw_probability: np.ndarray,
    centers: np.ndarray,
    shift_bins: float,
    log_scale: float,
    sigma_left_bins: float,
    sigma_right_bins: float,
    use_log_coordinate: bool,
) -> np.ndarray:
    raw = normalized_probability(raw_probability)
    if np.sum(raw) <= 0.0:
        return raw
    # endif

    if use_log_coordinate:
        positive_widths = np.diff(centers)
        epsilon = max(
            0.5 * float(np.median(positive_widths)),
            1.0e-6,
        )
        coordinate = np.log(np.maximum(centers + epsilon, epsilon))
        coordinate_step = float(np.median(np.diff(coordinate)))
        pivot = float(np.sum(raw * coordinate))
        mapped_coordinate = (
            pivot
            + math.exp(float(log_scale)) * (coordinate - pivot)
            + float(shift_bins) * coordinate_step
        )
        mapped_centers = np.exp(mapped_coordinate) - epsilon
    else:
        bin_width = float(np.median(np.diff(centers)))
        pivot = float(np.sum(raw * centers))
        mapped_centers = (
            pivot
            + math.exp(float(log_scale)) * (centers - pivot)
            + float(shift_bins) * bin_width
        )
    # endif

    deposited = deposit_probability(
        raw,
        mapped_centers,
        centers,
    )
    return asymmetric_gaussian_broadening(
        deposited,
        sigma_left_bins,
        sigma_right_bins,
    )


def simple_shift_smear_probability(
    raw_probability: np.ndarray,
    centers: np.ndarray,
    shift_bins: float,
    sigma_bins: float,
    use_log_coordinate: bool,
) -> np.ndarray:
    """
    Coherently shift and Gaussian-broaden an entire template.

    This mirrors the transport used by plot_exclusivity_data_dvcs_pi0_mc.py:
    every source-bin probability is transported into every target bin using
    integrated Gaussian probabilities. No interpolation cutoff is introduced.

    Parameters remain expressed in transformed-coordinate bin units so the
    existing --fit-max-shift-bins/--fit-max-smear-bins interface is preserved.
    """
    raw = normalized_probability(raw_probability)
    if np.sum(raw) <= 0.0 or centers.size < 2:
        return raw
    # endif

    linear_bin_width = float(np.median(np.diff(centers)))
    edges = np.concatenate(
        [
            np.asarray(
                [centers[0] - 0.5 * linear_bin_width],
                dtype=float,
            ),
            0.5 * (centers[:-1] + centers[1:]),
            np.asarray(
                [centers[-1] + 0.5 * linear_bin_width],
                dtype=float,
            ),
        ]
    )

    if use_log_coordinate:
        epsilon = max(
            0.25 * linear_bin_width,
            1.0e-8,
        )
        source_coordinate = np.log(
            np.maximum(centers + epsilon, 1.0e-15)
        )
        lower_coordinate = np.log(
            np.maximum(edges[:-1] + epsilon, 1.0e-15)
        )
        upper_coordinate = np.log(
            np.maximum(edges[1:] + epsilon, 1.0e-15)
        )
        coordinate_step = float(
            np.median(np.diff(source_coordinate))
        )
    else:
        source_coordinate = centers
        lower_coordinate = edges[:-1]
        upper_coordinate = edges[1:]
        coordinate_step = linear_bin_width
    # endif

    shift = float(shift_bins) * coordinate_step
    sigma = max(float(sigma_bins), 0.0) * coordinate_step

    if sigma <= 1.0e-12:
        mapped_coordinate = source_coordinate + shift

        if use_log_coordinate:
            epsilon = max(
                0.25 * linear_bin_width,
                1.0e-8,
            )
            mapped_centers = np.exp(mapped_coordinate) - epsilon
        else:
            mapped_centers = mapped_coordinate
        # endif

        target, _ = np.histogram(
            mapped_centers,
            bins=edges,
            weights=raw,
        )
    else:
        source_means = (source_coordinate + shift)[None, :]
        lower = lower_coordinate[:, None]
        upper = upper_coordinate[:, None]
        probabilities = (
            ndtr((upper - source_means) / sigma)
            - ndtr((lower - source_means) / sigma)
        )
        target = probabilities @ raw
    # endif

    target = np.clip(
        np.asarray(target, dtype=float),
        0.0,
        None,
    )
    total = float(np.sum(target))
    if total <= 0.0 or not np.isfinite(total):
        return raw
    # endif
    return target / total


def logit(value: float) -> float:
    value = float(np.clip(value, 1.0e-6, 1.0 - 1.0e-6))
    return math.log(value / (1.0 - value))


def logistic(value: float) -> float:
    if value >= 0.0:
        exponential = math.exp(-value)
        return 1.0 / (1.0 + exponential)
    # endif
    exponential = math.exp(value)
    return exponential / (1.0 + exponential)


def detect_two_component_divider(
    probability: np.ndarray,
    centers: np.ndarray,
    variable_key: str,
) -> Tuple[int, float]:
    """
    Find the valley between the two strongest smoothed template peaks.

    Delta_phi is always divided at pi. theta_gamma_gamma uses peak finding,
    with a weighted-quantile fallback when the second structure is weak.
    """
    probability = normalized_probability(probability)

    if variable_key == "Delta_phi":
        divider_index = int(np.searchsorted(centers, math.pi))
        divider_index = int(
            np.clip(divider_index, 1, centers.size - 2)
        )
        return divider_index, float(centers[divider_index])
    # endif

    smoothed = gaussian_filter1d(
        probability,
        sigma=1.5,
        mode="nearest",
    )
    prominence = max(
        0.015 * float(np.max(smoothed)),
        1.0e-8,
    )
    peaks, properties = find_peaks(
        smoothed,
        prominence=prominence,
        distance=max(3, centers.size // 20),
    )

    if peaks.size >= 2:
        peak_strength = smoothed[peaks]
        strongest = peaks[
            np.argsort(peak_strength)[-2:]
        ]
        left_peak, right_peak = sorted(
            int(index) for index in strongest
        )
        if right_peak - left_peak >= 2:
            valley_relative = int(
                np.argmin(smoothed[left_peak:right_peak + 1])
            )
            divider_index = left_peak + valley_relative
            divider_index = int(
                np.clip(divider_index, 1, centers.size - 2)
            )
            return divider_index, float(centers[divider_index])
        # endif
    # endif

    cumulative = np.cumsum(probability)
    divider_index = int(
        np.searchsorted(cumulative, 0.68)
    )
    divider_index = int(
        np.clip(divider_index, 1, centers.size - 2)
    )
    return divider_index, float(centers[divider_index])


def two_component_morph_probability(
    raw_probability: np.ndarray,
    centers: np.ndarray,
    divider_index: int,
    parameters: np.ndarray,
    use_log_coordinate: bool,
    component_cache: Optional[Mapping[str, object]] = None,
) -> Tuple[np.ndarray, Dict[str, float]]:
    """
    Morph two raw-template components independently.

    parameters:
      left shift, left smear,
      right shift, right smear,
      constrained component-weight logit change.
    """
    raw = normalized_probability(raw_probability)
    if component_cache:
        left_raw = np.asarray(
            component_cache["left_raw"],
            dtype=float,
        )
        right_raw = np.asarray(
            component_cache["right_raw"],
            dtype=float,
        )
        raw_left_weight = float(
            component_cache["raw_left_weight"]
        )
        raw_right_weight = float(
            component_cache["raw_right_weight"]
        )
    else:
        left_raw = raw.copy()
        right_raw = raw.copy()
        left_raw[divider_index + 1:] = 0.0
        right_raw[:divider_index + 1] = 0.0
        raw_left_weight = float(np.sum(left_raw))
        raw_right_weight = float(np.sum(right_raw))
    # endif
    if raw_left_weight <= 0.0 or raw_right_weight <= 0.0:
        fallback = simple_shift_smear_probability(
            raw,
            centers,
            float(parameters[0]),
            float(parameters[1]),
            use_log_coordinate,
        )
        return fallback, {
            "raw_left_weight": raw_left_weight,
            "fitted_left_weight": raw_left_weight,
        }
    # endif

    left_probability = simple_shift_smear_probability(
        left_raw,
        centers,
        float(parameters[0]),
        float(parameters[1]),
        use_log_coordinate,
    )
    right_probability = simple_shift_smear_probability(
        right_raw,
        centers,
        float(parameters[2]),
        float(parameters[3]),
        use_log_coordinate,
    )

    fitted_left_weight = logistic(
        logit(raw_left_weight) + float(parameters[4])
    )
    combined = (
        fitted_left_weight * left_probability
        + (1.0 - fitted_left_weight) * right_probability
    )
    return normalized_probability(combined), {
        "raw_left_weight": raw_left_weight,
        "fitted_left_weight": fitted_left_weight,
    }


def variable_model_spec(
    variable_key: str,
    multiplicity_key: str,
) -> Dict[str, str]:
    """
    Return category-specific template morph models.

    Exact-one deliberately uses only coherent whole-template transformations,
    matching the validated exclusivity-selection fit philosophy.

    Exact-two keeps the previously developed flexible models because the
    diphoton-mass selection is the primary topology discriminator there.
    """
    if multiplicity_key == "one_photon":
        if variable_key == "Delta_phi":
            return {
                "dvcs": "simple_linear",
                "pi0": "simple_linear",
                "description": (
                    "coherent additive shift+smear for complete DVCSGEN "
                    "and AAOGEN templates (exclusivity-style)"
                ),
            }
        # endif
        if variable_key == "theta_gamma_gamma":
            return {
                "dvcs": "simple_log",
                "pi0": "simple_log",
                "description": (
                    "coherent log(theta+epsilon) shift+smear for complete "
                    "DVCSGEN and AAOGEN templates (exclusivity-style)"
                ),
            }
        # endif
        if variable_key == "pTmiss":
            return {
                "dvcs": "simple_log",
                "pi0": "simple_log",
                "description": (
                    "coherent log(pTmiss+epsilon) shift+smear for complete "
                    "DVCSGEN and AAOGEN templates (exclusivity-style)"
                ),
            }
        # endif
        if variable_key == "Emiss2":
            return {
                "dvcs": "simple_linear",
                "pi0": "simple_linear",
                "description": (
                    "coherent additive shift+smear for complete DVCSGEN "
                    "and AAOGEN templates (exclusivity-style)"
                ),
            }
        # endif
    # endif

    # Exact-two categories retain the established flexible modeling.
    if variable_key == "Delta_phi":
        return {
            "dvcs": "simple_linear",
            "pi0": "two_component_linear",
            "description": (
                "DVCSGEN shift+smear; AAOGEN independent left/right "
                "Delta_phi lobes with constrained relative weight"
            ),
        }
    # endif
    if variable_key == "theta_gamma_gamma":
        return {
            "dvcs": "two_component_log",
            "pi0": "two_component_log",
            "description": (
                "Independent low/high reconstruction components in "
                "log(theta+epsilon)"
            ),
        }
    # endif
    if variable_key == "pTmiss":
        return {
            "dvcs": "affine_log_asymmetric",
            "pi0": "affine_log_asymmetric",
            "description": (
                "Log-space affine location/scale plus asymmetric broadening"
            ),
        }
    # endif
    if variable_key == "Emiss2":
        return {
            "dvcs": "affine_linear_asymmetric",
            "pi0": "affine_linear_asymmetric",
            "description": (
                "Linear affine location/scale plus asymmetric broadening "
                "with DVCSGEN/AAOGEN mean ordering"
            ),
        }
    # endif

    raise KeyError(
        f"No fit model configured for {multiplicity_key}/{variable_key}"
    )


def model_parameter_count(model_name: str) -> int:
    if model_name.startswith("simple_"):
        return 2
    # endif
    if model_name.startswith("two_component_"):
        return 5
    # endif
    if model_name.startswith("affine_"):
        return 4
    # endif
    raise KeyError(f"Unknown model name: {model_name}")


def model_initial_parameters(model_name: str) -> np.ndarray:
    return np.zeros(model_parameter_count(model_name), dtype=float)


def model_bounds(
    model_name: str,
    args: argparse.Namespace,
) -> List[Tuple[float, float]]:
    shift = float(args.fit_max_shift_bins)
    smear = float(args.fit_max_smear_bins)
    stretch = float(args.fit_max_log_stretch)

    if model_name.startswith("simple_"):
        return [(-shift, shift), (0.0, smear)]
    # endif
    if model_name.startswith("two_component_"):
        return [
            (-shift, shift),
            (0.0, smear),
            (-shift, shift),
            (0.0, smear),
            (-2.5, 2.5),
        ]
    # endif
    if model_name.startswith("affine_"):
        return [
            (-shift, shift),
            (-stretch, stretch),
            (0.0, smear),
            (0.0, smear),
        ]
    # endif
    raise KeyError(f"Unknown model name: {model_name}")


def nuisance_penalty(
    parameters: np.ndarray,
    model_name: str,
    args: argparse.Namespace,
    context: Optional[Mapping[str, object]] = None,
) -> float:
    """
    Gaussian nuisance penalties.

    For coherent simple models this follows the validated exclusivity fitter:
      * additive: widths are specified in histogram-bin units;
      * positive/log: transformed-space prior widths are 0.20 (shift) and
        0.40 (smearing).

    More flexible exact-two models retain their established penalties.
    """
    if args.disable_fit_nuisance_penalties:
        return 0.0
    # endif

    shift_prior = max(float(args.fit_shift_prior_bins), 1.0e-6)
    smear_prior = max(float(args.fit_smear_prior_bins), 1.0e-6)
    stretch_prior = max(float(args.fit_log_stretch_prior), 1.0e-6)
    weight_prior = max(
        float(args.fit_component_weight_prior),
        1.0e-6,
    )

    if model_name == "simple_linear":
        return float(
            0.5 * (parameters[0] / shift_prior) ** 2
            + 0.5 * (parameters[1] / smear_prior) ** 2
        )
    # endif

    if model_name == "simple_log":
        coordinate_step = (
            float(context["positive_coordinate_step"])
            if context is not None
            else 1.0
        )
        log_shift = float(parameters[0]) * coordinate_step
        log_sigma = float(parameters[1]) * coordinate_step
        return float(
            0.5 * (log_shift / 0.20) ** 2
            + 0.5 * (log_sigma / 0.40) ** 2
        )
    # endif

    if model_name.startswith("two_component_"):
        return float(
            (parameters[0] / shift_prior) ** 2
            + (parameters[1] / smear_prior) ** 2
            + (parameters[2] / shift_prior) ** 2
            + (parameters[3] / smear_prior) ** 2
            + (parameters[4] / weight_prior) ** 2
        )
    # endif

    if model_name.startswith("affine_"):
        return float(
            (parameters[0] / shift_prior) ** 2
            + (parameters[1] / stretch_prior) ** 2
            + (parameters[2] / smear_prior) ** 2
            + (parameters[3] / smear_prior) ** 2
        )
    # endif

    raise KeyError(f"Unknown model name: {model_name}")


def morph_model_probability(
    raw_probability: np.ndarray,
    centers: np.ndarray,
    model_name: str,
    parameters: np.ndarray,
    divider_index: int,
    component_cache: Optional[Mapping[str, object]] = None,
) -> Tuple[np.ndarray, Dict[str, float]]:
    if model_name == "simple_linear":
        return (
            simple_shift_smear_probability(
                raw_probability,
                centers,
                parameters[0],
                parameters[1],
                False,
            ),
            {},
        )
    # endif
    if model_name == "simple_log":
        return (
            simple_shift_smear_probability(
                raw_probability,
                centers,
                parameters[0],
                parameters[1],
                True,
            ),
            {},
        )
    # endif
    if model_name == "two_component_linear":
        return two_component_morph_probability(
            raw_probability,
            centers,
            divider_index,
            parameters,
            False,
        )
    # endif
    if model_name == "two_component_log":
        return two_component_morph_probability(
            raw_probability,
            centers,
            divider_index,
            parameters,
            True,
        )
    # endif
    if model_name == "affine_log_asymmetric":
        return (
            affine_morph_probability(
                raw_probability,
                centers,
                parameters[0],
                parameters[1],
                parameters[2],
                parameters[3],
                True,
            ),
            {},
        )
    # endif
    if model_name == "affine_linear_asymmetric":
        return (
            affine_morph_probability(
                raw_probability,
                centers,
                parameters[0],
                parameters[1],
                parameters[2],
                parameters[3],
                False,
            ),
            {},
        )
    # endif
    raise KeyError(f"Unknown model name: {model_name}")


def nuisance_dictionary(
    model_name: str,
    parameters: np.ndarray,
    metadata: Mapping[str, float],
) -> Dict[str, float]:
    if model_name.startswith("simple_"):
        output = {
            "shift_bins": float(parameters[0]),
            "smear_bins": float(parameters[1]),
        }
    elif model_name.startswith("two_component_"):
        output = {
            "left_shift_bins": float(parameters[0]),
            "left_smear_bins": float(parameters[1]),
            "right_shift_bins": float(parameters[2]),
            "right_smear_bins": float(parameters[3]),
            "component_weight_logit_delta": float(parameters[4]),
        }
    elif model_name.startswith("affine_"):
        output = {
            "shift_bins": float(parameters[0]),
            "log_scale": float(parameters[1]),
            "scale": float(math.exp(parameters[1])),
            "lower_smear_bins": float(parameters[2]),
            "upper_smear_bins": float(parameters[3]),
        }
    else:
        raise KeyError(f"Unknown model name: {model_name}")
    # endif

    output.update(
        {key: float(value) for key, value in metadata.items()}
    )
    return output


def prepare_variable_fit_context(
    variable: FitVariable,
    data_counts: np.ndarray,
    dvcs_counts: np.ndarray,
    pi0_counts: np.ndarray,
    multiplicity_key: str,
) -> Dict[str, object]:
    centers = np.linspace(
        variable.low,
        variable.high,
        variable.bins,
        endpoint=False,
        dtype=float,
    )
    bin_width = (variable.high - variable.low) / variable.bins
    centers += 0.5 * bin_width

    positive_epsilon = max(
        0.25 * float(bin_width),
        1.0e-8,
    )
    positive_coordinate = np.log(
        np.maximum(centers + positive_epsilon, 1.0e-15)
    )
    positive_coordinate_step = (
        float(np.median(np.diff(positive_coordinate)))
        if positive_coordinate.size > 1
        else 1.0
    )

    raw_dvcs = normalized_probability(dvcs_counts)
    raw_pi0 = normalized_probability(pi0_counts)
    model_spec = variable_model_spec(variable.key, multiplicity_key)

    dvcs_divider_index, dvcs_divider_value = (
        detect_two_component_divider(
            raw_dvcs,
            centers,
            variable.key,
        )
        if model_spec["dvcs"].startswith("two_component_")
        else (max(1, centers.size // 2), math.nan)
    )
    pi0_divider_index, pi0_divider_value = (
        detect_two_component_divider(
            raw_pi0,
            centers,
            variable.key,
        )
        if model_spec["pi0"].startswith("two_component_")
        else (max(1, centers.size // 2), math.nan)
    )

    def component_cache(
        raw_probability: np.ndarray,
        model_name: str,
        divider_index: int,
    ) -> Dict[str, object]:
        if not model_name.startswith("two_component_"):
            return {}
        # endif
        left_raw = raw_probability.copy()
        right_raw = raw_probability.copy()
        left_raw[divider_index + 1:] = 0.0
        right_raw[:divider_index + 1] = 0.0
        return {
            "left_raw": left_raw,
            "right_raw": right_raw,
            "raw_left_weight": float(np.sum(left_raw)),
            "raw_right_weight": float(np.sum(right_raw)),
        }

    return {
        "variable": variable,
        "multiplicity_key": multiplicity_key,
        "centers": centers,
        "bin_width": float(bin_width),
        "positive_coordinate_step": positive_coordinate_step,
        "data": np.asarray(data_counts, dtype=float),
        "data_total": float(np.sum(data_counts)),
        "raw_dvcs": raw_dvcs,
        "raw_pi0": raw_pi0,
        "dvcs_model": model_spec["dvcs"],
        "pi0_model": model_spec["pi0"],
        "model_description": model_spec["description"],
        "dvcs_divider_index": dvcs_divider_index,
        "dvcs_divider_value": dvcs_divider_value,
        "pi0_divider_index": pi0_divider_index,
        "pi0_divider_value": pi0_divider_value,
        "dvcs_component_cache": component_cache(
            raw_dvcs,
            model_spec["dvcs"],
            dvcs_divider_index,
        ),
        "pi0_component_cache": component_cache(
            raw_pi0,
            model_spec["pi0"],
            pi0_divider_index,
        ),
    }


def evaluate_variable_model(
    context: Mapping[str, object],
    f_pi0: float,
    dvcs_parameters: np.ndarray,
    pi0_parameters: np.ndarray,
    args: argparse.Namespace,
) -> Dict[str, object]:
    centers = context["centers"]
    dvcs_probability, dvcs_metadata = morph_model_probability(
        context["raw_dvcs"],
        centers,
        context["dvcs_model"],
        dvcs_parameters,
        context["dvcs_divider_index"],
        component_cache=context["dvcs_component_cache"],
    )
    pi0_probability, pi0_metadata = morph_model_probability(
        context["raw_pi0"],
        centers,
        context["pi0_model"],
        pi0_parameters,
        context["pi0_divider_index"],
        component_cache=context["pi0_component_cache"],
    )

    data_total = float(context["data_total"])
    f_pi0 = float(np.clip(f_pi0, 1.0e-6, 1.0 - 1.0e-6))
    dvcs_component = data_total * (1.0 - f_pi0) * dvcs_probability
    pi0_component = data_total * f_pi0 * pi0_probability
    model = dvcs_component + pi0_component

    deviance = poisson_deviance(
        context["data"],
        model,
    )
    penalty = (
        nuisance_penalty(
            dvcs_parameters,
            str(context["dvcs_model"]),
            args,
            context,
        )
        + nuisance_penalty(
            pi0_parameters,
            str(context["pi0_model"]),
            args,
            context,
        )
    )

    variable_key = str(context["variable"].key)
    if variable_key == "Emiss2":
        dvcs_mean = weighted_mean(dvcs_probability, centers)
        pi0_mean = weighted_mean(pi0_probability, centers)
        if (
            np.isfinite(dvcs_mean)
            and np.isfinite(pi0_mean)
            and dvcs_mean >= pi0_mean
        ):
            violation_bins = (
                (dvcs_mean - pi0_mean)
                / max(float(context["bin_width"]), 1.0e-12)
            )
            penalty += (
                float(args.fit_mean_order_penalty)
                * (1.0 + violation_bins ** 2)
            )
        # endif
    # endif

    return {
        "objective": float(deviance + penalty),
        "deviance": float(deviance),
        "penalty": float(penalty),
        "dvcs_probability": dvcs_probability,
        "pi0_probability": pi0_probability,
        "dvcs_component": dvcs_component,
        "pi0_component": pi0_component,
        "model": model,
        "dvcs_metadata": dvcs_metadata,
        "pi0_metadata": pi0_metadata,
    }


def optimize_variable_independent_fraction(
    context: Mapping[str, object],
    args: argparse.Namespace,
) -> Dict[str, object]:
    dvcs_model = str(context["dvcs_model"])
    pi0_model = str(context["pi0_model"])
    ndvcs = model_parameter_count(dvcs_model)
    npi0 = model_parameter_count(pi0_model)

    bounds = (
        [(1.0e-4, 1.0 - 1.0e-4)]
        + model_bounds(dvcs_model, args)
        + model_bounds(pi0_model, args)
    )

    requested_starts = max(
        1,
        int(args.fit_independent_multistarts),
    )
    candidate_fractions = np.linspace(
        0.15,
        0.95,
        requested_starts,
    )
    best_result = None

    for initial_fraction in candidate_fractions:
        initial = np.concatenate(
            [
                np.asarray([initial_fraction], dtype=float),
                model_initial_parameters(dvcs_model),
                model_initial_parameters(pi0_model),
            ]
        )

        def objective(parameters: np.ndarray) -> float:
            return float(
                evaluate_variable_model(
                    context,
                    float(parameters[0]),
                    parameters[1:1 + ndvcs],
                    parameters[1 + ndvcs:1 + ndvcs + npi0],
                    args,
                )["objective"]
            )

        result = minimize(
            objective,
            initial,
            method="L-BFGS-B",
            bounds=bounds,
            options={
                "maxiter": 250 if args.fast_fits else 700,
                "ftol": 1.0e-8 if args.fast_fits else 1.0e-10,
                "gtol": 1.0e-5 if args.fast_fits else 1.0e-6,
                "maxls": 25 if args.fast_fits else 40,
            },
        )
        if (
            best_result is None
            or float(result.fun) < float(best_result.fun)
        ):
            best_result = result
        # endif
    # endfor

    if best_result is None:
        raise RuntimeError("Independent variable fit produced no result.")
    # endif

    parameters = np.asarray(best_result.x, dtype=float)
    evaluation = evaluate_variable_model(
        context,
        float(parameters[0]),
        parameters[1:1 + ndvcs],
        parameters[1 + ndvcs:1 + ndvcs + npi0],
        args,
    )

    fraction_error = curvature_fraction_error(
        context,
        float(parameters[0]),
        parameters[1:1 + ndvcs],
        parameters[1 + ndvcs:1 + ndvcs + npi0],
        args,
    )

    return {
        "success": bool(best_result.success),
        "message": str(best_result.message),
        "f_pi0": float(parameters[0]),
        "f_pi0_err": fraction_error,
        "dvcs_parameters": parameters[1:1 + ndvcs].copy(),
        "pi0_parameters": parameters[
            1 + ndvcs:1 + ndvcs + npi0
        ].copy(),
        "evaluation": evaluation,
        "optimizer_objective": float(best_result.fun),
    }


def optimize_nuisance_fixed_fraction(
    context: Mapping[str, object],
    f_pi0: float,
    initial_dvcs: np.ndarray,
    initial_pi0: np.ndarray,
    args: argparse.Namespace,
) -> Dict[str, object]:
    dvcs_model = str(context["dvcs_model"])
    pi0_model = str(context["pi0_model"])
    ndvcs = model_parameter_count(dvcs_model)

    initial = np.concatenate(
        [
            np.asarray(initial_dvcs, dtype=float),
            np.asarray(initial_pi0, dtype=float),
        ]
    )
    bounds = (
        model_bounds(dvcs_model, args)
        + model_bounds(pi0_model, args)
    )

    def objective(parameters: np.ndarray) -> float:
        return float(
            evaluate_variable_model(
                context,
                f_pi0,
                parameters[:ndvcs],
                parameters[ndvcs:],
                args,
            )["objective"]
        )

    result = minimize(
        objective,
        initial,
        method="L-BFGS-B",
        bounds=bounds,
        options={
            "maxiter": 500,
            "ftol": 1.0e-10,
            "gtol": 1.0e-6,
            "maxls": 40,
        },
    )
    parameters = np.asarray(result.x, dtype=float)
    evaluation = evaluate_variable_model(
        context,
        f_pi0,
        parameters[:ndvcs],
        parameters[ndvcs:],
        args,
    )

    return {
        "success": bool(result.success),
        "message": str(result.message),
        "dvcs_parameters": parameters[:ndvcs].copy(),
        "pi0_parameters": parameters[ndvcs:].copy(),
        "evaluation": evaluation,
    }


def curvature_fraction_error(
    context: Mapping[str, object],
    f_pi0: float,
    dvcs_parameters: np.ndarray,
    pi0_parameters: np.ndarray,
    args: argparse.Namespace,
) -> float:
    step = max(2.0e-4, min(0.01, 0.05 * min(f_pi0, 1.0 - f_pi0)))
    low = max(1.0e-5, f_pi0 - step)
    high = min(1.0 - 1.0e-5, f_pi0 + step)
    if high <= f_pi0 or low >= f_pi0:
        return math.nan
    # endif

    center_value = float(
        evaluate_variable_model(
            context,
            f_pi0,
            dvcs_parameters,
            pi0_parameters,
            args,
        )["objective"]
    )
    low_value = float(
        evaluate_variable_model(
            context,
            low,
            dvcs_parameters,
            pi0_parameters,
            args,
        )["objective"]
    )
    high_value = float(
        evaluate_variable_model(
            context,
            high,
            dvcs_parameters,
            pi0_parameters,
            args,
        )["objective"]
    )

    effective_step = 0.5 * (high - low)
    second_derivative = (
        high_value - 2.0 * center_value + low_value
    ) / max(effective_step ** 2, 1.0e-12)
    if second_derivative <= 0.0:
        return math.nan
    # endif

    return float(math.sqrt(2.0 / second_derivative))


def shared_fraction_objective(
    f_pi0: float,
    contexts: Mapping[str, Mapping[str, object]],
    states: Mapping[str, Mapping[str, object]],
    args: argparse.Namespace,
) -> float:
    total = 0.0
    for variable_key, context in contexts.items():
        state = states[variable_key]
        total += float(
            evaluate_variable_model(
                context,
                f_pi0,
                np.asarray(state["dvcs_parameters"], dtype=float),
                np.asarray(state["pi0_parameters"], dtype=float),
                args,
            )["objective"]
        )
    # endfor
    return float(total)


def shared_fraction_curvature_error(
    f_pi0: float,
    contexts: Mapping[str, Mapping[str, object]],
    states: Mapping[str, Mapping[str, object]],
    args: argparse.Namespace,
) -> float:
    step = max(2.0e-4, min(0.006, 0.04 * min(f_pi0, 1.0 - f_pi0)))
    low = max(1.0e-5, f_pi0 - step)
    high = min(1.0 - 1.0e-5, f_pi0 + step)
    center = shared_fraction_objective(
        f_pi0,
        contexts,
        states,
        args,
    )
    low_value = shared_fraction_objective(
        low,
        contexts,
        states,
        args,
    )
    high_value = shared_fraction_objective(
        high,
        contexts,
        states,
        args,
    )
    effective_step = 0.5 * (high - low)
    second_derivative = (
        high_value - 2.0 * center + low_value
    ) / max(effective_step ** 2, 1.0e-12)
    if second_derivative <= 0.0:
        return math.nan
    # endif
    return float(math.sqrt(2.0 / second_derivative))


def fit_histogram_inputs(
    sample_payloads: Mapping[str, Mapping[str, object]],
    multiplicity_key: str,
) -> Tuple[
    Dict[str, np.ndarray],
    Dict[str, np.ndarray],
    Dict[str, np.ndarray],
]:
    fit_variables = fit_variables_for_category(multiplicity_key)

    def extract(sample_key: str) -> Dict[str, np.ndarray]:
        return {
            variable.key: np.asarray(
                sample_payloads[sample_key][multiplicity_key][
                    "after_histograms"
                ][variable.key]["counts"],
                dtype=np.float64,
            )
            for variable in fit_variables
        }

    return extract("data"), extract("dvcsgen"), extract("aaogen")



def raw_template_model_counts(
    context: Mapping[str, object],
    f_pi0: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Fraction-only model using the untouched DVCSGEN and AAOGEN templates.
    """
    f_pi0 = float(np.clip(f_pi0, 1.0e-8, 1.0 - 1.0e-8))
    data_total = float(context["data_total"])
    raw_dvcs = np.asarray(context["raw_dvcs"], dtype=float)
    raw_pi0 = np.asarray(context["raw_pi0"], dtype=float)

    dvcs_component = data_total * (1.0 - f_pi0) * raw_dvcs
    pi0_component = data_total * f_pi0 * raw_pi0
    model = dvcs_component + pi0_component
    return model, dvcs_component, pi0_component


def raw_template_objective(
    context: Mapping[str, object],
    f_pi0: float,
) -> float:
    model, _dvcs, _pi0 = raw_template_model_counts(
        context,
        f_pi0,
    )
    return poisson_deviance(
        np.asarray(context["data"], dtype=float),
        model,
    )


def scalar_fraction_curvature_error(
    objective,
    f_pi0: float,
) -> float:
    """Curvature error for a one-dimensional fraction objective."""
    step = max(
        2.0e-4,
        min(0.006, 0.04 * min(f_pi0, 1.0 - f_pi0)),
    )
    low = max(1.0e-5, f_pi0 - step)
    high = min(1.0 - 1.0e-5, f_pi0 + step)

    center_value = float(objective(f_pi0))
    low_value = float(objective(low))
    high_value = float(objective(high))
    effective_step = 0.5 * (high - low)

    second_derivative = (
        high_value - 2.0 * center_value + low_value
    ) / max(effective_step ** 2, 1.0e-12)

    if second_derivative <= 0.0 or not np.isfinite(second_derivative):
        return math.nan
    # endif

    return float(math.sqrt(2.0 / second_derivative))


def fit_raw_template_fraction(
    context: Mapping[str, object],
) -> Dict[str, float]:
    """Independent fraction-only raw-template fit for one projection."""
    result = minimize_scalar(
        lambda fraction: raw_template_objective(
            context,
            float(fraction),
        ),
        bounds=(1.0e-5, 1.0 - 1.0e-5),
        method="bounded",
        options={
            "xatol": 2.0e-7,
            "maxiter": 300,
        },
    )
    f_pi0 = float(result.x)
    error = scalar_fraction_curvature_error(
        lambda fraction: raw_template_objective(
            context,
            float(fraction),
        ),
        f_pi0,
    )
    return {
        "success": bool(result.success),
        "f_pi0": f_pi0,
        "f_pi0_err": error,
        "deviance": float(result.fun),
    }


def fit_shared_raw_template_fraction(
    contexts: Mapping[str, Mapping[str, object]],
) -> Dict[str, float]:
    """Shared fraction-only raw-template fit across all projections."""
    def objective(fraction: float) -> float:
        return float(
            sum(
                raw_template_objective(
                    context,
                    float(fraction),
                )
                for context in contexts.values()
            )
        )

    result = minimize_scalar(
        objective,
        bounds=(1.0e-5, 1.0 - 1.0e-5),
        method="bounded",
        options={
            "xatol": 2.0e-7,
            "maxiter": 400,
        },
    )
    f_pi0 = float(result.x)
    error = scalar_fraction_curvature_error(
        objective,
        f_pi0,
    )
    return {
        "success": bool(result.success),
        "f_pi0": f_pi0,
        "f_pi0_err": error,
        "deviance": float(result.fun),
    }

def run_one_shared_template_fit(
    period: PeriodConfig,
    multiplicity_key: str,
    multiplicity_label: str,
    sample_payloads: Mapping[str, Mapping[str, object]],
    args: argparse.Namespace,
) -> Tuple[AdvancedSharedFit, Dict[str, np.ndarray]]:
    data_histograms, dvcs_histograms, pi0_histograms = (
        fit_histogram_inputs(
            sample_payloads,
            multiplicity_key,
        )
    )

    contexts: Dict[str, Dict[str, object]] = {}
    independent_states: Dict[str, Dict[str, object]] = {}
    raw_independent_states: Dict[str, Dict[str, float]] = {}
    fit_variables = fit_variables_for_category(multiplicity_key)

    for variable in fit_variables:
        data_count = int(np.sum(data_histograms[variable.key]))
        dvcs_count = int(np.sum(dvcs_histograms[variable.key]))
        pi0_count = int(np.sum(pi0_histograms[variable.key]))
        if (
            data_count < args.fit_min_counts
            or dvcs_count <= 0
            or pi0_count <= 0
        ):
            raise RuntimeError(
                f"{period.label} {multiplicity_label} {variable.key}: "
                f"insufficient fit support: data={data_count}, "
                f"DVCSGEN={dvcs_count}, AAOGEN={pi0_count}"
            )
        # endif

        context = prepare_variable_fit_context(
            variable,
            data_histograms[variable.key],
            dvcs_histograms[variable.key],
            pi0_histograms[variable.key],
            multiplicity_key,
        )
        contexts[variable.key] = context

        raw_independent_states[variable.key] = (
            fit_raw_template_fraction(context)
        )
        independent_states[variable.key] = (
            optimize_variable_independent_fraction(
                context,
                args,
            )
        )

        raw_state = raw_independent_states[variable.key]
        morphed_state = independent_states[variable.key]
        log(
            f"{period.label} {multiplicity_label} {variable.key}: "
            f"raw f_pi0={raw_state['f_pi0']:.5f} +/- "
            f"{raw_state['f_pi0_err']:.5f}; "
            f"morphed f_pi0={morphed_state['f_pi0']:.5f} +/- "
            f"{morphed_state['f_pi0_err']:.5f}; "
            f"Delta={morphed_state['f_pi0'] - raw_state['f_pi0']:+.5f}"
        )
    # endfor

    raw_shared_state = fit_shared_raw_template_fraction(
        contexts
    )
    raw_shared_fraction = float(raw_shared_state["f_pi0"])
    raw_shared_error = float(raw_shared_state["f_pi0_err"])

    weights = []
    fractions = []
    for variable_key, state in independent_states.items():
        error = float(state["f_pi0_err"])
        weight = (
            1.0 / error ** 2
            if np.isfinite(error) and error > 0.0
            else 1.0
        )
        weights.append(weight)
        fractions.append(float(state["f_pi0"]))
    # endfor
    shared_fraction = float(
        np.average(fractions, weights=weights)
    )
    shared_fraction = float(
        np.clip(shared_fraction, 1.0e-4, 1.0 - 1.0e-4)
    )

    states: Dict[str, Dict[str, object]] = {
        variable_key: {
            "dvcs_parameters": np.asarray(
                state["dvcs_parameters"],
                dtype=float,
            ).copy(),
            "pi0_parameters": np.asarray(
                state["pi0_parameters"],
                dtype=float,
            ).copy(),
        }
        for variable_key, state in independent_states.items()
    }

    completed_iterations = 0
    coordinate_iterations = max(
        1,
        int(args.fit_coordinate_iterations),
    )
    if args.fast_fits:
        coordinate_iterations = min(coordinate_iterations, 2)
    # endif

    for iteration in range(coordinate_iterations):
        previous_fraction = shared_fraction

        for variable_key, context in contexts.items():
            optimized = optimize_nuisance_fixed_fraction(
                context,
                shared_fraction,
                np.asarray(
                    states[variable_key]["dvcs_parameters"],
                    dtype=float,
                ),
                np.asarray(
                    states[variable_key]["pi0_parameters"],
                    dtype=float,
                ),
                args,
            )
            states[variable_key] = {
                "dvcs_parameters": optimized["dvcs_parameters"],
                "pi0_parameters": optimized["pi0_parameters"],
            }
        # endfor

        fraction_result = minimize_scalar(
            lambda fraction: shared_fraction_objective(
                float(fraction),
                contexts,
                states,
                args,
            ),
            bounds=(1.0e-4, 1.0 - 1.0e-4),
            method="bounded",
            options={
                "xatol": 2.0e-5 if args.fast_fits else 2.0e-6,
                "maxiter": 250,
            },
        )
        shared_fraction = float(fraction_result.x)
        completed_iterations = iteration + 1

        if abs(shared_fraction - previous_fraction) < 2.0e-5:
            break
        # endif
    # endfor

    shared_error = shared_fraction_curvature_error(
        shared_fraction,
        contexts,
        states,
        args,
    )

    variable_results: Dict[str, AdvancedVariableFit] = {}
    total_deviance = 0.0
    total_bins = 0
    total_nuisance = 1

    for variable in fit_variables:
        context = contexts[variable.key]
        state = states[variable.key]
        evaluation = evaluate_variable_model(
            context,
            shared_fraction,
            np.asarray(state["dvcs_parameters"], dtype=float),
            np.asarray(state["pi0_parameters"], dtype=float),
            args,
        )
        independent = independent_states[variable.key]
        data_total = float(context["data_total"])
        raw_dvcs = np.asarray(context["raw_dvcs"], dtype=float)
        raw_pi0 = np.asarray(context["raw_pi0"], dtype=float)

        (
            raw_shared_model,
            raw_shared_dvcs_component,
            raw_shared_pi0_component,
        ) = raw_template_model_counts(
            context,
            raw_shared_fraction,
        )
        raw_shared_variable_deviance = poisson_deviance(
            np.asarray(context["data"], dtype=float),
            raw_shared_model,
        )
        raw_independent = raw_independent_states[variable.key]

        ndvcs = model_parameter_count(str(context["dvcs_model"]))
        npi0 = model_parameter_count(str(context["pi0_model"]))
        nuisance_count = ndvcs + npi0
        ndf = max(
            1,
            int(np.count_nonzero(np.asarray(context["data"]) >= 0.0))
            - nuisance_count
            - 1,
        )

        result = AdvancedVariableFit(
            success=True,
            message="advanced component-preserving morph fit",
            variable=variable.key,
            model_name=str(context["model_description"]),
            f_pi0=shared_fraction,
            f_pi0_err=shared_error,
            independent_f_pi0=float(independent["f_pi0"]),
            independent_f_pi0_err=float(independent["f_pi0_err"]),
            deviance=float(evaluation["deviance"]),
            ndf=ndf,
            nuisance_count=nuisance_count,
            model_counts=np.asarray(evaluation["model"], dtype=float),
            dvcs_component_counts=np.asarray(
                evaluation["dvcs_component"],
                dtype=float,
            ),
            pi0_component_counts=np.asarray(
                evaluation["pi0_component"],
                dtype=float,
            ),
            raw_dvcs_component_counts=(
                data_total * (1.0 - shared_fraction) * raw_dvcs
            ),
            raw_pi0_component_counts=(
                data_total * shared_fraction * raw_pi0
            ),
            morphed_dvcs_probability=np.asarray(
                evaluation["dvcs_probability"],
                dtype=float,
            ),
            morphed_pi0_probability=np.asarray(
                evaluation["pi0_probability"],
                dtype=float,
            ),
            raw_dvcs_probability=raw_dvcs,
            raw_pi0_probability=raw_pi0,
            dvcs_nuisance=nuisance_dictionary(
                str(context["dvcs_model"]),
                np.asarray(state["dvcs_parameters"], dtype=float),
                evaluation["dvcs_metadata"],
            ),
            pi0_nuisance=nuisance_dictionary(
                str(context["pi0_model"]),
                np.asarray(state["pi0_parameters"], dtype=float),
                evaluation["pi0_metadata"],
            ),
            divider_dvcs=float(context["dvcs_divider_value"]),
            divider_pi0=float(context["pi0_divider_value"]),
            fit_total=float(np.sum(evaluation["model"])),
            data_total=data_total,
            raw_shared_f_pi0=raw_shared_fraction,
            raw_shared_f_pi0_err=raw_shared_error,
            raw_independent_f_pi0=float(
                raw_independent["f_pi0"]
            ),
            raw_independent_f_pi0_err=float(
                raw_independent["f_pi0_err"]
            ),
            raw_shared_deviance=float(
                raw_shared_variable_deviance
            ),
            raw_independent_deviance=float(
                raw_independent["deviance"]
            ),
            raw_shared_model_counts=np.asarray(
                raw_shared_model,
                dtype=float,
            ),
            raw_shared_dvcs_component_counts=np.asarray(
                raw_shared_dvcs_component,
                dtype=float,
            ),
            raw_shared_pi0_component_counts=np.asarray(
                raw_shared_pi0_component,
                dtype=float,
            ),
        )
        variable_results[variable.key] = result
        total_deviance += result.deviance
        total_bins += variable.bins
        total_nuisance += nuisance_count
    # endfor

    global_ndf = max(1, total_bins - total_nuisance)
    summary = AdvancedSharedFit(
        success=True,
        message=(
            "coordinate-minimized shared fraction across four projections with "
            "variable-specific component-preserving morphs"
        ),
        f_pi0=shared_fraction,
        f_pi0_err=shared_error,
        deviance=float(total_deviance),
        ndf=int(global_ndf),
        variable_results=variable_results,
        coordinate_iterations=completed_iterations,
        multiplicity_key=multiplicity_key,
        raw_f_pi0=raw_shared_fraction,
        raw_f_pi0_err=raw_shared_error,
        raw_deviance=float(raw_shared_state["deviance"]),
        raw_ndf=max(
            1,
            int(sum(variable.bins for variable in fit_variables)) - 1,
        ),
    )
    return summary, data_histograms


def variable_fit_to_dict(
    result: AdvancedVariableFit,
) -> Dict[str, object]:
    return {
        "success": result.success,
        "message": result.message,
        "variable": result.variable,
        "model_name": result.model_name,
        "shared_f_pi0": result.f_pi0,
        "shared_f_pi0_error": result.f_pi0_err,
        "independent_f_pi0": result.independent_f_pi0,
        "independent_f_pi0_error": result.independent_f_pi0_err,
        "deviance": result.deviance,
        "ndf": result.ndf,
        "deviance_per_ndf": (
            result.deviance / result.ndf
            if result.ndf > 0
            else math.nan
        ),
        "raw_shared_f_pi0": result.raw_shared_f_pi0,
        "raw_shared_f_pi0_error": result.raw_shared_f_pi0_err,
        "raw_independent_f_pi0": result.raw_independent_f_pi0,
        "raw_independent_f_pi0_error": result.raw_independent_f_pi0_err,
        "delta_f_morph_shared": (
            result.f_pi0 - result.raw_shared_f_pi0
        ),
        "delta_f_morph_independent": (
            result.independent_f_pi0
            - result.raw_independent_f_pi0
        ),
        "raw_shared_deviance": result.raw_shared_deviance,
        "raw_independent_deviance": result.raw_independent_deviance,
        "nuisance_count": result.nuisance_count,
        "dvcs_nuisance": result.dvcs_nuisance,
        "pi0_nuisance": result.pi0_nuisance,
        "divider_dvcs": result.divider_dvcs,
        "divider_pi0": result.divider_pi0,
        "fit_total": result.fit_total,
        "data_total": result.data_total,
    }


def shared_fit_to_dict(
    summary: AdvancedSharedFit,
) -> Dict[str, object]:
    return {
        "success": summary.success,
        "message": summary.message,
        "f_pi0": summary.f_pi0,
        "f_pi0_error": summary.f_pi0_err,
        "deviance": summary.deviance,
        "ndf": summary.ndf,
        "deviance_per_ndf": (
            summary.deviance / summary.ndf
            if summary.ndf > 0
            else math.nan
        ),
        "coordinate_iterations": summary.coordinate_iterations,
        "raw_f_pi0": summary.raw_f_pi0,
        "raw_f_pi0_error": summary.raw_f_pi0_err,
        "raw_deviance": summary.raw_deviance,
        "raw_ndf": summary.raw_ndf,
        "raw_deviance_per_ndf": (
            summary.raw_deviance / summary.raw_ndf
            if summary.raw_ndf > 0
            else math.nan
        ),
        "delta_f_morph": (
            summary.f_pi0 - summary.raw_f_pi0
        ),
        "variables": {
            variable_key: variable_fit_to_dict(result)
            for variable_key, result
            in summary.variable_results.items()
        },
    }


def plot_shared_template_fit(
    path: Path,
    period_label: str,
    multiplicity_label: str,
    summary: AdvancedSharedFit,
    data_histograms: Mapping[str, np.ndarray],
    args: argparse.Namespace,
) -> None:
    """
    Compare the fraction-only raw-template fit with the constrained morphed fit.

    For exact-one events the raw-template result is the primary diagnostic and
    the morphed fit is a model-sensitivity comparison.
    """
    path.parent.mkdir(parents=True, exist_ok=True)

    sample_colors = {
        "data": "black",
        "dvcs_mc": "tab:blue",
        "pi0_mc": "tab:red",
        "morphed_fit": "tab:green",
        "raw_fit": "0.35",
    }

    fit_variables = fit_variables_for_category(
        summary.multiplicity_key
    )
    number_of_variables = len(fit_variables)
    fig, axes = plt.subplots(
        2,
        number_of_variables,
        figsize=(4.75 * number_of_variables, 9.5),
        sharex="col",
        squeeze=False,
        gridspec_kw={"height_ratios": [3.0, 1.15]},
    )

    for column, variable in enumerate(fit_variables):
        result = summary.variable_results[variable.key]
        top = axes[0, column]
        bottom = axes[1, column]

        data_counts = np.asarray(
            data_histograms[variable.key],
            dtype=float,
        )
        data_error = np.sqrt(np.maximum(data_counts, 0.0))
        edges = np.linspace(
            variable.low,
            variable.high,
            variable.bins + 1,
        )
        centers = 0.5 * (edges[:-1] + edges[1:])

        top.errorbar(
            centers,
            data_counts,
            yerr=data_error,
            fmt="o",
            markersize=2.4,
            linewidth=0.8,
            capsize=0.0,
            color=sample_colors["data"],
            label=r"$e'p'\gamma$ data",
            zorder=6,
        )

        # Raw-template components are plotted with thin dashed outlines and are
        # normalized using the shared RAW fraction, not the morphed fraction.
        top.stairs(
            result.raw_shared_dvcs_component_counts,
            edges,
            color=sample_colors["dvcs_mc"],
            linewidth=0.75,
            linestyle="--",
            label="raw DVCS MC component",
            zorder=1,
        )
        top.stairs(
            result.raw_shared_pi0_component_counts,
            edges,
            color=sample_colors["pi0_mc"],
            linewidth=0.75,
            linestyle="--",
            label=r"raw $e\pi^0$ MC component",
            zorder=1,
        )
        top.stairs(
            result.raw_shared_model_counts,
            edges,
            color=sample_colors["raw_fit"],
            linewidth=1.8,
            linestyle="-.",
            label="raw-template fraction-only fit",
            zorder=4,
        )

        # Morphed comparison.
        top.stairs(
            result.dvcs_component_counts,
            edges,
            color=sample_colors["dvcs_mc"],
            linewidth=1.45,
            label="morphed DVCS component",
            zorder=3,
        )
        top.stairs(
            result.pi0_component_counts,
            edges,
            color=sample_colors["pi0_mc"],
            linewidth=1.45,
            label=r"morphed $e\pi^0$ component",
            zorder=2,
        )
        top.stairs(
            result.model_counts,
            edges,
            color=sample_colors["morphed_fit"],
            linewidth=1.9,
            linestyle="--",
            label="constrained morphed fit",
            zorder=5,
        )

        quality = (
            rf"raw shared: $f_{{\pi^0}}={summary.raw_f_pi0:.3f}$"
            rf"$\pm{summary.raw_f_pi0_err:.3f}$"
            + "\n"
            + rf"raw independent: {result.raw_independent_f_pi0:.3f}"
            rf"$\pm{result.raw_independent_f_pi0_err:.3f}$"
            + "\n"
            + rf"morph shared: {summary.f_pi0:.3f}"
            rf"$\pm{summary.f_pi0_err:.3f}$"
            + "\n"
            + rf"morph independent: {result.independent_f_pi0:.3f}"
            rf"$\pm{result.independent_f_pi0_err:.3f}$"
            + "\n"
            + rf"$\Delta f_{{morph}}={summary.f_pi0-summary.raw_f_pi0:+.3f}$"
            + "\n"
            + rf"$D/ndf$ raw={result.raw_shared_deviance:.1f}/"
            f"{max(1, variable.bins - 1)}; "
            + f"morph={result.deviance:.1f}/{result.ndf}"
        )
        top.text(
            0.98,
            0.96,
            quality,
            transform=top.transAxes,
            ha="right",
            va="top",
            fontsize=7.7,
            bbox={
                "facecolor": "white",
                "alpha": 0.80,
                "edgecolor": "none",
            },
            zorder=10,
        )

        top.set_xlim(variable.low, variable.high)
        top.set_ylabel("events / bin")
        top.grid(axis="y", alpha=0.25)

        if args.shape_log_y:
            positive_arrays = (
                data_counts,
                result.raw_shared_model_counts,
                result.model_counts,
            )
            positive_parts = [
                np.asarray(array, dtype=float)[
                    np.asarray(array, dtype=float) > 0.0
                ]
                for array in positive_arrays
                if np.any(np.asarray(array, dtype=float) > 0.0)
            ]
            if positive_parts:
                positive_values = np.concatenate(positive_parts)
                top.set_yscale("log")
                top.set_ylim(
                    bottom=max(
                        0.5,
                        0.5 * float(np.min(positive_values)),
                    )
                )
            # endif
        else:
            top.set_ylim(bottom=0.0)
        # endif

        raw_model = np.asarray(
            result.raw_shared_model_counts,
            dtype=float,
        )
        morph_model = np.asarray(
            result.model_counts,
            dtype=float,
        )

        raw_valid = raw_model > 0.0
        morph_valid = morph_model > 0.0

        raw_ratio = np.full_like(data_counts, np.nan, dtype=float)
        raw_ratio_error = np.full_like(data_counts, np.nan, dtype=float)
        raw_ratio[raw_valid] = (
            data_counts[raw_valid] / raw_model[raw_valid]
        )
        raw_ratio_error[raw_valid] = (
            data_error[raw_valid] / raw_model[raw_valid]
        )

        morph_ratio = np.full_like(data_counts, np.nan, dtype=float)
        morph_ratio_error = np.full_like(data_counts, np.nan, dtype=float)
        morph_ratio[morph_valid] = (
            data_counts[morph_valid] / morph_model[morph_valid]
        )
        morph_ratio_error[morph_valid] = (
            data_error[morph_valid] / morph_model[morph_valid]
        )

        bin_width = float(edges[1] - edges[0])
        offset = 0.10 * bin_width

        bottom.errorbar(
            centers[raw_valid] - offset,
            raw_ratio[raw_valid],
            yerr=raw_ratio_error[raw_valid],
            fmt="o",
            markersize=2.2,
            markerfacecolor="none",
            linewidth=0.7,
            capsize=0.0,
            color=sample_colors["raw_fit"],
            label="data / raw fit",
            zorder=3,
        )
        bottom.errorbar(
            centers[morph_valid] + offset,
            morph_ratio[morph_valid],
            yerr=morph_ratio_error[morph_valid],
            fmt="o",
            markersize=2.2,
            linewidth=0.7,
            capsize=0.0,
            color=sample_colors["data"],
            label="data / morphed fit",
            zorder=4,
        )
        bottom.axhline(
            1.0,
            color=sample_colors["morphed_fit"],
            linewidth=1.2,
            linestyle="--",
            zorder=1,
        )
        bottom.set_xlim(variable.low, variable.high)
        bottom.set_ylim(0.0, 2.0)
        bottom.set_xlabel(variable.label)
        bottom.set_ylabel("data / fit")
        bottom.grid(axis="y", alpha=0.25)
    # endfor

    handles, labels = axes[0, 0].get_legend_handles_labels()
    ratio_handles, ratio_labels = axes[1, 0].get_legend_handles_labels()
    handles = [*handles, *ratio_handles]
    labels = [*labels, *ratio_labels]

    fraction_drivers = ", ".join(
        variable.key for variable in fit_variables
    )
    primary_note = (
        "raw fraction-only fit is primary diagnostic"
        if summary.multiplicity_key == "one_photon"
        else "raw fraction-only fit shown as diagnostic comparison"
    )

    fig.suptitle(
        "DVCS+pi0 raw-template versus constrained-morph fits: "
        f"{period_label}, {multiplicity_label}\n"
        f"fraction drivers: {fraction_drivers}; {primary_note}; "
        rf"raw shared $f_{{\pi^0}}={summary.raw_f_pi0:.4f}"
        rf"\pm{summary.raw_f_pi0_err:.4f}$; "
        rf"morphed={summary.f_pi0:.4f}\pm{summary.f_pi0_err:.4f}$; "
        rf"$\Delta f_{{morph}}={summary.f_pi0-summary.raw_f_pi0:+.4f}$",
        fontsize=14,
        y=0.992,
    )
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.875),
        ncol=5,
        frameon=False,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.80))
    fig.savefig(
        path,
        dpi=args.fit_dpi,
        bbox_inches="tight",
    )
    plt.close(fig)


def plot_fraction_consistency(
    path: Path,
    period_label: str,
    multiplicity_label: str,
    summary: AdvancedSharedFit,
) -> None:
    """Compare raw-template and morphed fractions projection by projection."""
    fit_variables = fit_variables_for_category(
        summary.multiplicity_key
    )
    labels = [variable.label for variable in fit_variables]

    raw_values = np.asarray(
        [
            summary.variable_results[
                variable.key
            ].raw_independent_f_pi0
            for variable in fit_variables
        ],
        dtype=float,
    )
    raw_errors = np.asarray(
        [
            summary.variable_results[
                variable.key
            ].raw_independent_f_pi0_err
            for variable in fit_variables
        ],
        dtype=float,
    )
    morphed_values = np.asarray(
        [
            summary.variable_results[
                variable.key
            ].independent_f_pi0
            for variable in fit_variables
        ],
        dtype=float,
    )
    morphed_errors = np.asarray(
        [
            summary.variable_results[
                variable.key
            ].independent_f_pi0_err
            for variable in fit_variables
        ],
        dtype=float,
    )

    x = np.arange(len(labels), dtype=float)
    offset = 0.08

    fig, axis = plt.subplots(figsize=(12, 6))
    axis.errorbar(
        x - offset,
        raw_values,
        yerr=raw_errors,
        fmt="o",
        capsize=3,
        label="Raw-template independent fits",
    )
    axis.errorbar(
        x + offset,
        morphed_values,
        yerr=morphed_errors,
        fmt="s",
        capsize=3,
        label="Morphed independent fits",
    )

    axis.axhline(
        summary.raw_f_pi0,
        linewidth=1.5,
        linestyle="-.",
        label=(
            "Raw shared fit "
            f"({summary.raw_f_pi0:.3f})"
        ),
    )
    axis.axhline(
        summary.f_pi0,
        linewidth=1.5,
        linestyle="--",
        label=(
            "Morphed shared fit "
            f"({summary.f_pi0:.3f})"
        ),
    )

    if np.isfinite(summary.raw_f_pi0_err):
        axis.axhspan(
            summary.raw_f_pi0 - summary.raw_f_pi0_err,
            summary.raw_f_pi0 + summary.raw_f_pi0_err,
            alpha=0.08,
        )
    # endif
    if np.isfinite(summary.f_pi0_err):
        axis.axhspan(
            summary.f_pi0 - summary.f_pi0_err,
            summary.f_pi0 + summary.f_pi0_err,
            alpha=0.08,
        )
    # endif

    axis.set_xticks(x, labels, rotation=20, ha="right")
    axis.set_ylabel(r"Fitted $\pi^0$ fraction")
    axis.set_ylim(0.0, 1.0)
    axis.grid(alpha=0.25)
    axis.legend()
    axis.set_title(
        f"{period_label}: {multiplicity_label} fraction consistency; "
        rf"$\Delta f_{{morph}}="
        f"{summary.f_pi0-summary.raw_f_pi0:+.4f}$"
    )
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def flatten_shared_fit_row(
    period: PeriodConfig,
    multiplicity_key: str,
    multiplicity_label: str,
    summary: AdvancedSharedFit,
) -> Dict[str, object]:
    row: Dict[str, object] = {
        "period": period.key,
        "period_label": period.label,
        "multiplicity": multiplicity_key,
        "multiplicity_label": multiplicity_label,
        "fit_success": summary.success,
        "fit_message": summary.message,
        "f_pi0": summary.f_pi0,
        "f_pi0_error": summary.f_pi0_err,
        "deviance": summary.deviance,
        "ndf": summary.ndf,
        "deviance_per_ndf": (
            summary.deviance / summary.ndf
            if summary.ndf > 0
            else math.nan
        ),
        "coordinate_iterations": summary.coordinate_iterations,
        "raw_f_pi0": summary.raw_f_pi0,
        "raw_f_pi0_error": summary.raw_f_pi0_err,
        "raw_deviance": summary.raw_deviance,
        "raw_ndf": summary.raw_ndf,
        "raw_deviance_per_ndf": (
            summary.raw_deviance / summary.raw_ndf
            if summary.raw_ndf > 0
            else math.nan
        ),
        "delta_f_morph": summary.f_pi0 - summary.raw_f_pi0,
        "primary_one_photon_fraction": (
            summary.raw_f_pi0
            if multiplicity_key == "one_photon"
            else summary.f_pi0
        ),
    }

    fit_variables = fit_variables_for_category(
        summary.multiplicity_key
    )

    for variable in fit_variables:
        result = summary.variable_results[variable.key]
        prefix = variable.key
        row[f"{prefix}_raw_independent_f_pi0"] = (
            result.raw_independent_f_pi0
        )
        row[f"{prefix}_raw_independent_f_pi0_error"] = (
            result.raw_independent_f_pi0_err
        )
        row[f"{prefix}_raw_independent_deviance"] = (
            result.raw_independent_deviance
        )
        row[f"{prefix}_raw_shared_deviance"] = (
            result.raw_shared_deviance
        )
        row[f"{prefix}_delta_f_morph_independent"] = (
            result.independent_f_pi0
            - result.raw_independent_f_pi0
        )
        row[f"{prefix}_independent_f_pi0"] = (
            result.independent_f_pi0
        )
        row[f"{prefix}_independent_f_pi0_error"] = (
            result.independent_f_pi0_err
        )
        row[f"{prefix}_deviance"] = result.deviance
        row[f"{prefix}_ndf"] = result.ndf
        row[f"{prefix}_deviance_per_ndf"] = (
            result.deviance / result.ndf
            if result.ndf > 0
            else math.nan
        )
        row[f"{prefix}_model"] = result.model_name
        row[f"{prefix}_dvcs_nuisance_json"] = json.dumps(
            result.dvcs_nuisance,
            sort_keys=True,
        )
        row[f"{prefix}_pi0_nuisance_json"] = json.dumps(
            result.pi0_nuisance,
            sort_keys=True,
        )
        row[f"{prefix}_dvcs_divider"] = result.divider_dvcs
        row[f"{prefix}_pi0_divider"] = result.divider_pi0
    # endfor

    return row


def plot_all_period_fitted_pi0_fractions(
    path: Path,
    rows: Sequence[Mapping[str, object]],
) -> None:
    period_labels = [period.label for period in PERIODS]
    x = np.arange(len(period_labels), dtype=float)
    offsets = {
        "one_photon": -0.08,
        "two_photon": 0.08,
    }
    display_labels = {
        "one_photon": "One-photon events",
        "two_photon": "Two-photon events",
    }

    fig, axis = plt.subplots(figsize=(12, 7))
    for multiplicity_key in ("one_photon", "two_photon"):
        selected = []
        for period in PERIODS:
            match = next(
                (
                    row for row in rows
                    if row["period"] == period.key
                    and row["multiplicity"] == multiplicity_key
                    and bool(row["fit_success"])
                ),
                None,
            )
            selected.append(match)
        # endfor

        values = np.asarray(
            [
                row["f_pi0"] if row is not None else math.nan
                for row in selected
            ],
            dtype=float,
        )
        errors = np.asarray(
            [
                row["f_pi0_error"] if row is not None else math.nan
                for row in selected
            ],
            dtype=float,
        )

        axis.errorbar(
            x + offsets[multiplicity_key],
            values,
            yerr=errors,
            fmt="o",
            capsize=3,
            label=display_labels[multiplicity_key],
        )
    # endfor

    axis.set_xticks(x, period_labels, rotation=20, ha="right")
    axis.set_ylabel(r"Fitted $\pi^0$ fraction")
    axis.set_ylim(0.0, 1.0)
    axis.grid(alpha=0.25)
    axis.legend()
    axis.set_title(
        "Advanced shared four-projection template fits"
    )
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_template_fit_outputs(
    fit_dir: Path,
    fit_rows: Sequence[Mapping[str, object]],
    fit_payloads: Mapping[str, object],
) -> None:
    fit_dir.mkdir(parents=True, exist_ok=True)

    if fit_rows:
        fieldnames = list(fit_rows[0].keys())
        with open(
            fit_dir / "template_fit_results.csv",
            "w",
            newline="",
            encoding="utf-8",
        ) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=fieldnames,
            )
            writer.writeheader()
            writer.writerows(fit_rows)
        # endwith
    # endif

    with open(
        fit_dir / "template_fit_results.json",
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(
            {
                "fit_model": {
                    "shared_fraction": True,
                    "independent_variable_fits": True,
                    "raw_fraction_only_fit": True,
                    "raw_fraction_only_has_shape_nuisances": False,
                    "one_photon_primary_diagnostic": "raw_fraction_only",
                    "morphed_fit_role": (
                        "secondary model-sensitivity comparison for one-photon; "
                        "existing fit result retained for exact-two categories"
                    ),
                    "data_population": "three energy-ranked photon-entry categories",
                    "templates": [
                        "DVCSGEN",
                        "AAOGEN-as-epgamma pi0",
                    ],
                    "variable_models": {
                        category_key: {
                            variable.key: variable_model_spec(
                                variable.key,
                                category_key,
                            )
                            for variable in fit_variables_for_category(
                                category_key
                            )
                        }
                        for category_key in (
                            "one_photon",
                            "two_photon_more_energetic",
                            "two_photon_less_energetic",
                        )
                    },
                    "projection_likelihood_note": (
                        "The four one-dimensional projections contain the "
                        "same events and are combined as a composite "
                        "likelihood. Projection correlations are not included "
                        "in the curvature uncertainty."
                    ),
                },
                "rows": list(fit_rows),
                "fits": fit_payloads,
            },
            handle,
            indent=2,
        )
    # endwith

    if fit_rows:
        plot_all_period_fitted_pi0_fractions(
            fit_dir / "all_periods_fitted_pi0_fraction.png",
            fit_rows,
        )
    # endif



def run_data_shape_stage(
    periods: Sequence[PeriodConfig],
    args: argparse.Namespace,
    workers: int,
) -> Dict[str, object]:
    """
    Run 15 independent file tasks with no nested multiprocessing.
    """
    data_root = Path(args.output_dir) / "data"
    shape_dir = data_root / "shape_comparison"
    shape_dir.mkdir(parents=True, exist_ok=True)

    tasks = [
        (period, sample_key, sample_label, period_attribute)
        for period in periods
        for sample_key, sample_label, period_attribute
        in RAW_SHAPE_SAMPLE_INFO
    ]
    file_workers = max(
        1,
        min(int(workers), MAX_WORKERS, len(tasks)),
    )

    log(
        f"DATA SHAPE + MULTIPLICITY STAGE: {len(tasks)} files, "
        f"{file_workers} worker(s). Each ROOT file is read once."
    )

    results: Dict[str, Dict[str, Dict[str, object]]] = {
        period.key: {} for period in periods
    }

    if file_workers == 1:
        for period, sample_key, sample_label, period_attribute in tasks:
            key, returned_sample, payload = process_data_shape_sample(
                period,
                sample_key,
                sample_label,
                period_attribute,
                vars(args),
            )
            results[key][returned_sample] = payload
        # endfor
    else:
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=file_workers
        ) as executor:
            future_map = {
                executor.submit(
                    process_data_shape_sample,
                    period,
                    sample_key,
                    sample_label,
                    period_attribute,
                    vars(args),
                ): (period.key, sample_key)
                for (
                    period,
                    sample_key,
                    sample_label,
                    period_attribute,
                ) in tasks
            }
            for future in concurrent.futures.as_completed(future_map):
                key, returned_sample, payload = future.result()
                results[key][returned_sample] = payload
                log(
                    f"Completed {payload['period_label']} "
                    f"{payload['sample_label']}."
                )
            # endfor
        # endwith
    # endif

    period_order = {
        period.key: index for index, period in enumerate(PERIODS)
    }
    sample_order = {
        sample_key: index
        for index, (sample_key, _, _) in enumerate(
            RAW_SHAPE_SAMPLE_INFO
        )
    }

    multiplicity_rows: List[Dict[str, object]] = []
    audit: Dict[str, object] = {}
    fit_rows: List[Dict[str, object]] = []
    fit_payloads: Dict[str, object] = {}
    fit_dir = data_root / "template_fits"
    if not args.skip_data_template_fits:
        fit_dir.mkdir(parents=True, exist_ok=True)
    # endif

    for period in periods:
        sample_payloads = results[period.key]
        missing_samples = [
            sample_key
            for sample_key, _, _ in RAW_SHAPE_SAMPLE_INFO
            if sample_key not in sample_payloads
        ]
        if missing_samples:
            raise RuntimeError(
                f"{period.label}: missing processed samples: "
                + ", ".join(missing_samples)
            )
        # endif

        category_specs = (
            (
                "one_photon",
                "exactly one reconstructed-photon entry",
            ),
            (
                "two_photon_more_energetic",
                "exactly two entries: more energetic photon",
            ),
            (
                "two_photon_less_energetic",
                "exactly two entries: less energetic photon",
            ),
        )

        shape_plot_paths: Dict[str, Path] = {}
        for category_key, category_label in category_specs:
            category_plot_path = (
                shape_dir
                / f"{period.key}_{category_key}_shape_comparison.png"
            )
            shape_plot_paths[category_key] = category_plot_path
            plot_before_after_shape_canvas(
                category_plot_path,
                period.label,
                category_label,
                sample_payloads,
                category_key,
                display_variables_for_category(category_key),
                "shape comparisons",
                args,
            )
        # endfor

        period_fit_payload: Dict[str, object] = {}
        fit_stage_seconds = 0.0
        if not args.skip_data_template_fits:
            template_fit_start = time.perf_counter()
            fit_categories = category_specs

            for multiplicity_key, multiplicity_label in fit_categories:
                try:
                    summary, data_histograms = (
                        run_one_shared_template_fit(
                            period,
                            multiplicity_key,
                            multiplicity_label,
                            sample_payloads,
                            args,
                        )
                    )

                    fit_plot_path = (
                        fit_dir
                        / f"{period.key}_{multiplicity_key}"
                        "_template_fit.png"
                    )
                    plot_shared_template_fit(
                        fit_plot_path,
                        period.label,
                        multiplicity_label,
                        summary,
                        data_histograms,
                        args,
                    )

                    consistency_plot_path = (
                        fit_dir
                        / f"{period.key}_{multiplicity_key}"
                        "_fraction_consistency.png"
                    )
                    plot_fraction_consistency(
                        consistency_plot_path,
                        period.label,
                        multiplicity_label,
                        summary,
                    )

                    fit_row = flatten_shared_fit_row(
                        period,
                        multiplicity_key,
                        multiplicity_label,
                        summary,
                    )
                    fit_rows.append(fit_row)
                    period_fit_payload[multiplicity_key] = {
                        "plot": str(fit_plot_path),
                        "fraction_consistency_plot": str(
                            consistency_plot_path
                        ),
                        **shared_fit_to_dict(summary),
                    }

                    log(
                        f"{period.label} {multiplicity_label}: "
                        f"raw f_pi0={summary.raw_f_pi0:.6f} +/- "
                        f"{summary.raw_f_pi0_err:.6f}; "
                        f"morphed f_pi0={summary.f_pi0:.6f} +/- "
                        f"{summary.f_pi0_err:.6f}; "
                        f"Delta={summary.f_pi0-summary.raw_f_pi0:+.6f}; "
                        f"raw D/ndf={summary.raw_deviance:.2f}/"
                        f"{summary.raw_ndf}; "
                        f"morph D/ndf={summary.deviance:.2f}/"
                        f"{summary.ndf}"
                    )
                except Exception as exc:
                    failure = {
                        "success": False,
                        "message": str(exc),
                    }
                    period_fit_payload[multiplicity_key] = failure
                    if not args.allow_data_fit_failures:
                        raise
                    # endif
                    log(
                        f"WARNING: {period.label} {multiplicity_label} "
                        f"template fit failed: {exc}"
                    )
                # endtry
            # endfor
        # endif

        fit_payloads[period.key] = period_fit_payload

        period_audit: Dict[str, object] = {
            "plots": {
                key: str(value)
                for key, value in shape_plot_paths.items()
            },
            "template_fits": period_fit_payload,
            "fit_stage_seconds": fit_stage_seconds,
            "samples": {},
        }
        for sample_key, _, _ in RAW_SHAPE_SAMPLE_INFO:
            payload = sample_payloads[sample_key]
            multiplicity_rows.append(
                payload["photon_multiplicity"]
            )
            period_audit["samples"][sample_key] = {
                "sample_label": payload["sample_label"],
                "path": payload["path"],
                "tree_entries": payload["tree_entries"],
                "entries_read": payload["entries_read"],
                "resolved_branches": payload["resolved_branches"],
                "multiplicity_audit": payload[
                    "multiplicity_audit"
                ],
                "cutflow": payload["cutflow"],
                "photon_multiplicity": payload[
                    "photon_multiplicity"
                ],
                "pi0_mass_diagnostics": payload.get(
                    "pi0_mass_diagnostics",
                    {},
                ),
                **{
                    category_key: {
                        "before_histograms": serializable_histograms(
                            payload[category_key]["before_histograms"]
                        ),
                        "after_histograms": serializable_histograms(
                            payload[category_key]["after_histograms"]
                        ),
                    }
                    for category_key, _ in category_specs
                },
            }
        # endfor
        audit[period.key] = period_audit
    # endfor

    multiplicity_rows.sort(
        key=lambda row: (
            period_order.get(str(row["period"]), 999),
            sample_order.get(str(row["sample"]), 999),
        )
    )
    write_photon_multiplicity_outputs(
        data_root,
        multiplicity_rows,
    )

    if not args.skip_data_template_fits:
        write_template_fit_outputs(
            fit_dir,
            fit_rows,
            fit_payloads,
        )
    # endif

    with open(
        data_root / "shape_comparison_audit.json",
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(audit, handle, indent=2)
    # endwith

    for row in multiplicity_rows:
        log(
            f"{row['period_label']} {row['sample_label']}: "
            f"events by photon multiplicity "
            f"1={row['events_with_1_photon']:,}, "
            f"2={row['events_with_2_photons']:,}, "
            f"3={row['events_with_3_photons']:,}, "
            f"4={row['events_with_4_photons']:,}, "
            f"5={row['events_with_5_photons']:,}, "
            f"6={row['events_with_6_photons']:,}, "
            f">6={row['events_with_more_than_6_photons']:,}"
        )
    # endfor

    return {
        "periods": audit,
        "photon_multiplicity_rows": multiplicity_rows,
        "template_fit_rows": fit_rows,
        "template_fit_payloads": fit_payloads,
    }



def process_period_mc_only(
    period: PeriodConfig,
    args_dict: Mapping[str, object],
) -> Tuple[str, List[Dict[str, object]], Dict[str, object]]:
    args = argparse.Namespace(**args_dict)
    period_dir = (
        Path(args.output_dir)
        / "mc"
        / period.key
    )
    period_dir.mkdir(parents=True, exist_ok=True)

    opportunities, _, cutflow = read_opportunities(
        period.pi0_epg_mc,
        period.beam_energy_GeV,
        "AAOGEN epgamma opportunities",
        args,
        deduplicate=True,
        collect_candidates=False,
    )
    emit_cutflow_diagnostics(
        period_dir,
        period.label,
        {"aaogen_epgamma": cutflow},
    )

    epgg, epgg_diagnostics = read_epgg(
        period.pi0_epgg_mc,
        period.beam_energy_GeV,
        "AAOGEN",
        args,
    )
    match = match_truth_partners(
        opportunities,
        epgg,
        "mc",
        args,
    )

    summary = match.summary
    log(
        f"{period.label} AAOGEN: opportunities="
        f"{summary['opportunities']:,}, identity="
        f"{summary['identity_matches']:,}, valid daughter solution="
        f"{summary['valid_solution_opportunities']:,}, tag matched="
        f"{summary['tag_matched']:,}, partner above "
        f"{args.found_probe_E_min:g} GeV="
        f"{summary['partner_above_threshold']:,}, reconstructed="
        f"{summary['matched_opportunities']:,}"
    )

    plot_expected_probe_diagnostics(
        period_dir / "selected_opportunity_kinematics.png",
        period.label,
        {"AAOGEN": opportunities},
    )
    plot_matching_residuals_mc_only(
        period_dir / "matching_residuals.png",
        period.label,
        match,
    )

    official_rows: List[MCEfficiencyRow] = []
    for topology in ("FT", "FD"):
        mask = topology_mask(opportunities, topology)
        row = efficiency_row(
            period,
            topology,
            mask,
            match.matched,
        )
        official_rows.append(row)
        log(
            f"{period.label} {topology}: "
            f"N_reconstructed={row.reconstructed:,}, "
            f"N_opportunities={row.opportunities:,}, "
            f"epsilon_MC={row.efficiency:.6f} +/- "
            f"{row.efficiency_error:.6f}"
        )
    # endfor

    official_dicts = [asdict(row) for row in official_rows]
    write_dict_rows(
        period_dir / "integrated_mc_efficiency.csv",
        official_dicts,
    )

    topology_masks = {
        topology: topology_mask(opportunities, topology)
        for topology in ("FT", "FD")
    }

    energy_edges = np.linspace(
        args.probe_E_min,
        args.probe_E_max,
        args.energy_bins + 1,
    )
    theta_edges = np.linspace(
        args.ft_theta_min,
        args.fd_theta_max,
        args.theta_bins + 1,
    )
    phi_edges = np.linspace(
        0.0,
        2.0 * math.pi,
        args.phi_bins + 1,
    )

    dependence_specs = {
        "probe_energy": (
            opportunities.probe_E,
            energy_edges,
            "Predicted probe energy (GeV)",
        ),
        "probe_theta": (
            np.degrees(opportunities.probe_theta),
            theta_edges,
            "Predicted probe theta (deg)",
        ),
        "probe_phi": (
            np.mod(opportunities.probe_phi, 2.0 * math.pi),
            phi_edges,
            "Predicted probe phi (rad)",
        ),
    }

    binned_outputs: Dict[str, Dict[str, List[Dict[str, object]]]] = {}
    for name, (values, edges, label) in dependence_specs.items():
        topology_rows: Dict[str, List[Dict[str, object]]] = {}
        for topology, mask in topology_masks.items():
            rows = binned_efficiency(
                values,
                mask,
                match.matched,
                edges,
            )
            topology_rows[topology] = rows
            write_dict_rows(
                period_dir / f"mc_efficiency_vs_{name}_{topology.lower()}.csv",
                rows,
            )
        # endfor

        binned_outputs[name] = topology_rows
        plot_efficiency_dependence(
            period_dir / f"mc_efficiency_vs_{name}.png",
            f"{period.label}: AAOGEN photon efficiency",
            label,
            topology_rows,
            args.minimum_bin_opportunities,
        )
    # endfor

    sector_rows = plot_fd_sector_diagnostics(
        period_dir / "fd_sector_diagnostic_efficiency.png",
        period.label,
        opportunities,
        match.matched,
    )
    write_dict_rows(
        period_dir / "fd_sector_diagnostic_efficiency.csv",
        sector_rows,
    )

    metadata = {
        "period": {
            "key": period.key,
            "label": period.label,
            "beam_energy_GeV": period.beam_energy_GeV,
            "pi0_epg_mc": period.pi0_epg_mc,
            "pi0_epgg_mc": period.pi0_epgg_mc,
        },
        "arguments": vars(args),
        "opportunity_cutflow": cutflow,
        "epgammagamma_reconstruction": epgg_diagnostics,
        "matching": summary,
        "integrated_efficiencies": official_dicts,
        "fd_sector_diagnostics": sector_rows,
        "binned_efficiencies": binned_outputs,
    }
    with open(
        period_dir / "mc_efficiency_metadata.json",
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(metadata, handle, indent=2)
    # endwith

    return period.key, official_dicts, metadata



def simple_data_event_key(
    runnum: np.ndarray,
    eventnum: np.ndarray,
) -> np.ndarray:
    """
    Build an exact structured (runnum, evnum) key for data entries.
    """
    runnum = np.asarray(runnum, dtype=np.int64)
    eventnum = np.asarray(eventnum, dtype=np.int64)
    keys = np.empty(
        runnum.size,
        dtype=[("runnum", "<i8"), ("eventnum", "<i8")],
    )
    keys["runnum"] = runnum
    keys["eventnum"] = eventnum
    return keys


def simple_binomial_error(
    matched: int,
    expected: int,
) -> float:
    if expected <= 0:
        return math.nan
    # endif
    efficiency = matched / expected
    return float(
        math.sqrt(
            max(efficiency * (1.0 - efficiency), 0.0)
            / expected
        )
    )



SIMPLE_DIAGNOSTIC_SPECS: Mapping[str, Tuple[int, float, float]] = {
    "low_energy": (100, 0.0, 2.0),
    "low_theta_deg": (100, 0.0, 40.0),
    "predicted_theta_deg": (100, 0.0, 40.0),
    "mgg_before": (120, 0.0, 0.40),
    "high_energy_before": (100, 2.0, 9.5),
    "high_energy_after": (100, 2.0, 9.5),
}


def simple_empty_diagnostic_histograms() -> Dict[str, Dict[str, Dict[str, np.ndarray]]]:
    result: Dict[str, Dict[str, Dict[str, np.ndarray]]] = {}
    for region in ("FT", "FD"):
        result[region] = {}
        for name, (bins, low, high) in SIMPLE_DIAGNOSTIC_SPECS.items():
            edges = np.linspace(low, high, bins + 1)
            result[region][name] = {
                "edges": edges,
                "counts": np.zeros(bins, dtype=np.int64),
            }
        # endfor
    # endfor
    return result


def simple_fill_diagnostic_histogram(
    histograms: Dict[str, Dict[str, Dict[str, np.ndarray]]],
    region: str,
    name: str,
    values: np.ndarray,
) -> None:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return
    # endif

    edges = np.asarray(
        histograms[region][name]["edges"],
        dtype=float,
    )
    counts, _ = np.histogram(values, bins=edges)
    histograms[region][name]["counts"] += counts.astype(
        np.int64,
        copy=False,
    )


def simple_finalize_diagnostic_histograms(
    histograms: Mapping[str, Mapping[str, Mapping[str, np.ndarray]]],
) -> Dict[str, object]:
    payload: Dict[str, object] = {}
    for region in ("FT", "FD"):
        payload[region] = {}
        for name in SIMPLE_DIAGNOSTIC_SPECS:
            edges = np.asarray(
                histograms[region][name]["edges"],
                dtype=float,
            )
            counts = np.asarray(
                histograms[region][name]["counts"],
                dtype=np.int64,
            )
            total = int(np.sum(counts))
            unit_area = (
                counts.astype(float) / total
                if total > 0
                else np.zeros_like(counts, dtype=float)
            )
            payload[region][name] = {
                "edges": edges,
                "counts": counts,
                "unit_area": unit_area,
                "in_range_entries": total,
            }
        # endfor
    # endfor
    return payload


def simple_region_key(detector_code: int) -> Optional[str]:
    if int(detector_code) == 0:
        return "FT"
    # endif
    if int(detector_code) == 1:
        return "FD"
    # endif
    return None

def process_period_simple_data_efficiency(
    period: PeriodConfig,
    args_dict: Mapping[str, object],
) -> Tuple[str, Dict[str, object]]:
    """
    Temporary direct data efficiency estimator requested by --simple.

    Denominator:
        every p2_p < 2 GeV photon entry is one expected pi0 opportunity.

    Numerator:
        the exact same (runnum, evnum) must contain at least one p2_p > 2 GeV
        photon that forms 0.11 < M_gamma_gamma < 0.16 GeV with the low-energy
        photon.

    Predicted detector region is still obtained from the missing-partner
    four-vector reconstructed from beam+target-e-p-low_gamma.
    """
    start_time = time.perf_counter()
    args = argparse.Namespace(**args_dict)
    path = period.epg_data

    total_entries, available_keys = require_tree(path)
    logical_names = (
        "runnum", "eventnum",
        "e_p", "e_theta", "e_phi",
        "p_p", "p_theta", "p_phi",
        "p2_p", "p2_theta", "p2_phi",
    )
    resolved = resolve_branches_from_keys(path, available_keys, logical_names)
    expressions = sorted(set(resolved.values()))

    entry_stop = (
        min(total_entries, int(args.max_events))
        if args.max_events is not None
        else total_entries
    )

    low_key_chunks: List[np.ndarray] = []
    low_detector_chunks: List[np.ndarray] = []
    low_sector_chunks: List[np.ndarray] = []
    low_probe_E_chunks: List[np.ndarray] = []
    low_probe_theta_deg_chunks: List[np.ndarray] = []
    low_E_chunks: List[np.ndarray] = []
    low_theta_chunks: List[np.ndarray] = []
    low_phi_chunks: List[np.ndarray] = []

    high_key_chunks: List[np.ndarray] = []
    high_E_chunks: List[np.ndarray] = []
    high_theta_chunks: List[np.ndarray] = []
    high_phi_chunks: List[np.ndarray] = []

    entries_read = 0
    finite_identity_entries = 0
    finite_low_kinematics = 0
    low_entries_total = 0
    high_entries_total = 0
    exactly_2GeV_entries = 0

    for arrays in uproot.iterate(
        f"{path}:{TREE_NAME}",
        expressions=expressions,
        step_size=args.step_size,
        entry_stop=entry_stop,
        library="np",
    ):
        if not arrays:
            continue
        # endif

        chunk_size = len(next(iter(arrays.values())))
        entries_read += chunk_size
        values = {
            name: np.asarray(arrays[resolved[name]], dtype=float)
            for name in logical_names
        }

        run_f = values["runnum"]
        event_f = values["eventnum"]
        photon_E = values["p2_p"]
        photon_theta = values["p2_theta"]
        photon_phi = values["p2_phi"]

        valid_identity = np.isfinite(run_f) & np.isfinite(event_f)
        finite_identity_entries += int(np.count_nonzero(valid_identity))

        finite_photon = (
            valid_identity
            & np.isfinite(photon_E)
            & np.isfinite(photon_theta)
            & np.isfinite(photon_phi)
        )
        low_mask = finite_photon & (photon_E < 2.0)
        high_mask = finite_photon & (photon_E > 2.0)
        exactly_mask = finite_photon & (photon_E == 2.0)

        low_entries_total += int(np.count_nonzero(low_mask))
        high_entries_total += int(np.count_nonzero(high_mask))
        exactly_2GeV_entries += int(np.count_nonzero(exactly_mask))

        if np.any(high_mask):
            high_key_chunks.append(simple_data_event_key(
                np.rint(run_f[high_mask]).astype(np.int64),
                np.rint(event_f[high_mask]).astype(np.int64),
            ))
            high_E_chunks.append(np.asarray(photon_E[high_mask], dtype=float))
            high_theta_chunks.append(np.asarray(photon_theta[high_mask], dtype=float))
            high_phi_chunks.append(np.asarray(photon_phi[high_mask], dtype=float))
        # endif

        if not np.any(low_mask):
            continue
        # endif

        low_required_finite = low_mask.copy()
        for name in ("e_p", "e_theta", "e_phi", "p_p", "p_theta", "p_phi"):
            low_required_finite &= np.isfinite(values[name])
        # endfor
        finite_low_kinematics += int(np.count_nonzero(low_required_finite))
        if not np.any(low_required_finite):
            continue
        # endif

        probe = reconstruct_probe(
            period.beam_energy_GeV,
            values["e_p"][low_required_finite],
            values["e_theta"][low_required_finite],
            values["e_phi"][low_required_finite],
            values["p_p"][low_required_finite],
            values["p_theta"][low_required_finite],
            values["p_phi"][low_required_finite],
            photon_E[low_required_finite],
            photon_theta[low_required_finite],
            photon_phi[low_required_finite],
        )
        probe_theta_deg = np.degrees(
            np.asarray(probe["theta"], dtype=float)
        )
        predicted_detector, predicted_sector = classify_probe(
            probe_theta_deg,
            np.asarray(probe["phi"], dtype=float),
            args,
        )

        low_key_chunks.append(simple_data_event_key(
            np.rint(run_f[low_required_finite]).astype(np.int64),
            np.rint(event_f[low_required_finite]).astype(np.int64),
        ))
        low_detector_chunks.append(np.asarray(predicted_detector, dtype=np.int16))
        low_sector_chunks.append(np.asarray(predicted_sector, dtype=np.int16))
        low_probe_E_chunks.append(
            np.asarray(probe["E"], dtype=float)
        )
        low_probe_theta_deg_chunks.append(
            np.asarray(probe_theta_deg, dtype=float)
        )
        low_E_chunks.append(
            np.asarray(photon_E[low_required_finite], dtype=float)
        )
        low_theta_chunks.append(np.asarray(photon_theta[low_required_finite], dtype=float))
        low_phi_chunks.append(np.asarray(photon_phi[low_required_finite], dtype=float))
    # endfor

    if low_key_chunks:
        low_keys = np.concatenate(low_key_chunks)
        low_detector = np.concatenate(low_detector_chunks)
        low_sector = np.concatenate(low_sector_chunks)
        low_probe_E = (
            np.concatenate(low_probe_E_chunks)
            if low_probe_E_chunks
            else np.empty(0, dtype=float)
        )
        low_probe_theta_deg = (
            np.concatenate(low_probe_theta_deg_chunks)
            if low_probe_theta_deg_chunks
            else np.empty(0, dtype=float)
        )
        low_E = np.concatenate(low_E_chunks)
        low_theta = np.concatenate(low_theta_chunks)
        low_phi = np.concatenate(low_phi_chunks)

        if low_probe_E.size != low_E.size:
            raise RuntimeError(
                f"{period.label} simple DATA diagnostic bookkeeping error: "
                f"probe-E entries={low_probe_E.size}, low tags={low_E.size}."
            )
        # endif
        if low_probe_theta_deg.size != low_E.size:
            raise RuntimeError(
                f"{period.label} simple DATA diagnostic bookkeeping error: "
                f"probe-theta entries={low_probe_theta_deg.size}, "
                f"low tags={low_E.size}."
            )
        # endif
    else:
        low_keys = np.empty(0, dtype=[("runnum", "<i8"), ("eventnum", "<i8")])
        low_detector = np.empty(0, dtype=np.int16)
        low_sector = np.empty(0, dtype=np.int16)
        low_probe_E = np.empty(0, dtype=float)
        low_probe_theta_deg = np.empty(0, dtype=float)
        low_E = np.empty(0, dtype=float)
        low_theta = np.empty(0, dtype=float)
        low_phi = np.empty(0, dtype=float)
    # endif

    if high_key_chunks:
        high_keys = np.concatenate(high_key_chunks)
        high_E = np.concatenate(high_E_chunks)
        high_theta = np.concatenate(high_theta_chunks)
        high_phi = np.concatenate(high_phi_chunks)
        order = np.argsort(high_keys, kind="stable")
        high_keys = high_keys[order]
        high_E = high_E[order]
        high_theta = high_theta[order]
        high_phi = high_phi[order]
    else:
        high_keys = np.empty(0, dtype=[("runnum", "<i8"), ("eventnum", "<i8")])
        high_E = np.empty(0, dtype=float)
        high_theta = np.empty(0, dtype=float)
        high_phi = np.empty(0, dtype=float)
    # endif

    diagnostic_histograms = simple_empty_diagnostic_histograms()
    pair_diagnostic_values: Dict[str, Dict[str, List[np.ndarray]]] = {
        region: {
            "mgg_before": [],
            "high_energy_before": [],
            "high_energy_after": [],
        }
        for region in ("FT", "FD")
    }

    for region, region_mask in (
        ("FT", low_detector == 0),
        ("FD", low_detector == 1),
    ):
        simple_fill_diagnostic_histogram(
            diagnostic_histograms,
            region,
            "low_energy",
            low_E[region_mask],
        )
        simple_fill_diagnostic_histogram(
            diagnostic_histograms,
            region,
            "low_theta_deg",
            np.degrees(low_theta[region_mask]),
        )
        simple_fill_diagnostic_histogram(
            diagnostic_histograms,
            region,
            "predicted_theta_deg",
            low_probe_theta_deg[region_mask],
        )
    # endfor

    same_event_high_exists = np.zeros(low_keys.size, dtype=bool)
    pi0_mass_partner_found = np.zeros(low_keys.size, dtype=bool)

    if low_keys.size > 0 and high_keys.size > 0:
        left = np.searchsorted(high_keys, low_keys, side="left")
        right = np.searchsorted(high_keys, low_keys, side="right")
        same_event_high_exists = right > left

        mlow2 = float(PI0_MASS_WINDOW[0]) ** 2
        mhigh2 = float(PI0_MASS_WINDOW[1]) ** 2

        for i in np.flatnonzero(same_event_high_exists):
            start = int(left[i])
            stop = int(right[i])
            hE = high_E[start:stop]
            htheta = high_theta[start:stop]
            hphi = high_phi[start:stop]

            cos_open = (
                math.cos(float(low_theta[i])) * np.cos(htheta)
                + math.sin(float(low_theta[i])) * np.sin(htheta)
                * np.cos(float(low_phi[i]) - hphi)
            )
            cos_open = np.clip(cos_open, -1.0, 1.0)
            m2 = 2.0 * float(low_E[i]) * hE * (1.0 - cos_open)
            finite_mass = np.isfinite(m2)
            passing_pair = (
                finite_mass
                & (m2 > mlow2)
                & (m2 < mhigh2)
            )

            region = simple_region_key(int(low_detector[i]))
            if region is not None:
                pair_diagnostic_values[region][
                    "mgg_before"
                ].append(
                    np.sqrt(
                        np.clip(
                            m2[finite_mass],
                            0.0,
                            None,
                        )
                    )
                )
                pair_diagnostic_values[region][
                    "high_energy_before"
                ].append(hE[finite_mass])
                if np.any(passing_pair):
                    pair_diagnostic_values[region][
                        "high_energy_after"
                    ].append(hE[passing_pair])
                # endif
            # endif

            if np.any(passing_pair):
                pi0_mass_partner_found[i] = True
            # endif
        # endfor
    # endif

    region_specs = [
        ("FT", 0, 0),
        *[(f"FD sector {sector}", 1, sector) for sector in range(1, 7)],
    ]

    for region in ("FT", "FD"):
        for name in (
            "mgg_before",
            "high_energy_before",
            "high_energy_after",
        ):
            chunks = pair_diagnostic_values[region][name]
            if chunks:
                values = np.concatenate(chunks)
                simple_fill_diagnostic_histogram(
                    diagnostic_histograms,
                    region,
                    name,
                    values,
                )
            # endif
        # endfor
    # endfor

    rows: List[Dict[str, object]] = []
    for region_label, detector_code, sector_code in region_specs:
        region_mask = low_detector == detector_code
        if detector_code == 1:
            region_mask &= low_sector == sector_code
        # endif

        expected = int(np.count_nonzero(region_mask))
        same_event = int(np.count_nonzero(region_mask & same_event_high_exists))
        found = int(np.count_nonzero(region_mask & pi0_mass_partner_found))
        efficiency = found / expected if expected > 0 else math.nan
        rows.append({
            "period": period.key,
            "period_label": period.label,
            "region": region_label,
            "detector_code": detector_code,
            "sector": sector_code,
            "expected_low_energy_tags": expected,
            "same_event_high_energy_partner": same_event,
            "matched_pi0_mass_partner": found,
            "efficiency": efficiency,
            "efficiency_error": simple_binomial_error(found, expected),
        })
    # endfor

    classified_mask = (
        (low_detector == 0)
        | ((low_detector == 1) & (low_sector >= 1) & (low_sector <= 6))
    )

    payload: Dict[str, object] = {
        "period": period.key,
        "period_label": period.label,
        "beam_energy_GeV": period.beam_energy_GeV,
        "path": path,
        "tree_entries": total_entries,
        "entries_read": entries_read,
        "threshold_GeV": 2.0,
        "pi0_mass_window_GeV": [float(PI0_MASS_WINDOW[0]), float(PI0_MASS_WINDOW[1])],
        "low_definition": "p2_p < 2 GeV",
        "high_definition": "p2_p > 2 GeV",
        "matching": "same (runnum,eventnum) and 0.11<M_gamma_gamma<0.16 GeV",
        "finite_identity_entries": finite_identity_entries,
        "low_entries_total": low_entries_total,
        "high_entries_total": high_entries_total,
        "exactly_2GeV_entries_excluded": exactly_2GeV_entries,
        "low_entries_with_finite_reconstruction_kinematics": finite_low_kinematics,
        "low_opportunities_built": int(low_keys.size),
        "high_energy_entries_stored": int(high_keys.size),
        "low_opportunities_with_same_event_high_partner": int(np.count_nonzero(same_event_high_exists)),
        "low_opportunities_with_pi0_mass_partner": int(np.count_nonzero(pi0_mass_partner_found)),
        "classified_low_opportunities": int(np.count_nonzero(classified_mask)),
        "unclassified_low_opportunities": int(np.count_nonzero(~classified_mask)),
        "low_opportunities_with_reconstructed_probe_E_gt_2": int(np.count_nonzero(np.isfinite(low_probe_E) & (low_probe_E > 2.0))),
        "low_opportunities_with_reconstructed_probe_E_le_2": int(np.count_nonzero(np.isfinite(low_probe_E) & (low_probe_E <= 2.0))),
        "rows": rows,
        "diagnostic_histograms": (
            simple_finalize_diagnostic_histograms(
                diagnostic_histograms
            )
        ),
        "elapsed_seconds": elapsed_seconds(start_time),
    }
    return period.key, payload


def print_simple_data_efficiency_block(
    title: str,
    rows: Sequence[Mapping[str, object]],
) -> None:
    log("=" * 72)
    log(f"SIMPLE DATA EFFICIENCY: {title}")
    log(
        "Definition: each p2_p<2 GeV data photon is one expected pi0 "
        "opportunity; success = same (runnum,evnum) contains p2_p>2 GeV "
        "with 0.11<Mgg<0.16 GeV."
    )
    log("Region is the PREDICTED missing-partner FT/FD sector.")
    for row in rows:
        expected = int(row["expected_low_energy_tags"])
        matched = int(row["matched_pi0_mass_partner"])
        same_event = int(row["same_event_high_energy_partner"])
        efficiency = float(row["efficiency"])
        error = float(row["efficiency_error"])
        if np.isfinite(efficiency):
            log(
                f"  {str(row['region']):<11s}: "
                f"{matched:,} / {expected:,} = "
                f"{100.0 * efficiency:7.3f}% +/- {100.0 * error:6.3f}% "
                f"(same-event >2 GeV before Mgg: {same_event:,})"
            )
        else:
            log(
                f"  {str(row['region']):<11s}: "
                f"{matched:,} / {expected:,} = undefined "
                f"(same-event >2 GeV before Mgg: {same_event:,})"
            )
        # endif
    # endfor
    log("=" * 72)



def write_simple_rows_csv(
    path: Path,
    rows: Sequence[Mapping[str, object]],
) -> None:
    """
    Write the compact --simple efficiency table without relying on any
    non-existent global CSV helper.
    """
    path.parent.mkdir(parents=True, exist_ok=True)

    if not rows:
        with open(path, "w", encoding="utf-8", newline=""):
            pass
        # endwith
        return
    # endif

    fieldnames: List[str] = []
    for row in rows:
        for key in row.keys():
            key_string = str(key)
            if key_string not in fieldnames:
                fieldnames.append(key_string)
            # endif
        # endfor
    # endfor

    with open(
        path,
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: row.get(key, "")
                    for key in fieldnames
                }
            )


def simple_json_default(value: object) -> object:
    """
    JSON serializer used only by the temporary --simple audit.

    Convert the NumPy scalar/array and Path types that can occur in the
    metadata payload to ordinary JSON-compatible Python objects.
    """
    if isinstance(value, np.integer):
        return int(value)
    # endif
    if isinstance(value, np.floating):
        return float(value)
    # endif
    if isinstance(value, np.bool_):
        return bool(value)
    # endif
    if isinstance(value, np.ndarray):
        return value.tolist()
    # endif
    if isinstance(value, Path):
        return str(value)
    # endif

    raise TypeError(
        f"Object of type {type(value).__name__} is not JSON serializable"
    )

        # endfor
    # endwith


def simple_mc_event_key(
    arrays_by_logical_name: Mapping[str, np.ndarray],
    mask: np.ndarray,
    decimals: int,
) -> np.ndarray:
    """
    Build the established reconstructed-MC event identity from electron and
    proton kinematics for the selected entries.
    """
    selected = {
        name: np.asarray(arrays_by_logical_name[name], dtype=float)[mask]
        for name in (
            "e_p",
            "e_theta",
            "e_phi",
            "p_p",
            "p_theta",
            "p_phi",
        )
    }
    signature, valid = mc_identity_signature(
        selected,
        decimals,
    )
    if not np.all(valid):
        raise RuntimeError(
            "simple_mc_event_key received non-finite electron/proton "
            "kinematics after its finite-input mask."
        )
    # endif
    return signature


def process_period_simple_mc_efficiency(
    period: PeriodConfig,
    args_dict: Mapping[str, object],
) -> Tuple[str, Dict[str, object]]:
    """
    Apply the exact same temporary reconstructed-photon algorithm used by
    --simple data to the reconstructed AAOGEN epgamma sample.

    The only unavoidable difference is event identity:
      data -> exact (runnum, evnum);
      MC   -> quantized electron+proton kinematic identity.
    """
    start_time = time.perf_counter()
    args = argparse.Namespace(**args_dict)
    path = period.pi0_epg_mc

    total_entries, available_keys = require_tree(path)
    logical_names = (
        "e_p", "e_theta", "e_phi",
        "p_p", "p_theta", "p_phi",
        "p2_p", "p2_theta", "p2_phi",
    )
    resolved = resolve_branches_from_keys(
        path,
        available_keys,
        logical_names,
    )
    expressions = sorted(set(resolved.values()))

    entry_stop = (
        min(total_entries, int(args.max_events))
        if args.max_events is not None
        else total_entries
    )

    low_key_chunks: List[np.ndarray] = []
    low_detector_chunks: List[np.ndarray] = []
    low_sector_chunks: List[np.ndarray] = []
    low_probe_E_chunks: List[np.ndarray] = []
    low_probe_theta_deg_chunks: List[np.ndarray] = []
    low_E_chunks: List[np.ndarray] = []
    low_theta_chunks: List[np.ndarray] = []
    low_phi_chunks: List[np.ndarray] = []

    high_key_chunks: List[np.ndarray] = []
    high_E_chunks: List[np.ndarray] = []
    high_theta_chunks: List[np.ndarray] = []
    high_phi_chunks: List[np.ndarray] = []

    entries_read = 0
    finite_identity_entries = 0
    finite_low_kinematics = 0
    low_entries_total = 0
    high_entries_total = 0
    exactly_2GeV_entries = 0

    for arrays in uproot.iterate(
        f"{path}:{TREE_NAME}",
        expressions=expressions,
        step_size=args.step_size,
        entry_stop=entry_stop,
        library="np",
    ):
        if not arrays:
            continue
        # endif

        chunk_size = len(next(iter(arrays.values())))
        entries_read += chunk_size
        values = {
            name: np.asarray(arrays[resolved[name]], dtype=float)
            for name in logical_names
        }

        identity_finite = np.ones(chunk_size, dtype=bool)
        for name in (
            "e_p", "e_theta", "e_phi",
            "p_p", "p_theta", "p_phi",
        ):
            identity_finite &= np.isfinite(values[name])
        # endfor
        finite_identity_entries += int(
            np.count_nonzero(identity_finite)
        )

        photon_E = values["p2_p"]
        photon_theta = values["p2_theta"]
        photon_phi = values["p2_phi"]

        finite_photon = (
            identity_finite
            & np.isfinite(photon_E)
            & np.isfinite(photon_theta)
            & np.isfinite(photon_phi)
        )
        low_mask = finite_photon & (photon_E < 2.0)
        high_mask = finite_photon & (photon_E > 2.0)
        exactly_mask = finite_photon & (photon_E == 2.0)

        low_entries_total += int(np.count_nonzero(low_mask))
        high_entries_total += int(np.count_nonzero(high_mask))
        exactly_2GeV_entries += int(np.count_nonzero(exactly_mask))

        if np.any(high_mask):
            high_key_chunks.append(
                simple_mc_event_key(
                    values,
                    high_mask,
                    args.mc_signature_decimals,
                )
            )
            high_E_chunks.append(
                np.asarray(photon_E[high_mask], dtype=float)
            )
            high_theta_chunks.append(
                np.asarray(photon_theta[high_mask], dtype=float)
            )
            high_phi_chunks.append(
                np.asarray(photon_phi[high_mask], dtype=float)
            )
        # endif

        if not np.any(low_mask):
            continue
        # endif

        low_required_finite = low_mask.copy()
        finite_low_kinematics += int(
            np.count_nonzero(low_required_finite)
        )

        probe = reconstruct_probe(
            period.beam_energy_GeV,
            values["e_p"][low_required_finite],
            values["e_theta"][low_required_finite],
            values["e_phi"][low_required_finite],
            values["p_p"][low_required_finite],
            values["p_theta"][low_required_finite],
            values["p_phi"][low_required_finite],
            photon_E[low_required_finite],
            photon_theta[low_required_finite],
            photon_phi[low_required_finite],
        )
        probe_theta_deg = np.degrees(
            np.asarray(probe["theta"], dtype=float)
        )
        predicted_detector, predicted_sector = classify_probe(
            probe_theta_deg,
            np.asarray(probe["phi"], dtype=float),
            args,
        )

        low_key_chunks.append(
            simple_mc_event_key(
                values,
                low_required_finite,
                args.mc_signature_decimals,
            )
        )
        low_detector_chunks.append(
            np.asarray(predicted_detector, dtype=np.int16)
        )
        low_sector_chunks.append(
            np.asarray(predicted_sector, dtype=np.int16)
        )
        low_probe_E_chunks.append(
            np.asarray(probe["E"], dtype=float)
        )
        low_probe_theta_deg_chunks.append(
            np.asarray(probe_theta_deg, dtype=float)
        )
        low_E_chunks.append(
            np.asarray(photon_E[low_required_finite], dtype=float)
        )
        low_theta_chunks.append(
            np.asarray(photon_theta[low_required_finite], dtype=float)
        )
        low_phi_chunks.append(
            np.asarray(photon_phi[low_required_finite], dtype=float)
        )
    # endfor

    key_dtype = np.dtype([("hash1", "<u8"), ("hash2", "<u8")])

    if low_key_chunks:
        low_keys = np.concatenate(low_key_chunks)
        low_detector = np.concatenate(low_detector_chunks)
        low_sector = np.concatenate(low_sector_chunks)
        low_probe_E = (
            np.concatenate(low_probe_E_chunks)
            if low_probe_E_chunks
            else np.empty(0, dtype=float)
        )
        low_probe_theta_deg = (
            np.concatenate(low_probe_theta_deg_chunks)
            if low_probe_theta_deg_chunks
            else np.empty(0, dtype=float)
        )
        low_E = np.concatenate(low_E_chunks)
        low_theta = np.concatenate(low_theta_chunks)
        low_phi = np.concatenate(low_phi_chunks)

        if low_probe_E.size != low_E.size:
            raise RuntimeError(
                f"{period.label} simple AAOGEN diagnostic bookkeeping error: "
                f"probe-E entries={low_probe_E.size}, low tags={low_E.size}."
            )
        # endif
        if low_probe_theta_deg.size != low_E.size:
            raise RuntimeError(
                f"{period.label} simple AAOGEN diagnostic bookkeeping error: "
                f"probe-theta entries={low_probe_theta_deg.size}, "
                f"low tags={low_E.size}."
            )
        # endif
    else:
        low_keys = np.empty(0, dtype=key_dtype)
        low_detector = np.empty(0, dtype=np.int16)
        low_sector = np.empty(0, dtype=np.int16)
        low_probe_E = np.empty(0, dtype=float)
        low_probe_theta_deg = np.empty(0, dtype=float)
        low_E = np.empty(0, dtype=float)
        low_theta = np.empty(0, dtype=float)
        low_phi = np.empty(0, dtype=float)
    # endif

    if high_key_chunks:
        high_keys = np.concatenate(high_key_chunks)
        high_E = np.concatenate(high_E_chunks)
        high_theta = np.concatenate(high_theta_chunks)
        high_phi = np.concatenate(high_phi_chunks)

        high_order = np.argsort(
            high_keys,
            kind="stable",
        )
        high_keys = high_keys[high_order]
        high_E = high_E[high_order]
        high_theta = high_theta[high_order]
        high_phi = high_phi[high_order]
    else:
        high_keys = np.empty(0, dtype=key_dtype)
        high_E = np.empty(0, dtype=float)
        high_theta = np.empty(0, dtype=float)
        high_phi = np.empty(0, dtype=float)
    # endif

    diagnostic_histograms = simple_empty_diagnostic_histograms()
    pair_diagnostic_values: Dict[str, Dict[str, List[np.ndarray]]] = {
        region: {
            "mgg_before": [],
            "high_energy_before": [],
            "high_energy_after": [],
        }
        for region in ("FT", "FD")
    }

    for region, region_mask in (
        ("FT", low_detector == 0),
        ("FD", low_detector == 1),
    ):
        simple_fill_diagnostic_histogram(
            diagnostic_histograms,
            region,
            "low_energy",
            low_E[region_mask],
        )
        simple_fill_diagnostic_histogram(
            diagnostic_histograms,
            region,
            "low_theta_deg",
            np.degrees(low_theta[region_mask]),
        )
        simple_fill_diagnostic_histogram(
            diagnostic_histograms,
            region,
            "predicted_theta_deg",
            low_probe_theta_deg[region_mask],
        )
    # endfor

    same_event_high_exists = np.zeros(
        low_keys.size,
        dtype=bool,
    )
    pi0_mass_partner_found = np.zeros(
        low_keys.size,
        dtype=bool,
    )

    if low_keys.size > 0 and high_keys.size > 0:
        left_indices = np.searchsorted(
            high_keys,
            low_keys,
            side="left",
        )
        right_indices = np.searchsorted(
            high_keys,
            low_keys,
            side="right",
        )
        same_event_high_exists = right_indices > left_indices

        candidate_low_indices = np.flatnonzero(
            same_event_high_exists
        )
        mass_low_sq = float(PI0_MASS_WINDOW[0]) ** 2
        mass_high_sq = float(PI0_MASS_WINDOW[1]) ** 2

        for low_index in candidate_low_indices:
            start = int(left_indices[low_index])
            stop = int(right_indices[low_index])

            candidate_E = high_E[start:stop]
            candidate_theta = high_theta[start:stop]
            candidate_phi = high_phi[start:stop]

            cosine_opening = (
                math.cos(float(low_theta[low_index]))
                * np.cos(candidate_theta)
                + math.sin(float(low_theta[low_index]))
                * np.sin(candidate_theta)
                * np.cos(
                    float(low_phi[low_index]) - candidate_phi
                )
            )
            cosine_opening = np.clip(
                cosine_opening,
                -1.0,
                1.0,
            )

            mass_squared = (
                2.0
                * float(low_E[low_index])
                * candidate_E
                * (1.0 - cosine_opening)
            )
            finite_mass = np.isfinite(mass_squared)
            passing_pair = (
                finite_mass
                & (mass_squared > mass_low_sq)
                & (mass_squared < mass_high_sq)
            )

            region = simple_region_key(
                int(low_detector[low_index])
            )
            if region is not None:
                pair_diagnostic_values[region][
                    "mgg_before"
                ].append(
                    np.sqrt(
                        np.clip(
                            mass_squared[finite_mass],
                            0.0,
                            None,
                        )
                    )
                )
                pair_diagnostic_values[region][
                    "high_energy_before"
                ].append(candidate_E[finite_mass])
                if np.any(passing_pair):
                    pair_diagnostic_values[region][
                        "high_energy_after"
                    ].append(candidate_E[passing_pair])
                # endif
            # endif

            if np.any(passing_pair):
                pi0_mass_partner_found[low_index] = True
            # endif
        # endfor
    # endif

    for region in ("FT", "FD"):
        for name in (
            "mgg_before",
            "high_energy_before",
            "high_energy_after",
        ):
            chunks = pair_diagnostic_values[region][name]
            if chunks:
                values = np.concatenate(chunks)
                simple_fill_diagnostic_histogram(
                    diagnostic_histograms,
                    region,
                    name,
                    values,
                )
            # endif
        # endfor
    # endfor

    rows: List[Dict[str, object]] = []
    for region_label, detector_code, sector_code in [
        ("FT", 0, 0),
        *[
            (f"FD sector {sector}", 1, sector)
            for sector in range(1, 7)
        ],
    ]:
        region_mask = low_detector == detector_code
        if detector_code == 1:
            region_mask &= low_sector == sector_code
        # endif

        expected = int(np.count_nonzero(region_mask))
        same_event_found = int(
            np.count_nonzero(
                region_mask & same_event_high_exists
            )
        )
        found = int(
            np.count_nonzero(
                region_mask & pi0_mass_partner_found
            )
        )
        efficiency = (
            found / expected if expected > 0 else math.nan
        )

        rows.append(
            {
                "period": period.key,
                "period_label": period.label,
                "region": region_label,
                "detector_code": detector_code,
                "sector": sector_code,
                "expected_low_energy_tags": expected,
                "same_event_high_energy_partner": same_event_found,
                "matched_pi0_mass_partner": found,
                "efficiency": efficiency,
                "efficiency_error": simple_binomial_error(
                    found,
                    expected,
                ),
            }
        )
    # endfor

    classified_mask = (
        (low_detector == 0)
        | (
            (low_detector == 1)
            & (low_sector >= 1)
            & (low_sector <= 6)
        )
    )

    payload: Dict[str, object] = {
        "period": period.key,
        "period_label": period.label,
        "sample": "reconstructed_AAOGEN_epgamma_simple",
        "path": path,
        "tree_entries": total_entries,
        "entries_read": entries_read,
        "threshold_GeV": 2.0,
        "pi0_mass_window_GeV": [
            float(PI0_MASS_WINDOW[0]),
            float(PI0_MASS_WINDOW[1]),
        ],
        "event_identity": (
            "quantized electron+proton reconstructed kinematics using "
            f"mc_signature_decimals={args.mc_signature_decimals}"
        ),
        "finite_identity_entries": finite_identity_entries,
        "low_entries_total": low_entries_total,
        "high_entries_total": high_entries_total,
        "exactly_2GeV_entries_excluded": exactly_2GeV_entries,
        "low_entries_with_finite_reconstruction_kinematics": (
            finite_low_kinematics
        ),
        "low_opportunities_built": int(low_keys.size),
        "high_energy_entries_stored": int(high_keys.size),
        "low_opportunities_with_same_event_high_partner": int(
            np.count_nonzero(same_event_high_exists)
        ),
        "low_opportunities_with_pi0_mass_partner": int(
            np.count_nonzero(pi0_mass_partner_found)
        ),
        "classified_low_opportunities": int(
            np.count_nonzero(classified_mask)
        ),
        "unclassified_low_opportunities": int(
            np.count_nonzero(~classified_mask)
        ),
        "low_opportunities_with_reconstructed_probe_E_gt_2": int(
            np.count_nonzero(
                np.isfinite(low_probe_E)
                & (low_probe_E > 2.0)
            )
        ),
        "rows": rows,
        "diagnostic_histograms": (
            simple_finalize_diagnostic_histograms(
                diagnostic_histograms
            )
        ),
        "elapsed_seconds": elapsed_seconds(start_time),
    }

    return period.key, payload


def aggregate_simple_fd(
    rows: Sequence[Mapping[str, object]],
) -> Dict[str, float]:
    """Combine the six FD sectors by summing numerators and denominators."""
    selected = [
        row
        for row in rows
        if str(row["region"]).startswith("FD sector ")
    ]
    expected = int(
        sum(int(row["expected_low_energy_tags"]) for row in selected)
    )
    matched = int(
        sum(int(row["matched_pi0_mass_partner"]) for row in selected)
    )
    efficiency = (
        matched / expected if expected > 0 else math.nan
    )
    return {
        "expected": expected,
        "matched": matched,
        "efficiency": efficiency,
        "efficiency_error": simple_binomial_error(
            matched,
            expected,
        ),
    }


def simple_ft_row(
    rows: Sequence[Mapping[str, object]],
) -> Mapping[str, object]:
    return next(row for row in rows if row["region"] == "FT")



def simple_histogram_arrays(
    payload: Mapping[str, object],
    region: str,
    name: str,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    histogram = payload["diagnostic_histograms"][region][name]
    edges = np.asarray(histogram["edges"], dtype=float)
    counts = np.asarray(histogram["counts"], dtype=float)
    unit_area = np.asarray(histogram["unit_area"], dtype=float)
    total = int(histogram["in_range_entries"])
    return edges, counts, unit_area, total


def plot_simple_diagnostic_overlays(
    path: Path,
    period_label: str,
    data_payload: Mapping[str, object],
    mc_payload: Mapping[str, object],
    variables: Sequence[Tuple[str, str]],
    title_suffix: str,
) -> None:
    """
    Plot predicted FT on the top row and aggregate predicted FD on the bottom.
    Each sample is independently normalized to unit integral in the displayed
    range.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(
        2,
        3,
        figsize=(18.5, 9.0),
        squeeze=False,
    )

    for row_index, region in enumerate(("FT", "FD")):
        for column_index, (name, xlabel) in enumerate(variables):
            axis = axes[row_index, column_index]

            data_edges, _data_counts, data_unit, data_total = (
                simple_histogram_arrays(
                    data_payload,
                    region,
                    name,
                )
            )
            mc_edges, _mc_counts, mc_unit, mc_total = (
                simple_histogram_arrays(
                    mc_payload,
                    region,
                    name,
                )
            )

            axis.stairs(
                data_unit,
                data_edges,
                linewidth=1.7,
                color="black",
                label=f"Data ({data_total:,})",
                zorder=4,
            )
            axis.stairs(
                mc_unit,
                mc_edges,
                linewidth=1.5,
                color="tab:red",
                label=f"AAOGEN ({mc_total:,})",
                zorder=3,
            )

            if name == "predicted_theta_deg":
                for boundary in (2.5, 5.0, 35.0):
                    axis.axvline(
                        boundary,
                        linewidth=1.0,
                        linestyle="--",
                        color="0.45",
                    )
                # endfor
            # endif

            if name == "mgg_before":
                axis.axvspan(
                    PI0_MASS_WINDOW[0],
                    PI0_MASS_WINDOW[1],
                    alpha=0.14,
                    color="0.5",
                    label=r"$0.11<M_{\gamma\gamma}<0.16$ GeV",
                )
                axis.axvline(
                    0.1349768,
                    linewidth=1.0,
                    linestyle="--",
                    color="0.35",
                )
            # endif

            axis.set_xlabel(xlabel)
            axis.set_ylabel("fraction / bin")
            axis.set_ylim(bottom=0.0)
            axis.grid(axis="y", alpha=0.25)
            axis.set_title(
                f"Predicted {region}: {xlabel}"
            )

            if row_index == 0 and column_index == 0:
                axis.legend(frameon=False)
            # endif
            if name == "mgg_before":
                axis.legend(frameon=False, fontsize=8)
            # endif
        # endfor
    # endfor

    fig.suptitle(
        f"{period_label}: --simple data vs reconstructed AAOGEN\n"
        f"{title_suffix}",
        fontsize=15,
        y=0.99,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    fig.savefig(
        path,
        dpi=180,
        bbox_inches="tight",
    )
    plt.close(fig)


def simple_conditional_rows(
    periods: Sequence[PeriodConfig],
    data_metadata: Mapping[str, object],
    mc_metadata: Mapping[str, object],
) -> List[Dict[str, object]]:
    """
    Summarize the two factors entering the simple efficiency:
      P(high-E same-event partner | low tag)
      P(pi0 mass | same-event high-E partner)
    """
    rows: List[Dict[str, object]] = []

    for period in periods:
        for sample_name, metadata in (
            ("data", data_metadata),
            ("reconstructed_aaogen", mc_metadata),
        ):
            period_rows = metadata["periods"][period.key]["rows"]
            for region in ("FT", "FD"):
                if region == "FT":
                    selected = [
                        row for row in period_rows
                        if row["region"] == "FT"
                    ]
                else:
                    selected = [
                        row for row in period_rows
                        if str(row["region"]).startswith("FD sector ")
                    ]
                # endif

                denominator = int(
                    sum(
                        int(row["expected_low_energy_tags"])
                        for row in selected
                    )
                )
                same_event = int(
                    sum(
                        int(row["same_event_high_energy_partner"])
                        for row in selected
                    )
                )
                pi0_pass = int(
                    sum(
                        int(row["matched_pi0_mass_partner"])
                        for row in selected
                    )
                )

                rows.append(
                    {
                        "period": period.key,
                        "period_label": period.label,
                        "sample": sample_name,
                        "region": region,
                        "low_tag_denominator": denominator,
                        "same_event_high_partner": same_event,
                        "pi0_mass_partner": pi0_pass,
                        "p_high_given_low": (
                            same_event / denominator
                            if denominator > 0
                            else math.nan
                        ),
                        "p_pi0_given_high": (
                            pi0_pass / same_event
                            if same_event > 0
                            else math.nan
                        ),
                        "simple_efficiency": (
                            pi0_pass / denominator
                            if denominator > 0
                            else math.nan
                        ),
                    }
                )
            # endfor
        # endfor
    # endfor

    return rows


def write_simple_diagnostic_package(
    periods: Sequence[PeriodConfig],
    data_metadata: Mapping[str, object],
    mc_metadata: Mapping[str, object],
    simple_parent: Path,
) -> None:
    diagnostic_root = simple_parent / "diagnostics"
    diagnostic_root.mkdir(parents=True, exist_ok=True)

    first_canvas_variables = (
        (
            "low_energy",
            r"$E_{\gamma,\mathrm{low}}$ (GeV)",
        ),
        (
            "predicted_theta_deg",
            r"predicted $\theta_X$ (deg)",
        ),
        (
            "mgg_before",
            r"$M_{\gamma\gamma}$ before mass cut (GeV)",
        ),
    )
    second_canvas_variables = (
        (
            "low_theta_deg",
            r"$\theta_{\gamma,\mathrm{low}}$ (deg)",
        ),
        (
            "high_energy_before",
            r"$E_{\gamma,\mathrm{high}}$ before mass cut (GeV)",
        ),
        (
            "high_energy_after",
            r"$E_{\gamma,\mathrm{high}}$ after mass cut (GeV)",
        ),
    )

    for period in periods:
        data_payload = data_metadata["periods"][period.key]
        mc_payload = mc_metadata["periods"][period.key]

        plot_simple_diagnostic_overlays(
            diagnostic_root
            / f"{period.key}_ft_fd_core_diagnostics.png",
            period.label,
            data_payload,
            mc_payload,
            first_canvas_variables,
            (
                "top row: predicted FT; bottom row: predicted FD; "
                "each sample independently unit normalized"
            ),
        )
        plot_simple_diagnostic_overlays(
            diagnostic_root
            / f"{period.key}_ft_fd_photon_diagnostics.png",
            period.label,
            data_payload,
            mc_payload,
            second_canvas_variables,
            (
                "low-tag and high-partner kinematics; "
                "each sample independently unit normalized"
            ),
        )
    # endfor

    conditional_rows = simple_conditional_rows(
        periods,
        data_metadata,
        mc_metadata,
    )
    write_simple_rows_csv(
        simple_parent / "simple_conditional_probabilities.csv",
        conditional_rows,
    )

    diagnostic_manifest = {
        "description": (
            "Data versus reconstructed-AAOGEN diagnostics for --simple."
        ),
        "regions": ["predicted FT", "aggregate predicted FD"],
        "pi0_mass_window_GeV": [
            float(PI0_MASS_WINDOW[0]),
            float(PI0_MASS_WINDOW[1]),
        ],
        "plots": {
            "core": [
                "low_energy",
                "predicted_theta_deg",
                "mgg_before",
            ],
            "photon": [
                "low_theta_deg",
                "high_energy_before",
                "high_energy_after",
            ],
        },
        "conditional_probability_csv": str(
            simple_parent / "simple_conditional_probabilities.csv"
        ),
    }
    with open(
        simple_parent / "simple_diagnostic_manifest.json",
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(
            diagnostic_manifest,
            handle,
            indent=2,
            default=simple_json_default,
        )
    # endwith

def print_simple_data_mc_final_summary(
    periods: Sequence[PeriodConfig],
    data_metadata: Mapping[str, object],
    simple_mc_metadata: Mapping[str, object],
) -> None:
    """
    Short final --simple comparison using identical reconstructed algorithms.
    """
    log("")
    log("=" * 86)
    log("FINAL --simple DATA / RECONSTRUCTED-AAOGEN SUMMARY")
    log(
        "same reconstructed algorithm on both samples: "
        "E_low<2 GeV, E_high>2 GeV, 0.11<Mgg<0.16 GeV"
    )
    log(
        f"{'Period':<12s} "
        f"{'FT data/MC':>12s} "
        f"{'FD data/MC':>12s} "
        f"{'FT data':>10s} "
        f"{'FT MC':>10s} "
        f"{'FD data':>10s} "
        f"{'FD MC':>10s}"
    )

    data_periods = data_metadata["periods"]
    mc_periods = simple_mc_metadata["periods"]

    for period in periods:
        data_rows = data_periods[period.key]["rows"]
        mc_rows = mc_periods[period.key]["rows"]

        data_ft = simple_ft_row(data_rows)
        mc_ft = simple_ft_row(mc_rows)
        data_fd = aggregate_simple_fd(data_rows)
        mc_fd = aggregate_simple_fd(mc_rows)

        data_ft_eff = float(data_ft["efficiency"])
        mc_ft_eff = float(mc_ft["efficiency"])
        data_fd_eff = float(data_fd["efficiency"])
        mc_fd_eff = float(mc_fd["efficiency"])

        ft_ratio = (
            data_ft_eff / mc_ft_eff
            if np.isfinite(data_ft_eff)
            and np.isfinite(mc_ft_eff)
            and mc_ft_eff > 0.0
            else math.nan
        )
        fd_ratio = (
            data_fd_eff / mc_fd_eff
            if np.isfinite(data_fd_eff)
            and np.isfinite(mc_fd_eff)
            and mc_fd_eff > 0.0
            else math.nan
        )

        log(
            f"{period.label:<12s} "
            f"{ft_ratio:12.3f} "
            f"{fd_ratio:12.3f} "
            f"{100.0*data_ft_eff:9.2f}% "
            f"{100.0*mc_ft_eff:9.2f}% "
            f"{100.0*data_fd_eff:9.2f}% "
            f"{100.0*mc_fd_eff:9.2f}%"
        )
    # endfor

    data_combined_rows = data_metadata["combined_rows"]
    mc_combined_rows = simple_mc_metadata["combined_rows"]

    data_ft = simple_ft_row(data_combined_rows)
    mc_ft = simple_ft_row(mc_combined_rows)
    data_fd = aggregate_simple_fd(data_combined_rows)
    mc_fd = aggregate_simple_fd(mc_combined_rows)

    data_ft_eff = float(data_ft["efficiency"])
    mc_ft_eff = float(mc_ft["efficiency"])
    data_fd_eff = float(data_fd["efficiency"])
    mc_fd_eff = float(mc_fd["efficiency"])

    ft_ratio = data_ft_eff / mc_ft_eff
    fd_ratio = data_fd_eff / mc_fd_eff

    log("-" * 86)
    log(
        f"{'COMBINED':<12s} "
        f"{ft_ratio:12.3f} "
        f"{fd_ratio:12.3f} "
        f"{100.0*data_ft_eff:9.2f}% "
        f"{100.0*mc_ft_eff:9.2f}% "
        f"{100.0*data_fd_eff:9.2f}% "
        f"{100.0*mc_fd_eff:9.2f}%"
    )
    log("=" * 86)

def run_simple_data_efficiency_stage(
    periods: Sequence[PeriodConfig],
    args: argparse.Namespace,
    workers: int,
    data_root: Path,
) -> Dict[str, object]:
    """
    Run the temporary --simple data-side override and persist a compact audit.
    """
    simple_root = data_root / "simple"
    simple_root.mkdir(parents=True, exist_ok=True)

    worker_count = max(
        1,
        min(int(workers), len(periods), MAX_WORKERS),
    )
    log(
        f"SIMPLE DATA OVERRIDE: {len(periods)} period(s), "
        f"{worker_count} worker(s); normal data shapes/template fits skipped."
    )

    payloads: Dict[str, Dict[str, object]] = {}

    if worker_count == 1:
        for period in periods:
            key, payload = process_period_simple_data_efficiency(
                period,
                vars(args),
            )
            payloads[key] = payload
        # endfor
    else:
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=worker_count
        ) as executor:
            future_map = {
                executor.submit(
                    process_period_simple_data_efficiency,
                    period,
                    vars(args),
                ): period.key
                for period in periods
            }
            for future in concurrent.futures.as_completed(
                future_map
            ):
                key, payload = future.result()
                payloads[key] = payload
            # endfor
        # endwith
    # endif

    ordered_rows: List[Dict[str, object]] = []
    for period in periods:
        payload = payloads[period.key]
        rows = payload["rows"]
        ordered_rows.extend(rows)
        print_simple_data_efficiency_block(
            period.label,
            rows,
        )
        log(
            f"{period.label} simple audit: "
            f"low tags={payload['low_entries_total']:,}, "
            f"high entries={payload['high_entries_total']:,}, "
            f"classified low opportunities="
            f"{payload['classified_low_opportunities']:,}, "
            f"unclassified="
            f"{payload['unclassified_low_opportunities']:,}."
        )
    # endfor

    # Aggregate exactly by numerator/denominator, not by averaging efficiencies.
    combined_rows: List[Dict[str, object]] = []
    for region_label, detector_code, sector_code in [
        ("FT", 0, 0),
        *[
            (f"FD sector {sector}", 1, sector)
            for sector in range(1, 7)
        ],
    ]:
        selected = [
            row for row in ordered_rows
            if row["region"] == region_label
        ]
        expected = int(
            sum(
                int(row["expected_low_energy_tags"])
                for row in selected
            )
        )
        same_event = int(
            sum(
                int(row["same_event_high_energy_partner"])
                for row in selected
            )
        )
        matched = int(
            sum(
                int(row["matched_pi0_mass_partner"])
                for row in selected
            )
        )
        efficiency = matched / expected if expected > 0 else math.nan
        combined_rows.append(
            {
                "period": "combined",
                "period_label": "Combined selected periods",
                "region": region_label,
                "detector_code": detector_code,
                "sector": sector_code,
                "expected_low_energy_tags": expected,
                "same_event_high_energy_partner": same_event,
                "matched_pi0_mass_partner": matched,
                "efficiency": efficiency,
                "efficiency_error": simple_binomial_error(matched, expected),
            }
        )
    # endfor

    print_simple_data_efficiency_block(
        "COMBINED SELECTED PERIODS",
        combined_rows,
    )

    all_rows = [*ordered_rows, *combined_rows]
    write_simple_rows_csv(
        simple_root / "simple_data_efficiency.csv",
        all_rows,
    )

    metadata: Dict[str, object] = {
        "mode": "simple_data_override",
        "threshold_GeV": 2.0,
        "periods": payloads,
        "combined_rows": combined_rows,
        "output_csv": str(
            simple_root / "simple_data_efficiency.csv"
        ),
    }
    with open(
        simple_root / "simple_data_efficiency.json",
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(
            metadata,
            handle,
            indent=2,
            default=simple_json_default,
        )
    # endwith

    return metadata


def run_simple_mc_efficiency_stage(
    periods: Sequence[PeriodConfig],
    args: argparse.Namespace,
    workers: int,
    data_root: Path,
) -> Dict[str, object]:
    """
    Run the identical reconstructed --simple algorithm on AAOGEN epgamma MC.
    """
    simple_root = data_root / "simple"
    simple_root.mkdir(parents=True, exist_ok=True)

    worker_count = max(
        1,
        min(int(workers), len(periods), MAX_WORKERS),
    )
    log(
        f"SIMPLE RECONSTRUCTED-AAOGEN COMPARISON: "
        f"{len(periods)} period(s), {worker_count} worker(s)."
    )

    payloads: Dict[str, Dict[str, object]] = {}

    if worker_count == 1:
        for period in periods:
            key, payload = process_period_simple_mc_efficiency(
                period,
                vars(args),
            )
            payloads[key] = payload
        # endfor
    else:
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=worker_count
        ) as executor:
            future_map = {
                executor.submit(
                    process_period_simple_mc_efficiency,
                    period,
                    vars(args),
                ): period.key
                for period in periods
            }
            for future in concurrent.futures.as_completed(
                future_map
            ):
                key, payload = future.result()
                payloads[key] = payload
            # endfor
        # endwith
    # endif

    ordered_rows: List[Dict[str, object]] = []
    for period in periods:
        payload = payloads[period.key]
        rows = payload["rows"]
        ordered_rows.extend(rows)
        print_simple_data_efficiency_block(
            f"{period.label} reconstructed AAOGEN",
            rows,
        )
    # endfor

    combined_rows: List[Dict[str, object]] = []
    for region_label, detector_code, sector_code in [
        ("FT", 0, 0),
        *[
            (f"FD sector {sector}", 1, sector)
            for sector in range(1, 7)
        ],
    ]:
        selected = [
            row for row in ordered_rows
            if row["region"] == region_label
        ]
        expected = int(
            sum(
                int(row["expected_low_energy_tags"])
                for row in selected
            )
        )
        same_event = int(
            sum(
                int(row["same_event_high_energy_partner"])
                for row in selected
            )
        )
        matched = int(
            sum(
                int(row["matched_pi0_mass_partner"])
                for row in selected
            )
        )
        efficiency = (
            matched / expected if expected > 0 else math.nan
        )
        combined_rows.append(
            {
                "period": "combined",
                "period_label": "Combined selected periods",
                "region": region_label,
                "detector_code": detector_code,
                "sector": sector_code,
                "expected_low_energy_tags": expected,
                "same_event_high_energy_partner": same_event,
                "matched_pi0_mass_partner": matched,
                "efficiency": efficiency,
                "efficiency_error": simple_binomial_error(
                    matched,
                    expected,
                ),
            }
        )
    # endfor

    all_rows = [*ordered_rows, *combined_rows]
    write_simple_rows_csv(
        simple_root / "simple_reconstructed_aaogen_efficiency.csv",
        all_rows,
    )

    metadata: Dict[str, object] = {
        "mode": "simple_reconstructed_AAOGEN_comparison",
        "threshold_GeV": 2.0,
        "pi0_mass_window_GeV": [
            float(PI0_MASS_WINDOW[0]),
            float(PI0_MASS_WINDOW[1]),
        ],
        "event_identity": (
            "quantized reconstructed electron+proton kinematics; "
            "no truth matching"
        ),
        "periods": payloads,
        "combined_rows": combined_rows,
        "output_csv": str(
            simple_root
            / "simple_reconstructed_aaogen_efficiency.csv"
        ),
    }
    with open(
        simple_root
        / "simple_reconstructed_aaogen_efficiency.json",
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(
            metadata,
            handle,
            indent=2,
            default=simple_json_default,
        )
    # endwith

    return metadata


def preflight_enabled_stages(
    periods: Sequence[PeriodConfig],
    args: argparse.Namespace,
) -> List[Dict[str, object]]:
    manifest: List[Dict[str, object]] = []

    for period in periods:
        period_entry: Dict[str, object] = {
            "key": period.key,
            "label": period.label,
            "beam_energy_GeV": period.beam_energy_GeV,
        }

        if args.simple:
            roles = ["epg_data", "pi0_epg_mc"]
        else:
            roles = ["epg_data", "dvcs_mc", "pi0_epg_mc"]
        # endif

        if not args.skip_mc_efficiency:
            roles.append("pi0_epgg_mc")
        # endif

        for role in roles:
            path = getattr(period, role)
            entries, branches = require_tree(path)
            period_entry[role] = path
            period_entry[f"{role}_entries"] = entries
            period_entry[f"{role}_branches"] = branches
            log(
                f"Preflight {period.label} {role}: {entries:,} entries"
            )
        # endfor

        manifest.append(period_entry)
    # endfor

    return manifest



def plot_all_period_integrated(
    path: Path,
    rows: Sequence[Mapping[str, object]],
) -> None:
    period_labels = [period.label for period in PERIODS]
    x = np.arange(len(period_labels), dtype=float)
    offsets = {"FT": -0.08, "FD": 0.08}

    fig, axis = plt.subplots(figsize=(12, 7))
    for topology in ("FT", "FD"):
        selected = []
        for label in period_labels:
            match = next(
                (
                    row for row in rows
                    if row["period_label"] == label
                    and row["topology"] == topology
                ),
                None,
            )
            selected.append(match)
        # endfor

        efficiencies = np.asarray(
            [
                row["efficiency"] if row is not None else math.nan
                for row in selected
            ],
            dtype=float,
        )
        errors = np.asarray(
            [
                row["efficiency_error"] if row is not None else math.nan
                for row in selected
            ],
            dtype=float,
        )
        axis.errorbar(
            x + offsets[topology],
            efficiencies,
            yerr=errors,
            fmt="o",
            capsize=3,
            label=topology,
        )
    # endfor

    axis.set_xticks(x, period_labels, rotation=20, ha="right")
    axis.set_ylabel(r"$\epsilon_{\mathrm{MC}}$")
    axis.set_ylim(0.0, 1.05)
    axis.grid(alpha=0.25)
    axis.legend()
    axis.set_title(
        "AAOGEN high-energy photon reconstruction efficiency"
    )
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    periods = selected_periods(args)
    output_root = Path(args.output_dir)
    data_root = output_root / "data"
    mc_root = output_root / "mc"
    simple_parent = output_root / "simple"

    output_root.mkdir(parents=True, exist_ok=True)
    data_root.mkdir(parents=True, exist_ok=True)
    mc_root.mkdir(parents=True, exist_ok=True)
    if args.simple:
        simple_parent.mkdir(parents=True, exist_ok=True)
    # endif

    workers = max(
        1,
        min(int(args.workers), MAX_WORKERS),
    )

    fitter_preflight: Dict[str, object] = {
        "enabled": (
            not args.simple
            and not args.skip_data_template_fits
        ),
        "implementation": "internal advanced component-preserving fitter",
        "scipy_optimizer": "L-BFGS-B plus bounded scalar coordinate update",
        "variables": {
            category_key: {
                variable.key: variable_model_spec(
                    variable.key,
                    category_key,
                )
                for variable in fit_variables_for_category(
                    category_key
                )
            }
            for category_key in (
                "one_photon",
                "two_photon_more_energetic",
                "two_photon_less_energetic",
            )
        },
        "shared_fraction": True,
        "independent_variable_fits": True,
    }

    log(
        "MC efficiency stage: "
        + (
            "DISABLED by --skip-mc-efficiency"
            if args.skip_mc_efficiency
            else "ENABLED (default)"
        )
    )

    manifest = {
        "script": Path(__file__).name,
        "mode": (
            "Fast data 2x5 shape comparison and photon multiplicity "
            "diagnostics with the complete AAOGEN MC efficiency enabled by default"
        ),
        "created_unix_time": time.time(),
        "arguments": vars(args),
        "resolved_stages": {
            "simple_data_override": bool(args.simple),
            "data_shape_and_multiplicity": not args.simple,
            "data_template_fits": (
                not args.simple
                and not args.skip_data_template_fits
            ),
            "mc_efficiency": not args.skip_mc_efficiency,
        },
        "fitter_preflight": fitter_preflight,
        "periods": preflight_enabled_stages(periods, args),
    }
    with open(
        output_root / "input_manifest.json",
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(manifest, handle, indent=2)
    # endwith

    if args.simple:
        data_metadata = run_simple_data_efficiency_stage(
            periods,
            args,
            workers,
            data_root,
        )
        simple_mc_metadata = run_simple_mc_efficiency_stage(
            periods,
            args,
            workers,
            data_root,
        )
    else:
        data_metadata = run_data_shape_stage(
            periods,
            args,
            workers,
        )
        simple_mc_metadata = {}
    # endif

    mc_rows: List[Dict[str, object]] = []
    mc_metadata: Dict[str, object] = {}

    if not args.skip_mc_efficiency:
        mc_workers = max(
            1,
            min(workers, len(periods)),
        )
        log(
            f"AAOGEN MC EFFICIENCY STAGE: "
            f"{len(periods)} period(s), {mc_workers} worker(s)."
        )

        if mc_workers == 1:
            for period in periods:
                key, rows, payload = process_period_mc_only(
                    period,
                    vars(args),
                )
                mc_metadata[key] = payload
                mc_rows.extend(rows)
            # endfor
        else:
            with concurrent.futures.ProcessPoolExecutor(
                max_workers=mc_workers
            ) as executor:
                future_map = {
                    executor.submit(
                        process_period_mc_only,
                        period,
                        vars(args),
                    ): period.key
                    for period in periods
                }
                for future in concurrent.futures.as_completed(future_map):
                    key, rows, payload = future.result()
                    mc_metadata[key] = payload
                    mc_rows.extend(rows)
                # endfor
            # endwith
        # endif

        period_order = {
            period.key: index
            for index, period in enumerate(PERIODS)
        }
        topology_order = {"FT": 0, "FD": 1}
        mc_rows.sort(
            key=lambda row: (
                period_order.get(str(row["period"]), 999),
                topology_order.get(str(row["topology"]), 999),
            )
        )

        write_dict_rows(
            mc_root / "photon_efficiency_mc.csv",
            mc_rows,
        )
        with open(
            mc_root / "photon_efficiency_mc.json",
            "w",
            encoding="utf-8",
        ) as handle:
            json.dump(
                {
                    "mode": "AAOGEN MC only",
                    "arguments": vars(args),
                    "rows": mc_rows,
                    "period_metadata": mc_metadata,
                },
                handle,
                indent=2,
            )
        # endwith

        plot_all_period_integrated(
            mc_root / "all_periods_integrated_mc_efficiency.png",
            mc_rows,
        )
    # endif

    with open(
        output_root / "stepwise_study_summary.json",
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(
            {
                "arguments": vars(args),
                "data_shape_and_multiplicity_stage": data_metadata,
                "simple_reconstructed_aaogen_stage": simple_mc_metadata,
                "mc_efficiency_rows": mc_rows,
                "mc_efficiency_metadata": mc_metadata,
            },
            handle,
            indent=2,
            default=simple_json_default,
        )
    # endwith

    if args.simple:
        write_simple_diagnostic_package(
            periods,
            data_metadata,
            simple_mc_metadata,
            simple_parent,
        )
        print_simple_data_mc_final_summary(
            periods,
            data_metadata,
            simple_mc_metadata,
        )
        log(
            "Wrote --simple FT/FD diagnostic plots to "
            f"{simple_parent / 'diagnostics'}"
        )
    # endif

    log(f"Wrote photon-efficiency study outputs to {output_root}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"FATAL ERROR: {exc}", file=sys.stderr, flush=True)
        raise
    # endtry
