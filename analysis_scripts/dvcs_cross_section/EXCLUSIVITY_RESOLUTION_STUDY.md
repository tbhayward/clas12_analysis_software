# DATA/MC exclusivity shape study

This is the first diagnostic before attempting particle-level MC smearing. It
compares conditional exclusivity-distribution widths, locations and tails. It
does not alter production events, cuts, acceptances, correction factors or output.

## Run

From `analysis_scripts/dvcs_cross_section`, with ROOT and Python numpy/matplotlib
available, use a completed, preserved nominal production run. That run must
contain `global_cuts.cpp`, `global_cuts.h`, `kinematic_bin_integration.h`, and
`output/jsons/{global_cuts_config,combined_cuts}.json`. The backend compiles
automatically against those saved sources; no `make` is required.

A small first check (one sorted ROOT part per period and sample):

```sh
python3 study_exclusivity_resolution.py --production-run "$PWD/test_runs/current_mc_fixed_001/proton_t_default_001" --periods fa18_inb fa18_out --workers 10 --max-files-per-sample 1 --output "$PWD/test_runs/exclusivity_resolution_fa18_smoke_001"
```

Full two-polarity study:

```sh
python3 study_exclusivity_resolution.py --production-run "$PWD/test_runs/current_mc_fixed_001/proton_t_default_001" --periods fa18_inb fa18_out --workers 10 --output "$PWD/test_runs/exclusivity_resolution_fa18_001"
```

`--input-root` defaults to `/volatile/clas12/thayward/dvcs_pass2_full_rebuild`.
All five production periods are supported; Fa18 inb/out are the defaults. Data
uses `data/PERIOD/dvcs_data/parts/*.root`; MC uses
`mc/PERIOD/dvcs_mc/CURRENT/parts/*.rec_only.root` (45 nA for Sp18 out,
50 nA otherwise). Missing compact MC inputs fail clearly rather than silently
reading the much larger generated population.

## What to open

Open `REPORT.html` first. It is self-contained (no internet/CDN), with selectors
for period, topology, variable, and momentum/angle interval. Choose **Inclusive**
first. Then examine proton momentum and polar angle, followed by electron and
photon subdivisions. The log-axis checkbox helps inspect tails. The full-range
checkbox reveals the complete histogram range; the initial zoom is only visual.

Each period/topology also has `inclusive_shapes.png` and `shape_summary.png`.
The latter shows width ratios, median shifts and excess tails across all six
kinematic axes; gray cells are flagged, not zero. Color scales saturate at their
labeled limits. `summary.csv` contains every numerical result and flag.
`REPORT.md` explains definitions and gives counts after all cuts.

The first question is where discrepancies consistently appear across different
variables. A large width ratio does not prove which particle is responsible.
These observations should inform a subsequent controlled smearing test, not be
converted directly into cross-section corrections.

## Selections and interpretation

The worker calls the saved production phase-space, proton-momentum and full
auxiliary/sector-aware global-cut functions, including the saved run blacklist.
It uses the saved period/topology-specific DATA or MC exclusivity windows.
In each projection it omits only the plotted variable's cut. If that cut is
absent, it retains all selected cuts. When plotting pTmiss, it also disables its
optional global pTmiss pre-cut; otherwise that pre-cut would truncate the plot.
All other global cuts, including physical t_min, remain enabled as configured.

The seven variables are Delta_phi, theta, theta_gamma_gamma, pTmiss, Emiss2,
Mx2 and Mx2_2. Missing required branches fail; Delta_phi alone has the same
absolute wrapped p1_phi−p2_phi fallback used in production. Branch units follow
the producer: theta_gamma_gamma is **degrees**, Emiss2 is **missing energy in
GeV**, theta is radians, and particle polar-angle slices convert radians to degrees.

Core spread is `(q84−q16)/2`; it equals approximately one sigma for a Gaussian,
but no Gaussian assumption or fit is made. The median shift is expressed in MC
core spreads. The optional reported histogram mode uses [1,2,1] smoothing;
it is descriptive, not a fitted peak position. Tail fractions use the MC central
95% interval and include finite underflow/overflow. Quantiles interpolate
within histogram bins using in-range counts. Results with less than 99% finite
histogram coverage, fewer than `--min-entries` (default 300) candidates, nonfinite
values, or a core spread smaller than two histogram bins are flagged. No
statistical significance is inferred from the width ratios.

All counts are unit weighted reconstructed DVCS candidates. No pi0 subtraction,
current correction, Neupane weighting or generator reweighting is applied.
Selection-dependent backgrounds, radiation, physics populations and different
other-variable cut windows may change these shapes. Six separate 1D projections
do not match the full joint kinematic population. Global phase space is used,
not a restriction to populated production CSV bins. A direct calibration would
need follow-up population/background controls; this report does not supply one.

## Cost, reproducibility and reuse

A compiled C++ worker streams active scalar branches once per file and fills
all variables and subdivisions together. Only fixed-size histograms reach
Python; no event arrays or Python event loop. A single pool of at most ten
workers covers all data, MC and periods, using batches balanced by file size.
Workers do not spawn ROOT thread pools. On failure, queued work is cancelled and
active workers are terminated; incomplete outputs lack a completed manifest.

The output directory must be new. Existing production and study outputs are
never overwritten. Each study has its own backend and batch logs, so independent
studies can run simultaneously without sharing mutable files. The manifest
records saved-cut/source hashes, input file paths/sizes/mtimes, entries, branch
fallbacks and histogram checksum. Changes to ROOT inputs during the scan fail.

To replot with a different low-statistics threshold, without ROOT reads:

```sh
python3 study_exclusivity_resolution.py --reuse-histograms "$PWD/test_runs/exclusivity_resolution_fa18_001" --min-entries 500 --output "$PWD/test_runs/exclusivity_resolution_fa18_replot_001"
```

Changing histogram ranges or particle-slice edges requires a new event scan.
Use `--config custom.json`, for example:

```json
{"variables":{"Mx2":{"bins":800,"low":-0.2,"high":0.2}},"axes":{"p1_theta":[0,35,45,50,55,60,65,180.000001]}}
```

Reuse validates the histogram checksum and uses saved selections/strata. It does
not need original ROOT files or the production directory. Standard tests:

```sh
python3 test_exclusivity_resolution.py
```

On a working ROOT installation, also run the tiny native end-to-end fixture:

```sh
python3 test_exclusivity_resolution.py --root-smoke --production-run "$PWD/test_runs/current_mc_fixed_001/proton_t_default_001"
```

That test writes tiny synthetic PhysicsEvents trees and checks both samples,
both periods, all topologies, native production global cuts, leave-one-out
behavior, optional global pT cut omission, blacklist, overflow and phi fallback.
It uses temporary files and does not read farm production ROOT data.

## Background-aware update (default)

The default now also reads two-photon pi0 data, nominal-current AAOgen two-photon
MC and pi0-as-DVCS leakage MC. Controls are pooled across topologies per period;
losing a photon can change topology. The ratio control DATA/control MC scales
the leakage sample. The residual fully selected DVCS-candidate yield fixes the
signal scale separately by topology. These scales remain fixed in every slice
and omitted-cut distribution; no fit to the displayed variable is performed.
Signal and leakage templates use the DATA exclusivity windows for a common
one-photon selection. Two-photon controls use their saved eppi0 DATA/MC windows.

This is a raw-candidate mixture diagnostic; it does not import current-corrected
yields or apply cross-section weights to raw histograms. A single period-wide
AAOgen scale is a first approximation, not a local generator reweighting.
`--aaogen-exposure-ratio` is generated exposure(control MC)/exposure(leakage MC),
default 1 following production's common-exposure convention. Different generated
exposures require the corresponding ratio; a filename filter cannot establish
that equality. Complete samples are required; file limits are supported only
with `--signal-only`. Controls must be signal-dominated for this normalization.

`normalizations.csv` records control counts, background scales, their counting-
only errors, selected background yields and signal scales. Invalid residual
signal yields stop the report rather than clipping the background. MC in the
main summary/heatmaps now means the SUM of signal and background; the original
signal-only comparisons are saved separately in `signal_only_summary.csv`.
The HTML and inclusive PNG overlays show DATA, signal, background and total.
All components share the DATA denominator so predicted integrals are preserved.

Avoid rereading the existing DATA events:

```sh
python3 study_exclusivity_resolution.py --production-run "$PWD/test_runs/current_mc_fixed_001/proton_t_default_001" --reuse-data-from "$PWD/test_runs/exclusivity_resolution_fa18_001" --periods fa18_inb fa18_out --workers 10 --output "$PWD/test_runs/exclusivity_resolution_fa18_pi0_001"
```

The seed DATA histograms are used only after checking cuts, source hashes,
strata, ranges, checksum and exact input-file inventory. Signal MC is reread
because its selection changes to DATA windows. Existing results are preserved.
Use `--signal-only` for the earlier no-background diagnostic. Version-1 studies
can seed DATA, but cannot be used as complete version-2 histogram replays.

### Recover an empty two-photon MC control (composite detector flag)

Direct pi0 MC uses `detector_gamma1` for topology assignment, matching
`total_counts.cpp`; its composite `detector2` can be -1. The ROOT smoke fixture
explicitly tests this representation. To recover a completed scan made before
this correction, use `--repair-controls-from OLD_STUDY` with a new `--output`.
This verifies/reuses DATA, DVCS MC and leakage histograms, rereading only the
two-photon DATA/MC controls. Do not use complete `--reuse-histograms` for that
repair: it would retain the original zero control counts.

## All-period joint kinematics (v3)

All five periods are now the default: Fa18 in/out, Sp18 in/out and Sp19 in.
The study still fills every topology, with one shared ten-worker budget.

Two new menus provide **joint** cells rather than independent projections:

- xB × Q² × proton momentum × proton polar angle (24 cells).
- xB × Q² × photon energy × photon polar angle (24 cells).

The initial bins are deliberately broad enough to retain statistics. Each
family partitions its supported domain without overlap. An event fills at most
one cell per joint family; these are separate comparisons, never summed across
families. xB splits at 0.2, Q² at 2.5 GeV²; proton momentum splits at 0.5/0.8
GeV and angle at 50 degrees; photon energy at 6/8 GeV and angle at 8 degrees.
Joint cells reduce physics/detector population mixing but do not match every
variable or the within-cell distributions. No weights are fitted to the
exclusivity variables; signal/background scales remain fixed by period/topology.
No particle smearing or new production correction is applied at this stage.

Missing mass squared epγ now covers [-1,+1] GeV², retaining the original
0.0004 GeV² bins over [-0.12,+0.12] and using wider bins in the tails.
Plots show density per variable unit so wider tail bins do not create false
shape features. Histogram quantiles integrate counts, not density. The mode is
estimated from smoothed density. The bin-resolution flag uses the local bin
width near the median. Underflow/overflow and range flags remain explicit.
CSV summaries separately report median/mode, core spread/asymmetry, and lower
and upper tail fractions outside the MC central 95% interval. These descriptive
quantities do not identify a particle resolution or supply statistical fits.

REPORT.html is now a small index linking separate period/topology interactive
reports. Zero histogram bins are stored sparsely in those pages to keep browser
memory and parsing costs down. joint_summary.csv contains the joint-cell
metrics; period_comparison.csv places each period's metrics/status alongside
one another for identical cells. Existing inclusive PNG/heatmaps remain compact;
the joint cells are accessed in the interactive menus or CSVs.

ROOT files are still streamed once each. Native mixed-radix lookup fills just
one cell per joint grid, instead of scanning every joint cell for each event.
Two-photon controls store counts only rather than empty diagnostic histograms.
This keeps the additional event-loop cost and memory modest.

A fresh scan is necessary: older aggregates cannot recover the new tails or
joint event correlations. Version-2 histogram replotting remains supported,
but seed/recovery options correctly reject incompatible new ranges/strata.

```sh
python3 test_exclusivity_resolution.py --root-smoke --production-run "$PWD/test_runs/current_mc_fixed_001/proton_t_default_001" && python3 study_exclusivity_resolution.py --production-run "$PWD/test_runs/current_mc_fixed_001/proton_t_default_001" --periods sp18_inb sp18_out fa18_inb fa18_out sp19_inb --workers 10 --output "$PWD/test_runs/exclusivity_resolution_all_periods_joint_001"
```

Read CD_FT first for response diagnostics, then the FD-photon topologies for
background sensitivity. Compare the same joint cells across periods, respecting
low-statistics/range flags. Choose which particle-level smearing variations to
test only after checking these more controlled shapes. The default AAOgen
common-exposure assumption and other background-transfer limitations still apply.
