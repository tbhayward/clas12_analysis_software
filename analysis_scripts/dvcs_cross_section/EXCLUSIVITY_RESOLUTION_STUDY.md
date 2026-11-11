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
