# Preliminary exclusivity selection

This analysis follows Maggie's selection sequence, with independently fitted cuts in every one of the 24 (xB, −t′) bins and each of the six periods. The revised workflow uses a linear background as the baseline, tests a quadratic that cannot curve upward, checks histogram sensitivity, and displays angular residuals in degrees.

## Run on one farm node

From the repository root:

```bash
python3 analysis_scripts/RGC_eppi-/channel_selection/exclusivity/run_selection.py --workers 6
```

All figures, fit tables, candidate caches and progress files are saved beside the scripts in `analysis_scripts/RGC_eppi-/channel_selection/exclusivity/output/`. Input ROOT files remain in the existing campaign on `/volatile`. The default destination is anchored to the script location, independent of your working directory. Six workers each handle one complete period. Avoid simultaneous launches into the same output directory.

Progress from another terminal:

```bash
python3 analysis_scripts/RGC_eppi-/channel_selection/exclusivity/run_selection.py --status
```

Dependencies are Python 3, NumPy, SciPy, uproot 5 and Matplotlib. ROOT itself is not required.

```bash
python3 -m pip install numpy scipy 'uproot>=5,<6' matplotlib
```

A limited smoke test uses a separate directory:

```bash
python3 analysis_scripts/RGC_eppi-/channel_selection/exclusivity/run_selection.py --workers 6 --max-entries 200000 --output analysis_scripts/RGC_eppi-/channel_selection/exclusivity/output/smoke_test
```

This reads the first stored candidates per period. It may not populate every bin and is marked PARTIAL SAMPLE.

## Procedure

1. Fit the deuteron-reference missing mass squared `(beam + stationary deuteron − electron − proton − pion)²` near the spectator-proton peak (mp² ≈ 0.88 GeV²). The nominal model is Gaussian plus linear background. Select fitted mean ± twice width.
2. Histogram the **full three-dimensional** missing momentum `|beam − electron − proton − pion|` after the nominal mass cut. Fit Gaussian plus linear background and select mean ± twice width. Clip the lower boundary at zero. This is distinct from transverse missing momentum.
3. Plot **theta_measured − theta_predicted in degrees** after both nominal cuts. The pion prediction uses a stationary neutron. A Gaussian is only a visual guide. No angular cut is imposed.

The default mass histogram has 0.010 GeV² bins and the momentum histogram has 0.005 GeV bins, twice their original widths. This retains all 24 analysis bins. Angular plots cover −6 to +6 degrees in 50 histogram bins. Angular guide means and widths in JSON and plot annotations are also degrees.

Plots contain one black data trace, the red nominal fit, a gray dotted background and green cut boundaries. Alternative fit curves are kept out of the primary figures. FD/CD counts after both cuts are recorded in the table. This pass shares cut windows between FD and CD protons within each bin and period.

## Background and stability checks

For each stage and bin, `stability.csv` records:

- Nominal Gaussian plus linear background.
- Gaussian plus a **nonnegative, concave-down quadratic** over the same fitted interval.
- The linear model fitted with another factor of two coarser histogram bins.

Both backgrounds are nonnegative throughout their fit interval. The quadratic is expressed as `b0*(1-u)^2 + 2*b1*u*(1-u) + b2*u^2`, with `b0,b2 >= 0` and `b1 = (b0+b2)/2 + c`, `c >= 0`. Its second derivative is nonpositive, so it cannot form an upward-curving valley under the signal. The linear baseline has zero curvature. The fits use the bin-integrated model, avoiding a bin-center approximation when comparing histogram widths.

Comparison rows include signed shifts of both cut boundaries. `cuts.csv` reports the maximum shift and flags a nominal usable cut if a comparison shifts a boundary by more than the larger of two nominal histogram bins or 10% of the nominal window width. This threshold is a review criterion, not a systematic uncertainty. If a comparison cannot resolve the peak, that gets a separate warning. Stability variations are conditional on the nominal upstream cut. They do not represent complete alternative selections or propagated downstream systematics.

### Fit-range checks

The initial mass fit range is 0.70–1.30 GeV². When the proposed cut goes outside it, the range is expanded with two histogram bins of margin and the fit is repeated. The mass range may expand as far as 0.40–1.50 GeV². The initial momentum range is 0–0.25 GeV and may expand to 0–0.45 GeV. Expansion stops once the cut is covered, a limit is reached, the fit becomes unresolved, or four fits have been attempted. Every interval actually fitted is saved in `range_history`. Remaining extrapolation is flagged. Plots expand to show the fitted interval and usable cut boundaries.

The mass-peak mean remains bounded to 0.78–1.05 GeV², with width 0.01–0.20 GeV², to avoid fitting a different large peak as the spectator. Momentum means are bounded to 0.005–0.18 GeV with widths 0.005–0.15 GeV. These are explicit starting assumptions for review.

Fits minimize binned Poisson deviance with three peak starting guesses. Fits with fewer than 100 entries, failed convergence, a mean/width at a bound, or less than 9 improvement in deviance over the same background-only model receive no usable cut. This is a practical peak-resolution screen, not a calibrated significance or a declaration that the bin contains no signal. Those bins do not contribute to downstream selected distributions and are labeled UNRESOLVED FIT. Empty panels after an unresolved upstream fit are not measured zero yields.

Pearson chi-square/ndf > 5 is flagged separately from background sensitivity, histogram sensitivity, unresolved comparisons and remaining range extrapolation. Warning-bearing nominal cuts are applied provisionally so their behavior can be inspected. Fitted statistical errors are local model-conditional approximations. Neither a small chi-square nor a narrow angular residual establishes purity. The nonnegative missing-momentum magnitude need not have a Gaussian signal shape. Nuclear backgrounds in ND3 can remain and still require the dilution study.

## Speed and cache behavior

The initial revised run needs **one ROOT scan per period**, replacing the earlier three scans. It selects the fiducial/kinematic sample, calculates only the three required observables and writes memory-mapped candidate chunks in `derived_cache/`. Each candidate retains double-precision mass, momentum and angular residual plus one-byte bin and detector labels: 26 bytes per combination (about 450 MiB for the six current samples combined). Only one source chunk is held in memory at a time.

The three selection stages then read this small cache. Changing the nominal mass cuts correctly rebuilds the downstream momentum and angular distributions from candidates, rather than reusing histograms selected with obsolete cuts. Later changes to fits, ranges or histogram widths can reuse the same candidate cache without ROOT decompression or new four-vector calculations. Input identity, kinematic/selection definitions and entry limit are checked independently of fit-model changes. Cache chunks are published for reuse only after the full source scan succeeds.

A normal rerun reuses complete matching results. To refit explicitly using a matching candidate cache:

```bash
python3 analysis_scripts/RGC_eppi-/channel_selection/exclusivity/run_selection.py --workers 6 --force
```

To force a new ROOT scan:

```bash
python3 analysis_scripts/RGC_eppi-/channel_selection/exclusivity/run_selection.py --workers 6 --rebuild-cache
```

An interrupted ROOT-cache build restarts that period's ROOT scan. Once that cache is complete, an interrupted fit run can restart from it. Old pre-revision histogram caches cannot supply changed downstream selections, so the first revised run builds the new candidate cache. No source data or final ROOT skim is modified.

## Inputs and accounting

The six default files are specified in `../../imports/datasets.json`. Corrected lab momenta and per-candidate beam energies are used, with p1 the proton and p2 the pion. Corrections are not applied again. Neutron-reference Mx2 is independently checked against the saved branch within 1e-6 GeV². The compact observable calculations are tested against the original shared diagnostic kinematics.

Default requirements remain finite required inputs, fiducial_status == 3, FD pion, FD or CD proton, Q² ≥ 1 GeV², W ≥ 2 GeV and y ≤ 0.80. The pi+ bin edges are unchanged: xB [0.10, 0.25, 0.35, 0.45, 0.60] and −t′ [0.05, 0.25, 0.45, 0.65, 0.85, 1.05, 1.25] GeV². Lower edges are included and upper edges excluded. Saved signed tprime is negated. Bin numbers are xB-major. Histograms count combinations, not unique events, with no charge or polarization weights.

## Outputs

Each period produces 12 PNG pages, each with six −t′ panels for one xB interval:

- `mass_xb1.png` through `mass_xb4.png`: deuteron missing mass, nominal fits and cuts.
- `momentum_xb1.png` through `momentum_xb4.png`: full missing momentum after the nominal mass cut.
- `angle_xb1.png` through `angle_xb4.png`: signed polar-angle residual in degrees after both nominal cuts.

There are 72 PNG pages for the six periods, with no PDF copies. `index.html` links the figures and tables. `cuts.csv` contains nominal parameters, cuts, quality and stability flags, and per-bin cut flows. `stability.csv` contains the comparison fits and boundary shifts. The fit JSON files preserve comparisons, ranges and errors. `angle_guides.json` explicitly declares degrees. `manifest.json` records schema version 2, input identity, cache provenance, full/partial status and per-bin counts including tails and undefined angular values. `histograms.npz` stores the final stage histograms. `progress.json` identifies ROOT-cache construction, cached histogram stages, fits or plotting.

Redraw figures without rerunning the fits:

```bash
python3 analysis_scripts/RGC_eppi-/channel_selection/exclusivity/run_selection.py --plots-only
```

The revised plotter expects schema-version-2 results. First rerun the analysis to upgrade old output. Use `--input-root PATH`, `--output PATH` and `--periods DATASET ...` to override defaults.

## Validation

```bash
python3 -m unittest discover -s analysis_scripts/RGC_eppi-/channel_selection/exclusivity -v
```

Tests check known-peak recovery, positive background curvature constraints, bin-integrated fits, adaptive range coverage, rebin accounting, degree-valued angular guides, cut boundaries, agreement with the original four-vector calculations and reuse of the derived cache without a second ROOT open. End-to-end synthetic runs check six-worker scanning, fitting, cut flow, rendering and cache-driven refits. Production cuts remain provisional until reviewed.
