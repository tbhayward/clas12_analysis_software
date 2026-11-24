# Preliminary exclusivity selection

This analysis follows the **sequence** in Maggie's RG-B preliminary note, Section 2.5. Cut parameters are fitted independently in each of the 24 (xB, −t′) bins and each of the six periods. It does not copy her period-integrated numerical cut limits.

## Run

From the analysis repository root, on one farm node:

```bash
python3 analysis_scripts/RGC_eppi-/channel_selection/exclusivity/run_selection.py --workers 6
```

Each of six workers handles one complete period, including its three scans, fits and plots. All outputs are saved under `analysis_scripts/RGC_eppi-/channel_selection/exclusivity/output/` beside this analysis. The output destination is anchored to the script location, independent of the working directory. The large input ROOT files are read from the existing campaign in `/volatile`. No new ROOT skim is written and the original diagnostics are preserved.

Check progress from another terminal:

```bash
python3 analysis_scripts/RGC_eppi-/channel_selection/exclusivity/run_selection.py --status
```

Status identifies the mass, momentum and angle passes separately. Each pass reports its own exact percentage of input rows read. Completed periods are cached and reused on a rerun. An interrupted period restarts its three passes. Run only one instance per output directory. A changed input, entry limit or analysis source invalidates its period cache.

Dependencies: Python 3, NumPy, SciPy, uproot 5 and Matplotlib. ROOT itself is not required.

```bash
python3 -m pip install numpy scipy 'uproot>=5,<6' matplotlib
```

An optional smoke test reads the first 200,000 candidate rows per period. It is deliberately separate from the full output and may lack enough events to fit many bins:

```bash
python3 analysis_scripts/RGC_eppi-/channel_selection/exclusivity/run_selection.py --workers 6 --max-entries 200000 --output analysis_scripts/RGC_eppi-/channel_selection/exclusivity/output/smoke_test
```

## Procedure

1. Construct the deuteron-reference missing mass squared, `(beam + stationary deuteron − electron − proton − pion)²`. Fit the spectator-proton peak locally with a Gaussian plus quadratic polynomial. Its expectation is near mp² ≈ 0.88 GeV². Use the fitted mean ± twice the fitted width as the provisional selection.
2. Rescan the file with the bin's mass selection. Histogram the **full three-dimensional** missing-momentum magnitude, `|beam − electron − proton − pion|`. Fit it with a Gaussian plus quadratic polynomial and select mean ± twice width. The lower limit is clipped at zero, as the magnitude is nonnegative. This is not the transverse missing momentum used in the initial diagnostic plots. Stationary neutron and stationary deuteron references give identical missing spatial momentum.
3. Rescan with both selections and show the signed pion polar-angle residual, `theta_measured − theta_predicted`, **in radians** to match the note. The predicted pion four-vector uses a stationary neutron. A Gaussian is drawn as a visual guide. No angular selection is applied and the guide's width is not used as a cut.

The primary plots contain one black data trace, a red fitted curve, a subtle dotted polynomial background and green cut boundaries. FD/CD overlays and correlation matrices are omitted. The numerical table retains the FD and CD counts after both cuts for detector comparisons. Cuts in this pass are shared by FD and CD protons within a bin and period. Detector-specific fits can be studied if the results warrant them.

### Fit choices and flags

The note specifies Gaussian plus polynomial without giving a polynomial order in its text. This implementation uses a quadratic background with nonnegative Bernstein coefficients. The model is positive over the fitted interval. Fits minimize the binned Poisson deviance and use three initial peak guesses. The table also reports Pearson chi-square divided by the number of degrees of freedom, which is a diagnostic rather than the minimized objective.

Mass plots show 0.65–1.35 GeV². The mass fit is restricted to 0.70–1.30 GeV² and the peak mean to 0.78–1.05 GeV², so a large broad peak at higher missing mass is not mistaken for the spectator peak. The allowed Gaussian width is 0.01–0.20 GeV². Histograms use 0.005 GeV² bins.

Momentum plots show 0–0.30 GeV. The fit range is 0–0.25 GeV, with peak mean 0.005–0.18 GeV and width 0.005–0.15 GeV. Histogram bins are 0.0025 GeV. The nonnegative distribution and its nuclear-motion tail need not be Gaussian. These limits are explicit starting assumptions for review, not established detector resolutions.

A bin gets no usable cut if it has fewer than 100 entries in the fit range, fails to converge, reaches a mean or width bound, or has less than 9 improvement in Poisson deviance over the same background-only model. This last requirement is a practical peak-resolution screen, not a calibrated significance. Such bins are marked clearly and do not contribute to downstream selected distributions. Do not interpret an empty final panel from such a bin as a measured zero yield.

Poor descriptions with chi-square/ndf > 5 and selections extending beyond the fit interval are flagged. Their fitted windows are still applied **provisionally**, following the note's approach, so their behavior can be examined. These results require review before adopting cuts for asymmetry extraction. Peak fit errors are local statistical approximations conditional on the selected model. They do not include uncertainty from the choice of background or signal shape.

No cuts are finalized automatically. This pass writes fit results and diagnostic cut flows, not analysis-selected event trees. Nitrogen contributions in ND3 can remain after the selections and still require the dilution study. A narrow angular residual does not by itself establish sample purity.

## Inputs, bins and selection

The default six input files and tree name come from `../../imports/datasets.json`. Shared corrected four-vector calculations are imported from `../diagnostics/kinematics.py`. Saved beam energy is used for each candidate. Proton is p1 and pion is p2. No momentum correction is applied again. Reconstructed neutron Mx2 must agree with the stored value to within 1e-6 GeV² or processing stops.

Default requirements remain fiducial_status == 3, FD pion, FD or CD proton, Q² ≥ 1 GeV², W ≥ 2 GeV, y ≤ 0.80, and finite required inputs. Bins are the same as pi+: xB edges [0.10, 0.25, 0.35, 0.45, 0.60] and −t′ edges [0.05, 0.25, 0.45, 0.65, 0.85, 1.05, 1.25] GeV². Intervals include the lower edge and exclude the upper edge. Signed saved tprime is negated. Bin numbers run in xB-major order. Histograms count candidate combinations, not unique events. No polarization or luminosity weights are applied.

## Outputs

Each period gets **12 PNG pages**:

- `mass_xb1.png` through `mass_xb4.png`: missing mass fits and cuts before exclusivity selections.
- `momentum_xb1.png` through `momentum_xb4.png`: missing momentum fits after the mass selection.
- `angle_xb1.png` through `angle_xb4.png`: signed polar-angle residual after both selections.

Each page has six −t′ panels for a single xB interval. The six periods produce 72 PNG pages in total, compared with 360 in the initial diagnostics. No PDF copies are generated.

`cuts.csv` records both stages' fitted means, widths, approximate errors, bounds, goodness-of-fit values, warnings and per-bin cut flow. `mass_fits.json`, `momentum_fits.json` and `angle_guides.json` preserve the complete fit results. `manifest.json` contains input provenance, source identity, code signature, full/partial sample flag and per-bin counts at each stage. The counts include underflow, overflow and undefined angular values. All retained candidates enter the cut flow even if outside the displayed range. `histograms.npz` caches the histograms and counts. `progress.json` tracks the current pass. Top-level `index.html` links the three stages and numerical tables.

Redraw figures from completed caches without rescanning ROOT data:

```bash
python3 analysis_scripts/RGC_eppi-/channel_selection/exclusivity/run_selection.py --plots-only
```

Use `--input-root PATH`, `--output PATH` or `--periods DATASET ...` to override inputs, destination or dataset selection. `--force` rebuilds all selected period caches.

## Validation

```bash
python3 -m unittest discover -s analysis_scripts/RGC_eppi-/channel_selection/exclusivity -v
```

Tests recover known mass and momentum peaks with backgrounds, reject missing peaks and insufficient statistics, verify cut boundaries and histogram accounting, and verify that the angular guide does not supply a cut. A separate synthetic ROOT integration checks the complete three-pass workflow and rendering. Production cuts must be inspected after the farm run.
