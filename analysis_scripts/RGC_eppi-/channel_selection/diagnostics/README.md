# Exclusive pi− exclusivity diagnostics

This standalone analysis reads the completed `PhysicsEvents` trees. It does not change the processing campaign or apply an exclusivity cut. The default inputs are the three RGC ND3 periods and the three RGB D2 periods listed in `../../imports/datasets.json`.

## Run on the farm

From the repository root, on **one node**, with Python 3, NumPy, uproot 5 and Matplotlib available:

```bash
python3 analysis_scripts/RGC_eppi-/channel_selection/diagnostics/run_diagnostics.py --workers 6
```

All generated PNG figures, histogram caches, tables, progress files and the HTML index are saved to `analysis_scripts/RGC_eppi-/channel_selection/diagnostics/output/` beside these scripts. This path is anchored to the script location, regardless of the working directory. On the farm it is therefore inside your home analysis repository, not on volatile storage. The input ROOT files continue to be read from the campaign on `/volatile`. The output folder is excluded from git. An explicit `--output PATH` can override the destination.

Run the command once. Six processes scan the six files concurrently, one process per period. Each file is read once in bounded chunks (default `100 MB` of decoded branch data). Plotting follows the scans, one period at a time. ROOT itself and the processing environment are not required. If needed, install the packages in your usual Python environment:

```bash
python3 -m pip install numpy 'uproot>=5,<6' matplotlib
```

For an initial end-to-end smoke test, add `--max-entries 200000 --output analysis_scripts/RGC_eppi-/channel_selection/diagnostics/output/smoke_test`. This reads the first 200,000 stored candidate rows per period, not a representative run sample. Figures are marked PARTIAL SAMPLE. Omit this option for the actual diagnostic distributions.

For a progress report in another terminal:

```bash
python3 analysis_scripts/RGC_eppi-/channel_selection/diagnostics/run_diagnostics.py --status
```

Workers print START immediately and report opening input, reading the first chunk and processing. Percentages are exact fractions of stored candidate rows read per period. The update age helps identify a stopped scan. After scanning, status distinguishes plots pending, plotting and complete. If interrupted, rerun the original command. Completed matching histogram caches are reused, while an interrupted period restarts its scan. Only run one scanner per period/output directory at a time. Previous results on `/volatile` are not moved or deleted by this change.

To redraw plots without reading the input ROOT files:

```bash
python3 analysis_scripts/RGC_eppi-/channel_selection/diagnostics/run_diagnostics.py --plots-only
```

Add `--log-y` for logarithmic one-dimensional count axes. `--input-root PATH` overrides the input ROOT directory. `--periods` can select any dataset listed in the dataset manifest, including carbon and the other control targets.

## Definitions

All four-vectors use the saved corrected lab momenta and each row's beam energy. Angles in the ROOT tree are radians. Proton is p1 and pion is p2. No momentum correction is applied again.

- Neutron reference: `Mx2_n = (k + Pn - electron - proton - pion)^2`, with a stationary neutron. This is checked against saved `Mx2` with a tolerance of 1e-6 GeV². A mismatch stops the scan.
- Deuteron reference: `Mx2_d = (k + Pd - electron - proton - pion)^2`, with a stationary deuteron. The line at mp² is the on-shell spectator-proton expectation for exclusive production on deuterium. This is also a reference for the deuterium contribution in ND3, not a nitrogen-target mass reconstruction.
- Predicted pion: `pi_pred = k + Pn - electron - proton`. Its energy is retained from this subtraction. It is not reset to an on-shell pion energy.
- Residuals are **measured minus predicted**: energy, momentum magnitude, Cartesian momentum components, polar angle and wrapped lab azimuth. The azimuth difference is wrapped to [−180°, 180°]. This is a lab track-angle residual, not the production angle `phi2`.
- Further distributions show the opening angle between measured and predicted pion directions, `pi_pred² - mpi²`, and the transverse missing momentum relative to the beam.

These observables share the same imbalance. In particular `Mx2_n = delta_E² - delta_px² - delta_py² - delta_pz²`. They are not independent measurements for a chi-square. The stationary-neutron prediction does not reconstruct the event's unmeasured spectator momentum. Fermi motion, binding, final-state interactions and nuclear backgrounds can broaden signal distributions. The deuteron missing mass and stationary-neutron missing mass contain different assumptions, so both are retained.

## Selection and bins

Default selection is `fiducial_status == 3`, an FD pion (`detector2 == 1`), and either an FD or CD proton (`detector1 == 1 or 2`). The existing Q² ≥ 1 GeV², W ≥ 2 GeV and y ≤ 0.80 requirements are retained. No further helicity, polarization, exclusivity or pion-angle selection is added. Nonfinite required inputs are counted and excluded. `--all-fiducials` includes candidates that fail the saved fiducial status for checks.

The pi+ bin edges are used unchanged, with half-open intervals including the lower edge and excluding the upper edge:

- xB: 0.10, 0.25, 0.35, 0.45, 0.60.
- −t′ (GeV²): 0.05, 0.25, 0.45, 0.65, 0.85, 1.05, 1.25.

Bin number is `6*x_index + t_index + 1`. The signed saved `tprime` is negated. Out-of-bin rows are counted. Candidate proton–pion combinations are counted, not unique events. No charge, polarization or acceptance weighting is applied. Raw yields across periods are not normalized comparisons.

## Outputs

Each of the six period directories contains:

- 48 PNG distribution pages: 12 observables × four xB intervals, with six −t′ panels per page. Each panel overlays all, FD-proton and CD-proton counts. FD + CD equals all. Dashed lines indicate the nominal reference, not an imposed cut.
- 12 PNG correlation pages: neutron missing mass versus energy residual, energy residual versus missing transverse momentum, and momentum residual versus polar-angle residual. Four xB pages per correlation show all proton detectors together, with logarithmic count colors.
- `statistics.csv`: finite counts, mean, standard deviation, underflow, overflow and undefined values in every bin, observable and proton category. Moments include the tails outside the displayed plotting range. The displayed upper histogram edge belongs to overflow.
- `histograms.npz`: compact histogram and moment cache. This permits later fitting or replotting without another ROOT scan.
- `manifest.json`: input identity, source entry count, selection configuration, bin edges, observable ranges, count accounting, full/partial flag and saved-versus-recomputed Mx2 agreement.
- `progress.json`: latest scan progress.
- `plot_status.json`: plotting or completion state.

The top-level `index.html` links the PNG pages. Copy the output directory locally to browse it. No PDF figures are produced. Default axis ranges are in `kinematics.py`. The per-panel outside-range count and CSV tails must be checked before adopting a plotting range or cut. Angular residuals undefined for zero predicted momentum are counted separately rather than removing the other observables for that candidate.

## Validation

```bash
python3 -m unittest discover -s analysis_scripts/RGC_eppi-/channel_selection/diagnostics -v
```

Tests construct conserved gamma* n → pi− p events in all 24 bins, verify the missing-mass and residual identities, check bin boundaries, scan six synthetic TTree files in chunks, check selection accounting, FD/CD sums and cache reuse, and reject an inconsistent saved Mx2 convention. Production data still need to be scanned on the farm.
