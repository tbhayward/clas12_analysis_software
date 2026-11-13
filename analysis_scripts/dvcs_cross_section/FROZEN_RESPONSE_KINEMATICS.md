# Frozen joint-response low/high diagnostics

```sh
python3 test_frozen_response_kinematics.py && python3 study_frozen_response_kinematics.py --joint-run "$PWD/test_runs/particle_response_joint_control_001" --workers 10 --scratch-root /scratch/thayward/dvcs_particle_response --output "$PWD/test_runs/frozen_response_kinematics_001"
```

Implements only the first agreed step: freeze parameters and locate kinematic
residuals. No new response parameters, fit, per-cell normalization or cut changes.
All five periods and all three topologies are retained. Predictions include
nominal MC and the fixed jointly fitted candidate.

Categories: FD electron, FD/CD proton, photon 1 FD/FT, photon 2 FD/FT. Momentum
and polar angle are split separately, not into four joint cells. Each numerical
boundary is the median of globally selected epgamma DATA half 0 before exclusivity
and before any agreement score is computed. The same boundary is used for DATA,
all MC, and both channels. Both photon roles use the SAME detector-specific edge;
photon 2 is displayed only for epgammagamma. The response law is shared between
photon roles; their independent random draws differ. Missing split support is
flagged, not optimized away. The soft second-photon population may have few
entries above a split defined by epgamma; low-statistics cells are flagged.

For each period/category/axis, REPORT.html links:

* scores: nominal-to-candidate discrepancy reduction on MC half 1, separately for
  channel/topology/low-high/observable. Gray cells lack support.
* residuals: frozen-candidate DATA/MC core-width ratios (target 1) and median
  differences in MC core-width units (target 0). Compare low with high and both
  channels; improvement alone does not imply agreement.
* examples: largest supported candidate residual in each channel/side. These are
  chosen for inspection only, not to choose split edges or fit parameters.

scores.csv contains every score for both MC halves, sample counts, width ratios,
median shifts and histogram-range coverage. splits.csv records exact thresholds
and their DATA source/support. photon_coverage.csv records unknown individual
photon coordinates. Unknown photon pairs retain their events for electron/proton
checks but are excluded from photon-specific subdivisions. Photon-role plots
share events and must not be summed. The composite pi0 is never used in place of
an individual photon.

All observables use full saved global selection and common DATA exclusivity
windows, omitting the displayed variable's cut. Epgamma predictions contain both
DVCSgen and pi0 leakage with the fitted normalization held fixed across groups;
epgammagamma uses AAOGEN. Each curve is displayed at unit area. Shapes and scores
are descriptive, not acceptance measurements or independent chi-square p-values.
DATA are shared between MC-half comparisons. Momentum, angles, and other kinematics
remain correlated, so these projections do not identify a unique particle cause.

The exporter adds an optional 28-column mode containing measured individual-photon
coordinates and tags. Existing 12/20-column formats are unchanged. MC use existing
uniform caches; DATA are read once per channel and period. Only globally selected
rows are saved. Candidate and nominal MC exports are independent deterministic
reevaluations of the same cached events. Large arrays remain on scratch, with
streamed checksums and a ten-worker total cap. Reused group masks avoid rebuilding
particle membership for every observable. Baseline DATA selection must reproduce
the frozen fit's recorded selected counts before evaluating its period.

To replot with a changed minimum-entry threshold without new event I/O:

```sh
python3 study_frozen_response_kinematics.py --joint-run "$PWD/test_runs/particle_response_joint_control_001" --reuse-exports "$PWD/test_runs/frozen_response_kinematics_001" --min-entries 200 --output "$PWD/test_runs/frozen_response_kinematics_002"
```

Requires matching input/candidate/export-source hashes. Export files must remain
on scratch. Transfer the report directory only. Completed production outputs and
joint-fit results are never overwritten.
