# FD–FD fixed-response background-transfer diagnostic

Run from dvcs_cross_section, after committing/pulling the new files:

```sh
python3 test_fd_background_transfer.py --root-smoke && python3 study_fd_background_transfer.py --common-study "$PWD/test_runs/particle_response_common_cuts_001" --workers 10 --scratch-root /scratch/thayward/dvcs_particle_response --output "$PWD/test_runs/fd_background_transfer_001"
```

The default covers all five periods individually, FD–FD only. Empty first-|t|
cells are explicitly listed. Use `--periods fa18_inb` for a restricted run.
Add-pT response parameters are frozen and applied to signal and leakage MC.
No scan, new response fit, extra histogram morph or production cut update occurs.
No nominal/original-training evaluation is run.

The completed common-cut study supplies the fixed period background scale and
FD–FD signal anchor, including the control-to-leakage exposure factor and saved
background multiplier. The two-photon control is **not** reread or refitted here.
This tests whether that control-derived transfer describes the FD–FD population;
it does not independently remeasure the transfer or fit away a failing budget.

Existing 12-column pre-cut exports lack particle coordinates. The new job enables
an optional 20-column export; previous exports and readers retain their 12-column
format. It exports only FD–FD/full-global events, with the eight extra coordinates
xB, Q², e p/theta, proton p/theta, photon p/theta, all after the fixed response in MC.
DATA are read once per period. MC use the verified existing sampled caches; no
new MC ROOT scan, two-photon reconstruction, or sampling is performed. Other
MC topologies are skipped before four-vector variation. At most 10 workers run
across DATA and both MC samples combined. Extended binaries stay under scratch;
only report files need transferring. Checksums are streamed with bounded memory.

Open REPORT.html. In each period and |t| bin, start with **shell_budget** and **inclusive**, then inspect
proton momentum/angle, photon momentum/angle, electron momentum/angle, xB and Q².
Each subdivision axis uses three DATA quantiles before exclusivity in that |t|
bin; the exact same edges are applied to MC. The axes are tested separately,
not as a joint multidimensional reweighting. Tables record edges and counts.

Top plots compare equal-area DATA (dashed) and the physical summed MC (solid),
with one color per kinematic subdivision. The same saved mixture anchors are
used throughout; equal-area display rescales the sum and its components together.
Bottom plots show R(c)/R(reference), where R = (DATA − predicted pi0)/DVCS MC.
One at the reference is guaranteed; flatness elsewhere is the test. Default
windows are the completed closure study's refitted DATA boundaries, applied
identically to all samples. `--window-family saved_DATA` repeats with the older
common DATA boundaries without changing or refitting either set.

`shell_budget.png` plots DATA versus predicted π⁰ in disjoint intervals, with
their difference below (negative residuals highlighted red). No DVCS signal
shape/normalization enters this budget. `shell_budgets.csv` is the central transfer check. For each disjoint missing-pT
interval (low,high] it records observed DATA, raw leakage MC, predicted background,
raw signal MC, background-subtracted DATA and its conditional statistical error.
Negative residuals are retained. The quantity DATA/predicted background is a point
budget diagnostic, not a measured transfer factor or confidence limit. Empty and
low-statistics cells cannot establish a failure. Conditional uncertainties exclude
background-transfer, response-fit and cut-window model errors.

`summary.csv` reports reference counts/background fractions and relative signal
scales at approximately 0.15 GeV missing pT and 0.01 GeV² mass half-width. Its
low-statistics flag requires at least 100 DATA and 50 sampled signal-MC events at
both test and reference thresholds. Other samples and subtraction uncertainty
must still be inspected. `scans.csv` supplies every threshold. `coverage.csv`
records the available global-selected FD–FD events, including empty cells.

Follow-up interpretation:

* An added interval with predicted background exceeding total DATA is a direct
  budget mismatch, conditional on transfer normalization and statistics.
* Drift disappearing within kinematic subdivisions suggests changing mixtures
  contribute; persistent local drift narrows the response/background problem.
* Neither outcome alone attributes an error to one track, validates absolute
  acceptance, nor yields a final cross-section correction.

To repeat the alternative windows without ROOT reads or response reevaluation:

```sh
python3 study_fd_background_transfer.py --common-study "$PWD/test_runs/particle_response_common_cuts_001" --reuse-exports "$PWD/test_runs/fd_background_transfer_001" --window-family saved_DATA --output "$PWD/test_runs/fd_background_transfer_saved_cuts_001"
```

The input periods, parameters, normalization and worker fingerprints must match.
Existing output directories are never overwritten.
