# Exclusive e′pπ− analysis

The development analysis now defaults to the relaxed campaign in
`/volatile/clas12/thayward/eppi-/loose/full_20261124/root/`, using the 19
`eppi_minus_loose_DATASET.root` files listed in `imports/datasets.json`.
Wait for **19/19 merged datasets**, not just 1076/1076 completed batches.

## First rerun on ifarm

Copy the updated development analysis directory to your farm checkout. No commit is required.
From the repository root, inside a manually opened tmux session:

```bash
python3 analysis_scripts/RGC_eppi-/run_relaxed_analysis.py --workers 6
```

This runs the existing channel calibration and 19-file export, selection audit,
kinematic plots, dilution calculation (including the 2,000-replica shared-period
fit), and the fixed-zero-UU polarized extraction, in that order. Any failed stage
stops the sequence. Python dependencies are numpy, scipy, uproot 5, matplotlib and
iminuit. threadpoolctl is optional. No ROOT Python installation is needed.
Run on one node. Do not launch two copies against the same output directory.

The first pass keeps all four xB bins and six −t′ bins for all five polarized
ratios. Both UU ratios remain fixed at zero. Decide target-sector coarsening after
inspecting the new errors. No theory inputs or old dilution factors are imported.

All new products are isolated under `analysis_scripts/RGC_eppi-/output_loose_20261124/`:

- `channel_selection/`: fresh fitted windows, selected ROOTs, cut flows and gallery.
- `kinematic_distributions/`: selected shapes and bin means.
- `dilution_factor/`: recalculated corrections, replica errors and pooling comparisons.
- `asymmetry_extraction/`: five polarized ratios, period checks and production PNGs.

Earlier `channel_selection/output`, `dilution_factor/output` and
`asymmetry_extraction/output_fixed_UU` results remain available. Matching selection
and extraction checkpoints are reused when restarting the wrapper. Dilution is
recomputed if its stage is rerun. Individual scripts retain their historical output
defaults, so use this wrapper for the isolated relaxed-data rerun.

Print the underlying commands without running anything:

```bash
python3 analysis_scripts/RGC_eppi-/run_relaxed_analysis.py --dry-run
```

Run an individual stage, for example after selection has finished:

```bash
python3 analysis_scripts/RGC_eppi-/run_relaxed_analysis.py --stage audit
python3 analysis_scripts/RGC_eppi-/run_relaxed_analysis.py --stage kinematics --workers 6
python3 analysis_scripts/RGC_eppi-/run_relaxed_analysis.py --stage dilution --workers 6
python3 analysis_scripts/RGC_eppi-/run_relaxed_analysis.py --stage asymmetry --workers 6
```

Selection progress:

```bash
python3 analysis_scripts/RGC_eppi-/channel_selection/run_channel_selection.py --output analysis_scripts/RGC_eppi-/output_loose_20261124/channel_selection --status
```

Use `--output-root PATH` on the wrapper to put all analysis products on a larger
filesystem. Repeat that option for subsequent stage commands. `--input-root PATH`
overrides the merged ROOT directory and `--run-info FILE` overrides charge/polarization
metadata. Metadata checks remain strict, including runs with zero selected candidates.
Do not replace absent charge information with selected candidate counts.

## Selection and conventions

This rerun retains the existing offline selection: finite required kinematics,
`fiducial_status == 3`, FD pion, FD/CD proton, original DIS and bin cuts, fitted
missing-mass and missing-momentum windows. RGB fall 2019 requires Q² > 1.5 GeV².
Other periods retain Q² ≥ 1 GeV². Windows are recalibrated from the new sample.
The old processing track cuts are **not** restored using the new `pass_*_legacy`
flags. The relaxed sample can therefore admit candidates previously rejected by
processing. This is a new selected sample, not a pure luminosity increase.
All original scalar branches, including the added diagnostic flags, survive export.

The detected proton is p1 and π− is p2. The production angle is phi2, not p2_phi.
RGB has Pt=0. Only RGC ND3 carries target polarization. Period-specific charge
normalization is preserved. Counts refer to candidate combinations, not unique events.
A larger raw entry count does not predict the final selected polarized precision.

See the component READMEs for calibration, dilution and likelihood details.

## Directory guide

The active code lives in `channel_selection`, `kinematic_distributions`,
`dilution_factor` and `asymmetry_extraction`. `imports` holds input configuration
and run metadata. `analysis_note` holds the paper source. Superseded selection
prototypes live under `archive/channel_selection`. Empty future-work placeholders
have been removed. Earlier numerical outputs and the UU sensitivity study are
retained for comparison.

Readable selection progress for this rerun:

```bash
python3 analysis_scripts/RGC_eppi-/run_relaxed_analysis.py --status
```

This reports scans, fitted groups, completed exports and per-dataset progress.
It does not estimate remaining wall time or report the later fitting stages.

## Compact review copy in home

The wrapper now updates `review_loose_20261124/` and
`review_loose_20261124.tar.gz` beside itself after every successful stage, even
when full outputs use `--output-root` on volatile. Override with `--review-output`.
Figures, fit JSONs, tables, small covariance/histogram NPZs and input metadata are
copied. ROOT files, calibration caches, NPY arrays and replica ensembles are omitted.
The review copy is for inspection, not a standalone input for refitting.

Export an already completed run without rerunning analysis:

```bash
python3 analysis_scripts/RGC_eppi-/export_review.py --source /volatile/clas12/thayward/eppi-/loose/full_20261124/analysis
```

Copy the resulting tar.gz to your desktop. The source files remain on volatile.

## Three-bin target-sector study

After the nominal relaxed run, use
`asymmetry_extraction/three_bin_study/run_three_bin_study.py --source /volatile/clas12/thayward/eppi-/loose/full_20261124/analysis --workers 12`
with Python. This recomputes dilution in 12 bins and compares free LL cosφ with
BM stage-2 and stage-3 fixed inputs. The six-bin LU result is retained. See the
study README for assumptions, outputs and the compact home review archive.
