# Standalone low-|t| migration and closure study

Run `study_low_t_closure.py` separately from `dvcs_analysis`. It defaults to
Fall-2018 inbending **and** outbending and uses at most ten native workers in
total, shared across the two periods. Results remain separate by period.

## Production inputs and selection

`--production-run` must identify the completed **nominal proton-t** run. The
study reads and snapshots its:

- `output/jsons/global_cuts_config.json`
- `output/jsons/combined_cuts.json`
- `output/dvcs_current_dependence/calibration/current_response_model.json`
- `output/csvs/dvcs_pass2_analysis.csv`
- saved `global_cuts.cpp`, `global_cuts.h` and `kinematic_bin_integration.h`

The small backend compiles against that run's **actual production global-cut
implementation**, including reconstructed electron/proton thresholds, physical
phase space, fiducials, sectors and any enabled ycol selection. It applies the
saved **MC** exclusivity windows for the appropriate period and reconstructed
topology, not the DATA windows and not newly optimized cuts. Window boundaries,
upper-quantile behavior and Delta_phi fallback match production counting.

Inputs are the original paired `MCEvents` trees in:

`/volatile/clas12/thayward/dvcs_pass2_full_rebuild/mc/<period>/dvcs_mc/50nA/parts/`

Companion `.rec_only.root` and incomplete files are excluded. No independent
truth/reconstruction files are joined: both records come from the same paired
entry. Required branches are checked explicitly. Particle-angle branches and
`phi2` are radians; production opening angle and saved exclusivity variables
use their existing branch conventions without an extra conversion.

Generated denominators use the production generated-channel selection, without
applying reconstructed detector cuts to truth. Counts are unit generated
counts, as in production. Both reconstructed halves receive the saved DVCS MC
current weights. DATA current slopes and Neupane are not applied to simulated
pseudo-data. Calibration parameters are held fixed in this conditional closure
test; their fit uncertainty is not propagated.

## Tests performed in the same input pass

Each file is read once. Entries are deterministically split into disjoint
response-training and pseudo-data halves. Train acceptance is computed from
reconstructed/generated yields in the actual **xB/Q²/|t|/phi bins**, including
reconstructed feed-in from outside the generated analysis bins. Pseudo-data
are corrected with that independently estimated acceptance and compared with
their known generated yields.

Five truth spectra are tested, including smooth exponential changes in |t| and
a localized first-bin enhancement. Weights vary **within bins** and apply to
outside-bin feed-in too. This avoids a closure test that merely reuses the
training spectrum's shape.

Three binning choices are accumulated simultaneously:

1. Nominal production bins.
2. A matched-coverage control, retaining the two low-|t| bins over the same total support as the merged case.
3. A merged 0.11–0.25 GeV² bin with the same coarsened phi partition.

The first two t bins sometimes have different phi partitions. Overlapping
intervals must be coarsened before merging; concatenating the original rows
would double count events. Comparing choices 2 and 3 isolates the t-merging
effect from phi/support changes. The control can contain more rows than the
nominal table because it includes the full merged support even where a
nominal low-t row was absent. Closure precision thresholds still apply. These are diagnostic binnings; the
production extraction is not changed.

The script also accumulates migration matrices, subthreshold feed-in,
reconstructed-topology purity/stability and truth-proton-angle slices.
Topology-conditional migration fractions are **not** absolute topology
efficiencies. Truth-angle closure is a conditional MC diagnostic, not a method
to assign unknown true angles to data.

No KM15 evaluation, radiative/centering-factor integration, data extraction or
exclusivity fitting is needed. Those reporting factors cancel in the closure
ratio. This study diagnoses correction/migration sensitivity, not absolute
DATA/MC detector-efficiency calibration or radiative-factor validity.

## Farm command

From `dvcs_cross_section`, this runs a small synthetic ROOT smoke test first,
then scans both production MC periods. No full-analysis rebuild is needed:

```sh
python3 test_low_t_closure.py --root-smoke && python3 study_low_t_closure.py --production-run "$PWD/test_runs/current_mc_fixed_001/proton_t_default_001" --periods fa18_inb fa18_out --workers 10 --output "$PWD/test_runs/low_t_closure_fa18_001"
```

The script automatically compiles its small worker once using `root-config`
and the C++ compiler, then shares it across all jobs. It needs the existing
ROOT development environment plus Python NumPy and Matplotlib. It does **not**
require uproot, PyROOT or a Python model package.

Use a new output directory. Existing destinations are refused. Input analysis
files are read only. A manifest records cut/source hashes, seed, precision
thresholds, spectra, file sizes and modification times. Files changing during
the scan cause an error. If a worker fails, queued work is cancelled and active
workers are stopped; details remain in `worker_results/job_*.log`.

For an input/schema check against a few real parts, use the same command with
`--max-files-per-period 1` and a separate output directory. Such a limited scan
must not be interpreted as the full study.

## Outputs and interpretation

- `REPORT.md`: definitions, first-bin summary and plots.
- `closure_by_cell.csv`: all period/binning/scenario/angle results, including
  cells below the precision threshold.
- `closure_summary.csv`: median, 16–84% cell spread and summed-yield closure.
- `matched_low_t_comparison.csv` and `.png`: compare the control and merged
  corrections in exactly the same x/Q/phi cells where both separate t bins
  and the merged bin pass the precision thresholds. The worker checks that
  generated and reconstructed sums match between the control and merged
  domains before making this comparison.
- `closure_nominal_vs_merged.png`: nominal/control/merged comparison for both periods.
- `closure_by_truth_proton_angle.png`: first-bin closure by generated angle.
- `migration_response.png`, `response_counts_<period>.csv`,
  `migration_metrics.csv`: true-to-reconstructed t response, purity/stability
  and feed-in below 0.11 GeV². These t projections cover the generated-channel
  sample after production global/exclusivity selection; closure itself uses
  the actual four-dimensional analysis bins.
- `migration_by_topology_and_truth_angle.csv`: conditional migration metrics.
- `cutflow.json`: input, validity and selection counts. Reconstructed entries
  without valid truth are recorded and excluded from this paired-truth study.
- `aggregates_<period>.npz`: compact reusable sufficient statistics.
- `production_inputs/`, `spec.json`, `manifest.json`: exact provenance.

Default usable-cell thresholds are 200 training generated events, 100 nominal
test generated events and 20 effective reconstructed events in each half, plus
the saved production valid-bin flag. All cells are still written. Error
estimates include same-cell truth/reconstruction covariance and uncertainty
in the trained acceptance. They omit inter-cell covariance and calibration
fit uncertainty; the 16–84% spread is **cell dispersion**, not an error band on
the median. A successful stress test is evidence of robustness for those
spectra, not a proof for every possible underlying distribution.

## Replot without another ROOT scan

```sh
python3 study_low_t_closure.py --production-run "$PWD/test_runs/current_mc_fixed_001/proton_t_default_001" --periods fa18_inb fa18_out --reuse-aggregates "$PWD/test_runs/low_t_closure_fa18_001" --output "$PWD/test_runs/low_t_closure_fa18_replot_001"
```

This skips compilation and ROOT I/O. Cut/binning/current-response/worker
provenance is checked before reuse. Plot settings and precision thresholds
can change without rescanning events; changing the event definitions or
spectrum stresses requires a new scan.

## Validation

`python3 test_low_t_closure.py` tests nonoverlapping merged/control grids,
closure/error algebra, actual C++ accumulation and exclusivity-window helpers,
production phase-space/proton-threshold helpers, deterministic splitting,
feed-in, two-period report rendering and reuse/overwrite guards without
initializing ROOT. `--root-smoke` additionally creates real paired ROOT files
and runs both polarities through the native worker, checking that the saved
MC windows rather than DATA windows govern selection.

Local helper/report tests and compilation pass. The local ROOT interpreter
has stale Cling headers and crashes even creating the synthetic input file,
so the ROOT end-to-end smoke test must be completed in the farm environment.
