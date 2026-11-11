# Optional low-|t| extraction tests

Nominal `./dvcs_analysis` behavior and imports are preserved. The new modes
require `--skip-systematics`, a full extraction, and the isolated launcher.
They cannot be combined with each other or the existing study modes.

Build once in the source directory (farm ROOT/GCC environment):

```sh
make -j10 dvcs_analysis
python3 test_t_overrides.py
```

Run separately, or on two nodes sharing the source directory:

```sh
python3 run_t_override.py --merge-low-t --skip-systematics --run-dir /path/to/runs/merged_t_001 --plot
python3 run_t_override.py --photon-t --skip-systematics --run-dir /path/to/runs/photon_t_001 --plot
```

Use a NEW absolute run directory for each invocation; existing destinations are
refused. Build before launching jobs. The launcher snapshots small code files,
executables and external Python helpers, links the read-only imports, and runs
with a private working directory. ROOT event samples stay at their existing
absolute shared paths. Each run owns its `output`, cut JSONs, model caches,
plots and metadata (`.dvcs_t_test_run`). No second checkout or ROOT-file copy is
needed. Keep shared imports fixed until both jobs complete.

The C++ stages retain their existing concurrency limits (up to ten workers;
five independent periods in the exclusivity stage). BLAS pools are limited to
one thread per worker; OpenMP defaults to ten threads. Allocate ten CPUs per
job. Two nodes still share storage bandwidth, so speedup depends on farm I/O.
There is no extra ROOT scan to implement either override.

## Merged-bin mode

The CSV initializer takes the union of valid 0.11–0.15 and 0.15–0.25 cells
with matching xB, Q2 and phi edges, yielding 0.11–0.25 cells. It selects the
valid upper-bin row as its pass-1 correction donor. Where no valid upper row
exists, it preserves the valid lower row and its factors. This is explicitly a
loose correction approximation, not an average of finished cross sections.
With the current import: 253 merged cells, 230 upper donors, 23 lower-only
fallbacks; total analysis rows decrease from 1445 to 1394. Donor counts are
printed during initialization. Unaffected rows retain their imported values.

Counts, means, acceptance, contamination, acceptance-corrected yields and
physical bin volumes are recalculated using the new bin edges. Imported
model/correction fields retain the donor approximation; the external final
plotter computes fresh BH/KM15 curves for the resulting cells.

## Photon mode

Bin membership and mean |t| use |(q-q')^2| computed from electron/photon
momenta and the beam energy for each period. Data and reconstructed MC use
measured momenta; paired generated MC uses gen_* momenta. All required
momentum branches are enabled and checked. The generated-count cache key
includes the production t definition. Existing proton-dependent event cuts
and detector efficiency weights retain their definitions, so this tests the
bin coordinate with the same event selection. Pass-1 correction factors are
still reused as a deliberate approximation. The original proton t branch is
not overwritten. Existing Sebastian diagnostics retain their separate flag.

## Final comparisons

`--plot` runs `external_scripts/plot_pass2_final_results_by_t.py` with ten
model workers, statistical errors, `--no-pass1` to focus on BH/KM15, and
explicit `--reference-kinematics means` for these approximate override studies.
Nominal production plots instead default to the imported Fbin reference.
Systematic columns are optional unless total-point-to-point errors are
requested. Models are cached within each run's output directory. The script
still requires the farm's usual dvcsgen/KM15 installations.

To plot later, from the private run directory:

```sh
python3 external_scripts/plot_pass2_final_results_by_t.py output/csvs/dvcs_pass2_analysis.csv --output-dir output/final_results_by_t --workers 10 --no-pass1 --reference-kinematics means
```

Omit `--xs-total-ptp` and `--bsa-total-ptp` for these runs. Comparison curves
are evaluated at representative bin kinematics; these are diagnostic tests,
not a complete resolution unfolding or a rederived photon-t correction set.
