# Cosine recovery at the observed phi coverage

`validate_cosine_coverage.py` scans injected `u1 = F_UU^cos(phi)/F_UU` and
`ll1 = F_LL,lab^cos(phi)/F_UU` independently. All seven physics coefficients and
three period dilution nuisances are free in each production-likelihood refit.
Nothing is written to the production results or the earlier validation directory.

## Defaults and scope

- Bins 4, 5, 6: the three highest -t-prime bins in 0.10 < xB < 0.25.
- Each cosine scans [-0.5, -0.25, 0, 0.25, 0.5], giving 25 pairs per bin.
- 100 attempts per physical pair: at most 7,500 main toy fits, plus constrained
  profile fits. Nonphysical inputs are explicitly recorded, never clipped or
  rescaled. Grid axes are configurable independently.
- Other injected coefficients and true period dilutions come from the saved
  nominal fit in `output/unpolarized_validation/checkpoints/bin_XX.json`.
- In `observed` mode every measured event location and period is retained,
  including the actual phi holes, kinematic correlations and selected counts.
- Run and helicity are jointly regenerated from the known conditional
  probabilities using measured charges and polarizations. Three new Gaussian
  auxiliary dilution measurements are generated per trial and profiled in fits.
- The fitter's five production starts are retained. Its initial u1 and ll1 are
  zero at every grid point; injected target coefficients are not passed as seeds.
  Other starting values are fixed across the grid at the baseline fit values.
  An accidental equality between a fixed start and a grid value is not seeding
  with truth. No coefficients are fixed during the main toy fit.

The experiment asks whether the available accepted coverage and counts allow
recovery under the assumed harmonic model and spin-independent acceptance. It
is not detector simulation and cannot fill an experimentally empty phi region.
It does not test omitted harmonics or spin-dependent acceptance.

## Run on both farm nodes

Copy the new script and this documentation into the repository on ifarm. The
script imports the existing `validate_unpolarized_modulations.py`, so keep that
file alongside it. The earlier validation needs to have a usable saved nominal
fit in each requested bin; it does not need to finish its other scenarios.

From the same shared RGC_enpi+ directory on **each** of ifarm2401 and ifarm2402:

```tcsh
env TMPDIR=/scratch/thayward OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 python3 asymmetry_extraction/validate_cosine_coverage.py --workers 10 --distributed --cache /scratch/thayward/rgc_enpi_thayward/selected_events_5774b50cd1143836.npz
```

The two invocations share grid cells through OS file locks: up to 20 concurrent
cells, rather than one whole bin per worker. Use the same shared output and input
paths on both nodes; the filesystem must support cross-node flock locks. Each
replica is saved atomically as soon as it finishes. Completed cells are reused;
partial cells resume without repeating saved replicas. Locks are released by the
OS on process exit; do not delete lock files while workers are active. Node-local
scratch caches can be rebuilt automatically if missing, using existing corrected
ROOT inputs. Do not rebuild a cache while another study is reading it.

Inputs, source code, baseline nominal fits, runtime versions, and grid settings
are fingerprinted. Changes require a new output directory; changing worker count
or using one/two nodes is allowed. The recognized previous migration-only
production-code upgrade can reuse its unchanged selection cache and baseline.
The new study does not alter the older study or discard its checkpoints.

For a smaller separate pilot (run the identical command on both nodes):

```tcsh
env TMPDIR=/scratch/thayward OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 python3 asymmetry_extraction/validate_cosine_coverage.py --bins 4 5 6 --u1-values -0.5 0 0.5 --ll1-values -0.5 0 0.5 --replicas 50 --workers 10 --distributed --cache /scratch/thayward/rgc_enpi_thayward/selected_events_5774b50cd1143836.npz --output-dir asymmetry_extraction/output/cosine_coverage_pilot
```

That is at most 1,350 main fits before physicality skips. Do not switch pilot
settings into the full-study directory. To scan the three highest-t-prime bins
in every xB range, supply `--bins 4 5 6 10 11 12 16 17 18 22 23 24` with a new
output directory.

## Optional matched-count forward-hole comparison

Add `--coverage-comparison` to both invocations, with a **different output directory**
if an observed-only study already exists there. It adds:

- `matched_uniform`: a uniform random subset of 70% of events in each period;
- `matched_forward_gap`: the same number per period, but no events within
  30 degrees of phi=0/360. Sampling uses no run or helicity labels.

Fractions and half-widths are configurable with `--retained-fraction` and
`--gap-half-width`. Insufficient eligible events are reported as
`coverage_unavailable`; counts are never silently reduced. The uniform and gap
controls have exactly equal period counts. Masks vary across replicas but are
independent of the injected coefficients, with common random ranks for the two
controls. Spin-generation and auxiliary seeds are paired between controls.
This separates losing forward-angle information from merely losing statistics.
It triples the maximum main-fit workload (22,500 attempts for the default grid).

The matched-comparison CSV uses replicas usable in both controls and reports
paired counts, the mean gap-minus-uniform recovery difference, and its paired
Monte Carlo standard error. Changes in fit-failure selection still need to be
considered alongside that paired comparison.

## Profile versus Hesse coverage

For each usable trial, the default also fixes each target coefficient separately
to its injected value and profiles the other six coefficients and dilution
nuisances. A warm start with zero target cosines and the fitted minimum provide
two deterministic starts for each constrained fit. It records

```
q = 2 * (NLL_profile_at_injected_value - NLL_best_fit).
```

The frequency of q <= 1 and q <= 3.84 measures coverage of the usual
one-parameter profile-LR confidence sets at approximately 68% and 95%. This
checks inclusion of the truth directly; it does not need to calculate both
endpoints of every asymmetric interval. These thresholds use the usual
asymptotic approximation, and the toys test its finite-sample performance.
Constrained fit failures and discoveries of a lower minimum are counted.
`--no-profile-coverage` skips these extra fits for a faster Hesse-only scan and
requires a separate output directory from a profile-enabled run.

The code also evaluates NLL at the injected physics and true dilution values.
If that point lies inside the fit bounds but has a lower NLL than the selected
minimum by more than 0.001, the trial is explicitly flagged as an optimizer
failure and excluded from usable-fit recovery. Such trials remain in the
all-attempt accounting. Truth is used for evaluation, not as an optimizer seed.

## Summarize, including while jobs run

From RGC_enpi+:

```tcsh
python3 asymmetry_extraction/validate_cosine_coverage.py --summarize-only
```

For the pilot add `--output-dir asymmetry_extraction/output/cosine_coverage_pilot`.
Only the summary command writes combined CSVs and PNGs during distributed runs;
run it from one terminal at a time. It can read the atomic checkpoints while
workers continue and needs no ROOT access for this step.

Outputs under `asymmetry_extraction/output/cosine_coverage_validation/`:

- `tables/cell_status.csv` and `.json`: completed/partial/nonphysical/unavailable
  cells, selected counts, start values, truth and reasons;
- `tables/injection_summary.csv`: mean recovery and bias with MC uncertainty,
  bias/data-stat ratio, empirical spread, reported errors, pull behavior, Hesse
  and profile coverage, all-attempt/usable/evaluated denominators and failures;
- `tables/toy_results.csv`: all seven fitted ratios for every attempted trial;
- `tables/matched_coverage_comparison.csv`, when controls are enabled;
- `plots/recovery_bin_XX_MODE.png`: recovered mean versus injected truth for each
  cosine, colored by the other injected cosine, with an identity line;
- `plots/bias_bin_XX_MODE.png`: the full two-dimensional bias grid;
- `plots/coverage_bin_XX_MODE.png`: Hesse and profile 68%/95% coverage;
- `plots/fit_quality_bin_XX_MODE.png`: usable fraction and empirical/fitted
  correlation between the two cosines;
- `plots/phi_coverage_bin_XX.png`: actual occupancy and a reproducible replica-0
  illustration of optional thinning;
- `plots/matched_coverage_difference_bin_XX.png`, when controls are enabled.

Recovery-plot error bars are the Monte Carlo standard error of the mean, not
one-experiment statistical errors. Gray heatmap cells are unsupported/missing or
have no usable fits; consult cell status. Hesse coverage among usable fits and
profile coverage among evaluated fits are shown separately from all-attempt
fractions. Do not treat good conditional means as evidence of stability if many
trials fail. No production correction or systematic uncertainty is assigned by
this diagnostic. PNG only; no PDFs are generated.
