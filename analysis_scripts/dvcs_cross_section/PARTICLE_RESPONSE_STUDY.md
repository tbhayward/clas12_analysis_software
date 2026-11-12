# Particle-level response test

This test varies reconstructed MC particle measurements, rather than adjusting
individual exclusivity distributions. It leaves the production analysis and
all existing outputs unchanged. Every run period/polarity is fitted independently:
Sp18 inb, Sp18 out, Fa18 inb, Fa18 out and Sp19 inb have separate parameter files.
No parameter is pooled across software versions or polarities.

## Run

From dvcs_cross_section after committing/pushing and pulling the change:

```sh
python3 test_particle_response.py --root-smoke --reference-study "$PWD/test_runs/exclusivity_resolution_all_periods_joint_001" && python3 study_particle_response.py --reference-study "$PWD/test_runs/exclusivity_resolution_all_periods_joint_001" --workers 10 --output "$PWD/test_runs/particle_response_all_periods_001"
```

The native worker compiles automatically; no full-analysis rebuild is required.
The reference supplies the saved global cuts/source snapshots, DATA histograms,
component scales, file inventory and exclusivity windows. DATA is not reread. Two-photon control MC is cached once on the first updated run, then reused. All five reference periods are processed; CD_FT is the
default fit topology. All three topologies are predicted for each period.

## Exactly what changes in a trial

For each of the electron (FD), proton (CD), and photon (FT), the default search
varies four quantities independently: extra fractional momentum/energy width,
extra polar-angle width, extra azimuthal-angle width, and fractional momentum
scale bias. Gaussian widths add independent random errors to reconstructed
measurements. They are not factors multiplying the existing resolution.

The particle mass is retained and energy follows its new momentum. Derived
Q2/xB/t, missing four-vector, transverse missing momentum, photon opening angle,
coplanarity, CM angle and electron-photon opening angle are recalculated together.
Stored quantities receive the recalculated *change* relative to original
four-vectors, preserving stored rounding/calibration offsets. A zero trial is
exactly the original event. The direct recomputation offsets are audited.
Detector identities stay fixed; measured fiducial angles/sectors are reevaluated.
Global selections, DATA exclusivity windows and joint kinematic cells are
reevaluated for every trial, with the plotted variable's cut omitted as before.
Signal, pi0 leakage, and two-photon control MC are varied consistently. Nominal plots retain reference scales; candidate scales follow the control/profile procedure below.

## Fitting and validation

Each period/polarity has its own bounded parameter-grid search, coordinate
updates and local refinement. Default tuning uses coplanarity, signed missing
energy and signed missing mass squared. Photon opening angle, pTmiss and missing
mass eγ are held-out observables. A deterministic independent half of the cached
MC provides an additional sampling check. Both halves compare to the same DATA;
this is not an independent DATA split. No other period is used as a surrogate
validation sample for a different detector/software configuration.

The fit score compares areas normalized within each projection to focus on
shape. It retains the signal/background ratio. Displayed predictions, all joint
cell comparisons use one candidate signal normalization per topology, fitted in the fully selected sample and held fixed across observables. Raw selection ratios are reported separately. The score
uses frozen nominal MC uncertainty and bin masks so increased sampling noise
cannot improve the objective simply by inflating its denominator. Every training
projection must be populated. This score is descriptive, not a formal independent
chi-square or a source of fit-confidence intervals.

The search is a candidate test, not a proof of a global optimum. Scan-boundary
parameters, degenerate template responses, failures on held-out observables and
remaining tails all require further investigation. Do not apply these parameters
to production without checking the outputs. The trial model is independent
Gaussian broadening and momentum bias: it cannot diagnose arbitrary correlated
errors, radiation, background-transfer mismodelling or reconstruction losses.

## Outputs

Open results/REPORT.html. Transfer only the results directory for review; keep
the larger event cache and evaluation files on the farm. Click one period, then CD_FT for the before/after overlays.
Each period directory contains:

- single_parameter_effects.png: which individual particle variations help
  training and held-out distributions.
- single_parameter_ranking.csv and single_parameter_projection_scores.csv:
  individual sensitivity and per-observable results.
- response_similarity.csv: similar finite template perturbations (cosine near
  +1 indicates ambiguous attribution, not a fitted parameter covariance).
- candidate_parameters.json/.csv: independent candidate values, units and
  scan-boundary flags.
- final_scores.csv: nominal/candidate training, independent-MC and held-out scores.
- selection_yields.csv: actual predicted changes after all selections.
- final_shapes.csv: nominal/candidate comparisons in inclusive and joint cells.
- nominal_replay.csv: sampled nominal yield versus the full reference with its
  sampling error. A gross mismatch stops the fit.

The root cache_manifest.json records all cached MC checksums/probabilities;
reconstruction_audit.csv records raw four-vector recomputation offsets. Backend
and per-evaluation logs are retained. An interrupted fit can reuse its completed
cache. Existing outputs are never overwritten.

## Computational design

The first stage takes a deterministic uniform sample BEFORE any analysis cuts,
including events which can migrate into selection. Default expected sample is
200,000 reconstructed events per period and MC category (signal/leakage); it is
not a cap on selected events. Every retained event is weighted by its inclusion
probability. Sampling limits are precision/runtime controls, not physics cuts.
Increase --expected-cache-events for sparse samples.

ROOT inputs are read once for this cache. Every later scan uses compact binary
records, common random numbers and native C++ calculations. Gaussian random draws
and nominal selection contributions are computed once per event and reused
across trials. Broad one-parameter/coordinate scans accumulate inclusive
histograms only; the final nominal/candidate pass also fills all joint cells.
No extra theory or cross-section model evaluations are required.

Independent period fits run concurrently but share ONE hard worker budget, ten
by default. Failed workers cancel queued/active work. DATA histograms are loaded
once rather than repeatedly decompressed inside the objective. Matplotlib runs
serially after fitting because it is not thread safe.

To repeat fitting without ROOT reads, keep the first run's cache and choose a
new output directory:

```sh
python3 study_particle_response.py --reference-study "$PWD/test_runs/exclusivity_resolution_all_periods_joint_001" --reuse-cache "$PWD/test_runs/particle_response_all_periods_001" --workers 10 --output "$PWD/test_runs/particle_response_all_periods_002"
```

--coordinate-rounds and --refinement-rounds control search cost. Defaults are one
each. --fit-topologies can include CD_FD/FD_FD; that adds FD photon/proton
parameters, but those fits are much more sensitive to the pi0 model. Default
zero values for those *unfitted* regions do not assert correct FD resolutions.

## Scratch storage

Large event caches and native evaluation binaries now live under
`/scratch/thayward/dvcs_particle_response` by default. Every run gets a unique
subdirectory, allowing independent runs without collisions. `--scratch-root`
overrides the parent. `scratch_location.json` records the actual directory.
Reports, CSVs and the compact results directory stay under --output.

Intermediate evaluation binaries are unmapped and deleted immediately after
their scores/sensitivities are recorded; final binaries are removed after plots
and summaries are produced. Event caches are retained for reuse. Their manifest
uses absolute paths, so --reuse-cache still points to the prior study output
even though its binary files are on scratch. Keep that scratch cache available
for later reuse. Pre-change relative cache manifests remain readable.

## Finish reports after a reporting failure

If all final evaluations finished but report generation failed, retain the study
and scratch directory. `--finish-reports` reads the saved final histograms and
candidate parameters, without ROOT reads or refitting:

```sh
python3 study_particle_response.py --finish-reports --reference-study "$PWD/test_runs/exclusivity_resolution_all_periods_joint_001" --output "$PWD/test_runs/particle_response_all_periods_001"
```

It validates the reference identity, final parameter agreement and binary sizes.
Completed report markers permit retries after some periods were already
rendered and their trial binaries deleted. Final binaries are removed only
after the corresponding period report is complete. Existing compact results
can be refreshed safely. CSV normalization values are explicitly converted
to numbers before reporting selected yields.

## Joint-topology model comparison (current default topology selection)

All three topologies now participate in each independent period/polarity fit.
The electron FD response is shared across all three; CD-proton response across
CD_FT/CD_FD; FD-photon response across CD_FD/FD_FD. FD protons and FT photons
have their own detector-specific responses. These same laws are applied to
signal and leakage MC: fractional momentum width gives sigma_p=a*p, with
constant angular widths and optional fractional momentum bias per detector
region. The identical law is used at identical particle kinematics in every
topology. This initial functional model is not a fine-grained momentum/angle
calibration; joint-cell validation can reveal a need for more dependence.

--compare-models runs broadening_only (all momentum biases exactly zero) and
with_bias (biases fitted). The second model reuses the first model's event cache
and includes its period-specific candidate as a starting point. Both searches
retain nominal and one-parameter starting candidates, then perform bounded
coordinate search/refinement. It is still not a proof of a global optimum.
Native worker budget stays ten TOTAL; the two model runs execute sequentially.
The all-period old cache can also seed the first model, avoiding all ROOT reads.

Scan limits are broadened: electron momentum sigma up to 3%, CD-proton up to
12%, FD-proton up to 6%, photon energy up to 6%; angular sigmas up to 12 mrad
(theta) and 60 mrad (phi); fitted momentum bias up to ±1.5%. These are exploratory
ADDITIONAL-error ranges, not measured detector resolutions. Boundary hits are
recorded and still require scrutiny.

One pi0 transfer multiplier per period is shared by all topologies and profiled
using training observables only. Its constraint combines the saved control
counting error with --background-transfer-uncertainty (default 0.10). That extra
10% is an explicit exploratory transfer uncertainty, NOT a measured systematic.
Set it to 0 for counting-only constraints, or rerun other assumptions from cache.
The grid considers 1, 1±sigma and 1±2sigma, with a squared standardized penalty.
A single multiplier cannot repair arbitrary background spectral mismodelling.
The same chosen multiplier is used for independent-MC and held-out validation.
The unsmeared nominal displayed model remains at its original background scale.

candidate_parameters.json records the multiplier, constraint width and whether
it hit a nuisance boundary. selection_yields.csv now separates signal and
background yields; a nuisance normalization change is not an acceptance change.
No per-variable/per-topology background normalization or refitting is performed.

Recommended command:

```sh
python3 study_particle_response.py --reference-study "$PWD/test_runs/exclusivity_resolution_all_periods_joint_001" --reuse-cache "$PWD/test_runs/particle_response_all_periods_001" --compare-models --workers 10 --scratch-root /scratch/thayward/dvcs_particle_response --output "$PWD/test_runs/particle_response_joint_topologies_001"
```

Open results/REPORT.html and compare the two models. model_comparison.csv lists
per-period training/held-out scores, nuisance multipliers, boundary counts and
separate signal/mixture yield changes. Transfer only results for assessment;
large event/trial files stay on scratch. Existing studies remain unchanged.
For report-only recovery of this paired run, add --compare-models --finish-reports
and use the same output/reference paths. The scalar response fit remains an
effective model; production acceptance and Frad are not modified.


## Consistent two-photon control response

Both individual photons are recovered from the saved composite three-momentum,
Trento azimuths, electron opening angles and measured gamma-gamma mass. Each is
varied with its own detector response and independent random numbers; the same
laws apply to signal and leakage photons. The composite is rebuilt with the
producer's nominal pion-mass convention. Unresolved pairs are retained unchanged;
the run stops if they exceed 1% of nominal selected control events. Inspect
`control_photon_reconstruction.csv` and each period's `control_transfer.csv`.

For each trial, alpha = alpha_reference * C_nominal/C_trial, where C is the
selected two-photon MC control count. A shared constrained transfer multiplier
is applied in addition. For each topology beta = (D-alpha*B_trial)/S_trial,
using only the common fully selected sample. Beta then stays fixed in every
omitted-cut projection and joint cell; there is no per-plot normalization fit.
Invalid negative signal residuals reject the trial. Counting diagnostics use
independent MC halves consistently with the shape checks.

Rate plots show this profiled prediction; `_shapes.png` plots normalize each
curve to unit area. `selection_yields.csv` includes `raw_signal_selection_ratio`:
with fixed generated sample this is the relative reconstructed selection
acceptance, before any fitted signal scale. A ratio below one means fewer
selected MC events and hence a larger acceptance-corrected data yield.

The producer already requires 0.11 < m(gamma gamma) < 0.16 GeV before saving
these control trees. The study reapplies saved analysis-stage cuts, records
migration out of that upstream mass window, and does not add another mass cut.
It cannot recover incoming candidates discarded upstream. These remain
conditional response tests; a wider upstream control skim is needed to establish
the full mass-window response before production adoption. Production analysis
and imported Frad/Fbin remain unchanged.


## Larger-cache wide scan and production cut redetermination

Use `--expected-cache-events 800000` for four times the original target per
period/MC sample (or every event in smaller samples). It is a pre-cut expected
count; MC halves each contain approximately half. Standard error should decrease
by roughly two where the target grows fourfold. Incompatible smaller reusable
caches now cause an explicit error. Omit `--reuse-cache` when increasing the target.
The paired second model automatically reuses the first model's larger cache.

`--scan-range wide` keeps the old grid and adds momentum biases to +/-10%,
electron fractional momentum broadening to 10%, CD proton broadening to 24%,
other momentum broadening to 16%, theta broadening to 40 mrad and phi broadening
to 160 mrad. These are exploratory search limits, not plausible calibration
errors or priors. Use `--coordinate-rounds 2 --refinement-rounds 3` to revisit
correlated coordinates and narrow local brackets; held-out scores never choose
parameters. Refinement retains newly inserted grid points so successive passes
actually narrow the search. A bounded coordinate scan still does not prove a
global optimum or determine parameter uncertainties.

The Mx2(epgamma) display uses robust core widths rather than extreme tails to
set its visible range. All histogram bins and tails remain in scores/tables.
`_shapes.png` normalizes each curve independently; ordinary plots retain a
single selected-sample normalization per topology. They can look similar by
construction. Neither alone validates the event rate.

`--redetermine-cuts` runs `study_response_cut_redetermination.py` after fitting.
It reads each DATA channel once at full statistics and exports pre-exclusivity
MC arrays from the shared cache, preserving aligned events across responses.
Large arrays stay in the unique scratch folder. It calls the existing production
`fit_single_template`, `run_dvcs_iterative_cuts`, and `run_pi0_iterative_cuts`
functions with their current default settings, including the direct-pi0 priors,
iterative order, stopping rule and separate DATA/MC 90/95/98% windows. Display-only
xF variables cannot affect any fitted cut and are not exported. Production cut
fitting uses its original minimal preselection; the subsequent relative-acceptance
counts additionally apply the saved full analysis global selection. No fixed old
exclusivity cut is imposed before determining new cuts. DATA themselves are not
morphed, but their fitted signal templates and boundaries can change.

A fresh nominal fit and both candidate fits are retained independently. Each
candidate is compared both with saved cuts and with the fresh nominal refit.
`relative_acceptance.csv` reports reconstructed-count ratios by period/topology
and |t|, with correlated MC counting errors. The generated denominator cancels
for a fixed generated sample and fixed reporting bin. Inverse acceptance alone
is not the cross-section change: DATA selection and background subtraction must
also be updated. Parameter/cut-fit uncertainties are not included in that error.
The upstream two-photon mass-window limitation remains unchanged.

Two-photon control figures now show DATA, nominal MC and candidate MC before
exclusivity cuts, with the newly fitted DATA/MC boundaries overlaid. They are
linked from `results/cut_redetermination/REPORT.html`; the same folder contains
per-topology production cut-development plots and trial cut JSONs. Nothing is
installed into production output/jsons. `DVCS_EXCLUSIVITY_ARRAY_ONLY=1` lets this
array-only caller import the production fitter without starting PyROOT.
