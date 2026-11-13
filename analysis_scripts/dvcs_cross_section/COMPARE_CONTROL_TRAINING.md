# Controlled comparison: epgamma-only versus joint shape training

```sh
python3 test_compare_control_training.py && python3 study_compare_control_training.py --joint-run "$PWD/test_runs/particle_response_joint_control_001" --workers 10 --output "$PWD/test_runs/compare_control_training_001"
```

Reuses the completed joint branch after checking its evaluator source hashes,
cache, local-refinement configuration and recorded starting candidate. Fits only
the missing epgamma-only branch using the same pre-joint starting candidate,
parameters, bounds, background normalization/nuisance prescription, cache size,
seed, train/reserved variables, topologies, and coordinate/refinement rounds.
Only the control contribution to the training objective is zero instead of 0.5.
Both branches still morph two-photon AAOGEN to calculate the common control
normalization. A missing/changed start or evaluator aborts instead of silently
presenting an unpaired comparison. Branches must both have selected that recorded
warm start over nominal. Equal schedules give equal maximum budgets; boundary
proposals can deduplicate differently, so exact unique evaluation counts can differ.

Every fit period stays independent. Existing MC caches and two-photon DATA exports
are verified/reused. No DATA ROOT read, new MC cache, extra broad parameter scan,
production cut update or control reweighting occurs. The ten-worker budget covers
all periods/samples of the new fit. A final common validation evaluation batches
nominal and both fitted candidates in one pass per sample, sequentially by period.
It reuses the newly compiled evaluator rather than compiling another backend.

Open REPORT.html. For each period:

* scores.png: side-by-side scores for epgamma-only and joint candidates, separately
  for epgamma and epgammagamma, separating training observables, other MC half and
  reserved observables. Lower is better.
* epgammagamma_TOPOLOGY.png: experimental two-photon DATA versus AAOGEN MC under
  BOTH fitted responses. These shapes are validation for the epgamma-only fit.
* epgamma_TOPOLOGY.png: DATA versus DVCS+pi0 leakage predictions from both fits.
* projection_comparison.csv: every channel/topology/observable and MC half.
* score_comparison.csv: channel means in each score category.

Both candidates use identical validation definitions and their frozen fitted
background multipliers. No new nuisance or response parameter is fitted on the
validation plots. MC half 1 and reserved observables never select parameters;
experimental DATA histograms are shared, so do not call this independent DATA
validation. Overlay plots combine MC halves for statistical precision; the score
CSVs distinguish the validation half. Variables/projections share events and
scores are descriptive, not independent chi-square p-values. Shape normalization
does not validate absolute efficiencies. The two-photon control may contain other
backgrounds; the comparison does not by itself establish a unique detector cause.

The feature is also accessible directly through study_particle_response.py:
--joint-control-study supplies the control shapes and checked warm start;
--joint-control-weight 0 trains on epgamma alone but still evaluates/reports the
two-photon channel, while 0.5 retains the current equal-channel objective.
