# Joint epgamma / epgammagamma response fit

Photon locations confirmed with the user: detector_gamma1 and detector_gamma2 are
independent labels, 1=FD, 0=FT. The control topology labels proton and photon 1;
photon 2 is not assigned photon 1's detector. The existing recovery cache stores
both tags. Electron/proton parameters are shared with epgamma in the same period;
each photon gets its own FD/FT response before the composite pi0 is rebuilt.
Unexpected detector tags fail rather than being silently classified as FT.

New opt-in --joint-control-study adds two-photon shape constraints. Default legacy
runs remain epgamma-only. Fits remain independent by running period, use shared
particle/detector parameters across topologies, and retain one control-constrained
background multiplier per period. No independent topology background multipliers,
control weights, new kinematic response parameters or production changes are added.

The objective is half the mean epgamma training-projection score plus half the
mean two-photon training-projection score, plus the existing background nuisance
penalty once. These are descriptive shape scores, not a joint likelihood or a
chi-square probability: variables share events. Coplanarity, missing pT, missing
energy and missing mass squared train the fit. Opening angle and missing mass
of e+neutral are reserved. Every required nominal projection must have enough MC.
The reserved MC half never selects fit parameters; DATA histograms are shared
between those comparisons, so this is not independent experimental-DATA validation.

Control shape histograms use full saved globals, the same saved DATA windows on
DATA and MC, and omit the displayed exclusivity variable. Their transfer/count
normalization keeps the existing separate MC control windows; normalization and
shape selections are explicitly different. Shapes are area-normalized; control
shape training does not independently validate absolute control efficiency.
No extra histogram morph is fitted. Unresolved photon recovery retains original
records as in existing response fitting; existing unresolved-fraction guard remains.

Uses verified pi0_control_transfer DATA exports without rereading DATA ROOT.
Candidate control histograms are accumulated during the existing MC control pass,
so no extra control MC scan is needed per candidate. Uniform MC caches are reused;
large caches/trials stay on scratch. The single shared budget caps workers at ten
across all periods and samples. Joint shapes are inclusive, not additional large
kinematic-grid histograms.

Example using the existing fresh-study paths:

```sh
python3 test_joint_control_fit.py && python3 study_particle_response.py --reference-study "$PWD/test_runs/exclusivity_resolution_fresh_001" --joint-control-study "$PWD/test_runs/pi0_control_transfer_001" --reuse-cache "$PWD/test_runs/particle_response_fresh_001" --expected-cache-events 800000 --bias-mode fit --scan-range wide --local-refinement --coordinate-rounds 1 --refinement-rounds 1 --workers 10 --scratch-root /scratch/thayward/dvcs_particle_response --output "$PWD/test_runs/particle_response_joint_control_001"
```

The reuse path accepts either the direct cache folder or its legacy with_bias
parent. The joint starting candidates are hash-checked against the completed
add-pT training study referenced by the control exports. --local-refinement skips
the broad initial single-parameter scan, compares nominal with the prior candidate,
and explores adjacent grid levels followed by midpoint refinement. It is a local
search, not a claim of a global optimum. Omit it to request a new broad scan.

Open the main REPORT.html then period reports. Existing epgamma plots remain;
control_plots/PERIOD/TOPOLOGY.png shows nominal and jointly fitted two-photon MC.
final_scores.csv separates epgamma_* and epgammagamma_* components. Projection
CSVs include a channel field. candidate_parameters.json records joint_control_fit,
training/reserved variables and parameter values. The old single-parameter ranking
is empty in local mode because that scan was deliberately omitted.

No production cuts/acceptances are updated by this command. Do not pass
--redetermine-cuts until the joint-fit shapes and component scores are reviewed.
