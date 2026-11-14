# Selection-first particle-response pipeline

This is an isolated diagnostic pipeline. It implements the agreed order of operations and leaves the production cross-section code/results untouched.

From `analysis_scripts/dvcs_cross_section`, after committing/transferring the update:

```sh
git pull && python3 test_selected_particle_response.py --native && python3 study_selected_particle_response.py --joint-run "$PWD/test_runs/particle_response_joint_control_001" --workers 10 --rounds 1 --scratch-root /scratch/thayward/dvcs_particle_response --output "$PWD/test_runs/selected_particle_response_001"
```

The joint run's original reference directory, DATA ROOT files and scratch MC caches must still exist. Report-only Desktop copies do not contain those inputs. One local-refinement round is the default for speed; increase `--rounds` only for a fresh run if refinement remains justified. Every period/polarity remains independent. The starting response is the saved joint fit; nominal and starting/fitted curves are all retained for comparison.

## Order of operations

1. Read DATA once and export full pre-exclusivity reconstructed MC from the checked existing caches. Preserve events that fail global/exclusivity selection; store their flags, do not discard them permanently. Only the main-analysis-supported groups are used.
2. Run the actual legacy template cut-development functions on nominal MC: three DVCS topologies and four unordered two-photon groups. Obtain initial separate DATA/MC 95% windows. For DVCS templates, pure FD/FD or FT/FT controls provide their corresponding calibration. Mixed-photon controls use the CD-FD optimizer policy on their own distributions and derive their own numerical windows.
3. Freeze the initial **DATA** windows for response fitting. Fit the bias+smearing model jointly to epgamma (DVCS signal + pi0 leakage) and epgammagamma (AAOGEN control). Every trial is applied to the entire pre-cut MC cache, then globals and fixed DATA windows are evaluated. Thus events may migrate into or out of the reference selection. DATA coordinates never change and cuts are not refitted during trial scoring.
4. Use N-minus-one distributions: omit the displayed variable's exclusivity cut while retaining the others. Training uses Delta_phi, pTmiss, Emiss2 and Mx2; control training also includes the pi0 mass shape. Opening angle and Mx2_2, plus the other MC half, are validation only. All DATA are shared, so this is not independent DATA validation. Projection support and reference variance are fixed from DATA/nominal MC, not selected by a response trial. Native proposals are batched in a single cache pass. Training scores alone select parameters.
5. After response fitting, run the legacy template cut procedure again with the modified signal, leakage and control MC. Produce separate DATA/MC 95% windows, then recount with these actual resulting selections. Record relative reconstructed acceptance, coarse pi0 transfer and changes to each window. This is one explicit iteration; inspect window changes before requesting another.

The native runtime uses one shared ten-process budget across periods and samples. Cut fits also cap their period concurrency; BLAS/OpenMP threads are restricted to one. Saved initial exports and completed fits can be reused. Large arrays, trial logs and export files stay under the requested scratch root. Intermediate trial histograms are deleted; only small final histograms are retained for reports. Reports are produced once, under this output tree.

## Supported control groups and selection

Control groups are FD proton + FD/FD photons, and CD proton + FD/FD, FD/FT or FT/FT photons. Group membership is invariant under exchanging photon roles. FD-proton groups involving FT photons are omitted.

For controls, the common global comparison domain requires both possible legacy detector-label assignments to pass, using the existing composite pi0 coordinates. This conservative intersection is inherited from the unordered diagnostic; it is not a new physical two-photon acceptance prescription. An upstream tree-processing rejection cannot be recovered by this analysis. Reconstructed MC caches contain events before analysis-global/exclusivity selection, not unreconstructed generator events.

## Mixture constraint

Use one period-wide selected two-photon normalization, the saved control-to-leakage MC exposure ratio and sampling probabilities. Recompute the control MC selected yield after every response trial. Signal normalization is then determined per DVCS topology inside the fixed selected DATA reference. A common background multiplier is profiled on training only, with the original exploratory transfer uncertainty plus control counting statistics. No independent topology background multipliers are introduced. An invalid nominal signal residual stops the fit, rather than clipping the background or inventing a signal scale.

## Pi0 mass definition

Other reconstructed observables already preserve stored offsets by adding only the change induced by particle modification. The isolated pipeline uses the analogous mass definition:

```
reported modified mass = saved mass + rebuilt modified mass - rebuilt original mass
```

For a zero response it preserves the stored mass exactly. For a nonzero response it retains the change of the actual two-photon invariant mass. This removes a possible numerical recovery offset from the mass objective. Both native scoring and exported modified masses use this definition. Production response headers are unchanged. Unresolved photon pairs retain stored composite quantities/mass; their count/coverage must be checked and their response cannot be claimed to be constrained by those unresolved events. The original upstream pi0 mass window still truncates the available sample.

## Outputs and interpretation

Open `REPORT.html`.

- `plots/PERIOD/scores.png`: training, other-MC-half and reserved scores, separated by channel; nominal versus previous response versus refitted response.
- `plots/PERIOD/epgamma_*.png` and `epgammagamma_*.png`: N-minus-one comparisons using the frozen initial DATA selection. Signal/leakage mixture and pi0 mass are shown.
- `postfit_common_selection/plots/PERIOD`: fitted response versus DATA with final DATA windows applied to all samples, omitting the displayed cut. Parameters and the fitted nuisance are held fixed; component normalization uses the new selected counts. These comparison windows differ from the separate MC windows used for the actual acceptance recount.
- `initial_cuts/PERIOD` and `final_cuts/PERIOD`: actual production-method template/cut-development plots and 90/95/98% window artifacts. The selected initial/final 95% windows are also saved together at the output root.
- `initial_selection_counts.csv`, `final_selection_counts.csv`: global/selected raw counts and retention using each sample's actual DATA/MC windows. These are not independent signal efficiencies in DATA.
- `relative_acceptance.csv`: candidate versus initial nominal **reconstructed-count ratio**, by DVCS topology and |t|. The unchanged generated denominator cancels. Error includes paired MC counting only, not response-fit or cut-fit uncertainty.
- `coarse_pi0_transfer.csv`: initial/final |t|-binned control normalization and leakage subtraction. This is deliberately labelled a **raw-count coarse diagnostic**: no production event current/proton weights or xB/Q2/phi binning. Unsupported control bins are flagged; missing controls never imply an invented zero background. The table uses the unprofiled control normalization, not the fitted exploratory nuisance.
- `window_changes.csv`: old/new boundaries, which can move even for unchanged DATA because the templates changed.
- `response/PERIOD/objective_support.csv`, `nominal_objective_projections.csv`, `fit_steps.csv`, `projection_scores.csv`, `evaluation_history.csv`: objective support, optimization history and individual checks.

No quoted inverse-acceptance factor alone is a cross-section prediction. Production adoption would need the same control definition, cuts and current calibration in the complete extraction, including weighted pi0 transfer. The pipeline does not silently adopt the fitted response or update production results.

## Recovery

Repeat the study command with `--resume` appended. Checked initial/candidate event exports, completed cut stages and completed response fits are reused. Candidate JSON alone is not a completion marker: its final histogram/score checksums must match. Source, fitter, input/settings and executable provenance must also match. A failed partial fit restarts that period's fit in a new scratch attempt; completed periods remain intact. Changed code/settings require a new output directory, rather than silently combining incompatible stages.

## Validation

`test_selected_particle_response.py --native` covers joint mass training, training/validation separation, fixed support, full driver ordering, checked resume, native MC migration into/out of fixed cuts, unordered groups, zero-mass preservation, modified-mass calibration and real Python/native evaluator integration. A separate local integration exercise ran the actual legacy fitter successfully for all three DVCS topologies and all four control groups on synthetic arrays. Farm numerical results require the real inputs.

The Mac's ROOT/Cling installation has a known DATA file-I/O problem. The cached-MC native tests do not need ROOT file reads. If desired, run `test_unordered_control_response.py --root-smoke --root-io` first on the farm to check the actual DATA ROOT reader path, as in the preceding unordered study.
