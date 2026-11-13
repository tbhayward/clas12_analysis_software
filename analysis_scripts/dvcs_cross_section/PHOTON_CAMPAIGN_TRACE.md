# Photon ordering: campaign-to-cross-section trace and unordered control study

Inspected processing branch: `dvcs_processing_branch`, commit `5a8b14f1a8a57f638079b99a0ea409fbd35a3c18`.
Inspected analysis branch: `development`, commit `940ce6e7cf893c6d52c1635a393fd1356a263091`.
Branch sources were read with `git show BRANCH:PATH`; the working checkout remained on development throughout. This examines both branches without disrupting the checkout or mixing their files.

## 1. The actual new campaign entry points

The relevant DATA path is not merely the older top-level script:

- `processing_scripts/rebuild_pilot/data_pilot_worker.py`: `SCRIPTS['pi0_data']` selects `processing_exclusive_pi0_stream.groovy`. It runs that script with `processing_classes/dist/processing_classes.jar`, then `convert_data_stream.cpp` writes `PhysicsEvents` using `reconstructed_data_schema.h`.
- `processing_scripts/rebuild_pilot/campaign_worker.py`: `MODE['pi0_mc']` selects the rebuild-pilot `processing_mc_rec_exclusive_pi0.groovy`, using the same JAR location. It streams paired generated/reconstructed rows through `convert_paired_stream.cpp`.
- `processing_scripts/rebuild_pilot/extract_reconstructed_mc.py`: copies valid `rec_*` fields into a companion `PhysicsEvents` tree, stripping the prefix. It does not swap photon labels, infer detector tags or reorder photons.

The DATA streaming script preserves the legacy ordered pair loop (lines 283–305): for two photons it visits `(0,1)` and then `(1,0)`. The last accepted permutation overwrites `detector_gamma1/2`, electron–photon opening angles and photon azimuths together. It explicitly assigns control `detector2 = detector_gamma1` (line 324).

The rebuilt MC script uses `choosePhotonPair` (lines 78–99): `j=i+1`, choose the pair closest to the pi0 mass. `pairedKinematics` gets both tags and both sets of angles from that same pair analyzer. Thus two-photon events use `(0,1)` instead of the DATA script's final `(1,0)`.

The branches have the same names and detector codes, but **photon 1 does not identify the same energy-ranked member**. The measured nearly complete role reversal agrees with this source behavior. No generator-detector assignment is needed to explain it.

## 2. Java classes and compiled artifacts

`ThreeParticles.java`:

- The constructor's `generated` argument skips reconstructed detector/fiducial assignment for truth.
- Reconstructed detector codes are FT=0, FD=1 and CD=2, obtained from reconstructed-bank detector-status tests.
- For the `(gamma,gamma)` analyzer, detector1/2 and phi1/2 refer to the two supplied photon occurrences.
- For the `(proton,pi0)` analyzer, particle 2 is synthetic PID 111; its detector field is not a meaningful photon detector.

`analysis_fitter.java`:

- Applies 2 GeV to the first accepted photon and 0.4 GeV to subsequent photons.
- Builds composite pi0 candidates with nested photon loops through `parent_hadron_creation.pi0_check`.
- `parent_hadron_creation` constructs the pi0 momentum as the photon-vector sum in a mass window.

I inspected committed bytecode with `javap`. The processing-branch JAR entries match the corresponding committed build classes byte-for-byte:

- `analyzers/ThreeParticles.class`: SHA256 `37dd69dc4df33670f8338acd0b9d4a7d01d55dc2dc500b63e79c33b91c834667`.
- `extended_kinematic_fitters/analysis_fitter.class`: SHA256 `22426de1d5c22137b4a08b320cb088dc28bbdb48ad95c87dc65563e4a6f438f2`.

This validates the checked-in compiled artifacts, not the identity of the historical JAR actually used on every farm node. The campaign `.ok` marker records input, period, mode and event counts but does not pin the JAR or script checksum. Existing tree-level ordering results provide independent evidence of the role difference.

## 3. What happens to generated detector fields

For pi0 truth, `pairedKinematics(..., generated=true)` writes NaN detector/fiducial fields explicitly. Generated DVCS also skips detector assignment; any default truth detector values are not detector observations.

In development `total_counts.cpp`:

- `BranchBinder` selects `gen_*` for generated MC and `rec_*` for reconstructed MC.
- `ready_for_generated_matching` (around line 1584) requires kinematics, not detector IDs.
- The event loop (around lines 2395–2438) applies topology, global and exclusivity selection only to DATA/reconstructed MC. Generated events use `passes_generated_phase_space_dispatch` instead.
- `acceptance.cpp` (around lines 723–767) sums reconstructed DVCS topologies and divides by the generated yield.

Therefore **unset truth detector IDs do not sort generated events into FT/FD acceptance bins in this extraction**. They do not cause the observed photon-order reversal.

## 4. Where development imports the reconstructed convention

`load_trees.cpp` opens paired `MCEvents` for generation and `.rec_only.root` `PhysicsEvents` for reconstructed MC.

The direct-pi0 reconstructed `detector2` is replaced by `detector_gamma1` (or `rec_detector_gamma1`) in:

- `total_counts.cpp`, `BranchBinder::physical` (around lines 1405–1413).
- `current_dependence.cpp`, regional MC binding (around lines 1363–1370).
- `eppi0_normalization.cpp`, branch binding (around lines 1041–1048).
- `plot_exclusivity_data_dvcs_pi0_mc.py`, branch-name resolution (around lines 760–763).
- The response/control workers use each individual photon tag for its own response but retain the historical photon-1 control topology.

Using a real photon tag instead of the synthetic pi0 detector is necessary. However, choosing photon 1 preserved the reconstructed DATA/MC role difference. This is distinct from the earlier repair of an invalid composite `rec_detector2` branch.

## 5. Can this affect the extracted cross section?

Yes, through the selected pi0 control/transfer, even if it does not change the primary DVCS acceptance directly.

`pi0_contamination.cpp`, `compute_contamination` (around line 483), computes per kinematic row and period:

```
c = N_misidentified_pi0_MC * N_two_photon_DATA
    / (N_two_photon_reconstructed_MC * N_DVCS_DATA)
```

All three reconstructed topology yields are summed. Generated pi0 counts do not appear in this formula. `pi0_corrected_counts.cpp` subtracts the resulting contamination from normalized DVCS DATA; `cross_sections.cpp` uses those corrected counts with the DVCS acceptance and other corrections.

Pure relabelling would cancel in a complete sum with identical cuts and weights. That condition is not met automatically:

- `total_counts.cpp` selects/rejects supported topologies before summing; FD-proton/FT-labelled events are outside its three-group list.
- It applies topology-dependent exclusivity windows.
- `global_cuts.cpp` applies detector-dependent auxiliary angle, momentum and sector conditions. For legacy direct-pi0 control these operate on the stored composite coordinates with the photon-1 detector tag.
- Regional current-response corrections use the FT/FD-sector assignment before summing.

Thus the numerator/denominator of the pi0 transfer can change when the control definition is harmonized. The audit's 21–29% DATA CD–FT membership exchange is not a 21–29% cross-section error: migration in/out largely balances, and that audit did not redo selections/current weights.

For fixed DVCS acceptance and other corrections, a change from contamination c to c' changes the signal cross section by `(1-c')/(1-c)`. Neither the size nor sign of c' is established by the grouping audit.

### Required production decision

The control grouping/selection should be made explicitly consistent before final precision extraction. This is **not** a reason to assign detector tags to generated tracks, flip tags without coordinates, regenerate detector simulation, or declare all existing DVCS acceptance values wrong.

The new unordered diagnostic is a comparison tool, not a drop-in production correction. A validated production adoption would require the common control definition in the exclusivity fitter, count loader, current calibration and pi0-transfer modules together, followed by recounting/recomputing the subtraction. Keep the current extraction as a baseline until the difference is quantified. No production physics code is changed by this update.

## 6. Two additional processing assumptions to audit

1. With multiple photon-pair candidates, DATA retains late pair metadata but analyzes composite PID 111 occurrence zero. Rebuilt MC selects the closest-mass photon pair while also analyzing composite occurrence zero. Those selection rules do not guarantee pair/composite agreement. Recovery flags some inconsistencies but does not prove every surviving pair is the intended one.
2. `ThreeParticles.getIndex` resolves PID occurrence in the original REC bank. The fitter can reject candidates before constructing `PhysicsEvent`, so filtered occurrence need not equal original-bank occurrence. The fitter exposes original-row tracking, but the inspected analyzer uses the raw occurrence lookup. This is an additional possible tag/fiducial association issue, not a demonstrated explanation of the measured DATA/MC ordering reversal. Quantifying it fully requires original HIPO row identity. The new report includes recovered photon angle distributions by tag and role as a useful, limited check.

For a future processing repair, use one shared pair selector for DATA and reconstructed MC, build the analyzed composite from that selected pair, carry its selected original REC indices, and save both photon four-vectors directly. Pin the script/JAR hashes in completion metadata. These are processing recommendations, not changes made to the historical files here.

## 7. Implemented unordered control comparison

`study_unordered_control_response.py` defines FD/CD proton × {FD–FD, FD–FT, FT–FT}: six disjoint groups. The pair label is invariant under photon exchange. Both photon roles retain the same detector response law.

The study reads two-photon DATA ROOT files once and uses the existing pi0 MC caches for nominal and frozen morphed MC. Fresh exports avoid inheriting the old global-only export's grouping-dependent losses. A checked isolated worker is compiled under this study output; production headers and the saved fit are untouched. A ten-worker total budget and scratch exports are retained.

The common diagnostic global selection requires **both** possible legacy detector-label assignments to pass, using the existing composite pi0 kinematics. This conservative intersection is explicitly order independent; it is not a newly derived two-photon acceptance prescription. It retains all six groups if allowed by the configured cuts. The original tree-processing losses cannot be recovered.

Compare:

- `global_only`: no exclusivity windows, to see what the frozen particle response does before truncation.
- `common_windows_Nminus1`: the same DATA-window intersection within each proton region, leaving out the displayed exclusivity variable. CD intersects the two saved CD window sets; FD uses its available FD–FD set. These comparison windows are not newly optimized 95% cuts.

Each panel overlays unit-area DATA, nominal AAOGEN and frozen morphed AAOGEN. Pi0 mass is included. Scores retain both MC halves, reference support/variance and insufficient-statistics flags. DATA are shared; this is not independent DATA validation. Counts remain raw diagnostic counts, not current-corrected contamination estimates.

Run from `analysis_scripts/dvcs_cross_section` after committing/transferring the update:

```sh
git pull && python3 test_unordered_control_response.py --root-smoke --root-io && python3 study_unordered_control_response.py --joint-run "$PWD/test_runs/particle_response_joint_control_001" --workers 10 --scratch-root /scratch/thayward/dvcs_particle_response --output "$PWD/test_runs/unordered_control_response_001"
```

Open `REPORT.html`; start with `plots/PERIOD/global_only/CDp_FT_FT.png`, `CDp_FD_FT.png` and `CDp_FD_FD.png`. CSVs provide scores, selected counts/unresolved coverage and photon angle quantiles by tag/role. Transfer the output report directory; retain its scratch arrays on the farm.

## Validation and local limitation

Local native cached-MC exporter tests pass, including all six groups and both photon permutations, as do the complete Python plot/score/window tests. The separate `--root-io` test creates a DATA ROOT fixture and exercises the actual branch-binding/recovery path. Its local execution is blocked by this Mac's ROOT/Cling module cache referencing missing macOS SDK headers (`experimental/__memory`, `_types.modulemap`), which crashes ROOT file creation. The farm command explicitly includes that DATA I/O test before starting the study. No ROOT installation or existing module cache was modified. Full campaign numerical results require the farm inputs.
