# Exclusive e' p pi- processing campaign

Dedicated data-only processing for the exclusive_pion_processing_branch. The default Groovy processor, Java fitter/analyzer, converter, and launchers are unchanged. No git commits are made by these scripts.

## Test launch on the farm

Transfer this directory to `processing_scripts/exclusive_pion_campaign/` in the farm checkout. All commands below run from that checkout's root. There is no SSH coordinator or startup barrier: each command starts workers only on the host where you are logged in.

First, prepare the shared test campaign **once**, on either node:

```sh
csh -f processing_scripts/exclusive_pion_campaign/environment.csh python3 processing_scripts/exclusive_pion_campaign/campaign.py prepare --repo "$PWD" --campaign /volatile/clas12/thayward/eppi-/test_20261122 --mode test
```

Wait for preparation to finish. Then log into **ifarm2401** and run:

```sh
tmux new -s eppi-minus-test 'csh -f processing_scripts/exclusive_pion_campaign/run_test.csh /volatile/clas12/thayward/eppi-/test_20261122'
```

Log into **ifarm2402** and run the same command there. Each host runs up to ten workers; twenty total when both are busy. Both commands must use the same campaign path. They can start independently and share the queue through atomic batch claims. A host can finish before the other. Whichever finishes the last batches automatically merges the complete campaign under a single finalization claim. An explicit `merge` command is also available for recovery.

Detach from tmux with Ctrl-b then d; reattach on that host with `tmux attach -t eppi-minus-test`.

The test selects one eligible run per period/target combination (19 runs), processes **all HIPO files in each selected run**, and produces 19 final files. The default selection is the eligible run with the smallest total HIPO size. `inventory.json` lists selections and all discovered runs. Preparation fails on missing/empty datasets or ambiguous filename run numbers. Processing checks the inferred run against RUN::config. Excluded runs are recorded separately.

For custom selections, add `--test-runs selections.json` to preparation (JSON object mapping dataset IDs such as `RGC_summer2022_ND3` to integer run numbers). For unusual filenames, add `--run-map paths.json` (exact absolute filenames mapped to integer run numbers).

Preparation builds the optimized converter, compiles current helper and campaign Java sources, checks real QADB availability and processor compilation, and runs ROOT conversion tests. It freezes source, dependencies, QA data, and the polarization CSV. Workers use the frozen physics implementation; the current repository worker handles launch-time settings and records its own checksum and Java settings in each batch runtime.json and completion record. Existing processing.csh and run-groovy launchers are bypassed because they modify shared files.

## Outputs and monitoring

Each campaign has `inventory.json`, `manifest.json`, `checksums.json`, and these directories:

- `build/`: frozen sources/dependencies, compiled classes/converter, QA snapshot, run metadata.
- `claims/`: atomic per-batch ownership; failed claims are retained.
- `hosts/`: one coordinator lock per host, preventing duplicate coordinators.
- `work/`: exact input lists and temporary outputs; large text output is removed after successful conversion.
- `chunks/`: validated intermediate ROOT files and completion records.
- `logs/`: per-batch logs and `ifarm2401.log` / `ifarm2402.log`.
- `reports/`: per-dataset final entry totals and output checksums.
- `root/`: the 19 final ROOT files below, each with one `PhysicsEvents` tree.

```text
eppi_minus_RGB_spring2019_D2.root
eppi_minus_RGB_fall2019_D2.root
eppi_minus_RGB_spring2020_D2.root
eppi_minus_RGC_summer2022_ND3.root
eppi_minus_RGC_summer2022_C.root
eppi_minus_RGC_summer2022_CH2.root
eppi_minus_RGC_summer2022_ET.root
eppi_minus_RGC_summer2022_He.root
eppi_minus_RGC_fall2022_ND3.root
eppi_minus_RGC_fall2022_C.root
eppi_minus_RGC_fall2022_CH2.root
eppi_minus_RGC_fall2022_ET.root
eppi_minus_RGC_fall2022_He.root
eppi_minus_RGC_spring2023_ND3.root
eppi_minus_RGC_spring2023_C.root
eppi_minus_RGC_spring2023_CH2.root
eppi_minus_RGC_spring2023_CD2.root
eppi_minus_RGC_spring2023_ET.root
eppi_minus_RGC_spring2023_He.root
```

Test and production use separate campaign roots. Their final filenames match, but test files contain only the selected runs. No files are published directly under the parent `/volatile/clas12/thayward/eppi-/root/`; this prevents the test replacing full-campaign outputs.

```sh
python3 processing_scripts/exclusive_pion_campaign/campaign.py status --campaign /volatile/clas12/thayward/eppi-/test_20261122
```

Batch logs include event, skipped-first-record, missing-bank, and accepted-row counts. A successful empty tree is permitted and flagged in its final report and merge output; inspect such cases before the full launch. A run entirely rejected by QA or channel cuts can legitimately yield no rows.

## Live progress and empty-run diagnostics

New campaigns print START, CONVERT and DONE messages, separate processing/conversion timings, and a combined two-host snapshot every 30 seconds while local jobs are running. Snapshots show completed/active/queued/failed batches, output rows, completed input GiB, and the host and current stage of each active batch. Processor progress files update about every 15 seconds (checked every 1,000 events), with an immediate startup and final update. Long file-open or QA-startup waits may leave an older timestamp visible.

From a separate terminal, monitor an existing campaign without stopping its workers:

```sh
python3 processing_scripts/exclusive_pion_campaign/campaign.py monitor --campaign /volatile/clas12/thayward/eppi-/test_20261122_v3 --interval 30
```

Use `--once` for one snapshot. Older campaigns without progress files still show batch completion, host ownership, elapsed claim time and totals. Their frozen processors do not acquire new event counters retroactively.

For newly prepared campaigns, each batch log ends with an ordered event selection summary: seen, skipped first records, missing banks, invalid/nonphysics records, explicit run exclusions, zero helicity, QADB rejection, QADB acceptance, no selected electron, no selected proton, no selected pion, particle-pair candidates, kinematic rejections, accepted events, and output rows. QA rejection counts apply after the nonzero-helicity check. Candidate-pair counts are distinct from event counts because an event can have more than one pair. Empty outputs print these counters directly. Completion JSON files retain the same diagnostics. Existing zero-entry trees alone cannot identify QADB as their cause.

## Recovery

Successful chunks are reused. Failed chunks never enter a merge. After the two host coordinators have exited, release failed claims and relaunch:

```sh
python3 processing_scripts/exclusive_pion_campaign/campaign.py retry-failed --campaign /volatile/clas12/thayward/eppi-/test_20261122
csh -f processing_scripts/exclusive_pion_campaign/run_test.csh /volatile/clas12/thayward/eppi-/test_20261122
```

Failed retry uses the same frozen software. If code/calibration changes are required, create a fresh campaign. Do not rerun preparation on an existing campaign root. If preparation fails midway, use a fresh directory after fixing the reported problem.

A crash or killed terminal session can leave host/batch ownership behind. Inspect `owner.json` and confirm the named coordinator and its children have stopped on that host before manually removing a stale claim. There is no timeout-based reassignment: a slow or disconnected live worker must not race a replacement.

## Recover conversion failures without repeating HIPO processing

After both node worker commands have exited, update the checkout and run once on either node:

```sh
csh -f processing_scripts/exclusive_pion_campaign/environment.csh python3 processing_scripts/exclusive_pion_campaign/campaign.py recover-conversion --campaign /volatile/clas12/thayward/eppi-/test_20261122_v2
```

This command uses the current checkout to build a repaired converter under `recoveries/`, records its source/binary checksums, and converts existing text only when the batch log proves the processor completed. It verifies the ROOT entry count against the processor's summary. It retains all text, logs and original frozen software, skips previously successful chunks, and merges when all batches are complete. Do not run worker commands during recovery. Converter failures print their actual messages. If a processor did not finish, that batch is left for inspection rather than rerunning it automatically.

The shared metadata CSV contains six-column RGC polarization rows and five-column legacy charge records. Only the requested ND3 run's six-column polarization row is used; unpolarized target polarization is explicitly zero. Future campaign completion records split processing and conversion timings. Existing failed campaigns have no recorded stage timings, so their HIPO processing speed cannot be inferred from converter exit messages alone.

## Full launch, after the test has been reviewed

Prepare once on either node:

```sh
csh -f processing_scripts/exclusive_pion_campaign/environment.csh python3 processing_scripts/exclusive_pion_campaign/campaign.py prepare --repo "$PWD" --campaign /volatile/clas12/thayward/eppi-/full_20261122 --mode full
```

Then run separately on **each** node:

```sh
tmux new -s eppi-minus-full 'csh -f processing_scripts/exclusive_pion_campaign/run_full.csh /volatile/clas12/thayward/eppi-/full_20261122'
```

Inputs are the three specified RGB cache paths, all ten specified RGC cache target/period paths, and the six ET/He directories. Production batches are capped at approximately 4 GiB of HIPO input (an individual larger file stays intact). Runs are divided into batches and merged back by dataset. Long batches are offered first to reduce the long tail. Workers perform processing then conversion sequentially within their ten-per-host budget. Merging is sequential after all batches finish; no extra parallel merge pool is created.

## Selection and corrections

- p1 = proton (2212), FD or CD; p2 = pi- (-211), FD only. Additional particles are allowed; every accepted p/pi- pair is written with original REC::Particle indices.
- Preserve the first-record skip, nonzero beam helicity requirement, existing QADB defects/Misc allowances, and existing explicit run exclusions. Stored helicity follows the existing ThreeParticles run-dependent sign mapping, including the RGB spring-2019 and low-energy fall-2019 sign reversals; RGC retains the raw sign.
- Reuse existing electron PID/fiducial cuts, pion `-3 < chi2pid < 3`, proton lower momentum thresholds, theta/vertex/fiducial cuts, and proton energy-loss correction. The **proton upper momentum cut is removed**. Existing proton chi2pid cuts remain disabled; diagnostic values are stored.
- Recompute Q2, W, x, y from corrected electron momentum. Require Q2 >= 1 GeV^2, W >= 2 GeV, y <= 0.80. Use the existing run-based BeamEnergy helper.
- Apply the pi+ analysis-note RGC **electron** momentum correction once, after detector/PID selection and before all kinematics. Scale all three momentum components, preserving direction, and recompute energy with the electron mass. Su22 and Sp23 use their respective six-sector tables; Fa22 selects the negative/positive-solenoid table from RUN::config.solenoid. Theta is in degrees and clamped to [8,15]. RGB receives no RGC correction. No pi+ momentum correction is applied to pi-.
- `electron_correction_reference.tex` preserves the tables and captions used from the locally exported pi+ note. The requested Overleaf URL could not be fetched directly in this environment.
- Stationary neutron reference throughout: no Fermi sampling. t = (q - p_pi)^2; tmin is the exact forward two-body gamma* n -> pi- p limit with finite pion and nucleon masses. tprime = t - tmin (normally nonpositive); plot -tprime. Do not clamp reconstructed tprime.
- Production phi1/phi2 retain the existing oriented-plane convention. No dihadron phiH/phiR/rest-frame angle calculations. Missing-system and background variables are retained without tight exclusivity cuts.
- Converter checks an exact schema header, exact field counts, finite numbers, integral IDs, expected run, nonzero helicity, and detector requirements. ND3 requires run-level target-polarization metadata. Unpolarized targets explicitly have target_pol=0. Beam polarization uses the existing converter table.

## Complete tree schema (57 branches)

Integer branches:

```text
runnum evnum helicity p1_index p2_index detector1 detector2
sector_e sector1 sector2 fiducial_status
```

Double branches:

```text
beam_energy beam_pol target_pol
e_p e_theta e_phi vz_e
p1_p p1_theta p1_phi vz_p1
p2_p p2_theta p2_phi vz_p2
chi2pid1 chi2pid2 beta1 beta2
Q2 W x y nu
t tmin tprime
Mx2 Mx2_1 Mx2_2 Emiss pTmiss Mh
z2 xF2 pT2
phi1 phi2 Delta_phi
open_angle_ep1 open_angle_ep2
DepA DepB DepC DepV DepW
```

Momenta/energies are GeV, squared invariants GeV^2, vertices cm. Lab and production angles are radians, azimuths in [0,2pi); opening angles are degrees. Delta_phi is phi2-phi1 without wrapping. detector1/2: FD=1, CD=2; sector values use the existing REC::Track lookup (CD may have no sector). fiducial_status follows the existing convention (3 all pass; 0/1/2 only electron/proton/pion fails; -1 multiple fail).

Mx2 is (q+n-p-pi)^2; Mx2_1 is (q+n-p)^2; Mx2_2 is (q+n-pi)^2. Emiss is the full missing four-vector's energy; pTmiss its transverse momentum relative to the beam. z2=E_pi/nu; xF2=2 p_pi_parallel*/W and pT2 is pion transverse momentum relative to q in the gamma*n CM. Mh is the p/pi invariant mass. The electron and proton kinematics stored are corrected, while pi- momenta retain the current reconstruction values.

## Validation

Local checks: campaign Java/helper compilation; Groovy processor compilation with a QA test double; physical two-body events checking conservation, phi sign, pion pT/xF and exact t/tmin/tprime; electron correction periods/sectors/endpoints; 19-dataset inventory; batching; 20 simultaneous contenders for a batch; failure isolation; C-shell syntax; exact Java/C++ schema agreement. The C++ converter compiles locally. The local ROOT runtime has a pre-existing Cling/macOS SDK incompatibility, so actual ROOT conversion tests are enforced during farm preparation. The real QADB and real HIPO data are checked by the farm test launch.

## Restore original memory settings on an existing campaign

The worker now uses the original processing launcher defaults (1 GiB initial, 2 GiB maximum), removing the campaign-specific 1.5 GiB maximum and ActiveProcessorCount=1. This restores the baseline; adequate memory still needs verification on representative files. Advanced overrides can use PION_JAVA_OPTS. Each launched batch records its settings in work/BATCH/runtime.json.

For a stalled old worker, stop that specific Java process normally with kill -TERM PID, wait for both node coordinator commands to return, and update the repository files on the farm. Release failed batches once with retry-failed, then run the existing run_test.csh command separately on each host against the SAME campaign directory. Completed chunks are skipped. Frozen physics sources and checksums are preserved, and the updated worker supports old processors accepting three arguments. Do not prepare a fresh campaign just to change these runtime settings. An exited process might have completed or failed; its completion record and log establish which.
