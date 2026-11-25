# Loose e′pπ− campaign

This is a separate campaign. The original `exclusive_pion_campaign` and its outputs are unchanged. Proton is p1 and π− is p2. The loose manager rejects original campaign manifests. It writes one final ROOT file per period and target, with the prefix `eppi_minus_loose_`.

## Selection

Hard requirements:
- Event Builder PID 11, 2212 and −211 for the electron, proton and pion.
- **Raw reconstructed electron momentum strictly greater than 2 GeV.** This is tested before electron corrections. A corrected electron may consequently have momentum below 2 GeV.
- Finite, nonzero reconstructed particle momenta.
- Required RUN::config, REC::Event and REC::Particle banks, a valid matching run, nonzero beam helicity and the original QADB mask and run exclusions. The first train-file record is skipped. The summer-2022 Hall C exclusion and spring-2023 outbending exclusion are retained.

No additional fiducial, chi2pid, vertex, hadron momentum, hadron angular, DIS or exclusivity cuts reject candidates. All electron × proton × pion combinations are written. Electrons need not be the particle at index zero. FD, CD, FT and other detector statuses are retained and tagged. In particular, CD pions are now retained so that their rejection can be made offline.

Calorimeter, trajectory, track and Cherenkov banks are optional. Missing diagnostic information does not remove a candidate. A helper error sets the affected flag to unavailable. Invalid momentum tracks are counted, as are electrons below the required threshold.

RGC electron momentum corrections and the original proton energy-loss helper are retained when evaluable. Raw Cartesian momenta are saved for all three particles, together with correction status. For an unavailable correction, the stored corrected vector defaults to its raw vector. The correction status indicates execution, not independent validation of the calibration for newly admitted tracks.

## Flags and fields

Selection flags use 1 = passes, 0 = fails and −1 = unavailable or not applicable. A missing bank must not be interpreted as passing a cut.

- `pass_e_vertex`, `pass_e_sampling`, `pass_e_diagonal`, `pass_e_pcal1`, `pass_e_pcal2`, `pass_e_dc`: original electron helper decisions. PCAL modes 1 and 2 are both retained.
- `pass_p1_momentum`: original lower proton thresholds, CD >0.3 GeV, FD >0.42 GeV for positive torus and >0.50 GeV for negative torus. No upper threshold.
- `pass_p1_theta`: proton polar angle <64.23 degrees.
- `pass_p1_vertex`, `pass_p1_dc`, `pass_p1_cvt1`, `pass_p1_cvt2`: proton vertex and fiducial helper decisions.
- `pass_p1_chi2pid`: optional |chi2pid|<3 diagnostic. This was not part of the original proton selection.
- `pass_p2_pid`: original pass-2 charged-hadron chi2pid helper. `pass_p2_vertex`, `pass_p2_dc`, `pass_p2_cvt1`, `pass_p2_cvt2` retain additional offline decisions.
- `pass_e_legacy`, `pass_p1_legacy`, `pass_p2_legacy`: conjunctions of the original track selections. The electron flag includes index zero and FD. The pion flag includes FD. These aggregate flags are zero when a required component is unavailable.
- `pass_dis`: Q²≥1 GeV², W≥2 GeV, y≤0.8 and ν>0. `pass_q2_1p5` records Q²>1.5 GeV² separately. Neither rejects candidates during processing.
- `kinematics_valid`: all derived kinematic values are finite. Outside their mathematical domain, derived quantities are saved as NaN rather than invented values. Use this flag and the desired DIS selection before analysis. Nonfinite optional vertex, chi2pid or beta measurements are also saved as NaN.
- `electron_correction`: 1 applied, 0 RGB (not required), −1 unavailable. `proton_correction`: 1 helper evaluated with finite output, −1 unavailable.
- `detector_e`, `detector1`, `detector2`: 1 FD, 2 CD, 3 FT, 0 other. Signed REC::Particle statuses are `status_e`, `status1`, `status2`.
- `e_index`, `p1_index`, `p2_index` and run/event identify each candidate combination. `e_raw_px/py/pz`, `p1_raw_px/py/pz`, `p2_raw_px/py/pz` are uncorrected components in GeV.

The original `fiducial_status` summary is retained. Value 3 means all three pass. Values 0, 1 and 2 identify a sole electron, proton or pion failure. Value −1 means multiple failures or unavailable information. Prefer explicit flags for new selections.

To approximate the previous selection, require all three `pass_*_legacy` flags, `pass_dis`, an available electron correction and `proton_correction==1`. There is still no upper proton momentum restriction. An offline comparison must also account for the old campaign requiring all optional banks to be present. This loose campaign retains such otherwise missing-bank events deliberately.

The target reference remains a stationary neutron. t=(q−pπ)² is the mesonic definition, tmin is the exact forward gamma*n→π−p value and tprime=t−tmin. `phi1` is the proton angle and `phi2` the pion angle. Momentum is in GeV, angles in radians except opening angles in degrees. Raw/corrected track momenta and these definitions preserve the existing analysis conventions.

Complete text schema, in column order (the ROOT tree adds `beam_pol` and `target_pol`):

```
runnum evnum helicity beam_energy p1_index p2_index detector1 detector2 sector_e sector1 sector2 fiducial_status e_p e_theta e_phi vz_e p1_p p1_theta p1_phi vz_p1 p2_p p2_theta p2_phi vz_p2 chi2pid1 chi2pid2 beta1 beta2 Q2 W x y nu t tmin tprime Mx2 Mx2_1 Mx2_2 Emiss pTmiss Mh z2 xF2 pT2 phi1 phi2 Delta_phi open_angle_ep1 open_angle_ep2 DepA DepB DepC DepV DepW e_index detector_e status_e status1 status2 e_raw_px e_raw_py e_raw_pz p1_raw_px p1_raw_py p1_raw_pz p2_raw_px p2_raw_py p2_raw_pz chi2pid_e beta_e pass_e_vertex pass_e_sampling pass_e_diagonal pass_e_pcal1 pass_e_pcal2 pass_e_dc pass_p1_momentum pass_p1_theta pass_p1_vertex pass_p1_dc pass_p1_cvt1 pass_p1_cvt2 pass_p1_chi2pid pass_p2_pid pass_p2_vertex pass_p2_dc pass_p2_cvt1 pass_p2_cvt2 pass_e_legacy pass_p1_legacy pass_p2_legacy electron_correction proton_correction pass_dis pass_q2_1p5 kinematics_valid
```

`run_info.csv` is the campaign-local copy of the current ND3 run-information table, including RGB charges and the previously missing run 17523. Preparation freezes it with checksums. Override it explicitly with `prepare --run-info /path/to/table.csv` if needed. ND3 polarization metadata remain mandatory even for empty files. RGB and auxiliary targets have zero target polarization.

## Launch

Copy this entire directory into `processing_scripts/exclusive_pion_loose_campaign/` in the farm repository on `exclusive_pion_processing_branch`. Commands below run from that repository root. Prepare once on either host, then launch workers separately on each host. There is no SSH coordinator. Ten workers per host give twenty total. Use a fresh campaign directory for every independently prepared test or full production.

Prepare a one-run-per-dataset smoke test once. The explicit `test_runs.json` defaults select known physics runs rather than tiny metadata-only files. Some input files are large, so this is an end-to-end test rather than a short event-limit test.

```csh
csh -f processing_scripts/exclusive_pion_loose_campaign/environment.csh python3 processing_scripts/exclusive_pion_loose_campaign/campaign.py prepare --repo "$PWD" --campaign /volatile/clas12/thayward/eppi-/loose/test_20261124 --mode test
```

On each farm node, open tmux manually, then run the processing command inside it:

```csh
tmux new -s eppi-loose-test
csh -f processing_scripts/exclusive_pion_loose_campaign/run_test.csh /volatile/clas12/thayward/eppi-/loose/test_20261124
```

Check progress from another terminal:

```csh
python3 processing_scripts/exclusive_pion_loose_campaign/campaign.py monitor --campaign /volatile/clas12/thayward/eppi-/loose/test_20261124 --once
```

When the test is finished, audit its ROOT files, schemas and completion records:

```csh
csh -f processing_scripts/exclusive_pion_loose_campaign/environment.csh python3 processing_scripts/exclusive_pion_loose_campaign/audit_campaign.py --campaign /volatile/clas12/thayward/eppi-/loose/test_20261124 --checksums
```

After a successful test, prepare full production once:

```csh
csh -f processing_scripts/exclusive_pion_loose_campaign/environment.csh python3 processing_scripts/exclusive_pion_loose_campaign/campaign.py prepare --repo "$PWD" --campaign /volatile/clas12/thayward/eppi-/loose/full_20261124 --mode full
```

Run separately on both hosts inside manually opened tmux sessions:

```csh
csh -f processing_scripts/exclusive_pion_loose_campaign/run_full.csh /volatile/clas12/thayward/eppi-/loose/full_20261124
```

Monitor once or continuously:

```csh
python3 processing_scripts/exclusive_pion_loose_campaign/campaign.py monitor --campaign /volatile/clas12/thayward/eppi-/loose/full_20261124 --once
python3 processing_scripts/exclusive_pion_loose_campaign/campaign.py monitor --campaign /volatile/clas12/thayward/eppi-/loose/full_20261124 --interval 30
```

The manager processes larger batches first across all datasets. Atomic claims prevent hosts duplicating work. It retains the original 2 GB Java heap and progress updates, and evaluates track flags/corrections once per track rather than for every combination. The looser selection can substantially increase output size and processing time.

Failed batches can be released without deleting completed results:

```csh
python3 processing_scripts/exclusive_pion_loose_campaign/campaign.py retry-failed --campaign /volatile/clas12/thayward/eppi-/loose/full_20261124
```

Then rerun `run_full.csh` separately on the desired hosts. Completed chunks are retained. Do not delete claims for active workers. When all batches finish successfully, final datasets merge automatically. Audit with:

```csh
csh -f processing_scripts/exclusive_pion_loose_campaign/environment.csh python3 processing_scripts/exclusive_pion_loose_campaign/audit_campaign.py --campaign /volatile/clas12/thayward/eppi-/loose/full_20261124 --checksums
```

## Final files

Full-production outputs go under `/volatile/clas12/thayward/eppi-/loose/full_20261124/root/`:

```
eppi_minus_loose_RGB_spring2019_D2.root
eppi_minus_loose_RGB_fall2019_D2.root
eppi_minus_loose_RGB_spring2020_D2.root
eppi_minus_loose_RGC_summer2022_ND3.root
eppi_minus_loose_RGC_summer2022_C.root
eppi_minus_loose_RGC_summer2022_CH2.root
eppi_minus_loose_RGC_summer2022_ET.root
eppi_minus_loose_RGC_summer2022_He.root
eppi_minus_loose_RGC_fall2022_ND3.root
eppi_minus_loose_RGC_fall2022_C.root
eppi_minus_loose_RGC_fall2022_CH2.root
eppi_minus_loose_RGC_fall2022_ET.root
eppi_minus_loose_RGC_fall2022_He.root
eppi_minus_loose_RGC_spring2023_ND3.root
eppi_minus_loose_RGC_spring2023_C.root
eppi_minus_loose_RGC_spring2023_CH2.root
eppi_minus_loose_RGC_spring2023_CD2.root
eppi_minus_loose_RGC_spring2023_ET.root
eppi_minus_loose_RGC_spring2023_He.root
```

## Validation

Preparation compiles the Java processor and converter, runs kinematic closure, loose-selection and repeated-metadata tests, exercises the converter's shared strict row parser, writes and reopens small ROOT trees, then resolves the installed QADB and compiles the Groovy processor. Preparation must succeed before launching workers. These checks supplement the real one-run-per-dataset smoke test, which is required before full production.

Local Java tests, C++ parser/metadata tests, Python manager/audit tests and Groovy syntax parsing pass. The ROOT converter compiles locally, but local ROOT 6.32 cannot complete its runtime IO test because its C++ standard-library modules fail to load. Installed QADB resolution and ROOT runtime IO remain to be verified by preparation on the farm. No farm jobs have been launched.
