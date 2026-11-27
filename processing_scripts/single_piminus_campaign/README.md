# Single pi-minus campaign

Run from the repository root on `exclusive_pion_processing_branch`. This campaign writes to its own `/volatile/clas12/thayward/epi-/single/` directory. It does not reuse the detected-proton campaign outputs.

## Selection and output

One row per accepted electron–pi-minus combination. The trigger electron is REC::Particle row zero. Electron and pion selection call the shared `analysis_fitter.electron_test` and `pion_test` used by the pi-plus processor. Electron momentum, vertex, sampling fraction, diagonal, PCAL2 and DC cuts apply. The pion pass2 chi2pid cut applies. There is no proton requirement, proton veto, pion fiducial cut or exclusive missing-mass window. Forward and central pions are retained and tagged (`detector`: 1=FD, 2=CD, 0=other). Pion DC fiducials are diagnostic only. `fiducial_status` retains the legacy combined diagnostic, including PCAL1. It is not another processing cut.

Require nonzero helicity, QADB acceptance, Q2 >= 1 GeV², W >= 2 GeV and y <= 0.8. Preserve the existing run exclusions and RGC electron correction. No random Fermi momentum is injected. Kinematics use a stationary neutron, an unobserved recoil proton, exact forward tmin and signed tprime=t-tmin. The proton missing-mass peak is near 0.880 GeV². `p_p`, `p_theta`, `p_phi` describe the detected pion. Missing momentum describes the unobserved recoil, not the spectator momentum of the detected-proton analysis.

All 19 target/period datasets are included. Outputs are `root/epi_minus_DATASET.root`, tree `PhysicsEvents`. Frozen code, dependencies, intermediate text, ROOT chunks, logs and reports all stay inside the campaign directory on volatile. These are processing inputs for a subsequent single-pion analysis, not inputs to the existing three-particle selection without adaptation.

## Test first

Pull the branch on the farm. Prepare once, on either host, outside or inside tmux:

```sh
git pull --ff-only origin exclusive_pion_processing_branch
csh -f processing_scripts/single_piminus_campaign/environment.csh python3 processing_scripts/single_piminus_campaign/campaign.py prepare --repo "$PWD" --mode test --campaign /volatile/clas12/thayward/epi-/single/test_20261126
```

Preparation freezes dependencies and metadata, compiles Java and ROOT conversion, and runs kinematic, topology, record, metadata and converter checks. It must succeed before starting workers. The test uses one designated run per dataset.

On **each** of ifarm2401 and ifarm2402, start a separate tmux session, then run:

```sh
tmux new -s single-piminus-test
csh -f processing_scripts/single_piminus_campaign/run_test.csh /volatile/clas12/thayward/epi-/single/test_20261126
```

Detach with Ctrl-b then d. The two hosts share the queue and claim distinct batches. Ten worker JVMs per host, 2 GB heap each. Completed batches are reused on restart. Do not start a second worker manager on the same host while one is running.

From another window:

```sh
python3 processing_scripts/single_piminus_campaign/campaign.py monitor --campaign /volatile/clas12/thayward/epi-/single/test_20261126
```

When the queue is finished and merges are complete, audit:

```sh
csh -f processing_scripts/single_piminus_campaign/environment.csh python3 processing_scripts/single_piminus_campaign/audit_campaign.py --campaign /volatile/clas12/thayward/epi-/single/test_20261126 --checksums
```

Check the audit passes, inspect yields and any unexpectedly empty datasets before full production. The frozen metadata snapshot includes the additional RGB fall-2019 charge rows from development.

## Full weekend campaign

After the test passes, prepare once:

```sh
csh -f processing_scripts/single_piminus_campaign/environment.csh python3 processing_scripts/single_piminus_campaign/campaign.py prepare --repo "$PWD" --mode full --campaign /volatile/clas12/thayward/epi-/single/full_20261126
```

On each host, in tmux:

```sh
tmux new -s single-piminus-full
csh -f processing_scripts/single_piminus_campaign/run_full.csh /volatile/clas12/thayward/epi-/single/full_20261126
```

Monitor and audit with the same commands above, substituting `full_20261126` for `test_20261126`. Use `status` instead of `monitor` for a single snapshot. Do not prepare an existing campaign directory again. No submission or SSH coordination happens automatically.

## Local validation limits

Java compilation, kinematic closure, correction, topology, metadata and Python bookkeeping checks pass locally. The ROOT converter compiles locally, but this Mac ROOT installation fails at runtime due to stale Cling/SDK modules. The farm preparation deliberately runs the ROOT converter checks again and stops if they fail. No farm processing has been launched by these changes.
