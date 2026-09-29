# High-luminosity suite path/archive fixes — 2026-09-29

- `prepare_pi0_gk_inputs.py`: default Story handoff changed to the actual `.zip`; accepts extracted directories, ZIP archives, and tar archives.
- `prepare_pi0_gk_stage2_covariance.py`: default Story handoff changed to the actual `.zip`; accepts directories/ZIP/tar; provenance hashing is only attempted for archive files, so extracted-directory input no longer crashes at summary writing.
- `run_pi0_gk_partons_stage3.py`: native campaign discovery now supports the actual timestamped ZIP in addition to unpacked directories and legacy tar archives.
- Verified `prepare_pi0_gk_inputs.py` and `prepare_pi0_gk_stage2_covariance.py` end-to-end against the supplied `high_luminosity/import/fa18_rosenbluth_inputs_20260924T165948Z.zip` with default arguments.
- Python syntax compilation passed for the Stage 1/2/3, projection, reach-study, and closure-test scripts.
