# CD-proton angle-dependent response test

Run from `analysis_scripts/dvcs_cross_section` after committing/pulling the changes:

```sh
git pull && python3 test_cd_proton_angle_response.py --root-smoke && python3 study_cd_proton_angle_response.py --joint-run "$PWD/test_runs/particle_response_joint_control_001" --frozen-study "$PWD/test_runs/frozen_response_kinematics_001" --workers 10 --scratch-root /scratch/thayward/dvcs_particle_response --output "$PWD/test_runs/cd_proton_angle_response_001"
```

The farm's original joint-response directory, reference studies and scratch caches/exports must still exist. Desktop report copies do not contain these large inputs. Provenance and checksums are checked before fitting. A new output directory is required; previous fits and production files are preserved.

## Hypothesis and limits

The frozen low/high diagnostic found excess DATA missing-transverse-momentum width at lower CD-proton polar angles in both channels. Test **one additional fractional momentum width**, independently for each period/polarity. All 20 saved particle-response parameters remain fixed. The background nuisance follows the same constrained profiling as the original joint fit.

For original reconstructed CD-proton polar angle theta:

```
extra_sigma_p / p = a * clamp((55 degrees - theta) / 20 degrees, 0, 1)
total_fractional_sigma = hypot(saved_fractional_sigma, extra_sigma_p / p)
p_new = p * (1 + saved_bias + total_fractional_sigma * Z)
```

The smooth ramp has full amplitude below 35 degrees and zero above 55 degrees. These are explicit diagnostic choices, motivated by the approximately 54-degree split in the previous study; neither edge is fitted. The same standardized Gaussian is reused across hypotheses to reduce comparison noise. Quadrature widths describe independent Gaussian contributions in distribution. The default amplitude scan is 0–0.08: this is additional fractional **particle momentum width**, not a percentage increase of a missing-mass width. Bounds can be changed with `--max-sigma`; a boundary optimum is flagged, not accepted automatically.

The proton law is identical for DVCS signal, pi0 leakage and two-photon control. Photon responses remain shared by detector, independent of photon role. The production 20-parameter headers remain unchanged: a checked isolated source copy is built for this study with a 21st parameter. Zero extension reproduces the original native response exactly. `angle_candidate.json` is deliberately a study-only artifact, not a production candidate file that another script might accidentally adopt.

## Selection and validation

- Joint epgamma (DVCSgen + AAOGEN leakage) and epgammagamma (AAOGEN) shape objective, equal channel weights; all three topologies.
- Training variables: Delta_phi, pTmiss, Emiss2 and Mx2 on MC half 0.
- Only training scores select the new amplitude; zero remains an eligible candidate.
- Other MC half, reserved theta_gamma_gamma/Mx2_2 and CD–FD regressions test the result. DATA histograms are shared, so this is not independent DATA validation.
- Inclusive scoring and low/high diagnostics retain their existing respective selection definitions. Low/high plots use full globals and common DATA exclusivity windows, omitting the plotted variable's exclusivity cut. No new windows are determined here.
- Low/high edges are the same DATA-defined per-period medians used in the frozen study. They are distinct from the fixed 35–55-degree ramp edges.
- No acceptance, contamination or cross-section changes are adopted. A shape improvement does not by itself establish the correct efficiency.

## Outputs

Open `REPORT.html`. Start with `plots/PERIOD/scores.png`, then `plots/PERIOD/angle_subdivisions/epgamma_CD_FT.png` and `epgammagamma_CD_FT.png`. Each low/high pair now uses the same observable. Check corresponding CD_FD plots before interpreting an improvement as useful. Inclusive overlays are under `plots/PERIOD/inclusive/`, including FD_FD as a track-response control (its mixture nuisance may still change).

`score_comparison.csv` contains channel-wise training, other-half and reserved scores. `PERIOD/scan.csv` records the scan; `angle_candidate.json` records amplitude and boundary status; `angle_residuals.csv` records low/high widths, shifts, counts and scores for both halves. Scores without adequate support are flagged. Low/high scores use the frozen-model reference variance for both hypotheses.

`photon_order_audit.csv` compares which photon is harder in both-FD control events after full globals, before exclusivity. DATA use recovered exported photon coordinates; MC use the frozen morphed control export. This is an ordering diagnostic, not an absolute-rate comparison. It does not establish whether a reversal originates in source processing or coordinate recovery. Role-specific photon-momentum fits remain out of scope.

## Cost and tests

No ROOT events are reread. Existing MC caches and DATA exports are reused. Seventeen coarse amplitudes are evaluated together, followed by a small refinement and final paired evaluation. Period fits share a ten-process cap; plotting is serialized because Matplotlib is not thread safe. Large temporary exports and logs remain under the chosen scratch root; only one report tree is produced.

`test_cd_proton_angle_response.py` compiles native core fixtures and verifies bit-for-bit zero-extension replay, angle taper, quadrature width, common proton response, untouched FD response, photon-role symmetry and training-only parameter selection. `--root-smoke` additionally builds real ROOT workers and exercises the extended control and signal evaluator paths and export replay. Full farm-data results require running the study command above.
