# ND3 dilution factor: preliminary Method 1

Run from the analysis repository on one node. This reads the completed nominal selection, not the original large campaign ROOT files. It should take seconds. No epoch subdivision or additional exclusivity fitting is performed.

`python3 analysis_scripts/RGC_eppi-/dilution_factor/determine_dilution_factor.py`

Defaults are `../channel_selection/output/selected_inputs.json`, `../imports/clas12_run_info_ND3_source.csv`, and `output/` beside this script. Override with `--selected-inputs`, `--run-info`, or `--output`. Requires numpy, uproot, matplotlib and the existing channel-selection dependencies. No changes to selected ROOT files are made.

## Calculation

Use the October 16, 2026 [RG-C dilution note](https://clasweb.jlab.org/wiki/images/7/72/Dilution_Factor_Determination_for_RG_C_10-16-26.pdf), Eqs. 7, 8, 10 and the **ND3** half of Table II. These updated coefficients differ from the older NH3 coefficients in the pi+ implementation. `model.py` records the exact numbers used. MT in the note is the helium-filled target, named `He` here. F is foils-only, named `ET` here.

For every target and each of the 24 nominal bins, n = N/(Q+ + Q-). Sum charge over unique runs present in the original pre-exclusivity ROOT export, including runs with no surviving exclusive candidates. Never infer luminosity from the selected run list or selected count. Missing or nonpositive source-run charge stops the calculation. The output records all included runs and charges, plus CSV runs omitted because they did not appear in that source ROOT sample, for audit.

For Sp23, use measured n_CD2. For Su22 and Fa22, follow Sec. III.3 item 2: measure R_b = n_CD2,Sp23 / n_CH2,Sp23 separately in each analysis bin, then use n_CD2,period = R_b n_CH2,period. This is the note's ratio of complete charge-normalized target rates. No additional ET or He subtraction is applied to the transfer ratio. Their period-specific contributions enter Method 1 itself. Transferring the ratio assumes its applicability across periods and is a preliminary approximation.

With the primed coefficients from Table II:

- S = n_CD2 + A' n_C + B' n_He + C' n_ET.
- D = n_CD2 + D' n_C + E' n_He + F' n_ET.
- DF = (n_ND3 - n_He) S / (n_ND3 D).
- PF_bath = G' (n_ND3 - n_He) / D.
- PF_cell = (5.86/5) PF_bath, following Eq. 10.

Packing fraction is only a cross-check. It is not averaged or inserted back into the dilution estimator. The integrated packing cross-check sums bin rates after applying the bin-specific CD2 transfer, using the same nominal selections. Different cuts across bins and changing target fills mean it is an effective period average, not an epoch-resolved physical measurement.

## Uncertainties and validity

20,000 independent Poisson count replicas propagate statistical uncertainties. Each Sp23 replica is reused in all three period calculations, preserving correlations from the shared CD2/CH2 transfer and between DF and PF inputs. Output includes standard deviations and central 16%-84% intervals. Rare denominator crossings can make standard deviations and covariance unstable, so inspect the intervals and crossing fraction too. No clipping, smoothing, pseudo-counts or replacement of problematic bins is performed.

A denominator sign-crossing fraction above 1%, an undefined estimate, a central DF or cell PF outside (0,1), or less than 99% finite replicas causes a review flag and `ready_for_asymmetry=false`. This status is a numerical screen, not final physics approval. Targets with fewer than 10 counts receive informational notes. A zero observed count also has zero variance under this Poisson plug-in calculation, so sparse ET bins require care. This matches the count-resampling convention, not an upper-limit treatment.

Uncertainties assume independent candidate-combination counts, as in pi+. They do not include event-level correlations from multiple retained combinations, charge uncertainty, target-material systematics, transfer model dependence, or radiative corrections. No 4% NH3 scale uncertainty is automatically assigned to ND3. RGB is not assigned an RGC dilution factor.

## Outputs

- `method1.csv`: 72 period/bin estimates, uncertainty intervals, stability flags.
- `method1.json`: downstream data, cut-set ID, coefficients, normalization and provenance.
- `cd2_ch2_transfer.csv`: 24 Sp23 transfer ratios and uncertainties.
- `packing_fraction_integrated.csv`: integrated cell packing cross-check per period.
- `statistical_covariance.npz`: 72-by-72 dilution and packing covariance matrices. Order is Su22 bins 1-24, Fa22 bins 1-24, Sp23 bins 1-24. Covariance uses replicas finite in all 72 components. Complete-replica counts are recorded in JSON. The file contains separate DF and PF matrices, not their joint cross-covariance.
- `figures/dilution.png`, `figures/packing_fraction.png`, `figures/cd2_ch2_transfer.png`.

Check `cutset_id` before combining these results with later asymmetry inputs. Review flags must be resolved before using the corresponding estimate as a nominal correction.

Checks: `python3 -m unittest discover -s analysis_scripts/RGC_eppi-/dilution_factor/tests -v`
