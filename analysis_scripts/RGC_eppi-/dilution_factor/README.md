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
- `figures/dilution.png`, `figures/packing_fraction.png`, `figures/cd2_ch2_transfer.png`: unconnected points with matching y-axis ranges across each canvas.
- `figures/cd2_ch2_distributions.png` and `cd2_ch2_distributions.csv`: measured Sp23 CD2 and CH2 yields divided by their respective good charges, versus minus-tprime in four xB panels. Uses 0.05 GeV² display bins. These do not change the 24-bin dilution calculation. Curves are not area normalized, and no areal-density conversion to nuclear luminosity is applied.

Check `cutset_id` before combining these results with later asymmetry inputs. Review flags must be resolved before using the corresponding estimate as a nominal correction.

Checks: `python3 -m unittest discover -s analysis_scripts/RGC_eppi-/dilution_factor/tests -v`

## xB-only comparison and ratio display

The CD2/CH2 distribution canvas now includes ratio panels in the six analysis tprime bins below the finer yield distributions. Yield labels are Counts / nC, using the source charge units without the former factor of one million. All four yield panels share a y range, as do the four ratio panels.

`method1_xb_comparison.csv`, `method1_xb_comparison.json`, and `figures/method1_xb_comparison.png` integrate each target count over the six tprime bins before the transfer ratio and Method 1 are calculated. These are comparison outputs, separate from the 24-bin nominal values. Statistical trials retain shared Sp23 inputs. Applying one value to all six tprime bins would introduce a common dilution uncertainty across those bins and assumes a tprime-independent correction within each xB interval.

## Pooling comparisons

Every execution also produces `method1_binning_comparison.csv/json`, `method1_tprime_half_differences.csv`, `binning_comparison_covariance.npz`, and the figures `method1_tprime_only.png` and `method1_two_tprime_halves.png`. Three schemes are compared: four xB bins integrated over tprime, six tprime bins integrated over xB, and four xB bins each split into [0.05,0.65) and [0.65,1.25) GeV². Charges are unchanged when summing event counts within each group. The CD2/CH2 transfer is recomputed at each scheme's binning, before Method 1.

All comparisons reuse the same original 24-bin Poisson trials. JSON records bin memberships and covariance ordering. Joint covariance includes all periods and both dilution and packing factors within each scheme. High-minus-low differences include central 68% and 95% replica intervals. They are uncertainty summaries, not formal hypothesis-test confidence intervals. Large replica standard deviations may be driven by rare denominator crossings. Such cases receive a heavy-tail warning in the difference table. Comparisons and the three periods share Sp23 control data, so their apparent trends are not independent evidence. None of these outputs replaces the 24-bin nominal correction automatically.

The calculation also writes `period_consistency_summary.csv`, `period_consistency_bins.csv` and `period_consistency.json`. These compare the three period dilution values in each one-dimensional projection using a common-value GLS fit with the full shared-control covariance. Global tests fit a separate common value in each projected bin. The xB-only test has 8 degrees of freedom, the tprime-only test 12. Probabilities are approximate and statistical-only. The two projections reuse the same events. Their agreement does not demonstrate equal dilution in each two-dimensional bin.

## Combined-period 24-bin result

`joint_period_method1.json/csv` and `figures/joint_period_method1.png` contain one shared dilution factor per original two-dimensional bin. A joint Poisson likelihood fits all 16 measured target/period counts, with independent target rates for each period and a shared Sp23 CD2/CH2 transfer. Three material-balance constraints impose the same DF in each period. This uses individual period charge normalizations without assuming identical accepted count rates. Each bin has 16 observations, 17 parameters and 3 constraints, leaving 2 degrees of freedom for an approximate consistency check. DF is bounded to (0,1). Zero-count auxiliary measurements enter the Poisson likelihood directly. The nominal statistical error is now symmetric: the sample standard deviation (ddof=1) of 2,000 combined fits to independent Poisson fluctuations of all target/period counts. This matches the pi+ reporting convention. The same Sp23 CD2 and CH2 draw is used coherently in all three period constraints. Bins are processed with six local workers by default. Use `--joint-replicas` and `--joint-workers` on the main command to change those settings. Boundary solutions are rechecked with multiple starting values. Failed-fit and boundary fractions, replica means and quantiles are retained in JSON. More than 1% failed fits stops the calculation. The physical DF bounds remain part of the joint estimator, so a symmetric error does not assert a Gaussian sampling distribution. `joint_period_replicas.npz` retains every successful fitted replica for downstream propagation and checks. They do not include charge, material, transfer-model or radiative systematics. Original separate-period results remain available for cross-checks.

`figures/period_agreement.png` presents both integrated consistency checks with matching scales. The projected checks use replica-based statistical intervals, while the combined result uses symmetric replica standard deviations. The common-DF assumption is checked by the joint fit deviance. It does not establish identical acceptance across periods or remove the preliminary assumptions of Method 1.
