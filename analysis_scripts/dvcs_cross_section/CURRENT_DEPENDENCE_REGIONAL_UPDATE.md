# Regional current-dependence update

This update changes the nominal current-dependent reconstruction correction from a period-wide post-binning scalar correction to an event-level correction using the reconstructed neutral-particle region (FT or FD sector 1--6) and the event run/current.

## Raw kinematic means and yield totals are unchanged

`bin_means.cpp` is unchanged. `total_counts.cpp` continues to accumulate and write all original raw DATA and unit-weight generated/reconstructed MC yield columns exactly as before. The existing raw-yield analysis-note tables/plots are made from those original accumulators. A second weighted accumulator fills the already-distinct `normalized raw yield` DATA columns and `reconstructed current corrected yield` MC columns. Thus the kinematic-means/raw-yield-total analysis-note section does not need to be regenerated because of this correction change.

## Calibration stage

`current_dependence.cpp` now runs before `total_counts.cpp` in `main.cpp` and writes:

- the existing period-integrated current-efficiency columns for bookkeeping and cross-checks;
- FT/S1--S6 current-dependence results for `ep->epg` and direct `ep->eppi0`;
- DATA `f(region,theta_e)` and `f(region,theta_gamma)` diagnostics;
- analysis-note-oriented current-response and regional-correction plots;
- `output/dvcs_current_dependence/calibration/current_response_model.json`.

The integrated DATA cross-check now evaluates the fitted response at the charge-weighted current, rather than a reconstructed-event-count-weighted current.

For Sp19 Inb, the Fa18 Inb transfer is applied only to the DATA current response. The independently determined Sp19 MC response is retained.

The previous inclusive `theta_e` production correction is disabled by default. The new regional-theta plots are diagnostic; they are intended to determine whether an additional polar-angle dependence remains after conditioning on detector region.

## Event-level DATA correction

For DATA in period `p`, reconstructed neutral-particle region `r`, and run current `I`, the calibration stores the relative fit slope

`R_data(p,r,I) = 1 + s_rel(p,r) I`, with `s_rel = m/b`.

Each accepted event contributes weight

`w_data = 1/R_data`.

The fit uncertainty on `s_rel` includes the covariance of the fitted linear slope/intercept. In a final analysis bin, its uncertainty is propagated coherently as a regional nuisance via the derivative of the complete weighted yield; it is not added as an independent uncertainty once per event.

## Event-level reconstructed-MC correction

For reconstructed MC, each accepted event is assigned FT/S1--S6 and weighted by

`w_mc = 1/f_mc(p,r)`.

Generated MC remains unweighted. For the misidentified `ep->eppi0->epg` reconstructed sample, the default convention uses the `ep->epg` regional MC response because the reconstructed topology is `epgamma`.

## Downstream flow

`pi0_contamination.cpp` now consumes the already current-corrected reconstructed `ep->eppi0->epg` columns and no longer applies an additional period-wide scalar MC current factor. Acceptance already consumes the reconstructed current-corrected DVCS MC columns. `yield_totals.cpp` (the later current-binned diagnostic, not the raw total-count note output) now also uses the regional run/current DATA response instead of the legacy period-wide scalar.

Automatic exclusivity-cut variations retain the same calibration but re-evaluate the event-level weighting for the varied event selection. They no longer overwrite corrected yields by rescaling from the nominal CSV, so changes in FT/sector composition are preserved.

## Analysis-note outputs

The nominal current stage produces an `analysis_note` directory under `output/dvcs_current_dependence/` containing current-response canvases, the regional correction relative to the previous integrated correction, and regional polar-angle diagnostics. The detailed developer plots remain under the per-channel diagnostic directories.

## Validation note

A full compile could not be completed in the ChatGPT execution container because ROOT and `root-config` are not installed there (`TFile.h`/other ROOT headers are unavailable). The modified source was structurally checked and the project should be compiled/run in the normal JLab environment before treating the numerical outputs as validated.
