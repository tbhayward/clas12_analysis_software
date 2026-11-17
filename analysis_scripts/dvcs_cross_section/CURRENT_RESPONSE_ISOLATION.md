# Nominal current-systematic response isolation

Nominal total counting explicitly exports its calibration responses. Other TotalCountsOptions callers now default to no export and the automatic cut-variation runner explicitly disables it. This prevents alternate-binning BSA variations and nominal cut variations from replacing the nominal response CSV.

The exporter writes an explicit nominal record for every bin/sample/period, including genuine zero counts. It also writes `<response CSV>.provenance.csv` with the eight bin boundaries and deterministic checksums of the final exclusivity cuts, current calibration and response CSV. These are identity checks, not cryptographic authentication.

The reader rejects missing metadata, mismatched row counts, changed bin boundaries and changed files before interpreting row-indexed responses. It independently checks every nominal component against the production CSV, including topology-summed current-corrected DATA yields. A wrong selection with identical binning therefore also fails. Old response files without metadata must be regenerated through nominal counting.

Explicitly zero selected background numerators contribute zero background. Missing records still fail validation. Nuisance excursions crossing zero signal retain their signed responses. Invalid acceptance/control denominator excursions fail explicitly rather than silently dropping an uncertainty contribution. Current-systematic columns are cleared before supported values are filled so that reevaluation cannot retain stale entries.

Run `python3 test_current_response_isolation.py`. The ROOT-free test compiles actual reader and propagation functions and tests export isolation, provenance checks, nominal-yield mismatches, zero versus missing controls, shared calibration and Sp19 transfer, signed zero crossings and invalid excursions. The full farm extraction remains to be exercised after maintenance.

This patch does not repair previously computed uncertainty columns. The nominal response sums must be recounted using the existing fitted detector response, final cuts and current calibration. The fit itself need not be repeated. Regenerate current uncertainties, the current-subtracted period residual, correlated-scale allocations, summaries and global fits using that correlated direction. Direct central-value ratios and statistical-plus-point-to-point pulls do not consume the affected allocation.

## Saved comparison plots

`external_scripts/plot_cross_section_ratio_pull_summary.py` creates Pass-1, BH and KM15 ratio/pull summaries from the saved comparison tables. All left axes are 0–2.5. Only Pass-1 has a ±31% normalization reference band around unity. The band is not added to each point's pull denominator. Top-edge triangles mark percentile ranges extending beyond the displayed ratio limit. The right-column range is common across the three rows of each figure. Pass-1 uses its 1289 matched points while the models use all 1445 supported production points.
