# Real-sample pi0 mass reconstruction closure

Existing native tests already exercise recovery of known photon pairs, rounded angles, photon-exchange consistency and response changes. This supplementary check uses the actual DATA and MC per-event exports from the unordered control study. It performs no fits or ROOT reads.

```sh
git pull && python3 test_control_mass_closure.py && python3 study_control_mass_closure.py --control-study "$PWD/test_runs/unordered_control_response_001" --output "$PWD/test_runs/control_mass_closure_001"
```

The farm scratch exports must still exist; the Desktop report alone does not contain them. Checksums, schema and inventory are verified before report creation. Previous outputs are preserved; use a fresh output path.

`mass_closure.csv` reports saved versus rebuilt mass for the identical resolved events, globally and inside the simultaneous common windows. Values include median, 1–99% range, RMS, maximum difference, fractions exceeding 0.1/1 MeV and saved/rebuilt peak widths. Those fractions are descriptive, not automatically chosen acceptance criteria. Plots use common-window events and show mass overlays and per-event mass differences.

- DATA / nominal MC: compare saved mass with the invariant mass reconstructed from their unmodified recovered photons.
- Frozen MC: compare the reported modified mass with the invariant mass from its exported modified photons. This checks the two output paths agree; it is not zero-response closure.
- Unresolved pairs remain counted explicitly and are omitted from both sides of the matched-event overlay.
- The stored mass was an input to photon recovery. This checks numerical consistency and possible reconstruction offsets; it is not independent validation of the physical mass model or pair identity.

Only main-analysis-supported control groups are plotted: FD proton with two FD photons, and CD proton with FD–FD, FD–FT or FT–FT photons. FD-proton groups involving FT photons are omitted. The unordered exporter retains its six-group encoding for compatibility and bookkeeping; its default comparison plots/scores now use the same four supported groups.

## Selection sequence for a later response refit

Use the legacy signal/background template method to obtain initial DATA selection boundaries. Keep those boundaries fixed while comparing response trials. Fit the response using signal-enriched, leave-one-out distributions with the appropriate background prediction. Retain independent reserved-observable and MC-half checks.

Always modify the full pre-cut reconstructed MC and then evaluate selection. Do not permanently discard initially failing MC before response fitting; some can enter after modification. DATA coordinates are unchanged. Broad inclusive tails are diagnostic background information, not automatically a target signal-resolution shape.

After the response fit, rerun the legacy template/cut procedure with the modified templates, retaining separate DATA/MC 95% windows if that remains the production convention. Reevaluate acceptance and pi0 transfer using the actual resulting selections, and iterate only if window changes materially affect the fitted response. Selected-core agreement alone does not validate efficiency.
