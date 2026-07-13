# Automatic cut-variation systematics

Set in `main.cpp`:

```cpp
cut_variation_opts.enabled = true;
```

A single execution first completes the nominal analysis and then automatically produces:

- `exclusivity_95.csv`: 2-sigma centered cuts and 95% upper-quantile cuts
- `exclusivity_99p99.csv`: 4-sigma centered cuts and 99.99% upper-quantile cuts
- `fiducial_loose.csv`: all angular lower boundaries -2 degrees and upper boundaries +2 degrees
- `fiducial_tight.csv`: all angular lower boundaries +2 degrees and upper boundaries -2 degrees

The variation runner clones the completed nominal CSV and recomputes only the cut-dependent chain:

1. exclusivity extraction
2. event and MC counts
3. fixed nominal current-correction application
4. pi0 contamination
5. pi0-subtracted signal yields
6. acceptance
7. unfolding
8. cross sections
9. normalized cross sections
10. Barlow and RMS systematic extraction

The exact effective nominal data and MC current corrections are inferred bin-by-bin from the completed nominal CSV and reapplied to each varied count. Therefore the variation study does not refit the current-dependence correction and does not mix that systematic with the cut systematic.

Barlow filtering uses `B >= 1`. Both raw and filtered RMS values are retained. The systematic columns and the point-to-point total are stored as absolute cross-section uncertainties in `nb/GeV^4`.

Outputs are under `output/cut_variation_systematics/`.
