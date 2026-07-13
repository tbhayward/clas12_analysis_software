# Automatic cut-variation total-count plot suppression

The automatic cut-variation runner now sets `TotalCountsOptions::make_plots = false`.
This suppresses the large per-period, per-topology and per-xB `total_counts_plots`
canvas production while retaining the complete numerical counting and CSV update.

The nominal workflow is unchanged because `TotalCountsOptions::make_plots` defaults
to `true`.

Old PNG files are not deleted automatically. Remove the old variation output before
rerunning if a clean output tree is desired:

```bash
rm -rf output/cut_variation_systematics/variations
```
