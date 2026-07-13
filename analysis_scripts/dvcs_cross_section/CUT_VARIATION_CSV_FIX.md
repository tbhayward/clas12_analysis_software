# Cut-variation CSV parsing and validation fix

## Problem

Some variation CSV files contained tuple cells with an extra quote layer, for example:

```csv
"""(2.483319049, 0.1106748048, 0)"""
```

After ordinary CSV parsing this left literal quotation marks around the tuple, so the
cross-section parser rejected the value. The systematic stage then treated the missing
variation as a zero difference and wrote zero raw and Barlow-filtered uncertainties.

## Changes

- `cut_variation_runner.cpp`
  - Normalizes complete extra wrapper-quote layers when reading CSV fields.
  - Rewrites fields with exactly one standards-compliant CSV escaping layer.
  - Makes tuple parsing tolerant of old malformed variation files.

- `cut_variation_systematics.cpp`
  - Applies the same wrapper-quote normalization to all input CSVs.
  - Performs a preflight validation of all four variation files.
  - Prints readable-bin counts for exclusivity loose/tight and fiducial loose/tight.
  - Requires every nominally populated cross-section bin to have a readable value in
    every variation.
  - Aborts before modifying the nominal CSV if validation fails.
  - Reports example missing or unreadable bin keys in the fatal error.
  - No longer substitutes zero differences for unreadable varied tuples.

## Expected startup summary

Before the RMS/Barlow calculation, output should include lines such as:

```text
[cut-systematics] exclusivity loose          valid bins: 1441 / 1441
[cut-systematics] exclusivity tight          valid bins: 1441 / 1441
[cut-systematics] fiducial loose             valid bins: 1441 / 1441
[cut-systematics] fiducial tight             valid bins: 1441 / 1441
```

The exact number is the number of nominally populated bins, so it may differ from 1441.
If any count is smaller than the nominal count, the main CSV is not modified.
