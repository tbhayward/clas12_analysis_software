# Additional proton-efficiency correction

Production now defaults to no additional Neupane proton DATA/MC efficiency correction. The particle-response fit, MC acceptance, current response, radiative corrections and bin-centering corrections remain active. The optional eppi0 efficiency normalization also remains off by default.

`--no-neupane-proton-correction` explicitly selects the default. `--neupane-proton-correction` restores the former weighting for comparison. Neupane and the eppi0 MC efficiency map cannot be enabled together.

The chosen Neupane mode is propagated to nominal DVCS and two-photon counts, Pass-1 and both Sebastian binning checks and all their automatic exclusivity/fiducial variations. The startup log states the mode. The per-bin proton-efficiency diagnostic factors are unity when disabled.

Use a fresh output directory for the new baseline. Reusing old corrected counts or systematic products is not a valid comparison. A full run uses `./dvcs_analysis --particle-response --no-neupane-proton-correction --all-bsa-binning-cross-checks --particle-response-workers 10`; omitting `--skip-systematics` runs all configured systematic studies. Large response intermediates belong under the configured scratch root.
