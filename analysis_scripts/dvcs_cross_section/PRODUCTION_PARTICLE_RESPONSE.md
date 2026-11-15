# Optional production particle response

The nominal invocation remains `./dvcs_analysis`. Particle response is enabled
explicitly with `--particle-response-run FITTED_RUN`. This first production
integration **imports a frozen fitted response**; it does not silently retune
particle parameters during extraction. The existing joint selected-response
study remains the fitting/validation workflow. No broadening-only fit is added.

## Farm build and first extraction

From `analysis_scripts/dvcs_cross_section`, after committing/pushing these changes:

```sh
git pull && make -j10 dvcs_analysis && python3 test_production_particle_response.py --native && python3 run_response_analysis.py --response-run "$PWD/test_runs/selected_particle_response_001" --run-dir "$PWD/test_runs/production_response_001" --scratch-root /scratch/thayward/dvcs_particle_response/production --workers 10 --skip-systematics --plot
```

The launcher creates a new private working directory, snapshots executables and
small source/model files, and shares imports read-only. It refuses an existing
run destination. Nominal outputs and earlier studies are preserved. Independent
run directories may run on different farm nodes using the same repository/build.
Build once before starting those jobs; running jobs use executable snapshots.
The scratch root is configurable and may be node-local. Same-node jobs use
per-file locks and can share complete caches without overwriting one another.

`--plot` runs the final-results-by-|t| script with Pass-1 autodetection enabled,
KM15/theory and measured reference means, using its usual phi sampling default.
To seed valid theory cache entries, add `--theory-cache-from /path/to/previous/run`.
Coordinates/settings still determine cache validity; changed means can require
new evaluations. Use a fresh `--run-dir` for each extraction.

The equivalent direct executable option is:

```sh
./dvcs_analysis --particle-response-run /absolute/path/to/fitted/run --particle-response-scratch /scratch/thayward/dvcs_particle_response/production --particle-response-workers 10 --skip-systematics
```

Run directly only in a private working directory with the usual scripts/imports.
The launcher is recommended. Reusing a response output without its response flag
is rejected, preventing nominal MC from being combined with response-derived
cuts/CSVs. Diagnostic-only executable modes are not accepted together with this
full-extraction option. Automatic downstream systematics inside a full response
run receive the modified maps normally; `--skip-systematics` omits those runs.

## Order in main.cpp

1. Determine the initial 90/95/98% production template windows using nominal MC.
   Save them under `output/particle_response/initial_cuts`. Initial detailed
   development canvases are suppressed to avoid duplicating final canvases.
2. Open the ordinary DATA, generated and reconstructed maps and write their
   explicit input inventory. Import the frozen response, checking the two shared
   response-core hashes and completed parameter checksum where available.
3. Materialize full-statistics reconstructed MC caches for the reference DVCS,
   pi0 leakage, direct two-photon control, all current scans, and any loaded
   reconstructed radiative inputs. Identical reference/current files are built
   only once. No event rows are removed. DATA and generated files/maps are
   unchanged, including generated denominators and MC weights.
4. Replace every reconstructed C++ map and install the same manifest for Python.
   Redetermine all three production cut containments from those modified inputs.
   Save final windows under `particle_response/final_cuts` and install the usual
   `output/jsons/combined_cuts*.json` consumed throughout production.
5. Generate the response note package, then continue through bin means, current
   calibration, production counts, pi0 subtraction, acceptance, normalization and
   all requested systematics. Current calibration therefore uses the same
   response as the acceptance numerator. The existing DATA proton-efficiency
   prescription is retained.
6. Associate the completed CSV with its model/manifest and final-cut hashes in
   `particle_response/extraction_manifest.json` and the CSV-adjacent
   `dvcs_pass2_analysis.response.json`.

There is one explicit response/cut pass, without an uncontrolled outer fit loop.
Imported Frad/Fbin conventions and the bin-volume calculation are unchanged.

## Event calculation and preservation

The same 20 fitted parameters describe electron FD, proton FD/CD, and photon
FD/FT. Momentum changes follow `p' = p*(1+bias+sigma_fraction*Z)` and angle changes
use independent Gaussian widths in radians. Bias and widths remain independent.
Each photon's own detector tag selects its response. Photon roles are not sorted
or exchanged; identical detector/kinematics use identical laws in either role,
with independent fluctuations.

A compiled persistent batch evaluator calls the existing fitted response laws.
A production extension updates xB, Q2, W, y, transfer variables, analysis azimuths
(including phi2), opening angles, missing quantities, pair/particle quantities,
depolarization factors and the other supported physics branches. Circular
quantities wrap correctly. Original stored calibration/rounding offsets are
preserved by adding the calculated change; zero response reproduces original
physics values exactly. The fitted exclusivity helper remains authoritative for
its original columns. An unknown physics branch causes a clear failure rather
than silently retaining a stale value.

For direct pi0, the individual photons are recovered using the same method as the
fit, modified independently, and the composite momentum is rebuilt. Composite
pi0 kinematics retain the fitter/analyzer's fixed-pi0-mass convention; reported
measured mgg is `saved_mgg + rebuilt_modified_mgg - rebuilt_original_mgg`.
Unresolved controls retain the **whole original event**, matching the fit; their
counts are recorded. Nonphysical modified events remain as rows and are marked
invalid with NaN matching kinematics, so selection rejects them consistently.

Event sequence/IDs, helicity, vertices, detector assignments, generator fields
and weights never undergo a floating-point round trip. They retain their source
types/values. Deterministic draws reproduce the fit's source-path/entry/seed keys,
independent of batching and worker count. Input files must not be renamed and
then assumed to have the same random realization.

Input trees already encode processing thresholds/fiducials and the two-photon
mass window. This stage cannot restore events absent from those trees or model
changes in reconstruction finding efficiency. It modifies the available
reconstructed response. Neither nominal physics inputs nor raw trees are edited.

## Selection and normalization conventions

The joint imported fit may use unordered two-photon validation groups. Production
continues to use its established three topology keys, with the first-photon tag
for direct pi0. Python final cut determination and C++ production counting use
that same convention; diagnostic four-group cut JSONs are never substituted into
three-key consumers. This distinction is explicit in the note report.

The exploratory fitted background nuisance is **not** installed as a new event
weight or efficiency factor. Production uses its existing binwise control/leakage
transfer and current/proton prescriptions, now with consistently modified MC and
final selections. Note overlays use raw common-selection control normalization,
sampling probabilities and the saved exposure ratio. They are shape diagnostics,
not substitutes for final weighted cross sections.

## Performance and files

Native calculations operate in bounded 16,384-event batches; ROOT files stream
through uproot with nested worker pools disabled and lossless LZ4 compression.
At most ten file workers run. No per-event Python response calculation is used.
Photon recovery, random keys and original derived values are reused within each
event. Unused Gaussian coordinates are not evaluated. Large event caches stay
on scratch, indexed by inputs, parameters, seed, evaluator and note-sampling settings. Global-cut
changes reuse the same modified four-vectors; saved diagnostic frames allow
new selections without rebuilding the response.
Changing the run directory alone does not invalidate a compatible cache.

Complete per-file outputs are installed atomically under a file lock and checked
by source identity, output identities and entry counts. Failures never install a
partial dataset. The generated data and reconstructed cache have separate paths;
no generated denominator is derived from the modified reconstructed files.

The initial cut finder caches unchanged DATA arrays on scratch. The final finder
and note stage reuse those arrays with matching sources/global cuts/reader code.
Downstream current/count/systematics stages read the cached modified ROOT trees
through their normal interfaces; DATA and generated trees remain original. Building modified MC
is a one-time full-statistics pass, which can be substantial on the first run;
compatible subsequent extractions reuse it. Moving a small frozen-model snapshot
does not force a rebuild; changing one period’s parameters does not rebuild other
periods’ caches.

## Analysis-note outputs

Open:

`test_runs/production_response_001/output/particle_response/analysis_note/REPORT.html`

The package contains vector PDF and 220-dpi PNG figures, grouped by period and
DVCS/two-photon channel, plus one combined `particle_response_figures.pdf`:

* Signed momentum biases and additional momentum/polar/azimuth widths.
* Conditional score influence of each of the 20 individual terms, obtained by
  setting one term to zero without refitting other parameters.
* Training-variable and other-MC-half scores plus reserved missing-mass checks.
* Four-panel DATA/nominal mixture/modified mixture comparisons of coplanarity,
  missing transverse momentum, missing energy and missing mass, including
  modified signal/background components. The displayed exclusivity cut is
  omitted, with common final DATA windows on the other variables.
* Direct-control counterparts and explicit two-photon mass comparisons.

The bounded note sample defaults to an expected **50,000 events per period and
reference MC category**, before selection. All samples/scenarios are paired. The
production extraction itself always uses full MC statistics. Bar scores depend
on note-sample precision and are not formal p-values or independent DATA tests.
Influences are conditional and not additive; no parameter uncertainty is implied.
Mass discrepancies remain visible rather than being concealed by the mixture.

Small tables and `plot_histograms.npz` preserve the parameters, sampled
normalizations, component histograms, scores and support. To change presentation
without reading any ROOT files or large event caches:

```sh
python3 plot_production_particle_response.py --replot --output "$PWD/test_runs/production_response_001/output/particle_response/analysis_note"
```

The existing production acceptance, contamination, current-response and final
cross-section note outputs are also produced by their normal downstream stages.
Transfer the private run's `output` for assessment, retaining large caches on the
farm. The private run also snapshots the small fitted response for reproducibility.

## Validation

`python3 test_production_particle_response.py --native` checks real ROOT streaming,
exact zero replay, independently calculated changed DIS/transfer/azimuth/missing
quantities, all note scenarios, batch invariance, uint64 identity and weight
preservation, own-detector photon response and role symmetry, unresolved records,
unknown physics branches, deduplication/cache reuse, moved-model reuse, shared
routing, DATA/GEN preservation, the actual production DATA reader/cache, note
rendering and histogram-only replot. Changed C++ translation units are checked
separately when developing. These tests do not claim a full farm-data extraction;
that is the first production command above.
