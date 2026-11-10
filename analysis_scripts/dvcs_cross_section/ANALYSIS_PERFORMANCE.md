# Result-preserving analysis speedups

These changes preserve event order, cuts, model grids, numerical precision,
weights, fitting settings and histogram/CSV/plot content. They complement the
existing fused current scans and shared generated-count scans. They do not
skip any analysis stage. The earlier user-requested 24-point plotting default
is unchanged by this update.

## 1. Batch KM15 model requests

`km15_cli.py --batch` accepts a JSON list of scalar argument lists on stdin
and returns a versioned JSON response. It imports gepard/KM15 once, then calls
the existing scalar evaluator for every request. Every request still uses
its original .12g input formatting, .12g scalar output formatting, helicity,
phi grid and BSA arithmetic. No interpolation, approximation or persistent
model-value cache is added.

The final-by-t plotter batches all requested KM15 helicities for one cell.
At 24 phi points, XS+BSA now uses one model startup instead of 72. The shared
Pass-2/Pass-1 comparison helper also batches its XS and BSA curve requests.
Older custom CLIs automatically fall back to their original scalar path.
Existing numerical model caches remain valid.

Local real-gepard check: all nine scalar/batch predictions matched exactly;
6.305 seconds scalar versus 0.661 seconds batch. This is a model-startup
benchmark, not an end-to-end farm speedup estimate. BH evaluation is unchanged.

For an exact farm comparison, `DVCS_VERIFY_KM15_BATCH=1` checks every batched
value against its original scalar call (intentionally slower). To force the
old plotting path, set `DVCS_KM15_SCALAR=1`.

## 2. Prepare current-calibration exclusivity windows once per scan

The current-dependence event loops previously built a topology/cut key and
looked up seven named windows for every event. Each branch binder now holds
an invocation-local prepared window set for each topology. The cache resets
on binding and on changes to cut-map identity, channel prefix or period.
Cut order, optional-cut behavior, nonfinite handling, fallback bounds,
upper-quantile handling, missing-branch diagnostics and pi0 angle mapping
are preserved. There is no cross-run cache and no change to accumulation
order or fit inputs.

Tests compare the actual prepared/reference dispatch functions across
100,000 events plus 10,000 randomized window configurations, boundary values,
nonfinite values, channel/topology changes and missing-branch failures.
A one-million-call local fixture took approximately 0.12 seconds with the
reference dispatch and 0.012 seconds prepared. This isolates cut evaluation;
it does not represent total calibration speedup.

Set `DVCS_VERIFY_CURRENT_CUTS=1` to additionally check every optimized decision
against the original implementation during a real extraction.

## 3. Restore bounded ROOT read-ahead for the bin-means pass

Bin means previously disabled ROOT's read cache. It now uses a 16 MiB cache
per tree worker, explicitly listing only the already-enabled branches and
stopping automatic cache learning. This batches reads on farm storage while
retaining every event and the original summation order. At five simultaneous
periods the cache budget is about 80 MiB. No extra scan is added.

Set `DVCS_BIN_MEANS_NO_CACHE=1` to restore uncached reads for a timing or
result comparison. The regression script builds three small ROOT files and
compares every branch value and accumulated sum across cached/uncached
TChain file transitions.

Both modified C++ units passed local ROOT-header syntax checks. The local
ROOT runtime has incompatible precompiled system headers and crashes before
the synthetic chain test can execute; that runtime test must be run in the
farm ROOT environment. No full farm extraction has been benchmarked locally.

## Build and check

From the suite directory, in the farm ROOT/GCC environment:

```sh
make -j10 dvcs_analysis && python3 test_analysis_performance.py --real-km15
```

`--real-km15` uses the Python environment running the test and requires its
gepard installation. Omit it to run the deterministic model fixture tests
without gepard. The default test also runs the ROOT-chain cache check;
`--skip-root` explicitly disables that check for environments without a
working ROOT runtime. No farm event files are needed for these checks.

Normal run commands are unchanged. Recompile once before launching jobs.
Previously launched isolated runs retain their snapshotted code and executables.
