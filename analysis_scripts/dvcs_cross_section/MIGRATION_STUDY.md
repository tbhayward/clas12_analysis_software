# Paired-MC migration study

The first milestone measures one-dimensional `|t|` migration using generated
and reconstructed quantities from the same `MCEvents` entry. It is diagnostic
only and does not alter the production cross-section correction.

The script uses `numpy`, `matplotlib`, and `uproot`; run it in the same Python
environment used by the production exclusivity-template script.

Run after the production exclusivity stage has written
`output/jsons/global_cuts_config.json` and `output/jsons/combined_cuts.json`:

```sh
python3 study_t_migration.py --workers 10 2>&1 | tee output/migration_study.log
```

For a quick branch/plot test using one paired file per period:

```sh
python3 study_t_migration.py --workers 5 --max-files-per-period 1
```

Outputs are written under `output/migration_study/t_response/`. Each period and
the combined result receives response-count, response-probability,
reconstructed-origin, quality-metric, residual, cutflow, reconstructed-topology,
and dedicated low-|t| feed-in plots, along with CSV and compressed numerical
arrays. `migration_response.root` stores response cells and metrics for
downstream ROOT inspection.

The next milestone will add the nominal multidimensional analysis-bin response
and independent train/test closure after these outputs have been checked.
