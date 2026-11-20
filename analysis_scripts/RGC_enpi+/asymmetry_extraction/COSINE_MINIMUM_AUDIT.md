# Saved-toy minimum audit

Run `audit_cosine_minima.py` from the RGC parent directory. Default: audit every existing bin-6 trial. Use `--audit-bins 4 5 6` to include all three bins. This does not alter production fits or original validation checkpoints. Finish the original study before taking the final audit snapshot; rerunning the audit will also pick up newly saved original trials.

On either/both farm nodes (shared original and audit output directories):

```sh
env TMPDIR=/scratch/thayward OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 python3 asymmetry_extraction/audit_cosine_minima.py --workers 10 --distributed --cache /scratch/thayward/rgc_enpi_thayward/selected_events_5774b50cd1143836.npz
```

After both nodes finish:

```sh
python3 asymmetry_extraction/audit_cosine_minima.py --summarize-only
```

Output: `asymmetry_extraction/output/cosine_minimum_audit/`. Same recovery/bias/coverage PNG families and CSV summaries as the original grid, using audited trials only. `tables/minimum_comparison.csv` provides paired old/new coefficients, usability, likelihood improvement and unresolved lower-profile flags. Seven rows per trial: do not sum likelihood improvements across parameter rows.

Reconstructs the exact original trial from saved spin seeds and deterministic thinning/auxiliary seeds; checks auxiliary centers and saved-fit likelihood before proceeding. Uses original inputs, limits and nuisance treatment. Original study/support/production hashes, runtime and inputs must still match. Audit settings and code are separately pinned; use a different audit output directory if changing them.

Search uses the saved fit, baseline nominal, all-zero physics and eight reproducible random physics starts independent of injected truth. All seven physics parameters and applicable dilution nuisances float. Each start runs MIGRAD, SIMPLEX, then strategy-2 MIGRAD. The lowest finite likelihood is retained even if its covariance fails quality checks; a higher minimum is not substituted just because its errors look better. HESSE is evaluated at the selected solution. Original lower solutions are retained if a retry is worse.

Both true-value profile checks use the stronger search. If a constrained search finds a better likelihood, its parameters seed a fully released fit and both profiles are repeated (at most two promotion rounds). Remaining inconsistencies are explicitly excluded and flagged. The injected value is used only in the diagnostic constrained profile, not in independent unrestricted starts. These checks improve the search but cannot prove a global minimum.

Default eight extra starts means substantially more work per trial; checkpoint locks are per replica so two nodes can share work. Completed replicas are reused on restart. No new toys are generated and no truth-dependent subset of trials is selected: all saved trials in the requested bins are audited. Recovery means remain conditional on usable fits; report unusable and profile-failure fractions alongside coverage.

Interpretation: if paired likelihood improvements remove bias/undercoverage, optimizer failures contributed. If recovery problems persist despite consistent minima, investigate finite-sample behavior, parameter bounds and weak information. This remains a conditional test at observed kinematics, not an independent detector-acceptance test.
