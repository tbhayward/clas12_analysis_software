# Photon-truth processing update

## Git branch

```bash
git checkout -b photon_efficiency
git push -u origin photon_efficiency
```

## New unified processors

The canonical scripts are:

- `processing_epgamma.groovy`
- `processing_epgammagamma.groovy`

The `processing_mc_*.groovy` copies are included only for compatibility with the
existing naming convention. All four scripts auto-detect MC event-by-event from
`MC::Particle`; separate data and MC logic is not required.

## Input arguments

Both new processors use:

```text
<hipo_dir> <output.txt> [n_files] [beam_energy] [run_override] [qadb_override] [max_output_rows]
```

`run_override <= 0` means use `RUN::config`.

`max_output_rows = 0` means unlimited.

## Recommended entry point

Use `processing.csh`, not the Groovy script directly. For the new photon
processors it automatically:

1. runs from the beginning until 1,000,000 candidate rows are written;
2. converts that text file immediately to `<output_base>_1M.root`;
3. restarts from the beginning;
4. processes all requested HIPO files;
5. converts the full file to `<output_base>.root`.

The 1M limit is **candidate rows**, i.e. ROOT-tree entries, rather than unique
HIPO event numbers. This is the most useful definition for getting a BDT-ready
sample quickly.

Example CLASDIS MC:

```bash
./processing.csh \
    processing_scripts/processing_epgamma.groovy \
    /path/to/clasdis/hipos \
    /scratch/$USER/clasdis_fa18_inb_epgamma_truth \
    0 10.604 11 1
```

Example data (`0` means no run override):

```bash
./processing.csh \
    processing_scripts/processing_epgamma.groovy \
    /path/to/data/hipos \
    /scratch/$USER/fa18_inb_epgamma_truth \
    0 10.604 0 0
```

The same calling convention applies to `processing_epgammagamma.groovy`.

## `processing_epgamma.groovy`

Each reconstructed photon candidate produces one row.

The tree preserves the reconstructed variables used in the existing DVCS
processing, and adds:

- `is_mc`
- reconstructed candidate indices / REC::Particle indices
- `matching_gamma_pid`
- `gamma_mcindex`
- `gamma_parent_index`, `gamma_parent_pid`
- `gamma_grandparent_index`, `gamma_grandparent_pid`
- analogous truth-match information for the electron and proton
- `gen_valid`
- generated counterparts of the reconstructed kinematic variables (`gen_*`)

Truth matching follows:

```text
REC::Particle pindex
    -> MC::RecMatch.mcindex
    -> MC::Particle.pid
    -> MC::Lund.parent
    -> parent MC::Lund.pid
```

The grandparent is saved as well. This is useful for a chain such as a gamma
whose immediate parent is a pi0 and whose pi0 itself came from a rho/omega/etc.

If `MC::Particle` is absent, `is_mc=0` and all truth/generated fields are
zero/sentinel/NaN as appropriate.

If `MC::Particle` exists but `MC::Lund` does not, the reconstructed-to-truth
PID can still be saved when `MC::RecMatch` exists, but parent/grandparent PID
remains zero.

## `processing_epgammagamma.groovy`

Every **unique** reconstructed photon pair is saved once:

```text
gamma2_index > gamma1_index
```

There is deliberately **no** `0.11 < M_gamma_gamma < 0.16 GeV` requirement.

The tree saves:

- electron, proton, gamma1, and gamma2 reconstructed kinematics;
- detector assignment of both photons;
- all FourParticles missing-mass combinations;
- `Mh_gammagamma` explicitly;
- single-photon `epgamma` topology variables separately for gamma1 and gamma2;
- truth PID, parent PID, and grandparent PID separately for both photons;
- generated FourParticles and single-photon topology variables when a complete
  truth match exists.

This makes the same file usable for pi0, eta, continuum, BDT, and tag-and-probe
photon-efficiency studies without rebuilding the two photons from a stored pi0.

## `convert_txt_to_root.cpp`

Legacy indices 0--8 are left on their existing hard-coded conversion paths.

New indices:

- `9` = epgamma truth schema
- `10` = epgammagamma truth schema

The new Groovy processors write a self-describing first line:

```text
#SCHEMA branch:I branch:D ...
```

For indices 9/10 only, the converter creates the ROOT branches from this
schema. This avoids adding hundreds of fragile positional `infile >> ...`
expressions and does not alter the old conversion behavior.

`beam_pol` and `target_pol` are still added by the converter using the same
run-information logic as the existing trees.

## Important caveats to validate on the first CLASDIS file

1. The parent lookup intentionally mirrors the existing
   `processing_mc_two_particles.groovy` convention that uses the
   `MC::RecMatch.mcindex` row to address `MC::Lund`. Check the first output's
   parent-PID spectrum before launching the full production.
2. A reconstructed calorimeter photon can truth-match to something other than
   PID 22 (for example an electron-associated shower). That is intentional and
   is exactly why `matching_gamma_pid` is saved separately from
   `gamma_parent_pid`.
3. `gen_valid=1` requires usable e, p, and gamma truth matches. A truth parent
   can therefore be available even for a row where the full generated
   Three/FourParticles kinematics could not be reconstructed.
4. AAOgen/DVCSgen may not carry useful `MC::Lund` ancestry. Those samples can
   still be processed with the same scripts, but parent PID should not be
   inferred when the bank is absent.
5. Before a day-scale production, run one small CLASDIS input and inspect:
   `matching_gamma_pid`, `gamma_parent_pid`, `gamma_grandparent_pid`,
   `gen_valid`, and the generated/reconstructed momentum closure.
