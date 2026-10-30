#!/usr/bin/env python3
"""Read-only, sampled readiness checks for paired CLAS12 DVCS and pi0 MC ROOT files.

Requires numpy and uproot. No ROOT files are modified.
"""
import argparse
from pathlib import Path
import sys

import numpy as np
import uproot

REQUIRED = (
    'event_sequence', 'gen_kinematics_valid', 'rec_kinematics_valid',
    'gen_Q2', 'rec_Q2', 'gen_x', 'rec_x', 'gen_t1', 'rec_t1',
)


def inspect_file(path, sample_size):
    print(f'\nFILE: {path}')
    with uproot.open(path) as root_file:
        if 'MCEvents' not in root_file:
            print('FAIL: missing MCEvents tree')
            return False
        tree = root_file['MCEvents']
        print(f'TREE: {tree.name}  ENTRIES: {tree.num_entries:,}  BRANCHES: {len(tree.keys())}')
        missing = sorted(set(REQUIRED) - set(tree.keys()))
        if missing:
            print('FAIL: missing branches:', ', '.join(missing))
            return False
        count = min(sample_size, tree.num_entries)
        if count == 0:
            print('WARNING: empty MCEvents tree')
            return False
        values = tree.arrays(list(REQUIRED), entry_stop=count, library='np')

    gen_raw = values['gen_kinematics_valid']
    rec_raw = values['rec_kinematics_valid']
    for label, raw in (('GEN', gen_raw), ('REC', rec_raw)):
        unique, counts = np.unique(raw, return_counts=True)
        print(f'{label} FLAG dtype={raw.dtype}, values={list(zip(unique.tolist(), counts.tolist()))}')
        if not np.all(np.isin(unique, [0, 1])):
            print(f'WARNING: {label} flag contains values other than 0 and 1')

    gen = np.isfinite(gen_raw) & (gen_raw == 1)
    rec = np.isfinite(rec_raw) & (rec_raw == 1)
    both = gen & rec
    print(f'SAMPLED: {count:,}  GEN VALID: {gen.sum():,}  REC VALID: {rec.sum():,}')
    print(f'BOTH VALID: {both.sum():,}  GEN-ONLY: {(gen & ~rec).sum():,}  REC-ONLY: {(rec & ~gen).sum():,}')

    seq = values['event_sequence']
    print('EVENT SEQUENCE first 10:', seq[:10].tolist())
    if not np.all(np.isfinite(seq)):
        print('FAIL: nonfinite event_sequence')
        return False
    diffs = np.diff(seq)
    print(f'SEQUENCE unique={len(np.unique(seq)) == len(seq)}, strictly_increasing={bool(np.all(diffs > 0))}, consecutive={bool(np.all(diffs == 1))}')
    if len(np.unique(seq)) != len(seq) or not np.all(diffs > 0):
        print('WARNING: event sequence not strictly increasing and unique in sample')

    for stem in ('Q2', 'x', 't1'):
        g = values['gen_' + stem]
        r = values['rec_' + stem]
        valid = both & np.isfinite(g) & np.isfinite(r)
        n = int(valid.sum())
        if n >= 3 and np.std(g[valid]) > 0 and np.std(r[valid]) > 0:
            corr = float(np.corrcoef(g[valid], r[valid])[0, 1])
            delta = r[valid] - g[valid]
            print(f'{stem}: matched_finite={n:,} correlation={corr:.6f} median(rec-gen)={np.median(delta):.6g} RMS(rec-gen)={np.sqrt(np.mean(delta**2)):.6g}')
        else:
            print(f'{stem}: insufficient jointly valid, finite entries ({n})')

    print('NOTE: correlations and sequence checks cannot prove correct event-level pairing; inspect converter pairing logic separately.')
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', default='/volatile/clas12/thayward/dvcs_pass2_full_rebuild/mc', help='MC campaign root directory')
    parser.add_argument('--sample', type=int, default=100000, help='Maximum entries per selected file')
    args = parser.parse_args()
    if args.sample < 1:
        parser.error('--sample must be positive')
    base = Path(args.root)
    paths = sorted(base.glob('**/parts/*.root'))
    print(f'MC ROOT FILES: {len(paths)}')
    if not paths:
        print('FAIL: no ROOT files found')
        return 1
    failures = 0
    for category in ('dvcs_mc', 'pi0_mc'):
        matches = [p for p in paths if category in p.parts]
        print(f'\nCATEGORY: {category}  FILES: {len(matches)}')
        if not matches:
            print('FAIL: no files found for category')
            failures += 1
            continue
        try:
            if not inspect_file(matches[0], args.sample):
                failures += 1
        except Exception as exc:
            print(f'FAIL: {type(exc).__name__}: {exc}')
            failures += 1
        #endif
    #endfor
    print(f'\nRESULT: {"PASS (sampled structural checks)" if failures == 0 else f"FAIL ({failures} categories)"}')
    return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(main())
