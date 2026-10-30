#!/usr/bin/env python3
"""Inspect reconstructed AAOgen topology assignments in paired MCEvents ROOT files."""
import argparse
import collections
import glob
from pathlib import Path
import numpy as np
import uproot

CURRENTS = {'sp18_inb': '50nA', 'sp18_out': '45nA', 'fa18_inb': '50nA', 'fa18_out': '50nA', 'sp19_inb': '50nA'}

def summarize(counter, title, limit=30):
    print(f'  {title}:', flush=True)
    for key, count in counter.most_common(limit):
        print(f'    {str(key):>24} : {count:,}', flush=True)
    if not counter:
        print('    (none)', flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--periods', nargs='+', default=['sp18_inb', 'sp18_out', 'fa18_inb'])
    p.add_argument('--root', default='/volatile/clas12/thayward/dvcs_pass2_full_rebuild')
    p.add_argument('--max-files', type=int, default=0, help='0 means all matching files')
    p.add_argument('--max-entries', type=int, default=5000000, help='per period; 0 means all')
    p.add_argument('--step-size', default='50 MB')
    args = p.parse_args()
    fields = ['rec_kinematics_valid', 'rec_detector1', 'rec_detector2', 'rec_detector_gamma1', 'rec_detector_gamma2']
    for period in args.periods:
        pattern = f'{args.root}/mc/{period}/pi0_mc/nominal/parts/*_{CURRENTS[period]}_*.root'
        paths = sorted(glob.glob(pattern))
        if args.max_files:
            paths = paths[:args.max_files]
        print(f'\n=== {period} === files={len(paths)}', flush=True)
        if not paths:
            print('  No matching files', flush=True)
            continue
        with uproot.open(paths[0]) as f:
            missing = [x for x in fields if x not in f['MCEvents'].keys()]
        if missing:
            print('  Missing branches:', missing, flush=True)
            continue
        combos = collections.Counter()
        pairs = collections.Counter()
        photons = collections.Counter()
        d1s = collections.Counter()
        d2s = collections.Counter()
        n_all = n_valid = n_pair_match = 0
        for chunk in uproot.iterate([f'{x}:MCEvents' for x in paths], expressions=fields, library='np', step_size=args.step_size):
            n = len(chunk[fields[0]])
            if args.max_entries and n_all + n > args.max_entries:
                n = args.max_entries - n_all
            if n <= 0:
                break
            n_all += n
            good = np.asarray(chunk['rec_kinematics_valid'][:n]) >= 0.5
            n_valid += int(np.count_nonzero(good))
            vals = [np.asarray(chunk[x][:n])[good] for x in fields[1:]]
            if not len(vals[0]):
                continue
            d1, d2, g1, g2 = [np.rint(v).astype(np.int32) for v in vals]
            n_pair_match += int(np.count_nonzero(d2 == g1))
            for value, count in zip(*np.unique(np.stack((d1, d2, g1, g2), axis=1), axis=0, return_counts=True)):
                combos[tuple(map(int, value))] += int(count)
            for a, b, target in [(d1, d2, pairs), (g1, g2, photons)]:
                for value, count in zip(*np.unique(np.stack((a, b), axis=1), axis=0, return_counts=True)):
                    target[tuple(map(int, value))] += int(count)
            for values, target in [(d1, d1s), (d2, d2s)]:
                for value, count in zip(*np.unique(values, return_counts=True)):
                    target[int(value)] += int(count)
            if args.max_entries and n_all >= args.max_entries:
                break
        print(f'  paired entries scanned: {n_all:,}; reconstructed valid: {n_valid:,} ({n_valid/max(n_all,1):.5%})', flush=True)
        print(f'  rec_detector2 == rec_detector_gamma1: {n_pair_match:,}/{n_valid:,}', flush=True)
        summarize(d1s, 'proton detector1')
        summarize(d2s, 'composite detector2')
        summarize(pairs, 'saved topology (detector1, detector2)')
        summarize(photons, 'individual photon detector pair (gamma1, gamma2)')
        summarize(combos, 'full (proton, composite, gamma1, gamma2) combinations')
        fd_fd = sum(v for (d, _, g1, g2), v in combos.items() if d == 1 and g1 == 1 and g2 == 1)
        cd_fd = sum(v for (d, _, g1, g2), v in combos.items() if d == 2 and g1 == 1 and g2 == 1)
        cd_ft = sum(v for (d, _, g1, g2), v in combos.items() if d == 2 and (g1 == 0 or g2 == 0))
        print(f'  photon-derived candidates: FD-FD={fd_fd:,} CD-FD={cd_fd:,} CD-FT(with >=1 FT photon)={cd_ft:,}', flush=True)

if __name__ == '__main__':
    main()
