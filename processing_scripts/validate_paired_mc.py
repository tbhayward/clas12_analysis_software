#!/usr/bin/env python3
"""Validate paired generated/reconstructed MCEvents trees from DVCS and pi0 pilots.

Requires PyROOT, numpy and matplotlib. Produces CSV summaries and PNG plots only.
"""
import argparse
import csv
import glob
import os
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CORE = ('Q2', 'x', 't1', 'phi2')
PI0_EXTRA = ('Mh_gammagamma', 'detector_gamma1', 'detector_gamma2',
             'open_angle_egamma1', 'open_angle_egamma2', 'gamma_phi1',
             'gamma_phi2', 'theta_pi0_pi0')


def read_columns(path, names, ROOT):
    frame = ROOT.RDataFrame('MCEvents', path)
    return {k: np.asarray(v, dtype=float) for k, v in frame.AsNumpy(list(names)).items()}


def save_plot(path, x, y, title, xlab, ylab, kind):
    fig, ax = plt.subplots(figsize=(6.8, 5.0))
    if kind == 'hist':
        x = x[np.isfinite(x)]
        y = y[np.isfinite(y)]
        if len(x) or len(y):
            combined = np.concatenate((x, y))
            lo, hi = np.nanpercentile(combined, [0.5, 99.5])
            if not np.isfinite(lo) or not np.isfinite(hi) or lo == hi:
                lo, hi = float(np.nanmin(combined))-1, float(np.nanmax(combined))+1
            bins = np.linspace(lo, hi, 61)
            ax.hist(x, bins=bins, histtype='step', label=f'generated ({len(x)})')
            ax.hist(y, bins=bins, histtype='step', label=f'reconstructed ({len(y)})')
            ax.legend()
    elif kind == 'migration':
        mask = np.isfinite(x) & np.isfinite(y)
        if np.any(mask):
            x, y = x[mask], y[mask]
            lo, hi = np.nanpercentile(np.concatenate((x, y)), [0.5, 99.5])
            if lo == hi: lo, hi = lo-1, hi+1
            ax.hist2d(x, y, bins=55, range=[[lo, hi], [lo, hi]], cmap='viridis')
            ax.plot([lo, hi], [lo, hi], '--', color='white', linewidth=1)
    elif kind == 'residual':
        z = (y-x)[np.isfinite(x) & np.isfinite(y)]
        if len(z):
            lo, hi = np.percentile(z, [0.5, 99.5])
            if lo == hi: lo, hi = lo-1, hi+1
            ax.hist(z, bins=np.linspace(lo, hi, 61), histtype='step')
    ax.set(title=title, xlabel=xlab, ylabel=ylab)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base', required=True, help='Pilot fa18_inb directory containing dvcs_rec and pi0_rec')
    p.add_argument('--output', default='paired_mc_validation_v2')
    args = p.parse_args()
    import ROOT
    ROOT.gROOT.SetBatch(True)
    os.makedirs(args.output, exist_ok=True)
    files = sorted(glob.glob(os.path.join(args.base, 'dvcs_rec', '**', '*.mc_events.root'), recursive=True) +
                   glob.glob(os.path.join(args.base, 'pi0_rec', '**', '*.mc_events.root'), recursive=True))
    if not files:
        p.error('No *.mc_events.root files found beneath --base')
    rows = []
    for i, path in enumerate(files):
        mode = 'pi0' if '/pi0_rec/' in path else 'dvcs'
        f = ROOT.TFile.Open(path)
        tree = f.Get('MCEvents') if f else None
        if not tree:
            print('ERROR: missing MCEvents:', path)
            continue
        names = {b.GetName() for b in tree.GetListOfBranches()}
        gen_fields = {n[4:] for n in names if n.startswith('gen_')}
        rec_fields = {n[4:] for n in names if n.startswith('rec_')}
        physics = gen_fields & rec_fields
        expected = 78 if mode == 'pi0' else 71
        missing = [f'{prefix}_{field}' for prefix in ('gen', 'rec') for field in CORE + (PI0_EXTRA if mode == 'pi0' else ()) if f'{prefix}_{field}' not in names]
        needed = {'gen_kinematics_valid', 'rec_kinematics_valid'} | {f'{side}_{k}' for side in ('gen', 'rec') for k in CORE}
        missing += sorted(needed - names)
        row = {'file': path, 'mode': mode, 'entries': int(tree.GetEntries()),
               'gen_prefixed': len(gen_fields), 'rec_prefixed': len(rec_fields),
               'paired_common_fields': len(physics), 'expected_physics_fields': expected,
               'missing_branches': ';'.join(sorted(set(missing))), 'status': 'OK'}
        if missing:
            row['status'] = 'MISSING_BRANCHES'
            print(f'[{i+1}/{len(files)}] {mode} MISSING BRANCHES: {path}: {missing}')
            rows.append(row)
            f.Close()
            continue
        # Include all paired physics variables, so we can test the complete schema.
        pairs = sorted((gen_fields & rec_fields) - {'kinematics_valid'})
        keys = {'gen_kinematics_valid', 'rec_kinematics_valid'} | {f'{side}_{k}' for side in ('gen', 'rec') for k in pairs}
        values = read_columns(path, keys, ROOT)
        gv = values['gen_kinematics_valid'] > 0.5
        rv = values['rec_kinematics_valid'] > 0.5
        both = gv & rv
        row.update(gen_valid=int(gv.sum()), rec_valid=int(rv.sum()), both_valid=int(both.sum()),
                   gen_only=int((gv & ~rv).sum()), rec_only=int((rv & ~gv).sum()),
                   neither=int((~gv & ~rv).sum()),
                   gen_valid_fraction=float(gv.mean()), rec_valid_fraction=float(rv.mean()))
        for side, mask in [('gen', gv), ('rec', rv)]:
            for k in CORE:
                v = values[f'{side}_{k}']
                row[f'{side}_{k}_finite_valid'] = int(np.isfinite(v[mask]).sum())
                row[f'{side}_{k}_nonfinite_valid'] = int((~np.isfinite(v[mask])).sum())
                row[f'{side}_{k}_finite_invalid'] = int(np.isfinite(v[~mask]).sum())
            row[f'{side}_any_nonfinite_valid_fields'] = sum(int((~np.isfinite(values[f'{side}_{k}'][mask])).sum()) for k in CORE)
        # Only numeric physics quantities should be tested for finite values; a few
        # bookkeeping fields are intentionally sentinel-valued at truth level.
        for k in CORE:
            gx = values[f'gen_{k}'][gv]
            rx = values[f'rec_{k}'][rv]
            prefix = os.path.join(args.output, f'{i:02d}_{mode}_{Path(path).stem.replace(".mc_events", "")}')
            save_plot(prefix + f'_{k}_distribution.png', gx, rx,
                      f'{mode.upper()} generated vs reconstructed: {k}', k, 'Events', 'hist')
        if np.any(both):
            for k in ('Q2', 't1'):
                gx = values[f'gen_{k}'][both]
                rx = values[f'rec_{k}'][both]
                row[f'{k}_median_rec_minus_gen'] = float(np.nanmedian(rx-gx))
                row[f'{k}_rms_rec_minus_gen'] = float(np.nanstd(rx-gx))
                prefix = os.path.join(args.output, f'{i:02d}_{mode}_{Path(path).stem.replace(".mc_events", "")}')
                save_plot(prefix + f'_{k}_migration.png', gx, rx,
                          f'{mode.upper()} migration: {k}', f'generated {k}', f'reconstructed {k}', 'migration')
                save_plot(prefix + f'_{k}_residual.png', gx, rx,
                          f'{mode.upper()} residual: {k}', f'reconstructed - generated {k}', 'Events', 'residual')
        # Legacy PhysicsEvents contains selected reconstructed candidates, which
        # need not map 1:1 to MCEvents entries. Report counts; do not assert an
        # eventwise match without an explicit common event index.
        legacy_path = path.replace('.mc_events.root', '.root')
        row['legacy_file_exists'] = os.path.isfile(legacy_path)
        row['legacy_entries'] = ''
        if row['legacy_file_exists']:
            lf = ROOT.TFile.Open(legacy_path)
            if lf:
                for key in lf.GetListOfKeys():
                    obj = lf.Get(key.GetName())
                    if obj and obj.InheritsFrom('TTree'):
                        row['legacy_entries'] = int(obj.GetEntries())
                        break
                lf.Close()
        if any(row[f'{side}_{k}_nonfinite_valid'] for side in ('gen','rec') for k in CORE):
            row['status'] = 'NONFINITE_VALID_CORE'
        if any(row[f'rec_{k}_finite_invalid'] for k in CORE):
            row['status'] = 'REC_VALUES_WHEN_INVALID'
        print(f'[{i+1:02d}/{len(files):02d}] {mode:4s} N={row["entries"]:6d} gen={row["gen_valid"]:6d} rec={row["rec_valid"]:6d} both={row["both_valid"]:6d} fields={len(physics):3d} legacy={row["legacy_entries"]} {row["status"]}')
        rows.append(row)
        f.Close()
    columns = sorted(set().union(*(r.keys() for r in rows)))
    with open(os.path.join(args.output, 'summary.csv'), 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    print('\nSummary:', os.path.join(args.output, 'summary.csv'))
    print('PNG plots:', args.output)
    print('NOTE: Legacy entry counts are reported, but eventwise comparison requires an explicit shared event index.')


if __name__ == '__main__':
    main()
