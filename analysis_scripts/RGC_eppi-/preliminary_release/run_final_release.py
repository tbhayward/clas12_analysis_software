#!/usr/bin/env python3
"""Preliminary detected-proton release from the completed loose analysis."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
import json
import os
import shutil
from pathlib import Path
import sys

for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
os.environ['MPLBACKEND'] = 'Agg'

import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
for folder in ('asymmetry_extraction', 'asymmetry_extraction/three_bin_study',
               'dilution_factor', 'channel_selection'):
    sys.path.insert(0, str(BASE/folder))

from inputs import load, PERIODS
from kinematics import XB, TP
from likelihood import fit, NAMES
from joint_period_fit import replica_bin

TP_GROUPS = ((0,), (1,), (2, 3), (4, 5))


def slices(groups, names, members):
    result = []
    for name in names:
        group = groups[name]
        events = group['events']
        keep = np.isin(events['bin'], members)
        if keep.any():
            result.append(dict(events={key: value[keep] for key, value in events.items()},
                               states=group['states']))
    if not result:
        raise ValueError(f'No events in {names}, fine bins {members}')
    return result


def binned_check(groups):
    """Fit 12 helicity-asymmetry points using their event-averaged sensitivity."""
    events = {key: np.concatenate([g['events'][key] for g in groups])
              for key in ('phi', 'pb', 'rw', 'h')}
    phi = events['phi']
    index = np.minimum((phi*12/(2*np.pi)).astype(int), 11)
    count = np.bincount(index, minlength=12)
    signed = np.bincount(index, weights=events['h'], minlength=12)
    sensitivity = np.bincount(index, weights=events['pb']*events['rw']*np.sin(phi),
                              minlength=12)
    use = (count > 1) & (np.abs(signed) < count)
    n = count[use]
    y, x = signed[use]/n, sensitivity[use]/n
    variance = (1-y*y)/n
    information = np.sum(x*x/variance)
    if information <= 0 or use.sum() < 2:
        return dict(status='invalid')
    value = np.sum(x*y/variance)/information
    return dict(status='valid', value=float(value), stat=float(information**-.5),
                chi2=float(np.sum((y-x*value)**2/variance)), ndf=int(use.sum()-1))


def save(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True, help='Completed loose analysis directory')
    p.add_argument('--output', type=Path, default=BASE/'output_preliminary_final')
    p.add_argument('--run-info', type=Path, default=BASE/'imports/clas12_run_info_ND3_source.csv')
    p.add_argument('--workers', type=int, default=4)
    p.add_argument('--replicas', type=int, default=2000)
    a = p.parse_args()
    if a.workers < 1 or a.replicas < 2:
        p.error('Require positive workers and at least two replicas')
    source = a.source.resolve()
    out = a.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    groups, _, audit, identity = load(
        source/'channel_selection/selected_inputs.json', a.run_info,
        source/'dilution_factor/joint_period_method1.json', BASE/'channel_selection')
    material = json.loads((source/'dilution_factor/method1.json').read_text())
    if material['cutset_id'] != identity['cutset_id']:
        raise ValueError('Material counts and selected events have different cut sets')
    xb = [list(range(i*6+1, i*6+7)) for i in range(4)]
    tp = [[i*6+j+1 for i in range(4) for j in columns] for columns in TP_GROUPS]
    projections = {'xB': (xb, XB.tolist()), 'minus_tprime': (tp, TP[[0,1,2,4,6]].tolist())}
    fine = []
    comparisons = []
    for b in range(1, 25):
        sample = slices(groups, PERIODS[:3], [b])
        result = fit(sample, (1., 0.), fixed_uu=(0., 0.))
        check = binned_check(sample)
        comparisons.append(dict(bin=b, likelihood=result, chi2=check))
        if result['status'] == 'valid':
            fine.append(dict(scope='RGB', axis='2D', bin=b, observable=NAMES[2],
                             value=result['values'][2], stat=result['errors'][2],
                             events=result['events']))
    internal = out/'internal'
    internal.mkdir(exist_ok=True)
    save(internal/'rgb_fit_comparison.json', comparisons)
    diagnostics = out/'exclusivity'
    diagnostics.mkdir(exist_ok=True)
    from selection import load_selected_inputs
    registry = load_selected_inputs(source/'channel_selection/selected_inputs.json')
    cut_file = Path(registry['cut_file'])
    shutil.copy2(cut_file, diagnostics/cut_file.name)
    for image in sorted((cut_file.parent/'fits').glob('*/*.png')):
        destination = diagnostics/image.parent.name
        destination.mkdir(exist_ok=True)
        shutil.copy2(image, destination/image.name)
    rows = list(fine)
    for axis, (members, edges) in projections.items():
        counts = [dict(item, selected_counts=[
            sum(item['selected_counts'][b-1] for b in ids) for ids in members])
                  for item in material['inputs']]
        jobs = [(dict(inputs=counts), i, a.replicas, 20261128) for i in range(4)]
        with ProcessPoolExecutor(max_workers=a.workers) as pool:
            results = list(pool.map(replica_bin, jobs))
        dilution = [row for row, _ in results]
        period, transfer = dilution_checks(counts, a.replicas)
        save(out/f'dilution_{axis}.json', dict(edges=edges, inputs=counts, bins=dilution,
                                              periods=period, transfer=transfer))
        fits = []
        for b, (ids, df) in enumerate(zip(members, dilution), 1):
            sample = slices(groups, PERIODS[3:], ids)
            result = fit(sample, (df['df'], df['stat_std']), fixed_uu=(0., 0.))
            fits.append(dict(bin=b, fine_bins=ids, **result))
            if result['status'] == 'valid':
                for j in range(2, 7):
                    rows.append(dict(scope='RGC', axis=axis, bin=b, observable=NAMES[j],
                                     value=result['values'][j], stat=result['errors'][j],
                                     events=result['events']))
        save(out/f'fits_{axis}.json', fits)
    save(out/'provenance.json', dict(identity=identity, samples=audit,
        xB_edges=XB.tolist(), fine_tprime_edges=TP.tolist(), projections=projections))
    with (out/'measurements.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=['scope','axis','bin','observable','value','stat','events'])
        writer.writeheader()
        writer.writerows(rows)
    plot(out, rows, projections)


def dilution_checks(inputs, replicas):
    from model import method1, transfer_ratio

    rng = np.random.default_rng(20261128)
    rates, draws = {}, {}
    for item in inputs:
        _, period, target = item['dataset'].split('_')
        counts = np.asarray(item['selected_counts'])
        charge = item['charge']
        rates.setdefault(period, {})[target] = counts/charge
        draws.setdefault(period, {})[target] = rng.poisson(counts, (replicas, len(counts)))/charge

    def summarize(values, samples):
        return [dict(bin=i+1, value=float(v) if np.isfinite(v) else None,
                     stat=float(np.std(valid, ddof=1)) if len(valid)>1 else None,
                     finite_replicas=len(valid))
                for i, (v, column) in enumerate(zip(values, samples.T))
                for valid in [column[np.isfinite(column)]]]

    transfer = transfer_ratio(rates['spring2023'])
    transfer_draws = transfer_ratio(draws['spring2023'])
    periods = {}
    for name in rates:
        measured = name == 'spring2023'
        nominal = method1(rates[name], None if measured else transfer)
        sampled = method1(draws[name], None if measured else transfer_draws)
        periods[name] = summarize(nominal['dilution'], sampled['dilution'])
    return periods, summarize(transfer, transfer_draws)


def plot(out, rows, projections):
    import matplotlib.pyplot as plt
    from plots import mark

    labels = [r'$F_{LU}^{\sin\phi}/F_{UU}$', r'$F_{UL,\mathrm{lab}}^{\sin\phi}/F_{UU}$',
              r'$F_{UL,\mathrm{lab}}^{\sin 2\phi}/F_{UU}$', r'$F_{LL,\mathrm{lab}}/F_{UU}$',
              r'$F_{LL,\mathrm{lab}}^{\cos\phi}/F_{UU}$']

    def points(ax, sample, centers):
        ax.axhline(0, color='.65', lw=1)
        ax.errorbar([centers[r['bin']] for r in sample], [r['value'] for r in sample],
                    yerr=[r['stat'] for r in sample], fmt='o', color='black', capsize=3)
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
    for i, ax in enumerate(axes.flat):
        centers = {i*6+j+1: (TP[j]+TP[j+1])/2 for j in range(6)}
        points(ax, [r for r in rows if r['scope']=='RGB' and r['bin'] in centers], centers)
        ax.set(ylim=(-.4,.4), ylabel=labels[0], title=rf'${XB[i]:.2f}<x_B<{XB[i+1]:.2f}$', xlabel=r"$-t'$ (GeV$^2$)")
    fig.suptitle(r'RGB deuterium · $F_{LU}^{\sin\phi}/F_{UU}$')
    mark(fig)
    fig.tight_layout()
    fig.savefig(out/'RGB_LU.png', dpi=300)
    plt.close(fig)

    for axis, (_, edges) in projections.items():
        centers = {i+1: (lo+hi)/2 for i, (lo, hi) in enumerate(zip(edges, edges[1:]))}
        fig, axes = plt.subplots(2, 3, figsize=(12, 7), sharex=True)
        for j, (ax, name) in enumerate(zip(axes.flat, NAMES[2:])):
            points(ax, [r for r in rows if r['axis']==axis and r['observable']==name], centers)
            ax.set(ylabel=labels[j], ylim=(-.4,.4) if j<3 else (-1,1), xlabel=r'$x_B$' if axis=='xB' else r"$-t'$ (GeV$^2$)")
        axes.flat[-1].set_visible(False)
        fig.suptitle(r'RGC polarized ND$_3$')
        mark(fig)
        fig.tight_layout()
        fig.savefig(out/f'RGC_{axis}.png', dpi=300)
        plt.close(fig)

        diagnostics = json.loads((out/f'dilution_{axis}.json').read_text())
        dilution = diagnostics['bins']
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.errorbar([centers[r['bin']] for r in dilution], [r['df'] for r in dilution],
                    yerr=[r['stat_std'] for r in dilution], fmt='o', color='black', capsize=3)
        ax.set(xlabel=r'$x_B$' if axis=='xB' else r"$-t'$ (GeV$^2$)",
               ylabel=r'Dilution factor $D_f$', ylim=(0, 1))
        fig.tight_layout()
        fig.savefig(out/f'dilution_{axis}.png', dpi=300)
        plt.close(fig)

        xlabel = r'$x_B$' if axis=='xB' else r"$-t'$ (GeV$^2$)"
        for suffix, series, ylabel in (
            ('periods', diagnostics['periods'], r'Dilution factor $D_f$'),
            ('transfer', {'Spring 2023': diagnostics['transfer']}, 'CD₂ / CH₂ charge-normalized yield'),
        ):
            fig, ax = plt.subplots(figsize=(6, 4), layout='constrained')
            for marker, (label, values) in zip(('o','s','^'), series.items()):
                valid = [r for r in values if r['value'] is not None and r['stat'] is not None]
                ax.errorbar([centers[r['bin']] for r in valid],
                            [r['value'] for r in valid], yerr=[r['stat'] for r in valid],
                            fmt=marker, capsize=3, label=label)
            ax.set(xlabel=xlabel, ylabel=ylabel)
            ax.legend()
            fig.savefig(out/f'dilution_{axis}_{suffix}.png', dpi=300)
            plt.close(fig)


if __name__ == '__main__':
    main()
