#!/usr/bin/env python3
"""Select missing-proton candidates and run the shared dilution and spin fits."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import uproot

HERE = Path(__file__).resolve().parent
CAMPAIGN = Path('/volatile/clas12/thayward/epi-/single/full_20261126')
sys.path.insert(0, str(HERE / 'channel_selection'))
from kinematics import bin_indices

BRANCHES = dict(runnum='int32', evnum='int32', helicity='int32',
                beam_pol='float64', target_pol='float64', DepA='float64',
                DepB='float64', DepC='float64', DepV='float64', DepW='float64',
                Mx2='float64', x='float64', tprime='float64', phi2='float64', analysis_bin='int32')


def select(job):
    name, source, output, cuts, cutset = job
    focused = cuts['mass_window'] is None
    edges = np.linspace(-1, 4, 251)
    spectra = np.zeros((24, len(edges)-1), dtype=np.int64)
    counts = np.zeros((24, 2), dtype=np.int64)
    runs = {}
    read = 0
    final = output / f'selected_{name}.root'
    temporary = final.with_suffix('.partial.root')
    columns = (set(BRANCHES) - {'analysis_bin', 'phi2'}) | {
        'phi', 'x', 'tprime', 'detector', 'pion_dc_fiducial'}
    with uproot.open(source) as incoming:
        tree = incoming['PhysicsEvents']
        entries = tree.num_entries
        from contextlib import nullcontext
        with (nullcontext() if focused else uproot.recreate(temporary)) as target:
            if target is not None:
                target.mktree('PhysicsEvents', BRANCHES)
            for a in tree.iterate(sorted(columns), library='np', step_size='64 MB'):
                read += len(a['x'])
                bins = bin_indices(a['x'], -a['tprime'])
                good = (bins >= 0) & (a['detector'] == 1) & np.isfinite(a['Mx2'])
                if cuts['pion_dc']:
                    good &= a['pion_dc_fiducial'] == 1
                if focused:
                    inside = good & (a['Mx2'] >= edges[0]) & (a['Mx2'] < edges[-1])
                    mass_bins = np.searchsorted(edges, a['Mx2'][inside], side='right')-1
                    spectra += np.bincount(bins[inside]*spectra.shape[1]+mass_bins,
                                           minlength=spectra.size).reshape(spectra.shape)
                    continue
                low, high = cuts['mass_window']
                keep = good & (a['Mx2'] >= low) & (a['Mx2'] < high)
                counts[:, 0] += np.bincount(bins[good], minlength=24)
                counts[:, 1] += np.bincount(bins[keep], minlength=24)
                for run, n in zip(*np.unique(a['runnum'], return_counts=True)):
                    runs.setdefault(str(run), {'source_combinations': 0})['source_combinations'] += int(n)
                if keep.any():
                    target['PhysicsEvents'].extend({
                        key: (bins[keep]+1 if key == 'analysis_bin'
                              else a['phi' if key == 'phi2' else key][keep])
                        for key in BRANCHES})
            if read != entries:
                raise ValueError(f'{name}: incomplete scan ({read}/{entries})')
            selected = int(counts[:, 1].sum())
            if target is not None:
                report = dict(dataset=name, fit_group='single_pion',
                              partial=False, source_entries=entries,
                              read_rows=read, selected_rows=selected,
                              counts=counts.tolist(), runs=runs)
                target['selection_manifest'] = json.dumps(report)
                target['selection_cutset_id'] = cutset
    if focused:
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(4, 6, figsize=(16, 9), sharex=True)
        for i, ax in enumerate(axes.flat):
            ax.stairs(spectra[i], edges)
            ax.axvline(.9382720813**2, color='tab:red', lw=1)
            ax.set_title(f'Bin {i+1}')
        fig.suptitle(name + r' — neutron-reference $M_X^2$ (GeV$^2$)')
        fig.tight_layout()
        fig.savefig(output / f'{name}_missing_mass.png', dpi=150)
        plt.close(fig)
        print(f'{name}: scanned {read:,} rows', flush=True)
        return None
    temporary.replace(final)
    print(f'{name}: selected {selected:,} of {read:,} rows', flush=True)
    return name, dict(path=final.name, selected_rows=selected)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-root', type=Path, default=CAMPAIGN/'root')
    parser.add_argument('--output-root', type=Path, default=CAMPAIGN/'analysis')
    parser.add_argument('--run-info', type=Path, default=HERE/'imports/clas12_run_info_ND3_source.csv')
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--pion-dc', action='store_true')
    parser.add_argument('--mass-window', type=float, nargs=2, metavar=('LOW', 'HIGH'))
    parser.add_argument('--stage', choices=['spectra','selection','dilution','asymmetry','all'],
                        default='spectra')
    args = parser.parse_args()
    if args.workers < 1:
        parser.error('--workers must be positive')
    if args.stage in {'selection', 'all'} and (
        args.mass_window is None or not np.isfinite(args.mass_window).all()
        or args.mass_window[0] >= args.mass_window[1]
    ):
        parser.error('Supply an increasing --mass-window LOW HIGH in GeV^2; inspect spectra first')
    out = args.output_root.resolve()
    selected = out / 'channel_selection'
    registry = selected / 'selected_inputs.json'
    dilution = out / 'dilution_factor/joint_period_method1.json'
    if args.stage in {'spectra', 'selection', 'all'}:
        datasets = json.loads((HERE/'imports/datasets.json').read_text())['files']
        sources = {name: args.input_root/f'epi_minus_{name}.root' for name in datasets}
        missing = [str(p) for p in sources.values() if not p.is_file()]
        if missing:
            parser.error('Dataset merges pending or missing:\n'+'\n'.join(missing))
        spectra = args.stage == 'spectra'
        output = out/'spectra' if spectra else selected
        output.mkdir(parents=True, exist_ok=True)
        cuts = dict(topology='single_pion', mass_window=None if spectra else args.mass_window,
                    pion_dc=args.pion_dc, detector=1)
        if not spectra:
            cuts['sources'] = {name: [str(path.resolve()), path.stat().st_size,
                                      path.stat().st_mtime_ns]
                               for name, path in sources.items()}
        cutset = hashlib.sha256(json.dumps(cuts, sort_keys=True).encode()).hexdigest()
        if not spectra and registry.exists():
            registry.unlink()
        jobs = [(name, path, output, cuts, cutset) for name, path in sources.items()]
        os.environ['MPLBACKEND'] = 'Agg'
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            results = list(pool.map(select, jobs))
        if not spectra:
            payload = dict(topology='single_pion', tree='PhysicsEvents', cuts=cuts,
                           cutset_id=cutset, files=dict(results))
            temp = registry.with_suffix('.tmp')
            temp.write_text(json.dumps(payload, indent=2)+'\n')
            temp.replace(registry)

    common = ['--selected-inputs', str(registry), '--run-info', str(args.run_info.resolve())]
    env = dict(os.environ, MPLBACKEND='Agg', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    if args.stage in {'all', 'dilution'}:
        subprocess.run([sys.executable, str(HERE/'dilution_factor/determine_dilution_factor.py'),
                        *common, '--output', str(dilution.parent),
                        '--joint-workers', str(args.workers)], check=True, env=env)
    if args.stage in {'all', 'asymmetry'}:
        subprocess.run([sys.executable, str(HERE/'asymmetry_extraction/extract_structure_function_ratios.py'),
                        *common, '--dilution', str(dilution), '--output',
                        str(out/'asymmetry_extraction'), '--workers', str(args.workers)],
                       check=True, env=env)


if __name__ == '__main__':
    main()
