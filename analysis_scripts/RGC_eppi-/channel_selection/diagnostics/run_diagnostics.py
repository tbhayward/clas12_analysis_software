#!/usr/bin/env python3
"""Stream six period ROOT trees into exclusivity histograms and PNG diagnostics."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import hashlib
import html
import json
import os
from pathlib import Path
import time
import numpy as np
from kinematics import BRANCHES, METRICS, CORRELATIONS, XB, TP, calculate, bin_indices

DEFAULT_OUTPUT = Path(__file__).resolve().parent/'output'
VERSION = 1
DEFAULT_PERIODS = ['RGB_spring2019_D2', 'RGB_fall2019_D2', 'RGB_spring2020_D2',
                   'RGC_summer2022_ND3', 'RGC_fall2022_ND3', 'RGC_spring2023_ND3']
EDGES = {k: np.linspace(v[1], v[2], 121) for k, v in METRICS.items()}
CATS = ['all', 'FD', 'CD']

def atomic_json(path, value):
    tmp = path.with_suffix('.tmp')
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
    tmp.replace(path)

def process(job):
    import uproot
    name, source, out, options = job
    source, out = Path(source), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    stat = source.stat()
    signature = dict(version=VERSION, source=str(source.resolve()), size=stat.st_size,
                     mtime_ns=stat.st_mtime_ns, max_entries=options['max_entries'],
                     all_fiducials=options['all_fiducials'], tree=options['tree'],
                     code_sha256=hashlib.sha256(Path(__file__).read_bytes()+Path(__file__).with_name('kinematics.py').read_bytes()).hexdigest())
    manifest = out/'manifest.json'
    if not options['force'] and manifest.exists() and (out/'histograms.npz').exists():
        previous = json.loads(manifest.read_text())
        if previous['signature'] == signature:
            print(f'CACHED {name}: {previous["counts"]["in_bins"]:,} combinations', flush=True)
            return name
    manifest.unlink(missing_ok=True)
    (out/'plot_status.json').unlink(missing_ok=True)
    h = {k: np.zeros((3, 24, 120), dtype=np.int64) for k in METRICS}
    # Per metric/category/bin: finite count, sum, sum of squares, underflow, overflow, nonfinite.
    stats = {k: np.zeros((3, 24, 6)) for k in METRICS}
    corr = {f'{x}__{y}': np.zeros((24, 60, 60), dtype=np.int64) for x, y in CORRELATIONS}
    counts = dict(read=0, nonfinite_input=0, quality_rejected=0, outside_bins=0, in_bins=0)
    max_mx2_difference = 0.
    started = time.monotonic()
    atomic_json(out/'progress.json', dict(stage='opening input', period=name, pid=os.getpid(), host=os.uname().nodename))
    print(f'START {name}: opening input (PID {os.getpid()})', flush=True)
    with uproot.open(source) as f:
        tree = f[options['tree']]
        absent = set(BRANCHES)-set(tree.keys())
        if absent:
            raise ValueError(f'{name}: missing branches {sorted(absent)}')
        total = min(tree.num_entries, options['max_entries'] or tree.num_entries)
        atomic_json(out/'progress.json', dict(stage='reading first chunk', period=name, entries_read=0, entries_total=total, percent=0.0))
        for a in tree.iterate(BRANCHES, library='np', step_size=options['step_size'], entry_stop=total):
            n = len(a['x'])
            counts['read'] += n
            finite = np.logical_and.reduce([np.isfinite(a[k]) for k in BRANCHES])
            counts['nonfinite_input'] += int((~finite).sum())
            quality = (a['detector2'] == 1) & np.isin(a['detector1'], [1, 2]) & (a['Q2'] >= 1) & (a['W'] >= 2) & (a['y'] <= .8)
            if not options['all_fiducials']:
                quality &= a['fiducial_status'] == 3
            counts['quality_rejected'] += int((finite & ~quality).sum())
            bi = bin_indices(a['x'], -a['tprime'])
            counts['outside_bins'] += int((finite & quality & (bi < 0)).sum())
            keep = finite & quality & (bi >= 0)
            a = {k: v[keep] for k, v in a.items()}
            bi = bi[keep]
            counts['in_bins'] += len(bi)
            values = calculate(a)
            if len(bi):
                diff = float(np.max(np.abs(values['mx2_n']-a['Mx2'])))
                max_mx2_difference = max(max_mx2_difference, diff)
                if diff > 1e-6:
                    raise ValueError(f'{name}: saved and neutron-reconstructed Mx2 disagree by {diff:g} GeV^2')
            for k, v in values.items():
                edges = EDGES[k]
                vi = np.searchsorted(edges, v, side='right')-1
                for ci in range(3):
                    sel = np.ones(len(bi), dtype=bool) if ci == 0 else a['detector1'] == ci
                    good = sel & np.isfinite(v)
                    inside = good & (vi >= 0) & (vi < 120)
                    h[k][ci] += np.bincount(bi[inside]*120+vi[inside], minlength=24*120).reshape(24, 120)
                    for j, mask, weight in [(0, good, None), (1, good, v), (2, good, v*v),
                                             (3, good & (v < edges[0]), None),
                                             (4, good & (v >= edges[-1]), None), (5, sel & ~np.isfinite(v), None)]:
                        stats[k][ci, :, j] += np.bincount(bi[mask], weights=None if weight is None else weight[mask], minlength=24)
            for x, y in CORRELATIONS:
                xe, ye = EDGES[x][::2], EDGES[y][::2]
                ix, iy = np.searchsorted(xe, values[x], side='right')-1, np.searchsorted(ye, values[y], side='right')-1
                ok = np.isfinite(values[x]) & np.isfinite(values[y]) & (ix >= 0) & (ix < 60) & (iy >= 0) & (iy < 60)
                corr[f'{x}__{y}'] += np.bincount(bi[ok]*3600+ix[ok]*60+iy[ok], minlength=24*3600).reshape(24, 60, 60)
            elapsed = time.monotonic()-started
            progress = dict(stage='processing', period=name, entries_read=counts['read'], entries_total=total,
                            percent=100*counts['read']/max(total, 1), selected_combinations=counts['in_bins'], elapsed_seconds=elapsed)
            atomic_json(out/'progress.json', progress)
            print(f'{name}: {progress["percent"]:5.1f}%  {counts["read"]:,}/{total:,} rows, {counts["in_bins"]:,} in bins ({elapsed:.0f}s)', flush=True)
        manifest_value = dict(signature=signature, counts=counts, source_entries=tree.num_entries,
                              partial=total < tree.num_entries, max_mx2_difference=max_mx2_difference,
                              xb_edges=XB.tolist(), minus_tprime_edges=TP.tolist(),
                              metric_definitions=METRICS, elapsed_seconds=time.monotonic()-started)
    with (out/'histograms.tmp').open('wb') as stream:
        np.savez_compressed(stream, **{f'h_{k}': v for k, v in h.items()},
                            **{f's_{k}': v for k, v in stats.items()}, **{f'c_{k}': v for k, v in corr.items()})
    (out/'histograms.tmp').replace(out/'histograms.npz')
    atomic_json(manifest, manifest_value)
    with (out/'statistics.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['metric', 'proton_detector', 'bin', 'xB_low', 'xB_high', 'minus_tprime_low', 'minus_tprime_high', 'finite', 'mean', 'std', 'underflow', 'overflow', 'nonfinite'])
        for k in METRICS:
            for ci, category in enumerate(CATS):
                for b in range(24):
                    n, s, ss, lo, hi, bad = stats[k][ci, b]
                    mean = s/n if n else 0
                    sd = np.sqrt(max(0, ss/n-mean*mean)) if n else 0
                    writer.writerow([k, category, b+1, XB[b//6], XB[b//6+1], TP[b%6], TP[b%6+1], int(n), mean if n else '', sd if n else '', int(lo), int(hi), int(bad)])
    print(f'DONE {name}: {counts["in_bins"]:,} combinations', flush=True)
    return name

def plot_period(out, log_y=False):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    meta = json.loads((out/'manifest.json').read_text())
    name = out.name
    atomic_json(out/'plot_status.json', dict(stage='plotting', started=time.time()))
    print(f'PLOTTING {name}', flush=True)
    note = 'PARTIAL SAMPLE' if meta['partial'] else 'Full sample'
    with np.load(out/'histograms.npz') as data:
        for xi in range(4):
            title = f'{name} | {XB[xi]:.2f} ≤ xB < {XB[xi+1]:.2f} | {note}'
            for k, definition in METRICS.items():
                fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
                for ti, ax in enumerate(axes.flat):
                    b = xi*6+ti
                    for ci, (color, label) in enumerate(zip(['black', 'tab:blue', 'tab:orange'], CATS)):
                        ax.stairs(data[f'h_{k}'][ci, b], EDGES[k], color=color, label=label, linewidth=1)
                    ax.axvline(definition[3], color='gray', linestyle='--', linewidth=1)
                    s = data[f's_{k}'][0, b]
                    ax.set_title(f'Bin {b+1}: {TP[ti]:.2f} ≤ −t′ < {TP[ti+1]:.2f}')
                    ax.text(.03, .96, f'N={int(s[0]):,}\noutside range={int(s[3]+s[4]):,}', transform=ax.transAxes, va='top', fontsize=9)
                    if log_y:
                        ax.set_yscale('log')
                        ax.set_ylim(bottom=.8)
                    else:
                        ax.set_ylim(bottom=0)
                    ax.set_xlabel(definition[0])
                    ax.set_ylabel('Candidate combinations / bin')
                axes.flat[0].legend(loc='upper right', fontsize=8)
                fig.suptitle(title)
                fig.savefig(out/f'{k}_xb{xi+1}.png', dpi=140)
                plt.close(fig)
            for x, y in CORRELATIONS:
                fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
                for ti, ax in enumerate(axes.flat):
                    z = data[f'c_{x}__{y}'][xi*6+ti]
                    if z.max() > 0:
                        artist = ax.pcolormesh(EDGES[x][::2], EDGES[y][::2], np.ma.masked_less(z.T, 1), norm=LogNorm(vmin=1, vmax=max(2, int(z.max()))), rasterized=True)
                        fig.colorbar(artist, ax=ax, label='Combinations')
                    else:
                        ax.text(.5, .5, 'No entries in displayed range', transform=ax.transAxes, ha='center')
                    ax.set_title(f'Bin {xi*6+ti+1}: {TP[ti]:.2f} ≤ −t′ < {TP[ti+1]:.2f}')
                    ax.set_xlabel(METRICS[x][0])
                    ax.set_ylabel(METRICS[y][0])
                fig.suptitle(title+' | all protons')
                fig.savefig(out/f'correlation_{x}_{y}_xb{xi+1}.png', dpi=140)
                plt.close(fig)
    atomic_json(out/'plot_status.json', dict(stage='complete', finished=time.time(), png_pages=60))
    print(f'PLOTTED {name}: 60 PNG pages', flush=True)

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--datasets', type=Path, default=Path(__file__).resolve().parents[2]/'imports/datasets.json')
    p.add_argument('--input-root', type=Path)
    p.add_argument('--output', type=Path, default=DEFAULT_OUTPUT, help='Output folder (default: output/ beside this script)')
    p.add_argument('--periods', nargs='+', default=DEFAULT_PERIODS)
    p.add_argument('--workers', type=int, default=6, help='Total processes, one input file per process')
    p.add_argument('--step-size', default='100 MB')
    p.add_argument('--max-entries', type=int, help='Read only first N candidate rows per period for a smoke test')
    p.add_argument('--all-fiducials', action='store_true', help='Include stored candidates failing fiducial cuts')
    p.add_argument('--force', action='store_true', help='Rebuild matching caches')
    p.add_argument('--plots-only', action='store_true', help='Replot existing complete caches without opening ROOT files')
    p.add_argument('--log-y', action='store_true')
    p.add_argument('--status', action='store_true', help='Show saved progress and exit')
    args = p.parse_args()
    if args.status:
        for name in args.periods:
            folder = args.output/name
            if (folder/'manifest.json').exists():
                m = json.loads((folder/'manifest.json').read_text())
                plot_status = folder/'plot_status.json'
                state = json.loads(plot_status.read_text())['stage'] if plot_status.exists() else 'scan complete, plots pending'
                print(f'{name}: {state}, {m["counts"]["in_bins"]:,} combinations' + (' (PARTIAL)' if m['partial'] else ''))
            elif (folder/'progress.json').exists():
                v = json.loads((folder/'progress.json').read_text())
                age = time.time()-(folder/'progress.json').stat().st_mtime
                detail = f', {v["percent"]:.1f}%, {v["entries_read"]:,}/{v["entries_total"]:,} rows' if 'percent' in v else ''
                print(f'{name}: {v.get("stage", "processing")}{detail}, update {age:.0f}s ago')
            else:
                print(f'{name}: not started')
        return
    if args.workers < 1 or (args.max_entries is not None and args.max_entries < 1):
        p.error('workers and max-entries must be positive')
    args.output.mkdir(parents=True, exist_ok=True)
    print(f'Output directory: {args.output.resolve()}', flush=True)
    if not args.plots_only:
        config = json.loads(args.datasets.read_text())
        root = args.input_root or Path(config['campaign'])/'root'
        jobs = []
        for name in args.periods:
            if name not in config['files']:
                p.error(f'Unknown dataset {name}')
            source = root/config['files'][name]
            if not source.is_file():
                p.error(f'Missing input {source}')
            options = dict(max_entries=args.max_entries, all_fiducials=args.all_fiducials,
                           tree=config['tree'], step_size=args.step_size, force=args.force)
            jobs.append((name, str(source), str(args.output/name), options))
        print(f'Starting {min(args.workers, len(jobs))} parallel scan workers for {len(jobs)} periods', flush=True)
        with ProcessPoolExecutor(max_workers=min(args.workers, len(jobs))) as pool:
            futures = [pool.submit(process, job) for job in jobs]
            for future in as_completed(futures):
                future.result()
    for name in args.periods:
        plot_period(args.output/name, args.log_y)
    links = []
    # Include previously completed periods when hosts process disjoint subsets.
    for name in sorted(d.name for d in args.output.iterdir() if d.is_dir() and (d/'manifest.json').exists()):
        links.append(f'<h2>{html.escape(name)}</h2>')
        links.extend(f'<a href="{name}/{f.name}">{f.stem}</a><br>' for f in sorted((args.output/name).glob('*.png')))
    index_tmp = args.output/f'index.{os.getpid()}.tmp'
    index_tmp.write_text('<!doctype html><meta charset="utf-8"><title>Exclusive pion diagnostics</title><h1>Exclusivity diagnostics</h1><p>Candidate combinations. No exclusivity cuts or fits.</p>'+''.join(links))
    index_tmp.replace(args.output/'index.html')
    print(f'Browse {args.output/"index.html"}', flush=True)

if __name__ == '__main__':
    main()
