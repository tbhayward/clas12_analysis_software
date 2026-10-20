#!/usr/bin/env python3
"""Pass-2 DVCS final-results canvases, grouped by |t| and gridded in Q2/xB.

Place alongside plot_pass2_vs_pass1_model_comparison.py.
Example:
 python3 plot_pass2_final_results_by_t.py imports/dvcs_pass2_analysis.csv \
   --output-dir output/final_results_by_t --phi-dense 73

Requires existing dvcsgen and km15_cli.py installations used by the comparison code.
A dedicated persistent JSON model cache avoids reevaluation on subsequent runs.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import csv
import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import plot_pass2_vs_pass1_model_comparison as base


def finite(x):
    return math.isfinite(x)


def model_cache_key(k, args, kind):
    return json.dumps([k.xb_min, k.xb_max, k.q2_min, k.q2_max,
                       k.t_min, k.t_max, args.e_beam, args.phi_dense,
                       args.dvcsgen_dir, args.km15_cli, args.python_km15, kind])


def compute_cell_models(job):
    """One process per kinematic cell; compute only missing observables."""
    k, cfg, kinds = job
    phi = [360.0 * i / (cfg.phi_dense - 1) for i in range(cfg.phi_dense)]
    xb = (k.xb_min + k.xb_max) / 2
    q2 = (k.q2_min + k.q2_max) / 2
    t = (k.t_min + k.t_max) / 2
    result = {}
    if 'xs' in kinds:
        result['xs'] = {'phi': phi, 'km15': [], 'bh': []}
    if 'bsa' in kinds:
        result['bsa'] = {'phi': phi, 'km15': [], 'bh': [0.0] * len(phi)}
    for angle in phi:
        if 'xs' in kinds:
            result['xs']['bh'].append(base.bh_xs(xb, q2, t, angle, cfg))
            result['xs']['km15'].append(base.km15_xs(xb, q2, t, angle, cfg, helicity=0))
        if 'bsa' in kinds:
            plus = base.km15_xs(xb, q2, t, angle, cfg, helicity=+1)
            minus = base.km15_xs(xb, q2, t, angle, cfg, helicity=-1)
            denom = plus + minus
            result['bsa']['km15'].append((plus - minus) / denom if math.isfinite(denom) and denom != 0 else float('nan'))
    return k, result


def save_cache_atomic(path, cache):
    temp = path.with_name(path.name + '.tmp')
    temp.write_text(json.dumps(cache, separators=(',', ':')))
    os.replace(temp, path)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('csv', type=Path, help='Pass-2 extracted 4D CSV')
    ap.add_argument('--output-dir', type=Path, default=Path('output/final_results_by_t'))
    ap.add_argument('--phi-dense', type=int, default=73)
    ap.add_argument('--e-beam', type=float, default=10.6)
    ap.add_argument('--dvcsgen-dir', default=base.DEFAULT_DVCSGEN_DIR)
    ap.add_argument('--km15-cli', default=base.DEFAULT_KM15_CLI)
    ap.add_argument('--python-km15', default=sys.executable)
    ap.add_argument('--xs-column', default=None)
    ap.add_argument('--xs-total-ptp', action='store_true', help='XS error bars: stat and point-to-point systematics in quadrature')
    ap.add_argument('--bsa-total-ptp', action='store_true', help='BSA error bars: stat and point-to-point systematics in quadrature, if available')
    ap.add_argument('--recompute-models', action='store_true')
    ap.add_argument('--dpi', type=int, default=180)
    ap.add_argument('--workers', type=int, default=8, help='Parallel model workers (1–8; default 8)')
    args = ap.parse_args()
    if not 1 <= args.workers <= 8:
        ap.error('--workers must be between 1 and 8')
    if args.phi_dense < 2:
        ap.error('--phi-dense must be >=2')
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with args.csv.open(newline='') as f:
        reader = csv.DictReader(f)
        fields = reader.fieldnames or []
        cols = base.detect_common_columns(fields, 'Pass-2')
        xs_col = args.xs_column or base.find_column(fields, base.PASS2_CENTRAL_XS_CANDIDATES)
        bsa_col = base.find_column(fields, [base.PASS2_BSA_TARGET_COLUMNS['10.6 GeV']])
        xs_ptp_col = base.find_column(fields, [base.PASS2_POINT_TO_POINT_TOTAL_COLUMN], required=False)
        bsa_ptp_col = base.find_column(fields, ['BSA, counts, 10.6 GeV, point-to-point total sys', 'BSA, counts, 10.6 GeV, point-to-point sys', 'BSA, counts, 10.6 GeV, total point-to-point sys'], required=False)
        if args.bsa_total_ptp and bsa_ptp_col is None:
            raise ValueError('No recognized BSA point-to-point systematic column; use statistical errors or extend column candidates.')
        panels = {}
        for row in reader:
            if not base.row_is_valid(row, cols['valid']):
                continue
            key = base.make_bin_key(row, cols)
            phi = base.row_phi(row, cols)
            if not finite(phi):
                continue
            item = panels.setdefault(key, {'xs': [], 'bsa': []})
            for observable, column, sys_col, enabled in [
                ('xs', xs_col, xs_ptp_col, args.xs_total_ptp),
                ('bsa', bsa_col, bsa_ptp_col, args.bsa_total_ptp),
            ]:
                v = base.parse_tuple_cell(row.get(column, ''))
                if not v.ok or not finite(v.value) or not finite(v.stat) or v.stat <= 0:
                    continue
                if observable == 'xs' and v.value <= 0:
                    continue
                err = v.stat
                if enabled and sys_col:
                    sys_err = base.to_float(row.get(sys_col, ''), default=0.0)
                    if finite(sys_err):
                        err = math.hypot(err, sys_err)
                item[observable].append((phi, v.value, err))
    panels = {k: v for k, v in panels.items() if v['xs'] or v['bsa']}
    if not panels:
        raise RuntimeError('No valid Pass-2 results found in CSV.')

    cfg = base.ModelConfig(e_beam=args.e_beam, dvcsgen_dir=args.dvcsgen_dir,
                           km15_cli=args.km15_cli, py_km15=args.python_km15,
                           phi_dense=args.phi_dense, use_cache=True,
                           allow_missing_models=False, verbose_worker_models=False)
    cache_path = args.output_dir / 'model_curves_cache.json'
    cache = json.loads(cache_path.read_text()) if cache_path.exists() and not args.recompute_models else {}
    groups = {}
    for k in panels:
        groups.setdefault((k.t_min, k.t_max), []).append(k)
    print(f'[results] {len(panels)} populated kinematic cells; {len(groups)} |t| bins', flush=True)
    # Calculate missing model curves once, in parallel across independent cells.
    # The plotting phase is deliberately serial: matplotlib is not process-safe
    # and parallel PNG rendering would waste memory and initialization time.
    jobs = []
    for k in panels:
        kinds = tuple(kind for kind in ('xs', 'bsa')
                      if panels[k][kind] and model_cache_key(k, args, kind) not in cache)
        if kinds:
            jobs.append((k, cfg, kinds))
    if jobs:
        print(f'[models] Computing {len(jobs)} cells using up to {args.workers} workers', flush=True)
        with ProcessPoolExecutor(max_workers=min(args.workers, len(jobs))) as pool:
            futures = [pool.submit(compute_cell_models, job) for job in jobs]
            try:
                for i, future in enumerate(as_completed(futures), 1):
                    k, curves = future.result()
                    for kind, entry in curves.items():
                        cache[model_cache_key(k, args, kind)] = entry
                    # Checkpoint occasionally, not after every curve/phi point.
                    if i % 8 == 0 or i == len(jobs):
                        save_cache_atomic(cache_path, cache)
                    print(f'[models] {i}/{len(jobs)} cells complete', flush=True)
            except BaseException:
                save_cache_atomic(cache_path, cache)
                raise
    else:
        print('[models] All required curves loaded from cache', flush=True)

    for it, (trange, keys) in enumerate(sorted(groups.items()), start=1):
        xbins = sorted(set((k.xb_min, k.xb_max) for k in keys))
        qbins = sorted(set((k.q2_min, k.q2_max) for k in keys))
        print(f'[results] |t| {trange}: {len(keys)} cells, {len(xbins)} xB columns, {len(qbins)} Q2 rows', flush=True)
        for kind in ('xs', 'bsa'):
            fig, axs = plt.subplots(len(qbins), len(xbins), figsize=(4.0 * len(xbins), 2.85 * len(qbins)), squeeze=False, sharex=True, sharey=(kind == 'bsa'))
            fig.subplots_adjust(left=0.075, right=0.985, bottom=0.10, top=0.87, hspace=0.20, wspace=0.20)
            fig.suptitle(rf'CLAS12 Pass-2 $ep\to ep\gamma$ — {"cross section" if kind == "xs" else "beam-spin asymmetry"}' + '\n' + rf'$|t|\in[{trange[0]:.3f},\,{trange[1]:.3f}]\ \mathrm{{GeV}}^2$, $E_{{beam}}={args.e_beam:g}$ GeV', fontsize=15, y=0.985)
            for iq, qb in enumerate(qbins):
                for ix, xb in enumerate(xbins):
                    ax = axs[iq, ix]
                    matches = [k for k in keys if (k.xb_min, k.xb_max) == xb and (k.q2_min, k.q2_max) == qb]
                    if not matches:
                        ax.axis('off')
                        continue
                    k = matches[0]
                    dat = sorted(panels[k][kind])
                    if not dat:
                        ax.axis('off')
                        continue
                    phi, yy, ee = np.array(dat).T
                    ax.errorbar(phi, yy, yerr=ee, fmt='o', ms=3.5, lw=0.9, capsize=1.5, color='#1f77b4', label='CLAS12 Pass-2', zorder=4)
                    cache_key = model_cache_key(k, args, kind)
                    curves = cache[cache_key]
                    for model, color, name in [('km15', '#d62728', 'KM15'), ('bh', '#2ca02c', 'BH')]:
                        vv = np.asarray(curves[model], dtype=float)
                        xx = np.asarray(curves['phi'], dtype=float)
                        valid = np.isfinite(vv) & np.isfinite(xx)
                        if np.any(valid):
                            ax.plot(xx[valid], vv[valid], '-', lw=1.65, color=color, label=name, zorder=2)
                    ax.set_xlim(0, 360)
                    ax.set_xticks([0, 90, 180, 270, 360])
                    ax.grid(alpha=0.18)
                    ax.text(0.035, 0.96, rf'$x_B=[{xb[0]:.3f},{xb[1]:.3f}]$' + '\n' + rf'$Q^2=[{qb[0]:.2f},{qb[1]:.2f}]\ \mathrm{{GeV}}^2$', transform=ax.transAxes, va='top', ha='left', fontsize=8, bbox={'facecolor': 'white', 'alpha': 0.78, 'edgecolor': 'none', 'pad': 2})
                    if kind == 'bsa':
                        ax.set_ylim(-0.65, 0.65)
                    if iq == len(qbins)-1:
                        ax.set_xlabel(r'$\phi$ (deg)')
                    if ix == 0:
                        ax.set_ylabel(r'$d^4\sigma$ (analysis units)' if kind == 'xs' else r'$A_{LU}$')
            from matplotlib.lines import Line2D
            handles = [Line2D([], [], color='#1f77b4', marker='o', ls='none', markersize=5, label='CLAS12 Pass-2'), Line2D([], [], color='#d62728', lw=2, label='KM15'), Line2D([], [], color='#2ca02c', lw=2, label='BH')]
            fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 0.914), ncol=3, frameon=False)
            output = args.output_dir / f'pass2_{kind}_tbin_{it:02d}_{trange[0]:.3f}_{trange[1]:.3f}.png'
            fig.savefig(output, dpi=args.dpi)
            plt.close(fig)
            print(f'[write] {output}', flush=True)
    print('[done] All requested PNG canvases generated.', flush=True)


if __name__ == '__main__':
    main()
