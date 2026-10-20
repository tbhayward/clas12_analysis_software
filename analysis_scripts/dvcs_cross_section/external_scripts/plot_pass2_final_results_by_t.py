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

# Avoid rewriting tracked .pyc files on the ifarm when importing local helpers.
# This must precede the local module import.
sys.dont_write_bytecode = True
import plot_pass2_vs_pass1_model_comparison as base


def finite(x):
    return math.isfinite(x)


def model_cache_key(k, args, kind):
    return json.dumps([k.xb_min, k.xb_max, k.q2_min, k.q2_max,
                       k.t_min, k.t_max, args.e_beam if args.e_beam is not None else getattr(args, "_active_energy", None), args.phi_dense,
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
    ap.add_argument('--e-beam', type=float, choices=(10.2, 10.6), default=None, help='Restrict to one beam energy (default: both)')
    ap.add_argument('--xs-column-10p2', default=None, help='Override the Sp19 10.2 GeV cross-section column')
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
        xs_col = args.xs_column or base.find_column(fields, base.PASS2_CENTRAL_XS_CANDIDATES, required=False)
        xs_10p2_candidates = [
            'normed cross sections, ep->epg, exp, Sp19 Inb, unpol',
            'cross sections, ep->epg, exp, Sp19 Inb, unpol',
            'normed cross sections, ep->epg, exp, 10.2 GeV, unpol',
            'cross sections, ep->epg, exp, 10.2 GeV, unpol',
        ]
        xs_10p2 = args.xs_column_10p2 or base.find_column(fields, xs_10p2_candidates, required=False)
        if xs_10p2 is not None and xs_10p2 not in fields:
            raise ValueError(f'10.2 GeV XS column not found: {xs_10p2}')
        if xs_col is None and (args.e_beam is None or args.e_beam == 10.6):
            raise ValueError('No 10.6 GeV cross-section column found')
        if xs_10p2 is None and (args.e_beam is None or args.e_beam == 10.2):
            print('[warning] No 10.2 GeV cross-section column detected; 10.2 GeV BSA will still be plotted. Use --xs-column-10p2 to specify the cross-section column.', flush=True)
        bsa_col = base.find_column(fields, [base.PASS2_BSA_TARGET_COLUMNS['10.6 GeV']], required=False)
        bsa_10p2 = base.find_column(fields, [base.PASS2_BSA_TARGET_COLUMNS['Sp19 Inb']], required=False)
        xs_ptp_col = base.find_column(fields, [base.PASS2_POINT_TO_POINT_TOTAL_COLUMN], required=False)
        bsa_ptp_col = base.find_column(fields, ['BSA, counts, 10.6 GeV, point-to-point total sys', 'BSA, counts, 10.6 GeV, point-to-point sys', 'BSA, counts, 10.6 GeV, total point-to-point sys'], required=False)
        if args.bsa_total_ptp and bsa_ptp_col is None:
            raise ValueError('No recognized BSA point-to-point systematic column; use statistical errors or extend column candidates.')
        panels_by_energy = {energy: {} for energy in (10.2, 10.6) if args.e_beam is None or args.e_beam == energy}
        for row in reader:
            if not base.row_is_valid(row, cols['valid']):
                continue
            key = base.make_bin_key(row, cols)
            phi = base.row_phi(row, cols)
            if not finite(phi):
                continue
            for energy, panels in panels_by_energy.items():
                item = panels.setdefault(key, {'xs': [], 'bsa': []})
                selected_xs = xs_col if energy == 10.6 else xs_10p2
                selected_bsa = bsa_col if energy == 10.6 else bsa_10p2
                for observable, column, sys_col, enabled in [
                    ('xs', selected_xs, xs_ptp_col if energy == 10.6 else None, args.xs_total_ptp),
                    ('bsa', selected_bsa, bsa_ptp_col if energy == 10.6 else None, args.bsa_total_ptp),
                ]:
                    if column is None:
                        continue
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
    panels_by_energy = {e: {k: v for k, v in p.items() if v['xs'] or v['bsa']} for e, p in panels_by_energy.items()}
    if not any(panels_by_energy.values()):
        raise RuntimeError('No valid Pass-2 results found in CSV.')

    for energy, panels in panels_by_energy.items():
        if not panels:
            continue
        args._active_energy = energy
        cfg = base.ModelConfig(e_beam=energy, dvcsgen_dir=args.dvcsgen_dir,
                               km15_cli=args.km15_cli, py_km15=args.python_km15,
                               phi_dense=args.phi_dense, use_cache=True,
                               allow_missing_models=False, verbose_worker_models=False)
        cache_path = args.output_dir / f'model_curves_cache_{energy:g}GeV.json'
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
                fig.subplots_adjust(left=0.075, right=0.985, bottom=0.10, top=0.87, hspace=0.0, wspace=0.0)
                fig.suptitle(rf'CLAS12 Pass-2 $ep\to ep\gamma$ — {"cross section" if kind == "xs" else "beam-spin asymmetry"}' + '\n' + rf'$|t|\in[{trange[0]:.3f},\,{trange[1]:.3f}]\ \mathrm{{GeV}}^2$, $E_{{beam}}={energy:g}$ GeV', fontsize=15, y=0.985)
                occupied = {(len(qbins) - 1 - iq, ix) for iq, qb in enumerate(qbins) for ix, xb in enumerate(xbins) if any((k.xb_min, k.xb_max) == xb and (k.q2_min, k.q2_max) == qb and panels[k][kind] for k in keys)}
                bottom_by_col = {ix: max(row for row, col in occupied if col == ix) for ix in range(len(xbins)) if any(col == ix for row, col in occupied)}
                left_by_row = {row: min(col for rr, col in occupied if rr == row) for row in range(len(qbins)) if any(rr == row for rr, col in occupied)}
                for iq, qb in enumerate(qbins):
                    for ix, xb in enumerate(xbins):
                        ax = axs[len(qbins) - 1 - iq, ix]
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
                        ax.errorbar(phi, yy, yerr=ee, fmt='o', ms=3.5, lw=0.9, capsize=1.5, color='black', label='CLAS12 Pass-2', zorder=4)
                        cache_key = model_cache_key(k, args, kind)
                        curves = cache[cache_key]
                        for model, color, name in [('km15', '#1f77b4', 'KM15'), ('bh', '#ff7f0e', 'BH')]:
                            vv = np.asarray(curves[model], dtype=float)
                            xx = np.asarray(curves['phi'], dtype=float)
                            valid = np.isfinite(vv) & np.isfinite(xx)
                            if np.any(valid):
                                ax.plot(xx[valid], vv[valid], '--' if model == 'bh' else '-', lw=1.65, color=color, label=name, zorder=2)
                        ax.set_xlim(0, 360)
                        ax.set_xticks([0, 90, 180, 270, 360])
                        ax.grid(alpha=0.18)
                        ax.text(0.035, 0.96, rf'$x_B=[{xb[0]:.3f},{xb[1]:.3f}]$' + '\n' + rf'$Q^2=[{qb[0]:.2f},{qb[1]:.2f}]\ \mathrm{{GeV}}^2$', transform=ax.transAxes, va='top', ha='left', fontsize=8, bbox={'facecolor': 'white', 'alpha': 0.78, 'edgecolor': 'none', 'pad': 2})
                        if kind == 'xs':
                            # Logarithmic cross-section axis; nonpositive points cannot be displayed.
                            ax.set_yscale('log')
                        else:
                            ax.set_ylim(-0.65, 0.65)
                        row = len(qbins) - 1 - iq
                        show_x = row == bottom_by_col.get(ix)
                        show_y = ix == left_by_row.get(row)
                        ax.tick_params(axis='x', labelbottom=show_x, bottom=show_x)
                        ax.tick_params(axis='y', labelleft=show_y, left=show_y)
                        if show_x:
                            ax.set_xlabel(r'$\phi$ (deg)')
                        if show_y:
                            ax.set_ylabel(r'$d^4\sigma$ (analysis units)' if kind == 'xs' else r'$A_{LU}$')
                if kind == 'xs':
                    # Each Q2 row has one common logarithmic range, including
                    # measured uncertainties and both model predictions.
                    # Do not share y across rows: their scales can differ greatly.
                    for row in range(len(qbins)):
                        positives = []
                        for col in range(len(xbins)):
                            if (row, col) not in occupied:
                                continue
                            qb = qbins[len(qbins) - 1 - row]
                            xb = xbins[col]
                            matches = [k for k in keys if (k.xb_min, k.xb_max) == xb and (k.q2_min, k.q2_max) == qb and panels[k][kind]]
                            if not matches:
                                continue
                            k = matches[0]
                            for _, val, err in panels[k][kind]:
                                for candidate in (val, val - err, val + err):
                                    if finite(candidate) and candidate > 0:
                                        positives.append(candidate)
                            curves = cache[model_cache_key(k, args, kind)]
                            for model in ('km15', 'bh'):
                                positives.extend(float(v) for v in curves[model] if finite(float(v)) and float(v) > 0)
                        #endfor
                        if positives:
                            lo, hi = min(positives), max(positives)
                            padding = max(0.10, 0.06 * math.log10(hi / lo)) if hi > lo else 0.15
                            limits = (10 ** (math.log10(lo) - padding), 10 ** (math.log10(hi) + padding))
                            for col in range(len(xbins)):
                                if (row, col) in occupied:
                                    axs[row, col].set_ylim(*limits)
                                    # Only the leftmost occupied panel in this row
                                    # displays y tick labels (including minor ticks).
                                    if col != left_by_row[row]:
                                        axs[row, col].tick_params(axis='y', which='both', left=False, right=False, labelleft=False, labelright=False)
                                        axs[row, col].yaxis.get_offset_text().set_visible(False)
                            #endfor
                        #endif
                    #endfor
                # Shared axes may hide tick labels on non-bottom rows; restore only
                # the lowest populated panel in each column (staircase layouts).
                for ix, row in bottom_by_col.items():
                    axs[row, ix].tick_params(axis='x', labelbottom=True, bottom=True)
                # Keep endpoint tick labels inside their own panel. At zero
                # inter-panel spacing, centered 0/360 labels cross the shared spine.
                for ix, row in bottom_by_col.items():
                    ax = axs[row, ix]
                    # Only the leftmost populated panel in this row labels 0°.
                    # Adjacent panels otherwise visually concatenate 360 and 0.
                    labels = ['0' if ix == left_by_row[row] else '',
                              '90', '180', '270', '360']
                    ax.set_xticks([0, 90, 180, 270, 360])
                    ax.set_xticklabels(labels)
                    tick_labels = ax.get_xticklabels()
                    tick_labels[0].set_horizontalalignment('left')
                    tick_labels[-1].set_horizontalalignment('right')
                from matplotlib.lines import Line2D
                handles = [Line2D([], [], color='black', marker='o', ls='none', markersize=5, label='CLAS12 Pass-2'), Line2D([], [], color='#1f77b4', lw=2, ls='-', label='KM15'), Line2D([], [], color='#ff7f0e', lw=2, ls='--', label='BH')]
                fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 0.914), ncol=3, frameon=False)
                output = args.output_dir / f'pass2_{energy:g}GeV_{kind}_tbin_{it:02d}_{trange[0]:.3f}_{trange[1]:.3f}.png'
                fig.savefig(output, dpi=args.dpi)
                plt.close(fig)
                print(f'[write] {output}', flush=True)
    print('[done] All requested PNG canvases generated.', flush=True)


if __name__ == '__main__':
    main()
