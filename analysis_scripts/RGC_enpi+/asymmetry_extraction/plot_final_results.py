#!/usr/bin/env python3
"""Regenerate the four final-result figures from the finished extraction CSV."""
import argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
OBSERVABLES = ('lu1', 'ul1', 'ul2', 'll0', 'll1')
LABELS = (r'$F_{LU}^{\sin\phi}/F_{UU}$',
          r'$F_{UL,\mathrm{lab}}^{\sin\phi}/F_{UU}$',
          r'$F_{UL,\mathrm{lab}}^{\sin2\phi}/F_{UU}$',
          r'$F_{LL,\mathrm{lab}}/F_{UU}$',
          r'$F_{LL,\mathrm{lab}}^{\cos\phi}/F_{UU}$')
TP_EDGES = np.array((.05, .25, .45, .65, .85, 1.05, 1.25))
XB_EDGES = (.10, .25, .35, .45, .60)


def plot_results(frame, output_dir):
    required = {'bin_number', 'x_index', 't_index', 'mean_minus_tprime_gev2'}
    required.update(o+s for o in OBSERVABLES for s in ('', '_stat', '_point_to_point_systematic'))
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f'Missing columns: {sorted(missing)}')
    frame = frame.sort_values('bin_number')
    if frame.bin_number.tolist() != list(range(1, 25)):
        raise ValueError('Expected unique analysis bins 1 through 24')
    if not np.isfinite(frame[sorted(required)].to_numpy(dtype=float)).all():
        raise ValueError('Nonfinite result values')
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for xi in range(4):
        g = frame[frame.x_index == xi].sort_values('t_index')
        if g.t_index.tolist() != list(range(6)):
            raise ValueError(f'Incomplete t-prime binning in xB group {xi}')
        fig, axes = plt.subplots(2, 3, figsize=(12, 7))
        axes = axes.ravel()
        for panel, (ax, obs, label) in enumerate(zip(axes, OBSERVABLES, LABELS)):
            syst = g[obs+'_point_to_point_systematic'].to_numpy()
            stat = g[obs+'_stat'].to_numpy()
            if (syst < 0).any() or (stat <= 0).any():
                raise ValueError(f'Invalid uncertainties for {obs}')
            # Systematic magnitude extends upward from zero across the full bin.
            ax.bar(TP_EDGES[:-1], syst, width=np.diff(TP_EDGES), align='edge',
                   bottom=0, color='0.65', alpha=.35, linewidth=0, zorder=1)
            ax.axhline(0, color='0.45', linewidth=.8, zorder=2)
            ax.errorbar(g.mean_minus_tprime_gev2, g[obs], yerr=stat,
                        fmt='o', color='black', ecolor='black', markersize=4,
                        capsize=3, linewidth=1.1, zorder=3)
            ax.set_xlim(0, 1.30)
            ax.set_ylim((-.5, .5) if panel < 3 else (-1, 1))
            ax.set_yticks((-.5, -.25, 0, .25, .5) if panel < 3 else (-1, -.5, 0, .5, 1))
            ax.set_xlabel(r'$-t^\prime$ (GeV$^2$)')
            ax.set_ylabel(label)
            ax.grid(alpha=.18, zorder=0)
        axes[5].axis('off')
        fig.suptitle(rf'${XB_EDGES[xi]:.2f}<x_B<{XB_EDGES[xi+1]:.2f}$')
        fig.tight_layout(rect=(0, 0, 1, .96))
        path = output_dir/f'final_results_xB_bin_{xi+1}.png'
        fig.savefig(path, dpi=200)
        paths.append(path)
        plt.close(fig)
    return paths


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-csv', type=Path, default=HERE/'output/asymmetry_extraction/nominal/tables/structure_function_ratios.csv')
    parser.add_argument('--output-dir', type=Path, default=HERE/'output/final_results')
    args = parser.parse_args()
    for path in plot_results(pd.read_csv(args.input_csv), args.output_dir):
        print(path)


if __name__ == '__main__':
    main()
