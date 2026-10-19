#!/usr/bin/env python3
"""DVCS Pass-2/Pass-1 analysis-note plots from existing comparison CSVs.

Direct Pass-2/Pass-1 ratios use the authoritative Pass-1 matching output.
The matched CSV does not include Pass-2 statistical errors: interval bars
represent the central 68% of the observed point distribution, NOT errors on
an estimated mean or median. No KM15 transport is applied to these plots.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 11, 'axes.labelsize': 12,
                     'axes.titlesize': 12, 'legend.fontsize': 9,
                     'savefig.dpi': 250,
                     'axes.linewidth': 1.0})
PASS1_NORM = 0.31  # Pass-1 pi0-model normalization scale, fully correlated
COLORS = {'Inclusive': '#222222', 'CD-FD': '#0072B2',
          'CD-FT': '#D55E00', 'FD-FD': '#009E73'}


def find_csv(folder, filename):
    exact = folder / filename
    if exact.is_file():
        return exact
    matches = sorted(folder.rglob(filename)) if folder.exists() else []
    if not matches:
        raise FileNotFoundError(f'Cannot find {filename} anywhere under {folder}')
    if len(matches) > 1:
        raise RuntimeError(f'Ambiguous {filename}: {matches}')
    return matches[0]


def numeric(frame, name):
    return pd.to_numeric(frame[name], errors='coerce')


def summarize(frame, xname, edges, ratio_name, label, rng=None):
    if rng is None:
        rng = np.random.default_rng(20261019)
    x = numeric(frame, xname).to_numpy(dtype=float)
    r = numeric(frame, ratio_name).to_numpy(dtype=float)
    rows = []
    for low, high in zip(edges[:-1], edges[1:]):
        mask = np.isfinite(x) & np.isfinite(r) & (r > 0) & (x >= low) & (x < high)
        if not mask.any():
            continue
        values = r[mask]
        # Resample matched ratios to estimate the sampling uncertainty of
        # the median. This is NOT the propagated cross-section stat error.
        if len(values) >= 2:
            draws = rng.choice(values, size=(2000, len(values)), replace=True)
            medians = np.median(draws, axis=1)
            ci16, ci84 = np.percentile(medians, [16, 84])
        else:
            ci16 = ci84 = float(values[0])
        rows.append({'selection': label, 'variable': xname,
                     'lower': low, 'upper': high, 'N': len(values),
                     'x_mean': np.mean(x[mask]), 'median': np.median(values),
                     'median_bootstrap_p16': ci16, 'median_bootstrap_p84': ci84,
                     'p16': np.percentile(values, 16),
                     'p84': np.percentile(values, 84),
                     'arithmetic_mean': np.mean(values)})
    return pd.DataFrame(rows)


def plot_series(ax, stats, label, color, marker='o', offset=0,
                show_counts=False):
    if stats.empty:
        return
    x = stats.x_mean.to_numpy(dtype=float) + offset
    y = stats['median'].to_numpy(dtype=float)
    # Gray: 16th–84th percentiles of individual matched ratios.
    spread_lo = np.maximum(y - stats.p16.to_numpy(dtype=float), 0)
    spread_hi = np.maximum(stats.p84.to_numpy(dtype=float) - y, 0)
    ax.errorbar(x, y, yerr=np.vstack([spread_lo, spread_hi]), fmt='none',
                ecolor='0.68', elinewidth=1.15, capsize=3, capthick=1,
                zorder=2)
    # Colored: bootstrap 68% confidence interval on the binned median.
    ci_lo = np.maximum(y - stats.median_bootstrap_p16.to_numpy(dtype=float), 0)
    ci_hi = np.maximum(stats.median_bootstrap_p84.to_numpy(dtype=float) - y, 0)
    ax.errorbar(x, y, yerr=np.vstack([ci_lo, ci_hi]), fmt=marker + '-',
                lw=1.5, ms=5.5, capsize=4, capthick=2.2, elinewidth=2.6,
                color=color, label=label, zorder=4)
    if show_counts:
        # Annotate above the upper distribution bar. Keep labels clear of data.
        for xi, yi, upper, count in zip(x, y, spread_hi, stats.N):
            ax.annotate(f'N={int(count)}', (xi, yi + upper),
                        xytext=(0, 6), textcoords='offset points',
                        ha='center', va='bottom', fontsize=8, color='0.32')


def save(fig, path):
    target = path.with_suffix('.png')
    fig.savefig(target, bbox_inches='tight', dpi=250)
    plt.close(fig)
    print('[write]', target)


def pass1_normalization_band(ax, label=False):
    """Pass-1 +/-31% scale, mapped exactly onto a Pass-2/Pass-1 ratio.

    Since Pass-1 is in the denominator, the resulting ratio interval is
    asymmetric: [1/(1+0.31), 1/(1-0.31)]. It is a single correlated
    normalization nuisance, NOT independent point-to-point uncertainty.
    """
    low, high = 1 / (1 + PASS1_NORM), 1 / (1 - PASS1_NORM)
    ax.axhspan(low, high, color='#B8C6D8', alpha=.24, lw=0,
               zorder=0, label='Pass-1 norm. (31%, correlated)' if label else None)
    ax.axhline(1, ls='--', lw=1, color='0.45', zorder=1)


def attach_inclusive_angles(inclusive, csv_path):
    """Join original Pass-2 combined mean angles to matched comparison rows.

    The matching study exports the 4D combined Pass-2 means but not angles.
    Use those four original values as a strict composite key; do not assign
    detector angles based on the published Pass-1 bin center.
    """
    if not csv_path.is_file():
        raise FileNotFoundError(f'Angle plots need the original inclusive Pass-2 CSV: {csv_path}')
    names = {'xB_pass2': 'xBavg, 10.6 GeV',
             'Q2_pass2': 'Q2avg, 10.6 GeV',
             't_pass2': 't_abs_avg, 10.6 GeV',
             'phi_pass2': 'phiavg, 10.6 GeV',
             'e_theta': 'e_theta, 10.6 GeV',
             'p_theta': 'p_theta, 10.6 GeV',
             'g_theta': 'g_theta, 10.6 GeV'}
    raw = pd.read_csv(csv_path, usecols=list(names.values()), low_memory=False)
    raw = raw.rename(columns={v: k for k, v in names.items()})
    left = inclusive.copy()
    keys = ['xB_pass2', 'Q2_pass2', 't_pass2', 'phi_pass2']
    # Both tables originate from the same Pass-2 CSV; 7-decimal rounding
    # accommodates CSV float formatting without permitting coarse bin matching.
    for key in keys:
        raw[key] = pd.to_numeric(raw[key], errors='coerce')
        left[key] = pd.to_numeric(left[key], errors='coerce')
        raw['_join_' + key] = raw[key].round(7)
        left['_join_' + key] = left[key].round(7)
    join_keys = ['_join_' + key for key in keys]
    if raw.duplicated(join_keys).any():
        raise ValueError('Original inclusive Pass-2 CSV has duplicate combined 4D mean kinematics; cannot assign angles uniquely')
    angles = raw[join_keys + ['e_theta', 'p_theta', 'g_theta']]
    merged = left.merge(angles, how='left', on=join_keys, validate='many_to_one')
    matched = merged['e_theta'].notna() & merged['p_theta'].notna() & merged['g_theta'].notna()
    print(f'[angles] Matched {matched.sum()}/{len(merged)} inclusive points to original Pass-2 angle means')
    if matched.sum() < .95 * len(merged):
        raise ValueError('Less than 95% of matched Pass-2 points recovered angles. Check inclusive CSV provenance and combined kinematics.')
    return merged.drop(columns=join_keys)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--topology-dir', type=Path, required=True)
    parser.add_argument('--common-dir', type=Path)
    parser.add_argument('--outdir', type=Path,
                        default=Path('output/world_cross_section_note_figures'))
    parser.add_argument('--bins', type=int, default=7)
    parser.add_argument('--inclusive-csv', type=Path, default=Path('output/csvs/dvcs_pass2_INCLUSIVE.csv'),
                        help='Original inclusive Pass-2 CSV supplying combined mean detector angles')
    args = parser.parse_args()
    if args.bins < 2:
        parser.error('--bins must be at least 2')
    args.outdir.mkdir(parents=True, exist_ok=True)

    match_path = find_csv(args.topology_dir,
                          'pass2_vs_authoritative_pass1_matched_points.csv')
    print('[input]', match_path)
    data = pd.read_csv(match_path)
    required = {'topology', 'xB_pass2', 'Q2_pass2', 't_pass2',
                'pass2_over_pass1'}
    missing = required - set(data.columns)
    if missing:
        raise ValueError(f'Missing required matched-table columns: {sorted(missing)}')
    data['abs_t_pass2'] = numeric(data, 't_pass2').abs()
    inclusive = data[data.topology.astype(str).str.casefold() == 'inclusive'].copy()
    if inclusive.empty:
        raise ValueError('Matched comparison table contains no Inclusive rows')

    # Main comparison: identical nominal binning; no model-based transport.
    specs = [('Q2_pass2', r'$Q^2$ (GeV$^2$)'),
             ('xB_pass2', r'$x_B$'),
             ('abs_t_pass2', r'$|t|$ (GeV$^2$)')]
    summaries = []
    fig, axes = plt.subplots(1, 3, figsize=(14.3, 4.2), layout='constrained')
    for ax, (key, xlabel) in zip(axes, specs):
        values = numeric(inclusive, key).to_numpy(dtype=float)
        values = values[np.isfinite(values)]
        edges = np.unique(np.quantile(values, np.linspace(0, 1, args.bins + 1)))
        if len(edges) < 3:
            raise ValueError(f'Insufficient distinct values for {key}')
        edges[-1] = np.nextafter(edges[-1], np.inf)
        stats = summarize(inclusive, key, edges, 'pass2_over_pass1', 'Inclusive')
        summaries.append(stats)
        plot_series(ax, stats, 'Inclusive', COLORS['Inclusive'], show_counts=True)
        pass1_normalization_band(ax)
        ax.set_xlabel(xlabel)
        ax.grid(alpha=0.18)
        ax.set_ylim(bottom=0)
    axes[0].set_ylabel('Pass-2 / Pass-1')
    fig.text(0.5, -0.035, 'Thick black: 68% bootstrap median CI   |   Thin gray: point spread (16–84%)   |   Shading: Pass-1 norm. (31%, correlated)   |   N: points',
             ha='center', va='top', fontsize=9, color='0.35')
    fig.suptitle('CLAS12 DVCS cross sections: direct Pass-2 / Pass-1', fontsize=13)
    save(fig, args.outdir / 'pass2_pass1_vs_kinematics')
    pd.concat(summaries, ignore_index=True).to_csv(
        args.outdir / 'pass2_pass1_binned_summary.csv', index=False)

    # Angular dependence from the original inclusive Pass-2 mean angles.
    with_angles = attach_inclusive_angles(inclusive, args.inclusive_csv)
    angle_specs = [('e_theta', r'$\theta_e$ (deg)'),
                   ('p_theta', r'$\theta_p$ (deg)'),
                   ('g_theta', r'$\theta_\gamma$ (deg)')]
    fig, axes = plt.subplots(1, 3, figsize=(14.3, 4.2), layout='constrained')
    angle_summaries = []
    for ax, (key, xlabel) in zip(axes, angle_specs):
        values = numeric(with_angles, key).to_numpy(dtype=float)
        values = values[np.isfinite(values)]
        edges = np.unique(np.quantile(values, np.linspace(0, 1, args.bins + 1)))
        if len(edges) < 3:
            raise ValueError(f'Insufficient distinct values for {key}')
        edges[-1] = np.nextafter(edges[-1], np.inf)
        stats = summarize(with_angles, key, edges, 'pass2_over_pass1', 'Inclusive')
        angle_summaries.append(stats)
        plot_series(ax, stats, 'Inclusive', COLORS['Inclusive'], show_counts=True)
        pass1_normalization_band(ax)
        ax.set_xlabel(xlabel)
        ax.grid(alpha=.18)
        ax.set_ylim(bottom=0)
    axes[0].set_ylabel('Pass-2 / Pass-1')
    fig.suptitle('CLAS12 DVCS cross sections: dependence on laboratory angles', fontsize=13)
    fig.text(0.5, -0.035,
             'Thick black: 68% bootstrap median CI   |   Thin gray: point spread (16–84%)   |   Shading: Pass-1 norm. (31%, correlated)',
             ha='center', va='top', fontsize=9, color='0.35')
    save(fig, args.outdir / 'pass2_pass1_vs_angles')
    pd.concat(angle_summaries, ignore_index=True).to_csv(
        args.outdir / 'pass2_pass1_angle_binned_summary.csv', index=False)

    # Topology comparison. These samples have different coverage; this is
    # descriptive, not an acceptance-matched topology-efficiency measurement.
    fig, ax = plt.subplots(figsize=(8.2, 5.0), layout='constrained')
    edges = np.array([0, .2, .3, .4, 1.0, np.inf])
    topology_stats = []
    for i, tag in enumerate(['Inclusive', 'CD-FD', 'CD-FT', 'FD-FD']):
        part = data[data.topology.astype(str).str.casefold() == tag.casefold()]
        if part.empty:
            continue
        stats = summarize(part, 'abs_t_pass2', edges, 'pass2_over_pass1', tag)
        topology_stats.append(stats)
        plot_series(ax, stats, f'{tag} (N={len(part)})', COLORS[tag],
                    offset=(i - 1.5) * .004)
    pass1_normalization_band(ax, label=True)
    ax.set(xlabel=r'$|t|$ (GeV$^2$)', ylabel='Pass-2 / Pass-1')
    ax.grid(alpha=.18)
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, ncol=2)
    fig.suptitle('CLAS12 DVCS cross sections by detector topology', fontsize=13)
    fig.text(0.5, -0.035, 'Thick colored: 68% bootstrap median CI   |   Thin gray: point spread (16–84%)   |   Shading: Pass-1 norm. (31%, correlated)',
             ha='center', va='top', fontsize=9, color='0.35')
    save(fig, args.outdir / 'pass2_pass1_by_topology')
    pd.concat(topology_stats, ignore_index=True).to_csv(
        args.outdir / 'topology_binned_summary.csv', index=False)

    if args.common_dir is not None:
        common_path = find_csv(args.common_dir,
                               'topology_common_kinematics_bin_by_bin.csv')
        print('[input]', common_path)
        common = pd.read_csv(common_path)
        needed = {'dxB_fd_minus_ft', 'dQ2_fd_minus_ft', 'R_raw', 't_common'}
        missing = needed - set(common.columns)
        if missing:
            raise ValueError(f'Missing common-kinematics columns: {sorted(missing)}')
        dx = numeric(common, 'dxB_fd_minus_ft').abs()
        dq = numeric(common, 'dQ2_fd_minus_ft').abs()
        raw = numeric(common, 'R_raw')
        valid = np.isfinite(raw) & (raw > 0)
        tight = (dx < .01) & (dq < .10)
        fig, ax = plt.subplots(figsize=(8.0, 5.0), layout='constrained')
        common_edges = np.array([0, .15, .25, .4, .6, np.inf])
        for label, mask, color, filename in [
            ('All common bins', valid, COLORS['CD-FD'], 'common_all_summary.csv'),
            (r'Tight: $|\Delta x_B|<0.01$, $|\Delta Q^2|<0.10$ GeV$^2$',
             valid & tight, COLORS['CD-FT'], 'common_tight_summary.csv')]:
            selected = common.loc[mask]
            stats = summarize(selected, 't_common', common_edges, 'R_raw', label)
            stats.to_csv(args.outdir / filename, index=False)
            print(f'[matching] {label}: {len(selected)} points')
            plot_series(ax, stats, f'{label} (N={len(selected)})', color)
        ax.axhline(1, ls='--', lw=1, color='0.45')
        ax.set(xlabel=r'$|t|$ (GeV$^2$)', ylabel='CD-FD / CD-FT')
        ax.set_ylim(bottom=0)
        ax.grid(alpha=.18)
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(frameon=False, fontsize=9)
        fig.suptitle('CD-FD / CD-FT: sensitivity to mean-kinematics matching',
                     fontsize=13)
        fig.text(0.5, -0.035, 'Thick colored: 68% bootstrap interval on median    |    Thin gray: 16th–84th percentile spread',
                 ha='center', va='top', fontsize=9, color='0.35')
        save(fig, args.outdir / 'cd_fd_cd_ft_tight_kinematics')
        scan = []
        for max_dx, max_dq in [(.04, .4), (.03, .3), (.02, .2),
                               (.01, .1), (.005, .05)]:
            sample = common.loc[valid & (dx < max_dx) & (dq < max_dq)]
            scan.append({'max_abs_dxB': max_dx, 'max_abs_dQ2': max_dq,
                         'N': len(sample),
                         'median_raw_ratio': sample.R_raw.median(),
                         'median_km15_residual': sample.R_residual_km15.median()
                         if 'R_residual_km15' in sample else np.nan})
        pd.DataFrame(scan).to_csv(
            args.outdir / 'kinematics_matching_threshold_scan.csv', index=False)
    print('[done]', args.outdir)
    print('[note] Thick bars: bootstrap 68% confidence intervals on binned medians; '
          'thin gray bars: 16th–84th percentiles of individual ratios. '
          'Neither is a propagated cross-section statistical error.')


if __name__ == '__main__':
    main()
