#!/usr/bin/env python3
"""BSA point-to-point systematic components vs xB, Q2 and |t|.

Uses the row indices in output/systematics/bsa_systematics_summary.csv to join
with output/csvs/dvcs_pass2_analysis.csv. The latter stores group-specific mean
kinematics as columns such as "xBavg, Sp18" (not plain "xBavg"). Each kinematic interval is summarized
by the median of the individual phi-bin absolute uncertainties, with a shaded
16th–84th percentile envelope. Sp18, Fa18 and Sp19 remain separate.
"""
import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

GROUPS = [('Sp18', 'Sp18', '#E69F00'), ('Fa18', 'Fa18', '#0072B2'),
          ('Sp19 Inb', 'Sp19 (10.2 GeV)', '#D55E00')]
COMPONENTS = [
    ('pi0_subtraction_sys', r'$\pi^0$ background subtraction'),
    ('exclusivity_sys', 'Exclusivity selections'),
    ('fiducial_sys', 'Fiducial selections'),
    ('point_to_point_sys', 'Total point-to-point systematic'),
]
AXES = [
    ('xB', 'xBmin', 'xBmax', 'xBavg', r'$x_B$'),
    ('Q2', 'Q2min', 'Q2max', 'Q2avg', r'$Q^2\ (\mathrm{GeV}^2)$'),
    ('t', 't_abs_min', 't_abs_max', 't_abs_avg', r'$|t|\ (\mathrm{GeV}^2)$'),
]


def read_csv(path):
    with path.open(newline='') as f:
        reader = csv.DictReader(f)
        return reader.fieldnames or [], list(reader)


def require_columns(columns, required, path):
    missing = set(required) - set(columns)
    if missing:
        raise ValueError(f'{path}: missing columns {sorted(missing)}')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--systematics', default='output/systematics/bsa_systematics_summary.csv')
    p.add_argument('--kinematics', default='output/csvs/dvcs_pass2_analysis.csv')
    p.add_argument('--output-dir', default='output/systematics/analysis_note/bsa')
    p.add_argument('--statistic', choices=('median', 'mean'), default='median',
                   help='Central tendency for each kinematic bin (default: median)')
    args = p.parse_args()
    sys_path, kin_path = Path(args.systematics), Path(args.kinematics)
    for path in (sys_path, kin_path):
        if not path.exists():
            p.error(f'Missing {path}')
    sys_cols, sys_rows = read_csv(sys_path)
    kin_cols, kin_rows = read_csv(kin_path)
    require_columns(sys_cols, ['group', 'row'] + [x[0] for x in COMPONENTS], sys_path)
    require_columns(kin_cols, [v for _, low, high, avg, _ in AXES for v in (low, high)], kin_path)
    # Mean kinematics are stored separately for each period, e.g. 'Q2avg, Fa18'.
    # Fall back to the geometric bin midpoint if a group-specific mean is absent.
    for group, _, _ in GROUPS:
        for _, _, _, avg, _ in AXES:
            col = f'{avg}, {group}'
            if col not in kin_cols:
                print(f'[warning] {kin_path}: no {col!r}; using bin midpoint')
    groups = {g for g, _, _ in GROUPS}
    joined = []
    for row in sys_rows:
        if row['group'].strip() not in groups:
            continue
        try:
            index = int(row['row'])
            if index < 0 or index >= len(kin_rows):
                raise IndexError(index)
            kin = kin_rows[index]
            vals = {field: float(row[field]) for field, _ in COMPONENTS}
            if not all(np.isfinite(v) and v >= 0 for v in vals.values()):
                continue
            coords = {}
            for name, low, high, avg, _ in AXES:
                l, h = (float(kin[key]) for key in (low, high))
                avg_key = f"{avg}, {row['group'].strip()}"
                try:
                    a = float(kin.get(avg_key, ''))
                except (ValueError, TypeError):
                    a = float('nan')
                if not np.isfinite(a):
                    a = 0.5 * (l + h)
                if not (np.isfinite(l) and np.isfinite(h) and np.isfinite(a)):
                    raise ValueError('nonfinite kinematic value')
                coords[name] = (l, h, a)
            joined.append((row['group'].strip(), coords, vals))
        except (ValueError, TypeError, KeyError, IndexError) as exc:
            raise ValueError(f'Cannot join systematics row {row.get("row")} ({row.get("group")}) '
                             f'to {kin_path}: {exc}') from exc
    if not joined:
        p.error('No valid rows for Sp18, Fa18 or Sp19')

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    summary = []
    for axis, _, _, _, xlabel in AXES:
        fig, axes = plt.subplots(2, 2, figsize=(12.8, 8.8), constrained_layout=True)
        for ax, (field, title) in zip(axes.flat, COMPONENTS):
            for group, label, color in GROUPS:
                bins = defaultdict(list)
                centers = defaultdict(list)
                for g, coords, values in joined:
                    if g != group:
                        continue
                    lo, hi, center = coords[axis]
                    key = (round(lo, 7), round(hi, 7))
                    bins[key].append(values[field])
                    centers[key].append(center)
                records = []
                for key, values in bins.items():
                    arr = np.asarray(values)
                    x = float(np.median(centers[key]))
                    center_y = float(np.median(arr) if args.statistic == 'median' else np.mean(arr))
                    p16, p84 = np.percentile(arr, [16, 84])
                    records.append((x, center_y, float(p16), float(p84), len(arr), *key))
                records.sort(key=lambda r: r[0])
                if not records:
                    continue
                x = np.asarray([r[0] for r in records])
                y = np.asarray([r[1] for r in records])
                lower = np.asarray([r[2] for r in records])
                upper = np.asarray([r[3] for r in records])
                ax.plot(x, y, '-o', color=color, label=label, linewidth=1.8, markersize=4.8)
                ax.fill_between(x, lower, upper, color=color, alpha=0.12, linewidth=0)
                for r in records:
                    summary.append(dict(variable=axis, component=field, group=group,
                                        bin_min=r[5], bin_max=r[6], bin_center=r[0],
                                        n=r[4], statistic=args.statistic, central=r[1],
                                        p16=r[2], p84=r[3]))
            ax.set_title(title, fontsize=13)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(r'Absolute BSA uncertainty $\delta A_{LU}$')
            ax.set_ylim(bottom=0)
            ax.grid(alpha=0.20)
            ax.legend(frameon=False, fontsize=9)
        path = out / f'bsa_point_to_point_systematics_vs_{axis}.png'
        fig.savefig(path, dpi=220)
        plt.close(fig)
        print('[write]', path)
    summary_path = out / 'bsa_systematics_kinematic_summary.csv'
    with summary_path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)
    print('[write]', summary_path)
    print(f'[info] Joined {len(joined)} group/row entries to the nominal kinematic CSV')
    print(f'[info] Shaded regions: 16th–84th percentiles (not uncertainties on the {args.statistic})')


if __name__ == '__main__':
    main()
