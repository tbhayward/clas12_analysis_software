#!/usr/bin/env python3
"""Plot final BSA point-to-point systematics for Sp18, Fa18 and Sp19.

Reads output/systematics/bsa_systematics_summary.csv created by the C++
--finalize-only systematics workflow. No event-level data or ROOT needed.
"""
import argparse
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

GROUPS = [('Sp18', '#E69F00'), ('Fa18', '#0072B2'), ('Sp19 Inb', '#D55E00')]
FIELDS = [
    ('pi0_subtraction_sys', r'$\pi^0$ background subtraction'),
    ('exclusivity_sys', 'Exclusivity selections'),
    ('fiducial_sys', 'Fiducial selections'),
    ('point_to_point_sys', 'Total point-to-point systematic'),
]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input', default='output/systematics/bsa_systematics_summary.csv')
    ap.add_argument('--output-dir', default='output/systematics/analysis_note')
    ap.add_argument('--bins', type=int, default=28)
    args = ap.parse_args()
    source = Path(args.input)
    if not source.exists():
        ap.error(f'Missing {source}. Run the C++ systematics finalization first.')
    data = {g: {f: [] for f, _ in FIELDS} for g, _ in GROUPS}
    with source.open(newline='') as fh:
        reader = csv.DictReader(fh)
        needed = {'group', *(f for f, _ in FIELDS)}
        missing = needed - set(reader.fieldnames or [])
        if missing:
            ap.error('Missing CSV columns: ' + ', '.join(sorted(missing)))
        for row in reader:
            group = row['group'].strip()
            if group not in data:
                continue
            for field, _ in FIELDS:
                try:
                    value = float(row[field])
                except (ValueError, TypeError):
                    continue
                if np.isfinite(value) and value >= 0:
                    data[group][field].append(value)
    for group, _ in GROUPS:
        if not data[group]['point_to_point_sys']:
            ap.error(f'No finite BSA systematics for {group} in {source}')

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    fig, axs = plt.subplots(2, 2, figsize=(12.5, 8.8), constrained_layout=True)
    rows = []
    for ax, (field, title) in zip(axs.flat, FIELDS):
        arrays = [np.asarray(data[g][field], dtype=float) for g, _ in GROUPS]
        maximum = max(float(a.max()) for a in arrays if len(a))
        upper = max(maximum * 1.035, 0.002)
        edges = np.linspace(0, upper, args.bins + 1)
        for (group, color), arr in zip(GROUPS, arrays):
            if not len(arr):
                continue
            label = ('Sp19 (10.2 GeV)' if group == 'Sp19 Inb' else group)
            ax.hist(arr, bins=edges, density=True, histtype='step', linewidth=2.0,
                    color=color, label=f'{label} (median {np.median(arr):.4f})')
            rows.append({
                'component': field, 'group': group, 'n': len(arr),
                'mean': float(np.mean(arr)), 'median': float(np.median(arr)),
                'p16': float(np.percentile(arr, 16)),
                'p84': float(np.percentile(arr, 84)),
                'maximum': float(np.max(arr)),
            })
        ax.set_title(title, fontsize=13)
        ax.set_xlabel(r'Absolute BSA uncertainty $\delta A_{LU}$')
        ax.set_ylabel('Probability density')
        ax.set_xlim(0, upper)
        ax.set_ylim(bottom=0)
        ax.grid(axis='y', alpha=0.2)
        ax.legend(fontsize=8.3, loc='best', frameon=False)
    fig.savefig(out / 'bsa_point_to_point_systematics.png', dpi=220)
    fig.savefig(out / 'bsa_point_to_point_systematics.pdf')
    plt.close(fig)
    with (out / 'bsa_systematics_statistics.csv').open('w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f'[write] {out / "bsa_point_to_point_systematics.png"}')
    print(f'[write] {out / "bsa_point_to_point_systematics.pdf"}')
    print(f'[write] {out / "bsa_systematics_statistics.csv"}')
    for row in rows:
        print(f"{row['group']:9s} {row['component']:35s} n={row['n']:4d} median={row['median']:.5f} p16-p84=({row['p16']:.5f}, {row['p84']:.5f})")

if __name__ == '__main__':
    main()
