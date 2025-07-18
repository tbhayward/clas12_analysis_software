#!/usr/bin/env python3

from collections import OrderedDict, defaultdict

def main():
    # Path to your CSV file
    filepath = '/u/home/thayward/clas12_analysis_software/analysis_scripts/asymmetry_extraction/imports/clas12_run_info.csv'

    # Exact group names to track
    groups = [
        'RGC Su22 NH3',
        'RGC Su22 C',
        'RGC Su22 CH2',
        'RGC Su22 He',
        'RGC Su22 ET',
        'RGC Fa22 NH3',
        'RGC Fa22 C',
        'RGC Fa22 CH2',
        'RGC Fa22 He',
        'RGC Fa22 ET',
        'RGC Sp23 NH3',
        'RGC Sp23 C',
        'RGC Sp23 CH2',
        'RGC Sp23 He',
        'RGC Sp23 ET',
    ]

    # Initialize totals
    group_totals = OrderedDict((g, 0.0) for g in groups)
    current_group = None

    # Read and accumulate
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith('#'):
                header = line.lstrip('#').strip()
                # Skip ND3 sections
                if 'ND3' in header:
                    current_group = None
                    continue
                # Strip trailing "runs" if present
                if header.lower().endswith('runs'):
                    header = header[:-len('runs')].strip()
                # Track only if in our list
                if header in group_totals:
                    current_group = header
                else:
                    current_group = None
                continue

            # Data line: accumulate if inside a tracked group
            if current_group:
                parts = [p.strip() for p in line.split(',')]
                if len(parts) < 2:
                    continue
                try:
                    charge = float(parts[1])
                    group_totals[current_group] += charge
                except ValueError:
                    pass

    # Compute totals per period
    period_totals = defaultdict(float)
    for grp, tot in group_totals.items():
        period = ' '.join(grp.split()[:2])  # "RGC Su22", "RGC Fa22", or "RGC Sp23"
        period_totals[period] += tot

    # Print results
    print(f"{'Target Type':<20}{'Total Charge':>15}{'Fraction (%)':>15}")
    for grp in groups:
        tot = group_totals[grp]
        period = ' '.join(grp.split()[:2])
        frac = (tot / period_totals[period] * 100) if period_totals[period] else 0
        print(f"{grp:<20}{tot:15.6f}{frac:15.2f}%")

if __name__ == '__main__':
    main()