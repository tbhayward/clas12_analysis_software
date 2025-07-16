#!/usr/bin/env python3

from collections import OrderedDict, defaultdict

def main():
    # Path to your CSV file
    filepath = '/u/home/thayward/clas12_analysis_software/analysis_scripts/asymmetry_extraction/imports/clas12_run_info.csv'

    # We'll accumulate charge only for the 10 groups of interest
    group_totals = OrderedDict()
    current_group = None

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # Header lines define a new group
            if line.startswith('#'):
                raw = line.lstrip('#').strip()
                tokens = raw.split()

                # Skip ND3 sections entirely
                if 'ND3' in tokens:
                    current_group = None
                    continue

                # Drop a trailing 'runs' token if present
                if tokens[-1].lower() == 'runs':
                    tokens = tokens[:-1]

                # Build period and target
                period = ' '.join(tokens[0:2])   # "RGC Su22" or "RGC Fa22"
                target = tokens[-1]              # e.g. "NH3", "C", "CH2", "He", or "ET"

                # Only track the ten groups we're interested in
                if period in ('RGC Su22', 'RGC Fa22') and target in ('NH3','C','CH2','He','ET'):
                    current_group = f"{period} {target}"
                    group_totals[current_group] = 0.0
                else:
                    current_group = None
                continue

            # Data lines: accumulate charge if we're inside one of the tracked groups
            if current_group is not None:
                parts = [p.strip() for p in line.split(',')]
                if len(parts) < 2:
                    continue
                try:
                    charge = float(parts[1])
                except ValueError:
                    continue
                group_totals[current_group] += charge

    # Compute total per run period
    period_totals = defaultdict(float)
    for grp, total in group_totals.items():
        period = ' '.join(grp.split()[0:2])
        period_totals[period] += total

    # Print out the ten groups with exact names and their fractions
    print(f"{'Target Type':<20}{'Total Charge':>15}{'Fraction (%)':>15}")
    for grp, total in group_totals.items():
        period = ' '.join(grp.split()[0:2])
        frac = (total / period_totals[period] * 100) if period_totals[period] > 0 else 0
        print(f"{grp:<20}{total:15.6f}{frac:15.2f}%")

if __name__ == '__main__':
    main()