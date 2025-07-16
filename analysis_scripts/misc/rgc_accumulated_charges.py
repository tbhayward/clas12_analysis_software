#!/usr/bin/env python3

import os
import csv
from collections import defaultdict, OrderedDict

def main():
    # Path to your CSV file
    filepath = '/u/home/thayward/clas12_analysis_software/analysis_scripts/asymmetry_extraction/imports/clas12_run_info.csv'

    # Keep track of sum of charges per target-type group
    group_totals = OrderedDict()
    current_group = None

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # Detect header lines (group names)
            if line.startswith('#'):
                header = line.lstrip('#').strip()
                # Skip the ND3 section entirely
                if 'ND3' in header:
                    current_group = None
                    continue
                # Remove trailing 'runs' if present
                if header.lower().endswith('runs'):
                    header = header[:-len('runs')].strip()
                current_group = header
                group_totals[current_group] = 0.0
            else:
                # Data line: first field is run number, second is accumulated charge
                if current_group is None:
                    continue
                parts = [p.strip() for p in line.split(',')]
                if len(parts) < 2:
                    continue
                try:
                    charge = float(parts[1])
                except ValueError:
                    continue
                group_totals[current_group] += charge
        #endfor

    # Now compute the total charge per run period (e.g. 'RGC Su22 Inb', 'RGC Fa22')
    period_totals = defaultdict(float)
    for group, total in group_totals.items():
        # period is first three tokens of the group name
        period = ' '.join(group.split()[:3])
        period_totals[period] += total
    #endfor

    # Print a summary: each target type, its total charge, and its fraction of the period
    print(f"{'Target Type':<30}{'Total Charge':>15}{'Fraction (%)':>15}")
    for group, total in group_totals.items():
        period = ' '.join(group.split()[:3])
        pct = (total / period_totals[period] * 100) if period_totals[period] > 0 else 0
        print(f"{group:<30}{total:15.6f}{pct:15.2f}%")
    #endfor

if __name__ == '__main__':
    main()
#endif