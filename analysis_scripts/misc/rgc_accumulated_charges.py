#!/usr/bin/env python3

from collections import OrderedDict, defaultdict

def main():
    # Path to your CSV file
    filepath = '/u/home/thayward/clas12_analysis_software/analysis_scripts/asymmetry_extraction/imports/clas12_run_info.csv'

    # Sum of accumulated charge per target‐type group
    group_totals = OrderedDict()
    current_group = None

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            #endif

            # Header lines define a new group
            if line.startswith('#'):
                header = line.lstrip('#').strip()
                # Skip any ND3 sections
                if 'ND3' in header:
                    current_group = None
                    continue
                #endif
                current_group = header
                group_totals[current_group] = 0.0
            #endif

            else:
                # Only accumulate if we're inside a valid group
                if current_group is None:
                    continue
                #endif

                parts = [p.strip() for p in line.split(',')]
                if len(parts) < 2:
                    continue
                #endif

                try:
                    charge = float(parts[1])
                except ValueError:
                    continue
                #endif

                group_totals[current_group] += charge
        #endfor

    # Only keep the 5 RGC Su22 and 5 RGC Fa22 groups
    filtered_totals = OrderedDict(
        (grp, total)
        for grp, total in group_totals.items()
        if grp in (
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
        )
    )

    # Compute total per run‐period ("RGC Su22" vs "RGC Fa22")
    period_totals = defaultdict(float)
    for grp, total in filtered_totals.items():
        period = ' '.join(grp.split()[:2])  # e.g. "RGC Su22" or "RGC Fa22"
        period_totals[period] += total
    #endfor

    # Print summary using the exact names provided
    print(f"{'Target Type':<25}{'Total Charge':>15}{'Fraction (%)':>15}")
    for grp, total in filtered_totals.items():
        period = ' '.join(grp.split()[:2])
        frac = (total / period_totals[period] * 100) if period_totals[period] > 0 else 0
        print(f"{grp:<25}{total:15.6f}{frac:15.2f}%")
    #endfor

if __name__ == '__main__':
    main()
#endif