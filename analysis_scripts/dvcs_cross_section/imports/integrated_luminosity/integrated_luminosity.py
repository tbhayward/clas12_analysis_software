#!/usr/bin/env python3

import csv
import os

def process_file(filename):
    rows = []
    # Read and parse full rows
    with open(filename, newline='') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            try:
                # ensure at least 4 columns parseable
                key = int(row[0])
                _ = float(row[1]); _ = float(row[2]); _ = float(row[3])
                rows.append((key, row))
            except (ValueError, IndexError):
                # skip malformed lines
                continue
        #endfor

    # Sort by the first column
    rows.sort(key=lambda x: x[0])
    sorted_rows = [r for _, r in rows]

    # Write sorted data back to a new file, rounding floats to 4 decimals
    base, ext = os.path.splitext(filename)
    outname = f"{base}_sorted{ext}"
    with open(outname, 'w', newline='') as outcsv:
        writer = csv.writer(outcsv)
        for r in sorted_rows:
            formatted = [ str(int(r[0])) ]
            for item in r[1:]:
                try:
                    num = float(item)
                    formatted.append(f"{num:.4f}")
                except ValueError:
                    formatted.append(item)
            writer.writerow(formatted)
        #endfor

    # Compute sums of columns 2, 3, and 4
    sum_col2 = sum(float(r[1]) for _, r in rows)
    sum_col3 = sum(float(r[2]) for _, r in rows)
    sum_col4 = sum(float(r[3]) for _, r in rows)
    sum_col34 = sum_col3 + sum_col4

    # Print results
    print(f"File: {filename}")
    print(f"  total charge: {sum_col2} (nC)")
    print(f"  sum of pos+neg helicity charge: {sum_col34} (nC)")
    print(f"  pos beam helicity charge: {sum_col3} (nC)")
    print(f"  neg beam helicity charge: {sum_col4} (nC)")
    print(f"  sorted file written to: {outname}")
    print()
#endfor

def main():
    files = [
        'rga_fa18_inb.txt',
        'rga_fa18_out.txt',
        'rga_sp19_inb.txt',
    ]
    for f in files:
        process_file(f)
    #endfor

if __name__ == '__main__':
    main()
#endif