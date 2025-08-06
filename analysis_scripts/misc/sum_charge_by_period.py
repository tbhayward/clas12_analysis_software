#!/usr/bin/env python3
import argparse
import sys
from collections import defaultdict

def parse_csv(path):
    data = defaultdict(lambda: defaultdict(float))
    periods = ["RGC Su22", "RGC Fa22", "RGC Sp23"]
    species_list = ["NH3", "C", "CH2", "He", "ET"]

    current_period = None
    current_species = None

    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if line.startswith("# RGC"):
                # e.g. "# RGC Su22 NH3"
                tokens = line.split()
                current_period = f"{tokens[1]} {tokens[2]}"
                current_species = tokens[3]
                if current_period not in periods or current_species not in species_list:
                    # skip any unexpected header
                    current_period = None
                    current_species = None
                continue
            if line.startswith("# RGA"):
                break
            if not current_period or not current_species: 
                continue
            if line.startswith("#") or not line:
                continue

            parts = line.split(",")
            try:
                pos = float(parts[2])
                neg = float(parts[3])
            except (IndexError, ValueError):
                continue

            data[current_period][current_species] += (pos + neg)

    return data, periods, species_list

def summarize(data, periods, species_list):
    for period in periods:
        species_data = data.get(period, {})
        total = sum(species_data[s] for s in species_list)
        print(f"{period}:")
        if total == 0:
            print("  (no data)")
            continue
        for s in species_list:
            val = species_data.get(s, 0.0)
            frac = 100.0 * val / total
            print(f"  {s:4s} : {val:12.4f}  ({frac:6.2f}%)")
        print(f"  {'Total':4s} : {total:12.4f}\n")

def main():
    parser = argparse.ArgumentParser(
        description="Sum pos+neg charge per species and fraction by run period"
    )
    parser.add_argument("csv_file", help="Path to clas12_run_info.csv")
    args = parser.parse_args()

    data, periods, species_list = parse_csv(args.csv_file)
    summarize(data, periods, species_list)

if __name__ == "__main__":
    main()