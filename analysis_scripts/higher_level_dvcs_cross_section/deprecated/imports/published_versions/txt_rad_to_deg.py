#!/usr/bin/env python
import math
import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Convert the first column (phi in radians) to degrees"
    )
    parser.add_argument("input_file", help="Path to the input text file")
    parser.add_argument("output_file", help="Path for the converted output file")
    args = parser.parse_args()

    with open(args.input_file, "r") as fin, open(args.output_file, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            # split on whitespace
            parts = line.strip().split()
            # convert first field from rad to deg
            phi_rad = float(parts[0])
            phi_deg = phi_rad * 180.0 / math.pi
            parts[0] = "{:.5f}".format(phi_deg)
            fout.write(" ".join(parts) + "\n")
        #endfor

if __name__ == "__main__":
    main()