#!/usr/bin/env python
import csv
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(
        description="Convert a CSV into a whitespace-delimited text file"
    )
    parser.add_argument("input_csv", help="Path to input CSV file")
    parser.add_argument("output_txt", help="Path to output text file")
    args = parser.parse_args()

    eb = 10.604  # beam energy in GeV

    with open(args.input_csv, "r") as csvfile:
        reader = csv.DictReader(csvfile)
        # build a map from lowercase/stripped names to the actual header
        field_map = {fn.strip().lower(): fn for fn in reader.fieldnames}

        # required columns (all lowercased)
        required = [
            "valid bin",
            "phiavg",
            "q2avg",
            "xbavg",
            "t_abs_avg",
            "cross sections, ep->epg, exp",
            "cross sections, ep->epg, exp, stat. unc."
        ]
        missing = [r for r in required if r not in field_map]
        if missing:
            sys.stderr.write("Error: missing columns: {}\n".format(missing))
            sys.stderr.write("Available columns: {}\n".format(reader.fieldnames))
            sys.exit(1)

        # map to actual CSV header names
        valid_col = field_map["valid bin"]
        phi_col   = field_map["phiavg"]
        q2_col    = field_map["q2avg"]
        xb_col    = field_map["xbavg"]
        t_col     = field_map["t_abs_avg"]
        A_col     = field_map["cross sections, ep->epg, exp"]
        sigA_col  = field_map["cross sections, ep->epg, exp, stat. unc."]

        with open(args.output_txt, "w") as fout:
            # header line
            fout.write("# phi(rad) q2(GeV2) xb t(GeV2) Eb(GeV) A sigA\n")

            # process rows
            for row in reader:
                if row[valid_col].strip() == "1":
                    phi  = float(row[phi_col])
                    q2   = float(row[q2_col])
                    xb   = float(row[xb_col])
                    t    = -float(row[t_col])
                    A    = float(row[A_col])
                    sigA = float(row[sigA_col])

                    # write out seven columns rounded to five decimal places
                    fout.write("{:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f}\n".format(
                        phi, q2, xb, t, eb, A, sigA
                    ))
                #endif
            #endfor

if __name__ == "__main__":
    main()