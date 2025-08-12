#!/usr/bin/env python3
import sys
import math

def parse_file(path, Eb_GeV, rows):
    print(f"[INFO] Reading: {path} (Eb={Eb_GeV} GeV)")
    with open(path, "r") as f:
        line_num = 0
        for raw in f:
            line_num += 1
            s = raw.strip()
            if not s:
                continue  #endif
            if s.startswith("#") or s.lower().startswith("x_b"):
                continue  #endif

            parts = s.split()
            if len(parts) < 7:
                print(f"[WARN] Skipping malformed line {line_num} in {path}: {s}")
                continue  #endif

            try:
                xb = float(parts[0])
                q2 = float(parts[1])
                t_mag = float(parts[2])   # input is "-t" (positive); we keep the magnitude
                phi = float(parts[3])
                A = float(parts[4])
                err_stat = float(parts[5])
                err_syst = float(parts[6])
            except ValueError:
                print(f"[WARN] Non-numeric line {line_num} in {path}: {s}")
                continue  #endif

            sigA = math.sqrt(err_stat * err_stat + err_syst * err_syst)
            rows.append((phi, q2, xb, t_mag, Eb_GeV, A, sigA))
        #endfor
    print(f"[INFO] Parsed {len(rows)} total rows so far.")
#enddef

def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <input_file_1> <input_file_2> <output_file>")
        sys.exit(1)
    #endif

    in1, in2, out_path = sys.argv[1], sys.argv[2], sys.argv[3]

    rows = []
    parse_file(in1, 6.535, rows)
    parse_file(in2, 7.546, rows)

    # # Optional: sort for stable ordering (phi, then q2, xb, t)
    # rows.sort(key=lambda r: (r[0], r[1], r[2], r[3]))

    with open(out_path, "w") as out:
        out.write("# phi(deg) q2(GeV2) xb t(GeV2) Eb(GeV) A sigA\n")
        for phi, q2, xb, t_mag, Eb, A, sigA in rows:
            out.write(
                "{:0.6f} {:0.6f} {:0.6f} {:0.6f} {:0.3f} {:0.6f} {:0.6f}\n".format(
                    phi, q2, xb, -t_mag, Eb, A, sigA
                )
            )
        #endfor

    print(f"[OK] Wrote {len(rows)} rows to: {out_path}")
#enddef

if __name__ == "__main__":
    main()
#endif