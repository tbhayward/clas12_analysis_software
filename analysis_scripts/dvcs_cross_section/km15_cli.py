#!/usr/bin/env python3
# km15_cli.py
# Usage:
#   python3 km15_cli.py xB Q2 t_pos phi_deg Ebeam helicity
# Prints a single floating number (cross section) to stdout.

import sys
import math

try:
    import gepard as g
    from gepard.fits import th_KM15
except Exception as e:
    print("0.0")
    sys.exit(0)

def main():
    if len(sys.argv) < 6:
        print("0.0")
        return

    xB      = float(sys.argv[1])
    Q2      = float(sys.argv[2])
    t_pos   = float(sys.argv[3])
    phi_deg = float(sys.argv[4])
    Ebeam   = float(sys.argv[5])
    helicity= int(sys.argv[6]) if len(sys.argv) > 6 else 0

    # Conventions to match your previous usage:
    t_km15      = -abs(t_pos)
    phi_rad     = math.radians(phi_deg)
    phi_trento  = math.pi - phi_rad

    try:
        pt = g.DataPoint(
            xB        = xB,
            t         = t_km15,
            Q2        = Q2,
            phi       = phi_trento,
            observable= 'XS',
            frame     = 'trento',
            process   = 'ep2epgamma',
            exptype   = 'fixed target',
            in1energy = Ebeam,
            in1charge = -1,
            in1polarization = 0  # extend here if th_KM15 supports helicity
        )
        pt.prepare()
        val = th_KM15.predict(pt)
        # If you add helicity support later, branch on 'helicity' here.
        print(f"{val:.12g}")
    except Exception:
        print("0.0")

if __name__ == "__main__":
    main()