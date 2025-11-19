#!/usr/bin/env python3
# km15_cli.py
# Usage:
#   python3 km15_cli.py xB Q2 t_pos phi_deg Ebeam helicity [mode]
#
# Args:
#   xB, Q2, t_pos, phi_deg, Ebeam: floats
#   helicity: -1, 0, or +1  (beam polarization state)
#   mode (optional): one of:
#       "XS"  -> total cross section for the given helicity (default)
#       "XSD" -> helicity difference sigma(+) - sigma(-) at same kinematics
#
# Prints a single float on stdout. On any error, prints "0.0".
#
# Notes:
# - We keep your historical conventions:
#     t = -abs(t_pos)
#     phi_trento = pi - phi_rad
# - If you request mode="XSD", we will compute sigma(+) - sigma(-)
#   by making two calls with in1polarization=+1 and -1.
# - For mode="XS" with helicity=0, you get the unpolarized prediction.

import sys
import math
import warnings

# Suppress the noisy RuntimeWarnings (invalid value in sqrt, etc.)
warnings.filterwarnings("ignore", category=RuntimeWarning)

try:
    import gepard as g
    from gepard.fits import th_KM15
except Exception:
    # If KM15 is not available, just return 0.0 so the caller can fall back.
    print("0.0")
    sys.exit(0)


def km15_value(xB, Q2, t_pos, phi_deg, Ebeam, helicity):
    """
    Safe wrapper around KM15 prediction.
    Returns 0.0 if KM15 gives a non-finite or negative value, or if
    the kinematics are obviously unphysical.
    """
    try:
        # Quick sanity checks to avoid clearly unphysical calls.
        if not (0.0 < xB < 1.0):
            return 0.0
        if Q2 <= 0.0:
            return 0.0
        if t_pos <= 0.0:
            return 0.0

        t_km15     = -abs(t_pos)
        phi_rad    = math.radians(phi_deg)
        phi_trento = math.pi - phi_rad

        pt = g.DataPoint(
            xB        = xB,
            t         = t_km15,
            Q2        = Q2,
            phi       = phi_trento,
            observable= "XS",
            frame     = "trento",
            process   = "ep2epgamma",
            exptype   = "fixed target",
            in1energy = Ebeam,
            in1charge = -1,
            in1polarization = helicity  # -1, 0, +1
        )
        pt.prepare()

        # Compute prediction
        val = th_KM15.predict(pt)

        # Guard against NaN, inf, or negative values
        if not math.isfinite(val) or val <= 0.0:
            return 0.0

        return val
    except Exception:
        return 0.0
#enddef


def main():
    if len(sys.argv) < 7:
        print("0.0")
        return

    try:
        xB      = float(sys.argv[1])
        Q2      = float(sys.argv[2])
        t_pos   = float(sys.argv[3])
        phi_deg = float(sys.argv[4])
        Ebeam   = float(sys.argv[5])
        helicity= int(sys.argv[6])   # -1, 0, +1
    except Exception:
        print("0.0")
        return

    mode = "XS"
    if len(sys.argv) >= 8 and sys.argv[7]:
        mode = str(sys.argv[7]).upper()

    try:
        if mode == "XSD":
            # sigma(+) - sigma(-)
            sp = km15_value(xB, Q2, t_pos, phi_deg, Ebeam, +1)
            sm = km15_value(xB, Q2, t_pos, phi_deg, Ebeam, -1)
            val = sp - sm
            # Even here, make sure final value is finite
            if not math.isfinite(val):
                val = 0.0
            print(f"{val:.12g}")
            return

        # Default: total cross section for requested helicity
        val = km15_value(xB, Q2, t_pos, phi_deg, Ebeam, helicity)
        if not math.isfinite(val):
            val = 0.0
        print(f"{val:.12g}")
    except Exception:
        print("0.0")
#enddef


if __name__ == "__main__":
    main()