#!/usr/bin/env python3
"""Fast KM15 bin-centering evaluator for the pass-2 DVCS analysis.

The expensive part of KM15 evaluation is importing/initializing Gepard.  This
script therefore processes the complete pass-2 CSV in one Python invocation
instead of spawning one Python process for every quadrature point.

The CSV convention used by the C++ extraction is the *multiplicative* bin-
centering correction

    C_bin = sigma_KM15(<xB>,<Q2>,<|t|>,<phi>) / <sigma_KM15>_bin,

so the bin-averaged measured cross section is multiplied by C_bin to report it
at the measured mean kinematics.  This is the inverse of the F_bin definition
(avg/center) sometimes used in prose in the pass-1 note, but is exactly the
factor needed by cross_sections.cpp, which multiplies the stored CSV value.
"""

import argparse
import csv
import math
import os
import sys
import tempfile
import warnings
from multiprocessing import get_context
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore", category=RuntimeWarning)

_G = None
_TH = None


def _init_worker():
    global _G, _TH
    import gepard as g
    from gepard.fits import th_KM15
    _G = g
    _TH = th_KM15


def _km15(point):
    xB, Q2, t_abs, phi_deg, ebeam = point
    try:
        if not (0.0 < xB < 1.0 and Q2 > 0.0 and t_abs > 0.0):
            return 0.0
        pt = _G.DataPoint(
            xB=xB,
            t=-abs(t_abs),
            Q2=Q2,
            phi=math.pi - math.radians(phi_deg),
            observable="XS",
            frame="trento",
            process="ep2epgamma",
            exptype="fixed target",
            in1energy=ebeam,
            in1charge=-1,
            in1polarization=0,
        )
        pt.prepare()
        val = float(_TH.predict(pt))
        return val if math.isfinite(val) and val > 0.0 else 0.0
    except Exception:
        return 0.0
#enddef


def _as_float(row, key):
    try:
        s = row.get(key, "")
        if s is None or str(s).strip() == "":
            return math.nan
        return float(s)
    except Exception:
        return math.nan
#enddef


def _fmt_triple(value, stat=0.0, syst=0.0):
    return f"({value:.8f}, {stat:.8f}, {syst:.8f})"
#enddef


def _quadrature_points(row, ebeam, mean_cols, order):
    xb = _as_float(row, mean_cols[0])
    q2 = _as_float(row, mean_cols[1])
    tt = _as_float(row, mean_cols[2])
    ph = _as_float(row, mean_cols[3])
    center = (xb, q2, tt, ph, ebeam)

    edges = [
        (_as_float(row, "xBmin"), _as_float(row, "xBmax")),
        (_as_float(row, "Q2min"), _as_float(row, "Q2max")),
        (_as_float(row, "t_abs_min"), _as_float(row, "t_abs_max")),
        (_as_float(row, "phimin"), _as_float(row, "phimax")),
    ]
    if not all(math.isfinite(v) for v in center):
        return None
    if not all(math.isfinite(lo) and math.isfinite(hi) and hi > lo for lo, hi in edges):
        return None

    nodes, weights = np.polynomial.legendre.leggauss(order)
    vals = []
    wts = []
    for ix, nx in enumerate(nodes):
        x = 0.5 * ((edges[0][1] - edges[0][0]) * nx + edges[0][1] + edges[0][0])
        for iq, nq in enumerate(nodes):
            q = 0.5 * ((edges[1][1] - edges[1][0]) * nq + edges[1][1] + edges[1][0])
            for it, nt in enumerate(nodes):
                t = 0.5 * ((edges[2][1] - edges[2][0]) * nt + edges[2][1] + edges[2][0])
                for ip, np_ in enumerate(nodes):
                    p = 0.5 * ((edges[3][1] - edges[3][0]) * np_ + edges[3][1] + edges[3][0])
                    vals.append((x, q, t, p, ebeam))
                    wts.append(weights[ix] * weights[iq] * weights[it] * weights[ip])
                #endfor
            #endfor
        #endfor
    #endfor
    return center, vals, np.asarray(wts, dtype=float)
#enddef


def _evaluate_jobs(jobs, workers):
    # Deduplicate points globally.  Repeated center/grid points arise when a
    # kinematic bin is represented by several bookkeeping rows.
    unique = {}
    ordered = []
    for point in jobs:
        key = tuple(round(float(x), 12) for x in point)
        if key not in unique:
            unique[key] = len(ordered)
            ordered.append(point)
        #endif
    #endfor

    if workers <= 1:
        _init_worker()
        values = [_km15(p) for p in ordered]
    else:
        ctx = get_context("spawn")
        chunksize = max(8, len(ordered) // max(1, workers * 32))
        with ctx.Pool(processes=workers, initializer=_init_worker) as pool:
            values = list(pool.imap(_km15, ordered, chunksize=chunksize))
        #endwith
    #endif
    return unique, np.asarray(values, dtype=float)
#enddef


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv_path")
    ap.add_argument("--order", type=int, default=3,
                    help="Gauss-Legendre order per dimension (default: 3 => 81 points/bin)")
    ap.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 1))
    args = ap.parse_args()

    if args.order < 2 or args.order > 6:
        raise SystemExit("--order must be in [2,6]")

    path = Path(args.csv_path)
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)
    #endwith

    required = {
        "xBmin", "xBmax", "Q2min", "Q2max", "t_abs_min", "t_abs_max",
        "phimin", "phimax", "xBavg, 10.6 GeV", "Q2avg, 10.6 GeV",
        "t_abs_avg, 10.6 GeV", "phiavg, 10.6 GeV", "xBavg, Sp19 Inb",
        "Q2avg, Sp19 Inb", "t_abs_avg, Sp19 Inb", "phiavg, Sp19 Inb",
        "Fbin, 10.6 GeV", "Fbin, 10.2 GeV",
    }
    missing = sorted(required.difference(fieldnames))
    if missing:
        raise SystemExit("Missing CSV columns: " + ", ".join(missing))

    groups = [
        (10.6,
         ("xBavg, 10.6 GeV", "Q2avg, 10.6 GeV", "t_abs_avg, 10.6 GeV", "phiavg, 10.6 GeV"),
         "Fbin, 10.6 GeV"),
        (10.2,
         ("xBavg, Sp19 Inb", "Q2avg, Sp19 Inb", "t_abs_avg, Sp19 Inb", "phiavg, Sp19 Inb"),
         "Fbin, 10.2 GeV"),
    ]

    descriptors = []
    all_points = []
    for ir, row in enumerate(rows):
        for ebeam, means, outcol in groups:
            desc = _quadrature_points(row, ebeam, means, args.order)
            if desc is None:
                continue
            center, grid, weights = desc
            descriptors.append((ir, outcol, center, grid, weights))
            all_points.append(center)
            all_points.extend(grid)
        #endfor
    #endfor

    print(f"[bincenter-batch] rows={len(rows)}, populated energy-rows={len(descriptors)}, "
          f"quadrature order={args.order}, raw evaluations={len(all_points)}, workers={args.workers}",
          flush=True)

    lookup, values = _evaluate_jobs(all_points, max(1, args.workers))
    print(f"[bincenter-batch] unique KM15 evaluations={len(values)}", flush=True)

    def value_for(point):
        key = tuple(round(float(x), 12) for x in point)
        return float(values[lookup[key]])
    #enddef

    written = {"Fbin, 10.6 GeV": 0, "Fbin, 10.2 GeV": 0}
    invalid = 0
    for ir, outcol, center, grid, weights in descriptors:
        cval = value_for(center)
        gvals = np.asarray([value_for(p) for p in grid], dtype=float)
        good = np.isfinite(gvals) & (gvals > 0.0) & np.isfinite(weights) & (weights > 0.0)
        if not (math.isfinite(cval) and cval > 0.0 and np.any(good)):
            rows[ir][outcol] = ""
            invalid += 1
            continue
        #endif
        avg = float(np.sum(weights[good] * gvals[good]) / np.sum(weights[good]))
        if not (math.isfinite(avg) and avg > 0.0):
            rows[ir][outcol] = ""
            invalid += 1
            continue
        #endif
        correction = cval / avg
        rows[ir][outcol] = _fmt_triple(correction, 0.0, 0.0)
        written[outcol] += 1
    #endfor

    fd, tmpname = tempfile.mkstemp(prefix=path.name + ".bincenter.", suffix=".tmp", dir=str(path.parent))
    os.close(fd)
    try:
        with open(tmpname, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        #endwith
        os.replace(tmpname, path)
    finally:
        if os.path.exists(tmpname):
            os.unlink(tmpname)
        #endif
    #endtry

    print(f"[bincenter-batch] wrote 10.6={written['Fbin, 10.6 GeV']}, "
          f"10.2={written['Fbin, 10.2 GeV']}, invalid={invalid}", flush=True)
#enddef


if __name__ == "__main__":
    main()
