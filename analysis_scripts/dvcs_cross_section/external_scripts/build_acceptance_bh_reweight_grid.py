#!/usr/bin/env python3
"""
Build the pure-BH subcell grid used by the pass-2 acceptance reweighting study.

The production 4D (xB,Q2,|t|,phi) bin is split 2x2x2x2.  The pure BH
cross section is evaluated at each subcell center with Gepard/KM15's BH term.
Only BH is evaluated; no VGG model call is made.

The absolute normalization is irrelevant.  The C++ study combines the BH value
with the subcell phase-space volume and the nominal generated-MC population to
construct a relative event weight.
"""

from __future__ import annotations
import argparse
import math
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
import numpy as np


def make_point(g, xB: float, Q2: float, t_abs: float,
               phi_deg: float, ebeam: float):
    # Match the convention already used by extract_emff_from_dvcs_bh.py.
    phi_rad = math.radians(float(phi_deg))
    phi_trento = math.pi - phi_rad
    pt = g.DataPoint(
        xB=float(xB),
        t=-abs(float(t_abs)),
        Q2=float(Q2),
        phi=float(phi_trento),
        observable="XS",
        frame="trento",
        process="ep2epgamma",
        exptype="fixed target",
        in1energy=float(ebeam),
        in1charge=-1,
        in1polarization=0,
        in2particle="p",
    )
    pt.prepare()
    return pt


def evaluate(task):
    energy_tag, row_index, subcell, xb, q2, tabs, phi, ebeam, volume = task
    try:
        import gepard as g
        from gepard.fits import th_KM15
    except Exception as exc:
        raise RuntimeError(
            "Could not import gepard/KM15. Run in the same Python environment "
            "used by the existing BH/EMFF scripts."
        ) from exc

    pt = make_point(g, xb, q2, tabs, phi, ebeam)
    th = th_KM15
    pref = float(th.PreFacSigma(pt))
    bh = pref * float(th.TBH2unp(pt))
    return {
        "energy_tag": energy_tag,
        "row": int(row_index),
        "subcell": int(subcell),
        "xB": float(xb),
        "Q2": float(q2),
        "t_abs": float(tabs),
        "phi_deg": float(phi),
        "ebeam": float(ebeam),
        "volume": float(volume),
        "bh_xs": float(bh),
    }


def truthy(v) -> bool:
    return str(v).strip().lower() in {"1", "1.0", "true"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--analysis-csv",
                    default="output/csvs/dvcs_pass2_analysis.csv")
    ap.add_argument("--output",
                    default="output/systematics/acceptance_reweighting/bh_subcell_grid.csv")
    ap.add_argument("--workers", type=int, default=7)
    ap.add_argument("--subdivisions", type=int, default=2)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    out = Path(args.output)
    if out.exists() and not args.force:
        print(f"[acceptance-bh-grid] Reusing existing cache: {out}")
        return 0

    df = pd.read_csv(args.analysis_csv, low_memory=False)
    required = ["xBmin","xBmax","Q2min","Q2max",
                "t_abs_min","t_abs_max","phimin","phimax","valid bin"]
    missing = [c for c in required if c not in df]
    if missing:
        raise RuntimeError(f"Missing required columns: {missing}")

    nsub = int(args.subdivisions)
    if nsub != 2:
        raise RuntimeError("This first implementation intentionally uses subdivisions=2.")

    tasks = []
    energy_specs = [
        ("10.6", 10.604),
        ("10.2", 10.200),
    ]

    for row_index, row in df.iterrows():
        if not truthy(row["valid bin"]):
            continue

        x0,x1 = float(row.xBmin), float(row.xBmax)
        q0,q1 = float(row.Q2min), float(row.Q2max)
        t0,t1 = float(row.t_abs_min), float(row.t_abs_max)
        p0,p1 = float(row.phimin), float(row.phimax)

        dx=(x1-x0)/2.0
        dq=(q1-q0)/2.0
        dt=(t1-t0)/2.0
        dp=(p1-p0)/2.0
        volume = dx*dq*dt*math.radians(dp)

        for energy_tag, ebeam in energy_specs:
            for ix in range(2):
                for iq in range(2):
                    for it in range(2):
                        for ip in range(2):
                            sub = (((ix*2)+iq)*2+it)*2+ip
                            xb=x0+(ix+0.5)*dx
                            q2=q0+(iq+0.5)*dq
                            tabs=t0+(it+0.5)*dt
                            phi=p0+(ip+0.5)*dp
                            tasks.append(
                                (energy_tag,row_index,sub,xb,q2,tabs,phi,
                                 ebeam,volume)
                            )

    print(f"[acceptance-bh-grid] Evaluating {len(tasks)} pure-BH subcell points "
          f"with {args.workers} worker(s).")

    workers=max(1,min(int(args.workers),16))
    with ProcessPoolExecutor(max_workers=workers) as ex:
        rows=list(ex.map(evaluate,tasks,chunksize=64))

    out.parent.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(rows).to_csv(out,index=False)
    print(f"[acceptance-bh-grid] Wrote {len(rows)} rows to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
