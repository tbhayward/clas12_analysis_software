#!/usr/bin/env python3
"""Evaluate BH and KM15 only at physically occupied generated-MC cells.

Input representative points come from generated MC itself.  Empty Cartesian
cells are never evaluated.  Invalid Gepard points are written explicitly as
invalid and are never silently assigned unit weight by the C++ study.
"""
from __future__ import annotations
import argparse, math, warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
import pandas as pd

def make_point(g,xB,Q2,t_abs,phi_deg,ebeam):
    pt=g.DataPoint(xB=float(xB),t=-abs(float(t_abs)),Q2=float(Q2),
        phi=math.pi-math.radians(float(phi_deg)),observable="XS",frame="trento",
        process="ep2epgamma",exptype="fixed target",in1energy=float(ebeam),
        in1charge=-1,in1polarization=0,in2particle="p")
    pt.prepare(); return pt

def evaluate(task):
    energy_tag,period,row,subcell,n_gen,sumw,xb,q2,tabs,phi=task
    ebeam=10.200 if str(energy_tag)=="10.2" else 10.604
    bh=km=np.nan; valid=False; err=""
    try:
        import gepard as g
        from gepard.fits import th_KM15
        with warnings.catch_warnings():
            warnings.simplefilter("ignore",RuntimeWarning)
            pt=make_point(g,xb,q2,tabs,phi,ebeam)
            pref=float(th_KM15.PreFacSigma(pt)); b=pref*float(th_KM15.TBH2unp(pt))
            k=b+pref*float(th_KM15.TINTunp(pt))+pref*float(th_KM15.TDVCS2unp(pt))
        if np.isfinite(b) and b>0 and np.isfinite(k) and k>0:
            bh,km,valid=float(b),float(k),True
        else: err="nonfinite_or_nonpositive_model"
    except Exception as exc: err=type(exc).__name__+":"+str(exc)[:120]
    return dict(energy_tag=energy_tag,period=period,row=int(row),subcell=int(subcell),
        n_gen=int(n_gen),sumw_gen=float(sumw),xB=float(xb),Q2=float(q2),t_abs=float(tabs),
        phi_deg=float(phi),ebeam=ebeam,bh_xs=bh,km15_xs=km,model_valid=int(valid),error=err)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--occupancy-csv',required=True); ap.add_argument('--output',required=True)
    ap.add_argument('--workers',type=int,default=8); ap.add_argument('--force',action='store_true'); a=ap.parse_args()
    out=Path(a.output)
    if out.exists() and not a.force: print(f"[acceptance-bh-grid] Reusing {out}"); return 0
    d=pd.read_csv(a.occupancy_csv)
    req=['energy_tag','period','row','subcell','n_gen','sumw_gen','xB_mean','Q2_mean','t_abs_mean','phi_deg_circular_mean']
    miss=[x for x in req if x not in d.columns]
    if miss: raise RuntimeError(f"Missing occupancy columns: {miss}")
    tasks=[tuple(x) for x in d[req].itertuples(index=False,name=None) if int(x[4])>0 and float(x[5])>0]
    print(f"[acceptance-bh-grid] Evaluating {len(tasks)} occupied generated-MC cells with {a.workers} worker(s).")
    with ProcessPoolExecutor(max_workers=max(1,min(a.workers,16))) as ex: rows=list(ex.map(evaluate,tasks,chunksize=32))
    z=pd.DataFrame(rows); out.parent.mkdir(parents=True,exist_ok=True); z.to_csv(out,index=False)
    nv=int(z.model_valid.sum()) if len(z) else 0
    print(f"[acceptance-bh-grid] Wrote {len(z)} occupied cells to {out}; valid model cells={nv}, invalid={len(z)-nv}.")
    return 0
if __name__=='__main__': raise SystemExit(main())
