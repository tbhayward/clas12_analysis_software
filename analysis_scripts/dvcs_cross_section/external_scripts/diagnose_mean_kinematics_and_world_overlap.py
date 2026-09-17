#!/usr/bin/env python3
"""Diagnose pass-1/pass-2 bin-mean shifts and CLAS6 overlap.

Primary questions:
  1) Did the mean (xB,Q2,|t|,phi) values move materially between the original
     pass-1 all_bin_v3 analysis and pass-2, despite identical bin definitions?
  2) How does the SAME-PERIOD Fa18 pass-2 result compare with pass-1?
  3) How does combined 10.6-GeV pass-2 compare with pass-1?
  4) Where does pass-2 overlap Jo 2015 and Saylor 2018, especially at low |t|?
  5) Are observed cross-section differences correlated with small changes in
     the mean kinematics, and what does KM15 predict for those finite shifts?

This script never changes analysis products. It writes CSV tables and PNGs to
output/mean_kinematics_overlap_study/.
"""
from __future__ import annotations
import argparse, math, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
if str(HERE) not in sys.path: sys.path.insert(0, str(HERE))

import extract_emff_from_dvcs_bh as emff
import compare_world_dvcs_cross_sections as worldmod

PASS2 = ROOT / "output" / "csvs" / "dvcs_pass2_analysis.csv"
PASS1_LEGACY = ROOT / "imports" / "all_bin_v3.csv"
SAYLOR = HERE / "import" / "saylor_CLAS6.txt"
OUT = ROOT / "output" / "mean_kinematics_overlap_study"

P2_XS = {
    "Fa18": "normed cross sections, ep->epg, exp, Fa18, unpol",
    "10.6 GeV": "normed cross sections, ep->epg, exp, 10.6 GeV, unpol",
}
P2_MEAN = lambda v, lab: f"{v}, {lab}"


def parse_tuple(v):
    try:
        s=str(v).strip()
        if not (s.startswith('(') and s.endswith(')')): return (np.nan,np.nan,np.nan)
        a=[float(x.strip()) for x in s[1:-1].split(',')]
        return tuple((a+[np.nan]*3)[:3])
    except Exception: return (np.nan,np.nan,np.nan)


def circdiff(a,b): return ((np.asarray(a)-np.asarray(b)+180.0)%360.0)-180.0


def finite_stats(x):
    a=np.asarray(x,float); a=a[np.isfinite(a)]
    if not len(a): return dict(N=0,median=np.nan,p16=np.nan,p84=np.nan,p95abs=np.nan,maxabs=np.nan)
    return dict(N=len(a), median=float(np.median(a)), p16=float(np.percentile(a,16)),
                p84=float(np.percentile(a,84)), p95abs=float(np.percentile(np.abs(a),95)),
                maxabs=float(np.max(np.abs(a))))


def load_pass1_legacy(path):
    d=pd.read_csv(path,low_memory=False)
    if "valid bin" in d: d=d[pd.to_numeric(d["valid bin"],errors="coerce")==1].copy()
    xs=[]; st=[]
    for v in d["cross sections, ep->epg, exp"]:
        # legacy file stores central xs as scalar
        try: xs.append(float(v))
        except: xs.append(np.nan)
    st=pd.to_numeric(d.get("cross sections, ep->epg, exp, stat. unc.",np.nan),errors="coerce")
    d["xs_p1"]=xs; d["stat_p1"]=st
    d["row_p1"]=d.index.astype(int)
    return d


def load_pass2(path,label):
    d=pd.read_csv(path,low_memory=False)
    req=["bin index","Bin Name","xBmin","xBmax","Q2min","Q2max","t_abs_min","t_abs_max","phimin","phimax",
         P2_XS[label],P2_MEAN("xBavg",label),P2_MEAN("Q2avg",label),P2_MEAN("t_abs_avg",label),P2_MEAN("phiavg",label)]
    miss=[c for c in req if c not in d]
    if miss: raise RuntimeError(f"Pass-2 CSV missing columns for {label}: {miss}")
    vals=[parse_tuple(v) for v in d[P2_XS[label]]]
    d=d.copy(); d["xs_p2"]=[v[0] for v in vals]; d["stat_p2"]=[v[1] for v in vals]
    d["xB_p2"]=pd.to_numeric(d[P2_MEAN("xBavg",label)],errors="coerce")
    d["Q2_p2"]=pd.to_numeric(d[P2_MEAN("Q2avg",label)],errors="coerce")
    d["t_p2"]=pd.to_numeric(d[P2_MEAN("t_abs_avg",label)],errors="coerce")
    d["phi_p2"]=np.mod(pd.to_numeric(d[P2_MEAN("phiavg",label)],errors="coerce"),360)
    return d[np.isfinite(d.xs_p2)&(d.xs_p2>0)].copy()


def exact_pass1_pass2(p1,p2,label):
    # all_bin_v3 and pass2 are intentionally based on the same bin definitions.
    # Match using Bin Name + exact bin edges; within a 3D bin choose same phi-bin edges.
    keys=["Bin Name","xBmin","xBmax","Q2min","Q2max","t_abs_min","t_abs_max","phimin","phimax"]
    a=p1.copy(); b=p2.copy()
    for c in keys:
        a[c]=pd.to_numeric(a[c],errors="coerce").round(10); b[c]=pd.to_numeric(b[c],errors="coerce").round(10)
    m=b.merge(a[keys+["row_p1","xBavg","Q2avg","t_abs_avg","phiavg","xs_p1","stat_p1"]],on=keys,how="inner",validate="one_to_one")
    m=m.rename(columns={"xBavg":"xB_p1","Q2avg":"Q2_p1","t_abs_avg":"t_p1","phiavg":"phi_p1"})
    m["dxB"]=m.xB_p2-m.xB_p1; m["dQ2"]=m.Q2_p2-m.Q2_p1; m["dt"]=m.t_p2-m.t_p1; m["dphi"]=circdiff(m.phi_p2,m.phi_p1)
    for v,lo,hi in [("xB","xBmin","xBmax"),("Q2","Q2min","Q2max"),("t","t_abs_min","t_abs_max")]:
        w=pd.to_numeric(m[hi],errors="coerce")-pd.to_numeric(m[lo],errors="coerce")
        m[f"d{v}_over_bin_width"]=m[{"xB":"dxB","Q2":"dQ2","t":"dt"}[v]]/w
    pw=(pd.to_numeric(m.phimax,errors="coerce")-pd.to_numeric(m.phimin,errors="coerce"))%360
    pw=np.where(pw==0,360,pw); m["dphi_over_bin_width"]=m.dphi/pw
    m["p2_over_p1_direct"]=m.xs_p2/m.xs_p1
    m["comparison"]=label
    return m


def km15_one(xB,Q2,t,phi,ebeam,dataset="pass2"):
    from types import SimpleNamespace
    r=SimpleNamespace(dataset=dataset,xB=float(xB),Q2=float(Q2),t_abs=float(t),phi_deg=float(phi))
    return float(worldmod.evaluate_one_km15(emff,r,float(ebeam))["km15_ep"])


def add_km15_exact(m, ebeam=10.604):
    p1=[];p2=[]
    for i,r in enumerate(m.itertuples(index=False)):
        if i%100==0: print(f"  KM15 exact-bin means {i}/{len(m)}",flush=True)
        try:
            p1.append(km15_one(r.xB_p1,r.Q2_p1,r.t_p1,r.phi_p1,ebeam,"lee2026"))
            p2.append(km15_one(r.xB_p2,r.Q2_p2,r.t_p2,r.phi_p2,ebeam,"pass2"))
        except Exception: p1.append(np.nan);p2.append(np.nan)
    m=m.copy(); m["km15_p1_mean"]=p1; m["km15_p2_mean"]=p2
    m["km15_p1_to_p2_factor"]=m.km15_p2_mean/m.km15_p1_mean
    m["km15_mean_shift_pct"]=100*(m.km15_p1_to_p2_factor-1)
    m["p2_over_p1_after_km15_mean_transport"]=m.xs_p2/(m.xs_p1*m.km15_p1_to_p2_factor)
    return m


def nearest_world(pass2, ref, name, ebeam_ref, max_norm_dist=2.0):
    # Purely geometric matching first. Scales are explicit and approximately the
    # pass-2 bin granularity; they are NOT model based.
    scales=dict(xB=0.03,Q2=0.30,t=0.10,phi=15.0)
    rows=[]; used=set()
    for _,p in pass2.iterrows():
        dx=(ref.xB.to_numpy(float)-p.xB_p2)/scales["xB"]
        dq=(ref.Q2.to_numpy(float)-p.Q2_p2)/scales["Q2"]
        dt=(ref.t_abs.to_numpy(float)-p.t_p2)/scales["t"]
        dp=circdiff(ref.phi_deg.to_numpy(float),p.phi_p2)/scales["phi"]
        dist=np.sqrt(dx*dx+dq*dq+dt*dt+dp*dp)
        order=np.argsort(dist)
        j=next((int(k) for k in order if int(k) not in used),None)
        if j is None or not np.isfinite(dist[j]) or dist[j]>max_norm_dist: continue
        used.add(j); q=ref.iloc[j]
        rows.append(dict(comparison=name,pass2_bin_index=p["bin index"],pass2_Bin_Name=p["Bin Name"],
            p2_xB=p.xB_p2,p2_Q2=p.Q2_p2,p2_t=p.t_p2,p2_phi=p.phi_p2,p2_xs=p.xs_p2,p2_stat=p.stat_p2,
            ref_id=q.get("point_id",q.get("source_row",j)),ref_xB=q.xB,ref_Q2=q.Q2,ref_t=q.t_abs,ref_phi=q.phi_deg,
            ref_xs=q.xs,ref_stat=q.stat_abs,ref_point_unc=q.point_unc_abs,
            dxB=p.xB_p2-q.xB,dQ2=p.Q2_p2-q.Q2,dt=p.t_p2-q.t_abs,dphi=float(circdiff(p.phi_p2,q.phi_deg)),geom_distance=float(dist[j]),
            p2_over_ref_direct=p.xs_p2/q.xs,ref_ebeam=ebeam_ref))
    out=pd.DataFrame(rows)
    if out.empty:return out
    kref=[];kp2=[]
    for i,r in enumerate(out.itertuples(index=False)):
        if i%100==0: print(f"  KM15 {name} {i}/{len(out)}",flush=True)
        try:
            kref.append(km15_one(r.ref_xB,r.ref_Q2,r.ref_t,r.ref_phi,r.ref_ebeam,"saylor2018" if "Saylor" in name else "jo2015"))
            kp2.append(km15_one(r.p2_xB,r.p2_Q2,r.p2_t,r.p2_phi,10.604,"pass2"))
        except Exception:kref.append(np.nan);kp2.append(np.nan)
    out["km15_ref"]=kref;out["km15_p2"]=kp2;out["km15_ref_to_p2_factor"]=out.km15_p2/out.km15_ref
    out["p2_over_ref_after_km15_transport"]=out.p2_xs/(out.ref_xs*out.km15_ref_to_p2_factor)
    return out


def summary_exact(m):
    rows=[]
    selections={"all":np.ones(len(m),bool),"|t|<0.15":m.t_p2<.15,"0.15<=|t|<0.25":(m.t_p2>=.15)&(m.t_p2<.25),
                "0.25<=|t|<0.40":(m.t_p2>=.25)&(m.t_p2<.40),"0.40<=|t|<0.60":(m.t_p2>=.40)&(m.t_p2<.60),"|t|>=0.60":m.t_p2>=.60}
    for name,mask in selections.items():
        d=m.loc[mask]
        if d.empty:continue
        rows.append(dict(comparison=d.comparison.iloc[0],selection=name,N=len(d),
            median_p2_over_p1=float(np.nanmedian(d.p2_over_p1_direct)),
            median_after_km15=float(np.nanmedian(d.p2_over_p1_after_km15_mean_transport)),
            median_abs_dxB=float(np.nanmedian(abs(d.dxB))),median_abs_dQ2=float(np.nanmedian(abs(d.dQ2))),
            median_abs_dt=float(np.nanmedian(abs(d.dt))),median_abs_dphi=float(np.nanmedian(abs(d.dphi))),
            p95_abs_dxB=float(np.nanpercentile(abs(d.dxB),95)),p95_abs_dQ2=float(np.nanpercentile(abs(d.dQ2),95)),
            p95_abs_dt=float(np.nanpercentile(abs(d.dt),95)),p95_abs_dphi=float(np.nanpercentile(abs(d.dphi),95)),
            median_abs_dxB_binfrac=float(np.nanmedian(abs(d.dxB_over_bin_width))),
            median_abs_dQ2_binfrac=float(np.nanmedian(abs(d.dQ2_over_bin_width))),
            median_abs_dt_binfrac=float(np.nanmedian(abs(d.dt_over_bin_width))),
            median_abs_dphi_binfrac=float(np.nanmedian(abs(d.dphi_over_bin_width))),
            median_abs_km15_shift_pct=float(np.nanmedian(abs(d.km15_mean_shift_pct)))))
    return pd.DataFrame(rows)


def plots_exact(m,outdir,tag):
    outdir.mkdir(parents=True,exist_ok=True)
    specs=[("dxB",r"$\Delta x_B$"),("dQ2",r"$\Delta Q^2$ (GeV$^2$)"),("dt",r"$\Delta |t|$ (GeV$^2$)"),("dphi",r"$\Delta\phi$ (deg)")]
    for col,xlab in specs:
        fig,ax=plt.subplots(figsize=(7.2,5.2)); ax.hist(m[col].dropna(),bins=50); ax.axvline(0,linewidth=1); ax.set_xlabel(xlab);ax.set_ylabel("Matched bins");ax.set_title(f"{tag}: pass-2 minus pass-1 mean kinematics");fig.tight_layout();fig.savefig(outdir/f"{tag}_{col}_distribution.png",dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7.2,5.2)); sc=ax.scatter(m.t_p2,m.p2_over_p1_direct,c=np.abs(m.km15_mean_shift_pct),s=12);ax.axhline(1,linewidth=1);ax.set_xlabel(r"pass-2 mean $|t|$ (GeV$^2$)");ax.set_ylabel("pass-2 / pass-1");cb=fig.colorbar(sc,ax=ax);cb.set_label("|KM15 mean-kinematic effect| (%)");fig.tight_layout();fig.savefig(outdir/f"{tag}_ratio_vs_t_colored_by_km15_shift.png",dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7.2,5.2));ax.scatter(m.km15_mean_shift_pct,100*(m.p2_over_p1_direct-1),s=12);ax.axhline(0,linewidth=1);ax.axvline(0,linewidth=1);ax.set_xlabel("KM15 predicted change from pass-1 mean to pass-2 mean (%)");ax.set_ylabel("Observed pass-2/pass-1 - 1 (%)");fig.tight_layout();fig.savefig(outdir/f"{tag}_observed_vs_predicted_mean_shift.png",dpi=180);plt.close(fig)


def main():
    ap=argparse.ArgumentParser();ap.add_argument("--pass2",type=Path,default=PASS2);ap.add_argument("--pass1",type=Path,default=PASS1_LEGACY);ap.add_argument("--out",type=Path,default=OUT);ap.add_argument("--skip-km15",action="store_true",help="Fast geometry-only pass")
    args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True);(args.out/"tables").mkdir(exist_ok=True);(args.out/"figures").mkdir(exist_ok=True)
    print(f"PASS1: {args.pass1}\nPASS2: {args.pass2}\nOUT: {args.out}",flush=True)
    p1=load_pass1_legacy(args.pass1)
    exacts=[]; sums=[]
    for lab,tag in [("Fa18","fa18_vs_pass1"),("10.6 GeV","combined10p6_vs_pass1")]:
        p2=load_pass2(args.pass2,lab); m=exact_pass1_pass2(p1,p2,lab)
        print(f"[{lab}] exact identical-bin-definition matches: {len(m)}",flush=True)
        if not args.skip_km15:m=add_km15_exact(m)
        else:m["km15_mean_shift_pct"]=np.nan;m["p2_over_p1_after_km15_mean_transport"]=np.nan
        m.to_csv(args.out/"tables"/f"{tag}_points.csv",index=False); exacts.append(m)
        s=summary_exact(m);s.to_csv(args.out/"tables"/f"{tag}_summary.csv",index=False);sums.append(s)
        plots_exact(m,args.out/"figures",tag)

    # Direct mean-shift QA across all exact pass1/pass2 bins, independent of xs.
    qa=pd.concat(sums,ignore_index=True);qa.to_csv(args.out/"tables"/"pass1_pass2_mean_kinematics_summary.csv",index=False)
    print("\nMEAN-KINEMATICS SUMMARY\n",qa.to_string(index=False),flush=True)

    # World comparisons use the same validated loaders as the production world script.
    if not args.skip_km15:
        p2=load_pass2(args.pass2,"10.6 GeV")
        jo=worldmod.canonicalize_jo(emff.load_clas6_gepard_dataset(),emff)
        say=worldmod.canonicalize_saylor(emff.load_saylor_supplement(str(SAYLOR)),emff)
        for ref,name,e in [(jo,"pass2_vs_jo2015",5.75),(say,"pass2_vs_saylor2018",5.75)]:
            q=nearest_world(p2,ref,name,e)
            q.to_csv(args.out/"tables"/f"{name}_matched_points.csv",index=False)
            if q.empty:continue
            cuts={"all":np.ones(len(q),bool),"low_t_<0.31":q.p2_t<.31,"high_t_>=0.31":q.p2_t>=.31,"very_low_t_<0.20":q.p2_t<.20}
            rows=[]
            for sel,mask in cuts.items():
                d=q.loc[mask]
                if d.empty:continue
                rows.append(dict(comparison=name,selection=sel,N=len(d),median_geom_distance=np.nanmedian(d.geom_distance),
                    median_direct_ratio=np.nanmedian(d.p2_over_ref_direct),median_transported_ratio=np.nanmedian(d.p2_over_ref_after_km15_transport),
                    median_abs_dxB=np.nanmedian(abs(d.dxB)),median_abs_dQ2=np.nanmedian(abs(d.dQ2)),median_abs_dt=np.nanmedian(abs(d.dt)),median_abs_dphi=np.nanmedian(abs(d.dphi))))
            pd.DataFrame(rows).to_csv(args.out/"tables"/f"{name}_summary.csv",index=False)
            fig,ax=plt.subplots(figsize=(7.2,5.2));ax.scatter(q.p2_t,q.p2_over_ref_after_km15_transport,s=12);ax.axhline(1,linewidth=1);ax.axvline(.31,linestyle='--',linewidth=1);ax.set_xlabel(r"pass-2 $|t|$ (GeV$^2$)");ax.set_ylabel(f"pass-2 / {name.split('_vs_')[-1]} after KM15 transport");fig.tight_layout();fig.savefig(args.out/"figures"/f"{name}_transported_ratio_vs_t.png",dpi=180);plt.close(fig)
    print("\nDone. No analysis inputs were modified.",flush=True)

if __name__=="__main__": main()
