#!/usr/bin/env python3
"""
Fold external pi0/VPK validation information through the *accepted* RGA AAOgen
population exported by eppi0_normalization.cpp.

This script is deliberately downstream of the C++ event selection: it never
reimplements the eppi0 cuts.  It consumes the sparse accepted-AAOgen grid and
external CLAS6/Hall-A validation outputs.

If either external validation output is absent, the corresponding fast
validation script is run automatically.

Primary outputs
---------------
  output/rga_coverage/rga_external_pi0_coverage_summary.csv
  output/rga_coverage/rga_vpk_fold_summary.csv
  output/rga_coverage/rga_candidate_photon_efficiency.csv   (when topology CSV exists)
  output/rga_coverage/png/*.png

Coverage convention
-------------------
"CLAS6 interpolation" means the AAOgen cell center lies inside the 3-D convex
hull of published CLAS6 sigma_U points in (xB,Q2,-t).  "Near CLAS6 boundary"
means it is outside that hull but within a normalized Euclidean distance 0.10
of the nearest CLAS6 point, where each coordinate is normalized by the CLAS6
published range.  Everything else is labeled "Outside CLAS6 support".

The 0.10 boundary band is a transparent diagnostic convention, not a physics
uncertainty.  Results are also written with the raw continuous nearest-point
distance so the threshold can be varied without rerunning AAOgen.

Model folding
-------------
The default folded correction M is a piecewise-linear interpolation of the
point-level sigma_U data/VPK residual in 3-D.  We report CLAS6-only and
CLAS6+Hall-A interpolants separately.  No extrapolation outside the respective
convex hull is performed.  Consequently M_acc is explicitly a supported-yield
average, accompanied by the supported fraction.  This prevents unsupported
high-Q2/edge AAOgen cells from silently acquiring a model correction.
"""
from __future__ import annotations

import argparse
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay, cKDTree, QhullError

HERE = Path(__file__).resolve().parent
DVCS_DIR = HERE.parent
EXT_OUTPUT = HERE / "output"
OUT = EXT_OUTPUT / "rga_coverage"
PNG = OUT / "png"

DEFAULT_AAO = DVCS_DIR / "output" / "data_mc_normalization" / "accepted_aao_population.csv"
DEFAULT_TOPOLOGY = DVCS_DIR / "output" / "data_mc_normalization" / "photon_topology" / "photon_topology_summary.csv"
CLAS6_POINTS = EXT_OUTPUT / "bedlinskiy2014_vpk_point_comparison.csv"
HALLA_POINTS = EXT_OUTPUT / "dlamini2021_vpk_point_comparison.csv"

PERIOD_ORDER = ["Sp18 Inb", "Sp18 Out", "Fa18 Inb", "Fa18 Out", "Sp19 Inb"]
TOPOLOGY_ORDER = ["ALL", "FT-FT", "FT-FD", "FD-FD"]
BOUNDARY_DISTANCE = 0.10
# Local smooth residual field. Distances are in coordinates normalized to the
# combined external-data ranges. The bandwidth is deliberately broad enough
# that the correction varies smoothly rather than following individual points.
LOCAL_BANDWIDTH = 0.22
LOCAL_K = 24
# Additional uncertainty assigned away from measured kinematics. This is not
# used to move the central correction; it makes extrapolation visibly less
# constraining. 0.20 means +20% absolute uncertainty per unit normalized
# nearest-data distance.
DISTANCE_UNCERTAINTY_SLOPE = 0.20


def ensure_validation_outputs() -> None:
    jobs = []
    if not HALLA_POINTS.exists():
        jobs.append(("Hall-A Dlamini", HERE / "validate_aao_pi0_model.py"))
    if not CLAS6_POINTS.exists():
        jobs.append(("CLAS6 Bedlinskiy", HERE / "validate_clas6_bedlinskiy_vpk.py"))
    for label, script in jobs:
        print(f"[rga-coverage] {label} output missing; running {script.name} ...")
        subprocess.run([sys.executable, str(script)], cwd=HERE, check=True)
    for path in (CLAS6_POINTS, HALLA_POINTS):
        if not path.exists():
            raise RuntimeError(f"Expected validation output was not produced: {path}")


def load_external_u_points() -> tuple[pd.DataFrame, pd.DataFrame]:
    c = pd.read_csv(CLAS6_POINTS)
    c = c[(c["observable"] == "U") & c["vpk_kine_valid"].astype(bool)].copy()
    c = c.rename(columns={
        "xB": "xB", "Q2_GeV2": "Q2", "minus_t_GeV2": "minus_t",
        "data_over_vpk": "ratio", "data_over_vpk_err": "ratio_err",
    })
    c = c[["xB", "Q2", "minus_t", "ratio", "ratio_err", "source_table"]]
    c["experiment"] = "CLAS6 Bedlinskiy"

    h = pd.read_csv(HALLA_POINTS)
    h = h[h["vpk_kine_valid"].astype(bool)].copy()
    h = h.rename(columns={
        "xB_mean": "xB", "Q2_GeV2": "Q2", "vpk_minus_t_GeV2": "minus_t",
        "data_over_vpk_U": "ratio", "data_over_vpk_U_err": "ratio_err",
    })
    h = h[["xB", "Q2", "minus_t", "ratio", "ratio_err"]]
    h["source_table"] = "Dlamini2021"
    h["experiment"] = "Hall A Dlamini"

    for d in (c, h):
        for col in ["xB", "Q2", "minus_t", "ratio", "ratio_err"]:
            d[col] = pd.to_numeric(d[col], errors="coerce")
        d.dropna(subset=["xB", "Q2", "minus_t", "ratio"], inplace=True)
    return c.reset_index(drop=True), h.reset_index(drop=True)


def load_aao(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(
            f"Accepted AAOgen population not found: {path}\n"
            "Run: ./dvcs_analysis --eppi0-normalization-only"
        )
    a = pd.read_csv(path)
    required = {
        "period", "photon_topology", "sum_weight", "raw_events",
        "xB_weighted_mean", "Q2_weighted_mean_GeV2", "minus_t_weighted_mean_GeV2",
    }
    missing = sorted(required - set(a.columns))
    if missing:
        raise RuntimeError("Accepted-AAOgen CSV missing columns: " + ", ".join(missing))
    a = a[a.sum_weight > 0].copy()
    a["xB"] = a.xB_weighted_mean
    a["Q2"] = a.Q2_weighted_mean_GeV2
    a["minus_t"] = a.minus_t_weighted_mean_GeV2
    return a


def unique_cloud(df: pd.DataFrame) -> pd.DataFrame:
    # Delaunay/LinearNDInterpolator require unique coordinates.  If duplicate
    # published coordinates occur, combine residuals with inverse-variance weights.
    rows = []
    for key, g in df.groupby(["xB", "Q2", "minus_t"], sort=False):
        e = g.ratio_err.to_numpy(float)
        r = g.ratio.to_numpy(float)
        good = np.isfinite(e) & (e > 0) & np.isfinite(r)
        if good.any():
            w = 1.0 / e[good]**2
            rr = float(np.sum(w*r[good])/np.sum(w))
            ee = float(math.sqrt(1.0/np.sum(w)))
        else:
            rr = float(np.nanmean(r)); ee = np.nan
        rows.append((*key, rr, ee))
    return pd.DataFrame(rows, columns=["xB", "Q2", "minus_t", "ratio", "ratio_err"])


def make_hull(points: np.ndarray):
    try:
        return Delaunay(points, qhull_options="QJ")
    except QhullError as exc:
        raise RuntimeError("Could not construct external-data 3-D convex hull") from exc


def add_clas6_coverage(a: pd.DataFrame, clas6: pd.DataFrame) -> pd.DataFrame:
    out = a.copy()
    xyz_c = clas6[["xB", "Q2", "minus_t"]].to_numpy(float)
    lo = xyz_c.min(axis=0); span = xyz_c.max(axis=0) - lo
    if np.any(span <= 0):
        raise RuntimeError("Degenerate CLAS6 support range")
    c_norm = (xyz_c - lo) / span
    q = out[["xB", "Q2", "minus_t"]].to_numpy(float)
    q_norm = (q - lo) / span
    hull = make_hull(c_norm)
    inside = hull.find_simplex(q_norm) >= 0
    tree = cKDTree(c_norm)
    dist, _ = tree.query(q_norm, k=1)
    out["clas6_nearest_distance_normalized"] = dist
    out["clas6_inside_convex_hull"] = inside
    out["clas6_support_class"] = np.where(
        inside, "CLAS6 interpolation",
        np.where(dist <= BOUNDARY_DISTANCE, "Near CLAS6 boundary", "Outside CLAS6 support")
    )
    return out


def build_interpolator(df: pd.DataFrame):
    d = unique_cloud(df)
    xyz = d[["xB", "Q2", "minus_t"]].to_numpy(float)
    vals = d.ratio.to_numpy(float)
    return LinearNDInterpolator(xyz, vals, fill_value=np.nan, rescale=True), d


def add_model_interpolations(a: pd.DataFrame, clas6: pd.DataFrame, halla: pd.DataFrame) -> pd.DataFrame:
    out = a.copy()
    q = out[["xB", "Q2", "minus_t"]].to_numpy(float)
    f_c, _ = build_interpolator(clas6)
    f_ch, _ = build_interpolator(pd.concat([clas6, halla], ignore_index=True))
    out["M_clas6"] = np.asarray(f_c(q), float)
    out["M_clas6_halla"] = np.asarray(f_ch(q), float)
    return out


def add_local_external_field(a: pd.DataFrame, clas6: pd.DataFrame, halla: pd.DataFrame) -> pd.DataFrame:
    """Smooth external data/VPK residual with continuous support diagnostics.

    This is intentionally not a hard-hull extrapolator.  For each AAOgen cell,
    use nearby CLAS6+Hall-A measurements in normalized (xB,Q2,-t), weighted by
    their quoted uncertainty and a Gaussian distance kernel.  The reported
    uncertainty contains the local quoted errors, local point-to-point scatter,
    and a term that grows continuously with distance from the nearest external
    measurement.
    """
    out=a.copy()
    ext=pd.concat([clas6,halla],ignore_index=True)
    xyz=ext[["xB","Q2","minus_t"]].to_numpy(float)
    q=out[["xB","Q2","minus_t"]].to_numpy(float)
    lo=xyz.min(axis=0); span=xyz.max(axis=0)-lo
    span=np.where(span>0,span,1.0)
    xn=(xyz-lo)/span; qn=(q-lo)/span
    tree=cKDTree(xn)
    k=min(LOCAL_K,len(ext))
    dist,idx=tree.query(qn,k=k)
    if k==1:
        dist=dist[:,None]; idx=idx[:,None]
    ratio=ext.ratio.to_numpy(float); err=ext.ratio_err.to_numpy(float)
    vals=ratio[idx]; errs=err[idx]
    # Avoid a point with tiny quoted uncertainty dominating the field.
    finite_err=err[np.isfinite(err)&(err>0)]
    err_floor=float(np.nanmedian(finite_err))*0.35 if len(finite_err) else 0.02
    ee=np.where(np.isfinite(errs)&(errs>0),np.maximum(errs,err_floor),err_floor)
    kw=np.exp(-0.5*(dist/LOCAL_BANDWIDTH)**2)/(ee**2)
    sw=np.sum(kw,axis=1)
    M=np.sum(kw*vals,axis=1)/sw
    # Error of local weighted mean plus observed local residual scatter.
    meas=np.sqrt(1.0/np.maximum(sw,1e-300))
    scatter=np.sqrt(np.sum(kw*(vals-M[:,None])**2,axis=1)/np.maximum(sw,1e-300))
    neff=(sw**2)/np.maximum(np.sum(kw**2,axis=1),1e-300)
    scatter_mean=scatter/np.sqrt(np.maximum(neff,1.0))
    dnear=dist[:,0]
    extrap=DISTANCE_UNCERTAINTY_SLOPE*dnear
    sigma=np.sqrt(meas**2+scatter_mean**2+extrap**2)
    out["external_nearest_distance_normalized"]=dnear
    out["external_local_M"]=M
    out["external_local_M_unc"]=sigma
    out["external_local_n_eff"]=neff
    return out


def local_fold_summary(a: pd.DataFrame) -> pd.DataFrame:
    rows=[]
    for (p,t),g in a.groupby(["period","photon_topology"],sort=False):
        w=g.sum_weight.to_numpy(float); W=float(w.sum())
        M=g.external_local_M.to_numpy(float); S=g.external_local_M_unc.to_numpy(float)
        d=g.external_nearest_distance_normalized.to_numpy(float)
        mean=float(np.sum(w*M)/W)
        # Conservative population-level model uncertainty: retain the weighted
        # cell uncertainty rather than allowing huge MC statistics to average it away.
        unc=float(np.sqrt(np.sum(w*S**2)/W))
        qd=weighted_quantile(d,w,[.16,.5,.84,.95])
        rows.append({"period":p,"photon_topology":t,"sum_weight":W,
                     "local_M_acc":mean,"local_M_acc_unc":unc,
                     "nearest_distance_p16":qd[0],"nearest_distance_median":qd[1],
                     "nearest_distance_p84":qd[2],"nearest_distance_p95":qd[3]})
    return pd.DataFrame(rows)


def robustness_efficiencies(local_fold: pd.DataFrame, topology_path: Path) -> pd.DataFrame:
    if not topology_path.exists(): return pd.DataFrame()
    top=pd.read_csv(topology_path)
    top=top[(top.comparison=="raw") & top.photon_topology.isin(["FT-FT","FD-FD"])].copy()
    z=top.merge(local_fold,on=["period","photon_topology"],how="left")
    z["candidate_single_photon_ratio_local"]=np.sqrt(z.data_over_mc/z.local_M_acc)
    # Model correction needed for no photon inefficiency: R = M*r^2 -> M_required=R at r=1.
    z["M_required_for_unit_photon_efficiency"]=z.data_over_mc
    z["required_M_shift_from_local"]=z.M_required_for_unit_photon_efficiency-z.local_M_acc
    z["required_shift_in_local_sigma"]=np.abs(z.required_M_shift_from_local)/z.local_M_acc_unc
    # How much VPK would have to overpredict data in the relevant accepted population.
    z["required_vpk_overprediction_percent"]=100.0*(1.0-z.M_required_for_unit_photon_efficiency)
    return z[["period","photon_topology","data_over_mc","stat_err","local_M_acc","local_M_acc_unc",
              "candidate_single_photon_ratio_local","M_required_for_unit_photon_efficiency",
              "required_M_shift_from_local","required_shift_in_local_sigma","required_vpk_overprediction_percent",
              "nearest_distance_median","nearest_distance_p84","nearest_distance_p95"]].sort_values(["period","photon_topology"])


def weighted_quantile(v, w, qs):
    v=np.asarray(v,float); w=np.asarray(w,float); qs=np.asarray(qs,float)
    m=np.isfinite(v)&np.isfinite(w)&(w>0)
    if not m.any(): return np.full(len(qs),np.nan)
    v=v[m]; w=w[m]; o=np.argsort(v); v=v[o]; w=w[o]
    c=np.cumsum(w)-0.5*w; c/=np.sum(w)
    return np.interp(qs,c,v)


def coverage_summary(a: pd.DataFrame) -> pd.DataFrame:
    rows=[]
    for (p,t),g in a.groupby(["period","photon_topology"],sort=False):
        W=float(g.sum_weight.sum())
        q2=weighted_quantile(g.Q2,g.sum_weight,[.16,.5,.84])
        xb=weighted_quantile(g.xB,g.sum_weight,[.16,.5,.84])
        mt=weighted_quantile(g.minus_t,g.sum_weight,[.16,.5,.84])
        row={"period":p,"photon_topology":t,"sum_weight":W,"raw_events":int(g.raw_events.sum()),
             "xB_p16":xb[0],"xB_median":xb[1],"xB_p84":xb[2],
             "Q2_p16_GeV2":q2[0],"Q2_median_GeV2":q2[1],"Q2_p84_GeV2":q2[2],
             "minus_t_p16_GeV2":mt[0],"minus_t_median_GeV2":mt[1],"minus_t_p84_GeV2":mt[2]}
        for label in ["CLAS6 interpolation","Near CLAS6 boundary","Outside CLAS6 support"]:
            key={"CLAS6 interpolation":"clas6_interpolation_fraction",
                 "Near CLAS6 boundary":"clas6_near_boundary_fraction",
                 "Outside CLAS6 support":"clas6_outside_fraction"}[label]
            row[key]=float(g.loc[g.clas6_support_class==label,"sum_weight"].sum()/W) if W>0 else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def fold_summary(a: pd.DataFrame) -> pd.DataFrame:
    rows=[]
    for (p,t),g in a.groupby(["period","photon_topology"],sort=False):
        W=float(g.sum_weight.sum())
        row={"period":p,"photon_topology":t,"sum_weight":W}
        for col,tag in [("M_clas6","clas6"),("M_clas6_halla","clas6_halla")]:
            m=np.isfinite(g[col].to_numpy(float)); ws=g.sum_weight.to_numpy(float)
            supported=float(ws[m].sum())
            row[f"{tag}_supported_weight_fraction"] = supported/W if W>0 else np.nan
            row[f"{tag}_M_acc_supported"] = float(np.sum(ws[m]*g.loc[m,col])/supported) if supported>0 else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def candidate_efficiencies(fold: pd.DataFrame, topology_path: Path) -> pd.DataFrame:
    if not topology_path.exists():
        print(f"[rga-coverage] topology summary not found; skipping candidate photon efficiencies: {topology_path}")
        return pd.DataFrame()
    top=pd.read_csv(topology_path)
    top=top[top.comparison=="raw"].copy()
    top=top[top.photon_topology.isin(["FT-FT","FD-FD"])].copy()
    z=top.merge(fold,on=["period","photon_topology"],how="left")
    for tag in ["clas6","clas6_halla"]:
        M=z[f"{tag}_M_acc_supported"]
        frac=z[f"{tag}_supported_weight_fraction"]
        z[f"candidate_single_photon_ratio_{tag}"]=np.sqrt(z.data_over_mc/M)
        # Flag rather than hide candidates based on incomplete model support.
        z[f"candidate_{tag}_support_ok_80pct"] = frac >= 0.80
    keep=["period","photon_topology","data_over_mc","stat_err",
          "clas6_supported_weight_fraction","clas6_M_acc_supported","candidate_single_photon_ratio_clas6","candidate_clas6_support_ok_80pct",
          "clas6_halla_supported_weight_fraction","clas6_halla_M_acc_supported","candidate_single_photon_ratio_clas6_halla","candidate_clas6_halla_support_ok_80pct"]
    return z[keep].sort_values(["period","photon_topology"])


def interference_closure(clas6_path: Path) -> pd.DataFrame:
    """Summarize U/TT/LT closure without treating unstable TT/LT ratios as normalizations."""
    p = pd.read_csv(clas6_path)
    rows=[]
    for obs in ["U","TT","LT"]:
        g=p[(p.observable==obs)&p.vpk_kine_valid].copy()
        err=pd.to_numeric(g.data_totalerr,errors="coerce").to_numpy(float)
        delta=(pd.to_numeric(g.data_value,errors="coerce")-pd.to_numeric(g.vpk_value,errors="coerce")).to_numpy(float)
        good=np.isfinite(delta)&np.isfinite(err)&(err>0)
        pull=delta[good]/err[good]
        ratio=pd.to_numeric(g.data_over_vpk,errors="coerce").to_numpy(float)
        ratio=ratio[np.isfinite(ratio)]
        rows.append({"observable":obs,"n_points":int(good.sum()),
                     "pull_mean":float(np.mean(pull)),"pull_rms":float(np.sqrt(np.mean(pull**2))),
                     "data_over_vpk_median":float(np.median(ratio)) if len(ratio) else np.nan,
                     "data_over_vpk_p16":float(np.percentile(ratio,16)) if len(ratio) else np.nan,
                     "data_over_vpk_p84":float(np.percentile(ratio,84)) if len(ratio) else np.nan})
    return pd.DataFrame(rows)


def plot_interference_pulls(clas6_path: Path):
    p=pd.read_csv(clas6_path)
    fig,ax=plt.subplots(figsize=(9,5.8))
    for obs,marker in [("U","o"),("TT","s"),("LT","^")]:
        g=p[(p.observable==obs)&p.vpk_kine_valid].copy()
        pull=(g.data_value-g.vpk_value)/g.data_totalerr
        ax.scatter(g.Q2_GeV2,pull,s=22,marker=marker,alpha=.75,label=obs)
    ax.axhline(0,linewidth=1); ax.axhline(2,linewidth=.8,linestyle="--"); ax.axhline(-2,linewidth=.8,linestyle="--")
    ax.set_xlabel(r"$Q^2$ (GeV$^2$)"); ax.set_ylabel(r"(CLAS6 - VPK) / published total uncertainty")
    ax.set_title(r"CLAS6 VPK closure: $σ_U$, $σ_{TT}$ and $σ_{LT}$")
    ax.grid(alpha=.2); ax.legend()
    fig.tight_layout(); fig.savefig(PNG/"clas6_U_TT_LT_closure_pulls.png",dpi=220); plt.close(fig)

def plot_population_projection(a, ext_c, ext_h, xcol, ycol, xlabel, ylabel, stem, topology="ALL"):
    # Five period panels for one accepted photon topology.
    periods=[p for p in PERIOD_ORDER if p in set(a.period)]
    fig,axes=plt.subplots(2,3,figsize=(14,9),sharex=True,sharey=True)
    axes=axes.ravel()
    for ax,p in zip(axes,periods):
        g=a[(a.period==p)&(a.photon_topology==topology)]
        hb=ax.hexbin(g[xcol],g[ycol],C=g.sum_weight,reduce_C_function=np.sum,
                     gridsize=35,mincnt=1,norm=LogNorm())
        ax.scatter(ext_c[xcol],ext_c[ycol],marker="o",facecolors="white",edgecolors="black",linewidths=.8,s=38,label="CLAS6",zorder=5)
        ax.scatter(ext_h[xcol],ext_h[ycol],marker="s",facecolors="none",edgecolors="red",linewidths=1.2,s=48,label="Hall A",zorder=6)
        ax.set_title(p); ax.grid(alpha=.15)
    for ax in axes[len(periods):]: ax.axis("off")
    for ax in axes[-3:]:
        if ax.axison: ax.set_xlabel(xlabel)
    axes[0].set_ylabel(ylabel); axes[3].set_ylabel(ylabel)
    axes[0].legend(fontsize=8)
    fig.suptitle(f"Accepted RGA AAOgen {topology} population and external π⁰ structure-function coverage",y=.995)
    # Reserve a dedicated strip at far right for the shared colorbar.  Do not
    # let matplotlib place it over the upper-right physics panel.
    fig.subplots_adjust(left=.07,right=.86,bottom=.08,top=.92,wspace=.08,hspace=.16)
    active=[ax for ax in axes if ax.axison]
    if active:
        cax=fig.add_axes([.885,.16,.018,.68])
        cb=fig.colorbar(hb,cax=cax); cb.set_label("Accepted AAOgen weighted population / hex")
    fig.savefig(PNG/f"accepted_aao_{topology.replace(chr(45),chr(95))}_{stem}_external_coverage.png",dpi=220)
    plt.close(fig)


def plot_1d_periods(a):
    specs=[("Q2",r"$Q^2$ (GeV$^2$)","accepted_aao_Q2_by_period.png"),
           ("xB",r"$x_B$","accepted_aao_xB_by_period.png"),
           ("minus_t",r"$-t$ (GeV$^2$)","accepted_aao_minus_t_by_period.png")]
    for col,xlab,name in specs:
        fig,ax=plt.subplots(figsize=(8.6,5.8))
        allg=a[a.photon_topology=="ALL"]
        lo,hi=np.nanpercentile(allg[col],[.2,99.8]); bins=np.linspace(lo,hi,55)
        for p in PERIOD_ORDER:
            g=allg[allg.period==p]
            if g.empty: continue
            hist,edges=np.histogram(g[col],bins=bins,weights=g.sum_weight,density=False)
            if hist.sum()>0: hist=hist/hist.sum()
            ctr=.5*(edges[:-1]+edges[1:])
            ax.step(ctr,hist,where="mid",label=p)
        ax.set_xlabel(xlab); ax.set_ylabel("Fraction of accepted AAOgen weight / bin")
        ax.set_title("Accepted AAOgen kinematics by RGA period")
        ax.grid(alpha=.2); ax.legend()
        fig.tight_layout(); fig.savefig(PNG/name,dpi=220); plt.close(fig)


def plot_topology_coverage(summary):
    s=summary[summary.photon_topology.isin(["FT-FT","FD-FD"])].copy()
    if s.empty: return
    labels=[f"{p}\n{t}" for p,t in zip(s.period,s.photon_topology)]
    x=np.arange(len(s)); w=.75
    fig,ax=plt.subplots(figsize=(12,6))
    a=s.clas6_interpolation_fraction.to_numpy()*100
    b=s.clas6_near_boundary_fraction.to_numpy()*100
    c=s.clas6_outside_fraction.to_numpy()*100
    ax.bar(x,a,w,label="CLAS6 interpolation")
    ax.bar(x,b,w,bottom=a,label="Near boundary")
    ax.bar(x,c,w,bottom=a+b,label="Outside support")
    ax.set_xticks(x,labels,rotation=45,ha="right")
    ax.set_ylabel("Accepted AAOgen weight (%)"); ax.set_ylim(0,100)
    ax.set_title("External-model support for photon-pair topology samples")
    ax.legend(); ax.grid(axis="y",alpha=.2)
    fig.tight_layout(); fig.savefig(PNG/"clas6_support_fraction_by_topology.png",dpi=220); plt.close(fig)


def plot_folded_M(fold):
    s=fold[fold.photon_topology.isin(["ALL","FT-FT","FD-FD"])].copy()
    if s.empty:return
    labels=[f"{p}\n{t}" for p,t in zip(s.period,s.photon_topology)]
    x=np.arange(len(s)); width=.38
    fig,ax=plt.subplots(figsize=(13,6))
    ax.bar(x-width/2,s.clas6_M_acc_supported,width,label="CLAS6-supported")
    ax.bar(x+width/2,s.clas6_halla_M_acc_supported,width,label="CLAS6 + Hall A supported")
    ax.axhline(1,linewidth=1)
    ax.set_xticks(x,labels,rotation=45,ha="right"); ax.set_ylabel(r"Folded $M_{acc}$ = data / VPK")
    ax.set_title("External π⁰ model residual folded through accepted AAOgen population")
    ax.legend(); ax.grid(axis="y",alpha=.2)
    fig.tight_layout(); fig.savefig(PNG/"folded_external_model_correction_by_period_topology.png",dpi=220); plt.close(fig)


def plot_external_distance(a: pd.DataFrame, topology="FT-FT"):
    periods=[p for p in PERIOD_ORDER if p in set(a.period)]
    fig,axes=plt.subplots(2,3,figsize=(14,9),sharex=True,sharey=True); axes=axes.ravel()
    vmax=float(np.nanpercentile(a.loc[a.photon_topology==topology,"external_nearest_distance_normalized"],99))
    sc=None
    for ax,p in zip(axes,periods):
        g=a[(a.period==p)&(a.photon_topology==topology)]
        # Cell marker size follows accepted weight weakly; color is the actual support distance.
        size=8+30*np.sqrt(g.sum_weight/np.nanmax(g.sum_weight))
        sc=ax.scatter(g.xB,g.Q2,c=g.external_nearest_distance_normalized,s=size,
                      vmin=0,vmax=vmax,alpha=.75,rasterized=True)
        ax.set_title(p); ax.grid(alpha=.15)
    for ax in axes[len(periods):]: ax.axis("off")
    for ax in axes[-3:]:
        if ax.axison: ax.set_xlabel(r"$x_B$")
    axes[0].set_ylabel(r"$Q^2$ (GeV$^2$)"); axes[3].set_ylabel(r"$Q^2$ (GeV$^2$)")
    fig.suptitle(f"Accepted {topology} AAOgen: continuous distance to nearest external π⁰ measurement",y=.995)
    fig.subplots_adjust(left=.07,right=.86,bottom=.08,top=.92,wspace=.08,hspace=.16)
    if sc is not None:
        cax=fig.add_axes([.885,.16,.018,.68]); cb=fig.colorbar(sc,cax=cax)
        cb.set_label("Nearest-data distance (normalized 3D kinematics)")
    fig.savefig(PNG/f"accepted_aao_{topology.replace('-','_')}_external_distance.png",dpi=220); plt.close(fig)


def plot_local_correction(a: pd.DataFrame, topology="FT-FT"):
    periods=[p for p in PERIOD_ORDER if p in set(a.period)]
    fig,axes=plt.subplots(2,3,figsize=(14,9),sharex=True,sharey=True); axes=axes.ravel()
    vals=a.loc[a.photon_topology==topology,"external_local_M"]
    lo,hi=np.nanpercentile(vals,[1,99]); d=max(abs(lo-1),abs(hi-1)); lo,hi=1-d,1+d
    sc=None
    for ax,p in zip(axes,periods):
        g=a[(a.period==p)&(a.photon_topology==topology)]
        size=8+30*np.sqrt(g.sum_weight/np.nanmax(g.sum_weight))
        sc=ax.scatter(g.xB,g.Q2,c=g.external_local_M,s=size,vmin=lo,vmax=hi,alpha=.75,rasterized=True)
        ax.set_title(p); ax.grid(alpha=.15)
    for ax in axes[len(periods):]: ax.axis("off")
    for ax in axes[-3:]:
        if ax.axison: ax.set_xlabel(r"$x_B$")
    axes[0].set_ylabel(r"$Q^2$ (GeV$^2$)"); axes[3].set_ylabel(r"$Q^2$ (GeV$^2$)")
    fig.suptitle(f"Accepted {topology} AAOgen: smooth external data/VPK correction field",y=.995)
    fig.subplots_adjust(left=.07,right=.86,bottom=.08,top=.92,wspace=.08,hspace=.16)
    if sc is not None:
        cax=fig.add_axes([.885,.16,.018,.68]); cb=fig.colorbar(sc,cax=cax); cb.set_label(r"Local $M=$ data/VPK")
    fig.savefig(PNG/f"accepted_aao_{topology.replace('-','_')}_local_model_correction.png",dpi=220); plt.close(fig)


def plot_robustness(r: pd.DataFrame):
    if r.empty:return
    s=r[r.photon_topology=="FT-FT"].copy()
    x=np.arange(len(s)); fig,ax=plt.subplots(figsize=(9.5,5.8))
    ax.errorbar(x,s.local_M_acc,yerr=s.local_M_acc_unc,fmt="o",capsize=4,label="Smooth external constraint")
    ax.scatter(x,s.M_required_for_unit_photon_efficiency,marker="s",s=55,label=r"Required if $r_{FT}=1$")
    ax.axhline(1,linewidth=1,linestyle="--")
    ax.set_xticks(x,s.period,rotation=25,ha="right"); ax.set_ylabel(r"Accepted-population model ratio $M=$ data/VPK")
    ax.set_title("FT robustness: external π⁰ model constraint vs correction required to remove photon inefficiency")
    ax.grid(axis="y",alpha=.2); ax.legend(); fig.tight_layout()
    fig.savefig(PNG/"ft_model_correction_robustness.png",dpi=220); plt.close(fig)


def plot_distance_cdf(a: pd.DataFrame):
    fig,ax=plt.subplots(figsize=(9,6))
    for topology in ["FD-FD","FT-FT"]:
        g=a[a.photon_topology==topology]
        d=g.external_nearest_distance_normalized.to_numpy(float); w=g.sum_weight.to_numpy(float)
        o=np.argsort(d); d=d[o]; w=w[o]; c=np.cumsum(w)/np.sum(w)
        ax.plot(d,c,label=topology)
    ax.set_xlabel("Nearest external-data distance (normalized 3D kinematics)"); ax.set_ylabel("Cumulative accepted AAOgen weight")
    ax.set_ylim(0,1.01); ax.grid(alpha=.2); ax.legend(); ax.set_title("How far are accepted photon topologies from measured π⁰ kinematics?")
    fig.tight_layout(); fig.savefig(PNG/"external_distance_cdf_by_topology.png",dpi=220); plt.close(fig)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--accepted-aao",type=Path,default=DEFAULT_AAO)
    ap.add_argument("--topology-summary",type=Path,default=DEFAULT_TOPOLOGY)
    args=ap.parse_args()

    OUT.mkdir(parents=True,exist_ok=True); PNG.mkdir(parents=True,exist_ok=True)
    ensure_validation_outputs()
    clas6,halla=load_external_u_points()
    a=load_aao(args.accepted_aao)
    a=add_clas6_coverage(a,clas6)
    a=add_model_interpolations(a,clas6,halla)
    a=add_local_external_field(a,clas6,halla)

    cov=coverage_summary(a); fold=fold_summary(a); local_fold=local_fold_summary(a)
    cov.to_csv(OUT/"rga_external_pi0_coverage_summary.csv",index=False)
    interference=interference_closure(CLAS6_POINTS)
    interference.to_csv(OUT/"clas6_U_TT_LT_closure_summary.csv",index=False)
    fold.to_csv(OUT/"rga_vpk_fold_summary.csv",index=False)
    local_fold.to_csv(OUT/"rga_local_external_model_fold.csv",index=False)
    a.to_csv(OUT/"accepted_aao_population_with_external_support.csv",index=False)

    cand=candidate_efficiencies(fold,args.topology_summary)
    if not cand.empty:
        cand.to_csv(OUT/"rga_candidate_photon_efficiency.csv",index=False)
    robust=robustness_efficiencies(local_fold,args.topology_summary)
    if not robust.empty:
        robust.to_csv(OUT/"rga_photon_efficiency_robustness.csv",index=False)

    for topology in ["ALL","FD-FD","FT-FT"]:
        plot_population_projection(a,clas6,halla,"xB","Q2",r"$x_B$",r"$Q^2$ (GeV$^2$)","xB_Q2",topology)
        plot_population_projection(a,clas6,halla,"minus_t","Q2",r"$-t$ (GeV$^2$)",r"$Q^2$ (GeV$^2$)","minus_t_Q2",topology)
    plot_interference_pulls(CLAS6_POINTS)
    plot_1d_periods(a)
    plot_topology_coverage(cov)
    plot_folded_M(fold)
    for topology in ["FD-FD","FT-FT"]:
        plot_external_distance(a,topology)
        plot_local_correction(a,topology)
    plot_distance_cdf(a)
    plot_robustness(robust)

    print("\n=== RGA accepted-AAOgen external pi0 coverage ===")
    show=cov[cov.photon_topology=="ALL"][["period","raw_events","Q2_median_GeV2","clas6_interpolation_fraction","clas6_near_boundary_fraction","clas6_outside_fraction"]]
    print(show.to_string(index=False,formatters={c:"{:.3f}".format for c in show.columns if c not in ["period","raw_events"]}))
    print("\n=== CLAS6 U/TT/LT closure (TT/LT ratios can be unstable near zero) ===")
    print(interference.to_string(index=False))
    print("\n=== Supported-yield VPK residual fold (ALL / FT-FT / FD-FD) ===")
    print(fold[fold.photon_topology.isin(["ALL","FT-FT","FD-FD"])].to_string(index=False,formatters={
        "clas6_supported_weight_fraction":"{:.3f}".format,"clas6_M_acc_supported":"{:.3f}".format,
        "clas6_halla_supported_weight_fraction":"{:.3f}".format,"clas6_halla_M_acc_supported":"{:.3f}".format}))
    if not robust.empty:
        print("\n=== Smooth-field robustness test (all accepted kinematics) ===")
        print(robust.to_string(index=False))
    if not cand.empty:
        print("\n=== Candidate single-photon residuals (diagnostic; support fraction shown) ===")
        print(cand.to_string(index=False))
    print(f"\nWrote tables under: {OUT}")
    print(f"Wrote PNG figures under: {PNG}")

if __name__=="__main__":
    main()
