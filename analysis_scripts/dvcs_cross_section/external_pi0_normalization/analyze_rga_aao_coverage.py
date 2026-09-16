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
# The transport kernel is no longer fixed by hand.  Candidate bandwidths and
# neighbor counts are tested by withholding entire published kinematic settings.
# The selected values are those minimizing the setting-level closure RMSE.
TRANSPORT_BANDWIDTH_SCAN = [0.08, 0.12, 0.16, 0.20, 0.25, 0.32, 0.40, 0.55]
TRANSPORT_K_SCAN = [4, 6, 8, 12, 16, 24, 32, 48]
MIN_DISTANCE_BIN_COUNT = 8
UNCERTAINTY_K_SCAN = [2, 3, 4, 5, 6, 8, 10, 12]


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
    # New eppi0 exporter columns.  epsilon_data/epsilon_MC is the proton
    # efficiency ratio itself; the production DATA correction is its inverse.
    if "proton_efficiency_ratio_mean" not in a.columns:
        a["proton_efficiency_ratio_mean"] = np.nan
    if "proton_data_weight_mean" not in a.columns:
        a["proton_data_weight_mean"] = np.nan
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


def _external_scale(ext: pd.DataFrame):
    """Return a transparent dimensionless coordinate scaling.

    The scale itself is not used to claim an uncertainty: it only defines a
    distance coordinate.  All transport accuracy versus that distance is then
    measured by withheld-data closure below.
    """
    xyz=ext[["xB","Q2","minus_t"]].to_numpy(float)
    lo=xyz.min(axis=0); span=xyz.max(axis=0)-lo
    span=np.where(span>0,span,1.0)
    return lo,span


def _predict_kernel(train: pd.DataFrame, query_xyz: np.ndarray, lo: np.ndarray,
                    span: np.ndarray, bandwidth: float, k: int):
    """Kernel prediction with no uncertainty floor and no imposed distance error."""
    xyz=train[["xB","Q2","minus_t"]].to_numpy(float)
    xn=(xyz-lo)/span; qn=(np.asarray(query_xyz,float)-lo)/span
    tree=cKDTree(xn)
    kk=min(int(k),len(train))
    dist,idx=tree.query(qn,k=kk)
    if kk==1:
        dist=dist[:,None]; idx=idx[:,None]
    vals=train.ratio.to_numpy(float)[idx]
    errs=train.ratio_err.to_numpy(float)[idx]
    # Published errors are used as published.  Invalid/missing errors simply
    # remove inverse-variance preference rather than introducing an arbitrary floor.
    good=np.isfinite(errs)&(errs>0)
    finite=train.ratio_err.to_numpy(float)
    finite=finite[np.isfinite(finite)&(finite>0)]
    fallback=float(np.median(finite)) if len(finite) else 1.0
    ee=np.where(good,errs,fallback)
    kw=np.exp(-0.5*(dist/float(bandwidth))**2)/(ee**2)
    sw=np.sum(kw,axis=1)
    pred=np.sum(kw*vals,axis=1)/np.maximum(sw,1e-300)
    return pred,dist[:,0]


def _setting_ids(ext: pd.DataFrame) -> pd.Series:
    # CLAS6 source_table is one published (Q2,xB) setting.  Hall A has four t'
    # points per (Q2,xB) setting, so construct the corresponding group explicitly.
    out=[]
    for _,r in ext.iterrows():
        if str(r["experiment"]).startswith("CLAS6"):
            out.append("CLAS6:"+str(r["source_table"]))
        else:
            out.append(f"HallA:Q2={float(r.Q2):.3f}:xB={float(r.xB):.3f}")
    return pd.Series(out,index=ext.index,dtype=str)


def _cross_validate_transport(ext: pd.DataFrame, bandwidth: float, k: int,
                              mode: str) -> pd.DataFrame:
    """Withhold points or whole kinematic settings and predict them."""
    lo,span=_external_scale(ext)
    rows=[]
    if mode=="point":
        groups=[[i] for i in ext.index]
    elif mode=="setting":
        sid=_setting_ids(ext)
        groups=[list(ix) for _,ix in sid.groupby(sid).groups.items()]
    else:
        raise ValueError(mode)
    for hold in groups:
        test=ext.loc[hold]
        train=ext.drop(index=hold)
        if len(train)<4: continue
        pred,dnear=_predict_kernel(train,test[["xB","Q2","minus_t"]].to_numpy(float),lo,span,bandwidth,k)
        for pos,(idx,r) in enumerate(test.iterrows()):
            rows.append({"cv_mode":mode,"held_index":int(idx),"experiment":r.experiment,
                         "setting_id":_setting_ids(test).loc[idx],"xB":r.xB,"Q2":r.Q2,
                         "minus_t":r.minus_t,"measured_M":r.ratio,"measured_M_err":r.ratio_err,
                         "predicted_M":pred[pos],"residual":pred[pos]-r.ratio,
                         "abs_residual":abs(pred[pos]-r.ratio),"nearest_training_distance":dnear[pos],
                         "bandwidth":bandwidth,"k":k})
    return pd.DataFrame(rows)


def calibrate_external_transport(clas6: pd.DataFrame, halla: pd.DataFrame):
    """Choose transport settings from withheld-setting closure, not from RGA yields.

    Hyperparameters are selected using CLAS6 leave-one-setting-out closure only.
    Hall A is retained as a genuinely independent transport test because VPK was
    fitted to the CLAS6 data.  Point-LOO is reported as a less stringent diagnostic.
    """
    ext=pd.concat([clas6,halla],ignore_index=True)
    scans=[]
    # Tune only on CLAS6 whole-setting holdout.  This avoids choosing parameters
    # because they happen to make the RGA or Hall-A result look favorable.
    for bw in TRANSPORT_BANDWIDTH_SCAN:
        for k in TRANSPORT_K_SCAN:
            cv=_cross_validate_transport(clas6,bw,k,"setting")
            if cv.empty: continue
            rmse=float(np.sqrt(np.mean(cv.residual.to_numpy(float)**2)))
            mae=float(np.mean(cv.abs_residual))
            bias=float(np.mean(cv.residual))
            scans.append({"bandwidth":bw,"k":k,"clas6_setting_rmse":rmse,
                          "clas6_setting_mae":mae,"clas6_setting_bias":bias,"n":len(cv)})
    scan=pd.DataFrame(scans).sort_values(["clas6_setting_rmse","clas6_setting_mae","k","bandwidth"]).reset_index(drop=True)
    if scan.empty: raise RuntimeError("External transport scan produced no valid closure tests")
    best=scan.iloc[0]
    bw=float(best.bandwidth); k=int(best.k)

    # Final closure products.  Combined-data setting holdout measures local
    # interpolation behavior; CLAS6->HallA is the independent extrapolation test.
    point=_cross_validate_transport(ext,bw,k,"point")
    setting=_cross_validate_transport(ext,bw,k,"setting")
    lo,span=_external_scale(ext)
    hall_pred,hall_dist=_predict_kernel(clas6,halla[["xB","Q2","minus_t"]].to_numpy(float),lo,span,bw,k)
    independent=halla.copy()
    independent["cv_mode"]="CLAS6_to_HallA"
    independent["held_index"]=np.arange(len(independent))
    independent["setting_id"]=_setting_ids(independent).values
    independent["measured_M"]=independent.ratio
    independent["measured_M_err"]=independent.ratio_err
    independent["predicted_M"]=hall_pred
    independent["residual"]=hall_pred-independent.ratio.to_numpy(float)
    independent["abs_residual"]=np.abs(independent.residual)
    independent["nearest_training_distance"]=hall_dist
    independent["bandwidth"]=bw; independent["k"]=k
    independent=independent[[c for c in point.columns if c in independent.columns]]
    closure=pd.concat([point,setting,independent],ignore_index=True,sort=False)

    # Empirical transport uncertainty versus distance comes from the stringent
    # whole-setting holdout.  Use expanding distance thresholds so every estimate
    # is based on enough actually-withheld measurements.  RMS residual is the
    # directly observed prediction error; no a*d term is imposed.
    base=setting.copy().sort_values("nearest_training_distance")
    dvals=np.unique(np.quantile(base.nearest_training_distance,[0,.15,.30,.45,.60,.75,.90,1]))
    erows=[]
    for dmax in dvals[1:]:
        g=base[base.nearest_training_distance<=dmax]
        if len(g)<MIN_DISTANCE_BIN_COUNT: continue
        erows.append({"distance_max":float(dmax),"n":len(g),
                      "bias":float(np.mean(g.residual)),
                      "rmse":float(np.sqrt(np.mean(g.residual**2))),
                      "mae":float(np.mean(g.abs_residual)),
                      "p68_abs_residual":float(np.quantile(g.abs_residual,.68)),
                      "p95_abs_residual":float(np.quantile(g.abs_residual,.95))})
    empirical=pd.DataFrame(erows)
    return scan,closure,empirical,bw,k


def setting_closure_summary(closure: pd.DataFrame) -> pd.DataFrame:
    """Collapse whole-setting holdout residuals to one empirical observation per setting.

    The RMS here is the error actually made when an entire published kinematic
    setting was absent from the training sample.  These setting-level errors,
    rather than individual t points, calibrate the phase-space uncertainty.
    """
    g=closure[closure.cv_mode=="setting"].copy()
    rows=[]
    for sid,q in g.groupby("setting_id",sort=False):
        rows.append({
            "setting_id":sid,
            "experiment":str(q.experiment.iloc[0]),
            "xB_mean":float(q.xB.mean()),
            "Q2_mean_GeV2":float(q.Q2.mean()),
            "minus_t_mean_GeV2":float(q.minus_t.mean()),
            "minus_t_min_GeV2":float(q.minus_t.min()),
            "minus_t_max_GeV2":float(q.minus_t.max()),
            "n_points":int(len(q)),
            "mean_signed_residual":float(q.residual.mean()),
            "rmse_residual":float(np.sqrt(np.mean(q.residual.to_numpy(float)**2))),
            "mae_residual":float(q.abs_residual.mean()),
            "p68_abs_residual":float(np.quantile(q.abs_residual,.68)),
            "nearest_training_distance_median":float(q.nearest_training_distance.median()),
            "nearest_training_distance_max":float(q.nearest_training_distance.max()),
        })
    return pd.DataFrame(rows)


def calibrate_phase_space_uncertainty(setting_summary: pd.DataFrame,
                                      ext: pd.DataFrame):
    """Calibrate a local uncertainty from held-out *setting* failures.

    No analytic Q2/xB/t dependence is assumed.  For each held-out setting, its
    uncertainty is predicted from the RMS failures of the k nearest OTHER
    held-out settings in normalized (xB,Q2,-t) space.  k is selected by a
    second leave-one-setting-out layer using Gaussian log score.  This asks a
    direct question: did the uncertainty inferred from neighboring closure
    failures assign sensible probability to the unseen setting's actual error?

    Beyond the phase-space separation exercised by this second withholding
    layer, the method falls back to the global setting-holdout RMS.
    """
    if setting_summary.empty:
        return pd.DataFrame(),None,np.nan,np.nan
    lo,span=_external_scale(ext)
    X=setting_summary[["xB_mean","Q2_mean_GeV2","minus_t_mean_GeV2"]].to_numpy(float)
    Xn=(X-lo)/span
    e=setting_summary.rmse_residual.to_numpy(float)
    global_sigma=float(np.sqrt(np.mean(e**2)))
    rows=[]
    for k0 in UNCERTAINTY_K_SCAN:
        k=min(int(k0),len(setting_summary)-1)
        if k<1: continue
        pred=[]; dnearest=[]
        for i in range(len(setting_summary)):
            mask=np.arange(len(setting_summary))!=i
            dist=np.linalg.norm(Xn[mask]-Xn[i],axis=1)
            order=np.argsort(dist)[:k]
            # Equal weights: no additional bandwidth or distance-decay constant.
            sig=float(np.sqrt(np.mean(e[mask][order]**2)))
            pred.append(sig); dnearest.append(float(dist[order[0]]))
        pred=np.asarray(pred); dnearest=np.asarray(dnearest)
        # Proper Gaussian log score (constant omitted). Lower is better.
        score=float(np.mean(np.log(pred)+0.5*(e/pred)**2))
        z=e/pred
        rows.append({"k_uncertainty_neighbors":k0,"effective_k":k,
                     "nested_log_score":score,
                     "standardized_error_rms":float(np.sqrt(np.mean(z**2))),
                     "fraction_error_below_1sigma":float(np.mean(e<=pred)),
                     "global_setting_rmse":global_sigma,
                     "max_tested_nearest_setting_distance":float(np.max(dnearest))})
    scan=pd.DataFrame(rows).sort_values(["nested_log_score","effective_k"]).reset_index(drop=True)
    best=int(scan.iloc[0].effective_k)
    maxdist=float(scan.iloc[0].max_tested_nearest_setting_distance)
    return scan,best,global_sigma,maxdist


def _phase_space_sigma(query_xyz: np.ndarray, setting_summary: pd.DataFrame,
                       ext: pd.DataFrame, k: int, global_sigma: float,
                       max_valid_distance: float):
    lo,span=_external_scale(ext)
    S=setting_summary[["xB_mean","Q2_mean_GeV2","minus_t_mean_GeV2"]].to_numpy(float)
    Sn=(S-lo)/span; Qn=(np.asarray(query_xyz,float)-lo)/span
    tree=cKDTree(Sn)
    kk=min(int(k),len(Sn))
    dist,idx=tree.query(Qn,k=kk)
    if kk==1:
        dist=np.asarray(dist)[:,None]; idx=np.asarray(idx)[:,None]
    err=setting_summary.rmse_residual.to_numpy(float)
    local=np.sqrt(np.mean(err[idx]**2,axis=1))
    nearest=dist[:,0]
    fallback=nearest>max_valid_distance
    sigma=local.copy(); sigma[fallback]=global_sigma
    return sigma,nearest,fallback


def add_local_external_field(a: pd.DataFrame, clas6: pd.DataFrame, halla: pd.DataFrame,
                             bandwidth: float, k: int,
                             setting_summary: pd.DataFrame,
                             uncertainty_k: int, global_sigma: float,
                             max_valid_uncertainty_distance: float) -> pd.DataFrame:
    """Apply central M transport plus independently calibrated local uncertainty."""
    out=a.copy(); ext=pd.concat([clas6,halla],ignore_index=True)
    lo,span=_external_scale(ext)
    xyz=out[["xB","Q2","minus_t"]].to_numpy(float)
    pred,dnear=_predict_kernel(ext,xyz,lo,span,bandwidth,k)
    sigma,dclosure,fallback=_phase_space_sigma(
        xyz,setting_summary,ext,uncertainty_k,global_sigma,max_valid_uncertainty_distance)
    out["external_nearest_distance_normalized"]=dnear
    out["external_local_M"]=pred
    out["external_local_M_unc"]=sigma
    out["closure_setting_nearest_distance_normalized"]=dclosure
    out["external_uncertainty_used_global_fallback"]=fallback
    out["external_local_n_eff"]=np.nan
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
                     "nearest_distance_p84":qd[2],"nearest_distance_p95":qd[3],
                     "global_uncertainty_fallback_weight_fraction":float(np.sum(w*g.external_uncertainty_used_global_fallback.to_numpy(bool))/W)})
    return pd.DataFrame(rows)


def proton_fold_summary(a: pd.DataFrame) -> pd.DataFrame:
    """Fold Krishna Neupane epsilon_data/epsilon_MC through each accepted topology.

    The C++ exporter evaluates exactly the production parameterization on each
    accepted reconstructed-AAOgen proton.  No Krishna weight is applied to the
    eppi0 yield itself; this table only decomposes the observed DATA/MC ratio.
    """
    rows=[]
    for (p,t),g in a.groupby(["period","photon_topology"],sort=False):
        w=g.sum_weight.to_numpy(float); W=float(w.sum())
        pr=pd.to_numeric(g.proton_efficiency_ratio_mean,errors="coerce").to_numpy(float)
        pw=pd.to_numeric(g.proton_data_weight_mean,errors="coerce").to_numpy(float)
        good=np.isfinite(pr)&(pr>0)&np.isfinite(w)&(w>0)
        if not good.any():
            mean_pr=mean_pw=np.nan
        else:
            mean_pr=float(np.sum(w[good]*pr[good])/np.sum(w[good]))
            goodw=np.isfinite(pw)&good
            mean_pw=float(np.sum(w[goodw]*pw[goodw])/np.sum(w[goodw])) if goodw.any() else np.nan
        rows.append({"period":p,"photon_topology":t,"sum_weight":W,
                     "krishna_epsilon_data_over_mc":mean_pr,
                     "krishna_production_data_weight":mean_pw})
    return pd.DataFrame(rows)


def topology_model_ratio_jackknife(a: pd.DataFrame, clas6: pd.DataFrame, halla: pd.DataFrame,
                                   bandwidth: float, k: int) -> tuple[pd.DataFrame,pd.DataFrame]:
    """Propagate common external-model variations jointly through FT-FT and FD-FD.

    Each published external setting is removed in turn, the central M field is
    rebuilt from all remaining external measurements, and that SAME varied field
    is folded through both photon topologies.  The spread of M_FT/M_FD therefore
    measures the part of external-model transport sensitivity that survives the
    common-mode cancellation.  This is intentionally reported separately from
    the absolute phase-space transport uncertainty; it is not obtained by
    pretending the two ~15% absolute uncertainties are independent.
    """
    ext=pd.concat([clas6,halla],ignore_index=True).copy()
    ext["setting_id"]=_setting_ids(ext).values
    lo,span=_external_scale(ext)
    xyz=a[["xB","Q2","minus_t"]].to_numpy(float)
    nominal,_=_predict_kernel(ext,xyz,lo,span,bandwidth,k)
    base=a[["period","photon_topology","sum_weight"]].copy(); base["M"]=nominal
    def folds(frame):
        rr=[]
        for (p,t),g in frame.groupby(["period","photon_topology"],sort=False):
            W=float(g.sum_weight.sum())
            rr.append((p,t,float(np.sum(g.sum_weight*g.M)/W)))
        return pd.DataFrame(rr,columns=["period","photon_topology","M"])
    nom=folds(base)
    nomp=nom[nom.photon_topology.isin(["FT-FT","FD-FD"])].pivot(index="period",columns="photon_topology",values="M")
    rows=[]
    for sid in ext.setting_id.drop_duplicates():
        train=ext[ext.setting_id!=sid]
        pred,_=_predict_kernel(train,xyz,lo,span,bandwidth,k)
        q=a[["period","photon_topology","sum_weight"]].copy(); q["M"]=pred
        f=folds(q); piv=f[f.photon_topology.isin(["FT-FT","FD-FD"])].pivot(index="period",columns="photon_topology",values="M")
        for period in piv.index:
            if not {"FT-FT","FD-FD"}.issubset(piv.columns): continue
            ratio=float(piv.loc[period,"FT-FT"]/piv.loc[period,"FD-FD"])
            rows.append({"held_setting":sid,"period":period,
                         "M_FT_FT":float(piv.loc[period,"FT-FT"]),
                         "M_FD_FD":float(piv.loc[period,"FD-FD"]),
                         "M_FT_over_FD":ratio})
    variants=pd.DataFrame(rows)
    summary=[]
    for period,g in variants.groupby("period",sort=False):
        r0=float(nomp.loc[period,"FT-FT"]/nomp.loc[period,"FD-FD"])
        d=g.M_FT_over_FD.to_numpy(float)-r0
        summary.append({"period":period,
                        "nominal_M_FT_FT":float(nomp.loc[period,"FT-FT"]),
                        "nominal_M_FD_FD":float(nomp.loc[period,"FD-FD"]),
                        "nominal_M_FT_over_FD":r0,
                        "n_setting_omissions":int(len(g)),
                        "loso_ratio_rms":float(np.sqrt(np.mean(d*d))),
                        "loso_ratio_max_abs_shift":float(np.max(np.abs(d))),
                        "loso_ratio_p68_abs_shift":float(np.quantile(np.abs(d),.68))})
    return pd.DataFrame(summary),variants


def relative_photon_efficiency(topology_path: Path, model_ratio: pd.DataFrame,
                               proton_fold: pd.DataFrame) -> pd.DataFrame:
    if not topology_path.exists() or model_ratio.empty: return pd.DataFrame()
    top=pd.read_csv(topology_path)
    top=top[(top.comparison=="raw") & top.photon_topology.isin(["FT-FT","FD-FD"])].copy()
    r=top.pivot(index="period",columns="photon_topology",values="data_over_mc")
    e=top.pivot(index="period",columns="photon_topology",values="stat_err")
    pp=proton_fold[proton_fold.photon_topology.isin(["FT-FT","FD-FD"])].pivot(index="period",columns="photon_topology",values="krishna_epsilon_data_over_mc")
    rows=[]
    for _,m in model_ratio.iterrows():
        p=m.period
        if p not in r.index or p not in pp.index: continue
        Rft=float(r.loc[p,"FT-FT"]); Rfd=float(r.loc[p,"FD-FD"])
        Pft=float(pp.loc[p,"FT-FT"]); Pfd=float(pp.loc[p,"FD-FD"])
        Mr=float(m.nominal_M_FT_over_FD)
        rel=math.sqrt((Rft/Rfd)/(Mr*(Pft/Pfd)))
        # Statistical uncertainty from the two raw topology ratios only.
        sft=float(e.loc[p,"FT-FT"]); sfd=float(e.loc[p,"FD-FD"])
        stat=0.5*rel*math.sqrt((sft/Rft)**2+(sfd/Rfd)**2)
        # Convert the empirically observed LOSO M-ratio sensitivity to rFT/rFD.
        model=0.5*rel*float(m.loso_ratio_rms)/Mr
        rows.append({"period":p,"R_FT_over_FD":Rft/Rfd,
                     "M_FT_over_FD":Mr,"M_FT_over_FD_loso_rms":float(m.loso_ratio_rms),
                     "proton_eff_FT_over_FD":Pft/Pfd,
                     "r_FT_over_r_FD":rel,"stat_err":stat,
                     "external_model_ratio_sensitivity":model})
    return pd.DataFrame(rows)


def absolute_photon_with_proton(local_fold: pd.DataFrame, topology_path: Path,
                                proton_fold: pd.DataFrame) -> pd.DataFrame:
    if not topology_path.exists(): return pd.DataFrame()
    top=pd.read_csv(topology_path)
    top=top[(top.comparison=="raw") & top.photon_topology.isin(["FT-FT","FD-FD"])].copy()
    z=top.merge(local_fold,on=["period","photon_topology"],how="left").merge(
        proton_fold,on=["period","photon_topology"],how="left",suffixes=("","_p"))
    P=z.krishna_epsilon_data_over_mc
    z["candidate_single_photon_ratio_after_krishna"]=np.sqrt(z.data_over_mc/(z.local_M_acc*P))
    z["candidate_single_photon_ratio_before_krishna"]=np.sqrt(z.data_over_mc/z.local_M_acc)
    return z[["period","photon_topology","data_over_mc","stat_err","local_M_acc","local_M_acc_unc",
              "krishna_epsilon_data_over_mc","krishna_production_data_weight",
              "candidate_single_photon_ratio_before_krishna","candidate_single_photon_ratio_after_krishna"]]


def plot_topology_model_cancellation(model_ratio: pd.DataFrame) -> None:
    if model_ratio.empty: return
    q=model_ratio.copy(); q["period"]=pd.Categorical(q.period,PERIOD_ORDER,ordered=True); q=q.sort_values("period")
    x=np.arange(len(q))
    fig,ax=plt.subplots(figsize=(9.2,5.6))
    y=q.nominal_M_FT_over_FD.to_numpy(float); e=q.loso_ratio_rms.to_numpy(float)
    ax.errorbar(x,y,yerr=e,fmt="o",capsize=4,label="M(FT-FT) / M(FD-FD)")
    ax.axhline(1.0,ls="--",lw=1)
    ax.set_xticks(x,q.period.astype(str),rotation=25,ha="right")
    ax.set_ylabel("External-model topology ratio")
    ax.set_title("Common VPK transport cancellation between photon topologies")
    ax.legend(frameon=False); fig.tight_layout()
    fig.savefig(PNG/"topology_model_ratio_common_mode_cancellation.png",dpi=220); plt.close(fig)


def plot_krishna_and_relative(proton_fold: pd.DataFrame, relative: pd.DataFrame) -> None:
    if not proton_fold.empty:
        q=proton_fold[proton_fold.photon_topology.isin(["FT-FT","FD-FD"])].copy()
        q["period"]=pd.Categorical(q.period,PERIOD_ORDER,ordered=True); q=q.sort_values(["period","photon_topology"])
        fig,ax=plt.subplots(figsize=(9.4,5.7))
        for topo,dx in [("FD-FD",-0.10),("FT-FT",0.10)]:
            g=q[q.photon_topology==topo]; xx=np.arange(len(g))+dx
            ax.plot(xx,g.krishna_epsilon_data_over_mc,"o",label=topo)
        ax.axhline(1.0,ls="--",lw=1); ax.set_xticks(np.arange(len(PERIOD_ORDER)),PERIOD_ORDER,rotation=25,ha="right")
        ax.set_ylabel("Folded proton efficiency ratio, data / MC")
        ax.set_title("Krishna Neupane proton efficiency in accepted eπ⁰ topologies")
        ax.legend(frameon=False); fig.tight_layout()
        fig.savefig(PNG/"krishna_proton_efficiency_by_photon_topology.png",dpi=220); plt.close(fig)
    if not relative.empty:
        q=relative.copy(); q["period"]=pd.Categorical(q.period,PERIOD_ORDER,ordered=True); q=q.sort_values("period")
        x=np.arange(len(q)); y=q.r_FT_over_r_FD.to_numpy(float)
        estat=q.stat_err.to_numpy(float); emod=q.external_model_ratio_sensitivity.to_numpy(float)
        fig,ax=plt.subplots(figsize=(9.4,5.7))
        ax.errorbar(x,y,yerr=np.sqrt(estat**2+emod**2),fmt="o",capsize=4,label="stat ⊕ external-model ratio sensitivity")
        ax.errorbar(x,y,yerr=estat,fmt="none",capsize=3,label="statistical only")
        ax.axhline(1.0,ls="--",lw=1); ax.set_xticks(x,q.period.astype(str),rotation=25,ha="right")
        ax.set_ylabel("Relative photon efficiency, FT / FD")
        ax.set_title("eπ⁰ relative FT-to-FD photon efficiency after proton correction")
        ax.legend(frameon=False); fig.tight_layout()
        fig.savefig(PNG/"relative_ft_fd_photon_efficiency.png",dpi=220); plt.close(fig)



def krishna_uncertainty_summary(proton_fold: pd.DataFrame) -> pd.DataFrame:
    """Attach the *existing production* Krishna Neupane normalization uncertainty.

    The DVCS production analysis currently treats Krishna's proton-efficiency
    uncertainty as an overall normalization source: 2.83% for Fa18 and Sp19,
    and a deliberately conservative 11.32% transfer uncertainty for Sp18.

    Because that published/production uncertainty is an overall scale on the
    proton correction, it is 100% common between FT-FT and FD-FD populations
    of a given period and therefore cancels in the FT/FD *relative* ratio.
    This does NOT claim that the differential p/theta/phi shape of Krishna's
    map is known perfectly.  No covariance/variation set for that map is
    available in this repository, so a topology-dependent Krishna shape
    uncertainty cannot honestly be manufactured here.
    """
    if proton_fold.empty:
        return pd.DataFrame()
    def frac(period: str) -> float:
        return 0.1132 if period in ("Sp18 Inb", "Sp18 Out") else 0.0283
    q=proton_fold.copy()
    q["krishna_overall_norm_sys_fraction"]=[frac(str(x)) for x in q.period]
    q["krishna_relative_ft_fd_sys_fraction_from_overall_scale"]=0.0
    q["krishna_topology_shape_sys_available"]=False
    return q


def pass2_normalization_summary(local_fold: pd.DataFrame, topology_path: Path,
                                proton_fold: pd.DataFrame) -> pd.DataFrame:
    """Build the normalization-focused result that motivated this study.

    For same-detector pi0 topologies,

      R_pi0 = M_VPK * P_proton * r_gamma^2,

    so the inferred single-photon data/MC efficiency is

      r_gamma = sqrt(R_pi0 / (M_VPK P_proton)),

    and the candidate one-photon correction for DVCS is C_gamma=1/r_gamma.

    The table propagates only uncertainty components that are currently
    quantified and defensible in this repository: pi0 counting statistics,
    the empirically calibrated absolute VPK transport uncertainty, and the
    existing production Krishna Neupane overall-normalization uncertainty.
    The VPK and proton uncertainties are halved in relative size because the
    pi0 control channel contains two same-detector photons whereas DVCS contains
    one.  The result is explicitly a candidate normalization uncertainty, not
    yet the final production assignment: photon sector/momentum dependence and
    any differential Krishna-map transfer uncertainty remain validation items.
    """
    if not topology_path.exists() or local_fold.empty or proton_fold.empty:
        return pd.DataFrame()
    top=pd.read_csv(topology_path)
    top=top[(top.comparison=="raw") & top.photon_topology.isin(["FT-FT","FD-FD"])].copy()
    pf=krishna_uncertainty_summary(proton_fold)
    z=top.merge(local_fold,on=["period","photon_topology"],how="left").merge(
        pf,on=["period","photon_topology"],how="left",suffixes=("","_p"))
    R=pd.to_numeric(z.data_over_mc,errors="coerce")
    sR=pd.to_numeric(z.stat_err,errors="coerce")
    M=pd.to_numeric(z.local_M_acc,errors="coerce")
    sM=pd.to_numeric(z.local_M_acc_unc,errors="coerce")
    P=pd.to_numeric(z.krishna_epsilon_data_over_mc,errors="coerce")
    sPfrac=pd.to_numeric(z.krishna_overall_norm_sys_fraction,errors="coerce")
    r=np.sqrt(R/(M*P))
    C=1.0/r
    stat_frac=0.5*sR/R
    vpk_frac=0.5*sM/M
    proton_frac=0.5*sPfrac
    quantified=np.sqrt(stat_frac**2+vpk_frac**2+proton_frac**2)
    out=pd.DataFrame({
        "period":z.period,
        "photon_detector":np.where(z.photon_topology.eq("FT-FT"),"FT","FD"),
        "control_topology":z.photon_topology,
        "raw_eppi0_data_over_mc":R,
        "vpk_M_acc":M,
        "vpk_transport_unc_abs":sM,
        "krishna_proton_eff_data_over_mc":P,
        "krishna_overall_norm_sys_fraction":sPfrac,
        "single_photon_eff_data_over_mc":r,
        "candidate_dvcs_photon_correction":C,
        "candidate_norm_stat_fraction":stat_frac,
        "candidate_norm_vpk_fraction":vpk_frac,
        "candidate_norm_krishna_fraction":proton_frac,
        "candidate_norm_quantified_total_fraction":quantified,
        "pass1_reference_norm_fraction":0.31,
        "improvement_factor_vs_31pct":0.31/quantified,
        "final_production_ready":False,
        "remaining_unquantified":"photon sector/momentum dependence; differential Krishna-map transfer; phi-dependent VPK acceptance check",
    })
    return out


def plot_pass2_normalization_summary(summary: pd.DataFrame) -> None:
    if summary.empty: return
    q=summary.copy(); q["period"]=pd.Categorical(q.period,PERIOD_ORDER,ordered=True); q=q.sort_values(["period","photon_detector"])
    # Candidate one-photon correction.
    fig,ax=plt.subplots(figsize=(10.0,5.8))
    for det,dx in [("FD",-0.10),("FT",0.10)]:
        g=q[q.photon_detector==det]; x=np.arange(len(g))+dx
        y=g.candidate_dvcs_photon_correction.to_numpy(float)
        e=y*g.candidate_norm_quantified_total_fraction.to_numpy(float)
        ax.errorbar(x,y,yerr=e,fmt="o",capsize=4,label=det)
    ax.axhline(1.0,ls="--",lw=1)
    ax.set_xticks(np.arange(len(PERIOD_ORDER)),PERIOD_ORDER,rotation=25,ha="right")
    ax.set_ylabel("Candidate one-photon correction to DVCS yield")
    ax.set_title("eπ⁰-based pass-2 photon normalization candidate")
    ax.legend(frameon=False); fig.tight_layout()
    fig.savefig(PNG/"pass2_candidate_dvcs_photon_normalization.png",dpi=220); plt.close(fig)

    # Uncertainty budget, kept explicit rather than hiding the 31% comparison.
    fig,ax=plt.subplots(figsize=(10.0,5.8))
    labels=[]; x=[]; total=[]; vpk=[]; kp=[]
    ii=0
    for period in PERIOD_ORDER:
        for det in ("FD","FT"):
            g=q[(q.period==period)&(q.photon_detector==det)]
            if g.empty: continue
            r=g.iloc[0]; labels.append(f"{period} {det}"); x.append(ii); ii+=1
            total.append(100*r.candidate_norm_quantified_total_fraction)
            vpk.append(100*r.candidate_norm_vpk_fraction); kp.append(100*r.candidate_norm_krishna_fraction)
    ax.plot(x,total,"o",label="currently quantified total")
    ax.plot(x,vpk,"s",label="VPK transport contribution")
    ax.plot(x,kp,"^",label="Krishna Neupane contribution")
    ax.axhline(31.0,ls="--",lw=1,label="pass-1 31% reference")
    ax.set_xticks(x,labels,rotation=35,ha="right"); ax.set_ylabel("Relative normalization uncertainty (%)")
    ax.set_title("Pass-2 eπ⁰ normalization: currently quantified uncertainty budget")
    ax.legend(frameon=False,ncol=2); fig.tight_layout()
    fig.savefig(PNG/"pass2_candidate_normalization_uncertainty_budget.png",dpi=220); plt.close(fig)

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


def plot_transport_closure(scan: pd.DataFrame, closure: pd.DataFrame, empirical: pd.DataFrame):
    if scan.empty or closure.empty: return
    # Method scan: minimum RMSE over k for each bandwidth and vice versa.
    fig,ax=plt.subplots(figsize=(8.4,5.8))
    for k,g in scan.groupby("k"):
        ax.plot(g.bandwidth,g.clas6_setting_rmse,marker="o",ms=3,lw=1,label=f"k={k}")
    ax.set_xlabel("Gaussian bandwidth in normalized (xB,Q²,-t) distance")
    ax.set_ylabel("CLAS6 leave-setting-out RMS of predicted − measured M")
    ax.set_title("Transport-method selection from withheld CLAS6 settings")
    ax.grid(alpha=.25); ax.legend(ncol=2,fontsize=8)
    fig.tight_layout(); fig.savefig(PNG/"external_transport_method_scan.png",dpi=220); plt.close(fig)

    fig,ax=plt.subplots(figsize=(8.4,5.8))
    marks={"point":"o","setting":"s","CLAS6_to_HallA":"^"}
    for mode,g in closure.groupby("cv_mode"):
        ax.scatter(g.nearest_training_distance,g.residual,s=18,alpha=.65,marker=marks.get(mode,"o"),label=mode)
    ax.axhline(0,lw=1)
    if not empirical.empty:
        ax.plot(empirical.distance_max,empirical.rmse,"k--",lw=1.5,label="setting-holdout RMS")
        ax.plot(empirical.distance_max,-empirical.rmse,"k--",lw=1.5)
    ax.set_xlabel("Distance to nearest retained external measurement")
    ax.set_ylabel("Predicted M − measured M")
    ax.set_title("External π⁰ transport closure versus actual kinematic separation")
    ax.grid(alpha=.25); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(PNG/"external_transport_closure_vs_distance.png",dpi=220); plt.close(fig)

    ind=closure[closure.cv_mode=="CLAS6_to_HallA"]
    if not ind.empty:
        fig,ax=plt.subplots(figsize=(6.4,6.1))
        ax.errorbar(ind.measured_M,ind.predicted_M,xerr=ind.measured_M_err,fmt="o",ms=4,capsize=2)
        lo=min(ind.measured_M.min(),ind.predicted_M.min()); hi=max(ind.measured_M.max(),ind.predicted_M.max())
        ax.plot([lo,hi],[lo,hi],"k--",lw=1)
        ax.set_xlabel("Measured Hall-A data/VPK M")
        ax.set_ylabel("Prediction from CLAS6 residuals only")
        ax.set_title("Independent CLAS6 → Hall-A transport test")
        ax.grid(alpha=.25); fig.tight_layout()
        fig.savefig(PNG/"clas6_to_halla_transport_closure.png",dpi=220); plt.close(fig)


def plot_phase_space_uncertainty(setting_summary: pd.DataFrame, uscan: pd.DataFrame,
                                 a: pd.DataFrame):
    if setting_summary.empty:return
    # Setting-level closure errors versus the principal physics coordinates.
    for x,col,label in [("Q2_mean_GeV2","Q2_mean_GeV2",r"$Q^2$ (GeV$^2$)"),
                        ("xB_mean","xB_mean",r"$x_B$"),
                        ("minus_t_mean_GeV2","minus_t_mean_GeV2",r"$-t$ (GeV$^2$)")]:
        fig,ax=plt.subplots(figsize=(8.2,5.6))
        for exp,g in setting_summary.groupby("experiment"):
            ax.scatter(g[col],g.rmse_residual,s=48,label=exp)
        ax.set_xlabel(label); ax.set_ylabel("Whole-setting holdout RMS in data/VPK M")
        ax.set_title("Measured transport failure by external-data phase space")
        ax.grid(alpha=.25); ax.legend(fontsize=8); fig.tight_layout()
        fig.savefig(PNG/f"external_setting_closure_rmse_vs_{x.replace('_mean_GeV2','').replace('_mean','')}.png",dpi=220); plt.close(fig)

    fig,ax=plt.subplots(figsize=(8.4,5.8))
    sc=ax.scatter(setting_summary.xB_mean,setting_summary.Q2_mean_GeV2,
                  c=setting_summary.rmse_residual,s=70)
    ax.set_xlabel(r"$x_B$"); ax.set_ylabel(r"$Q^2$ (GeV$^2$)")
    ax.set_title("Whole-setting holdout error across measured π⁰ phase space")
    cbar=fig.colorbar(sc,ax=ax,pad=.03); cbar.set_label("RMS(predicted − measured M)")
    fig.tight_layout(); fig.savefig(PNG/"external_setting_closure_rmse_xB_Q2.png",dpi=220); plt.close(fig)

    if not uscan.empty:
        fig,ax=plt.subplots(figsize=(8.2,5.6))
        ax.plot(uscan.effective_k,uscan.nested_log_score,marker="o")
        ax.set_xlabel("Number of neighboring withheld settings used for local uncertainty")
        ax.set_ylabel("Nested leave-setting-out Gaussian log score (lower is better)")
        ax.set_title("Phase-space uncertainty complexity selected by withheld data")
        ax.grid(alpha=.25); fig.tight_layout()
        fig.savefig(PNG/"external_uncertainty_neighbor_scan.png",dpi=220); plt.close(fig)

    # What uncertainty does the final calibration assign to actual RGA cells?
    for topology in ["FD-FD","FT-FT"]:
        g=a[a.photon_topology==topology]
        fig,ax=plt.subplots(figsize=(8.4,5.8))
        sc=ax.scatter(g.xB,g.Q2,c=g.external_local_M_unc,s=8,alpha=.55)
        ax.set_xlabel(r"$x_B$"); ax.set_ylabel(r"$Q^2$ (GeV$^2$)")
        ax.set_title(f"{topology}: empirically calibrated local VPK-transport uncertainty")
        cbar=fig.colorbar(sc,ax=ax,pad=.03); cbar.set_label(r"$\delta M_{transport}$")
        fig.tight_layout(); fig.savefig(PNG/f"accepted_aao_{topology.replace('-','_')}_phase_space_uncertainty.png",dpi=220); plt.close(fig)

def _period_ebeam(period: str) -> float:
    return {"Sp18 Inb":10.594,"Sp18 Out":10.594,"Fa18 Inb":10.604,"Fa18 Out":10.604,"Sp19 Inb":10.200}[period]


def _angular_residual_training(clas6_path: Path) -> dict[str,pd.DataFrame]:
    """Dimensionless additive LT/TT residuals, scaled to VPK sigma_U.

    Ratios data/VPK are deliberately avoided for LT and TT because either
    interference structure function can cross zero.  The transported quantity
    is (data-model)/sigma_U(model) at the same published kinematics.
    """
    from validate_aao_pi0_model import epsilon, xsigma_t, xsigma_l
    p=pd.read_csv(clas6_path)
    out={}
    for obs in ["LT","TT"]:
        g=p[(p.observable==obs)&p.vpk_kine_valid].copy()
        su=[]
        for r in g.itertuples(index=False):
            ep=epsilon(float(r.xB),float(r.Q2_GeV2),float(r.Ebeam_GeV))
            su.append(xsigma_t(float(r.t_GeV2),float(r.xB),float(r.Q2_GeV2),float(r.Ebeam_GeV))+ep*xsigma_l(float(r.t_GeV2),float(r.xB),float(r.Q2_GeV2),float(r.Ebeam_GeV)))
        su=np.asarray(su,float)
        g["angular_residual_over_U"]=(pd.to_numeric(g.data_value,errors="coerce").to_numpy(float)-pd.to_numeric(g.vpk_value,errors="coerce").to_numpy(float))/su
        g["angular_residual_over_U_err"]=pd.to_numeric(g.data_totalerr,errors="coerce").to_numpy(float)/np.abs(su)
        g=g[np.isfinite(g.angular_residual_over_U)&np.isfinite(g.angular_residual_over_U_err)&(g.angular_residual_over_U_err>0)].copy()
        out[obs]=g
    return out


def phi_dependent_vpk_acceptance(a: pd.DataFrame, phi_path: Path, clas6_path: Path,
                                 bandwidth: float, k: int) -> pd.DataFrame:
    """Fold observed CLAS6 LT/TT residuals through the accepted AAOgen phi population.

    This is a diagnostic of the change relative to the present U-only external
    correction.  It uses the same externally selected transport bandwidth/k;
    no additional smoothing choice is introduced here.
    """
    if not phi_path.exists():
        print(f"[rga-coverage] accepted AAOgen phi population missing; skipping phi-dependent VPK diagnostic: {phi_path}")
        return pd.DataFrame()
    from validate_aao_pi0_model import epsilon, xsigma_t, xsigma_l, xsigma_tt, xsigma_lt
    ph=pd.read_csv(phi_path)
    base=a[["period","photon_topology","ix","iq","it","xB","Q2","minus_t","external_local_M"]].copy()
    z=ph.merge(base,on=["period","photon_topology","ix","iq","it"],how="inner")
    trains=_angular_residual_training(clas6_path)
    ext=pd.concat(list(trains.values()),ignore_index=True).rename(columns={"Q2_GeV2":"Q2","minus_t_GeV2":"minus_t"})
    lo,scale=_external_scale(ext)
    for obs in ["LT","TT"]:
        tr=trains[obs].copy(); xyz=tr[["xB","Q2_GeV2","minus_t_GeV2"]].to_numpy(float)
        tr2=pd.DataFrame({"xB":tr.xB,"Q2":tr.Q2_GeV2,"minus_t":tr.minus_t_GeV2,
                          "ratio":tr.angular_residual_over_U,"ratio_err":tr.angular_residual_over_U_err})
        pred=[]
        for q in z[["xB","Q2","minus_t"]].to_numpy(float):
            v,_=_predict_kernel(tr2,np.asarray(q,float)[None,:],lo,scale,bandwidth,k)
            pred.append(float(v[0]))
        z[f"delta_{obs}_over_U"]=pred
    ratios=[]
    for r in z.itertuples(index=False):
        x=float(r.xB); q2=float(r.Q2); t=-float(r.minus_t); e=_period_ebeam(str(r.period)); phi=math.radians(0.5*(float(r.phi_low_deg)+float(r.phi_high_deg)))
        ep=epsilon(x,q2,e); U=xsigma_t(t,x,q2,e)+ep*xsigma_l(t,x,q2,e); TT=xsigma_tt(t,x,q2,e); LT=xsigma_lt(t,x,q2,e)
        A=math.sqrt(max(0.0,2.0*ep*(1.0+ep)))
        nom=U + A*LT*math.cos(phi) + ep*TT*math.cos(2.0*phi)
        cor=float(r.external_local_M)*U + A*(LT+float(r.delta_LT_over_U)*U)*math.cos(phi) + ep*(TT+float(r.delta_TT_over_U)*U)*math.cos(2.0*phi)
        ratios.append(cor/nom if np.isfinite(nom) and abs(nom)>1e-12 else np.nan)
    z["full_angular_over_vpk"]=ratios
    z["u_only_over_vpk"]=z.external_local_M
    z["angular_relative_to_u_only"]=z.full_angular_over_vpk/z.u_only_over_vpk
    rows=[]
    for (period,topo),g in z.groupby(["period","photon_topology"]):
        good=g[np.isfinite(g.angular_relative_to_u_only)&(g.sum_weight>0)]
        if good.empty: continue
        w=good.sum_weight.to_numpy(float); rr=good.angular_relative_to_u_only.to_numpy(float)
        mean=float(np.sum(w*rr)/np.sum(w))
        rms=float(np.sqrt(np.sum(w*(rr-mean)**2)/np.sum(w)))
        rows.append({"period":period,"photon_topology":topo,"angular_over_u_only":mean,
                     "angular_shift_fraction":mean-1.0,"weighted_cell_phi_rms":rms,
                     "n_cell_phi_bins":len(good)})
    return pd.DataFrame(rows).sort_values(["period","photon_topology"])


def corrected_kinematic_closure(a: pd.DataFrame, grid_dir: Path) -> tuple[pd.DataFrame,pd.DataFrame]:
    """Correct coarse eppi0 DATA/MC cells for local VPK and proton terms, then test residual shape.

    The residual is normalized to its own period/topology weighted mean.  Thus
    this is a shape/closure test, not another determination of the absolute
    normalization.
    """
    files=sorted(grid_dir.glob("*_kinematic_closure_grid.csv"))
    if not files:
        print(f"[rga-coverage] no kinematic closure grids found under {grid_dir}; skipping corrected closure")
        return pd.DataFrame(),pd.DataFrame()
    g=pd.concat([pd.read_csv(f) for f in files],ignore_index=True)
    fine=a.copy()
    fine["cix"]=(fine.xB/0.1).astype(int).clip(0,6); fine["ciq"]=((fine.Q2-1.0)/1.0).astype(int).clip(0,5); fine["cit"]=(fine.minus_t/0.2).astype(int).clip(0,4)
    rows=[]
    for key,h in fine.groupby(["period","photon_topology","cix","ciq","cit"]):
        w=h.sum_weight.to_numpy(float); sw=w.sum()
        if sw<=0: continue
        rows.append({"period":key[0],"photon_topology":key[1],"ix":key[2],"iq":key[3],"it":key[4],
                     "local_M":float(np.sum(w*h.external_local_M.to_numpy(float))/sw),
                     "proton_ratio":float(np.sum(w*h.proton_efficiency_ratio_mean.to_numpy(float))/sw)})
    f=pd.DataFrame(rows)
    z=g.merge(f,on=["period","photon_topology","ix","iq","it"],how="left")
    z["corrected_two_photon_residual"]=z.data_over_mc/(z.local_M*z.proton_ratio)
    summaries=[]
    z["closure_normalized_residual"]=np.nan
    for (period,topo),idx in z.groupby(["period","photon_topology"]).groups.items():
        h=z.loc[idx].copy(); good=np.isfinite(h.corrected_two_photon_residual)&np.isfinite(h.stat_err)&(h.stat_err>0)
        h=h[good]
        if h.empty: continue
        # MC-yield weighting asks whether the accepted control-channel population closes as a whole,
        # while retaining cell-level residuals for trend plots.
        w=h.mc_sum_weight.to_numpy(float); rr=h.corrected_two_photon_residual.to_numpy(float)
        cen=float(np.sum(w*rr)/np.sum(w)); nr=rr/cen
        z.loc[h.index,"closure_normalized_residual"]=nr
        rms=float(np.sqrt(np.sum(w*(nr-1.0)**2)/np.sum(w)))
        p90=float(weighted_quantile(np.abs(nr-1.0),w,[0.90])[0])
        summaries.append({"period":period,"photon_topology":topo,"n_cells":len(h),"central_residual":cen,
                          "weighted_rms_fraction":rms,"weighted_p90_abs_fraction":p90,
                          "max_abs_fraction":float(np.max(np.abs(nr-1.0)))})
    return z,pd.DataFrame(summaries).sort_values(["period","photon_topology"])


def plot_phi_acceptance_diagnostic(df: pd.DataFrame) -> None:
    if df.empty:return
    q=df[df.photon_topology.isin(["FT-FT","FD-FD"])].copy()
    fig,ax=plt.subplots(figsize=(9.2,5.8))
    x=np.arange(len(PERIOD_ORDER)); width=.34
    for j,topo in enumerate(["FT-FT","FD-FD"]):
        vals=[]
        for p in PERIOD_ORDER:
            h=q[(q.period==p)&(q.photon_topology==topo)]
            vals.append(100*float(h.angular_shift_fraction.iloc[0]) if len(h) else np.nan)
        ax.bar(x+(j-.5)*width,vals,width,label=topo)
    ax.axhline(0,linewidth=1); ax.set_xticks(x,PERIOD_ORDER,rotation=20,ha="right")
    ax.set_ylabel("Change from U-only external correction (%)")
    ax.set_title("LT/TT discrepancies folded through accepted AAOgen phi")
    ax.legend(); ax.grid(axis="y",alpha=.2); fig.tight_layout(); fig.savefig(PNG/"phi_dependent_vpk_acceptance_shift.png",dpi=220); plt.close(fig)


def plot_corrected_closure(z: pd.DataFrame) -> None:
    if z.empty:return
    for var,label in [("xB_low",r"$x_B$"),("Q2_low_GeV2",r"$Q^2$ (GeV$^2$)"),("minus_t_low_GeV2",r"$-t$ (GeV$^2$)")]:
        fig,axes=plt.subplots(2,3,figsize=(14,8.5),sharey=True); axes=axes.ravel()
        for ax,p in zip(axes,PERIOD_ORDER):
            h=z[(z.period==p)&z.photon_topology.isin(["FT-FT","FD-FD"])&np.isfinite(z.closure_normalized_residual)]
            for topo,marker in [("FT-FT","o"),("FD-FD","s")]:
                q=h[h.photon_topology==topo]
                if q.empty:continue
                # weighted mean in each displayed coordinate bin
                pts=[]
                for xv,g in q.groupby(var):
                    w=g.mc_sum_weight.to_numpy(float); y=g.closure_normalized_residual.to_numpy(float)
                    pts.append((xv,float(np.sum(w*y)/np.sum(w))))
                if pts:
                    xx,yy=zip(*pts); ax.plot(xx,yy,marker=marker,label=topo)
            ax.axhline(1,linestyle="--",linewidth=1); ax.set_title(p); ax.grid(alpha=.2)
        for ax in axes[len(PERIOD_ORDER):]:ax.axis("off")
        axes[0].legend(); fig.supxlabel(label); fig.supylabel("Corrected eπ0 residual / period-topology mean")
        fig.suptitle("Corrected eπ0 kinematic closure",y=.995); fig.tight_layout(rect=[0,0,1,.97]); fig.savefig(PNG/f"corrected_eppi0_closure_vs_{var}.png",dpi=220); plt.close(fig)

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
    scan,closure,empirical,best_bw,best_k=calibrate_external_transport(clas6,halla)
    scan.to_csv(OUT/"external_transport_method_scan.csv",index=False)
    closure.to_csv(OUT/"external_transport_closure_predictions.csv",index=False)
    empirical.to_csv(OUT/"external_transport_distance_uncertainty.csv",index=False)
    setting_summary=setting_closure_summary(closure)
    setting_summary.to_csv(OUT/"external_transport_setting_closure_summary.csv",index=False)
    ext_all=pd.concat([clas6,halla],ignore_index=True)
    uscan,unc_k,global_sigma,max_unc_dist=calibrate_phase_space_uncertainty(setting_summary,ext_all)
    uscan.to_csv(OUT/"external_transport_uncertainty_model_scan.csv",index=False)
    print(f"[rga-coverage] central transport selected from CLAS6 setting closure: bandwidth={best_bw:g}, k={best_k}")
    print(f"[rga-coverage] phase-space uncertainty selected by nested setting closure: k={unc_k}, global RMS={global_sigma:.4f}, fallback beyond d={max_unc_dist:.4f}")
    a=add_local_external_field(a,clas6,halla,best_bw,best_k,setting_summary,unc_k,global_sigma,max_unc_dist)

    cov=coverage_summary(a); fold=fold_summary(a); local_fold=local_fold_summary(a)
    proton_fold=proton_fold_summary(a)
    proton_fold.to_csv(OUT/"rga_krishna_proton_efficiency_fold.csv",index=False)
    model_ratio,model_variants=topology_model_ratio_jackknife(a,clas6,halla,best_bw,best_k)
    model_ratio.to_csv(OUT/"rga_topology_model_ratio_uncertainty.csv",index=False)
    model_variants.to_csv(OUT/"rga_topology_model_ratio_setting_omissions.csv",index=False)
    relative=relative_photon_efficiency(args.topology_summary,model_ratio,proton_fold)
    if not relative.empty:
        relative.to_csv(OUT/"rga_relative_ft_fd_photon_efficiency.csv",index=False)
    absolute_p=absolute_photon_with_proton(local_fold,args.topology_summary,proton_fold)
    if not absolute_p.empty:
        absolute_p.to_csv(OUT/"rga_absolute_photon_efficiency_with_krishna.csv",index=False)
    krishna_unc=krishna_uncertainty_summary(proton_fold)
    if not krishna_unc.empty:
        krishna_unc.to_csv(OUT/"rga_krishna_proton_efficiency_with_uncertainty.csv",index=False)
    norm_summary=pass2_normalization_summary(local_fold,args.topology_summary,proton_fold)
    if not norm_summary.empty:
        norm_summary.to_csv(OUT/"pass2_eppi0_normalization_candidate_summary.csv",index=False)
    cov.to_csv(OUT/"rga_external_pi0_coverage_summary.csv",index=False)
    interference=interference_closure(CLAS6_POINTS)
    interference.to_csv(OUT/"clas6_U_TT_LT_closure_summary.csv",index=False)
    phi_diag=phi_dependent_vpk_acceptance(a, DVCS_DIR/"output"/"data_mc_normalization"/"accepted_aao_phi_population.csv", CLAS6_POINTS, best_bw, best_k)
    if not phi_diag.empty:
        phi_diag.to_csv(OUT/"phi_dependent_vpk_acceptance_diagnostic.csv",index=False)
    closure_cells,closure_summary=corrected_kinematic_closure(a, DVCS_DIR/"output"/"data_mc_normalization"/"photon_topology")
    if not closure_cells.empty:
        closure_cells.to_csv(OUT/"corrected_eppi0_kinematic_closure_cells.csv",index=False)
        closure_summary.to_csv(OUT/"corrected_eppi0_kinematic_closure_summary.csv",index=False)
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
    plot_transport_closure(scan,closure,empirical)
    plot_phase_space_uncertainty(setting_summary,uscan,a)
    plot_topology_model_cancellation(model_ratio)
    plot_krishna_and_relative(proton_fold,relative)
    plot_pass2_normalization_summary(norm_summary)
    plot_phi_acceptance_diagnostic(phi_diag)
    plot_corrected_closure(closure_cells)

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
        print("\n=== Empirically calibrated transport robustness test (all accepted kinematics) ===")
        print(robust.to_string(index=False))
    if not cand.empty:
        print("\n=== Candidate single-photon residuals (diagnostic; support fraction shown) ===")
        print(cand.to_string(index=False))
    if not proton_fold.empty:
        print("\n=== Krishna proton efficiency folded through accepted AAOgen ===")
        print(proton_fold[proton_fold.photon_topology.isin(["FT-FT","FD-FD"])].to_string(index=False))
    if not model_ratio.empty:
        print("\n=== Common external-model cancellation in M_FT/M_FD ===")
        print(model_ratio.to_string(index=False))
    if not relative.empty:
        print("\n=== Relative photon efficiency after Krishna: r_FT/r_FD ===")
        print(relative.to_string(index=False))
    if not norm_summary.empty:
        print("\n=== Pass-2 eppi0 normalization candidate (NOT final production assignment) ===")
        cols=["period","photon_detector","single_photon_eff_data_over_mc","candidate_dvcs_photon_correction",
              "candidate_norm_vpk_fraction","candidate_norm_krishna_fraction","candidate_norm_quantified_total_fraction","improvement_factor_vs_31pct"]
        print(norm_summary[cols].to_string(index=False))
        print("NOTE: quantified total excludes photon sector/momentum dependence, differential Krishna-map transfer, and phi-dependent VPK acceptance effects.")
    if not phi_diag.empty:
        print("\n=== Phi-dependent VPK acceptance diagnostic: LT/TT relative to U-only ===")
        print(phi_diag[phi_diag.photon_topology.isin(["FT-FT","FD-FD"])].to_string(index=False))
    if not closure_summary.empty:
        print("\n=== Corrected eppi0 kinematic closure ===")
        print(closure_summary[closure_summary.photon_topology.isin(["FT-FT","FD-FD"])].to_string(index=False))
    print(f"\nWrote tables under: {OUT}")
    print(f"Wrote PNG figures under: {PNG}")

if __name__=="__main__":
    main()
