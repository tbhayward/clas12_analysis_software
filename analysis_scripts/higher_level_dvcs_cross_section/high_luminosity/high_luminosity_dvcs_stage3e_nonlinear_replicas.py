#!/usr/bin/env python3
"""
Stage 3E: nonlinear replica validation of the high-luminosity D-term and
mechanical-structure projection.

Purpose
-------
Stage 3C used a local covariance around the Volker/Burkert generalized
multipole truth

    d1^Q(t) = -0.9 C (1 - t/M^2)^(-alpha),

and found that M^2 and alpha are individually weakly determined and almost
perfectly correlated.  Stage 3D then propagated that local Gaussian covariance
through the Fourier transform.  This script checks the result without assuming
that (C,M^2,alpha) are Gaussian.

Each pseudo-experiment:
  1. fluctuates every XS and BSA point with its projected statistical + PTP
     uncertainty and with common XS-normalization / BSA-polarization shifts;
  2. refits C, M^2 and alpha *nonlinearly*, simultaneously with the independent
     active KM15 ImH directions and the two correlated experimental nuisances;
  3. stores d1(t) from the fitted nonlinear parameters;
  4. separates a valid observable refit from a safe unrestricted high-q
     extrapolation. Fits that converge at a broad M^2/alpha diagnostic bound
     are retained for measured-range d1(t), rather than silently discarded;
  5. Fourier transforms every valid replica only up to the controlled CLAS12
     endpoint q_data=sqrt(|t|max). These are the primary, data-anchored
     mechanical bands;
  6. treats the unrestricted generalized-multipole transform as a separate
     model-continuation diagnostic, using only fits away from shape bounds;
  7. provides a simple cumulative-q comparison (data endpoint, 2 GeV, full).

Performance strategy
--------------------
A full Gepard evaluation of ~900 points for every optimizer iteration would be
prohibitively expensive.  Stage 3C already computed exact finite derivatives
of every observable.  We use those to build a fast response surrogate:

  * the D-term is nonlinear exactly in (C,M^2,alpha);
  * at each measured point, its arbitrary subtraction-function change DeltaS(t)
    is mapped to XS/BSA using dO/dS inferred from Stage-3C dO/dC;
  * the remaining active KM15 ImH directions use their validated local finite
    derivatives.

Thus the weakly constrained D-term shape is treated nonlinearly -- the issue we
need to validate -- while the expensive detector-kinematics response is cached.
A dedicated surrogate-closure table compares the nonlinear surrogate derivative
at the truth to the original Stage-3C derivatives.  This script does NOT claim
that the surrogate replaces an eventual smaller exact-Gepard spot check.

Output is standardized under output/stage3e/{tables,figures}.
"""
from pathlib import Path
import argparse, math, os
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from scipy.special import spherical_jn
from scipy.integrate import simpson

LUMI_FACTORS=(1,2,5,10)
SCENARIOS=("statistics_only","ptp_half","baseline")
DTERM_NAMES=("C","dM2","dalpha")
IMH_NAMES=("Nsea","Nv","alv","alpv","mv2")
FIT_NAMES=DTERM_NAMES+IMH_NAMES+("beta_xs_norm","beta_bsa_pol")
MP=0.9382720813
HBARC=0.1973269804
GEV4_TO_GEV_FM3=1.0/HBARC**3

# Volker/Burkert central values used in Stage 3C.
C0=2.04/0.9
M20=1.02
ALPHA0=2.76


def scenario_input(d, scenario):
    q=d.copy()
    if scenario=="baseline": return q
    if scenario=="ptp_half":
        q["xs_ptp_sys_pseudo_abs"]*=0.5; q["bsa_ptp_sys_pseudo_abs"]*=0.5
        return q
    if scenario=="statistics_only":
        q["xs_ptp_sys_pseudo_abs"]=0.0; q["bsa_ptp_sys_pseudo_abs"]=0.0
        q["xs_scale_frac"]=0.0; q["bsa_scale_frac"]=0.0
        return q
    raise ValueError(scenario)


def subtraction(tabs,C,M2,alpha):
    return C*(1.0+tabs/M2)**(-alpha)


def prepare_config(data, deriv, scenario, central_imh):
    d=scenario_input(data,scenario).merge(
        deriv,on=["point_id","bin","xB","Q2","t_abs","phi_deg"],
        how="inner",validate="one_to_one")
    if len(d)!=len(data): raise RuntimeError("Stage-2 / Stage-3C derivative merge lost points")
    tabs=d.t_abs.to_numpy(float); xs0=d.xs0.to_numpy(float); bsa0=d.bsa0.to_numpy(float)
    S0=subtraction(tabs,C0,M20,ALPHA0)
    # Since S=C*f(t), dS/dC=f=S/C.  Hence dO/dS=(dO/dC)/(S/C).
    f0=S0/C0
    dxs_dS=d.d_xs_d_C.to_numpy(float)/f0
    dbsa_dS=d.d_bsa_d_C.to_numpy(float)/f0
    DX=np.column_stack([d[f"d_xs_d_{p}"].to_numpy(float) for p in IMH_NAMES])
    DA=np.column_stack([d[f"d_bsa_d_{p}"].to_numpy(float) for p in IMH_NAMES])
    sx=np.hypot(d.xs_stat_pseudo_abs.to_numpy(float),d.xs_ptp_sys_pseudo_abs.to_numpy(float))
    sa=np.hypot(d.bsa_stat_pseudo_abs.to_numpy(float),d.bsa_ptp_sys_pseudo_abs.to_numpy(float))
    return dict(t=tabs,xs0=xs0,bsa0=bsa0,S0=S0,dxs_dS=dxs_dS,dbsa_dS=dbsa_dS,
                DX=DX,DA=DA,sx=sx,sa=sa,xscale=d.xs_scale_frac.to_numpy(float),
                ascale=d.bsa_scale_frac.to_numpy(float),central_imh=np.array([central_imh[p] for p in IMH_NAMES]))


def unpack_z(z,cfg):
    # Positive D-term shape parameters are optimized in log coordinates.
    C=C0*np.exp(z[0]); M2=M20*np.exp(z[1]); alpha=ALPHA0*np.exp(z[2])
    # Dimensionless nuisance coordinates: one unit = natural central parameter scale.
    scales=np.maximum(np.abs(cfg["central_imh"]),1.0)
    imh=cfg["central_imh"]+z[3:8]*scales
    return C,M2,alpha,imh,z[8],z[9]


def predict(z,cfg):
    C,M2,alpha,imh,bx,ba=unpack_z(z,cfg)
    dS=subtraction(cfg["t"],C,M2,alpha)-cfg["S0"]
    dimh=imh-cfg["central_imh"]
    xs=cfg["xs0"]+cfg["dxs_dS"]*dS+cfg["DX"]@dimh
    bsa=cfg["bsa0"]+cfg["dbsa_dS"]*dS+cfg["DA"]@dimh
    # Match Stage-3C correlated-nuisance convention (linear about central prediction).
    xs=xs+bx*cfg["xscale"]*cfg["xs0"]
    bsa=bsa+ba*cfg["ascale"]*cfg["bsa0"]
    return xs,bsa


def residual(z,cfg,yx,ya):
    px,pa=predict(z,cfg)
    r=np.r_[(px-yx)/cfg["sx"],(pa-ya)/cfg["sa"]]
    # Unit Gaussian priors on correlated scale nuisances when they are active.
    if np.any(cfg["xscale"]): r=np.r_[r,z[8]]
    if np.any(cfg["ascale"]): r=np.r_[r,z[9]]
    return r


def one_replica(payload):
    seed,cfg,max_nfev=payload
    rng=np.random.default_rng(seed)
    bx_true=rng.normal() if np.any(cfg["xscale"]) else 0.0
    ba_true=rng.normal() if np.any(cfg["ascale"]) else 0.0
    yx=cfg["xs0"]*(1.0+bx_true*cfg["xscale"])+rng.normal(0,cfg["sx"])
    ya=cfg["bsa0"]*(1.0+ba_true*cfg["ascale"])+rng.normal(0,cfg["sa"])
    # Broad log-bounds prevent unphysical M2/alpha while allowing the large non-Gaussian tails.
    lo=np.array([math.log(.05/C0),math.log(.05/M20),math.log(.20/ALPHA0),-8,-8,-8,-8,-8,-6,-6])
    hi=np.array([math.log(12/C0),math.log(12/M20),math.log(12/ALPHA0), 8, 8, 8, 8, 8, 6, 6])
    fit=least_squares(residual,np.zeros(10),args=(cfg,yx,ya),bounds=(lo,hi),
                      method="trf",xtol=2e-7,ftol=2e-7,gtol=2e-7,max_nfev=max_nfev,x_scale="jac")
    C,M2,alpha,imh,bx,ba=unpack_z(fit.x,cfg)
    at_bound=bool(np.any((fit.x-lo)<2e-4)|np.any((hi-fit.x)<2e-4))
    return (seed,fit.success,fit.cost*2,fit.nfev,at_bound,C,M2,alpha,*imh,bx,ba,bx_true,ba_true)


def percentile_band(x,axis=0):
    return np.nanpercentile(x,[16,50,84],axis=axis)


def d1_curve(C,M2,alpha,t):
    return -0.9*C*(1+t/M2)**(-alpha)


def mechanics_batch(C,M2,alpha,r,qmax=40.0,nq=5001,batch=128):
    """Fast vectorized Simpson transform for many replicas."""
    q=np.linspace(0,qmax,nq); rq=np.outer(r/HBARC,q)
    common=q**4
    Kp=spherical_jn(0,rq)*common
    Ks=spherical_jn(2,rq)*common
    P=[]; S=[]
    for i in range(0,len(C),batch):
        sl=slice(i,min(i+batch,len(C)))
        d1=-0.9*C[sl,None]*(1+q[None,:]**2/M2[sl,None])**(-alpha[sl,None])
        Ip=simpson(d1[:,None,:]*Kp[None,:,:],x=q,axis=2)
        Is=simpson(d1[:,None,:]*Ks[None,:,:],x=q,axis=2)
        P.append(-Ip/(15*np.pi**2*MP)*GEV4_TO_GEV_FM3)
        S.append(-Is/(10*np.pi**2*MP)*GEV4_TO_GEV_FM3)
    return np.vstack(P),np.vstack(S)


def cumulative_central(C,M2,alpha,r,qcuts,nq_per_GeV=600):
    rows=[]
    for qc in qcuts:
        nq=max(401,int(qc*nq_per_GeV)+1); q=np.linspace(0,qc,nq)
        d1=-0.9*C*(1+q*q/M2)**(-alpha); rr=np.outer(r/HBARC,q); w=q**4*d1
        p=-simpson(spherical_jn(0,rr)*w[None,:],x=q,axis=1)/(15*np.pi**2*MP)*GEV4_TO_GEV_FM3
        s=-simpson(spherical_jn(2,rr)*w[None,:],x=q,axis=1)/(10*np.pi**2*MP)*GEV4_TO_GEV_FM3
        for k,rv in enumerate(r): rows.append(dict(qmax=qc,r_fm=rv,pressure=p[k],shear=s[k]))
    return pd.DataFrame(rows)


def savefig(fig,p): fig.tight_layout(); fig.savefig(p,dpi=250,bbox_inches="tight"); plt.close(fig)


def make_plots(rep,bands,mech,cum,figdir,qdata):
    labels={"statistics_only":"Statistics only",
            "ptp_half":"Point-to-point systematics / 2",
            "baseline":"Current point-to-point systematics"}

    fig,ax=plt.subplots(figsize=(8.8,6))
    for L in LUMI_FACTORS:
        d=bands[(bands.scenario=="statistics_only")&(bands.luminosity_factor==L)]
        rel=50*(d.q84-d.q16)/np.maximum(np.abs(d.q50),1e-12)
        ax.plot(d.t_abs,rel,label=f"{L}x")
    ax.set(xlabel=r"$|t|$ (GeV$^2$)",
           ylabel=r"68% uncertainty on $d_1^Q(t)$ (%)",
           title="Nonlinear-replica D-term precision in the controlled CLAS12 range")
    ax.grid(alpha=.2); ax.legend(title="Pass-2 exposure")
    savefig(fig,figdir/"01_d1_replica_relative_precision.png")

    fig,ax=plt.subplots(figsize=(8.8,6))
    for sc in SCENARIOS:
        d=bands[(bands.scenario==sc)&(bands.luminosity_factor==10)]
        rel=50*(d.q84-d.q16)/np.maximum(np.abs(d.q50),1e-12)
        ax.plot(d.t_abs,rel,label=labels[sc])
    ax.set(xlabel=r"$|t|$ (GeV$^2$)",
           ylabel=r"68% uncertainty on $d_1^Q(t)$ (%)",
           title="Systematic limitation of the 10x D-term projection")
    ax.grid(alpha=.2); ax.legend()
    savefig(fig,figdir/"02_d1_replica_systematics_10x.png")

    # Backup only: demonstrates that M2 and alpha are poor coordinates for
    # describing the experimentally constrained D-term curve.
    d=rep[(rep.scenario=="statistics_only")&(rep.luminosity_factor==10)&rep.observable_valid]
    fig,ax=plt.subplots(figsize=(7,6))
    interior=d[~d.at_bound]; edge=d[d.at_bound]
    ax.scatter(interior.dM2,interior.dalpha,s=8,alpha=.22,label="interior fits")
    if len(edge):
        ax.scatter(edge.dM2,edge.dalpha,s=14,alpha=.45,marker="x",label="diagnostic-bound fits")
    ax.scatter([M20],[ALPHA0],marker="*",s=100,label="truth")
    ax.set(xlabel=r"$M^2$ (GeV$^2$)",ylabel=r"$\alpha$",
           title="Backup diagnostic: nonlinear shape-parameter degeneracy")
    ax.grid(alpha=.2); ax.legend()
    savefig(fig,figdir/"03_M2_alpha_replica_scatter_backup.png")

    if len(mech):
        # Workshop headline mechanics: full Nature/Volker-style generalized-
        # multipole transform. These are projections conditional on that
        # functional continuation, not model-independent pressure/shear.
        for quantity,title,ylabel,fname in [
            ("r2p_full","Pressure projection from the full fitted D-term",
             "$r^2p(r)$ (GeV fm$^{-1}$)","04_pressure_full_replicas.png"),
            ("r2s_full","Shear-stress projection from the full fitted D-term",
             "$r^2s(r)$ (GeV fm$^{-1}$)","05_shear_full_replicas.png")]:
            fig,ax=plt.subplots(figsize=(8.8,6))
            for L in LUMI_FACTORS:
                d=mech[(mech.scenario=="statistics_only")&
                       (mech.luminosity_factor==L)&(mech.quantity==quantity)]
                if not len(d): continue
                ax.fill_between(d.r_fm,d.q16,d.q84,alpha=.13)
                ax.plot(d.r_fm,d.q50,label=f"{L}x")
            ax.axhline(0,lw=.8,alpha=.5)
            ax.set(xlabel="r (fm)",ylabel=ylabel,title=title)
            ax.grid(alpha=.2); ax.legend(title="Pass-2 exposure")
            savefig(fig,figdir/fname)

        # Systematics comparison at 10x, again for the full fitted form.
        for quantity,title,ylabel,fname in [
            ("r2p_full","10x pressure projection: impact of point-to-point systematics",
             "$r^2p(r)$ (GeV fm$^{-1}$)","06_pressure_full_systematics_10x.png"),
            ("r2s_full","10x shear-stress projection: impact of point-to-point systematics",
             "$r^2s(r)$ (GeV fm$^{-1}$)","07_shear_full_systematics_10x.png")]:
            fig,ax=plt.subplots(figsize=(8.8,6))
            for sc in SCENARIOS:
                d=mech[(mech.scenario==sc)&(mech.luminosity_factor==10)&
                       (mech.quantity==quantity)]
                if not len(d): continue
                ax.fill_between(d.r_fm,d.q16,d.q84,alpha=.11)
                ax.plot(d.r_fm,d.q50,label=labels[sc])
            ax.axhline(0,lw=.8,alpha=.5)
            ax.set(xlabel="r (fm)",ylabel=ylabel,title=title)
            ax.grid(alpha=.2); ax.legend()
            savefig(fig,figdir/fname)

        # Compact stability diagnostics: fractions of full-transform replicas
        # with Nature-like pressure/shear topology. These are diagnostics, not
        # hard cuts on the workshop projection.
        stab=mech[mech.quantity=="stability_summary"].copy()
        if len(stab):
            fig,ax=plt.subplots(figsize=(8.8,5.8))
            d=stab[stab.scenario=="statistics_only"]
            ax.plot(d.luminosity_factor,100*d.pressure_naturelike_fraction,marker="o",
                    label="pressure: one + to - crossing")
            ax.plot(d.luminosity_factor,100*d.shear_positive_fraction,marker="o",
                    label="shear: non-negative")
            ax.set(xlabel="Pass-2 exposure factor",ylabel="Replica fraction (%)",
                   title="Mechanical-shape diagnostic for the full fitted D-term")
            ax.set_xticks(LUMI_FACTORS); ax.set_ylim(0,105)
            ax.grid(alpha=.2); ax.legend()
            savefig(fig,figdir/"08_mechanical_stability_diagnostic.png")

    # Simple central-value momentum-support diagnostic: only 3 curves.
    if len(cum):
        qvals=sorted(cum.qmax.unique())
        qfull=max(qvals)
        qmid=min(qvals,key=lambda x:abs(x-2.0))
        qdat=min(qvals,key=lambda x:abs(x-qdata))
        selected=[qdat,qmid,qfull]
        for col,title,fname in [
            ("pressure","Pressure: measured-range contribution versus model continuation",
             "09_pressure_momentum_support_diagnostic.png"),
            ("shear","Shear: measured-range contribution versus model continuation",
             "10_shear_momentum_support_diagnostic.png")]:
            fig,ax=plt.subplots(figsize=(8.8,6))
            for qc in selected:
                d=cum[np.isclose(cum.qmax,qc)]
                if np.isclose(qc,qdat):
                    lab=f"Controlled CLAS12 contribution (q < {qc:.3f} GeV)"
                elif np.isclose(qc,qmid):
                    lab="Extended to q < 2 GeV"
                else:
                    lab=f"Full assumed multipole (q < {qc:g} GeV)"
                ax.plot(d.r_fm,d[col],label=lab)
            ax.axhline(0,lw=.8,alpha=.5)
            ax.set(xlabel="r (fm)",
                   ylabel=("p(r) (GeV fm$^{-3}$)" if col=="pressure" else "s(r) (GeV fm$^{-3}$)"),
                   title=title)
            ax.grid(alpha=.2); ax.legend()
            savefig(fig,figdir/fname)


def main():
    here=Path(__file__).resolve().parent
    ap=argparse.ArgumentParser()
    ap.add_argument("--stage2",default=str(here/"output"/"stage2"/"tables"))
    ap.add_argument("--stage3c",default=str(here/"output"/"stage3c"/"tables"))
    ap.add_argument("--outdir",default=str(here/"output"/"stage3e"))
    ap.add_argument("--replicas",type=int,default=300,help="replicas per luminosity/scenario (default 300)")
    ap.add_argument("--workers",type=int,default=min(8,os.cpu_count() or 1))
    ap.add_argument("--seed",type=int,default=381902)
    ap.add_argument("--max-nfev",type=int,default=350)
    ap.add_argument("--qmax",type=float,default=40.0)
    ap.add_argument("--mechanics-nq",type=int,default=5001)
    ap.add_argument("--run-high-t-ablation",action="store_true",
                    help="reserved backup diagnostic; not part of the main luminosity projection")
    args=ap.parse_args()
    s2=Path(args.stage2)
    if not s2.exists():
        legacy=here/"output_stage2_pass2"/"tables"
        if legacy.exists(): print(f"[input] using legacy Stage-2 directory {legacy}"); s2=legacy
    s3=Path(args.stage3c); out=Path(args.outdir); tab=out/"tables"; fig=out/"figures"; tab.mkdir(parents=True,exist_ok=True); fig.mkdir(parents=True,exist_ok=True)
    deriv=pd.read_csv(s3/"volker_observable_parameter_derivatives.csv")
    try:
        from gepard.fits import th_KM15
        central_imh={p:float(th_KM15.parameters[p]) for p in IMH_NAMES}
    except Exception as e: raise RuntimeError("Stage 3E needs the same Gepard environment as Stage 3C to read KM15 central parameters") from e
    print("="*96); print("STAGE 3E: NONLINEAR D-TERM REFIT REPLICAS + MECHANICAL-STRUCTURE VALIDATION"); print("="*96)
    print(f"replicas/config={args.replicas}; workers={args.workers}; 12 configurations")
    print("nonlinear parameters: C, M^2, alpha; cached local response: active KM15 ImH directions")
    print("correlated nuisances: XS normalization and BSA polarization generated + refitted")
    print("mechanics: empirical replica percentiles; cumulative-q diagnostic separates data reach from extrapolation")
    rng=np.random.default_rng(args.seed); allrows=[]
    for L in LUMI_FACTORS:
        data=pd.read_csv(s2/f"joint_fit_input_km15_{L}x.csv"); data["bin"]=data["bin"].astype(int)
        for sc in SCENARIOS:
            cfg=prepare_config(data,deriv,sc,central_imh)
            seeds=rng.integers(1,2**63-1,size=args.replicas,dtype=np.int64)
            payload=[(int(x),cfg,args.max_nfev) for x in seeds]
            if args.workers>1:
                ctx=mp.get_context("spawn")
                with ProcessPoolExecutor(max_workers=args.workers,mp_context=ctx) as ex: vals=list(ex.map(one_replica,payload,chunksize=max(1,args.replicas//(args.workers*4))))
            else: vals=[one_replica(x) for x in payload]
            for v in vals:
                row=dict(zip(["seed","success","chi2","nfev","at_bound","C","dM2","dalpha",*IMH_NAMES,"beta_xs_norm","beta_bsa_pol","beta_xs_true","beta_bsa_true"],v))
                row.update(scenario=sc,luminosity_factor=L)
                # Hitting a deliberately broad shape bound is not an observable
                # fit failure. Keep it for d1(t) inside the controlled range,
                # but do not use it for unrestricted high-q continuation.
                row["observable_valid"]=bool(row["success"] and np.isfinite(row["chi2"]))
                row["extrapolation_safe"]=bool(row["observable_valid"] and not row["at_bound"])
                row["accepted"]=row["observable_valid"]  # compatibility with old tables
                allrows.append(row)
            subset=[r for r in allrows if r["scenario"]==sc and r["luminosity_factor"]==L]
            nvalid=sum(r["observable_valid"] for r in subset)
            nsafe=sum(r["extrapolation_safe"] for r in subset)
            print(f"[replicas] {sc:15s} {L:2d}x: observable-valid {nvalid}/{args.replicas}; "
                  f"interior/high-q-safe {nsafe}/{args.replicas}")
    rep=pd.DataFrame(allrows); rep.to_csv(tab/"replica_fit_results.csv",index=False)

    # D-term bands: restrict the displayed/result grid to the actually controlled
    # CLAS12 t range. All converged observable fits contribute, including fits
    # whose arbitrary M2/alpha coordinates land at a diagnostic bound.
    tmin=float(deriv.t_abs.min()); tmax=float(deriv.t_abs.max())
    tg=np.linspace(tmin,tmax,181); brows=[]
    for (sc,L),g in rep[rep.observable_valid].groupby(["scenario","luminosity_factor"]):
        vals=np.array([d1_curve(r.C,r.dM2,r.dalpha,tg) for r in g.itertuples(index=False)])
        q16,q50,q84=percentile_band(vals)
        for k,t in enumerate(tg):
            brows.append(dict(scenario=sc,luminosity_factor=L,t_abs=t,
                              q16=q16[k],q50=q50[k],q84=q84[k],
                              n_replicas=len(g),controlled_tmin=tmin,controlled_tmax=tmax))
    bands=pd.DataFrame(brows); bands.to_csv(tab/"d1_replica_bands.csv",index=False)

    sums=[]
    for (sc,L),g in rep[rep.observable_valid].groupby(["scenario","luminosity_factor"]):
        for p,truth in [("C",C0),("dM2",M20),("dalpha",ALPHA0)]:
            q=np.percentile(g[p],[16,50,84])
            sums.append(dict(scenario=sc,luminosity_factor=L,parameter=p,truth=truth,
                             q16=q[0],median=q[1],q84=q[2],n=len(g),
                             bound_hit_fraction=float(g.at_bound.mean())))
    pd.DataFrame(sums).to_csv(tab/"replica_parameter_summaries.csv",index=False)

    # Workshop mechanics projection: follow the Nature/Volker procedure and
    # transform the complete fitted generalized-multipole form.  We retain all
    # observable-valid replicas for d1(t), but for an unrestricted Fourier
    # transform we use only replicas that stay away from the deliberately broad
    # M2/alpha diagnostic bounds.  This is a model-conditional projection.
    qdata=float(np.sqrt(tmax)); rgrid=np.linspace(.05,2.0,99); mrows=[]
    wanted={("statistics_only",L) for L in LUMI_FACTORS}|{(sc,10) for sc in SCENARIOS}
    for sc,L in sorted(wanted,key=lambda x:(x[1],x[0])):
        g=rep[(rep.scenario==sc)&(rep.luminosity_factor==L)&rep.extrapolation_safe]
        if len(g)<20: continue
        P,S=mechanics_batch(g.C.to_numpy(),g.dM2.to_numpy(),g.dalpha.to_numpy(),
                            rgrid,args.qmax,args.mechanics_nq)
        R2P=P*rgrid[None,:]**2
        R2S=S*rgrid[None,:]**2
        for name,A in [("r2p_full",R2P),("r2s_full",R2S),
                       ("pressure_full",P),("shear_full",S)]:
            q16,q50,q84=percentile_band(A)
            for k,rv in enumerate(rgrid):
                mrows.append(dict(scenario=sc,luminosity_factor=L,quantity=name,
                                  r_fm=rv,q16=q16[k],q50=q50[k],q84=q84[k],
                                  n_replicas=len(g),qmax_GeV=args.qmax,
                                  interpretation="full generalized-multipole projection"))

        # Mechanical-shape sanity checks. Nature-like pressure means positive
        # at small r and exactly one + -> - crossing before 2 fm. Shear should
        # remain non-negative over the displayed radial interval. We report,
        # rather than impose, these conditions.
        pressure_ok=[]; shear_ok=[]; vonlaue=[]
        dr=np.gradient(rgrid)
        for pp,ss in zip(P,S):
            sig=np.sign(pp)
            # Ignore numerical zeros; count robust sign changes.
            cross=np.where(sig[:-1]*sig[1:]<0)[0]
            naturelike=(pp[0]>0 and len(cross)==1 and pp[-1]<=0)
            pressure_ok.append(naturelike)
            shear_ok.append(bool(np.nanmin(ss)>=-1e-8))
            # Finite-r diagnostic only; full von-Laue check is also written
            # separately below on an extended r grid.
            vonlaue.append(float(np.sum(rgrid**2*pp*dr)))
        mrows.append(dict(scenario=sc,luminosity_factor=L,quantity="stability_summary",
                          r_fm=np.nan,q16=np.nan,q50=np.nan,q84=np.nan,
                          n_replicas=len(g),qmax_GeV=args.qmax,
                          interpretation="mechanical-shape diagnostic",
                          pressure_naturelike_fraction=float(np.mean(pressure_ok)),
                          shear_positive_fraction=float(np.mean(shear_ok)),
                          median_finite_range_vonlaue=float(np.median(vonlaue))))

    mech=pd.DataFrame(mrows)
    mech.to_csv(tab/"mechanics_replica_bands.csv",index=False)

    # A more direct central-truth sanity check, including a wider r interval for
    # the von-Laue integral. This verifies the implementation against the
    # expected Nature-like topology without using it as a fit constraint.
    rcheck=np.linspace(.01,8.0,800)
    P0,S0=mechanics_batch(np.array([C0]),np.array([M20]),np.array([ALPHA0]),
                          rcheck,args.qmax,args.mechanics_nq)
    P0=P0[0]; S0=S0[0]
    cross=np.where(np.sign(P0[:-1])*np.sign(P0[1:])<0)[0]
    zeros=[]
    for i in cross:
        x1,x2=rcheck[i],rcheck[i+1]; y1,y2=P0[i],P0[i+1]
        zeros.append(float(x1-y1*(x2-x1)/(y2-y1)))
    central_diag=pd.DataFrame([dict(
        pressure_zero_crossings=len(zeros),
        first_pressure_zero_fm=(zeros[0] if zeros else np.nan),
        shear_min_GeV_fm3=float(np.min(S0)),
        shear_nonnegative=bool(np.min(S0)>=-1e-8),
        von_laue_integral_GeV=float(np.trapz(rcheck**2*P0,rcheck)),
        qmax_GeV=args.qmax,
        note="full fitted generalized-multipole central truth")])
    central_diag.to_csv(tab/"mechanics_central_stability_check.csv",index=False)

    # Simplified support diagnostic: controlled endpoint, 2 GeV, and effectively full.
    qcuts=sorted(set([qdata,2.0,args.qmax]))
    cum=cumulative_central(C0,M20,ALPHA0,rgrid,qcuts)
    cum["q_data_endpoint_GeV"]=qdata
    cum.to_csv(tab/"mechanics_cumulative_q_support.csv",index=False)
    make_plots(rep,bands,mech,cum,fig,qdata)

    diag=rep.groupby(["scenario","luminosity_factor"]).agg(
        attempted=("observable_valid","size"),
        observable_valid=("observable_valid","sum"),
        extrapolation_safe=("extrapolation_safe","sum"),
        median_nfev=("nfev","median"),
        bound_hits=("at_bound","sum")).reset_index()
    diag["observable_valid_fraction"]=diag.observable_valid/diag.attempted
    diag["extrapolation_safe_fraction"]=diag.extrapolation_safe/diag.attempted
    diag.to_csv(tab/"replica_diagnostics.csv",index=False)
    print("\nReplica diagnostics:"); print(diag.to_string(index=False))
    print(f"\n[data/Fourier] controlled endpoint q_data=sqrt(tmax)={qdata:.3f} GeV; transforms also evaluated above this to expose model continuation")
    print(f"[output] {out}")
    print("[interpretation] d1 bands use every converged observable fit inside the controlled t range.")
    print("[interpretation] mechanics headline follows Nature/Volker: full fitted generalized-multipole transform, explicitly model-conditional.")
    print("[interpretation] hard q_data truncation is used only as a momentum-support diagnostic; it is not interpreted as a physical pressure/shear distribution.")
    if args.run_high_t_ablation:
        print("[high-t ablation] not executed in the main workflow: a real luminosity upgrade improves the full accepted kinematic range.")
        print("[high-t ablation] flag retained only so a future targeted diagnostic can be added without changing the main physics projection.")

if __name__=="__main__": main()
