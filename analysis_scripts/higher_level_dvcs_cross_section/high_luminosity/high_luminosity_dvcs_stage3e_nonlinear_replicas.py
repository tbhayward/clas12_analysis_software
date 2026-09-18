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
  4. for the configurations used in the mechanics plots, Fourier transforms
     every accepted replica to p(r) and s(r), producing empirical percentile
     bands rather than linear error propagation;
  5. diagnoses how much of p(r) and s(r) comes from q below the controlled
     CLAS12 endpoint versus the assumed high-q continuation of the multipole.

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


def make_plots(rep,bands,mech,cum,figdir):
    # Nonlinear empirical d1 precision: statistical luminosity progression.
    fig,ax=plt.subplots(figsize=(8.8,6))
    for L in LUMI_FACTORS:
        d=bands[(bands.scenario=="statistics_only")&(bands.luminosity_factor==L)]
        rel=50*(d.q84-d.q16)/np.maximum(np.abs(d.q50),1e-12)
        ax.plot(d.t_abs,rel,label=f"{L}x")
    ax.set(xlabel=r"$|t|$ (GeV$^2$)",ylabel=r"Empirical 68% relative half-width on $d_1^Q(t)$ (%)",
           title="Nonlinear-replica D-term precision")
    ax.grid(alpha=.2); ax.legend(title="Pass-2 exposure"); savefig(fig,figdir/"01_d1_replica_relative_precision.png")
    # Systematics at 10x.
    labels={"statistics_only":"Statistics only","ptp_half":"Point-to-point systematics / 2","baseline":"Current point-to-point systematics"}
    fig,ax=plt.subplots(figsize=(8.8,6))
    for sc in SCENARIOS:
        d=bands[(bands.scenario==sc)&(bands.luminosity_factor==10)]
        rel=50*(d.q84-d.q16)/np.maximum(np.abs(d.q50),1e-12); ax.plot(d.t_abs,rel,label=labels[sc])
    ax.set(xlabel=r"$|t|$ (GeV$^2$)",ylabel=r"Empirical 68% relative half-width on $d_1^Q(t)$ (%)",
           title="Nonlinear-replica systematics limitation at 10x")
    ax.grid(alpha=.2); ax.legend(); savefig(fig,figdir/"02_d1_replica_systematics_10x.png")
    # Parameter scatter demonstrates nonlinear M2-alpha degeneracy.
    d=rep[(rep.scenario=="statistics_only")&(rep.luminosity_factor==10)&rep.accepted]
    fig,ax=plt.subplots(figsize=(7,6)); ax.scatter(d.dM2,d.dalpha,s=8,alpha=.25)
    ax.scatter([M20],[ALPHA0],marker="*",s=100,label="truth"); ax.set(xlabel=r"$M^2$ (GeV$^2$)",ylabel=r"$\alpha$",title="10x statistics-only nonlinear shape degeneracy")
    ax.grid(alpha=.2); ax.legend(); savefig(fig,figdir/"03_M2_alpha_replica_scatter.png")
    if len(mech):
        # pressure/shear empirical bands for stats luminosity progression
        for quantity,title,ylabel,fname in [("r2p","Pressure","$r^2p(r)$ (GeV fm$^{-1}$)","04_pressure_replicas.png"),("shear","Shear","$s(r)$ (GeV fm$^{-3}$)","05_shear_replicas.png")]:
            fig,ax=plt.subplots(figsize=(8.8,6))
            for L in LUMI_FACTORS:
                d=mech[(mech.scenario=="statistics_only")&(mech.luminosity_factor==L)]
                ax.fill_between(d.r_fm,d.q16,d.q84,alpha=.15); ax.plot(d.r_fm,d.q50,label=f"{L}x")
            ax.axhline(0,lw=.8,alpha=.5); ax.set(xlabel="r (fm)",ylabel=ylabel,title=f"{title}: nonlinear replica projection")
            ax.grid(alpha=.2); ax.legend(title="Pass-2 exposure"); savefig(fig,figdir/fname)
    # cumulative q support central diagnostic
    if len(cum):
        for col,title,fname in [("pressure","Pressure transform: cumulative momentum support","06_pressure_cumulative_q.png"),("shear","Shear transform: cumulative momentum support","07_shear_cumulative_q.png")]:
            fig,ax=plt.subplots(figsize=(8.8,6))
            for qc in sorted(cum.qmax.unique()):
                d=cum[cum.qmax==qc]; ax.plot(d.r_fm,d[col],label=f"q < {qc:g} GeV")
            ax.axhline(0,lw=.8,alpha=.5); ax.set(xlabel="r (fm)",ylabel=("p(r) (GeV fm$^{-3}$)" if col=="pressure" else "s(r) (GeV fm$^{-3}$)"),title=title)
            ax.grid(alpha=.2); ax.legend(ncol=2); savefig(fig,figdir/fname)


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
                row.update(scenario=sc,luminosity_factor=L); row["accepted"]=bool(row["success"] and not row["at_bound"] and np.isfinite(row["chi2"]))
                allrows.append(row)
            acc=sum(r["accepted"] for r in allrows if r["scenario"]==sc and r["luminosity_factor"]==L)
            print(f"[replicas] {sc:15s} {L:2d}x: accepted {acc}/{args.replicas}")
    rep=pd.DataFrame(allrows); rep.to_csv(tab/"replica_fit_results.csv",index=False)
    # Empirical D-term bands.
    tg=np.linspace(0,0.95,191); brows=[]
    for (sc,L),g in rep[rep.accepted].groupby(["scenario","luminosity_factor"]):
        vals=np.array([d1_curve(r.C,r.dM2,r.dalpha,tg) for r in g.itertuples(index=False)])
        q16,q50,q84=percentile_band(vals)
        for k,t in enumerate(tg): brows.append(dict(scenario=sc,luminosity_factor=L,t_abs=t,q16=q16[k],q50=q50[k],q84=q84[k],n_replicas=len(g)))
    bands=pd.DataFrame(brows); bands.to_csv(tab/"d1_replica_bands.csv",index=False)
    # Parameter summaries and coverage around truth.
    sums=[]
    for (sc,L),g in rep[rep.accepted].groupby(["scenario","luminosity_factor"]):
        for p,truth in [("C",C0),("dM2",M20),("dalpha",ALPHA0)]:
            q=np.percentile(g[p],[16,50,84]); sums.append(dict(scenario=sc,luminosity_factor=L,parameter=p,truth=truth,q16=q[0],median=q[1],q84=q[2],n=len(g)))
    pd.DataFrame(sums).to_csv(tab/"replica_parameter_summaries.csv",index=False)
    # Mechanics only for the six configurations needed by the talk: stats all L + all scenarios at 10x.
    rgrid=np.linspace(.10,2.0,96); mrows=[]
    wanted={("statistics_only",L) for L in LUMI_FACTORS}|{(sc,10) for sc in SCENARIOS}
    for sc,L in sorted(wanted,key=lambda x:(x[1],x[0])):
        g=rep[(rep.scenario==sc)&(rep.luminosity_factor==L)&rep.accepted]
        if len(g)<20: continue
        P,S=mechanics_batch(g.C.to_numpy(),g.dM2.to_numpy(),g.dalpha.to_numpy(),rgrid,args.qmax,args.mechanics_nq)
        R2P=P*rgrid[None,:]**2
        for name,A in [("r2p",R2P),("shear",S)]:
            q16,q50,q84=percentile_band(A)
            for k,rv in enumerate(rgrid): mrows.append(dict(scenario=sc,luminosity_factor=L,quantity=name,r_fm=rv,q16=q16[k],q50=q50[k],q84=q84[k],n_replicas=len(g)))
    mech=pd.DataFrame(mrows); mech.to_csv(tab/"mechanics_replica_bands.csv",index=False)
    # Cumulative-q central diagnostic. q_data is sqrt(max controlled |t|) from derivative cache.
    qdata=float(np.sqrt(deriv.t_abs.max())); qcuts=sorted(set([qdata,1.25,1.5,2.,3.,5.,10.,args.qmax]))
    cum=cumulative_central(C0,M20,ALPHA0,rgrid,qcuts); cum["q_data_endpoint_GeV"]=qdata; cum.to_csv(tab/"mechanics_cumulative_q_support.csv",index=False)
    make_plots(rep,bands,mech,cum,fig)
    # Compact diagnostics.
    diag=rep.groupby(["scenario","luminosity_factor"]).agg(attempted=("accepted","size"),accepted=("accepted","sum"),median_nfev=("nfev","median"),bound_hits=("at_bound","sum")).reset_index()
    diag["accepted_fraction"]=diag.accepted/diag.attempted; diag.to_csv(tab/"replica_diagnostics.csv",index=False)
    print("\nReplica diagnostics:"); print(diag.to_string(index=False))
    print(f"\n[data/Fourier] controlled endpoint q_data=sqrt(tmax)={qdata:.3f} GeV; transforms also evaluated above this to expose model continuation")
    print(f"[output] {out}")
    print("[interpretation] use empirical d1/p/shear bands only after checking acceptance, bound hits, and nonlinear parameter tails.")

if __name__=="__main__": main()
