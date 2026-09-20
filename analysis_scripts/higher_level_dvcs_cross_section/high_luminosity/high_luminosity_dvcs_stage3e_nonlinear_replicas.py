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
from scipy.special import kv, gamma
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


def full_analytic_mechanics(C,M2,alpha,r_fm):
    """
    Analytic 3D Breit-frame transform for the generalized multipole.

    d1(-q^2) = -0.9 C (1 + q^2/M^2)^(-alpha)
    D_Q(t)   = (4/5) d1^Q(t)

    For the full transform,
      Dtilde(r) = D0 * M^3 /[(2pi)^(3/2) 2^(alpha-1) Gamma(alpha)]
                  * x^(alpha-3/2) K_(alpha-3/2)(x),  x=M r.

    Pressure and shear are obtained analytically from radial derivatives of
    Dtilde. This function has no q cutoff and no numerical Fourier integration:
    it is the full Nature/Volker-style generalized-multipole transform.
    """
    C=np.asarray(C,float); M2=np.asarray(M2,float); alpha=np.asarray(alpha,float)
    r_fm=np.asarray(r_fm,float)
    rr=r_fm/HBARC                         # fm -> GeV^-1
    M=np.sqrt(M2)
    D0=(4.0/5.0)*(-0.9*C)

    # Broadcast: replicas x radii.
    nu=alpha[:,None]-1.5
    x=M[:,None]*rr[None,:]
    pref=(D0*M**3/((2*np.pi)**1.5 * 2.0**(alpha-1.0) * gamma(alpha)))[:,None]

    # F(x)=x^nu K_nu(x)
    # F'=-x^nu K_(nu-1)
    # F''=x^nu K_(nu-2)-x^(nu-1)K_(nu-1)
    Fp = -np.power(x,nu)*kv(nu-1.0,x)
    Fpp= np.power(x,nu)*kv(nu-2.0,x) - np.power(x,nu-1.0)*kv(nu-1.0,x)

    # d/dr and d2/dr2 with r in GeV^-1.
    d1r=pref*M[:,None]*Fp
    d2r=pref*(M[:,None]**2)*Fpp

    # p = (D'' + 2D'/r)/(6 m_p)
    # s = -(D'' - D'/r)/(4 m_p)
    P_nat=(d2r + 2.0*d1r/rr[None,:])/(6.0*MP)
    S_nat=-(d2r - d1r/rr[None,:])/(4.0*MP)

    # GeV^4 -> GeV/fm^3.
    conv=1.0/(HBARC**3)
    return P_nat*conv, S_nat*conv


def finite_q_support_diagnostic(C,M2,alpha,r_fm,qmax,nq=5001):
    """Finite-q diagnostic only. A hard q cutoff is not a physical pressure distribution."""
    C=np.asarray(C,float); M2=np.asarray(M2,float); alpha=np.asarray(alpha,float)
    q=np.linspace(0,float(qmax),int(nq))
    if len(q)%2==0: q=q[:-1]
    rr=np.asarray(r_fm,float)/HBARC
    qr=np.outer(q,rr)
    j0=np.sinc(qr/np.pi)
    j2=spherical_jn(2,qr)
    w=simpson_weights(len(q),q[1]-q[0])
    d1=(-0.9*C[:,None])*(1+q[None,:]**2/M2[:,None])**(-alpha[:,None])
    base=d1*(w*q**4)[None,:]
    P=-(base@j0)/(15*np.pi**2*MP)
    S=-(base@j2)/(10*np.pi**2*MP)
    conv=1/(HBARC**3)
    return P*conv,S*conv


def robust_zero_crossings(r,y,relative_floor=1e-8):
    """Zero crossings after discarding the numerically negligible far tail."""
    r=np.asarray(r,float); y=np.asarray(y,float)
    floor=relative_floor*max(float(np.nanmax(np.abs(y))),1e-300)
    active=np.abs(y)>floor
    if not np.any(active):
        return []
    last=np.where(active)[0][-1]
    rr=r[:last+1]; yy=y[:last+1]
    cross=np.where(np.sign(yy[:-1])*np.sign(yy[1:])<0)[0]
    zeros=[]
    for i in cross:
        x1,x2=rr[i],rr[i+1]; y1,y2=yy[i],yy[i+1]
        zeros.append(float(x1-y1*(x2-x1)/(y2-y1)))
    return zeros


def validate_analytic_mechanics():
    """
    Cheap preflight run before any replicas.

    Checks:
      1) analytic transform normalization using the alpha=1 Yukawa limit;
      2) Volker/Nature central truth has one robust pressure crossing;
      3) central shear stays positive over the workshop plotting range.
    """
    # Transform normalization test.
    rtest_fm=np.array([0.2,0.5,1.0])
    rr=rtest_fm/HBARC
    M=1.0; alpha=1.0; D0=1.0
    nu=alpha-1.5; x=M*rr
    analytic=(D0*M**3/((2*np.pi)**1.5*2**(alpha-1)*gamma(alpha))
              *x**nu*kv(nu,x))
    yukawa=D0*M**2*np.exp(-M*rr)/(4*np.pi*rr)
    rel=np.max(np.abs((analytic-yukawa)/yukawa))

    # Central mechanical topology over the range actually shown in the talk.
    r=np.linspace(0.01,2.0,500)
    P,S=full_analytic_mechanics(np.array([C0]),np.array([M20]),
                                np.array([ALPHA0]),r)
    P=P[0]; S=S[0]
    zeros=robust_zero_crossings(r,P)
    shear_floor=-1e-10*max(float(np.max(np.abs(S))),1.0)
    shear_ok=bool(np.min(S)>=shear_floor)
    pressure_ok=bool(len(zeros)==1 and P[0]>0 and P[-1]<0)

    result=dict(
        yukawa_max_relative_difference=float(rel),
        pressure_zero_crossings=len(zeros),
        first_pressure_zero_fm=(zeros[0] if zeros else np.nan),
        pressure_naturelike=pressure_ok,
        shear_min_GeV_fm3=float(np.min(S)),
        shear_nonnegative=shear_ok,
    )
    ok=bool(rel<1e-10 and pressure_ok and shear_ok)
    return ok,result


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
    ap.add_argument("--replicas",type=int,default=100000)
    ap.add_argument("--workers",type=int,default=min(8,os.cpu_count() or 1))
    ap.add_argument("--seed",type=int,default=381902)
    ap.add_argument("--max-nfev",type=int,default=350)
    ap.add_argument("--qmax",type=float,default=40.0)
    ap.add_argument("--mechanics-nq",type=int,default=5001)
    ap.add_argument("--summary-sample",type=int,default=100000,
                    help="max replicas/config used for d1/mechanics percentile bands")
    ap.add_argument("--scatter-sample",type=int,default=20000)
    ap.add_argument("--save-replica-sample",type=int,default=0,
                    help="save at most N fitted replicas/config; default 0 = no raw replica table")
    args=ap.parse_args()

    s2=Path(args.stage2); s3=Path(args.stage3c)
    out=Path(args.outdir); tab=out/"tables"; fig=out/"figures"
    tab.mkdir(parents=True,exist_ok=True); fig.mkdir(parents=True,exist_ok=True)
    if not s2.exists(): raise FileNotFoundError(f"Missing Stage-2 input: {s2}")
    if not s3.exists(): raise FileNotFoundError(f"Missing Stage-3C input: {s3}")
    deriv=pd.read_csv(s3/"volker_observable_parameter_derivatives.csv")
    from gepard.fits import th_KM15
    central_imh={p:float(th_KM15.parameters[p]) for p in IMH_NAMES}

    print("="*96); print("STAGE 3E: STREAMED NONLINEAR REPLICA PROJECTION"); print("="*96)
    print(f"replicas/config={args.replicas:,}; workers={args.workers}; 12 configurations")
    print("disk policy: no per-replica CSV by default; only compact summaries and figures")

    ok,pre=validate_analytic_mechanics()
    pd.DataFrame([pre]).to_csv(tab/"mechanics_analytic_preflight.csv",index=False)
    print(f"[preflight] Yukawa rel.diff={pre['yukawa_max_relative_difference']:.3e}; "
          f"pressure crossings={pre['pressure_zero_crossings']}; "
          f"zero={pre['first_pressure_zero_fm']:.3f} fm; "
          f"topology={'PASS' if pre['pressure_naturelike'] else 'FAIL'}; "
          f"shear={'PASS' if pre['shear_nonnegative'] else 'FAIL'}")
    if not ok: raise RuntimeError("Analytic mechanics preflight failed.")
    for f in fig.glob("*.png"): f.unlink()

    tmin=float(deriv.t_abs.min()); tmax=float(deriv.t_abs.max())
    tg=np.linspace(tmin,tmax,181); qdata=float(np.sqrt(tmax)); rgrid=np.linspace(.05,2.0,99)
    rng=np.random.default_rng(args.seed)
    brows=[]; sums=[]; mrows=[]; diagrows=[]; plotrows=[]; saverows=[]
    wanted={("statistics_only",L) for L in LUMI_FACTORS}|{(sc,10) for sc in SCENARIOS}

    for L in LUMI_FACTORS:
        data=pd.read_csv(s2/f"joint_fit_input_km15_{L}x.csv"); data["bin"]=data["bin"].astype(int)
        for sc in SCENARIOS:
            cfg=prepare_config(data,deriv,sc,central_imh)
            seeds=rng.integers(1,2**63-1,size=args.replicas,dtype=np.int64)
            C=np.empty(args.replicas); M2=np.empty(args.replicas); A=np.empty(args.replicas)
            valid=np.zeros(args.replicas,dtype=bool); safe=np.zeros(args.replicas,dtype=bool)
            nfev=np.empty(args.replicas,dtype=np.int16)
            payload=((int(seed),cfg,args.max_nfev) for seed in seeds)
            if args.workers>1:
                ctx=mp.get_context("spawn")
                with ProcessPoolExecutor(max_workers=args.workers,mp_context=ctx) as ex:
                    it=ex.map(one_replica,payload,chunksize=max(1,min(5000,args.replicas//(args.workers*8))))
                    for i,v in enumerate(it):
                        _,success,chi2,nf,at_bound,c,m2,a,*_=v
                        good=bool(success and np.isfinite(chi2))
                        C[i]=c; M2[i]=m2; A[i]=a; valid[i]=good; safe[i]=good and not at_bound; nfev[i]=nf
            else:
                for i,x in enumerate(payload):
                    v=one_replica(x); _,success,chi2,nf,at_bound,c,m2,a,*_=v
                    good=bool(success and np.isfinite(chi2))
                    C[i]=c; M2[i]=m2; A[i]=a; valid[i]=good; safe[i]=good and not at_bound; nfev[i]=nf

            iv=np.flatnonzero(valid); isafe=np.flatnonzero(safe)
            diagrows.append(dict(scenario=sc,luminosity_factor=L,attempted=args.replicas,
                observable_valid=len(iv),extrapolation_safe=len(isafe),median_nfev=float(np.median(nfev)),
                bound_hits=len(iv)-len(isafe),observable_valid_fraction=len(iv)/args.replicas,
                extrapolation_safe_fraction=len(isafe)/args.replicas))
            print(f"[replicas] {sc:15s} {L:2d}x: valid {len(iv):,}/{args.replicas:,}; high-q-safe {len(isafe):,}")

            for name,arr,truth in [("C",C,C0),("dM2",M2,M20),("dalpha",A,ALPHA0)]:
                q=np.percentile(arr[iv],[16,50,84])
                sums.append(dict(scenario=sc,luminosity_factor=L,parameter=name,truth=truth,
                    q16=q[0],median=q[1],q84=q[2],n=len(iv),bound_hit_fraction=1-len(isafe)/len(iv)))

            take=min(args.summary_sample,len(iv))
            sel=iv if take==len(iv) else rng.choice(iv,take,replace=False)
            V=d1_curve(C[sel,None],M2[sel,None],A[sel,None],tg[None,:])
            q16,q50,q84=percentile_band(V)
            for k,t in enumerate(tg):
                brows.append(dict(scenario=sc,luminosity_factor=L,t_abs=t,q16=q16[k],q50=q50[k],q84=q84[k],
                    n_replicas=take,controlled_tmin=tmin,controlled_tmax=tmax))
            del V

            nplot=min(args.scatter_sample,len(iv)); psel=iv if nplot==len(iv) else rng.choice(iv,nplot,replace=False)
            plotrows.extend(dict(scenario=sc,luminosity_factor=L,C=C[j],dM2=M2[j],dalpha=A[j],
                                 observable_valid=True,at_bound=not safe[j]) for j in psel)
            if args.save_replica_sample:
                nk=min(args.save_replica_sample,len(iv)); ksel=iv if nk==len(iv) else rng.choice(iv,nk,replace=False)
                saverows.extend(dict(scenario=sc,luminosity_factor=L,C=C[j],dM2=M2[j],dalpha=A[j],
                                     extrapolation_safe=bool(safe[j])) for j in ksel)

            if (sc,L) in wanted and len(isafe)>=20:
                nm=min(args.summary_sample,len(isafe)); msel=isafe if nm==len(isafe) else rng.choice(isafe,nm,replace=False)
                P,S=full_analytic_mechanics(C[msel],M2[msel],A[msel],rgrid)
                for name,X in [("r2p_full",P*rgrid[None,:]**2),("r2s_full",S*rgrid[None,:]**2),
                               ("pressure_full",P),("shear_full",S)]:
                    q16,q50,q84=percentile_band(X)
                    for k,r in enumerate(rgrid):
                        mrows.append(dict(scenario=sc,luminosity_factor=L,quantity=name,r_fm=r,
                            q16=q16[k],q50=q50[k],q84=q84[k],n_replicas=nm,qmax_GeV=np.inf,
                            interpretation="full generalized-multipole projection"))
                pok=[]; sok=[]
                for pp,ss in zip(P,S):
                    z=robust_zero_crossings(rgrid,pp)
                    pok.append(pp[0]>0 and len(z)==1 and pp[-1]<=0)
                    sok.append(np.nanmin(ss)>=-1e-10*max(float(np.nanmax(np.abs(ss))),1.0))
                mrows.append(dict(scenario=sc,luminosity_factor=L,quantity="stability_summary",
                    r_fm=np.nan,q16=np.nan,q50=np.nan,q84=np.nan,n_replicas=nm,qmax_GeV=np.inf,
                    interpretation="mechanical-shape diagnostic",
                    pressure_naturelike_fraction=float(np.mean(pok)),shear_positive_fraction=float(np.mean(sok))))
                del P,S
            del C,M2,A,valid,safe,nfev

    bands=pd.DataFrame(brows); pars=pd.DataFrame(sums); mech=pd.DataFrame(mrows); diag=pd.DataFrame(diagrows)
    bands.to_csv(tab/"d1_replica_bands.csv",index=False)
    pars.to_csv(tab/"replica_parameter_summaries.csv",index=False)
    mech.to_csv(tab/"mechanics_replica_bands.csv",index=False)
    diag.to_csv(tab/"replica_diagnostics.csv",index=False)
    if args.save_replica_sample: pd.DataFrame(saverows).to_csv(tab/"replica_fit_sample.csv",index=False)

    rcheck=np.linspace(.01,8.0,800)
    P0,S0=full_analytic_mechanics(np.array([C0]),np.array([M20]),np.array([ALPHA0]),rcheck); P0=P0[0]; S0=S0[0]
    zeros=robust_zero_crossings(rcheck,P0)
    pd.DataFrame([dict(pressure_zero_crossings=len(zeros),first_pressure_zero_fm=zeros[0] if zeros else np.nan,
        shear_min_GeV_fm3=float(np.min(S0)),shear_nonnegative=bool(np.min(S0)>=-1e-10*max(float(np.max(np.abs(S0))),1.0)),
        von_laue_integral_GeV=float(np.trapezoid(rcheck**2*P0,rcheck)),qmax_GeV=np.inf,
        note="analytic full generalized-multipole central truth")]).to_csv(tab/"mechanics_central_stability_check.csv",index=False)

    qcuts=sorted(set([qdata,2.0,args.qmax])); cum=cumulative_central(C0,M20,ALPHA0,rgrid,qcuts)
    cum["q_data_endpoint_GeV"]=qdata; cum.to_csv(tab/"mechanics_cumulative_q_support.csv",index=False)
    make_plots(pd.DataFrame(plotrows),bands,mech,cum,fig,qdata)
    print("\nReplica diagnostics:"); print(diag.to_string(index=False))
    print(f"\n[output] {out}")
    print("[disk] no raw per-replica table written." if not args.save_replica_sample else
          f"[disk] compact sample only: <= {args.save_replica_sample:,} replicas/config.")


if __name__=="__main__": main()
