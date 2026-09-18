#!/usr/bin/env python3
"""
Stage 3D: propagate the Stage-3C generalized-multipole quark D-term into
radial pressure and shear-force distributions.

Conventions follow Burkert, Elouadrhiri & Girod, Nature 557, 396 (2018):
  d1^Q(t) = d1^Q(0) (1 - t/M^2)^(-alpha)
  D_Q(t)  = (4/5) d1^Q(t)

For a spherically symmetric Breit-frame transform,
  Dtilde(r) = int d^3q/(2pi)^3 exp(i q.r) D_Q(-q^2)
  p(r) = [1/(6 m_p)] (1/r^2) d/dr [r^2 d Dtilde/dr]
  s(r) = -[1/(4 m_p)] r d/dr [(1/r) d Dtilde/dr]

Equivalent radial integrals used here are
  p(r) = -1/(15 pi^2 m_p) int dq q^4 j0(qr) d1^Q(-q^2)
  s(r) = -1/(10 pi^2 m_p) int dq q^4 j2(qr) d1^Q(-q^2).

The Stage-3C full (C, M^2, alpha) covariance is propagated coherently.
This is a *linearized* uncertainty projection.  Because M^2 and alpha are
very strongly correlated and individually weakly constrained, final talk-level
mechanics bands should be validated with nonlinear pseudo-experiment/refit
replicas (planned next), not interpreted as a final confidence interval.
"""
from pathlib import Path
import argparse, math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import quad_vec
from scipy.special import spherical_jn

LUMI_FACTORS=(1,2,5,10)
SCENARIOS=("statistics_only","ptp_half","baseline")
MP=0.9382720813       # GeV
HBARC=0.1973269804    # GeV fm
GEV4_TO_GEV_FM3=1.0/HBARC**3
GEV_FM3_TO_PA=1.602176634e35


def load_covariance(df, scenario, lumi):
    names=("C","dM2","dalpha")
    sub=df[(df.scenario==scenario)&(df.luminosity_factor==lumi)]
    cov=np.zeros((3,3))
    for i,a in enumerate(names):
        for j,b in enumerate(names):
            z=sub[(sub.parameter_i==a)&(sub.parameter_j==b)]
            if len(z)!=1: raise RuntimeError(f"Missing covariance {scenario} {lumi}x {a},{b}")
            cov[i,j]=float(z.iloc[0].covariance)
    return cov


def d1_q(q, C, M2, alpha):
    # Stage 3C convention: C_H=-C(...), d1^Q=0.9 C_H.
    return -0.9*C*(1.0+q*q/M2)**(-alpha)


def mechanics_profile(C,M2,alpha,r_fm,qmax=80.0):
    r=np.asarray(r_fm)/HBARC  # GeV^-1
    def fp(q): return q**4*d1_q(q,C,M2,alpha)*spherical_jn(0,q*r)
    def fs(q): return q**4*d1_q(q,C,M2,alpha)*spherical_jn(2,q*r)
    Ip=quad_vec(fp,0.0,qmax,epsabs=1e-8,epsrel=2e-5,limit=800)[0]
    Is=quad_vec(fs,0.0,qmax,epsabs=1e-8,epsrel=2e-5,limit=800)[0]
    p_nat=-Ip/(15.0*np.pi**2*MP)
    s_nat=-Is/(10.0*np.pi**2*MP)
    return p_nat*GEV4_TO_GEV_FM3, s_nat*GEV4_TO_GEV_FM3


def profile_and_gradient(C,M2,alpha,r_fm,qmax):
    p,s=mechanics_profile(C,M2,alpha,r_fm,qmax)
    vals=np.array([C,M2,alpha],float)
    fracs=np.array([1e-3,1e-3,1e-3])
    gp=np.empty((len(r_fm),3)); gs=np.empty_like(gp)
    for k in range(3):
        h=fracs[k]*max(abs(vals[k]),1.0)
        vp=vals.copy(); vm=vals.copy(); vp[k]+=h; vm[k]-=h
        pp,sp=mechanics_profile(*vp,r_fm,qmax)
        pm,sm=mechanics_profile(*vm,r_fm,qmax)
        gp[:,k]=(pp-pm)/(2*h); gs[:,k]=(sp-sm)/(2*h)
    return p,s,gp,gs


def propagate(p,s,gp,gs,cov):
    vp=np.einsum('ri,ij,rj->r',gp,cov,gp)
    vs=np.einsum('ri,ij,rj->r',gs,cov,gs)
    return np.sqrt(np.maximum(vp,0)),np.sqrt(np.maximum(vs,0))


def savefig(fig,path):
    fig.tight_layout(); fig.savefig(path,dpi=250,bbox_inches='tight'); plt.close(fig)


def make_plots(df,figdir):
    labels={"baseline":"Current point-to-point systematics",
            "ptp_half":"Point-to-point systematics / 2",
            "statistics_only":"Statistics only"}
    # Statistical progression: pressure, shown as r^2 p(r) like Nature 2018.
    fig,ax=plt.subplots(figsize=(8.8,6.0))
    for L in LUMI_FACTORS:
        d=df[(df.scenario=='statistics_only')&(df.luminosity_factor==L)]
        y=d.r_fm**2*d.pressure_GeV_fm3; e=d.r_fm**2*d.sigma_pressure_GeV_fm3
        ax.fill_between(d.r_fm,y-e,y+e,alpha=.15); ax.plot(d.r_fm,y,label=f'{L}x')
    ax.axhline(0,lw=.8,alpha=.5); ax.set_xlabel('r (fm)'); ax.set_ylabel(r'$r^2p(r)$ (GeV fm$^{-1}$)')
    ax.set_title('Statistical potential for the proton pressure distribution'); ax.grid(alpha=.2); ax.legend(title='Pass-2 exposure')
    savefig(fig,figdir/'01_pressure_statistics_only_luminosity.png')

    fig,ax=plt.subplots(figsize=(8.8,6.0))
    for L in LUMI_FACTORS:
        d=df[(df.scenario=='statistics_only')&(df.luminosity_factor==L)]
        ax.fill_between(d.r_fm,d.shear_GeV_fm3-d.sigma_shear_GeV_fm3,d.shear_GeV_fm3+d.sigma_shear_GeV_fm3,alpha=.15)
        ax.plot(d.r_fm,d.shear_GeV_fm3,label=f'{L}x')
    ax.axhline(0,lw=.8,alpha=.5); ax.set_xlabel('r (fm)'); ax.set_ylabel(r'$s(r)$ (GeV fm$^{-3}$)')
    ax.set_title('Statistical potential for the proton shear-force distribution'); ax.grid(alpha=.2); ax.legend(title='Pass-2 exposure')
    savefig(fig,figdir/'02_shear_statistics_only_luminosity.png')

    # Systematics at 10x.
    fig,ax=plt.subplots(figsize=(8.8,6.0))
    for sc in SCENARIOS:
        d=df[(df.scenario==sc)&(df.luminosity_factor==10)]
        y=d.r_fm**2*d.pressure_GeV_fm3; e=d.r_fm**2*d.sigma_pressure_GeV_fm3
        ax.fill_between(d.r_fm,y-e,y+e,alpha=.15); ax.plot(d.r_fm,y,label=labels[sc])
    ax.axhline(0,lw=.8,alpha=.5); ax.set_xlabel('r (fm)'); ax.set_ylabel(r'$r^2p(r)$ (GeV fm$^{-1}$)')
    ax.set_title('Systematic limitation of the 10x pressure projection'); ax.grid(alpha=.2); ax.legend()
    savefig(fig,figdir/'03_pressure_systematics_10x.png')

    fig,ax=plt.subplots(figsize=(8.8,6.0))
    for sc in SCENARIOS:
        d=df[(df.scenario==sc)&(df.luminosity_factor==10)]
        ax.fill_between(d.r_fm,d.shear_GeV_fm3-d.sigma_shear_GeV_fm3,d.shear_GeV_fm3+d.sigma_shear_GeV_fm3,alpha=.15)
        ax.plot(d.r_fm,d.shear_GeV_fm3,label=labels[sc])
    ax.axhline(0,lw=.8,alpha=.5); ax.set_xlabel('r (fm)'); ax.set_ylabel(r'$s(r)$ (GeV fm$^{-3}$)')
    ax.set_title('Systematic limitation of the 10x shear projection'); ax.grid(alpha=.2); ax.legend()
    savefig(fig,figdir/'04_shear_systematics_10x.png')


def main():
    here=Path(__file__).resolve().parent
    ap=argparse.ArgumentParser()
    ap.add_argument('--stage3c',default=str(here/'output'/'stage3c'/'tables'))
    ap.add_argument('--outdir',default=str(here/'output'/'stage3d'))
    ap.add_argument('--qmax',type=float,default=80.0)
    args=ap.parse_args()
    inp=Path(args.stage3c); out=Path(args.outdir); tab=out/'tables'; fig=out/'figures'; tab.mkdir(parents=True,exist_ok=True); fig.mkdir(parents=True,exist_ok=True)
    pars=pd.read_csv(inp/'dr_parameter_uncertainties.csv'); covdf=pd.read_csv(inp/'dterm_parameter_covariances.csv')
    # Avoid r=0: short-distance transform is most sensitive to unmeasured large |t|.
    r=np.linspace(0.10,2.00,96)
    C=float(pars.C.iloc[0]); M2=float(pars.dM2.iloc[0]); alpha=float(pars.dalpha.iloc[0])
    p,s,gp,gs=profile_and_gradient(C,M2,alpha,r,args.qmax)
    rows=[]
    for L in LUMI_FACTORS:
        for sc in SCENARIOS:
            cov=load_covariance(covdf,sc,L); ep,es=propagate(p,s,gp,gs,cov)
            for k,rr in enumerate(r):
                rows.append({'scenario':sc,'luminosity_factor':L,'r_fm':rr,
                             'pressure_GeV_fm3':p[k],'sigma_pressure_GeV_fm3':ep[k],
                             'pressure_Pa':p[k]*GEV_FM3_TO_PA,'sigma_pressure_Pa':ep[k]*GEV_FM3_TO_PA,
                             'shear_GeV_fm3':s[k],'sigma_shear_GeV_fm3':es[k]})
    df=pd.DataFrame(rows); df.to_csv(tab/'pressure_shear_bands_linearized.csv',index=False)
    # Fourier-tail convergence diagnostic: central profiles at qmax/2 and qmax.
    p2,s2=mechanics_profile(C,M2,alpha,r,args.qmax/2)
    conv=pd.DataFrame({'r_fm':r,'pressure_qmax':p,'pressure_qmax_half':p2,
                       'pressure_rel_change':np.abs(p-p2)/np.maximum(np.abs(p),1e-12),
                       'shear_qmax':s,'shear_qmax_half':s2,
                       'shear_rel_change':np.abs(s-s2)/np.maximum(np.abs(s),1e-12)})
    conv.to_csv(tab/'fourier_qmax_convergence.csv',index=False)
    make_plots(df,fig)
    print('='*96); print('STAGE 3D: D-TERM -> PRESSURE AND SHEAR')
    print(f'input  : {inp}'); print(f'output : {out}')
    print('relation: D_Q(t)=4/5 d1^Q(t); 3D Breit-frame Fourier transform')
    print('warning : bands are linear covariance propagation; validate with nonlinear refit replicas before final use')
    print(f'qmax convergence diagnostic written for {args.qmax/2:g} vs {args.qmax:g} GeV')

if __name__=='__main__': main()
