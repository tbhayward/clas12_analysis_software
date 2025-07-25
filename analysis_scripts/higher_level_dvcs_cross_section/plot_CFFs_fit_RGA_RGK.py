#!/usr/bin/env python3
"""
plot_CFFs_fit_RGA_RGK.py

Usage:
    python plot_CFFs_fit_RGA_RGK.py <fit_results1.txt> [<fit_results2.txt>]

Reads which CFFs were fit from the header of each results file, then for each
enabled Im CFF makes two sets of figures:
  1) “single” — Im CFF vs. ξ and vs. −t for the first fit only
  2) “comparison” — same plots overlaying fit1 (blue) and fit2 (red)

Always shows median curve and 95% CI bands. Y-axes run from 0 to 20 with ticks
at [0,5,10,15,20], but the top‐row plots omit the “0” tick so it doesn’t
overlap the bottom‐row “20”. And the bottom middle/right panels drop the “0.0”
x-tick so it doesn’t clip the neighbor’s label.
"""
import os, sys, re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# ─── Command‐line ────────────────────────────────────────────────────────────────
if not (2 < len(sys.argv) <= 3):
    print("Usage: python plot_CFFs_fit_RGA_RGK.py <fit1.txt> [<fit2.txt>]")
    sys.exit(1)
fit1file = sys.argv[1]
fit2file = sys.argv[2] if len(sys.argv) == 3 else None

def parse_timestamp(fn):
    m = re.search(r'fit_results_(\d{8}_\d{6})\.txt$', fn)
    return m.group(1) if m else os.path.splitext(fn)[0]

ts1 = parse_timestamp(fit1file)
ts2 = parse_timestamp(fit2file) if fit2file else None

# ─── Parse results ──────────────────────────────────────────────────────────────
def parse_fit_results(fname):
    with open(fname) as f:
        L = [l.strip() for l in f if l.strip()]
    fl = next(l for l in L if l.startswith("H "))
    toks = fl.split()
    flags = { toks[i]: int(toks[i+1]) for i in range(0,len(toks),2) }
    pnames=[]
    for i,l in enumerate(L):
        if l.startswith("# parameters"):
            pnames=L[i].split()[2:]
            break
    vals=errs=None
    chi2=ndf=None
    for i,l in enumerate(L):
        if l.startswith("# values"): vals=list(map(float,L[i+1].split()))
        if l.startswith("# errors"): errs=list(map(float,L[i+1].split()))
        if l.startswith("# chi2"):
            a=L[i+1].split(); chi2,ndf=float(a[0]),int(a[1])
    return flags,pnames,np.array(vals),np.array(errs),chi2,ndf

flags1,pn1,vals1,errs1,chi21,ndf1 = parse_fit_results(fit1file)
if fit2file:
    flags2,pn2,vals2,errs2,chi22,ndf2 = parse_fit_results(fit2file)

def get_idx(pn,name): return pn.index(name)
def build_dicts(flags,pn,vals,errs):
    ren=vals[get_idx(pn,"renormImag")]
    renerr=errs[get_idx(pn,"renormImag")]
    fp,fe={},{}
    for cff in ("H","Ht","E","Et"):
        if flags[cff]:
            ks=["r","alpha0","alpha1","n","b","M2","P"]
            idxs={k:get_idx(pn,f"{k}_{cff}") for k in ks}
            fp[cff]={k:vals[i] for k,i in idxs.items()}
            fe[cff]={k:errs[i] for k,i in idxs.items()}
    return ren,renerr,fp,fe

ren1,ren1err,fitp1,fite1 = build_dicts(flags1,pn1,vals1,errs1)
if fit2file:
    ren2,ren2err,fitp2,fite2 = build_dicts(flags2,pn2,vals2,errs2)

defaults={
 "H":  dict(r=0.9,alpha0=0.43,alpha1=0.85,n=1.35,b=0.4, M2=0.64,P=1.0,factor=2.0),
 "Ht": dict(r=7.0,alpha0=0.43,alpha1=0.85,n=0.6, b=2.0, M2=0.8, P=1.0,factor=0.4),
 "E":  dict(r=0.9,alpha0=0.43,alpha1=0.85,n=1.35,b=0.4,M2=0.64,P=1.0,factor=1.0),
 "Et": dict(r=0.9,alpha0=0.43,alpha1=0.85,n=0.0, b=0.4, M2=0.64,P=1.0,factor=1.0),
}

def make_Im_func(cff,params,renorm):
    d=defaults[cff]
    def Im(xi,t):
        r    = params.get("r",d["r"])
        a0   = params.get("alpha0",d["alpha0"])
        a1   = params.get("alpha1",d["alpha1"])
        nval = params.get("n",d["n"])
        bval = params.get("b",d["b"])
        M2   = params.get("M2",d["M2"])
        Pval = params.get("P",d["P"])
        fac  = d["factor"]
        alpha= a0 + a1*t
        pref = np.pi*5/9 * nval*r/(1+xi)
        xfac = (2*xi/(1+xi))**(-alpha)
        yfac = ((1-xi)/(1+xi))**(bval)
        tfac = (1 - ((1-xi)/(1+xi))*t/M2)**(-Pval)
        return renorm*pref*xfac*yfac*tfac*fac
    return Im

def generate_reps(cent,err,n=200):
    out=[]
    for _ in range(n):
        d={}
        for k,v in cent.items():
            d[k]=np.random.normal(v,err[k]/1.96)
        out.append(d)
    return out

def band(cff,xi_vals,t_vals,ren,cent,err,n=200):
    xi=np.atleast_1d(xi_vals)
    t =np.atleast_1d(t_vals)
    xi_b,t_b = np.broadcast_arrays(xi,-t)
    reps=generate_reps(cent[cff],err[cff],n)
    ren_reps=np.random.normal(ren,ren/np.sqrt(len(xi_b)),n)
    arr=[]
    for i in range(n):
        f=make_Im_func(cff,reps[i],ren_reps[i])
        arr.append(f(xi_b,t_b))
    A=np.array(arr)
    return (np.median(A,0), np.percentile(A,2.5,0), np.percentile(A,97.5,0))

# ─── Plot setup ────────────────────────────────────────────────────────────────
os.makedirs("output/plots",exist_ok=True)
t_fixed=np.linspace(0.1,1.0,6)
xi_fixed=np.linspace(0.05,0.50,6)
xi_range=np.linspace(0.0,0.5,200)
t_range =np.linspace(0.0,1.0,200)

plt.style.use('classic')
plt.rcParams.update({'font.size':14,'font.family':'serif'})

style1={'color':'tab:red','ls':'--','lw':2.5}
band1 ={'color':'tab:red','alpha':0.2}
style2={'color':'tab:blue','ls':'-','lw':2.5}
band2 ={'color':'tab:blue','alpha':0.2}
zero  ={'color':'gray','ls':'--','lw':1}

lbl1="RGK preliminary"
lbl2="RGA pass-1"
tex_map={"H":"H","Ht":r"\tilde H","E":"E","Et":r"\tilde E"}

# ─── Single (RGK only) ──────────────────────────────────────────────────────────
for cff in fitp1:
    tex=tex_map[cff]
    # Im vs ξ
    fig,axes=plt.subplots(2,3,figsize=(12,8),sharex=True,sharey=True)
    axes=axes.flatten()
    fig.suptitle(rf"$\mathrm{{Im}}\,{tex}$",fontsize=16,y=0.98)
    for ax,t in zip(axes,t_fixed):
        m,lo,up=band(cff,xi_range,t,ren1,fitp1,fite1)
        ax.plot(xi_range,m,**style1)
        ax.fill_between(xi_range,lo,up,**band1)
        ax.axhline(0,**zero)
        ax.text(0.6,0.75,rf"$-t={t:.2f}$",transform=ax.transAxes)
        ax.set_xlim(0,0.5); ax.set_ylim(0,20); ax.set_yticks([0,5,10,15,20])
    for ax in axes[:3]:
        ax.set_yticks([5,10,15,20])
    # drop 0.0 x‐tick on bottom middle/right
    for ax in (axes[4],axes[5]):
        for lbl in ax.get_xticklabels():
            if lbl.get_text() in ('0','0.0'): lbl.set_visible(False)
    axes[2].legend([Line2D([0],[0],**style1,label=lbl1)],loc='upper right',fontsize=10)
    fig.text(0.06,0.5,rf"$\mathrm{{Im}}\,{tex}(\xi,\,-t)$",va='center',ha='center',rotation='vertical')
    fig.text(0.5,0.02,r"$\xi$",ha='center')
    fig.subplots_adjust(left=0.10,right=0.98,top=0.92,bottom=0.08,wspace=0,hspace=0)
    out=f"output/plots/RGK_Im{cff}_vs_xi_{ts1}.pdf"
    fig.savefig(out,bbox_inches='tight'); print("Saved:",out)
    plt.close(fig)

    # Im vs −t
    fig,axes=plt.subplots(2,3,figsize=(12,8),sharex=True,sharey=True)
    axes=axes.flatten()
    fig.suptitle(rf"$\mathrm{{Im}}\,{tex}$",fontsize=16,y=0.98)
    for ax,xi0 in zip(axes,xi_fixed):
        m,lo,up=band(cff,xi0,t_range,ren1,fitp1,fite1)
        ax.plot(t_range,m,**style1)
        ax.fill_between(t_range,lo,up,**band1)
        ax.axhline(0,**zero)
        ax.text(0.6,0.75,rf"$\xi={xi0:.2f}$",transform=ax.transAxes)
        ax.set_xlim(0,1.0); ax.set_ylim(0,20); ax.set_yticks([0,5,10,15,20])
    for ax in axes[:3]:
        ax.set_yticks([5,10,15,20])
    for ax in (axes[4],axes[5]):
        for lbl in ax.get_xticklabels():
            if lbl.get_text() in ('0','0.0'): lbl.set_visible(False)
    axes[2].legend([Line2D([0],[0],**style1,label=lbl1)],loc='upper right',fontsize=10)
    fig.text(0.06,0.5,rf"$\mathrm{{Im}}\,{tex}(\xi,\,-t)$",va='center',ha='center',rotation='vertical')
    fig.text(0.5,0.02,r"$-t\,(\mathrm{GeV^2})$",ha='center')
    fig.subplots_adjust(left=0.10,right=0.98,top=0.92,bottom=0.08,wspace=0,hspace=0)
    out=f"output/plots/RGK_Im{cff}_vs_t_{ts1}.pdf"
    fig.savefig(out,bbox_inches='tight'); print("Saved:",out)
    plt.close(fig)

# ─── Comparison ─────────────────────────────────────────────────────────────────
if fit2file:
    for cff in fitp1:
        tex=tex_map[cff]
        # Im vs ξ
        fig,axes=plt.subplots(2,3,figsize=(12,8),sharex=True,sharey=True)
        axes=axes.flatten()
        fig.suptitle(rf"$\mathrm{{Im}}\,{tex}$",fontsize=16,y=0.98)
        for ax,t in zip(axes,t_fixed):
            m1,l1,u1=band(cff,xi_range,t,ren1,fitp1,fite1)
            m2,l2,u2=band(cff,xi_range,t,ren2,fitp2,fite2)
            ax.plot(xi_range,m2,**style2)
            ax.fill_between(xi_range,l2,u2,**band2)
            ax.plot(xi_range,m1,**style1)
            ax.fill_between(xi_range,l1,u1,**band1)
            ax.axhline(0,**zero)
            ax.text(0.6,0.75,rf"$-t={t:.2f}$",transform=ax.transAxes)
            ax.set_xlim(0,0.5); ax.set_ylim(0,20); ax.set_yticks([0,5,10,15,20])
        for ax in axes[:3]:
            ax.set_yticks([5,10,15,20])
        for ax in (axes[4],axes[5]):
            for lbl in ax.get_xticklabels():
                if lbl.get_text() in ('0','0.0'): lbl.set_visible(False)
        handles=[
          Line2D([0],[0],**style2,label=lbl2),
          Line2D([0],[0],**style1,label=lbl1),
          Patch(**band2,label="RGA 95% CI"),
          Patch(**band1,label="RGK 95% CI")
        ]
        axes[2].legend(handles=handles,loc='upper right',fontsize=10)
        fig.text(0.06,0.5,rf"$\mathrm{{Im}}\,{tex}(\xi,\,-t)$",va='center',ha='center',rotation='vertical')
        fig.text(0.5,0.02,r"$\xi$",ha='center')
        fig.subplots_adjust(left=0.10,right=0.98,top=0.92,bottom=0.08,wspace=0,hspace=0)
        out=f"output/plots/RGA_vs_RGK_Im{cff}_vs_xi_{ts2}_{ts1}.pdf"
        fig.savefig(out,bbox_inches='tight'); print("Saved:",out)
        plt.close(fig)

        # Im vs −t
        fig,axes=plt.subplots(2,3,figsize=(12,8),sharex=True,sharey=True)
        axes=axes.flatten()
        fig.suptitle(rf"$\mathrm{{Im}}\,{tex}$",fontsize=16,y=0.98)
        for ax,xi0 in zip(axes,xi_fixed):
            m1,l1,u1=band(cff,xi0,t_range,ren1,fitp1,fite1)
            m2,l2,u2=band(cff,xi0,t_range,ren2,fitp2,fite2)
            ax.plot(t_range,m2,**style2)
            ax.fill_between(t_range,l2,u2,**band2)
            ax.plot(t_range,m1,**style1)
            ax.fill_between(t_range,l1,u1,**band1)
            ax.axhline(0,**zero)
            ax.text(0.6,0.75,rf"$\xi={xi0:.2f}$",transform=ax.transAxes)
            ax.set_xlim(0,1.0); ax.set_ylim(0,20); ax.set_yticks([0,5,10,15,20])
        for ax in axes[:3]:
            ax.set_yticks([5,10,15,20])
        for ax in (axes[4],axes[5]):
            for lbl in ax.get_xticklabels():
                if lbl.get_text() in ('0','0.0'): lbl.set_visible(False)
        axes[2].legend(handles=handles,loc='upper right',fontsize=10)
        fig.text(0.06,0.5,rf"$\mathrm{{Im}}\,{tex}(\xi,\,-t)$",va='center',ha='center',rotation='vertical')
        fig.text(0.5,0.02,r"$-t\,(\mathrm{GeV^2})$",ha='center')
        fig.subplots_adjust(left=0.10,right=0.98,top=0.92,bottom=0.08,wspace=0,hspace=0)
        out=f"output/plots/RGA_vs_RGK_Im{cff}_vs_t_{ts2}_{ts1}.pdf"
        fig.savefig(out,bbox_inches='tight'); print("Saved:",out)
        plt.close(fig)

print("All done.")