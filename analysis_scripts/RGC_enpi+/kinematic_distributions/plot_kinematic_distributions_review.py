#!/usr/bin/env python3
"""Publication-style RGC e n pi+ kinematic-distribution plots for reviewer response.

Reads nominal NH3 paper-version ROOT trees for Su22/Fa22/Sp23 and applies the
period- and (xB,-t')-dependent nominal +/-2 sigma Mx2 windows from channel
selection. No asymmetry-extraction quantities are fitted or modified.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import uproot
from matplotlib.colors import LogNorm

PERIODS = ("su22", "fa22", "sp23")
XB_BINS = ((0.10,0.25),(0.25,0.35),(0.35,0.45),(0.45,0.60))
TP_BINS = ((0.05,0.25),(0.25,0.45),(0.45,0.65),(0.65,0.85),(0.85,1.05),(1.05,1.25))
PAPER_DIR = Path("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/paper_versions")
DEFAULT_INPUTS = {p: PAPER_DIR / f"rgc_{p}_inb_NH3_epi+_mom_corrections.root" for p in PERIODS}
DEFAULT_CUTS = Path("../channel_selection/output/channel_selection_mx2_fit_stability/final_carbon_assisted_cuts/tables/final_carbon_assisted_mx2_cuts.json")
DEFAULT_OUT = Path("output/kinematic_distributions_review")
TREE, CHUNK = "PhysicsEvents", "250 MB"
CMAP = "turbo"  # blue -> cyan/green -> yellow -> red; requested replacement for viridis
ALIASES = {
    "xB": ("xB","x","xb","x_b"), "Q2": ("Q2","q2"),
    "tprime": ("tprime","t_prime","tp","tPrime"),
    "Mx2": ("Mx2","mx2","Mx2_epi","Mx2_epip","missing_mass_squared","missing_mass2"),
    "phi": ("phi","phi1","phi_h","trento_phi"), "W": ("W","w"), "y": ("y","inelasticity")
}

def resolve(tree, logical, required=True):
    names=set(tree.keys())
    for c in ALIASES[logical]:
        if c in names: return c
    if required: raise KeyError(f"Could not resolve {logical}; tried {ALIASES[logical]}")
    return None

def phi_degrees(values, branch):
    values=np.asarray(values,dtype=float); finite=values[np.isfinite(values)]
    if not finite.size: return values
    if np.nanmax(np.abs(finite)) <= 2*np.pi+0.25: values=np.degrees(values)
    return np.mod(values,360.0)

def bin_index(xb, mtp):
    out=np.full(xb.shape,-1,dtype=np.int16)
    for ix,(xl,xh) in enumerate(XB_BINS):
        for it,(tl,th) in enumerate(TP_BINS):
            out[(xb>=xl)&(xb<xh)&(mtp>=tl)&(mtp<th)] = ix*len(TP_BINS)+it+1
    return out

def load_nominal_windows(path):
    payload=json.loads(path.read_text())
    if float(payload.get("nominal_sigma_multiple",np.nan)) != 2.0:
        raise RuntimeError(f"Expected nominal 2-sigma windows in {path}")
    windows={}
    for p in PERIODS:
        windows[p]={int(r["bin_number"]):tuple(map(float,r["nominal"])) for r in payload["periods"][p]}
        if len(windows[p]) != 24: raise RuntimeError(f"Expected 24 nominal windows for {p}")
    return windows

def collect(period,path,windows):
    print(f"[load] {period}: {path}",flush=True)
    with uproot.open(path) as f:
        tree=f[TREE]
        branches={k:resolve(tree,k,required=(k not in {"W","y"})) for k in ALIASES}
        expressions=[v for v in branches.values() if v is not None]
        acc={k:[] for k in ("xB","Q2","minus_tprime","Mx2","phi","W","y")}
        seen=selected=0
        for ichunk,arrays in enumerate(tree.iterate(expressions=expressions,step_size=CHUNK,library="np"),1):
            xb=np.asarray(arrays[branches["xB"]],float); q2=np.asarray(arrays[branches["Q2"]],float)
            mtp=-np.asarray(arrays[branches["tprime"]],float); mx2=np.asarray(arrays[branches["Mx2"]],float)
            phi=phi_degrees(arrays[branches["phi"]],branches["phi"])
            w=np.asarray(arrays[branches["W"]],float) if branches["W"] else np.full(xb.shape,np.nan)
            y=np.asarray(arrays[branches["y"]],float) if branches["y"] else np.full(xb.shape,np.nan)
            b=bin_index(xb,mtp)
            keep=(b>0)&np.isfinite(xb)&np.isfinite(q2)&np.isfinite(mtp)&np.isfinite(mx2)&np.isfinite(phi)
            for bn in range(1,25):
                lo,hi=windows[bn]; m=b==bn; keep[m]&=(mx2[m]>=lo)&(mx2[m]<=hi)
            keep &= q2>1.0
            if branches["W"]: keep &= w>2.0
            if branches["y"]: keep &= y<0.80
            seen += xb.size; selected += int(np.count_nonzero(keep))
            vals={"xB":xb,"Q2":q2,"minus_tprime":mtp,"Mx2":mx2,"phi":phi,"W":w,"y":y}
            for k,v in vals.items(): acc[k].append(v[keep])
            print(f"[load] {period}: chunk {ichunk}; seen={seen:,}; selected={selected:,}",flush=True)
    result={k:np.concatenate(v) if v else np.empty(0) for k,v in acc.items()}
    result["period"]=np.full(result["xB"].size,period)
    print(f"[load] {period}: DONE; selected={result['xB'].size:,}",flush=True)
    return result

def finite(a): return a[np.isfinite(a)]

def hist2d(ax,x,y,bins,ranges,xlabel,ylabel):
    good=np.isfinite(x)&np.isfinite(y); h,xe,ye=np.histogram2d(x[good],y[good],bins=bins,range=ranges)
    pos=h[h>0]; norm=LogNorm(vmin=max(1.,float(pos.min())) if pos.size else 1.,vmax=max(2.,float(h.max())))
    mesh=ax.pcolormesh(xe,ye,h.T,shading="auto",norm=norm,cmap=CMAP)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); return mesh

def save(fig,path):
    fig.savefig(path,dpi=220,bbox_inches="tight"); plt.close(fig); print(f"[plot] wrote {path}",flush=True)

def analysis_bin_title(ix,it,n):
    xl,xh=XB_BINS[ix]; tl,th=TP_BINS[it]
    return rf"${xl:.2f}\leq x_B<{xh:.2f}$, ${tl:.2f}\leq -t'<{th:.2f}$"+f"\nN={n:,}"

def plot_variable_by_analysis_bin(data,key,edges,xlabel,path):
    """4x6 grid matching the actual xB (rows) x -t' (columns) analysis binning."""
    fig,axes=plt.subplots(4,6,figsize=(17.5,10.5),sharex=True,sharey=False)
    for ix in range(4):
        for it in range(6):
            ax=axes[ix,it]; xl,xh=XB_BINS[ix]; tl,th=TP_BINS[it]
            m=(data["xB"]>=xl)&(data["xB"]<xh)&(data["minus_tprime"]>=tl)&(data["minus_tprime"]<th)&np.isfinite(data[key])
            vals=data[key][m]
            ax.hist(vals,bins=edges,histtype="step",linewidth=1.25)
            ax.set_title(analysis_bin_title(ix,it,vals.size),fontsize=8)
            ax.tick_params(direction="in",top=True,right=True,labelsize=8)
            if ix==3: ax.set_xlabel(xlabel,fontsize=9)
            if it==0: ax.set_ylabel("Events",fontsize=9)
    fig.suptitle(rf"{xlabel} distributions in the 24 analysis bins",fontsize=15)
    fig.tight_layout(rect=(0,0,1,0.965)); save(fig,path)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--cuts",type=Path,default=DEFAULT_CUTS); ap.add_argument("--output-dir",type=Path,default=DEFAULT_OUT)
    for p in PERIODS: ap.add_argument(f"--{p}",type=Path,default=DEFAULT_INPUTS[p])
    args=ap.parse_args(); args.output_dir.mkdir(parents=True,exist_ok=True)
    windows=load_nominal_windows(args.cuts); datasets=[collect(p,getattr(args,p),windows[p]) for p in PERIODS]
    data={k:np.concatenate([d[k] for d in datasets]) for k in datasets[0]}; n=data["xB"].size
    print(f"[summary] combined selected sample: {n:,} events",flush=True)
    q2max=max(8.0,np.nanpercentile(data["Q2"],99.7)); wmax=max(4.5,np.nanpercentile(data["W"],99.7)) if np.isfinite(data["W"]).any() else 4.5

    one_d=[("xB",np.linspace(.10,.60,51),r"$x_B$"),("Q2",np.linspace(1,q2max,55),r"$Q^2$ (GeV$^2$)"),("minus_tprime",np.linspace(.05,1.25,49),r"$-t'$ (GeV$^2$)"),("phi",np.linspace(0,360,49),r"$\phi$ (deg)")]
    if np.isfinite(data["W"]).any(): one_d.append(("W",np.linspace(2,wmax,50),r"$W$ (GeV)"))
    if np.isfinite(data["y"]).any(): one_d.append(("y",np.linspace(0,.8,49),r"$y$"))
    fig,axes=plt.subplots(2,3,figsize=(11.5,7.0))
    for ax,(key,bins,label) in zip(axes.flat,one_d): ax.hist(finite(data[key]),bins=bins,histtype="step",linewidth=1.5); ax.set_xlabel(label); ax.set_ylabel("Events"); ax.tick_params(direction="in",top=True,right=True)
    for ax in axes.flat[len(one_d):]: ax.axis("off")
    fig.suptitle(r"RGC $e n \pi^+$ kinematics after nominal exclusivity selection"); fig.tight_layout(rect=(0,0,1,.96)); save(fig,args.output_dir/"01_kinematic_1d.png")

    fig,axes=plt.subplots(1,3,figsize=(14,4.1)); specs=[
      (data["xB"],data["Q2"],(50,55),((.10,.60),(1,q2max)),r"$x_B$",r"$Q^2$ (GeV$^2$)"),
      (data["xB"],data["minus_tprime"],(50,48),((.10,.60),(.05,1.25)),r"$x_B$",r"$-t'$ (GeV$^2$)"),
      (data["Q2"],data["minus_tprime"],(55,48),((1,q2max),(.05,1.25)),r"$Q^2$ (GeV$^2$)",r"$-t'$ (GeV$^2$)")]
    for ax,s in zip(axes,specs): mesh=hist2d(ax,*s); fig.colorbar(mesh,ax=ax,label="Events")
    fig.suptitle("Non-redundant kinematic correlations"); fig.tight_layout(rect=(0,0,1,.94)); save(fig,args.output_dir/"02_core_correlations.png")

    fig,axes=plt.subplots(3,1,figsize=(9,10),sharex=True); phi_edges=np.linspace(0,360,49)
    ys=[(data["Q2"],np.linspace(1,q2max,55),r"$Q^2$ (GeV$^2$)"),(data["xB"],np.linspace(.10,.60,51),r"$x_B$"),(data["minus_tprime"],np.linspace(.05,1.25,49),r"$-t'$ (GeV$^2$)")]
    for ax,(yy,yedges,ylabel) in zip(axes,ys):
        good=np.isfinite(data["phi"])&np.isfinite(yy); h,xe,ye=np.histogram2d(data["phi"][good],yy[good],bins=(phi_edges,yedges)); pos=h[h>0]
        mesh=ax.pcolormesh(xe,ye,h.T,shading="auto",norm=LogNorm(vmin=max(1.,float(pos.min())) if pos.size else 1.,vmax=max(2.,float(h.max()))),cmap=CMAP)
        ax.set_ylabel(ylabel); fig.colorbar(mesh,ax=ax,label="Events")
    axes[-1].set_xlabel(r"$\phi$ (deg)"); axes[-1].set_xlim(0,360); fig.suptitle(r"Common-event-sample correlations with $\phi$"); fig.tight_layout(rect=(0,0,1,.96)); save(fig,args.output_dir/"03_phi_correlations.png")

    fig,axes=plt.subplots(2,2,figsize=(10.5,7.5),sharex=True,sharey=True)
    for ax,(xl,xh) in zip(axes.flat,XB_BINS):
        m=(data["xB"]>=xl)&(data["xB"]<xh); counts,edges=np.histogram(data["phi"][m],bins=phi_edges); centers=.5*(edges[:-1]+edges[1:]); density=counts/counts.sum() if counts.sum() else counts.astype(float)
        ax.step(centers,density,where="mid",linewidth=1.5); ax.set_title(fr"${xl:.2f} \leq x_B < {xh:.2f}$  ($N={counts.sum():,}$)"); ax.set_xlim(0,360); ax.tick_params(direction="in",top=True,right=True)
    for ax in axes[-1,:]: ax.set_xlabel(r"$\phi$ (deg)")
    for ax in axes[:,0]: ax.set_ylabel("Fraction of events / bin")
    fig.suptitle(r"$\phi$ acceptance shape in the four analysis $x_B$ intervals"); fig.tight_layout(rect=(0,0,1,.95)); save(fig,args.output_dir/"04_phi_projections.png")

    # New: directly show the Q2 and W populations entering every extracted analysis bin.
    plot_variable_by_analysis_bin(data,"Q2",np.linspace(1,q2max,45),r"$Q^2$ (GeV$^2$)",args.output_dir/"05_Q2_by_analysis_bin.png")
    if np.isfinite(data["W"]).any():
        plot_variable_by_analysis_bin(data,"W",np.linspace(2,wmax,45),r"$W$ (GeV)",args.output_dir/"06_W_by_analysis_bin.png")

    # New: Q2-W correlation in each xB row. This compactly exposes how DIS phase space evolves across xB.
    if np.isfinite(data["W"]).any():
        fig,axes=plt.subplots(2,2,figsize=(10.5,8.0),sharex=True,sharey=True)
        for ax,(xl,xh) in zip(axes.flat,XB_BINS):
            m=(data["xB"]>=xl)&(data["xB"]<xh); mesh=hist2d(ax,data["Q2"][m],data["W"][m],(50,45),((1,q2max),(2,wmax)),r"$Q^2$ (GeV$^2$)",r"$W$ (GeV)"); ax.set_title(fr"${xl:.2f}\leq x_B<{xh:.2f}$  ($N={np.count_nonzero(m):,}$)"); fig.colorbar(mesh,ax=ax,label="Events")
        fig.suptitle(r"$Q^2$--$W$ phase space in the four analysis $x_B$ intervals"); fig.tight_layout(rect=(0,0,1,.95)); save(fig,args.output_dir/"07_Q2_W_by_xB.png")

    summary={"cuts_json":str(args.cuts.resolve()),"nominal_sigma_multiple":2.0,"total_selected_events":int(n),"selected_by_period":{p:int(np.count_nonzero(data["period"]==p)) for p in PERIODS},"inputs":{p:str(getattr(args,p).resolve()) for p in PERIODS},"explicit_dis_cuts":{"Q2_min_gev2":1.0,"W_min_gev_if_branch_present":2.0,"y_max_if_branch_present":.80},"analysis_binning":{"xB":XB_BINS,"minus_tprime_gev2":TP_BINS},"plot_colormap_2d":CMAP}
    (args.output_dir/"kinematic_summary.json").write_text(json.dumps(summary,indent=2)); print(f"[summary] wrote {args.output_dir/'kinematic_summary.json'}",flush=True)

if __name__ == "__main__": main()
