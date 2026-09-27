#!/usr/bin/env python3
"""
Chunked analysis of pid_exclusivity_skim.groovy worker text files.

Designed for very large and still-growing files.  Nothing proportional to the
full input size is retained in memory: cut-flow counters and histograms are
accumulated chunk by chunk.

Usage:
  python analyze_pid_exclusivity.py /scratch/thayward/test_su22.txt
  python analyze_pid_exclusivity.py worker1.txt worker2.txt --output-dir out
"""
import argparse, time
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

COLS = "run event period assigned_pid charge status fd beta chi2pid p px py pz theta phi vz e_p e_px e_py e_pz e_theta e_phi e_vz Ebeam Q2 W xB y t tmin minus_tprime Mx2 analysis_bin mx2_mu mx2_sigma mx2_nsigma pass_qa pass_fd pass_p pass_chi2pid pass_Q2 pass_W pass_y pass_phase_space pass_mx2 pass_final".split()

# Only columns needed by this diagnostic are parsed.
USE = ["assigned_pid","beta","chi2pid","p","analysis_bin",
       "pass_fd","pass_p","pass_chi2pid","pass_Q2","pass_W","pass_y",
       "pass_phase_space","pass_mx2"]

PIDS = [211,321,2212]
PID_NAMES = ["pi","K","p"]
PID_LABELS = {211:r"REC $\pi^+$",321:r"REC $K^+$",2212:r"REC $p$"}
MASSES = {211:0.13957039,321:0.493677,2212:0.9382720813}

STAGES = [
    "Positive FD",
    "+ 0.5 < p < 5 GeV",
    "+ DIS (Q2,W,y)",
    "+ analysis (xB,-t') bin",
    "+ nominal Mx2 exclusivity",
    "+ |chi2PID| < 3.5",
]

P_EDGES = np.linspace(0.5,5.0,181)
BETA_EDGES = np.linspace(0.2,1.05,151)
CHI_EDGES = np.linspace(-12,12,193)

def stage_masks(d):
    finite=np.isfinite(d["p"]) & np.isfinite(d["beta"])
    m=[]
    m.append(finite & (d.pass_fd==1))
    m.append(m[-1] & (d.pass_p==1))
    m.append(m[-1] & (d.pass_Q2==1) & (d.pass_W==1) & (d.pass_y==1))
    m.append(m[-1] & (d.pass_phase_space==1))
    m.append(m[-1] & (d.pass_mx2==1))
    m.append(m[-1] & (d.pass_chi2pid==1))
    return m

def beta_expected(p,m):
    return p/np.sqrt(p*p+m*m)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("inputs",nargs="+")
    ap.add_argument("--output-dir",default="pid_exclusivity_diagnostics")
    ap.add_argument("--chunksize",type=int,default=250000)
    args=ap.parse_args()
    out=Path(args.output_dir); out.mkdir(parents=True,exist_ok=True)

    nstage=len(STAGES)
    counts=np.zeros(nstage,dtype=np.int64)
    composition=np.zeros((nstage,4),dtype=np.int64) # pi,K,p,other
    hbeta=np.zeros((nstage,len(P_EDGES)-1,len(BETA_EDGES)-1),dtype=np.int64)
    # pre/post exclusivity split by assigned species
    hspecies=np.zeros((2,3,len(P_EDGES)-1,len(BETA_EDGES)-1),dtype=np.int64)
    hchi=np.zeros((3,len(CHI_EDGES)-1),dtype=np.int64)
    bybin=np.zeros((2,24,4),dtype=np.int64)

    total=0; nchunk=0
    t0=time.time(); last=t0

    print(f"[start] chunk size = {args.chunksize:,} rows",flush=True)

    for fn in args.inputs:
        print(f"[open] {fn}",flush=True)
        # comment='#' ignores header; names supplies schema; usecols avoids parsing/storing unused columns.
        reader=pd.read_csv(fn,sep=r"\s+",comment="#",names=COLS,usecols=USE,
                           chunksize=args.chunksize,engine="c",
                           on_bad_lines="skip")
        for d in reader:
            nchunk+=1
            n=len(d); total+=n
            masks=stage_masks(d)

            for isg,sel in enumerate(masks):
                q=d.loc[sel]
                counts[isg]+=len(q)
                for j,pid in enumerate(PIDS):
                    composition[isg,j]+=int((q.assigned_pid==pid).sum())
                composition[isg,3]+=int((~q.assigned_pid.isin(PIDS)).sum())
                if len(q):
                    h,_,_=np.histogram2d(q.p,q.beta,bins=(P_EDGES,BETA_EDGES))
                    hbeta[isg]+=h.astype(np.int64)

            pre=masks[4]; post=masks[5]
            for ip,sel in enumerate((pre,post)):
                for j,pid in enumerate(PIDS):
                    q=d.loc[sel & (d.assigned_pid==pid)]
                    if len(q):
                        h,_,_=np.histogram2d(q.p,q.beta,bins=(P_EDGES,BETA_EDGES))
                        hspecies[ip,j]+=h.astype(np.int64)
                for ib in range(1,25):
                    qb=d.loc[sel & (d.analysis_bin==ib)]
                    for j,pid in enumerate(PIDS):
                        bybin[ip,ib-1,j]+=int((qb.assigned_pid==pid).sum())
                    bybin[ip,ib-1,3]+=int((~qb.assigned_pid.isin(PIDS)).sum())

            for j,pid in enumerate(PIDS):
                z=d.loc[pre & (d.assigned_pid==pid),"chi2pid"].to_numpy()
                z=z[np.isfinite(z)]
                hchi[j]+=np.histogram(z,bins=CHI_EDGES)[0]

            now=time.time()
            dt=now-last; elapsed=now-t0
            rate=total/elapsed if elapsed else 0
            print(f"[read] {total:,} candidates | chunk {nchunk:,} | "
                  f"{elapsed:.1f} s elapsed | {rate:,.0f} rows/s | "
                  f"final-so-far {counts[-1]:,}",flush=True)
            last=now

    # Tables
    rows=[]
    prev=total
    for i,name in enumerate(STAGES):
        rows.append(dict(stage=name,candidates=int(counts[i]),
                         fraction_previous=counts[i]/prev if prev else np.nan,
                         assigned_pi=int(composition[i,0]),
                         assigned_K=int(composition[i,1]),
                         assigned_p=int(composition[i,2]),
                         assigned_other=int(composition[i,3])))
        prev=counts[i]
    cf=pd.DataFrame(rows)
    cf.to_csv(out/"cutflow.csv",index=False)
    print("\n"+cf.to_string(index=False),flush=True)

    surv=[]
    for j,pid in enumerate(PIDS):
        npre=composition[4,j]; npost=composition[5,j]
        surv.append([PID_NAMES[j],npre,npost,npost/npre if npre else np.nan])
    pd.DataFrame(surv,columns=["assigned_species","before_chi2pid","after_chi2pid",
                               "survival_fraction"]).to_csv(
        out/"pid_survival_after_exclusivity.csv",index=False)

    # Plot accumulated beta(p) maps.
    pp=np.linspace(.5,5,500)
    extent=[P_EDGES[0],P_EDGES[-1],BETA_EDGES[0],BETA_EDGES[-1]]
    for i,name in enumerate(STAGES):
        H=hbeta[i].T
        fig,ax=plt.subplots(figsize=(8,6))
        positive=H[H>0]
        if positive.size:
            im=ax.imshow(np.ma.masked_where(H<=0,H),origin="lower",aspect="auto",
                         extent=extent,norm=LogNorm(vmin=1,vmax=max(2,H.max())))
            fig.colorbar(im,ax=ax,label="Candidates")
        for pid in PIDS:
            ax.plot(pp,beta_expected(pp,MASSES[pid]),lw=1.4,label=PID_LABELS[pid]+" mass")
        ax.set(xlabel="p (GeV)",ylabel=r"$\beta$",xlim=(.5,5),ylim=(.2,1.05),
               title=f"{name}\nN={counts[i]:,}")
        ax.legend(loc="lower right",fontsize=9); fig.tight_layout()
        fig.savefig(out/f"beta_vs_p_stage_{i+1:02d}.png",dpi=180); plt.close(fig)

    # Species-separated pre/post PID maps.
    for ip,tag in enumerate(("pre_pid","post_pid")):
        fig,axs=plt.subplots(1,3,figsize=(15,4.8),sharex=True,sharey=True)
        for j,(ax,pid) in enumerate(zip(axs,PIDS)):
            H=hspecies[ip,j].T
            if np.any(H):
                ax.imshow(np.ma.masked_where(H<=0,H),origin="lower",aspect="auto",
                          extent=extent,norm=LogNorm(vmin=1,vmax=max(2,H.max())))
            for hyp in PIDS: ax.plot(pp,beta_expected(pp,MASSES[hyp]),lw=1.0)
            ax.set_title(f"{PID_LABELS[pid]}\nN={H.sum():,}"); ax.set_xlabel("p (GeV)")
        axs[0].set_ylabel(r"$\beta$"); axs[0].set_xlim(.5,5); axs[0].set_ylim(.2,1.05)
        fig.suptitle("After nominal exclusivity, "+("before PID cut" if ip==0 else "after PID cut"))
        fig.tight_layout(); fig.savefig(out/f"assigned_species_beta_vs_p_{tag}.png",dpi=180); plt.close(fig)

    # chi2PID
    centers=.5*(CHI_EDGES[:-1]+CHI_EDGES[1:])
    fig,ax=plt.subplots(figsize=(8,5.6))
    for j,pid in enumerate(PIDS):
        ax.step(centers,hchi[j],where="mid",lw=1.6,label=f"{PID_LABELS[pid]} (N={hchi[j].sum():,})")
    ax.axvline(-3.5,ls="--"); ax.axvline(3.5,ls="--"); ax.set_yscale("log")
    ax.set(xlabel=r"$\chi^2_{\rm PID}$",ylabel="Candidates",
           title="After all non-PID channel/exclusivity cuts")
    ax.legend(); fig.tight_layout(); fig.savefig(out/"chi2pid_after_exclusivity.png",dpi=180); plt.close(fig)

    # 24-bin tables and plots.
    for ip,tag in enumerate(("pre_pid","post_pid")):
        t=pd.DataFrame({"bin":np.arange(1,25),"pi":bybin[ip,:,0],"K":bybin[ip,:,1],
                        "p":bybin[ip,:,2],"other":bybin[ip,:,3]})
        t.to_csv(out/f"assigned_pid_by_analysis_bin_{tag}.csv",index=False)
        fig,ax=plt.subplots(figsize=(9,5.3)); bottom=np.zeros(24)
        for c in ("pi","K","p","other"):
            ax.bar(t["bin"],t[c],bottom=bottom,label=c); bottom+=t[c].to_numpy()
        ax.set_yscale("log"); ax.set(xlabel="Analysis bin",ylabel="Candidates",
                                    title=f"REC-assigned PID composition: {tag}")
        ax.legend(); fig.tight_layout(); fig.savefig(out/f"assigned_pid_by_analysis_bin_{tag}.png",dpi=180); plt.close(fig)

    print(f"\n[done] {total:,} candidates in {time.time()-t0:.1f} s; outputs: {out}",flush=True)

if __name__=="__main__":
    main()
