#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

COLS="run event period assigned_pid charge status fd beta chi2pid p px py pz theta phi vz e_p e_px e_py e_pz e_theta e_phi e_vz Ebeam Q2 W xB y t tmin minus_tprime Mx2 analysis_bin mx2_mu mx2_sigma mx2_nsigma pass_qa pass_fd pass_p pass_chi2pid pass_Q2 pass_W pass_y pass_phase_space pass_mx2 pass_final".split()
MASSES={211:0.13957039,321:0.493677,2212:0.9382720813}
LABELS={211:"REC pi+",321:"REC K+",2212:"REC p"}

def read_file(path):
    rows=[]; bad=0
    with open(path,errors="replace") as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            z=line.split()
            if len(z)!=len(COLS): bad+=1; continue
            rows.append(z)
    d=pd.DataFrame(rows,columns=COLS)
    if d.empty: raise RuntimeError(f"No complete rows in {path}")
    for c in COLS:
        if c!="period": d[c]=pd.to_numeric(d[c],errors="coerce")
    print(f"[read] {path}: {len(d):,} rows; {bad:,} incomplete/malformed skipped")
    return d

def masks(d):
    finite=np.isfinite(d.p)&np.isfinite(d.beta)&np.isfinite(d.Q2)&np.isfinite(d.W)&np.isfinite(d.y)&np.isfinite(d.Mx2)
    m={}
    m["Positive FD"]=finite&(d.pass_fd==1)
    m["+ 0.5 < p < 5 GeV"]=m["Positive FD"]&(d.pass_p==1)
    m["+ DIS (Q2,W,y)"]=m["+ 0.5 < p < 5 GeV"]&(d.pass_Q2==1)&(d.pass_W==1)&(d.pass_y==1)
    m["+ analysis (xB,-t') bin"]=m["+ DIS (Q2,W,y)"]&(d.pass_phase_space==1)
    m["+ nominal Mx2 exclusivity"]=m["+ analysis (xB,-t') bin"]&(d.pass_mx2==1)
    m["+ |chi2PID| < 3.5"]=m["+ nominal Mx2 exclusivity"]&(d.pass_chi2pid==1)
    return m

def beta(p,m): return p/np.sqrt(p*p+m*m)

def main():
    a=argparse.ArgumentParser(description="Analyze RGC PID/exclusivity skim, including a still-growing worker file.")
    a.add_argument("inputs",nargs="+")
    a.add_argument("--output-dir",default="pid_exclusivity_diagnostics")
    x=a.parse_args(); out=Path(x.output_dir); out.mkdir(parents=True,exist_ok=True)
    d=pd.concat([read_file(p) for p in x.inputs],ignore_index=True)
    m=masks(d)
    print(f"[total] {len(d):,} candidate rows")

    # Sequential cut flow and REC-assigned composition.
    rows=[]; prev=len(d)
    for name,sel in m.items():
        q=d.loc[sel]; n=len(q)
        rows.append(dict(stage=name,candidates=n,fraction_previous=n/prev if prev else np.nan,
                         assigned_pi=int((q.assigned_pid==211).sum()),
                         assigned_K=int((q.assigned_pid==321).sum()),
                         assigned_p=int((q.assigned_pid==2212).sum()),
                         assigned_other=int((~q.assigned_pid.isin([211,321,2212])).sum())))
        prev=n
    cf=pd.DataFrame(rows); cf.to_csv(out/"cutflow.csv",index=False)
    print("\n"+cf.to_string(index=False))

    # beta(p) at every cut stage.
    pp=np.linspace(.5,5,500)
    for i,(name,sel) in enumerate(m.items(),1):
        q=d.loc[sel,["p","beta"]]
        if q.empty: continue
        fig,ax=plt.subplots(figsize=(8,6))
        h=ax.hist2d(q.p,q.beta,bins=(180,150),range=((.5,5),(.2,1.05)),norm=LogNorm())
        fig.colorbar(h[3],ax=ax,label="Candidates")
        for pid in (211,321,2212): ax.plot(pp,beta(pp,MASSES[pid]),lw=1.4,label=LABELS[pid]+" mass")
        ax.set(xlabel="p (GeV)",ylabel=r"$\beta$",xlim=(.5,5),ylim=(.2,1.05),title=f"{name}\nN={len(q):,}")
        ax.legend(loc="lower right",fontsize=9); fig.tight_layout()
        fig.savefig(out/f"beta_vs_p_stage_{i:02d}.png",dpi=180); plt.close(fig)

    # Critical pre-PID population split by REC assignment.
    pre=m["+ nominal Mx2 exclusivity"]; post=m["+ |chi2PID| < 3.5"]
    for tag,sel in [("pre_pid",pre),("post_pid",post)]:
        fig,axs=plt.subplots(1,3,figsize=(15,4.8),sharex=True,sharey=True)
        for ax,pid in zip(axs,(211,321,2212)):
            q=d.loc[sel&(d.assigned_pid==pid),["p","beta"]]
            if len(q): ax.hist2d(q.p,q.beta,bins=(120,110),range=((.5,5),(.2,1.05)),norm=LogNorm())
            for hyp in (211,321,2212): ax.plot(pp,beta(pp,MASSES[hyp]),lw=1.0)
            ax.set_title(f"{LABELS[pid]}\nN={len(q):,}"); ax.set_xlabel("p (GeV)")
        axs[0].set_ylabel(r"$\beta$"); axs[0].set_xlim(.5,5); axs[0].set_ylim(.2,1.05)
        fig.suptitle("After exclusivity, "+("before PID cut" if tag=="pre_pid" else "after PID cut"))
        fig.tight_layout(); fig.savefig(out/f"assigned_species_beta_vs_p_{tag}.png",dpi=180); plt.close(fig)

    # chi2PID after all non-PID cuts.
    fig,ax=plt.subplots(figsize=(8,5.6)); bins=np.linspace(-12,12,193)
    surv=[]
    for pid in (211,321,2212):
        z=d.loc[pre&(d.assigned_pid==pid),"chi2pid"].dropna()
        zz=d.loc[post&(d.assigned_pid==pid),"chi2pid"].dropna()
        surv.append((LABELS[pid],len(z),len(zz),len(zz)/len(z) if len(z) else np.nan))
        if len(z): ax.hist(z,bins=bins,histtype="step",lw=1.6,label=f"{LABELS[pid]} (N={len(z):,})")
    ax.axvline(-3.5,ls="--"); ax.axvline(3.5,ls="--"); ax.set_yscale("log")
    ax.set(xlabel=r"$\chi^2_{\rm PID}$",ylabel="Candidates",title="After all non-PID channel/exclusivity cuts")
    ax.legend(); fig.tight_layout(); fig.savefig(out/"chi2pid_after_exclusivity.png",dpi=180); plt.close(fig)
    pd.DataFrame(surv,columns=["assigned_species","before_chi2pid","after_chi2pid","survival_fraction"]).to_csv(out/"pid_survival_after_exclusivity.csv",index=False)

    # 24-bin composition, before/after PID.
    for tag,sel in [("pre_pid",pre),("post_pid",post)]:
        rr=[]
        for ib in range(1,25):
            q=d.loc[sel&(d.analysis_bin==ib)]
            rr.append([ib,*[int((q.assigned_pid==p).sum()) for p in (211,321,2212)],int((~q.assigned_pid.isin([211,321,2212])).sum())])
        t=pd.DataFrame(rr,columns=["bin","pi","K","p","other"]); t.to_csv(out/f"assigned_pid_by_analysis_bin_{tag}.csv",index=False)
        fig,ax=plt.subplots(figsize=(9,5.3)); bottom=np.zeros(24)
        for c in ["pi","K","p","other"]:
            ax.bar(t["bin"],t[c],bottom=bottom,label=c); bottom+=t[c].to_numpy()
        ax.set_yscale("log"); ax.set(xlabel="Analysis bin",ylabel="Candidates",title=f"REC-assigned PID composition: {tag}")
        ax.legend(); fig.tight_layout(); fig.savefig(out/f"assigned_pid_by_analysis_bin_{tag}.png",dpi=180); plt.close(fig)

    print(f"\nDone: {out}")

if __name__=="__main__": main()
