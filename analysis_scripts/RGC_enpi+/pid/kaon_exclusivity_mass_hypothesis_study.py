#!/usr/bin/env python3
"""Kaon-mass-hypothesis exclusivity diagnostic for RGC e pi+ analysis.

Starts from the nominal NH3 paper-version ROOT trees.  It reproduces the final
24-bin + nominal Mx2 selection using the pion hypothesis, then keeps the measured
hadron three-momentum fixed, replaces m_pi -> m_K, and recomputes Mx2.

This is a diagnostic of the *kinematic discriminating power* of the neutron
exclusivity requirement.  The kaon-hypothesis survival fraction is NOT itself a
kaon contamination estimate because the selected sample is overwhelmingly real
pions.
"""

from pathlib import Path
import json
import math
import numpy as np
import pandas as pd
import uproot
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
TREE_NAME = "PhysicsEvents"
CHUNK = "250 MB"
OUTDIR = Path("output/kaon_mass_hypothesis")

PAPER_DIR = Path("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/paper_versions")
INPUTS = {
    "su22": PAPER_DIR / "rgc_su22_inb_NH3_epi+_mom_corrections.root",
    "fa22": PAPER_DIR / "rgc_fa22_inb_NH3_epi+_mom_corrections.root",
    "sp23": PAPER_DIR / "rgc_sp23_inb_NH3_epi+_mom_corrections.root",
}
BEAM_E = {"su22": 10.5473, "fa22": 10.5563, "sp23": 10.5593}

CUT_JSON = Path("../channel_selection/output/channel_selection_mx2_fit_stability/final_carbon_assisted_cuts/tables/final_carbon_assisted_mx2_cuts.json")

MP = 0.9382720813
ME = 0.00051099895
MPI = 0.13957039
MK = 0.493677

XB_BINS = ((0.10,0.25),(0.25,0.35),(0.35,0.45),(0.45,0.60))
MTP_BINS = ((0.05,0.25),(0.25,0.45),(0.45,0.65),(0.65,0.85),(0.85,1.05),(1.05,1.25))

# Broad where useful; merge the sparse highest-momentum region.
P_RANGES = ((0.5,5.0),(5.0,5.5),(5.5,6.0),(6.0,6.5),(6.5,7.25),(7.25,8.0),(8.0,10.5))

ALIASES = {
    "xB": ("xB","x","xb","x_b"),
    "Q2": ("Q2","q2"),
    "tprime": ("tprime","t_prime","tp","tPrime"),
    "Mx2": ("Mx2","mx2","Mx2_epi","Mx2_epip","missing_mass_squared","missing_mass2"),
    "W": ("W","w"),
    "y": ("y","inelasticity"),
    "e_p": ("e_p","p_e","electron_p"),
    "e_theta": ("e_theta","theta_e","electron_theta"),
    "e_phi": ("e_phi","phi_e","electron_phi"),
    "p_p": ("p_p","p_pi","p_pip","pion_p"),
    "p_theta": ("p_theta","theta_p","theta_pi","pion_theta"),
    "p_phi": ("p_phi","phi_p","phi_pi","pion_phi"),
}


def resolve(tree, logical, required=True):
    names = set(tree.keys())
    for c in ALIASES[logical]:
        if c in names:
            return c
    if required:
        raise KeyError(f"Could not resolve {logical}; tried {ALIASES[logical]}")
    return None


def angle_to_rad(a):
    a = np.asarray(a, dtype=float)
    finite = a[np.isfinite(a)]
    if len(finite) and np.nanpercentile(np.abs(finite), 99) > 2.0*math.pi + 0.2:
        return np.deg2rad(a)
    return a


def components(p, th, ph):
    st = np.sin(th)
    return p*st*np.cos(ph), p*st*np.sin(ph), p*np.cos(th)


def bin_number(xb, minus_tp):
    out = np.zeros(len(xb), dtype=np.int16)
    for ix,(xl,xh) in enumerate(XB_BINS):
        xm = (xb >= xl) & (xb < xh if ix < len(XB_BINS)-1 else xb <= xh)
        for it,(tl,th) in enumerate(MTP_BINS):
            tm = (minus_tp >= tl) & (minus_tp < th if it < len(MTP_BINS)-1 else minus_tp <= th)
            out[xm & tm] = ix*len(MTP_BINS) + it + 1
    return out


def load_windows(path):
    payload = json.loads(path.read_text())
    pp = payload["periods"]
    out = {}
    for period in INPUTS:
        out[period] = {}
        for row in pp[period]:
            b = int(row.get("bin_number", int(row["x_index"])*6 + int(row["t_index"]) + 1))
            lo,hi = row["nominal"]
            out[period][b] = (float(lo),float(hi))
    return out


def alternative_mx2(beam_e, ep, eth, eph, hp, hth, hph, mass):
    eE = np.sqrt(ep*ep + ME*ME)
    epx,epy,epz = components(ep,eth,eph)
    hpx,hpy,hpz = components(hp,hth,hph)
    beam_p = math.sqrt(max(beam_e*beam_e - ME*ME, 0.0))
    qE = beam_e - eE
    qx,qy,qz = -epx,-epy,beam_p-epz
    hE = np.sqrt(hp*hp + mass*mass)
    Em = MP + qE - hE
    pxm,pym,pzm = qx-hpx,qy-hpy,qz-hpz
    return Em*Em - pxm*pxm - pym*pym - pzm*pzm


def load_period(period, windows):
    path = INPUTS[period]
    print(f"[load] {period}: {path}", flush=True)
    with uproot.open(path) as f:
        tree = f[TREE_NAME]
        br = {k:resolve(tree,k,required=(k not in {"W","y"})) for k in ALIASES}
        expr = [x for x in br.values() if x]
        acc = {k:[] for k in ("bin","p","mx_pi","mx_k","k_survives")}
        nseen=nsel=0
        for arrays in tree.iterate(expressions=expr, step_size=CHUNK, library="np"):
            xb=np.asarray(arrays[br["xB"]],float); q2=np.asarray(arrays[br["Q2"]],float)
            mtp=-np.asarray(arrays[br["tprime"]],float); mx=np.asarray(arrays[br["Mx2"]],float)
            hp=np.asarray(arrays[br["p_p"]],float)
            ep=np.asarray(arrays[br["e_p"]],float)
            eth=angle_to_rad(arrays[br["e_theta"]]); eph=angle_to_rad(arrays[br["e_phi"]])
            hth=angle_to_rad(arrays[br["p_theta"]]); hph=angle_to_rad(arrays[br["p_phi"]])
            b=bin_number(xb,mtp)
            keep=np.isfinite(xb)&np.isfinite(q2)&np.isfinite(mtp)&np.isfinite(mx)&np.isfinite(hp)&np.isfinite(ep)&np.isfinite(eth)&np.isfinite(eph)&np.isfinite(hth)&np.isfinite(hph)&(b>0)&(q2>1.0)
            if br["W"]: keep &= np.asarray(arrays[br["W"]],float)>2.0
            if br["y"]: keep &= np.asarray(arrays[br["y"]],float)<0.8
            for ib in range(1,25):
                m=b==ib
                if np.any(m):
                    lo,hi=windows[ib]
                    keep[m] &= (mx[m]>=lo)&(mx[m]<=hi)
            idx=np.flatnonzero(keep)
            nseen += len(xb); nsel += len(idx)
            if not len(idx): continue
            mxk=alternative_mx2(BEAM_E[period],ep[idx],eth[idx],eph[idx],hp[idx],hth[idx],hph[idx],MK)
            surv=np.zeros(len(idx),dtype=bool)
            for ib in range(1,25):
                m=b[idx]==ib
                if np.any(m):
                    lo,hi=windows[ib]
                    surv[m]=(mxk[m]>=lo)&(mxk[m]<=hi)
            acc["bin"].append(b[idx]); acc["p"].append(hp[idx]); acc["mx_pi"].append(mx[idx]); acc["mx_k"].append(mxk); acc["k_survives"].append(surv)
        print(f"[load] {period}: selected {nsel:,}/{nseen:,}", flush=True)
    return {k:np.concatenate(v) if v else np.array([]) for k,v in acc.items()}


def make_summary(data):
    rows=[]
    for period,d in data.items():
        for ib in range(1,25):
            m=d["bin"]==ib
            if not np.any(m): continue
            rows.append({"period":period,"analysis_bin":ib,"N":int(m.sum()),"kaon_hypothesis_survival_fraction":float(np.mean(d["k_survives"][m])),"mean_delta_Mx2_GeV2":float(np.mean(d["mx_k"][m]-d["mx_pi"][m]))})
        for lo,hi in P_RANGES:
            m=(d["p"]>=lo)&(d["p"]<hi)
            if not np.any(m): continue
            rows.append({"period":period,"analysis_bin":0,"p_min_GeV":lo,"p_max_GeV":hi,"N":int(m.sum()),"kaon_hypothesis_survival_fraction":float(np.mean(d["k_survives"][m])),"mean_delta_Mx2_GeV2":float(np.mean(d["mx_k"][m]-d["mx_pi"][m]))})
    return pd.DataFrame(rows)


def plots(data):
    # Combined arrays
    p=np.concatenate([d["p"] for d in data.values()]); mpi=np.concatenate([d["mx_pi"] for d in data.values()]); mk=np.concatenate([d["mx_k"] for d in data.values()]); surv=np.concatenate([d["k_survives"] for d in data.values()]); bins=np.concatenate([d["bin"] for d in data.values()])

    # 1) Mx2 distributions in momentum slices.
    fig,axs=plt.subplots(3,3,figsize=(15,11)); axs=axs.flat
    for ia,(lo,hi) in enumerate(P_RANGES):
        ax=axs[ia]; m=(p>=lo)&(p<hi); vals=np.concatenate((mpi[m],mk[m])) if np.any(m) else np.array([])
        if len(vals):
            qlo,qhi=np.nanpercentile(vals,[0.5,99.5]); qlo=max(-0.5,qlo); qhi=min(2.0,qhi)
            edges=np.linspace(qlo,qhi,100)
            ax.hist(mpi[m],bins=edges,histtype="step",density=True,label=r"$M_{X,\pi}^2$")
            ax.hist(mk[m],bins=edges,histtype="step",density=True,label=r"$M_{X,K}^2$")
        ax.set_title(f"{lo:.2f} < p < {hi:.2f} GeV\nN = {m.sum():,}, K-hyp. survival = {100*np.mean(surv[m]) if np.any(m) else np.nan:.1f}%")
        ax.set_xlabel(r"$M_X^2$ (GeV$^2$)"); ax.set_ylabel("Normalized entries"); ax.legend(fontsize=8)
    for ax in list(axs)[len(P_RANGES):]: ax.axis("off")
    fig.suptitle(r"Final-selected events under $\pi^+$ and $K^+$ mass hypotheses")
    fig.tight_layout(rect=(0,0,1,0.97)); fig.savefig(OUTDIR/"mx2_pion_vs_kaon_hypothesis_by_momentum.png",dpi=180); plt.close(fig)

    # 2) 2D migration plot.
    fig,ax=plt.subplots(figsize=(8,7)); good=(p<8.0)&np.isfinite(mpi)&np.isfinite(mk)
    h=ax.hist2d(mpi[good],mk[good],bins=180,norm=LogNorm())
    limlo=max(-0.2,min(np.nanpercentile(mpi[good],0.2),np.nanpercentile(mk[good],0.2))); limhi=min(1.6,max(np.nanpercentile(mpi[good],99.8),np.nanpercentile(mk[good],99.8)))
    ax.plot([limlo,limhi],[limlo,limhi],"--",linewidth=1,label="No mass-hypothesis shift")
    ax.set_xlim(limlo,limhi); ax.set_ylim(limlo,limhi); ax.set_xlabel(r"$M_{X,\pi}^2$ (GeV$^2$)"); ax.set_ylabel(r"$M_{X,K}^2$ (GeV$^2$)"); ax.legend(); fig.colorbar(h[3],ax=ax,label="Events")
    fig.tight_layout(); fig.savefig(OUTDIR/"mx2_pion_vs_kaon_hypothesis_2d.png",dpi=180); plt.close(fig)

    # 3) survival versus momentum.
    edges=np.arange(0.5,8.01,0.25); centers=0.5*(edges[:-1]+edges[1:]); frac=[]; err=[]; nvals=[]
    for lo,hi in zip(edges[:-1],edges[1:]):
        m=(p>=lo)&(p<hi); n=m.sum(); f=np.mean(surv[m]) if n else np.nan
        frac.append(f); err.append(math.sqrt(f*(1-f)/n) if n and np.isfinite(f) else np.nan); nvals.append(n)
    fig,ax=plt.subplots(figsize=(10,6)); ax.errorbar(centers,100*np.asarray(frac),yerr=100*np.asarray(err),fmt="o-"); ax.set_xlabel(r"$p_{\pi^+}$ (GeV)"); ax.set_ylabel(r"Events also inside nominal $M_X^2$ window under $K^+$ hypothesis (%)"); ax.set_ylim(0,105); ax.grid(alpha=0.25); fig.tight_layout(); fig.savefig(OUTDIR/"kaon_hypothesis_survival_vs_momentum.png",dpi=180); plt.close(fig)
    pd.DataFrame({"p_min_GeV":edges[:-1],"p_max_GeV":edges[1:],"N":nvals,"survival_fraction":frac,"binomial_error":err}).to_csv(OUTDIR/"kaon_hypothesis_survival_vs_momentum.csv",index=False)

    # 4) survival in the 24 bins, combined periods.
    xs=np.arange(1,25); f=[]; e=[]; ns=[]
    for ib in xs:
        m=bins==ib; n=m.sum(); ff=np.mean(surv[m]) if n else np.nan
        f.append(ff); e.append(math.sqrt(ff*(1-ff)/n) if n and np.isfinite(ff) else np.nan); ns.append(n)
    fig,ax=plt.subplots(figsize=(12,5.5)); ax.errorbar(xs,100*np.asarray(f),yerr=100*np.asarray(e),fmt="o-"); ax.set_xticks(xs); ax.set_xlabel("Analysis bin"); ax.set_ylabel(r"$K^+$-hypothesis exclusivity survival (%)"); ax.set_ylim(0,105); ax.grid(alpha=0.25); fig.tight_layout(); fig.savefig(OUTDIR/"kaon_hypothesis_survival_24bins.png",dpi=180); plt.close(fig)
    pd.DataFrame({"analysis_bin":xs,"N":ns,"survival_fraction":f,"binomial_error":e}).to_csv(OUTDIR/"kaon_hypothesis_survival_24bins_combined.csv",index=False)


def main():
    OUTDIR.mkdir(parents=True,exist_ok=True)
    if not CUT_JSON.is_file():
        raise FileNotFoundError(f"Missing nominal channel-selection JSON: {CUT_JSON}\nRun from the pid directory, or edit CUT_JSON.")
    windows=load_windows(CUT_JSON)
    data={p:load_period(p,windows[p]) for p in INPUTS}
    make_summary(data).to_csv(OUTDIR/"kaon_mass_hypothesis_summary.csv",index=False)
    plots(data)
    print("\nKaon mass-hypothesis diagnostic complete")
    print(f"Output: {OUTDIR}")
    print("IMPORTANT: survival fraction is a kinematic diagnostic, not a kaon contamination estimate.")

if __name__ == "__main__":
    main()
