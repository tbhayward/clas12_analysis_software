#!/usr/bin/env python3
"""Truth-anchored proton/photon |t| comparison in paired DVCS MCEvents."""
import argparse
from array import array
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from compare_t_definitions import (BEAM, MC_CURRENT, binning, output_path,
                                   exclusivity_expression)

TGAMMA_CPP = r'''
#include <cmath>
double truth_study_tgamma(double ep,double eth,double eph,
                          double gp,double gth,double gph,double E) {
  constexpr double me=0.00051099895;
  double ex=ep*std::sin(eth)*std::cos(eph), ey=ep*std::sin(eth)*std::sin(eph);
  double ez=ep*std::cos(eth);
  double gx=gp*std::sin(gth)*std::cos(gph), gy=gp*std::sin(gth)*std::sin(gph);
  double gz=gp*std::cos(gth);
  double q0=E-std::sqrt(ep*ep+me*me)-gp;
  double qx=-ex-gx, qy=-ey-gy, qz=E-ez-gz;
  return q0*q0-qx*qx-qy*qy-qz*qz;
}'''

def matrix(handle,nx,ny):
    h=handle.GetValue()
    return np.array([[h.GetBinContent(i+1,j+1) for j in range(ny)]
                     for i in range(nx)],dtype=np.int64)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input-root",type=Path,
                    default=Path("/volatile/clas12/thayward/dvcs_pass2_full_rebuild"))
    ap.add_argument("--period",choices=BEAM,default="fa18_inb")
    ap.add_argument("--max-files",type=int,default=0)
    ap.add_argument("--workers",type=int,default=10)
    ap.add_argument("--cuts-json",type=Path,
                    default=Path(__file__).resolve().parent/"output/jsons/combined_cuts.json")
    ap.add_argument("--output-dir",type=Path,default=Path("output/t_definition_quick_test"))
    ap.add_argument("--merge-low-t",action="store_true",
                    help="Combine the two lowest |t| bins for this run")
    args=ap.parse_args()
    EDGES,LABELS,MIGRATION_EDGES,MIGRATION_LABELS=binning(args.merge_low_t)
    n=len(LABELS)
    nm=len(MIGRATION_LABELS)
    if not 1<=args.workers<=10:
        ap.error("--workers must be 1–10")
    import ROOT
    ROOT.EnableImplicitMT(args.workers)
    directory=args.input_root/"mc"/args.period/"dvcs_mc"/MC_CURRENT[args.period]/"parts"
    paths=sorted(p for p in directory.glob("*.root")
                 if not p.name.endswith(".rec_only.root") and ".partial" not in p.name)
    if args.max_files>0:
        paths=paths[:args.max_files]
    if not paths:
        ap.error(f"No paired MCEvents files in {directory}")
    chain=ROOT.TChain("MCEvents")
    for path in paths:
        chain.Add(str(path))
    branches={b.GetName() for b in chain.GetListOfBranches()}
    gen_required={"gen_kinematics_valid","gen_t1","gen_x","gen_Q2"}
    rec_required={"rec_kinematics_valid","rec_t1","rec_Q2","rec_x",
                  "rec_e_p","rec_e_theta","rec_e_phi","rec_p2_p",
                  "rec_p2_theta","rec_p2_phi","rec_detector1","rec_detector2"}
    missing=(gen_required|rec_required)-branches
    if missing:
        ap.error(f"Missing paired branches: {sorted(missing)}")
    logical={name[4:] for name in branches if name.startswith("rec_")}
    exclusive=exclusivity_expression(args.cuts_json,args.period,logical,"mc")
    ROOT.gInterpreter.Declare(TGAMMA_CPP)
    beam=BEAM[args.period]
    frame=ROOT.RDataFrame(chain)
    generated=frame.Filter(
        "gen_kinematics_valid>0.5 && std::isfinite(gen_t1) && gen_t1<0 && "
        "std::isfinite(gen_Q2) && gen_Q2>1 && std::isfinite(gen_x) && "
        "gen_x>0 && gen_x<1 && "
        "0.9382720813*0.9382720813+gen_Q2*(1/gen_x-1)>4")
    generated=generated.Define("truth_t","-gen_t1")
    truth_edges=array("d",EDGES)
    outcome_edges=array("d",MIGRATION_EDGES)
    ngen=generated.Histo1D(ROOT.RDF.TH1DModel("ngen","",n,truth_edges),"truth_t")
    selected=generated.Filter("rec_kinematics_valid>0.5")
    needed_aliases={"t1","Q2","x","e_p","e_theta","e_phi","p2_p",
                    "p2_theta","p2_phi","detector1","detector2",
                    "Mx2","pTmiss","theta","theta_gamma_gamma",
                    "Delta_phi","Emiss2","Mx2_2","p1_phi"}
    for name in sorted(needed_aliases & logical):
        selected=selected.Define(name,"rec_"+name)
    selected=selected.Filter(
        "std::isfinite(t1) && t1<0 && -t1<1 && "
        "std::isfinite(Q2) && Q2>1 && std::isfinite(x) && x>0 && x<1 && "
        "0.9382720813*0.9382720813+Q2*(1/x-1)>4 && "
        "std::isfinite(e_p) && e_p>2 && std::isfinite(p2_p) && p2_p>0 && "
        "(detector1==1 || detector1==2) && (detector2==0 || detector2==1)")
    selected=selected.Filter(exclusive)
    selected=selected.Define("rec_tp","-t1")
    selected=selected.Define("rec_tgamma_signed",
        f"truth_study_tgamma(e_p,e_theta,e_phi,p2_p,p2_theta,p2_phi,{beam:.17g})")
    selected=selected.Filter("std::isfinite(rec_tgamma_signed)")
    selected=selected.Define("rec_tg","-rec_tgamma_signed")
    selected=selected.Define("proton_error","rec_tp-truth_t")
    selected=selected.Define("photon_error","rec_tg-truth_t")
    selected=selected.Define("proton_abs_error","std::abs(proton_error)")
    selected=selected.Define("photon_abs_error","std::abs(photon_error)")
    proton=selected.Histo2D(ROOT.RDF.TH2DModel(
        "hproton","",n,truth_edges,nm,outcome_edges),"truth_t","rec_tp")
    photon=selected.Histo2D(ROOT.RDF.TH2DModel(
        "hphoton","",n,truth_edges,nm,outcome_edges),"truth_t","rec_tg")
    recp=selected.Histo1D(ROOT.RDF.TH1DModel("hrecp","",n,truth_edges),"rec_tp")
    recg=selected.Histo1D(ROOT.RDF.TH1DModel("hrecg","",n,truth_edges),"rec_tg")
    selected_truth=selected.Histo1D(ROOT.RDF.TH1DModel(
        "hselected_truth","",n,truth_edges),"truth_t")
    abs_p=selected.Histo1D(ROOT.RDF.TH1DModel(
        "habs_p","",n,truth_edges),"truth_t","proton_abs_error")
    abs_g=selected.Histo1D(ROOT.RDF.TH1DModel(
        "habs_g","",n,truth_edges),"truth_t","photon_abs_error")
    bias_p=selected.Histo1D(ROOT.RDF.TH1DModel(
        "hbias_p","",n,truth_edges),"truth_t","proton_error")
    bias_g=selected.Histo1D(ROOT.RDF.TH1DModel(
        "hbias_g","",n,truth_edges),"truth_t","photon_error")
    near_p=selected.Filter("proton_abs_error<0.05").Histo1D(ROOT.RDF.TH1DModel(
        "hnear_p","",n,truth_edges),"truth_t")
    near_g=selected.Filter("photon_abs_error<0.05").Histo1D(ROOT.RDF.TH1DModel(
        "hnear_g","",n,truth_edges),"truth_t")
    residual_edges=np.linspace(-2.0,5.0,351)
    residual_bins=array("d",residual_edges)
    residual_p=selected.Histo2D(ROOT.RDF.TH2DModel(
        "hresidual_p","",n,truth_edges,350,residual_bins),"truth_t","proton_error")
    residual_g=selected.Histo2D(ROOT.RDF.TH2DModel(
        "hresidual_g","",n,truth_edges,350,residual_bins),"truth_t","photon_error")
    nselected=selected.Count()
    ROOT.RDF.RunGraphs([ngen,proton,photon,recp,recg,nselected,
                        selected_truth,abs_p,abs_g,bias_p,bias_g,near_p,near_g,
                        residual_p,residual_g])
    gen=np.array([ngen.GetValue().GetBinContent(i+1) for i in range(n)],dtype=np.int64)
    pmat=matrix(proton,n,nm)
    gmat=matrix(photon,n,nm)
    rp=np.array([recp.GetValue().GetBinContent(i+1) for i in range(n)],dtype=np.int64)
    rg=np.array([recg.GetValue().GetBinContent(i+1) for i in range(n)],dtype=np.int64)
    def bins(handle):
        h=handle.GetValue()
        return np.array([h.GetBinContent(i+1) for i in range(n)],dtype=float)
    selected_n=bins(selected_truth)
    mae_p=np.divide(bins(abs_p),selected_n,out=np.zeros(n),where=selected_n>0)
    mae_g=np.divide(bins(abs_g),selected_n,out=np.zeros(n),where=selected_n>0)
    mean_p=np.divide(bins(bias_p),selected_n,out=np.zeros(n),where=selected_n>0)
    mean_g=np.divide(bins(bias_g),selected_n,out=np.zeros(n),where=selected_n>0)
    close_p=np.divide(bins(near_p),selected_n,out=np.zeros(n),where=selected_n>0)
    close_g=np.divide(bins(near_g),selected_n,out=np.zeros(n),where=selected_n>0)
    residual_counts_p=matrix(residual_p,n,350)
    residual_counts_g=matrix(residual_g,n,350)
    out=output_path(args.output_dir,args.period,args.merge_low_t)
    out.mkdir(parents=True,exist_ok=True)
    with (out/"truth_response.csv").open("w",newline="") as f:
        w=csv.writer(f)
        w.writerow(["truth_bin","reconstructed_bin","proton_count","photon_count"])
        for i in range(n):
            for j in range(nm):
                w.writerow([LABELS[i],MIGRATION_LABELS[j],pmat[i,j],gmat[i,j]])
    with (out/"truth_acceptance_diagnostic.csv").open("w",newline="") as f:
        w=csv.writer(f)
        w.writerow(["t_min_GeV2","t_max_GeV2","generated","reconstructed_proton",
                    "reconstructed_photon","A_proton","A_photon","A_photon_over_A_proton"])
        for i in range(n):
            w.writerow([EDGES[i],EDGES[i+1],gen[i],rp[i],rg[i],
                        rp[i]/gen[i] if gen[i] else "",rg[i]/gen[i] if gen[i] else "",
                        rg[i]/rp[i] if rp[i] else ""])
    with (out/"truth_resolution.csv").open("w",newline="") as f:
        w=csv.writer(f)
        w.writerow(["t_min_GeV2","t_max_GeV2","selected_matched",
                    "mean_error_proton_GeV2","mean_error_photon_GeV2",
                    "mean_absolute_error_proton_GeV2","mean_absolute_error_photon_GeV2",
                    "fraction_within_0p05_proton","fraction_within_0p05_photon"])
        for i in range(n):
            w.writerow([EDGES[i],EDGES[i+1],int(selected_n[i]),
                        mean_p[i],mean_g[i],mae_p[i],mae_g[i],close_p[i],close_g[i]])
    with (out/"truth_residual_hist.csv").open("w",newline="") as f:
        w=csv.writer(f)
        w.writerow(["residual_min_GeV2","residual_max_GeV2",
                    *[f"proton_truth_bin_{i+1}" for i in range(n)],
                    *[f"photon_truth_bin_{i+1}" for i in range(n)]])
        for j in range(350):
            w.writerow([residual_edges[j],residual_edges[j+1],
                        *residual_counts_p[:,j],*residual_counts_g[:,j]])
    fig,axes=plt.subplots(2,2,figsize=(13,10),constrained_layout=True)
    for ax,counts,title in ((axes[0,0],pmat,"Reconstructed proton |t|"),
                            (axes[0,1],gmat,"Reconstructed photon |t|")):
        fraction=np.divide(counts,gen[:,None],out=np.zeros_like(counts,dtype=float),
                           where=gen[:,None]>0)*100
        image=ax.imshow(fraction,origin="lower",aspect="auto",cmap="Blues",vmin=0)
        ax.set(xticks=np.arange(nm),yticks=np.arange(n),
               xticklabels=MIGRATION_LABELS,yticklabels=LABELS,
               xlabel="Reconstructed |t| bin (GeV²)",ylabel="Generated |t| bin (GeV²)",
               title=title+" / generated (%)")
        ax.tick_params(axis="x",rotation=35)
        fig.colorbar(image,ax=ax,label="Percent of generated truth bin")
    xx=np.arange(n)
    axes[1,0].plot(xx,np.divide(rp,gen,out=np.zeros(n),where=gen>0),"o-",label="Proton definition")
    axes[1,0].plot(xx,np.divide(rg,gen,out=np.zeros(n),where=gen>0),"s-",label="Photon definition")
    axes[1,0].set(xticks=xx,xticklabels=LABELS,xlabel="|t| bin (GeV²)",
                  ylabel="Selected reconstructed / generated",title="Acceptance diagnostic")
    axes[1,0].tick_params(axis="x",rotation=30)
    axes[1,0].legend()
    axes[1,0].grid(alpha=.2)
    axes[1,1].plot(xx,mae_p,"o-",label="Proton definition")
    axes[1,1].plot(xx,mae_g,"s-",label="Photon definition")
    axes[1,1].set(xticks=xx,xticklabels=LABELS,xlabel="Generated |t| bin (GeV²)",
                  ylabel="Mean |reconstructed − generated| (GeV²)",
                  title="Which reconstruction is closer to truth? Lower is better")
    axes[1,1].tick_params(axis="x",rotation=30)
    axes[1,1].legend()
    axes[1,1].grid(alpha=.2)
    fig.suptitle(f"{args.period}: paired MC truth comparison; {len(paths)} files")
    target=out/"truth_diagnostics.png"
    fig.savefig(target,dpi=160)
    plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(15,8.5),sharex=True,sharey=True,
                          constrained_layout=True)
    for i,ax in enumerate(axes.flat):
        if i>=n:
            ax.axis("off")
            continue
        for counts,color,label in ((residual_counts_p[i],"#2673a8","Proton reconstruction"),
                                   (residual_counts_g[i],"#d97706","Photon reconstruction")):
            total=counts.sum()
            if total:
                ax.stairs(counts/total,residual_edges,color=color,label=label,
                          linewidth=1.7)
        ax.axvline(0,color="black",linestyle="--",linewidth=.8)
        ax.set(xlim=(-.5,1.0),title=f"Truth |t| = {LABELS[i]} GeV²; N={int(selected_n[i]):,}",
               xlabel="Reconstructed |t| − generated |t| (GeV²)",
               ylabel="Fraction per 0.02 GeV²")
        ax.text(.97,.94,
                f"Mean absolute error: p={mae_p[i]:.3f}, γ={mae_g[i]:.3f} GeV²\n"
                f"Within 0.05 GeV²: p={close_p[i]:.0%}, γ={close_g[i]:.0%}",
                transform=ax.transAxes,ha="right",va="top",fontsize=8,
                bbox=dict(facecolor="white",alpha=.88,edgecolor="none"))
        ax.grid(alpha=.15)
    axes[0,0].legend(fontsize=8,loc="upper left")
    fig.suptitle(f"{args.period}: reconstructed proton versus photon |t| errors, same selected MC events")
    residual_target=out/"truth_reconstruction_errors.png"
    fig.savefig(residual_target,dpi=160)
    plt.close(fig)
    print(f"[TRUTH] files={len(paths)} selected={int(nselected.GetValue()):,}")
    print(f"[OUTPUT] {target}")
    print(f"[OUTPUT] {residual_target}")

if __name__=="__main__":
    main()
